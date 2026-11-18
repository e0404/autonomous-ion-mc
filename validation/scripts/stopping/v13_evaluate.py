# ruff: noqa: E501
"""Sealed evaluator of the V3-013 rows V-C1, V-C1-A, V-C2, V-C2-D and V-C3 (plan
``validation/plans/v3-013-acceptance.md`` sections A.6, D steps 3-6 and F).

Launched ONLY through the sealed launcher (it refuses any other launch)::

    .venv/bin/python -I -S -B -X pycache_prefix="$(mktemp -d)" \
        validation/scripts/stopping/sealed_launch.py validation/scripts/stopping/v13_evaluate.py \
        --seal-sha <sha> --code-sha <sha> --records <dir> --out <dir>

It never modifies the evidence registry and writes only ``<out>/v13_result.json``. It fails closed
(exit status 2, no result file) in this order:

(a) ``seal-review.json`` (and ``samples-review.json``) are codex-worker-signed review records,
    verified against ``infrastructure/keys/trusted_attestation_keys.json`` as committed at the seal
    SHA, passed, ``head_sha`` equal to the pinned SHA; the plan carries exactly one line
    ``seal review <id> at <sha>; manifest sha256 <digest>`` and one ``samples review <id> at <sha>``
    naming the same records, and the manifest digest equals the on-disk manifest's sha256;
(b) every path of ``v13-frozen-manifest.json`` (and the manifest, this evaluator and
    ``sealed_launch.py``): on-disk sha256 = manifest hash = sha256 of the blob at the seal SHA
    (``git cat-file -p``); the ``src/ionmc`` tree at the seal equals the manifest's pinned tree;
(c) the launch record of ``sealed_launch.py`` (the launcher and this file are the sealed bytes) and
    the environment digests recomputed with the sealed ``env_digest.py`` equal
    ``seal-environment.json``.

Every evidence run is then checked and marked invalid (never silently dropped) on any failure:
(d) sealed blobs and the ``src/ionmc`` tree at the run's ``code_sha`` identical to the seal's, the
    seal an ancestor of ``code_sha``; ionmc runs: signed host-run record, published outputs,
    step document bound to the run, identity block and table identities equal the pins of
    ``vc-predictions.json``, environment digests equal the seal, ``warp_modules`` identical across
    all runs of one backend, start after the seal review finished; TOPAS runs (signed
    reference-runner records): start after both reviews finished (the seal review and the
    ionmc-samples review), (e) executed ``input.txt`` / ``case.json`` equal the sealed bundle, the
    r2 run-validity checks (stderr, histories, UrbanMsc, frozen physics lines, modules).
Then the science: V-C1 per case (TD-EM (a)-(e), TOST, repeat rule), V-C1-A (report-only), V-C2
(paired delete-one-batch jackknife), V-C2-D (gates (i) and (ii)), V-C3 and the integration state.

Records directory (``--records``) layout, fixed here::

    index.json                                  {"ionmc": [{run_id, row, case, model, backend, k,
                                                 step_stdout}], "topas": [{bundle, run_id}]}
    host-records/<RUN-id>.json                  signed host-runner record of an ionmc run
    published/validation/generated/<RUN-id>/..  the run's published outputs (step stdout, npz);
                                                ``step_stdout`` is the path inside that directory
    topas/<REF-id>/{record.json, artifacts.json, request.json, stdout.txt, stderr.txt,
                    work/idd_dose.csv}          signed reference-runner record and outputs
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]


def _load(name: str, path: Path) -> ModuleType:
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# First statement after the stdlib imports, before numpy or ionmc is imported: capture the launch
# state (enforced by evaluate() before anything else).
env_seal = _load("env_seal", HERE / "env_seal.py")  # stdlib only at import
LAUNCH: dict[str, Any] = env_seal.capture_launch()

import numpy as np  # noqa: E402

rio = _load("r2_io", HERE.parent / "ion_em" / "r2_io.py")
st = _load("r2_statistics", HERE.parent / "ion_em" / "r2_statistics.py")
r2e = _load("r2_evaluate", HERE.parent / "ion_em" / "r2_evaluate.py")
envd = _load("env_digest", HERE / "env_digest.py")

RecordError = rio.RecordError
TASK_ID = "V3-013"
RESULT_SCHEMA = "ionmc-v13-result-1"
RESULT_NAME = "v13_result.json"
MANIFEST_SCHEMA = "ionmc-v13-frozen-manifest-1"
PREDICTIONS_SCHEMA = "ionmc-v13-vc-predictions-1"
RESULTS_DIR = "validation/results/stopping/v3-013"
SEAL_REVIEW = f"{RESULTS_DIR}/seal-review.json"
SAMPLES_REVIEW = f"{RESULTS_DIR}/samples-review.json"
ENV_SEAL = f"{RESULTS_DIR}/seal-environment.json"
PREDICTIONS = f"{RESULTS_DIR}/vc-predictions.json"
VC3 = f"{RESULTS_DIR}/v-c3.json"
MANIFEST = "validation/reference_cases/topas/v13-frozen-manifest.json"
BUNDLE_DIR = "validation/reference_cases/topas"
PLAN_PATH = "validation/plans/v3-013-acceptance.md"
TRUST_PATH = "infrastructure/keys/trusted_attestation_keys.json"
OWN_REL = "validation/scripts/stopping/v13_evaluate.py"
LAUNCHER_REL = env_seal.LAUNCHER_REL
ENV_DIGEST_REL = "validation/scripts/stopping/env_digest.py"
STEPS_REL = "validation/scripts/stopping/steps_v13.py"
#: paths the manifest must list, besides every TOPAS bundle file and the import closure
REQUIRED_SEALED = (
    OWN_REL, LAUNCHER_REL, ENV_DIGEST_REL, STEPS_REL, PREDICTIONS, TRUST_PATH, ENV_SEAL,
    "pyproject.toml", "uv.lock",
)  # fmt: skip
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
SEAL_NOTE_RE = re.compile(
    r"seal review (?P<id>REVIEW-[0-9a-f]{32}) at (?P<sha>[0-9a-f]{40}); "
    r"manifest sha256 (?P<digest>[0-9a-f]{64})(?![0-9A-Za-z])"
)
SAMPLES_NOTE_RE = re.compile(
    r"samples review (?P<id>REVIEW-[0-9a-f]{32}) at (?P<sha>[0-9a-f]{40})(?![0-9A-Za-z])"
)

CASES = (1, 2, 3)  # He 75 / 125 / 200 MeV/u
MODEL_ICRU90 = "icru90-water-v1"
MODEL_DEFAULT = "bethe-ionmc-v1"
MODELS = (MODEL_ICRU90, MODEL_DEFAULT)
BACKENDS = ("python", "warp-cpu")
ION13_BASE = 20612004
TOPAS13_BASE = 20610000
EVIDENCE_PARTS = (0, 1, 2, 3)
REPEAT_OFFSET = 4
N_HISTORIES = 1_000_000
N_BATCHES = 40
N_HISTORIES_D = 2000
N_BATCHES_D = 20
TOPAS_HISTORIES = 25000
TOL_FLOOR_MM = 0.02  # V-C2 tolerance max(0.02 mm, 0.10 |dR_pred|)
TOL_FRACTION = 0.10
ROUND_TRIP = 1e-4  # V-C2-D (i) and (ii)
Z_TOST = st.Z95
EVIDENCE_PASS = "topas-agreement-shared-lineage:V-C1:PASS"
EVIDENCE_FAIL = "failed:V-C1"


# ---------------------------------------------------------------------------------------------
# git (read-only)
# ---------------------------------------------------------------------------------------------


def _git_env() -> dict[str, str]:
    """The environment for git calls on an explicit repository: inherited ``GIT_*`` variables
    (the controlled host runner exports ``GIT_WORK_TREE=/workspace``) would redirect ``git -C``
    to another work tree, so they are dropped."""
    return {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}


def git(repo: Path | str, *args: str) -> bytes:
    proc = subprocess.run(  # noqa: S603 - fixed argv, no shell
        ["git", "-C", str(repo), *args], capture_output=True, check=False, env=_git_env()
    )
    if proc.returncode != 0:
        raise RecordError(f"git {' '.join(args)} failed: {proc.stderr.decode().strip()}")
    return proc.stdout


def blob_at(repo: Path | str, sha: str, rel: str) -> bytes:
    """The blob of ``rel`` at commit ``sha`` (``git cat-file -p``)."""
    return git(repo, "cat-file", "-p", f"{sha}:{rel}")


def blob_ids(repo: Path | str, sha: str, rels: Sequence[str]) -> dict[str, str]:
    """{path: blob id} of ``rels`` at commit ``sha`` (absent paths are absent)."""
    out: dict[str, str] = {}
    for ln in git(repo, "ls-tree", "-r", "--full-tree", sha, "--", *rels).decode().splitlines():
        meta, path = ln.split("\t", 1)
        out[path] = meta.split()[2]
    return out


def tree_id(repo: Path | str, sha: str, rel: str = "src/ionmc") -> str:
    return git(repo, "rev-parse", f"{sha}:{rel}").decode().strip()


def sha256_bytes(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def require_clean(repo: Path | str) -> None:
    if git(repo, "status", "--porcelain", "--untracked-files=no").strip():
        raise RecordError("the evaluator worktree is not clean (modified tracked files)")


def _time(value: Any, what: str) -> Any:
    return rio._parse_time(value, what)


def _read_json(path: Path, what: str) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RecordError(f"{what}: cannot read {path}: {exc}") from exc


# ---------------------------------------------------------------------------------------------
# (a) review records, plan pins, trust map
# ---------------------------------------------------------------------------------------------


def review_records_module() -> Any:
    rio.attest_module()  # puts the checkout root on sys.path
    return importlib.import_module("infrastructure.validation.review_records")


def host_records_module() -> Any:
    rio.attest_module()
    return importlib.import_module("infrastructure.validation.host_run_records")


def parse_trust(blob: bytes, where: str) -> dict[str, Any]:
    try:
        data = json.loads(blob.decode("utf-8"))
        return {s: {e["key_id"]: e["public_key"] for e in entries} for s, entries in data.items()}
    except (ValueError, KeyError, TypeError, AttributeError) as exc:
        raise RecordError(f"{where} is malformed: {exc}") from exc


def check_review(doc: Any, what: str, trust: Mapping[str, Any], head_sha: str) -> dict[str, Any]:
    """A codex-worker-signed review record, verified against ``trust`` (task V3-013), passed, and
    reviewing exactly ``head_sha``; anything else raises RecordError."""
    if not isinstance(doc, dict):
        raise RecordError(f"{what}: not a review record object")
    status, why = review_records_module().verify_review_record(
        doc, dict(trust), expected_task_id=TASK_ID
    )
    if status != "verified":
        raise RecordError(f"{what}: codex-worker review record {status}: {why}")
    report = doc.get("report")
    if (
        doc.get("status") != "passed"
        or not isinstance(report, dict)
        or report.get("verdict") != "pass"
    ):
        raise RecordError(f"{what}: the review did not pass")
    if doc.get("head_sha") != head_sha:
        raise RecordError(f"{what}: head_sha {doc.get('head_sha')} is not the pinned {head_sha}")
    return dict(doc)


def parse_pins(plan_text: str) -> tuple[dict[str, str], dict[str, str]]:
    """The one seal-pin note and the one samples-pin note of the plan (zero or several: fail)."""
    seal = list(SEAL_NOTE_RE.finditer(plan_text))
    samples = list(SAMPLES_NOTE_RE.finditer(plan_text))
    if len(seal) != 1:
        raise RecordError(
            f"{PLAN_PATH}: expected exactly one seal review pin note, found {len(seal)}"
        )
    if len(samples) != 1:
        raise RecordError(
            f"{PLAN_PATH}: expected exactly one samples review pin note, found {len(samples)}"
        )
    return seal[0].groupdict(), samples[0].groupdict()


# ---------------------------------------------------------------------------------------------
# (b) the sealed closure
# ---------------------------------------------------------------------------------------------


@dataclass
class Seal:
    repo: Path
    seal_sha: str
    head: str
    manifest_sha256: str
    files: dict[str, str]  # sealed path -> sha256 (manifest)
    tree: str
    blob_ids: dict[str, str]  # sealed path (and the manifest) -> blob id at the seal
    trust: dict[str, Any]
    trust_sha256: str
    seal_review: dict[str, Any]
    samples_review: dict[str, Any]
    notes: dict[str, dict[str, str]]
    predictions: dict[str, Any]
    environment: dict[str, Any]
    finished: Any = None  # seal review finished_at
    samples_finished: Any = None
    samples_started: Any = None


def load_manifest(repo: Path, seal_sha: str) -> tuple[dict[str, Any], str, dict[str, str]]:
    """(b): the manifest, its sha256 and {sealed path: sha256}; every sealed path (and the
    manifest, the evaluator and the launcher) byte-identical on disk, in the manifest and in the
    blob at the seal SHA; the ``src/ionmc`` tree equal to the pinned one."""
    disk = (repo / MANIFEST).read_bytes() if (repo / MANIFEST).is_file() else None
    if disk is None:
        raise RecordError(f"{MANIFEST} does not exist")
    if blob_at(repo, seal_sha, MANIFEST) != disk:
        raise RecordError(f"{MANIFEST} differs between the seal {seal_sha} and the working copy")
    try:
        doc = json.loads(disk.decode("utf-8"))
    except ValueError as exc:
        raise RecordError(f"{MANIFEST} is not JSON: {exc}") from exc
    files = doc.get("files") if isinstance(doc, dict) else None
    if (
        not isinstance(doc, dict)
        or doc.get("schema") != MANIFEST_SCHEMA
        or not isinstance(files, dict)
    ):
        raise RecordError(f"{MANIFEST}: schema is not {MANIFEST_SCHEMA} or files missing")
    if MANIFEST in files:
        raise RecordError("the manifest must not list itself")
    missing = [p for p in REQUIRED_SEALED if p not in files]
    if missing:
        raise RecordError(f"the manifest does not list required sealed paths {missing}")
    if not any(p.startswith(f"{BUNDLE_DIR}/") and p.endswith("/input.txt") for p in files):
        raise RecordError("the manifest lists no TOPAS bundle")
    for rel, want in sorted(files.items()):
        path = repo / rel
        if path.is_symlink() or not path.is_file():
            raise RecordError(f"sealed path {rel} is not a regular file")
        on_disk = sha256_bytes(path.read_bytes())
        at_seal = sha256_bytes(blob_at(repo, seal_sha, rel))
        if not (on_disk == want == at_seal):
            raise RecordError(
                f"sealed path {rel}: sha256 differs (disk {on_disk[:12]}, manifest {str(want)[:12]}, "
                f"seal blob {at_seal[:12]})"
            )
    tree = tree_id(repo, seal_sha)
    if doc.get("src_ionmc_tree") != tree:
        raise RecordError(
            f"src/ionmc tree at the seal {tree} differs from the manifest's {doc.get('src_ionmc_tree')}"
        )
    return doc, sha256_bytes(disk), dict(files)


def sealed_problems(seal: Seal, sha: str, what: str) -> list[str]:
    """(d): problems of commit ``sha`` against the seal: the seal an ancestor of ``sha``, the blob
    of every sealed path and the manifest, and the ``src/ionmc`` tree identical to the seal's."""
    out: list[str] = []
    try:
        if not SHA_RE.match(sha):
            return [f"{what}: {sha!r} is not a commit SHA"]
        if not rio.is_ancestor(seal.repo, seal.seal_sha, sha):
            out.append(f"{what}: the seal {seal.seal_sha[:12]} is not an ancestor of {sha[:12]}")
        ids = blob_ids(seal.repo, sha, sorted(seal.blob_ids))
        bad = sorted(p for p, b in seal.blob_ids.items() if ids.get(p) != b)
        if bad:
            out.append(
                f"{what}: sealed blobs differ at {sha[:12]}: {bad[:5]}{' ...' if len(bad) > 5 else ''}"
            )
        if tree_id(seal.repo, sha) != seal.tree:
            out.append(f"{what}: the src/ionmc tree at {sha[:12]} differs from the seal's")
    except RecordError as exc:
        out.append(f"{what}: {exc}")
    return out


def load_seal(
    repo: Path,
    seal_sha: str,
    environment: Callable[[], dict[str, Any]] | None = None,
    launch_check: bool = True,
) -> Seal:
    """Steps (a)-(c); any failure raises RecordError (evaluation stops)."""
    if not SHA_RE.match(seal_sha):
        raise RecordError(f"--seal-sha {seal_sha!r} is not a 40-hex commit SHA")
    head = git(repo, "rev-parse", "HEAD").decode().strip()
    if not rio.is_ancestor(repo, seal_sha, head):
        raise RecordError(f"the seal {seal_sha} is not an ancestor of HEAD {head}")
    require_clean(repo)
    trust_blob = blob_at(repo, seal_sha, TRUST_PATH)
    trust = parse_trust(trust_blob, f"{TRUST_PATH} at the seal")
    manifest, manifest_sha, files = load_manifest(repo, seal_sha)
    # (a) review records and the plan pins
    plan = (repo / PLAN_PATH).read_text(encoding="utf-8")
    if blob_at(repo, head, PLAN_PATH).decode("utf-8") != plan:
        raise RecordError(f"{PLAN_PATH}: the working copy differs from HEAD")
    seal_note, samples_note = parse_pins(plan)
    if seal_note["sha"] != seal_sha:
        raise RecordError(f"the plan pins the seal {seal_note['sha']}, not --seal-sha {seal_sha}")
    if seal_note["digest"] != manifest_sha:
        raise RecordError("the plan pin note's manifest sha256 differs from the on-disk manifest")
    seal_review = check_review(
        _read_json(repo / SEAL_REVIEW, "seal review record") if (repo / SEAL_REVIEW).is_file() else None,
        SEAL_REVIEW, trust, seal_sha,
    )  # fmt: skip
    if seal_review["review_id"] != seal_note["id"]:
        raise RecordError(f"{SEAL_REVIEW} is not the review {seal_note['id']} named by the plan")
    samples_sha = samples_note["sha"]
    if not rio.is_ancestor(repo, seal_sha, samples_sha) or not rio.is_ancestor(
        repo, samples_sha, head
    ):
        raise RecordError("the samples commit is not between the seal and HEAD")
    samples_review = check_review(
        _read_json(repo / SAMPLES_REVIEW, "samples review record") if (repo / SAMPLES_REVIEW).is_file() else None,
        SAMPLES_REVIEW, trust, samples_sha,
    )  # fmt: skip
    if samples_review["review_id"] != samples_note["id"]:
        raise RecordError(
            f"{SAMPLES_REVIEW} is not the review {samples_note['id']} named by the plan"
        )
    # (c) launch record and environment
    if launch_check:
        try:
            launch = env_seal.check_launch(LAUNCH, Path(__file__).resolve())
            env_seal.numpy_modules_check()
        except env_seal.EnvironmentSealError as exc:
            raise RecordError(f"launch: {exc}") from exc
        if launch["launcher_record"].get("launcher_sha256") != files[LAUNCHER_REL]:
            raise RecordError(
                "the launcher that started this run is not the sealed sealed_launch.py"
            )
        if launch["launcher_record"].get("target_sha256") != files[OWN_REL]:
            raise RecordError("the launched target differs from the sealed evaluator")
    seal_env = _read_json(repo / ENV_SEAL, "seal-environment.json")
    now = (environment or envd.environment_digests)()
    diff = envd.differences(now, seal_env) if isinstance(seal_env, dict) else ["unreadable"]
    if diff:
        raise RecordError(f"the environment differs from {ENV_SEAL} in {diff}")
    predictions = _read_json(repo / PREDICTIONS, "vc-predictions.json")
    if predictions.get("schema") != PREDICTIONS_SCHEMA:
        raise RecordError(f"{PREDICTIONS}: schema is not {PREDICTIONS_SCHEMA}")
    if predictions.get("src_ionmc_tree") != manifest["src_ionmc_tree"]:
        raise RecordError(f"{PREDICTIONS}: src_ionmc_tree differs from the manifest's")
    if predictions.get("script", {}).get("sha256") != files[STEPS_REL]:
        raise RecordError(
            f"{PREDICTIONS}: the pinned steps_v13.py hash differs from the manifest's"
        )
    ids = blob_ids(repo, seal_sha, [*files, MANIFEST])
    if len(ids) != len(files) + 1:
        raise RecordError("some sealed path has no blob at the seal")
    return Seal(
        repo=repo, seal_sha=seal_sha, head=head, manifest_sha256=manifest_sha, files=files,
        tree=manifest["src_ionmc_tree"], blob_ids=ids, trust=trust,
        trust_sha256=sha256_bytes(trust_blob), seal_review=seal_review,
        samples_review=samples_review, notes={"seal": seal_note, "samples": samples_note},
        predictions=predictions, environment=seal_env,
        finished=_time(seal_review["finished_at"], "seal review finished_at"),
        samples_finished=_time(samples_review["finished_at"], "samples review finished_at"),
        samples_started=_time(samples_review["started_at"], "samples review started_at"),
    )  # fmt: skip


# ---------------------------------------------------------------------------------------------
# estimators (pure functions; TD-EM (a)-(e) as r2_evaluate.td_statistics, EM branch)
# ---------------------------------------------------------------------------------------------


def td_em(
    ion_batches: np.ndarray,
    topas_parts: list[tuple[np.ndarray, int]],
    hz: float,
    r_plan_mm: float,
) -> dict[str, Any]:
    """TD-EM (r2-topas-target.md section 5.1): every criterion (a)-(e) gating, (b) only where
    W_T >= 1.0 mm. Same formulas and helpers as ``r2_evaluate.td_statistics`` (which cannot be called
    here: it looks the case up in the V3-008 table); ``ion_batches`` (B, nz) energy per primary per
    batch [MeV], ``topas_parts`` K = 4 (run-total energy per bin [MeV], histories)."""
    et = sum(p[0] for p in topas_parts) / float(sum(p[1] for p in topas_parts)) / hz
    et_p = [p[0] / float(p[1]) / hz for p in topas_parts]
    ei_b = [np.asarray(b) / hz for b in ion_batches]
    ei = np.mean(np.asarray(ion_batches), axis=0) / hz
    qi, qt = r2e._quantities(ei, hz, r_plan_mm), r2e._quantities(et, hz, r_plan_mm)
    qib = [r2e._quantities(b, hz, r_plan_mm) for b in ei_b]
    qtp = [r2e._quantities(p, hz, r_plan_mm) for p in et_p]
    crit = [r2e._criterion("(a)", "r80", qi, qib, qt, qtp, r2e._r80_tol)]
    w_t = qt.get("w")
    crit.append(
        r2e._criterion(
            "(b)", "w", qi, qib, qt, qtp, lambda r: 0.10 * r, gated=(w_t is None or w_t >= 1.0)
        )
    )
    crit.append(r2e._criterion("(c)", "z_peak", qi, qib, qt, qtp, lambda r: 0.5))
    crit.append(r2e._criterion("(d)", "pe", qi, qib, qt, qtp, lambda r: 0.03 * r))
    r80_t = qt.get("r80")
    nan = math.nan
    if r80_t is None:
        res = st.tost(nan, nan, nan, nan, nan)
        res.update(
            criterion="(e)", observable="plateau", gated=True, note="invalid: R80_T undefined"
        )
        crit.append(res)
    else:
        groups = [
            r2e._safe(lambda e=e: st.plateau_groups(e, hz, r_plan_mm, r80_t))
            for e in [ei, et, *ei_b, *et_p]
        ]
        if any(g[0] is None for g in groups):
            res = st.tost(nan, nan, nan, nan, nan)
            res.update(
                criterion="(e)", observable="plateau", gated=True,
                note="invalid: " + str(groups[0][1] or groups[1][1]),
            )  # fmt: skip
            crit.append(res)
        else:
            gi, gt = groups[0][0], groups[1][0]
            gib = np.array([g[0] for g in groups[2 : 2 + len(ei_b)]])
            gtp = np.array([g[0] for g in groups[2 + len(ei_b) :]])
            per = [
                st.tost(
                    float(gi[j]),
                    float(gt[j]),
                    st.batch_se(gib[:, j]),
                    st.topas_se(gtp[:, j]),
                    0.02 * float(gt[j]),
                )
                for j in range(gi.size)
            ]
            worst = max(per, key=lambda p: (p["statistic"] or math.inf) - p["tol"])
            cls = st.range_class(r_plan_mm)
            crit.append(
                {
                    "criterion": "(e)", "observable": "plateau", "gated": True, "groups": len(per),
                    "group_mm": cls["group_mm"], "start_mm": cls["ent_lo"], "end_mm": 0.8 * r80_t,
                    "failing_groups": [j for j, p in enumerate(per) if p["verdict"] != "PASS"],
                    "worst_group": worst, "per_group": per,
                    "verdict": "PASS" if all(p["verdict"] == "PASS" for p in per) else "FAIL",
                }
            )  # fmt: skip
    report = {
        "r90": r2e._criterion("R90", "r90", qi, qib, qt, qtp, lambda r: math.inf, gated=False)
    }
    for c in crit:
        if c.get("gated") is False:
            c["verdict"] = "REPORT-ONLY"
    gating = [c for c in crit if c["verdict"] != "REPORT-ONLY"]
    status = "PASS" if gating and all(c["verdict"] == "PASS" for c in gating) else "FAIL"
    pooled = {
        k: {"ionmc": qi.get(k), "topas": qt.get(k)}
        for k in ("z_peak", "d_peak", "r80", "r20", "w", "d_ent", "pe")
    }
    return {"status": status, "criteria": crit, "pooled": pooled, "report_only": report}


def jackknife_delta_r80(
    a: np.ndarray, b: np.ndarray, hz: float, r_plan_mm: float
) -> dict[str, Any]:
    """dR80 = R80(a) - R80(b) of the pooled curves of two paired (B, nz) batch arrays and its
    delete-one-batch jackknife sigma over the B paired batches."""
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.shape != b.shape or a.ndim != 2 or a.shape[0] < 2:
        raise RecordError("paired batch arrays must share a (B >= 2, nz) shape")

    def r80(x: np.ndarray) -> float | None:
        return r2e._quantities(np.mean(x, axis=0) / hz, hz, r_plan_mm)["r80"]

    full_a, full_b = r80(a), r80(b)
    if full_a is None or full_b is None:
        return {"delta_mm": None, "sigma_mm": None, "note": "R80 undefined"}
    n = a.shape[0]
    deltas = []
    for i in range(n):
        keep = np.arange(n) != i
        ra, rb = r80(a[keep]), r80(b[keep])
        if ra is None or rb is None:
            return {
                "delta_mm": None,
                "sigma_mm": None,
                "note": "R80 undefined in a jackknife replicate",
            }
        deltas.append(ra - rb)
    d = np.asarray(deltas)
    sigma = math.sqrt((n - 1) / n * float(np.sum((d - d.mean()) ** 2)))
    return {"delta_mm": full_a - full_b, "sigma_mm": sigma, "r80_a_mm": full_a, "r80_b_mm": full_b,
            "n_batches": n, "note": None}  # fmt: skip


def vc2_verdict(delta: float | None, sigma: float | None, pred: float) -> dict[str, Any]:
    """PASS iff |dR80 - dR_pred| + 1.645 sigma <= max(0.02 mm, 0.10 |dR_pred|)."""
    tol = max(TOL_FLOOR_MM, TOL_FRACTION * abs(pred))
    out: dict[str, Any] = {
        "delta_r80_mm": delta,
        "sigma_mm": sigma,
        "delta_r_pred_mm": pred,
        "tol_mm": tol,
    }
    if delta is None or sigma is None or not (math.isfinite(delta) and math.isfinite(sigma)):
        return out | {"statistic_mm": None, "verdict": "FAIL"}
    stat = abs(delta - pred) + Z_TOST * sigma
    return out | {"statistic_mm": stat, "verdict": "PASS" if stat <= tol else "FAIL"}


def vc2d_verdict(
    z_icru: float, r_icru: float, z_def: float, r_def: float, pred_d: float,
    tol_icru: float, tol_def: float, tol_diff: float,
) -> dict[str, Any]:  # fmt: skip
    """(i) |Z_m - R_m(E0)| <= 1e-4 R_m(E0) for both models, (ii) |(Z_icru - Z_def) - dR_pred,D| <=
    1e-4 R_CSDA,default(E0); the tolerances are the sealed ones of vc-predictions.json."""
    e_i, e_d = abs(z_icru - r_icru), abs(z_def - r_def)
    e_diff = abs((z_icru - z_def) - pred_d)
    ok_i = e_i <= tol_icru and e_d <= tol_def
    ok_ii = e_diff <= tol_diff
    return {
        "round_trip_icru90_mm": e_i, "round_trip_default_mm": e_d,
        "round_trip_tol_icru90_mm": tol_icru, "round_trip_tol_default_mm": tol_def,
        "difference_error_mm": e_diff, "difference_tol_mm": tol_diff,
        "gate_i": "PASS" if ok_i else "FAIL", "gate_ii": "PASS" if ok_ii else "FAIL",
        "verdict": "PASS" if ok_i and ok_ii else "FAIL",
    }  # fmt: skip


def warp_modules_problems(runs: Sequence[Mapping[str, Any]]) -> dict[str, str]:
    """{run_id: problem} for every ionmc run of a backend whose ``warp_modules`` are not identical
    to those of all other runs of that backend (python: all None; warp-cpu: identical non-empty
    maps). When the runs of a backend disagree every one of them is marked (the reference cannot
    be known)."""
    out: dict[str, str] = {}
    for backend in sorted({r["backend"] for r in runs}):
        group = [r for r in runs if r["backend"] == backend]
        mods = [r.get("warp_modules") for r in group]
        if backend == "python":
            bad = [r for r in group if r.get("warp_modules") is not None]
            note = "python run records warp_modules"
        else:
            same = all(isinstance(m, dict) and m and m == mods[0] for m in mods)
            bad = [] if same else group
            note = "warp_modules are empty or not identical across the runs of the backend"
        out.update({r["run_id"]: note for r in bad})
    return out


# ---------------------------------------------------------------------------------------------
# evidence runs
# ---------------------------------------------------------------------------------------------


def case_info(pred: Mapping[str, Any], case: int) -> dict[str, Any]:
    for e in pred["energies"]:
        if e["case"] == case:
            return e
    raise RecordError(f"vc-predictions.json has no case {case}")


def bundle_name(energy: float, part: int) -> str:
    return f"he-water-{energy:g}mevu-v13-idd-em-part{part}"


def expected_ionmc_runs() -> list[tuple[str, int, str, str]]:
    runs = [("vc-sample", c, m, "warp-cpu") for c in CASES for m in MODELS]
    runs += [("vc2d-sample", c, m, b) for c in CASES for m in MODELS for b in BACKENDS]
    return runs


def step_expectation(row: str, case: int, backend: str) -> dict[str, Any]:
    common = {"case": case, "ion": "he", "k": 0, "backend": backend, "precision": "float64",
              "valid": True, "pass": True}  # fmt: skip
    if row == "vc-sample":
        return common | {
            "row": row, "seed": ION13_BASE + 100 * (case - 1), "r_index": case - 1,
            "histories": N_HISTORIES, "frozen_histories": N_HISTORIES, "reduced": False,
            "n_batches": N_BATCHES,
        }  # fmt: skip
    return common | {
        "row": row, "seed": ION13_BASE + 1500, "r_index": 15, "histories": N_HISTORIES_D,
        "n_batches": N_BATCHES_D, "end_codes_all_cutoff": True, "e_end_all_at_or_below_cut": True,
    }  # fmt: skip


def identity_problems(
    doc: Mapping[str, Any], pred: Mapping[str, Any], model: str, case: int
) -> list[str]:
    """The run's identity block and table identities equal the seal pins exactly."""
    eff = doc.get("effective_config")
    if not isinstance(eff, dict) or not isinstance(eff.get("stopping"), dict):
        return ["the step document has no effective_config.stopping identity block"]
    out = []
    stop, tables = eff["stopping"], eff.get("tables") if isinstance(eff.get("tables"), dict) else {}
    try:
        want = pred["identity_blocks"][model][str(case)]
        pins = pred["models"][model]
    except KeyError:
        return [f"vc-predictions.json has no pins for {model} case {case}"]
    if stop != want:
        out.append("the identity block differs from the sealed identity block")
    if doc.get("stopping_identity") != stop:
        out.append("stopping_identity differs from effective_config.stopping")
    for key in ("model_id_effective", "version", "spec_sha256", "content_sha256"):
        if stop.get(key) != pins.get(key):
            out.append(f"identity {key} differs from the seal pin")
    if tables.get("sha256") != pins["tables_sha256"].get(str(case)):
        out.append("the transport tables sha256 differs from the seal pin")
    if tables.get("materials") != pins["tables_materials"].get(str(case)):
        out.append("the transport tables materials differ from the seal pins")
    return out


@dataclass
class IonmcRun:
    run_id: str
    row: str
    case: int
    model: str
    backend: str
    problems: list[str] = field(default_factory=list)
    doc: dict[str, Any] = field(default_factory=dict)
    code_sha: str = ""
    batches: np.ndarray | None = None
    edges_z: np.ndarray | None = None

    @property
    def valid(self) -> bool:
        return not self.problems

    def summary(self) -> dict[str, Any]:
        return {"run_id": self.run_id, "row": self.row, "case": self.case, "model": self.model,
                "backend": self.backend, "code_sha": self.code_sha, "valid": self.valid,
                "problems": list(self.problems)}  # fmt: skip


def check_ionmc_run(seal: Seal, records: Path, entry: Mapping[str, Any]) -> IonmcRun:
    """One ionmc evidence run; every failure is a problem of the run (it is then invalid)."""
    run = IonmcRun(str(entry.get("run_id")), str(entry.get("row")), int(entry["case"]),
                   str(entry.get("model")), str(entry.get("backend")))  # fmt: skip
    try:
        _check_ionmc(seal, records, entry, run)
    except (RecordError, OSError, ValueError, KeyError, TypeError) as exc:
        run.problems.append(f"{type(exc).__name__}: {exc}")
    return run


def _check_ionmc(seal: Seal, records: Path, entry: Mapping[str, Any], run: IonmcRun) -> None:
    rid = run.run_id
    if not re.fullmatch(r"RUN-[0-9A-Za-z-]+", rid):
        raise RecordError(f"run id {rid!r} is malformed")
    rec = _read_json(records / "host-records" / f"{rid}.json", f"host-run record {rid}")
    sha = rec.get("sha")
    if not isinstance(sha, str) or not SHA_RE.match(sha):
        raise RecordError("the host-run record has no commit sha")
    status, why = host_records_module().verify_run_record(rec, TASK_ID, sha, rid, seal.trust)
    if status != "verified":
        raise RecordError(f"host-run record {status}: {why}")
    rio._bind(rid, "signed", rec, {
        "succeeded": True, "exit_code": 0, "timed_out": False, "dirty_before": False,
        "dirty_after": False, "committed_state_unchanged": True, "output_publication_error": None,
    })  # fmt: skip
    run.code_sha = sha
    run.problems += sealed_problems(seal, sha, "code_sha")
    started = _time(rec["started_at"], "host record started_at")
    if not started > seal.finished:
        run.problems.append("the run started before the seal review finished")
    published = records / "published"
    digest, files = rio.published_manifest(published, rid)
    if digest != rec.get("output_manifest_sha256"):
        raise RecordError(
            "the published outputs do not reproduce the signed output_manifest_sha256"
        )
    rel = str(entry["step_stdout"])
    base = published / "validation" / "generated" / rid
    key = f"validation/generated/{Path(rel).as_posix()}"
    blob = (base / rel).read_bytes()
    if files.get(key) != sha256_bytes(blob):
        raise RecordError("the step stdout is not a signed published output")
    doc = rio.step_document(blob.decode("utf-8"), rid)
    run.doc = doc
    rio._bind(rid, "step document", doc, {**step_expectation(run.row, run.case, run.backend),
                                          "model": run.model})  # fmt: skip
    counters = doc.get("counters")
    if (
        not isinstance(counters, dict)
        or not counters
        or not all(rio._exact(v, 0) for v in counters.values())
    ):
        run.problems.append("counters are not all zero")
    run.problems += identity_problems(doc, seal.predictions, run.model, run.case)
    env = doc.get("environment_digests")
    diff = envd.differences(env, seal.environment) if isinstance(env, dict) else ["missing"]
    if diff:
        run.problems.append(f"environment digests differ from the seal in {diff}")
    if run.row == "vc-sample":
        npz_rel = str(Path(rel).parent / str(doc["npz"]))
        npz_key = f"validation/generated/{Path(npz_rel).as_posix()}"
        npz_path = base / npz_rel
        if (
            files.get(npz_key) != doc.get("npz_sha256")
            or sha256_bytes(npz_path.read_bytes()) != doc["npz_sha256"]
        ):
            raise RecordError("the npz is not a signed published output with the recorded sha256")
        with np.load(npz_path, allow_pickle=False) as z:
            arr = np.asarray(z["idd_batch_energy_mev"], dtype=np.float64)
            edges = np.asarray(z["idd_edges2_mm"], dtype=np.float64)
        if (
            arr.ndim != 4
            or arr.shape[0] != N_BATCHES
            or arr.shape[1:3] != (1, 1)
            or not np.all(np.isfinite(arr))
        ):
            raise RecordError("idd_batch_energy_mev is not a finite (40, 1, 1, nz) array")
        if edges.shape != (arr.shape[3] + 1,):
            raise RecordError("idd_edges2_mm does not match the depth grid")
        run.batches, run.edges_z = arr[:, 0, 0, :], edges
    else:
        for k in ("z_m_mm", "r_e0_mm"):
            if not isinstance(doc.get(k), float) or not math.isfinite(doc[k]):
                raise RecordError(f"step document {k} is not a finite number")


@dataclass
class TopasRun:
    bundle: str
    part: int
    case: int
    run_id: str
    problems: list[str] = field(default_factory=list)
    energy: np.ndarray | None = None
    hz: float = 0.0
    histories: int = 0
    run: Any = None
    code_sha: str = ""
    started: str = ""

    @property
    def valid(self) -> bool:
        return not self.problems

    def summary(self) -> dict[str, Any]:
        return {"bundle": self.bundle, "part": self.part, "case": self.case, "run_id": self.run_id,
                "code_sha": self.code_sha, "reference_started": self.started,
                "valid": self.valid, "problems": list(self.problems)}  # fmt: skip


def sealed_bundle(seal: Seal, bundle: str) -> dict[str, dict[str, Any]]:
    out = {}
    for name in ("case.json", "input.txt"):
        rel = f"{BUNDLE_DIR}/{bundle}/{name}"
        if rel not in seal.files:
            raise RecordError(f"bundle {bundle} is not sealed ({rel})")
        blob = blob_at(seal.repo, seal.seal_sha, rel)
        out[name] = {"sha256": sha256_bytes(blob), "bytes": len(blob)}
    return out


def check_topas_run(
    seal: Seal, records: Path, entry: Mapping[str, Any], case: int, part: int
) -> TopasRun:
    tr = TopasRun(str(entry["bundle"]), part, case, str(entry["run_id"]))
    try:
        _check_topas(seal, records, tr)
    except (RecordError, OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        tr.problems.append(f"{type(exc).__name__}: {exc}")
    return tr


def _check_topas(seal: Seal, records: Path, tr: TopasRun) -> None:
    rid, bundle = tr.run_id, tr.bundle
    rdir = records / "topas" / rid
    result = _read_json(rdir / "record.json", f"signed run record {rid}")
    status, why = rio.attest_module().classify_record(result, seal.trust, rio.REFERENCE_SIGNER)
    if status != "verified":
        raise RecordError(f"record.json attestation {status}: {why}")
    for key in rio.RESULT_FIELDS:
        rio.require(result, key, f"result record {rid}")
    if (
        result["schema_version"] != 3
        or result.get("timed_out") is not False
        or result.get("task_id") != TASK_ID
    ):
        raise RecordError("schema_version / timed_out / task_id are not 3 / false / V3-013")
    if (
        not re.fullmatch(
            r"REF-" + re.escape(str(result["configuration_id"])[:20]) + r"-[0-9a-f]{8}", rid
        )
        or result["run_id"] != rid
    ):
        raise RecordError("run id is not bound to the record / configuration")
    if result["engine"] != "topas" or not (
        result["succeeded"] is True and result["exit_code"] == 0 and result["missing_outputs"] == []
    ):
        raise RecordError("not a successful topas run")
    code_sha = str(result["code_sha"])
    tr.code_sha, tr.started = code_sha, str(result["started_at"])
    tr.problems += sealed_problems(seal, code_sha, "code_sha")
    started, finished = (
        _time(result["started_at"], "started_at"),
        _time(result["finished_at"], "finished_at"),
    )
    if started > finished:
        raise RecordError("started_at is later than finished_at")
    for what, t in (("seal review", seal.finished), ("samples review", seal.samples_finished)):
        if not started > t:
            tr.problems.append(f"reference_started is not after the {what} finished_at")
    sizes: dict[str, Any] = {}
    listing = rio._listing(
        _read_json(rdir / "artifacts.json", "artifact listing"), rid, code_sha, sizes
    )
    work = {
        p[5:]: {"sha256": d, "bytes": sizes.get(p)}
        for p, d in listing.items()
        if p.startswith("work/")
    }
    if rio.canonical_digest(work) != result["output_files_sha256"] or result["output_count"] != len(
        work
    ):
        raise RecordError("listed work/ outputs do not reproduce the signed digest")
    for name in ("request.json", "stdout.txt", "stderr.txt", "record.json"):
        if name not in listing:
            raise RecordError(f"artifact listing lacks {name}")
    if sha256_bytes((rdir / "record.json").read_bytes()) != listing["record.json"]:
        raise RecordError("record.json differs from the listed record.json")
    run = rio.TopasRun(
        run_id=rid, case=tr.case, config="idd-em", part=tr.part, bundle=bundle,
        configuration_id=str(result["configuration_id"]), code_sha=code_sha,
        finished_at=str(result["finished_at"]), seconds=float(result["seconds"]),
        reference_started=tr.started, histories=0, expected_outputs=[], listing=listing,
        outputs_dir=rdir,
    )  # fmt: skip
    for name, key in (("stdout.txt", "stdout_sha256"), ("stderr.txt", "stderr_sha256")):
        rio.read_verified(run, name)
        if result[key] != listing[name]:
            raise RecordError(f"signed {key} differs from the listing")
    request_blob = rio.read_verified(run, "request.json")
    if sha256_bytes(request_blob) != result["request_sha256"]:
        raise RecordError("request.json does not hash to the signed request_sha256")
    request = json.loads(request_blob.decode("utf-8"))
    if rio.canonical_digest(request) != run.configuration_id or request.get("code_sha") != code_sha:
        raise RecordError("request.json is not bound to the record")
    sealed = sealed_bundle(seal, bundle)
    if request.get("input_files") != sealed:  # (e) executed input hashes = sealed bundle hashes
        raise RecordError(
            f"the executed input.txt/case.json of {bundle} differ from the sealed bundle"
        )
    case_doc = json.loads(blob_at(seal.repo, seal.seal_sha, f"{BUNDLE_DIR}/{bundle}/case.json"))
    if request.get("case") != case_doc:
        raise RecordError("request case differs from the sealed case.json")
    if case_doc.get("seeds") != [TOPAS13_BASE + 100 * tr.part + tr.case]:
        raise RecordError(f"seed {case_doc.get('seeds')} breaks the topas13 formula")
    hist = case_doc.get("histories")
    if hist != TOPAS_HISTORIES:
        raise RecordError(f"histories {hist!r} is not {TOPAS_HISTORIES}")
    run.histories, run.case_doc, run.request = hist, case_doc, request
    run.expected_outputs = list(case_doc["expected_outputs"])
    run.native = blob_at(seal.repo, seal.seal_sha, f"{BUNDLE_DIR}/{bundle}/input.txt").decode(
        "utf-8", "replace"
    )
    tr.problems += rio.run_checks(run)["problems"]
    e, hz = rio.topas_idd(run)
    tr.energy, tr.hz, tr.histories, tr.run = e, float(hz), int(hist), run


def load_index(
    records: Path, pred: Mapping[str, Any]
) -> tuple[list[dict], dict[tuple[int, int], dict]]:
    """The run index: exactly the expected ionmc runs, TOPAS evidence/repeat parts at most once."""
    idx = _read_json(records / "index.json", "run index")
    if (
        not isinstance(idx, dict)
        or not isinstance(idx.get("ionmc"), list)
        or not isinstance(idx.get("topas"), list)
    ):
        raise RecordError("index.json must hold the lists 'ionmc' and 'topas'")
    ion: dict[tuple, dict] = {}
    for e in idx["ionmc"]:
        key = (e.get("row"), e.get("case"), e.get("model"), e.get("backend"))
        if key in ion or key not in expected_ionmc_runs() or e.get("k") != 0:
            raise RecordError(f"unexpected or duplicate ionmc run {key}")
        ion[key] = e
    missing = [k for k in expected_ionmc_runs() if k not in ion]
    if missing:
        raise RecordError(f"the index lacks ionmc runs {missing}")
    names = {
        bundle_name(case_info(pred, c)["e_mev_per_u"], p): (c, p)
        for c in CASES
        for p in (*EVIDENCE_PARTS, *(q + REPEAT_OFFSET for q in EVIDENCE_PARTS))
    }
    top: dict[tuple[int, int], dict] = {}
    for e in idx["topas"]:
        if e.get("bundle") not in names:
            raise RecordError(
                f"index lists {e.get('bundle')!r}, which is not a V-C1 evidence or repeat bundle"
            )
        if names[e["bundle"]] in top:
            raise RecordError(f"index lists bundle {e['bundle']} twice")
        top[names[e["bundle"]]] = e
    return list(ion.values()), top


def select_parts(
    case: int, runs: Mapping[int, TopasRun]
) -> tuple[list[TopasRun] | None, list[str]]:
    """Repeat rule: an invalid evidence part is replaced once by its frozen repeat (part + 4); a
    second failure, a missing repeat or a repeat without trigger is a FAIL of the case (None)."""
    used: list[TopasRun] = []
    notes: list[str] = []
    for p in EVIDENCE_PARTS:
        orig, rep = runs.get(p), runs.get(p + REPEAT_OFFSET)
        if orig is None:
            return None, [f"evidence part {p} has no run"]
        if orig.valid:
            if rep is not None:
                return None, [
                    f"repeat part {p + REPEAT_OFFSET} was run without a validity failure of part {p}"
                ]
            used.append(orig)
        elif rep is None:
            return None, [
                f"part {p} is invalid ({orig.problems[:1]}) and its repeat part {p + REPEAT_OFFSET} is missing"
            ]
        elif rep.valid:
            used.append(rep)
            notes.append(f"part {p} invalid; repeat part {p + REPEAT_OFFSET} used")
        else:
            return None, [f"part {p} and its repeat part {p + REPEAT_OFFSET} are both invalid"]
    return used, notes


# ---------------------------------------------------------------------------------------------
# the science
# ---------------------------------------------------------------------------------------------


def _safe_json(v: Any) -> Any:
    if isinstance(v, dict):
        return {str(k): _safe_json(x) for k, x in v.items()}
    if isinstance(v, list | tuple):
        return [_safe_json(x) for x in v]
    if isinstance(v, np.ndarray):
        return _safe_json(v.tolist())
    if isinstance(v, np.generic):
        return _safe_json(v.item())
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v


def _ion(runs: Mapping[tuple, IonmcRun], row: str, case: int, model: str, backend: str) -> IonmcRun:
    return runs[(row, case, model, backend)]


def vc1_case(
    case: int, info: Mapping[str, Any], icru: IonmcRun, topas: Mapping[int, TopasRun], what: str
) -> dict[str, Any]:
    """V-C1 (what = 'V-C1') or V-C1-A (the default's samples) of one case."""
    out: dict[str, Any] = {"case": case, "model": icru.model, "row": what, "reasons": []}
    parts, notes = select_parts(case, topas)
    out["repeat_notes"] = notes
    if parts is None:
        return out | {"status": "FAIL", "reasons": notes, "td": None}
    out["topas_parts"] = [p.part for p in parts]
    reasons = []
    if not icru.valid or icru.batches is None:
        reasons.append(f"ionmc run {icru.run_id} is invalid: {icru.problems[:2]}")
    gate = rio.replicate_gate([p.run for p in parts])
    if gate:
        reasons.append(f"replicate gate: {gate}")
    hzs = {p.hz for p in parts}
    nzs = {p.energy.size for p in parts}  # type: ignore[union-attr]
    if len(hzs) != 1 or len(nzs) != 1:
        reasons.append("TOPAS parts differ in their depth grid")
    if reasons:
        return out | {"status": "FAIL", "reasons": reasons, "td": None}
    hz, nz = parts[0].hz, parts[0].energy.size  # type: ignore[union-attr]
    edges = np.arange(nz + 1) * hz
    try:
        rio.same_edges(icru.edges_z, edges, f"case {case}")  # type: ignore[arg-type]
        td = td_em(
            icru.batches, [(p.energy, p.histories) for p in parts], hz, float(info["r_plan_mm"])
        )  # type: ignore[arg-type]
    except (RecordError, st.StatisticsError, ValueError) as exc:
        return out | {"status": "FAIL", "reasons": [str(exc)], "td": None}
    return out | {"status": td["status"], "td": td}


def evaluate_science(
    seal: Seal, ion: Mapping[tuple, IonmcRun], topas: Mapping[int, Mapping[int, TopasRun]]
) -> dict[str, Any]:
    pred = seal.predictions
    res: dict[str, Any] = {"vc1": {}, "vc1a": {}, "vc2": {}, "vc2d": {}}
    for c in CASES:
        info = case_info(pred, c)
        icru = _ion(ion, "vc-sample", c, MODEL_ICRU90, "warp-cpu")
        dflt = _ion(ion, "vc-sample", c, MODEL_DEFAULT, "warp-cpu")
        res["vc1"][str(c)] = vc1_case(c, info, icru, topas[c], "V-C1")
        a = vc1_case(c, info, dflt, topas[c], "V-C1-A")
        a["gating"] = False  # report-only: never changes a verdict, the default or a spec
        res["vc1a"][str(c)] = a
        # V-C2
        gated = bool(info["vc2_precondition_holds"])
        entry: dict[str, Any] = {"case": c, "gating": gated, "precondition_holds": gated}
        if icru.valid and dflt.valid and icru.batches is not None and dflt.batches is not None:
            try:
                rio.same_edges(icru.edges_z, dflt.edges_z, f"V-C2 case {c}")  # type: ignore[arg-type]
                hz = float(icru.edges_z[1] - icru.edges_z[0])  # type: ignore[index]
                jk = jackknife_delta_r80(icru.batches, dflt.batches, hz, float(info["r_plan_mm"]))
                entry |= jk | vc2_verdict(
                    jk["delta_mm"], jk["sigma_mm"], float(info["delta_r_pred_mm"])
                )
            except (RecordError, st.StatisticsError, ValueError) as exc:
                entry |= {"verdict": "FAIL", "note": str(exc)}
        else:
            entry |= {"verdict": "FAIL", "note": "an ionmc run of the pair is invalid"}
        res["vc2"][str(c)] = entry
        # V-C2-D
        informative = bool(info["vc2d_informative"])
        for b in BACKENDS:
            ri = _ion(ion, "vc2d-sample", c, MODEL_ICRU90, b)
            rd = _ion(ion, "vc2d-sample", c, MODEL_DEFAULT, b)
            e2: dict[str, Any] = {"case": c, "backend": b, "gating": informative}
            if not informative:
                e2["gate"] = info["vc2d_gate"]
            if ri.valid and rd.valid:
                tols = info["vc2d_round_trip_tolerance_mm"]
                e2 |= vc2d_verdict(
                    ri.doc["z_m_mm"], ri.doc["r_e0_mm"], rd.doc["z_m_mm"], rd.doc["r_e0_mm"],
                    float(info["delta_r_pred_d_mm"]), float(tols[MODEL_ICRU90]),
                    float(tols[MODEL_DEFAULT]), float(info["vc2d_difference_tolerance_mm"]),
                )  # fmt: skip
            else:
                e2 |= {"verdict": "FAIL", "note": "a run of the pair is invalid"}
            res["vc2d"][f"{c}/{b}"] = e2
    return res


def integration_state(sci: Mapping[str, Any]) -> dict[str, Any]:
    """Section F, the part decided by this evaluator: V-C2-D PASS on both backends at every energy
    where it is informative (substitute gate elsewhere), V-C2 PASS where its precondition holds,
    V-C1 graded (never omitted)."""
    blocking = [
        f"V-C2-D {k}" for k, v in sci["vc2d"].items() if v["gating"] and v["verdict"] != "PASS"
    ]
    blocking += [
        f"V-C2 case {k}" for k, v in sci["vc2"].items() if v["gating"] and v["verdict"] != "PASS"
    ]
    graded = all(v["status"] in ("PASS", "FAIL") for v in sci["vc1"].values())
    if not graded:
        blocking.append("V-C1 not graded")
    return {
        "science_gates": "PASS" if not blocking else "BLOCKED",
        "blocking": blocking,
        "substitute_gate_energies": sorted({k.split("/")[0] for k, v in sci["vc2d"].items() if not v["gating"]}),
        "not_evaluated_here": ["A.1-A.5 rows", "the CC-T/CC-ID/PER-1 substitute gate", "step S6 documentation", "Codex review, validation record, CI, merge"],
    }  # fmt: skip


def evaluate(
    repo: Path | str,
    seal_sha: str,
    code_sha: str,
    records: Path | str,
    out: Path | str,
    *,
    launch_check: bool = True,
    environment: Callable[[], dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Run the whole evaluation and write ``<out>/v13_result.json``; RecordError stops it.
    ``launch_check`` and ``environment`` exist for the unit tests; the CLI always checks."""
    repo, records, out = Path(repo).resolve(), Path(records).resolve(), Path(out)
    if launch_check:  # first: refuse any launch that is not the sealed one
        try:
            env_seal.check_launch(LAUNCH, Path(__file__).resolve())
        except env_seal.EnvironmentSealError as exc:
            raise RecordError(f"launch: {exc}") from exc
    seal = load_seal(repo, seal_sha, environment, launch_check)
    if code_sha != seal.head:
        raise RecordError(f"--code-sha {code_sha} is not the evaluated HEAD {seal.head}")
    bad = sealed_problems(seal, code_sha, "evaluator commit")
    if bad:
        raise RecordError("; ".join(bad))
    ion_entries, top_entries = load_index(records, seal.predictions)
    ion_runs = [check_ionmc_run(seal, records, e) for e in ion_entries]
    for r in ion_runs:
        r.problems = list(dict.fromkeys(r.problems))
    mods = warp_modules_problems(
        [
            {"run_id": r.run_id, "backend": r.backend, "warp_modules": r.doc.get("warp_modules")}
            for r in ion_runs
        ]
    )
    for r in ion_runs:
        if r.run_id in mods:
            r.problems.append(mods[r.run_id])
    ion = {(r.row, r.case, r.model, r.backend): r for r in ion_runs}
    topas: dict[int, dict[int, TopasRun]] = {c: {} for c in CASES}
    for (c, p), e in sorted(top_entries.items()):
        topas[c][p] = check_topas_run(seal, records, e, c, p)
    sci = evaluate_science(seal, ion, topas)
    vc3 = _read_json(repo / VC3, "v-c3.json") if (repo / VC3).is_file() else None
    all_pass = all(v["status"] == "PASS" for v in sci["vc1"].values())
    result = {
        "schema": RESULT_SCHEMA,
        "seal": {
            "seal_sha": seal.seal_sha, "evaluated_head": seal.head, "manifest_sha256": seal.manifest_sha256,
            "sealed_path_count": len(seal.files), "src_ionmc_tree": seal.tree,
            "seal_review": {"id": seal.seal_review["review_id"], "finished_at": seal.seal_review["finished_at"]},
            "samples_review": {"id": seal.samples_review["review_id"], "finished_at": seal.samples_review["finished_at"]},
            "trust_sha256": seal.trust_sha256, "environment": seal.environment,
        },  # fmt: skip
        "runs": {"ionmc": [r.summary() for r in ion_runs],
                 "topas": [t.summary() for c in CASES for _, t in sorted(topas[c].items())]},  # fmt: skip
        "warp_modules_runs": {b: sorted(r.run_id for r in ion_runs if r.backend == b) for b in BACKENDS},
        **sci,
        "vc3": {"note": "report-only; committed aggregates", "source": VC3, "content": vc3},
        "evidence_record": EVIDENCE_PASS if all_pass else EVIDENCE_FAIL,
        "integration": integration_state(sci),
    }  # fmt: skip
    out.mkdir(parents=True, exist_ok=True)
    (out / RESULT_NAME).write_text(
        json.dumps(_safe_json(result), indent=1, sort_keys=True, allow_nan=False) + "\n"
    )
    return result


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Sealed V3-013 evaluator (launch via sealed_launch.py)"
    )
    ap.add_argument("--seal-sha", required=True)
    ap.add_argument("--code-sha", required=True)
    ap.add_argument("--records", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--repo", default=str(REPO_ROOT))
    a = ap.parse_args(argv)
    try:
        res = evaluate(a.repo, a.seal_sha, a.code_sha, a.records, a.out)
    except RecordError as exc:
        print(f"FAIL CLOSED: {exc}", file=sys.stderr)
        return 2
    print(f"{res['evidence_record']}; science gates {res['integration']['science_gates']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
