"""V3-008 revision 2: fail-closed loaders of the TOPAS-agreement evaluator (numpy + stdlib + the
pure-numpy TOPAS parsers of ``ionmc.reference.parsers``).

Three inputs, each bound to a controlled record before any number is used:

1. **Plan and seal review.** ``seal_sha`` / ``seal_review_id`` come from the single seal-note
   JSON block of ``validation/plans/v3-008-acceptance.md`` as committed at HEAD (the one fenced
   ```json block containing the exact substring ``"seal_sha"``, closing quote included).
   **Held-out seal (amendment R2-STOP3 item 3):** the held-out cases 10-12 are graded against a
   second, held-out seal. Its note is the one fenced ```json block of the plan containing
   ``"seal_sha_heldout"`` and exactly the keys ``seal_sha_heldout`` / ``seal_review_id_heldout``
   (``"seal_sha_heldout"`` does not contain ``"seal_sha"``, so the main search is unaffected);
   its review record is ``<records>/review-heldout.json`` (:data:`REVIEW_RECORD_HELDOUT` at the
   default records directory; same shape as ``review.json``). The seal kind of a case is
   :func:`seal_kind_of` (``CASES`` set "held-out" -> ``heldout``, else ``main``) and selects note,
   review record, seal record, seed base and npz subdirectory everywhere below (:data:`SEAL_KINDS`);
   held-out grading and ``check-run`` of a held-out bundle never fall back to the main seal.
   The review record ``<records>/review.json`` is read as committed at HEAD (fields used:
   ``review_id, head_sha, status, finished_at, report.verdict``; all required;
   ``review_id``/``head_sha`` must match the seal note). It is a verbatim copy of the controlled
   codex-worker response (its ``report_sha256`` is the service's report digest); a signed review
   record is a future infrastructure item (V3-INFRA-005), so the binding today is to the
   reviewed, committed bytes.
   All other records (index, signed run records, listings, seal, trust file) are likewise read
   through :class:`CommittedReader`.
   **Chronology pinning (Codex REVIEW-096af6a3):** the chronology of each run uses the seal note
   and ``review.json`` (held-out case: the held-out note and ``review-heldout.json``) as committed
   at the run's SIGNED ``code_sha`` (``git show <code_sha>:...``;
   they predate the run start and cannot be altered later); ``code_sha`` must be an ancestor of
   HEAD and the HEAD seal note must equal the pinned one; the report records the review digests at
   ``code_sha`` and at HEAD. **Bundle binding:** the frozen bundle manifest
   ``validation/results/transport/v8r2-topas/frozen-bundles.json`` must be byte-identical at the
   run's ``code_sha`` and at HEAD (amendment rule: a bundle change after a run invalidates that
   run's grading; for a SMOKE checked by ``check-run`` only the smoke bundle's own manifest entry
   must be identical, see :func:`frozen_manifest`), and the run's ``request.input_files`` sha256
   of ``case.json`` / ``input.txt``
   must equal the manifest entry of the bundle named for the index's case/config/part (else a
   run-validity FAIL); the frozen module/parameter-line checks remain as an additional layer.
2. **ionmc samples.** The sealed record ``sample-seal.json`` (``--seal``) must be byte-identical to
   the file committed at the plan's ``seal_sha`` (``git show <seal_sha>:<SEAL_RECORD_PATH>``).
   Each chunk npz (``<samples-root>/validation/generated/<run_id>/r2-samples/case<c>/<npz>``) is
   read only after its size and sha256 equal the sealed ``npz_bytes`` / ``npz_sha256``; the chunk
   set of a case must equal the sealed chunk plan (k = 1..n, seed formula, histories summing
   to 1e6).
   **Held-out samples (R2-STOP3):** ``sample-seal-heldout-r2stop3.json`` (schema
   ``ionmc-v8r2-sample-seal-heldout-1``), byte-identical to the file committed at the plan's
   ``seal_sha_heldout``; it carries ``seed_base`` = 20621004 (seed = seed_base + 100 (29 + case) +
   k, ``r_index`` = 29 + case), ``npz_subdir`` = ``r2-samples-heldout`` (chunks at
   ``<samples-root>/validation/generated/<run_id>/r2-samples-heldout/case<c>/<npz>``) and an
   ``implementation_sha`` that must be an ancestor of (or equal to) ``seal_sha_heldout``.
   **Host-run binding (Codex REVIEW-2f0b4eec finding 1), both seal kinds:** before any number of a
   chunk is read, :func:`verify_host_records` requires the complete signed host-run record of the
   chunk's run, ``validation/results/transport/v8r2-samples/host-records/<run_id>.json``
   (:data:`HOST_RECORDS_DIR`; verbatim record of ``infrastructure/host_runner/host_runner.py``) as
   committed at HEAD, ``verified`` for signer ``host-runner`` against the committed trust file; its
   attestation must be the chunk's sealed ``host_runner_attestation``, its signed ``sha`` the
   seal's ``implementation_sha``, its run a clean success of task V3-008, its times / stdout digest
   / ``output_manifest_sha256`` those of the seal (chunk and ``host_runs``), and the signed output
   manifest must be reproduced from the published run directory
   ``<samples-root>/validation/generated/<run_id>/`` and hold exactly the npz, step stdout and
   stderr of the run's sealed chunks (npz with the sealed sha256, ``case<c>-k<k>.stdout.txt`` with
   the sealed ``step_stdout_sha256``; no unsealed output); that stdout's one step document must
   state exactly the chunk (:func:`_step_bound`, Codex REVIEW-e844de8a). The signed ``argv`` must be
   ``["python", "-c", <source>]`` with a source byte-identical to one of the FOUR reviewed wrapper
   templates of ``r2_wrapper_templates.py`` (T1/T2/T3 main, T4 held-out; each runs ``steps_v8.py
   td-tl-sample`` once per k with the chunk's case/k/scale, writes the subprocess stdout/stderr
   verbatim to the published files and propagates failure through its exit code) rendered for the
   run's case, sealed ks and scale (:func:`sample_invocation`, Codex REVIEW-84382c07); any other
   signed source fails closed. The report lists the verified record per chunk
   (``ionmc_chunks[].host_record``).
3. **TOPAS runs.** Per run the orchestrator materializes into ``<records>/runs/<run_id>/``:
   ``record.json`` = the SIGNED run record of the reference runner (V3-INFRA-004, schema_version 3;
   fields used: ``run_id, configuration_id, code_sha, task_id`` (may be null), ``engine,
   request_sha256, started_at, finished_at, seconds, exit_code, timed_out, missing_outputs,
   succeeded, stdout_sha256, stderr_sha256`` and the ``attestation`` block; the runner's unsigned
   ``result.json`` archive index is never used as evidence) and
   ``artifacts.json`` = the runner's artifact listing (one page or a list of pages of
   ``{"run_id", "files": [{"path", "sha256", "bytes"?}]}``). ``<records>/
   index.json`` binds ``{"case", "config", "part", "run_id"}`` (one run per expected part). The
   output files themselves are read from ``<outputs>/<run_id>/<path>`` (the materialize layout,
   default ``<repo>/.ionmc-cache/reference-runs``), each verified against the listing's sha256.
   Binding of a run to its bundle: ``request.json`` (sha256 from the listing) must hash
   (canonical JSON) to the record's ``configuration_id``, carry the record's ``code_sha``, the
   seed ``20619000 + 100 part + case`` and ``input_files`` equal to the bundle
   ``validation/reference_cases/topas/<bundle>/`` as committed at ``code_sha``.

``reference_started`` = the record's ``started_at`` (ISO-8601 with offset), captured by the runner
at its controlled ``reference_started`` event, and nothing else (``finished_at`` and ``seconds`` are
never used for the chronology; Codex REVIEW-01b6acd6, orchestrator decision 2026-11-09). The run
record must carry an Ed25519 attestation (``ionmc-attestation-1``, signer ``reference-runner``,
``infrastructure/experiment_v3/attest.py``) that ``attest.classify_record`` classifies as
``verified`` against the trust file (default: the committed
``infrastructure/keys/trusted_attestation_keys.json``; no bootstrap tolerance: ``unverified`` and
``invalid`` both fail closed). Its signed ``run_id`` (equal to the index entry and of the
form ``REF-<configuration_id[:20]>-<8 hex>``), ``configuration_id`` (canonical hash of
``request.json``), ``request_sha256`` (sha256 of the request.json bytes read), ``code_sha`` (the
request's), ``stdout_sha256`` / ``stderr_sha256`` (the listing's and the hashed files) must all
agree, ``started_at <= finished_at``, ``record.json`` itself must be in the artifact listing with
the same sha256 (replay protection, as ``host_run_records.verify_run_record``).
The signed ``output_files_sha256`` (and ``output_count``) must be reproduced from the listing's
``work/`` entries (``{path: {"sha256", "bytes"}}``, canonical JSON digest, as the runner computes
it), and every output file read must hash to its listed sha256: the outputs are thereby bound to the
signed record. The trust map is the committed trust file at HEAD (:func:`committed_trust`).
Until the runner signs its records and records ``started_at`` (infra task V3-INFRA-004) every grade
fails closed at this step, by design.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

from ionmc.reference.parsers import ParseError, parse_topas_binary, parse_topas_csv

SEAL_RECORD_PATH = "validation/results/transport/v8r2-samples/sample-seal.json"
SEAL_SCHEMA = "ionmc-v8r2-sample-seal-1"
# amendment R2-STOP3 item 3 (iv): the held-out samples, re-run at base 20621004, own seal record
SEAL_RECORD_PATH_HELDOUT = (
    "validation/results/transport/v8r2-samples/sample-seal-heldout-r2stop3.json"
)
SEAL_SCHEMA_HELDOUT = "ionmc-v8r2-sample-seal-heldout-1"
HELDOUT_SEED_BASE = 20621004  # R2-STOP3 item 3 (iv): re-run base of the held-out samples
NPZ_SUBDIR, NPZ_SUBDIR_HELDOUT = "r2-samples", "r2-samples-heldout"
BUNDLE_ROOT = "validation/reference_cases/topas"
N_BATCHES = 40  # section 6: ionmc 40 batches
N_HISTORIES = 1_000_000  # section 6: ionmc 1e6 histories per case
LV8_BASE, R_INDEX_OFFSET = 20581004, 29  # section 7 (as steps_v8.py)
TOPAS_SEED_BASE = 20619000  # section 7: seed = 20619000 + 100 part + case
FIELD_MM = 120.0  # section 2: water box 120 x 120 mm
WATER_G_PER_MM3 = 1.0e-3  # section 4.1: "slab mass (120 x 120 mm x h_z at 1.0 g/cm3)"
J_PER_MEV = 1.602176634e-13  # exact SI definition of the electronvolt (unit conversion only)
EDGE_RTOL = 1e-9
# binary output stem -> TOPAS scorer name (bundles' input.txt, Sc/<name>/OutputFile)
SCORER_OF = {
    "lat_x": "LatX",
    "lat_y": "LatY",
    "lat_x_primary": "LatXPrimary",
    "lat_y_primary": "LatYPrimary",
}
_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0}

# r2-topas-target.md section 3: case -> (ion, MeV/u, R_plan mm, box L mm, h_z mm, h_x mm, set)
CASES: dict[int, tuple[str, float, float, float, float, float, str]] = {
    1: ("he", 50.0, 22.5, 50.0, 0.02, 0.01, "main"),
    2: ("he", 150.0, 159.3, 250.0, 0.1, 0.02, "main"),
    3: ("he", 250.0, 382.6, 500.0, 0.1, 0.02, "main"),
    4: ("carbon", 80.0, 17.5, 50.0, 0.02, 0.01, "main"),
    5: ("carbon", 290.0, 163.2, 250.0, 0.1, 0.02, "main"),
    6: ("carbon", 430.0, 309.3, 450.0, 0.1, 0.02, "main"),
    7: ("oxygen", 100.0, 19.6, 50.0, 0.02, 0.01, "main"),
    8: ("oxygen", 350.0, 166.6, 250.0, 0.1, 0.02, "main"),
    9: ("oxygen", 430.0, 231.9, 350.0, 0.1, 0.02, "main"),
    10: ("he", 100.0, 78.0, 150.0, 0.05, 0.02, "held-out"),
    11: ("carbon", 150.0, 53.1, 100.0, 0.05, 0.02, "held-out"),
    12: ("oxygen", 200.0, 65.5, 100.0, 0.05, 0.02, "held-out"),
}
# section 3: (low end, standard point, upper end) of each ion in the main set
ION_ENDS: dict[str, tuple[int, int, int]] = {
    "he": (1, 2, 3),
    "carbon": (4, 5, 6),
    "oxygen": (7, 8, 9),
}
# section 7: part blocks per configuration (parts 16-19, the FULL fallback, need an amendment)
CONFIG_PARTS: dict[str, tuple[int, ...]] = {
    "idd-em": (0, 1, 2, 3),
    "lat-em": (4, 5, 6, 7),
    "idd-full": (8, 9, 10, 11),
    "lat-full": (12, 13, 14, 15),
}
REFERENCE_SIGNER = "reference-runner"
DEFAULT_TRUST = "infrastructure/keys/trusted_attestation_keys.json"
RESULT_FIELDS = (
    "schema_version",
    "started_at",
    "request_sha256",
    "stdout_sha256",
    "stderr_sha256",
    "output_files_sha256",
    "output_count",
    "run_id",
    "configuration_id",
    "code_sha",
    "engine",
    "finished_at",
    "seconds",
    "exit_code",
    "succeeded",
    "missing_outputs",
)
REVIEW_FIELDS = ("review_id", "head_sha", "status", "finished_at", "report")


class RecordError(Exception):
    """A record or input file is missing, malformed or not bound to the graded run: fail closed."""


# ---------------------------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------------------------


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def canonical_digest(value: Any) -> str:
    """sha256 of canonical JSON (sorted keys, compact separators): the runner's configuration id."""
    blob = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(blob).hexdigest()


def read_json(path: Path, what: str) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RecordError(f"{what}: cannot read {path}: {exc}") from exc


def git_show(repo: Path | str, sha: str, rel: str) -> bytes:
    proc = subprocess.run(  # noqa: S603 - fixed argv
        ["git", "-C", str(repo), "show", f"{sha}:{rel}"], capture_output=True, check=False
    )
    if proc.returncode != 0:
        raise RecordError(f"{rel} is not committed at {sha}: {proc.stderr.decode().strip()}")
    return proc.stdout


def git_ls(repo: Path | str, sha: str, rel_dir: str) -> list[str]:
    proc = subprocess.run(  # noqa: S603 - fixed argv
        ["git", "-C", str(repo), "ls-tree", "-r", "--name-only", sha, "--", rel_dir + "/"],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RecordError(f"git ls-tree {sha} {rel_dir}: {proc.stderr.strip()}")
    return [ln for ln in proc.stdout.splitlines() if ln]


def require(record: Any, key: str, what: str) -> Any:
    if not isinstance(record, dict) or key not in record or record[key] is None:
        raise RecordError(f"missing record field: {what} ({key})")
    if isinstance(record[key], str) and not record[key]:
        raise RecordError(f"empty record field: {what} ({key})")
    return record[key]


def case_info(case: int) -> dict[str, Any]:
    if case not in CASES:
        raise RecordError(f"unknown case {case}; cases are 1..12")
    ion, e_u, r_plan, length, hz, hx, cset = CASES[case]
    return {
        "case": case,
        "ion": ion,
        "mev_u": e_u,
        "r_plan_mm": r_plan,
        "length_mm": length,
        "hz_mm": hz,
        "hx_mm": hx,
        "set": cset,
    }


def bundle_name(case: int, config: str, part: int) -> str:
    """``<ion>-water-<E>mevu-r2-<config>-part<p>`` (section 7)."""
    ion, e_u = CASES[case][0], CASES[case][1]
    return f"{ion}-water-{int(e_u)}mevu-r2-{config}-part{part}"


def topas_seed(case: int, part: int) -> int:
    return TOPAS_SEED_BASE + 100 * part + case


_CONFIG_RE = "(idd-em|lat-em|idd-full|lat-full)"
_PART_BUNDLE = re.compile(rf"^(he|carbon|oxygen)-water-(\d+)mevu-r2-{_CONFIG_RE}-part(\d+)$")
_SMOKE_BUNDLE = re.compile(
    rf"^(he|carbon|oxygen)-water-(\d+)mevu-r2-{_CONFIG_RE}(-ecut1km|-lp)?"
    r"-(timing|steps-verbose)-smoke$"
)


def parse_bundle(bundle: str) -> dict[str, Any]:
    """``{"kind": "part"|"smoke", "case", "config", "part"}`` of an r2 bundle directory name.
    Evidence parts: ``<ion>-water-<E>mevu-r2-<config>-part<p>`` (section 7; ``part`` must belong
    to ``config``). Smokes: ``<ion>-water-<E>mevu-r2-<config>[-ecut1km]-(timing|steps-verbose)-
    smoke`` and the ListProcesses variants (R2-SMOKE-LP) ``<ion>-water-<E>mevu-r2-<config>-lp-
    (timing|steps-verbose)-smoke`` (``-lp-`` never combines with ``-ecut1km-``); the configuration
    (frozen module set and lines of :func:`run_checks`) is the name's ``<config>`` (the
    ``-ecut1km-`` smokes are LAT-FULL under R2-SIZE-1), ``part`` is None. The
    case is the section-3 case with that ion and energy; anything else raises RecordError."""
    m = _PART_BUNDLE.match(bundle)
    kind = "part"
    if m is None:
        m, kind = _SMOKE_BUNDLE.match(bundle), "smoke"
    if m is None:
        raise RecordError(f"{bundle!r} is not an r2 evidence-part or smoke bundle name")
    ion, energy, config = m.group(1), int(m.group(2)), m.group(3)
    cases = [c for c, v in CASES.items() if v[0] == ion and int(v[1]) == energy]
    if len(cases) != 1:
        raise RecordError(f"{bundle}: no unique section-3 case for {ion} {energy} MeV/u")
    part = int(m.group(4)) if kind == "part" else None
    if part is not None and part not in CONFIG_PARTS[config]:
        raise RecordError(f"{bundle}: part {part} is not in the {config} block")
    return {"kind": kind, "case": cases[0], "config": config, "part": part}


# ---------------------------------------------------------------------------------------------
# plan, review, index
# ---------------------------------------------------------------------------------------------


PLAN_PATH = "validation/plans/v3-008-acceptance.md"
DEFAULT_RECORDS = "validation/results/transport/v8r2-topas/records"
REVIEW_NAME, REVIEW_NAME_HELDOUT = "review.json", "review-heldout.json"
REVIEW_RECORD_HELDOUT = f"{DEFAULT_RECORDS}/{REVIEW_NAME_HELDOUT}"
HELDOUT_NOTE_KEYS = ("seal_sha_heldout", "seal_review_id_heldout")


@dataclass(frozen=True)
class SealKind:
    """Everything that differs between the main seal and the held-out seal (R2-STOP3 item 3)."""

    name: str  # "main" | "heldout"
    record_path: str
    schema: str
    seed_base: int
    npz_subdir: str
    review_name: str  # file name inside the records directory


SEAL_KINDS = {
    "main": SealKind("main", SEAL_RECORD_PATH, SEAL_SCHEMA, LV8_BASE, NPZ_SUBDIR, REVIEW_NAME),
    "heldout": SealKind(
        "heldout", SEAL_RECORD_PATH_HELDOUT, SEAL_SCHEMA_HELDOUT, HELDOUT_SEED_BASE,
        NPZ_SUBDIR_HELDOUT, REVIEW_NAME_HELDOUT,
    ),
}  # fmt: skip


def seal_kind_of(case: int) -> str:
    """``heldout`` for a section-3 held-out case (10-12), ``main`` otherwise (incl. case 0)."""
    spec = CASES.get(case)
    return "heldout" if spec is not None and spec[6] == "held-out" else "main"


class CommittedReader:
    """Reads repository files only as committed at one resolved HEAD and fails closed when the
    working copy differs (Codex REVIEW-a50374b1: no caller-controlled input). ``reads`` records
    every path read with its sha256 for the report."""

    def __init__(self, repo: Path | str) -> None:
        self.repo = Path(repo).resolve()
        head = subprocess.run(  # noqa: S603 - fixed argv
            ["git", "-C", str(self.repo), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=False,
        )  # fmt: skip
        if head.returncode != 0:
            raise RecordError(f"cannot resolve HEAD of {self.repo}: {head.stderr.strip()}")
        self.head = head.stdout.strip()
        self.reads: dict[str, str] = {}

    def rel(self, path: Path | str) -> str:
        p = Path(path)
        p = (p if p.is_absolute() else self.repo / p).resolve()
        try:
            return p.relative_to(self.repo).as_posix()
        except ValueError as exc:
            raise RecordError(f"{path} is outside the evaluator repository") from exc

    def read(self, path: Path | str) -> bytes:
        rel = self.rel(path)
        blob = git_show(self.repo, self.head, rel)
        work = self.repo / rel
        if work.is_symlink() or not work.is_file() or work.read_bytes() != blob:
            raise RecordError(f"the working copy of {rel} differs from HEAD {self.head}")
        self.reads[rel] = hashlib.sha256(blob).hexdigest()
        return blob

    def json(self, path: Path | str, what: str) -> Any:
        blob = self.read(path)
        try:
            return json.loads(blob.decode("utf-8"))
        except ValueError as exc:
            raise RecordError(f"{what}: {path} is not JSON: {exc}") from exc


def load_plan(reader: CommittedReader, kind: str = "main") -> dict[str, Any]:
    """``seal_sha`` / ``seal_review_id`` from the one fenced ```json block of the committed
    acceptance plan that contains "seal_sha" (zero or several blocks: fail closed); ``kind``
    "heldout": the held-out note (:func:`parse_heldout_seal_note`), returned in the same
    normalized shape ``{"seal_sha", "seal_review_id"}``."""
    return parse_note(reader.read(PLAN_PATH).decode("utf-8"), kind)


def parse_note(text: str, kind: str) -> dict[str, Any]:
    if kind == "main":
        return parse_seal_note(text)
    if kind == "heldout":
        return parse_heldout_seal_note(text)
    raise RecordError(f"unknown seal kind {kind!r}")


def parse_seal_note(text: str) -> dict[str, Any]:
    blocks = [b for b in re.findall(r"```json[ \t]*\n(.*?)```", text, re.S) if '"seal_sha"' in b]
    if len(blocks) != 1:
        raise RecordError(
            f"{PLAN_PATH}: expected exactly one seal-note JSON block, found {len(blocks)}"
        )
    try:
        plan = json.loads(blocks[0])
    except ValueError as exc:
        raise RecordError(f"{PLAN_PATH}: the seal-note block is not JSON: {exc}") from exc
    for key in ("seal_sha", "seal_review_id"):
        require(plan, key, "plan seal note")
    return {"seal_sha": plan["seal_sha"], "seal_review_id": plan["seal_review_id"]}


def parse_heldout_seal_note(text: str) -> dict[str, Any]:
    """The dated held-out seal note (R2-STOP3 item 3): the one fenced ```json block containing
    ``"seal_sha_heldout"``, with exactly the keys :data:`HELDOUT_NOTE_KEYS` (zero or several
    blocks, other keys or a non-object: fail closed). Returned normalized as ``{"seal_sha",
    "seal_review_id"}`` so that the chronology checks are the same functions for both seals."""
    blocks = [
        b for b in re.findall(r"```json[ \t]*\n(.*?)```", text, re.S) if '"seal_sha_heldout"' in b
    ]
    if len(blocks) != 1:
        raise RecordError(
            f"{PLAN_PATH}: expected exactly one held-out seal-note JSON block, found {len(blocks)}"
        )
    try:
        note = json.loads(blocks[0])
    except ValueError as exc:
        raise RecordError(f"{PLAN_PATH}: the held-out seal-note block is not JSON: {exc}") from exc
    if not isinstance(note, dict) or sorted(note) != sorted(HELDOUT_NOTE_KEYS):
        raise RecordError(
            f"{PLAN_PATH}: the held-out seal note must hold exactly the keys {HELDOUT_NOTE_KEYS}"
        )
    return {"seal_sha": note["seal_sha_heldout"], "seal_review_id": note["seal_review_id_heldout"]}


def load_review(
    reader: CommittedReader, records: Path | str, plan: dict[str, Any], kind: str = "main"
) -> dict[str, Any]:
    """The seal review record as committed at HEAD (``<records>/review.json``; reported only, the
    chronology uses the copy at each run's code_sha, :func:`load_pinned_chronology`). It is a
    verbatim copy of the controlled codex-worker response (its ``report_sha256`` is the
    service's report digest); a signed review record is a future infrastructure item
    (V3-INFRA-005), so the binding today is to the reviewed, committed bytes. ``review_id`` and
    ``head_sha`` must equal the plan's seal note. ``kind`` "heldout": ``review-heldout.json``
    against the (normalized) held-out note."""
    name = SEAL_KINDS[kind].review_name
    return check_review(reader.json(Path(records) / name, f"seal review record {name}"), plan)


def check_review(review: Any, plan: dict[str, Any]) -> dict[str, Any]:
    for key in REVIEW_FIELDS:
        require(review, key, "seal review record")
    require(review["report"], "verdict", "seal review report")
    if review["review_id"] != plan["seal_review_id"] or review["head_sha"] != plan["seal_sha"]:
        raise RecordError("the committed review record is not the seal review named by the plan")
    return review


def load_pinned_chronology(
    repo: Path | str, code_sha: str, records_rel: str, kind: str = "main"
) -> dict[str, Any]:
    """The seal note and the seal review record as committed at the run's signed ``code_sha``
    (Codex REVIEW-096af6a3): the runner executed exactly that commit, so these bytes existed before
    the run started and cannot be changed afterwards. A later commit (e.g. a forged earlier
    ``finished_at``) does not affect the chronology of this run. Missing at ``code_sha``: fail
    closed. ``kind`` "heldout" (a run of a held-out case): the held-out note and
    ``review-heldout.json`` at ``code_sha`` (the main note is not consulted)."""
    plan_blob = git_show(repo, code_sha, PLAN_PATH)
    plan = parse_note(plan_blob.decode("utf-8"), kind)
    name = SEAL_KINDS[kind].review_name
    rev_blob = git_show(repo, code_sha, f"{records_rel}/{name}")
    try:
        review = json.loads(rev_blob.decode("utf-8"))
    except ValueError as exc:
        raise RecordError(f"{name} at {code_sha} is not JSON: {exc}") from exc
    check_review(review, plan)
    return {
        "seal_kind": kind, "plan": plan, "review": review,
        "plan_sha256": hashlib.sha256(plan_blob).hexdigest(),
        "review_sha256": hashlib.sha256(rev_blob).hexdigest(),
    }  # fmt: skip


def is_ancestor(repo: Path | str, older: str, newer: str) -> bool:
    proc = subprocess.run(  # noqa: S603 - fixed argv
        ["git", "-C", str(repo), "merge-base", "--is-ancestor", older, newer],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    if proc.returncode not in (0, 1):
        raise RecordError(f"git merge-base failed: {proc.stderr.strip()}")
    return proc.returncode == 0


FROZEN_MANIFEST = "validation/results/transport/v8r2-topas/frozen-bundles.json"
MANIFEST_SCHEMA = "ionmc-v8r2-frozen-bundles-1"


def _manifest_doc(blob: bytes, where: str) -> dict[str, Any]:
    try:
        doc = json.loads(blob.decode("utf-8"))
    except ValueError as exc:
        raise RecordError(f"{FROZEN_MANIFEST} at {where} is not JSON: {exc}") from exc
    if require(doc, "schema", "frozen manifest") != MANIFEST_SCHEMA:
        raise RecordError(f"{FROZEN_MANIFEST}: schema {doc['schema']!r} != {MANIFEST_SCHEMA}")
    require(doc, "bundle_set_digest", "frozen manifest")
    if not isinstance(require(doc, "bundles", "frozen manifest"), dict):
        raise RecordError(f"{FROZEN_MANIFEST}: bundles must be an object")
    return doc


def frozen_manifest(
    reader: CommittedReader, code_sha: str, smoke_bundle: str | None = None
) -> dict[str, Any]:
    """The frozen bundle manifest as committed at HEAD, bound to the run's ``code_sha`` by run kind
    (amendment rule: a bundle change after a run invalidates that run's grading):

    - evidence part (``smoke_bundle`` None): the manifest must be byte-identical at ``code_sha``
      and at HEAD. An evidence row is graded against the complete frozen bundle set, so any bundle
      change after the run voids it.
    - smoke (``smoke_bundle`` = its bundle name; Codex REVIEW-8286be8c allowed "a reviewed
      versioned/per-bundle amendment that preserves the intended binding"): only the smoke
      bundle's own entry (its ``case.json`` / ``input.txt`` sha256 pair) must be present and
      identical at ``code_sha`` and at HEAD; the whole-manifest digest may differ and is only
      reported. A smoke is an execution-only sizing record bound to its own frozen bundle, so a
      later re-sizing of other (evidence) bundles does not change what that smoke executed.

    A violation of either rule raises RecordError (ERROR, fail closed). The returned document is
    HEAD's, with ``_binding`` ("set" | "entry") and ``_bundle_set_digest_at_code_sha`` added."""
    at_run = git_show(reader.repo, code_sha, FROZEN_MANIFEST)
    at_head = reader.read(FROZEN_MANIFEST)
    if smoke_bundle is None and at_run != at_head:
        raise RecordError(
            f"{FROZEN_MANIFEST} differs between the run's code_sha {code_sha} and HEAD "
            f"{reader.head}: amendment rule, bundles changed after the run invalidate its grading"
        )
    doc = _manifest_doc(at_head, f"HEAD {reader.head}")
    run_doc = _manifest_doc(at_run, f"code_sha {code_sha}")
    if smoke_bundle is not None:
        old, new = run_doc["bundles"].get(smoke_bundle), doc["bundles"].get(smoke_bundle)
        pair = ("case.json", "input.txt")
        if not isinstance(old, dict) or not isinstance(new, dict):
            raise RecordError(
                f"smoke bundle {smoke_bundle} is not in {FROZEN_MANIFEST} at both code_sha "
                f"{code_sha} and HEAD {reader.head}"
            )
        if any(old.get(k) is None or old.get(k) != new.get(k) for k in pair):
            raise RecordError(
                f"the frozen-manifest entry of smoke bundle {smoke_bundle} differs between the "
                f"run's code_sha {code_sha} and HEAD {reader.head}: the smoke's bundle changed "
                "after the run"
            )
    doc["_binding"] = "set" if smoke_bundle is None else "entry"
    doc["_bundle_set_digest_at_code_sha"] = run_doc["bundle_set_digest"]
    return doc


def load_index(reader: CommittedReader, records: Path | str) -> list[dict[str, Any]]:
    index = reader.json(Path(records) / "index.json", "run index")
    runs = require(index, "runs", "run index")
    if not isinstance(runs, list):
        raise RecordError("run index: runs must be a list")
    for e in runs:
        for key in ("case", "config", "part", "run_id"):
            require(e, key, "run index entry")
    return runs


def runs_for(index: list[dict[str, Any]], case: int, config: str) -> list[tuple[int, str]]:
    """(part, run_id) for the K = 4 parts of ``config`` (exactly one run per part)."""
    parts = CONFIG_PARTS[config]
    found: dict[int, str] = {}
    for e in index:
        if e["case"] == case and e["config"] == config:
            if e["part"] not in parts:
                raise RecordError(f"case {case} {config}: part {e['part']} is not in {parts}")
            if e["part"] in found:
                raise RecordError(f"case {case} {config}: two runs listed for part {e['part']}")
            found[e["part"]] = str(e["run_id"])
    missing = [p for p in parts if p not in found]
    if missing:
        raise RecordError(f"case {case} {config}: no run record for parts {missing}")
    return [(p, found[p]) for p in parts]


# ---------------------------------------------------------------------------------------------
# TOPAS runs
# ---------------------------------------------------------------------------------------------


@dataclass
class TopasRun:
    run_id: str
    case: int
    config: str
    part: int | None
    bundle: str
    configuration_id: str
    code_sha: str
    finished_at: str
    seconds: float
    reference_started: str
    histories: int
    expected_outputs: list[str]
    listing: dict[str, str]
    outputs_dir: Path
    read: dict[str, str] = field(default_factory=dict)
    case_doc: dict[str, Any] = field(default_factory=dict)
    request: dict[str, Any] = field(default_factory=dict)
    native: str = ""
    pinned: dict[str, Any] = field(default_factory=dict)
    bundle_problems: list[str] = field(default_factory=list)

    def chronology_record(self) -> dict[str, Any]:
        """The run record handed to ``check_run_validity``."""
        return {
            "run_id": self.run_id,
            "code_sha": self.code_sha,
            "reference_started": self.reference_started,
        }

    def summary(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "config": self.config,
            "part": self.part,
            "bundle": self.bundle,
            "configuration_id": self.configuration_id,
            "code_sha": self.code_sha,
            "finished_at": self.finished_at,
            "seconds": self.seconds,
            "reference_started": self.reference_started,
            "histories": self.histories,
            "files_sha256": dict(sorted(self.read.items())),
        }


def _listing(
    raw: Any, run_id: str, code_sha: str, sizes: dict[str, Any] | None = None
) -> dict[str, str]:
    """path -> sha256 of the artifact listing; ``sizes`` (if given) receives path -> bytes."""
    pages = raw if isinstance(raw, list) else [raw]
    files: dict[str, str] = {}
    for page in pages:
        if require(page, "run_id", "artifact listing") != run_id:
            raise RecordError(f"artifact listing belongs to {page['run_id']!r}, not {run_id!r}")
        if "sha" in page and page["sha"] != code_sha:
            raise RecordError("artifact listing sha differs from the run record code_sha")
        entries = require(page, "files", "artifact listing")
        if not isinstance(entries, list):
            raise RecordError("artifact listing: files must be a list")
        for e in entries:
            p, d = require(e, "path", "artifact entry"), require(e, "sha256", "artifact entry")
            if p in files and files[p] != d:
                raise RecordError(f"artifact listing: conflicting entries for {p}")
            files[p] = d
            if sizes is not None:
                sizes[p] = e.get("bytes")
    return files


def _parse_time(value: Any, what: str) -> Any:
    from datetime import UTC, datetime

    if not isinstance(value, str):
        raise RecordError(f"{what} is not a timestamp string: {value!r}")
    try:
        t = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise RecordError(f"{what} is not ISO-8601: {value!r}") from exc
    if t.tzinfo is None:
        raise RecordError(f"{what} has no UTC offset: {value!r}")
    return t.astimezone(UTC)


def load_topas_run(
    repo: Path | str,
    records: Path | str,
    outputs: Path | str,
    case: int,
    config: str,
    part: int,
    run_id: str,
    trust: dict[str, Any],
    reader: CommittedReader,
) -> TopasRun:
    """Load and bind one TOPAS evidence run (see the module docstring); any gap raises
    RecordError. ``trust`` is ``attest.load_trusted(<trust file>)``. The record lives under
    ``<records>/runs/<run_id>/`` and the bundle is ``bundle_name(case, config, part)`` with the
    seed formula of section 7."""
    return load_signed_run(
        repo, Path(records) / "runs" / run_id, outputs, run_id, bundle_name(case, config, part),
        config, trust, reader, review_records=records, case=case, part=part,
        expected_seeds=[topas_seed(case, part)],
    )  # fmt: skip


def load_signed_run(
    repo: Path | str,
    rdir: Path | str,
    outputs: Path | str,
    run_id: str,
    bundle: str,
    config: str,
    trust: dict[str, Any],
    reader: CommittedReader,
    *,
    review_records: Path | str,
    case: int = 0,
    part: int | None = None,
    expected_seeds: list[int] | None = None,
    manifest_binding: str = "set",
) -> TopasRun:
    """Load and bind one signed TOPAS run of ``bundle`` (evidence part or smoke) whose
    ``record.json`` / ``artifacts.json`` are in ``rdir`` (see the module docstring for every
    binding); the seal note and ``<review_records>/review.json`` are pinned at the run's
    ``code_sha``. ``expected_seeds``: the seed formula of an evidence part; None for a smoke, whose
    seed is bound only through the committed bundle at ``code_sha`` (case.json ``seeds``).
    ``manifest_binding``: "set" (evidence parts: whole frozen manifest identical at ``code_sha``
    and HEAD) or "entry" (smokes: only the bundle's own manifest entry; :func:`frozen_manifest`)."""
    if manifest_binding not in ("set", "entry"):
        raise RecordError(f"unknown manifest binding {manifest_binding!r}")
    rdir = Path(rdir)
    record_path = rdir / "record.json"
    result = reader.json(record_path, f"signed run record {run_id}")
    status, reason = attest_module().classify_record(result, trust, REFERENCE_SIGNER)
    if status != "verified":
        raise RecordError(
            f"{run_id}: record.json is not attested by a trusted {REFERENCE_SIGNER!r} key "
            f"(status {status}: {reason}); signed run records with started_at are required "
            "(gap until the reference runner signs them, V3-INFRA-004)"
        )
    for key in RESULT_FIELDS:
        require(result, key, f"result record {run_id}")
    if result["schema_version"] != 3:  # V3-INFRA-004: signed record with started_at
        raise RecordError(f"{run_id}: result schema_version {result['schema_version']!r} != 3")
    if "task_id" not in result:
        raise RecordError(f"missing record field: run record {run_id} (task_id)")
    if result.get("timed_out") is not False:
        raise RecordError(f"{run_id}: timed_out is {result.get('timed_out')!r}, not false")
    if not re.fullmatch(
        r"REF-" + re.escape(str(result["configuration_id"])[:20]) + r"-[0-9a-f]{8}", run_id
    ):
        raise RecordError(f"{run_id}: run id is not derived from the configuration_id")
    if result["run_id"] != run_id:
        raise RecordError(f"result record {result['run_id']!r} is not bound to run {run_id!r}")
    if result["engine"] != "topas":
        raise RecordError(f"{run_id}: engine {result['engine']!r} is not topas")
    if not (
        result["succeeded"] is True and result["exit_code"] == 0 and result["missing_outputs"] == []
    ):
        raise RecordError(
            f"{run_id}: the run did not succeed (succeeded={result['succeeded']!r}, exit_code="
            f"{result['exit_code']!r}, missing_outputs={result['missing_outputs']!r})"
        )
    seconds = result["seconds"]
    if not isinstance(seconds, int | float) or not math.isfinite(seconds) or seconds < 0:
        raise RecordError(f"{run_id}: seconds is not a non-negative number: {seconds!r}")
    finished = _parse_time(result["finished_at"], f"{run_id} finished_at")
    started = _parse_time(result["started_at"], f"{run_id} started_at")
    if started > finished:
        raise RecordError(f"{run_id}: started_at is later than finished_at")
    started_iso = str(result["started_at"])
    code_sha = str(result["code_sha"])
    sizes: dict[str, Any] = {}
    listing = _listing(
        reader.json(rdir / "artifacts.json", "artifact listing"), run_id, code_sha, sizes
    )
    # Bind the outputs to the signed record: the runner signs output_files_sha256 =
    # common.digest(inventory(run_dir / "work", allow_internal_links=True)), i.e. the canonical
    # JSON digest of {path relative to work/: {"sha256", "bytes"}} (V3-INFRA-004 reference.py).
    # The unsigned listing's work/ entries must reproduce it exactly (and output_count).
    work = {}
    for p, d in listing.items():
        if p.startswith("work/"):
            b = sizes.get(p)
            if not isinstance(b, int) or isinstance(b, bool) or b < 0:
                raise RecordError(f"{run_id}: artifact listing lacks the byte size of {p}")
            work[p[len("work/") :]] = {"sha256": d, "bytes": b}
    if canonical_digest(work) != result.get("output_files_sha256"):
        raise RecordError(
            f"{run_id}: the listed work/ outputs do not reproduce the signed output_files_sha256"
        )
    if result.get("output_count") != len(work):
        raise RecordError(f"{run_id}: signed output_count differs from the listed outputs")
    for name in ("request.json", "stdout.txt", "stderr.txt", "record.json"):
        if name not in listing:
            raise RecordError(f"{run_id}: artifact listing lacks {name}")
    if reader.reads[reader.rel(record_path)] != listing["record.json"]:
        raise RecordError(f"{run_id}: record.json differs from the listed record.json")
    for name, key in (("stdout.txt", "stdout_sha256"), ("stderr.txt", "stderr_sha256")):
        if result[key] != listing[name]:
            raise RecordError(f"{run_id}: signed {key} differs from the listing's {name}")
    run = TopasRun(
        run_id=run_id,
        case=case,
        config=config,
        part=part,
        bundle=bundle,
        configuration_id=str(result["configuration_id"]),
        code_sha=code_sha,
        finished_at=str(result["finished_at"]),
        seconds=float(seconds),
        reference_started=started_iso,
        histories=0,
        expected_outputs=[],
        listing=listing,
        outputs_dir=Path(outputs) / run_id,
    )
    for name in ("stdout.txt", "stderr.txt"):
        read_verified(run, name)  # the files themselves hash to the listing (= signed) digests
    request_blob = read_verified(run, "request.json")
    if hashlib.sha256(request_blob).hexdigest() != result["request_sha256"]:
        raise RecordError(f"{run_id}: request.json does not hash to the signed request_sha256")
    try:
        request = json.loads(request_blob.decode("utf-8"))
    except ValueError as exc:
        raise RecordError(f"{run_id}: request.json is not JSON: {exc}") from exc
    if canonical_digest(request) != run.configuration_id:
        raise RecordError(f"{run_id}: request.json does not hash to the configuration_id")
    if require(request, "code_sha", "request") != code_sha:
        raise RecordError(f"{run_id}: request code_sha differs from the result record")
    case_doc = require(request, "case", "request")
    input_files = require(request, "input_files", "request")
    bundle_dir = f"{BUNDLE_ROOT}/{bundle}"
    committed: dict[str, dict[str, Any]] = {}
    for rel in git_ls(repo, code_sha, bundle_dir):
        blob = git_show(repo, code_sha, rel)
        committed[rel[len(bundle_dir) + 1 :]] = {
            "sha256": hashlib.sha256(blob).hexdigest(),
            "bytes": len(blob),
        }
    if not committed:
        raise RecordError(f"{run_id}: bundle {bundle} is not committed at {code_sha}")
    if not is_ancestor(repo, code_sha, reader.head):
        raise RecordError(f"{run_id}: code_sha {code_sha} is not an ancestor of HEAD {reader.head}")
    manifest = frozen_manifest(reader, code_sha, bundle if manifest_binding == "entry" else None)
    entry = manifest["bundles"].get(bundle)
    if not isinstance(entry, dict):
        run.bundle_problems.append(f"bundle {bundle} is not in the frozen manifest")
    else:
        for name in ("case.json", "input.txt"):
            got = input_files.get(name, {}).get("sha256") if isinstance(input_files, dict) else None
            if got is None or got != entry.get(name):
                run.bundle_problems.append(f"{bundle}/{name} differs from the frozen manifest")
    run.pinned = load_pinned_chronology(  # held-out case: held-out note + review (R2-STOP3)
        repo, code_sha, reader.rel(review_records), seal_kind_of(case)
    )
    run.pinned["manifest_bundle_set_digest"] = manifest["bundle_set_digest"]
    run.pinned["manifest_bundle_set_digest_at_code_sha"] = manifest[
        "_bundle_set_digest_at_code_sha"
    ]
    run.pinned["manifest_binding"] = manifest["_binding"]
    if input_files != committed:
        raise RecordError(
            f"{run_id}: request input_files differ from the committed bundle {bundle}"
        )
    case_json = json.loads(git_show(repo, code_sha, f"{bundle_dir}/case.json").decode("utf-8"))
    if case_doc != case_json:
        raise RecordError(f"{run_id}: request case differs from the committed case.json")
    seeds = require(case_doc, "seeds", "case.json")
    if expected_seeds is not None and seeds != expected_seeds:
        raise RecordError(f"{run_id}: seed {seeds!r} != formula {expected_seeds}")
    hist = require(case_doc, "histories", "case.json")
    if not isinstance(hist, int) or hist <= 0:
        raise RecordError(f"{run_id}: histories {hist!r} is not a positive integer")
    run.histories = hist
    run.expected_outputs = list(require(case_doc, "expected_outputs", "case.json"))
    run.case_doc, run.request = case_doc, request
    native_name = require(case_doc, "input", "case.json")
    run.native = git_show(repo, code_sha, f"{bundle_dir}/{native_name}").decode("utf-8", "replace")
    return run


# ---------------------------------------------------------------------------------------------
# run validity of each TOPAS part (section 6 "Run validity (each TOPAS part; a violation
# invalidates the part, which is re-run at a fresh re-run base declared by amendment, the invalid
# part preserved)")
# ---------------------------------------------------------------------------------------------

_BEGIN = re.compile(r"Begin processing", re.IGNORECASE)
_KILLED = re.compile(r"killed|stuck", re.IGNORECASE)
_HISTORIES = re.compile(r"Particle source Beam: Total number of histories: (\d+)")
_MSC_HEAD = re.compile(r"^msc:\s+for\s+(\S+)\s+SubType", re.M)
_MODEL = re.compile(r"^\s*(\S+)\s*:\s*Emin=", re.M)
MSC_PARTICLES = ("GenericIon", "alpha")


def attest_module() -> Any:
    """``infrastructure/experiment_v3/attest.py`` of this checkout (pure stdlib)."""
    import sys

    root = str(Path(__file__).resolve().parents[3])
    if root not in sys.path:
        sys.path.insert(0, root)
    from infrastructure.experiment_v3 import attest

    return attest


def load_trust(path: Path | str) -> dict[str, Any]:
    """Trust map of a trust FILE (tests only; production uses :func:`committed_trust`)."""
    try:
        return dict(attest_module().load_trusted(path))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise RecordError(f"trust file {path} unreadable: {exc}") from exc


def committed_trust(reader: CommittedReader) -> tuple[dict[str, Any], dict[str, Any]]:
    """The trust map of ``infrastructure/keys/trusted_attestation_keys.json`` as committed at the
    evaluator repository's HEAD (Codex REVIEW-6b4731a7 / REVIEW-a50374b1); a working copy that
    differs fails closed. Returns (trust map, provenance)."""
    blob = reader.read(DEFAULT_TRUST)
    try:
        data = json.loads(blob.decode("utf-8"))
        trust = {s: {e["key_id"]: e["public_key"] for e in entries} for s, entries in data.items()}
    except (ValueError, KeyError, TypeError, AttributeError) as exc:
        raise RecordError(f"{DEFAULT_TRUST} at {reader.head} is malformed: {exc}") from exc
    info = {"path": DEFAULT_TRUST, "head_sha": reader.head, "sha256": reader.reads[DEFAULT_TRUST]}
    return trust, info


# Module -> printout marker (Geant4 11.4.2 / OpenTOPAS 4.3 stdout). Each marker was checked, format
# only, against materialized r2 smoke logs under .ionmc-cache/reference-runs (FULL: e.g.
# REF-171fec5d..., EM opt4: e.g. REF-188fa1d8..., EM SS: REF-285075c5... as the opt4 negative).
_MODULES_LINE = re.compile(r"^[ \t]*sv:Ph/Default/Modules[ \t]*=[ \t]*(\d+)[ \t]+(.*)$", re.M)
_HAD_HEAD = re.compile(r"^\s*Hadronic Processes for (\S+)\s*$")
_REGISTERED = re.compile(r"^Registered Physics Processes:\s*$", re.M)
HADRONIC_MARKERS = {  # forbidden in an EM-only (no hadronic module) run
    "hInelastic": re.compile(r"\bhInelastic\b"),
    "QMDModel": re.compile(r"QMDModel"),
    "Radioactive Decay Physics Parameters": re.compile(r"Radioactive Decay Physics Parameters"),
    "HADRONIC PROCESSES SUMMARY": re.compile(r"HADRONIC PROCESSES SUMMARY"),
}


def registered_processes(stdout: str) -> list[str] | None:
    """Process names of the TOPAS ``Ph/ListProcesses`` printout, or None if it is absent. Format
    (confirmed on materialized EM smoke logs, e.g. REF-188fa1d8...): a line ``Registered Physics
    Processes:`` followed by comma-separated process names over several lines, ended by a blank
    line. The printout lists names only, without the particle they are attached to."""
    m = _REGISTERED.search(stdout)
    if m is None:
        return None
    block = stdout[m.end() :].lstrip("\n").split("\n\n", 1)[0]
    return [t for t in re.split(r"[,\s]+", block) if t]


def hadronic_block(stdout: str, particle: str) -> str:
    """Lines of the "Hadronic Processes for <particle>" block of the hadronic summary (up to the
    next dashed separator); all such blocks concatenated."""
    out, inside = [], False
    for ln in stdout.splitlines():
        m = _HAD_HEAD.match(ln)
        if m:
            inside = m.group(1) == particle
            continue
        if inside and ln.lstrip().startswith("-----"):
            inside = False
        if inside:
            out.append(ln)
    return "\n".join(out)


def _process_models(block: str, process: str) -> list[str]:
    """Model lines under ``Process: <process>`` of one hadronic block."""
    out, inside = [], False
    for ln in block.splitlines():
        if ln.strip().startswith("Process:"):
            inside = ln.strip() == f"Process: {process}"
            continue
        if inside:
            out.append(ln)
    return out


def _marker(module: str, stdout: str) -> bool:
    if module == "g4em-standard_opt4":  # ICRU90 stopping on and GS msc for e- (option4 only)
        icru = re.search(r"^Use ICRU90 data\s+1\s*$", stdout, re.M) is not None
        return icru and "GoudsmitSaunderson" in msc_models(stdout, "e-")
    if module == "g4h-phy_QGSP_BIC_HP":
        return re.search(r"^\s*hInelastic QGSP_BIC_HP\b", stdout, re.M) is not None
    if module == "g4ion-QMD":
        models = _process_models(hadronic_block(stdout, "GenericIon"), "ionInelastic")
        return any("QMDModel" in ln for ln in models)
    if module == "g4radioactivedecay":
        return "Radioactive Decay Physics Parameters" in stdout
    if module == "g4h-elastic_HP":
        models = _process_models(hadronic_block(stdout, "neutron"), "hadElastic")
        return any(re.search(r"Model:\s*NeutronHPElastic\b", ln) for ln in models)
    if module == "g4stopping":
        return "Process: muMinusCaptureAtRest" in hadronic_block(stdout, "mu-")
    if module == "g4decay":
        # Requested rule: a process "Decay" listed for pi+/pi-/mu+/mu-/kaon+/neutron. The confirmed
        # ListProcesses format carries no particle association, so the check is: the name "Decay"
        # (G4Decay's process name; exact token, not RadioactiveDecay) is registered. Presence of
        # the token in a g4decay run is UNCONFIRMED until the first FULL ListProcesses smoke.
        names = registered_processes(stdout)
        return names is not None and "Decay" in names
    raise KeyError(module)


VERIFIABLE_MODULES = (
    "g4em-standard_opt4", "g4h-phy_QGSP_BIC_HP", "g4ion-QMD", "g4radioactivedecay",
    "g4h-elastic_HP", "g4stopping", "g4decay",
)  # fmt: skip


def bundle_modules(native: str) -> list[str]:
    """The ``sv:Ph/Default/Modules`` list of the committed native input (exactly one line, count
    equal to the number of quoted names)."""
    lines = _MODULES_LINE.findall(native)
    if len(lines) != 1:
        raise RecordError(f"expected exactly one sv:Ph/Default/Modules line, found {len(lines)}")
    n, rest = int(lines[0][0]), lines[0][1]
    names = re.findall(r'"([^"]+)"', rest)
    if len(names) != n:
        raise RecordError(f"Modules line declares {n} modules but names {names}")
    return names


EM_MODULES = ("g4em-standard_opt4",)
FULL_MODULES = (
    "g4em-standard_opt4", "g4h-phy_QGSP_BIC_HP", "g4decay", "g4h-elastic_HP", "g4stopping",
    "g4ion-QMD", "g4radioactivedecay",
)  # fmt: skip
FROZEN_MODULES = {"idd-em": EM_MODULES, "lat-em": EM_MODULES, "idd-full": FULL_MODULES,
                  "lat-full": FULL_MODULES}  # fmt: skip
# frozen native lines per configuration (section 2; R2-SIZE-1: LAT-FULL electron cut 1 km;
# ListProcesses amendment: process list in every evidence bundle); exactly once each
_COMMON_LINES = ('b:Ph/ListProcesses = "True"', "d:Ph/Default/CutForAllParticles = 0.05 mm")
_ECUT = "d:Ph/Default/CutForElectron = 1 km"
FROZEN_LINES = {"idd-em": (*_COMMON_LINES, _ECUT), "lat-em": (*_COMMON_LINES, _ECUT),
                "idd-full": _COMMON_LINES, "lat-full": (*_COMMON_LINES, _ECUT)}  # fmt: skip


def frozen_lines_check(native: str, config: str) -> list[str]:
    """The frozen physics lines of ``config`` appear exactly once (whitespace-normalised, comment
    lines ignored); IDD-FULL must carry no CutForElectron line."""
    lines = [" ".join(ln.split()) for ln in native.splitlines() if not ln.lstrip().startswith("#")]
    problems = [f"frozen line {w!r} missing or repeated" for w in FROZEN_LINES[config]
                if lines.count(w) != 1]  # fmt: skip
    if config == "idd-full" and any("CutForElectron" in ln for ln in lines):
        problems.append("IDD-FULL bundle sets CutForElectron (not frozen for this configuration)")
    return problems


def module_check(native: str, stdout: str, config: str) -> tuple[list[str], list[str]]:
    """(problems, unverifiable modules) of the executed module list. Interim rule (Codex
    REVIEW-6b4731a7, orchestrator 2026-11-09): every bundled module must show its runtime marker
    (see :func:`_marker`); a bundled module without a marker is a run-validity FAIL as well (it is
    also returned as unverifiable); a run whose bundle lists only ``g4em-*`` modules must show no
    hadronic marker and no registered ``Decay`` process. The bundle's module list must equal the
    frozen set of its configuration (``FROZEN_MODULES``; order-insensitive, Codex REVIEW-a50374b1);
    the runtime markers are checked for exactly that frozen set."""
    try:
        modules = bundle_modules(native)
    except RecordError as exc:
        return [str(exc)], []
    problems, unverifiable = [], []
    frozen = FROZEN_MODULES[config]
    if sorted(modules) != sorted(frozen):
        problems.append(f"bundle modules {modules} != the frozen {config} set {list(frozen)}")
        modules = list(frozen)  # still check the runtime markers of the frozen set
    for mod in modules:
        if mod in VERIFIABLE_MODULES:
            if not _marker(mod, stdout):
                problems.append(f"module {mod}: printout marker missing")
        else:
            unverifiable.append(mod)
            problems.append(f"module {mod}: no runtime evidence marker is defined (fail closed)")
    if all(m.startswith("g4em-") for m in modules):
        for name, rx in HADRONIC_MARKERS.items():
            if rx.search(stdout):
                problems.append(f"EM-only run shows the hadronic marker {name!r}")
        if "Decay" in (registered_processes(stdout) or []):
            problems.append("EM-only run registers the Decay process")
    return problems, unverifiable


def msc_models(stdout: str, particle: str) -> list[str]:
    """Model names of every ``msc:  for <particle>`` block of the Geant4 EM printout (the block
    ends at the first blank line)."""
    models: list[str] = []
    for m in _MSC_HEAD.finditer(stdout):
        if m.group(1) != particle:
            continue
        block = stdout[m.end() :].split("\n\n", 1)[0]
        models.extend(_MODEL.findall(block))
    return models


def run_checks(run: TopasRun) -> dict[str, list[str]]:
    """Problems of one part (empty = valid): stderr must be empty ("zero killed or stuck tracks in
    stderr", applied strictly: any stderr output invalidates); no stdout line outside the
    "Begin processing" progress lines mentions killed or stuck tracks; "Particle source Beam: Total
    number of histories: N" present with every N equal to the bundle's histories
    (NumberOfHistoriesInRun completed); the msc printout of GenericIon and alpha lists UrbanMsc
    only; the executed module list matches the bundle (:func:`module_check`). Returns
    ``{"problems": [...], "module_markers_unverifiable": [...]}``."""
    problems = list(run.bundle_problems)  # frozen-manifest binding (REVIEW-096af6a3)
    stderr = read_verified(run, "stderr.txt")
    stdout = read_verified(run, "stdout.txt").decode("utf-8", "replace")
    if stderr.strip():
        problems.append("stderr is not empty")
    bad = [ln for ln in stdout.splitlines() if _KILLED.search(ln) and not _BEGIN.search(ln)]
    if bad:
        problems.append(f"stdout reports killed/stuck tracks ({len(bad)} line(s))")
    found = [int(v) for v in _HISTORIES.findall(stdout)]
    if not found or any(v != run.histories for v in found):
        problems.append(f"histories completed {found} != NumberOfHistoriesInRun {run.histories}")
    for p in MSC_PARTICLES:
        models = msc_models(stdout, p)
        if not models or any(m != "UrbanMsc" for m in models):
            problems.append(f"msc models for {p} are {models}, not UrbanMsc")
    problems += frozen_lines_check(run.native, run.config)
    mod_problems, unverifiable = module_check(run.native, stdout, run.config)
    return {"problems": problems + mod_problems, "module_markers_unverifiable": unverifiable}


def replicate_gate(runs: list[TopasRun]) -> list[str]:
    """Section 6: "the replicate gate of compare_batches.py (identical configuration except the
    seed) holds within each group" (group = the K parts of one case and configuration). Reuses the
    seed-normalisation constants of ``compare_batches.py`` (its ``i:Ts/Seed`` line pattern, the
    case fields allowed to differ ``seeds``/``rationale``, the request identity keys) on the
    committed bundle bound to each run. Deviation from ``compare_batches.run_fingerprint``:
    full-line ``#`` comments of the native input are also removed (TOPAS ignores them; every r2
    bundle names its part in its first comment line, so the unmodified fingerprint would reject
    every group).
    Its other protocol conditions (>= 1e5 histories, depth bin <= 0.5 mm, an IDD in every run,
    equal source commit) are not part of the plan's definition and are not applied."""
    cb = _compare_batches()
    problems = []
    seeds = [r.case_doc.get("seeds") for r in runs]
    if len({json.dumps(x) for x in seeds}) != len(runs):
        problems.append(f"seeds are not distinct: {seeds}")
    prints = set()
    for r in runs:
        nat_seed = re.findall(r"^[ \t]*i:Ts/Seed[ \t]*=[ \t]*(\d+)[ \t]*$", r.native, re.M)
        if [int(v) for v in nat_seed] != list(r.case_doc.get("seeds") or []):
            problems.append(f"{r.run_id}: native Ts/Seed {nat_seed} != case seeds")
        nat_hist = [int(v) for v in cb._NATIVE_HISTORIES["topas"].findall(r.native)]
        if nat_hist != [r.histories]:
            problems.append(f"{r.run_id}: native histories {nat_hist} != {r.histories}")
        native = cb._SEED_LINE["topas"].sub("", r.native)
        native = "\n".join(ln for ln in native.splitlines() if not ln.lstrip().startswith("#"))
        case = {k: v for k, v in r.case_doc.items() if k not in cb._CASE_SEED_FIELDS}
        identity = {k: r.request.get(k) for k in cb._IDENTITY_KEYS}
        if identity.get("engine") is None or identity.get("dirty") is not False:
            problems.append(f"{r.run_id}: request lacks the engine identity or is not clean")
        prints.add(canonical_digest({"native": native, "case": case, "identity": identity}))
    if len(prints) != 1:
        problems.append("parts differ in the executed configuration beyond the seed")
    return problems


def _compare_batches() -> Any:
    import importlib.util
    import sys

    path = Path(__file__).resolve().parents[1] / "reference" / "compare_batches.py"
    name = "r2_compare_batches"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RecordError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def read_verified(run: TopasRun, rel: str) -> bytes:
    """Bytes of ``<outputs>/<run_id>/<rel>`` after the sha256 equals the artifact listing."""
    if rel not in run.listing:
        raise RecordError(f"{run.run_id}: {rel} is not in the artifact listing")
    path = run.outputs_dir / rel
    if path.is_symlink() or not path.is_file():
        raise RecordError(f"{run.run_id}: {rel} is not materialized as a regular file")
    blob = path.read_bytes()
    digest = hashlib.sha256(blob).hexdigest()
    if digest != run.listing[rel]:
        raise RecordError(f"{run.run_id}: {rel} sha256 mismatch against the artifact listing")
    run.read[rel] = digest
    return blob


def _mm(width: float, unit: str, what: str) -> float:
    if unit not in _TO_MM:
        raise RecordError(f"{what}: unknown length unit {unit!r}")
    return width * _TO_MM[unit]


def topas_idd(run: TopasRun, name: str = "idd_dose.csv") -> tuple[np.ndarray, float]:
    """Run-total energy per depth bin [MeV] of a TOPAS IDD (DoseToMedium, Gy, 1 x 1 x nz over the
    full 120 x 120 mm field): dose x slab mass (120 x 120 mm x h_z at 1.0 g/cm3) / (J per MeV),
    section 4.1. Returns (energy per bin, h_z in mm). Missing bins fail closed."""
    if name not in run.expected_outputs:
        raise RecordError(f"{run.run_id}: {name} is not an expected output of {run.bundle}")
    rel = f"work/{name}"
    try:
        s = parse_topas_csv(read_verified(run, rel).decode("utf-8"), rel)
    except ParseError as exc:
        raise RecordError(f"{run.run_id}: {exc}") from exc
    if s.bins[:2] != (1, 1) or "Sum" not in s.values:
        raise RecordError(f"{run.run_id}: {rel} must be a 1 x 1 x N column with a Sum")
    if s.unit != "Gy":
        raise RecordError(f"{run.run_id}: {rel} unit {s.unit!r} is not Gy")
    for ax in (0, 1):
        if abs(_mm(s.bin_width[ax], s.bin_unit[ax], rel) - FIELD_MM) > 1e-6:
            raise RecordError(f"{run.run_id}: {rel} lateral extent is not {FIELD_MM} mm")
    hz = _mm(s.bin_width[2], s.bin_unit[2], rel)
    dose = s.values["Sum"][0, 0, :]
    if not np.all(np.isfinite(dose)):
        raise RecordError(f"{run.run_id}: {rel} has missing (unwritten) bins")
    mass_kg = FIELD_MM * FIELD_MM * hz * WATER_G_PER_MM3 * 1e-3
    return dose * mass_kg / J_PER_MEV, hz


def topas_lateral(run: TopasRun, stem: str, axis: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Projected lateral profile of a TOPAS binary scorer ``<stem>.bin`` + ``.binheader`` (dose,
    run total). Format: the header layout of the real r2 smoke (``# TOPAS Version``, ``# Parameter
    File``, ``# Results for scorer: LatX``, ``# Scored in component``, ``# X in 6000 bins of
    0.002 cm``, ``# Y in 1 bin  of 12 cm``, ``# Z in 500 bins of 0.1 cm``, ``# DoseToMedium ( Gy )
    : Sum``, ``# Binary file: lat_x.bin``); the scorer name and the ``Binary file`` line are bound
    to the stem; values float64 little-endian (TOPAS default; float32 is rejected here) in x-fastest
    (Fortran) order ``index = ix + nx (iy + ny iz)``, the order ``ionmc.reference.parsers``
    documents as verified on a real OpenTOPAS 4.3 run (REF-95ef7515...).
    ``axis`` 0 (x profile, n x 1 x nz) or 1 (y profile, 1 x n x nz). Returns (values (n, nz),
    lateral edges (n + 1,) centred on the axis, depth bin width mm)."""
    for ext in (".bin", ".binheader"):
        if stem + ext not in run.expected_outputs:
            raise RecordError(f"{run.run_id}: {stem}{ext} is not an expected output")
    header = read_verified(run, f"work/{stem}.binheader").decode("utf-8")
    blob = read_verified(run, f"work/{stem}.bin")
    try:
        s = parse_topas_binary(header, blob, f"work/{stem}.bin", order="F")
    except ParseError as exc:
        raise RecordError(f"{run.run_id}: {exc}") from exc
    if "Sum" not in s.values:
        raise RecordError(f"{run.run_id}: {stem} has no Sum")
    if s.meta.get("dtype") != "<f8":
        raise RecordError(f"{run.run_id}: {stem}.bin is not float64 (TOPAS default)")
    if s.meta.get("scorer") != SCORER_OF.get(stem):
        raise RecordError(f"{run.run_id}: {stem}.binheader scorer {s.meta.get('scorer')!r}")
    if f"# Binary file: {stem}.bin" not in [ln.strip() for ln in header.splitlines()]:
        raise RecordError(f"{run.run_id}: {stem}.binheader does not name {stem}.bin")
    if s.unit != "Gy":
        raise RecordError(f"{run.run_id}: {stem} unit {s.unit!r} is not Gy")
    other = 1 - axis
    if (
        s.bins[other] != 1
        or abs(_mm(s.bin_width[other], s.bin_unit[other], stem) - FIELD_MM) > 1e-6
    ):
        raise RecordError(f"{run.run_id}: {stem} is not projected over the full field")
    n = s.bins[axis]
    h = _mm(s.bin_width[axis], s.bin_unit[axis], stem)
    if abs(n * h - FIELD_MM) > 1e-6 * FIELD_MM:
        raise RecordError(f"{run.run_id}: {stem} lateral bins do not span {FIELD_MM} mm")
    vals = s.values["Sum"]
    vals = vals[:, 0, :] if axis == 0 else vals[0, :, :]
    edges = -FIELD_MM / 2.0 + h * np.arange(n + 1, dtype=np.float64)
    return vals, edges, _mm(s.bin_width[2], s.bin_unit[2], stem)


# ---------------------------------------------------------------------------------------------
# ionmc samples
# ---------------------------------------------------------------------------------------------


def load_seal(
    reader: CommittedReader, seal_path: Path | str, seal_sha: str, kind: str = "main"
) -> tuple[dict[str, Any], str]:
    """The sealed sample record of ``kind`` as committed at HEAD (path ``SealKind.record_path``),
    byte-identical to the one committed at the seal commit, with the kind's schema. Held-out
    (R2-STOP3 item 3) additionally: ``seed_base`` == 20621004, ``npz_subdir`` ==
    ``r2-samples-heldout`` and ``implementation_sha`` an ancestor of (or equal to) ``seal_sha``."""
    sk = SEAL_KINDS[kind]
    if reader.rel(seal_path) != sk.record_path:
        raise RecordError(f"the {kind} seal record must be {sk.record_path}, not {seal_path}")
    blob = reader.read(seal_path)
    if blob != git_show(reader.repo, seal_sha, sk.record_path):
        raise RecordError(f"{seal_path} differs from {sk.record_path} at the seal {seal_sha}")
    record = json.loads(blob.decode("utf-8"))
    if not isinstance(record, dict):
        raise RecordError(f"{sk.record_path} is not a JSON object")
    if require(record, "schema", "sample seal") != sk.schema:
        raise RecordError(f"sample seal schema {record['schema']!r} != {sk.schema}")
    if kind == "heldout":
        if require(record, "seed_base", "held-out seal") != sk.seed_base:
            raise RecordError(
                f"held-out seal seed_base {record['seed_base']!r} != {sk.seed_base} (R2-STOP3)"
            )
        if require(record, "npz_subdir", "held-out seal") != sk.npz_subdir:
            raise RecordError(
                f"held-out seal npz_subdir {record['npz_subdir']!r} != {sk.npz_subdir!r}"
            )
        impl = require(record, "implementation_sha", "held-out seal")
        if not isinstance(impl, str) or not re.fullmatch(r"[0-9a-f]{40}", impl):
            raise RecordError(f"held-out seal implementation_sha {impl!r} is not a commit SHA")
        if not is_ancestor(reader.repo, impl, seal_sha):
            raise RecordError(
                f"held-out seal implementation_sha {impl} is not an ancestor of the held-out "
                f"seal commit {seal_sha}"
            )
    return record, hashlib.sha256(blob).hexdigest()


def case_chunks(record: dict[str, Any], case: int, kind: str = "main") -> list[dict[str, Any]]:
    """The sealed chunks of ``case``, checked against the sealed chunk plan (k = 1..n, the lv8
    seed formula with the kind's base -- main 20581004, held-out 20621004 -- every chunk
    counter-clean, histories summing to 1e6; held-out also ``r_index`` == 29 + case). The case
    must belong to the seal ``kind`` (a held-out case is never read from the main seal)."""
    if seal_kind_of(case) != kind:
        raise RecordError(f"case {case} is not graded against the {kind} seal")
    base = SEAL_KINDS[kind].seed_base
    plan_n = require(require(record, "chunk_plan", "sample seal"), "n", "chunk plan")
    n = plan_n.get(str(case))
    chunks = [c for c in require(record, "chunks", "sample seal") if c.get("case") == case]
    if not isinstance(n, int) or n < 1:
        raise RecordError(f"sample seal: no chunk plan for case {case}")
    if sorted(c.get("k") for c in chunks) != list(range(1, n + 1)):
        raise RecordError(f"sample seal: case {case} chunks are not k = 1..{n}")
    for c in chunks:
        for key in ("k", "seed", "histories", "run_id", "npz", "npz_sha256", "npz_bytes"):
            require(c, key, f"sample seal chunk case {case}")
        if c["seed"] != base + 100 * (R_INDEX_OFFSET + case) + c["k"]:
            raise RecordError(f"sample seal: case {case} k {c['k']} seed breaks the formula")
        if kind == "heldout" and c.get("r_index") != R_INDEX_OFFSET + case:
            raise RecordError(f"sample seal: case {case} k {c['k']} r_index != {29 + case}")
        if c.get("counters_all_zero") is not True:
            raise RecordError(f"sample seal: case {case} k {c['k']} is not counter-clean")
        if c["npz"] != f"ionmc-case{case}-k{c['k']}.npz":
            raise RecordError(f"sample seal: case {case} k {c['k']} npz name {c['npz']!r}")
        if c["histories"] % N_BATCHES:
            raise RecordError(f"sample seal: case {case} k {c['k']} histories not a multiple of 40")
    if sum(c["histories"] for c in chunks) != N_HISTORIES:
        raise RecordError(f"sample seal: case {case} histories do not sum to {N_HISTORIES}")
    return sorted(chunks, key=lambda c: c["k"])


def npz_path(
    samples_root: Path | str, case: int, chunk: dict[str, Any], subdir: str = NPZ_SUBDIR
) -> Path:
    if subdir not in (NPZ_SUBDIR, NPZ_SUBDIR_HELDOUT):
        raise RecordError(f"unknown chunk subdirectory {subdir!r}")
    return (
        Path(samples_root)
        / "validation"
        / "generated"
        / str(chunk["run_id"])
        / subdir
        / f"case{case}"
        / str(chunk["npz"])
    )


def verify_npz(path: Path, chunk: dict[str, Any]) -> None:
    if path.is_symlink() or not path.is_file():
        raise RecordError(f"sealed chunk not found as a regular file: {path}")
    if path.stat().st_size != chunk["npz_bytes"]:
        raise RecordError(f"{path.name}: size differs from the sealed record")
    if sha256_file(path) != chunk["npz_sha256"]:
        raise RecordError(f"{path.name}: sha256 differs from the sealed record")


# Host-run records of the sample chunks (Codex REVIEW-2f0b4eec finding 1): every chunk of a graded
# case must be proven to be an output of the trusted host runner at the sealed implementation SHA.
HOST_RECORDS_DIR = "validation/results/transport/v8r2-samples/host-records"
HOST_SIGNER = "host-runner"
HOST_TASK_ID = "V3-008"
_ATTESTATION_CORE = ("signer", "key_id", "signed_at_utc", "payload_sha256", "signature")
STEP_ROW = "td-tl-sample"


def _exact(got: Any, want: Any) -> bool:
    """Equality with exact JSON types (bools and None by identity; 0 is not False)."""
    if isinstance(want, bool) or want is None:
        return got is want
    return type(got) is type(want) and got == want


def _bind(run_id: str, what: str, doc: dict[str, Any], expected: dict[str, Any]) -> None:
    for key, want in expected.items():
        if key not in doc:
            raise RecordError(f"host run {run_id}: {what} lacks {key!r}")
        if not _exact(doc[key], want):
            raise RecordError(f"host run {run_id}: {what} {key} {doc[key]!r} != {want!r}")


def wrapper_templates() -> Any:
    """``r2_wrapper_templates`` (next to this file): the four reviewed chunk-wrapper sources."""
    import importlib.util
    import sys

    name = "r2_wrapper_templates"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            name, Path(__file__).resolve().with_name(f"{name}.py")
        )
        assert spec is not None and spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
    mod = sys.modules[name]
    mod.check_digests()
    return mod


def sample_invocation(argv: Any, kind: str, case: int, ks: list[int], scale: str) -> dict[str, Any]:
    """Bind the signed host-run ``argv`` of one sample run (Codex REVIEW-84382c07).

    Accepted is ONLY ``["python", "-c", <source>]`` whose source is byte-identical to one of the
    four reviewed wrapper templates of ``r2_wrapper_templates`` allowed for the seal ``kind`` (T1,
    T2, T3 main; T4 held-out), rendered for this run's sealed chunks: ``case``, the sorted ``ks``
    of every sealed chunk of the run and the sealed ``scale``. Each template runs ``steps_v8.py
    td-tl-sample`` once per k with exactly that case/k/scale (T4: ``--seed-base 20621004``),
    writes the subprocess stdout/stderr verbatim to the published ``case<c>-k<k>.stdout.txt`` /
    ``.stderr.txt`` and propagates a failed step through its exit code; any other signed source
    (extra statements, a changed character, other k, another step) fails closed. Returns
    ``{"template", "template_sha256", "case", "ks", "scale", "seed_base", "out_dir"}``."""
    w = wrapper_templates()
    if (
        not isinstance(argv, list)
        or len(argv) != 3
        or not all(isinstance(t, str) for t in argv)
        or argv[0] != "python"
        or argv[1] != "-c"
    ):
        raise RecordError('signed argv is not ["python", "-c", <reviewed wrapper source>]')
    matches = []
    for t in w.TEMPLATES:
        if t.kind != kind:
            continue
        try:
            if w.render(t, case, ks, scale) == argv[2]:
                matches.append(t)
        except ValueError:
            continue
    if len(matches) != 1:
        raise RecordError(
            f"signed argv: the wrapper source is not a reviewed {kind} template rendered for case "
            f"{case}, k {ks}, scale {scale!r}"
        )
    t = matches[0]
    out_dir = f"validation/generated/{t.subdir}/case{case}"
    if t.seed_base != SEAL_KINDS[kind].seed_base or t.subdir != SEAL_KINDS[kind].npz_subdir:
        raise RecordError(f"wrapper template {t.id} does not belong to the {kind} seal")
    return {
        "template": t.id, "template_sha256": t.sha256, "case": case, "ks": sorted(ks),
        "scale": scale, "seed_base": t.seed_base, "out_dir": out_dir,
    }  # fmt: skip


def step_document(text: str, run_id: str) -> dict[str, Any]:
    """The one JSON document of a step stdout, between the lines ``#JSON-BEGIN`` and
    ``#JSON-END`` (``steps.emit``); zero or several documents fail closed."""
    lines = text.splitlines()
    begins = [i for i, x in enumerate(lines) if x == "#JSON-BEGIN"]
    ends = [i for i, x in enumerate(lines) if x == "#JSON-END"]
    if len(begins) != 1 or len(ends) != 1 or ends[0] <= begins[0]:
        raise RecordError(f"host run {run_id}: step stdout does not hold exactly one JSON document")
    try:
        doc = json.loads("\n".join(lines[begins[0] + 1 : ends[0]]))
    except ValueError as exc:
        raise RecordError(f"host run {run_id}: step document is not JSON: {exc}") from exc
    if not isinstance(doc, dict):
        raise RecordError(f"host run {run_id}: step document is not a JSON object")
    return doc


_PUBLISHED_ROOTS = ("validation/generated", "benchmarks/generated")  # host_runner.publish_outputs


def published_manifest(samples_root: Path | str, run_id: str) -> tuple[str, dict[str, str]]:
    """Reproduce the host runner's ``output_manifest_sha256`` of ``run_id`` from its published
    outputs (``infrastructure/host_runner/host_runner.py``: ``output_manifest_sha256`` over
    ``<run_dir>/outputs``, published by ``publish_outputs`` to ``<worktree>/<root>/<run_id>/`` for
    ``root`` in ``validation/generated`` and ``benchmarks/generated``). The digest is sha256 over
    the sorted lines ``f"{root}/{path relative to the published run directory}\0{file sha256}\n"``
    of every regular file. Any symlink or non-regular entry fails closed (the runner never
    publishes links). Returns (digest, {manifest path: sha256})."""
    base = Path(samples_root)
    entries: dict[str, str] = {}
    found = False
    for root in _PUBLISHED_ROOTS:
        top = base / root / run_id
        if top.is_symlink():
            raise RecordError(f"host run {run_id}: published directory {top} is a symlink")
        if not top.exists():
            continue
        if not top.is_dir():
            raise RecordError(f"host run {run_id}: published output {top} is not a directory")
        found = True
        for p in sorted(top.rglob("*")):
            if p.is_symlink():
                raise RecordError(f"host run {run_id}: published output {p} is a symlink")
            if p.is_dir():
                continue
            if not p.is_file():
                raise RecordError(f"host run {run_id}: published output {p} is not a regular file")
            entries[f"{root}/{p.relative_to(top).as_posix()}"] = sha256_file(p)
    if not found or not entries:
        raise RecordError(f"host run {run_id}: no published outputs under {base}")
    lines = "".join(f"{rel}\0{entries[rel]}\n" for rel in sorted(entries))
    return hashlib.sha256(lines.encode("utf-8")).hexdigest(), entries


def _attestation_bound(record_block: Any, chunk_block: Any, strict: bool, run_id: str) -> None:
    """The chunk's sealed ``host_runner_attestation`` must be the signed record's attestation. The
    held-out seal records the complete block (exact equality); the main seal (sealed before the
    block was recorded in full) carries the subset signer/key_id/signed_at_utc/payload_sha256/
    signature, each of which must be present and equal (``payload_sha256`` and ``signature`` bind
    the whole signed payload)."""
    if not isinstance(chunk_block, dict) or not isinstance(record_block, dict):
        raise RecordError(f"host run {run_id}: sealed host_runner_attestation is not an object")
    if strict:
        if chunk_block != record_block:
            raise RecordError(
                f"host run {run_id}: the signed record's attestation differs from the sealed "
                "host_runner_attestation"
            )
        return
    missing = [k for k in _ATTESTATION_CORE if k not in chunk_block]
    if missing:
        raise RecordError(f"host run {run_id}: sealed host_runner_attestation lacks {missing}")
    if any(k not in record_block or record_block[k] != v for k, v in chunk_block.items()):
        raise RecordError(
            f"host run {run_id}: the signed record's attestation differs from the sealed "
            "host_runner_attestation"
        )


def verify_host_records(
    reader: CommittedReader | None,
    trust: dict[str, Any],
    seal: dict[str, Any],
    case: int,
    chunks: list[dict[str, Any]],
    samples_root: Path | str,
    subdir: str,
) -> list[dict[str, Any]]:
    """Bind every sealed chunk of ``case`` to the complete signed host-run record that produced it
    (Codex REVIEW-2f0b4eec finding 1), before any number of the chunk is read. Fails closed.

    Per chunk: ``<HOST_RECORDS_DIR>/<run_id>.json`` as committed at HEAD (:class:`CommittedReader`)
    must be classified ``verified`` by ``attest.classify_record`` against the committed trust map
    with signer ``host-runner``; its attestation must be the chunk's sealed
    ``host_runner_attestation`` (:func:`_attestation_bound`); the signed fields must state a clean
    successful run of this task at the seal's ``implementation_sha`` (``run_id``, ``sha``,
    ``task_id`` V3-008, ``succeeded`` true, ``exit_code`` 0, ``timed_out`` false,
    ``dirty_before`` / ``dirty_after`` false, ``committed_state_unchanged`` true,
    ``output_publication_error`` null) and agree with the seal (``started_at`` / ``finished_at`` =
    ``run_started_at`` / ``run_finished_at``, ``stdout_sha256`` = ``wrapper_stdout_sha256``,
    ``output_manifest_sha256`` = the chunk's and the seal's single ``host_runs`` entry of the run).
    The npz is checked against the sealed size/sha256 (:func:`verify_npz`), the signed output
    manifest is reproduced from the published run directory (:func:`published_manifest`) and must
    contain the entry ``validation/generated/<subdir>/case<c>/<npz>`` with the sealed sha256. The
    step stdout and its document are bound by :func:`_step_bound` (Codex REVIEW-e844de8a). Once
    per run: the sealed chunks of the run share case and ``scale`` (1e6 * scale = histories), the
    signed ``argv`` is one of the four reviewed wrapper sources rendered for the run's case,
    sorted sealed ks and scale (:func:`sample_invocation`, Codex REVIEW-84382c07), and the signed
    output manifest holds exactly the npz, stdout and stderr of the run's sealed chunks (no
    unsealed output).
    Returns one report block per chunk."""
    if reader is None:
        raise RecordError("host-run records must be read through the committed reader")
    impl = require(seal, "implementation_sha", "sample seal")
    if not isinstance(impl, str) or not re.fullmatch(r"[0-9a-f]{40}", impl):
        raise RecordError(f"sample seal implementation_sha {impl!r} is not a commit SHA")
    strict = seal.get("schema") == SEAL_SCHEMA_HELDOUT
    kind = "heldout" if strict else "main"
    if SEAL_KINDS[kind].npz_subdir != subdir:
        raise RecordError(f"chunk subdirectory {subdir!r} is not the {kind} seal's")
    host_runs = require(seal, "host_runs", "sample seal")
    if not isinstance(host_runs, list):
        raise RecordError("sample seal host_runs is not a list")
    sealed = require(seal, "chunks", "sample seal")
    attest = attest_module()
    manifests: dict[str, tuple[str, dict[str, str]]] = {}
    invocations: dict[str, dict[str, Any]] = {}
    out = []
    for c in chunks:
        run_id = require(c, "run_id", f"sample seal chunk case {case}")
        if not isinstance(run_id, str) or not re.fullmatch(r"RUN-[0-9A-Za-z-]+", run_id):
            raise RecordError(f"sample seal: case {case} k {c.get('k')} run_id {run_id!r}")
        for key in ("run_started_at", "run_finished_at", "output_manifest_sha256",
                    "wrapper_stdout_sha256", "host_runner_attestation"):  # fmt: skip
            require(c, key, f"sample seal chunk case {case} k {c.get('k')}")
        rel = f"{HOST_RECORDS_DIR}/{run_id}.json"
        try:
            blob = reader.read(rel)
        except RecordError as exc:
            raise RecordError(
                f"host run {run_id}: no committed signed host-run record: {exc}"
            ) from exc
        try:
            rec = json.loads(blob.decode("utf-8"))
        except ValueError as exc:
            raise RecordError(f"host run {run_id}: {rel} is not JSON: {exc}") from exc
        if not isinstance(rec, dict):
            raise RecordError(f"host run {run_id}: {rel} is not a JSON object")
        status, why = attest.classify_record(rec, trust, HOST_SIGNER)
        if status != "verified":
            raise RecordError(f"host run {run_id}: host-run record attestation {status}: {why}")
        _attestation_bound(rec["attestation"], c["host_runner_attestation"], strict, run_id)
        expected = {
            "run_id": run_id, "sha": impl, "task_id": HOST_TASK_ID, "succeeded": True,
            "exit_code": 0, "timed_out": False, "dirty_before": False, "dirty_after": False,
            "committed_state_unchanged": True, "output_publication_error": None,
            "started_at": c["run_started_at"], "finished_at": c["run_finished_at"],
            "stdout_sha256": c["wrapper_stdout_sha256"],
            "output_manifest_sha256": c["output_manifest_sha256"],
        }  # fmt: skip
        _bind(run_id, "signed", rec, expected)
        signed_manifest = rec["output_manifest_sha256"]
        entries = [h for h in host_runs if isinstance(h, dict) and h.get("run_id") == run_id]
        if len(entries) != 1 or entries[0].get("output_manifest_sha256") != signed_manifest:
            raise RecordError(
                f"host run {run_id}: not exactly one seal host_runs entry with the signed "
                "output_manifest_sha256"
            )
        path = npz_path(samples_root, case, c, subdir)
        verify_npz(path, c)  # sealed size and sha256 (unchanged check)
        if run_id not in invocations:  # once per run: argv and the run's complete output set
            peers = [x for x in sealed if isinstance(x, dict) and x.get("run_id") == run_id]
            scale = c.get("scale")
            if {x.get("case") for x in peers} != {case} or {x.get("scale") for x in peers} != {
                scale
            }:
                raise RecordError(f"host run {run_id}: sealed chunks differ in case or scale")
            if not isinstance(scale, str) or round(N_HISTORIES * Fraction(scale)) != c["histories"]:
                raise RecordError(f"host run {run_id}: sealed scale {scale!r} != histories / 1e6")
            ks = sorted(x["k"] for x in peers)
            inv = sample_invocation(rec.get("argv"), kind, case, ks, scale)
            _bind(run_id, "signed argv", inv, {
                "case": case, "ks": ks, "scale": scale, "seed_base": SEAL_KINDS[kind].seed_base,
                "out_dir": f"validation/generated/{subdir}/case{case}",
            })  # fmt: skip
            digest, files = published_manifest(samples_root, run_id)
            if digest != signed_manifest:
                raise RecordError(
                    f"host run {run_id}: the published outputs do not reproduce the signed "
                    "output_manifest_sha256"
                )
            base = f"validation/generated/{subdir}/case{case}"
            expected_files = {
                f"{base}/{name}"
                for x in peers
                for name in (x["npz"], f"case{case}-k{x['k']}.stdout.txt",
                             f"case{case}-k{x['k']}.stderr.txt")
            }  # fmt: skip
            if set(files) != expected_files:
                extra = sorted(set(files) - expected_files)
                raise RecordError(
                    f"host run {run_id}: the signed output manifest is not exactly the npz, stdout "
                    f"and stderr of the sealed chunks k {ks} (unsealed or missing: "
                    f"{extra or sorted(expected_files - set(files))})"
                )
            invocations[run_id], manifests[run_id] = inv, (digest, files)
        inv = invocations[run_id]
        digest, files = manifests[run_id]
        entry = f"validation/generated/{subdir}/case{case}/{c['npz']}"
        if files.get(entry) != c["npz_sha256"]:
            raise RecordError(
                f"host run {run_id}: the signed output manifest has no entry {entry} with the "
                "sealed npz_sha256"
            )
        step = _step_bound(samples_root, run_id, case, c, files, subdir, kind, impl)
        block = rec["attestation"]
        out.append({
            "case": case, "k": c["k"], "run_id": run_id, "record_path": rel,
            "record_sha256": reader.reads[rel], "signer": HOST_SIGNER,
            "key_id": block["key_id"], "payload_sha256": block["payload_sha256"],
            "signed_at_utc": block["signed_at_utc"], "sha": rec["sha"],
            "output_manifest_sha256": signed_manifest, "manifest_entry": entry,
            "attestation_verified": True, "manifest_verified": True,
            "step_document_verified": True, "step_stdout_sha256": step,
            "invocation": {key: inv[key] for key in ("template", "template_sha256", "ks", "scale",
                                                     "seed_base", "out_dir")},
            "argv_bound": True,
        })  # fmt: skip
    return out


def _step_bound(
    samples_root: Path | str, run_id: str, case: int, c: dict[str, Any], files: dict[str, str],
    subdir: str, kind: str, impl: str,
) -> str:  # fmt: skip
    """The step stdout ``case<c>-k<k>.stdout.txt`` published next to the npz must be in the signed
    output manifest with the sealed ``step_stdout_sha256``, and its one JSON document must state
    exactly this chunk (row, suite, backend, case, k, seed, seed base, r_index, histories, npz and
    its sha256, git SHA = implementation_sha, pass/valid, every counter 0, effective n_histories /
    seed / single CPU worker). Returns the verified sha256."""
    want_sha = require(c, "step_stdout_sha256", f"sample seal chunk case {case} k {c['k']}")
    name = f"case{case}-k{c['k']}.stdout.txt"
    entry = f"validation/generated/{subdir}/case{case}/{name}"
    if files.get(entry) != want_sha:
        raise RecordError(
            f"host run {run_id}: the signed output manifest has no entry {entry} with the sealed "
            "step_stdout_sha256"
        )
    path = Path(samples_root) / "validation" / "generated" / run_id / subdir / f"case{case}" / name
    if path.is_symlink() or not path.is_file():
        raise RecordError(f"host run {run_id}: step stdout {path} is not a regular file")
    blob = path.read_bytes()
    if hashlib.sha256(blob).hexdigest() != want_sha:
        raise RecordError(f"host run {run_id}: step stdout differs from step_stdout_sha256")
    doc = step_document(blob.decode("utf-8"), run_id)
    _bind(run_id, "step document", doc, {
        "row": STEP_ROW, "suite": "lv8", "backend": "warp-cpu", "case": case, "k": c["k"],
        "seed": c["seed"], "seed_base": SEAL_KINDS[kind].seed_base,
        "r_index": R_INDEX_OFFSET + case, "histories": c["histories"],
        "frozen_histories": N_HISTORIES, "npz": c["npz"], "npz_sha256": c["npz_sha256"],
        "git_sha": impl, "pass": True, "valid": True,
    })  # fmt: skip
    counters = doc.get("counters")
    if (
        not isinstance(counters, dict)
        or not counters
        or not all(_exact(v, 0) for v in counters.values())
    ):
        raise RecordError(f"host run {run_id}: step document counters are not all zero")
    eff = doc.get("effective_config")
    if not isinstance(eff, dict):
        raise RecordError(f"host run {run_id}: step document lacks effective_config")
    _bind(run_id, "step document effective_config", eff,
          {"n_histories": c["histories"], "seed": c["seed"], "cpu_workers": 1})  # fmt: skip
    return want_sha


@dataclass
class IonmcGrid:
    """Combined per-batch energy per primary [MeV] of one grid, selected and squeezed."""

    batches: np.ndarray  # (B, ...) per primary
    edges: list[np.ndarray]  # the three edge arrays of the grid [mm]


def combine_chunks(arrays: list[tuple[int, np.ndarray]]) -> np.ndarray:
    """Per-batch combination of chunks (section 7, "per-batch grids are combined by the
    evaluator"): batch b of the case is the union of batch b of every chunk, i.e. the
    history-weighted mean of the chunks' per-primary batch arrays,
    ``E_b = sum_c h_c e_{c,b} / sum_c h_c`` with h_c the histories per batch of chunk c. This keeps
    B = 40 independent batches of N/40 histories each, exactly the batches a single 1e6-history run
    with 40 batches would have had (implementation decision, the plan says "40 batches")."""
    if not arrays:
        raise RecordError("no chunks to combine")
    shape = arrays[0][1].shape
    acc = np.zeros(shape, dtype=np.float64)
    total = 0
    for h, a in arrays:
        if a.shape != shape:
            raise RecordError(f"chunk arrays differ in shape: {a.shape} != {shape}")
        acc += float(h) * a
        total += h
    return acc / float(total)


def load_ionmc(
    samples_root: Path | str,
    case: int,
    chunks: list[dict[str, Any]],
    grid: str,
    slabs: list[int] | None = None,
    subdir: str = NPZ_SUBDIR,
) -> IonmcGrid:
    """Load grid ``idd`` (-> (B, nz)), ``latx`` (-> (B, nx, len(slabs))) or ``laty`` (-> (B, ny,
    len(slabs))) of every sealed chunk of ``case`` (sha256 verified first) and combine per batch."""
    parts: list[tuple[int, np.ndarray]] = []
    edges: list[np.ndarray] | None = None
    for c in chunks:
        path = npz_path(samples_root, case, c, subdir)
        verify_npz(path, c)
        try:
            with np.load(path, allow_pickle=False) as z:
                arr = np.asarray(z[f"{grid}_batch_energy_mev"], dtype=np.float64)
                e = [np.asarray(z[f"{grid}_edges{ax}_mm"], dtype=np.float64) for ax in range(3)]
        except (OSError, KeyError, ValueError) as exc:
            raise RecordError(f"{path.name}: cannot read grid {grid}: {exc}") from exc
        if arr.ndim != 4 or arr.shape[0] != N_BATCHES or not np.all(np.isfinite(arr)):
            raise RecordError(
                f"{path.name}: {grid} is not a finite ({N_BATCHES}, nx, ny, nz) array"
            )
        if edges is None:
            edges = e
        elif any(
            a.shape != b.shape or not np.array_equal(a, b) for a, b in zip(edges, e, strict=True)
        ):
            raise RecordError(f"{path.name}: {grid} edges differ between chunks")
        if grid == "idd":
            sel = arr[:, 0, 0, :]
        elif grid == "latx":
            sel = arr[:, :, 0, :][:, :, slabs] if slabs is not None else arr[:, :, 0, :]
        elif grid == "laty":
            sel = arr[:, 0, :, :][:, :, slabs] if slabs is not None else arr[:, 0, :, :]
        else:
            raise RecordError(f"unknown grid {grid}")
        parts.append((c["histories"] // N_BATCHES, np.array(sel)))
        del arr
    assert edges is not None
    return IonmcGrid(combine_chunks(parts), edges)


def same_edges(a: np.ndarray, b: np.ndarray, what: str) -> None:
    """Identical bin edges of both engines (section 4) or fail closed."""
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    scale = max(1.0, float(np.max(np.abs(a))) if a.size else 1.0)
    if a.shape != b.shape or not np.all(np.abs(a - b) <= EDGE_RTOL * scale):
        raise RecordError(f"{what}: ionmc and TOPAS bin edges differ")
