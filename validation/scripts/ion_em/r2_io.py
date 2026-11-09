"""V3-008 revision 2: fail-closed loaders of the TOPAS-agreement evaluator (numpy + stdlib + the
pure-numpy TOPAS parsers of ``ionmc.reference.parsers``).

Three inputs, each bound to a controlled record before any number is used:

1. **Plan and seal review.** ``--plan`` is the JSON block of the plan's dated seal note
   (``{"seal_sha", "seal_review_id"}``). The review record is the codex-worker record as returned by
   the controlled service, materialized by the orchestrator as ``<records>/review.json`` (fields
   used: ``review_id, head_sha, status, finished_at, report.verdict``; all required).
2. **ionmc samples.** The sealed record ``sample-seal.json`` (``--seal``) must be byte-identical to
   the file committed at the plan's ``seal_sha`` (``git show <seal_sha>:<SEAL_RECORD_PATH>``).
   Each chunk npz (``<samples-root>/validation/generated/<run_id>/r2-samples/case<c>/<npz>``) is
   read only after its size and sha256 equal the sealed ``npz_bytes`` / ``npz_sha256``; the chunk
   set of a case must equal the sealed chunk plan (k = 1..n, seed formula, histories summing
   to 1e6).
3. **TOPAS runs.** Per run the orchestrator materializes into ``<records>/runs/<run_id>/``:
   ``result.json`` = the controlled result record of the reference runner (schema_version 2:
   ``run_id, configuration_id, code_sha, engine, finished_at, seconds, exit_code, succeeded,
   missing_outputs``; all required) and ``artifacts.json`` = the runner's artifact listing (one page
   or a list of pages of ``{"run_id", "files": [{"path", "sha256", "bytes"?}]}``). ``<records>/
   index.json`` binds ``{"case", "config", "part", "run_id"}`` (one run per expected part). The
   output files themselves are read from ``<outputs>/<run_id>/<path>`` (the materialize layout,
   default ``<repo>/.ionmc-cache/reference-runs``), each verified against the listing's sha256.
   Binding of a run to its bundle: ``request.json`` (sha256 from the listing) must hash
   (canonical JSON) to the record's ``configuration_id``, carry the record's ``code_sha``, the
   seed ``20619000 + 100 part + case`` and ``input_files`` equal to the bundle
   ``validation/reference_cases/topas/<bundle>/`` as committed at ``code_sha``.

``reference_started`` = ``finished_at - seconds`` of the controlled result record (both written by
the controlled runner). CAVEAT: this satisfies the plan's "the runner's own start record" only until
a signed start record exists: the runner emits its ``reference_started`` event before it starts the
``seconds`` clock and stamps ``finished_at`` just before stopping it, so the derived start differs
from the event time by the runner's bookkeeping time (milliseconds); a review finishing within that
interval of a run start cannot be ordered by these records.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from dataclasses import dataclass, field
from datetime import timedelta
from pathlib import Path
from typing import Any

import numpy as np

from ionmc.reference.parsers import ParseError, parse_topas_binary, parse_topas_csv

SEAL_RECORD_PATH = "validation/results/transport/v8r2-samples/sample-seal.json"
SEAL_SCHEMA = "ionmc-v8r2-sample-seal-1"
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
RESULT_FIELDS = (
    "schema_version",
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


# ---------------------------------------------------------------------------------------------
# plan, review, index
# ---------------------------------------------------------------------------------------------


def load_plan(path: Path | str) -> dict[str, Any]:
    plan = read_json(Path(path), "plan seal note")
    for key in ("seal_sha", "seal_review_id"):
        require(plan, key, "plan seal note")
    return {"seal_sha": plan["seal_sha"], "seal_review_id": plan["seal_review_id"]}


def load_review(records: Path | str) -> dict[str, Any]:
    review = read_json(Path(records) / "review.json", "seal review record")
    for key in REVIEW_FIELDS:
        require(review, key, "seal review record")
    require(review["report"], "verdict", "seal review report")
    return review


def load_index(records: Path | str) -> list[dict[str, Any]]:
    index = read_json(Path(records) / "index.json", "run index")
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
    part: int
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


def _listing(path: Path, run_id: str, code_sha: str) -> dict[str, str]:
    raw = read_json(path, "artifact listing")
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
) -> TopasRun:
    """Load and bind one TOPAS run (see the module docstring); any gap raises RecordError."""
    rdir = Path(records) / "runs" / run_id
    result = read_json(rdir / "result.json", f"result record {run_id}")
    for key in RESULT_FIELDS:
        require(result, key, f"result record {run_id}")
    if result["schema_version"] != 2:
        raise RecordError(f"{run_id}: result schema_version {result['schema_version']!r} != 2")
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
    started = (finished - timedelta(seconds=float(seconds))).isoformat()
    code_sha = str(result["code_sha"])
    listing = _listing(rdir / "artifacts.json", run_id, code_sha)
    for name in ("request.json", "stdout.txt", "stderr.txt"):
        if name not in listing:
            raise RecordError(f"{run_id}: artifact listing lacks {name}")
    bundle = bundle_name(case, config, part)
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
        reference_started=started,
        histories=0,
        expected_outputs=[],
        listing=listing,
        outputs_dir=Path(outputs) / run_id,
    )
    try:
        request = json.loads(read_verified(run, "request.json").decode("utf-8"))
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
    if input_files != committed:
        raise RecordError(
            f"{run_id}: request input_files differ from the committed bundle {bundle}"
        )
    case_json = json.loads(git_show(repo, code_sha, f"{bundle_dir}/case.json").decode("utf-8"))
    if case_doc != case_json:
        raise RecordError(f"{run_id}: request case differs from the committed case.json")
    seeds = require(case_doc, "seeds", "case.json")
    if seeds != [topas_seed(case, part)]:
        raise RecordError(f"{run_id}: seed {seeds!r} != formula {topas_seed(case, part)}")
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


def run_checks(run: TopasRun) -> list[str]:
    """Problems of one part (empty = valid): stderr must be empty ("zero killed or stuck tracks in
    stderr", applied strictly: any stderr output invalidates); no stdout line outside the
    "Begin processing" progress lines mentions killed or stuck tracks; "Particle source Beam: Total
    number of histories: N" present with every N equal to the bundle's histories
    (NumberOfHistoriesInRun completed); the msc printout of GenericIon and alpha lists UrbanMsc
    only. The module-list check of the plan is NOT made here (see the evaluator docstring)."""
    problems = []
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
    return problems


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


def load_seal(repo: Path | str, seal_path: Path | str, seal_sha: str) -> tuple[dict[str, Any], str]:
    """The sealed sample record, byte-identical to the one committed at the seal commit."""
    blob = Path(seal_path).read_bytes()
    if blob != git_show(repo, seal_sha, SEAL_RECORD_PATH):
        raise RecordError(f"{seal_path} differs from {SEAL_RECORD_PATH} at the seal {seal_sha}")
    record = json.loads(blob.decode("utf-8"))
    if require(record, "schema", "sample seal") != SEAL_SCHEMA:
        raise RecordError(f"sample seal schema {record['schema']!r} != {SEAL_SCHEMA}")
    return record, hashlib.sha256(blob).hexdigest()


def case_chunks(record: dict[str, Any], case: int) -> list[dict[str, Any]]:
    """The sealed chunks of ``case``, checked against the sealed chunk plan (k = 1..n, the lv8
    seed formula, every chunk counter-clean, histories summing to 1e6)."""
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
        if c["seed"] != LV8_BASE + 100 * (R_INDEX_OFFSET + case) + c["k"]:
            raise RecordError(f"sample seal: case {case} k {c['k']} seed breaks the formula")
        if c.get("counters_all_zero") is not True:
            raise RecordError(f"sample seal: case {case} k {c['k']} is not counter-clean")
        if c["npz"] != f"ionmc-case{case}-k{c['k']}.npz":
            raise RecordError(f"sample seal: case {case} k {c['k']} npz name {c['npz']!r}")
        if c["histories"] % N_BATCHES:
            raise RecordError(f"sample seal: case {case} k {c['k']} histories not a multiple of 40")
    if sum(c["histories"] for c in chunks) != N_HISTORIES:
        raise RecordError(f"sample seal: case {case} histories do not sum to {N_HISTORIES}")
    return sorted(chunks, key=lambda c: c["k"])


def npz_path(samples_root: Path | str, case: int, chunk: dict[str, Any]) -> Path:
    return (
        Path(samples_root)
        / "validation"
        / "generated"
        / str(chunk["run_id"])
        / "r2-samples"
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
) -> IonmcGrid:
    """Load grid ``idd`` (-> (B, nz)), ``latx`` (-> (B, nx, len(slabs))) or ``laty`` (-> (B, ny,
    len(slabs))) of every sealed chunk of ``case`` (sha256 verified first) and combine per batch."""
    parts: list[tuple[int, np.ndarray]] = []
    edges: list[np.ndarray] | None = None
    for c in chunks:
        path = npz_path(samples_root, case, c)
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
