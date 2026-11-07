"""Provenance-bound measured calibration-part rate of the V7-R2 rep-based plans (plan R-4b, Amendment 20 (f)8;
Codex REVIEW-fe47f290 finding 2, REVIEW-61f4b0e3). Standard library only: ``run_suite.py`` (preflight) and
``steps_v5d.py`` (throughput gate) import it without ``ionmc``.

The rep-based V7-R2 plans (``v7r.r2_rep_plans``: 3200 repetitions over 39 parts per surrogate, the kappa parts) rest on
a slow-host seconds-per-repetition rate whose default (``v7r.R2_CAL_S_PER_REP`` = 12.2 s) is EXTRAPOLATED, and the 39-part
layout has 1.6 % headroom to the 1100 s planning target. The lv5d suite therefore accepts the rep plans only from a
MEASURED rate record, a JSON file named by the environment variable ``IONMC_V7R2_CAL_RATE_RECORD``::

    {"schema": "ionmc-v7r2-cal-rate-1",
     "s_per_rep": <wall_mc_s / n of the measured part>,   "fixed_s": <per-part fixed cost, >= the 70 s default>,
     "reps": <repetitions of the measured part (n)>,       "wall_mc_s": <its Monte Carlo wall time>,
     "wall_total_s": <pytest wall time of the part>,
     "host": {"node": ..., "platform": ..., "machine": ..., "cpu_count": ..., "python": ...},  # = the sealed part record's
     "host_run_id": <controlled host-runner run id, RUN-YYYYMMDDTHHMMSSZ-xxxxxxxx>,
     "host_run_token": <the token the wrapper set as IONMC_HOST_RUN_ID inside the run; = the sealed part record's>,
     "validation_id": <exact-SHA local validation record id, VAL-YYYYMMDD-HHMMSS-XXXXXX>,
     "validation_record": <file name (in log_dir) of the verbatim copy of that validation-manager record>,
     "validation_record_sha256": <64 hex>,
     "git_sha": <40 hex of the measured tree>,             "log_dir": <directory with the pytest log and part record>,
     "part_record": <file name of the sealed calibration part record>, "part_record_sha256": <64 hex>,
     "part_log": <file name of the pytest log>,            "part_log_sha256": <64 hex>}

``python v7r2_cal_rate.py --part-record ... --part-log ... --log-dir ... --git-sha ... --host-run-id RUN-...
--host-run-token ... --validation-id ... --validation-record ... --out ...`` builds the record from the sealed part
record (``wall_mc_s``, ``n``, ``host``, the sealed token), the pytest log (``1 passed in X s``) and the validation-record
copy; the rate is the Monte Carlo wall time per repetition and the fixed cost is ``max(70 s, wall_total - wall_mc)``
(conservative). No measurement is made here.

Binding (Codex REVIEW-ffe77aac). The record carries exactly the fields above and no other: every planning field is
DERIVED and recomputed by the loader from the verified artifacts, ``s_per_rep = wall_mc_s / reps`` and
``fixed_s = max(MIN_FIXED_S, wall_total_s - wall_mc_s)`` (stored values must agree, relative 1e-12); an undeclared
field is refused. ``host`` is bound to the SEALED part record: the part writer
(``tests/ionmc/test_v7r2_calibration.py``) seals a ``host`` block {node, platform, machine, cpu_count, python} of the
executing process and the token read from ``IONMC_HOST_RUN_ID`` (its sealed field is named ``host_run_id``; "unknown"
when unset, which is refused). The host-runner exposes no run id inside the process, so the wrapper sets that variable
to a token (``R4B-<timestamp>-<sha8>``); the controlled runner id (``RUN-...``) is a different string, and the
validation-manager record is what binds the two.

What is verified (fail closed; Codex REVIEW-6519d4ff, REVIEW-61f4b0e3). A record is only accepted, by
:func:`load_record` and before :func:`build_record` writes one, after the controlled-run artifacts have been re-opened
and checked (:func:`verify_artifacts`); every value of a built record is taken from them:

* ``log_dir`` (a relative path is resolved against the repository root) holds ``part_record``, ``part_log`` and
  ``validation_record``, and their sha256 equal ``part_record_sha256`` / ``part_log_sha256`` /
  ``validation_record_sha256``;
* the part record is a SEALED V7-R2 calibration part record: ``content_sha256`` equals the sha256 of its canonical JSON
  (sorted keys, compact separators) without that field (``steps_v5.content_digest``, re-implemented here); schema
  ``ionmc-v7r2-calibration-record-1``; the file name ``v7r2-calibration-<surrogate>-p<NN>.json`` agrees with its
  ``surrogate`` / ``part`` (39 parts per surrogate); layout ``{replicates 10800, reps_per_shard 225, parts 39}`` and
  ``reps_env`` 3200 (``v7r.r2_part_layout``); ``suite`` is ``lv5d``; ``run_sha`` and ``git_sha`` both equal the rate
  record's ``git_sha``; ``n`` is the planned size of that part (83 for parts 0-1, 82 otherwise, 3200 = 2 x 83 + 37 x
  82) and equals ``reps``; ``wall_mc_s`` equals the record's ``wall_mc_s``; its sealed token equals ``host_run_token``;
* the pytest log has exactly ONE pytest summary line, which reports exactly one passed test and no other outcome
  except deselected / warnings (no failed, error, skipped, xfailed, xpassed, rerun), no ``FAILED`` / ``ERROR`` line,
  and exactly one ``PARTIAL <part_record name> <content_sha256>`` line, which matches the part record (run the part
  with ``-s`` so the line is not captured); ``wall_total_s`` equals the "in X s" time of the summary (relative
  tolerance 1e-9: both are the same decimal text);
* the validation-record copy (the JSON the validation manager returned, copied verbatim) has ``schema_version``,
  ``validation_id`` = the record's, ``task_id`` "V3-005B", ``status`` "passed", ``git_sha`` (and ``details.head_sha``
  if present) = the record's ``git_sha``, a parseable ISO-8601 ``recorded_at_utc``, and ``details``: ``host_run_id`` =
  the record's runner id, ``host_run_token`` = the record's token (= the sealed one), ``part_record`` /
  ``part_record_sha256`` = the part record's name / SEAL digest (``content_sha256``, the digest of the ``PARTIAL`` line;
  the validation record names the seal, not the file bytes, whose sha256 is the rate record's ``part_record_sha256``), ``n`` = ``reps``, ``wall_mc_s`` (relative 1e-9),
  ``wall_total_s`` and ``fixed_s`` equal the record's, ``s_per_rep`` equal to the record's within the rounding of the
  copy (relative 1e-4). Time ordering against the part log is not checked (the log carries no wall-clock time).

A record that fails any of these raises ``ValueError`` naming the field.

What is NOT verifiable by repository code. The protected originals (the host-runner run record with its argv, exit
status and timestamps; the validation-manager store) cannot be read by repository code, and the unkeyed part seal
and the copy's JSON can be recomputed by anyone. The loader therefore cannot show that the copy is a faithful copy of
the protected record, nor that the run exited 0 on the slow host. That fidelity is attested by the reviewed commit
that adds the copy and by the validation id (the orchestrator copied the record verbatim and a reviewer can compare it
with ``inspect_local_validation``); the loader binds everything it CAN check (the copy's content, its digest, and its
agreement with the sealed artifacts) to that attestation through the plan pin.

Plan pin (Codex REVIEW-323f4c12, REVIEW-61f4b0e3). The checks above are recomputable by anyone holding synthetic
artifacts, so they do not by themselves establish the measurement. The controlled attestation is a PIN in the
committed acceptance plan ``validation/plans/v3-005-acceptance.md`` (a member of the hashed-source set of the suite,
changed only by reviewed commits). A measured rate record is CONFORMANT only if the plan carries a pin line (a whole
line, fixed grammar)::

    RATE-RECORD-OF-RECORD: sha256=<64 hex> git_sha=<40 hex> host_run_id=<RUN-...> validation_id=<VAL-...>
    host_run_token=<token> validation_record_sha256=<64 hex> host_node=<node>        (ONE line, single spaces)

The pin names the record (``sha256`` = :func:`record_digest`, canonical JSON), the controlled runner id, the validation
id, the in-run token and the sha256 of the validation-record copy, so a self-generated record with synthetic artifacts
is conformant only if a reviewed commit pins it together with a validation-record copy that binds it. :func:`load_record`
refuses the record unless the LAST pin line of the plan (a later pin supersedes every earlier one; a pin of an older
grammar is no pin) equals the record on every field, and the copy hashes to ``validation_record_sha256``. The plan
path is the repository plan, resolved from this module's location; another plan may be given ONLY as the ``plan_path``
function argument, never through the environment. The builder never writes the plan: ``--print-pin-line`` prints the
line to copy. Without a pin the record is nonconformant. Consumers (``steps_v5d.declared_plans``,
``rep_plan_table``, ``run_suite.rehearsal_preflight``) accept only the result of :func:`load_record`:
:func:`require_validated` re-loads the named file and demands an identical result, so a hand-made dictionary is
refused."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

CAL_RATE_ENV = "IONMC_V7R2_CAL_RATE_RECORD"
CAL_RATE_SCHEMA = "ionmc-v7r2-cal-rate-1"
MIN_FIXED_S = 70.0  # v7r.R2_CAL_FIXED_S (the planning value; the record may only be more conservative)
NOT_MEASURED = "not ok: calibration rate not measured"
_HEX40 = re.compile(r"[0-9a-f]{40}")
_HEX64 = re.compile(r"[0-9a-f]{64}")
_RUN_ID = re.compile(r"RUN-\d{8}T\d{6}Z-[0-9a-f]{8}")
_TOKEN = re.compile(r"[A-Za-z0-9_.:-]{8,}")
_VALIDATION_ID = re.compile(r"VAL-\d{8}-\d{6}-[0-9A-F]{6}")
HOST_FIELDS = ("node", "platform", "machine", "cpu_count", "python")
RECORD_FIELDS = frozenset((
    "schema", "s_per_rep", "fixed_s", "reps", "wall_mc_s", "wall_total_s", "host", "host_run_id", "host_run_token",
    "validation_id", "validation_record", "validation_record_sha256", "git_sha", "log_dir", "part_record",
    "part_record_sha256", "part_log", "part_log_sha256"))  # fmt: skip
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ACCEPTANCE_PLAN = REPO / "validation" / "plans" / "v3-005-acceptance.md"
PIN_PREFIX = "RATE-RECORD-OF-RECORD:"
_PIN = re.compile(
    r"^RATE-RECORD-OF-RECORD: sha256=([0-9a-f]{64}) git_sha=([0-9a-f]{40}) "
    r"host_run_id=(RUN-\d{8}T\d{6}Z-[0-9a-f]{8}) validation_id=(VAL-\d{8}-\d{6}-[0-9A-F]{6}) "
    r"host_run_token=([A-Za-z0-9_.:-]{8,}) validation_record_sha256=([0-9a-f]{64}) host_node=(\S+)[ \t]*$",
    re.MULTILINE)
PIN_KEYS = ("sha256", "git_sha", "host_run_id", "validation_id", "host_run_token", "validation_record_sha256",
            "host_node")  # the order of the pin line
NOT_PINNED = "not pinned in the acceptance plan (nonconformant until the orchestrator pins it)"
PART_SCHEMA = "ionmc-v7r2-calibration-record-1"  # v7r.CAL2_SCHEMA
PART_SUITE = "lv5d"
PART_LAYOUT = {"replicates": 10800, "reps_per_shard": 225, "parts": 39}  # v7r.r2_part_layout()
PART_REPS_ENV = 3200  # v7r.V7R_CAL_REPS
PART_SURROGATES = ("archive_fitted", "rehearsal_block_sums")
_PART_NAME = re.compile(r"v7r2-calibration-(archive_fitted|rehearsal_block_sums)-p(\d\d)\.json")
_SUMMARY = re.compile(r"^(?:=+ )?((?:\d+ [a-z]+(?:, )?)+) in ([0-9]+(?:\.[0-9]+)?)s(?: \(\d+:\d\d:\d\d\))?(?: =+)?$")
_COUNT = re.compile(r"(\d+) ([a-z]+)")
_PARTIAL = re.compile(r"^PARTIAL (\S+) ([0-9a-f]{64})$")


def _number(rec: dict[str, Any], key: str, *, positive: bool) -> float:
    v = rec.get(key)
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise ValueError(f"{key} is not a finite number")
    if (v <= 0.0) if positive else (v < 0.0):
        raise ValueError(f"{key} must be {'> 0' if positive else '>= 0'}")
    return float(v)


def validate_host(host: Any) -> dict[str, Any]:
    """The host block {node, platform, machine, cpu_count, python} (exactly these fields) or ``ValueError``."""
    if not (isinstance(host, dict) and set(host) == set(HOST_FIELDS)):
        raise ValueError(f"host must be an object with exactly the fields {list(HOST_FIELDS)}")
    for k in ("node", "platform", "machine", "python"):
        if not (isinstance(host[k], str) and host[k]):
            raise ValueError(f"host.{k} must be a non-empty string")
    cpu = host["cpu_count"]
    if isinstance(cpu, bool) or not isinstance(cpu, int) or cpu < 1:
        raise ValueError("host.cpu_count must be a positive integer")
    return host


def validate_token(token: Any) -> None:
    """Format of the in-run token (``IONMC_HOST_RUN_ID`` of the wrapper; "unknown" = unset is refused)."""
    if not (isinstance(token, str) and token != "unknown" and _TOKEN.fullmatch(token)):
        raise ValueError("host_run_token must be the in-run host run token (>= 8 chars of [A-Za-z0-9_.:-], "
                         "not 'unknown')")  # fmt: skip


def validate_ids(host_run_id: Any, validation_id: Any) -> None:
    """Format of the two protected-store ids (presence and format only: the stores are not readable here)."""
    if not (isinstance(host_run_id, str) and _RUN_ID.fullmatch(host_run_id)):
        raise ValueError("host_run_id must be a controlled host-runner run id RUN-YYYYMMDDTHHMMSSZ-xxxxxxxx")
    if not (isinstance(validation_id, str) and _VALIDATION_ID.fullmatch(validation_id)):
        raise ValueError("validation_id must be a validation record id VAL-YYYYMMDD-HHMMSS-XXXXXX")


def validate_record(rec: Any) -> dict[str, Any]:
    """The record itself or ``ValueError`` (every field of the module docstring is required and checked)."""
    if not isinstance(rec, dict):
        raise ValueError("the record is not a JSON object")
    if rec.get("schema") != CAL_RATE_SCHEMA:
        raise ValueError(f"schema is {rec.get('schema')!r}, expected {CAL_RATE_SCHEMA!r}")
    s = _number(rec, "s_per_rep", positive=True)
    fixed = _number(rec, "fixed_s", positive=False)
    wall = _number(rec, "wall_mc_s", positive=True)
    _number(rec, "wall_total_s", positive=True)
    reps = rec.get("reps")
    if isinstance(reps, bool) or not isinstance(reps, int) or reps < 1:
        raise ValueError("reps must be a positive integer")
    extra = set(rec) - RECORD_FIELDS
    if extra or RECORD_FIELDS - set(rec):
        raise ValueError(f"record fields differ from the declared set (undeclared {sorted(extra)}, "
                         f"missing {sorted(RECORD_FIELDS - set(rec))})")  # fmt: skip
    if fixed < MIN_FIXED_S:
        raise ValueError(f"fixed_s {fixed} is below the planning value {MIN_FIXED_S}")
    validate_host(rec.get("host"))
    validate_ids(rec.get("host_run_id"), rec.get("validation_id"))
    validate_token(rec.get("host_run_token"))
    if not (isinstance(rec.get("git_sha"), str) and _HEX40.fullmatch(rec["git_sha"])):
        raise ValueError("git_sha must be 40 lowercase hex digits")
    if not (isinstance(rec.get("log_dir"), str) and rec["log_dir"]):
        raise ValueError("log_dir must be a non-empty string")
    for key in ("part_record", "part_log", "validation_record"):
        if not (isinstance(rec.get(key), str) and rec[key]):
            raise ValueError(f"{key} must be a non-empty file name")
        if not (isinstance(rec.get(f"{key}_sha256"), str) and _HEX64.fullmatch(rec[f"{key}_sha256"])):
            raise ValueError(f"{key}_sha256 must be 64 lowercase hex digits")
    return rec


def _sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _resolve_dir(log_dir: str) -> Path:
    p = Path(log_dir)
    return p if p.is_absolute() else REPO / p


def part_size(part: int) -> int:
    """Planned repetitions of calibration part ``part`` (``v7r.cal_part_sizes(3200, 39)``)."""
    q, r = divmod(PART_REPS_ENV, PART_LAYOUT["parts"])
    return q + (1 if part < r else 0)


def seal_digest(doc: dict[str, Any]) -> str:
    """sha256 of the canonical JSON (sorted keys, compact separators) of ``doc`` without ``content_sha256``, after one
    JSON round trip (``steps_v5.content_digest``)."""
    norm = json.loads(json.dumps(doc))
    norm.pop("content_sha256", None)
    return hashlib.sha256(json.dumps(norm, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def verify_part_record(path: Path, git_sha: str) -> dict[str, Any]:
    """The sealed V7-R2 calibration part record at ``path`` (rules of the module docstring) or ``ValueError``."""
    try:
        doc = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f"part record {path.name} unreadable: {exc}") from exc
    if not isinstance(doc, dict):
        raise ValueError(f"part record {path.name} is not a JSON object")
    if doc.get("content_sha256") != seal_digest(doc):
        raise ValueError(f"part record {path.name} is not sealed: content_sha256 does not match the document")
    if doc.get("schema") != PART_SCHEMA:
        raise ValueError(f"part record {path.name}: schema {doc.get('schema')!r}, expected {PART_SCHEMA!r}")
    m = _PART_NAME.fullmatch(path.name)
    part = doc.get("part")
    if not m or isinstance(part, bool) or doc.get("surrogate") != m.group(1) or part != int(m.group(2)):
        raise ValueError(f"part record {path.name}: file name and surrogate / part of the record disagree")
    if not 0 <= part < PART_LAYOUT["parts"]:
        raise ValueError(f"part record {path.name}: part {part} outside 0..{PART_LAYOUT['parts'] - 1}")
    if doc.get("layout") != PART_LAYOUT or doc.get("reps_env") != PART_REPS_ENV:
        raise ValueError(f"part record {path.name}: layout / reps_env are not the lv5d V7-R2 calibration layout")
    if doc.get("suite") != PART_SUITE:
        raise ValueError(f"part record {path.name}: suite {doc.get('suite')!r}, expected {PART_SUITE!r}")
    if doc.get("run_sha") != git_sha or doc.get("git_sha", git_sha) != git_sha:
        raise ValueError(f"part record {path.name}: run_sha / git_sha differ from the record's git_sha {git_sha}")
    n, wall = doc.get("n"), doc.get("wall_mc_s")
    if isinstance(n, bool) or n != part_size(part):
        raise ValueError(f"part record {path.name}: n {n!r} is not the planned size {part_size(part)} of part {part}")
    if isinstance(wall, bool) or not isinstance(wall, (int, float)) or not math.isfinite(wall) or wall <= 0.0:
        raise ValueError(f"part record {path.name}: wall_mc_s is not a positive finite number")
    try:
        validate_host(doc.get("host"))
        validate_token(doc.get("host_run_id"))  # the sealed field of the part record holds the in-run token
    except ValueError as exc:
        raise ValueError(f"part record {path.name}: sealed host identity unusable: {exc}") from exc
    return doc


def verify_log(text: str, part_name: str, part_digest: str) -> float:
    """The pytest wall time ("in X s") of a log of exactly one passing calibration part, or ``ValueError``."""
    lines = [ln.strip() for ln in text.splitlines()]
    if any(ln.startswith(("FAILED ", "ERROR ")) for ln in lines):
        raise ValueError("the pytest log contains a FAILED / ERROR line")
    summaries = [m for ln in lines if (m := _SUMMARY.match(ln))]
    if len(summaries) != 1:
        raise ValueError(f"the pytest log has {len(summaries)} summary lines (need exactly 1)")
    counts: dict[str, int] = {}
    for num, word in _COUNT.findall(summaries[0].group(1)):
        counts[word] = counts.get(word, 0) + int(num)
    extra = set(counts) - {"passed", "deselected", "warning", "warnings"}
    if counts.get("passed") != 1 or extra:
        raise ValueError(f"the pytest summary is not exactly one passed test ({summaries[0].group(1)!r})")
    partials = [m.groups() for ln in lines if (m := _PARTIAL.match(ln))]
    if partials != [(part_name, part_digest)]:
        raise ValueError(f"the pytest log needs exactly one 'PARTIAL {part_name} {part_digest}' line, found {partials}")
    return float(summaries[0].group(2))


def _close(a: Any, b: Any, rel: float) -> bool:
    return (isinstance(a, (int, float)) and not isinstance(a, bool) and isinstance(b, (int, float))
            and math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=rel))  # fmt: skip


def verify_validation_record(path: Path, rec: dict[str, Any], part_seal: str) -> None:
    """The verbatim validation-manager record copy at ``path`` against the (shape-validated) rate record ``rec`` and
    the seal ``part_seal`` (``content_sha256``) of the part record (module docstring); ``ValueError`` naming the first field that disagrees."""
    try:
        doc = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f"validation record {path.name} unreadable: {exc}") from exc
    if not isinstance(doc, dict):
        raise ValueError(f"validation record {path.name} is not a JSON object")
    if "schema_version" not in doc:
        raise ValueError(f"validation record {path.name}: schema_version missing")
    for key, want in (("validation_id", rec["validation_id"]), ("task_id", "V3-005B"), ("status", "passed"),
                      ("git_sha", rec["git_sha"])):  # fmt: skip
        if doc.get(key) != want:
            raise ValueError(f"validation record {path.name}: {key} {doc.get(key)!r} differs from the expected "
                             f"{want!r}")  # fmt: skip
    at = doc.get("recorded_at_utc")
    try:
        if not isinstance(at, str):
            raise ValueError("not a string")
        datetime.fromisoformat(at)
    except ValueError as exc:
        raise ValueError(f"validation record {path.name}: recorded_at_utc {at!r} is not ISO-8601") from exc
    det = doc.get("details")
    if not isinstance(det, dict):
        raise ValueError(f"validation record {path.name}: details is not an object")
    for key, want in (("host_run_id", rec["host_run_id"]), ("host_run_token", rec["host_run_token"]),
                      ("part_record", rec["part_record"]), ("part_record_sha256", part_seal),
                      ("n", rec["reps"]), ("fixed_s", rec["fixed_s"])):  # fmt: skip
        if det.get(key) != want or isinstance(det.get(key), bool):
            raise ValueError(f"validation record {path.name}: details.{key} {det.get(key)!r} differs from {want!r}")
    if "head_sha" in det and det["head_sha"] != rec["git_sha"]:
        raise ValueError(f"validation record {path.name}: details.head_sha {det['head_sha']!r} differs from the "
                         f"record's git_sha")  # fmt: skip
    for key, rel in (("wall_mc_s", 1e-9), ("wall_total_s", 1e-9), ("s_per_rep", 1e-4)):
        if not _close(det.get(key), rec[key], rel):
            raise ValueError(f"validation record {path.name}: details.{key} {det.get(key)!r} differs from "
                             f"{rec[key]!r}")  # fmt: skip


def verify_artifacts(rec: dict[str, Any]) -> None:
    """Re-open the part record, the pytest log and the validation-record copy named by the (shape-validated) ``rec``
    and check them against it (module docstring); ``ValueError`` on the first mismatch."""
    d = _resolve_dir(rec["log_dir"])
    files = {}
    for key in ("part_record", "part_log", "validation_record"):
        name = rec[key]
        if Path(name).name != name:
            raise ValueError(f"{key} {name!r} is not a plain file name")
        p = d / name
        if not p.is_file():
            raise ValueError(f"{key} {p} is not a file")
        if _sha256(p) != rec[f"{key}_sha256"]:
            raise ValueError(f"{key} {p.name}: sha256 differs from the record")
        files[key] = p
    doc = verify_part_record(files["part_record"], rec["git_sha"])
    if doc["n"] != rec["reps"]:
        raise ValueError(f"reps {rec['reps']} differs from n {doc['n']} of the part record")
    if doc["wall_mc_s"] != rec["wall_mc_s"]:
        raise ValueError(f"wall_mc_s {rec['wall_mc_s']} differs from the part record's {doc['wall_mc_s']}")
    if rec["host"] != doc["host"]:
        raise ValueError(f"host {rec['host']} differs from the sealed part record's {doc['host']}")
    if rec["host_run_token"] != doc["host_run_id"]:
        raise ValueError(f"host_run_token {rec['host_run_token']!r} differs from the sealed part record's "
                         f"{doc['host_run_id']!r}")  # fmt: skip
    total = verify_log(files["part_log"].read_text(errors="replace"), files["part_record"].name,
                       doc["content_sha256"])  # fmt: skip
    if not math.isclose(total, rec["wall_total_s"], rel_tol=1e-9):
        raise ValueError(f"wall_total_s {rec['wall_total_s']} differs from the pytest summary time {total}")
    # derived planning fields, recomputed from the verified artifacts (Codex REVIEW-ffe77aac)
    want_s = doc["wall_mc_s"] / doc["n"]
    if not math.isclose(rec["s_per_rep"], want_s, rel_tol=1e-12):
        raise ValueError(f"s_per_rep {rec['s_per_rep']} differs from wall_mc_s / reps = {want_s}")
    want_fixed = max(MIN_FIXED_S, total - doc["wall_mc_s"])
    if not math.isclose(rec["fixed_s"], want_fixed, rel_tol=1e-12):
        raise ValueError(f"fixed_s {rec['fixed_s']} differs from max({MIN_FIXED_S}, wall_total_s - wall_mc_s) = "
                         f"{want_fixed}")  # fmt: skip
    verify_validation_record(files["validation_record"], rec, doc["content_sha256"])  # Codex REVIEW-61f4b0e3


def record_digest(rec: dict[str, Any]) -> str:
    """sha256 of the canonical JSON (sorted keys, compact separators) of the rate record: the digest a pin names."""
    return hashlib.sha256(json.dumps(rec, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _pin_fields(rec: dict[str, Any]) -> dict[str, str]:
    """The pin fields of the (validated) record, in the order of the pin line."""
    return {"sha256": record_digest(rec), "git_sha": rec["git_sha"], "host_run_id": rec["host_run_id"],
            "validation_id": rec["validation_id"], "host_run_token": rec["host_run_token"],
            "validation_record_sha256": rec["validation_record_sha256"], "host_node": rec["host"]["node"]}  # fmt: skip


def pin_line(rec: dict[str, Any]) -> str:
    """The exact pin line of the (validated) record, as the orchestrator writes it into the acceptance plan."""
    return PIN_PREFIX + "".join(f" {k}={v}" for k, v in _pin_fields(rec).items())


def plan_pin(plan_path: str | os.PathLike[str] | None = None) -> dict[str, str] | None:
    """The LAST pin line of the plan (default: the repository acceptance plan) as a dict, or ``None`` if there is none.
    A pin is a whole line of the exact grammar of the module docstring; a placeholder with angle brackets, and a pin of
    the older grammar (without the runner id / token / validation-record digest), is no pin."""
    p = Path(plan_path) if plan_path is not None else ACCEPTANCE_PLAN
    try:
        text = p.read_text(encoding="utf-8")
    except (OSError, ValueError) as exc:
        raise ValueError(f"acceptance plan {p} unreadable: {exc}") from exc
    pins = _PIN.findall(text)
    if not pins:
        return None
    return dict(zip(PIN_KEYS, pins[-1], strict=True))


def verify_pin(rec: dict[str, Any], plan_path: str | os.PathLike[str] | None = None) -> None:
    """``ValueError`` unless the last pin of the plan names exactly this record (digest, runner id, validation id,
    token, validation-record digest, git sha and host node)."""
    pin = plan_pin(plan_path)
    if pin is None:
        raise ValueError(f"the record is {NOT_PINNED}")
    bad = [k for k, v in _pin_fields(rec).items() if pin[k] != v]
    if bad:
        raise ValueError(f"the last pin of the acceptance plan does not match the record ({', '.join(bad)} differ); "
                         f"the record is {NOT_PINNED}")  # fmt: skip


def _load_verified(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Shape and artifact verification only (no pin): ``{"record", "sha256", "path"}`` or ``ValueError``."""
    p = Path(path)
    if not p.is_file():
        raise ValueError(f"{p} is not a file")
    raw = p.read_bytes()
    try:
        rec = validate_record(json.loads(raw))
        verify_artifacts(rec)
    except ValueError as exc:  # JSON errors are ValueErrors as well
        raise ValueError(f"calibration rate record {p} invalid: {exc}") from exc
    return {"record": rec, "sha256": hashlib.sha256(raw).hexdigest(), "path": str(p)}


def load_record(path: str | os.PathLike[str], plan_path: str | os.PathLike[str] | None = None) -> dict[str, Any]:
    """``{"record", "sha256", "path"}`` of the rate record file, or ``ValueError`` (absent, unreadable, malformed,
    whose controlled-run artifacts fail :func:`verify_artifacts`, or which is not pinned in the acceptance plan: the
    measured rate is conformant only when pinned). ``plan_path`` (a function argument only) replaces the repository
    plan."""
    got = _load_verified(path)
    try:
        verify_pin(got["record"], plan_path)
    except ValueError as exc:
        raise ValueError(f"calibration rate record {got['path']} invalid: {exc}") from exc
    return got


def record_from_env() -> dict[str, Any] | None:
    """:func:`load_record` (repository plan) of ``IONMC_V7R2_CAL_RATE_RECORD`` (``None`` if unset or empty)."""
    path = os.environ.get(CAL_RATE_ENV, "")
    return load_record(path) if path else None


def build_record(part_record: Path, part_log: Path, log_dir: Path, git_sha: str, host_run_id: str,
                 host_run_token: str, validation_id: str, validation_record: Path) -> dict[str, Any]:  # fmt: skip
    """The record of one measured calibration part (see the module docstring). Every value comes from the verified
    sealed part record, pytest log and validation-record copy (all in ``log_dir``); ``host`` and the token are copied
    from the sealed record (a differing ``host_run_token`` argument is refused); the ids are required; the result
    passes :func:`verify_artifacts` (so the copy's bindings must agree with the artifacts) or ``ValueError`` is
    raised."""
    validate_ids(host_run_id, validation_id)
    validate_token(host_run_token)
    if not (isinstance(git_sha, str) and _HEX40.fullmatch(git_sha)):
        raise ValueError("git_sha must be 40 lowercase hex digits")
    if not (part_record.is_file() and part_log.is_file() and validation_record.is_file()):
        raise ValueError("the part record, the pytest log and the validation record must be files")
    ld = log_dir if log_dir.is_absolute() else Path.cwd() / log_dir
    if any(f.resolve().parent != ld.resolve() for f in (part_record, part_log, validation_record)):
        raise ValueError("the part record, the pytest log and the validation record must lie in log_dir")
    doc = verify_part_record(part_record, git_sha)
    if doc["host_run_id"] != host_run_token:
        raise ValueError(f"host_run_token {host_run_token!r} differs from the sealed part record's "
                         f"{doc['host_run_id']!r}")  # fmt: skip
    total = verify_log(part_log.read_text(errors="replace"), part_record.name, doc["content_sha256"])
    n, wall = doc["n"], float(doc["wall_mc_s"])
    try:
        stored = str(ld.resolve().relative_to(REPO))
    except ValueError:
        stored = str(ld.resolve())
    rec = {
        "schema": CAL_RATE_SCHEMA, "s_per_rep": wall / n, "fixed_s": max(MIN_FIXED_S, total - wall),
        "reps": n, "wall_mc_s": wall, "wall_total_s": total,
        "host": {k: doc["host"][k] for k in HOST_FIELDS}, "host_run_id": host_run_id,
        "host_run_token": doc["host_run_id"], "validation_id": validation_id,
        "validation_record": validation_record.name, "validation_record_sha256": _sha256(validation_record),
        "git_sha": git_sha, "log_dir": stored,
        "part_record": part_record.name, "part_record_sha256": _sha256(part_record),
        "part_log": part_log.name, "part_log_sha256": _sha256(part_log),
    }  # fmt: skip
    validate_record(rec)
    verify_artifacts(rec)
    return rec


def require_validated(got: Any, plan_path: str | os.PathLike[str] | None = None) -> dict[str, Any]:
    """``got`` if it is exactly the result of :func:`load_record` of the file it names, else ``ValueError``: the file
    is re-loaded (all artifact checks and the plan pin again, against ``plan_path`` or the repository plan) and the
    result must be identical; a hand-made dictionary, a record whose file changed or whose pin was removed, or
    ``None`` is refused."""
    if not (isinstance(got, dict) and isinstance(got.get("path"), str) and isinstance(got.get("sha256"), str)):
        raise ValueError("not a calibration rate record loaded by v7r2_cal_rate.load_record")
    fresh = load_record(got["path"], plan_path)
    if fresh != got:
        raise ValueError("the calibration rate record differs from the file it names (not the result of load_record)")
    return fresh


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument(
        "--part-record",
        required=True,
        help="sealed v7r2-calibration-<surrogate>-p<NN>.json of the part",
    )
    ap.add_argument("--part-log", required=True, help="the logged pytest output of that single part")
    ap.add_argument("--log-dir", required=True, help="directory holding both files (recorded)")
    ap.add_argument("--git-sha", required=True, help="full SHA of the tree the part was measured at")
    ap.add_argument("--host-run-id", required=True, help="controlled host-runner run id RUN-... of the measurement")
    ap.add_argument("--host-run-token", required=True,
                    help="the in-run token (IONMC_HOST_RUN_ID of the wrapper); must equal the sealed part record's")  # fmt: skip
    ap.add_argument(
        "--validation-id",
        required=True,
        help="exact-SHA validation record id (VAL-YYYYMMDD-HHMMSS-XXXXXX)",
    )
    ap.add_argument("--validation-record", required=True,
                    help="verbatim copy of the validation-manager record (inside --log-dir)")  # fmt: skip
    ap.add_argument("--out", required=True, help="rate record to write (refuses to overwrite)")
    ap.add_argument("--print-pin-line", action="store_true",
                    help="also print the pin line of the record for the acceptance plan (the orchestrator writes it "
                         "after verifying the protected records; this tool never writes the plan)")  # fmt: skip
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists")
    try:
        rec = build_record(
            Path(a.part_record),
            Path(a.part_log),
            Path(a.log_dir),
            a.git_sha,
            a.host_run_id,
            a.host_run_token,
            a.validation_id,
            Path(a.validation_record),
        )
    except ValueError as exc:
        raise SystemExit(f"not ok: {exc}") from exc
    out.write_text(json.dumps(rec, sort_keys=True, indent=1) + "\n")
    print(f"{out}: s_per_rep {rec['s_per_rep']:.4f}, fixed_s {rec['fixed_s']:.1f}, "
          f"sha256 {_load_verified(out)['sha256']}")  # fmt: skip
    if a.print_pin_line:
        print(pin_line(rec))
    return 0


if __name__ == "__main__":
    sys.exit(main())
