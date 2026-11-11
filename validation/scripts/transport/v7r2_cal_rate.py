"""Provenance-bound measured calibration-part rate of the V7-R2 rep-based plans (plan R-4b, Amendment 20 (f)8;
Codex REVIEW-fe47f290 finding 2, REVIEW-61f4b0e3, REVIEW-cf3f019e). Standard library plus the repository verifier
``attest_verify`` (same directory) only: ``run_suite.py`` (preflight) and ``steps_v5d.py`` (throughput gate) import it
without ``ionmc``.

The rep-based V7-R2 plans (``v7r.r2_rep_plans``: 3200 repetitions over 39 parts per surrogate, the kappa parts) rest on
a slow-host seconds-per-repetition rate whose default (``v7r.R2_CAL_S_PER_REP`` = 12.2 s) is EXTRAPOLATED, and the 39-part
layout has 1.6 % headroom to the 1100 s planning target. The lv5d suite therefore accepts the rep plans only from a
MEASURED rate record, a JSON file named by the environment variable ``IONMC_V7R2_CAL_RATE_RECORD``::

    {"schema": "ionmc-v7r2-cal-rate-2",
     "s_per_rep": <wall_mc_s / n of the measured part>,   "fixed_s": <per-part fixed cost, >= the 70 s default>,
     "reps": <repetitions of the measured part (n)>,       "wall_mc_s": <its Monte Carlo wall time>,
     "wall_total_s": <pytest wall time of the part>,
     "host": {"node": ..., "platform": ..., "machine": ..., "cpu_count": ..., "python": ...},  # = the sealed part record's
     "git_sha": <40 hex of the measured tree>,             "log_dir": <directory holding the four files below>,
     "part_record": <file name of the sealed calibration part record>, "part_record_sha256": <64 hex>,
     "part_log": <file name of the pytest log of the part>, "part_log_sha256": <64 hex>,
     "host_run_id": <RUN-YYYYMMDDTHHMMSSZ-xxxxxxxx>,       "host_run_record": <file name of the signed RUN- copy>,
     "host_run_record_sha256": <64 hex>,                   "host_run_signature": {"key_id", "payload_sha256", "signature"},
     "validation_id": <VAL-YYYYMMDD-HHMMSS-XXXXXX>,        "validation_record": <file name of the signed VAL- copy>,
     "validation_record_sha256": <64 hex>,                 "validation_signature": {"key_id", "payload_sha256", "signature"},
     "trust_file_sha256": <sha256 of the trust file the signatures were verified against>}

``python v7r2_cal_rate.py --part-record ... --part-log ... --log-dir ... --git-sha ... --host-run-record RUN-....json
--validation-record VAL-....json [--trust-file ...] --out ...`` builds the record; the rate is the Monte Carlo wall
time per repetition and the fixed cost is ``max(70 s, wall_total - wall_mc)`` (conservative). No measurement is made
here. Schema 1 (the token-attested record of 2026-11-07, ``host_run_token`` and a repository-authored validation copy;
rejected by Codex REVIEW-61f4b0e3 / REVIEW-cf3f019e) is refused with an explicit "legacy" message.

Signed attestation chain (V3-INFRA-002/003; ``experiment/v3/SETUP.md``, "Signed attestations"). The two protected
services sign their records with Ed25519 keys whose public halves are committed in the trust file
``infrastructure/keys/trusted_attestation_keys.json`` (:data:`TRUST_FILE`; a test passes another file). The rate record
names VERBATIM copies of the two signed records, both in ``log_dir``; every load re-verifies them (fail closed):

* the host-runner record copy (``RUN-...``): ``attest_verify.require_verified`` with signer ``host-runner``, task
  ``V3-005B``, git SHA = the record's ``git_sha``, run id = its own ``run_id`` (= the record's ``host_run_id``) and
  output manifest = the ``output_manifest_sha256`` the validation record claims for that run (64 hex, never null);
  status ``verified`` only (untrusted key = unverified = refused; altered payload or bad signature = invalid =
  refused); and the run succeeded: ``succeeded`` true, ``exit_code`` 0, ``timed_out`` false,
  ``committed_state_unchanged`` true, ``dirty_after`` false (all present), ``output_publication_error`` null;
* the validation-manager record copy (``VAL-...``): ``require_verified`` with signer ``validation-manager``, task
  ``V3-005B`` and the record's ``git_sha``; ``validation_id`` = the record's; ``status`` ``passed``;
  ``details.host_runs`` is exactly ``[{"run_id": <host_run_id>, "output_manifest_sha256": <64 hex>}]``, with neither
  the legacy ``details.host_run_id`` / ``host_run_ids`` nor ``host_run_manifest_lookup`` "missing";
* the record's ``host_run_signature`` / ``validation_signature`` equal the copies' attestation ``key_id``,
  ``payload_sha256`` and ``signature``, ``host_run_record_sha256`` / ``validation_record_sha256`` the copies' file
  sha256, and ``trust_file_sha256`` the sha256 of the trust file used for the verification (a changed trust file
  invalidates the record until it is rebuilt and re-pinned).

Binding of the measurement to the signed run. The calibration part record and its pytest log are unsigned outputs;
they are bound to the run through its SIGNED stdout (the ``stdout`` field of the RUN record, the inline text the
signature covers; when it is not truncated, its sha256 must equal the signed ``stdout_sha256``):

* the in-run token sealed in the part record (field ``host_run_id``; the wrapper sets ``IONMC_HOST_RUN_ID`` =
  ``R4B-<timestamp>-<sha8>``, the runner exposes no run id inside the process) appears in an element of the signed
  ``argv`` or as a whole stdout line ``token <token>`` (:func:`token_line`);
* the signed stdout carries exactly one line ``CAL-RATE-ARTIFACT part_record <name> <sha256>`` and exactly one line
  ``CAL-RATE-ARTIFACT part_log <name> <sha256>`` (:func:`artifact_lines`, printed by the wrapper after the part, i.e.
  at the END of stdout, which survives the runner's tail truncation of long streams), naming exactly the record's
  ``part_record`` / ``part_log`` and their file sha256.

Beyond that, every value of the record is derived from the verified artifacts (Codex REVIEW-ffe77aac, REVIEW-6519d4ff):

* the part record is a SEALED V7-R2 calibration part record: ``content_sha256`` equals the sha256 of its canonical JSON
  (sorted keys, compact separators) without that field (``steps_v5.content_digest``, re-implemented here); schema
  ``ionmc-v7r2-calibration-record-1``; the file name ``v7r2-calibration-<surrogate>-p<NN>.json`` agrees with its
  ``surrogate`` / ``part`` (39 parts per surrogate); layout ``{replicates 10800, reps_per_shard 225, parts 39}`` and
  ``reps_env`` 3200 (``v7r.r2_part_layout``); ``suite`` is ``lv5d``; ``run_sha`` (and ``git_sha`` if present) equal the
  rate record's ``git_sha``; ``n`` is the planned size of that part (83 for parts 0-1, 82 otherwise) and equals
  ``reps``; ``wall_mc_s`` equals the record's ``wall_mc_s``; ``host`` equals the record's ``host``;
* the pytest log has exactly ONE pytest summary line, which reports exactly one passed test and no other outcome
  except deselected / warnings, no ``FAILED`` / ``ERROR`` line, and exactly one ``PARTIAL <part_record name>
  <content_sha256>`` line, which matches the part record (run the part with ``-s``); ``wall_total_s`` equals the
  "in X s" time of the summary (relative 1e-9);
* ``s_per_rep = wall_mc_s / reps`` and ``fixed_s = max(MIN_FIXED_S, wall_total_s - wall_mc_s)`` (relative 1e-12); the
  record carries exactly the declared fields (an undeclared field is refused).

A record that fails any of these raises ``ValueError`` naming the field.

Plan pin (Codex REVIEW-323f4c12, REVIEW-61f4b0e3). A measured rate record is CONFORMANT only if the committed acceptance
plan ``validation/plans/v3-005-acceptance.md`` carries a pin line (a whole line, fixed grammar, single spaces)::

    RATE-RECORD-OF-RECORD: sha256=<64 hex> git_sha=<40 hex> host_run_id=<RUN-...> validation_id=<VAL-...>
    host_run_record_sha256=<64 hex> validation_record_sha256=<64 hex> trust_file_sha256=<64 hex>      (ONE line)

``sha256`` is :func:`record_digest` (canonical JSON of the record). :func:`load_record` refuses the record unless the
LAST pin line of the plan (a later pin supersedes every earlier one; a pin of an older grammar, e.g. with
``host_run_token`` / ``host_node``, is no pin) equals the record on every field. The plan path is the repository plan,
resolved from this module's location; another plan may be given ONLY as the ``plan_path`` function argument, never
through the environment. The builder never writes the plan: ``--print-pin-line`` prints the line to copy. Consumers
(``steps_v5d.declared_plans``, ``rep_plan_table``, ``run_suite.rehearsal_preflight``) accept only the result of
:func:`load_record`: :func:`require_validated` re-loads the named file (signatures included) and demands an identical
result, so a hand-made dictionary is refused."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any

try:
    import attest_verify
except ImportError:  # imported by path (tests) without the scripts directory on sys.path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import attest_verify

CAL_RATE_ENV = "IONMC_V7R2_CAL_RATE_RECORD"
CAL_RATE_SCHEMA = "ionmc-v7r2-cal-rate-2"
LEGACY_SCHEMAS = ("ionmc-v7r2-cal-rate-1",)
LEGACY_REFUSED = ("legacy token-attested calibration rate record (schema ionmc-v7r2-cal-rate-1 / host_run_token; "
                  "Codex REVIEW-61f4b0e3, REVIEW-cf3f019e) refused: rebuild it from the signed RUN- and VAL- record "
                  "copies")  # fmt: skip
TASK_ID = "V3-005B"
RUN_SIGNER = "host-runner"
VAL_SIGNER = "validation-manager"
MIN_FIXED_S = 70.0  # v7r.R2_CAL_FIXED_S (the planning value; the record may only be more conservative)
NOT_MEASURED = "not ok: calibration rate not measured"
_HEX40 = re.compile(r"[0-9a-f]{40}")
_HEX64 = re.compile(r"[0-9a-f]{64}")
_KEY_ID = re.compile(r"[0-9a-f]{16}")
_SIG = re.compile(r"[0-9a-f]{128}")
_RUN_ID = re.compile(r"RUN-\d{8}T\d{6}Z-[0-9a-f]{8}")
_TOKEN = re.compile(r"[A-Za-z0-9_.:-]{8,}")
_VALIDATION_ID = re.compile(r"VAL-\d{8}-\d{6}-[0-9A-F]{6}")
HOST_FIELDS = ("node", "platform", "machine", "cpu_count", "python")
SIGNATURE_FIELDS = ("key_id", "payload_sha256", "signature")
FILE_KEYS = ("part_record", "part_log", "host_run_record", "validation_record")
RECORD_FIELDS = frozenset((
    "schema", "s_per_rep", "fixed_s", "reps", "wall_mc_s", "wall_total_s", "host", "git_sha", "log_dir",
    *FILE_KEYS, *(f"{k}_sha256" for k in FILE_KEYS), "host_run_id", "host_run_signature", "validation_id",
    "validation_signature", "trust_file_sha256"))  # fmt: skip
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ACCEPTANCE_PLAN = REPO / "validation" / "plans" / "v3-005-acceptance.md"
TRUST_FILE = Path(attest_verify.TRUSTED_KEYS_PATH)
PIN_PREFIX = "RATE-RECORD-OF-RECORD:"
_PIN = re.compile(
    r"^RATE-RECORD-OF-RECORD: sha256=([0-9a-f]{64}) git_sha=([0-9a-f]{40}) "
    r"host_run_id=(RUN-\d{8}T\d{6}Z-[0-9a-f]{8}) validation_id=(VAL-\d{8}-\d{6}-[0-9A-F]{6}) "
    r"host_run_record_sha256=([0-9a-f]{64}) validation_record_sha256=([0-9a-f]{64}) "
    r"trust_file_sha256=([0-9a-f]{64})[ \t]*$",
    re.MULTILINE)
PIN_KEYS = ("sha256", "git_sha", "host_run_id", "validation_id", "host_run_record_sha256",
            "validation_record_sha256", "trust_file_sha256")  # the order of the pin line
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
ARTIFACT_PREFIX = "CAL-RATE-ARTIFACT"
_ARTIFACT = re.compile(r"^CAL-RATE-ARTIFACT (part_record|part_log) (\S+) ([0-9a-f]{64})$")
TOKEN_PREFIX = "token"


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
    """Format of the in-run token sealed in the part record (``IONMC_HOST_RUN_ID`` of the wrapper; "unknown" = unset
    is refused)."""
    if not (isinstance(token, str) and token != "unknown" and _TOKEN.fullmatch(token)):
        raise ValueError("the sealed in-run token (part record host_run_id) must be >= 8 chars of [A-Za-z0-9_.:-], "
                         "not 'unknown'")  # fmt: skip


def validate_ids(host_run_id: Any, validation_id: Any) -> None:
    """Format of the two protected-store ids."""
    if not (isinstance(host_run_id, str) and _RUN_ID.fullmatch(host_run_id)):
        raise ValueError("host_run_id must be a controlled host-runner run id RUN-YYYYMMDDTHHMMSSZ-xxxxxxxx")
    if not (isinstance(validation_id, str) and _VALIDATION_ID.fullmatch(validation_id)):
        raise ValueError("validation_id must be a validation record id VAL-YYYYMMDD-HHMMSS-XXXXXX")


def _validate_signature(sig: Any, key: str) -> None:
    if not (isinstance(sig, dict) and set(sig) == set(SIGNATURE_FIELDS)):
        raise ValueError(f"{key} must be an object with exactly the fields {list(SIGNATURE_FIELDS)}")
    for field, pat in (("key_id", _KEY_ID), ("payload_sha256", _HEX64), ("signature", _SIG)):
        if not (isinstance(sig[field], str) and pat.fullmatch(sig[field])):
            raise ValueError(f"{key}.{field} is malformed")


def validate_record(rec: Any) -> dict[str, Any]:
    """The record itself or ``ValueError`` (every field of the module docstring is required and checked; a legacy
    schema-1 record is refused with :data:`LEGACY_REFUSED`)."""
    if not isinstance(rec, dict):
        raise ValueError("the record is not a JSON object")
    if rec.get("schema") in LEGACY_SCHEMAS or "host_run_token" in rec:
        raise ValueError(LEGACY_REFUSED)
    if rec.get("schema") != CAL_RATE_SCHEMA:
        raise ValueError(f"schema is {rec.get('schema')!r}, expected {CAL_RATE_SCHEMA!r}")
    _number(rec, "s_per_rep", positive=True)
    fixed = _number(rec, "fixed_s", positive=False)
    _number(rec, "wall_mc_s", positive=True)
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
    if not (isinstance(rec.get("git_sha"), str) and _HEX40.fullmatch(rec["git_sha"])):
        raise ValueError("git_sha must be 40 lowercase hex digits")
    if not (isinstance(rec.get("log_dir"), str) and rec["log_dir"]):
        raise ValueError("log_dir must be a non-empty string")
    for key in FILE_KEYS:
        if not (isinstance(rec.get(key), str) and rec[key]):
            raise ValueError(f"{key} must be a non-empty file name")
        if not (isinstance(rec.get(f"{key}_sha256"), str) and _HEX64.fullmatch(rec[f"{key}_sha256"])):
            raise ValueError(f"{key}_sha256 must be 64 lowercase hex digits")
    if not (isinstance(rec.get("trust_file_sha256"), str) and _HEX64.fullmatch(rec["trust_file_sha256"])):
        raise ValueError("trust_file_sha256 must be 64 lowercase hex digits")
    _validate_signature(rec.get("host_run_signature"), "host_run_signature")
    _validate_signature(rec.get("validation_signature"), "validation_signature")
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


# -- the signed attestation chain ---------------------------------------------------------------------------------
def token_line(token: str) -> str:
    """The stdout line by which the chain wrapper announces its in-run token (``IONMC_HOST_RUN_ID``)."""
    return f"{TOKEN_PREFIX} {token}"


def artifact_lines(part_record: Path, part_log: Path) -> list[str]:
    """The two stdout lines the chain wrapper prints AFTER the calibration part (end of stdout), binding the file
    bytes of the part record and its pytest log to the signed host-runner record."""
    return [f"{ARTIFACT_PREFIX} part_record {part_record.name} {_sha256(part_record)}",
            f"{ARTIFACT_PREFIX} part_log {part_log.name} {_sha256(part_log)}"]  # fmt: skip


def _trust(trust_path: str | os.PathLike[str] | None) -> tuple[dict[str, dict[str, str]], str]:
    """``(trust map, sha256 of the trust file)`` of ``trust_path`` (default :data:`TRUST_FILE`) or ``ValueError``."""
    p = Path(trust_path) if trust_path is not None else TRUST_FILE
    try:
        raw = p.read_bytes()
        trusted = attest_verify.load_trusted(p)
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f"trust file {p} unreadable: {exc!r}") from exc
    return trusted, hashlib.sha256(raw).hexdigest()


def _read_json(path: Path, what: str) -> dict[str, Any]:
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"{what} {path.name} unreadable: {exc}") from exc
    if not isinstance(doc, dict):
        raise ValueError(f"{what} {path.name} is not a JSON object")
    return doc


def _signature_of(doc: dict[str, Any]) -> dict[str, Any]:
    block = doc.get("attestation")
    return {k: block.get(k) for k in SIGNATURE_FIELDS} if isinstance(block, dict) else {}


def verify_validation_record(doc: dict[str, Any], git_sha: str, trusted: dict[str, Any], name: str) -> tuple[str, str]:
    """The signed validation-manager record copy ``doc`` (module docstring): ``(host run id, claimed output manifest
    sha256)`` of its single ``details.host_runs`` entry, or ``ValueError``."""
    try:
        attest_verify.require_verified(doc, VAL_SIGNER, trusted, expected_task_id=TASK_ID, expected_git_sha=git_sha)
    except attest_verify.AttestationRequired as exc:
        raise ValueError(f"validation record {name}: {exc}") from exc
    vid = doc.get("validation_id")
    if not (isinstance(vid, str) and _VALIDATION_ID.fullmatch(vid)):
        raise ValueError(f"validation record {name}: validation_id {vid!r} is not a VAL- id")
    if doc.get("status") != "passed":
        raise ValueError(f"validation record {name}: status {doc.get('status')!r} is not 'passed'")
    det = doc.get("details")
    if not isinstance(det, dict):
        raise ValueError(f"validation record {name}: details is not an object")
    if "host_run_id" in det or "host_run_ids" in det:
        raise ValueError(f"validation record {name}: legacy details.host_run_id / host_run_ids (no manifest binding)")
    if det.get("host_run_manifest_lookup") == "missing":
        raise ValueError(f"validation record {name}: details.host_run_manifest_lookup is 'missing'")
    runs = det.get("host_runs")
    if not (isinstance(runs, list) and len(runs) == 1 and isinstance(runs[0], dict)):
        raise ValueError(f"validation record {name}: details.host_runs must reference exactly one host run")
    run_id, claim = runs[0].get("run_id"), runs[0].get("output_manifest_sha256")
    if not (isinstance(run_id, str) and _RUN_ID.fullmatch(run_id)):
        raise ValueError(f"validation record {name}: details.host_runs[0].run_id {run_id!r} is not a RUN- id")
    if not (isinstance(claim, str) and _HEX64.fullmatch(claim)):
        raise ValueError(f"validation record {name}: details.host_runs[0].output_manifest_sha256 {claim!r} is not "
                         f"64 hex (the measurement run publishes its outputs)")  # fmt: skip
    return run_id, claim


def verify_host_run_record(doc: dict[str, Any], git_sha: str, run_id: str, manifest: str, trusted: dict[str, Any],
                           name: str) -> None:  # fmt: skip
    """The signed host-runner record copy ``doc`` of run ``run_id`` with output manifest ``manifest`` (module
    docstring), or ``ValueError``."""
    if doc.get("run_id") != run_id:
        raise ValueError(f"host run record {name}: run_id {doc.get('run_id')!r} is not the run {run_id!r} the "
                         f"validation record references")  # fmt: skip
    try:
        attest_verify.require_verified(doc, RUN_SIGNER, trusted, expected_task_id=TASK_ID, expected_git_sha=git_sha,
                                       expected_run_id=run_id, expected_output_manifest_sha256=manifest)  # fmt: skip
    except attest_verify.AttestationRequired as exc:
        raise ValueError(f"host run record {name}: {exc}") from exc
    for key, want in (("succeeded", True), ("timed_out", False), ("committed_state_unchanged", True),
                      ("dirty_after", False)):  # fmt: skip
        if doc.get(key) is not want:
            raise ValueError(f"host run record {name}: {key} is {doc.get(key)!r}, required {want!r} (run did not "
                             f"succeed)")  # fmt: skip
    if isinstance(doc.get("exit_code"), bool) or doc.get("exit_code") != 0:
        raise ValueError(f"host run record {name}: exit_code is {doc.get('exit_code')!r}, required 0")
    if "output_publication_error" not in doc or doc["output_publication_error"] is not None:
        raise ValueError(f"host run record {name}: output_publication_error is "
                         f"{doc.get('output_publication_error', '<absent>')!r}, required null")  # fmt: skip
    out, argv = doc.get("stdout"), doc.get("argv")
    if not isinstance(out, str) or not (isinstance(argv, list) and all(isinstance(a, str) for a in argv)):
        raise ValueError(f"host run record {name}: stdout / argv missing or malformed")
    if doc.get("stdout_truncated") is False and hashlib.sha256(out.encode("utf-8")).hexdigest() != doc.get(
            "stdout_sha256"):  # fmt: skip
        raise ValueError(f"host run record {name}: stdout does not hash to the signed stdout_sha256")


def verify_run_binding(run: dict[str, Any], token: str, files: dict[str, Path], name: str) -> None:
    """The sealed token and the file bytes of the part record and the part log are bound to the signed run (module
    docstring), or ``ValueError``."""
    lines = [ln.rstrip("\r") for ln in run["stdout"].splitlines()]
    if not (any(token in a for a in run["argv"]) or token_line(token) in lines):
        raise ValueError(f"host run record {name}: the sealed in-run token {token!r} appears neither in the signed "
                         f"argv nor as a stdout line '{token_line(token)}'")  # fmt: skip
    found = [m.groups() for ln in lines if (m := _ARTIFACT.match(ln))]
    for kind in ("part_record", "part_log"):
        want = (kind, files[kind].name, _sha256(files[kind]))
        mine = [g for g in found if g[0] == kind]
        if mine != [want]:
            raise ValueError(f"host run record {name}: the signed stdout needs exactly one line "
                             f"'{' '.join((ARTIFACT_PREFIX, *want))}', found {mine}")  # fmt: skip


def verify_chain(git_sha: str, files: dict[str, Path], token: str,
                 trust_path: str | os.PathLike[str] | None = None) -> dict[str, Any]:  # fmt: skip
    """Verify the signed VAL- and RUN- copies (``files["validation_record"]`` / ``files["host_run_record"]``) and the
    run binding of the part files; returns the attestation fields of the rate record, or ``ValueError``."""
    trusted, trust_sha = _trust(trust_path)
    vpath, rpath = files["validation_record"], files["host_run_record"]
    val = _read_json(vpath, "validation record")
    run_id, manifest = verify_validation_record(val, git_sha, trusted, vpath.name)
    run = _read_json(rpath, "host run record")
    verify_host_run_record(run, git_sha, run_id, manifest, trusted, rpath.name)
    verify_run_binding(run, token, files, rpath.name)
    return {"host_run_id": run_id, "host_run_signature": _signature_of(run),
            "validation_id": val["validation_id"], "validation_signature": _signature_of(val),
            "trust_file_sha256": trust_sha}  # fmt: skip


def verify_artifacts(rec: dict[str, Any], trust_path: str | os.PathLike[str] | None = None) -> None:
    """Re-open the four files named by the (shape-validated) ``rec``, re-verify the signatures against the trust file
    (``trust_path``, default :data:`TRUST_FILE`) and check everything against ``rec`` (module docstring);
    ``ValueError`` on the first mismatch."""
    d = _resolve_dir(rec["log_dir"])
    files = {}
    for key in FILE_KEYS:
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
    chain = verify_chain(rec["git_sha"], files, doc["host_run_id"], trust_path)  # V3-INFRA-002 signatures
    for key, got in chain.items():
        if rec[key] != got:
            what = "trust file" if key == "trust_file_sha256" else "signed record copy"
            raise ValueError(f"{key} of the record differs from the {what} ({rec[key]!r} != {got!r})")


def record_digest(rec: dict[str, Any]) -> str:
    """sha256 of the canonical JSON (sorted keys, compact separators) of the rate record: the digest a pin names."""
    return hashlib.sha256(json.dumps(rec, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _pin_fields(rec: dict[str, Any]) -> dict[str, str]:
    """The pin fields of the (validated) record, in the order of the pin line."""
    return {"sha256": record_digest(rec), "git_sha": rec["git_sha"], "host_run_id": rec["host_run_id"],
            "validation_id": rec["validation_id"], "host_run_record_sha256": rec["host_run_record_sha256"],
            "validation_record_sha256": rec["validation_record_sha256"],
            "trust_file_sha256": rec["trust_file_sha256"]}  # fmt: skip


def pin_line(rec: dict[str, Any]) -> str:
    """The exact pin line of the (validated) record, as the orchestrator writes it into the acceptance plan."""
    return PIN_PREFIX + "".join(f" {k}={v}" for k, v in _pin_fields(rec).items())


def plan_pin(plan_path: str | os.PathLike[str] | None = None) -> dict[str, str] | None:
    """The LAST pin line of the plan (default: the repository acceptance plan) as a dict, or ``None`` if there is none.
    A pin is a whole line of the exact grammar of the module docstring; a placeholder with angle brackets, and a pin of
    an older grammar (``host_run_token`` / ``host_node``), is no pin."""
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
    """``ValueError`` unless the last pin of the plan names exactly this record (digest, git sha, runner id,
    validation id, the two record-copy digests and the trust-file digest)."""
    pin = plan_pin(plan_path)
    if pin is None:
        raise ValueError(f"the record is {NOT_PINNED}")
    bad = [k for k, v in _pin_fields(rec).items() if pin[k] != v]
    if bad:
        raise ValueError(f"the last pin of the acceptance plan does not match the record ({', '.join(bad)} differ); "
                         f"the record is {NOT_PINNED}")  # fmt: skip


def _load_verified(path: str | os.PathLike[str], trust_path: str | os.PathLike[str] | None = None) -> dict[str, Any]:
    """Shape, artifact and signature verification only (no pin): ``{"record", "sha256", "path"}`` or ``ValueError``."""
    p = Path(path)
    if not p.is_file():
        raise ValueError(f"{p} is not a file")
    raw = p.read_bytes()
    try:
        rec = validate_record(json.loads(raw))
        verify_artifacts(rec, trust_path)
    except ValueError as exc:  # JSON errors are ValueErrors as well
        raise ValueError(f"calibration rate record {p} invalid: {exc}") from exc
    return {"record": rec, "sha256": hashlib.sha256(raw).hexdigest(), "path": str(p)}


def load_record(path: str | os.PathLike[str], plan_path: str | os.PathLike[str] | None = None,
                trust_path: str | os.PathLike[str] | None = None) -> dict[str, Any]:  # fmt: skip
    """``{"record", "sha256", "path"}`` of the rate record file, or ``ValueError`` (absent, unreadable, malformed,
    legacy, whose artifacts or signed record copies fail :func:`verify_artifacts`, or which is not pinned in the
    acceptance plan). ``plan_path`` / ``trust_path`` (function arguments only) replace the repository plan / trust
    file."""
    got = _load_verified(path, trust_path)
    try:
        verify_pin(got["record"], plan_path)
    except ValueError as exc:
        raise ValueError(f"calibration rate record {got['path']} invalid: {exc}") from exc
    return got


def record_from_env() -> dict[str, Any] | None:
    """:func:`load_record` (repository plan and trust file) of ``IONMC_V7R2_CAL_RATE_RECORD`` (``None`` if unset or
    empty)."""
    path = os.environ.get(CAL_RATE_ENV, "")
    return load_record(path) if path else None


def build_record(part_record: Path, part_log: Path, log_dir: Path, git_sha: str, host_run_record: Path,
                 validation_record: Path, trust_path: str | os.PathLike[str] | None = None) -> dict[str, Any]:  # fmt: skip
    """The record of one measured calibration part (module docstring). Every value comes from the verified sealed
    part record, pytest log and the signed RUN- / VAL- record copies (all four in ``log_dir``); the result passes
    :func:`verify_artifacts` or ``ValueError`` is raised."""
    if not (isinstance(git_sha, str) and _HEX40.fullmatch(git_sha)):
        raise ValueError("git_sha must be 40 lowercase hex digits")
    files = {"part_record": part_record, "part_log": part_log, "host_run_record": host_run_record,
             "validation_record": validation_record}  # fmt: skip
    if not all(f.is_file() for f in files.values()):
        raise ValueError("the part record, the pytest log and the host-run and validation record copies must be files")
    ld = log_dir if log_dir.is_absolute() else Path.cwd() / log_dir
    if any(f.resolve().parent != ld.resolve() for f in files.values()):
        raise ValueError("the part record, the pytest log and the two record copies must lie in log_dir")
    doc = verify_part_record(part_record, git_sha)
    total = verify_log(part_log.read_text(errors="replace"), part_record.name, doc["content_sha256"])
    chain = verify_chain(git_sha, files, doc["host_run_id"], trust_path)
    n, wall = doc["n"], float(doc["wall_mc_s"])
    try:
        stored = str(ld.resolve().relative_to(REPO))
    except ValueError:
        stored = str(ld.resolve())
    rec: dict[str, Any] = {
        "schema": CAL_RATE_SCHEMA, "s_per_rep": wall / n, "fixed_s": max(MIN_FIXED_S, total - wall),
        "reps": n, "wall_mc_s": wall, "wall_total_s": total,
        "host": {k: doc["host"][k] for k in HOST_FIELDS}, "git_sha": git_sha, "log_dir": stored, **chain,
    }  # fmt: skip
    for key, f in files.items():
        rec[key], rec[f"{key}_sha256"] = f.name, _sha256(f)
    validate_record(rec)
    verify_artifacts(rec, trust_path)
    return rec


def require_validated(got: Any, plan_path: str | os.PathLike[str] | None = None,
                      trust_path: str | os.PathLike[str] | None = None) -> dict[str, Any]:  # fmt: skip
    """``got`` if it is exactly the result of :func:`load_record` of the file it names, else ``ValueError``: the file
    is re-loaded (all artifact and signature checks and the plan pin again) and the result must be identical; a
    hand-made dictionary, a record whose file changed or whose pin was removed, or ``None`` is refused."""
    if not (isinstance(got, dict) and isinstance(got.get("path"), str) and isinstance(got.get("sha256"), str)):
        raise ValueError("not a calibration rate record loaded by v7r2_cal_rate.load_record")
    fresh = load_record(got["path"], plan_path, trust_path)
    if fresh != got:
        raise ValueError("the calibration rate record differs from the file it names (not the result of load_record)")
    return fresh


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--part-record", required=True, help="sealed v7r2-calibration-<surrogate>-p<NN>.json of the part")
    ap.add_argument("--part-log", required=True, help="the logged pytest output of that single part")
    ap.add_argument("--log-dir", required=True, help="directory holding the four files (recorded)")
    ap.add_argument("--git-sha", required=True, help="full SHA of the tree the part was measured at")
    ap.add_argument("--host-run-record", required=True,
                    help="verbatim copy of the signed host-runner record RUN-....json (inside --log-dir)")  # fmt: skip
    ap.add_argument("--validation-record", required=True,
                    help="verbatim copy of the signed validation-manager record VAL-....json (inside --log-dir)")  # fmt: skip
    ap.add_argument("--trust-file", default=None,
                    help=f"trust file of the attestation keys (default: the committed {TRUST_FILE.name})")  # fmt: skip
    ap.add_argument("--out", required=True, help="rate record to write (refuses to overwrite)")
    ap.add_argument("--print-pin-line", action="store_true",
                    help="also print the pin line of the record for the acceptance plan (the orchestrator writes it "
                         "in a reviewed commit; this tool never writes the plan)")  # fmt: skip
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists")
    try:
        rec = build_record(Path(a.part_record), Path(a.part_log), Path(a.log_dir), a.git_sha,
                           Path(a.host_run_record), Path(a.validation_record), a.trust_file)  # fmt: skip
    except ValueError as exc:
        raise SystemExit(f"not ok: {exc}") from exc
    out.write_text(json.dumps(rec, sort_keys=True, indent=1) + "\n")
    print(f"{out}: s_per_rep {rec['s_per_rep']:.4f}, fixed_s {rec['fixed_s']:.1f}, "
          f"sha256 {_load_verified(out, a.trust_file)['sha256']}")  # fmt: skip
    if a.print_pin_line:
        print(pin_line(rec))
    return 0


if __name__ == "__main__":
    sys.exit(main())
