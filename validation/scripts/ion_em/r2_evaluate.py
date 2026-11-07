"""V3-008 revision 2: evaluator skeleton of the TOPAS-agreement rows (TD-EM, TD-FULL, TL-EM, TL-FULL).

Implemented now (design record ``validation/plans/v3-008/r2-topas-target.md`` section 6, "Run validity
against the seal"): the two chronology checks, failing closed on any missing record.

1. Ancestry: ``git merge-base --is-ancestor <seal> <code_sha>`` succeeds (the seal itself counts).
2. Review before start: the controlled review record proves that the seal commit itself passed its
   independent review (``head_sha`` equals the plan's ``seal_sha``, ``status`` is ``passed``, report
   verdict ``pass`` when present) and its ``finished_at`` is strictly earlier than the run's
   ``reference_started`` (records of the controlled services, never Git timestamps; the run's own
   ``finished_at`` is not used).

A run failing either check invalidates the TD/TL rows of its case (FAIL of the blinding control,
never re-graded). The TD/TL statistics are stubs that raise ``NotImplementedError`` (fail closed): they
are implemented with the rows, after the sample seal.

Record shapes (dicts; ISO-8601 timestamps with an explicit UTC offset or ``Z``):

* ``plan``: ``{"seal_sha": str, "seal_review_id": str}``;
* ``review``: the controlled review record of the codex-worker service, fields used:
  ``{"review_id", "head_sha", "status", "finished_at"}`` and ``report.verdict`` when a report is present;
* ``run``: ``{"code_sha": str, "reference_started": str}`` plus ``run_id`` when the caller binds it.

Reference-run loader contract (the loader is still a stub and not part of this module): the ``run``
record handed to :func:`check_run_validity` must be the run's own signed/controlled record (the
reference-runner record bound to the run id and the executed code SHA); ``reference_started`` and
``code_sha`` are taken from that one record, never from an unlinked timestamp or a separately supplied
SHA. Pass ``expected_run_id`` to bind the record to the run being graded.
"""

from __future__ import annotations

import re
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

SHA_RE = re.compile(r"^[0-9a-f]{40}$")


class ChronologyError(Exception):
    """A chronology record is missing or malformed: the check fails closed (never passes)."""


def _require(record: dict[str, Any] | None, key: str, what: str) -> Any:
    if not isinstance(record, dict) or record.get(key) in (None, ""):
        raise ChronologyError(f"missing record: {what} ({key})")
    return record[key]


def _sha(value: Any, what: str) -> str:
    if not isinstance(value, str) or not SHA_RE.match(value):
        raise ChronologyError(f"{what} is not a full 40-hex commit SHA: {value!r}")
    return value


def _time(value: Any, what: str) -> datetime:
    if not isinstance(value, str):
        raise ChronologyError(f"{what} is not a timestamp string: {value!r}")
    try:
        t = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ChronologyError(f"{what} is not an ISO-8601 timestamp: {value!r}") from exc
    if t.tzinfo is None:
        raise ChronologyError(f"{what} has no UTC offset: {value!r}")
    return t.astimezone(UTC)


def check_ancestry(repo: Path | str, seal_sha: Any, code_sha: Any) -> bool:
    """True iff ``seal_sha`` is an ancestor of (or equal to) ``code_sha`` in ``repo``; raises
    :class:`ChronologyError` for a missing/malformed SHA or a git failure other than "not ancestor"."""
    seal = _sha(seal_sha, "seal SHA")
    code = _sha(code_sha, "run code_sha")
    proc = subprocess.run(  # noqa: S603 - fixed argv, validated hex SHAs
        ["git", "-C", str(repo), "merge-base", "--is-ancestor", seal, code],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode == 0:
        return True
    if proc.returncode == 1:
        return False
    raise ChronologyError(f"git merge-base failed (rc={proc.returncode}): {proc.stderr.strip()}")


def check_review_before_start(
    review: dict[str, Any] | None, plan: dict[str, Any] | None, run: dict[str, Any] | None
) -> bool:
    """True iff the seal review recorded in ``plan`` passed for the seal commit and finished strictly
    before the run started. Any missing, malformed, mismatched or non-passing field raises
    :class:`ChronologyError`."""
    review_id = _require(plan, "seal_review_id", "plan seal note")
    seal = _sha(_require(plan, "seal_sha", "plan seal note"), "plan seal_sha")
    if _require(review, "review_id", "seal review record") != review_id:
        raise ChronologyError("the review record is not the review named by the plan's seal note")
    head = _sha(_require(review, "head_sha", "seal review head_sha"), "review head_sha")
    if head != seal:
        raise ChronologyError(
            f"the review covers commit {head}, not the plan's seal commit {seal}"
        )
    status = _require(review, "status", "seal review status")
    if status != "passed":
        raise ChronologyError(f"the seal review did not pass (status={status!r})")
    report = review.get("report") if isinstance(review, dict) else None
    if report is not None:
        if not isinstance(report, dict):
            raise ChronologyError(f"the review report is malformed: {report!r}")
        if "verdict" in report and report["verdict"] != "pass":
            raise ChronologyError(f"the seal review verdict is not pass: {report['verdict']!r}")
    finished = _time(_require(review, "finished_at", "seal review finished_at"), "finished_at")
    started = _time(_require(run, "reference_started", "run reference_started"), "reference_started")
    return finished < started


def check_run_validity(
    repo: Path | str,
    plan: dict[str, Any] | None,
    review: dict[str, Any] | None,
    run: dict[str, Any] | None,
    expected_run_id: str | None = None,
) -> dict[str, Any]:
    """Both chronology checks of one TOPAS run. ``valid`` is True only if both hold; any missing
    record raises :class:`ChronologyError` (fail closed). ``run`` must be the run's own controlled
    record (see the module docstring); with ``expected_run_id`` its ``run_id`` must match."""
    if expected_run_id is not None and _require(run, "run_id", "run record") != expected_run_id:
        raise ChronologyError("the run record is not bound to the run being evaluated")
    seal = _require(plan, "seal_sha", "plan seal note")
    ancestry = check_ancestry(repo, seal, _require(run, "code_sha", "run record"))
    reviewed = check_review_before_start(review, plan, run)
    return {"ancestry": ancestry, "review_before_start": reviewed, "valid": ancestry and reviewed}


def td_statistics(*_args: Any, **_kwargs: Any) -> Any:
    """TD-EM / TD-FULL statistics (R80, W, z_peak, PE, plateau groups; TOST with the pre-registered
    tolerances). Not implemented: fails closed."""
    raise NotImplementedError("TD statistics are implemented with the rows after the sample seal")


def tl_statistics(*_args: Any, **_kwargs: Any) -> Any:
    """TL-EM / TL-FULL statistics (core-width fits, ratio band, D_ext). Not implemented: fails
    closed."""
    raise NotImplementedError("TL statistics are implemented with the rows after the sample seal")
