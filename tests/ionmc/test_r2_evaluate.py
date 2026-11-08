"""V3-008 revision 2: chronology checks of ``validation/scripts/ion_em/r2_evaluate.py`` (synthetic
records, temporary git repository)."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

_PATH = Path(__file__).resolve().parents[2] / "validation" / "scripts" / "ion_em" / "r2_evaluate.py"
_spec = importlib.util.spec_from_file_location("r2_evaluate_under_test", _PATH)
assert _spec is not None and _spec.loader is not None
R2 = importlib.util.module_from_spec(_spec)
sys.modules["r2_evaluate_under_test"] = R2
_spec.loader.exec_module(R2)


def _git(repo: Path, *args: str) -> str:
    env = {
        "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@x", "GIT_COMMITTER_NAME": "t",
        "GIT_COMMITTER_EMAIL": "t@x", "PATH": "/usr/bin:/bin:/usr/local/bin", "HOME": str(repo),
    }  # fmt: skip
    out = subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=True, env=env
    )
    return out.stdout.strip()


@pytest.fixture
def repo(tmp_path: Path) -> dict[str, Any]:
    _git(tmp_path, "init", "-q", "-b", "main")
    shas = {}
    for name in ("a", "seal", "after"):
        (tmp_path / "f.txt").write_text(name)
        _git(tmp_path, "add", "f.txt")
        _git(tmp_path, "-c", "commit.gpgsign=false", "commit", "-q", "-m", name)
        shas[name] = _git(tmp_path, "rev-parse", "HEAD")
    _git(tmp_path, "checkout", "-q", "-b", "side", shas["a"])
    (tmp_path / "g.txt").write_text("side")
    _git(tmp_path, "add", "g.txt")
    _git(tmp_path, "-c", "commit.gpgsign=false", "commit", "-q", "-m", "side")
    shas["side"] = _git(tmp_path, "rev-parse", "HEAD")
    return {"path": tmp_path, **shas}


def test_ancestry(repo: dict[str, Any]) -> None:
    p, seal = repo["path"], repo["seal"]
    assert R2.check_ancestry(p, seal, repo["after"]) is True
    assert R2.check_ancestry(p, seal, seal) is True  # the seal commit itself counts
    assert R2.check_ancestry(p, seal, repo["a"]) is False  # predates the seal
    assert R2.check_ancestry(p, seal, repo["side"]) is False  # divergent branch
    for bad in (None, "", "abc", "z" * 40, repo["after"][:12]):
        with pytest.raises(R2.ChronologyError):
            R2.check_ancestry(p, bad, repo["after"])
        with pytest.raises(R2.ChronologyError):
            R2.check_ancestry(p, seal, bad)
    with pytest.raises(R2.ChronologyError):  # unknown (well-formed) commit: git error, not "no"
        R2.check_ancestry(p, "0" * 40, repo["after"])


SEAL = "a" * 40
T_OK = {"reference_started": "2026-11-08T11:00:00Z"}
PLAN = {"seal_sha": SEAL, "seal_review_id": "REVIEW-1"}
REVIEW = {
    "review_id": "REVIEW-1",
    "head_sha": SEAL,
    "status": "passed",
    "finished_at": "2026-11-08T10:00:00Z",
    "report": {"verdict": "pass"},
}


def test_review_before_start() -> None:
    ok = {"reference_started": "2026-11-08T10:00:01+00:00"}
    assert R2.check_review_before_start(REVIEW, PLAN, ok) is True
    # offsets are compared as instants
    assert R2.check_review_before_start(
        REVIEW, PLAN, {"reference_started": "2026-11-08T12:00:30+02:00"}
    )
    for t in ("2026-11-08T10:00:00Z", "2026-11-08T09:59:59Z"):  # equal or earlier start: fails
        assert R2.check_review_before_start(REVIEW, PLAN, {"reference_started": t}) is False
    # a run begun before the review but finished after must not pass (finished_at is ignored)
    late = {"reference_started": "2026-11-08T09:00:00Z", "finished_at": "2026-11-08T11:00:00Z"}
    assert R2.check_review_before_start(REVIEW, PLAN, late) is False


@pytest.mark.parametrize(
    ("review", "plan", "run"),
    [
        (None, PLAN, T_OK),
        ({"review_id": "REVIEW-1"}, PLAN, T_OK),
        (REVIEW, {"seal_sha": SEAL}, {"reference_started": "2026-11-08T11:00:00Z"}),
        (REVIEW, None, {"reference_started": "2026-11-08T11:00:00Z"}),
        (REVIEW, PLAN, {}),
        (REVIEW, PLAN, None),
        (REVIEW, PLAN, {"reference_started": "yesterday"}),
        (REVIEW, PLAN, {"reference_started": "2026-11-08T11:00:00"}),  # no UTC offset
        ({**REVIEW, "review_id": "OTHER"}, PLAN, T_OK),
        ({**REVIEW, "head_sha": "b" * 40}, PLAN, T_OK),
        ({**REVIEW, "head_sha": SEAL[:12]}, PLAN, T_OK),
        ({**REVIEW, "status": "changes_required"}, PLAN, T_OK),
        ({**REVIEW, "status": "running"}, PLAN, T_OK),
        ({**REVIEW, "report": {"verdict": "fail"}}, PLAN, T_OK),
        ({**REVIEW, "report": "pass"}, PLAN, T_OK),
        ({k: v for k, v in REVIEW.items() if k != "report"}, PLAN, T_OK),  # no report
        ({**REVIEW, "report": None}, PLAN, T_OK),
        ({**REVIEW, "report": {}}, PLAN, T_OK),  # report without verdict
        ({**REVIEW, "report": {"summary": "ok"}}, PLAN, T_OK),
        ({**REVIEW, "report": {"verdict": None}}, PLAN, T_OK),
        ({**REVIEW, "report": {"verdict": "PASS"}}, PLAN, T_OK),
        ({**REVIEW, "finished_at": None}, PLAN, T_OK),
        ({**REVIEW, "finished_at": "soon"}, PLAN, T_OK),
        ({**REVIEW, "finished_at": "2026-11-08T10:00:00"}, PLAN,
         {"reference_started": "2026-11-08T11:00:00Z"}),
        ({k: v for k, v in REVIEW.items() if k != "head_sha"}, PLAN,
         {"reference_started": "2026-11-08T11:00:00Z"}),
        ({k: v for k, v in REVIEW.items() if k != "status"}, PLAN,
         {"reference_started": "2026-11-08T11:00:00Z"}),
    ],
)  # fmt: skip
def test_review_check_fails_closed(review: Any, plan: Any, run: Any) -> None:
    with pytest.raises(R2.ChronologyError):
        R2.check_review_before_start(review, plan, run)


def test_review_finished_after_start() -> None:
    # finished_at after reference_started: not before start (False, the run is invalid)
    late = {**REVIEW, "finished_at": "2026-11-08T12:00:00Z"}
    assert not R2.check_review_before_start(
        late, PLAN, {"reference_started": "2026-11-08T11:00:00Z"}
    )


def _good_run(repo: dict[str, Any], **kw: Any) -> dict[str, Any]:
    run = {"run_id": "r1", "code_sha": repo["after"], "reference_started": "2026-11-08T11:00:00Z"}
    run.update(kw)
    return run


def test_run_record_binding(repo: dict[str, Any]) -> None:
    plan = {"seal_sha": repo["seal"], "seal_review_id": "REVIEW-1"}
    review = {**REVIEW, "head_sha": repo["seal"]}
    run = _good_run(repo)
    assert R2.check_run_validity(repo["path"], plan, review, run, expected_run_id="r1")["valid"]
    with pytest.raises(R2.ChronologyError):  # mismatched run id
        R2.check_run_validity(repo["path"], plan, review, run, expected_run_id="r2")
    for key in ("run_id", "code_sha", "reference_started"):  # every field is mandatory
        partial = {k: v for k, v in run.items() if k != key}
        with pytest.raises(R2.ChronologyError):
            R2.check_run_validity(repo["path"], plan, review, partial, expected_run_id="r1")
    for bad in (None, ""):  # the binding is mandatory: no unlinked default
        with pytest.raises(R2.ChronologyError):
            R2.check_run_validity(repo["path"], plan, review, run, expected_run_id=bad)
    with pytest.raises(TypeError):  # no default value for expected_run_id
        R2.check_run_validity(repo["path"], plan, review, run)  # type: ignore[call-arg]
    with pytest.raises(R2.ChronologyError):
        R2.check_run_validity(repo["path"], plan, review, None, expected_run_id="r1")


def test_run_validity_combined(repo: dict[str, Any]) -> None:
    plan = {"seal_sha": repo["seal"], "seal_review_id": "REVIEW-1"}
    review = {**REVIEW, "head_sha": repo["seal"]}
    good = _good_run(repo)
    assert R2.check_run_validity(repo["path"], plan, review, good, "r1")["valid"] is True
    early = _good_run(repo, reference_started="2026-11-08T09:00:00Z")
    v = R2.check_run_validity(repo["path"], plan, review, early, "r1")
    assert v == {"ancestry": True, "review_before_start": False, "valid": False}
    old = _good_run(repo, code_sha=repo["a"])
    assert R2.check_run_validity(repo["path"], plan, review, old, "r1")["valid"] is False
    with pytest.raises(R2.ChronologyError):
        R2.check_run_validity(repo["path"], {"seal_review_id": "REVIEW-1"}, review, good, "r1")
    for rep in (None, {}, {"verdict": "fail"}):  # the seal review must carry verdict pass
        bad = {**review, "report": rep}
        with pytest.raises(R2.ChronologyError):
            R2.check_run_validity(repo["path"], plan, bad, good, "r1")


def test_statistics_stubs_fail_closed() -> None:
    with pytest.raises(NotImplementedError):
        R2.td_statistics()
    with pytest.raises(NotImplementedError):
        R2.tl_statistics()
