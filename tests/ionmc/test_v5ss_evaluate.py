"""V5-SS evaluator ``validation/scripts/ion_em/ss_evaluate.py`` (synthetic data, temporary git
repositories; CI-fast). The phase-space fixture ``tests/ionmc/data/v5ss/`` is in the ASSUMED
OpenTOPAS ASCII layout (TO VERIFY on the scorer smoke)."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
_PATH = ROOT / "validation" / "scripts" / "ion_em" / "ss_evaluate.py"
_spec = importlib.util.spec_from_file_location("ss_evaluate_under_test", _PATH)
assert _spec is not None and _spec.loader is not None
SS = importlib.util.module_from_spec(_spec)
sys.modules["ss_evaluate_under_test"] = SS
_spec.loader.exec_module(SS)
IO, BUDGET = SS.io, SS.ss_budget
FIX = ROOT / "tests" / "ionmc" / "data" / "v5ss"
KEY_SEED = bytes(range(32))
SHA = "1" * 40


# ---------------------------------------------------------------------------------- statistics
def test_windowed_rms_of_a_gaussian_matches_the_closed_form() -> None:
    rng = np.random.default_rng(1)
    sig, n = 1e-3, 200_000
    tx, ty = rng.normal(0, sig, n), rng.normal(0, sig, n)
    got, used = SS.windowed_rms(tx, ty, 2.5 * sig)
    want = BUDGET.gauss_windowed(sig, 2.5 * sig)
    p, m2, m4 = BUDGET.gauss_window_moments(sig, 2.5 * sig)
    se = BUDGET.rms_rel_se(p, m2, m4, n, joint=(p * p, m2 * p, m2 * m2)) * want
    assert abs(got - want) <= 3 * se
    assert used == pytest.approx(2 * n * p, rel=0.01)
    assert SS.windowed_rms(tx, ty, None)[0] == pytest.approx(sig, rel=0.01)


def test_bootstrap_sigma_matches_the_delta_method() -> None:
    rng = np.random.default_rng(2)
    sig, n = 1.0, 20_000
    tx, ty = rng.normal(0, sig, n), rng.normal(0, sig, n)
    s, c = SS.contributions(tx, ty, 2.5 * sig)
    reps = SS.boot_rms(s, c, [n // 2, n - n // 2], 600, np.random.default_rng(3))
    rel = np.std(reps, ddof=1) / SS.windowed_rms(tx, ty, 2.5 * sig)[0]
    p, m2, m4 = BUDGET.gauss_window_moments(sig, 2.5 * sig)
    want = BUDGET.rms_rel_se(p, m2, m4, n, joint=(p * p, m2 * p, m2 * m2))
    assert rel == pytest.approx(want, rel=0.10)


def _samples(widths: dict[int, float], n_ref: int, n_model: int, seed: int, ref_scale=None):
    rng = np.random.default_rng(seed)
    ref_scale = ref_scale or {}
    ref, model, windows = {}, {}, {}
    for c, w in widths.items():
        model[c] = (rng.normal(0, w, n_model), rng.normal(0, w, n_model))
        wr = w * ref_scale.get(c, 1.0)
        h = n_ref // 2
        ref[c] = (
            [rng.normal(0, wr, h), rng.normal(0, wr, h)],
            [rng.normal(0, wr, h), rng.normal(0, wr, h)],
        )
        windows[c] = SS.WINDOW_K * SS.windowed_rms(*model[c], None)[0]
    return ref, model, windows


WIDTHS = {
    22: 2e-3,
    25: 1.5e-3,
    26: 1.4e-3,
    27: 1.0e-3,
    28: 1.3e-3,
    29: 0.7e-3,
    30: 0.65e-3,
    31: 0.9e-3,
}


def test_identical_samples_give_zero_dev_and_pass() -> None:
    ref, model, windows = _samples(WIDTHS, 20_000, 20_000, 4)
    same = {c: (np.concatenate(ref[c][0]), np.concatenate(ref[c][1])) for c in ref}
    res = SS.grade_pairs(ref, same, windows, n_boot=50)
    assert all(abs(p["dev_pct"]) < 1e-12 for p in res["pairs"])
    assert res["row_verdict"] == "PASS"
    assert [p["expected_dev_pct"] for p in res["pairs"]] == [-0.49, -0.69, -0.38, -0.57, -0.17]


def test_planted_two_percent_width_change() -> None:
    ref, model, windows = _samples(WIDTHS, 60_000, 60_000, 5, ref_scale={22: 1.02, 28: 1.02})
    res = SS.grade_pairs(ref, model, windows, n_boot=100)
    by = {p["pair"]: p for p in res["pairs"]}
    for name in ("He/C 150", "He/O 150", "He/C 250", "He/O 250"):
        assert by[name]["dev_pct"] == pytest.approx(-1.96, abs=0.8), name
    assert abs(by["C/O 430"]["dev_pct"]) < 0.8
    assert 0.1 < by["He/C 150"]["sigma_pct"] < 1.0


def test_bootstrap_is_reproducible() -> None:
    ref, model, windows = _samples(WIDTHS, 4_000, 4_000, 6)
    a = SS.grade_pairs(ref, model, windows, n_boot=20)
    b = SS.grade_pairs(ref, model, windows, n_boot=20)
    assert a == b and a["analysis_seed"] == 20261113


def test_gate_examples_and_row() -> None:
    assert BUDGET.gate(2.99, 0.9) == "FAIL"
    assert BUDGET.gate(-0.69, 0.30) == "PASS"
    assert BUDGET.gate(0.0, 1.01) == "INCONCLUSIVE"
    ok = [{"verdict": "PASS"}] * 5
    assert SS.row_verdict(ok) == "PASS"
    assert SS.row_verdict(ok[:4] + [{"verdict": "INCONCLUSIVE"}]) == "FAIL"
    assert SS.row_verdict(ok[:4] + [{"verdict": "FAIL"}]) == "FAIL"
    assert SS.row_verdict(ok[:4]) == "FAIL"


# ---------------------------------------------------------------------------------- phase space
def _fixture() -> tuple[str, str]:
    return (FIX / "ss_exit.header").read_text(), (FIX / "ss_exit.phsp").read_text()


def test_parser_fixture() -> None:
    h, rows = _fixture()
    assert SS.parse_header(h)["columns"] == list(SS.EXPECTED_COLUMNS)
    out = SS.parse_phasespace(h, rows, 3, 1000020040)
    u, v = 0.0012345, -0.0006789
    w = math.sqrt(1 - u * u - v * v)
    assert out["theta_x"][0] == pytest.approx(math.atan2(u, w), rel=1e-12)
    assert out["theta_y"][0] == pytest.approx(math.atan2(v, w), rel=1e-12)


def _row_edit(rows: str, col: int, value: str) -> str:
    lines = rows.splitlines()
    f = lines[1].split()
    f[col] = value
    lines[1] = "  ".join(f)
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize(
    "edit",
    [
        lambda h, r: (h.replace("Weight", "Weight X"), r),
        lambda h, r: (h.replace("10: Flag", "11: Flag"), r),
        lambda h, r: (h.replace("Scored Particles: 3", "Scored Particles: 2"), r),
        lambda h, r: (h.replace("Original Histories: 3", "Original Histories: 4"), r),
        lambda h, r: (h, "\n".join(r.splitlines()[:2]) + "\n"),
        lambda h, r: (h, r + r.splitlines()[0] + "\n"),
        lambda h, r: (h, _row_edit(r, 7, "1000060120")),
        lambda h, r: (h, _row_edit(r, 8, "1")),
        lambda h, r: (h, _row_edit(r, 9, "0")),
        lambda h, r: (h, _row_edit(r, 2, "0.05")),
        lambda h, r: (h, _row_edit(r, 3, "nan")),
        lambda h, r: (h, _row_edit(r, 3, "1.2")),
        lambda h, r: (h, _row_edit(r, 4, "abc")),
    ],
)
def test_parser_fails_closed(edit: Any) -> None:
    h, r = edit(*_fixture())
    with pytest.raises(SS.RecordError):
        SS.parse_phasespace(h, r, 3, 1000020040)


def test_unexpected_columns_are_named() -> None:
    h, r = _fixture()
    with pytest.raises(SS.RecordError, match="observed column list .*Weight X"):
        SS.parse_header(h.replace("Weight", "Weight X"))


# ---------------------------------------------------------------------------------- log checks
def _stdout(n: int, particle: str = "alpha", model: str = "IonCoulombScattering") -> str:
    return (
        "Registered Physics Processes:\nTransportation, CoulombScat, ionIoni, eIoni\n\n"
        f"CoulombScat:  for {particle} XStype:1 SubType=1 BuildTable=1\n"
        "      ===== EM models for the G4Region  DefaultRegionForTheWorld ======\n"
        f"  {model} : Emin=    0 eV  Emax=   10 GeV\n\n"
        " Range cuts        :  gamma  50 um     e-  1 km      e+  50 um  proton 50 um\n"
        " Energy thresholds :  gamma  990 eV     e-  9.55 GeV    e+  57 keV proton 5 keV\n"
        f"Particle source Beam: Total number of histories: {n}\n"
    )


def test_log_checks() -> None:
    assert SS.log_checks(_stdout(5), b"", 5, "alpha") == []
    assert SS.log_checks(_stdout(5, "GenericIon"), b"", 5, "GenericIon") == []
    # a global msc table entry without an msc block for the beam particle passes
    glob = _stdout(5).replace("ionIoni", "msc")
    assert SS.log_checks(glob, b"", 5, "alpha") == []
    assert SS.msc_info(glob) == {"msc_blocks_for": [], "registered_msc_entry": True}
    # an msc block for e- is recorded, not a failure
    ele = _stdout(5) + "msc:  for e-  SubType= 10\n"
    assert SS.log_checks(ele, b"", 5, "alpha") == []
    assert SS.msc_info(ele)["msc_blocks_for"] == ["e-"]
    assert SS.msc_info(_stdout(5))["registered_msc_entry"] is False
    bad = [
        SS.log_checks(_stdout(5), b"x", 5, "alpha"),
        SS.log_checks(_stdout(4), b"", 5, "alpha"),
        SS.log_checks(_stdout(5, model="eCoulombScattering"), b"", 5, "alpha"),
        SS.log_checks(_stdout(5) + "msc:  for alpha  SubType= 10\n", b"", 5, "alpha"),
        SS.log_checks(_stdout(5).replace("CoulombScat,", "ionIoni,"), b"", 5, "alpha"),
        SS.log_checks(_stdout(5) + "msc:  for GenericIon  SubType= 10\n", b"", 5, "GenericIon"),
        SS.log_checks(_stdout(5).replace("9.55 GeV", "990 keV"), b"", 5, "alpha"),
        SS.log_checks(_stdout(5) + "track killed\n", b"", 5, "alpha"),
        SS.log_checks(_stdout(5), b"", 5, "GenericIon"),
    ]
    assert all(bad), bad


# ---------------------------------------------------------------------------------- seal
def test_names_and_seeds() -> None:
    assert SS.bundle_name(22, 1) == "v5ss-he-water-150mevu-case22-part1"
    assert [SS.ionmc_seed(c) for c in SS.SS_CASES] == list(range(20631004, 20631705, 100))
    assert SS.topas_seed(31, 1) == 20619131
    assert SS.parse_bundle(SS.SMOKE_BUNDLE)["kind"] == "smoke"
    with pytest.raises(SS.RecordError):
        SS.parse_bundle("v5ss-he-water-150mevu-case23-part0")


def _sample(case: int, tmp: Path, seed: int = 7) -> tuple[dict[str, Any], Path]:
    rng = np.random.default_rng(seed)
    n = SS.SS_HISTORIES
    tx, ty = rng.normal(0, 1e-3, n), rng.normal(0, 1e-3, n)
    anchor = SS.windowed_rms(tx, ty, None)[0]
    npz = tmp / f"ionmc-ss-case{case}-k0.npz"
    np.savez_compressed(
        npz,
        theta_x=tx,
        theta_y=ty,
        z_exit_mm=np.ones(n),
        escape_history=np.arange(n),
        theta_rms_untruncated=np.float64(anchor),
    )
    seed_ = SS.ionmc_seed(case)
    doc = {
        "row": "v5ss",
        "case": case,
        "r_index": SS.SS_R_INDEX[case],
        "k": 0,
        "seed": seed_,
        "seed_base": 20631004,
        "histories": n,
        "frozen_histories": n,
        "reduced": False,
        "backend": "warp-cpu",
        "precision": "float64",
        "valid": True,
        "pass": True,
        "git_sha": SHA,
        "suite": "v5ss",
        "n_escaped": n,
        "npz": npz.name,
        "npz_sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
        "counters": {"stall": 0, "energy_inversion": 0},
        "theta_anchor_rad": anchor,
        "effective_config": {"cpu_workers": 1, "seed": seed_, "n_histories": n},
    }
    return doc, npz


@pytest.fixture(scope="module")
def sample22(tmp_path_factory: pytest.TempPathFactory) -> tuple[dict[str, Any], Path]:
    return _sample(22, tmp_path_factory.mktemp("s22"))


def test_case_entry_window_and_rejections(sample22: tuple[dict[str, Any], Path]) -> None:
    doc, npz = sample22
    e = SS.case_entry(22, doc, npz, SHA)
    assert e["window_rad"] == 2.5 * e["theta_anchor_rad"]
    assert e["theta_w_rad"] < e["theta_anchor_rad"] and 0.97 < e["p_window"] < 0.99
    for change in (
        {"seed": doc["seed"] + 1},
        {"histories": 399_960},
        {"npz_sha256": "0" * 64},
        {"git_sha": "2" * 40},
        {"pass": False},
        {"counters": {"stall": 1}},
        {"theta_anchor_rad": doc["theta_anchor_rad"] * 1.001},
        {"effective_config": {**doc["effective_config"], "cpu_workers": 2}},
    ):
        with pytest.raises(SS.RecordError):
            SS.case_entry(22, {**doc, **change}, npz, SHA)
    with pytest.raises(SS.RecordError):
        SS.case_entry(25, doc, npz, SHA)  # another case's seed


def test_build_seal_rejections(sample22: tuple[dict[str, Any], Path]) -> None:
    doc, npz = sample22
    e = SS.case_entry(22, doc, npz, SHA) | {"run_id": "RUN-x", "output_manifest_sha256": "m"}
    entries = [{**e, "case": c, "seed": SS.ionmc_seed(c)} for c in SS.SS_CASES]
    seal = SS.build_seal(entries, SHA)
    assert seal["schema"] == "ionmc-v8-ss-sample-seal-1" and len(seal["cases"]) == 8
    for bad in (
        entries[:7],
        [*entries, entries[0]],
        [{**entries[0], "seed": 1}, *entries[1:]],
        [{**entries[0], "histories": 1000}, *entries[1:]],
        [{**entries[0], "window_rad": 3 * entries[0]["theta_anchor_rad"]}, *entries[1:]],
    ):
        with pytest.raises(SS.RecordError):
            SS.build_seal(bad, SHA)
    with pytest.raises(SS.RecordError):
        SS.build_seal(entries, "abc")


def test_wrapper_render_round_trip() -> None:
    src = SS.render_wrapper(22, SHA, SS.NPZ_SUBDIR, "1")
    assert (
        src.split("\n")[1]
        == f"case, scale, expected_sha, subdir = 22, '1', '{SHA}', 'v5ss-samples'"
    )
    assert src.replace(src.split("\n")[1], SS.PLACEHOLDER) == SS.SS_WRAPPER_TEMPLATE
    assert hashlib.sha256(SS.SS_WRAPPER_TEMPLATE.encode()).hexdigest() == SS.SS_WRAPPER_SHA256
    assert SS.SS_WRAPPER_SHA256 == PINNED_WRAPPER_SHA256
    compile(src, "<wrapper>", "exec")
    assert "'ss-sample'" in src and "'--seed-base', '20631004'" in src
    for args in (
        (22.0, SHA, "v5ss-samples", "1"),
        (22, "abc", "v5ss-samples", "1"),
        (22, SHA, "../x", "1"),
        (22, SHA, "v5ss-samples", "1; rm"),
    ):
        with pytest.raises(ValueError):
            SS.render_wrapper(*args)


PINNED_WRAPPER_SHA256 = "021088be16f79292e3a20a420160eff6e463dd599842255114d8b946d01f06eb"


def _sign(record: dict[str, Any], signer: str) -> dict[str, Any]:
    at = IO.attest_module()
    old = at.SIGNERS
    at.SIGNERS = (*old, "reference-runner")
    try:
        return at.sign_record(record, at.keypair_from_seed(KEY_SEED), signer)
    finally:
        at.SIGNERS = old


def _trust() -> dict[str, Any]:
    kp = IO.attest_module().keypair_from_seed(KEY_SEED)
    return {
        "host-runner": {kp.key_id: kp.public_key.hex()},
        "reference-runner": {kp.key_id: kp.public_key.hex()},
    }


def _host_world(tmp: Path, sample: tuple[dict[str, Any], Path], argv_case: int = 22) -> dict:
    doc, npz = sample
    run_id = "RUN-20261114T000000Z-00000001"
    top = tmp / "validation" / "generated" / run_id / SS.NPZ_SUBDIR / "case22"
    top.mkdir(parents=True)
    shutil.copy(npz, top / npz.name)
    text = f"#JSON-BEGIN\n{json.dumps(doc, sort_keys=True)}\n#JSON-END\n"
    (top / "case22-k0.stdout.txt").write_text(text)
    (top / "case22-k0.stderr.txt").write_text("")
    digest, _ = IO.published_manifest(tmp, run_id)
    rec = {
        "schema_version": 1,
        "run_id": run_id,
        "task_id": "V3-008",
        "sha": SHA,
        "dirty_before": False,
        "dirty_after": False,
        "committed_state_unchanged": True,
        "argv": ["python", "-c", SS.render_wrapper(argv_case, SHA, SS.NPZ_SUBDIR, "1")],
        "started_at": "2026-11-14T00:00:00+00:00",
        "finished_at": "2026-11-14T00:10:00+00:00",
        "exit_code": 0,
        "timed_out": False,
        "succeeded": True,
        "output_publication_error": None,
        "stdout_sha256": "a" * 64,
        "output_manifest_sha256": digest,
    }
    return _sign(rec, "host-runner")


def test_verify_sample_run(tmp_path: Path, sample22: tuple[dict[str, Any], Path]) -> None:
    rec = _host_world(tmp_path / "ok", sample22)
    e = SS.verify_sample_run(rec, _trust(), SHA, tmp_path / "ok", "d" * 64)
    assert e["case"] == 22 and e["run_id"] == rec["run_id"] and e["host_record_sha256"] == "d" * 64
    with pytest.raises(SS.RecordError):  # wrong implementation sha
        SS.verify_sample_run(rec, _trust(), "2" * 40, tmp_path / "ok", "d")
    with pytest.raises(SS.RecordError):  # tampered after signing
        SS.verify_sample_run({**rec, "exit_code": 1}, _trust(), SHA, tmp_path / "ok", "d")
    bad = _host_world(tmp_path / "argv", sample22, argv_case=25)  # argv of another case
    with pytest.raises(SS.RecordError):
        SS.verify_sample_run(bad, _trust(), SHA, tmp_path / "argv", "d")
    rec2 = _host_world(tmp_path / "extra", sample22)
    (tmp_path / "extra" / "validation" / "generated" / rec2["run_id"] / "x.txt").write_text("x")
    with pytest.raises(SS.RecordError):
        SS.verify_sample_run(rec2, _trust(), SHA, tmp_path / "extra", "d")


# ---------------------------------------------------------------------------------- plan note
def test_ss_note_parsing() -> None:
    note = f'```json\n{{"ss_seal_sha": "{SHA}", "ss_seal_review_id": "REVIEW-s"}}\n```\n'
    assert SS.parse_ss_note("x\n" + note) == {"seal_sha": SHA, "seal_review_id": "REVIEW-s"}
    for text in (
        "no note",
        note + note,
        note.replace("REVIEW-s", 'REVIEW-s", "x": "1'),
        note.replace(SHA, "abc"),
    ):
        with pytest.raises(SS.RecordError):
            SS.parse_ss_note(text)
    plan = (ROOT / "validation" / "plans" / "v3-008-acceptance.md").read_text()
    committed = SS.parse_ss_note(plan)  # the dated SS seal note of 2026-11-13
    assert committed == {
        "seal_sha": "7fe6ab1d7b7eb18ab377156e0551057264fba40b",
        "seal_review_id": "REVIEW-5d9e8df442e94669aef5fe609794b608",
    }
    assert IO.parse_seal_note(plan)["seal_sha"].startswith("0c68b78a")  # r2 note unaffected


class _Reader:
    def __init__(self, doc: Any) -> None:
        self.doc = doc

    def json(self, path: Any, what: str) -> Any:
        return self.doc


def test_index_fails_closed() -> None:
    runs = [{"case": c, "part": p, "run_id": f"REF-{c}-{p}"} for c in SS.SS_CASES for p in (0, 1)]
    assert len(SS.load_index(_Reader({"runs": runs}))) == 16
    for bad in (
        runs[:15],
        [*runs[:15], runs[0]],
        [*runs[:15], {**runs[15], "case": 23}],
        [*runs[:15], {**runs[15], "part": 2}],
    ):
        with pytest.raises(SS.RecordError):
            SS.load_index(_Reader({"runs": bad}))


# ---------------------------------------------------------------------------------- check-run
def _git(repo: Path, *args: str) -> str:
    env = {
        "GIT_AUTHOR_NAME": "t",
        "GIT_AUTHOR_EMAIL": "t@x",
        "GIT_COMMITTER_NAME": "t",
        "GIT_COMMITTER_EMAIL": "t@x",
        "PATH": "/usr/bin:/bin:/usr/local/bin",
        "HOME": str(repo),
    }
    return subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=True, env=env
    ).stdout.strip()


def _commit(repo: Path, msg: str) -> str:
    _git(repo, "add", "-A")
    _git(repo, "-c", "commit.gpgsign=false", "commit", "-q", "-m", msg)
    return _git(repo, "rev-parse", "HEAD")


def _world(
    tmp: Path,
    *,
    started: str = "2026-11-15T00:00:00+00:00",
    seal_first: bool = True,
    phsp_edit: Any = None,
    bundle_case_edit: Any = None,
    manifest_late_edit: bool = False,
    stdout: str | None = None,
) -> dict[str, Any]:
    repo = tmp / "repo"
    repo.mkdir(parents=True)
    _git(repo, "init", "-q", "-b", "main")
    case, part = 28, 0
    bundle = SS.bundle_name(case, part)
    src = ROOT / "validation" / "reference_cases" / "topas" / bundle
    dst = repo / IO.BUNDLE_ROOT / bundle
    shutil.copytree(src, dst)
    if bundle_case_edit:
        cj = json.loads((dst / "case.json").read_text())
        (dst / "case.json").write_text(json.dumps(bundle_case_edit(cj)))
    manifest = {
        "schema": SS.MANIFEST_SCHEMA,
        "bundle_set_digest": "x",
        "bundles": {
            bundle: {
                f: hashlib.sha256((dst / f).read_bytes()).hexdigest()
                for f in ("case.json", "input.txt")
            }
        },
    }
    (repo / SS.MANIFEST_PATH).write_text(json.dumps(manifest))
    kp = IO.attest_module().keypair_from_seed(KEY_SEED)
    trust = repo / IO.DEFAULT_TRUST
    trust.parent.mkdir(parents=True)
    key = [{"key_id": kp.key_id, "public_key": kp.public_key.hex()}]
    trust.write_text(json.dumps({"reference-runner": key, "host-runner": key}))
    (repo / "validation" / "plans").mkdir(parents=True)
    plan = repo / IO.PLAN_PATH
    plan.write_text("# plan\n")
    impl = _commit(repo, "implementation")
    (repo / SS.SS_DIR).mkdir(parents=True, exist_ok=True)
    (repo / SS.SEAL_PATH).write_text("{}")
    seal_sha = _commit(repo, "seal")
    review = {
        "review_id": "REVIEW-s",
        "head_sha": seal_sha,
        "status": "passed",
        "finished_at": "2026-11-14T00:00:00+00:00",
        "report": {"verdict": "pass"},
    }
    (repo / SS.REVIEW_PATH).write_text(json.dumps(review))
    note = json.dumps({"ss_seal_sha": seal_sha, "ss_seal_review_id": "REVIEW-s"})
    plan.write_text(f"# plan\n\n```json\n{note}\n```\n")
    code_sha = _commit(repo, "seal note") if seal_first else impl
    outputs = tmp / "outputs"
    files = {}
    for name in ("case.json", "input.txt"):
        blob = (repo / IO.BUNDLE_ROOT / bundle / name).read_bytes()
        files[name] = {"sha256": hashlib.sha256(blob).hexdigest(), "bytes": len(blob)}
    case_doc = json.loads((repo / IO.BUNDLE_ROOT / bundle / "case.json").read_text())
    request = {
        "schema_version": 2,
        "code_sha": code_sha,
        "input_files": files,
        "case": case_doc,
        "dirty": False,
        "engine": {"engine": "topas"},
    }
    cfg_id = IO.canonical_digest(request)
    run_id = f"REF-{cfg_id[:20]}-0000abcd"
    odir = outputs / run_id
    (odir / "work").mkdir(parents=True)
    (odir / "request.json").write_text(json.dumps(request))
    (odir / "stdout.txt").write_text(stdout if stdout is not None else _stdout(30000))
    (odir / "stderr.txt").write_text("")
    h, _ = _fixture()
    h = h.replace(": 3\n", ": 30000\n")
    rng = np.random.default_rng(9)
    u, v = rng.normal(0, 1e-3, 30000), rng.normal(0, 1e-3, 30000)
    rows = "".join(
        f"0 0 0.1 {float(a)!r} {float(b)!r} 999 1 1000020040 0 1\n"
        for a, b in zip(u, v, strict=True)
    )
    if phsp_edit:
        h, rows = phsp_edit(h, rows)
    (odir / "work" / "ss_exit.header").write_text(h)
    (odir / "work" / "ss_exit.phsp").write_text(rows)
    listing = [
        {
            "path": p.relative_to(odir).as_posix(),
            "bytes": p.stat().st_size,
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        }
        for p in sorted(odir.rglob("*"))
        if p.is_file()
    ]
    sha = {e["path"]: e["sha256"] for e in listing}
    work = {
        e["path"][5:]: {"sha256": e["sha256"], "bytes": e["bytes"]}
        for e in listing
        if e["path"].startswith("work/")
    }
    result = _sign(
        {
            "schema_version": 3,
            "run_id": run_id,
            "configuration_id": cfg_id,
            "code_sha": code_sha,
            "task_id": None,
            "engine": "topas",
            "request_sha256": sha["request.json"],
            "started_at": started,
            "finished_at": "2026-11-15T01:00:00+00:00",
            "seconds": 600.0,
            "exit_code": 0,
            "timed_out": False,
            "missing_outputs": [],
            "succeeded": True,
            "output_count": len(work),
            "output_files_sha256": IO.canonical_digest(work),
            "stdout_sha256": sha["stdout.txt"],
            "stderr_sha256": sha["stderr.txt"],
        },
        "reference-runner",
    )
    rdir = repo / SS.TOPAS_RECORDS / run_id
    rdir.mkdir(parents=True)
    (rdir / "record.json").write_text(json.dumps(result))
    rb = (rdir / "record.json").read_bytes()
    listing.append(
        {"path": "record.json", "sha256": hashlib.sha256(rb).hexdigest(), "bytes": len(rb)}
    )
    (rdir / "artifacts.json").write_text(json.dumps({"run_id": run_id, "files": listing}))
    if manifest_late_edit:
        manifest["bundles"][bundle]["input.txt"] = "0" * 64
        (repo / SS.MANIFEST_PATH).write_text(json.dumps(manifest))
    _commit(repo, "records")
    return {"repo": repo, "outputs": outputs, "run_id": run_id, "bundle": bundle}


def _check(w: dict[str, Any]) -> dict[str, Any]:
    return SS.check_run(w["run_id"], w["bundle"], None, w["outputs"], repo=w["repo"])


def test_check_run_valid(tmp_path: Path) -> None:
    rep = _check(_world(tmp_path))
    assert rep["verdict"] == "VALID", rep["problems"]
    assert rep["chronology"]["ancestry"] and rep["chronology"]["review_before_start"]
    SS._json_safe(rep)  # no arrays in the report


def test_check_run_blinding_failures(tmp_path: Path) -> None:
    rep = _check(_world(tmp_path / "a", started="2026-11-13T00:00:00+00:00"))
    assert rep["verdict"] == "BLINDING FAIL" and not rep["chronology"]["review_before_start"]


def test_check_run_before_seal_is_error(tmp_path: Path) -> None:
    w = _world(tmp_path, seal_first=False)  # code_sha predates the seal note and review
    with pytest.raises(IO.RecordError):
        _check(w)


@pytest.mark.parametrize(
    "edit",
    [
        lambda h, r: (h, r.replace("1000020040", "1000060120", 1)),
        lambda h, r: (h, "\n".join(r.splitlines()[1:]) + "\n"),
        lambda h, r: (h.replace("Weight", "Weight X"), r),
    ],
)
def test_check_run_phase_space_problems(tmp_path: Path, edit: Any) -> None:
    rep = _check(_world(tmp_path, phsp_edit=edit))
    assert rep["verdict"] == "RUN-VALIDITY FAIL" and rep["problems"]


def test_check_run_log_problem(tmp_path: Path) -> None:
    rep = _check(_world(tmp_path, stdout=_stdout(30000, model="eCoulombScattering")))
    assert rep["verdict"] == "RUN-VALIDITY FAIL"


@pytest.mark.parametrize(
    "kw",
    [
        {"bundle_case_edit": lambda c: {**c, "seeds": [20619128]}},
        {"bundle_case_edit": lambda c: {**c, "histories": 29999}},
        {"manifest_late_edit": True},
    ],
)
def test_check_run_record_errors(tmp_path: Path, kw: dict[str, Any]) -> None:
    with pytest.raises(IO.RecordError):
        _check(_world(tmp_path, **kw))


def test_cli_error_report(tmp_path: Path) -> None:
    out = tmp_path / "r.json"
    code = SS.main(["check-run", "--run-id", "REF-x", "--bundle", "nope", "--out", str(out)])
    assert code == 2 and json.loads(out.read_text())["status"] == "ERROR"


def test_grade_fails_closed_without_the_run_index(tmp_path: Path) -> None:
    w = _world(tmp_path)
    with pytest.raises(IO.RecordError):
        SS.grade(w["outputs"], tmp_path, repo=w["repo"])
