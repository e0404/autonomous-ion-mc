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


# =============================================================================================
# TD/TL statistics, loaders and CLI (synthetic data only: no TOPAS output, no physics reference
# value). The TOPAS file formats are written from the header lines of the TOPAS CSV/binary
# documentation as parsed by ``ionmc.reference.parsers``; no r2 TOPAS output was materialized
# when these tests were written, and no number of any materialized output is used.
# =============================================================================================

import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402

import numpy as np  # noqa: E402

ST, IO = R2.st, R2.io


def _phi(x: np.ndarray) -> np.ndarray:
    return 0.5 * np.array([math.erfc(-v / math.sqrt(2.0)) for v in np.ravel(x)]).reshape(
        np.shape(x)
    )


def _gauss_bins(edges: np.ndarray, sigma: float, n: float = 1.0) -> np.ndarray:
    return n * np.diff(_phi(edges / sigma))


def test_tost_pass_fail_inconclusive() -> None:
    assert ST.tost(10.0, 10.02, 0.01, 0.01, 0.05)["verdict"] == "PASS"
    inc = ST.tost(10.0, 10.04, 0.01, 0.01, 0.05)  # |d| <= tol but margin exceeds it
    assert inc["verdict"] == "FAIL" and inc["note"] == "inconclusive"
    out = ST.tost(10.0, 10.1, 0.0, 0.0, 0.05)
    assert out["verdict"] == "FAIL" and out["note"] == "outside tolerance"
    assert ST.tost(math.nan, 1.0, 0.0, 0.0, 1.0)["verdict"] == "FAIL"
    r = ST.tost(1.0, 1.0, 0.03, 0.04, 1.0)
    assert r["sigma_comb"] == pytest.approx(0.05) and r["statistic"] == pytest.approx(1.645 * 0.05)
    assert ST.band(1.1, 0.01, -0.05, 0.40)["verdict"] == "PASS"
    assert ST.band(0.96, 0.01, -0.05, 0.40)["note"] == "inconclusive"
    assert ST.band(1.5, 0.0, -0.05, 0.40)["verdict"] == "FAIL"


def test_uncertainties() -> None:
    v = [1.0, 2.0, 3.0, 4.0]
    assert ST.topas_se(v) == pytest.approx(1.430 * np.std(v, ddof=1) / 2.0)
    assert ST.batch_se(v) == pytest.approx(np.std(v, ddof=1) / 2.0)
    for bad in ([1.0, 2.0, 3.0], [1.0] * 6, [1.0, 2.0, 3.0, math.nan]):
        with pytest.raises(ST.StatisticsError):
            ST.topas_se(bad)


def test_peak_parabola_recovers_vertex() -> None:
    hz, z0 = 0.02, 3.05  # z0 at a group centre: group means of a quadratic are a quadratic
    zc = (np.arange(500) + 0.5) * hz
    e = 10.0 - 0.5 * (zc - z0) ** 2
    zp, dp = ST.peak(e, hz)
    # mean of 5 bin-centre samples (offsets -2h..2h) of 10 - 0.5 u^2: 10 - 0.5 * 2 h^2
    assert zp == pytest.approx(z0, abs=1e-9) and dp == pytest.approx(10.0 - hz * hz)
    zp2, _ = ST.peak(10.0 - 0.5 * (zc - 3.0) ** 2, hz)  # off-centre vertex: exact as well
    assert zp2 == pytest.approx(3.0, abs=1e-9)
    with pytest.raises(ST.StatisticsError):
        ST.peak(zc, hz)  # monotonic: maximal group at the end


def test_distal_levels_on_a_linear_falloff() -> None:
    hz = 0.05
    zc = (np.arange(400) + 0.5) * hz
    e = np.where(zc < 10.0, 1.0, np.clip(1.0 - (zc - 10.0) / 4.0, 0.0, None))
    assert ST.distal_level(e, hz, 9.0, 0.8) == pytest.approx(10.8)
    assert ST.distal_level(e, hz, 9.0, 0.2) == pytest.approx(13.2)
    with pytest.raises(ST.StatisticsError):
        ST.distal_level(np.ones(400), hz, 9.0, 0.8)


def _bragg(zc: np.ndarray, shift: float = 0.0) -> np.ndarray:
    z = zc - shift
    tail = 0.5 * np.array([math.erfc((v - 7.0) / (math.sqrt(2.0) * 0.15)) for v in z])
    return (1.0 + 0.05 * z) * tail + 2.0 * np.exp(-((z - 6.9) ** 2) / (2 * 0.25**2)) * tail


def test_depth_quantities_on_synthetic_bragg_curve() -> None:
    hz = 0.05
    zc = (np.arange(200) + 0.5) * hz
    q = ST.depth_quantities(_bragg(zc), hz, r_plan_mm=22.5)
    assert 6.5 < q["z_peak"] < q["r80"] < q["r20"] < 7.6
    assert q["w"] == pytest.approx(q["r20"] - q["r80"])
    assert q["d_ent"] == pytest.approx(np.mean(_bragg(zc)[10:30]))  # [0.5, 1.5] mm
    assert q["pe"] == pytest.approx(q["d_peak"] / q["d_ent"])
    g = ST.plateau_groups(_bragg(zc), hz, 22.5, q["r80"])
    assert g.size == int((0.8 * q["r80"] - 0.5) // 1.0)  # 1 mm groups from 0.5 mm
    assert g[0] == pytest.approx(np.mean(_bragg(zc)[10:30]))
    shifted = ST.depth_quantities(_bragg(zc, 0.3), hz, 22.5)
    assert shifted["r80"] - q["r80"] == pytest.approx(0.3, abs=0.01)
    assert ST.range_class(159.3)["group_mm"] == 5.0 and ST.range_class(53.1)["group_mm"] == 2.0


def test_core_sigma_recovers_gaussian() -> None:
    edges = -60.0 + 0.02 * np.arange(6001)
    for sigma in (0.07, 0.5, 3.0):
        fit = ST.core_sigma(_gauss_bins(edges, sigma, 1e6), edges)
        assert fit.sigma == pytest.approx(sigma, rel=1e-3)
    # robust to a broad low pedestal (the halo): the +-20 mm moment is not
    prof = _gauss_bins(edges, 0.5, 1e6) + _gauss_bins(edges, 8.0, 5e4)
    assert ST.core_sigma(prof, edges).sigma == pytest.approx(0.5, rel=0.01)
    assert ST.second_moment_sigma(prof, edges) > 0.7
    assert ST.halo_fraction(prof, edges, 0.5) > 0.04
    rng = np.random.default_rng(5)
    noisy = rng.poisson(_gauss_bins(edges, 0.4, 2e5)).astype(float)
    assert ST.core_sigma(noisy, edges).sigma == pytest.approx(0.4, rel=0.02)
    with pytest.raises(ST.StatisticsError):  # flat profile: no finite ML scale in the window
        ST.core_sigma(np.ones(6000), edges)


def test_chunk_combination_equals_single_sample() -> None:
    rng = np.random.default_rng(1)
    sums = [rng.random((40, 7)) * h for h in (250, 250, 500)]  # raw per-batch sums of 3 chunks
    hpb = [250, 250, 500]  # histories per batch
    combined = IO.combine_chunks([(h, s / h) for h, s in zip(hpb, sums, strict=True)])
    single = sum(sums) / sum(hpb)  # one run with the same 40 batches of 1000 histories
    np.testing.assert_allclose(combined, single, rtol=1e-14)
    with pytest.raises(IO.RecordError):
        IO.combine_chunks([(1, np.zeros(3)), (1, np.zeros(4))])


def test_verify_npz_fails_closed(tmp_path: Path) -> None:
    p = tmp_path / "a.npz"
    np.savez_compressed(p, a=np.zeros(3))
    good = {"npz_sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "npz_bytes": p.stat().st_size}
    IO.verify_npz(p, good)
    with pytest.raises(IO.RecordError, match="sha256"):
        IO.verify_npz(p, {**good, "npz_sha256": "0" * 64})
    with pytest.raises(IO.RecordError, match="size"):
        IO.verify_npz(p, {**good, "npz_bytes": 1})
    with pytest.raises(IO.RecordError):
        IO.verify_npz(tmp_path / "missing.npz", good)


# ---- synthetic evidence world (git repo, records, outputs, samples) --------------------------

HZ, HX, L = 0.05, 0.1, 10.0
T_LATE, T_EARLY = "2026-11-09T12:00:00+00:00", "2026-11-09T08:59:00Z"  # review ends 09:00
SMALL = {1: ("he", 50.0, 22.5, L, HZ, HX, "main"), 3: ("he", 250.0, 22.5, L, HZ, HX, "main")}
J_MEV = 1.602176634e-13
SIG = [0.5 + 0.05 * k for k in range(int(L))]  # lateral sigma per 1 mm slab (>= 3 h_x: gated)
OUTS = {
    "idd-em": ["idd_dose.csv"],
    "idd-full": ["idd_dose.csv", "idd_primary_dose.csv", "primary_fluence.csv"],
    "lat-em": ["lat_x.bin", "lat_x.binheader", "lat_y.bin", "lat_y.binheader"],
    "lat-full": ["lat_x.bin", "lat_x.binheader", "lat_y.bin", "lat_y.binheader",
                 "lat_x_primary.bin", "lat_x_primary.binheader", "lat_y_primary.bin",
                 "lat_y_primary.binheader"],
}  # fmt: skip
# Synthetic printout in the Geant4 11.4 / OpenTOPAS 4.3 layout (format only, checked against the
# materialized r2 smoke logs; no number of those logs is used).
MSC = "msc:  for {p}  SubType= 10\n      ===== EM models for the G4Region  DefaultRegionForTheWorld ======\n            {m} : Emin=    0 eV  Emax=   10 GeV\n          StepLim=Minimal Rfact=0.2 Gfact=2.5\n\n"  # noqa: E501
REGISTERED = (  # TOPAS Ph/ListProcesses layout, as in the materialized EM smoke logs (names only)
    "Registered Physics Processes:\n     Transportation,                msc,              eIoni\n"
    "{extra}            ionIoni,              hIoni\n\n"
)
STDOUT_EM = (
    REGISTERED.format(extra="")
    + "Use ICRU90 data                                    1\n"
    + MSC.format(p="e-", m="GoudsmitSaunderson")
    + MSC.format(p="GenericIon", m="UrbanMsc") + MSC.format(p="alpha", m="UrbanMsc")
    + "G4WT0 > Begin processing for Run: 0, History: 0\n"
    + "Particle source Beam: Total number of histories: 1000\n"
)  # fmt: skip
SEP = "-----------------------------------------------------------------------\n"
STDOUT_FULL = (
    REGISTERED.format(extra="              Decay,     Radioactivation,\n")
    + " hInelastic QGSP_BIC_HP  Thresholds: \n"
    "======          Radioactive Decay Physics Parameters           =======\n"
    "Threshold for very long decay time at rest        1 y  \n\n"
    + STDOUT_EM
    + "                  HADRONIC PROCESSES SUMMARY (verbose level 1)\n" + SEP
    + "                           Hadronic Processes for GenericIon\n"
    + "  Process: ionInelastic\n"
    + "        Model:                  QMDModel: 100 MeV/n ---> 6 GeV/n\n"
    + "  Process: RadioactiveDecay\n" + SEP
    + "                           Hadronic Processes for mu-\n"
    + "  Process: muMinusCaptureAtRest\n" + SEP
    + "                           Hadronic Processes for neutron\n"
    + "  Process: hadElastic\n"
    + "        Model:          NeutronHPElastic: 0 eV  ---> 20 MeV\n" + SEP
)  # fmt: skip
MODULES = {
    "em": 'sv:Ph/Default/Modules = 1 "g4em-standard_opt4"\n',
    "full": 'sv:Ph/Default/Modules = 7 "g4em-standard_opt4" "g4h-phy_QGSP_BIC_HP" "g4decay" '
    '"g4h-elastic_HP" "g4stopping" "g4ion-QMD" "g4radioactivedecay"\n',
}
KEY_SEED = bytes(range(32))


def _attest() -> Any:
    return IO.attest_module()


def _sign(record: dict[str, Any], signer: str = "reference-runner") -> dict[str, Any]:
    at = _attest()
    old = at.SIGNERS
    at.SIGNERS = (*old, "reference-runner")
    try:
        return at.sign_record(record, at.keypair_from_seed(KEY_SEED), signer)
    finally:
        at.SIGNERS = old


def _resign(path: Path, changes: dict[str, Any]) -> None:
    """Re-sign a record.json with ``changes`` and keep its artifact-listing entry consistent."""
    rec = json.loads(path.read_text())
    rec.pop("attestation")
    path.write_text(json.dumps(_sign(rec | changes)))
    _relist(path)


def _relist(path: Path) -> None:
    lst = path.parent / "artifacts.json"
    doc = json.loads(lst.read_text())
    for e in doc["files"]:
        if e["path"] == "record.json":
            e["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            e["bytes"] = path.stat().st_size
    lst.write_text(json.dumps(doc))


def _csv(energy_total: np.ndarray, scorer: str = "IDD") -> str:
    mass = 120.0 * 120.0 * HZ * 1e-6
    dose = energy_total * J_MEV / mass
    head = (
        f"# TOPAS Version: 4.3\n# Parameter File: input.txt\n# Results for scorer: {scorer}\n"
        "# Scored in component: Water\n# X in 1 bin  of 12 cm\n# Y in 1 bin  of 12 cm\n"
        f"# Z in {dose.size} bins of {HZ / 10} cm\n"
        "# DoseToMedium ( Gy ) : Sum   Standard_Deviation\n"
    )
    return head + "".join(f"0, 0, {i}, {float(v)!r}, 0\n" for i, v in enumerate(dose))


def _binheader(stem: str, axis: int, n: int) -> str:
    """The real r2 .binheader layout (He 250 LAT-EM smoke, format only)."""
    ax = [f"# X in {n} bins of {HX / 10} cm", "# Y in 1 bin  of 12 cm"]
    if axis == 1:
        ax = ["# X in 1 bin  of 12 cm", f"# Y in {n} bins of {HX / 10} cm"]
    return "\n".join([
        "# TOPAS Version: 4.3", "# Parameter File: input.txt",
        f"# Results for scorer: {IO.SCORER_OF[stem]}", "# Scored in component: Water", *ax,
        f"# Z in {int(L)} bins of 0.1 cm", "# DoseToMedium ( Gy ) : Sum",
        f"# Binary file: {stem}.bin", "",
    ])  # fmt: skip


def _commit(repo: Path, msg: str) -> str:
    _git(repo, "add", "-A")
    _git(repo, "-c", "commit.gpgsign=false", "commit", "-q", "-m", msg)
    return _git(repo, "rev-parse", "HEAD")


def _world(
    tmp: Path, *, cases: dict[int, dict[str, Any]] | None = None,
    configs: tuple[str, ...] = ("idd-em", "lat-em"), started_ok: bool = True,
    tamper: dict[str, Any] | None = None,
) -> dict[str, Any]:  # fmt: skip
    """A complete synthetic evidence world. ``cases``: case -> {shift (TOPAS EM IDD), full_shift,
    sig (TOPAS sigma per slab), ion_scale (ionmc sigma factor)}; ``tamper``: run-validity faults
    (stderr / stdout callable / input_extra) applied to part 1 of the first config of case 1."""
    cases = cases or {1: {}}
    tamper = tamper or {}
    rng = np.random.default_rng(11)
    repo = tmp / "repo"
    repo.mkdir(parents=True)
    _git(repo, "init", "-q", "-b", "main")
    nz, nl = round(L / HZ), round(120.0 / HX)
    zc = (np.arange(nz) + 0.5) * HZ
    ledges = -60.0 + HX * np.arange(nl + 1)
    first = (1, configs[0], IO.CONFIG_PARTS[configs[0]][1])

    def lat_of(sig: list[float]) -> np.ndarray:
        return np.stack([_gauss_bins(ledges, s) for s in sig], axis=1)  # (nl, L) per primary

    for case in cases:  # bundles (committed before the seal)
        for cfg in configs:
            for part in IO.CONFIG_PARTS[cfg]:
                d = repo / IO.BUNDLE_ROOT / IO.bundle_name(case, cfg, part)
                d.mkdir(parents=True)
                seed = IO.topas_seed(case, part)
                case_doc = {"histories": 1000, "seeds": [seed], "expected_outputs": OUTS[cfg],
                            "input": "input.txt", "rationale": f"part {part}"}  # fmt: skip
                (d / "case.json").write_text(json.dumps(case_doc))
                extra = tamper.get("input_extra", "") if (case, cfg, part) == first else ""
                mods = MODULES["full" if cfg.endswith("full") else "em"]
                if (case, cfg, part) == first and "modules" in tamper:
                    mods = tamper["modules"]
                (d / "input.txt").write_text(
                    f"# case {case} part {part}\ni:Ts/Seed = {seed}\n"
                    f"i:So/Beam/NumberOfHistoriesInRun = 1000\n{mods}{extra}"
                )
    samples = tmp / "samples"  # ionmc chunks + seal
    chunks = []
    for case, spec in cases.items():
        sig_i = [s * spec.get("ion_scale", 1.0) for s in spec.get("sig", SIG)]
        for k in (1, 2):
            run_id = f"RUN-2026110{k}T00000{case}Z-0000000{k}"
            out = samples / "validation" / "generated" / run_id / "r2-samples" / f"case{case}"
            out.mkdir(parents=True)
            idd = (_bragg(zc) * HZ)[None, :] * (1.0 + 1e-3 * rng.standard_normal((40, 1)))
            lx = lat_of(sig_i)[None] * (1.0 + 1e-3 * rng.standard_normal((40, 1, 1)))
            arrays = {
                "idd_batch_energy_mev": idd[:, None, None, :],
                "idd_edges0_mm": np.array([-60.0, 60.0]), "idd_edges1_mm": np.array([-60.0, 60.0]),
                "idd_edges2_mm": HZ * np.arange(nz + 1),
                "latx_batch_energy_mev": lx[:, :, None, :], "latx_edges0_mm": ledges,
                "latx_edges1_mm": np.array([-60.0, 60.0]),
                "latx_edges2_mm": np.arange(int(L) + 1.0),
                "laty_batch_energy_mev": lx[:, None, :, :],
                "laty_edges0_mm": np.array([-60.0, 60.0]), "laty_edges1_mm": ledges,
                "laty_edges2_mm": np.arange(int(L) + 1.0),
            }  # fmt: skip
            npz = out / f"ionmc-case{case}-k{k}.npz"
            np.savez_compressed(npz, **arrays)
            chunks.append({
                "case": case, "k": k, "seed": IO.LV8_BASE + 100 * (29 + case) + k,
                "histories": 500000, "run_id": run_id, "npz": npz.name,
                "npz_bytes": npz.stat().st_size,
                "npz_sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
                "counters_all_zero": True,
            })  # fmt: skip
    seal_rec = {"schema": IO.SEAL_SCHEMA, "chunk_plan": {"n": {str(c): 2 for c in cases}},
                "chunks": chunks}  # fmt: skip
    seal_file = repo / IO.SEAL_RECORD_PATH
    seal_file.parent.mkdir(parents=True)
    seal_file.write_text(json.dumps(seal_rec))
    seal_sha = _commit(repo, "seal")
    (repo / "later.txt").write_text("after the seal")
    kp = _attest().keypair_from_seed(KEY_SEED)
    trust = repo / IO.DEFAULT_TRUST  # committed: the evaluator reads it from HEAD
    trust.parent.mkdir(parents=True)
    trust.write_text(json.dumps({"reference-runner": [
        {"key_id": kp.key_id, "public_key": kp.public_key.hex()}]}))  # fmt: skip
    code_sha = _commit(repo, "after")
    records, outputs = tmp / "records", tmp / "outputs"  # TOPAS records and outputs
    index, ids = [], {}
    for case, spec in cases.items():
        sig_t = spec.get("sig", SIG)
        for cfg in configs:
            for part in IO.CONFIG_PARTS[cfg]:
                bdir = f"{IO.BUNDLE_ROOT}/{IO.bundle_name(case, cfg, part)}"
                files = {}
                for name in ("case.json", "input.txt"):
                    blob = (repo / bdir / name).read_bytes()
                    files[name] = {"sha256": hashlib.sha256(blob).hexdigest(), "bytes": len(blob)}
                request = {
                    "schema_version": 2, "code_sha": code_sha, "input_files": files,
                    "case": json.loads((repo / bdir / "case.json").read_text()), "dirty": False,
                    "engine": {"engine": "topas"}, "os_release_sha256": "1" * 64,
                    "runner_source_sha256": "2" * 64, "sandbox_binary_sha256": "3" * 64,
                    "gpu": False,
                }  # fmt: skip
                cfg_id = IO.canonical_digest(request)
                run_id = f"REF-{cfg_id[:20]}-{case:02x}{part:02x}0000"  # runner's id derivation
                ids[(case, cfg, part)] = run_id
                odir = outputs / run_id
                (odir / "work").mkdir(parents=True)
                (odir / "request.json").write_text(json.dumps(request, indent=2))
                hit = (case, cfg, part) == first
                base_out = STDOUT_FULL if cfg.endswith("full") else STDOUT_EM
                stdout = tamper["stdout"](base_out) if hit and "stdout" in tamper else base_out
                (odir / "stdout.txt").write_text(stdout)
                (odir / "stderr.txt").write_text(tamper.get("stderr", "") if hit else "")
                pn = 1000 * (1.0 + 1e-3 * rng.standard_normal())
                w = odir / "work"
                if cfg.startswith("idd"):
                    sh = spec.get("shift", 0.0) if cfg == "idd-em" else spec.get("full_shift", 0.0)
                    (w / "idd_dose.csv").write_text(_csv(_bragg(zc, sh) * HZ * pn))
                    if cfg == "idd-full":
                        (w / "idd_primary_dose.csv").write_text(
                            _csv(0.9 * _bragg(zc, sh) * HZ * pn, "IDDPrimary")
                        )
                        (w / "primary_fluence.csv").write_text("# synthetic, never read\n")
                else:
                    profiles = {"lat_x": lat_of(sig_t), "lat_y": lat_of(sig_t)}
                    if cfg == "lat-full":
                        broad = lat_of([1.3 * s for s in sig_t])
                        profiles = {"lat_x": broad, "lat_y": broad,
                                    "lat_x_primary": lat_of(sig_t),
                                    "lat_y_primary": lat_of(sig_t)}  # fmt: skip
                    for stem, prof in profiles.items():
                        axis = 1 if stem.startswith("lat_y") else 0
                        m = (prof * pn * J_MEV / (HX * 120 * 1e-6)).astype("<f8")
                        shape = (nl, 1, int(L)) if axis == 0 else (1, nl, int(L))
                        (w / f"{stem}.bin").write_bytes(m.reshape(shape).tobytes("F"))
                        (w / f"{stem}.binheader").write_text(_binheader(stem, axis, nl))
                listing = [
                    {"path": p.relative_to(odir).as_posix(), "bytes": p.stat().st_size,
                     "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                    for p in sorted(odir.rglob("*")) if p.is_file()
                ]  # fmt: skip
                rdir = records / "runs" / run_id
                rdir.mkdir(parents=True)
                sha = {e["path"]: e["sha256"] for e in listing}
                work_inv = {  # the runner's common.digest(inventory(run_dir / "work"))
                    e["path"][len("work/") :]: {"sha256": e["sha256"], "bytes": e["bytes"]}
                    for e in listing if e["path"].startswith("work/")
                }  # fmt: skip
                result = {  # V3-INFRA-004 signed run record, schema_version 3
                    "schema_version": 3, "run_id": run_id, "configuration_id": cfg_id,
                    "code_sha": code_sha, "task_id": None, "engine": "topas",
                    "request_sha256": sha["request.json"],
                    "started_at": T_LATE if started_ok else T_EARLY,
                    "finished_at": "2026-11-09T13:00:00+00:00", "seconds": 100.0,
                    "exit_code": 0, "timed_out": False, "error": None, "missing_outputs": [],
                    "succeeded": True, "scientifically_validated": False,
                    "metadata_status": "declared; inspect native logs/outputs",
                    "output_count": len(work_inv),
                    "output_files_sha256": IO.canonical_digest(work_inv),
                    "stdout": stdout[-200:], "stderr": "",
                    "stdout_sha256": sha["stdout.txt"], "stderr_sha256": sha["stderr.txt"],
                }  # fmt: skip
                if not (hit and tamper.get("unsigned")):
                    result = _sign(result, tamper.get("signer", "reference-runner") if hit
                                   else "reference-runner")  # fmt: skip
                if hit and "started_after_signing" in tamper:
                    result["started_at"] = tamper["started_after_signing"]
                (rdir / "record.json").write_text(json.dumps(result))
                listing.append({"path": "record.json", "sha256": "", "bytes": 0})
                (rdir / "artifacts.json").write_text(
                    json.dumps({"run_id": run_id, "files": listing})
                )
                _relist(rdir / "record.json")
                index.append({"case": case, "config": cfg, "part": part, "run_id": run_id})
    (records / "index.json").write_text(json.dumps({"runs": index}))
    review = {"review_id": "REVIEW-x", "head_sha": seal_sha, "status": "passed",
              "finished_at": "2026-11-09T09:00:00Z", "report": {"verdict": "pass"}}  # fmt: skip
    (records / "review.json").write_text(json.dumps(review))
    plan = tmp / "plan.json"
    plan.write_text(json.dumps({"seal_sha": seal_sha, "seal_review_id": "REVIEW-x"}))
    return {"repo": repo, "records": records, "outputs": outputs, "samples": samples,
            "plan": plan, "seal": seal_file, "trust": trust, "ids": ids}  # fmt: skip


def _cli(w: dict[str, Any], row: str, out: Path, case: int = 1) -> tuple[int, dict[str, Any]]:
    code = R2.main([
        "grade", "--plan", str(w["plan"]), "--records", str(w["records"]), "--seal",
        str(w["seal"]), "--case", str(case), "--row", row, "--outputs", str(w["outputs"]),
        "--samples-root", str(w["samples"]), "--repo", str(w["repo"]), "--out", str(out),
    ])  # fmt: skip
    return code, json.loads(out.read_text())


@pytest.fixture
def small_case(monkeypatch: pytest.MonkeyPatch) -> None:
    for case, spec in SMALL.items():
        monkeypatch.setitem(IO.CASES, case, spec)


def _no_arrays(obj: Any) -> None:
    if isinstance(obj, dict):
        for v in obj.values():
            _no_arrays(v)
    elif isinstance(obj, list):
        assert len(obj) < 100, "report must not carry raw curves"
        for v in obj:
            _no_arrays(v)


def test_cli_td_em_end_to_end(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path)
    code, rep = _cli(w, "TD-EM", tmp_path / "td.json")
    assert rep["status"] == "PASS" and code == 0, rep
    names = {c["criterion"]: c for c in rep["criteria"]}
    assert set(names) == {"(a)", "(b)", "(c)", "(d)", "(e)"}
    assert names["(b)"]["verdict"] == "REPORT-ONLY"  # W_T < 1.0 mm: report-only
    assert all(v["valid"] for v in rep["run_validity"]) and len(rep["topas_runs"]) == 4
    assert {c["k"] for c in rep["ionmc_chunks"]} == {1, 2}
    assert all(f.startswith("work/") or f in ("request.json", "stdout.txt", "stderr.txt")
               for r in rep["topas_runs"] for f in r["files_sha256"])  # fmt: skip
    _no_arrays(rep)


def test_cli_td_em_detects_a_range_shift(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path, cases={1: {"shift": 0.3}})
    code, rep = _cli(w, "TD-EM", tmp_path / "td.json")
    assert code == 1 and rep["status"] == "FAIL"
    assert next(c for c in rep["criteria"] if c["criterion"] == "(a)")["verdict"] == "FAIL"


def test_cli_tl_em_end_to_end(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path)
    code, rep = _cli(w, "TL-EM", tmp_path / "tl.json")
    assert rep["status"] == "PASS" and code == 0, rep
    assert len(rep["points"]) == 4 and all(p["gated"] for p in rep["points"])
    for p in rep["points"]:
        assert p["a"]["r"] == pytest.approx(1.0, abs=0.01)
        assert p["sigma_topas"] == pytest.approx(SIG[p["slab"]], rel=2e-3)
    assert rep["criterion_b"] is None and "upper-end" in rep["criterion_b_note"]


def test_cli_chronology_failure_is_a_fail(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path, started_ok=False)  # runs started before the seal review finished
    code, rep = _cli(w, "TD-EM", tmp_path / "td.json")
    assert code == 1 and rep["status"] == "FAIL" and "blinding control" in rep["reason"]
    assert rep["category"] == "blinding control" and rep["rerun_allowed"] is False
    assert "criteria" not in rep


def test_cli_fails_closed_on_tampering(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path)
    npz = next(w["samples"].rglob("*.npz"))
    npz.write_bytes(npz.read_bytes() + b"x")
    code, rep = _cli(w, "TD-EM", tmp_path / "a.json")
    assert code == 2 and rep["status"] == "ERROR" and "sealed record" in rep["error"]
    w2 = _world(tmp_path / "b")
    csv = next(w2["outputs"].rglob("idd_dose.csv"))
    csv.write_text(csv.read_text() + "\n")
    code, rep = _cli(w2, "TD-EM", tmp_path / "b.json")
    assert code == 2 and "sha256 mismatch" in rep["error"]


@pytest.mark.parametrize("field", [*R2.io.RESULT_FIELDS, "task_id", "timed_out"])
def test_missing_result_field_fails_closed(tmp_path: Path, small_case: None, field: str) -> None:
    w = _world(tmp_path)
    rec = w["records"] / "runs" / w["ids"][(1, "idd-em", 0)] / "record.json"
    data = json.loads(rec.read_text())
    data.pop("attestation")
    del data[field]
    rec.write_text(json.dumps(_sign(data)))  # validly signed record without the field
    _relist(rec)
    code, rep = _cli(w, "TD-EM", tmp_path / "r.json")
    assert code == 2 and rep["status"] == "ERROR"


def test_unbound_or_failed_runs_fail_closed(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path)
    rec = w["records"] / "runs" / w["ids"][(1, "idd-em", 1)] / "record.json"
    good = rec.read_text()
    for bad in ({"succeeded": False}, {"configuration_id": "0" * 64}, {"run_id": "REF-other"},
                {"code_sha": "b" * 40}, {"seconds": -1}, {"started_at": "yesterday"},
                {"started_at": "2026-11-09T12:00:00"}, {"request_sha256": "0" * 64},
                {"stdout_sha256": "0" * 64}, {"stderr_sha256": "0" * 64}, {"timed_out": True},
                {"started_at": "2026-11-09T14:00:00Z"}, {"exit_code": 1}):  # fmt: skip
        rec.write_text(good)
        _resign(rec, bad)  # validly signed, but not bound / not successful
        code, rep = _cli(w, "TD-EM", tmp_path / "r.json")
        assert code == 2 and rep["status"] == "ERROR", bad
    rec.write_text(good)
    _relist(rec)
    assert _cli(w, "TD-EM", tmp_path / "r.json")[0] == 0
    rec.write_text(good.replace('"task_id": null', '"task_id": "V3-008"'))  # unsigned change
    _relist(rec)
    assert _cli(w, "TD-EM", tmp_path / "r.json")[0] == 2
    rec.write_text(good)
    lst = rec.parent / "artifacts.json"
    doc = json.loads(lst.read_text())
    doc["files"] = [e for e in doc["files"] if e["path"] != "record.json"]
    lst.write_text(json.dumps(doc))
    assert _cli(w, "TD-EM", tmp_path / "r.json")[0] == 2  # record.json must be listed
    _relist(rec)
    doc["files"].append({"path": "record.json", "sha256": "0" * 64})
    lst.write_text(json.dumps(doc))
    assert _cli(w, "TD-EM", tmp_path / "r.json")[0] == 2  # listed with another digest
    (w["records"] / "review.json").write_text(json.dumps({"review_id": "REVIEW-x"}))
    assert _cli(w, "TD-EM", tmp_path / "r.json")[0] == 2


def test_d_ext_rules() -> None:
    def pt(f: float, r: float, gated: bool = True) -> dict[str, Any]:
        return {"f": f, "gated": gated, "sigma_ionmc": r, "sigma_topas": 1.0,
                "_batches": [r * (1 + 1e-3 * ((b % 3) - 1)) for b in range(40)],
                "_parts": [1.0 + 1e-3 * (j - 1.5) for j in range(4)]}  # fmt: skip

    low = [pt(0.75, 1.10), pt(0.9, 1.12)]
    assert R2.d_ext(low, [pt(0.75, 1.05), pt(0.9, 1.07)])["verdict"] == "PASS"
    assert R2.d_ext(low, [pt(0.75, 1.30), pt(0.9, 1.32)])["verdict"] == "FAIL"
    none = [pt(0.75, 1.1, False), pt(0.9, 1.1, False)]
    assert R2.d_ext(none, low)["verdict"] == "NOT EVALUABLE"
    res = R2.tl_statistics([{"f": 0.25, "gated": False, "a": {"verdict": "REPORT-ONLY"}}], None)
    assert res["status"] == "NOT EVALUABLE"


def test_case_table_matches_steps_v8() -> None:
    sys.path.insert(0, str(_PATH.parents[1] / "transport"))
    try:
        steps_v8 = pytest.importorskip("steps_v8")
    finally:
        sys.path.pop(0)
    for case, (ion, eu, length, hz, hx, cset) in steps_v8.CASES.items():
        assert IO.CASES[case][:2] == (ion, eu) and IO.CASES[case][3:] == (length, hz, hx, cset)


def test_cli_td_full_uses_a_frag(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path, cases={1: {"full_shift": 0.062}}, configs=("idd-em", "idd-full"))
    code, rep = _cli(w, "TD-FULL", tmp_path / "f.json")
    a = next(c for c in rep["criteria"] if c["criterion"] == "(a)")
    r80 = rep["pooled"]["r80"]["topas"]
    assert a["tol"] == pytest.approx(0.05 + 0.0015 * r80)  # TD-EM floor + A_frag (He)
    assert a["verdict"] == "PASS" and abs(a["x"] - a["r"]) > 0.05  # passes only with A_frag
    assert rep["status"] == "PASS" and code == 0, rep
    gating = {c["criterion"] for c in rep["criteria"] if c["verdict"] != "REPORT-ONLY"}
    assert gating == {"(a)", "(c)"}
    ro = rep["report_only"]
    assert ro["delta_frag_mm"]["value"] == pytest.approx(0.054, abs=0.002)
    assert ro["r80_topas_iddprimary_mm"] == pytest.approx(r80, abs=1e-6)
    assert "idd_primary_dose.csv" in json.dumps(rep["topas_runs"])
    assert rep["not_computed"]
    w2 = _world(tmp_path / "b", cases={1: {"full_shift": 0.08}}, configs=("idd-em", "idd-full"))
    code, rep = _cli(w2, "TD-FULL", tmp_path / "f2.json")
    assert code == 1 and rep["status"] == "FAIL"


def test_cli_tl_full_end_to_end(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path, configs=("idd-em", "lat-em", "lat-full"))
    code, rep = _cli(w, "TL-FULL", tmp_path / "t.json")
    assert rep["status"] == "PASS" and code == 0, rep
    for p in rep["points"]:  # graded on the generation-0 profiles, gated by TOPAS-EM
        assert p["sigma_topas"] == pytest.approx(SIG[p["slab"]], rel=2e-3)
        assert p["gate_sigma_topas_em"] == pytest.approx(SIG[p["slab"]], rel=2e-3)
    ro = rep["report_only"]
    for p, s_all, ratio in zip(rep["points"], ro["sigma_topas_full_all_particle"],
                               ro["sigma_full_primary_over_em"], strict=True):  # fmt: skip
        assert s_all == pytest.approx(1.3 * SIG[p["slab"]], rel=0.02)
        assert ratio == pytest.approx(1.0, abs=2e-3)
    assert "lat_x_primary.bin" in json.dumps(rep["topas_runs"])


@pytest.mark.parametrize(("scale", "verdict"), [(1.05, "PASS"), (1.2, "FAIL")])
def test_cli_tl_b_through_grade(
    tmp_path: Path, small_case: None, scale: float, verdict: str
) -> None:
    w = _world(tmp_path, cases={1: {}, 3: {"ion_scale": scale}})
    code, rep = _cli(w, "TL-EM", tmp_path / "b.json", case=3)
    b = rep["criterion_b"]
    assert (b["low_case"], b["upper_case"]) == (1, 3) and b["depths"] == [0.75, 0.9]
    assert b["x"] == pytest.approx(scale - 1.0, abs=0.01)
    assert b["verdict"] == verdict and rep["criterion_a"] == "PASS"
    assert rep["status"] == verdict and code == (0 if verdict == "PASS" else 1)
    assert {c["case"] for c in rep["ionmc_chunks"]} == {1, 3}


def test_cli_tl_b_not_evaluable(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path, cases={1: {"sig": [0.2] * int(L)}, 3: {}})  # low end: nothing gated
    code, rep = _cli(w, "TL-EM", tmp_path / "n.json", case=3)
    assert rep["criterion_b"]["verdict"] == "NOT EVALUABLE"
    assert rep["status"] == "NOT EVALUABLE" and code == 1


@pytest.mark.parametrize(
    "tamper",
    [
        {"stderr": "G4Exception: something\n"},
        {"stdout": lambda s: s + "Track stuck or not moving: killed\n"},
        {"stdout": lambda s: s.replace("histories: 1000", "histories: 999")},
        {"stdout": lambda s: s.replace("histories: 1000\n", "")},
        {"stdout": lambda s: s.replace("UrbanMsc", "WentzelVIUni", 2)},
        {"input_extra": "d:Ph/Default/CutForAllParticles = 0.1 mm\n"},
    ],
)
def test_cli_run_validity_fail(tmp_path: Path, small_case: None, tamper: dict[str, Any]) -> None:
    w = _world(tmp_path, tamper=tamper)
    code, rep = _cli(w, "TD-EM", tmp_path / "v.json")
    assert code == 1 and rep["status"] == "FAIL", rep
    assert rep["category"] == "run validity" and rep["rerun_allowed"] is True
    assert "criteria" not in rep
    probs = [p for x in rep["part_validity"] + rep["replicate_gate"] for p in x["problems"]]
    assert probs


def test_begin_processing_lines_are_not_killed_tracks(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path, tamper={"stdout": lambda s: s + "Begin processing (killed=0)\n"})
    assert _cli(w, "TD-EM", tmp_path / "v.json")[1]["status"] == "PASS"


@pytest.mark.parametrize(
    ("tamper", "needle"),
    [
        ({"unsigned": True}, "status unverified"),
        ({"signer": "host-runner"}, "status invalid"),
        ({"started_after_signing": "2026-11-09T12:30:00+00:00"}, "status invalid"),
    ],
)
def test_result_attestation_fails_closed(
    tmp_path: Path, small_case: None, tamper: dict[str, Any], needle: str
) -> None:
    w = _world(tmp_path, tamper=tamper)
    code, rep = _cli(w, "TD-EM", tmp_path / "s.json")
    assert code == 2 and rep["status"] == "ERROR" and needle in rep["error"], rep
    assert "V3-INFRA-004" in rep["error"]


def test_untrusted_key_fails_closed(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path)
    good = w["trust"].read_text()
    w["trust"].write_text(json.dumps({"reference-runner": []}))  # uncommitted change
    code, rep = _cli(w, "TD-EM", tmp_path / "s.json")
    assert code == 2 and "differs from HEAD" in rep["error"]
    _commit(w["repo"], "drop the key")  # committed trust without the key
    code, rep = _cli(w, "TD-EM", tmp_path / "s.json")
    assert code == 2 and "status unverified" in rep["error"]
    w["trust"].write_text(good)
    _commit(w["repo"], "restore the key")
    code, rep = _cli(w, "TD-EM", tmp_path / "s.json")
    assert code == 0 and rep["trust_file"]["head_sha"] == _git(w["repo"], "rev-parse", "HEAD")


def test_cli_has_no_trust_option(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path)
    with pytest.raises(SystemExit):
        R2.main(["grade", "--plan", str(w["plan"]), "--records", str(w["records"]), "--seal",
                 str(w["seal"]), "--case", "1", "--row", "TD-EM", "--trust", "x.json"])  # fmt: skip


def test_coordinated_output_tampering_fails_closed(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path)
    rid = w["ids"][(1, "idd-em", 0)]
    csv = w["outputs"] / rid / "work" / "idd_dose.csv"
    lines = csv.read_text().splitlines()
    lines[-1] = lines[-1].rsplit(",", 2)[0] + ", 1.0, 0"  # change one value
    csv.write_text("\n".join(lines) + "\n")
    lst = w["records"] / "runs" / rid / "artifacts.json"
    doc = json.loads(lst.read_text())
    for e in doc["files"]:  # the unsigned listing is updated consistently
        if e["path"] == "work/idd_dose.csv":
            e["sha256"] = hashlib.sha256(csv.read_bytes()).hexdigest()
            e["bytes"] = csv.stat().st_size
    lst.write_text(json.dumps(doc))
    code, rep = _cli(w, "TD-EM", tmp_path / "t.json")
    assert code == 2 and "output_files_sha256" in rep["error"], rep


def test_signed_start_before_review_is_blinding_fail(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path)
    assert _cli(w, "TD-EM", tmp_path / "a.json")[1]["status"] == "PASS"
    rec = w["records"] / "runs" / w["ids"][(1, "idd-em", 2)] / "record.json"
    _resign(rec, {"started_at": "2026-11-09T08:59:59Z"})  # finished_at stays late: unused
    code, rep = _cli(w, "TD-EM", tmp_path / "b.json")
    assert code == 1 and rep["category"] == "blinding control"
    bad = [v["run_id"] for v in rep["run_validity"] if not v["valid"]]
    assert bad == [w["ids"][(1, "idd-em", 2)]]


def _drop_block(text: str, start: str) -> str:
    """Remove the hadronic-summary block that starts with ``start`` (up to the next separator)."""
    i = text.index(start)
    j = text.index(SEP, i)
    return text[:i] + text[j + len(SEP) :]


@pytest.mark.parametrize(
    "edit",
    [
        lambda s: s.replace(
            "        Model:                  QMDModel: 100 MeV/n ---> 6 GeV/n\n", ""
        ),
        lambda s: s.replace("Radioactive Decay Physics Parameters", "Decay Parameters"),
        lambda s: s.replace(" hInelastic QGSP_BIC_HP  Thresholds: \n", ""),
        lambda s: _drop_block(s, "                           Hadronic Processes for mu-"),
        lambda s: s.replace("NeutronHPElastic", "hElasticCHIPS"),
        lambda s: s.replace("Use ICRU90 data                                    1", "x"),
    ],
)
def test_full_module_markers_missing_fail(tmp_path: Path, small_case: None, edit: Any) -> None:
    w = _world(tmp_path, configs=("idd-full",), tamper={"stdout": edit})
    code, rep = _cli(w, "TD-FULL", tmp_path / "m.json")
    assert code == 1 and rep["category"] == "run validity", rep
    assert any("module" in p for x in rep["part_validity"] for p in x["problems"])


def test_full_logs_intact_pass_and_list_unverifiable(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path, configs=("idd-full",))
    code, rep = _cli(w, "TD-FULL", tmp_path / "m.json")
    assert code == 0 and rep["status"] == "PASS", rep
    assert rep["module_markers_unverifiable"] == []


@pytest.mark.parametrize(
    "inject",
    [" hInelastic QGSP_BIC_HP  Thresholds: \n", "        Model:  QMDModel: x\n",
     "======          Radioactive Decay Physics Parameters           =======\n"],
)  # fmt: skip
def test_em_log_with_hadronic_marker_fails(tmp_path: Path, small_case: None, inject: str) -> None:
    w = _world(tmp_path, tamper={"stdout": lambda s: inject + s})
    code, rep = _cli(w, "TD-EM", tmp_path / "e.json")
    assert code == 1 and rep["category"] == "run validity"
    assert any("hadronic marker" in p for x in rep["part_validity"] for p in x["problems"])


def test_decay_marker_and_unknown_modules(tmp_path: Path, small_case: None) -> None:
    w = _world(tmp_path, configs=("idd-full",),
               tamper={"stdout": lambda s: s.replace("Decay,     ", "", 1)})  # fmt: skip
    code, rep = _cli(w, "TD-FULL", tmp_path / "a.json")
    assert code == 1 and rep["category"] == "run validity"
    assert any("g4decay" in p for x in rep["part_validity"] for p in x["problems"])
    w2 = _world(tmp_path / "b", tamper={"stdout": lambda s: s.replace("eIoni\n", "eIoni, Decay\n")})
    code, rep = _cli(w2, "TD-EM", tmp_path / "b.json")
    assert code == 1 and any("Decay" in p for x in rep["part_validity"] for p in x["problems"])
    mods = 'sv:Ph/Default/Modules = 2 "g4em-standard_opt4" "g4foo"\n'
    w3 = _world(tmp_path / "c", tamper={"modules": mods})
    code, rep = _cli(w3, "TD-EM", tmp_path / "c.json")
    assert code == 1 and rep["category"] == "run validity"
    assert rep["module_markers_unverifiable"] == ["g4foo"]
