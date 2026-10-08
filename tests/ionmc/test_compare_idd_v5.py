"""compare_idd_v5.py (row V5): synthetic curves with known R80 and plateau, synthetic runs."""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pytest

from tests.ionmc.test_compare_batches import _write

SCRIPT = Path(__file__).resolve().parents[2] / "validation/scripts/reference/compare_idd_v5.py"
NZ = 348
PLATEAU, PEAK, R80 = 2.0, 6.0, 122.25


def _load():
    spec = importlib.util.spec_from_file_location("compare_idd_v5", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load()


def curve(shift_mm: float = 0.0, plateau_scale: float = 1.0, extra: float = 0.0) -> np.ndarray:
    """Plateau 2.0 to 100 mm, linear rise to 6.0 at 120.25 mm, linear fall to 0 at 130.25 mm."""
    z = M.depth_centres(NZ) - shift_mm
    up = np.interp(z, [100.0, 120.25], [PLATEAU, PEAK])
    down = np.clip(PEAK * (130.25 - z) / 10.0, 0.0, None)
    y = np.where(z <= 120.25, np.where(z < 100.0, PLATEAU, up), down) * plateau_scale
    return y + extra * (M.depth_centres(NZ) < 100.0)


def noisy(base: np.ndarray, n: int, rel: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return base * (1.0 + rel * rng.standard_normal((n, base.size)))


def test_synthetic_curve_metrics_are_known() -> None:
    m = M.curve_metrics(curve())
    assert m["plateau"] == pytest.approx(PLATEAU)
    assert m["r80_mm"] == pytest.approx(R80, abs=1e-9)
    assert m["peak_over_plateau"] == pytest.approx(3.0)
    assert m["total_deposit"] == pytest.approx(curve().sum() * 0.05)
    n = int(round((60.0 - 20.0) / 0.5))
    assert m["plateau_integral_mev"] == pytest.approx(PLATEAU * n * 0.05)


def test_t_quantile_and_tost() -> None:
    assert M.t95(2) == pytest.approx(2.92)
    assert M.t95(19) == pytest.approx(1.7291)
    assert M.t95(500) == pytest.approx(1.6973)  # conservative beyond the table
    ok = M.tost(M.Est(0.004, 0.003, 19), 0.02)  # CI +-0.0052
    assert ok["pass"] and ok["ci90"][1] < 0.02
    bad = M.tost(M.Est(0.015, 0.004, 19), 0.02)  # upper 0.0219 > tol although |diff| < tol
    assert not bad["pass"]
    inconclusive = M.tost(M.Est(0.0, 0.05, 2), 0.02)
    assert not inconclusive["pass"]


def test_jackknife_se_matches_analytic_mean() -> None:
    b = noisy(curve(), 20, 0.05, 1)
    est = M.jackknife_metrics(b)["total_deposit"]
    per = b.sum(axis=1) * 0.05
    assert est.value == pytest.approx(per.mean())
    assert est.se == pytest.approx(per.std(ddof=1) / math.sqrt(20), rel=1e-9)
    assert est.df == 19


def _evaluate(ion_scale=1.0, shift=0.0, delta_ion=0.2, delta_ref=0.2, ion_peak_scale=1.0):
    def ion_curve(nuc: bool, e: int) -> np.ndarray:
        c = curve(shift_mm=shift, plateau_scale=ion_scale, extra=delta_ion if nuc else 0.0)
        return noisy(c, 20, 0.03, 10 + e + nuc)

    def eng_curve(nuc: bool, e: int, s: int) -> np.ndarray:
        return noisy(curve(extra=delta_ref if nuc else 0.0), 3, 0.0067, 100 + e + nuc + s)

    ion = {(e, n): ion_curve(n, e) for e in (150, 200) for n in (True, False)}
    eng = {
        "topas": {(e, n): eng_curve(n, e, 0) for e in (150, 200) for n in (True, False)},
        "mcsquare": {(e, n): eng_curve(n, e, 7) for e in (150, 200) for n in (True, False)},
    }
    return M.evaluate(ion, eng)


def test_identical_curves_pass_all_criteria() -> None:
    doc = _evaluate()
    assert doc["pass"] and doc["kappa_rule"]["decision"] == "kappa-not-used"
    assert doc["engines"]["topas"]["role"] == "gating"
    assert doc["engines"]["mcsquare"]["role"] == "report-only"
    for e in ("150", "200"):
        rows = doc["engines"]["topas"]["energies"][e]
        assert rows["r80_mm"]["ionmc"]["value"] == pytest.approx(R80, abs=0.3)
        assert rows["plateau"]["reference"]["value"] == pytest.approx(PLATEAU + 0.2, rel=0.01)
        assert set(M.TOLERANCES) <= set(rows)


def test_plateau_shift_fails_plateau_only() -> None:
    doc = _evaluate(ion_scale=1.03)
    rows = doc["engines"]["topas"]["energies"]["150"]
    assert not rows["plateau"]["pass"] and not doc["pass"]
    assert rows["r80_mm"]["pass"]  # R80 does not depend on the amplitude


def test_r80_shift_fails() -> None:
    rows = _evaluate(shift=1.0)["engines"]["topas"]["energies"]["150"]
    assert not rows["r80_mm"]["pass"]
    assert rows["r80_mm"]["diff"] == pytest.approx(1.0, abs=0.2)
    assert _evaluate(shift=0.25)["engines"]["topas"]["energies"]["150"]["r80_mm"]["pass"]


def test_delta_idd_criterion() -> None:
    assert _evaluate(delta_ion=0.2, delta_ref=0.2)["engines"]["topas"]["energies"]["150"][
        "delta_idd_plateau_integral"
    ]["pass"]
    bad = _evaluate(delta_ion=0.3, delta_ref=0.2)["engines"]["topas"]["energies"]["150"][
        "delta_idd_plateau_integral"
    ]
    assert not bad["pass"] and bad["diff"] == pytest.approx(0.5, abs=0.05)
    assert bad["reference_delta_significant"]


def test_total_deposit_and_kappa_rule() -> None:
    doc = _evaluate(ion_scale=1.012)  # all bins +1.2 %: total fails with D > 0, |D| <= 1.5 %
    t = doc["engines"]["topas"]["energies"]["150"]["total_deposit"]
    assert not t["pass"] and t["diff"] == pytest.approx(0.012 * 0.9, abs=0.004) or t["diff"] > 0.008
    assert doc["kappa_rule"]["decision"] == "build-kappa"
    assert _evaluate(ion_scale=0.985)["kappa_rule"]["decision"] == "stop-kappa-cannot-be-the-remedy"
    assert _evaluate(ion_scale=1.03)["kappa_rule"]["decision"] == "stop-kappa-cannot-be-the-remedy"


def test_kappa_rule_cases() -> None:
    p = lambda d, ok: {"diff": d, "pass": ok}  # noqa: E731
    assert M.kappa_rule({150: p(0.001, True), 200: p(-0.002, True)})["decision"] == "kappa-not-used"
    assert M.kappa_rule({150: p(0.012, False), 200: p(0.001, True)})["decision"] == "build-kappa"
    assert (
        M.kappa_rule({150: p(0.016, False), 200: p(0.0, True)})["decision"]
        == "stop-kappa-cannot-be-the-remedy"
    )
    mixed = M.kappa_rule({150: p(0.012, False), 200: p(-0.012, False)})
    assert mixed["decision"] == "stop-kappa-cannot-be-the-remedy" and mixed["blocked_energies"] == [
        200
    ]


# -- synthetic runs and the CLI ---------------------------------------------------------------
def _scaled(e: int, nuc: bool) -> np.ndarray:
    off = curve()
    return curve(extra=0.2 if nuc else 0.0) * (0.90 * e / (off.sum() * 0.05))  # common scale


def _topas_input(e: int, nuc: bool, seed: int, extra: str = "") -> bytes:
    mods = (
        'sv:Ph/Default/Modules = 2 "g4em-standard_opt4" "g4h-phy_QGSP_BIC_HP"'
        if nuc
        else 'sv:Ph/Default/Modules = 1 "g4em-standard_opt4"'
    )
    return (
        f"i:Ts/Seed = {seed}\nd:So/Beam/BeamEnergy = {e} MeV\n{mods}\n"
        f"i:So/Beam/NumberOfHistoriesInRun = 100000\n{extra}"
    ).encode()


def _mc_input(nuc: bool, seed: int, extra: str = "") -> bytes:
    return (
        f"RNG_Seed {seed}\nNum_Primaries 100000\n"
        f"Simulate_Nuclear_Interactions {'True' if nuc else 'False'}\n{extra}"
    ).encode()


def _topas_run(
    tmp: Path, e: int, nuc: bool, seed: int, noise: float, label: tuple[int, bool] | None = None,
    extra: str = "",
) -> Path:  # fmt: skip
    idd = _scaled(e, nuc) * (1.0 + noise)
    gy = idd * 1e5 * 0.05 * M.MEV_J / 0.08  # 400 x 400 x 0.5 mm = 80 cm3 = 0.08 kg
    rows = "".join(f"0,0,{i},{float(v)!r},0.0\n" for i, v in enumerate(gy))
    csv = (
        f"# X in 1 bin of 400 mm\n# Y in 1 bin of 400 mm\n# Z in {NZ} bins of 0.5 mm\n"
        f"# DoseToMedium ( Gy ) : Sum Standard_Deviation\n{rows}"
    )
    case = {
        "histories": 100000,
        "seeds": [seed],
        "input": "input.txt",
        "v5": {
            "row": "V5",
            "energy_mev": (label or (e, nuc))[0],
            "nuclear": (label or (e, nuc))[1],
        },
    }
    files = {
        "work/idd_dose.csv": csv.encode(),
        "inputs/input.txt": _topas_input(e, nuc, seed, extra),
    }
    out = "Particle source Beam: Total number of histories: 100000\n"
    return _write(tmp, "topas", files, case, f"T-{seed}", stdout=out)


def _mc_run(
    tmp: Path, e: int, nuc: bool, seed: int, noise: float, label: tuple[int, bool] | None = None,
    extra: str = "",
) -> Path:  # fmt: skip
    idd = _scaled(e, nuc) * (1.0 + noise)
    vol = np.zeros((5, NZ, 5), dtype="<f4")  # (z, y, x); beam in the centre voxel
    vol[2, ::-1, 2] = idd / 64.0 * 1e6  # MeV/g/primary x 1e6 -> eV/g; 8 x 8 cm2 voxel face
    mhd = (
        f"NDims = 3\nDimSize = 5 {NZ} 5\nElementSpacing = 80 0.5 80\nOffset = 0 0 0\n"
        "ElementType = MET_FLOAT\nElementByteOrderMSB = False\nElementDataFile = Dose.raw\n"
    )
    case = {
        "histories": 100000,
        "seeds": [seed],
        "input": "config.txt",
        "v5": {
            "row": "V5",
            "energy_mev": (label or (e, nuc))[0],
            "nuclear": (label or (e, nuc))[1],
        },
    }
    files = {
        "work/Outputs/Dose.mhd": mhd.encode(),
        "work/Outputs/Dose.raw": vol.tobytes(),
        "inputs/config.txt": _mc_input(nuc, seed, extra),
        "inputs/Plan.txt": f"####Energy (MeV)\n{e}\n####NbOfScannedSpots\n".encode(),
    }
    return _write(
        tmp, "mcsquare", files, case, f"M-{seed}", stdout="Nbr primaries simulated: 100000\n"
    )


def test_engine_loaders_reproduce_the_idd(tmp_path: Path) -> None:
    runs_t = [_topas_run(tmp_path / f"t{s}", 150, True, 1 + s, 0.0) for s in range(3)]
    runs_m = [_mc_run(tmp_path / f"m{s}", 150, True, 1 + s, 0.0) for s in range(3)]
    for engine, runs in (("topas", runs_t), ("mcsquare", runs_m)):
        g = M.load_engine_groups(runs, engine, {150: NZ})[(150, True)]
        assert g["curves"].shape == (3, NZ) and g["seeds"] == [1, 2, 3]
        np.testing.assert_allclose(g["curves"][0], _scaled(150, True), rtol=2e-6)


def test_loader_failures(tmp_path: Path) -> None:
    runs = [_topas_run(tmp_path / f"t{s}", 150, True, 1 + s, 0.0) for s in range(2)]
    with pytest.raises(M.IddError, match="need >= 3"):
        M.load_engine_groups(runs, "topas", {150: NZ})
    dup = [_topas_run(tmp_path / f"d{s}", 150, True, 7, 0.0) for s in range(3)]
    with pytest.raises(M.IddError, match="duplicate seed"):
        M.load_engine_groups(dup, "topas", {150: NZ})
    wrong_unit = _topas_run(tmp_path / "u", 150, True, 9, 5.0)  # 6x the beam energy deposited
    with pytest.raises(M.IddError, match="unit or geometry"):
        M.load_engine_groups([wrong_unit], "topas", {150: NZ})
    with pytest.raises(M.IddError, match="engine"):
        M.load_engine_groups([_mc_run(tmp_path / "x", 150, True, 3, 0.0)], "topas", {150: NZ})


def _partial(path: Path, e: int, nuc: bool, batches: np.ndarray) -> None:
    doc = {
        "row": f"v5-{e}-{'on' if nuc else 'off'}",
        "energy": float(e),
        "nuclear": nuc,
        "seed": 1,
        "n": 100000,
        "n_batches": 20,
        "bin_mm": 0.5,
        "unit": M.UNIT,
        "density_g_cm3": 1.0,
        "idd_batches": batches.tolist(),
        "valid": True,
        "reduced": False,
    }
    doc["content_sha256"] = M.content_digest(doc)
    path.write_text(json.dumps(doc, sort_keys=True))


def test_cli_end_to_end_and_partial_checks(tmp_path: Path) -> None:
    ion = tmp_path / "ion"
    ion.mkdir()
    for e in (150, 200):
        for nuc in (True, False):
            _partial(
                ion / f"v5-{e}-{'on' if nuc else 'off'}.json",
                e,
                nuc,
                noisy(_scaled(e, nuc), 20, 0.03, e + nuc),
            )
    topas = [
        _topas_run(tmp_path / f"t{e}{n}{s}", e, n, 100 * e + 10 * n + s + 1, 0.001 * (s - 1))
        for e in (150, 200)
        for n in (True, False)
        for s in range(3)
    ]
    mc = [
        _mc_run(tmp_path / f"m{e}{n}{s}", e, n, 100 * e + 10 * n + s + 1, 0.001 * (s - 1))
        for e in (150, 200)
        for n in (True, False)
        for s in range(3)
    ]
    out = tmp_path / "out" / "v5.json"
    argv = [
        "--ionmc-dir",
        str(ion),
        "--topas-runs",
        *map(str, topas),
        "--mcsquare-runs",
        *map(str, mc),
        "--output",
        str(out),
    ]
    assert M.main(argv) == 0
    doc = json.loads(out.read_text())
    assert (
        doc["row"] == "V5"
        and doc["pass"] is True
        and doc["kappa_rule"]["decision"] == "kappa-not-used"
    )
    assert set(doc["engines"]) == {"topas", "mcsquare"}
    assert doc["engine_inputs"]["topas"]["150-on"]["seeds"] == [15011, 15012, 15013]
    # tampered partial: content hash mismatch fails closed
    p = ion / "v5-150-on.json"
    d = json.loads(p.read_text())
    d["n"] = 99999
    p.write_text(json.dumps(d))
    with pytest.raises(M.IddError, match="content_sha256"):
        M.main(argv)


def test_lineage_mislabeled_case_json_is_refused(tmp_path: Path) -> None:
    """C19 F2: energy and configuration come from the verified native input; a case.json label
    that disagrees (150 MeV EM-only labelled as 200 MeV full) is refused, not trusted."""
    for engine, maker in (("topas", _topas_run), ("mcsquare", _mc_run)):
        bad = maker(tmp_path / f"{engine}-e", 150, False, 5, 0.0, label=(200, False))
        with pytest.raises(M.IddError, match="mislabeled"):
            M.load_engine_groups([bad], engine, {150: NZ, 200: NZ})
        bad = maker(tmp_path / f"{engine}-n", 150, False, 5, 0.0, label=(150, True))
        with pytest.raises(M.IddError, match="mislabeled"):
            M.load_engine_groups([bad], engine, {150: NZ, 200: NZ})


def test_lineage_derives_the_group_from_the_native_input(tmp_path: Path) -> None:
    runs = [_topas_run(tmp_path / f"t{s}", 200, False, 1 + s, 0.0) for s in range(3)]
    g = M.load_engine_groups(runs, "topas", {150: NZ, 200: NZ})
    assert set(g) == {(200, False)}
    fps = g[(200, False)]["fingerprints"]
    assert len(fps) == 3 and {f["group"] for f in fps.values()} == {"topas-emonly"}
    assert len({f["config_sha256"] for f in fps.values()}) == 1


def test_lineage_failed_fingerprint_and_differing_replicate_are_refused(tmp_path: Path) -> None:
    ok = [_topas_run(tmp_path / f"t{s}", 150, True, 1 + s, 0.0) for s in range(2)]
    odd = _topas_run(
        tmp_path / "odd", 150, True, 3, 0.0, extra="d:Ph/Default/CutForProton = 1 mm\n"
    )
    with pytest.raises(M.IddError, match="replicates differ in config_sha256"):
        M.load_engine_groups([*ok, odd], "topas", {150: NZ})
    # a run whose engine summary contradicts the declared histories fails the fingerprint
    bad = _mc_run(tmp_path / "m", 150, True, 3, 0.0)
    (bad / "stdout.txt").write_text("Nbr primaries simulated: 5\n")
    mf = json.loads((bad / "transfer-manifest.json").read_text())
    import hashlib

    blob = (bad / "stdout.txt").read_bytes()
    mf["files"]["stdout.txt"] = {"bytes": len(blob), "sha256": hashlib.sha256(blob).hexdigest()}
    (bad / "transfer-manifest.json").write_text(json.dumps(mf))
    with pytest.raises(M.IddError, match="summary histories"):
        M.load_engine_groups([bad], "mcsquare", {150: NZ})


def test_cli_exits_nonzero_when_the_gating_verdict_fails(tmp_path: Path) -> None:
    ion = tmp_path / "ion"
    ion.mkdir()
    for e in (150, 200):
        for nuc in (True, False):
            b = noisy(_scaled(e, nuc) * 0.9, 20, 0.03, e + nuc)  # plateau 10 % low: V5 fails
            _partial(ion / f"v5-{e}-{'on' if nuc else 'off'}.json", e, nuc, b)
    topas = [
        _topas_run(tmp_path / f"t{e}{n}{s}", e, n, 100 * e + 10 * n + s + 1, 0.001 * (s - 1))
        for e in (150, 200)
        for n in (True, False)
        for s in range(3)
    ]
    out = tmp_path / "out.json"
    argv = ["--ionmc-dir", str(ion), "--topas-runs", *map(str, topas), "--output", str(out)]
    assert M.main(argv) == 1
    doc = json.loads(out.read_text())
    assert doc["pass"] is False
    fp = doc["engine_inputs"]["topas"]["150-on"]["fingerprints"]
    assert len(fp) == 3 and all(len(v["config_sha256"]) == 64 for v in fp.values())
