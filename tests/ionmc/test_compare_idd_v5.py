"""compare_idd_v5.py (row V5): synthetic curves with known R80 and plateau, synthetic runs."""

# ruff: noqa: E501

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pytest

from ionmc.reference.runs import load_run, read_verified
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


FULL_MODULES = (
    'sv:Ph/Default/Modules = 6 "g4em-standard_opt4" "g4h-phy_QGSP_BIC_HP" "g4h-elastic_HP" '
    '"g4stopping" "g4ion-binarycascade" "g4decay"'
)
EMONLY_MODULES = 'sv:Ph/Default/Modules = 1 "g4em-standard_opt4"'


def _topas_input(e: int, nuc: bool, seed: int, extra: str = "") -> bytes:
    mods = FULL_MODULES if nuc else EMONLY_MODULES
    return (
        f"i:Ts/Seed = {seed}\nd:So/Beam/BeamEnergy = {e} MeV\nu:So/Beam/BeamEnergySpread = 0\n"
        f"{mods}\nd:Ph/Default/CutForAllParticles = 0.05 mm\n"
        f"i:So/Beam/NumberOfHistoriesInRun = 100000\n{extra}"
    ).encode()


def _mc_input(nuc: bool, seed: int, extra: str = "") -> bytes:
    flag = "True" if nuc else "False"
    return (
        f"RNG_Seed {seed}\nNum_Primaries 100000\nSimulate_Nuclear_Interactions {flag}\n"
        f"Simulate_Secondary_Protons {flag}\nSimulate_Secondary_Deuterons {flag}\n"
        f"Simulate_Secondary_Alphas {flag}\n{extra}"
    ).encode()


def _topas_run(
    tmp: Path, e: int, nuc: bool, seed: int, noise: float, label: tuple[int, bool] | None = None,
    extra: str = "", input_bytes: bytes | None = None,
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
            "mode": "full" if (label or (e, nuc))[1] else "emonly",
        },
    }
    files = {
        "work/idd_dose.csv": csv.encode(),
        "inputs/input.txt": input_bytes or _topas_input(e, nuc, seed, extra),
    }
    out = "Particle source Beam: Total number of histories: 100000\n"
    return _write(tmp, "topas", files, case, f"T-{seed}", stdout=out)


def _mc_run(
    tmp: Path, e: int, nuc: bool, seed: int, noise: float, label: tuple[int, bool] | None = None,
    extra: str = "", input_bytes: bytes | None = None,
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
            "mode": "nuclear-on" if (label or (e, nuc))[1] else "nuclear-off",
        },
    }
    files = {
        "work/Outputs/Dose.mhd": mhd.encode(),
        "work/Outputs/Dose.raw": vol.tobytes(),
        "inputs/config.txt": input_bytes or _mc_input(nuc, seed, extra),
        "inputs/Plan.txt": f"####Energy (MeV)\n{e}\n####NbOfScannedSpots\n".encode(),
    }
    return _write(
        tmp, "mcsquare", files, case, f"M-{seed}", stdout="Nbr primaries simulated: 100000\n"
    )


def _freeze(root: Path, runs: list[Path], strip_v5: bool = False) -> Path:
    """A committed-cases dir ``root/<engine>/<family>-seed<k>`` of byte-identical copies of the run inputs."""
    for rd in runs:
        run = load_run(rd)
        v5 = run.case.get("v5", {})
        fam = f"syn-{v5.get('energy_mev')}-{v5.get('mode')}"
        base = root / run.engine / f"{fam}-seed{run.case['seeds'][0]}"
        d, k = base, 0
        while d.exists():  # a second run of the same seed (negative tests)
            k += 1
            d = base.with_name(f"{base.name}x{k}")
        for rel in sorted(run.files):
            if rel.startswith("inputs/"):
                f = d / rel[len("inputs/") :]
                f.parent.mkdir(parents=True, exist_ok=True)
                f.write_bytes(read_verified(run, rel))
        if strip_v5:
            case = json.loads((d / "case.json").read_text())
            case.pop("v5")
            (d / "case.json").write_text(json.dumps(case))
    for eng in ("topas", "mcsquare"):
        (root / eng).mkdir(parents=True, exist_ok=True)
    return root


def _load_groups(tmp: Path, runs: list[Path], engine: str, **kw):
    cases = _freeze(tmp / "cases", kw.pop("frozen", runs))
    return M.load_engine_groups(runs, engine, {150: NZ, 200: NZ}, cases)


def test_engine_loaders_reproduce_the_idd(tmp_path: Path) -> None:
    runs_t = [_topas_run(tmp_path / f"t{s}", 150, True, 1 + s, 0.0) for s in range(3)]
    runs_m = [_mc_run(tmp_path / f"m{s}", 150, True, 1 + s, 0.0) for s in range(3)]
    cases = _freeze(tmp_path / "cases", [*runs_t, *runs_m])
    for engine, runs in (("topas", runs_t), ("mcsquare", runs_m)):
        g = M.load_engine_groups(runs, engine, {150: NZ}, cases)[(150, True)]
        assert g["curves"].shape == (3, NZ) and g["seeds"] == [1, 2, 3]
        np.testing.assert_allclose(g["curves"][0], _scaled(150, True), rtol=2e-6)
        assert g["family"] == "syn-150-" + ("full" if engine == "topas" else "nuclear-on")
        b = g["bound_cases"]["T-1" if engine == "topas" else "M-1"]
        assert (
            b["bound_case"].startswith(f"{engine}/syn-150-")
            and "case.json" in b["committed_sha256"]
        )


def test_loader_failures(tmp_path: Path) -> None:
    runs = [_topas_run(tmp_path / f"t{s}", 150, True, 1 + s, 0.0) for s in range(2)]
    with pytest.raises(M.IddError, match="need >= 3"):
        _load_groups(tmp_path / "a", runs, "topas")
    # same seed, otherwise distinguishable inputs (each bound to its own committed case)
    dup = [_topas_run(tmp_path / f"d{s}", 150, True, 7, 0.0, extra=f"# {s}\n") for s in range(3)]
    with pytest.raises(M.IddError, match="duplicate seed"):
        _load_groups(tmp_path / "b", dup, "topas")
    wrong_unit = _topas_run(tmp_path / "u", 150, True, 9, 5.0)  # 6x the beam energy deposited
    with pytest.raises(M.IddError, match="unit or geometry"):
        _load_groups(tmp_path / "c", [wrong_unit], "topas")
    with pytest.raises(M.IddError, match="engine"):
        _load_groups(tmp_path / "d", [_mc_run(tmp_path / "x", 150, True, 3, 0.0)], "topas")


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
    cases = _freeze(tmp_path / "cases", [*topas, *mc])
    argv = [
        "--ionmc-dir",
        str(ion),
        "--cases-dir",
        str(cases),
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
    for eng, fam in (("topas", "syn-150-full"), ("mcsquare", "syn-150-nuclear-on")):
        grp = doc["engine_inputs"][eng]["150-on"]
        assert grp["family"] == fam and len(grp["binding"]) == 3
        b = next(iter(grp["binding"].values()))
        assert (
            b["bound_case"].startswith(f"{eng}/{fam}-seed")
            and len(b["committed_sha256"]["case.json"]) == 64
        )
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
            _load_groups(tmp_path / f"c{engine}e", [bad], engine)
        bad = maker(tmp_path / f"{engine}-n", 150, False, 5, 0.0, label=(150, True))
        with pytest.raises(M.IddError, match="mislabeled"):
            _load_groups(tmp_path / f"c{engine}n", [bad], engine)


def test_lineage_derives_the_group_from_the_native_input(tmp_path: Path) -> None:
    runs = [_topas_run(tmp_path / f"t{s}", 200, False, 1 + s, 0.0) for s in range(3)]
    g = _load_groups(tmp_path, runs, "topas")
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
        _load_groups(tmp_path / "c1", [*ok, odd], "topas")
    # a run whose engine summary contradicts the declared histories fails the fingerprint
    bad = _mc_run(tmp_path / "m", 150, True, 3, 0.0)
    (bad / "stdout.txt").write_text("Nbr primaries simulated: 5\n")
    mf = json.loads((bad / "transfer-manifest.json").read_text())
    import hashlib

    blob = (bad / "stdout.txt").read_bytes()
    mf["files"]["stdout.txt"] = {"bytes": len(blob), "sha256": hashlib.sha256(blob).hexdigest()}
    (bad / "transfer-manifest.json").write_text(json.dumps(mf))
    with pytest.raises(M.IddError, match="summary histories"):
        _load_groups(tmp_path / "c2", [bad], "mcsquare")


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
    mc = [
        _mc_run(tmp_path / f"m{e}{n}{s}", e, n, 100 * e + 10 * n + s + 1, 0.001 * (s - 1))
        for e in (150, 200)
        for n in (True, False)
        for s in range(3)
    ]
    out = tmp_path / "out.json"
    cases = _freeze(tmp_path / "cases", [*topas, *mc])
    argv = ["--ionmc-dir", str(ion), "--topas-runs", *map(str, topas), "--mcsquare-runs",
            *map(str, mc), "--cases-dir", str(cases), "--output", str(out)]  # fmt: skip
    assert M.main(argv) == 1
    doc = json.loads(out.read_text())
    assert doc["pass"] is False
    fp = doc["engine_inputs"]["topas"]["150-on"]["fingerprints"]
    assert len(fp) == 3 and all(len(v["config_sha256"]) == 64 for v in fp.values())


# -- binding to the frozen committed cases (C20 G1) -------------------------------------------------
def _three(tmp: Path, maker, e: int, nuc: bool, **kw0) -> list[Path]:
    return [maker(tmp / f"r{s}", e, nuc, 1 + s, 0.0, **(kw0 if s == 0 else {})) for s in range(3)]


def _with_modules(nuc: bool, line: str) -> bytes:
    base = _topas_input(150, nuc, 1)
    return base.replace((FULL_MODULES if nuc else EMONLY_MODULES).encode(), line.encode())


@pytest.mark.parametrize(
    ("maker", "engine", "nuc", "bad_input", "reason"),
    [
        (  # five modules: the frozen full set minus g4decay
            _topas_run, "topas", True,
            _with_modules(True, FULL_MODULES.replace(' "g4decay"', "").replace("= 6", "= 5")),
            "neither the frozen full set",
        ),
        (  # opt4 plus one hadronic module, labelled emonly
            _topas_run, "topas", False,
            _with_modules(False, 'sv:Ph/Default/Modules = 2 "g4em-standard_opt4" "g4h-elastic_HP"'),
            "neither the frozen full set",
        ),
        (  # a different hadronic physics list
            _topas_run, "topas", True,
            _with_modules(True, FULL_MODULES.replace("QGSP_BIC_HP", "QGSP_BIC")),
            "neither the frozen full set",
        ),
        (_topas_run, "topas", True, _topas_input(150, True, 1).replace(b"Spread = 0", b"Spread = 0.5"),
         "BeamEnergySpread"),
        (_topas_run, "topas", True,
         _topas_input(150, True, 1).replace(b"CutForAllParticles", b"CutForGamma"), "CutForAllParticles"),
        (  # one secondary flag False in an "on" case
            _mc_run, "mcsquare", True,
            _mc_input(True, 1).replace(b"Simulate_Secondary_Alphas True", b"Simulate_Secondary_Alphas False"),
            "not all equal",
        ),
        (_mc_run, "mcsquare", True, _mc_input(True, 1).replace(b"Simulate_Secondary_Protons True\n", b""),
         "exactly one Simulate_Secondary_Protons"),
    ],
)  # fmt: skip
def test_binding_refuses_wrong_frozen_physics(
    tmp_path, maker, engine, nuc, bad_input, reason
) -> None:
    """The run is byte-identical to a (equally wrong) committed case, but the frozen physics is not."""
    bad = maker(tmp_path / "bad", 150, nuc, 1, 0.0, input_bytes=bad_input)
    with pytest.raises(M.IddError, match=reason):
        _load_groups(tmp_path / "c", [bad], engine)


def test_binding_refuses_non_integral_or_unequal_energy(tmp_path: Path) -> None:
    bad = _topas_run(tmp_path / "t", 150, True, 1, 0.0,
                     input_bytes=_topas_input(150, True, 1).replace(b"= 150 MeV", b"= 150.5 MeV"))  # fmt: skip
    with pytest.raises(M.IddError, match="not unique/integral"):
        _load_groups(tmp_path / "ct", [bad], "topas")
    bad = _mc_run(tmp_path / "m", 150, True, 1, 0.0)
    plan = (bad / "inputs" / "Plan.txt").read_bytes().replace(b"\n150\n", b"\n150.5\n")
    (bad / "inputs" / "Plan.txt").write_bytes(plan)
    mf = json.loads((bad / "transfer-manifest.json").read_text())
    import hashlib

    mf["files"]["inputs/Plan.txt"] = {
        "bytes": len(plan),
        "sha256": hashlib.sha256(plan).hexdigest(),
    }
    (bad / "transfer-manifest.json").write_text(json.dumps(mf))
    with pytest.raises(M.IddError):
        _load_groups(tmp_path / "cm", [bad], "mcsquare")


@pytest.mark.parametrize(
    ("maker", "engine", "old", "new"),
    [
        (_topas_run, "topas", b"i:So/Beam/NumberOfHistoriesInRun", b"# moved\ni:So/Beam/NumberOfHistoriesInRun"),
        (_topas_run, "topas", b"d:So/Beam/BeamEnergy = 150 MeV", b"d:So/Beam/BeamEnergy = 150 MeV\ns:So/Beam/Particle = \"proton\""),
        (_topas_run, "topas", b"CutForAllParticles = 0.05 mm", b"CutForAllParticles = 0.5 mm"),
        (_mc_run, "mcsquare", b"Num_Primaries 100000", b"Num_Primaries 100000\nE_Cut_Pro 5"),
    ],
)  # fmt: skip
def test_binding_refuses_a_modified_input_line(tmp_path, maker, engine, old, new) -> None:
    """Source/geometry/cut/scorer lines changed in the run but not in the committed case: no match."""
    runs = _three(tmp_path / "ok", maker, 150, True)
    cases = _freeze(tmp_path / "cases", runs)
    base = (runs[0] / "inputs" / ("input.txt" if engine == "topas" else "config.txt")).read_bytes()
    assert old in base
    bad = maker(tmp_path / "bad", 150, True, 1, 0.0, input_bytes=base.replace(old, new))
    with pytest.raises(M.IddError, match="match 0 committed"):
        M.load_engine_groups([bad, *runs[1:]], engine, {150: NZ}, cases)


def test_binding_refuses_extra_and_missing_aux_files(tmp_path: Path) -> None:
    runs = _three(tmp_path / "ok", _mc_run, 150, True)
    cases = _freeze(tmp_path / "cases", runs)
    # extra aux input in the run (manifested), absent from the committed case
    extra = runs[0] / "inputs" / "BDL_mono.txt"
    extra.write_bytes(b"bdl")
    mf = json.loads((runs[0] / "transfer-manifest.json").read_text())
    import hashlib

    mf["files"]["inputs/BDL_mono.txt"] = {"bytes": 3, "sha256": hashlib.sha256(b"bdl").hexdigest()}
    (runs[0] / "transfer-manifest.json").write_text(json.dumps(mf))
    with pytest.raises(M.IddError, match="match 0 committed"):
        M.load_engine_groups(runs, "mcsquare", {150: NZ}, cases)
    # committed case carries an aux file the run lacks
    runs2 = _three(tmp_path / "ok2", _mc_run, 150, True)
    cases2 = _freeze(tmp_path / "cases2", runs2)
    next(p for p in (cases2 / "mcsquare").iterdir() if p.name.endswith("seed1")).joinpath(
        "CT.mhd"
    ).write_bytes(b"x")
    with pytest.raises(M.IddError, match="match 0 committed"):
        M.load_engine_groups(runs2, "mcsquare", {150: NZ}, cases2)


def test_binding_refuses_missing_and_ambiguous_cases(tmp_path: Path) -> None:
    runs = _three(tmp_path / "ok", _topas_run, 150, True)
    empty = _freeze(tmp_path / "empty", [])
    with pytest.raises(M.IddError, match="match 0 committed"):
        M.load_engine_groups(runs, "topas", {150: NZ}, empty)
    cases = _freeze(tmp_path / "cases", runs)
    src = next(p for p in (cases / "topas").iterdir() if p.name.endswith("seed1"))
    import shutil

    shutil.copytree(src, cases / "topas" / "syn-150-full-copy-seed9")
    with pytest.raises(M.IddError, match="match 2 committed"):
        M.load_engine_groups(runs, "topas", {150: NZ}, cases)


def test_binding_refuses_one_case_bound_twice_and_missing_v5_block(tmp_path: Path) -> None:
    runs = _three(tmp_path / "ok", _topas_run, 150, True)
    cases = _freeze(tmp_path / "cases", runs)
    with pytest.raises(M.IddError, match="already bound"):
        M.load_engine_groups([runs[0], runs[0], runs[1]], "topas", {150: NZ}, cases)
    no_v5 = _freeze(tmp_path / "nov5", runs, strip_v5=True)
    # the run's own case.json still has the block, so it no longer matches the stripped committed case
    with pytest.raises(M.IddError, match="match 0 committed"):
        M.load_engine_groups(runs, "topas", {150: NZ}, no_v5)
    # run frozen WITHOUT a V5 block (committed == run): refused for the missing block
    bare = []
    for s in range(3):
        rd = _topas_run(tmp_path / f"bare{s}", 150, True, 1 + s, 0.0)
        case = json.loads((rd / "inputs" / "case.json").read_text())
        case.pop("v5")
        blob = json.dumps(case).encode()
        (rd / "inputs" / "case.json").write_bytes(blob)
        import hashlib

        mf = json.loads((rd / "transfer-manifest.json").read_text())
        mf["files"]["inputs/case.json"] = {
            "bytes": len(blob),
            "sha256": hashlib.sha256(blob).hexdigest(),
        }
        (rd / "transfer-manifest.json").write_text(json.dumps(mf))
        bare.append(rd)
    with pytest.raises(M.IddError, match="no V5 block"):
        _load_groups(tmp_path / "cbare", bare, "topas")


def test_binding_requires_one_family_per_group(tmp_path: Path) -> None:
    runs = _three(tmp_path / "ok", _topas_run, 150, True)
    cases = _freeze(tmp_path / "cases", runs)
    d = next(p for p in (cases / "topas").iterdir() if p.name.endswith("seed1"))
    d.rename(d.with_name("other-family-seed1"))
    with pytest.raises(M.IddError, match="several families"):
        M.load_engine_groups(runs, "topas", {150: NZ}, cases)
    d = cases / "topas" / "other-family-seed1"
    d.rename(d.with_name("not-a-seed-variant"))
    with pytest.raises(M.IddError, match="seed<k>"):
        M.load_engine_groups(runs, "topas", {150: NZ}, cases)


def test_build_verdict_without_engine_runs_is_refused(tmp_path: Path) -> None:
    ion = tmp_path / "ion"
    ion.mkdir()
    for e in (150, 200):
        for nuc in (True, False):
            _partial(
                ion / f"v5-{e}-{'on' if nuc else 'off'}.json",
                e,
                nuc,
                noisy(_scaled(e, nuc), 20, 0.03, e),
            )
    topas = [_topas_run(tmp_path / f"t{e}{n}{s}", e, n, 100 * e + 10 * n + s + 1, 0.0)
             for e in (150, 200) for n in (True, False) for s in range(3)]  # fmt: skip
    cases = _freeze(tmp_path / "cases", topas)
    with pytest.raises(M.IddError, match="no mcsquare reference runs"):
        M.build_verdict(ion, topas, [], cases)
    with pytest.raises(M.IddError, match="no topas reference runs"):
        M.build_verdict(ion, [], [], cases)
