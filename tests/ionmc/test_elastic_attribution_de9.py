"""D-E9 investigation harness (exploratory, report-only): B0 independent p-p level reconstruction, B3 variant family, B2 genealogy split,
the de9 record and the seed-override guard, on synthetic inputs."""

# ruff: noqa: E501
from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from validation.scripts.elastic_attribution import common as c
from validation.scripts.elastic_attribution import de9_pp_level as b0


# ---------------------------------------------------------------------------------------------------------------------------------
# B0
# ---------------------------------------------------------------------------------------------------------------------------------
def ef(x: float) -> str:
    """ENDF-6 float field (11 characters)."""
    m, e = f"{x:.6e}".split("e")
    return f"{m}{int(e):+d}".rjust(11)


def _line(vals: list[Any], mat: int, mf: int, mt: int, seq: int) -> str:
    f = "".join(
        v if isinstance(v, str) else (str(v).rjust(11) if isinstance(v, int) else ef(v))
        for v in vals
    )
    return f"{f:<66}{mat:4d}{mf:2d}{mt:3d}{seq:5d}"


def _synthetic_endf(energies_ev: list[float], b0s: list[float], a0: complex = 0j) -> str:
    """H-1 like MF1/MT451 header (AWI) and an MF6/MT2 LAW=5 LTP=1 LIDP=1 section with NL = 0 records (b_0, Re a_0, Im a_0)."""
    mat, lines = 125, []
    lines.append(_line([1001.0, 0.99917, 0, 0, 0, 0], mat, 1, 451, 1))
    lines.append(_line([0.0, 0.0, 0, 0, 0, 0], mat, 1, 451, 2))
    lines.append(_line([0.99862, 0.0, 0, 0, 0, 0], mat, 1, 451, 3))
    seq = [0]

    def add(vals: list[Any]) -> None:
        seq[0] += 1
        lines.append(_line(vals, mat, 6, 2, seq[0]))

    def block(vals: list[float]) -> None:
        for i in range(0, len(vals), 6):
            chunk: list[Any] = list(vals[i : i + 6])
            add(chunk + [0.0] * (6 - len(chunk)))

    add([1001.0, 0.99917, 0, 0, 1, 0])
    add([1.0, 0.99662, 0, 5, 1, 2])
    block([2.0, 2.0])
    block([1.0e3, 1.0, 1.5e8, 1.0])
    add([0.5, 0.0, 1, 0, 1, len(energies_ev)])
    block([float(len(energies_ev)), 2.0])
    for e, b in zip(energies_ev, b0s, strict=True):
        add([0.0, e, 1, 0, 3, 0])
        block([b, a0.real, a0.imag])
    return "\n".join(lines) + "\n"


def _sec(
    b0s: list[float] | None = None, energies: list[float] | None = None, a0: complex = 0j
) -> dict[str, Any]:
    return b0.parse_h1_mt2(
        _synthetic_endf(energies or [1.0e7, 1.0e8, 1.5e8], b0s or [2.0e-3, 2.0e-3, 2.0e-3], a0)
    )


def test_endf_fields_and_section_parser() -> None:
    assert (
        b0.endf_float(" 1.001000+3") == 1001.0
        and b0.endf_float("-5.000000-1") == -0.5
        and b0.endf_float("           ") == 0.0
        and b0.endf_float("1.5e8") == 1.5e8
    )
    sec = _sec()
    assert (
        sec["spi"] == 0.5 and sec["lidp"] == 1 and sec["awr"] == 0.99917 and sec["awi"] == 0.99862
    )
    assert [r["e_ev"] for r in sec["records"]] == [1.0e7, 1.0e8, 1.5e8] and sec["records"][0][
        "nl"
    ] == 0
    with pytest.raises(ValueError, match="LTP"):
        b0.parse_h1_mt2(
            _synthetic_endf([1e7, 1e8], [1e-3, 1e-3]).replace(
                "        1          0          3          0",
                "       12          0          3          0",
                1,
            )
        )


def test_nuclear_only_record_integrates_analytically_at_a_node() -> None:
    """a_l = 0: NI = (1/2) b_0 P_0 exactly, so sigma = 2 pi (0.96) (b_0 / 2) at a node energy (the ratio interpolation is exact there)."""
    sec = _sec()
    assert b0.sigma_ni_barn(sec, 100.0) == pytest.approx(
        2.0 * math.pi * 0.96 * 0.5 * 2.0e-3, rel=1e-12
    )
    assert b0.sigma_ni_barn(sec, 10.0) == pytest.approx(
        2.0 * math.pi * 0.96 * 0.5 * 2.0e-3, rel=1e-12
    )
    assert b0.sigma_ni_barn(sec, 100.0, n=256) == pytest.approx(
        b0.sigma_ni_barn(sec, 100.0, n=1024), rel=1e-12
    )


def test_ratio_interpolation_between_nodes_and_the_coulomb_symmetry() -> None:
    sec = _sec([2.0e-3, 4.0e-3, 4.0e-3])
    mu = np.linspace(0.0, 0.96, 7)
    e = 5.5e7
    r1 = b0.ratio_at(sec, e, mu)
    k1, eta1 = b0.coulomb_k_eta(1.0e7, sec["awr"], sec["awi"])
    k2, eta2 = b0.coulomb_k_eta(1.0e8, sec["awr"], sec["awi"])
    f = (e - 1.0e7) / (1.0e8 - 1.0e7)
    want = (1 - f) * (1 + 0.5 * 2.0e-3 / b0.sigma_coulomb(mu, k1, eta1)) + f * (
        1 + 0.5 * 4.0e-3 / b0.sigma_coulomb(mu, k2, eta2)
    )
    assert np.allclose(r1, want, rtol=1e-12)
    x = np.array([0.1, 0.5, 0.9])
    assert np.allclose(
        b0.sigma_coulomb(x, k1, eta1), b0.sigma_coulomb(-x, k1, eta1), rtol=1e-12
    )  # identical particles: symmetric in mu
    with pytest.raises(ValueError, match="outside"):
        b0.ratio_at(sec, 2.0e8, mu)


def test_interference_term_follows_equation_6_14() -> None:
    """With a_0 = (0.01, 0) and b_0 = 0: sigma_e - sigma_c = -2 eta/(1-mu^2) Re{ [(1+mu)e1 + (1-mu)e2] a_0/2 } (P_0 = 1)."""
    sec = _sec([0.0, 0.0, 0.0], a0=0.01 + 0.0j)
    rec = sec["records"][1]
    k, eta = b0.coulomb_k_eta(rec["e_ev"], sec["awr"], sec["awi"])
    mu = np.array([0.2, 0.7])
    got = b0.sigma_e(mu, rec, k, eta) - b0.sigma_coulomb(mu, k, eta)
    e1, e2 = np.exp(1j * eta * np.log((1 - mu) / 2)), np.exp(1j * eta * np.log((1 + mu) / 2))
    want = -2 * eta / (1 - mu**2) * np.real(((1 + mu) * e1 + (1 - mu) * e2) * 0.01 / 2)
    assert np.allclose(got, want, rtol=1e-12, atol=0.0)


def test_dip_rule_and_the_pre_declared_reading() -> None:
    e = np.arange(100.0, 151.0)
    flat = np.linspace(1.0, 0.9, e.size)
    assert not b0.dip_test(e, flat)["persists"]
    d = b0.dip_test(e, np.where((e > 120) & (e < 135), 0.8, np.linspace(1.0, 0.95, e.size)))
    assert (
        d["persists"]
        and 120 < d["min_energy_mev"] < 135
        and d["min_over_lower_endpoint_minus_1"] < -0.01
    )
    shallow = np.where(
        (e > 120) & (e < 135),
        0.995 * np.linspace(1.0, 0.99, e.size),
        np.linspace(1.0, 0.99, e.size),
    )
    assert not b0.dip_test(e, shallow)["persists"]  # less than 1 % below the endpoints


def _table_for(
    sec: dict[str, Any], grid: np.ndarray, scale: dict[float, float] | None = None
) -> np.ndarray:
    return np.array(
        [
            (b0.sigma_ni_barn(sec, float(x)) if x <= 150.0 else 1.0)
            * (scale or {}).get(float(x), 1.0)
            for x in grid
        ]
    )


def test_compare_match_feature_and_defect_readings(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        b0, "N_GAUSS", 96
    )  # synthetic integrands are smooth; the production value is 1024
    grid = np.arange(
        90.0, 176.0, 1.0
    )  # the builder's 1 MeV grid: linear midpoints of the table match the reconstruction
    dip_b0 = [4.0e-3, 4.0e-3, 4.0e-3, 4.0e-3, 4.0e-3]
    sec_dip = _sec(
        [4.0e-3, 4.0e-3, 1.0e-3, 1.0e-3, 4.0e-3], energies=[1.0e7, 1.0e8, 1.25e8, 1.3e8, 1.5e8]
    )
    sec_flat = _sec(dip_b0, energies=[1.0e7, 1.0e8, 1.25e8, 1.3e8, 1.5e8])
    res = b0.compare(sec_dip, grid, _table_for(sec_dip, grid))
    assert (
        res["match"]["all_within"]
        and res["dip"]["persists"]
        and res["reading"].startswith("evaluation feature")
    )
    assert (
        set(res["named_nodes"]) == {"100", "125", "130", "140", "150"}
        and "100.5" in res["midpoints"]
        and "175" not in res["nodes"]
    )
    assert res["levels_over_125mev"]["125"] == {"rec": 1.0, "table": 1.0}
    bad = b0.compare(
        sec_dip, grid, _table_for(sec_dip, grid, {140.0: 1.01})
    )  # a 1 % node mismatch beyond 1e-3
    assert (
        not bad["match"]["all_within"]
        and bad["reading"].startswith("builder defect")
        and not bad["nodes"]["140"]["within"]
    )
    nodip = b0.compare(sec_flat, grid, _table_for(sec_flat, grid))
    assert (
        nodip["match"]["all_within"]
        and not nodip["dip"]["persists"]
        and nodip["reading"].startswith("builder defect")
    )


def test_projection_has_no_absolute_cross_section_and_the_leak_check_fires(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(b0, "N_GAUSS", 96)
    grid = np.arange(100.0, 151.0, 5.0)
    sec = _sec(
        [4.0e-3, 4.0e-3, 1.0e-3, 1.0e-3, 4.0e-3], energies=[1.0e7, 1.0e8, 1.25e8, 1.3e8, 1.5e8]
    )
    res = b0.compare(sec, grid, _table_for(sec, grid))
    absolutes = res.pop("_absolutes")
    doc = {
        "banner": c.BANNER,
        "disclosure": b0.DISCLOSURE,
        **res,
        "provenance": {"quadrature": "x"},
    }
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), b0.GIT_ALLOW)
    assert b0.leak_violations(proj, absolutes) == []
    assert "rec_over_table" in json.dumps(proj) and len(absolutes) > 20
    leaked = json.loads(json.dumps(proj))
    leaked["nodes"]["125"]["rec_over_table"] = absolutes[0] * 1e3  # barn -> mb of a node value
    assert b0.leak_violations(leaked, absolutes)
    assert "Reading (pre-declared)" in b0.render_report(proj)


def test_b0_is_independent_of_the_builder_endf6_and_law5_code() -> None:
    tree = ast.parse(Path(b0.__file__).read_text())
    imported: list[str] = []
    for n in ast.walk(tree):
        if isinstance(n, ast.ImportFrom) and n.module:
            imported += [f"{n.module}.{a.name}" for a in n.names]
        elif isinstance(n, ast.Import):
            imported += [a.name for a in n.names]
    ionmc = [m for m in imported if "ionmc" in m]
    assert ionmc == ["ionmc.nuclear.bgg"], (
        ionmc
    )  # only the Geant4 p-p S(E) formula (carry-up check); no builder, endf6 or LAW=5 code
    assert not [m for m in imported if any(w in m for w in ("elastic_build", "endf6", "law5"))]


def test_real_data_parses_and_reproduces_the_100_mev_node() -> None:
    """Parser regression on the raw LA150 H-1 file: skipped when the pinned cache is absent (no verdict is asserted here)."""
    import zipfile

    try:
        zp, _ = b0.locate_raw(c.CACHE_DIR)
        z = np.load(
            c.CACHE_DIR / "derived" / f"elastic-proton-{c.TABLE_ID}.npz", allow_pickle=False
        )
    except (FileNotFoundError, KeyError):
        pytest.skip("pinned cache / LA150 object absent")
    sec = b0.parse_h1_mt2(zipfile.ZipFile(zp).read(b0.MEMBER).decode("ascii"))
    assert len(sec["records"]) == 131 and sec["records"][-1]["e_ev"] == 1.5e8
    k = int(np.argmin(np.abs(z["grid_e_mev"] - 100.0)))
    assert b0.sigma_ni_barn(sec, 100.0) == pytest.approx(float(z["sigma_barn"][0][k]), rel=1e-9)


# ---------------------------------------------------------------------------------------------------------------------------------
# B3 variant family
# ---------------------------------------------------------------------------------------------------------------------------------
def _base_arrays() -> dict[str, np.ndarray]:
    grid = np.array([100.0, 125.0, 130.0, 140.0, 150.0, 175.0, 200.0, 250.0])
    sig = np.array([[5.0, 4.6, 3.8, 4.0, 4.5, 4.1, 3.7, 3.4], [1.0] * 8])
    return {
        "grid_e_mev": grid,
        "sigma_barn": sig,
        "edges_mu": np.zeros((2, 8, 5)),
        "valid": np.ones((2, 8), dtype=np.int8),
    }


def test_pp_level_scale_variants_scale_only_the_h1_nodes_at_or_above_150_mev() -> None:
    from validation.scripts.elastic_attribution import make_variant_tables as mv

    base = _base_arrays()
    for name, s in (("pp-level-0.90-ge150", 0.90), ("pp-level-0.95-ge150", 0.95)):
        a, stats = mv._build_arrays(name, base, ("H-1", "C-12"))
        want = base["sigma_barn"][0].copy()
        want[4:] *= s
        assert (
            np.allclose(a["sigma_barn"][0], want, rtol=0, atol=0)
            and np.array_equal(a["sigma_barn"][1], base["sigma_barn"][1])
            and stats == {}
        )
        assert np.array_equal(a["edges_mu"], base["edges_mu"]) and np.array_equal(
            a["valid"], base["valid"]
        )  # shape unchanged
    assert mv.PP_LEVEL_SCALES == {
        "pp-level-0.90-ge150": 0.90,
        "pp-level-0.95-ge150": 0.95,
    }  # literal factors, never selected


def test_monotone_variant_carries_the_125_mev_level_by_the_tables_own_s_of_e() -> None:
    from ionmc.nuclear import bgg
    from validation.scripts.elastic_attribution import make_variant_tables as mv

    base = _base_arrays()
    a, _ = mv._build_arrays(mv.PP_LEVEL_MONOTONE, base, ("H-1", "C-12"))
    s = a["sigma_barn"][0]
    assert np.array_equal(s[:2], base["sigma_barn"][0][:2])  # nodes <= 125 MeV unchanged
    for k in range(2, 8):
        assert s[k] == pytest.approx(
            base["sigma_barn"][0][1]
            * bgg.pp_bgg_mb(float(base["grid_e_mev"][k]))
            / bgg.pp_bgg_mb(125.0),
            rel=1e-14,
        )
    assert (
        np.all(np.diff(s[1:]) <= 0.0) and base["sigma_barn"][0][4] > base["sigma_barn"][0][2]
    )  # the input has the upturn, the variant does not
    assert np.array_equal(a["edges_mu"], base["edges_mu"]) and np.array_equal(
        a["sigma_barn"][1], base["sigma_barn"][1]
    )
    bad = {**base, "grid_e_mev": base["grid_e_mev"] + 0.5}
    with pytest.raises(AssertionError, match="no 125 MeV node"):
        mv._build_arrays(mv.PP_LEVEL_MONOTONE, bad, ("H-1", "C-12"))


def test_variants_are_registered_with_parameters_and_the_no_inversion_rule() -> None:
    import inspect

    from validation.scripts.elastic_attribution import make_variant_tables as mv

    assert set(mv.PP_LEVEL_VARIANTS) <= set(mv.VARIANTS)
    for v in mv.PP_LEVEL_VARIANTS:
        assert mv._parameters(v)
    src = inspect.getsource(mv._build_arrays) + inspect.getsource(mv._parameters)
    assert (
        "residual" not in src and "TOPAS" not in src and "closure" not in src
    )  # no level or scale is derived from a residual or a reference
    import validation.scripts.elastic_attribution.de9_record as rec

    assert not [
        n for n in dir(rec) if any(w in n.lower() for w in ("invert", "closing", "required_scale"))
    ]


# ---------------------------------------------------------------------------------------------------------------------------------
# seed override, genealogy split
# ---------------------------------------------------------------------------------------------------------------------------------
def test_seed_override_is_guarded_to_the_registered_block_and_once_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from validation.scripts.elastic_attribution import de3_probes as d3

    for e, seeds in ((200, [20992000, 20992007, 20992101]), (150, [20992010, 20992013, 20992100])):
        for s in seeds:
            d3.check_seed_override(e, s)
    for e, s in (
        (200, 20992008),
        (200, 20992010),
        (150, 20992000),
        (150, 20992014),
        (200, 20990101),
        (150, 20992101),
        (200, 20992200),
    ):
        with pytest.raises(SystemExit, match="not registered"):
            d3.check_seed_override(e, s)
    assert (
        d3.run_name(200, "pinned", 3) == "de3-200-pinned-r3.json"
        and d3.run_name(150, "all") == "de3-150-all.json"
    )
    out, att = tmp_path / "de9", tmp_path / "attempts.json"
    out.mkdir()
    calls: list[dict[str, Any]] = []

    def fake(
        variant: str, energy: int, n: int, batches: int, lateral: bool, out: Path, **kw: Any
    ) -> dict[str, Any]:
        calls.append(kw)
        (out / d3.run_name(energy, variant, kw.get("realisation"))).write_text("{}")
        return {"valid": True}

    monkeypatch.setattr(d3, "_run_one", fake)
    with pytest.raises(SystemExit, match="not registered"):
        d3.run_one(
            "pinned", 200, 100000, 20, False, out, att, seed_override=20990101, realisation=0
        )  # refused before any attempt record
    assert not att.exists() and not calls
    d3.run_one(
        "pinned",
        200,
        100000,
        20,
        False,
        out,
        att,
        seed_override=20992000,
        realisation=0,
        genealogy_split=True,
    )
    assert calls == [{"seed_override": 20992000, "realisation": 0, "genealogy_split": True}]
    log = json.loads(att.read_text())
    assert (
        [r["outcome"] for r in log] == ["started", "published-valid"]
        and log[0]["output"].endswith("de3-200-pinned-r0.json")
        and log[0]["args"]["seed_override"] == 20992000
    )
    with pytest.raises(SystemExit, match="exists"):
        d3.run_one(
            "pinned", 200, 100000, 20, False, out, att, seed_override=20992000, realisation=0
        )  # a realisation is never repeated
    d3.run_one(
        "pinned", 200, 100000, 20, False, out, att, seed_override=20992001, realisation=1
    )  # the next realisation is a different target


def test_genealogy_split_uses_existing_species_and_generation_channels() -> None:
    from ionmc.scoring import TallyRequest

    sec = TallyRequest(
        "pp_secondary_proton", "dose", "edep", species=("proton",), generation="secondary"
    )
    pri = TallyRequest("primary_proton", "dose", "edep", species=("proton",), generation="primary")
    assert (sec.species, sec.generation, pri.generation) == (("proton",), "secondary", "primary")


# ---------------------------------------------------------------------------------------------------------------------------------
# de9 record: B1, B3, A1, B2, Lead A sensitivity table
# ---------------------------------------------------------------------------------------------------------------------------------
NZ = 400
PL_EM = 6.0


def _curves(
    plateau: float, peak_over_plateau: float, jitter: float = 0.0, scale_b: float = 1.0
) -> list[list[float]]:
    rows = []
    for i in range(20):
        r = np.full(NZ, plateau * (1.0 + 4e-4 * math.cos(i + jitter)))
        r[300] = plateau * peak_over_plateau * (1.0 + 1e-3 * math.sin(2 * i + jitter))
        rows.append((r * scale_b).tolist())
    return rows


REUSED_V = ("all", "pp-only", "pA-only", "emonly")


def _expected() -> dict[str, dict[str, str]]:
    return {
        v: {"table_id": f"T-{v}", "npz_sha256": f"S-{v}"}
        for v in (
            "all",
            "pp-only",
            "pA-only",
            "emonly",
            "pinned",
            "m2prime",
            "pp-level-0.90-ge150",
            "pp-level-0.95-ge150",
            "pp-level-monotone",
        )
    }


def _doc(energy: int, variant: str, plateau: float, seed: int, **extra: Any) -> dict[str, Any]:
    d: dict[str, Any] = {
        "energy": energy, "variant": variant, "seed": seed, "idd_batches": _curves(plateau, 5.0, jitter=seed % 7), "valid": True, "relative_residual": 1e-16,
        "n": 100000, "n_batches": 20, "mode": "emonly" if variant == "emonly" else "emel", "backend": "warp-cpu", "precision": "float64", "elastic_table_id": f"T-{variant}",
    }  # fmt: skip
    if (
        "realisation" in extra or variant not in REUSED_V
    ):  # reused D-E3 documents predate the npz field
        d["elastic_table_npz_sha256"] = f"S-{variant}"
    return {**d, **extra}


GAINS = {
    "all": 0.025,
    "pp-only": 0.020,
    "pA-only": 0.005,
    "m2prime": 0.0213,
    "pp-level-0.90-ge150": 0.0225,
    "pp-level-0.95-ge150": 0.0238,
    "pp-level-monotone": 0.0242,
}


def _write_inputs(
    tmp: Path, a1_gains: dict[int, list[float]] | None = None, with_b2: bool = True
) -> dict[str, Path]:
    de3, de9 = tmp / "de3", tmp / "de9"
    de3.mkdir()
    de9.mkdir()
    seeds = {150: 20990100, 200: 20990101}
    for e in (150, 200):
        (de3 / f"de3-{e}-emonly.json").write_text(json.dumps(_doc(e, "emonly", PL_EM, seeds[e])))
        for v, g in GAINS.items():
            (de3 / f"de3-{e}-{v}.json").write_text(
                json.dumps(_doc(e, v, PL_EM * (1 + g), seeds[e]))
            )
    rng = np.random.default_rng(1)
    for e, k in ((200, 8), (150, 4)):
        base = {200: 20992000, 150: 20992010}[e]
        gs = (a1_gains or {}).get(e) or list(0.0255 + 0.001 * rng.standard_normal(k))
        for r in range(k):
            for v, pl in (("pinned", PL_EM * (1 + gs[r])), ("emonly", PL_EM)):
                (de9 / f"de3-{e}-{v}-r{r}.json").write_text(
                    json.dumps(_doc(e, v, pl, base + r, realisation=r))
                )
        if with_b2:
            pl = PL_EM * 1.025
            d = _doc(e, "pinned", pl, {150: 20992100, 200: 20992101}[e], realisation=100)
            d["recoil_idd_batches"] = _curves(0.002, 1.0)
            d["genealogy_split"] = {
                "pp_secondary_proton": _curves(0.1, 1.0),
                "primary_proton": _curves(5.0, 1.0),
            }
            (de9 / f"de3-{e}-pinned-r100.json").write_text(json.dumps(d))
    steps = {
        "16-v11-compare": {
            "verdict": {
                "energies": {
                    str(e): {
                        "plateau_gain": {
                            "ionmc": {"value": 0.02672217441583613, "se": 0.00077},
                            "reference": {"value": 0.018419142339446415},
                        }
                    }
                    for e in (150, 200)
                }
            }
        }
    }
    (tmp / "archive.json").write_text(json.dumps({"steps": steps}))
    hashes = tmp / "hashes.md"
    hashes.write_text(
        "Input files (sha256):\n\n"
        + "".join(
            f"- `x/{p.name}`  {hashlib.sha256(p.read_bytes()).hexdigest()}\n"
            for v in REUSED_V
            for p in sorted(de3.glob(f"de3-*-{v}.json"))
        )
    )
    b0 = tmp / "de9-b0.json"
    b0.write_text(
        json.dumps(
            {
                "match": {"all_within": True},
                "dip": {"persists": True},
                "reading": "r",
                "named_nodes": {},
                "provenance": {"git_head": "a" * 40, "git_clean": True},
            }
        )
    )
    return {"de3": de3, "de9": de9, "archive": tmp / "archive.json", "hashes": hashes, "b0": b0}


def _build(inp: dict[str, Path], monkeypatch: pytest.MonkeyPatch, **kw: Any) -> dict[str, Any]:
    from validation.scripts.elastic_attribution import de9_record as rec

    monkeypatch.setattr(
        c, "repo_state", lambda repo=None: {"git_head": "a" * 40, "git_clean": True}
    )
    return rec.build_document(
        inp["de3"],
        inp["de9"],
        kw["b0"] if "b0" in kw else inp["b0"],
        archive=inp["archive"],
        absolutes=kw.get("absolutes"),
        expected_tables=kw.get("expected", _expected()),
        hash_sources=(inp["hashes"],),
        allow_incomplete=kw.get("allow_incomplete", False),
        freeze=kw.get("freeze"),
    )


def test_b1_b3_and_the_gain_definitions(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    doc = _build(_write_inputs(tmp_path), monkeypatch)
    r = doc["b1"]["200"]
    assert (
        r["gain_all"] == pytest.approx(0.025, abs=2e-4)
        and r["gain_pp_only"] == pytest.approx(0.020, abs=2e-4)
        and r["pp_share_of_sum"] == pytest.approx(0.8, abs=0.02)
    )
    assert r["additivity_gap"] == pytest.approx(0.0, abs=4e-4)
    b3 = doc["b3"]
    assert b3["pp-level-0.90-ge150"]["200"]["delta_gain"] == pytest.approx(0.0225 - 0.025, abs=2e-4)
    assert b3["pp-level-monotone"]["150"]["delta_gain"] == pytest.approx(0.0242 - 0.025, abs=2e-4)
    s = b3["slope"]["200"]
    assert s["dg_dln_s_0.90_to_0.95"] == pytest.approx(
        (0.0238 - 0.0225) / math.log(0.95 / 0.90), rel=0.05
    ) and s["dg_dln_s_0.95_to_1.00"] == pytest.approx((0.025 - 0.0238) / -math.log(0.95), rel=0.05)
    text = json.dumps(doc).lower()
    assert (
        "closing" not in text and "required" not in text and "residual" not in text
    )  # the record never reports a closing scale or level


def test_a1_reading_rule_is_applied_as_pre_registered(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    g200 = [0.0245, 0.0250, 0.0255, 0.0260, 0.0265, 0.0258, 0.0252, 0.0248]
    g150 = [0.0200, 0.0210, 0.0205, 0.0195]
    doc = _build(_write_inputs(tmp_path, {200: g200, 150: g150}), monkeypatch)
    a = doc["a1"]["200"]
    mu, sd = float(np.mean(g200)), float(np.std(g200, ddof=1))
    assert (
        a["k"] == 8
        and a["mean_gain"] == pytest.approx(mu, abs=2e-4)
        and a["sd"] == pytest.approx(sd, abs=2e-4)
        and a["se_mean"] == pytest.approx(a["sd"] / math.sqrt(8), rel=1e-12)
    )
    assert (a["rho_range_lo"], a["rho_range_hi"]) == (0.556, 1.418) and (
        doc["a1"]["150"]["rho_range_lo"],
        doc["a1"]["150"]["rho_range_hi"],
    ) == (0.342, 1.614)
    assert a["z_archived"] == pytest.approx(
        (0.02672217441583613 - a["mean_gain"]) / math.hypot(0.00077, a["se_mean"]), rel=1e-12
    )
    assert a["archived_reading"] in (
        "consistent with the realisation distribution",
        "outlier at 2 sigma",
    ) and a["per_run_se_calibration"].startswith("per-run SE")
    assert (abs(a["z_archived"]) <= 2) == (
        a["archived_reading"] == "consistent with the realisation distribution"
    )


def test_lead_a_sensitivity_table_is_labelled_and_never_compared_with_a_band(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    doc = _build(_write_inputs(tmp_path), monkeypatch)
    t = doc["a1_sensitivity"]
    delta = GAINS["m2prime"] - GAINS["all"]
    frozen = 0.02672217441583613 + delta - 0.018419142339446415
    assert (
        t["frozen"]["d"] == pytest.approx(frozen, abs=3e-4)
        and t["frozen"]["label"] == "the frozen value; the STOP stands"
    )
    label = f"sensitivity; frozen d = {t['frozen']['d']:.5f} stands"
    assert t["archived_minus_2se"]["d"] == pytest.approx(
        t["frozen"]["d"] - 2 * 0.00077, abs=1e-12
    ) and t["archived_plus_2se"]["d"] == pytest.approx(t["frozen"]["d"] + 2 * 0.00077, abs=1e-12)
    for k, row in t.items():
        assert (
            row["compared_with_band"] is False
            and "holds" not in row
            and "pass" not in json.dumps(row).lower()
        )
        if k != "frozen":
            assert row["label"] == label
    assert set(t) >= {
        "frozen",
        "archived_minus_2se",
        "archived_plus_2se",
        "de3_baseline",
        "pooled_baseline",
        "a1_mean",
    }
    w = 1 / 0.00077**2, 1 / doc["b1"]["200"]["gain_all_se"] ** 2
    assert t["pooled_baseline"]["se"] == pytest.approx(1 / math.sqrt(sum(w)), rel=1e-9)


def test_b2_share_and_the_dropped_case(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    doc = _build(_write_inputs(tmp_path), monkeypatch)
    r = doc["b2"]["200"]
    excess = PL_EM * 0.025
    assert (
        r["status"] == "recorded"
        and r["secondary_proton_share"] == pytest.approx(0.1 / excess, rel=1e-2)
        and r["recoil_share"] == pytest.approx(0.002 / excess, rel=1e-2)
    )
    assert r["remainder_share"] == pytest.approx(1 - (0.1 + 0.002) / excess, rel=1e-2)
    dropped = _build(
        _write_inputs(
            tmp_path / "x" if (tmp_path / "x").mkdir() is None else tmp_path, with_b2=False
        ),
        monkeypatch,
    )
    assert dropped["b2"]["150"]["status"].startswith("dropped")


def test_record_refuses_wrong_seeds_missing_documents_and_overwrites(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from validation.scripts.elastic_attribution import de9_record as rec

    inp = _write_inputs(tmp_path)
    p = inp["de9"] / "de3-200-pinned-r3.json"
    d = json.loads(p.read_text())
    d["seed"] = 20990101
    p.write_text(json.dumps(d))
    with pytest.raises(SystemExit, match="differs from the frozen configuration"):
        _build(inp, monkeypatch)
    p.unlink()
    with pytest.raises(SystemExit, match="is missing"):
        _build(inp, monkeypatch)
    inp2 = _write_inputs(tmp_path / "w" if (tmp_path / "w").mkdir() is None else tmp_path)
    monkeypatch.setattr(
        c, "repo_state", lambda repo=None: {"git_head": "a" * 40, "git_clean": True}
    )
    monkeypatch.setattr(rec, "check_freeze", lambda *a, **k: {})
    monkeypatch.setattr(rec, "hashes_at", lambda sha, repo=None: {})
    (tmp_path / "exp.json").write_text(json.dumps(_expected()))
    out = tmp_path / "out"
    (out / "git").mkdir(parents=True)
    (out / "git" / "de9.json").write_text("{}")
    with pytest.raises(SystemExit, match="already published"):
        rec.run(
            inp2["de3"],
            inp2["de9"],
            out,
            inp2["b0"],
            "REVIEW-0123abcd",
            "b" * 40,
            None,
            tmp_path / "exp.json",
        )


def test_projection_leak_check_and_report(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from validation.scripts.elastic_attribution import de9_record as rec

    doc = _build(_write_inputs(tmp_path), monkeypatch, absolutes=[123.456])
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), rec.GIT_ALLOW)
    assert c.git_safety_violations(proj, rec.GIT_ALLOW) == []
    md = rec.render_report(proj)
    assert (
        "not a re-evaluation" in md
        and "never compared with the band" in md
        and "sensitivity; frozen d =" in md
    )
    leaked = json.loads(json.dumps(proj))
    leaked["b1"]["150"]["gain_all"] = 123.456
    assert rec.leak_violations(leaked, [123.456])
    with pytest.raises(SystemExit, match="leak check failed"):
        _build(
            _write_inputs(tmp_path / "l" if (tmp_path / "l").mkdir() is None else tmp_path),
            monkeypatch,
            absolutes=[doc["b1"]["150"]["gain_pp_only"]],
        )


# ---------------------------------------------------------------------------------------------------------------------------------
# freeze binding of the record, registered seeds, B0 disclosure
# ---------------------------------------------------------------------------------------------------------------------------------
def _repo9(tmp: Path) -> tuple[Path, str, str]:
    import subprocess

    from validation.scripts.elastic_attribution import de9_record as rec

    repo = tmp / "repo"
    root = repo / rec.SCRIPT_REL
    root.mkdir(parents=True)
    for n in rec.DE9_CODE_FILES:
        (root / n).write_text(f"# {n}\n")

    def git(*a: str) -> str:
        return subprocess.run(
            ["git", "-c", "user.email=t@t", "-c", "user.name=t", *a],
            cwd=repo,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

    git("init", "-q")
    git("add", "-A")
    git("commit", "-q", "-m", "freeze")
    first = git("rev-parse", "HEAD")
    (repo / "other.txt").write_text("x\n")
    git("add", "-A")
    git("commit", "-q", "-m", "second")
    return repo, first, git("rev-parse", "HEAD")


def test_record_freeze_binding_in_a_temporary_git_repo(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de9_record as rec

    repo, first, head = _repo9(tmp_path)
    man = tmp_path / "de9-manifest.json"
    clean = {"git_head": head, "git_clean": True}
    assert set(rec.DE9_CODE_FILES) == {
        "de9_pp_level.py",
        "de9_record.py",
        "de3_probes.py",
        "make_variant_tables.py",
    }
    assert rec.hashes_at(first, repo) == rec.working_hashes(repo)
    assert rec.check_freeze("REVIEW-0123abcd", first, clean, man, repo) == {}
    with pytest.raises(SystemExit, match="not clean"):
        rec.check_freeze("REVIEW-0123abcd", first, {**clean, "git_clean": False}, man, repo)
    for bad in (None, "", "review-1", "REVIEW-0123abcdef"):
        with pytest.raises(SystemExit, match="review-id"):
            rec.check_freeze(bad, first, clean, man, repo)
    for bad_sha in (None, "abc", "g" * 40):
        with pytest.raises(SystemExit, match="reviewed-sha"):
            rec.check_freeze("REVIEW-0123abcd", bad_sha, clean, man, repo)
    with pytest.raises(SystemExit, match="not an ancestor"):
        rec.check_freeze("REVIEW-0123abcd", head, {**clean, "git_head": first}, man, repo)
    (repo / rec.SCRIPT_REL / "de9_record.py").write_text("# changed\n")
    with pytest.raises(SystemExit, match="differ between the reviewed commit and the working tree"):
        rec.check_freeze("REVIEW-0123abcd", first, clean, man, repo)
    (repo / rec.SCRIPT_REL / "de9_record.py").write_text("# de9_record.py\n")
    man.write_text(json.dumps({"reviews": {first: "REVIEW-0123abcd"}}))
    with pytest.raises(SystemExit, match="binds REVIEW-0123abcd"):
        rec.check_freeze("REVIEW-89abcdef", first, clean, man, repo)


def test_documents_must_come_from_the_evaluator_head_and_the_reviewed_code() -> None:
    from validation.scripts.elastic_attribution import de9_record as rec

    freeze = {
        "sha": "b" * 40,
        "head": "a" * 40,
        "scripts": {"de3_probes.py": "1" * 64},
        "ancestor": lambda s, h: True,
    }
    ok = {
        "code": {
            "git_head": "a" * 40,
            "git_clean": True,
            "scripts_sha256": {"de3_probes.py": "1" * 64, "de5_sampling.py": "9" * 64},
        }
    }
    rec.check_document_code(ok, "doc", freeze)
    for patch, msg in (
        ({"git_head": "c" * 40}, "differ from the evaluator HEAD"),
        ({"git_clean": False}, "clean flag"),
        ({"scripts_sha256": {"de3_probes.py": "2" * 64}}, "differs from the reviewed commit"),
    ):
        with pytest.raises(SystemExit, match=msg):
            rec.check_document_code({"code": {**ok["code"], **patch}}, "doc", freeze)
    with pytest.raises(SystemExit, match="does not contain the reviewed commit"):
        rec.check_document_code(ok, "doc", {**freeze, "ancestor": lambda s, h: False})
    b0doc = {
        "provenance": {"git_head": "a" * 40, "git_clean": True}
    }  # the projected B0 record carries its HEAD in provenance
    rec.check_document_code(b0doc, "de9-b0.json", freeze)
    with pytest.raises(SystemExit, match="differ from the evaluator HEAD"):
        rec.check_document_code(
            {"provenance": {"git_head": "d" * 40, "git_clean": True}}, "de9-b0.json", freeze
        )


def test_record_run_binds_the_review_writes_the_manifest_and_never_overwrites(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from validation.scripts.elastic_attribution import de9_record as rec

    inp = _write_inputs(tmp_path / "in" if (tmp_path / "in").mkdir() is None else tmp_path)
    for p in [
        *inp["de9"].glob("de3-*.json"),
        *[
            q
            for q in inp["de3"].glob("de3-*.json")
            if not any(q.name.endswith(f"-{v}.json") for v in REUSED_V)
        ],
    ]:  # the harness documents record the evaluator HEAD and the reviewed code
        d = json.loads(p.read_text())
        d["code"] = {
            "git_head": "a" * 40,
            "git_clean": True,
            "scripts_sha256": {"de3_probes.py": "1" * 64},
        }
        p.write_text(json.dumps(d))
    monkeypatch.setattr(
        c, "repo_state", lambda repo=None: {"git_head": "a" * 40, "git_clean": True}
    )
    monkeypatch.setattr(rec, "ARCHIVE_JSON", inp["archive"])
    monkeypatch.setattr(rec, "hashes_at", lambda sha, repo=None: {"de3_probes.py": "1" * 64})
    monkeypatch.setattr(rec, "working_hashes", lambda repo=None: {"de3_probes.py": "1" * 64})
    monkeypatch.setattr(rec, "is_ancestor", lambda s, h, repo=None: True)
    monkeypatch.setattr(rec, "DE9_CODE_FILES", ("de3_probes.py",))
    monkeypatch.setattr(rec, "HASH_SOURCES", (inp["hashes"],))
    exp = tmp_path / "expected.json"
    exp.write_text(json.dumps(_expected()))
    out = tmp_path / "out"
    doc = rec.run(
        inp["de3"],
        inp["de9"],
        out,
        inp["b0"],
        "REVIEW-0123abcd",
        "b" * 40,
        "REVIEW-0123abcd0123456789abcdef01234567",
        exp,
    )
    man = json.loads((out / "de9-manifest.json").read_text())
    assert (
        man["reviews"] == {"b" * 40: "REVIEW-0123abcd"}
        and man["review_id_full_text"].startswith("REVIEW-0123abcd0123")
        and man["reviewed_sha"] == "b" * 40
    )
    assert doc["provenance"]["reviewed_sha"] == "b" * 40 and (out / "git" / "de9.json").is_file()
    before = {str(p): p.read_bytes() for p in sorted(out.rglob("*")) if p.is_file()}
    with pytest.raises(SystemExit, match="already published"):
        rec.run(inp["de3"], inp["de9"], out, inp["b0"], "REVIEW-0123abcd", "b" * 40, None, exp)
    assert {str(p): p.read_bytes() for p in sorted(out.rglob("*")) if p.is_file()} == before
    out2 = tmp_path / "out2"
    p = inp["de9"] / "de3-200-pinned-r0.json"
    d = json.loads(p.read_text())
    d["code"]["git_head"] = "e" * 40  # a document from another commit
    p.write_text(json.dumps(d))
    with pytest.raises(SystemExit, match="differ from the evaluator HEAD"):
        rec.run(inp["de3"], inp["de9"], out2, inp["b0"], "REVIEW-0123abcd", "b" * 40, None, exp)


def test_seed_override_guard_mirrors_the_registered_de9_seeds() -> None:
    import importlib.util

    from validation.scripts.elastic_attribution import de3_probes as d3

    spec = importlib.util.spec_from_file_location(
        "seed_families", Path(__file__).with_name("test_v3_005d_seed_families.py")
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert d3.REGISTERED_20992[200] == frozenset(
        mod.DE9_A1_200 | {20992101}
    ) and d3.REGISTERED_20992[150] == frozenset(mod.DE9_A1_150 | {20992100})
    assert d3.REGISTERED_20992[200] | d3.REGISTERED_20992[150] == mod.DE9_USED and mod.DE9_B2 == {
        20992100,
        20992101,
    }
    from validation.scripts.elastic_attribution import de9_record as rec

    assert {rec.A1_SEED_BASE[200] + k for k in range(rec.A1_K[200])} == mod.DE9_A1_200 and {
        rec.A1_SEED_BASE[150] + k for k in range(rec.A1_K[150])
    } == mod.DE9_A1_150
    assert set(rec.B2_SEEDS.values()) == mod.DE9_B2


def test_b0_provenance_discloses_the_single_pre_run_look() -> None:
    assert (
        "100 MeV node only" in b0.PRE_RUN_LOOKS
        and "2e-12" in b0.PRE_RUN_LOOKS
        and "no other node" in b0.PRE_RUN_LOOKS
        and "dip test" in b0.PRE_RUN_LOOKS
    )
    src = Path(b0.__file__).read_text()
    assert (
        '"pre_run_looks": PRE_RUN_LOOKS' in src and "**c.repo_state()" in src
    )  # fixed provenance fields of the B0 record


# ---------------------------------------------------------------------------------------------------------------------------------
# review round REVIEW-a701b2c8: once-only variants, fail-closed inputs, incomplete records, carry-up check
# ---------------------------------------------------------------------------------------------------------------------------------
def test_every_investigation_variant_is_once_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from validation.scripts.elastic_attribution import de3_probes as d3

    out, att = tmp_path / "de3", tmp_path / "attempts.json"
    out.mkdir()
    ran: list[str] = []

    def fake(
        variant: str, energy: int, n: int, batches: int, lateral: bool, out: Path, **kw: Any
    ) -> dict[str, Any]:
        ran.append(variant)
        (out / d3.run_name(energy, variant)).write_text("{}")
        return {"valid": True}

    monkeypatch.setattr(d3, "_run_one", fake)
    for v in (
        "pp-level-0.90-ge150",
        "pp-level-0.95-ge150",
        "pp-level-monotone",
        "pA-s0.5",
        "m2prime",
        "pp-only",
    ):
        d3.run_one(v, 150, 100000, 20, False, out, att)
        with pytest.raises(SystemExit, match="exists"):
            d3.run_one(
                v, 150, 100000, 20, False, out, att
            )  # the existing target is never overwritten
    assert len(ran) == 6
    log = json.loads(att.read_text())
    assert [r["outcome"] for r in log if "attempt" in r] == ["published-valid"] * 6 and sum(
        r["outcome"] == "started" for r in log
    ) == 6
    n = len(log)
    for v in d3.UNGUARDED_VARIANTS:  # the pinned table, its base copy and the EM-only run keep the plain path (no attempt record)
        d3.run_one(v, 200, 100000, 20, False, out, att)
    assert len(json.loads(att.read_text())) == n and set(d3.UNGUARDED_VARIANTS) == {
        "pinned",
        "all",
        "emonly",
    }


def _freeze() -> dict[str, Any]:
    return {
        "sha": "b" * 40,
        "head": "a" * 40,
        "scripts": {"de3_probes.py": "1" * 64},
        "ancestor": lambda s, h: True,
    }


def _code_everywhere(inp: dict[str, Path]) -> None:
    for p in [
        *inp["de9"].glob("de3-*.json"),
        *[
            q
            for q in inp["de3"].glob("de3-*.json")
            if not any(q.name.endswith(f"-{v}.json") for v in REUSED_V)
        ],
    ]:
        d = json.loads(p.read_text())
        d["code"] = {
            "git_head": "a" * 40,
            "git_clean": True,
            "scripts_sha256": {"de3_probes.py": "1" * 64},
        }
        p.write_text(json.dumps(d))
    b = json.loads(inp["b0"].read_text())
    inp["b0"].write_text(json.dumps(b))


def _patch(path: Path, **kv: Any) -> None:
    d = json.loads(path.read_text())
    d.update(kv)
    path.write_text(json.dumps(d))


def test_every_document_is_validated_fail_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from validation.scripts.elastic_attribution import de9_record as rec

    inp = _write_inputs(tmp_path)
    _code_everywhere(inp)
    assert (
        _build(inp, monkeypatch, freeze=_freeze())["complete"] is True
    )  # the pristine set passes with the code binding as well
    b3 = inp["de3"] / "de3-200-pp-level-monotone.json"
    orig = b3.read_text()
    cases = [({"energy": 150}, "energy"), ({"variant": "all"}, "variant"), ({"n": 50000}, "n"), ({"n_batches": 10}, "n_batches"), ({"mode": "emonly"}, "mode"),
             ({"backend": "warp-cuda"}, "backend"), ({"precision": "float32"}, "precision"), ({"seed": 20990100}, "seed"), ({"elastic_table_id": "x"}, "elastic table id"),
             ({"elastic_table_npz_sha256": "x"}, "npz sha256"), ({"elastic_table_npz_sha256": None}, "npz sha256")]  # fmt: skip
    for patch, msg in cases:
        b3.write_text(orig)
        _patch(b3, **patch)
        with pytest.raises(SystemExit, match=msg):
            _build(inp, monkeypatch)
    b3.write_text(orig)
    d = json.loads(orig)
    del d["elastic_table_npz_sha256"]  # a new document must carry the npz hash
    b3.write_text(json.dumps(d))
    with pytest.raises(SystemExit, match="npz sha256"):
        _build(inp, monkeypatch)
    b3.write_text(orig)
    d = json.loads(orig)
    d["code"]["git_head"] = "e" * 40  # B3 documents go through the code binding
    b3.write_text(json.dumps(d))
    with pytest.raises(SystemExit, match="differ from the evaluator HEAD"):
        _build(inp, monkeypatch, freeze=_freeze())
    b3.write_text(orig)
    with pytest.raises(SystemExit, match="no expected table"):
        _build(
            inp,
            monkeypatch,
            expected={k: v for k, v in _expected().items() if k != "pp-level-monotone"},
        )
    # reused D-E3 documents: bound to the committed hashes and checked for identity
    base = inp["de3"] / "de3-150-all.json"
    keep = base.read_text()
    _patch(base, note="edited")  # different bytes -> hash differs from the committed list
    with pytest.raises(SystemExit, match="sha256 differs"):
        _build(inp, monkeypatch)
    base.write_text(keep)
    _patch(base, elastic_table_npz_sha256="wrong")
    inp["hashes"].write_text(
        inp["hashes"]
        .read_text()
        .replace(
            hashlib.sha256(keep.encode()).hexdigest(), hashlib.sha256(base.read_bytes()).hexdigest()
        )
    )
    with pytest.raises(SystemExit, match="npz sha256"):
        _build(inp, monkeypatch)
    assert rec.committed_hashes((inp["hashes"],))["de3-150-all.json"]
    two = tmp_path / "two.md"
    two.write_text("- `a/de3-150-all.json`  " + "0" * 64 + "\n")
    with pytest.raises(SystemExit, match="different sha256"):
        rec.committed_hashes((inp["hashes"], two))


def test_a_published_invalid_run_is_recorded_and_never_synthesised_as_valid(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from validation.scripts.elastic_attribution import de9_record as rec

    inp = _write_inputs(tmp_path)
    _patch(inp["de9"] / "de3-200-pinned-r2.json", valid=False)
    doc = _build(inp, monkeypatch)
    assert (
        doc["complete"] is False
        and doc["invalid_documents"] == {"de3-200-pinned-r2.json": "published-invalid"}
        and "b1" not in doc
        and "a1" not in doc
        and "a1_sensitivity" not in doc
    )
    md = rec.render_report(
        c.git_projection(json.loads(json.dumps(doc, default=float)), rec.GIT_ALLOW)
    )
    assert "INCOMPLETE — not the D-E9 synthesis" in md and "published-invalid" in md
    inp2 = _write_inputs(tmp_path / "r" if (tmp_path / "r").mkdir() is None else tmp_path)
    _patch(
        inp2["de9"] / "de3-150-emonly-r0.json", relative_residual=1e-3
    )  # the energy balance does not close
    assert _build(inp2, monkeypatch)["invalid_documents"]


def test_missing_m2prime_or_b0_is_required_unless_incomplete_is_allowed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from validation.scripts.elastic_attribution import de9_record as rec

    inp = _write_inputs(tmp_path)
    (inp["de3"] / "de3-200-m2prime.json").unlink()
    with pytest.raises(SystemExit, match="de3-200-m2prime.json is missing"):
        _build(inp, monkeypatch)
    doc = _build(inp, monkeypatch, allow_incomplete=True)
    assert (
        doc["complete"] is False
        and "a1_sensitivity" not in doc
        and "b1" in doc
        and any("m2prime" in r for r in doc["incomplete_reasons"])
    )
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), rec.GIT_ALLOW)
    assert "INCOMPLETE — not the D-E9 synthesis" in rec.render_report(
        proj
    ) and "D-E9: exploratory" in rec.render_report(proj)
    inp2 = _write_inputs(tmp_path / "b" if (tmp_path / "b").mkdir() is None else tmp_path)
    with pytest.raises(SystemExit, match="B0 result is a required input"):
        _build(inp2, monkeypatch, b0=None)
    with pytest.raises(SystemExit, match="B0 result is a required input"):
        _build(inp2, monkeypatch, b0=tmp_path / "absent.json")
    d2 = _build(inp2, monkeypatch, b0=None, allow_incomplete=True)
    assert d2["complete"] is False and "b0" not in d2 and "a1_sensitivity" in d2


def test_incomplete_run_is_labelled_in_the_header_and_the_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from validation.scripts.elastic_attribution import de9_record as rec

    inp = _write_inputs(tmp_path / "in" if (tmp_path / "in").mkdir() is None else tmp_path)
    (inp["de3"] / "de3-200-m2prime.json").unlink()
    _code_everywhere(inp)
    monkeypatch.setattr(
        c, "repo_state", lambda repo=None: {"git_head": "a" * 40, "git_clean": True}
    )
    monkeypatch.setattr(rec, "ARCHIVE_JSON", inp["archive"])
    monkeypatch.setattr(rec, "hashes_at", lambda sha, repo=None: {"de3_probes.py": "1" * 64})
    monkeypatch.setattr(rec, "working_hashes", lambda repo=None: {"de3_probes.py": "1" * 64})
    monkeypatch.setattr(rec, "is_ancestor", lambda s, h, repo=None: True)
    monkeypatch.setattr(rec, "DE9_CODE_FILES", ("de3_probes.py",))
    monkeypatch.setattr(rec, "HASH_SOURCES", (inp["hashes"],))
    exp = tmp_path / "expected.json"
    exp.write_text(json.dumps(_expected()))
    with pytest.raises(SystemExit, match="is missing"):
        rec.run(
            inp["de3"],
            inp["de9"],
            tmp_path / "o1",
            inp["b0"],
            "REVIEW-0123abcd",
            "b" * 40,
            None,
            exp,
        )
    rec.run(
        inp["de3"],
        inp["de9"],
        tmp_path / "o2",
        inp["b0"],
        "REVIEW-0123abcd",
        "b" * 40,
        None,
        exp,
        allow_incomplete=True,
    )
    assert json.loads((tmp_path / "o2" / "de9-manifest.json").read_text())["complete"] is False
    md = (tmp_path / "o2" / "git" / "de9-report.md").read_text()
    assert "INCOMPLETE — not the D-E9 synthesis" in md.split(c.BODY_MARKER)[0]


def test_b0_carry_up_check_reconstructs_the_level_above_150_mev(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from ionmc.nuclear import bgg

    monkeypatch.setattr(b0, "N_GAUSS", 64)
    sec = _sec([4.0e-3, 4.0e-3, 4.0e-3], energies=[1.0e7, 1.0e8, 1.5e8])
    grid = np.array([100.0, 150.0, 175.0, 200.0, 225.0, 250.0])
    rec150 = b0.sigma_ni_barn(sec, 150.0)
    table = np.array(
        [
            b0.sigma_ni_barn(sec, 100.0),
            rec150,
            *[rec150 * bgg.pp_bgg_mb(float(e)) / bgg.pp_bgg_mb(150.0) for e in grid[2:]],
        ]
    )
    res = b0.carry_up_check(sec, grid, table)
    absolutes = res.pop("_absolutes")
    assert (
        set(res["named_nodes"]) == {"175", "200", "250"}
        and res["match"]["all_within"]
        and "100" not in res["nodes"]
        and "150" not in res["nodes"]
    )
    g = res["grid_150_200_over_150"]
    assert (
        len(g) == 51
        and g["150"] == 1.0
        and all(g[f"{e}"] >= g[f"{e + 1}"] for e in range(150, 200))
    )
    bad = table.copy()
    bad[3] *= 1.01
    cu = b0.carry_up_check(sec, grid, bad)
    assert (
        not cu["match"]["all_within"]
        and not cu["named_nodes"]["200"]["within"]
        and cu["named_nodes"]["175"]["within"]
    )
    proj = c.git_projection(
        json.loads(json.dumps({"banner": c.BANNER, "carry_up_check": res}, default=float)),
        b0.GIT_ALLOW,
    )
    assert b0.leak_violations(proj, absolutes) == [] and "carried_over_table" in json.dumps(proj)
    assert "Carry-up check above 150 MeV" in b0.render_report({"disclosure": "d", "match": {"tolerance": 1e-3, "all_within": True, "max_abs_rel_diff": 0.0},
                                                               "dip": {"margin": 0.01, "persists": True, "min_energy_mev": 130.0, "min_over_lower_endpoint_minus_1": -0.1, "min_over_upper_endpoint_minus_1": -0.1},
                                                               "reading": "r", "named_nodes": {}, "levels_over_125mev": {}, "carry_up_check": {**res}})  # fmt: skip
