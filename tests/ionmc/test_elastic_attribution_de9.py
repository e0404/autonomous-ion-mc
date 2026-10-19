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


# ---------------------------------------------------------------------------------------------------------------------------------
# de9 record: B1, B3, A1, B2, Lead A table; fail-closed inputs (REVIEW-a701b2c8, REVIEW-6b2120c4)
# ---------------------------------------------------------------------------------------------------------------------------------
from validation.scripts.elastic_attribution import de9_record as rec  # noqa: E402

REAL_GEOMETRY = rec.expected_geometry
NZ = 400
PL_EM = 6.0
REUSED_V = ("all", "pp-only", "pA-only", "emonly")
B3_V = ("pp-level-0.90-ge150", "pp-level-0.95-ge150", "pp-level-monotone")
TABLE_V = ("all", "pp-only", "pA-only", "pinned", "m2prime", "m2", "rb", *B3_V)
GAINS = {
    "all": 0.025,
    "pp-only": 0.020,
    "pA-only": 0.005,
    "m2prime": 0.0213,
    "pp-level-0.90-ge150": 0.0225,
    "pp-level-0.95-ge150": 0.0238,
    "pp-level-monotone": 0.0242,
}
MP = pytest.MonkeyPatch()
CTX: dict[str, Any] = {}


@pytest.fixture(autouse=True)
def _undo_patches() -> Any:
    yield
    MP.undo()
    CTX.clear()


def _geometry(energy: int) -> dict[str, Any]:
    return {
        "origin_mm": [-200.0, -200.0, 0.0],
        "spacing_mm": [400.0, 400.0, 0.5],
        "shape": [1, 1, NZ],
        "name": "dose",
        "nz": NZ,
        "bin_mm": 0.5,
    }


def _curves(plateau: float, peak_over_plateau: float, jitter: float = 0.0) -> list[list[float]]:
    rows = []
    for i in range(20):
        r = np.full(NZ, plateau * (1.0 + 4e-4 * math.cos(i + jitter)))
        r[300] = plateau * peak_over_plateau * (1.0 + 1e-3 * math.sin(2 * i + jitter))
        rows.append(r.tolist())
    return rows


def _sidecar_for(v: str, pin_tid: str, sha: str) -> dict[str, Any]:
    """A synthetic table sidecar that DEFINES its variant the way the real builders do (``exploratory`` block)."""
    from validation.scripts.elastic_attribution import make_variant_tables as mv

    side: dict[str, Any] = {"label": v, "npz_sha256": sha}
    if v in rec.BUILDER_IDENTITY:
        side["exploratory"] = {
            **rec.BUILDER_CONFIG[v],
            "change": "p + A shapes",
            "purpose": "diagnostic",
        }
    elif v != "pinned":
        side["exploratory"] = {
            "d_e3_variant": v,
            "parameters": mv._parameters(v),
            "source_table_id": pin_tid,
        }
    return side


def _install_table(scratch: Path, v: str, side: dict[str, Any], payload: bytes) -> dict[str, str]:
    from ionmc.nuclear.build import table_id

    side = {**side, "npz_sha256": hashlib.sha256(payload).hexdigest()}
    tid = table_id(side)
    side["table_id"] = tid
    (scratch / "derived" / f"elastic-proton-{tid}.npz").write_bytes(payload)
    (scratch / "derived" / f"elastic-proton-{tid}.json").write_text(json.dumps(side))
    return {"table_id": tid, "npz_sha256": side["npz_sha256"]}


def _scratch_cache(tmp: Path) -> tuple[Path, dict[str, dict[str, str]]]:
    """A synthetic scratch cache (variants.json, builders, derived files with consistent ids) and the matching expected-tables map."""
    scratch = tmp / "scratch"
    (scratch / "derived").mkdir(parents=True)
    exp: dict[str, dict[str, str]] = {
        "pinned": _install_table(scratch, "pinned", _sidecar_for("pinned", "", ""), b"npz-pinned")
    }
    pin = dict(exp["pinned"])
    MP.setattr(c, "TABLE_ID", pin["table_id"])
    variants: dict[str, Any] = {}
    builders: dict[str, Any] = {}
    for v in TABLE_V:
        if v == "pinned":
            continue
        exp[v] = _install_table(
            scratch, v, _sidecar_for(v, pin["table_id"], ""), f"npz-{v}".encode()
        )
        (builders if v in rec.BUILDER_IDENTITY else variants)[v] = {"table_id": exp[v]["table_id"]}
    (scratch / "variants.json").write_text(json.dumps({"variants": variants, "builders": builders}))
    (tmp / "pin.json").write_text(json.dumps(pin))
    MP.setattr(rec, "PIN_JSON", tmp / "pin.json")
    MP.setattr(rec, "expected_geometry", _geometry)
    return scratch, exp


def _doc(
    energy: int,
    variant: str,
    plateau: float,
    seed: int,
    exp: dict[str, dict[str, str]],
    **extra: Any,
) -> dict[str, Any]:
    """A producer-shaped harness document (the fields ``de3_probes.run_one`` writes), synthetic numbers."""
    new = (
        "realisation" in extra or variant not in REUSED_V
    )  # the reused D-E3 documents predate the marker, the npz field and the code provenance
    emonly = variant == "emonly"
    d: dict[str, Any] = {
        "energy": energy, "variant": variant, "seed": seed, "idd_batches": _curves(plateau, 5.0, jitter=seed % 7), "valid": True, "relative_residual": 1e-16,
        "n": 100000, "n_batches": 20, "mode": "emonly" if emonly else "emel", "backend": "warp-cpu", "precision": "float64", "bin_mm": 0.5, "unit": rec.UNIT,
        "scoring": {"origin_mm": [-200.0, -200.0, 0.0], "spacing_mm": [400.0, 400.0, 0.5], "shape": [1, 1, NZ], "name": "dose"},
        "tallies": [] if emonly else ["recoil"], "elastic_table_id": exp["pinned" if emonly else variant]["table_id"],
    }  # fmt: skip
    if not emonly:
        d["recoil_idd_batches"] = _curves(0.002, 1.0)
    if new:
        d["elastic"] = not emonly
        d["code"] = {
            "git_head": "a" * 40,
            "git_clean": True,
            "scripts_sha256": {"de3_probes.py": "1" * 64},
        }
        if not emonly:
            d["elastic_table_npz_sha256"] = exp[variant]["npz_sha256"]
    d.update(extra)
    if "genealogy_split" in extra:
        d["tallies"] = ["recoil", "pp_secondary_proton", "primary_proton"]
    return d


def _is_new(p: Path) -> bool:
    return (
        p.name == "de9-b0.json"
        or p.parent.name == "de9"
        or any(p.name.endswith(f"-{v}.json") for v in (*B3_V, "m2prime"))
    )


def _relog(inp: dict[str, Path]) -> None:
    """Rebuild the attempt log, the committed hash list and the D-E8 record from the documents on disk (the state a correct producer run leaves)."""
    log: list[dict[str, Any]] = []
    m2p: dict[str, Any] = {}
    for p in sorted([*inp["de9"].glob("de3-*.json"), *inp["de3"].glob("de3-*.json"), inp["b0"]]):
        if not _is_new(p):
            continue
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
        if p.name.endswith("-m2prime.json"):
            m2p[p.name.split("-")[1]] = {"sha256": sha}
            continue
        n = len(log)
        log += [
            {"id": n, "outcome": "started", "output": str(p)},
            {"id": n + 1, "attempt": n, "outcome": "published-valid", "document_sha256": sha},
        ]
    inp["attempts"].write_text(json.dumps(log))
    inp["de8"].write_text(json.dumps({"runs": {"arms": {"m2prime": m2p}}}))
    inp["hashes"].write_text(
        "Input files (sha256):\n\n"
        + "".join(
            f"- `x/{p.name}`  {hashlib.sha256(p.read_bytes()).hexdigest()}\n"
            for v in REUSED_V
            for p in sorted(inp["de3"].glob(f"de3-*-{v}.json"))
        )
    )


def _put(path: Path, **kv: Any) -> None:
    """Edit a document and keep every authentication source consistent (so that the validator, not the authentication, is what refuses)."""
    d = json.loads(path.read_text())
    d.update(kv)
    path.write_text(json.dumps(d))
    _relog(CTX["inp"])


_B0_CACHE: dict[str, Any] = {}


def _synthetic_b0(tmp: Path, exp: dict[str, dict[str, str]], **override: Any) -> Path:
    """The projected B0 record as ``de9_pp_level.assemble`` produces it (synthetic LA150-like section and table values; built once per process)."""
    MP.setattr(c, "repo_state", lambda repo=None: {"git_head": "a" * 40, "git_clean": True})
    if "proj" not in _B0_CACHE:
        MP.setattr(b0, "N_GAUSS", 64)
        sec = _sec([4.0e-3, 1.0e-3, 4.0e-3], energies=[1.0e7, 1.25e8, 1.5e8])
        grid = np.arange(90.0, 251.0, 10.0)
        table = np.array([b0.sigma_ni_barn(sec, float(e)) if e <= 150.0 else 1.0 for e in grid])
        ident = {"source_sha256": rec.LA150_ZIP_SHA256, "table_id": "x", "npz_sha256": "y"}
        doc, _ = b0.assemble(sec, grid, table, ident, rec.LA150_H1_MEMBER_SHA256)
        _B0_CACHE["proj"] = c.git_projection(
            json.loads(json.dumps(doc, default=float)), b0.GIT_ALLOW
        )
    proj = json.loads(json.dumps(_B0_CACHE["proj"]))
    proj["provenance"].update(
        {
            "table_id": exp["pinned"]["table_id"],
            "npz_sha256": exp["pinned"]["npz_sha256"],
            "code_sha256": c.sha256_file(Path(b0.__file__)),
        }
    )
    proj["provenance"].update(override.pop("provenance", {}))
    proj.update(override)
    p = tmp / "de9-b0.json"
    p.write_text(json.dumps(proj))
    return p


def _write_inputs(
    tmp: Path, a1_gains: dict[int, list[float]] | None = None, with_b2: bool = True
) -> dict[str, Any]:
    scratch, exp = _scratch_cache(tmp)
    de3, de9 = tmp / "de3", tmp / "de9"
    de3.mkdir()
    de9.mkdir()
    seeds = {150: 20990100, 200: 20990101}
    for e in (150, 200):
        (de3 / f"de3-{e}-emonly.json").write_text(
            json.dumps(_doc(e, "emonly", PL_EM, seeds[e], exp))
        )
        for v, g in GAINS.items():
            (de3 / f"de3-{e}-{v}.json").write_text(
                json.dumps(_doc(e, v, PL_EM * (1 + g), seeds[e], exp))
            )
    rng = np.random.default_rng(1)
    for e, k in ((200, 8), (150, 4)):
        base = {200: 20992000, 150: 20992010}[e]
        gs = (a1_gains or {}).get(e) or list(0.0255 + 0.001 * rng.standard_normal(k))
        for r in range(k):
            for v, pl in (("pinned", PL_EM * (1 + gs[r])), ("emonly", PL_EM)):
                (de9 / f"de3-{e}-{v}-r{r}.json").write_text(
                    json.dumps(_doc(e, v, pl, base + r, exp, realisation=r))
                )
        if with_b2:
            d = _doc(e, "pinned", PL_EM * 1.025, {150: 20992100, 200: 20992101}[e], exp, realisation=100,
                     genealogy_split={"pp_secondary_proton": _curves(0.1, 1.0), "primary_proton": _curves(5.0, 1.0)})  # fmt: skip
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
    b0doc = _synthetic_b0(tmp, exp)
    de7 = tmp / "de7.json"
    de7.write_text(
        json.dumps(
            {"variants": {rec.DE7_NAME[v]: {"identity": dict(exp[v])} for v in rec.DE7_NAME}}
        )
    )
    inp: dict[str, Any] = {"de7": de7, "de3": de3, "de9": de9, "archive": tmp / "archive.json", "hashes": tmp / "hashes.md", "b0": b0doc, "expected": exp, "scratch": scratch,
                           "attempts": tmp / "attempts.json", "de8": tmp / "de8.json", "tmp": tmp}  # fmt: skip
    CTX["inp"] = inp
    _relog(inp)
    return inp


def _build(inp: dict[str, Any], **kw: Any) -> dict[str, Any]:
    MP.setattr(c, "repo_state", lambda repo=None: {"git_head": "a" * 40, "git_clean": True})
    return rec.build_document(
        inp["de3"], inp["de9"], kw["b0"] if "b0" in kw else inp["b0"], archive=inp["archive"], absolutes=kw.get("absolutes"), expected_tables=kw.get("expected", inp["expected"]),
        hash_sources=(inp["hashes"],), allow_incomplete=kw.get("allow_incomplete", False), freeze=kw.get("freeze"), scratch_cache=kw.get("scratch", inp["scratch"]),
        attempts_path=kw.get("attempts", inp["attempts"]), de8_record=inp["de8"], de7_record=inp["de7"],
    )  # fmt: skip


def _freeze() -> dict[str, Any]:
    return {
        "sha": "b" * 40,
        "head": "a" * 40,
        "scripts": {"de3_probes.py": "1" * 64},
        "ancestor": lambda s, h: True,
    }


def test_b1_b3_and_the_gain_definitions(tmp_path: Path) -> None:
    doc = _build(_write_inputs(tmp_path))
    r = doc["b1"]["200"]
    assert (
        r["gain_all"] == pytest.approx(0.025, abs=2e-4)
        and r["gain_pp_only"] == pytest.approx(0.020, abs=2e-4)
        and r["pp_share_of_sum"] == pytest.approx(0.8, abs=0.02)
    )
    assert r["additivity_gap"] == pytest.approx(0.0, abs=4e-4)
    b3 = doc["b3"]
    assert b3["pp-level-0.90-ge150"]["200"]["delta_gain"] == pytest.approx(
        0.0225 - 0.025, abs=2e-4
    ) and b3["pp-level-monotone"]["150"]["delta_gain"] == pytest.approx(0.0242 - 0.025, abs=2e-4)
    s = b3["slope"]["200"]
    assert s["dg_dln_s_0.90_to_0.95"] == pytest.approx(
        (0.0238 - 0.0225) / math.log(0.95 / 0.90), rel=0.05
    ) and s["dg_dln_s_0.95_to_1.00"] == pytest.approx((0.025 - 0.0238) / -math.log(0.95), rel=0.05)
    text = json.dumps(doc).lower()
    assert (
        "closing" not in text
        and "required" not in text
        and "residual" not in text
        and doc["complete"] is True
    )


def test_a1_reading_rule_is_applied_as_pre_registered(tmp_path: Path) -> None:
    g200 = [0.0245, 0.0250, 0.0255, 0.0260, 0.0265, 0.0258, 0.0252, 0.0248]
    doc = _build(_write_inputs(tmp_path, {200: g200, 150: [0.0200, 0.0210, 0.0205, 0.0195]}))
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
    assert (abs(a["z_archived"]) <= 2) == (
        a["archived_reading"] == "consistent with the realisation distribution"
    )


def test_lead_a_sensitivity_table_is_labelled_and_never_compared_with_a_band(
    tmp_path: Path,
) -> None:
    doc = _build(_write_inputs(tmp_path))
    t = doc["a1_sensitivity"]
    frozen = 0.02672217441583613 + (GAINS["m2prime"] - GAINS["all"]) - 0.018419142339446415
    assert (
        t["frozen"]["d"] == pytest.approx(frozen, abs=3e-4)
        and t["frozen"]["label"] == "the frozen value; the STOP stands"
    )
    label = f"sensitivity; frozen d = {t['frozen']['d']:.5f} stands"
    assert t["archived_minus_2se"]["d"] == pytest.approx(t["frozen"]["d"] - 2 * 0.00077, abs=1e-12)
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


def test_b2_share_and_the_dropped_case(tmp_path: Path) -> None:
    doc = _build(_write_inputs(tmp_path))
    r = doc["b2"]["200"]
    excess = PL_EM * 0.025
    assert (
        r["status"] == "recorded"
        and r["secondary_proton_share"] == pytest.approx(0.1 / excess, rel=1e-2)
        and r["recoil_share"] == pytest.approx(0.002 / excess, rel=1e-2)
    )
    dropped = _build(
        _write_inputs(
            tmp_path / "x" if (tmp_path / "x").mkdir() is None else tmp_path, with_b2=False
        )
    )
    assert dropped["b2"]["150"]["status"].startswith("dropped")


def test_projection_leak_check_and_report(tmp_path: Path) -> None:
    doc = _build(_write_inputs(tmp_path), absolutes=[123.456])
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
    inp = _write_inputs(tmp_path / "l" if (tmp_path / "l").mkdir() is None else tmp_path)
    with pytest.raises(SystemExit, match="leak check failed"):
        _build(inp, absolutes=[doc["b1"]["150"]["gain_pp_only"]])


def test_missing_documents_and_overwrites_are_refused(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path)
    (inp["de9"] / "de3-200-pinned-r3.json").unlink()
    _relog(inp)
    with pytest.raises(SystemExit, match="is missing"):
        _build(inp)
    inp2 = _write_inputs(tmp_path / "w" if (tmp_path / "w").mkdir() is None else tmp_path)
    MP.setattr(rec, "check_freeze", lambda *a, **k: {})
    MP.setattr(rec, "hashes_at", lambda sha, repo=None: {})
    (tmp_path / "exp.json").write_text(json.dumps(inp2["expected"]))
    out = tmp_path / "out"
    (out / "git").mkdir(parents=True)
    (out / "git" / "de9.json").write_text("{}")
    MP.setattr(c, "repo_state", lambda repo=None: {"git_head": "a" * 40, "git_clean": True})
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
            scratch_cache=inp2["scratch"],
        )


def test_every_document_is_validated_fail_closed(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path)
    assert (
        _build(inp, freeze=_freeze())["complete"] is True
    )  # the pristine set passes, with the code binding as well
    b3 = inp["de3"] / "de3-200-pp-level-monotone.json"
    orig = json.loads(b3.read_text())
    cases = [({"energy": 150}, "energy"), ({"variant": "all"}, "variant"), ({"n": 50000}, "n "), ({"n_batches": 10}, "n_batches"), ({"mode": "emonly"}, "mode"),
             ({"backend": "warp-cuda"}, "backend"), ({"precision": "float32"}, "precision"), ({"seed": 20990100}, "seed"), ({"elastic_table_id": "x"}, "elastic table id"),
             ({"elastic_table_npz_sha256": "x"}, "npz sha256"), ({"elastic_table_npz_sha256": None}, "npz sha256"), ({"elastic": False}, "elastic = true")]  # fmt: skip
    for patch, msg in cases:
        b3.write_text(json.dumps(orig))
        _put(b3, **patch)
        with pytest.raises(SystemExit, match=msg):
            _build(inp)
    d = dict(orig)
    del d["elastic_table_npz_sha256"]
    b3.write_text(json.dumps(d))
    _relog(inp)
    with pytest.raises(SystemExit, match="npz sha256"):
        _build(inp)
    for bad_code in (
        {"git_head": "e" * 40, "git_clean": True, "scripts_sha256": {"de3_probes.py": "1" * 64}},
        {"git_head": "a" * 40, "git_clean": True},
        {"git_head": "a" * 40, "git_clean": True, "scripts_sha256": {}},
        {"git_head": "a" * 40, "git_clean": True, "scripts_sha256": {"de5_sampling.py": "2" * 64}},
        {"git_head": "a" * 40, "git_clean": True, "scripts_sha256": {"de3_probes.py": "2" * 64}},
    ):  # B3 documents go through the code binding; hashes are REQUIRED
        b3.write_text(json.dumps({**orig, "code": bad_code}))
        _relog(inp)
        with pytest.raises(SystemExit, match="evaluator HEAD|scripts_sha256|sha256 of de3_probes"):
            _build(inp, freeze=_freeze())
    b3.write_text(json.dumps(orig))
    _relog(inp)
    with pytest.raises(SystemExit, match="no expected table"):
        _build(inp, expected={k: v for k, v in inp["expected"].items() if k != "pp-level-monotone"})
    base = (
        inp["de3"] / "de3-150-all.json"
    )  # reused D-E3 documents: bound to the committed hashes and checked for identity
    keep = json.loads(base.read_text())
    base.write_text(
        json.dumps({**keep, "note": "edited"})
    )  # different bytes: the committed hash no longer matches
    with pytest.raises(SystemExit, match="sha256 differs"):
        _build(inp)
    base.write_text(json.dumps({**keep, "elastic_table_npz_sha256": "wrong"}))
    _relog(inp)  # the list is made consistent: only the identity check can refuse now
    with pytest.raises(SystemExit, match="npz sha256"):
        _build(inp)
    two = tmp_path / "two.md"
    two.write_text("- `a/de3-150-all.json`  " + "0" * 64 + "\n")
    with pytest.raises(SystemExit, match="different sha256"):
        rec.committed_hashes((inp["hashes"], two))


def test_generated_documents_are_authenticated_against_the_attempt_log_and_the_de8_record(
    tmp_path: Path,
) -> None:
    inp = _write_inputs(tmp_path)
    p = inp["de9"] / "de3-200-pinned-r1.json"
    keep = p.read_text()
    p.write_text(
        json.dumps({**json.loads(keep), "note": "edited after publication"})
    )  # the log keeps the published sha
    with pytest.raises(SystemExit, match="attempt-log entry"):
        _build(inp)
    p.write_text(keep)
    log = json.loads(inp["attempts"].read_text())
    inp["attempts"].write_text(
        json.dumps(
            [
                r
                for r in log
                if r.get("output") != str(p) and log[r.get("attempt", -1)].get("output") != str(p)
            ]
        )
    )
    with pytest.raises(SystemExit, match="no \\(unique\\) published entry"):
        _build(inp)
    inp["attempts"].write_text(
        json.dumps(
            [{**r, "document_sha256": "0" * 64} if "document_sha256" in r else r for r in log]
        )
    )
    with pytest.raises(SystemExit, match="attempt-log entry"):
        _build(inp)
    inp["attempts"].write_text(
        json.dumps([{k: v for k, v in r.items() if k != "document_sha256"} for r in log])
    )  # a log without the document hash
    with pytest.raises(SystemExit, match="no \\(unique\\) published entry"):
        _build(inp)
    _relog(inp)
    inp["attempts"].unlink()
    with pytest.raises(SystemExit, match="attempt log is missing"):
        _build(inp)
    _relog(inp)
    m2p = inp["de3"] / "de3-200-m2prime.json"
    m2p.write_text(json.dumps({**json.loads(m2p.read_text()), "note": "x"}))
    with pytest.raises(SystemExit, match="committed D-E8 record"):
        _build(inp)


def test_emonly_contract_is_explicit(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path)
    p = inp["de9"] / "de3-150-emonly-r0.json"
    d = json.loads(p.read_text())
    assert (
        d["elastic"] is False and "elastic_table_npz_sha256" not in d
    )  # the producer contract: no npz hash for EM-only runs
    for k in ("elastic_table_id",):
        d2 = {kk: v for kk, v in d.items() if kk != k}
        p.write_text(json.dumps(d2))
        _relog(inp)
        assert _build(inp)["complete"] is True  # exempt from the table identity
    for marker in (True, None):
        p.write_text(
            json.dumps(
                {
                    **{k: v for k, v in d.items() if k != "elastic"},
                    **({"elastic": marker} if marker is not None else {}),
                }
            )
        )
        _relog(inp)
        with pytest.raises(SystemExit, match="elastic = false"):
            _build(inp)
    p.write_text(json.dumps(d))
    _relog(inp)
    reused = (
        inp["de3"] / "de3-150-emonly.json"
    )  # a reused (hash-bound, pre-marker) EM-only document needs no marker
    assert "elastic" not in json.loads(reused.read_text()) and _build(inp)["complete"] is True


def test_expected_tables_are_verified_against_the_scratch_cache(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path)
    assert _build(inp)["complete"] is True
    exp = json.loads(json.dumps(inp["expected"]))
    bad = json.loads(json.dumps(exp))
    bad["pp-level-monotone"]["npz_sha256"] = "0" * 64
    with pytest.raises(SystemExit, match="not consistent with the scratch cache"):
        _build(inp, expected=bad)
    bad = json.loads(json.dumps(exp))
    bad["pp-level-0.90-ge150"] = exp[
        "pp-level-0.95-ge150"
    ]  # id and npz of another variant: variant name <-> id mismatch
    with pytest.raises(SystemExit, match="variants.json with this table id"):
        _build(inp, expected=bad)
    bad = {**exp, "unknown-variant": exp["all"]}
    with pytest.raises(SystemExit, match="not in the scratch cache variants.json"):
        _build(inp, expected=bad)
    side_p = inp["scratch"] / "derived" / f"elastic-proton-{exp['all']['table_id']}.json"
    side = json.loads(side_p.read_text())
    side["variant"] = "tampered"  # the sidecar no longer reproduces its table id
    side_p.write_text(json.dumps(side))
    with pytest.raises(SystemExit, match="not consistent with the scratch cache"):
        _build(inp)
    side["variant"] = "all"
    side_p.write_text(json.dumps(side))
    (inp["tmp"] / "pin.json").write_text(
        json.dumps({"table_id": exp["pinned"]["table_id"], "npz_sha256": "1" * 64})
    )
    with pytest.raises(SystemExit, match="pinned"):
        _build(inp)
    with pytest.raises(SystemExit, match="scratch cache .* is required"):
        rec.build_document(
            inp["de3"],
            inp["de9"],
            inp["b0"],
            archive=inp["archive"],
            expected_tables=exp,
            hash_sources=(inp["hashes"],),
        )


def test_frozen_configuration_of_the_arrays_is_checked(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path)
    b3 = inp["de3"] / "de3-150-pp-level-monotone.json"
    orig = json.loads(b3.read_text())
    short = [row[:-1] for row in orig["idd_batches"]]
    nan_rows = json.loads(json.dumps(orig["idd_batches"]))
    nan_rows[3][5] = float("nan")
    cases = [({"scoring": {**orig["scoring"], "origin_mm": [-100.0, -100.0, 0.0]}}, "scoring origin_mm"), ({"scoring": {**orig["scoring"], "shape": [1, 1, 349]}}, "scoring shape"),
             ({"scoring": {**orig["scoring"], "spacing_mm": [400.0, 400.0, 1.0]}}, "scoring spacing_mm"), ({"bin_mm": 1.0}, "bin size"), ({"unit": "Gy"}, "bin size / unit"),
             ({"tallies": []}, "tally names"), ({"tallies": ["recoil", "extra"]}, "tally names"), ({"idd_batches": short}, "array idd_batches"),
             ({"idd_batches": orig["idd_batches"][:10]}, "array idd_batches"), ({"idd_batches": nan_rows}, "not finite"), ({"recoil_idd_batches": None}, "array recoil_idd_batches")]  # fmt: skip
    for patch, msg in cases:
        b3.write_text(json.dumps(orig))
        _put(b3, **patch)
        with pytest.raises(SystemExit, match=msg):
            _build(inp)
    b3.write_text(json.dumps(orig))
    _relog(inp)
    g = {
        **json.loads(b3.read_text()),
        "genealogy_split": {"pp_secondary_proton": _curves(0.1, 1.0)},
    }  # a split without its tally names and both arrays
    _put(b3, genealogy_split=g["genealogy_split"])
    with pytest.raises(SystemExit, match="tally names"):
        _build(inp)


def test_a_published_invalid_run_is_recorded_and_never_synthesised_as_valid(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path)
    _put(inp["de9"] / "de3-200-pinned-r2.json", valid=False)
    doc = _build(inp)
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
    _put(
        inp2["de9"] / "de3-150-emonly-r0.json", relative_residual=1e-3
    )  # the energy balance does not close
    assert _build(inp2)["invalid_documents"]


def test_missing_m2prime_or_b0_is_required_unless_incomplete_is_allowed(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path)
    (inp["de3"] / "de3-200-m2prime.json").unlink()
    _relog(inp)
    with pytest.raises(SystemExit, match="de3-200-m2prime.json is missing"):
        _build(inp)
    doc = _build(inp, allow_incomplete=True)
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
    for missing in (None, tmp_path / "absent.json"):
        with pytest.raises(SystemExit, match="B0 result is a required input"):
            _build(inp2, b0=missing)
    d2 = _build(inp2, b0=None, allow_incomplete=True)
    assert d2["complete"] is False and "b0" not in d2 and "a1_sensitivity" in d2


def _run_patched(inp: dict[str, Any]) -> None:
    MP.setattr(c, "repo_state", lambda repo=None: {"git_head": "a" * 40, "git_clean": True})
    MP.setattr(rec, "ARCHIVE_JSON", inp["archive"])
    MP.setattr(rec, "ATTEMPTS_DEFAULT", inp["attempts"])
    MP.setattr(rec, "DE8_RECORD", inp["de8"])
    MP.setattr(rec, "DE7_RECORD", inp["de7"])
    MP.setattr(rec, "hashes_at", lambda sha, repo=None: {"de3_probes.py": "1" * 64})
    MP.setattr(rec, "working_hashes", lambda repo=None: {"de3_probes.py": "1" * 64})
    MP.setattr(rec, "is_ancestor", lambda s, h, repo=None: True)
    MP.setattr(rec, "DE9_CODE_FILES", ("de3_probes.py",))
    MP.setattr(rec, "HASH_SOURCES", (inp["hashes"],))


def test_record_run_binds_the_review_writes_the_manifest_and_never_overwrites(
    tmp_path: Path,
) -> None:
    inp = _write_inputs(tmp_path / "in" if (tmp_path / "in").mkdir() is None else tmp_path)
    _run_patched(inp)
    exp = tmp_path / "expected.json"
    exp.write_text(json.dumps(inp["expected"]))
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
        scratch_cache=inp["scratch"],
    )
    man = json.loads((out / "de9-manifest.json").read_text())
    assert (
        man["reviews"] == {"b" * 40: "REVIEW-0123abcd"}
        and man["review_id_full_text"].startswith("REVIEW-0123abcd0123")
        and man["reviewed_sha"] == "b" * 40
        and man["complete"] is True
    )
    assert doc["provenance"]["reviewed_sha"] == "b" * 40 and (out / "git" / "de9.json").is_file()
    assert (
        doc["provenance"]["expected_tables_sha256"] == hashlib.sha256(exp.read_bytes()).hexdigest()
    )  # the map is recorded by its hash
    before = {str(p): p.read_bytes() for p in sorted(out.rglob("*")) if p.is_file()}
    with pytest.raises(SystemExit, match="already published"):
        rec.run(
            inp["de3"],
            inp["de9"],
            out,
            inp["b0"],
            "REVIEW-0123abcd",
            "b" * 40,
            None,
            exp,
            scratch_cache=inp["scratch"],
        )
    assert {str(p): p.read_bytes() for p in sorted(out.rglob("*")) if p.is_file()} == before
    p = inp["de9"] / "de3-200-pinned-r0.json"
    d = json.loads(p.read_text())
    d["code"]["git_head"] = "e" * 40  # a document from another commit
    _put(p, code=d["code"])
    with pytest.raises(SystemExit, match="differ from the evaluator HEAD"):
        rec.run(
            inp["de3"],
            inp["de9"],
            tmp_path / "out2",
            inp["b0"],
            "REVIEW-0123abcd",
            "b" * 40,
            None,
            exp,
            scratch_cache=inp["scratch"],
        )


def test_incomplete_run_is_labelled_in_the_header_and_the_manifest(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path / "in" if (tmp_path / "in").mkdir() is None else tmp_path)
    (inp["de3"] / "de3-200-m2prime.json").unlink()
    _relog(inp)
    _run_patched(inp)
    exp = tmp_path / "expected.json"
    exp.write_text(json.dumps(inp["expected"]))
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
            scratch_cache=inp["scratch"],
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
        scratch_cache=inp["scratch"],
    )
    assert json.loads((tmp_path / "o2" / "de9-manifest.json").read_text())["complete"] is False
    md = (tmp_path / "o2" / "git" / "de9-report.md").read_text()
    assert "INCOMPLETE — not the D-E9 synthesis" in md.split(c.BODY_MARKER)[0]


# ---------------------------------------------------------------------------------------------------------------------------------
# documents generated by the REAL producers (tiny settings): the producer and the validator contracts agree
# ---------------------------------------------------------------------------------------------------------------------------------
def test_validator_accepts_documents_generated_by_the_real_producer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from validation.scripts.elastic_attribution import de3_probes as d3

    pinned_npz = c.CACHE_DIR / "derived" / f"elastic-proton-{c.TABLE_ID}.npz"
    if not pinned_npz.is_file():
        pytest.skip("pinned elastic table absent from the data cache")
    state = {"git_head": "a" * 40, "git_clean": True}
    monkeypatch.setattr(c, "repo_state", lambda repo=None: dict(state))
    monkeypatch.setenv("IONMC_CACHE_DIR", str(c.CACHE_DIR))
    monkeypatch.setattr(rec, "N_HIST", 200)
    monkeypatch.setattr(
        rec, "N_BATCH", 2
    )  # tiny settings: the validator's frozen configuration is patched to match, nothing else
    out, att = tmp_path / "de9", tmp_path / "attempts.json"
    out.mkdir()
    pin = json.loads((c.REPO / "src" / "ionmc" / "data" / "elastic_table_pin.json").read_text())
    expected = {"pinned": {"table_id": pin["table_id"], "npz_sha256": pin["npz_sha256"]}}
    scripts = c.code_provenance()["scripts_sha256"]
    freeze = {
        "sha": "b" * 40,
        "head": "a" * 40,
        "scripts": {"de3_probes.py": scripts["de3_probes.py"]},
        "ancestor": lambda s, h: True,
    }
    for e, seed_a1, seed_b2 in ((150, 20992010, 20992100), (200, 20992000, 20992101)):
        em = d3.run_one("emonly", e, 200, 2, False, out, att, seed_override=seed_a1, realisation=0)
        el = d3.run_one("pinned", e, 200, 2, False, out, att, seed_override=seed_a1, realisation=0)
        b2 = d3.run_one(
            "pinned",
            e,
            200,
            2,
            False,
            out,
            att,
            seed_override=seed_b2,
            realisation=100,
            genealogy_split=True,
        )
        monkeypatch.setattr(rec, "expected_geometry", REAL_GEOMETRY)  # the real V5 geometry
        for doc, name, variant, seed in (
            (em, f"de3-{e}-emonly-r0.json", "emonly", seed_a1),
            (el, f"de3-{e}-pinned-r0.json", "pinned", seed_a1),
            (b2, f"de3-{e}-pinned-r100.json", "pinned", seed_b2),
        ):
            path = out / name
            assert (
                rec.validate_document(
                    doc,
                    name,
                    energy=e,
                    variant=variant,
                    seed=seed,
                    expected=expected,
                    freeze=freeze,
                )
                is True
            )
            assert (
                rec.attempt_log_sha(json.loads(att.read_text()), name)
                == hashlib.sha256(path.read_bytes()).hexdigest()
            )  # the log authenticates the published file
        assert (
            em["elastic"] is False
            and "elastic_table_npz_sha256" not in em
            and el["elastic"] is True
            and el["elastic_table_npz_sha256"] == pin["npz_sha256"]
        )
        assert b2["tallies"] == ["recoil", "pp_secondary_proton", "primary_proton"] and set(
            b2["genealogy_split"]
        ) == {"pp_secondary_proton", "primary_proton"}
    # the same documents are refused when one frozen property is violated
    bad = dict(el)
    bad["scoring"] = {**el["scoring"], "shape": [1, 1, 3]}
    with pytest.raises(SystemExit, match="scoring shape"):
        rec.validate_document(
            bad, "x", energy=200, variant="pinned", seed=20992000, expected=expected, freeze=freeze
        )


def test_b0_quadrature_convergence_block_and_the_corrected_docstring() -> None:
    sec = _sec([4.0e-3, 2.0e-3, 4.0e-3], energies=[1.0e7, 1.25e8, 1.5e8])
    q = b0.quadrature_convergence(sec, nodes=(100.0, 125.0, 150.0))
    absolutes = q.pop("_absolutes")
    assert set(q["nodes"]) == {"100", "125", "150"} and all(
        set(v) == {"gl1024_vs_gl2048", "simpson4001_vs_gl2048", "simpson20001_vs_gl2048"}
        for v in q["nodes"].values()
    )
    assert (
        q["acceptance"] == 1e-6
        and q["accepted"] is True
        and q["max_abs_rel_diff"] < 1e-6
        and len(absolutes) == 12
    )
    assert b0.simpson_uniform(np.linspace(0, 1, 5) ** 2, np.linspace(0, 1, 5)) == pytest.approx(
        1 / 3, rel=1e-12
    )  # exact for a quadratic
    doc = b0.__doc__ or ""
    assert (
        "4001-point uniform Simpson grid" in doc
        and "20001 points" not in doc
        and "quadrature_convergence" in doc
    )
    proj = c.git_projection(
        json.loads(json.dumps({"banner": c.BANNER, "quadrature_convergence": q}, default=float)),
        b0.GIT_ALLOW,
    )
    assert "gl1024_vs_gl2048" in json.dumps(proj) and b0.leak_violations(proj, absolutes) == []


def _retable(inp: dict[str, Any], v: str, mutate: Any) -> dict[str, dict[str, str]]:
    """Rewrite variant ``v``'s sidecar consistently (new id, npz hash, variants.json, expected map): only the SEMANTIC check can refuse."""
    old = inp["expected"][v]
    side = json.loads(
        (inp["scratch"] / "derived" / f"elastic-proton-{old['table_id']}.json").read_text()
    )
    mutate(side)
    new = _install_table(inp["scratch"], v, side, f"npz-{v}".encode())
    vj = json.loads((inp["scratch"] / "variants.json").read_text())
    for blk in ("variants", "builders"):
        if v in vj[blk]:
            vj[blk][v] = {"table_id": new["table_id"]}
    (inp["scratch"] / "variants.json").write_text(json.dumps(vj))
    if (
        v in rec.DE7_NAME
    ):  # keep the committed D-E7 record consistent: only the semantic check can refuse
        d7 = json.loads(inp["de7"].read_text())
        d7["variants"][rec.DE7_NAME[v]]["identity"] = dict(new)
        inp["de7"].write_text(json.dumps(d7))
    return {**inp["expected"], v: new}


def test_variant_identity_is_authenticated_from_the_sidecar_not_from_variants_json(
    tmp_path: Path,
) -> None:
    inp = _write_inputs(tmp_path)
    exp = inp["expected"]
    # relabel: the pp-only table presented as `all`, with a consistent variants.json (the hint) -- the sidecar says pp-only
    vj = json.loads((inp["scratch"] / "variants.json").read_text())
    vj["variants"]["all"] = {"table_id": exp["pp-only"]["table_id"]}
    (inp["scratch"] / "variants.json").write_text(json.dumps(vj))
    with pytest.raises(SystemExit, match="does not define this variant"):
        _build(inp, expected={**exp, "all": exp["pp-only"]})
    # relabel: an rb table presented as m2prime
    inp2 = _write_inputs(tmp_path / "b" if (tmp_path / "b").mkdir() is None else tmp_path)
    vj = json.loads((inp2["scratch"] / "variants.json").read_text())
    vj["builders"]["m2prime"] = {"table_id": inp2["expected"]["rb"]["table_id"]}
    (inp2["scratch"] / "variants.json").write_text(json.dumps(vj))
    d7 = json.loads(inp2["de7"].read_text())
    d7["variants"]["M2prime"]["identity"] = dict(
        inp2["expected"]["rb"]
    )  # even a D-E7 record relabelled the same way cannot make the sidecar define m2prime
    inp2["de7"].write_text(json.dumps(d7))
    with pytest.raises(SystemExit, match="does not define the m2prime construction"):
        _build(inp2, expected={**inp2["expected"], "m2prime": inp2["expected"]["rb"]})
    # each defining field, with everything else consistent
    cases = [("pp-level-0.90-ge150", lambda s: s["exploratory"].update(d_e3_variant="pp-level-0.95-ge150"), "does not define this variant"),
             ("pp-level-0.90-ge150", lambda s: s["exploratory"].update(parameters={"sigma_pp_scale_at_nodes_ge_150_mev": 0.95, "shape": "unchanged"}), "does not define this variant"),
             ("pp-level-monotone", lambda s: s["exploratory"].update(source_table_id="0" * 64), "does not define this variant"),
             ("pp-only", lambda s: s.pop("exploratory"), "does not define this variant"),
             ("m2prime", lambda s: s["exploratory"].update(pauli_variant="bdc-pn"), "does not define the m2prime construction"),
             ("m2prime", lambda s: s["exploratory"].update(pauli_emin_mev=5.0), "does not define the m2prime construction"),
             ("m2prime", lambda s: s["exploratory"].update(model="m2-optical-glauber"), "does not define the m2prime construction"),
             ("rb", lambda s: s["exploratory"].update(pauli_variant="clementel-villi"), "does not define the rb construction"),
             ("m2prime", lambda s: s["exploratory"].update(beta="horiuchi-table"), "does not define the m2prime construction"),
             ("m2prime", lambda s: s["exploratory"].update(pp_level="said-anchored-ge150"), "does not define the m2prime construction"),
             ("m2prime", lambda s: s["exploratory"].update(unfold="formfactor"), "does not define the m2prime construction"),
             ("m2prime", lambda s: s["exploratory"].update(inputs_unpinned=True), "does not define the m2prime construction"),
             ("m2prime", lambda s: s["exploratory"].update(variant="V-b"), "does not define the m2prime construction"),
             ("m2prime", lambda s: s["exploratory"].update(pauli_emin_mev=80), "does not define the m2prime construction"),
             ("m2prime", lambda s: s["exploratory"].pop("beta"), "does not define the m2prime construction"),
             ("m2", lambda s: s["exploratory"].update(pauli_variant="clementel-villi", pauli_emin_mev=80.0), "does not define the m2 construction"),
             ("rb", lambda s: s["exploratory"].update(pp_level="said-anchored-ge150"), "does not define the rb construction"),
             ("pinned", lambda s: s.update(exploratory={"d_e3_variant": "all"}), "pinned table differs")]  # fmt: skip
    for i, (v, mutate, msg) in enumerate(cases):
        inp3 = _write_inputs(
            tmp_path / f"case{i}" if (tmp_path / f"case{i}").mkdir() is None else tmp_path
        )
        MP.setattr(c, "TABLE_ID", inp3["expected"]["pinned"]["table_id"])
        new_exp = _retable(inp3, v, mutate)
        if v == "pinned":
            (inp3["tmp"] / "pin.json").write_text(json.dumps(new_exp["pinned"]))
            MP.setattr(c, "TABLE_ID", new_exp["pinned"]["table_id"])
        with pytest.raises(SystemExit, match=msg):
            rec.verify_expected_tables(new_exp, inp3["scratch"], inp3["de7"])


def test_b0_record_is_bound_authenticated_and_schema_checked(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path)
    assert _build(inp)["complete"] is True
    good = json.loads(inp["b0"].read_text())

    def put(**kv: Any) -> None:
        inp["b0"].write_text(json.dumps({**good, **kv}))
        _relog(inp)

    for key in (
        "match",
        "dip",
        "reading",
        "carry_up_check",
        "quadrature_convergence",
        "provenance",
        "named_nodes",
    ):
        inp["b0"].write_text(json.dumps({k: v for k, v in good.items() if k != key}))
        _relog(inp)
        with pytest.raises(SystemExit, match="lacks"):
            _build(inp)
    qc = good["quadrature_convergence"]
    nan_nodes = json.loads(json.dumps(qc["nodes"]))
    nan_nodes["125"]["gl1024_vs_gl2048"] = float("nan")
    big_nodes = json.loads(json.dumps(qc["nodes"]))
    big_nodes["150"]["simpson4001_vs_gl2048"] = 1e-5
    small = json.loads(json.dumps(qc["nodes"]))
    small["100"]["simpson20001_vs_gl2048"] = 4e-7
    for qbad, msg in (
        ({"accepted": True}, "lacks the acceptance threshold"),
        ({**qc, "accepted": "yes"}, "disagrees with the acceptance"),
        ({**qc, "acceptance": 1e-3}, "lacks the acceptance threshold"),
        ({**qc, "nodes": nan_nodes}, "not a finite number"),
        (
            {**qc, "nodes": big_nodes, "max_abs_rel_diff": 1e-5, "accepted": True},
            "disagrees with the acceptance",
        ),
        ({**qc, "nodes": small}, "not the maximum"),
        (
            {**qc, "nodes": {k: v for k, v in qc["nodes"].items() if k != "125"}},
            "lacks the acceptance threshold",
        ),
        ({**qc, "nodes": {**qc["nodes"], "100": {"gl1024_vs_gl2048": 0.0}}}, "lacks the GL 1024"),
        ({**qc, "accepted": not qc["accepted"]}, "disagrees with the acceptance"),
        ({k: v for k, v in qc.items() if k != "acceptance"}, "lacks the acceptance threshold"),
    ):
        put(quadrature_convergence=qbad)
        with pytest.raises(SystemExit, match=msg):
            _build(inp)
    put(match={"all_within": "maybe"})
    with pytest.raises(SystemExit, match="malformed"):
        _build(inp)
    for k, v, msg in (
        ("code_sha256", "", "lacks code_sha256"),
        ("table_id", "0" * 64, "differs from the expected pinned table"),
        ("npz_sha256", "0" * 64, "differs from the expected pinned table"),
        ("source_sha256", "0" * 64, "frozen pins"),
        ("member_sha256", "0" * 64, "frozen pins"),
        ("pre_run_looks", "edited", "frozen text"),
        ("member_sha256", "", "lacks member_sha256"),
    ):
        put(provenance={**good["provenance"], k: v})
        with pytest.raises(SystemExit, match=msg):
            _build(inp)
    put()  # restore
    assert _build(inp)["complete"] is True
    # authentication against the attempt log
    inp["b0"].write_text(json.dumps({**good, "reading": "edited after publication"}))
    with pytest.raises(SystemExit, match="attempt-log entry"):
        _build(inp)
    inp["attempts"].write_text(
        json.dumps(
            [
                r
                for r in json.loads(inp["attempts"].read_text())
                if "de9-b0.json" not in str(r.get("output"))
            ]
        )
    )
    with pytest.raises(SystemExit, match="no \\(unique\\) published entry"):
        _build(inp)
    # the code binding: the B0 record must carry the reviewed de9_pp_level.py hash
    _relog(inp)
    inp["b0"].write_text(json.dumps(good))
    _relog(inp)
    fz = {**_freeze(), "scripts": {"de3_probes.py": "1" * 64, "de9_pp_level.py": "2" * 64}}
    with pytest.raises(SystemExit, match="code_sha256 is missing or differs"):
        _build(inp, freeze=fz)


def test_unaccepted_quadrature_convergence_makes_the_record_incomplete(tmp_path: Path) -> None:
    inp = _write_inputs(tmp_path)
    good = json.loads(inp["b0"].read_text())
    qc = json.loads(json.dumps(good["quadrature_convergence"]))
    qc["nodes"]["125"]["simpson4001_vs_gl2048"] = (
        2.0e-6  # a consistent, genuinely failing record: above the pre-declared 1e-6
    )
    qc["max_abs_rel_diff"], qc["accepted"] = 2.0e-6, False
    inp["b0"].write_text(json.dumps({**good, "quadrature_convergence": qc}))
    _relog(inp)
    doc = _build(inp)
    assert (
        doc["complete"] is False
        and any("quadrature convergence not accepted" in r for r in doc["incomplete_reasons"])
        and doc["b0"]["quadrature_accepted"] is False
    )
    md = rec.render_report(
        c.git_projection(json.loads(json.dumps(doc, default=float)), rec.GIT_ALLOW)
    )
    assert (
        "INCOMPLETE — not the D-E9 synthesis" in md and "quadrature convergence not accepted" in md
    )


def _b0_cache(tmp: Path, member: str | None = None) -> Path:
    """A synthetic data cache that ``de9_pp_level.run`` can read: the LA150-like zip object, the pinned sidecar and the pinned node values."""
    import zipfile

    MP.setattr(b0, "N_GAUSS", 64)
    MP.setattr(c, "TABLE_ID", "t" * 64)
    cd = tmp / "cache"
    (cd / "derived").mkdir(parents=True)
    (cd / "objects").mkdir()
    sec = _sec([4.0e-3, 1.0e-3, 4.0e-3], energies=[1.0e7, 1.25e8, 1.5e8])
    zp = tmp / "la150.zip"
    with zipfile.ZipFile(zp, "w") as z:
        z.writestr(
            b0.MEMBER,
            member
            if member is not None
            else _synthetic_endf([1.0e7, 1.25e8, 1.5e8], [4.0e-3, 1.0e-3, 4.0e-3]),
        )
    zsha = hashlib.sha256(zp.read_bytes()).hexdigest()
    (cd / "objects" / zsha).write_bytes(zp.read_bytes())
    grid = np.arange(100.0, 251.0, 10.0)
    sig = np.zeros((2, grid.size))
    from ionmc.nuclear import bgg

    rec150 = b0.sigma_ni_barn(sec, 150.0)
    sig[0] = [
        b0.sigma_ni_barn(sec, float(e))
        if e <= 150.0
        else rec150 * bgg.pp_bgg_mb(float(e)) / bgg.pp_bgg_mb(150.0)
        for e in grid
    ]
    np.savez(cd / "derived" / f"elastic-proton-{'t' * 64}.npz", grid_e_mev=grid, sigma_barn=sig)
    side = {
        "table_id": "t" * 64,
        "npz_sha256": c.sha256_file(cd / "derived" / f"elastic-proton-{'t' * 64}.npz"),
        "sources": {b0.SOURCE_KEY: {"sha256": zsha}},
    }
    (cd / "derived" / f"elastic-proton-{'t' * 64}.json").write_text(json.dumps(side))
    return cd


def test_b0_run_is_once_only_and_the_attempt_precedes_any_work(tmp_path: Path) -> None:
    cd = _b0_cache(tmp_path)
    out, att = tmp_path / "out", tmp_path / "att.json"
    seen: list[bool] = []
    real_locate = b0.locate_raw

    def spy(cache_dir: Path) -> Any:
        seen.append(
            att.is_file() and json.loads(att.read_text())[-1]["outcome"] == "started"
        )  # the attempt record exists BEFORE the inputs are located
        return real_locate(cache_dir)

    MP.setattr(b0, "locate_raw", spy)
    doc = b0.run(cd, out, att)
    assert seen == [True] and doc["quadrature_convergence"]["accepted"] is True
    log = json.loads(att.read_text())
    assert [r["outcome"] for r in log] == ["started", "published-valid"] and log[0][
        "producer"
    ] == "de9_pp_level.run"
    assert (
        rec.attempt_log_sha(log, "de9-b0.json")
        == hashlib.sha256((out / "git" / "de9-b0.json").read_bytes()).hexdigest()
    )
    MP.setattr(
        b0,
        "locate_raw",
        lambda cache_dir: (_ for _ in ()).throw(
            AssertionError("no work may precede the published-target refusal")
        ),
    )
    with pytest.raises(SystemExit, match="exists"):
        b0.run(
            cd, out, att
        )  # a published record is refused before any work, never overwritten, never re-run
    assert len(json.loads(att.read_text())) == 2


def test_b0_run_failures_are_recorded_and_bounded_by_two_attempts(tmp_path: Path) -> None:
    cd = _b0_cache(tmp_path, member="not an ENDF section\n")
    cd2 = _b0_cache(
        tmp_path / "second" if (tmp_path / "second").mkdir() is None else tmp_path
    )  # built before parse_h1_mt2 is patched
    out, att = tmp_path / "out", tmp_path / "att.json"
    for _ in range(2):
        with pytest.raises(Exception) as exc:
            b0.run(cd, out, att)
        assert not isinstance(exc.value, SystemExit)  # a parse failure, not a refusal
        assert not (out / "git" / "de9-b0.json").exists()
    assert [r["outcome"] for r in json.loads(att.read_text())] == [
        "started",
        "failed",
        "started",
        "failed",
    ]
    parsed: list[str] = []
    with (
        pytest.MonkeyPatch.context() as m,
        pytest.raises(SystemExit, match="two attempts without a published document"),
    ):
        m.setattr(b0, "parse_h1_mt2", lambda text: parsed.append(text) or {})
        b0.run(cd, out, att)
    assert not parsed and len(json.loads(att.read_text())) == 4  # no further attempt, no work
    # a leak-check failure after the calculation is a failed attempt as well (no document)
    MP.setattr(b0, "leak_violations", lambda proj, absolutes: ["forced"])
    with pytest.raises(SystemExit, match="leak check failed"):
        b0.run(cd2, tmp_path / "out2", tmp_path / "att2.json")
    assert [r["outcome"] for r in json.loads((tmp_path / "att2.json").read_text())] == [
        "started",
        "failed",
    ] and not (tmp_path / "out2" / "git" / "de9-b0.json").exists()
