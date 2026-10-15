"""Exploratory D-E0 kernels (Amendment 19, report-only): analytic cases."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validation.scripts.elastic_attribution import de0_kernels as de0  # noqa: E402
from validation.scripts.elastic_attribution import de1_shapes as de1  # noqa: E402


def _lin_range(e: object) -> object:
    return np.asarray(e, dtype=float)  # toy range R(E) = E mm


def _ev(mu: np.ndarray, m_t: float = 1.0e12, e: float = 100.0) -> dict[str, np.ndarray]:
    from ionmc.nuclear.elastic_kin import two_body

    tb = two_body(e, 938.272, m_t, mu)
    return {"mu": mu, "t1": tb.t1, "cos": tb.cos_lab1, "t_other": tb.t2}


def test_isotropic_infinite_target_exceedance() -> None:
    """Isotropic CM, infinite target: Delta = R (1 - mu), P(Delta > w) = 1 - w / (2 R)."""
    u = (np.arange(4096) + 0.5) / 4096
    ev = _ev(2.0 * u - 1.0)
    d = de0.deficit_mm(100.0, ev, _lin_range)
    for w in (1.0, 2.0, 5.0):
        assert np.mean(d > w) == pytest.approx(1.0 - w / 200.0, abs=1e-3)


def test_forward_delta_has_no_deficit() -> None:
    ev = _ev(np.full(64, 1.0))
    assert np.all(de0.deficit_mm(100.0, ev, _lin_range) < 1e-6)
    assert de0.f_pred(0.0, 0.0) == 1.0


def test_threshold_angle() -> None:
    th = de0.theta_w_deg(150.0, 2.0, _lin_range)
    assert 150.0 * (1.0 - np.cos(np.radians(th))) == pytest.approx(2.0)


def test_disk_family_normalised_and_quantiles() -> None:
    """Disk form: 1 at q = 0, first zero at x = 3.8317; the families are 1 at q = 0."""
    assert float(de1.disk_a(np.array([0.0]))[0]) == 1.0
    assert float(de1.disk_a(np.array([de1.J1_FIRST_ZERO]))[0]) < 1e-9
    assert float(de1.family_weight("s0.5", np.array([0.0]), 2.6)[0]) == 1.0
    assert float(de1.family_weight("gauss", np.array([0.0]), 2.6)[0]) == 1.0


@pytest.mark.requires_tables("elastic")
def test_kernels_on_pinned_table() -> None:
    from validation.scripts.elastic_attribution import common as c

    tab = c.load_table()
    sh, water = de1.Shapes(tab), c.Water(tab)
    k = de0.kernels(sh, water, 150, "table", ("O-16",), False, n_z=40)["O-16"]
    assert k["P_el"] > k["K_1mm"] >= k["K_2mm"] >= k["K_3mm"] >= k["K_5mm"] > 0.0
    soft = de0.kernels(sh, water, 150, "s0.9", ("O-16",), False, n_z=40)["O-16"]
    assert soft["K_2mm"] < k["K_2mm"]


def test_kstar_sign_pinned_on_asymmetric_idd() -> None:
    """Synthetic IDD (slow rise, sharp fall): shortened range gives d(z_pk) = IDD(z_pk + Delta)
    (shifted, physical); the as-frozen z_pk - Delta samples the rise."""

    def idd(z: object) -> object:
        z = np.asarray(z, dtype=float)
        return np.where(z <= 100.0, z / 100.0, np.clip(1.0 - (z - 100.0) / 10.0, 0.0, None))

    lit, shifted = de0.kstar_node(idd, 100.0, np.array([5.0, 5.0]))
    assert shifted == pytest.approx(1.0 - 0.5)  # IDD(105) = 0.5
    assert lit == pytest.approx(1.0 - 0.95)  # IDD(95) = 0.95
    assert shifted > lit


# --- Git projection of the exploratory records (decision 0041: explicit allowlist) ---
import fnmatch  # noqa: E402
import importlib  # noqa: E402
import json  # noqa: E402

from validation.scripts.elastic_attribution import common as de_common  # noqa: E402

RESULTS = Path(__file__).resolve().parents[2] / "validation" / "results" / "elastic_attribution"
STEP_MODULES = {
    "de0": "de0_kernels",
    "de0h5": "de0_h5",
    "de1": "de1_shapes",
    "de2": "de2_pp",
    "de5": "de5_sampling",
}  # steps with a rendered ``<step>-report.md``
OWNERS = {
    **STEP_MODULES,
    "de3": "de3_probes",
}  # every committed JSON directory and the module that owns its allowlist


def test_allowlist_projection_semantics() -> None:
    doc = {
        "banner": "x",
        "note": "n",
        "per_energy": {
            "100": {
                "ratio": 1.1,
                "sigma_el_mb": 99.0,
                "ci": [1.0, 2.0],
                "hist": [1, 2, 3],
                "kernel": [1, 2, 3],
            }
        },
        "files": {"a/b/c.json": "h"},
        "other": 3.0,
        "xs_mb_sr": 1.0,
    }
    allow = (
        "per_energy/*/ratio",
        "per_energy/*/ci",
        "per_energy/*/hist",
        "per_energy/*/kernel",
        "files/**",
        "xs_mb_sr",
    )
    proj = de_common.git_projection(doc, allow, allow_long=("per_energy/*/kernel",))
    assert proj == {
        "banner": "x",
        "note": "n",
        "per_energy": {"100": {"ratio": 1.1, "ci": [1.0, 2.0], "kernel": [1, 2, 3]}},
        "files": {"a/b/c.json": "h"},
    }  # absolute cross-section names, the long list "hist" and the unlisted "other" are dropped
    assert de_common.git_projection(proj, allow, ("per_energy/*/kernel",)) == proj
    assert de_common.git_safety_violations(proj, allow, ("per_energy/*/kernel",)) == []
    bad = de_common.git_safety_violations(doc, allow)
    assert (
        any("sigma_el_mb" in b for b in bad)
        and any("hist" in b for b in bad)
        and any("other" in b for b in bad)
    )
    assert de_common.path_matches("a/b/c", ("a/*/c",)) and not de_common.path_matches(
        "a/b", ("a/*/c", "a/b/**")
    )
    assert not de_common.path_matches("a/b/c", ("a/b",))


def _allow_for(step: str, name: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    mod = importlib.import_module(f"validation.scripts.elastic_attribution.{OWNERS[step]}")
    for pat, allow in mod.GIT_FILES.items():
        if fnmatch.fnmatch(name, pat):
            return allow  # type: ignore[no-any-return]
    raise AssertionError(f"{OWNERS[step]} has no allowlist for {name}")


@pytest.mark.parametrize(
    "path", sorted(RESULTS.glob("*/*.json")), ids=lambda p: str(p.relative_to(RESULTS))
)
def test_committed_records_match_the_module_allowlist(path: Path) -> None:
    """Every committed leaf matches the module ``GIT_ALLOW``; no forbidden key, cross-section
    name or long numeric list, and the projection leaves it unchanged."""
    allow, allow_long = _allow_for(path.parent.name, path.name)
    doc = json.loads(path.read_text())
    assert de_common.git_safety_violations(doc, allow, allow_long) == []
    proj = de_common.git_projection(doc, allow, allow_long)
    assert doc == proj


@pytest.mark.parametrize("step", sorted(STEP_MODULES))
def test_committed_markdown_is_rendered_from_the_projection(step: str) -> None:
    """The committed report body equals the module rendering of the committed (projected) JSON."""
    mod = importlib.import_module(f"validation.scripts.elastic_attribution.{STEP_MODULES[step]}")
    (json_path,) = [
        p for p in sorted((RESULTS / step).glob("*.json")) if p.name != "quarantine-data-seen.json"
    ]
    md = (RESULTS / step / f"{step}-report.md").read_text()
    assert md.count(de_common.BODY_MARKER) == 1
    body = md.split(de_common.BODY_MARKER + "\n", 1)[1]
    assert body == mod.render_report(json.loads(json_path.read_text()))


def _render_pairs() -> list[tuple[str, str, str, str]]:
    """(step dir, json, md, "module:function") of every committed record with a report."""
    out = []
    for step, mod in sorted(STEP_MODULES.items()):
        (jp,) = [
            p for p in (RESULTS / step).glob("*.json") if p.name != "quarantine-data-seen.json"
        ]
        out.append((step, jp.name, f"{step}-report.md", f"{mod}:render_report"))
    for p in sorted((RESULTS / "de3").glob("*.json")):
        if p.name == "de3-summary.json":
            out.append(("de3", p.name, "de3-summary.md", "de3_probes:render_summary"))
        elif p.name.startswith("check-base-"):
            out.append(
                (
                    "de3",
                    p.name,
                    p.name.replace(".json", "-report.md"),
                    "de3_probes:render_check_base",
                )
            )
    for p in sorted((RESULTS / "de4").glob("*.json")) if (RESULTS / "de4").is_dir() else []:
        if p.name == "de4.json":
            out.append(("de4", p.name, "de4-report.md", "de4_parse:render_report"))
        else:
            out.append(
                ("de4", p.name, p.name.replace(".json", "-report.md"), "de4_parse:render_case")
            )
    return out


@pytest.mark.parametrize("pair", _render_pairs(), ids=lambda t: f"{t[0]}/{t[1]}")
def test_every_committed_report_is_rendered_from_its_projection(
    pair: tuple[str, str, str, str],
) -> None:
    """The committed Markdown body equals the module rendering of the committed projected JSON;
    the header lists input files with sha256."""
    step, jname, mname, target = pair
    modname, fn = target.split(":")
    mod = importlib.import_module(f"validation.scripts.elastic_attribution.{modname}")
    md = (RESULTS / step / mname).read_text()
    assert md.count(de_common.BODY_MARKER) == 1
    head, body = md.split(de_common.BODY_MARKER + "\n", 1)
    assert body == getattr(mod, fn)(json.loads((RESULTS / step / jname).read_text()))
    assert "Input files (sha256):" in head


def _run_file(
    variant: str, energy: int, lateral_valid: bool | None, residual: float | None = 0.0
) -> dict:
    rng = np.random.default_rng(7 if variant == "emonly" else 8)
    nz = 200
    base = np.concatenate([np.full(120, 5.0), np.linspace(5, 30, 30), np.linspace(30, 0, 50)])[:nz]
    idd = (
        base[None, :]
        * (1.0 + 0.01 * rng.standard_normal((20, nz)))
        * (0.9 if variant != "emonly" else 1.0)
    )
    mode = "emonly" if variant == "emonly" else "emel"
    doc = {"variant": variant, "mode": mode, "energy": energy,
           "seed": 1, "n": 1000, "n_batches": 20, "elastic_table_id": "x",
           "idd_batches": idd.tolist(), "valid": True, "wall_s": 1.0,
           "metrics": {"r80_mm": 100.0}, "recoil_total_mev_per_primary": 0.1}  # fmt: skip
    if residual is not None:
        doc["relative_residual"] = residual
    if lateral_valid is not None:
        doc["lateral"] = {
            "valid": lateral_valid,
            "relative_residual": 0.0,
            "rms_mm_batches": [1.0 + 0.01 * i for i in range(20)],
            "wall_s": 1.0,
        }
    return doc


def test_summary_validity_is_fail_closed(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de3_probes

    def rows(**kw: object) -> dict:
        out = tmp_path / str(len(list(tmp_path.iterdir())))
        out.mkdir()
        (out / "de3-150-emonly.json").write_text(json.dumps(_run_file("emonly", 150, True)))
        (out / "de3-150-synthetic.json").write_text(json.dumps(_run_file("synthetic", 150, **kw)))  # type: ignore[arg-type]
        return de3_probes.summarize(out)["rows"]["150/synthetic"]

    assert rows(lateral_valid=True)["valid"] is True
    assert rows(lateral_valid=False)["valid"] is False  # lateral run invalid -> row invalid
    assert (
        rows(lateral_valid=True, residual=1e-6)["valid"] is False
    )  # energy balance of the main run does not close
    assert (
        rows(lateral_valid=True, residual=None)["valid"] is False
    )  # residual not recorded -> invalid
    assert rows(lateral_valid=None)["valid"] is True  # a run without lateral part stays valid


def test_contrast_jackknife_of_a_multiplicative_composition() -> None:
    from validation.scripts.elastic_attribution import de3_probes

    cmp = de3_probes._cmp()
    rng = np.random.default_rng(3)
    nz = 200
    base = np.concatenate([np.full(120, 5.0), np.linspace(5, 30, 30), np.linspace(30, 0, 50)])[:nz]
    noise = 1.0 + 0.01 * rng.standard_normal((20, nz))
    den = base[None, :] * noise

    # peak and plateau scaled independently: F is multiplicative, the removal is not additive
    def var(peak: float, plat: float) -> np.ndarray:
        out = den.copy()
        out[:, :120] *= plat
        out[:, 120:] *= peak
        return out

    c_all, c_pa, c_pp = contrasts_args = (var(0.8 * 0.9, 1.0), var(0.8, 1.0), var(0.9, 1.0))
    res = de3_probes.contrasts(cmp, c_all, c_pa, c_pp, den)
    assert abs(res["F_all_minus_product"]) < 1e-12 and res["consistent_with_multiplicative"] is True
    assert res["removal_additivity"] == pytest.approx(
        0.2 + 0.1 - (1 - 0.72), abs=1e-9
    )  # 0.02: removal is not additive
    assert contrasts_args[0].shape == den.shape


def test_kin_curve_matches_two_body() -> None:
    from ionmc.nuclear.elastic_kin import two_body
    from validation.scripts.elastic_attribution import de5_sampling as de5

    mu = np.linspace(-0.96, 0.99, 400)
    for m_t in (938.272, 14903.0):
        tb = two_body(150.0, 938.272, m_t, mu)
        sel = tb.cos_lab1 > 0.3  # near 90 deg the root is ill-conditioned
        curve = de5.kin_curve_t(150.0, m_t, tb.cos_lab1[sel])
        assert np.max(np.abs(curve - tb.t1[sel]) / tb.t1[sel]) < 1e-5


def test_classification_rule_on_synthetic_tracks() -> None:
    from validation.scripts.elastic_attribution import de5_sampling as de5

    e0, m_h, m_a = 150.0, 938.272, 14903.0
    # conflict track (theta_lab 0.5 deg, on both curves within 1e-3): follows the secondary flag
    cos = np.cos(np.radians(np.full(2, 0.5)))
    t1 = de5.kin_curve_t(e0, m_a, cos)
    is_pp, conflict, off = de5.classify(cos, t1, e0, m_h, m_a, np.array([True, False]))
    assert conflict.all() and not off.any() and list(is_pp) == [True, False]
    # on-curve p-p track at 30 deg (far from the p+A curve), with and without the flag
    cos = np.cos(np.radians(np.array([30.0])))
    t1 = de5.kin_curve_t(e0, m_h, cos)
    is_pp, conflict, off = de5.classify(cos, t1, e0, m_h, m_a, True)
    assert is_pp[0] and not conflict[0] and not off[0]
    # on-curve p+A track at 30 deg
    t1 = de5.kin_curve_t(e0, m_a, cos)
    is_pp, conflict, off = de5.classify(cos, t1, e0, m_h, m_a, False)
    assert not is_pp[0] and not conflict[0] and not off[0]
    # off-curve track: counted as off-curve, classified by the smaller residual
    is_pp, conflict, off = de5.classify(cos, t1 * 0.9, e0, m_h, m_a, False)
    assert off[0] and not conflict[0]
    # crossover of the curves is about 1.8 deg at 150 MeV (p+O)
    assert 1.0 < de5.curve_crossover_deg(e0, m_h, m_a) < 3.0


# --- D-E4 parser and the frozen track selection (synthetic phase spaces) ---
HEADER = """TOPAS Binary Phase Space

Number of Original Histories: {hist}
Number of Original Histories that Reached Phase Space: {hist}
Number of Scored Particles: {n}
Number of Bytes per Particle: 46

Byte order of each record is as follows:
f4: Position X [cm]
f4: Position Y [cm]
f4: Position Z [cm]
f4: Direction Cosine X
f4: Direction Cosine Y
f4: Energy [MeV]
f4: Weight
i4: Particle Type (in PDG Format)
b1: Flag to tell if Third Direction Cosine is Negative (1 means true)
b1: Flag to tell if this is the First Scored Particle from this History (1 means true)
i4: Event ID
i4: Track ID
i4: Parent ID

Number of proton: {n}

Minimum Kinetic Energy of proton: 0.5 MeV

Maximum Kinetic Energy of proton: {emax} MeV
"""
INPUT = """i:Ts/Seed = 1
s:Ge/Target/Material = "{mat}"
d:Ge/Target/HL = {hl} mm
d:So/Beam/BeamEnergy = {e0} MeV
i:So/Beam/NumberOfHistoriesInRun = {hist}
"""


def _write_ps(d: Path, recs: list[tuple], hist: int, e0: float, mat: str, hl_mm: float) -> Path:
    """Write a TOPAS binary phase space. ``recs``: (event, track, parent, E, cos, w_negative)."""
    d.mkdir(parents=True, exist_ok=True)
    f4, i4 = "<f4", "<i4"
    names = [(n, f4) for n in ("X", "Y", "Z", "U", "V", "E", "w")]
    dt = np.dtype(
        [*names, ("pdg", i4), ("neg", "u1"), ("first", "u1"), ("ev", i4), ("tr", i4), ("pa", i4)]
    )
    arr = np.zeros(len(recs), dtype=dt)
    for i, (ev, tr, pa, e, cos, neg) in enumerate(recs):
        arr[i] = (
            0,
            0,
            0,
            math.sqrt(max(1 - cos * cos, 0.0)),
            0,
            e,
            1,
            2212,
            neg,
            int(pa == 0),
            ev,
            tr,
            pa,
        )
    assert dt.itemsize == 46
    (d / "ps.phsp").write_bytes(arr.tobytes())
    (d / "ps.header").write_text(HEADER.format(hist=hist, n=len(recs), emax=e0))
    (d / "input.txt").write_text(INPUT.format(mat=mat, hl=hl_mm, e0=e0, hist=hist))
    return d


def _tracks(
    e0: float, m_a: float, ev0: int, n_pa: int = 6, n_pp: int = 4, bad: bool = True
) -> list[tuple]:
    """Synthetic events: undeflected primaries, on-curve p+A and p-p events, and (``bad``) one
    off-curve p+A primary and one off-curve (plural) p-p event."""
    import validation.scripts.elastic_attribution.de5_sampling as de5
    from ionmc.nuclear.elastic_kin import two_body

    recs, ev = [], ev0
    for _ in range(5):
        recs.append((ev, 1, 0, e0, 1.0, 0))
        ev += 1
    for k in range(n_pa):
        cos = math.cos(math.radians(6.0 + 2.0 * k))
        recs.append((ev, 1, 0, float(de5.kin_curve_t(e0, m_a, np.array([cos]))[0]), cos, 0))
        ev += 1
    for k in range(n_pp):  # equal masses: both protons from two_body
        tb = two_body(e0, 938.272, 938.272, np.array([0.1 + 0.15 * k]))
        cos2 = float(tb.pz2[0] / np.hypot(tb.px2[0], tb.pz2[0]))
        recs += [
            (ev, 1, 0, float(tb.t1[0]), float(tb.cos_lab1[0]), 0),
            (ev, 2, 1, float(tb.t2[0]), cos2, 0),
        ]
        ev += 1
    if bad:
        cos = math.cos(math.radians(10.0))
        recs.append(
            (ev, 1, 0, 0.9 * float(de5.kin_curve_t(e0, m_a, np.array([cos]))[0]), cos, 0)
        )  # off both curves
        ev += 1
        tb = two_body(e0, 938.272, 938.272, np.array([0.2]))
        cos2 = float(tb.pz2[0] / np.hypot(tb.px2[0], tb.pz2[0]))
        recs += [
            (ev, 1, 0, 0.9 * float(tb.t1[0]), float(tb.cos_lab1[0]), 0),
            (ev, 2, 1, float(tb.t2[0]), cos2, 0),
        ]  # plural
    return recs


def _mass(name: str) -> float:
    from validation.scripts.elastic_attribution import common as c

    tab = c.load_table()
    return float(tab.arrays["target_mass_mev"][list(tab.target_names).index(name)])


@pytest.mark.requires_tables("elastic")
def test_de4_water_phase_space_is_identified_and_off_curve_tracks_are_excluded(
    tmp_path: Path,
) -> None:
    from validation.scripts.elastic_attribution import de4_parse as de4

    d = _write_ps(
        tmp_path / "w", _tracks(150.0, _mass("O-16"), 1), 4_000_000, 150.0, "G4_WATER", 5.0
    )
    doc = de4.analyse([d])
    assert (doc["case"], doc["row"], doc["material"], doc["target"], doc["length_cm"]) == (
        3,
        17,
        "water",
        "O-16",
        1.0,
    )
    assert doc["ke_filter"]["removed_primaries"] == 5
    assert doc["classes"]["O-16"]["events"] == 6  # the off-curve primary is excluded, not counted
    assert doc["off_curve"]["primaries_excluded"] == 1 and doc["pp_events_off_curve"] == 1
    assert (
        doc["classes"]["H-1"]["events"] + doc["pp_events_removed_by_cm_cut"] == 4
    )  # the plural p-p event is excluded
    assert doc["classes"]["O-16"]["rate_per_primary_after_rule"] == pytest.approx(6 / 4_000_000)
    assert set(doc["geant4_summary"]["3"]) == {"O-16", "H-1"}


@pytest.mark.requires_tables("elastic")
def test_de4_graphite_phase_space_has_only_the_carbon_class(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de4_parse as de4

    d = _write_ps(
        tmp_path / "g",
        _tracks(150.0, _mass("C-12"), 1, n_pp=0, bad=False),
        4_000_000,
        150.0,
        "G4_C",
        5.0,
    )
    doc = de4.analyse([d])
    assert (doc["case"], doc["material"], doc["target"]) == (8, "graphite", "C-12")
    assert list(doc["classes"]) == ["C-12"] and doc["classes"]["C-12"]["events"] == 6


@pytest.mark.requires_tables("elastic")
def test_de4_split_case_parts_are_concatenated_and_validated(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de4_parse as de4

    hist = (6_666_667, 6_666_667, 6_666_666)
    parts = [
        _write_ps(
            tmp_path / f"p{k}",
            _tracks(30.0, _mass("O-16"), 1, n_pa=3, n_pp=2, bad=False),
            h,
            30.0,
            "G4_WATER",
            1.0,
        )
        for k, h in enumerate(hist)
    ]  # identical event ids in every part
    doc = de4.analyse(parts)
    assert (
        doc["case"] == 0 and doc["primaries"] == 20_000_000 and doc["identification"]["parts"] == 3
    )
    assert doc["ke_filter"]["removed_primaries"] == 15 and doc["classes"]["O-16"]["events"] == 9
    assert doc["pp_events_all"] == 6  # event ids of different parts do not collide
    with pytest.raises(AssertionError):
        de4.analyse(parts[:2])  # the primaries of two parts do not add up to the case
    with pytest.raises(AssertionError):
        de4.analyse(parts, case=1)


def test_de4_header_layout_is_parsed_not_hard_coded(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de4_parse as de4

    d = _write_ps(tmp_path / "h", [(1, 1, 0, 150.0, 1.0, 0)], 1, 150.0, "G4_WATER", 5.0)
    rec, h = de4.read_phsp(d)
    assert (
        h["bytes"] == 46
        and rec.dtype.names is not None
        and "parent" in rec.dtype.names
        and float(rec["E"][0]) == 150.0
    )


def test_track_selection_excludes_off_curve_tracks_on_both_sides() -> None:
    """KE filter + off-curve exclusion of the shared selection (table side and Geant4 parser)."""
    from validation.scripts.elastic_attribution import de5_sampling as de5

    e0, m_h, m_a = 150.0, 938.272, 14903.0
    cos = np.cos(np.radians(np.array([0.0, 10.0, 10.0, 30.0])))
    curve_a, curve_h = de5.kin_curve_t(e0, m_a, cos), de5.kin_curve_t(e0, m_h, cos)
    t1 = np.array([e0, curve_a[1], 0.9 * curve_a[2], curve_h[3]])
    sel = de5.select_tracks(cos, t1, e0, m_h, m_a, np.array([False, False, False, True]))
    assert list(sel["ke_removed"]) == [True, False, False, False]
    assert list(sel["off"]) == [False, False, True, False]  # off both curves: counted and excluded
    assert list(sel["keep"]) == [False, True, False, True]
    assert list(sel["is_pp"][[1, 3]]) == [False, True]


@pytest.mark.requires_tables("elastic")
def test_table_side_applies_the_same_selection() -> None:
    from validation.scripts.elastic_attribution import common as c
    from validation.scripts.elastic_attribution import de5_sampling as de5

    tab = c.load_table()
    cs = de5.CASES[3]
    de5.N_EVENTS = 20_000
    try:
        doc = de5.sample_case(tab, c.Water(tab), cs)
    finally:
        de5.N_EVENTS = 1_000_000
    for v in doc[
        "classes"
    ].values():  # table events lie on a curve by construction: nothing is excluded
        assert v["off_curve_fraction_excluded"] == 0.0
        assert v["rate_per_primary_after_rule"] <= v["rate_per_primary_after_ke_cut"] * (1 + 1e-12)
    assert [cs["case"] for cs in de5.CASES] == list(range(10)) and de5.CASES[9]["row"] == 23
