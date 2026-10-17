"""Exploratory D-E0 kernels (Amendment 19, report-only): analytic cases."""
# ruff: noqa: E501  (parametrised mutation tables and synthetic fixtures)

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
import copy  # noqa: E402
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
    "de4": "de4_parse",
    "de6": "de6_swap",
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
    for p in (
        sorted((RESULTS / "de6").glob("de6-summary.json")) if (RESULTS / "de6").is_dir() else []
    ):
        out.append(("de6", p.name, "de6-report.md", "de6_swap:render_report"))
    for p in sorted((RESULTS / "de4").glob("*.json")) if (RESULTS / "de4").is_dir() else []:
        if p.name == "de4-ratios.json":
            out.append(("de4", p.name, "de4-ratios-report.md", "de4_parse:render_ratios"))
        elif p.name == "de4.json":
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
    assert md.endswith("\n") and not md.endswith("\n\n")  # exactly one final newline
    rendered = getattr(mod, fn)(json.loads((RESULTS / step / jname).read_text()))
    assert body.rstrip("\n") == rendered.rstrip("\n")
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


@pytest.mark.requires_tables("elastic")
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
BUNDLES = Path(__file__).resolve().parents[2] / "validation" / "reference_cases" / "topas"


def _native_text(
    mat: str,
    e0: float,
    hl_mm: float,
    seed: int,
    input_hist: int,
    modules: str,
    rmax: float,
    sphere: float,
    edits: tuple[tuple[str, str], ...],
) -> str:
    """The committed bundle input with seed / histories substituted and the mutations."""
    import re

    from validation.scripts.elastic_attribution import de4_parse as de4
    from validation.scripts.elastic_attribution import de5_sampling as de5

    cands = [cs for cs in de5.CASES if cs["topas_material"] == mat and cs["energy_mev"] == e0]
    cs = next((x for x in cands if abs(x["length_cm"] * 5.0 - hl_mm) < 1e-6), cands[0])
    j = seed - cs["topas_seed_base"]
    j = j if 1 <= j <= len(cs["part_histories"]) else 1
    text = (BUNDLES / de4.bundle_dirname(cs, j) / "input.txt").read_text()

    def sub(name: str, value: str) -> None:
        nonlocal text
        text, n = re.subn(
            rf"^({re.escape(name)}\s*=\s*).*$", lambda m: m.group(1) + value, text, flags=re.M
        )
        assert n == 1, name

    sub("i:Ts/Seed", str(seed))
    sub("i:So/Beam/NumberOfHistoriesInRun", str(input_hist))
    sub("sv:Ph/Default/Modules", modules)
    sub("d:Ge/Target/RMax", f"{rmax:g} mm")
    sub("d:Ge/Sphere/RMax", f"{sphere:g} mm")
    sub("d:Ge/Target/HL", f"{hl_mm:g} mm")
    for name, value in edits:
        sub(name, value)
    return text


def _write_ps(
    d: Path,
    recs: list[tuple],
    hist: int,
    e0: float,
    mat: str,
    hl_mm: float,
    seed: int = 20287001,
    input_hist: int | None = None,
    modules: str = '1 "g4h-elastic_HP"',
    rmax: float = 10.0,
    sphere: float = 500.0,
    edits: tuple[tuple[str, str], ...] = (),
) -> Path:
    """Write a TOPAS phase space and input.txt; ``recs``: (event, track, parent, E, cos, w_neg)."""
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
    ih = hist if input_hist is None else input_hist
    (d / "input.txt").write_text(
        _native_text(mat, e0, hl_mm, seed, ih, modules, rmax, sphere, edits)
    )
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
    assert set(doc["geant4_summary"]["3"]) == {"O-16", "H-1", "_selection_rule"}
    assert doc["selection_rule"] == "fast-proton-on-curve-v2"


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
        seed=20292001,
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
            seed=20284001 + k,
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


def _water_part(tmp_path: Path, name: str = "w", **kw: object) -> Path:
    return _write_ps(
        tmp_path / name,
        _tracks(150.0, _mass("O-16"), 1, bad=False),
        4_000_000,
        150.0,
        "G4_WATER",
        5.0,
        **kw,
    )  # type: ignore[arg-type]


@pytest.mark.requires_tables("elastic")
@pytest.mark.parametrize(
    "kw",
    [
        {
            "input_hist": 4_000_001
        },  # input.txt histories differ from the header (and the frozen part)
        {"seed": 20287002},  # not 20270000 + 1000 row + j
        {"modules": '2 "g4h-elastic_HP" "g4em-standard_opt4"'},  # extra module
        {"rmax": 12.0},  # wrong target radius
        {"sphere": 400.0},  # wrong sphere radius
    ],
    ids=["histories", "seed", "modules", "target-radius", "sphere-radius"],
)
def test_de4_identification_binds_seed_modules_geometry_histories(tmp_path: Path, kw: dict) -> None:
    from validation.scripts.elastic_attribution import de4_parse as de4

    de4.analyse([_water_part(tmp_path)])  # the unmodified part is accepted
    with pytest.raises(AssertionError):
        de4.analyse([_water_part(tmp_path, "bad", **kw)])


@pytest.mark.requires_tables("elastic")
def test_de4_identification_rejects_wrong_header_and_wrong_thickness(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de4_parse as de4

    d = _water_part(tmp_path, "hdr")
    hdr = (
        (d / "ps.header")
        .read_text()
        .replace(
            "Number of Original Histories: 4000000", "Number of Original Histories: 3999999", 1
        )
    )
    (d / "ps.header").write_text(hdr)
    with pytest.raises(AssertionError):
        de4.analyse([d])  # header histories differ from input.txt
    d2 = _write_ps(
        tmp_path / "hl",
        _tracks(150.0, _mass("O-16"), 1, bad=False),
        4_000_000,
        150.0,
        "G4_WATER",
        4.0,
    )
    with pytest.raises(AssertionError):
        de4.analyse([d2])  # thickness 8 mm matches no frozen case


@pytest.mark.requires_tables("elastic")
def test_de4_cdf_output_layout_and_equiprobable_edges(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de4_parse as de4
    from validation.scripts.elastic_attribution import de5_sampling as de5

    m_o = _mass("O-16")  # 300 on-curve p+A tracks at distinct angles
    cos = np.cos(np.radians(3.0 + 0.05 * np.arange(300)))
    recs = [
        (i + 1, 1, 0, float(de5.kin_curve_t(150.0, m_o, cos[i : i + 1])[0]), float(cos[i]), 0)
        for i in range(300)
    ]
    d = _write_ps(tmp_path / "w", recs, 4_000_000, 150.0, "G4_WATER", 5.0)
    doc = de4.run([d], tmp_path / "out", cdf_out=tmp_path / "cdf")
    z = np.load(tmp_path / "cdf" / "de4-case3-cdf.npz")
    side = json.loads((tmp_path / "cdf" / "de4-case3-cdf.json").read_text())
    assert side["npz_sha256"] == de4.c.sha256_file(tmp_path / "cdf" / "de4-case3-cdf.npz")
    assert (
        side["classes"]["O-16"]["n"] == doc["classes"]["O-16"]["events"] == z["O-16_mu_sorted"].size
    )
    mu, edges, q = z["O-16_mu_sorted"], z["O-16_edges_mu"], z["O-16_edges_q_mev"]
    assert mu.dtype == np.float64 and np.all(np.diff(mu) >= 0) and edges.size == q.size == 257
    assert (
        edges[0] == mu[0] and edges[-1] == mu[-1] and np.all(np.diff(edges) >= 0)
    )  # end edges = sample min / max
    assert (
        side["classes"]["O-16"]["sample_min"] == mu[0]
        and side["classes"]["O-16"]["tail_edge_gap_lo"] == edges[1] - edges[0]
    )
    counts = np.histogram(mu, bins=edges)[
        0
    ]  # 256 equiprobable bins (the ends are the range limits)
    assert counts.sum() == mu.size and counts.max() - counts.min() <= 2
    p_cm = side["classes"]["O-16"]["p_cm_mev"]
    assert np.all(np.diff(q) <= 1e-9) and q[0] == pytest.approx(
        2.0 * p_cm * np.sin(np.arccos(mu[0]) / 2.0)
    )
    assert (
        "H-1_edges_mu" not in z.files or side["classes"]["H-1"]["symmetrised"] is True
    )  # p-p: both signs of mu


def _case_record(cs: dict, **over: object) -> dict:
    rec = {"case": cs["case"], "row": cs["row"], "material": cs["material"], "target": cs["target"],
           "energy_mev": float(cs["energy_mev"]),
           "length_cm": cs["length_cm"], "primaries": cs["primaries"],
           "identification": {"geometry_bound": True, "native_input_bound": True},
           "selection_rule": "fast-proton-on-curve-v2", "off_curve": {}, "cdf": {"npz_sha256": "ab" * 32}}  # fmt: skip
    rec.update(over)
    return rec


def test_de4_combine_requires_exactly_the_ten_validated_cases(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de4_parse as de4
    from validation.scripts.elastic_attribution import de5_sampling as de5

    def write(n: int, rec: dict) -> None:
        (tmp_path / f"de4-case{n}.json").write_text(json.dumps(rec))

    for cs in de5.CASES:
        write(cs["case"], _case_record(cs))
    assert sorted(de4.validated_case_records(tmp_path), key=int) == [str(i) for i in range(10)]
    (tmp_path / "de4-case7.json").unlink()
    with pytest.raises(SystemExit, match="case 7 .* missing"):
        de4.validated_case_records(tmp_path)  # missing
    write(7, _case_record(de5.CASES[3]))  # a record of case 3 filed as case 7
    with pytest.raises(SystemExit, match="differ from the frozen case"):
        de4.validated_case_records(tmp_path)  # mismatched case fields
    write(7, _case_record(de5.CASES[7], identification={}))
    with pytest.raises(SystemExit, match="validated case identification"):
        de4.validated_case_records(tmp_path)  # not identified against the frozen matrix
    write(7, _case_record(de5.CASES[7], selection_rule="exactly-two-on-curve-v1"))
    with pytest.raises(SystemExit, match="legacy selection schema"):
        de4.validated_case_records(tmp_path)  # a record of the superseded p-p selection rule
    write(7, _case_record(de5.CASES[7], off_curve={"pp_events_not_two_records": 3}))
    with pytest.raises(SystemExit, match="legacy selection schema"):
        de4.validated_case_records(tmp_path)
    write(7, _case_record(de5.CASES[7], cdf={}))
    with pytest.raises(SystemExit, match="no CDF hash"):
        de4.validated_case_records(tmp_path)  # a record from before the CDF lineage binding
    write(7, _case_record(de5.CASES[7]))
    (tmp_path / "de4-case7-copy.json").write_text("{}")
    with pytest.raises(SystemExit, match="unexpected file"):
        de4.validated_case_records(tmp_path)  # duplicate / stray record
    (tmp_path / "de4-case7-copy.json").unlink()
    write(10, _case_record(de5.CASES[0]))
    with pytest.raises(SystemExit, match="not a case of the frozen matrix"):
        de4.validated_case_records(tmp_path)


def test_signature_ratios_raise_on_missing_class_or_signature() -> None:
    from validation.scripts.elastic_attribution import de5_sampling as de5

    sigs = {s: 1.0 for s in de5.SIGNATURES}
    table = {"cases": {"3": {"classes": {"O-16": dict(sigs), "H-1": dict(sigs)}}}}
    rule = {"_selection_rule": de5.SELECTION_RULE}
    g4 = {
        "3": {
            "O-16": {s: 2.0 for s in de5.SIGNATURES},
            "H-1": {s: 2.0 for s in de5.SIGNATURES},
            **rule,
        }
    }
    assert de5.signature_ratios(table, g4)["3"]["O-16"]["T_tail_moment"] == pytest.approx(0.5)
    with pytest.raises(KeyError, match="no case 3"):
        de5.signature_ratios(table, {})
    with pytest.raises(KeyError, match="no class H-1"):
        de5.signature_ratios(table, {"3": {"O-16": g4["3"]["O-16"], **rule}})
    with pytest.raises(ValueError, match="selection rule"):  # a legacy Geant4 summary is refused
        de5.signature_ratios(
            table, {"3": {k: v for k, v in g4["3"].items() if k != "_selection_rule"}}
        )
    bad = {
        "3": {
            "O-16": g4["3"]["O-16"],
            "H-1": {k: v for k, v in g4["3"]["H-1"].items() if k != "qslope"},
            **rule,
        }
    }
    with pytest.raises(KeyError, match="no signature qslope"):
        de5.signature_ratios(table, bad)


@pytest.mark.requires_tables("elastic")
def test_de4_pp_event_with_a_second_interaction_is_accepted(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de4_parse as de4

    recs = _tracks(150.0, _mass("O-16"), 1, n_pa=6, n_pp=2, bad=False)
    pp_event = recs[-1][
        0
    ]  # the second p-p event; a further secondary proton from a second interaction in the target
    recs.append((pp_event, 3, 1, 20.0, 0.5, 0))
    d = _write_ps(tmp_path / "w", recs, 4_000_000, 150.0, "G4_WATER", 5.0)
    doc = de4.analyse([d])
    assert (
        doc["pp_events_all"] == 2 and doc["pp_events_off_curve"] == 0
    )  # the fast proton of the first pair is on the p+H curve
    assert doc["classes"]["H-1"]["events"] + doc["pp_events_removed_by_cm_cut"] == 2


# --- D-E6 swap tables and summary arithmetic (toy sources, synthetic batches) ---
def _toy_source(
    energies: tuple[float, ...], p_cm: tuple[float, ...], lo: float = -1.0, hi: float = 1.0
):
    from validation.scripts.elastic_attribution import de6_swap as d6

    rng = np.random.default_rng(5)
    mus = [
        np.sort(1.0 - 2.0 * rng.beta(1.0, 6.0 * (1.0 + 0.5 * i), 20000) ** 0.5)
        for i in range(len(energies))
    ]
    mus = [np.clip(m, lo, hi) for m in mus]
    q = np.array([d6._edges_q(m, pc) for m, pc in zip(mus, p_cm, strict=False)])
    return d6.ClassSource(np.array(energies), mus, np.array(p_cm), lo, hi, q)


def test_swap_edges_reproduce_the_source_cdf_and_interpolate_linearly_in_q() -> None:
    from validation.scripts.elastic_attribution import de6_swap as d6

    src = _toy_source((30.0, 50.0), (300.0, 400.0))
    for i, (e, pc) in enumerate(((30.0, 300.0), (50.0, 400.0))):
        edges = d6.swap_edges(
            e, src, pc
        )  # at a source energy with its own p_CM the sample edges come back
        np.testing.assert_allclose(edges, np.quantile(src.mu[i], np.arange(257) / 256), atol=1e-9)
        assert edges[0] == pytest.approx(src.mu[i].min()) and edges[-1] == pytest.approx(
            src.mu[i].max()
        )
    mid = d6.swap_edges(
        40.0, src, 350.0
    )  # intermediate energy: q at fixed level is the average of the two sources
    q_mid = 0.5 * (src.edges_q[0] + src.edges_q[1])
    np.testing.assert_allclose(mid, 1.0 - 2.0 * np.clip(q_mid / 700.0, 0, 1) ** 2, atol=1e-12)
    assert (
        np.all(np.diff(mid) > 0.0) and mid.size == 257
    )  # normalised: 256 equiprobable bins, strictly ascending


def test_swap_arrays_touch_only_the_swapped_class_nodes_and_keep_sigma() -> None:
    from validation.scripts.elastic_attribution import de6_swap as d6

    grid = np.array([10.0, 30.0, 60.0, 150.0, 250.0])
    names = ("H-1", "O-16")
    p_cm = np.array([[200.0] * 5, [250.0, 300.0, 350.0, 450.0, 500.0]])
    base = {"grid_e_mev": grid, "valid": np.ones((2, 5), dtype=np.int8), "p_cm_mev": p_cm,
            "sigma_barn": np.arange(10, dtype=float).reshape(2, 5) + 1.0,
            "edges_mu": np.tile(np.linspace(-1.0, 1.0, 257), (2, 5, 1))}  # fmt: skip
    src = {"O-16": _toy_source((30.0, 200.0), (300.0, 480.0))}
    arr, stats = d6.swap_arrays(base, names, src, ("O-16",))
    assert np.array_equal(arr["edges_mu"][0], base["edges_mu"][0])  # H-1 untouched
    assert np.array_equal(
        arr["edges_mu"][1, 0], base["edges_mu"][1, 0]
    )  # below 30 MeV: ionmc shape kept
    assert (
        not np.array_equal(arr["edges_mu"][1, 1:], base["edges_mu"][1, 1:])
        and stats["O-16"]["e_last_mev"] == 250.0
    )
    assert np.array_equal(arr["sigma_barn"], base["sigma_barn"])  # sigma unchanged
    sig = {"energies_mev": [30.0, 100.0], "factors": [0.8, 0.9]}
    arr2, _ = d6.swap_arrays(base, names, src, (), sig)
    assert arr2["sigma_barn"][0, 0] == base["sigma_barn"][0, 0]  # factor 1 below 30 MeV
    assert arr2["sigma_barn"][0, 1] == pytest.approx(0.8 * base["sigma_barn"][0, 1])
    assert arr2["sigma_barn"][0, 2] == pytest.approx(
        (0.8 + (60 - 30) / 70 * 0.1) * base["sigma_barn"][0, 2]
    )
    assert arr2["sigma_barn"][0, 4] == pytest.approx(
        0.9 * base["sigma_barn"][0, 4]
    ) and np.array_equal(arr2["sigma_barn"][1], base["sigma_barn"][1])


def test_de6_closure_interaction_and_extended_arithmetic() -> None:
    from validation.scripts.elastic_attribution import de6_swap as d6

    a = d6.arithmetic(
        {"pA": 0.020, "pp": 0.008, "both": 0.029},
        se_both=0.001,
        se_i=0.001,
        obs=0.035,
        se_obs=0.0026,
        se_cdf=0.001,
    )
    assert a["interaction"] == pytest.approx(0.001) and a["closure_residual"] == pytest.approx(
        0.006
    )
    assert (
        a["se_r"] == pytest.approx(math.sqrt(0.0026**2 + 0.001**2 + 0.001**2))
        and a["closure_threshold"] == 0.010
        and a["clause2_holds"]
    )
    assert a["clause3"]["named_component"] is None  # 0.020 / 0.029 = 0.69 < 70 %
    e = d6.extended_arithmetic(
        {"pA": 0.020, "pp": 0.008, "ppS": 0.006, "all": 0.036}, 0.001, 0.001, 0.035, 0.0026, 0.001
    )
    assert (
        e["interaction_all"] == pytest.approx(0.002)
        and e["closure_residual_all"] == pytest.approx(-0.001)
        and e["clause2_holds_extended"]
    )
    assert d6.outcome_label(False, True).startswith(
        "closure reached only with the post-hoc p-p sigma arm"
    )
    assert d6.outcome_label(True, True).startswith(
        "closure holds on the frozen"
    ) and d6.outcome_label(False, False).startswith("closure not reached")
    assert (
        d6.clause3_status(True, False) == "inconclusive (metric-dependent)"
        and d6.clause3_status(True, True) == "holds under both metrics"
    )


def test_de6_h4_recoil_share() -> None:
    from validation.scripts.elastic_attribution import de6_swap as d6

    nz = 200
    base = np.full((20, nz), 0.010)
    arm = np.full((20, nz), 0.004)
    h = d6.recoil_h4({"all": base, "swap_all": arm}, den_plateau=6.0, gain_excess=0.0058)
    assert h["swap_all"]["reduction_vs_base"] == pytest.approx((0.010 - 0.004) / 6.0)
    assert (
        h["swap_all"]["share_of_gain_excess"] == pytest.approx(0.001 / 0.0058)
        and not h["swap_all"]["ge_20pct_bound"]
    )


def _point_generated_paths_at_an_empty_dir(monkeypatch: pytest.MonkeyPatch, empty: Path) -> None:
    """Emulate a clean checkout: validation/generated and the transport-sample cache are git-ignored and absent."""
    empty.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(de_common, "GEN", empty)
    monkeypatch.setattr(de_common, "RUN_V11", empty / "16-v11-compare.txt")
    monkeypatch.setattr(de_common, "RUN_V10A", empty / "03-v10-a.txt")
    monkeypatch.setattr(de_common, "SAMPLES", empty / "samples")
    monkeypatch.setattr(
        de_common,
        "V11_PARTIALS",
        {k: empty / f"{k[0]}-{k[1]}.json" for k in de_common.V11_PARTIALS},
    )


@pytest.fixture(autouse=True)
def _clean_checkout_emulation(
    request: pytest.FixtureRequest,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path_factory: pytest.TempPathFactory,
):  # type: ignore[no-untyped-def]
    """``IONMC_TEST_NO_ARCHIVE=1`` runs the whole file with the generated/archive paths empty (CI proof)."""
    import os

    if os.environ.get("IONMC_TEST_NO_ARCHIVE") == "1":
        _point_generated_paths_at_an_empty_dir(monkeypatch, tmp_path_factory.mktemp("noarchive"))


@pytest.fixture(autouse=True)
def _restore_de6_pinned_constants():  # type: ignore[no-untyped-def]
    from validation.scripts.elastic_attribution import de6_swap as d6

    old = (d6.PINNED_ID, d6.PIN_PATH, dict(d6.ARCHIVED_V11_PIN))
    yield
    d6.PINNED_ID, d6.PIN_PATH = old[:2]
    d6.ARCHIVED_V11_PIN.clear()
    d6.ARCHIVED_V11_PIN.update(old[2])


def _fake_cdf_dir(d: Path) -> Path:
    """Five water CDF files of synthetic samples: O-16 plain, H-1 as the mirrored (+mu, -mu) pairs."""
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(3)
    for n, e in enumerate((30.0, 50.0, 100.0, 150.0, 200.0)):
        o = np.sort(np.cos(np.radians(rng.exponential(10.0 + n, 300) + 0.5)))
        phys = np.cos(np.radians(rng.uniform(16.26, 90.0, 120)))
        z = {"O-16_mu_sorted": o, "H-1_mu_sorted": np.sort(np.concatenate([phys, -phys]))}
        npz = d / f"de4-case{n}-cdf.npz"
        np.savez(npz, **z)
        side = {"case": n, "row": 14 + n, "material": "water", "energy_mev": e, "primaries": 4_000_000, "npz": npz.name,
                "npz_sha256": de_common.sha256_file(npz),
                "classes": {"O-16": {"p_cm_mev": 400.0 + 20 * n, "symmetrised": False},
                            "H-1": {"p_cm_mev": 250.0 + 10 * n, "symmetrised": True}}}  # fmt: skip
        (d / f"de4-case{n}-cdf.json").write_text(json.dumps(side))
    return d


def _toy_table(der: Path, idx: int, exploratory: dict, arrays: dict | None = None) -> str:
    """A real (tiny) table file pair whose id is the recomputable ``table_id`` of its sidecar."""
    from ionmc.nuclear.build import table_id

    tmp = der / "tmp.npz"
    np.savez(tmp, **(arrays if arrays is not None else {"x": np.full(4, float(idx))}))
    side = {
        "schema": "toy",
        "npz_sha256": de_common.sha256_file(tmp),
        "exploratory": exploratory,
    }
    side["table_id"] = table_id(side)
    tmp.rename(der / f"elastic-proton-{side['table_id']}.npz")
    (der / f"elastic-proton-{side['table_id']}.json").write_text(json.dumps(side))
    return side["table_id"]


# Archived V11 numbers as recorded in the committed H5 record (block ``archived_v11``); the fixture's toy h5.json carries them so that the
# D-E6 summary needs no local archive (validation/generated is git-ignored).
_ARCHIVED_ENERGIES = {
    "150": {"F_ionmc": 0.8567840337079885, "F_ionmc_se": 0.0015196783488711665, "F_topas": 0.891584942125469, "F_topas_se": 0.0021600982500448217,
            "gain_ionmc": 0.0220987766133649, "gain_topas": 0.016252044281034106, "F_relative_diff_se": 0.002885429972284705},
    "200": {"F_ionmc": 0.8440660646539855, "F_ionmc_se": 0.0025408116155115545, "F_topas": 0.8873703022321767, "F_topas_se": 0.0023837942663498973,
            "gain_ionmc": 0.02672217441583613, "gain_topas": 0.018419142339446415, "F_relative_diff_se": 0.0038376917410936247},
}  # fmt: skip


_ARCHIVED_V11 = {"source": {**{"run_id": "RUN-20261012T142957Z-a8e41216", "step_path": "lv5c-dd976c5-16/16-v11-compare.txt",
                              "evidence_sha": "dd976c5a69850b9597ac4da938792d4273a71719",
                              "compare_sha256": "1f81c774d5209761fce23f757e6f6c2d4f731c533be9b31d0e017a253175f6e2"},
                            "input_sha256": [{"path": "x", "sha256": "0" * 64}]},
                 "energies": _ARCHIVED_ENERGIES}  # fmt: skip


def _de6_fixture(
    tmp_path: Path,
    mutate: object = None,
    post: object = None,
    edit: object = None,
    tamper: object = None,
    pinned_sub: object = None,
    sidecar_tweak: object = None,
) -> dict:
    """Frozen D-E6 inputs on real files: CDF files, toy tables with recomputable ids, the lineage manifest written by ``write_lineage``,
    D-E4 case records, the CDF error record and run files with the real V5 geometry. ``mutate(ctx)`` edits the in-memory documents before they
    are written, ``post(ctx)`` tampers with the files afterwards."""
    from types import SimpleNamespace

    from validation.scripts.elastic_attribution import common as c
    from validation.scripts.elastic_attribution import de6_swap as d6

    rng = np.random.default_rng(11)
    cdf_dir = _fake_cdf_dir(tmp_path / "cdf")
    sources, prov = d6.load_sources(cdf_dir)
    prov_committed = prov  # the hashes the committed D-E4 records carry
    if (
        sidecar_tweak is not None
    ):  # a changed sidecar (e.g. p_CM): tables and manifest are rebuilt consistently from it
        sp = cdf_dir / "de4-case0-cdf.json"
        side = json.loads(sp.read_text())
        sidecar_tweak(side)
        sp.write_text(json.dumps(side))
        sources, prov = d6.load_sources(cdf_dir)
    scratch = tmp_path / "cache"
    (scratch / "derived").mkdir(parents=True)
    peak = {
        "emonly": 1.0,
        "all": 0.86,
        "swap-pA": 0.90,
        "swap-pp": 0.93,
        "swap-both": 0.97,
        "swap-pp-sigma": 0.89,
        "swap-all": 1.0,
    }
    ratios = (
        tmp_path / "de4-ratios.json"
    )  # the committed ratios record (sigma source of the p-p sigma arms)
    from validation.scripts.elastic_attribution import de5_sampling as de5

    ratios.write_text(json.dumps({"cases": {str(n): {cl: {"signatures": {sg: {"ratio": 1.0 + 0.01 * (n + 1), "ratio_se": 0.01} for sg in de5.SIGNATURES}}
                                                     for cl in (("O-16", "H-1") if n < 5 else ("C-12",))} for n in range(10)}}))  # fmt: skip
    sigma = d6.load_sigma_factors(ratios)
    cites = [
        {
            k: q[k]
            for k in ("case", "row", "energy_mev", "primaries", "npz_sha256", "sidecar_sha256")
        }
        for q in prov
    ]
    names = ("H-1", "O-16")
    grid = np.array([10.0, 30.0, 60.0, 150.0, 250.0])
    base = {"grid_e_mev": grid, "valid": np.ones((2, 5), dtype=np.int8), "sigma_barn": np.arange(10, dtype=float).reshape(2, 5) + 1.0,
            "p_cm_mev": np.array([[200.0] * 5, [250.0, 300.0, 350.0, 450.0, 500.0]]), "edges_mu": np.tile(np.linspace(-1.0, 1.0, 257), (2, 5, 1))}  # fmt: skip
    from ionmc.nuclear.build import table_id

    pinned = tmp_path / "pinned"  # the toy pinned table with its pin
    pinned.mkdir()

    def write_pinned(arrs: dict, pid: str | None) -> tuple[str, str]:
        tmp = pinned / "tmp.npz"
        np.savez(tmp, **arrs)
        sha = c.sha256_file(tmp)
        side = {"schema": "toy", "target_names": list(names), "npz_sha256": sha}
        tid = pid or table_id(side)
        side["table_id"] = tid
        tmp.replace(pinned / f"elastic-proton-{tid}.npz")
        (pinned / f"elastic-proton-{tid}.json").write_text(json.dumps(side))
        return tid, sha

    pid, sha0 = write_pinned(base, None)
    pin_path = tmp_path / "pin.json"
    pin_path.write_text(
        json.dumps({"table_id": pid, "npz_sha256": sha0})
    )  # the pin of the ORIGINAL pinned table
    if (
        pinned_sub is not None
    ):  # the pinned files are substituted consistently (recomputable sidecar), the pin stays
        base = pinned_sub(base)
        write_pinned(base, pid)
    d6.PINNED_ID, d6.PIN_PATH = pid, pin_path
    ids = {}
    for i, v in enumerate(d6.LINEAGE_VARIANTS):
        if v == "all":
            ex = {"d_e3_variant": "all", "source_table_id": pid, "parameters": {}}
        else:
            ex = {
                "d_e6_variant": v,
                "source_table_id": pid,
                "parameters": d6.expected_parameters(v, sigma),
                "cdf_sources": [dict(q) for q in cites],
            }
        if edit is not None:
            edit(v, ex)
        arrays = d6.build_variant_arrays(
            v, (base, names), sources, sigma
        )  # the builder of --build, from the toy pinned table
        if tamper is not None:
            tamper(v, arrays)
        ids[v] = _toy_table(scratch / "derived", i, ex, arrays)
    (scratch / "variants.json").write_text(
        json.dumps({"variants": {v: {"table_id": t} for v, t in ids.items()}})
    )
    out = tmp_path / "out"
    d6.write_lineage(scratch, prov, cdf_dir, out, ratios_path=ratios)
    manifest = out / d6.LINEAGE_FILE
    de4 = {n: {"case": n, "cdf": {"npz_sha256": prov_committed[n]["npz_sha256"], "sidecar_sha256": prov_committed[n]["sidecar_sha256"]},
               "classes": {"O-16": {"events": prov[n]["class_events"]["O-16"]}, "H-1": {"events": prov[n]["class_events"]["H-1"] // 2}}}
           for n in range(5)}  # fmt: skip
    terms = {
        "bootstrap": {"F_sd": 0.0002},
        "tail_edge": {"dF": 1e-4},
        "leave_one_energy_out_at_150": {"dF": 1e-4},
    }
    one = {**terms, "F_total_error": math.sqrt(0.0002**2 + 1e-4**2 + 1e-4**2)}
    err = {"150": one, "200": one, "n_boot": 200, "seed_base": d6.CDF_SEED, "seeds": {"150": d6.CDF_SEED + 150, "200": d6.CDF_SEED + 200},
           "n_z": 100, "gain_dz_mm": 2.0,
           "source_energies_mev": [30.0, 50.0, 100.0, 150.0, 200.0],
           "events_per_energy": {"O-16": [q["class_events"]["O-16"] for q in prov], "H-1": [q["class_events"]["H-1"] for q in prov]},
           "lineage": {"cdf_sources": [{"case": q["case"], "energy_mev": q["energy_mev"], "npz_sha256": q["npz_sha256"]} for q in prov]}}  # fmt: skip
    files: dict = {}
    for e in (150, 200):
        geo = d6.frozen_geometry(e)
        nz = geo["nz"]
        shape = np.full(nz, 5.0)
        shape[120:150] = np.linspace(5, 30, 30)
        shape[150:200] = np.linspace(30, 0.0, 50)
        shape[200:] = 0.0
        noise = 1.0 + 0.003 * rng.standard_normal((20, nz))  # shared batch noise (CRN)
        for v, pf in peak.items():
            idd = shape[None, :] * noise
            idd[:, 120:] *= pf
            idd[:, :120] *= (
                1.0
                if v == "emonly"
                else 1.025 - 0.012 * (v == "swap-all") - 0.006 * (v == "swap-pA")
            )
            rec = 0.01 * noise * (0.3 if v == "swap-all" else 1.0)
            files[(e, v)] = {
                "variant": v, "mode": "emonly" if v == "emonly" else "emel", "energy": e, "seed": 20990100 + (e == 200),
                "n": 100000, "n_batches": 20, "backend": "warp-cpu", "precision": "float64",
                "elastic_table_id": pid if v == "emonly" else ids[v], "bin_mm": geo["bin_mm"], "unit": geo["unit"],
                "scoring": geo["scoring"], "tallies": [] if v == "emonly" else ["recoil"], "valid": True,
                "relative_residual": 0.0, "idd_batches": idd.tolist(), "recoil_idd_batches": rec.tolist(),
                "recoil_total_mev_per_primary": 0.3,
            }  # fmt: skip
    ctx = SimpleNamespace(pid=pid, files=files, err=err, de4=de4, scratch=scratch, cdf_dir=cdf_dir, ids=ids, out=out, manifest=manifest,
                          runs=tmp_path / "runs", de4_dir=tmp_path / "de4", pinned=pinned, cdf_error=tmp_path / "cdf-error.json", prov=prov, tmp=tmp_path,
                          ratios=ratios, err_true=copy.deepcopy(err), h5=tmp_path / "h5.json", archive=tmp_path / "no-archive" / "16-v11-compare.txt",
                          archived=copy.deepcopy(_ARCHIVED_V11))  # fmt: skip
    if mutate is not None:
        mutate(ctx)
    ctx.runs.mkdir()
    ctx.de4_dir.mkdir()
    for (e, v), d in ctx.files.items():
        (ctx.runs / f"de3-{e}-{v}.json").write_text(json.dumps(d))
    for e in (150, 200):
        (ctx.runs / f"check-base-{e}.json").write_text(json.dumps({"bit_for_bit": True}))
    for n, r in ctx.de4.items():
        (ctx.de4_dir / f"de4-case{n}.json").write_text(json.dumps(r))
    ctx.cdf_error.write_text(json.dumps(ctx.err))
    ctx.h5.write_text(json.dumps({"archived_v11": ctx.archived}))
    if post is not None:
        post(ctx)
    return ctx


def _summarize(ctx):  # type: ignore[no-untyped-def]
    from validation.scripts.elastic_attribution import de6_swap as d6

    return d6.summarize(ctx.runs, ctx.out, ctx.cdf_error, ratios_path=ctx.ratios, scratch=ctx.scratch,
                        manifest_path=ctx.manifest, cdf_dir=ctx.cdf_dir, de4_dir=ctx.de4_dir,
                        recompute=lambda d: ctx.err_true, pinned_dir=ctx.pinned, h5_path=ctx.h5, archive_path=ctx.archive)  # fmt: skip


def test_de6_summary_end_to_end_on_real_files_of_the_frozen_geometry(tmp_path: Path) -> None:
    """Toy tables with recomputable ids, a real lineage manifest and the real V5 geometry (identity is covered; the physics is synthetic)."""
    from validation.scripts.elastic_attribution import de6_swap as d6

    ctx = _de6_fixture(tmp_path)
    doc = _summarize(ctx)
    b = doc["energies"]["150"]
    assert (
        b["one_minus_F"]["delta"]["both"] > b["one_minus_F"]["delta"]["pA"] > 0
        and "closure_residual_all" in b["one_minus_F_extended"]
    )
    assert b["outcome_label"] in {
        d6.outcome_label(f, x) for f in (True, False) for x in (True, False)
    }
    assert (
        set(b["h4"]) >= {"all", "swap_all", "ppS"} and b["h4"]["swap_all"]["reduction_vs_base"] > 0
    )
    proj = json.loads((ctx.out / "git" / "de6-summary.json").read_text())
    assert de_common.git_safety_violations(proj, d6.GIT_ALLOW, d6.GIT_ALLOW_LONG) == []
    md = (ctx.out / "git" / "de6-report.md").read_text()
    assert "Extended set" in md and "H4" in md and md.endswith("\n") and not md.endswith("\n\n")
    lin = json.loads(
        (ctx.out / "git" / "de6-lineage.json").read_text()
    )  # the manifest through the projection
    assert de_common.git_safety_violations(
        lin, d6.GIT_ALLOW_LINEAGE, d6.GIT_ALLOW_LONG_LINEAGE
    ) == [] and set(lin["variants"]) == set(d6.LINEAGE_VARIANTS)
    assert (
        len(lin["sigma_source"]["ratios_sha256"]) == 64
        and len(lin["sigma_source"]["sigma"]["factors"]) == 5
    )  # survives the projection


def _set(key: str, value: object, variant: str = "swap-pA", energy: int = 150):
    return lambda ctx: ctx.files[(energy, variant)].__setitem__(key, value)


def _del(key: str, variant: str = "all"):
    return lambda ctx: ctx.files[(150, variant)].pop(key)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (_set("seed", 20990101), "frozen"), (_set("n", 50000), "frozen"), (_set("n_batches", 10), "frozen"),
        (_set("elastic_table_id", "0" * 64), "frozen"), (_set("backend", "python"), "frozen"), (_set("precision", "float32"), "frozen"),
        (_set("variant", "swap-pp"), "frozen"), (_set("energy", 200), "frozen"), (_set("valid", False), "not valid"),
        (_set("relative_residual", 1e-6), "not valid"), (_set("idd_batches", [[1.0] * 150] * 19), "shape"),
        (_set("bin_mm", 1.0), "frozen"), (_set("unit", "MeV"), "frozen"), (_set("tallies", []), "frozen"),
        (_set("scoring", {"origin_mm": [0, 0, 0], "spacing_mm": [1, 1, 1], "shape": [1, 1, 5], "name": "dose"}), "frozen"),
        (_del("backend"), "frozen"), (_del("relative_residual"), "not valid"), (_del("scoring"), "frozen"), (_del("tallies"), "frozen"),
    ],
    ids=["seed", "n", "batches", "table-id", "backend", "precision", "variant", "energy", "valid", "residual", "shape", "bin", "unit", "tallies",
         "scoring", "no-backend", "no-residual", "no-scoring", "no-tallies"],
)  # fmt: skip
def test_de6_summary_rejects_run_files_not_bound_to_the_frozen_configuration(
    tmp_path: Path, mutation: object, message: str
) -> None:
    ctx = _de6_fixture(tmp_path, mutation)
    with pytest.raises(SystemExit, match=message):
        _summarize(ctx)


def _err(fn):  # type: ignore[no-untyped-def]
    return lambda ctx: fn(ctx.err)


@pytest.mark.parametrize(
    "mutation",
    [
        _err(lambda e: e.__setitem__("n_boot", 100)),
        _err(lambda e: e.__setitem__("seed_base", 1)),
        _err(lambda e: e["seeds"].__setitem__("150", 20990200)),  # the per-energy bootstrap seeds are bound
        _err(lambda e: e.__setitem__("source_energies_mev", [30.0, 50.0, 100.0, 150.0])),
        _err(lambda e: e.pop("lineage")),
        _err(lambda e: e["lineage"]["cdf_sources"][0].__setitem__("npz_sha256", "f" * 64)),
        _err(lambda e: e["150"].__setitem__("F_total_error", 0.5)),  # inflated error
        _err(lambda e: e["200"].__setitem__("F_total_error", 1e-9)),  # not the quadrature of its terms
        _err(lambda e: e["150"]["bootstrap"].__setitem__("F_sd", float("nan"))),  # non-finite term
        _err(lambda e: e["events_per_energy"]["O-16"].__setitem__(2, 1)),  # wrong event count
        _err(lambda e: e["events_per_energy"]["H-1"].__setitem__(0, e["events_per_energy"]["H-1"][0] // 2)),
    ],
    ids=["n-boot", "seed", "energy-seed", "energies", "no-lineage", "sha-mismatch", "inflated", "quadrature", "nan", "counts-O", "counts-H"],
)  # fmt: skip
def test_de6_cdf_error_gate_checks_numbers_not_labels(tmp_path: Path, mutation: object) -> None:
    ctx = _de6_fixture(tmp_path, mutation)
    with pytest.raises(SystemExit, match="CDF error record not bound"):
        _summarize(ctx)


def _tamper_file(path_of):  # type: ignore[no-untyped-def]
    def post(ctx):  # type: ignore[no-untyped-def]
        p = path_of(ctx)
        p.write_bytes(p.read_bytes() + b"x")

    return post


@pytest.mark.parametrize(
    ("post", "message"),
    [
        (_tamper_file(lambda c: c.scratch / "derived" / f"elastic-proton-{c.ids['swap-all']}.npz"), "differs from the manifest"),
        (lambda c: (c.scratch / "derived" / f"elastic-proton-{c.ids['swap-pA']}.npz").unlink(), "missing"),
        (_tamper_file(lambda c: c.cdf_dir / "de4-case3-cdf.npz"), "CDF file"),
        (lambda c: (c.cdf_dir / "de4-case1-cdf.npz").unlink(), "CDF file"),
        (lambda c: (c.de4_dir / "de4-case2.json").write_text(json.dumps({**c.de4[2], "cdf": {"npz_sha256": "0" * 64}})), "committed D-E4 record"),
        (lambda c: (c.de4_dir / "de4-case0.json").write_text(json.dumps({**c.de4[0], "classes": {"O-16": {"events": 1}, "H-1": {"events": 1}}})), "event counts"),
    ],
    ids=["scratch-npz-bytes", "scratch-npz-missing", "cdf-bytes", "cdf-missing", "de4-record-hash", "de4-record-counts"],
)  # fmt: skip
def test_de6_identity_binds_scratch_tables_cdf_files_and_de4_records(
    tmp_path: Path, post: object, message: str
) -> None:
    ctx = _de6_fixture(tmp_path, None, post)
    with pytest.raises(SystemExit, match=message):
        _summarize(ctx)


def test_load_sources_reads_the_symmetrised_flag_and_the_bootstrap_mirrors(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de6_swap as d6

    cdf_dir = _fake_cdf_dir(tmp_path / "cdf")
    src, _ = d6.load_sources(cdf_dir)
    assert src["H-1"].symmetrised is True and src["O-16"].symmetrised is False
    boot = d6._bootstrap(src, np.random.default_rng(4), ("O-16", "H-1"))
    for m, orig in zip(boot["H-1"].mu, src["H-1"].mu, strict=True):
        assert np.allclose(m, -m[::-1]) and m.size == orig.size  # mirrored
        assert (
            np.unique(np.abs(m)).size <= orig.size // 2
        )  # effective sample size: N physical events
    assert not np.allclose(boot["O-16"].mu[0], -boot["O-16"].mu[0][::-1])
    side = json.loads((cdf_dir / "de4-case0-cdf.json").read_text())
    side["classes"]["H-1"]["symmetrised"] = False
    (cdf_dir / "de4-case0-cdf.json").write_text(json.dumps(side))
    with pytest.raises(AssertionError, match="symmetrised"):
        d6.load_sources(cdf_dir)


def test_run_log_provenance_binds_executed_bytes_only_for_clean_runs(tmp_path: Path) -> None:
    import subprocess

    from validation.scripts.elastic_attribution import de3_probes as d3

    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=d3.c.REPO
    ).stdout.strip()
    root = tmp_path / "test-logs"
    runs = (("T1-de3-swap-pA-150", "de3-swap-pA-150", "", None), ("T2-de6c-swap-all-200", "de6c-swap-all-200", " M x.py", None),
            ("T3-de6b-swap-pp-150", "de6b-swap-pp-150", "", {"de3_probes.py": "ab" * 32}))  # fmt: skip
    for name, label, st, sha in runs:
        (root / name).mkdir(parents=True)
        meta = {
            "label": label,
            "commit": head,
            "git_status": st,
            "start_utc": "2026-10-15T11:00:00Z",
            "end_utc": "2026-10-15T11:01:00Z",
            "exit_code": 0,
        }
        if sha:
            meta["scripts_sha256"] = sha
        (root / name / "meta.json").write_text(json.dumps(meta))
    info = d3.read_logs([r[0] for r in runs], tmp_path)
    assert [(r["variant"], r["energy"], r["status_clean"]) for r in info["run_logs"]] == [
        ("swap-pA", 150, True),
        ("swap-all", 200, False),
        ("swap-pp", 150, True),
    ]
    clean, dirty, hashed = info["run_logs"]
    assert clean["binding"].startswith("bound: clean commit") and clean["script_blob"]
    assert (
        dirty["binding"] == "unbound (dirty working tree)" and dirty["script_blob"] is None
    )  # never a blob for a dirty run
    assert hashed["binding"].startswith("bound: scripts_sha256")
    lines = "\n".join(d3._run_provenance_lines(info))
    assert (
        "bound: clean commit" in lines
        and "unbound (dirty working tree)" in lines
        and "M x.py" in lines
        and head in lines
    )


def _edit(fn):  # type: ignore[no-untyped-def]
    return lambda v, ex: fn(ex) if v == "swap-pA" else None


@pytest.mark.parametrize(
    "edit",
    [
        lambda v, ex: ex.__setitem__("d_e6_variant", "swap-pp") if v == "swap-pA" else None,  # the table claims another variant (relabelled)
        _edit(lambda ex: ex.__setitem__("source_table_id", "0" * 64)),  # not derived from the pinned table
        _edit(lambda ex: ex["parameters"].__setitem__("classes", ["H-1"])),  # wrong classes
        _edit(lambda ex: ex["parameters"].__setitem__("below_mev_ionmc_shape", 20.0)),  # wrong floor
        _edit(lambda ex: ex["parameters"].__setitem__("end_edges", "-1 / +1")),  # wrong edge rule
        _edit(lambda ex: ex["parameters"].__setitem__("interpolation", "linear in E at fixed mu")),
        _edit(lambda ex: ex["cdf_sources"][0].__setitem__("npz_sha256", "e" * 64)),  # cites another CDF input
        lambda v, ex: ex["parameters"]["sigma"]["factors"].__setitem__(0, 2.0) if v == "swap-all" else None,  # wrong sigma source numbers
    ],
    ids=["relabelled", "not-pinned", "classes", "floor", "edges", "interpolation", "cdf-sources", "sigma"],
)  # fmt: skip
def test_de6_lineage_checks_table_semantics_not_only_ids(tmp_path: Path, edit: object) -> None:
    ctx = _de6_fixture(
        tmp_path, None, None, edit
    )  # real sidecars with a valid recomputable id but wrong semantics
    with pytest.raises(SystemExit, match="lineage not bound"):
        _summarize(ctx)


def test_de6_lineage_rejects_a_changed_ratios_record(tmp_path: Path) -> None:
    def post(ctx):  # type: ignore[no-untyped-def]
        doc = json.loads(ctx.ratios.read_text())
        doc["cases"]["3"]["H-1"]["signatures"]["rate_per_primary_after_rule"]["ratio"] = 1.5
        ctx.ratios.write_text(json.dumps(doc))

    ctx = _de6_fixture(tmp_path, None, post)
    with pytest.raises(SystemExit, match="sigma source"):
        _summarize(ctx)


def test_de6_cdf_error_is_recomputed_and_must_equal_the_record(tmp_path: Path) -> None:
    def both(ctx):  # type: ignore[no-untyped-def]
        x = ctx.err[
            "150"
        ]  # a consistent forgery: F_sd and the total change together (the quadrature still holds)
        x["bootstrap"]["F_sd"] = 0.0003
        x["F_total_error"] = math.sqrt(
            0.0003**2 + x["tail_edge"]["dF"] ** 2 + x["leave_one_energy_out_at_150"]["dF"] ** 2
        )

    ctx = _de6_fixture(tmp_path, both)
    with pytest.raises(SystemExit, match="not reproduced"):
        _summarize(ctx)  # the stub recomputation returns the true record


def test_error_recomputation_comparison_is_exact_to_1e12() -> None:
    from validation.scripts.elastic_attribution import de6_swap as d6

    rec = {"150": {"bootstrap": {"F_sd": 0.0002}, "F_total_error": 0.0003}, "n_boot": 200}
    d6.verify_error_recomputed(rec, copy.deepcopy(rec))
    off = copy.deepcopy(rec)
    off["150"]["bootstrap"]["F_sd"] += 1e-9
    with pytest.raises(SystemExit, match="not reproduced"):
        d6.verify_error_recomputed(rec, off)
    off2 = copy.deepcopy(rec)
    off2["n_boot"] = 100
    with pytest.raises(SystemExit, match="not reproduced"):
        d6.verify_error_recomputed(rec, off2)


@pytest.mark.parametrize(
    "tamper",
    [
        lambda v, a: a["edges_mu"].__setitem__((1, 2, 100), a["edges_mu"][1, 2, 100] + 1e-6) if v == "swap-pA" else None,  # one edge
        lambda v, a: a["sigma_barn"].__setitem__((0, 3), a["sigma_barn"][0, 3] * 1.01) if v == "swap-all" else None,  # sigma content
        lambda v, a: a["sigma_barn"].__setitem__((1, 1), a["sigma_barn"][1, 1] * 1.5) if v == "all" else None,  # the base copy altered
    ],
    ids=["edge", "sigma", "base-copy"],
)  # fmt: skip
def test_de6_lineage_rebuilds_table_contents(tmp_path: Path, tamper: object) -> None:
    """A table whose npz has different array content but consistent metadata (valid recomputable id, manifest written from it) is rejected."""
    ctx = _de6_fixture(tmp_path, tamper=tamper)
    with pytest.raises(SystemExit, match="table contents differ"):
        _summarize(ctx)


def test_build_variant_arrays_is_pure_and_matches_the_arms() -> None:
    from validation.scripts.elastic_attribution import de6_swap as d6

    src = {
        "O-16": _toy_source((30.0, 200.0), (300.0, 480.0)),
        "H-1": _toy_source((30.0, 200.0), (150.0, 250.0), -0.96, 0.96),
    }
    grid = np.array([10.0, 30.0, 60.0, 150.0, 250.0])
    base = {"grid_e_mev": grid, "valid": np.ones((2, 5), dtype=np.int8), "sigma_barn": np.ones((2, 5)),
            "p_cm_mev": np.array([[200.0] * 5, [250.0, 300.0, 350.0, 450.0, 500.0]]), "edges_mu": np.tile(np.linspace(-0.9, 0.9, 257), (2, 5, 1))}  # fmt: skip
    before = {k: v.copy() for k, v in base.items()}
    sig = {"energies_mev": [30.0, 200.0], "factors": [0.8, 0.9]}
    a = {
        v: d6.build_variant_arrays(v, (base, ("H-1", "O-16")), src, sig)
        for v in d6.LINEAGE_VARIANTS
    }
    assert all(
        np.array_equal(base[k], before[k]) for k in base
    )  # the pinned arrays are not modified
    assert all(np.array_equal(a["all"][k], base[k]) for k in base)
    assert np.array_equal(a["swap-pp-sigma"]["edges_mu"], base["edges_mu"]) and not np.array_equal(
        a["swap-pp-sigma"]["sigma_barn"], base["sigma_barn"]
    )
    assert not np.array_equal(a["swap-pA"]["edges_mu"][1], base["edges_mu"][1]) and np.array_equal(
        a["swap-pA"]["edges_mu"][0], base["edges_mu"][0]
    )
    assert np.array_equal(a["swap-all"]["edges_mu"], a["swap-both"]["edges_mu"]) and np.array_equal(
        a["swap-all"]["sigma_barn"], a["swap-pp-sigma"]["sigma_barn"]
    )


def test_de6_content_rebuild_authenticates_the_pinned_table_against_the_pin(tmp_path: Path) -> None:
    """The pinned arrays are substituted consistently (tables and manifest rebuilt from the substitute, valid recomputable sidecar): only the pin catches it."""

    def sub(base: dict) -> dict:
        b = {k: v.copy() for k, v in base.items()}
        b["sigma_barn"][1, 3] *= 1.01
        return b

    ctx = _de6_fixture(tmp_path, pinned_sub=sub)
    with pytest.raises(SystemExit, match="pinned table does not match the pin"):
        _summarize(ctx)


def test_de6_content_rebuild_binds_the_cdf_sidecar_to_the_committed_hash(tmp_path: Path) -> None:
    """A changed CDF sidecar field (p_CM supplies the transformation) with tables and manifest rebuilt consistently is rejected."""
    ctx = _de6_fixture(
        tmp_path, sidecar_tweak=lambda side: side["classes"]["O-16"].__setitem__("p_cm_mev", 777.0)
    )
    with pytest.raises(SystemExit, match="committed D-E4 record"):
        _summarize(ctx)


def test_de6_summary_needs_only_the_committed_h5_record_when_the_archive_is_absent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Clean checkout: every generated/archive path points at an empty directory; the summary runs on the committed H5 block alone."""
    _point_generated_paths_at_an_empty_dir(monkeypatch, tmp_path / "empty")
    ctx = _de6_fixture(tmp_path / "fx")
    ctx.archive = None  # the module default (c.RUN_V11) is now the missing path
    doc = _summarize(ctx)
    assert (
        doc["archived_v11"]["status"].startswith("archive absent")
        and len(doc["archived_v11"]["h5_sha256"]) == 64
    )
    assert doc["energies"]["150"]["base"]["F_ionmc_archived"] == pytest.approx(
        _ARCHIVED_ENERGIES["150"]["F_ionmc"], abs=1e-15
    )


def _fake_archive(ctx, path: Path) -> dict:  # type: ignore[no-untyped-def]
    arch = {e: {"F": {"ionmc": {"value": v["F_ionmc"], "se": v["F_ionmc_se"]}, "reference": {"value": v["F_topas"], "se": v["F_topas_se"]},
                      "se": v["F_relative_diff_se"]},
               "plateau_gain": {"ionmc": {"value": v["gain_ionmc"]}, "reference": {"value": v["gain_topas"]}}} for e, v in ctx.archived["energies"].items()}  # fmt: skip
    path.write_text(
        "#JSON-BEGIN\n"
        + json.dumps({"git_sha": "x", "verdict": {"energies": arch}})
        + "\n#JSON-END\n"
    )
    return arch


def test_de6_archived_v11_cross_check_against_a_present_archive(tmp_path: Path) -> None:
    import hashlib

    from validation.scripts.elastic_attribution import de6_swap as d6

    ctx = _de6_fixture(tmp_path)
    p = tmp_path / "arch.txt"
    arch = _fake_archive(ctx, p)
    with pytest.raises(
        SystemExit, match="frozen pin"
    ):  # a different archive file than the pinned one
        d6.archived_v11(ctx.h5, p)
    d6.ARCHIVED_V11_PIN["compare_sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
    ctx.archived["source"]["compare_sha256"] = d6.ARCHIVED_V11_PIN["compare_sha256"]
    ctx.h5.write_text(json.dumps({"archived_v11": ctx.archived}))
    got, prov = d6.archived_v11(ctx.h5, p)
    assert (
        prov["status"].startswith("archive cross-checked")
        and got[150]["F_ionmc"] == ctx.archived["energies"]["150"]["F_ionmc"]
    )
    arch["150"]["F"]["ionmc"]["value"] += 1e-9
    _ = p.write_text(
        "#JSON-BEGIN\n"
        + json.dumps({"git_sha": "x", "verdict": {"energies": arch}})
        + "\n#JSON-END\n"
    )
    d6.ARCHIVED_V11_PIN["compare_sha256"] = ctx.archived["source"]["compare_sha256"] = (
        hashlib.sha256(p.read_bytes()).hexdigest()
    )
    ctx.h5.write_text(json.dumps({"archived_v11": ctx.archived}))
    with pytest.raises(SystemExit, match="differs from the archive"):
        d6.archived_v11(ctx.h5, p)
    ctx.h5.write_text("{}")
    with pytest.raises(SystemExit, match="regenerate"):
        d6.archived_v11(ctx.h5, None)


def _tamper_h5(ctx, fn) -> None:  # type: ignore[no-untyped-def]
    blk = copy.deepcopy(ctx.archived)
    fn(blk)
    ctx.h5.write_text(json.dumps({"archived_v11": blk}))


def test_de6_archived_v11_block_is_validated_and_pinned(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import de6_swap as d6

    ctx = _de6_fixture(tmp_path)
    nan = float("nan")
    cases = {
        "nan": lambda b: b["energies"]["150"].update(F_relative_diff_se=nan),
        "inf": lambda b: b["energies"]["200"].update(gain_ionmc=float("inf")),
        "nonpositive se": lambda b: b["energies"]["150"].update(F_topas_se=0.0),
        "missing key": lambda b: b["energies"]["150"].pop("gain_topas"),
        "extra key": lambda b: b["energies"]["150"].update(extra=1.0),
        "missing energy": lambda b: b["energies"].pop("200"),
        "extra energy": lambda b: b["energies"].update({"100": b["energies"]["150"]}),
        "extra block key": lambda b: b.update(more=1),
        "extra source key": lambda b: b["source"].update(more="x"),
        "string number": lambda b: b["energies"]["150"].update(F_ionmc="0.85"),
    }
    for fn in cases.values():
        _tamper_h5(ctx, fn)
        with pytest.raises(SystemExit):
            d6.archived_v11(ctx.h5, None)

    # a consistent substitution of the values together with a recomputed summary is stopped by the source pin
    def sub(b):  # type: ignore[no-untyped-def]
        b["energies"]["150"]["F_ionmc"] += 0.01
        b["source"]["compare_sha256"] = "f" * 64

    _tamper_h5(ctx, sub)
    with pytest.raises(SystemExit, match="frozen pin"):
        _summarize(ctx)
    _tamper_h5(ctx, lambda b: b["source"].update(run_id="RUN-other"))
    with pytest.raises(SystemExit, match="frozen pin"):
        d6.archived_v11(ctx.h5, None)
