"""Exploratory D-E0 kernels (Amendment 19, report-only): analytic cases."""

from __future__ import annotations

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
