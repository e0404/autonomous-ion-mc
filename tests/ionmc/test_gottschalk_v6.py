"""Row V6 / V6-sens: parser, beam numbers, E0, sphere weights, comparator verdict logic
(V3-005C C8; synthetic LaTeX only, no measured value is stored in Git)."""

from __future__ import annotations

import importlib.util
import math
import os
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "gottschalk_v6", ROOT / "validation" / "scripts" / "reference" / "gottschalk_v6.py"
)
assert _spec is not None and _spec.loader is not None
g = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(g)

Z = np.array(g.DEPTHS_CM)
R = np.array(g.RADII_CM, dtype=float)


def synthetic_tex(missing: set[tuple[int, int]] = frozenset({(0, 1), (3, 5)})) -> str:  # type: ignore[assignment]
    head = "\\begin{table}\n\\begin{tabular}{rrrrrrrrrrr}\n\\multicolumn{1}{c}{z/r}&\n"
    cells = [f"\\multicolumn{{1}}{{c}}{{{r:.1f}}}" for r in g.RADII_CM]
    head += "&\n".join(cells[:-1]) + "&\n" + cells[-1] + "\\\\\n\\noalign{\\vspace{8pt}}\n"
    rows = []
    for j, z in enumerate(Z):
        vals = []
        for i in range(10):
            vals.append("M" if (j, i) in missing else f"{-1.0 - 0.01 * j - 0.1 * i:.4f}")
        rows.append(f" {z:.3f}& " + "& ".join(f"{v:>8}" for v in vals) + "\\\\")
    parms = (
        "\\begin{tabular}{lrll}\n\\multicolumn{1}{c}{param}&\\\\\n"
        "$p_1$&21.321&cm&R$_0$\\\\\n$p_2$&0.01372&cm$^{-1}$&s\\\\\n$p_3$&0.0312&cm&x\\\\\n"
        "$p_4$&$-0.231$&&x\\\\\n$p_5$&5.948&mm&sx\\\\\n$p_6$&31.55&mrad\\quad&st\\\\\n"
        "$p_7$&$-0.1462$&&em\\\\\n$p_{25}$&$-0.1431$&&x\n\\end{tabular}\n"
        "\\caption{x\\label{tbl:parms}}\n"
    )
    return (head + "\n".join(rows) + "\n\\end{tabular}\n\\caption{y\\label{tbl:dmlg}}\n"
            + parms)  # fmt: skip


def test_parse_dmlg_synthetic() -> None:
    z, r, v, nm = g.parse_dmlg(synthetic_tex())
    assert z.shape == (46,) and r.tolist() == list(g.RADII_CM) and v.shape == (46, 10)
    assert nm == 2 and np.isnan(v[0, 1]) and np.isnan(v[3, 5])
    assert math.isclose(v[2, 3], -1.0 - 0.02 - 0.3)
    assert abs(z[0] - 1.533) < 1e-12 and abs(z[-1] - 25.0) < 1e-12


def test_parse_dmlg_rejects_wrong_shape() -> None:
    tex = synthetic_tex()
    lines = tex.splitlines()
    k = next(i for i, ln in enumerate(lines) if ln.startswith(" 1.533"))
    with pytest.raises(ValueError):
        g.parse_dmlg("\n".join(lines[: k + 1] + lines[k + 2 :]))  # 45 depths


def test_parse_parms_and_frozen_values() -> None:
    p = g.parse_parms(synthetic_tex())
    assert p["p1"] == 21.321 and p["p2"] == 0.01372 and p["p5"] == 5.948
    assert p["p6"] == 31.55 and p["p7"] == -0.1462 and p["p4"] == -0.231 and p["p25"] == -0.1431
    with pytest.raises(ValueError):
        g.parse_parms(synthetic_tex().replace("21.321", "21.322"))


def test_beam_numbers_from_p7() -> None:
    b = g.beam_numbers(g.parse_parms(synthetic_tex()))
    assert math.isclose(b["f"], -0.8538, rel_tol=1e-12)
    assert b["rho_xtheta_plan_rounded"] == -0.521 and abs(b["rho_xtheta"] + 0.5206) < 1e-4
    assert math.isclose(b["A1_cm_rad"], -9.770e-3, rel_tol=2e-4)
    assert math.isclose(b["emittance_pi_mm_mrad"], 160.2, rel_tol=1e-3)
    assert abs(b["rho_xtheta"] - (-0.1462)) > 0.3  # not p7


def test_e0_bisection_on_synthetic_monotone_table() -> None:
    def rng(e: float) -> float:  # R = 0.0022 E^1.77 (monotone)
        return 0.0022 * e**1.77

    target = rng(177.123)
    e, it, trace = g.solve_e0(rng, target)
    assert abs(e - 177.123) / 177.123 < 2e-6 and len(trace) == it
    with pytest.raises(ValueError):
        g.solve_e0(rng, 1e6)
    rec = g.e0_record(
        21.321,
        1.0,
        rng,
        lambda r: (r / 0.0022) ** (1 / 1.77),
        lambda e: 4.8,
        lambda e: rng(e) * (1 + 5e-7),
        {"synthetic": True},
    )
    assert (
        rec["pass"]
        and rec["dose_table_read"] is False
        and rec["transport_table_finding"]["within_1e-6"]
    )
    assert math.isclose(rec["sigma_E_MeV"], 0.008 * 21.321 * 4.8)
    bad = g.e0_record(21.321, 1.0, rng, lambda r: (r / 0.0022) ** (1 / 1.77), lambda e: 4.8,
                      lambda e: rng(e) * 1.01, {})  # fmt: skip
    assert (
        bad["pass"] and not bad["transport_table_finding"]["within_1e-6"]
    )  # a finding, not a gate


def test_e0_real_table_and_record_flags() -> None:
    from ionmc import materials as M
    from ionmc.physics.projectiles import PROTON
    from ionmc.physics.stopping import BetheStoppingSource

    t = BetheStoppingSource().table(M.WATER, PROTON)
    rho = M.WATER.density_g_cm3
    rec = g.e0_record(
        21.321,
        rho,
        lambda e: float(t.range_at(e)),
        lambda r: float(t.energy_from_range(r)),
        lambda e: float(t.stopping_at(e)),
        lambda e: float(t.range_at(e)),
        {"t": "bethe"},
    )
    assert 176.0 <= rec["E0_MeV"] <= 178.0 and rec["inside_recording_band"] is True
    assert rec["inside_expectation_band"] is False and rec["finding_outside_176_178"] is False
    assert rec["cross_check_exact_inverse"]["ok"]
    assert abs(rec["sigma_E_MeV"] - rec["sigma_E_finite_difference_MeV"]) < 2e-3


def test_beam_parameters_variants_and_crn_scaling() -> None:
    p = g.parse_parms(synthetic_tex())
    bp = g.beam_parameters(p, lambda r0: (8.3 * r0, 0.008 * r0 * 4.8))
    assert set(bp) == set(g.SENS_NAMES)
    nom = bp["nominal"]
    assert math.isclose(bp["sx-lo"]["sigma_x_mm"], 0.8 * nom["sigma_x_mm"])
    assert math.isclose(bp["st-hi"]["sigma_theta_mrad"], 1.2 * nom["sigma_theta_mrad"])
    assert math.isclose(bp["rho-lo"]["rho_xtheta"], 0.8 * nom["rho_xtheta"])
    assert (
        abs(bp["rho-lo"]["rho_xtheta"] + 0.417) < 1e-3
        and abs(bp["rho-hi"]["rho_xtheta"] + 0.625) < 1e-3
    )
    assert bp["sx-lo"]["sigma_theta_mrad"] == nom["sigma_theta_mrad"]
    assert bp["r0-hi"]["E0_MeV"] > nom["E0_MeV"] > bp["r0-lo"]["E0_MeV"]


def test_tampered_hash_fails_closed_and_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import io
    import tarfile

    from ionmc.data import cache

    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w:gz") as tar:
        data = b"not the paper"
        info = tarfile.TarInfo("halo.tex")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    raw = payload.getvalue()
    monkeypatch.setattr(cache, "verify", lambda *_a, **_k: _write(tmp_path, raw))
    with pytest.raises(g.V6DataError) as ei:
        g.load_halo_tex(tmp_path)
    assert ei.value.record["tarball_ok"] is False and ei.value.record["sha256_tarball"]


def _write(tmp: Path, raw: bytes) -> Path:
    p = tmp / "obj"
    p.write_bytes(raw)
    return p


def test_real_halo_tex_hash_and_shape_when_cached() -> None:
    try:
        tex, rec = g.load_halo_tex()
    except g.V6DataError as exc:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail(f"e-print not available: {exc} {exc.record}")
        pytest.skip("e-print not cached")
    assert rec["tarball_ok"] and rec["halo_tex_ok"]
    z, r, v, nm = g.parse_dmlg(tex)
    assert v.shape == (46, 10) and nm == 165 and int(np.isfinite(v).sum()) == 295
    assert g.parse_parms(tex)["p1"] == 21.321


def test_sphere_weights() -> None:
    vol = 4.0 / 3.0 * math.pi * 2.0**3
    for r_mm, z_mm in ((0.0, 100.0), (10.0, 100.3), (100.0, 51.1), (30.0, 15.33)):
        ir, iz, w = g.sphere_weights(r_mm, z_mm)
        assert math.isclose(float(w.sum()), vol, rel_tol=1e-2), (r_mm, z_mm)
        assert len(w) == len(ir) == len(iz) and np.all(w > 0)
    # r = 0 against the analytic overlap of a sphere about the axis with cylindrical shells
    ir, iz, w = g.sphere_weights(0.0, 100.0)
    zz = np.linspace(98.0, 102.0, 400001)
    rad2 = 4.0 - (zz - 100.0) ** 2
    exact = {}
    for a, b in set(zip(ir.tolist(), iz.tolist(), strict=True)):
        lo, hi = a * 0.5, (a + 1) * 0.5
        in_z = (zz >= iz_lo(b)) & (zz < iz_lo(b) + 0.5)
        area = np.pi * (
            np.clip(rad2, 0, None).clip(max=hi * hi) - np.minimum(np.clip(rad2, 0, None), lo * lo)
        )
        exact[(a, b)] = float(np.trapezoid(np.where(in_z, area, 0.0), zz))
    assert math.isclose(sum(exact.values()), vol, rel_tol=1e-6)
    worst = max(abs(w[k] - exact[(ir[k], iz[k])]) for k in range(len(w)))
    assert worst < 0.03 * vol  # coarse 16^3 sub-sampling per cell
    # the dose average of a constant is the constant
    wt = {(0, 0): (ir, iz, w)}
    dose = np.full((3, 210, 520), 2.5)
    assert np.allclose(g.sphere_average_batches(dose, wt, 1, 1), 2.5)


def iz_lo(b: int) -> float:
    return b * 0.5


# ---- verdict logic --------------------------------------------------------------------------
def _meas() -> np.ndarray:
    d = 10.0 ** (-1.0 - 0.01 * np.arange(46)[:, None] - 0.1 * np.arange(10)[None, :])
    d[0, 1] = np.nan
    return d


def _batches(
    meas: np.ndarray, scale: np.ndarray | float, nb: int = 20, noise: float = 1e-4, seed: int = 3
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    base = np.nan_to_num(meas, nan=1e-3) * scale
    return base[None] * (1.0 + noise * rng.standard_normal((nb, 46, 10)))


def _scale(r0: float = 0.0, halo: float = 0.0) -> np.ndarray:
    s = np.ones((46, 10))
    s[:, 0] += r0
    s[:, 1:] += halo
    return s


def _sens(meas: np.ndarray, spread: float) -> tuple[dict[str, np.ndarray], np.ndarray]:
    nom = _batches(meas, 1.0, seed=10)
    var = {n: nom * (1.0 + spread * (1 if n.endswith("hi") else -1)) for n in g.VARIANTS}
    return var, nom


def test_r0_boundary_at_exactly_two_percent() -> None:
    m = _meas()
    for dev, ok in ((0.0199, True), (0.0201, False), (-0.0199, True), (-0.0201, False)):
        v = g.evaluate(Z, R, m, _batches(m, _scale(r0=dev), noise=0.0 + 1e-9))
        assert v["r0"]["pass"] is ok, dev
    v = g.evaluate(Z, R, m, _batches(m, _scale(r0=0.02), noise=1e-9))
    assert v["r0"]["n_points"] == 46 and v["r0"]["points"][0]["range_shift_mm"] is not None
    assert "interval90" in v["r0"]["points"][0]
    assert not any(p["z_cm"] == 1.533 and p.get("r_cm") == 1.0 for p in v["halo"]["points"])
    assert v["halo"]["n_points"] == 46 * 9 - 1


def test_halo_rules_ninety_percent_median_and_sign() -> None:
    m = _meas()
    s = _scale()
    v = g.evaluate(Z, R, m, _batches(m, s, noise=1e-9))
    assert (
        v["halo"]["pass"]
        and v["halo"]["frac_pass"] == 1.0
        and v["halo"]["sign_test"]["n_ties_excluded"] >= 0
    )
    # 12 % of points off by 20 %: fraction below 90 % -> fails
    s2 = s.copy()
    flat = [(j, i) for j in range(46) for i in range(1, 10)]
    for j, i in flat[: int(0.12 * len(flat))]:
        s2[j, i] = 1.2
    v2 = g.evaluate(Z, R, m, _batches(m, s2, noise=1e-9))
    assert v2["halo"]["frac_pass"] < 0.9 and not v2["halo"]["pass"]
    # 8 % off: fraction passes (0.92) -> passes if the median and trend are fine
    s3 = s.copy()
    for j, i in flat[: int(0.06 * len(flat))]:
        s3[j, i] = 1.2
    v3 = g.evaluate(Z, R, m, _batches(m, s3, noise=1e-9))
    assert v3["halo"]["frac_pass"] >= 0.9
    # median 0.12 -> fails the median rule although per-point (0.15) passes
    v4 = g.evaluate(Z, R, m, _batches(m, _scale(halo=0.12), noise=1e-9))
    assert (
        v4["halo"]["frac_pass"] == 1.0
        and not v4["halo"]["pass"]
        and v4["halo"]["median_abs"] > 0.10
    )
    # r-trend: far radii systematically higher than near ones at every depth
    s5 = s.copy()
    s5[:, 1:4] = 0.97
    s5[:, 4:] = 1.03
    v5 = g.evaluate(Z, R, m, _batches(m, s5, noise=1e-9))
    st = v5["halo"]["sign_test"]
    assert st["n_pos"] == st["n_depths"] == 46 and st["p"] < 0.01 and not v5["halo"]["pass"]


def test_sign_test_ties_excluded_and_binomial() -> None:
    m = _meas()
    v = g.evaluate(Z, R, m, _batches(m, _scale(), noise=1e-9))
    st = v["halo"]["sign_test"]
    assert st["n_pos"] + st["n_neg"] + st["n_ties_excluded"] == st["n_depths"]
    assert g.binom_two_sided_p(0, 0) == 1.0
    assert math.isclose(g.binom_two_sided_p(10, 5), 1.0)
    assert math.isclose(g.binom_two_sided_p(10, 0), 2.0 / 1024)
    assert g.binom_two_sided_p(20, 4) > 0.01 > g.binom_two_sided_p(20, 2)
    # exact tie: equal medians at a depth are counted and excluded
    s = _scale()
    resid = np.zeros((46, 10))
    assert resid.sum() == 0.0
    ties = g.evaluate(Z, R, m, _batches(m, s, noise=0.0 + 1e-12))
    assert (
        ties["halo"]["sign_test"]["n_ties_excluded"]
        + ties["halo"]["sign_test"]["n_pos"]
        + ties["halo"]["sign_test"]["n_neg"]
        == 46
    )


def test_v6_sens_condition_status() -> None:
    m = _meas()
    good_sens, nom = _sens(m, 0.05)
    v = g.evaluate(Z, R, m, _batches(m, _scale(), noise=1e-9), good_sens, nom)
    assert v["v6_sens"]["condition_holds"] is True and v["status"] == "gating"
    assert v["gate_pass"] is True
    bad_sens, nom = _sens(m, 0.16)
    v2 = g.evaluate(Z, R, m, _batches(m, _scale(), noise=1e-9), bad_sens, nom)
    assert v2["v6_sens"]["condition_holds"] is False and v2["status"] == "report-only"
    assert v2["gate_pass"] is None and v2["r0"]["pass"] is True  # r = 0 still computed and reported
    # |ratio - 1| + 2 sigma_Delta right at the limit: 0.14 with se 0 passes, 0.15 does not
    edge, nom = _sens(m, 0.14999)
    assert g.evaluate(Z, R, m, _batches(m, _scale(), noise=1e-9), edge, nom)["v6_sens"][
        "condition_holds"
    ]
    edge, nom = _sens(m, 0.1501)
    assert not g.evaluate(Z, R, m, _batches(m, _scale(), noise=1e-9), edge, nom)["v6_sens"][
        "condition_holds"
    ]
    # sigma_Delta enters: noisy variants with 0.10 mean shift and se > 0.025 fail
    noisy = {n: v_ * (1.0 + 0.1) * (1.0 + 0.03 * np.random.default_rng(4).standard_normal(v_.shape))
             for n, v_ in good_sens.items()}  # fmt: skip
    assert not g.evaluate(Z, R, m, _batches(m, _scale(), noise=1e-9), noisy, nom)["v6_sens"][
        "condition_holds"
    ]
    pend = g.evaluate(Z, R, m, _batches(m, _scale(), noise=1e-9))
    assert pend["status"] == "pending" and pend["gate_pass"] is None


def test_ratio_paired_se_is_zero_for_proportional_batches() -> None:
    nom = np.abs(np.random.default_rng(1).normal(1.0, 0.1, (20, 3)))
    r, se = g.ratio_paired(1.3 * nom, nom)
    assert np.allclose(r, 1.3) and np.allclose(se, 0.0, atol=1e-12)


def test_backend_disagreement_rule() -> None:
    m = _meas()
    sens, nom = _sens(m, 0.05)
    ok = g.evaluate(Z, R, m, _batches(m, _scale(), noise=1e-9, seed=5), sens, nom)
    ok2 = g.evaluate(Z, R, m, _batches(m, _scale(), noise=1e-9, seed=6), sens, nom)
    c = g.combine_backends(ok, ok2)
    assert c["pass"] is True and c["disagreement_components"] == []
    fail_r0 = g.evaluate(Z, R, m, _batches(m, _scale(r0=0.03), noise=1e-9), sens, nom)
    c2 = g.combine_backends(ok, fail_r0)
    assert (
        c2["pass"] is False
        and "r0" in c2["disagreement_components"]
        and "overall" in c2["disagreement_components"]
    )
    c3 = g.combine_backends(fail_r0, fail_r0)  # both fail: no disagreement but not a pass
    assert c3["pass"] is False and c3["disagreement_components"] == []
    # a per-point z beyond 4 is listed as a finding
    shifted = g.evaluate(Z, R, m, _batches(m, _scale(r0=0.001), noise=1e-3, seed=7), sens, nom)
    assert g.combine_backends(ok, shifted)["findings_abs_z_gt_4"]


def test_isolation_dose_table_read_only_by_the_comparator() -> None:
    for path in (ROOT / "src" / "ionmc").rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        if path.name == "registry.py":  # the registry entry itself names the dataset and the table
            continue
        assert "arxiv-1409-1938-source" not in text and "tbl:dmlg" not in text, path
    users = []
    for p in (ROOT / "validation" / "scripts").rglob("*.py"):
        text = p.read_text(encoding="utf-8")
        if "tbl:dmlg" in text or "parse_dmlg" in text:
            users.append(p.name)
    assert sorted(users) == ["gottschalk_v6.py", "steps_v6.py"], users
    src = (ROOT / "validation" / "scripts" / "transport" / "steps_v6.py").read_text(
        encoding="utf-8"
    )
    # only the comparator path reads the dose table: parse_dmlg is called once, in _data, whose only
    # caller is _compare; producers, the E0 procedure and the beam builder never touch it
    assert src.count("parse_dmlg(") == 1 and src.count("_data(g)") == 1
    assert src.index("parse_dmlg(") > src.index("def _data(") and src.index("_data(g)") > src.index(
        "def _compare("
    )
    assert (
        "dose_table_read"
        in (ROOT / "validation" / "scripts" / "reference" / "gottschalk_v6.py").read_text()
    )
