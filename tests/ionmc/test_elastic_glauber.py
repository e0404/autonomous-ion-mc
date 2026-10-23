"""V3-005D M2 (optical-limit Glauber) elastic shape builder: synthetic unit tests. The one
data-backed build test needs the source cache AND SAID text outputs (``IONMC_M2_SAID_DIR``) and
skips cleanly without them."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from ionmc.nuclear import elastic_build as EB
from ionmc.nuclear import elastic_glauber as G

SNIPPET = """<hr>
<PRE>

    SM16 0-3.0 GEV PP=50276/25348 NP=22832/12938 RAA (147)  5/03/07
     PP    SGT  at ACM =    0.00    10/17/26
     TLAB  Obs,Err
    10.000  0.3496E+03  0.0000E+00
    15.000  0.2252E+03  0.0000E+00
    20.000  0.1630E+03  0.0000E+00

 DELETED from analysis  10/13/86
      0.00  0.23730E+02 -0.50000E+00

</PRE>
"""


def _series(ch: str, ob: str, t: list[float], v: list[float]) -> G.SaidSeries:
    return G.SaidSeries(ch, ob, np.asarray(t), np.asarray(v), ("h",), ("0" * 64,))


def _const_nn(sig_mb: float = 25.72, al: float = 1.623) -> G.NNInputs:
    t = [1.0, 300.0]
    s = {("PP", "SGT"): _series("PP", "SGT", t, [sig_mb] * 2),
         ("NP", "SGT"): _series("NP", "SGT", t, [52.14] * 2),
         ("PP", "ALFA"): _series("PP", "ALFA", t, [al] * 2),
         ("NP", "ALFA"): _series("NP", "ALFA", t, [0.9435] * 2)}  # fmt: skip
    return G.NNInputs(s)


def test_said_reader_prediction_block_only(tmp_path: Path) -> None:
    h, ch, ob, t, v = G.parse_said_prediction(SNIPPET)
    assert h.startswith("SM16 0-3.0 GEV") and (ch, ob) == ("PP", "SGT")
    assert t.tolist() == [10.0, 15.0, 20.0] and v[0] == pytest.approx(349.6)
    p = tmp_path / "a.txt"
    p.write_text(SNIPPET)
    s = G.read_said_series([p])
    assert s.at(12.5) == pytest.approx(0.5 * (349.6 + 225.2)) and len(s.sha256[0]) == 64
    with pytest.raises(G.M2Error):
        s.at(5.0)  # no extrapolation


def test_said_reader_merges_files_and_rejects_conflicts(tmp_path: Path) -> None:
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    a.write_text(SNIPPET)
    b.write_text(SNIPPET.replace("10.000  0.3496E+03", "12.000  0.3000E+03"))
    s = G.read_said_series([a, b])
    assert s.t_mev.tolist() == [10.0, 12.0, 15.0, 20.0]
    b.write_text(SNIPPET.replace("0.3496E+03", "0.3500E+03"))
    with pytest.raises(G.M2Error):
        G.read_said_series([a, b])


def test_beta_formula_at_150_mev() -> None:
    nn = _const_nn()
    assert nn.beta_fm2("PP", 150.0) == pytest.approx(0.186, abs=5e-4)
    assert nn.beta_fm2("NP", 150.0) == pytest.approx(0.196, abs=5e-4)
    assert nn.sigma_fm2("PP", 150.0) == pytest.approx(2.572)


def test_ho_rms_formula_and_unfolding() -> None:
    d = G.HODensity(1.082, 1.692)
    r = np.linspace(0.0, 14.0, 140001)
    rho = d.density(r)
    assert np.trapezoid(4 * np.pi * r**2 * rho, r) == pytest.approx(1.0, abs=1e-9)
    rms2_num = np.trapezoid(4 * np.pi * r**4 * rho, r)
    assert rms2_num == pytest.approx(d.rms2(), rel=1e-9)
    assert d.rms2() == pytest.approx(1.692**2 * (6 + 15 * 1.082) / (4 + 6 * 1.082))
    pt = d.unfolded(G.R_PROTON_FM)
    assert pt.rms2() == pytest.approx(d.rms2() - G.R_PROTON_FM**2, rel=1e-12)
    q = np.linspace(0.0, 6.0, 25)
    np.testing.assert_allclose(
        pt.form_factor(q), d.form_factor(q) / np.exp(-(q**2) * G.R_PROTON_FM**2 / 6.0), rtol=1e-12
    )
    # the form factor is the Fourier transform of the density
    qq = 1.3
    ft = np.trapezoid(4 * np.pi * r**2 * rho * np.sin(qq * r) / (qq * r + 1e-300), r)
    assert ft == pytest.approx(float(d.form_factor(np.array([qq]))[0]), abs=1e-8)


def test_gaussian_folded_ho_thickness_against_brute_force() -> None:
    pt = G.HODensity(1.082, 1.692).unfolded()
    beta = 0.19
    # T(b) from the z integral of the density
    z = np.linspace(-14.0, 14.0, 5601)
    for b0 in (0.0, 1.7, 3.1):
        tb = np.trapezoid(pt.density(np.sqrt(b0**2 + z**2)), z)
        assert tb == pytest.approx(float(pt.gauss_folded_thickness(np.array([b0]), 0.0)[0]),
                                   rel=1e-9)  # fmt: skip
    # Cartesian 2D convolution with the Gaussian profile (independent quadrature)
    h = 0.04
    x = np.arange(-9.0, 9.0 + h / 2, h)
    xx, yy = np.meshgrid(x, x, indexing="ij")
    bp = np.sqrt(xx**2 + yy**2)
    zz = np.linspace(-12.0, 12.0, 961)
    tb_grid = np.trapezoid(pt.density(np.sqrt(bp[..., None] ** 2 + zz**2)), zz, axis=-1)
    for b0 in (0.0, 1.5, 3.0):
        f = np.exp(-((xx - b0) ** 2 + yy**2) / (2 * beta)) / (2 * np.pi * beta)
        num = float(np.sum(tb_grid * f) * h * h)
        ana = float(pt.gauss_folded_thickness(np.array([b0]), beta)[0])
        assert num == pytest.approx(ana, rel=1e-6)


def test_j0_accuracy() -> None:
    x = np.array([0.0, 0.4, 2.9, 3.1, 10.0, 47.3, 120.0])
    t = np.linspace(0.0, np.pi, 4001)[:-1]
    ref = np.array([np.mean(np.cos(xi * np.sin(t))) for xi in x])
    np.testing.assert_allclose(G.bessel_j0(x), ref, atol=2e-7)


def test_amplitude_black_disk_limit_and_optical_theorem() -> None:
    assert G.black_disk_limit_check()["max_rel_dev"] < 1e-3
    r_d, k = 4.0, 2.5
    b = np.arange(0.0, r_d + 1e-9, 0.005)
    chi = np.full(b.size, 60.0j)
    s = np.linspace(0.0, 1.0, 513)
    f = G.eikonal_amplitude(b, chi, k, s)
    x = 2 * k * s * r_d
    ana = np.where(
        x < 1e-9,
        k * r_d**2 / 2,
        k * r_d**2 * EB.bessel_j1(np.maximum(x, 1e-9)) / np.maximum(x, 1e-9),
    )
    np.testing.assert_allclose(f.imag, ana, atol=2e-4 * ana[0], rtol=0)
    assert abs(f.real).max() < 2e-4 * ana[0]
    res = G.m2_shape(b, chi, k, 2049)
    assert res.sigma_tot_ot_fm2 == pytest.approx(2 * math.pi * r_d**2, rel=2e-3)
    assert res.sigma_react_fm2 == pytest.approx(math.pi * r_d**2, rel=1e-3)


def test_m2_unitarity_and_normalisation_on_synthetic_nn() -> None:
    d = G.HODensity(1.082, 1.692)
    th = G.Thickness(d, G.b_grid(d))
    chi = G.eikonal_chi(th, _const_nn(), 150.0, 6, 6)
    k = 560.0 / G.HBARC
    res = G.m2_shape(th.b, chi, k, 2049)
    assert res.sigma_tot_ot_fm2 == pytest.approx(res.sigma_tot_b_fm2, rel=1e-9)
    assert res.sigma_el_fm2 + res.sigma_react_fm2 == pytest.approx(res.sigma_tot_ot_fm2, rel=1e-3)
    assert res.density.min() >= 0.0
    edges, dev = G.quantiles_mu(res.s, res.density, 256)
    assert dev < 1e-12 and edges[0] == -1.0 and edges[-1] == 1.0 and np.all(np.diff(edges) >= 0)


def test_quantile_convention_equals_pinned_disk_quantiles() -> None:
    r_fm, pcm, n_s, n_q = 3.4, 520.0, 4097, 256
    s = np.linspace(0.0, 1.0, n_s)
    a = EB.disk_amplitude_sq(2.0 * pcm / EB.HBARC * s * r_fm)
    e_new, _ = G.quantiles_mu(s, a, n_q)
    e_old, _ = EB.disk_quantiles(r_fm, pcm, n_s, n_q)
    np.testing.assert_allclose(e_new, e_old, atol=1e-13)
    assert G.HBARC == EB.HBARC and G.MP == EB.MP


def test_fermi_form_factor_small_q_and_unfolding_modes() -> None:
    f = G.FermiDensity(3.106, 0.542)
    q = np.array([1e-9, 0.1])
    ff = f.form_factor(q)
    assert ff[0] == pytest.approx(1.0, abs=1e-9)
    assert ff[1] == pytest.approx(1 - 0.01 * f.rms2() / 6, abs=3e-4)
    b = G.b_grid(f)
    th = G.Thickness(f, b, unfold="rms-matched")
    c = th.conv(0.3)
    assert np.trapezoid(2 * np.pi * b * c, b) == pytest.approx(
        1.0, abs=8e-3
    )  # cusp tail of the 2pF form factor
    assert th.point_fermi.rms2() == pytest.approx(f.rms2() - G.R_PROTON_FM**2, rel=1e-9)
    # plan-literal form-factor unfolding is exactly well behaved at large beta, diverges at small
    thf = G.Thickness(f, b, unfold="formfactor")
    thf.conv(0.7)
    assert thf.tail < 1e-20
    thf.conv(0.05)
    assert thf.tail > 1.0


def test_mott_endf_form_and_nonrelativistic_limit() -> None:
    from ionmc.nuclear import law5

    th = np.array([20.0, 45.0, 90.0, 135.0])
    m = G.mott_pp_endf_mb_sr(150.0, th, 0.99862)
    assert np.all(m > 0) and m[1] > m[2]
    np.testing.assert_allclose(G.mott_pp_endf_mb_sr(150.0, th, 0.99862)[1], m[1])
    # symmetric under theta -> 180 - theta (identical particles)
    np.testing.assert_allclose(G.mott_pp_endf_mb_sr(150.0, 180.0 - th, 0.99862), m, rtol=1e-12)
    # explicit (eta/2k)^2 [1/s^4 + 1/c^4 - cos(eta ln tan^2)/(s^2 c^2)] form
    k, eta = law5.coulomb_params(150e6, 0.99862, 0.99862, 1, 1)
    s2, c2 = np.sin(np.radians(th) / 2) ** 2, np.cos(np.radians(th) / 2) ** 2
    expl = (eta / (2 * k)) ** 2 * (1 / s2**2 + 1 / c2**2
                                   - np.cos(eta * np.log(s2 / c2)) / (s2 * c2)) * 1e3  # fmt: skip
    np.testing.assert_allclose(m, expl, rtol=1e-12)
    # relativistic Born (Moller) cross reference -> ENDF form at 90 deg, low energy
    low = G.mott_pp_moller_mb_sr(0.05, np.array([90.0]))[0]
    assert low == pytest.approx(G.mott_pp_endf_mb_sr(0.05, np.array([90.0]), 0.99862)[0], rel=2e-3)


def test_dsg_anchor_half_sphere_on_synthetic_scan() -> None:
    a = np.arange(1.0, 180.0)
    t, awi = 175.0, 0.99862
    const = 7.5  # mb/sr
    dsg = G.mott_pp_endf_mb_sr(t, a, awi) + const
    r = G.dsg_anchor(t, a, dsg, 24.0, awi)
    assert r["usable"] and r["sigma_anchor_mb"] == pytest.approx(
        2 * math.pi * const * 0.96, rel=1e-6
    )
    assert r["ratio_to_sgt"] == pytest.approx(r["sigma_anchor_mb"] / 24.0)
    bad = G.dsg_anchor(t, a, G.mott_pp_endf_mb_sr(t, a, awi) * 0.5, 24.0, awi)
    assert not bad["usable"]


def test_dsg_segments_merge() -> None:
    def seg(lo: int, hi: int) -> str:
        return "header\n" + "\n".join(
            f"{x:8.2f} {100.0 / x:.4E} 0.0 {x:8.2f} 1.0" for x in range(lo, hi + 1)
        )

    a, d = G.merge_dsg_segments([seg(1, 60), seg(61, 120), seg(121, 179)])
    assert a.size == 179 and d[0] == pytest.approx(100.0)
    with pytest.raises(G.M2Error):
        G.merge_dsg_segments([seg(1, 60), seg(61, 120)])


def _rows(ratios: dict[float, float], sgt: float = 0.9) -> list[dict[str, object]]:
    return [{"t_mev": t, "usable": True, "ratio_to_sgt": sgt, "ratio_to_table": r,
             "sigma_anchor_mb": 20.0 + 0.01 * t,
             "sigma_table_mb": 20.0} for t, r in ratios.items()]  # fmt: skip


def test_pp_reanchor_ramp_and_blocking_checks() -> None:
    ts = (100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0)
    ok = G.PPReanchor.build(_rows({t: 1.01 for t in ts}))  # type: ignore[arg-type]
    assert ok.enabled
    assert ok.sigma_mb(120.0, 21.0) == 21.0
    assert ok.sigma_mb(137.5, 20.0) == pytest.approx(20.0 * (1 + 0.01 * 0.5))
    assert ok.sigma_mb(200.0, 99.0) == pytest.approx(22.0)
    assert ok.sigma_mb(225.0, 99.0) == pytest.approx(22.0 + 0.5 * (22.5 - 22.0))
    bad = G.PPReanchor.build(_rows({t: (1.068 if t == 125.0 else 1.0) for t in ts}))  # type: ignore[arg-type]
    assert not bad.enabled and "continuity" in bad.reason
    assert bad.sigma_mb(200.0, 25.0) == 25.0
    unit = G.PPReanchor.build(_rows({t: 1.0 for t in ts}, sgt=2.0))  # type: ignore[arg-type]
    assert not unit.enabled and "SGT" in unit.reason


def test_provenance_json_and_densities() -> None:
    inp = G.load_inputs()
    assert inp["iaea_rms_charge_radius_fm"]["acquisition_id"].startswith("8dd77914")
    assert inp["iaea_rms_charge_radius_fm"]["tolerance_fm"] == 0.03
    for name in ("C-12", "O-16", "Si-28", "Ca-40"):
        d, ent = G.target_density(inp, name)
        rec = G.rms_check(d, name, inp)
        assert not ent["surrogate"] and rec["pass"], (name, rec)
    ref = inp["iaea_rms_charge_radius_fm"]["values"]
    for name in ("N-14", "Al-27", "P-31"):
        d, ent = G.target_density(inp, name)
        assert ent["surrogate"] and math.sqrt(d.rms2()) == pytest.approx(ref[name], abs=1e-9)
    n14, _ = G.target_density(inp, "N-14")
    assert n14.alpha == pytest.approx(1.313) and n14.a_fm == pytest.approx(1.7393, abs=2e-4)
    al, _ = G.target_density(inp, "Al-27")
    assert al.c_fm == pytest.approx(2.9751, abs=3e-3) and al.z_fm == 0.542  # numeric rms
    al_s, _ = G.target_density(inp, "Al-27", {"Al-27": 0.586})
    assert al_s.z_fm == 0.586 and al_s.c_fm < al.c_fm
    ca, _ = G.target_density(inp, "Ca-40")
    assert isinstance(ca, G.FermiDensity) and ca.w == -0.161


def test_horiuchi_table_file_reproduces_eq10_and_holds_constant() -> None:
    doc = json.loads(G.HORIUCHI_JSON.read_text())
    assert doc["acquisition_id"].startswith("f56c546a")
    e, s, a, b = (
        np.asarray(doc["table"][k]) for k in ("energy_mev", "sigma_fm2", "alpha", "beta_fm2")
    )
    sel = e <= 200.0
    np.testing.assert_allclose(b[sel], (1 + a[sel] ** 2) * s[sel] / (16 * math.pi), atol=5e-4)
    h = G.HoriuchiBeta.load()
    assert h.beta(10.0) == h.beta(30.0) == pytest.approx(0.685)
    assert h.beta(150.0) == pytest.approx(0.195)
    nn = _const_nn()
    nn.horiuchi = h
    assert nn.beta_fm2("PP", 150.0) == nn.beta_fm2("NP", 150.0) == pytest.approx(0.195)


def test_horiuchi_empty_file_fails_closed(tmp_path: Path) -> None:
    p = tmp_path / "h.json"
    p.write_text(json.dumps({"table": {"energy_mev": [], "beta_fm2": []}}))
    with pytest.raises(G.M2Error):
        G.HoriuchiBeta.load(p)


def test_m2_config_variants() -> None:
    assert G.M2Config(Path(".")).variant() == "V-b"
    assert G.M2Config(Path("."), beta="horiuchi-table").variant() == "V-c"
    with pytest.raises(G.M2Error):
        G.M2Config(Path("."), pp_level="said-anchored-ge150")  # needs the DSG scans


@pytest.mark.requires_tables("elastic")
def test_m2_build_into_tmp(tmp_path: Path) -> None:
    sdir = os.environ.get("IONMC_M2_SAID_DIR")
    if not sdir or not Path(sdir).is_dir():
        pytest.skip("IONMC_M2_SAID_DIR (SAID SGT/ALFA text outputs) not set")
    from ionmc.data import cache

    cfg = G.M2Config(Path(sdir), unfold="rms-matched", n_s=1025, allow_unpinned_inputs=True)
    res = EB.build_elastic_proton(cache.resolve_cache_dir(None), m2=cfg, out_dir=tmp_path)
    assert res.info["exploratory"]["model"] == "m2-optical-glauber"
    assert res.info["m2"]["pp_reanchor"] is None
    assert res.info["shape_normalisation_max_deviation"] <= 1e-12
    t2 = res.info["m2"]["trigger_2_convergence"]
    assert t2["passed"] and set(t2["component_max"]) == {
        "b_step",
        "s_grid",
        "radial_step",
        "q_step",
        "q_cutoff",
    }
    tg = {(c["target"]) for c in t2["candidate_energies"]}
    assert {"C-12", "O-16"} <= tg  # evaluated at the nodes nearest 100/150/200 MeV
    assert all(v is not None for v in t2["component_max"].values())  # Fermi targets present
    kinds = {x["name"]: x["kind"] for x in res.info["targets"]}
    assert kinds["C-12"] == "m2-glauber" and kinds["Ca-40"] == "m2-glauber"
    assert "model_omitted_reason" in res.info["xendf"]


def test_thickness_error_handling_and_model_pin_checks(tmp_path: Path) -> None:
    f = G.FermiDensity(3.106, 0.542)
    with pytest.raises(G.M2Error):
        G.Thickness(f, G.b_grid(f), unfold="nonsense")
    inp = G.load_inputs()
    inp["densities"]["N-14"]["status"] = "pending"
    with pytest.raises(G.M2Error):
        G.target_density(inp, "N-14")
    with pytest.raises(G.M2Error):
        G.target_density(inp, "S-32")
    # empty pins: refused unless explicitly allowed (synthetic SAID dir)
    for tag, ch, ob in (("a", "PP", "SGT"), ("b", "PP", "ALFA"), ("c", "NP", "SGT"),
                        ("d", "NP", "ALFA")):  # fmt: skip
        (tmp_path / f"{tag}.txt").write_text(
            f"HDR\n {ch}    {ob}  at ACM =    0.00\n     TLAB  Obs,Err\n"
            "    5.000  0.1000E+02  0.0\n  300.000  0.2000E+02  0.0\n\n"
        )
    empty = G.load_inputs()
    empty["said"]["file_sha256"], empty["pp_dsg"]["file_sha256"] = {}, {}
    ej = tmp_path / "empty_pins.json"
    ej.write_text(json.dumps(empty))
    with pytest.raises(G.M2Error, match="not pinned"):
        G.M2Model(G.M2Config(tmp_path, inputs_json=ej), [("C-12", 6, 12)])
    cfg = G.M2Config(tmp_path, allow_unpinned_inputs=True, inputs_json=ej)
    m = G.M2Model(cfg, [("C-12", 6, 12)])
    assert m.unpinned and m.provenance()["inputs_unpinned"] is True
    # a pinned hash list that does not contain the read files fails closed
    pj = tmp_path / "inp.json"
    doc = G.load_inputs()
    doc["said"]["file_sha256"] = {"x": "ab" * 32}
    pj.write_text(json.dumps(doc))
    with pytest.raises(G.M2Error, match="unpinned sha256"):
        G.M2Model(G.M2Config(tmp_path, inputs_json=pj), [("C-12", 6, 12)])


def test_surrogate_sensitivity_changes_size_not_rms() -> None:
    inp = G.load_inputs()
    rms = inp["iaea_rms_charge_radius_fm"]["values"]["N-14"]
    lo, _ = G.target_density(inp, "N-14", {"N-14": 1.082})
    hi, _ = G.target_density(inp, "N-14", {"N-14": 1.544})
    assert lo.a_fm == pytest.approx(1.7575, abs=3e-4) and hi.a_fm == pytest.approx(1.7254, abs=3e-4)
    assert math.sqrt(lo.rms2()) == pytest.approx(rms, abs=1e-9)
    assert math.sqrt(hi.rms2()) == pytest.approx(rms, abs=1e-9)


def test_pp_reanchor_disabled_leaves_level_and_records_reason() -> None:
    ts = (100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0)
    rea = G.PPReanchor.build(_rows({t: (1.068 if t == 125.0 else 1.0) for t in ts}))  # type: ignore[arg-type]
    assert not rea.enabled
    for e in (100.0, 130.0, 150.0, 220.0):
        assert rea.sigma_mb(e, 23.4) == 23.4
    assert "125" in rea.reason


def test_cli_exposes_m2_flags() -> None:
    from ionmc.cli import build_parser

    a = build_parser().parse_args(
        ["data", "build", "elastic-proton", "--model", "m2", "--beta", "horiuchi-table",
         "--pp-level", "said-anchored-ge150", "--unfold", "rms-matched", "--said-dir", "x",
         "--dsg-dir", "y", "--out-dir", "z", "--allow-unpinned-inputs"]
    )  # fmt: skip
    assert a.model == "m2" and a.allow_unpinned_inputs and a.pp_level == "said-anchored-ge150"
    assert G.M2Config(Path("."), pp_level="unchanged").pp_level == "none"


def test_helm_q_h_matches_verified_values() -> None:
    c = G.helm_q_h(G.HODensity(1.082, 1.692))
    o = G.helm_q_h(G.HODensity(1.544, 1.833))
    assert c["q_h_fm_inv"] == pytest.approx(1.912, abs=2e-3)
    assert o["q_h_fm_inv"] == pytest.approx(1.661, abs=2e-3)
    assert c["s_pt_fm"] == pytest.approx(0.825, abs=2e-3)
    assert c["r_h_fm"] == pytest.approx(2.350, abs=2e-3)


def _bdc(e: float, rp: float, rn: float) -> float:
    """Independent scalar transcription of the Bertulani & De Conti closed form."""
    rho0 = 0.17
    rt = (rp + rn) / rho0
    asym = 1.0 / (1.0 + 1.892 * (abs(rp - rn) / (rt * rho0)) ** 2.75)
    t = rt ** (2.0 / 3.0)
    base = 1.0 - 37.02 * t / e if e > 46.27 * t else e / (231.38 * t)
    return asym * base


def test_pauli_factor_closed_form_points() -> None:
    for e, rp, rn in ((100.0, 0.08, 0.08), (150.0, 0.04, 0.05), (30.0, 0.085, 0.085),
                      (23.0, 0.085, 0.085), (200.0, 0.01, 0.0)):  # fmt: skip
        got = float(G.pauli_factor(e, np.array([rp]), np.array([rn]), "bdc-pn")[0])
        assert got == pytest.approx(_bdc(e, rp, rn), rel=1e-12)
    # N = Z: asymmetry factor 1; branches meet at 0.2 (E = 46.27 rho~^(2/3))
    rt23 = (0.34 / 0.17) ** (2 / 3)
    lo = float(
        G.pauli_factor(46.27 * rt23 * (1 - 1e-9), np.array([0.17]), np.array([0.17]), "bdc-pn")[0]
    )
    hi = float(
        G.pauli_factor(46.27 * rt23 * (1 + 1e-9), np.array([0.17]), np.array([0.17]), "bdc-pn")[0]
    )
    assert lo == pytest.approx(0.2, abs=2e-4) and hi == pytest.approx(0.2, abs=2e-4)
    # rho -> 0 gives the free cross section
    z = np.array([1e-14])
    assert float(G.pauli_factor(150.0, z, z, "bdc-pn")[0]) == pytest.approx(1.0, abs=1e-6)
    # strong asymmetry reduces the factor: |N-Z|/A = 0.25 with the same total density
    sym = float(G.pauli_factor(150.0, np.array([0.04]), np.array([0.04]), "bdc-pn")[0])
    asy = float(G.pauli_factor(150.0, np.array([0.03]), np.array([0.05]), "bdc-pn")[0])
    assert asy < sym


def test_m2prime_reduces_to_m2_when_pauli_term_vanishes() -> None:
    nn = _const_nn()
    for dens in (G.HODensity(1.082, 1.692), G.FermiDensity(3.106, 0.542)):
        d_pt = dens
        b = G.b_grid(d_pt)
        th2 = G.Thickness(d_pt, b, unfold="rms-matched") if isinstance(
            d_pt, G.FermiDensity) else G.Thickness(d_pt, b)  # fmt: skip
        pt = G.PauliThickness(d_pt, b, 6, 6, variant="bdc-pn", pauli_scale=1e-12)
        chi_pt = G.eikonal_chi(pt, nn, 150.0, 6, 6) if not pt.set_energy(150.0) else None
        chi_m2 = G.eikonal_chi(th2, nn, 150.0, 6, 6)
        assert chi_pt is not None
        np.testing.assert_allclose(chi_pt, chi_m2, rtol=0, atol=3e-4 * np.abs(chi_m2).max())
    # the real Pauli factor lowers the phase (reduces sigma_tot)
    d = G.HODensity(1.082, 1.692)
    b = G.b_grid(d)
    full = G.PauliThickness(d, b, 6, 6, variant="bdc-pn")
    full.set_energy(100.0)
    assert 0.2 < full.pauli_peak < 0.9
    chi_full = G.eikonal_chi(full, nn, 100.0, 6, 6)
    chi_free = G.eikonal_chi(G.Thickness(d, b), nn, 100.0, 6, 6)
    assert np.abs(chi_full).max() < np.abs(chi_free).max()


def test_m2prime_config_and_measured_table(tmp_path: Path) -> None:
    kw = {"pauli": True, "pauli_variant": "bdc-pn", "pauli_emin_mev": 80.0}
    assert G.M2Config(Path("."), **kw).variant() == "M2prime"  # type: ignore[arg-type]
    cfg_c = G.M2Config(Path("."), beta="horiuchi-table", **kw)  # type: ignore[arg-type]
    assert cfg_c.variant() == "M2prime-c"
    assert G.M2Config(Path(".")).variant() == "V-b"  # m2 unchanged
    assert G.load_measured_sigma_tot(None) == {}
    p = tmp_path / "m.json"
    p.write_text(json.dumps({"C-12": {"100": 400.0, "500": 1.0}}))
    nodes = {"C-12": [{"e_mev": 50.0, "sigma_tot_m2_mb": 500.0},
                      {"e_mev": 150.0, "sigma_tot_m2_mb": 300.0}]}  # fmt: skip
    cmp = G.sigma_tot_comparison(G.load_measured_sigma_tot(p), nodes)
    assert cmp["C-12"]["100_mev"]["model_mb"] == pytest.approx(400.0)
    assert cmp["C-12"]["100_mev"]["model_over_measured"] == pytest.approx(1.0)
    assert cmp["C-12"]["500_mev"] is None
    f = G.FermiDensity(3.106, 0.542)
    inp = G.load_inputs()
    assert inp["densities"]["Si-28"]["form"] == "2pF" and f.density(np.array([0.0]))[0] > 0
    r = np.linspace(0, 12, 120001)
    assert np.trapezoid(4 * np.pi * r**2 * f.density(r), r) == pytest.approx(1.0, abs=1e-6)


def test_pauli_variants_required_and_clementel_villi() -> None:
    with pytest.raises(G.M2Error, match="pauli_variant"):
        G.M2Config(Path("."), pauli=True)  # no default for the reading
    with pytest.raises(G.M2Error, match="pauli_emin_mev"):
        G.M2Config(Path("."), pauli=True, pauli_variant="clementel-villi")  # no default emin
    with pytest.raises(TypeError):
        G.pauli_factor(100.0, np.array([0.08]), np.array([0.08]))  # type: ignore[call-arg]
    for bad in ("nonsense", "bdc-single"):  # bdc-single is deliberately not implemented
        with pytest.raises(G.M2Error):
            G.pauli_factor(100.0, np.array([0.08]), np.array([0.08]), bad)
    hb2m = 197.3269804**2 / (2 * 938.27208816)
    ef0 = hb2m * (1.5 * math.pi**2 * 0.17) ** (2 / 3)
    assert ef0 == pytest.approx(38.39, abs=0.01)
    assert float(G.fermi_energy_mev(np.array([0.17]))[0]) == pytest.approx(ef0)
    rho = np.array([0.17])
    # high branch: E >= 2 E_F
    e = 150.0
    got = float(G.pauli_factor(e, rho / 2, rho / 2, "clementel-villi")[0])
    assert got == pytest.approx(1 - 7 * ef0 / (5 * e), rel=1e-12)
    # value and continuity at E = 2 E_F: both branches give 0.3
    for eps in (1e-9, -1e-9):
        v = float(G.pauli_factor(2 * ef0 * (1 + eps), rho / 2, rho / 2, "clementel-villi")[0])
        assert v == pytest.approx(0.3, abs=1e-6)
    # low branch: E_F < E < 2 E_F
    e = 1.5 * ef0
    x = ef0 / e
    ref = 1 - 1.4 * x + 0.4 * x * (2 - 1 / x) ** 2.5
    assert float(G.pauli_factor(e, rho / 2, rho / 2, "clementel-villi")[0]) == pytest.approx(ref)
    assert 0.0 < ref < 0.3
    # blocked below the Fermi energy; free limit at low density; independent of p/n split
    assert float(G.pauli_factor(0.9 * ef0, rho / 2, rho / 2, "clementel-villi")[0]) == 0.0
    z = np.array([1e-14])
    assert float(G.pauli_factor(150.0, z, z, "clementel-villi")[0]) == pytest.approx(1.0, abs=1e-6)
    a = G.pauli_factor(150.0, np.array([0.03]), np.array([0.05]), "clementel-villi")
    b = G.pauli_factor(150.0, np.array([0.04]), np.array([0.04]), "clementel-villi")
    np.testing.assert_allclose(a, b)
    # monotone in density and energy
    dens = np.linspace(0.01, 0.17, 9)
    f = G.pauli_factor(150.0, dens / 2, dens / 2, "clementel-villi")
    assert np.all(np.diff(f) < 0)
    assert float(G.pauli_factor(200.0, rho / 2, rho / 2, "clementel-villi")[0]) > got


def test_m2prime_clementel_villi_against_researcher_scratch() -> None:
    sdir = os.environ.get("IONMC_M2_SAID_DIR")
    if not sdir or not Path(sdir).is_dir():
        pytest.skip("IONMC_M2_SAID_DIR (SAID SGT/ALFA text outputs) not set")
    ref = {  # researcher's scratch script (m2prime2.py, R-C): sigma_el [mb] at 100/150/200 MeV
        "C-12": (12, 6, (150.6, 122.0, 94.8)),
        "O-16": (16, 8, (209.6, 168.9, 131.2)),
    }
    cfg = G.M2Config(Path(sdir), pauli=True, pauli_variant="clementel-villi", pauli_emin_mev=80.0,
                     allow_unpinned_inputs=True)  # fmt: skip
    m = G.M2Model(cfg, [("C-12", 6, 12), ("O-16", 8, 16)])
    for name, (a, _, vals) in ref.items():
        mass = a * 931.494
        for e, want in zip((100.0, 150.0, 200.0), vals, strict=True):
            e1 = e + 938.272
            p = math.sqrt(e1 * e1 - 938.272**2)
            pcm = p * mass / math.sqrt(938.272**2 + mass**2 + 2 * e1 * mass)
            _, rec = m.node(name, e, pcm, 256)
            assert rec["sigma_el_m2_mb"] == pytest.approx(want, rel=0.01), (name, e)


def test_clementel_villi_table_values_and_pinned_convention() -> None:
    hb = 197.3269804
    mp = 938.27208816
    ef = hb**2 * (1.5 * math.pi**2 * 0.17) ** (2 / 3) / (2 * mp)
    assert ef == pytest.approx(38.392, abs=5e-4)
    r = np.array([0.17])
    for ratio, want in ((1.2, 0.0241), (1.5, 0.1138), (3.0, 0.5333), (6.0, 0.7667), (2.0, 0.3)):
        got = float(G.pauli_factor(ratio * ef, r / 2, r / 2, "clementel-villi")[0])
        assert got == pytest.approx(want, abs=6e-5), ratio
    assert float(G.pauli_factor(ef, r / 2, r / 2, "clementel-villi")[0]) == 0.0


@pytest.mark.requires_tables("elastic")
def test_m2prime_build_emin_split_into_tmp(tmp_path: Path) -> None:
    sdir = os.environ.get("IONMC_M2_SAID_DIR")
    if not sdir or not Path(sdir).is_dir():
        pytest.skip("IONMC_M2_SAID_DIR (SAID SGT/ALFA text outputs) not set")
    from ionmc.data import cache

    cfg = G.M2Config(Path(sdir), pauli=True, pauli_variant="clementel-villi",
                     pauli_emin_mev=80.0, n_s=1025, allow_unpinned_inputs=True)  # fmt: skip
    res = EB.build_elastic_proton(cache.resolve_cache_dir(None), m2=cfg, out_dir=tmp_path)
    pz = res.info["m2"]["pauli"]
    assert pz["pauli_variant"] == "clementel-villi" and pz["pauli_emin_mev"] == 80.0
    assert res.info["exploratory"]["pauli_emin_mev"] == 80.0
    nb = pz["nodes_below_emin_black_disk_per_target"]["C-12"]
    na = pz["nodes_applied_per_target"]["C-12"]
    assert nb > 0 and na > 0
    assert all(r["pauli_applied"] is False for r in pz["nodes_below_emin_black_disk"]["C-12"])
    assert min(r["e_mev"] for r in res.info["m2"]["targets"]["C-12"]["nodes"]) >= 80.0


def _synthetic_said_and_pins(tmp: Path) -> tuple[Path, Path]:
    d = tmp / "said_conv"
    d.mkdir()
    vals = {("PP", "SGT"): (31.0, 24.0), ("NP", "SGT"): (70.0, 40.0),
            ("PP", "ALFA"): (1.8, 1.2), ("NP", "ALFA"): (0.95, 0.75)}  # fmt: skip
    for i, ((ch, ob), (v0, v1)) in enumerate(vals.items()):
        rows = "\n".join(
            f"{t:10.3f}  {v0 + (v1 - v0) * (t - 5) / 295:.4E}  0.0" for t in (5.0, 100.0, 300.0)
        )
        (d / f"s{i}.txt").write_text(
            f"HDR\n {ch}    {ob}  at ACM =    0.00\n TLAB Obs,Err\n{rows}\n\n"
        )
    inp = G.load_inputs()
    inp["said"]["file_sha256"], inp["pp_dsg"]["file_sha256"] = {}, {}
    pj = tmp / "pins_empty.json"
    pj.write_text(json.dumps(inp))
    return d, pj


def test_convergence_components_reported_separately(tmp_path: Path) -> None:
    d, pj = _synthetic_said_and_pins(tmp_path)
    comps = ("b_step", "s_grid", "radial_step", "q_step", "q_cutoff")
    cases = [
        ("C-12", 6, 12, {}),  # HO: radial and q components do not exist (analytic fold)
        ("Si-28", 14, 28, {}),  # 2pF: every component exists
        (
            "C-12",
            6,
            12,
            {"pauli": True, "pauli_variant": "clementel-villi", "pauli_emin_mev": 80.0},
        ),  # fmt: skip
    ]
    for name, z, a, extra in cases:
        cfg = G.M2Config(d, inputs_json=pj, allow_unpinned_inputs=True, n_s=1025, **extra)
        m = G.M2Model(cfg, [(name, z, a)])
        mass, e = a * 931.494, 150.0
        e1 = e + 938.272
        p = math.sqrt(e1 * e1 - 938.272**2)
        pcm = p * mass / math.sqrt(938.272**2 + mass**2 + 2 * e1 * mass)
        cv = m.convergence(name, e, pcm)
        assert set(cv["components"]) == set(comps) and cv["tol"] == 1e-3
        analytic = name == "C-12" and not extra
        for c in comps:
            if analytic and c in ("radial_step", "q_step", "q_cutoff"):
                assert cv["components"][c] is None
            else:
                v = cv["components"][c]
                assert v is not None and v >= 0.0 and math.isfinite(v), (name, extra, c)
        vals = [v for v in cv["components"].values() if v is not None]
        assert cv["max"] == max(vals) and cv["passed"] == (cv["max"] <= 1e-3)
        assert cv["max"] <= 1e-3, (name, extra, cv)
