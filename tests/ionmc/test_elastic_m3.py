"""V3-005D M3 (ENDF/B-VIII.0 neutron MF4/MT2 mirror) shape: reader and builder tests on synthetic
ENDF-6 fixtures (no real evaluation is stored in Git)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from ionmc.data import endf6
from ionmc.nuclear import elastic_build as EB
from ionmc.nuclear import elastic_m3 as M3
from ionmc.nuclear.elastic_glauber import quantiles_mu


def ef(x: float) -> str:
    if x == 0.0:
        return "0.000000+0"
    m, e = f"{x:.6e}".split("e")
    return f"{m}{int(e):+d}".rjust(11)


def ei(n: int) -> str:
    return str(n).rjust(11)


class Writer:
    def __init__(self, mat: int) -> None:
        self.mat, self.out, self.ns = mat, [], 0

    def rec(self, mf: int, mt: int, fields: list[str]) -> None:
        self.ns += 1
        body = "".join(f.rjust(11) for f in fields).ljust(66)
        self.out.append(f"{body}{self.mat:4d}{mf:2d}{mt:3d}{self.ns:5d}")

    def text(self, mf: int, mt: int, txt: str) -> None:
        self.ns += 1
        self.out.append(f"{txt.ljust(66)[:66]}{self.mat:4d}{mf:2d}{mt:3d}{self.ns:5d}")

    def cont(self, mf: int, mt: int, c1: float, c2: float, l1: int, l2: int, n1: int, n2: int):
        self.rec(mf, mt, [ef(c1), ef(c2), ei(l1), ei(l2), ei(n1), ei(n2)])

    def floats(self, mf: int, mt: int, vals: list[float]) -> None:
        for i in range(0, len(vals), 6):
            self.rec(mf, mt, [ef(v) for v in vals[i : i + 6]])

    def ints(self, mf: int, mt: int, vals: list[int]) -> None:
        for i in range(0, len(vals), 6):
            self.rec(mf, mt, [ei(v) for v in vals[i : i + 6]])

    def send(self, mf: int) -> None:
        self.out.append(f"{'':66s}{self.mat:4d}{mf:2d}{0:3d}{99999:5d}")


def material(
    za: int, awr: float, mat: int, *, ltt: int = 1, lct: int = 2, li: int = 0,
    leg: list[tuple[float, list[float]]] | None = None,
    tab: list[tuple[float, list[tuple[float, float]]]] | None = None,
    nsub: int = 10, nver: int = 8, lt: int = 0,
) -> str:  # fmt: skip
    w = Writer(mat)
    w.cont(1, 451, za, awr, 1, 0, 0, 0)
    w.cont(1, 451, 0.0, 0.0, 0, 0, 0, 7)
    w.cont(1, 451, 1.0, 2.0e7, 0, 0, nsub, nver)
    w.cont(1, 451, 0.0, 0.0, 0, 0, 2, 0)
    w.text(1, 451, "synthetic fixture")
    w.text(1, 451, "ENDF/B-VIII.0 test")
    w.send(1)
    w.cont(4, 2, za, awr, 0, ltt, 0, 0)
    w.cont(4, 2, 0.0, awr, li, lct, 0, 0)
    if li == 0:
        if ltt in (1, 3):
            assert leg is not None
            w.cont(4, 2, 0.0, 0.0, 0, 0, 1, len(leg))
            w.ints(4, 2, [len(leg), 2])
            for e, a in leg:
                w.cont(4, 2, 0.0, e, lt, 0, len(a), 0)
                w.floats(4, 2, a)
        if ltt in (2, 3):
            assert tab is not None
            w.cont(4, 2, 0.0, 0.0, 0, 0, 1, len(tab))
            w.ints(4, 2, [len(tab), 2])
            for e, pts in tab:
                w.cont(4, 2, 0.0, e, lt, 0, 1, len(pts))
                w.ints(4, 2, [len(pts), 2])
                w.floats(4, 2, [v for p in pts for v in p])
    w.send(4)
    return "\n".join(w.out) + "\n"


def parsed(text: str) -> M3.Mf4Angular:
    return M3.parse_mf4(endf6.parse_endf(text).sections[(4, 2)])


LEG = [(1.0e6, [0.3, 0.1]), (3.0e6, [0.4, 0.15, 0.03]), (2.0e7, [0.5, 0.2, 0.05, 0.01])]


def test_legendre_reconstruction_normalised_and_values() -> None:
    m = parsed(material(8016, 15.857, 825, leg=LEG))
    mu = np.linspace(-1.0, 1.0, 20001)
    for e in (1.0e6, 3.0e6, 1.0e7, 2.0e7):
        f = m.f_native(mu, e)
        assert np.trapezoid(f, mu) == pytest.approx(1.0, abs=1e-6)
    f1 = m.f_native(np.array([0.5]), 1.0e6)[0]
    assert f1 == pytest.approx(0.5 + 1.5 * 0.3 * 0.5 + 2.5 * 0.1 * (1.5 * 0.25 - 0.5))
    assert m.emax_ev == 2.0e7 and m.ltt == 1 and m.lct == 2


def test_energy_interpolation_linlin_of_coefficients() -> None:
    m = parsed(material(8016, 15.857, 825, leg=LEG))
    mu = np.linspace(-1, 1, 11)
    mid = m.f_native(mu, 2.0e6)
    avg = 0.5 * (m.f_native(mu, 1.0e6) + m.f_native(mu, 3.0e6))
    np.testing.assert_allclose(mid, avg, atol=1e-12)
    with pytest.raises(M3.M3Error):
        m.f_native(mu, 3.0e7)


def test_tabulated_and_mixed_ltt() -> None:
    tab = [(1.0e7, [(-1.0, 0.2), (0.0, 0.5), (1.0, 0.8)]),
           (2.0e7, [(-1.0, 0.4), (0.0, 0.5), (1.0, 0.6)])]  # fmt: skip
    m = parsed(material(6012, 11.896, 625, ltt=2, tab=tab))
    assert m.f_native(np.array([0.5]), 1.0e7)[0] == pytest.approx(0.65)
    assert m.f_native(np.array([0.5]), 1.5e7)[0] == pytest.approx(0.5 * (0.65 + 0.55))
    m3 = parsed(material(6012, 11.896, 625, ltt=3, leg=LEG[:2], tab=tab))
    assert m3.f_native(np.array([0.0]), 2.0e6)[0] == pytest.approx(
        0.5 * (m3.f_native(np.array([0.0]), 1.0e6)[0] + m3.f_native(np.array([0.0]), 3.0e6)[0])
    )
    assert m3.f_native(np.array([1.0]), 1.0e7)[0] == pytest.approx(0.8)  # tabulated block


def test_isotropic_and_lab_frame_conversion() -> None:
    m = parsed(material(6012, 11.896, 625, li=1, ltt=0))
    assert m.f_native(np.array([0.3]), 1e6)[0] == 0.5
    lab = parsed(material(6012, 11.896, 625, lct=1, leg=[(1.0e6, [0.0]), (2.0e7, [0.0])]))
    mu = np.linspace(-1.0, 1.0, 40001)
    f = lab.f_cm(mu, 5.0e6)
    assert np.trapezoid(f, mu) == pytest.approx(1.0, abs=1e-6)  # isotropic in lab, normalised in CM
    a = 11.896
    s = math.sqrt(a * a + 2 * a * 0.4 + 1)
    assert lab.f_cm(np.array([0.4]), 5.0e6)[0] == pytest.approx(0.5 * a * a * (a + 0.4) / s**3)


def test_unsupported_variants_fail_closed() -> None:
    with pytest.raises(endf6.UnsupportedEndfError):
        parsed(material(8016, 15.857, 825, leg=LEG, lt=1))
    with pytest.raises(endf6.UnsupportedEndfError):
        parsed(material(8016, 15.857, 825, leg=LEG, lct=3))


def _write(d: Path, name: str, text: str) -> None:
    (d / name).write_text(text)


def _inputs(tmp: Path, pins: dict[str, str] | None = None) -> Path:
    doc = json.loads(M3.M3_INPUTS_JSON.read_text())
    doc["file_sha256"] = pins or {}
    p = tmp / "m3_inputs.json"
    p.write_text(json.dumps(doc))
    return p


def _fixture_dir(tmp: Path, name: str = "endf", leg: list | None = None) -> Path:  # type: ignore[type-arg]
    """C-12 and a genuinely different O-16 (x0.9 coefficients, mixed toward isotropic; x1.1 gave a
    negative density at 45 MeV that the positivity gate rejects) with shared records above 20 MeV
    (default LEGHI); pass ``leg=LEG`` for the <= 20 MeV files that must fail the copy check."""
    leg = LEGHI if leg is None else leg
    d = tmp / name
    d.mkdir()
    _write(d, "c12.endf", material(6012, 11.896, 625, leg=leg))
    leg_o = [(e, [x * 0.9 for x in a]) for e, a in leg]  # 0.9 f + 0.05: positive if f is
    _write(d, "o16.endf", material(8016, 15.857, 825, leg=leg_o))
    return d


def test_mf1_copy_check_and_unpinned_refusal(tmp_path: Path) -> None:
    d = _fixture_dir(tmp_path)
    inp = _inputs(tmp_path)
    with pytest.raises(M3.M3Error, match="pinned"):
        M3.M3Model(M3.M3Config(d, inputs_json=inp))
    model = M3.M3Model(M3.M3Config(d, inputs_json=inp, allow_unpinned_inputs=True))
    pv = model.provenance()
    assert pv["inputs_unpinned"] and pv["files"]["O-16"]["mat"] == 825
    assert pv["files"]["C-12"]["mf1_check"]["checks"]["nsub"] is True
    assert len(pv["files"]["O-16"]["sha256"]) == 64
    bad = tmp_path / "bad"
    bad.mkdir()
    _write(bad, "c12.endf", material(6012, 11.896, 626, leg=LEG))  # wrong MAT
    _write(bad, "o16.endf", material(8016, 15.857, 825, leg=LEG))
    with pytest.raises(M3.M3Error, match="catalogue"):
        M3.M3Model(M3.M3Config(bad, inputs_json=inp, allow_unpinned_inputs=True))
    nsub = tmp_path / "nsub"
    nsub.mkdir()
    _write(nsub, "c12.endf", material(6012, 11.896, 625, leg=LEG, nsub=0))
    _write(nsub, "o16.endf", material(8016, 15.857, 825, leg=LEG))
    with pytest.raises(M3.M3Error, match="catalogue"):
        M3.M3Model(M3.M3Config(nsub, inputs_json=inp, allow_unpinned_inputs=True))


LEGHI = [*LEG, (3.0e7, [0.55, 0.25, 0.06, 0.02]), (5.0e7, [0.6, 0.3, 0.08, 0.03])]


def _dir2(tmp: Path, name: str, c_text: str, o_text: str) -> Path:
    d = tmp / name
    d.mkdir()
    _write(d, "c12.endf", c_text)
    _write(d, "o16.endf", o_text)
    return d


def _model(d: Path, tmp: Path) -> M3.M3Model:
    return M3.M3Model(M3.M3Config(d, inputs_json=_inputs(tmp), allow_unpinned_inputs=True))


def test_copy_gate_distribution_level(tmp_path: Path) -> None:
    c = material(6012, 11.896, 625, leg=LEGHI)
    # identical Legendre data above 20 MeV: a copy
    with pytest.raises(M3.M3Error, match="copy"):
        _model(_dir2(tmp_path, "cp", c, material(8016, 15.857, 825, leg=LEGHI)), tmp_path)
    # a copy given in a different representation (tabulated 2001-point f of the C-12 record)
    cm = parsed(c)
    mu = np.linspace(-1.0, 1.0, 2001)
    tab = [(e, [(float(m), float(v)) for m, v in zip(mu, cm.f_cm(mu, e), strict=True)])
           for e in (3.0e7, 5.0e7)]  # fmt: skip
    o_tab = material(8016, 15.857, 825, ltt=3, leg=LEG, tab=tab)
    with pytest.raises(M3.M3Error, match="copy"):
        _model(_dir2(tmp_path, "cp2", c, o_tab), tmp_path)
    # a genuinely different O-16 passes and records the evidence
    leg_o = [(e, [x * 0.9 for x in a]) for e, a in LEGHI]
    ok = _model(_dir2(tmp_path, "ok", c, material(8016, 15.857, 825, leg=leg_o)), tmp_path)
    cc = ok.copy_check
    assert cc["evaluated"] and cc["n_mu_points"] == 181 and cc["shared_energies"] == 2
    assert cc["energies_that_differ"] == 2
    # differs at the points of only some energies -> still a copy when <= 50 % of the energies
    mixed = [(e, a if i < 4 else [x * 1.1 for x in a]) for i, (e, a) in enumerate(LEGHI)]
    with pytest.raises(M3.M3Error, match="copy"):
        _model(_dir2(tmp_path, "mx", c, material(8016, 15.857, 825, leg=mixed)), tmp_path)
    # no shared energy above 20 MeV: the criterion cannot be established -> fail closed
    with pytest.raises(M3.M3Error, match="share no record energy above 20 MeV"):
        _model(_fixture_dir(tmp_path, "low", LEG), tmp_path)
    # O-16 file missing: also fail closed
    only_c = tmp_path / "only_c"
    only_c.mkdir()
    _write(only_c, "c12.endf", c)
    with pytest.raises(M3.M3Error):
        _model(only_c, tmp_path)


def test_normalisation_gate(tmp_path: Path) -> None:
    mu_c = [(-1.0, 0.2), (0.0, 0.5), (1.0, 0.8)]  # integral 1
    mu_o = [(-1.0, 0.3), (0.0, 0.5), (1.0, 0.7)]  # integral 1, differs from C-12
    energies = (1.0e6, 2.0e7, 3.0e7, 5.0e7)  # records above 20 MeV are shared with C-12
    ok = [(e, mu_o) for e in energies]
    bad = [(1.0e6, [(m, v * 1.0001) for m, v in mu_o]), *[(e, mu_o) for e in energies[1:]]]
    c = material(6012, 11.896, 625, ltt=2, tab=[(e, mu_c) for e in energies])
    good = _model(_dir2(tmp_path, "n_ok", c, material(8016, 15.857, 825, ltt=2, tab=ok)), tmp_path)
    good.node("C-12", 5.0, 256)
    m = _model(_dir2(tmp_path, "n_bad", c, material(8016, 15.857, 825, ltt=2, tab=bad)), tmp_path)
    with pytest.raises(M3.M3Error, match="normalisation"):
        m.node("O-16", 1.0, 256)
    m.node("O-16", 20.0, 256)  # the other record is normalised


def test_positivity_gates_no_clipping(tmp_path: Path) -> None:
    neg = [(1.0e6, [0.3, 0.1]), (3.0e6, [0.95, 0.9, 0.7]), *[(e, [x * 1.1 for x in a])
                                                            for e, a in LEGHI[2:]]]  # fmt: skip
    c = material(6012, 11.896, 625, leg=LEGHI)
    m = _model(_dir2(tmp_path, "neg", c, material(8016, 15.857, 825, leg=neg)), tmp_path)
    m.node("O-16", 1.0, 256)  # positive record
    with pytest.raises(M3.M3Error, match="negative"):
        m.node("O-16", 3.0, 256)  # negative record: refused, not clipped
    with pytest.raises(M3.M3Error, match="negative"):
        m.check_midpoint("O-16", 2.5)  # a midpoint between a positive and a negative record
    m.check_midpoint("O-16", 1.1)  # positive midpoint passes
    ok = _model(_fixture_dir(tmp_path), tmp_path)
    for e in (1.5, 2.0, 7.5, 12.5):
        ok.check_midpoint("C-12", e)


def test_moment_and_cdf_gate_records(tmp_path: Path) -> None:
    m = _model(_fixture_dir(tmp_path), tmp_path)
    _, rec = m.node("C-12", 2.0, 256)
    assert rec["cdf_monotone"] is True and 0.0 < rec["first_moment_1_minus_mu"] < 2.0
    iso = parsed(material(6012, 11.896, 625, li=1, ltt=0))
    m.mf4["C-12"] = iso
    _, rec = m.node("C-12", 2.0, 256)
    assert rec["first_moment_1_minus_mu"] == pytest.approx(1.0, abs=1e-6)  # isotropic <1 - mu>


def test_hold_last_is_fixed_in_momentum_transfer(tmp_path: Path) -> None:
    m = _model(_fixture_dir(tmp_path), tmp_path)
    p_top, p_node = 400.0, 520.0  # CM momenta at E_top and at the node [MeV/c]
    e_top, _ = m.node("C-12", 50.0, 256)
    with pytest.raises(M3.M3Error):
        m.node("C-12", 100.0, 256)  # needs the CM momenta
    with pytest.raises(M3.M3Error):
        m.node("C-12", 100.0, 256, 300.0, 400.0)  # node momentum below p_top
    e_node, rec = m.node("C-12", 100.0, 256, p_node, p_top)
    assert rec["q_scaled"] is True and rec["e_used_mev"] == pytest.approx(50.0)
    ratio = p_top / p_node
    assert rec["p_ratio_top_over_node"] == pytest.approx(ratio)
    # identity: q = 2 p sin(theta/2) of every edge is unchanged
    q_top = 2.0 * p_top * np.sqrt((1.0 - e_top) / 2.0)
    q_node = 2.0 * p_node * np.sqrt(np.maximum((1.0 - e_node) / 2.0, 0.0))
    np.testing.assert_allclose(q_node, q_top, rtol=1e-12, atol=1e-9)
    # support ends at the E_top kinematic limit: zero density below mu_min, edges monotone
    assert e_node[0] == pytest.approx(1.0 - 2.0 * ratio**2) and e_node[-1] == 1.0
    assert e_node[0] > -1.0 and np.all(np.diff(e_node) >= 0.0)
    # p_node == p_top reproduces the E_top edges exactly
    same, _ = m.node("C-12", 100.0, 256, p_top, p_top)
    np.testing.assert_allclose(same, e_top, atol=1e-12)


def test_pinned_hash_must_match(tmp_path: Path) -> None:
    d = _fixture_dir(tmp_path)
    with pytest.raises(M3.M3Error, match="not pinned"):
        M3.M3Model(M3.M3Config(d, inputs_json=_inputs(tmp_path, {"x": "ab" * 32})))


def test_node_edges_equiprobable_nonnegative_and_conventions(tmp_path: Path) -> None:
    d = _fixture_dir(tmp_path)
    model = M3.M3Model(M3.M3Config(d, inputs_json=_inputs(tmp_path), allow_unpinned_inputs=True))
    edges, rec = model.node("C-12", 5.0, 256)
    assert edges.shape == (257,) and edges[0] == -1.0 and edges[-1] == 1.0
    assert np.all(np.diff(edges) >= 0.0)
    assert rec["norm_int_f_dmu"] == pytest.approx(1.0, abs=1e-6) and rec["norm_deviation"] < 1e-12
    # equiprobable under f: CDF at the edges is i/256
    mu = np.linspace(-1, 1, 200001)
    f = model.mf4["C-12"].f_cm(mu, 5.0e6)
    cdf = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(mu))))
    np.testing.assert_allclose(np.interp(edges, mu, cdf), np.arange(257) / 256, atol=5e-5)
    # same machinery as the pinned convention: isotropic input gives uniform edges in mu
    s = np.linspace(0.0, 1.0, 4097)
    e_iso, _ = quantiles_mu(s, np.ones_like(s), 256)
    np.testing.assert_allclose(e_iso, np.linspace(-1, 1, 257), atol=1e-6)
    assert EB.N_Q == 256


def test_above_emax_policy_validation() -> None:
    with pytest.raises(M3.M3Error):
        M3.M3Config(Path("."), above_emax="extrapolate")


@pytest.mark.requires_tables("elastic")
def test_m3_build_into_tmp(tmp_path: Path) -> None:
    from ionmc.data import cache

    d = _fixture_dir(tmp_path)
    cfg = M3.M3Config(d, above_emax="hold-last", inputs_json=_inputs(tmp_path),
                      allow_unpinned_inputs=True)  # fmt: skip
    out = tmp_path / "out"
    res = EB.build_elastic_proton(cache.resolve_cache_dir(None), m3=cfg, out_dir=out)
    ex = res.info["exploratory"]
    assert ex["model"] == "m3-endf-mirror" and ex["inputs_unpinned"] is True
    t = res.info["m3"]["targets"]
    assert t["O-16"]["shape"] == "endf-mf4-mt2" and t["O-16"]["n_nodes_above_emax"] > 0
    assert t["O-16"]["n_nodes_q_scaled"] == t["O-16"]["n_nodes_above_emax"]
    assert t["O-16"]["n_midpoints_checked_positive"] > 0 and t["O-16"]["coverage_limitation"]
    assert t["Ca-40"]["shape"].startswith("black-disk")
    kinds = {x["name"]: x["kind"] for x in res.info["targets"]}
    assert kinds["C-12"] == kinds["O-16"] == "m3-endf-mirror" and kinds["Ca-40"] == "bgg-black-disk"
    assert "model_omitted_reason" in res.info["xendf"]
    with np.load(res.npz_path) as z:
        assert np.all(np.isnan(z["xendf_model_mb_sr"]))
    with pytest.raises(EB.BuildError, match="above the ENDF evaluation range"):
        EB.build_elastic_proton(
            cache.resolve_cache_dir(None),
            m3=M3.M3Config(d, inputs_json=_inputs(tmp_path), allow_unpinned_inputs=True),
            out_dir=out,
        )
