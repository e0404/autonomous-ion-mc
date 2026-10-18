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


def _fixture_dir(tmp: Path) -> Path:
    d = tmp / "endf"
    d.mkdir()
    _write(d, "c12.endf", material(6012, 11.896, 625, leg=LEG))
    leg_o = [(e, [x * 1.1 for x in a]) for e, a in LEG]
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


def test_o16_copy_of_c12_is_refused(tmp_path: Path) -> None:
    d = tmp_path / "copy"
    d.mkdir()
    _write(d, "c12.endf", material(6012, 11.896, 625, leg=LEG))
    _write(d, "o16.endf", material(8016, 15.857, 825, leg=LEG))
    with pytest.raises(M3.M3Error, match="copy"):
        M3.M3Model(M3.M3Config(d, inputs_json=_inputs(tmp_path), allow_unpinned_inputs=True))


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
    cfg = M3.M3Config(d, above_emax="black-disk", inputs_json=_inputs(tmp_path),
                      allow_unpinned_inputs=True)  # fmt: skip
    out = tmp_path / "out"
    res = EB.build_elastic_proton(cache.resolve_cache_dir(None), m3=cfg, out_dir=out)
    ex = res.info["exploratory"]
    assert ex["model"] == "m3-endf-mirror" and ex["inputs_unpinned"] is True
    t = res.info["m3"]["targets"]
    assert t["O-16"]["shape"] == "endf-mf4-mt2" and t["O-16"]["n_nodes_above_emax"] > 0
    assert t["Ca-40"]["shape"].startswith("black-disk")
    with pytest.raises(EB.BuildError, match="above the ENDF evaluation range"):
        EB.build_elastic_proton(
            cache.resolve_cache_dir(None),
            m3=M3.M3Config(d, inputs_json=_inputs(tmp_path), allow_unpinned_inputs=True),
            out_dir=out,
        )
