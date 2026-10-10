"""V3-005C P6 (1)-(4): MF6 LAW=5 parser, ENDF-102 eqs 6.9/6.10/6.14 reconstruction, LA150 values,
ratio interpolation. Fixture cases are hand-written ENDF records; data cases use the LA150 H-1 and
O-16 files of the cache and fail under ``IONMC_REQUIRE_DATA=1`` when it is missing."""

from __future__ import annotations

import math
import os
from typing import Any

import numpy as np
import pytest

from ionmc.data import cache, endf6
from ionmc.data.endf6 import EndfSection, UnsupportedEndfError
from ionmc.nuclear import law5


def _f(x: float) -> str:
    return f"{x:11.4E}"


def _i(n: int) -> str:
    return f"{n:11d}"


def _rec(c1: float, c2: float, l1: int, l2: int, n1: int, n2: int) -> str:
    return _f(c1) + _f(c2) + _i(l1) + _i(l2) + _i(n1) + _i(n2)


def _floats(vals: list[float]) -> list[str]:
    out = []
    for s in range(0, len(vals), 6):
        out.append("".join(_f(v) for v in vals[s : s + 6]))
    return out


def _section(
    lidp: int,
    records: list[tuple[float, int, int, list[float]]],
    *,
    law: int = 5,
    tab2_int: int = 2,
) -> EndfSection:
    ne = len(records)
    lines = [
        _rec(1001.0, 0.99917, 0, 2, 1, 0),
        _rec(1001.0, 0.99662, 0, law, 1, 2),
        _i(2) + _i(2) + " " * 44,
        _f(1e3) + _f(1.0) + _f(1.5e8) + _f(1.0) + " " * 22,
        _rec(0.5, 0.0, lidp, 0, 1, ne),
        _i(ne) + _i(tab2_int) + " " * 44,
    ]
    for e, ltp, nl, a in records:
        lines.append(_rec(0.0, e, ltp, 0, len(a), nl))
        lines += _floats(a)
    return EndfSection(6, 2, tuple(lines))


A_LIDP1 = [0.1, 0.02, 0.5, 0.1, -0.3, 0.2]  # NL=1: b0, b1, Re/Im a0, Re/Im a1 (NW = 6)
A_LIDP0 = [0.1, 0.02, 0.01, 0.5, 0.1, -0.3, 0.2]  # NL=1: b0..b2, a0, a1 (NW = 4NL+3 = 7)


def test_parse_ltp1_lidp1_and_lidp0_layouts() -> None:
    s1 = law5.parse_law5(_section(1, [(1.0e6, 1, 1, A_LIDP1), (2.0e6, 1, 1, A_LIDP1)]))
    assert s1.lidp == 1 and s1.spi == 0.5 and len(s1.records) == 2 and s1.zap == 1001
    assert s1.records[0].data.size == 6 and s1.records[1].energy_ev == 2.0e6
    s0 = law5.parse_law5(_section(0, [(1.0e6, 1, 1, A_LIDP0 + [0.0] * 0)]))
    assert s0.lidp == 0 and s0.records[0].data.size == 7
    ltp12 = law5.parse_law5(_section(0, [(1.0e6, 12, 2, [-1.0, 0.1, 0.9, 2.0])]))
    assert ltp12.records[0].ltp == 12 and ltp12.records[0].data.size == 4  # NW = 2 NL


@pytest.mark.parametrize(
    "bad",
    [
        _section(1, [(1.0e6, 1, 1, A_LIDP1[:5])]),  # NW != 3NL+3
        _section(0, [(1.0e6, 1, 1, A_LIDP1)]),  # NW != 4NL+3 for LIDP 0
        _section(1, [(1.0e6, 2, 1, [0.1, 0.2])]),  # LTP=2 unsupported
        _section(1, [(1.0e6, 14, 1, [0.1, 0.2])]),  # LTP=14 unsupported
        _section(1, [(1.0e6, 1, 1, A_LIDP1)], law=1),  # not LAW=5
        _section(1, [(1.0e6, 1, 1, A_LIDP1)], tab2_int=5),  # energy interpolation not lin-lin
        _section(1, [(2.0e6, 1, 1, A_LIDP1), (1.0e6, 1, 1, A_LIDP1)]),  # energies not ascending
    ],
)
def test_parser_fails_closed(bad: EndfSection) -> None:
    with pytest.raises(UnsupportedEndfError):
        law5.parse_law5(bad)


def test_parser_floats_are_bitwise() -> None:
    s = law5.parse_law5(_section(1, [(1.0e6, 1, 1, A_LIDP1)]))
    expected = np.array([float(f"{v:11.4E}") for v in A_LIDP1])
    assert np.array_equal(s.records[0].data, expected)


def test_distinguishable_coulomb_is_rutherford() -> None:
    """Eq 6.9 equals the Rutherford cross section (Z1 Z2 alpha hbarc / (4 T_cm))^2 / sin^4(theta/2)
    with T_cm = E A/(1+A) (non-relativistic, A = AWR/AWI), in b/sr."""
    e_ev, awr, awi, z2 = 20.0e6, 15.85751, 0.99862, 8
    k, eta = law5.coulomb_params(e_ev, awr, awi, 1, z2)
    a = awr / awi
    t_cm_mev = 20.0 * a / (1.0 + a)
    for th in (5.0, 20.0, 60.0, 120.0):
        mu = np.array([math.cos(math.radians(th))])
        ruth_fm2 = (z2 * law5.ALPHA * 197.3269804 / (4.0 * t_cm_mev)) ** 2 / math.sin(
            math.radians(th) / 2
        ) ** 4
        assert law5.sigma_c(mu, k, eta, 0, 0.0)[0] == pytest.approx(ruth_fm2 * 0.01, rel=2e-6)


def test_identical_coulomb_by_hand_at_90_degrees() -> None:
    """mu = 0, s = 1/2: (2 eta^2/k^2) [1 - cos(0)/2] = eta^2/k^2 (cos coefficient -1/2)."""
    k, eta = law5.coulomb_params(10.0e6, 0.99917, 0.99862, 1, 1)
    assert law5.sigma_c(np.array([0.0]), k, eta, 1, 0.5)[0] == pytest.approx(
        eta**2 / k**2, rel=1e-14
    )
    # symmetric under mu -> -mu (identical particles)
    mu = np.array([0.3, 0.7])
    assert np.allclose(
        law5.sigma_c(mu, k, eta, 1, 0.5), law5.sigma_c(-mu, k, eta, 1, 0.5), rtol=1e-13
    )


def test_ltp1_pure_nuclear_limit_and_symmetry() -> None:
    """With a_l = 0 the cross section is Coulomb + sum (4l+1)/2 b_l P_2l (even in mu)."""
    rec = law5.Law5Energy(1.0e7, 1, 1, np.array([0.1, 0.02, 0.0, 0.0, 0.0, 0.0]))
    k, eta = law5.coulomb_params(1.0e7, 0.99917, 0.99862, 1, 1)
    mu = np.array([0.0, 0.4, 0.9])
    se, nuc, intf = law5.sigma_e_ltp1(mu, rec, 1, k, eta, 0.5)
    p2 = 0.5 * (3 * mu**2 - 1)
    assert np.allclose(nuc, 0.5 * 0.1 + 2.5 * 0.02 * p2, rtol=1e-13) and np.all(intf == 0.0)
    assert np.allclose(se, law5.sigma_c(mu, k, eta, 1, 0.5) + nuc, rtol=1e-13)


def _material(member: str) -> endf6.EndfMaterial:
    cdir = cache.resolve_cache_dir(None)
    try:
        zpath = cache.verify("endf-b8.0-protons", cdir)
    except FileNotFoundError:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail("ENDF/B-VIII.0 proton sublibrary is not cached")
        pytest.skip("LA150 not cached")
    return endf6.parse_endf(endf6.read_member(zpath, f"ENDF-B-VIII.0_protons/{member}.endf"))


def test_h1_ni_cross_sections_match_la150_reconstruction() -> None:
    """P6 (3): H-1 NI over theta_CM >= 16.26 deg 49.8 / 27.3 / 26.3 mb at 50 / 100 / 150 MeV
    (1e-3 relative), and above 10 deg 49.4 / 27.4 / 26.9 mb (research reconstruction)."""
    m = _material("p-001_H_001")
    sec, awi = law5.parse_law5(m.sections[(6, 2)]), law5.projectile_awi(m)
    assert sec.lidp == 1 and {r.ltp for r in sec.records} == {1}
    cos10 = math.cos(math.radians(10.0))
    for t, ni1626, ni10 in ((50.0, 49.8, 49.4), (100.0, 27.3, 27.4), (150.0, 26.3, 26.9)):
        got = law5.ni_cross_section_b(sec, t * 1e6, law5.MU_CUT_PP, awi, 1, 1) * 1e3
        assert got == pytest.approx(ni1626, rel=1e-3 + 0.05 / ni1626)  # quoted to 0.1 mb
        got10 = law5.ni_cross_section_b(sec, t * 1e6, cos10, awi, 1, 1) * 1e3
        assert got10 == pytest.approx(ni10, abs=0.06)


def test_o16_100mev_ltp12_reconstruction() -> None:
    """P6 (3): O-16 100 MeV NI(theta >= 5 deg) 123.8 mb and sigma_e(theta >= 5 deg) 185.6 mb."""
    m = _material("p-008_O_016")
    sec, awi = law5.parse_law5(m.sections[(6, 2)]), law5.projectile_awi(m)
    idx = [r.energy_ev for r in sec.records].index(100.0e6)
    sig = float(m.cross_section(2).interpolate(100.0e6)[0])
    out = law5.ltp12_cross_sections_b(sec, idx, sig, awi, 1, 8)
    assert out["ni_b"] * 1e3 == pytest.approx(123.8, rel=1e-3)
    assert out["total_b"] * 1e3 == pytest.approx(185.6, rel=1e-3)
    assert out["int_p"] == pytest.approx(1.0, abs=1e-5)


def test_h1_ratio_interpolation_nodes_midpoints_and_positivity() -> None:
    """P6 (4): at a node the ratio is the node's own; between nodes it is the linear blend of the
    two node ratios; the transported density is non-negative at every node and node midpoint from
    20 MeV up over |mu| <= 0.96 (below that the builder's E_min_pp applies)."""
    m = _material("p-001_H_001")
    sec, awi = law5.parse_law5(m.sections[(6, 2)]), law5.projectile_awi(m)
    mu = np.linspace(0.0, law5.MU_CUT_PP, 2001)
    en = sec.energies_ev
    j = int(np.searchsorted(en, 60.0e6))
    k, eta = law5.coulomb_params(en[j], sec.awr, awi, 1, 1)
    se, _, _ = law5.sigma_e_ltp1(mu, sec.records[j], 1, k, eta, 0.5)
    r_node = se / law5.sigma_c(mu, k, eta, 1, 0.5)
    assert np.allclose(law5.ratio_ltp1(sec, en[j], mu, awi, 1, 1), r_node, rtol=1e-12)
    e_mid = 0.5 * (en[j] + en[j + 1])
    r_mid = law5.ratio_ltp1(sec, e_mid, mu, awi, 1, 1)
    r_hi = law5.ratio_ltp1(sec, en[j + 1], mu, awi, 1, 1)
    assert np.allclose(r_mid, 0.5 * (r_node + r_hi), rtol=1e-12)
    sel: list[Any] = [e for e in en if 20.0e6 <= e <= 150.0e6]
    worst = 0.0
    for e0, e1 in zip(sel[:-1], sel[1:], strict=True):
        for e in (e0, 0.5 * (e0 + e1)):
            worst = min(worst, float(law5.ni_density_ltp1(sec, e, mu, awi, 1, 1).min()))
    assert worst >= 0.0
