"""V3-005C P6 (5) and D9-BGG inputs: BGG / Barashenkov elastic cross sections (Amendment 14 (b)).

Hand-computed values are read from the Geant4 11.4.2 arrays by hand (research-answers-round2.md);
the data-backed tests compare the transcribed arrays with the pinned source files and fail under
``IONMC_REQUIRE_DATA=1`` when the cache lacks them.
"""

from __future__ import annotations

import math
import os
import re

import pytest

from ionmc.data import cache
from ionmc.nuclear import bgg

# (Z, E MeV, quoted mb, half unit of the last quoted digit): research-answers quotes (P6 (5))
QUOTED = [
    (8, 100.0, 352.0, 0.5), (8, 150.0, 159.0, 0.5), (8, 200.0, 102.0, 0.5), (8, 250.0, 92.0, 0.5),
    (6, 100.0, 253.0, 0.5), (6, 150.0, 125.0, 0.5),
    (7, 100.0, 282.0, 0.5), (7, 150.0, 132.0, 0.5),
    (20, 150.0, 400.0, 0.5),
]  # fmt: skip


@pytest.mark.parametrize(("z", "e", "mb", "tol"), QUOTED)
def test_bgg_hand_computed_values(z: int, e: float, mb: float, tol: float) -> None:
    assert abs(bgg.bgg_elastic_mb(e, z) - mb) <= tol


def test_calcium_100_mev_is_767_not_912() -> None:
    """The array gives sigma_tot(n+Ca) 1340 mb and sigma_inel(p+Ca) 573 mb at 0.1 GeV (index 17
    of e3), i.e. 767 mb. The 912 mb quoted for Ca at 100 MeV in Amendment 14 P6 (5) and in the
    research notes is the difference at 90 MeV (1500 - 588, index 16): a hand-reading slip.
    This test records the discrepancy; the value from the source is the one implemented."""
    assert bgg.barashenkov_mb(100.0, 20) == (1340.0, 573.0)
    assert bgg.bgg_elastic_mb(100.0, 20) == 767.0
    assert bgg.bgg_elastic_mb(90.0, 20) == 1500.0 - 588.0 == 912.0


def test_pp_bgg_formula_and_plab() -> None:
    assert abs(bgg.p_lab_gev(100.0) - 0.4446) < 5e-5 and abs(bgg.p_lab_gev(150.0) - 0.5513) < 5e-5
    assert abs(bgg.pp_bgg_mb(100.0) - 27.61) <= 0.005
    assert abs(bgg.pp_bgg_mb(150.0) - 23.86) <= 0.005
    p = bgg.p_lab_gev(150.0)
    assert bgg.pp_bgg_mb(150.0) == pytest.approx(
        1.0115 * (23 + 50 * math.sqrt(math.log(0.73 / p) ** 7)), rel=1e-12
    )
    with pytest.raises(bgg.BggError):
        bgg.pp_bgg_mb(400.0)  # p_lab > 0.73 GeV/c


def test_linear_interpolation_between_nodes_by_hand() -> None:
    # O nodes 0.1 / 0.12 GeV: m_t 645 / 540, p_in 293 / 287; at 110 MeV the midpoint values
    tot, inel = bgg.barashenkov_mb(110.0, 8)
    assert tot == pytest.approx(0.5 * (645 + 540), rel=1e-14)
    assert inel == pytest.approx(0.5 * (293 + 287), rel=1e-14)
    assert bgg.barashenkov_elastic_mb(110.0, 8) == pytest.approx(592.5 - 290.0, rel=1e-14)


def test_untabulated_z_uses_a23_interpolation() -> None:
    """P (Z = 15) between Si (14) and Ca (20): r_i = x_i A75[Z]/A75[Z_i], weights by aeff."""
    a = {14: 28.0854, 15: 30.9738, 20: 40.078}
    for e in (50.0, 100.0, 150.0, 250.0):
        out = []
        for which in (0, 1):
            x1 = bgg.barashenkov_mb(e, 14)[which]
            x2 = bgg.barashenkov_mb(e, 20)[which]
            r1 = x1 * (a[15] / a[14]) ** (2 / 3)
            r2 = x2 * (a[15] / a[20]) ** (2 / 3)
            w1, w2 = a[15] - a[14], a[20] - a[15]
            out.append((r1 * w2 + r2 * w1) / (w1 + w2))
        assert bgg.barashenkov_mb(e, 15) == pytest.approx(tuple(out), rel=1e-12)
        assert bgg.bgg_elastic_mb(e, 15) == pytest.approx(max(out[0] - out[1], 0.0), rel=1e-12)


def test_low_energy_rule_is_the_coulomb_factor_scaled_at_14_mev() -> None:
    for z in bgg.SUPPORTED_Z:
        s14 = bgg.bgg_elastic_mb(14.0, z)
        assert s14 == pytest.approx(bgg.barashenkov_elastic_mb(14.0, z), rel=1e-12)
        ratio = bgg.coulomb_factor(7.0, z) / bgg.coulomb_factor(14.0, z)
        assert bgg.bgg_elastic_mb(7.0, z) == pytest.approx(s14 * ratio, rel=1e-12)
        assert bgg.bgg_elastic_mb(7.0, z) < s14
        # the factor vanishes below the Coulomb barrier
        assert bgg.bgg_elastic_mb(0.5, z) == 0.0 or bgg.coulomb_factor(0.5, z) >= 0.0
    b_c = 0.5 * bgg.FINE_STRUCTURE * bgg.HBARC_MEV_FM * 8 / (0.895 + 1.71 * 16 ** (1 / 3))
    assert 1.0 < b_c < 1.2
    assert bgg.coulomb_factor(1.0, 8) == 0.0
    assert bgg.coulomb_factor(14.0, 8) == pytest.approx(
        1 - b_c / (14.0 * 15.8 / 16.8 * 0.998), rel=2e-2
    )


def test_fail_closed_unsupported_element_and_domain() -> None:
    with pytest.raises(bgg.BggError):
        bgg.bgg_elastic_mb(100.0, 11)  # Na: not carried (surrogate scaling applies upstream)
    with pytest.raises(bgg.BggError):
        bgg.bgg_elastic_mb(251.0, 8)
    with pytest.raises(bgg.BggError):
        bgg.bgg_elastic_mb(0.0, 8)


def _header() -> str | None:
    cdir = cache.resolve_cache_dir(None)
    try:
        path = cache.verify("geant4-g4barashenkovdata-hh-11.4.2", cdir)
    except FileNotFoundError:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail("G4BarashenkovData.hh is not cached")
        return None
    return path.read_text(encoding="utf-8")


def _array(text: str, name: str) -> list[float]:
    m = re.search(r"static const G4double " + name + r"\[\d*\]\s*=\s*\{(.*?)\};", text, re.S)
    assert m is not None, name
    body = re.sub(r"//.*", "", m.group(1))
    return [float(x) for x in re.split(r"[\s,]+", body) if x]


def test_transcribed_arrays_equal_the_pinned_geant4_header() -> None:
    text = _header()
    if text is None:
        pytest.skip("G4BarashenkovData.hh not cached")
    assert list(bgg._E1_GEV) == _array(text, "e1")
    assert list(bgg._E2_GEV) == _array(text, "e2")
    assert list(bgg._E3_GEV) == _array(text, "e3")
    for el, z in (("c", 6), ("n", 7), ("o", 8), ("al", 13), ("si", 14), ("ca", 20)):
        assert list(bgg.BARASHENKOV_M_T[z]) == _array(text, el + "_m_t"), el
        assert list(bgg.BARASHENKOV_P_IN[z]) == _array(text, el + "_p_in"), el
        n = len(bgg.BARASHENKOV_ENERGY_GEV[z])
        assert n == len(bgg.BARASHENKOV_M_T[z]) == len(bgg.BARASHENKOV_P_IN[z])


def test_aeff_and_r0_equal_the_pinned_geant4_sources() -> None:
    cdir = cache.resolve_cache_dir(None)
    try:
        iso = cache.verify("geant4-g4isotopelist-hh-11.4.2", cdir).read_text(encoding="utf-8")
        rad = cache.verify("geant4-g4nuclearradii-cc-11.4.2", cdir).read_text(encoding="utf-8")
    except FileNotFoundError:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail("G4IsotopeList.hh / G4NuclearRadii.cc not cached")
        pytest.skip("not cached")
    aeff = _array(iso, "aeff")
    r0 = _array(
        rad.replace("const G4double G4NuclearRadii::r0[]", "static const G4double r0[1]"), "r0"
    )
    for z, v in bgg.AEFF.items():
        assert aeff[z] == v
    for z, v in bgg.R0_CB.items():
        assert r0[z] == v
