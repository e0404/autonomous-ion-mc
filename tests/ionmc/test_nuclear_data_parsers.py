"""Tests for the ENDF-6, AME2020, EXFOR and geant-val parsers (V3-005A, check P1).

Fixtures are hand-written synthetic text with made-up numbers in the record layouts of the
formats; no data of the real files is stored in Git. The data-backed tests read the verified
datasets from the cache (``IONMC_CACHE_DIR``) and are skipped when a dataset is not cached;
``IONMC_REQUIRE_DATA=1`` turns a missing dataset into a failure.
"""

from __future__ import annotations

import json
import os
import zipfile
from pathlib import Path

import numpy as np
import pytest

from ionmc.data import ame, cache, endf6, exfor, geant_val

# --------------------------------------------------------------------------------------------
# Synthetic ENDF-6 text
# --------------------------------------------------------------------------------------------

MAT = 9999


def ef(x: float) -> str:
    """Format a float the ENDF way (no ``E``): 1.500000e+06 -> '1.500000+6'."""
    mant, exp = f"{x:.6e}".split("e")
    return f"{mant}{int(exp):+d}"


def rec(mf: int, mt: int, ns: int, *fields: str) -> str:
    body = "".join(f.rjust(11) for f in fields).ljust(66)
    return f"{body}{MAT:4d}{mf:2d}{mt:3d}{ns:5d}"


class Sec:
    """Accumulates the records of one section and numbers the lines."""

    def __init__(self, mf: int, mt: int) -> None:
        self.mf, self.mt, self.lines = mf, mt, []

    def cont(self, c1: str, c2: str, l1: int, l2: int, n1: int, n2: int) -> None:
        self.lines.append(rec(self.mf, self.mt, len(self.lines) + 1, c1, c2, str(l1), str(l2), str(n1), str(n2)))

    def values(self, vals: list[str]) -> None:
        for i in range(0, len(vals), 6):
            self.lines.append(rec(self.mf, self.mt, len(self.lines) + 1, *vals[i : i + 6]))

    def tab1(self, c1: str, c2: str, l1: int, l2: int, nbt_int: list[tuple[int, int]],
             xy: list[tuple[float, float]]) -> None:
        self.cont(c1, c2, l1, l2, len(nbt_int), len(xy))
        self.values([str(v) for pair in nbt_int for v in pair])
        self.values([ef(v) for pair in xy for v in pair])

    def tab2(self, c1: str, c2: str, l1: int, l2: int, nbt_int: list[tuple[int, int]], nz: int) -> None:
        self.cont(c1, c2, l1, l2, len(nbt_int), nz)
        self.values([str(v) for pair in nbt_int for v in pair])

    def lst(self, c1: str, c2: str, l1: int, l2: int, n2: int, body: list[float]) -> None:
        self.cont(c1, c2, l1, l2, len(body), n2)
        self.values([ef(v) for v in body])

    def text(self) -> str:
        send = rec(self.mf, 0, 99999, *[""] * 6)
        return "\n".join([*self.lines, send.replace(f"{self.mf:2d}{0:3d}", f"{self.mf:2d}  0", 1)])


def endf_text(sections: list[Sec]) -> str:
    tpid = " synthetic fixture, not evaluated data".ljust(66) + f"{1:4d}{0:2d}{0:3d}{0:5d}"
    fend = " " * 66 + f"{MAT:4d}{0:2d}{0:3d}{0:5d}"
    mend = " " * 66 + f"{0:4d}{0:2d}{0:3d}{0:5d}"
    tend = " " * 66 + f"{-1:4d}{0:2d}{0:3d}{0:5d}"
    return "\n".join([tpid, *(s.text() for s in sections), fend, mend, tend]) + "\n"


def mf1_451() -> Sec:
    s = Sec(1, 451)
    s.cont("6.012000+3", "1.190780+1", 0, 0, 0, 0)  # HEAD ZA AWR ...
    s.cont("0.000000+0", "0.000000+0", 0, 0, 0, 6)
    s.cont("9.986200-1", "1.500000+8", 0, 0, 10010, 8)  # AWI EMAX
    s.cont("0.000000+0", "0.000000+0", 0, 0, 1, 5)
    s.lines.append(rec(1, 451, len(s.lines) + 1, "synthetic text line"))
    return s


def mf3_mt5() -> Sec:
    s = Sec(3, 5)
    s.cont("6.012000+3", "1.190780+1", 0, 0, 0, 0)
    pts = [(1.0e6, 1.0), (2.0e6, 3.0), (4.0e6, 2.0), (8.0e6, 4.0), (1.6e7, 16.0)]
    s.tab1("0.000000+0", " 0.000000+0", 0, 0, [(3, 2), (5, 5)], pts)
    return s


def mf6_mt5() -> Sec:
    s = Sec(6, 5)
    s.cont("6.012000+3", "1.190780+1", 0, 3, 3, 0)  # NK = 3, LCT = 3
    # product 1: neutron, Kalbach-Mann, LEP = 1
    s.tab1("1.000000+0", "1.000000+0", 0, 1, [(2, 2)], [(1.0e6, 0.0), (1.0e8, 1.5)])
    s.tab2("0.000000+0", "0.000000+0", 2, 1, [(2, 2)], 2)
    s.lst("0.000000+0", ef(5.0e6), 0, 1, 2, [1.0e5, 2.0e-6, 0.25, 2.0e5, 1.0e-5, 0.5])
    s.lst("0.000000+0", ef(1.0e7), 0, 1, 3,
          [1.0e5, 1.0e-6, 0.125, 2.0e5, 2.0e-6, 0.25, 3.0e5, 3.0e-6, 0.75])
    # product 2: alpha, Legendre NA = 2, LEP = 2
    s.tab1("2.004000+3", "3.972600+0", 1, 1, [(2, 4)], [(1.0e6, 0.0), (1.0e8, 0.5)])
    s.tab2("0.000000+0", "0.000000+0", 1, 2, [(1, 2)], 1)
    s.lst("0.000000+0", ef(2.0e7), 2, 2, 2, [1.0e5, 1.0, 0.5, 0.25, 2.0e5, 2.0, 1.0, 0.75])
    # product 3: heavy residue (LANG=1, NA=0)
    s.tab1("6.011000+3", "1.000000+1", 0, 1, [(2, 2)], [(1.0e6, 0.0), (1.0e8, 1.0)])
    s.tab2("0.000000+0", "0.000000+0", 1, 1, [(1, 2)], 1)
    s.lst("0.000000+0", ef(2.0e7), 0, 0, 1, [1.0e3, 1.0])
    return s


def mf6_law5() -> Sec:
    s = Sec(6, 2)
    s.cont("1.001000+3", "9.986200-1", 0, 3, 2, 0)
    s.tab1("1.001000+3", "9.986200-1", 0, 5, [(2, 2)], [(1.0e6, 1.0), (2.0e8, 1.0)])
    s.tab2("5.000000-1", "0.000000+0", 0, 0, [(2, 2)], 2)
    s.lst("0.000000+0", ef(1.0e6), 1, 0, 3, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
    s.lst("0.000000+0", ef(2.0e6), 1, 0, 3, [0.1, 0.2, 0.3])
    s.tab1("1.001000+3", "9.986200-1", 0, 1, [(2, 2)], [(1.0e6, 1.0), (2.0e8, 1.0)])
    s.tab2("0.000000+0", "0.000000+0", 1, 1, [(1, 2)], 1)
    s.lst("0.000000+0", ef(2.0e7), 0, 0, 1, [1.0e3, 1.0])
    return s


def test_float_and_int_fields() -> None:
    assert endf6.parse_endf_float(" 8.016000+3") == 8016.0
    assert endf6.parse_endf_float(" 9.986200-1") == 0.9986200
    assert endf6.parse_endf_float(" 3.319869-4") == 3.319869e-4
    assert endf6.parse_endf_float("1.2345E+03") == 1234.5
    assert endf6.parse_endf_float(" 0.000000+0") == 0.0
    assert endf6.parse_endf_float("-1.500000+2") == -150.0
    assert endf6.parse_endf_float("          ") == 0.0
    assert endf6.parse_endf_float("       0.25") == 0.25
    assert endf6.parse_endf_int("         21") == 21
    assert endf6.parse_endf_int("         -1") == -1
    for bad in ("1.5.2+3", "abc"):
        with pytest.raises(endf6.EndfError):
            endf6.parse_endf_float(bad)
    with pytest.raises(endf6.EndfError):
        endf6.parse_endf_int("1.5")


def test_material_header_and_mf1() -> None:
    mat = endf6.parse_endf(endf_text([mf1_451(), mf3_mt5()]))
    assert mat.za == 6012.0
    assert mat.awr == 11.9078
    assert mat.emax == 1.5e8
    assert mat.mat == MAT
    assert set(mat.sections) == {(1, 451), (3, 5)}


def test_mf3_tab1_regions_and_values() -> None:
    mat = endf6.parse_endf(endf_text([mf1_451(), mf3_mt5()]))
    tab = mat.cross_section(5)
    assert tab.nbt.tolist() == [3, 5]
    assert tab.interp.tolist() == [2, 5]
    assert tab.x.tolist() == [1.0e6, 2.0e6, 4.0e6, 8.0e6, 1.6e7]
    assert tab.y.tolist() == [1.0, 3.0, 2.0, 4.0, 16.0]
    # region 1 (lin-lin): exact binary values
    assert tab.interpolate(1.5e6)[0] == 2.0
    assert tab.interpolate(3.0e6)[0] == 2.5
    # nodes, including the region boundary shared by both laws
    assert tab.interpolate(np.array([1.0e6, 2.0e6, 4.0e6, 8.0e6, 1.6e7])).tolist() == [
        1.0, 3.0, 2.0, 4.0, 16.0,
    ]
    # region 2 (log-log): y = 2 (x/4e6)^1 on (4e6, 8e6) and 4 (x/8e6)^2 on (8e6, 1.6e7)
    assert tab.interpolate(6.0e6)[0] == pytest.approx(3.0, rel=1e-14)
    assert tab.interpolate(1.2e7)[0] == pytest.approx(9.0, rel=1e-14)
    with pytest.raises(ValueError, match="outside"):
        tab.interpolate(9.0e5)
    with pytest.raises(ValueError, match="outside"):
        tab.interpolate(1.7e7)


def _one_region(law: int, x: list[float], y: list[float]) -> endf6.Tab1:
    return endf6.Tab1(
        nbt=np.array([len(x)]), interp=np.array([law]), x=np.array(x), y=np.array(y)
    )


def test_all_five_interpolation_laws() -> None:
    # 1 histogram: constant y_i over [x_i, x_{i+1})
    t1 = _one_region(1, [1.0, 2.0, 4.0], [10.0, 20.0, 30.0])
    assert t1.interpolate(np.array([1.0, 1.5, 1.999, 2.0, 3.0])).tolist() == [
        10.0, 10.0, 10.0, 20.0, 20.0,
    ]
    # 2 lin-lin
    t2 = _one_region(2, [1.0, 3.0], [10.0, 20.0])
    assert t2.interpolate(2.0)[0] == 15.0
    assert t2.interpolate(2.5)[0] == 17.5
    # 3 lin-log: y linear in ln x; x = 2 is the geometric mean of 1 and 4
    t3 = _one_region(3, [1.0, 4.0], [0.0, 10.0])
    assert t3.interpolate(2.0)[0] == pytest.approx(5.0, rel=1e-15)
    # 4 log-lin: ln y linear in x; y(1) is the geometric mean of 1 and 100
    t4 = _one_region(4, [0.0, 2.0], [1.0, 100.0])
    assert t4.interpolate(1.0)[0] == pytest.approx(10.0, rel=1e-15)
    # 5 log-log: y = x^2
    t5 = _one_region(5, [1.0, 100.0], [1.0, 1.0e4])
    assert t5.interpolate(10.0)[0] == pytest.approx(100.0, rel=1e-14)
    # log laws reject non-positive data on the used interval
    with pytest.raises(ValueError, match="positive"):
        _one_region(5, [1.0, 2.0], [0.0, 1.0]).interpolate(1.5)
    with pytest.raises(ValueError, match="positive"):
        _one_region(4, [1.0, 2.0], [1.0, 0.0]).interpolate(1.5)
    with pytest.raises(ValueError, match="positive"):
        _one_region(3, [0.0, 2.0], [1.0, 2.0]).interpolate(1.0)


def test_discontinuity_takes_upper_branch() -> None:
    tab = endf6.Tab1(
        nbt=np.array([2, 4]),
        interp=np.array([2, 2]),
        x=np.array([1.0, 2.0, 2.0, 3.0]),
        y=np.array([1.0, 2.0, 5.0, 6.0]),
    )
    assert tab.interpolate(np.array([1.5, 2.0, 2.5])).tolist() == [1.5, 5.0, 5.5]


def test_tab1_validation() -> None:
    with pytest.raises(endf6.EndfError):
        endf6.Tab1(np.array([2]), np.array([6]), np.array([1.0, 2.0]), np.array([1.0, 2.0]))
    with pytest.raises(endf6.EndfError):
        endf6.Tab1(np.array([3]), np.array([2]), np.array([1.0, 2.0]), np.array([1.0, 2.0]))
    with pytest.raises(endf6.EndfError):
        endf6.Tab1(np.array([2]), np.array([2]), np.array([2.0, 1.0]), np.array([1.0, 2.0]))


def test_mf6_law1_lang2_and_lang1() -> None:
    mat = endf6.parse_endf(endf_text([mf6_mt5()]))
    sec = mat.products(5)
    assert sec.lct == 3
    assert sec.za == 6012.0
    n, a, r = sec.products
    assert (n.zap, n.awp, n.lip, n.law, n.lang, n.lep, n.nr, n.ne) == (1, 1.0, 0, 1, 2, 1, 1, 2)
    assert n.yield_.y.tolist() == [0.0, 1.5]
    assert n.energy_nbt.tolist() == [2] and n.energy_interp.tolist() == [2]  # type: ignore[union-attr]
    d0, d1 = n.distributions
    assert (d0.energy, d0.nd, d0.na, d0.nw, d0.nep) == (5.0e6, 0, 1, 6, 2)
    assert d0.rows.tolist() == [[1.0e5, 2.0e-6, 0.25], [2.0e5, 1.0e-5, 0.5]]  # E', f0, r
    assert (d1.energy, d1.nw, d1.nep) == (1.0e7, 9, 3)
    assert d1.rows.shape == (3, 3)
    assert d1.rows[2].tolist() == [3.0e5, 3.0e-6, 0.75]
    assert (a.zap, a.awp, a.lip, a.lang, a.lep, a.ne) == (2004, 3.9726, 1, 1, 2, 1)
    assert a.yield_.interp.tolist() == [4]
    (da,) = a.distributions
    assert (da.nd, da.na, da.nw, da.nep) == (2, 2, 8, 2)
    assert da.rows.tolist() == [[1.0e5, 1.0, 0.5, 0.25], [2.0e5, 2.0, 1.0, 0.75]]
    assert (r.zap, r.distributions[0].na, r.distributions[0].rows.tolist()) == (
        6011, 0, [[1.0e3, 1.0]]
    )


def test_mf6_law5_is_opaque_and_following_product_is_parsed() -> None:
    sec = endf6.parse_endf(endf_text([mf6_law5()])).products(2)
    first, second = sec.products
    assert first.law == 5 and first.lang is None and first.distributions == ()
    # TAB2 (1 line + 1 interpolation line) + LIST(1 + 2 lines) + LIST(1 + 1 line)
    assert len(first.opaque_records) == 2 + 3 + 2
    assert second.law == 1 and second.distributions[0].rows.tolist() == [[1.0e3, 1.0]]


def _replace_in(sec: Sec, old: str, new: str) -> Sec:
    sec.lines = [ln.replace(old, new) for ln in sec.lines]
    return sec


def test_unsupported_and_malformed_mf6_fail_closed() -> None:
    # a product with LAW = 2
    s = Sec(6, 5)
    s.cont("6.012000+3", "1.190780+1", 0, 3, 1, 0)
    s.tab1("1.000000+0", "1.000000+0", 0, 2, [(2, 2)], [(1.0e6, 0.0), (1.0e8, 1.5)])
    with pytest.raises(endf6.UnsupportedEndfError, match="LAW=2"):
        endf6.parse_endf(endf_text([s])).products(5)
    # LAW = 1 with LANG = 14 (tabulated angular distribution)
    s = Sec(6, 5)
    s.cont("6.012000+3", "1.190780+1", 0, 3, 1, 0)
    s.tab1("1.000000+0", "1.000000+0", 0, 1, [(2, 2)], [(1.0e6, 0.0), (1.0e8, 1.5)])
    s.tab2("0.000000+0", "0.000000+0", 14, 1, [(2, 2)], 1)
    with pytest.raises(endf6.UnsupportedEndfError, match="LANG=14"):
        endf6.parse_endf(endf_text([s])).products(5)
    # Kalbach-Mann with NA = 2
    s = Sec(6, 5)
    s.cont("6.012000+3", "1.190780+1", 0, 3, 1, 0)
    s.tab1("1.000000+0", "1.000000+0", 0, 1, [(2, 2)], [(1.0e6, 0.0), (1.0e8, 1.5)])
    s.tab2("0.000000+0", "0.000000+0", 2, 1, [(1, 2)], 1)
    s.lst("0.000000+0", ef(1.0e6), 0, 2, 1, [1.0, 2.0, 3.0, 4.0])
    with pytest.raises(endf6.UnsupportedEndfError, match="NA=1"):
        endf6.parse_endf(endf_text([s])).products(5)
    # inconsistent NW
    s = Sec(6, 5)
    s.cont("6.012000+3", "1.190780+1", 0, 3, 1, 0)
    s.tab1("1.000000+0", "1.000000+0", 0, 1, [(2, 2)], [(1.0e6, 0.0), (1.0e8, 1.5)])
    s.tab2("0.000000+0", "0.000000+0", 1, 1, [(1, 2)], 1)
    s.cont("0.000000+0", ef(1.0e6), 0, 1, 5, 2)  # NW = 5 but NEP*(NA+2) = 6
    s.values([ef(v) for v in (1.0, 2.0, 3.0, 4.0, 5.0)])
    with pytest.raises(endf6.EndfError, match="NW"):
        endf6.parse_endf(endf_text([s])).products(5)
    # truncated section
    s = Sec(6, 5)
    s.cont("6.012000+3", "1.190780+1", 0, 3, 2, 0)
    s.tab1("1.000000+0", "1.000000+0", 0, 5, [(2, 2)], [(1.0e6, 0.0), (1.0e8, 1.5)])
    with pytest.raises(endf6.EndfError):
        endf6.parse_endf(endf_text([s])).products(5)


def test_parse_endf_rejects_bad_text() -> None:
    with pytest.raises(endf6.EndfError, match="no data"):
        endf6.parse_endf("")
    bad = endf_text([mf3_mt5()]).replace("9999", "xxxx", 1)
    with pytest.raises(endf6.EndfError):
        endf6.parse_endf(bad)
    two = endf_text([mf3_mt5()]) + endf_text([mf3_mt5()]).replace("9999", "8888")
    with pytest.raises(endf6.EndfError, match="one MAT"):
        endf6.parse_endf(two)


def test_zip_members_validated(tmp_path: Path) -> None:
    zp = tmp_path / "sub.zip"
    good = "ENDF-B-VIII.0_protons/p-006_C_012.endf"
    with zipfile.ZipFile(zp, "w") as zf:
        zf.writestr(good, endf_text([mf3_mt5()]))
        zf.writestr("ENDF-B-VIII.0_protons/readme.txt", "x")
        zf.writestr("../evil.endf", "x")
    assert endf6.list_zip_members(zp) == [good]
    assert endf6.parse_endf(endf6.read_member(zp, good)).za == 6012.0
    for name in ("../evil.endf", "ENDF-B-VIII.0_protons/readme.txt", "/etc/passwd",
                 "ENDF-B-VIII.0_protons/../p-006_C_012.endf"):
        with pytest.raises(ValueError, match="valid ENDF"):
            endf6.read_member(zp, name)
    with pytest.raises(KeyError):
        endf6.read_member(zp, "ENDF-B-VIII.0_protons/p-008_O_016.endf")


# --------------------------------------------------------------------------------------------
# Synthetic AME2020
# --------------------------------------------------------------------------------------------

AME_HEADER = """1    synthetic fixture, not AME2020 data
0    header line one
0    header line two      col 1 : Fortran character control
"""


def ame_line(cc: str, n: int, z: int, a: int, el: str, mass: str, unc: str,
             hi: int, lo: str, lounc: str, o: str = "") -> str:
    return (f"{cc}{z - n:3d}{n:5d}{z:5d}{a:5d} {el:>3}{o:>4} {mass:>14}{unc:>12}"
            f"{'0.0':>13} {'0.0':>10} B-{'*':>13}{'':>11} {hi:3d} {lo:>13}{lounc:>12}")


AME_TEXT = (
    AME_HEADER
    + ame_line("0", 1, 1, 2, "H", "100.5", "0.25", 2, "014000.5", "0.5") + "\n"
    + ame_line("0", 8, 8, 16, "O", "-2000.25", "0.5", 15, "990000.25", "0.75") + "\n"
    + ame_line("-", 2, 3, 5, "Li", "28667#", "2000#", 5, "030775#", "2147#", " -pp") + "\n"
)


def test_ame_parse_synthetic_lines() -> None:
    table = ame.load_ame2020(AME_TEXT)
    assert set(table) == {(1, 2), (8, 16), (3, 5)}
    o = table[(8, 16)]
    assert (o.element, o.mass_excess_kev, o.err_kev, o.estimated) == ("O", -2000.25, 0.5, False)
    assert o.atomic_mass_u == pytest.approx(15.99000025, abs=1e-12)
    assert o.err_u == pytest.approx(0.75e-6, abs=1e-18)
    li = table[(3, 5)]
    assert li.estimated is True
    assert (li.mass_excess_kev, li.err_kev) == (28667.0, 2000.0)  # '#' stands for the decimal point
    assert li.atomic_mass_u == pytest.approx(5.030775, abs=1e-12)
    assert li.err_u == pytest.approx(2147.0e-6, rel=1e-12)
    assert not table[(1, 2)].estimated
    assert table[(1, 2)].atomic_mass_u == pytest.approx(2.0140005, abs=1e-12)


def test_ame_nuclear_mass_and_errors() -> None:
    table = ame.load_ame2020(AME_TEXT)
    expect = table[(8, 16)].atomic_mass_u * ame.U_MEV - 8 * ame.ELECTRON_MASS_MEV
    assert ame.nuclear_mass_mev(table, 8, 16) == expect
    assert (ame.U_MEV, ame.ELECTRON_MASS_MEV) == (931.49410242, 0.51099895)
    with pytest.raises(KeyError):
        ame.nuclear_mass_mev(table, 9, 99)
    with pytest.raises(ame.AmeError, match="no AME2020"):
        ame.load_ame2020(AME_HEADER)
    with pytest.raises(ame.AmeError, match="duplicate"):
        ame.load_ame2020(AME_TEXT + ame_line("0", 8, 8, 16, "O", "1.0", "1.0", 15, "9.0", "1.0") + "\n")
    with pytest.raises(ame.AmeError, match="cannot parse"):
        ame.load_ame2020(AME_TEXT + "garbage line\n" + "x" * 130 + "\n")


# --------------------------------------------------------------------------------------------
# Synthetic EXFOR
# --------------------------------------------------------------------------------------------


def xl(key: str, text: str = "", ptr: str = " ") -> str:
    return f"{key:<10}{ptr}{text}".ljust(66)


def xrow(*vals: str | None) -> str:
    return "".join((v or "").rjust(11) for v in vals).rstrip()


def xhead(*vals: str) -> str:
    return "".join(v.ljust(11) for v in vals).rstrip()


EXFOR_TEXT = "\n".join([
    xl("ENTRY", "      X0001   19990101").replace("ENTRY     ", "ENTRY     ", 1),
    xl("SUBENT", "     X0001001   19990101"),
    xl("BIB", "         3          5"),
    xl("TITLE", "Synthetic cross sections for"),
    xl("", "an invented target"),
    xl("AUTHOR", "(A.Smith,B.Jones,"),
    xl("", "C.Brown)"),
    xl("REFERENCE", "(J,XX/Y,1,2,1999)"),
    xl("ENDBIB", "         5          0"),
    xl("NOCOMMON", "         0          0"),
    xl("ENDSUBENT", "         9          0"),
    xl("SUBENT", "     X0001002   19990101"),
    xl("BIB", "         1          1"),
    xl("REACTION", "(6-C-12(P,NON),,SIG) invented"),
    xl("ENDBIB", "         1          0"),
    xl("COMMON", "         1          3"),
    xl("ERR-SYS"),
    xl("PER-CENT"),
    xl("2.5"),
    xl("ENDCOMMON", "         3          0"),
    xl("DATA", "         3          3"),
    xhead("EN", "DATA", "ERR-S"),
    xhead("MEV", "MB", "MB"),
    xrow("100.", "275.", "21."),
    xrow("150.", "1.5+2", None),
    xrow("2.E+2", None, "3."),
    xl("ENDDATA", "         5          0"),
    xl("ENDSUBENT", "         9          0"),
    xl("SUBENT", "     X0001003   19990101"),
    xl("BIB", "         1          2"),
    xl("REACTION", "(8-O-16(P,NON),,SIG)", ptr="1"),
    xl("", "(8-O-16(P,NON),,SIG) cont", ptr="2"),
    xl("ENDBIB", "         2          0"),
    xl("NOCOMMON", "         0          0"),
    xl("DATA", "         2          2"),
    "EN".ljust(11) + "DATA".ljust(10) + "1",
    xhead("MEV/A", "B"),
    xrow("10.", "0.25"),
    xrow("20.", "0.5"),
    xl("ENDDATA", "         4          0"),
    xl("ENDSUBENT", "         9          0"),
    xl("ENDENTRY", "         3          0"),
])


def test_exfor_synthetic_entry() -> None:
    entry = exfor.parse_entry(EXFOR_TEXT)
    assert entry.entry_id == "X0001"
    assert [s.subentry_id for s in entry.subentries] == ["X0001001", "X0001002", "X0001003"]
    assert entry.first_author == "A.Smith"
    assert entry.reference_year == 1999
    assert entry.title == "Synthetic cross sections for an invented target"
    s1, s2, s3 = entry.subentries
    assert s1.common is None and s1.data is None
    assert s2.reactions == (exfor.Reaction("", "(6-C-12(P,NON),,SIG) invented"),)
    assert s2.common == exfor.CommonBlock(("ERR-SYS",), ("PER-CENT",), (2.5,))
    assert s2.data is not None
    assert s2.data.heads == ("EN", "DATA", "ERR-S")
    assert s2.data.units == ("MEV", "MB", "MB")
    assert s2.data.rows == ((100.0, 275.0, 21.0), (150.0, 150.0, None), (200.0, None, 3.0))
    col = s2.data.column(2)
    assert col[0] == 21.0 and np.isnan(col[1]) and col[2] == 3.0
    # pointers on REACTION and on a head field; MEV/A and B units
    assert [(r.pointer, r.text) for r in s3.reactions] == [
        ("1", "(8-O-16(P,NON),,SIG)"), ("2", "(8-O-16(P,NON),,SIG) cont"),
    ]
    assert s3.data is not None
    assert s3.data.units == ("MEV/A", "B") and s3.data.pointers == ("", "1")


def test_exfor_unit_helpers() -> None:
    e, per_a = exfor.energy_to_mev(np.array([100.0, 250.0]), "KEV")
    assert e.tolist() == [0.1, 0.25] and not per_a
    e, per_a = exfor.energy_to_mev(np.array([10.0]), "MEV/A")
    assert e.tolist() == [10.0] and per_a
    assert exfor.energy_total_mev(np.array([10.0]), "MEV/A", 12).tolist() == [120.0]
    assert exfor.energy_total_mev(np.array([10.0]), "MEV", None).tolist() == [10.0]
    with pytest.raises(exfor.ExforError):
        exfor.energy_total_mev(np.array([10.0]), "MEV/A", None)
    with pytest.raises(exfor.ExforError):
        exfor.energy_to_mev(np.array([1.0]), "GEV")
    assert exfor.xs_to_mb(np.array([0.25, 0.5]), "B").tolist() == [250.0, 500.0]
    assert exfor.xs_to_mb(np.array([275.0]), "MB").tolist() == [275.0]
    assert exfor.percent_to_fraction(np.array([2.5]), "PER-CENT").tolist() == [0.025]
    with pytest.raises(exfor.ExforError):
        exfor.percent_to_fraction(np.array([2.5]), "MB")
    assert exfor.parse_exfor_number("1.5+3") == 1500.0
    assert exfor.parse_exfor_number("  -2.E-2 ") == -0.02
    assert exfor.parse_exfor_number("     ") is None
    with pytest.raises(exfor.ExforError):
        exfor.parse_exfor_number("1.2.3")


def test_exfor_malformed_text_fails_closed() -> None:
    with pytest.raises(exfor.ExforError, match="ENTRY"):
        exfor.parse_entry("SUBENT X 1\n")
    with pytest.raises(exfor.ExforError, match="ENDSUBENT"):
        exfor.parse_entry("\n".join(EXFOR_TEXT.splitlines()[:14]))
    wide = EXFOR_TEXT.replace(xl("DATA", "         3          3"), xl("DATA", "         7          3"))
    with pytest.raises(exfor.ExforError, match="six fields"):
        exfor.parse_entry(wide)


# --------------------------------------------------------------------------------------------
# geant-val JSON
# --------------------------------------------------------------------------------------------


def _gv_record(rid: int, beam: str) -> dict[str, object]:
    return {
        "id": rid,
        "metadata": {"observableName": "inelastic cross section", "targetName": "X",
                     "beamParticle": beam},
        "chart": {"xAxisName": "E (MeV)", "yAxisName": "cross section, barn",
                  "xValues": [1.0, 2.0], "yValues": [0.5, 0.75],
                  "yStatErrorsPlus": [0.0, 0.0], "yStatErrorsMinus": [0.0, 0.0],
                  "ySysErrorsPlus": [0.125, 0.25], "ySysErrorsMinus": [0.0625, 0.5]},
    }


def test_geant_val_synthetic() -> None:
    curves = geant_val.parse_geant_val(json.dumps([_gv_record(1, "proton"), _gv_record(2, "C12")]))
    p, c = curves
    assert (p.record_id, p.target, p.beam, p.observable) == (1, "X", "proton", "inelastic cross section")
    assert p.x.tolist() == [1.0, 2.0] and p.y.tolist() == [0.5, 0.75]
    assert p.y_sys_plus.tolist() == [0.125, 0.25] and p.y_sys_minus.tolist() == [0.0625, 0.5]
    assert p.energy_caveat is None
    assert c.beam == "C12" and c.energy_caveat is not None and "per-nucleon" in c.energy_caveat
    bad = _gv_record(3, "proton")
    bad["chart"]["yValues"] = [1.0]  # type: ignore[index]
    with pytest.raises(geant_val.GeantValError, match="lengths"):
        geant_val.parse_geant_val(json.dumps([bad]))
    with pytest.raises(geant_val.GeantValError):
        geant_val.parse_geant_val(json.dumps({"not": "a list"}))
    with pytest.raises(geant_val.GeantValError, match="malformed"):
        geant_val.parse_geant_val(json.dumps([{"id": 1}]))


# --------------------------------------------------------------------------------------------
# Data-backed tests (cache; skipped when absent, failing under IONMC_REQUIRE_DATA=1)
# --------------------------------------------------------------------------------------------


def _cached(dataset_id: str) -> Path:
    """Verified cache object of ``dataset_id``; skip (or fail under IONMC_REQUIRE_DATA=1) if it
    is not cached. Integrity failures propagate and fail the test."""
    try:
        return cache.verify(dataset_id)
    except FileNotFoundError:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail(f"IONMC_REQUIRE_DATA=1 but dataset {dataset_id!r} is not cached")
        pytest.skip(f"dataset {dataset_id!r} is not cached")


def _member(zp: Path, element: str) -> endf6.EndfMaterial:
    return endf6.parse_endf(endf6.read_member(zp, f"ENDF-B-VIII.0_protons/p-{element}.endf"))


def test_real_endf_o16_structure_and_cross_section() -> None:
    zp = _cached("endf-b8.0-protons")
    members = endf6.list_zip_members(zp)
    for tag in ("p-001_H_001", "p-006_C_012", "p-007_N_014", "p-008_O_016", "p-013_Al_027",
                "p-014_Si_028", "p-015_P_031", "p-020_Ca_040"):
        assert f"ENDF-B-VIII.0_protons/{tag}.endf" in members
    o16 = _member(zp, "008_O_016")
    assert o16.za == 8016.0 and o16.emax == 1.5e8
    assert o16.cross_section(5).interpolate(1.0e8)[0] == pytest.approx(0.297, rel=0.02)
    sec = o16.products(5)
    assert sec.lct == 3 and len(sec.products) == 30
    light = [p for p in sec.products if p.zap in (1, 1001, 1002, 2004)]
    assert sorted(p.zap for p in light) == [1, 1001, 1002, 2004]
    for p in light:
        assert p.law == 1 and p.lang == 2
        assert p.distributions and all(d.na == 1 and d.rows.shape[1] == 3 for d in p.distributions)
    # the elastic section carries LAW=5, parsed opaque
    assert [p.law for p in o16.products(2).products] == [5]


def test_real_endf_c12_cross_section_and_h1() -> None:
    zp = _cached("endf-b8.0-protons")
    c12 = _member(zp, "006_C_012")
    assert c12.cross_section(5).interpolate(1.0e8)[0] == pytest.approx(0.227, rel=0.02)
    assert c12.products(5).lct == 3
    h1 = _member(zp, "001_H_001")
    assert (3, 5) not in h1.sections and (6, 5) not in h1.sections
    assert (3, 2) in h1.sections


def test_real_endf_targets_of_interest_parse() -> None:
    zp = _cached("endf-b8.0-protons")
    for tag in ("007_N_014", "013_Al_027", "014_Si_028", "015_P_031", "020_Ca_040"):
        mat = _member(zp, tag)
        tab = mat.cross_section(5)
        assert tab.x[-1] == mat.emax == 1.5e8
        assert 0.0 < tab.interpolate(1.0e8)[0] < 2.0
        assert all(p.law == 1 for p in mat.products(5).products)


def test_real_ame_o16() -> None:
    table = ame.load_ame2020(_cached("ame2020-mass").read_text(encoding="ascii"))
    o16 = table[(8, 16)]
    assert o16.mass_excess_kev == pytest.approx(-4737.00217, abs=1e-5)
    assert not o16.estimated
    assert any(e.estimated for e in table.values())
    assert 14895.0 < ame.nuclear_mass_mev(table, 8, 16) < 14895.2
    assert ame.nuclear_mass_mev(table, 1, 1) == pytest.approx(938.272, abs=1e-3)


def test_real_exfor_d0356_and_c1862() -> None:
    d = exfor.parse_entry(_cached("exfor-d0356").read_text(encoding="ascii"))
    assert d.entry_id == "D0356" and d.first_author == "A.Auce" and d.reference_year == 2005
    carbon = [s for s in d.subentries if any("6-C-12(P,NON)" in r.text for r in s.reactions)]
    assert len(carbon) == 1 and carbon[0].data is not None
    assert len(carbon[0].data.rows) == 6 and carbon[0].data.units == ("MEV", "MB", "MB")
    assert carbon[0].common is not None and carbon[0].common.units == ("PER-CENT",)
    c = exfor.parse_entry(_cached("exfor-c1862").read_text(encoding="ascii"))
    assert c.entry_id == "C1862" and c.reference_year == 1975
    oxygen = [s for s in c.subentries if any("8-O-16(P,NON)" in r.text for r in s.reactions)]
    assert len(oxygen) == 1 and oxygen[0].data is not None and oxygen[0].data.rows


def test_real_geant_val() -> None:
    curves = geant_val.parse_geant_val(_cached("geant-val-exfor-inelastic-7").read_text())
    assert len(curves) == 6
    assert {c.target for c in curves} == {"Al", "C", "Ca", "O"}
    flagged = [c for c in curves if c.energy_caveat]
    assert sorted(c.beam for c in flagged) == ["C12", "C12"]
    assert all(c.x.size == c.y.size > 0 for c in curves)
