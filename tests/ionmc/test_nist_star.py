"""Tests for NIST PSTAR/ASTAR access."""

from __future__ import annotations

import urllib.error
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import pytest

from ionmc.data import DataCache, DatasetSpec, nist_star

FIXTURE = """PSTAR: Stopping Powers and Range Tables for Protons

WATER, LIQUID

Kinetic   Electron. Nuclear   Total     CSDA      Projected Detour
Energy    Stp. Pow. Stp. Pow. Stp. Pow. Range     Range     Factor
MeV       MeV cm2/g MeV cm2/g MeV cm2/g g/cm2     g/cm2

1.000E-03 1.337E+02 4.315E+01 1.769E+02 6.319E-06 2.878E-06 0.4555
1.500E-03 1.638E+02 3.460E+01 1.984E+02 8.969E-06 4.400E-06 0.4906
1.000E+01 4.564E+01 1.0E-01 4.574E+01 1.230E-01 1.2E-01 0.9990
5.000E+01 1.239E+01 5.0E-03 1.240E+01 2.227E+00 2.2E+00 0.9994
8.000E+01 9.559E+00 4.0E-03 9.563E+00 5.184E+00 5.2E+00 0.9995
1.000E+02 7.286E+00 3.5E-03 7.290E+00 7.718E+00 7.7E+00 0.9995
"""


def _cache(tmp_path: Path) -> DataCache:
    return DataCache(tmp_path, download=lambda s: FIXTURE.encode())


def test_parse_skips_headers() -> None:
    t = nist_star.parse_star_table(FIXTURE, "PSTAR", "water_liquid")
    assert t.energy_mev.shape == (6,)
    assert t.energy_mev[0] == 1e-3
    assert t.electronic[-1] == 7.286
    assert t.nuclear[0] == 43.15
    assert t.total[1] == 198.4
    assert t.detour_factor[0] == 0.4555
    assert t.csda_range_g_cm2[-1] == pytest.approx(7.718)
    assert t.projected_range_g_cm2[-1] == pytest.approx(7.7)
    assert t.electronic.dtype == np.float64


def test_parse_empty_raises() -> None:
    with pytest.raises(ValueError):
        nist_star.parse_star_table("nothing here\n")


def test_astar_energy_per_nucleon() -> None:
    a = nist_star.parse_star_table(FIXTURE, "ASTAR", "water_liquid")
    p = nist_star.parse_star_table(FIXTURE, "PSTAR", "water_liquid")
    np.testing.assert_allclose(a.energy_per_nucleon_mev(), a.energy_mev / 4)
    np.testing.assert_array_equal(p.energy_per_nucleon_mev(), p.energy_mev)


def test_interpolation() -> None:
    t = nist_star.parse_star_table(FIXTURE, "PSTAR", "water_liquid")
    assert t.electronic_at(100.0) == pytest.approx(7.286)
    mid = t.electronic_at(np.sqrt(10.0 * 50.0))
    assert mid == pytest.approx(np.sqrt(45.64 * 12.39))
    with pytest.raises(ValueError):
        t.electronic_at(1000.0)


def test_spec() -> None:
    s = nist_star.spec("pstar", "water_liquid")
    assert s.dataset_id == "nist-pstar/water_liquid"
    assert s.version == "ICRU49"
    assert s.method == "POST"
    assert s.post_fields is not None
    assert s.post_fields["matno"] == "276"
    assert s.post_fields["prog"] == "PSTAR"
    assert "prog=PSTAR" in urlencode(s.post_fields)
    with pytest.raises(KeyError):
        nist_star.spec("PSTAR", "unobtainium")
    with pytest.raises(ValueError):
        nist_star.spec("XSTAR", "water_liquid")


def test_load_star_table_cached(tmp_path: Path) -> None:
    calls: list[DatasetSpec] = []

    def dl(s: DatasetSpec) -> bytes:
        calls.append(s)
        return FIXTURE.encode()

    cache = DataCache(tmp_path, download=dl)
    _, rec = nist_star.load_star_table(cache, "PSTAR", "water_liquid")
    assert rec.source == "network"
    t, rec2 = nist_star.load_star_table(cache, "PSTAR", "water_liquid", offline=True)
    assert rec2.source == "cache"
    assert len(calls) == 1
    assert t.material == "water_liquid"


@pytest.mark.local
def test_real_nist_pstar_water(tmp_path: Path) -> None:
    cache = DataCache(tmp_path)
    try:
        table, rec = nist_star.load_star_table(cache, "PSTAR", "water_liquid")
    except urllib.error.URLError as exc:
        pytest.skip(f"network unavailable: {exc}")
    assert rec.source == "network"
    assert table.electronic_at(100.0) == pytest.approx(7.286, rel=1e-3)
