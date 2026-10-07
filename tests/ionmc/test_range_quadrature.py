"""V3-003D acceptance rows D1, D1b, D2, D2b, D3 and D4 (plan validation/plans/v3-003d-acceptance.md,
criteria frozen before any result): the CSDA range is the exact integral of the log-log
interpolated stopping power (decision 0039, ``range_construction = exact-loglog-quadrature-v1``).

The negative controls replace ``ionmc.physics.stopping.exact_loglog_range_increments`` by the
trapezoid rule in ln E of ``a E / S`` (the construction before V3-003D) for the duration of the
table build; the tables are data, so one patch covers both backends. Analytic Bethe water unless
stated otherwise, ``E_cut`` = 2 MeV. CI tier, fixed small seeds (the runs are deterministic).
"""

from __future__ import annotations

import math
import os
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from unittest import mock

import numpy as np
import pytest

from ionmc.config import DiagnosticsOptions, SimulationConfig
from ionmc.materials import ALUMINIUM, WATER, Material
from ionmc.physics import stopping
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource, StoppingSource, StoppingTable
from ionmc.simulation import Simulation
from ionmc.transport.tables import TransportTables
from tests.ionmc.test_scoring_transport import (
    A4B_GRID,
    _a4b_worst,
    _csda_run,
    _cut_depth,
    _reqs,
    _slab,
)

MakeConfig = Callable[..., SimulationConfig]
E_CUT = 2.0
S_MAX = (1.0, 0.5, 0.25, 0.1)

_GL_X, _GL_W = np.polynomial.legendre.leggauss(32)


def _trapezoid_increments(
    ln_e: np.ndarray, ln_s: np.ndarray, a: float
) -> np.ndarray:  # the construction before V3-003D (control)
    f = a * np.exp(ln_e - ln_s)
    return 0.5 * (f[1:] + f[:-1]) * np.diff(ln_e)


@contextmanager
def trapezoid_tables() -> Iterator[None]:
    with mock.patch.object(stopping, "exact_loglog_range_increments", _trapezoid_increments):
        yield


def _gl_increment(t: TransportTables, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """32-point Gauss-Legendre quadrature in ``u = ln E`` of ``E / S_interp(E)`` from ``lo`` to
    ``hi`` (arrays of ln E inside one bin each); ``S`` from ``stopping_mass`` (independent of the
    closed form)."""
    out = np.empty(lo.size)
    for k in range(lo.size):
        u = 0.5 * (hi[k] - lo[k]) * _GL_X + 0.5 * (hi[k] + lo[k])
        s = np.array([t.stopping_mass(0, math.exp(x)) for x in u])
        out[k] = 0.5 * (hi[k] - lo[k]) * float(np.sum(_GL_W * np.exp(u) / s))
    return out


def _grid(t: TransportTables) -> np.ndarray:
    return t.ln_e0[0] + np.arange(t.n_e) / t.inv_dln_e[0]


def _nist_source() -> StoppingSource | None:
    """The NIST PSTAR water source from the local cache, or None when the dataset is not cached.

    Only the cache-missing condition (:class:`~ionmc.data.acquire.OfflineError`, offline fetch of
    an object that is not cached) returns None; checksum, parsing, registry and loader failures
    propagate and fail the test. ``IONMC_REQUIRE_NIST=1`` (set by the ``lv`` suite) makes a missing
    cache a failure, so the NIST-water case of D1 cannot be skipped in the controlled run."""
    from ionmc.data.acquire import OfflineError, fetch
    from ionmc.data.nist_star import load_star_table

    try:
        star = load_star_table(fetch("nist-pstar-water-2005", None, offline=True))
    except OfflineError:
        if os.environ.get("IONMC_REQUIRE_NIST") == "1":
            pytest.fail("IONMC_REQUIRE_NIST=1 but the NIST PSTAR water table is not cached")
        return None
    return stopping.NistStarStoppingSource(star)


def _source_tables(name: str) -> tuple[StoppingTable, TransportTables]:
    if name == "nist-water":
        src = _nist_source()
        if src is None:
            pytest.skip("NIST PSTAR water table not cached (D1 on NIST water runs in LV)")
        table = src.table(WATER, PROTON)
    else:
        material: Material = WATER if name == "water" else ALUMINIUM
        table = BetheStoppingSource().table(material, PROTON)
    return table, TransportTables.from_stopping_tables([table])


# --------------------------------------------------------------------------- D1


@pytest.mark.parametrize("name", ["water", "aluminium", "nist-water"])
def test_d1_table_equals_independent_quadrature(name: str) -> None:
    table, t = _source_tables(name)
    grid = _grid(t)
    # (a) every nodal increment equals the Gauss-Legendre quadrature of E / S_interp
    inc = np.diff(t.r_mass[0])
    gl = _gl_increment(t, grid[:-1], grid[1:])
    worst_a = float(np.max(np.abs(inc / gl - 1.0)))
    # (b) range_g_cm2 at 1e4 log-uniform energies = R_0 + quadrature from E_min
    rng = np.random.default_rng(20390001)
    e = np.exp(rng.uniform(grid[0], grid[-1], 10_000))
    ln_e = np.log(e)
    i = np.clip(np.floor((ln_e - grid[0]) * t.inv_dln_e[0]).astype(int), 0, t.n_e - 2)
    cum = np.concatenate(([0.0], np.cumsum(gl)))
    part = _gl_increment(t, grid[i], ln_e)
    ref = t.r_min_g_cm2[0] + cum[i] + part
    got = np.array([t.range_g_cm2(0, x) for x in e])
    worst_b = float(np.max(np.abs(got / ref - 1.0)))
    # (c) the StoppingTable nodes are the TransportTables nodes (same grid)
    assert table.csda_range_g_cm2.size == t.n_e
    worst_c = float(np.max(np.abs(table.csda_range_g_cm2 / t.r_mass[0] - 1.0)))
    print(f"D1 {name}: (a) {worst_a:.2e} (b) {worst_b:.2e} (c) {worst_c:.2e}")
    assert worst_a <= 1e-12
    assert worst_b <= 1e-12
    assert worst_c <= 1e-12
    assert table.metadata["range_construction"] == "exact-loglog-quadrature-v1"
    assert t.identity[0]["range_construction"] == "exact-loglog-quadrature-v1"
    # (d) negative control: the trapezoid increments are off by >= 1e-5 in >= 90 % of the intervals
    with trapezoid_tables():
        _, tc = _source_tables(name)
    dev = np.abs(np.diff(tc.r_mass[0]) / gl - 1.0)
    frac = float(np.mean(dev >= 1e-5))
    print(f"D1 {name}: trapezoid control, fraction of intervals off by >= 1e-5: {frac:.3f}")
    assert frac >= 0.9
    assert tc.sha256 != t.sha256


# --------------------------------------------------------------------------- D1b


def test_d1b_round_trip_and_monotonicity() -> None:
    _, t = _source_tables("water")
    energies = np.geomspace(t.e_min_mev[0], t.e_max_mev[0], 20001)
    ranges = np.array([t.range_g_cm2(0, e) for e in energies])
    assert np.all(np.diff(ranges) > 0.0)
    back = np.array([t.energy_from_range(0, r) for r in ranges])
    assert np.all(np.diff(back) > 0.0)
    eta_e = float(np.max(np.abs(back / energies - 1.0)))
    r_grid = np.geomspace(t.r_min_g_cm2[0], t.r_max_g_cm2[0], 20001)
    e_of_r = np.array([t.energy_from_range(0, r) for r in r_grid])
    assert np.all(np.diff(e_of_r) > 0.0)
    eta_r = float(np.max(np.abs(np.array([t.range_g_cm2(0, e) for e in e_of_r]) / r_grid - 1.0)))
    print(f"D1b: max |Rinv(R(E))/E - 1| = {eta_e:.2e}, max |R(Rinv(r))/r - 1| = {eta_r:.2e}")
    assert eta_e <= 1e-6
    assert eta_r <= 2e-6
    assert eta_e <= 1e-5 and eta_r <= 1e-5  # U2 as frozen


# --------------------------------------------------------------------------- D2 / D2b / D4

_CACHE: dict[tuple[str, float, bool], dict[str, float]] = {}


def _deterministic_run(
    make_config: MakeConfig, backend: str, s_max: float, control: bool
) -> dict[str, float]:
    """Boundary consistency and projected end depth of the deterministic 150 MeV history."""
    key = (backend, s_max, control)
    if key in _CACHE:
        return _CACHE[key]
    nz = int(round(170.0 / s_max))
    t0 = time.perf_counter()
    cfg = make_config(
        energy=150.0, n=2, n_batches=2, straggling=False, mcs=False, max_step=s_max,
        geometry=_slab(s_max, nz), position=(1.0, 1.0, 0.0), backend=backend,
        diagnostics=DiagnosticsOptions(trace_histories=1),
    )  # fmt: skip
    if control:
        with trapezoid_tables():
            res = Simulation(cfg).run()
    else:
        res = Simulation(cfg).run()
    t = res.effective_config.tables
    tr = res.diagnostics["trace"]
    sel = (tr["history"] == 0) & (tr["energy_mev"] > E_CUT)
    z_b, e_b = tr["z_mm"][sel], tr["energy_mev"][sel]
    r_b = np.array([t.range_g_cm2(0, e) for e in e_b]) * 10.0  # mm of water (rho = 1)
    resid = np.abs(r_b - (t.range_g_cm2(0, 150.0) * 10.0 - z_b))
    z_end = float(z_b[-1] + r_b[-1] - t.range_g_cm2(0, E_CUT) * 10.0)
    out = {
        "max_resid_um": float(np.max(resid)) * 1000.0,
        "z_end_mm": z_end,
        "r150_mm": t.range_g_cm2(0, 150.0) * 10.0,
        "r100_mm": t.range_g_cm2(0, 100.0) * 10.0,
        "r200_mm": t.range_g_cm2(0, 200.0) * 10.0,
        "seconds": time.perf_counter() - t0,
    }
    _CACHE[key] = out
    return out


@pytest.mark.parametrize("s_max", S_MAX)
@pytest.mark.parametrize("backend", ["python", "warp-cpu"])
def test_d2_branch_consistency(make_config: MakeConfig, backend: str, s_max: float) -> None:
    new = _deterministic_run(make_config, backend, s_max, control=False)
    print(f"D2 {backend} s_max {s_max}: max |R(E_b) - (R(E0) - z_b)| = {new['max_resid_um']:.3f} um"
          f" ({new['seconds']:.1f} s)")  # fmt: skip
    assert new["max_resid_um"] <= 0.5
    if backend == "python":
        ctl = _deterministic_run(make_config, backend, s_max, control=True)
        print(f"D2 control (trapezoid) s_max {s_max}: {ctl['max_resid_um']:.3f} um")
        assert ctl["max_resid_um"] >= 1.0


def test_d2b_end_depth_independent_of_step_size(make_config: MakeConfig) -> None:
    ends = [_deterministic_run(make_config, "python", s, False)["z_end_mm"] for s in S_MAX]
    ctrl = [_deterministic_run(make_config, "python", s, True)["z_end_mm"] for s in S_MAX]
    spread = (max(ends) - min(ends)) * 1000.0
    spread_c = (max(ctrl) - min(ctrl)) * 1000.0
    print(f"D2b: end-depth spread {spread:.3f} um (trapezoid control {spread_c:.3f} um)")
    assert spread <= 0.5
    assert spread_c >= 2.0


# --------------------------------------------------------------------------- D3


@pytest.mark.parametrize("max_step", [1.0, 0.5, 0.25])
def test_d3_a4b_against_the_table_csda_energy(make_config: MakeConfig, max_step: float) -> None:
    r = _csda_run(make_config, _reqs("let_t"), 150.0, max_step, A4B_GRID)
    z_cut = _cut_depth(r)
    worst_table = _a4b_worst(r, z_cut, 100.0, 1e9, table_csda=True)
    worst_transported = _a4b_worst(r, z_cut, 100.0, 1e9)
    print(f"D3 s_max {max_step}: A4b vs table CSDA {worst_table:.2e}, vs transported "
          f"{worst_transported:.2e}")  # fmt: skip
    assert worst_table <= 1e-4
    assert worst_transported <= 1e-4
    with trapezoid_tables():
        rc = _csda_run(make_config, _reqs("let_t"), 150.0, max_step, A4B_GRID)
        worst_c = _a4b_worst(rc, _cut_depth(rc), 100.0, 1e9, table_csda=True)
    print(f"D3 control (trapezoid) s_max {max_step}: {worst_c:.2e}")
    assert worst_c > 1e-4


# --------------------------------------------------------------------------- D4

# (deterministic transported end-depth shift [um] predicted from the mock, plan D4)
_END_SHIFT_UM = {1.0: -3.41, 0.5: -1.74, 0.25: -0.88, 0.1: -0.35}


def test_d4_range_shift_against_the_trapezoid_construction(make_config: MakeConfig) -> None:
    new = _deterministic_run(make_config, "python", 1.0, False)
    old = _deterministic_run(make_config, "python", 1.0, True)
    print(f"D4: R(150 MeV) new {new['r150_mm']:.5f} mm, trapezoid {old['r150_mm']:.5f} mm, "
          f"shift {(new['r150_mm'] - old['r150_mm']) * 1000.0:.2f} um")  # fmt: skip
    assert new["r150_mm"] == pytest.approx(158.6248, abs=0.0002)
    for k in ("r100_mm", "r200_mm"):
        shift = (new[k] - old[k]) * 1000.0
        print(f"D4: {k} shift {shift:.2f} um")
        assert abs(shift) <= 10.0
    for s_max, predicted in _END_SHIFT_UM.items():
        n = _deterministic_run(make_config, "python", s_max, False)["z_end_mm"]
        o = _deterministic_run(make_config, "python", s_max, True)["z_end_mm"]
        shift = (n - o) * 1000.0
        print(f"D4: end-depth shift at s_max {s_max}: {shift:.2f} um (predicted {predicted})")
        assert abs(shift - predicted) <= 0.5


def test_range_construction_identity_enters_the_hash() -> None:
    _, t = _source_tables("water")
    assert t.identity[0]["range_construction"] == stopping.RANGE_CONSTRUCTION


# --------------------------------------------------------------------------- NIST loader policy


def test_nist_source_missing_cache_skips_or_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Only a missing cache is tolerated (None), and ``IONMC_REQUIRE_NIST=1`` makes it a failure;
    every other loader error propagates."""
    from ionmc.data import acquire

    def missing(*_a: object, **_k: object) -> None:
        raise acquire.OfflineError("not cached")

    monkeypatch.setattr(acquire, "fetch", missing)
    monkeypatch.delenv("IONMC_REQUIRE_NIST", raising=False)
    assert _nist_source() is None
    monkeypatch.setenv("IONMC_REQUIRE_NIST", "1")
    with pytest.raises(pytest.fail.Exception, match="IONMC_REQUIRE_NIST"):
        _nist_source()

    from ionmc.data.cache import IntegrityError

    def corrupt(*_a: object, **_k: object) -> None:
        raise IntegrityError("checksum")

    monkeypatch.setattr(acquire, "fetch", corrupt)
    monkeypatch.delenv("IONMC_REQUIRE_NIST", raising=False)
    with pytest.raises(IntegrityError):
        _nist_source()
