"""Tests of the Stage-1 CSDA depth-dose transport (decision 0009).

Reference Python path (float64 oracle); the Warp paths are covered by
``@pytest.mark.warp`` here and by ``validation/v1_depth_dose_csda.py`` on the
host runner. The NIST/tabulated CSDA range is the physics reference.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc import materials, particles
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import (
    DepthDoseGrid,
    ParticleState,
    PencilBeamSource,
    TransportEngine,
    WaterSlab,
)
from ionmc.transport.state import Species, Status

ENERGIES = (100.0, 150.0, 200.0)
R80_TOL = 0.003  # decision 0009


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


@pytest.fixture(scope="module")
def engine(table):
    # deterministic (no straggling): the DEV-004 CSDA foundation tests
    return TransportEngine(
        table,
        WaterSlab(400.0, materials.WATER),
        DepthDoseGrid(400.0, 800),
        straggling=False,
    )


def _csda_range_mm(table, e0: float) -> float:
    tab = TabulatedStoppingPower(
        table, materials.WATER, particles.PROTON, path="python"
    )
    return float(tab.csda_range(e0)[0]) * 10.0 / materials.WATER.density_g_per_cm3


def test_state_layout_carries_identity_fields() -> None:
    state = ParticleState.allocate(4)
    assert state.size == 4
    assert state.species.dtype == np.int32 and np.all(
        state.species == int(Species.PROTON)
    )
    assert state.beamlet.shape == (4,) and state.rng_state.dtype == np.uint32
    assert np.all(state.status == int(Status.ALIVE))
    assert state.position_mm.shape == (4, 3) and state.direction.shape == (4, 3)
    np.testing.assert_array_equal(state.direction[0], [0.0, 0.0, 1.0])


def test_source_makes_distinct_per_history_rng_streams() -> None:
    state = PencilBeamSource(150.0).sample(8, seed=7)
    assert np.all(state.energy_mev == 150.0)
    assert len(set(state.rng_state.tolist())) == 8  # distinct streams
    with pytest.raises(ValueError):
        PencilBeamSource(-1.0)


def test_energy_is_conserved_on_reference_path(engine) -> None:
    res = engine.run(PencilBeamSource(150.0), n_histories=3, path="python")
    assert res.truncated == 0
    assert abs(res.energy_balance) <= 1e-9
    # energy deposited equals the energy started (all of it stops in the slab)
    assert res.energy_deposited_mev == pytest.approx(150.0 * 3, rel=1e-12)


@pytest.mark.parametrize("e0", ENERGIES)
def test_r80_matches_csda_range(engine, table, e0) -> None:
    res = engine.run(PencilBeamSource(e0), n_histories=1, path="python")
    r80 = res.r80_mm()
    r_csda = _csda_range_mm(table, e0)
    assert abs(r80 / r_csda - 1.0) <= R80_TOL, (e0, r80, r_csda)


def test_depth_dose_scales_linearly_with_histories(engine) -> None:
    one = engine.run(PencilBeamSource(120.0), n_histories=1, path="python")
    five = engine.run(PencilBeamSource(120.0), n_histories=5, path="python")
    np.testing.assert_allclose(five.edep_mev, 5.0 * one.edep_mev, rtol=1e-12)


def test_step_size_convergence_on_fixed_grid(table) -> None:
    # Step convergence is measured on the mean stopping depth (the transport
    # physics), which is grid-independent; R80 additionally depends on how a
    # step's energy is smeared across bins once deposition is distributed along
    # the step, so it is not the right step-convergence metric here.
    grid = DepthDoseGrid(400.0, 8000)
    coarse = TransportEngine(
        table,
        WaterSlab(400.0, materials.WATER),
        grid,
        max_fraction=0.04,
        straggling=False,
    )
    fine = TransportEngine(
        table,
        WaterSlab(400.0, materials.WATER),
        grid,
        max_fraction=0.002,
        straggling=False,
    )
    for e0 in ENERGIES:
        r_c = coarse.run(PencilBeamSource(e0), 1, path="python").range_mean_mm
        r_f = fine.run(PencilBeamSource(e0), 1, path="python").range_mean_mm
        assert abs(r_c / r_f - 1.0) <= 5e-4, (e0, r_c, r_f)


def test_range_cross_check_against_csda_table(engine, table) -> None:
    # the on-axis path length (stopping depth) equals the CSDA range from the
    # cut energy, i.e. R(E0) measured from the 0.5 MeV table floor.
    for e0 in ENERGIES:
        res = engine.run(PencilBeamSource(e0), 1, path="python")
        # last non-empty bin centre is the stopping depth to within a bin
        last = np.max(np.nonzero(res.edep_mev)[0])
        stop_depth = res.grid.centers_mm[last]
        r_csda = _csda_range_mm(table, e0)
        assert abs(stop_depth / r_csda - 1.0) <= 0.002, (e0, stop_depth, r_csda)


def test_escaped_history_when_slab_too_thin(table) -> None:
    thin = TransportEngine(
        table, WaterSlab(50.0, materials.WATER), DepthDoseGrid(50.0, 200)
    )
    res = thin.run(PencilBeamSource(150.0), n_histories=1, path="python")
    # a 150 MeV proton ranges ~158 mm, so it leaves a 50 mm slab: only part of
    # its energy is deposited.
    assert res.energy_deposited_mev < 150.0
    assert res.energy_balance < 0.0


def test_scoring_grid_shorter_than_slab_does_not_stop_transport(table) -> None:
    # Scoring grid (100 mm) is shorter than the slab (400 mm); a 150 MeV proton
    # (~158 mm range) must transport its full range in the medium and not be
    # terminated at the grid edge - only the first 100 mm is scored.
    slab = WaterSlab(400.0, materials.WATER)
    short_grid = DepthDoseGrid(100.0, 200)
    eng = TransportEngine(table, slab, short_grid, straggling=False)
    res = eng.run(PencilBeamSource(150.0), n_histories=1, path="python")
    # energy scored is only the part deposited within the first 100 mm
    assert 0.0 < res.energy_deposited_mev < 150.0
    # the proton stops in the slab (does not escape a 400 mm medium)
    full_grid = DepthDoseGrid(400.0, 800)
    full = TransportEngine(table, slab, full_grid, straggling=False).run(
        PencilBeamSource(150.0), 1, path="python"
    )
    assert abs(full.energy_balance) <= 1e-9  # full grid: all energy scored


def test_scoring_grid_longer_than_slab_marks_escape(table) -> None:
    # Slab (100 mm) shorter than the grid (400 mm): the proton leaves the medium
    # at 100 mm with residual energy, so less than its full energy is deposited.
    eng = TransportEngine(
        table,
        WaterSlab(100.0, materials.WATER),
        DepthDoseGrid(400.0, 800),
        straggling=False,
    )
    res = eng.run(PencilBeamSource(150.0), n_histories=1, path="python")
    assert res.energy_deposited_mev < 150.0
    # nothing scored beyond the 100 mm slab
    beyond = res.grid.centers_mm > 100.0
    assert np.all(res.edep_mev[beyond] == 0.0)


def test_truncation_counter_increments_on_step_cap(table) -> None:
    # A tiny step cap forces the history to hit max_steps before stopping; it
    # must be counted, not silently dropped.
    eng = TransportEngine(
        table,
        WaterSlab(400.0, materials.WATER),
        DepthDoseGrid(400.0, 800),
        max_steps=5,
        straggling=False,
    )
    res = eng.run(PencilBeamSource(150.0), n_histories=2, path="python")
    assert res.truncated == 2
    # only a few steps' worth of energy was deposited before truncation
    assert res.energy_deposited_mev < 150.0 * 2


# -- Warp paths (sandbox: CPU; CUDA on the host runner) -----------------------


@pytest.mark.warp
@pytest.mark.parametrize("e0", ENERGIES)
def test_warp_cpu_matches_reference_cumulative(warp_module, engine, e0) -> None:
    ref = engine.run(PencilBeamSource(e0), 1, path="python")
    warp = engine.run(PencilBeamSource(e0), 1, path="warp", device="cpu")
    total = np.sum(ref.edep_mev)
    cum_diff = (
        np.max(np.abs(np.cumsum(warp.edep_mev) - np.cumsum(ref.edep_mev))) / total
    )
    assert cum_diff <= 1e-4, (e0, cum_diff)
    # float32 transport conserves energy to ~1e-5 (the de values and the
    # along-step deposition fractions are float32); the accumulator is float64.
    assert abs(warp.energy_balance) <= 5e-5
    # R80 within one bin
    assert abs(warp.r80_mm() - ref.r80_mm()) <= ref.grid.bin_width_mm
    # edge-aware per-bin normalized difference
    peak = ref.edep_mev.max()
    nd = np.max(
        np.abs(warp.edep_mev - ref.edep_mev)
        / (1e-3 * peak + 2e-3 * np.abs(ref.edep_mev))
    )
    assert nd <= 1.0, (e0, nd)


@pytest.mark.cuda
def test_warp_cuda_matches_cpu(warp_module, cuda_available, engine) -> None:
    cpu = engine.run(PencilBeamSource(150.0), 1, path="warp", device="cpu")
    cuda = engine.run(PencilBeamSource(150.0), 1, path="warp", device="cuda:0")
    total = np.sum(cpu.edep_mev)
    cum_diff = (
        np.max(np.abs(np.cumsum(cuda.edep_mev) - np.cumsum(cpu.edep_mev))) / total
    )
    assert cum_diff <= 1e-5, cum_diff
