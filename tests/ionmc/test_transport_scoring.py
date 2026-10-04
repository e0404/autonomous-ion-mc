"""Track-length scoring, fixed-point accumulators and the T9-CI step-independence check."""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np
import pytest

from ionmc.config import SimulationConfig, validate
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.transport.funcs import make_transport_funcs
from ionmc.transport.tally import MAX_QUANTA, QUANTUM_MEV, PartialTransport, merge_partials

MakeConfig = Callable[..., SimulationConfig]


def _idd(make_config: MakeConfig, s_max: float, bin_mm: float, n: int = 20000) -> np.ndarray:
    """IDD [MeV per primary] of 150 MeV protons in a single-voxel water box, straggling and MCS
    off (a deterministic path), bins of ``bin_mm``; warp-cpu float32."""
    depth = 1.3 * 158.6
    nz = int(math.ceil(depth / bin_mm))
    cfg = make_config(
        energy=150.0,
        n=n,
        n_batches=20,
        geometry=BoxPhantom((-50.0, -50.0, 0.0), (100.0, 100.0, depth), WATER),
        scoring=(ScoringGrid((-40.0, -40.0, 0.0), (80.0, 80.0, bin_mm), (1, 1, nz)),),
        mcs=False,
        straggling=False,
        max_step=s_max,
        backend="warp-cpu",
        precision="float32",
        run_kwargs={"cpu_workers": 4},
    )
    res = Simulation(cfg).run()
    assert res.valid and res.energy_balance.grid_relative_residual(0) < 1e-5
    return np.asarray(res.grids[0].energy_mev).reshape(-1)


@pytest.mark.parametrize(
    ("bin_mm", "steps"), [(1.0, (0.33, 0.9, 1.0, 0.5)), (2.0, (2.0, 1.0, 0.9))], ids=["1mm", "2mm"]
)
def test_t9_ci_idd_step_independence_deterministic(
    make_config: MakeConfig, bin_mm: float, steps: tuple[float, ...]
) -> None:
    """T9-CI: single-voxel water box, 150 MeV, straggling and MCS off, 2e4 histories: the IDD of
    every ``s_max`` agrees with the 0.1 mm-step run within 2e-3 (20-120 mm) and 1e-2
    (125-140 mm); for 2 mm bins ``s_max`` = 2.0 is included.

    Pre-fix contrary evidence (midpoint scoring, same setup, 1 mm bins, deviations from the 0.1 mm
    run): 62 % (125-140 mm) for s_max = 1.0, 32 % (20-120 mm) for 0.33, 80 % for 0.9; with 2 mm
    bins and s_max = 2.0 the plateau deviated by 66 %. The aliasing of point deposits with the bin
    edges is removed by apportioning each step along its path."""
    ref = _idd(make_config, 0.1, bin_mm)
    z = (np.arange(len(ref)) + 0.5) * bin_mm
    keep = ref > 0.01 * ref.max()
    for s in steps:
        a = _idd(make_config, s, bin_mm)
        dev = np.abs(a / np.where(ref > 0, ref, 1.0) - 1.0)
        plateau = keep & (z >= 20.0) & (z <= 120.0)
        peak = keep & (z >= 125.0) & (z <= 140.0)
        assert dev[plateau].max() <= 2e-3, (s, bin_mm, dev[plateau].max())
        assert dev[peak].max() <= 1e-2, (s, bin_mm, dev[peak].max())


def test_deposits_are_apportioned_over_voxels_and_close_the_balance(
    make_config: MakeConfig,
) -> None:
    """Deterministic 20 MeV oblique beam through 1 mm scoring voxels with 1 mm steps: the grid
    sum plus the rounding
    residual plus outside deposits equals the deposited energy to 1e-12 (float64) on both
    backends, and the grids are integers times the quantum."""
    results = {}
    for backend in ("python", "warp-cpu"):
        cfg = make_config(
            energy=20.0,
            n=2,
            n_batches=2,
            position=(0.37, 0.11, 0.0),
            direction=(0.3, 0.2, 1.0),
            geometry=BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 12.0), WATER),
            scoring=(ScoringGrid((-10.0, -10.0, 0.0), (1.0, 1.0, 1.0), (20, 20, 12)),),
            mcs=False,
            straggling=False,
            max_step=1.0,
            backend=backend,
        )
        res = Simulation(cfg).run()
        eb = res.energy_balance
        assert eb.grid_relative_residual(0) <= 1e-12 and eb.relative_residual <= 1e-12
        assert abs(eb.quantization_mev[0]) <= 0.5 * QUANTUM_MEV * 1e5
        quanta = np.asarray(res.grids[0].batch_energy_mev) * 2 / QUANTUM_MEV
        assert np.allclose(quanta, np.rint(quanta), atol=1e-3)
        results[backend] = res
    a, b = (
        results["python"].grids[0].batch_energy_mev,
        results["warp-cpu"].grids[0].batch_energy_mev,
    )
    assert np.array_equal(a, b)  # float64 python and warp-cpu: bit-identical fixed-point grids


def test_seg_piece_walks_a_segment_through_planes() -> None:
    import warp as wp

    f = make_transport_funcs(wp.float64)
    v = f.vec3
    r = wp.float64
    origin, spacing = v(r(0.0), r(0.0), r(0.0)), v(r(1.0), r(1.0), r(1.0))
    # along +x from x = 0.5 in voxel 0: the plane x = 1 ends the piece after 0.5 mm
    piece, axis = f.seg_piece(v(r(0.5), r(0.5), r(0.5)), v(r(1.0), r(0.0), r(0.0)), 0, 0, 0,
                              origin, spacing, r(0.8))  # fmt: skip
    assert (float(piece), int(axis)) == (0.5, 0)
    # the segment ends before the plane: the whole remaining length, no axis
    piece, axis = f.seg_piece(v(r(0.5), r(0.5), r(0.5)), v(r(1.0), r(0.0), r(0.0)), 0, 0, 0,
                              origin, spacing, r(0.3))  # fmt: skip
    assert (float(piece), int(axis)) == (0.3, -1)
    # starting on a plane against the travel direction: zero-length piece, axis of that plane
    piece, axis = f.seg_piece(v(r(1.0), r(0.5), r(0.5)), v(r(-1.0), r(0.0), r(0.0)), 1, 0, 0,
                              origin, spacing, r(0.4))  # fmt: skip
    assert (float(piece), int(axis)) == (0.0, 0)


def test_capacity_guard_and_overflow_counter(make_config: MakeConfig) -> None:
    """A batch whose energy could exceed the int64 capacity is rejected before transport; a
    voxel at capacity at the end of a run invalidates the result (counter)."""
    n_big = 2**31 * 2  # 4.3e9 histories in 2 batches of 150 MeV each exceed 2**62 quanta
    with pytest.raises(UnsupportedCombinationError, match="capacity"):
        validate(make_config(energy=150.0, n=n_big, n_batches=2))
    ok = validate(make_config(energy=150.0, n=2**20, n_batches=2))
    assert ok.requested.run.n_histories == 2**20

    wrapped = PartialTransport(
        0, 2, [[1.0]] * 8, [0] * 9,
        [np.array([[MAX_QUANTA, 0]], dtype=np.int64).repeat(2, axis=0)], None,
    )  # fmt: skip
    from ionmc.transport.tally import rows_to_partial

    part = rows_to_partial(
        0, 2, np.ones((2, 8)), np.zeros((2, 9), dtype=np.int32),
        [np.array([[MAX_QUANTA, 0]], dtype=np.int64)], None,
    )  # fmt: skip
    raw = merge_partials([part], 2, 1)
    assert raw.counters["accumulator_overflow"] == 1
    assert wrapped.h1 == 2


@pytest.mark.parametrize("backend", ["python", "warp-cpu"])
def test_scoring_piece_bound_overflow_invalidates_the_run(
    make_config: MakeConfig, backend: str
) -> None:
    """A leg that needs more voxel pieces than the validated bound is never truncated silently:
    the remainder is deposited (energy is conserved), the counter is raised, the result invalid."""
    from dataclasses import replace

    cfg = make_config(
        energy=20.0,
        n=2,
        n_batches=2,
        direction=(0.3, 0.2, 1.0),
        geometry=BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 12.0), WATER),
        scoring=(ScoringGrid((-10.0, -10.0, 0.0), (0.5, 0.5, 0.5), (40, 40, 24)),),
        mcs=False,
        straggling=False,
        max_step=2.0,
        backend=backend,
        allow_invalid=True,
    )
    sim = Simulation(cfg)
    assert sim.effective.scoring_pieces == 3 * 4 + 4
    assert sim.run().valid  # the computed bound suffices
    sim.effective = replace(sim.effective, scoring_pieces=1)
    res = sim.run()
    assert not res.valid and res.counters.scoring_pieces_overflow > 0
    assert res.energy_balance.grid_relative_residual(0) <= 1e-12 or backend == "warp-cpu"


def test_steps_longer_than_the_scoring_bins_are_apportioned(make_config: MakeConfig) -> None:
    """Steps of 2.5 mm over 1 mm bins (formerly rejected) are spread over the voxels crossed. The
    IDD agrees with the 0.1 mm-step run to 1e-2 on the plateau (the deposit is uniform along a
    step: a shape error of the step length, not a scoring artefact) and the run is valid."""
    ref = _idd(make_config, 0.1, 1.0)
    a = _idd(make_config, 2.5, 1.0)
    z = np.arange(len(ref)) + 0.5
    keep = (ref > 0.01 * ref.max()) & (z >= 20.0) & (z <= 120.0)
    assert np.abs(a[keep] / ref[keep] - 1.0).max() <= 1e-2
