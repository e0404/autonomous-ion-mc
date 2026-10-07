"""The world clip plane ``VoxelGeometry.z_exit_mm`` (the far face of the world independent of the
voxel grid, used by T14 so that a shifted grid overhangs instead of lengthening the slab)."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest

from ionmc._wpfunc import python_twin
from ionmc.config import DiagnosticsOptions, SimulationConfig
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import VoxelGeometry
from ionmc.materials import WATER
from ionmc.scoring import ScoringGrid, voxel_mass_g
from ionmc.simulation import Simulation
from ionmc.transport.funcs import BIG_LENGTH_MM, make_transport_funcs
from ionmc.transport.parity import compare_traces

MakeConfig = Callable[..., SimulationConfig]


def _geo(z_exit: float | None, z0: float = -1.0) -> VoxelGeometry:
    shape = (12, 12, 14)
    return VoxelGeometry(
        (-12.0, -12.0, z0), (2.0, 2.0, 2.0), shape, (WATER,), np.zeros(shape, dtype=np.int32),
        z_exit_mm=z_exit,
    )  # fmt: skip


def test_clip_plane_is_validated() -> None:
    _geo(20.0)
    for bad in (-1.0, 100.0, float("nan")):
        with pytest.raises(UnsupportedCombinationError):
            _geo(bad)


@pytest.mark.parametrize("z0", [-1.0, 0.0])
def test_particles_leave_exactly_at_the_clip_plane_on_every_backend(
    make_config: MakeConfig, z0: float
) -> None:
    """A half-voxel-shifted grid (z0 = -1 mm of 2 mm voxels) and an aligned one: the world ends at
    z_exit either way; escaping particles are recorded exactly there, python and warp-cpu float64
    follow the same trajectories (T1) and the energy balance closes."""
    results = {}
    for backend in ("python", "warp-cpu"):
        cfg = make_config(
            energy=60.0, n=6, n_batches=2, lateral_sigma=0.4, geometry=_geo(20.0, z0),
            scoring=(ScoringGrid((-12.0, -12.0, 0.0), (2.0, 2.0, 2.0), (12, 12, 10)),),
            max_step=1.0, backend=backend,
            diagnostics=DiagnosticsOptions(escape_records=True, track_end_positions=True,
                                           trace_histories=6),
        )  # fmt: skip
        results[backend] = Simulation(cfg).run()
    py, wr = results["python"], results["warp-cpu"]
    esc = py.diagnostics["escape_position_mm"]
    assert len(esc) == 6 and np.all(esc[:, 2] == 20.0)  # all leave through the clip plane
    assert np.array_equal(esc, wr.diagnostics["escape_position_mm"])
    assert compare_traces(py.diagnostics, wr.diagnostics)["pass"]
    assert py.energy_balance.relative_residual <= 1e-12
    assert py.energy_balance.escaped_mev > 0.0 and py.valid and wr.valid


def test_default_geometry_has_no_clip_and_masses_respect_it() -> None:
    assert _geo(None).world_upper_mm == _geo(None).upper_mm
    assert _geo(20.0).world_upper_mm[2] == 20.0
    grid = ScoringGrid((-12.0, -12.0, 0.0), (12.0, 12.0, 4.0), (2, 2, 7))  # 0..28 mm
    full, clipped = voxel_mass_g(grid, _geo(None)), voxel_mass_g(grid, _geo(20.0))
    assert clipped[:, :, :5].sum() == pytest.approx(full[:, :, :5].sum())  # below 20 mm
    assert clipped[:, :, 5:].sum() == 0.0 and full[:, :, 5:].sum() > 0.0


def test_dda_next_clip_picks_the_clip_plane_and_reduces_to_dda_next() -> None:
    f = python_twin(make_transport_funcs)  # pure-Python twin, no Warp call
    v, r = f.vec3, float
    org, sp = v(r(0.0), r(0.0), r(0.0)), v(r(1.0), r(1.0), r(1.0))
    p, u = v(r(0.5), r(0.5), r(0.5)), v(r(0.0), r(0.0), r(1.0))
    d, ax = f.dda_next_clip(p, u, 0, 0, 0, org, sp, r(BIG_LENGTH_MM))
    d0, ax0 = f.dda_next(p, u, 0, 0, 0, org, sp)
    assert (float(d), int(ax)) == (float(d0), int(ax0)) == (0.5, 2)
    d, ax = f.dda_next_clip(
        p, u, 0, 0, 0, org, sp, r(0.8)
    )  # clip plane nearer than the voxel plane
    assert (float(d), int(ax)) == (pytest.approx(0.3), 3)
    d, ax = f.dda_next_clip(p, u, 0, 0, 0, org, sp, r(1.0))  # a tie goes to the clip plane
    assert (float(d), int(ax)) == (0.5, 3)
    d, ax = f.dda_next_clip(p, v(r(0.0), r(0.0), r(-1.0)), 0, 0, 0, org, sp, r(0.8))
    assert int(ax) == 2  # moving away from the clip plane: it is ignored
    leg, axis = f.leg2_limit_clip(p, u, 0, 0, 0, org, sp, r(2.0), r(0.8))
    assert (float(leg), int(axis)) == (pytest.approx(0.3), 3)
