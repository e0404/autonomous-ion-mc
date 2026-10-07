"""The truncated-hinge diagnostic (T14 negative control): off by default and then the engine
of the stored baseline trace (see the provenance in the test); on, python and warp-cpu float64
still follow the same trajectory (T1) and the trajectory differs from the default."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
)
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import VoxelGeometry
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport.parity import compare_traces
from ionmc.transport.tally import TRACE_COLUMNS

BASELINE = Path(__file__).parent / "data" / "trace_baseline_20mev.npz"


def _config(
    backend: str = "python", *, diagnostic: bool = False, legacy: bool = False
) -> SimulationConfig:
    shape = (16, 16, 16)
    geo = VoxelGeometry(
        (-12.0, -12.0, 0.0), (1.5, 1.5, 1.5), shape, (WATER,), np.zeros(shape, dtype=np.int32)
    )
    return SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.1, 0.05, 1.0), 20.0, 0.0, 0.3),
        geo,
        (ScoringGrid((-12.0, -12.0, 0.0), (2.0, 2.0, 2.0), (12, 12, 12)),),
        PhysicsOptions(
            nuclear=False,
            stopping=BetheStoppingSource(),
            max_step_mm=1.0,
            truncated_hinge_diagnostic=diagnostic,
            **(  # the physics defaults of the commit that generated the baseline trace fixture
                {"straggling_model": "bohr_gauss_clamped_gamma_v1", "short_step_fraction": 1e-3}
                if legacy
                else {}
            ),
        ),
        RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision="float64",
            seed=7,
            n_histories=4,
            n_batches=2,
        ),
        DiagnosticsOptions(trace_histories=4, track_end_positions=True),
    )


def test_default_trace_is_unchanged_from_the_stored_baseline() -> None:
    """Provenance of ``data/trace_baseline_20mev.npz``: generated at commit 6d58480 (before the
    truncated-hinge diagnostic existed) by the python reference with the physics defaults of that
    commit, i.e. straggling model ``bohr_gauss_clamped_gamma_v1`` and ``short_step_fraction`` 1e-3
    (both pinned here; the defaults have since changed). It guards that, with the diagnostic off,
    the engine reproduces that trajectory: discrete columns exactly, continuous columns to
    rtol 1e-9 and atol 2e-10. One documented, intended deviation exists: the midpoint rule of the
    linear short-step energy-loss branch (V3-003B, after 6d58480) changes the loss of the one step
    of this run that takes that branch by about 8e-11 MeV (the trace differs by at most 8e-11 in
    energy and 2e-11 mm in position). On the generating host python and warp-cpu otherwise agree
    bit for bit."""
    base = np.load(BASELINE)
    for backend in ("python", "warp-cpu"):
        res = Simulation(_config(backend, legacy=True)).run()
        tr = res.diagnostics["trace"]
        for name in TRACE_COLUMNS:
            if name in ("history", "step", "ix", "iy", "iz", "reason", "blocks", "attempts"):
                assert np.array_equal(tr[name], base[name]), (backend, name)
            else:
                # atol 2e-10 covers the documented 8e-11 MeV midpoint-rule change of one step
                np.testing.assert_allclose(tr[name], base[name], rtol=1e-9, atol=2e-10)


def test_diagnostic_on_keeps_python_warp_parity_and_changes_trajectories() -> None:
    py = Simulation(_config("python", diagnostic=True)).run()
    wp = Simulation(_config("warp-cpu", diagnostic=True)).run()
    v = compare_traces(py.diagnostics, wp.diagnostics)
    assert v["pass"], v
    base = np.load(BASELINE)
    ux = py.diagnostics["trace"]["ux"]
    assert len(ux) != len(base["ux"]) or not np.allclose(ux, base["ux"], rtol=0, atol=1e-12)
    assert py.effective_config.summary()["physics"]["truncated_hinge_diagnostic"] is True
    default = Simulation(_config("python")).run()
    assert default.effective_config.summary()["physics"]["truncated_hinge_diagnostic"] is False


def test_diagnostic_flag_is_validated() -> None:
    cfg = _config()
    with pytest.raises(UnsupportedCombinationError, match="bool"):
        replace(cfg, physics=replace(cfg.physics, truncated_hinge_diagnostic=1))  # type: ignore[arg-type]


@pytest.mark.parametrize("diagnostic", [False, True])
def test_steps_that_change_voxel_end_exactly_on_the_crossed_plane(diagnostic: bool) -> None:
    """Whatever the direction used for the second leg (default or the control's re-sampled one),
    a step that moves to another transport voxel ends exactly on the plane it crossed."""
    cfg = _config("python", diagnostic=diagnostic)
    res = Simulation(cfg).run()
    tr = res.diagnostics["trace"]
    origin, spacing = (-12.0, -12.0, 0.0), 1.5
    idx = np.stack([tr["ix"], tr["iy"], tr["iz"]], axis=1).astype(int)
    pos = np.stack([tr["x_mm"], tr["y_mm"], tr["z_mm"]], axis=1)
    hist = tr["history"].astype(int)
    n_checked = 0
    for k in range(1, len(hist)):
        if hist[k] != hist[k - 1]:
            continue
        d = idx[k] - idx[k - 1]
        if not d.any():
            continue
        assert np.count_nonzero(d) == 1 and abs(d.sum()) == 1
        a = int(np.nonzero(d)[0][0])
        plane = origin[a] + (idx[k, a] if d[a] > 0 else idx[k, a] + 1) * spacing
        assert pos[k, a] == pytest.approx(plane, abs=1e-12)
        n_checked += 1
    assert n_checked > 10


def test_truncate_first_control_ends_cut_steps_on_the_plane_with_bounded_snap() -> None:
    """Truncate-first control: the planned step already ends at the first plane the straight line
    reaches; leg 2 is not cut again and the end point is snapped onto that plane. The snap
    displacement relative to the step, recorded per history (python and kernel identical), is
    recorded and small; every geometry-limited step changes the voxel (it ends on its
    plane); the travelled hinge path equals the step exactly; the default records nothing."""
    py = Simulation(_config("python", diagnostic=True)).run()
    wp = Simulation(_config("warp-cpu", diagnostic=True)).run()
    assert np.array_equal(py.diagnostics["control_residual"], wp.diagnostics["control_residual"])
    tr = py.diagnostics["trace"]
    idx = np.stack([tr["ix"], tr["iy"], tr["iz"]], axis=1)
    hist = tr["history"]
    n_cut = 0
    for k in range(1, len(hist)):
        if hist[k] == hist[k - 1] and int(tr["reason"][k]) == 0:
            assert (idx[k] != idx[k - 1]).any()  # ended on its plane: a new voxel
            n_cut += 1
    assert n_cut > 10
    assert np.array_equal(
        py.diagnostics["control_displacement"], wp.diagnostics["control_displacement"]
    )
    assert py.diagnostics["control_residual"].max() > 0.0
    assert np.abs(py.diagnostics["control_displacement"]).max() > 0.0
    off = Simulation(_config("warp-cpu")).run().diagnostics["control_residual"]
    assert not off.any()
    assert not Simulation(_config("warp-cpu")).run().diagnostics["control_displacement"].any()


def test_control_direction_after_the_hinge_is_the_sampled_direction_exactly() -> None:
    """The snap moves only the end position, never the direction: in control mode the direction
    after every step equals, bit for bit, the direction returned by the scattering rotation."""
    from ionmc.config import validate
    from ionmc.transport.reference import _Reference

    eff = validate(_config("python", diagnostic=True))
    ref = _Reference(eff)
    sampled: list[tuple[float, float, float]] = []
    original = ref.EM.rotate_dir

    def spy(u, theta, phi):  # type: ignore[no-untyped-def]
        out = original(u, theta, phi)
        sampled.append((float(out[0]), float(out[1]), float(out[2])))
        return out

    ref.EM.rotate_dir = spy
    try:
        part = ref.run_range(0, 4)
    finally:
        ref.EM.rotate_dir = original  # the namespace is shared by every reference run
    tr = part.diagnostics.trace_float  # columns x y z ux uy uz E deposit step
    assert len(sampled) == len(tr) > 100
    for k, d in enumerate(sampled):
        assert (tr[k, 3], tr[k, 4], tr[k, 5]) == d
