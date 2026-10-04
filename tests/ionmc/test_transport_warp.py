"""Warp CPU backend: trajectory parity with the Python reference (T1) and the frozen
deterministic and conservation checks (T2-T4, T11) on the production float32 kernel.

The Python reference and the float64 Warp kernel draw the same Philox blocks in the same
order and execute the same shared Warp functions, so their per-step traces are compared
step by step: discrete columns exactly, continuous columns within rtol = atol = 1e-10.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np
import pytest

from ionmc.config import DiagnosticsOptions, SimulationConfig
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.materials import PMMA, WATER
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.transport.tally import COUNTER_NAMES, TRACE_COLUMNS, TRACE_N_DISCRETE

MakeConfig = Callable[..., SimulationConfig]

DISCRETE = TRACE_COLUMNS[:TRACE_N_DISCRETE]
CONTINUOUS = TRACE_COLUMNS[TRACE_N_DISCRETE:]


def _water_voxels(voxel_mm: float = 5.0, depth_mm: float = 160.0) -> VoxelGeometry:
    shape = (12, 12, int(depth_mm / voxel_mm))
    shape = (int(60.0 / voxel_mm),) * 2 + (shape[2],)
    return VoxelGeometry(
        origin_mm=(-30.0, -30.0, 0.0),
        spacing_mm=(voxel_mm,) * 3,
        shape=shape,
        materials=(WATER,),
        material_index=np.zeros(shape, dtype=np.int32),
    )


def _scoring(depth_mm: float = 160.0) -> tuple[ScoringGrid, ...]:
    return (ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, int(depth_mm / 2.0))),)


def test_t1_trajectory_parity_python_vs_warp_cpu_f64(make_config: MakeConfig) -> None:
    """K = 16 histories of 100 MeV protons in a water voxel box, all physics on: the
    python reference and the warp-cpu float64 kernel agree on every discrete trace column
    exactly and on every continuous column within rtol = atol = 1e-10; any branch flip
    fails."""
    k = 16
    diag = DiagnosticsOptions(track_end_positions=True, escape_records=True, trace_histories=k)

    def run(backend: str):  # type: ignore[no-untyped-def]
        cfg = make_config(
            energy=100.0,
            n=k,
            n_batches=2,
            seed=20261004,
            geometry=_water_voxels(),
            scoring=_scoring(),
            diagnostics=diag,
            backend=backend,
            lateral_sigma=1.0,
            energy_sigma=0.5,
        )
        return Simulation(cfg).run()

    ref, wrp = run("python"), run("warp-cpu")
    tr_r, tr_w = ref.diagnostics["trace"], wrp.diagnostics["trace"]
    assert len(tr_r["step"]) == len(tr_w["step"]) > 1000
    for name in DISCRETE:
        assert np.array_equal(tr_r[name], tr_w[name]), name
    for name in CONTINUOUS:
        np.testing.assert_allclose(tr_w[name], tr_r[name], rtol=1e-10, atol=1e-10, err_msg=name)
    for key in ("trace_end_history", "trace_end_code", "end_code"):
        assert np.array_equal(ref.diagnostics[key], wrp.diagnostics[key]), key
    np.testing.assert_allclose(
        wrp.diagnostics["trace_end_energy_mev"],
        ref.diagnostics["trace_end_energy_mev"],
        rtol=1e-10,
        atol=1e-10,
    )
    # the trace has power: voxel crossings, energy-loss and maximum-step limited steps
    assert {0, 1, 3}.issubset(set(tr_r["reason"].astype(int)))
    assert tr_r["attempts"].max() >= 1
    # counters and tallies of the two backends are identical per history sums
    assert wrp.counters == ref.counters and not wrp.counters.any_nonzero
    assert wrp.energy_balance.relative_residual <= 1e-12
    for field in ("initial_mev", "step_deposit_mev", "cutoff_mev", "escaped_mev"):
        assert getattr(wrp.energy_balance, field) == pytest.approx(
            getattr(ref.energy_balance, field), rel=1e-12
        )
    np.testing.assert_allclose(
        wrp.grids[0].batch_energy_mev, ref.grids[0].batch_energy_mev, rtol=1e-9, atol=1e-12
    )
    assert wrp.backend == "warp-cpu" and wrp.timings["compile"] >= 0.0


def _water_box(depth: float = 200.0) -> BoxPhantom:
    return BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, depth), WATER)


def test_t2_deterministic_csda_end_depth_f32(make_config: MakeConfig) -> None:
    """100 MeV, straggling and MCS off, float32 on warp-cpu: the end depth satisfies the frozen
    bounds ``R(E0) - R(e_cut) - 1e-4 R <= z_end <= R(E0) + 1e-4 R``."""
    cfg = make_config(
        energy=100.0,
        n=4,
        geometry=_water_box(),
        scoring=(ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 100)),),
        mcs=False,
        straggling=False,
        diagnostics=DiagnosticsOptions(track_end_positions=True),
        backend="warp-cpu",
        precision="float32",
    )
    sim = Simulation(cfg)
    res = sim.run()
    tab = sim.effective.tables
    r0 = tab.range_g_cm2(0, 100.0) * 10.0 / WATER.density_g_cm3
    r_cut = tab.range_g_cm2(0, cfg.physics.e_cut_mev) * 10.0 / WATER.density_g_cm3
    end = res.diagnostics["end_position_mm"]
    assert np.all(res.diagnostics["end_code"] == 0)
    z = end[:, 2]
    assert np.allclose(end[:, :2], 0.0, atol=1e-6)
    assert np.all(z >= r0 - r_cut - 1e-4 * r0) and np.all(z <= r0 + 1e-4 * r0)
    assert res.valid and res.energy_balance.relative_residual < 1e-5
    assert res.effective_config.production is True


def test_t3_finite_slab_escape_energy_f32(make_config: MakeConfig) -> None:
    """100 MeV through 30 mm of water, deterministic, float32: escaped energy per primary is
    Rinv(R(100) - 3.0 g/cm2) within 1e-4 relative."""
    cfg = make_config(
        energy=100.0,
        n=4,
        position=(0.0, 0.0, -5.0),
        geometry=BoxPhantom((-20.0, -20.0, 0.0), (40.0, 40.0, 30.0), WATER),
        scoring=(ScoringGrid((-20.0, -20.0, 0.0), (2.0, 2.0, 2.0), (20, 20, 15)),),
        mcs=False,
        straggling=False,
        diagnostics=DiagnosticsOptions(escape_records=True),
        backend="warp-cpu",
        precision="float32",
    )
    sim = Simulation(cfg)
    res = sim.run()
    tab = sim.effective.tables
    expected = tab.energy_from_range(0, tab.range_g_cm2(0, 100.0) - 3.0)
    assert res.energy_balance.escaped_mev / 4 == pytest.approx(expected, rel=1e-4)
    esc = res.diagnostics
    assert len(esc["escape_history"]) == 4
    assert np.allclose(esc["escape_position_mm"][:, 2], 30.0, atol=1e-4)
    assert np.allclose(esc["escape_energy_mev"], expected, rtol=1e-4, atol=0)
    assert res.energy_balance.relative_residual < 1e-5


def _two_material_geometry() -> VoxelGeometry:
    idx = np.zeros((6, 6, 6), dtype=np.int32)
    idx[:, :, 3:] = 1
    dens = np.where(idx == 1, 1.19, 1.0).astype(np.float64)
    dens[0, 0, 0] = 0.9  # a density override inside a material
    return VoxelGeometry(
        origin_mm=(0.0, 0.0, 0.0),
        spacing_mm=(2.0, 2.0, 2.0),
        shape=(6, 6, 6),
        materials=(WATER, PMMA),
        material_index=idx,
        density_g_cm3=dens,
    )


def test_t4_energy_balance_every_destination_f32_and_f64(make_config: MakeConfig) -> None:
    """Initial energy equals the independently tallied destinations: 1e-5 (float32) and 1e-12
    (float64) on warp-cpu; oblique beam with energy and lateral spread through a two-material
    density-modulated cube, two scoring grids (one smaller than the geometry)."""
    for precision, tol in (("float32", 1e-5), ("float64", 1e-12)):
        small = ScoringGrid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (6, 12, 12), name="half")
        wide = ScoringGrid((-3.0, -3.0, -3.0), (1.5, 1.5, 1.5), (14, 14, 14), name="wide")
        cfg = make_config(
            energy=40.0,
            n=24,
            n_batches=4,
            position=(1.0, 2.0, -3.0),
            direction=(1.0, 0.4, 2.0),
            energy_sigma=0.4,
            lateral_sigma=0.8,
            geometry=_two_material_geometry(),
            scoring=(small, wide),
            max_step=1.0,
            backend="warp-cpu",
            precision=precision,
        )
        res = Simulation(cfg).run()
        eb = res.energy_balance
        assert res.valid and not res.counters.any_nonzero
        assert eb.relative_residual <= tol, precision
        for g in range(2):
            assert eb.grid_relative_residual(g) <= tol, precision
        assert eb.outside_mev[0] > 0.0 and eb.escaped_mev > 0.0 and eb.cutoff_mev > 0.0
        assert eb.truncated_mev == 0.0 and eb.unaccounted_mev == 0.0
        for g, grid in enumerate(res.grids):
            assert grid.energy_mev.sum() * cfg.run.n_histories == pytest.approx(
                eb.in_grid_mev[g], rel=1e-6 if precision == "float32" else 1e-12
            )


S3 = 1.0 / math.sqrt(3.0)
ADVERSARIAL = [
    ("corner_diagonal", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), False),
    ("corner_diagonal_mcs", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), True),
    ("interior_vertex_up", (10.0, 10.0, 10.0), (1.0, 1.0, 1.0), True),
    ("interior_vertex_down", (10.0, 10.0, 10.0), (-1.0, -1.0, -1.0), False),
    ("along_edge_zero_components", (10.0, 10.0, 4.0), (0.0, 0.0, 1.0), False),
    ("along_edge_mcs", (10.0, 10.0, 4.0), (0.0, 0.0, 1.0), True),
    ("on_plane_zero_component", (0.0, 6.0, 6.0), (1.0, 0.0, 0.0), False),
    ("grazing_entry_on_face", (-5.0, 0.0, 10.0), (1.0, 0.0, 0.0), True),
    ("grazing_entry_far_face", (-5.0, 20.0, 9.0), (1.0, 0.0, 0.0), False),
    ("entry_through_corner", (-3.0, -3.0, -3.0), (1.0, 1.0, 1.0), True),
    ("exit_through_edge", (16.0, 16.0, 10.0), (1.0, 1.0, 0.0), False),
    ("oblique_on_plane", (4.0, 10.0, 10.0), (0.6, 0.8, 0.0), True),
    ("misses_world", (-5.0, -5.0, 30.0), (1.0, 0.0, 0.0), True),
]


@pytest.mark.parametrize("precision", ["float32", "float64"])
@pytest.mark.parametrize(
    ("name", "position", "direction", "mcs"), ADVERSARIAL, ids=[a[0] for a in ADVERSARIAL]
)
def test_t11_dda_adversarial_cases_warp_cpu(
    make_config: MakeConfig,
    name: str,
    position: tuple[float, float, float],
    direction: tuple[float, float, float],
    mcs: bool,
    precision: str,
) -> None:
    """Sources on planes, edges and corners, zero direction components, grazing entry: no stall
    or truncation, energy balance closes (1e-5 float32, 1e-12 float64)."""
    geo = VoxelGeometry(
        origin_mm=(0.0, 0.0, 0.0),
        spacing_mm=(2.0, 2.0, 2.0),
        shape=(10, 10, 10),
        materials=(WATER,),
        material_index=np.zeros((10, 10, 10), dtype=np.int32),
    )
    cfg = make_config(
        energy=30.0,
        n=8,
        n_batches=2,
        position=position,
        direction=direction,
        geometry=geo,
        scoring=(ScoringGrid((0.0, 0.0, 0.0), (2.0, 2.0, 2.0), (10, 10, 10)),),
        mcs=mcs,
        straggling=mcs,
        max_step=2.0,
        backend="warp-cpu",
        precision=precision,
    )
    res = Simulation(cfg).run()
    assert res.counters.as_dict() == dict.fromkeys(COUNTER_NAMES, 0), name
    tol = 1e-5 if precision == "float32" else 1e-12
    assert res.valid and res.energy_balance.relative_residual <= tol
    assert res.energy_balance.grid_relative_residual(0) <= tol
    if name == "misses_world":
        assert res.energy_balance.escaped_mev == res.energy_balance.initial_mev
        assert res.energy_balance.step_deposit_mev == 0.0


def test_trace_is_float64_only(make_config: MakeConfig) -> None:
    cfg = make_config(
        backend="warp-cpu",
        precision="float32",
        diagnostics=DiagnosticsOptions(trace_histories=1),
    )
    with pytest.raises(UnsupportedCombinationError, match="float64"):
        Simulation(cfg)


def test_float32_kernel_keeps_energy_and_range_bookkeeping_in_double(
    make_config: MakeConfig,
) -> None:
    """The float32 kernel holds positions and directions in float32 but energy, range, mean loss,
    straggling and deposits in float64: a deterministic (straggling and MCS off) 100 MeV track
    ends at the float64 end depth to well below the old float32 bookkeeping bias (about 2e-5 of
    the range), and every depth bin of the depth-dose agrees to 1e-5."""
    out = {}
    for precision in ("float32", "float64"):
        cfg = make_config(
            energy=100.0, n=4, geometry=_water_box(),
            scoring=(ScoringGrid((-30.0, -30.0, 0.0), (60.0, 60.0, 1.0), (1, 1, 100)),),
            mcs=False, straggling=False, max_step=1.0, backend="warp-cpu", precision=precision,
            diagnostics=DiagnosticsOptions(track_end_positions=True),
        )  # fmt: skip
        res = Simulation(cfg).run()
        out[precision] = (
            float(res.diagnostics["end_position_mm"][:, 2].mean()),
            np.asarray(res.grids[0].energy_mev).reshape(-1),
        )
    assert abs(out["float32"][0] - out["float64"][0]) < 5e-4  # mm (the old bias was ~1.5e-3 mm)
    a, b = out["float32"][1], out["float64"][1]
    keep = b > 0.01 * b.max()
    assert np.abs(a[keep] / b[keep] - 1.0).max() < 1e-5


def _path_loss(make_config: MakeConfig, s_max: float, fraction: float) -> float:
    """Deterministic energy lost by a 150 MeV proton in 40 mm of water (float64 warp-cpu)."""
    from dataclasses import replace

    cfg = make_config(
        energy=150.0, n=2, n_batches=2, position=(0.0, 0.0, -1.0),
        geometry=BoxPhantom((-20.0, -20.0, 0.0), (40.0, 40.0, 40.0), WATER),
        scoring=(ScoringGrid((-20.0, -20.0, 0.0), (40.0, 40.0, 40.0), (1, 1, 1)),),
        mcs=False, straggling=False, max_step=s_max, backend="warp-cpu", precision="float64",
        diagnostics=DiagnosticsOptions(escape_records=True),
    )  # fmt: skip
    cfg = replace(cfg, physics=replace(cfg.physics, short_step_fraction=fraction))
    res = Simulation(cfg).run()
    return 150.0 - float(res.diagnostics["escape_energy_mev"][0])


def test_short_step_branch_has_no_first_order_step_dependence(make_config: MakeConfig) -> None:
    """The linear short-step branch uses the stopping power at the midpoint energy
    ``E - S(E) t / 2`` (second order in the step) and is the default for ``t < 1e-2 R``. With the
    default fraction the loss over the same path agrees between 1 mm, 0.1 mm and 0.01 mm steps
    to the second-order level (1 vs 0.1 mm: 5e-6 relative; 0.1 vs 0.01 mm: 1e-7); a first-order
    branch would differ by about 2e-4 between 0.1 and 0.01 mm, and the telescoping form alone
    (fraction 1e-5) by about 3e-4 from the table round trip that enters once per step."""
    from ionmc.config import PhysicsOptions

    assert PhysicsOptions.__dataclass_fields__["short_step_fraction"].default == 1e-2
    loss = {s: _path_loss(make_config, s, 1e-2) for s in (1.0, 0.1, 0.01)}
    assert abs(loss[0.1] / loss[1.0] - 1.0) < 5e-6
    assert abs(loss[0.01] / loss[0.1] - 1.0) < 1e-7
    telescoping = {s: _path_loss(make_config, s, 1e-5) for s in (0.1, 0.01)}
    assert abs(telescoping[0.01] / telescoping[0.1] - 1.0) > 1e-4  # the bias the branch removes
