"""T2-T4, T11 and result/diagnostic checks of the python (reference) backend."""

import math
from collections.abc import Callable

import numpy as np
import pytest

from ionmc.config import DiagnosticsOptions, SimulationConfig
from ionmc.errors import TransportLimitError
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.materials import PMMA, WATER
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import MEV_PER_G_TO_GY, ScoringGrid
from ionmc.simulation import Simulation, capabilities

MakeConfig = Callable[..., SimulationConfig]


def _water_box(depth: float = 200.0) -> BoxPhantom:
    return BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, depth), WATER)


def test_t2_deterministic_csda_end_depth(make_config: MakeConfig) -> None:
    """100 MeV, straggling and MCS off: the track ends where the CSDA range says.

    The frozen text asks for 200 histories; the deterministic setup makes every history
    identical, so 2 are run (measured cost of the python backend: about 0.25 s per history at
    100 MeV with 2 mm steps, about 4 ms per step).
    """
    cfg = make_config(
        energy=100.0,
        n=2,
        geometry=_water_box(),
        scoring=(ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 100)),),
        mcs=False,
        straggling=False,
        diagnostics=DiagnosticsOptions(track_end_positions=True),
    )
    sim = Simulation(cfg)
    res = sim.run()
    tab = sim.effective.tables
    r0 = tab.range_g_cm2(0, 100.0) * 10.0 / WATER.density_g_cm3
    r_cut = tab.range_g_cm2(0, cfg.physics.e_cut_mev) * 10.0 / WATER.density_g_cm3
    end = res.diagnostics["end_position_mm"]
    assert np.all(res.diagnostics["end_code"] == 0)  # ended by the cutoff
    z = end[:, 2]
    assert np.allclose(end[:, :2], 0.0, atol=1e-12)  # no scattering: stays on the axis
    assert np.all(z >= r0 - r_cut - 1e-4 * r0) and np.all(z <= r0 + 1e-4 * r0)
    assert z[0] == pytest.approx(z[1], rel=1e-12)  # the hinge fraction only changes rounding
    assert res.valid and not res.counters.any_nonzero
    assert res.energy_balance.relative_residual < 1e-12


def test_t3_finite_slab_escape_energy(make_config: MakeConfig) -> None:
    """100 MeV through 30 mm (3.0 g/cm2) of water, deterministic: escaped energy is
    Rinv(R(100) - 3.0 g/cm2). The beam starts in vacuum 5 mm before the slab."""
    cfg = make_config(
        energy=100.0,
        n=2,
        position=(0.0, 0.0, -5.0),
        geometry=BoxPhantom((-20.0, -20.0, 0.0), (40.0, 40.0, 30.0), WATER),
        scoring=(ScoringGrid((-20.0, -20.0, 0.0), (2.0, 2.0, 2.0), (20, 20, 15)),),
        mcs=False,
        straggling=False,
        diagnostics=DiagnosticsOptions(escape_records=True),
    )
    sim = Simulation(cfg)
    res = sim.run()
    tab = sim.effective.tables
    expected = tab.energy_from_range(0, tab.range_g_cm2(0, 100.0) - 3.0)
    per_primary = res.energy_balance.escaped_mev / 2
    assert per_primary == pytest.approx(expected, rel=1e-4)
    assert res.energy_balance.cutoff_mev == 0.0
    esc = res.diagnostics
    assert np.allclose(esc["escape_position_mm"][:, 2], 30.0, atol=1e-9)
    assert np.allclose(esc["escape_direction"], [[0, 0, 1]] * 2)
    assert np.allclose(esc["escape_energy_mev"], expected, rtol=1e-4)
    # deposited energy per primary + escaped energy = initial
    assert res.energy_balance.relative_residual < 1e-12


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


def test_t4_energy_balance_with_every_destination(make_config: MakeConfig) -> None:
    """Initial energy equals the independently tallied destinations (1e-12, float64): grid
    deposits, deposits outside a smaller grid, cutoff, escape; an oblique beam with energy and
    lateral spread through a two-material, density-modulated 12 mm cube."""
    geo = _two_material_geometry()
    small = ScoringGrid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (6, 12, 12), name="half")
    wide = ScoringGrid((-3.0, -3.0, -3.0), (1.5, 1.5, 1.5), (14, 14, 14), name="wide")
    cfg = make_config(
        energy=40.0,
        n=6,
        n_batches=2,
        position=(1.0, 2.0, -3.0),
        direction=(1.0, 0.4, 2.0),
        energy_sigma=0.4,
        lateral_sigma=0.8,
        geometry=geo,
        scoring=(small, wide),
        max_step=1.0,
        diagnostics=DiagnosticsOptions(track_end_positions=True),
    )
    res = Simulation(cfg).run()
    eb = res.energy_balance
    assert res.valid and not res.counters.any_nonzero
    assert eb.relative_residual <= 1e-12
    for g in range(2):
        assert eb.grid_relative_residual(g) <= 1e-12
    assert eb.outside_mev[0] > 0.0  # the half-size grid misses part of the deposit
    assert eb.escaped_mev > 0.0 and eb.cutoff_mev > 0.0 and eb.truncated_mev == 0.0
    assert eb.unaccounted_mev == 0.0
    # per-primary grid energies are consistent with the raw tallies
    for g, grid in enumerate(res.grids):
        assert grid.energy_mev.sum() * cfg.run.n_histories == pytest.approx(
            eb.in_grid_mev[g], rel=1e-12
        )
    codes = set(res.diagnostics["end_code"].tolist())
    assert codes <= {0, 1}
    assert eb.initial_mev == pytest.approx(
        sum(40.0 for _ in range(6)), abs=6 * 5 * 0.4
    )  # sampled energies scatter around the mean


def _grid_20() -> tuple[VoxelGeometry, ScoringGrid]:
    geo = VoxelGeometry(
        origin_mm=(0.0, 0.0, 0.0),
        spacing_mm=(2.0, 2.0, 2.0),
        shape=(10, 10, 10),
        materials=(WATER,),
        material_index=np.zeros((10, 10, 10), dtype=np.int32),
    )
    return geo, ScoringGrid((0.0, 0.0, 0.0), (2.0, 2.0, 2.0), (10, 10, 10))


S3 = 1.0 / math.sqrt(3.0)
ADVERSARIAL = [
    # name, position, direction, mcs
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


@pytest.mark.parametrize(
    ("name", "position", "direction", "mcs"), ADVERSARIAL, ids=[a[0] for a in ADVERSARIAL]
)
def test_t11_dda_adversarial_cases(
    make_config: MakeConfig,
    name: str,
    position: tuple[float, float, float],
    direction: tuple[float, float, float],
    mcs: bool,
) -> None:
    """Sources on planes, edges and corners, zero direction components, grazing entry: no
    stall or truncation, energy balance closes."""
    geo, grid = _grid_20()
    cfg = make_config(
        energy=30.0,
        n=2,
        position=position,
        direction=direction,
        geometry=geo,
        scoring=(grid,),
        mcs=mcs,
        straggling=mcs,  # straight deterministic tracks where scattering is off
        max_step=2.0,
        diagnostics=DiagnosticsOptions(track_end_positions=True, trace_histories=2),
    )
    res = Simulation(cfg).run()
    assert res.counters.as_dict() == dict.fromkeys(res.counters.as_dict(), 0), name
    assert res.valid
    assert res.energy_balance.relative_residual <= 1e-12
    assert res.energy_balance.grid_relative_residual(0) <= 1e-12
    if name == "misses_world":
        assert res.energy_balance.escaped_mev == res.energy_balance.initial_mev
        assert res.energy_balance.step_deposit_mev == 0.0
        assert set(res.diagnostics["end_code"].tolist()) == {4}
    elif not mcs and name in ("corner_diagonal", "along_edge_zero_components"):
        # straight-line path: the track end lies on the ray at the CSDA distance, however
        # many corner crossings the voxel grid forces
        d = np.array(direction) / np.linalg.norm(direction)
        end = res.diagnostics["end_position_mm"]
        along = (end - np.array(position)) @ d
        perp = (end - np.array(position)) - np.outer(along, d)
        assert np.allclose(perp, 0.0, atol=1e-9)
        tab = res.effective_config.tables
        r_mm = tab.range_g_cm2(0, 30.0) * 10.0
        r_cut = tab.range_g_cm2(0, 2.0) * 10.0
        assert np.all(along <= r_mm + 1e-6) and np.all(along >= r_mm - r_cut - 1e-6)
    # trace: indices stay in range, no step leaves the world unrecorded
    tr = res.diagnostics["trace"]
    if len(tr["step"]):
        assert np.all(tr["ix"] >= -1) and np.all(tr["ix"] <= 10)


def test_invalid_result_fails_closed_unless_allowed(make_config: MakeConfig) -> None:
    cfg = make_config(energy=30.0, n=2, max_steps=3)
    with pytest.raises(TransportLimitError) as err:
        Simulation(cfg).run()
    result = err.value.result
    assert result is not None and not result.valid
    assert result.counters.step_truncation == 2
    assert result.energy_balance.truncated_mev > 0.0
    assert result.energy_balance.relative_residual <= 1e-12  # truncated energy is tallied
    assert all(not g.valid for g in result.grids)
    # the truncated energy was never scored
    assert result.energy_balance.step_deposit_mev + result.energy_balance.truncated_mev == (
        pytest.approx(result.energy_balance.initial_mev, rel=1e-12)
    )
    lenient = Simulation(make_config(energy=30.0, n=2, max_steps=3, allow_invalid=True)).run()
    assert not lenient.valid and lenient.counters.step_truncation == 2


def test_source_energy_out_of_range_is_counted(make_config: MakeConfig) -> None:
    """A sampled energy above the table maximum is never transported: counter, invalid."""
    cfg = make_config(energy=20.0, energy_sigma=1.0, n=4, allow_invalid=True)
    cfg = SimulationConfigWithTable.apply(cfg, e_max=20.0)
    res = Simulation(cfg).run()
    assert res.counters.source_energy_out_of_range > 0 and not res.valid
    assert res.energy_balance.relative_residual <= 1e-12


class SimulationConfigWithTable:
    """Helper: a configuration whose stopping tables end at ``e_max`` MeV."""

    @staticmethod
    def apply(cfg: SimulationConfig, e_max: float) -> SimulationConfig:
        from dataclasses import replace

        physics = replace(cfg.physics, stopping=BetheStoppingSource(e_max_per_u=e_max))
        return replace(cfg, physics=physics)


def test_result_statistics_and_dose(make_config: MakeConfig) -> None:
    geo = BoxPhantom((0.0, 0.0, 0.0), (8.0, 8.0, 12.0), WATER)
    grid = ScoringGrid((-2.0, 0.0, 0.0), (2.0, 2.0, 2.0), (6, 4, 6), name="dose")
    cfg = make_config(
        energy=20.0,
        n=6,
        n_batches=3,
        position=(4.0, 4.0, 0.0),
        geometry=geo,
        scoring=(grid,),
        seed=11,
    )
    res = Simulation(cfg).run()
    g = res.grid("dose")
    b = res.n_batches
    # mean and variance of the mean follow the batch method
    per_batch = g.batch_energy_mev
    assert per_batch.shape == (3, 6, 4, 6)
    assert np.allclose(per_batch.mean(axis=0), g.energy_mev, rtol=1e-12, atol=0)
    var = ((per_batch - g.energy_mev) ** 2).sum(axis=0) / (b * (b - 1))
    assert np.allclose(np.sqrt(var), g.energy_std_mev, rtol=1e-12, atol=0)
    # voxel mass: partial coverage of the grid (x from -2 to 10 vs box 0..8)
    assert g.mass_g[0, 0, 0] == 0.0 and g.mass_g[2, 1, 1] == pytest.approx(8e-3)
    assert np.isnan(g.dose_gy[0, 0, 0]) and not g.defined_mask[0, 0, 0]
    d = g.defined_mask
    assert d.any()
    assert np.allclose(
        g.dose_gy[d], g.energy_mev[d] * MEV_PER_G_TO_GY / g.mass_g[d], rtol=1e-12, atol=0
    )
    # undefined uncertainty: NaN where nothing was deposited, never 0
    empty = g.energy_mev == 0.0
    assert empty.any() and np.all(np.isnan(g.relative_uncertainty[empty]))
    assert np.all(g.relative_uncertainty[d] >= 0.0)
    assert np.array_equal(g.n_nonzero_batches > 0, g.energy_mev > 0.0)
    assert g.low_batch_count  # 3 batches
    assert np.all(g.sparse_mask <= d)
    assert res.requested_config is cfg and res.effective_config.requested is cfg
    assert res.effective_config.summary()["tables"]["sha256"] == res.effective_config.tables.sha256
    assert set(res.timings) == {"setup", "compile", "transport", "reduce", "total"}
    assert res.environment["python_version"] and res.rng["generator"] == "philox4x32-10"
    assert res.rng["u01"] == "((w >> 8) + 0.5) * 2**-24"
    with pytest.raises(KeyError):
        res.grid("missing")


def test_diagnostics_trace_and_end_records(make_config: MakeConfig) -> None:
    diag = DiagnosticsOptions(track_end_positions=True, escape_records=True, trace_histories=2)
    res = Simulation(make_config(energy=15.0, n=4, diagnostics=diag)).run()
    d = res.diagnostics
    assert d["end_position_mm"].shape == (4, 3) and d["end_code"].shape == (4,)
    tr = d["trace"]
    assert set(d["trace_columns"]) == set(tr)
    assert set(tr["history"].astype(int)) == {0, 1}
    for h in (0, 1):
        sel = tr["history"] == h
        assert np.all(np.diff(tr["step"][sel]) == 1)  # consecutive steps
        assert np.all(np.diff(tr["blocks"][sel]) >= 2)  # block A plus at least one block B
        assert np.isclose(tr["deposit_mev"][sel].sum() + d["trace_end_energy_mev"][h], 15.0)
    assert d["escape_energy_mev"].shape == (0,)
    off = Simulation(make_config(energy=15.0, n=2)).run()
    assert off.diagnostics == {}


def test_capabilities_report() -> None:
    cap = capabilities()
    assert cap["species"] == ["proton"] and cap["physics"]["nuclear"] is False
    assert cap["backends"]["python"].startswith("available")
    assert cap["backends"]["warp-cpu"].startswith("available")
    assert cap["chunk_histories"]["default"] == 2**18
    assert "warp-cuda" in cap["backend_names"]
    import ionmc

    assert ionmc.capabilities() == cap


def test_getting_started_example_runs() -> None:
    """The minimal example of docs/getting-started/index.md runs and is valid."""
    import re
    from pathlib import Path

    text = (Path(__file__).resolve().parents[2] / "docs/getting-started/index.md").read_text()
    blocks = [b for b in re.findall(r"```python\n(.*?)```", text, re.S) if "Simulation(" in b]
    assert len(blocks) == 1
    namespace: dict[str, object] = {}
    exec(compile(blocks[0], "getting-started", "exec"), namespace)
    result = namespace["result"]
    assert result.valid and result.energy_balance.relative_residual <= 1e-12  # type: ignore[attr-defined]
    assert result.grid("dose").dose_gy.shape == (10, 10, 20)  # type: ignore[attr-defined]


def test_reproducibility_and_seed_dependence(make_config: MakeConfig) -> None:
    def run(seed: int) -> np.ndarray:
        cfg = make_config(energy=15.0, n=4, seed=seed)
        res = Simulation(cfg).run()
        return np.asarray(res.grids[0].batch_energy_mev)

    a, b, c = run(5), run(5), run(6)
    assert np.array_equal(a, b)  # identical streams, identical result
    assert not np.array_equal(a, c)


def test_history_belongs_to_batch_history_mod_b(make_config: MakeConfig) -> None:
    """With one voxel covering everything, batch b holds the energies of histories h % B == b."""
    n, b = 6, 3
    cfg = make_config(
        energy=15.0,
        n=n,
        n_batches=b,
        energy_sigma=0.5,
        scoring=(ScoringGrid((-30.0, -30.0, 0.0), (60.0, 60.0, 60.0), (1, 1, 1)),),
        max_step=2.0,
        diagnostics=DiagnosticsOptions(trace_histories=n),
    )
    res = Simulation(cfg).run()
    d = res.diagnostics
    energy = np.zeros(n)  # sampled energy of every history = all deposits + cutoff energy
    for h in range(n):
        energy[h] = d["trace"]["deposit_mev"][d["trace"]["history"] == h].sum()
        energy[h] += d["trace_end_energy_mev"][d["trace_end_history"] == h][0]
    assert len(set(np.round(energy, 9))) == n  # distinct per-history energies
    per_batch = res.grids[0].batch_energy_mev[:, 0, 0, 0] * (n // b)
    for k in range(b):
        assert per_batch[k] == pytest.approx(energy[k::b].sum(), rel=1e-12)
