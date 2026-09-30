"""Warp CPU/CUDA backends: statistical parity with the reference path, determinism, contract.

CPU tests run everywhere (Warp compiles CPU kernels without a GPU). Tests
marked ``cuda`` require a CUDA device and run through the controlled host
runner; they are skipped elsewhere. Parity uses per-voxel z-scores between
independent runs (different backends use different random streams), as
decision 0042 prescribes for stochastic output.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from ionmc.config import (
    PhysicsConfig,
    SimulationConfig,
    UnsupportedConfigurationError,
    capabilities,
)
from ionmc.geometry import homogeneous_box
from ionmc.scoring import ScoringGrid
from ionmc.simulation import build_tables, run
from ionmc.sources import PencilBeam


def depth_profile(result):
    z, prof = result.depth_dose(2)
    err = np.sqrt(np.nansum(np.nan_to_num(result.energy_stderr) ** 2, axis=(0, 1)))
    return z, prof, err


def zscores(a, ea, b, eb, floor=0.02):
    m = (a > floor * a.max()) & (b > floor * b.max())
    return (a[m] - b[m]) / np.sqrt(ea[m] ** 2 + eb[m] ** 2)


@pytest.fixture(scope="module")
def case():
    geo = homogeneous_box((40.0, 40.0, 100.0), 1.0, "water")
    beam = PencilBeam("proton", 100.0, (0.0, 0.0, -0.5))
    tables = build_tables(SimulationConfig(beam, geo, 1, 1), offline=True)
    grid = ScoringGrid.coarse(geo, (40.0, 40.0, 1.0))
    return geo, beam, tables, grid


def run_case(case, backend, precision, histories, seed, batches=8, physics=None):
    geo, beam, tables, grid = case
    cfg = SimulationConfig(
        beam,
        geo,
        histories=histories,
        batches=batches,
        seed=seed,
        backend=backend,
        precision=precision,
        scoring=grid,
        physics=physics or PhysicsConfig(),
    )
    return run(cfg, tables=tables, offline=True)


def r80(z, prof):
    peak = prof.max()
    i = int(np.where(prof >= 0.8 * peak)[0][-1])
    return z[i] + (0.8 * peak - prof[i]) / (prof[i + 1] - prof[i]) * (z[i + 1] - z[i])


def test_warp_cpu_float32_matches_reference_statistically(case):
    # Small CI case: depth bins of one run are correlated (each proton feeds all
    # bins), so the mean z is only loosely bounded here; the strict criterion of
    # decision 0042 is applied in the local/CUDA tests with many histories.
    ref = run_case(case, "python", "float64", 500, 11)
    cpu = run_case(case, "warp-cpu", "float32", 4000, 12)
    zc, pr, er = depth_profile(ref)
    _, pc, ec = depth_profile(cpu)
    z = zscores(pr, er, pc, ec)
    assert abs(z.mean()) < 1.0 and z.std() < 2.0 and np.abs(z).max() < 5.0, (
        z.mean(),
        z.std(),
    )
    assert abs(r80(zc, pr) - r80(zc, pc)) < 0.3
    assert pr.sum() == pytest.approx(
        pc.sum(), rel=1e-4
    )  # same initial energy, all deposited
    a_ref = ref.metadata["energy_accounting_per_primary"]
    a_cpu = cpu.metadata["energy_accounting_per_primary"]
    assert a_cpu["steps_per_history"] == pytest.approx(
        a_ref["steps_per_history"], rel=0.05
    )
    assert a_cpu["deposited_cutoff_mev"] == pytest.approx(
        a_ref["deposited_cutoff_mev"], rel=0.05
    )


def test_warp_cpu_float64_matches_float32(case):
    f32 = run_case(case, "warp-cpu", "float32", 1500, 21)
    f64 = run_case(case, "warp-cpu", "float64", 1500, 22)
    _, p32, e32 = depth_profile(f32)
    _, p64, e64 = depth_profile(f64)
    z = zscores(p32, e32, p64, e64)
    assert abs(z.mean()) < 1.0 and z.std() < 2.0
    assert p32.sum() == pytest.approx(p64.sum(), rel=1e-5)
    assert f64.metadata["effective"]["precision"] == "float64"


def test_warp_cpu_is_deterministic_for_a_seed(case):
    a = run_case(case, "warp-cpu", "float32", 400, 5)
    b = run_case(case, "warp-cpu", "float32", 400, 5)
    c = run_case(case, "warp-cpu", "float32", 400, 6)
    assert np.array_equal(a.energy_mean, b.energy_mean)
    assert not np.array_equal(a.energy_mean, c.energy_mean)


def test_warp_energy_is_accounted_exactly(case):
    res = run_case(case, "warp-cpu", "float32", 200, 3)
    acc = res.metadata["energy_accounting_per_primary"]
    total = (
        acc["deposited_continuous_mev"]
        + acc["deposited_cutoff_mev"]
        + acc["escaped_mev"]
        + acc["truncated_max_steps_mev"]
    )
    assert total == pytest.approx(acc["initial_mev"], rel=1e-6)
    assert res.energy_mean.sum() == pytest.approx(
        acc["deposited_continuous_mev"] + acc["deposited_cutoff_mev"], rel=1e-5
    )
    assert res.metadata["timing"]["kernel_seconds"] > 0


def test_warp_backend_contract_and_capabilities(case):
    geo, beam, tables, grid = case
    caps = capabilities()
    assert set(caps["warp-cpu"]["precision"]) == {"float32", "float64"}
    with pytest.raises(UnsupportedConfigurationError, match="nuclear"):
        run(
            SimulationConfig(
                beam,
                geo,
                histories=2,
                batches=1,
                backend="warp-cpu",
                physics=PhysicsConfig(nuclear=True),
            ),
            tables=tables,
            offline=True,
        )
    if not wp.is_cuda_available():
        with pytest.raises(UnsupportedConfigurationError, match="CUDA"):
            run(
                SimulationConfig(
                    beam,
                    geo,
                    histories=2,
                    batches=1,
                    backend="warp-cuda",
                    precision="float32",
                ),
                tables=tables,
                offline=True,
            )


def test_warp_cpu_physics_switches_are_discriminating(case):
    geo, beam, tables, _ = case
    lateral_grid = ScoringGrid.coarse(geo, (1.0, 1.0, 100.0))
    full = run(
        SimulationConfig(
            beam,
            geo,
            histories=600,
            batches=4,
            seed=31,
            backend="warp-cpu",
            precision="float32",
            scoring=lateral_grid,
        ),
        tables=tables,
        offline=True,
    )
    no_mcs = run(
        SimulationConfig(
            beam,
            geo,
            histories=600,
            batches=4,
            seed=32,
            backend="warp-cpu",
            precision="float32",
            scoring=lateral_grid,
            physics=PhysicsConfig(multiple_scattering=False),
        ),
        tables=tables,
        offline=True,
    )
    lateral_full = full.energy_mean.sum(axis=2)
    lateral_off = no_mcs.energy_mean.sum(axis=2)
    # without MCS all energy stays in the central column
    assert lateral_off.max() / lateral_off.sum() > 0.999
    assert lateral_full.max() / lateral_full.sum() < 0.95


cuda = pytest.mark.cuda


@cuda
@pytest.mark.skipif(not wp.is_cuda_available(), reason="needs a CUDA device")
def test_warp_cuda_matches_warp_cpu_statistically(case):
    cpu = run_case(case, "warp-cpu", "float32", 20000, 41, batches=10)
    gpu = run_case(case, "warp-cuda", "float32", 100000, 42, batches=10)
    _, pc, ec = depth_profile(cpu)
    _, pg, eg = depth_profile(gpu)
    z = zscores(pc, ec, pg, eg)
    assert abs(z.mean()) < 0.6 and z.std() < 1.3 and np.abs(z).max() < 4.5, (
        z.mean(),
        z.std(),
        np.abs(z).max(),
    )
    assert pc.sum() == pytest.approx(pg.sum(), rel=1e-4)
    # Depth bins are correlated within a run, so the mean z has a standard
    # deviation well above 1/sqrt(n_bins); decision 0042 bounds it at 0.6 and
    # relies on integral quantities for tight checks.
    zc_axis, _, _ = depth_profile(cpu)
    assert abs(r80(zc_axis, pc) - r80(zc_axis, pg)) < 0.15
    plateau = (zc_axis > 10) & (zc_axis < 50)
    assert pc[plateau].mean() == pytest.approx(pg[plateau].mean(), rel=0.005)
    assert gpu.metadata["timing"]["device"].startswith("cuda")


@cuda
@pytest.mark.skipif(not wp.is_cuda_available(), reason="needs a CUDA device")
def test_warp_cuda_float64_matches_float32_and_is_deterministic(case):
    a = run_case(case, "warp-cuda", "float32", 50000, 51, batches=10)
    b = run_case(case, "warp-cuda", "float64", 50000, 52, batches=10)
    _, pa, ea = depth_profile(a)
    _, pb, eb = depth_profile(b)
    z = zscores(pa, ea, pb, eb)
    assert abs(z.mean()) < 0.6 and z.std() < 1.3
    again = run_case(case, "warp-cuda", "float32", 50000, 51, batches=10)
    # atomics make float32 device sums order-dependent at the ulp level; compare to 1e-5
    assert np.allclose(a.energy_mean, again.energy_mean, rtol=1e-5, atol=1e-9)
