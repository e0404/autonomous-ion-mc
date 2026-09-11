"""Tests of dose-averaged LET (LET_d) scoring (decision 0023).

LET_d is scored on a grid co-registered with ``DoseGrid3D`` as
``Sum eps_i * L_i / Sum eps_i`` ("Method C": L_i is the tabulated unrestricted
electronic linear stopping power at the step-mean energy, MeV/mm == keV/um). The
tests cover the ratio/masking helper, the dose-grid requirement, the thin-voxel
analytic limit against NIST PSTAR water, denominator/dose consistency, the
distal-LET depth trend, and reference/Warp CPU/CUDA agreement.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.backend import reference
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import (
    DepthDoseGrid,
    DoseGrid3D,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
)
from ionmc.transport.depth_dose import DepthLateralGrid


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(depth, int(depth * 2), 30.0, 400)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (21, 21, 400), (1.0, 1.0, 0.5), 1.0, origin_mm=(-10.5, -10.5, 0.0)
    )


def _s_lin(table, energy: float, density: float = 1.0) -> float:
    """Reference (numpy-bound) linear stopping power [MeV/mm] at ``energy``."""
    tp = reference.load_bound_module(
        "ionmc.physics.transport",
        "python",
        rebind_dependencies=["ionmc.physics.tabulated"],
    )
    return float(
        tp.linear_stopping_power(
            energy,
            density,
            table.energy_mev,
            table.stopping_mev_cm2_per_g,
            table.slope,
            table.size,
            table.bisection_steps,
        )
    )


def test_let_d_ratio_helper() -> None:
    """``let_d_kev_um`` is num/den where dose>0, with a low-dose mask."""
    g = DoseGrid3D(shape=(2, 1, 2), spacing_mm=(1.0, 1.0, 1.0))
    num = np.array([[[2.0, 0.0]], [[9.0, 100.0]]])  # MeV/mm
    den = np.array([[[1.0, 0.0]], [[3.0, 100.0]]])  # MeV (dose energy)
    letd = g.let_d_kev_um(num, den)
    assert letd[0, 0, 0] == pytest.approx(2.0)  # 2/1
    assert letd[0, 0, 1] == 0.0  # zero dose -> 0
    assert letd[1, 0, 0] == pytest.approx(3.0)  # 9/3
    assert letd[1, 0, 1] == pytest.approx(1.0)  # 100/100
    # low-dose mask: keep only voxels above 5% of the peak dose (=100)
    masked = g.let_d_kev_um(num, den, min_dose_frac=0.05)
    assert masked[0, 0, 0] == 0.0  # dose 1 < 5 -> masked
    assert masked[1, 0, 0] == 0.0  # dose 3 < 5 -> masked
    assert masked[1, 0, 1] == pytest.approx(1.0)  # dose 100 kept


def test_let_requires_dose_grid(table) -> None:
    """``score_let`` without a ``DoseGrid3D`` (the denominator) is rejected."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(200.0, 10))
    with pytest.raises(ValueError, match="LET_d scoring needs a DoseGrid3D"):
        eng.run_scattering(
            PencilBeamSource(100.0), _lat(200.0), 1, path="python", score_let=True
        )


@pytest.mark.parametrize("energy", [250.0, 150.0, 100.0])
def test_thin_voxel_analytic_let(table, energy: float) -> None:
    """Thin-voxel limit: a monoenergetic proton's LET_d -> S_el(E)/rho (decision
    0023). Short steps and a 1 mm dose voxel at entrance keep E ~ constant so the
    dose-averaged LET equals the tabulated linear stopping power."""
    eng = TransportEngine(
        table,
        _box(),
        DepthDoseGrid(200.0, 10),
        straggling=False,
        scattering=False,
        max_step_mm=0.1,
    )
    # a single 1 mm-thick dose voxel centred on the beam axis (x=y=0)
    dose = DoseGrid3D(
        shape=(1, 1, 1), origin_mm=(-1.5, -1.5, 0.0), spacing_mm=(3.0, 3.0, 1.0)
    )
    r = eng.run_scattering(
        PencilBeamSource(energy),
        _lat(200.0),
        1,
        seed=1,
        path="python",
        dose_grid=dose,
        score_let=True,
    )
    letd = dose.let_d_kev_um(r.let3d_num_mev_per_mm, r.dose3d_mev)[0, 0, 0]
    s_lin = _s_lin(table, energy)  # MeV/mm == keV/um
    # slightly above S(E): the step-mean energy is below E and the voxel spans a
    # little energy loss, both raising LET_d; within 1 % at these energies.
    assert letd == pytest.approx(s_lin, rel=0.01)
    assert letd >= s_lin  # dose-averaging over a falling-E voxel never lowers LET


def test_let_denominator_equals_dose_and_positive(table) -> None:
    """The LET_d denominator is the scored dose energy (same eps_i), and the
    numerator is strictly positive wherever dose is deposited."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    dose = DoseGrid3D(
        shape=(21, 21, 100), origin_mm=(-10.5, -10.5, 0.0), spacing_mm=(1.0, 1.0, 2.0)
    )
    r = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(200.0),
        1,
        seed=3,
        path="python",
        dose_grid=dose,
        score_let=True,
    )
    # the denominator used by LET_d IS dose3d_mev; num > 0 exactly where dose > 0
    pos = r.dose3d_mev > 0.0
    assert np.all(r.let3d_num_mev_per_mm[pos] > 0.0)
    assert np.all(r.let3d_num_mev_per_mm[~pos] == 0.0)


def test_distal_let_trend(table) -> None:
    """LET_d rises with depth; its peak lies at or distal to the Bragg dose peak,
    the entrance value is ~0.5 keV/um for 150 MeV, and the distal region reaches
    several keV/um (decision 0023 analytic limits)."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    dose = DoseGrid3D(
        shape=(21, 21, 200), origin_mm=(-10.5, -10.5, 0.0), spacing_mm=(1.0, 1.0, 1.0)
    )
    r = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(200.0),
        1,
        seed=5,
        path="python",
        dose_grid=dose,
        score_let=True,
    )
    zc, ddz = dose.axis_marginals(r.dose3d_mev)
    _, dnz = dose.axis_marginals(r.let3d_num_mev_per_mm)
    letz = np.divide(dnz, ddz, out=np.zeros_like(ddz), where=ddz > 0)
    sig = ddz > 0.01 * ddz.max()  # voxels with appreciable dose
    k_dose = int(np.argmax(ddz))
    k_let = int(np.argmax(np.where(sig, letz, 0.0)))
    assert zc[k_let] >= zc[k_dose]  # LET peak at or distal to the dose peak
    first = int(np.argmax(sig))
    assert letz[first] == pytest.approx(0.5445, rel=0.05)  # entrance ~ PSTAR 150 MeV
    assert letz[k_let] > 3.0  # distal region several keV/um
    # monotone-ish rise: near-peak LET exceeds entrance by a large factor
    assert letz[k_let] > 5.0 * letz[first]


@pytest.mark.warp
def test_warp_cpu_let_matches_reference_deterministic(warp_module, table) -> None:
    """Deterministic +z: the Warp CPU LET_d grid matches the reference per voxel
    (float32 round-off), the tight cross-backend spatial gate (decision 0023)."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    dose = DoseGrid3D(
        shape=(21, 21, 100), origin_mm=(-10.5, -10.5, 0.0), spacing_mm=(1.0, 1.0, 2.0)
    )
    src = PencilBeamSource(150.0)
    lat = _lat(200.0)
    ref = eng.run_scattering(
        src, lat, 1, seed=4, path="python", dose_grid=dose, score_let=True
    )
    cpu = eng.run_scattering(
        src, lat, 1, seed=4, path="warp", device="cpu", dose_grid=dose, score_let=True
    )
    letd_r = dose.let_d_kev_um(ref.let3d_num_mev_per_mm, ref.dose3d_mev)
    letd_c = dose.let_d_kev_um(cpu.let3d_num_mev_per_mm, cpu.dose3d_mev)
    mask = letd_r > 0.0
    max_rel = np.max(np.abs(letd_c[mask] - letd_r[mask]) / letd_r[mask])
    assert max_rel < 5e-3  # float32 budget


@pytest.mark.cuda
def test_warp_cuda_let_matches_cpu_deterministic(
    warp_module, cuda_available, table
) -> None:
    """Deterministic +z: CUDA vs CPU LET_d per voxel to the float32 budget (same
    kernel, non-chaotic straight track -> spatially identical, decision 0023)."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    dose = DoseGrid3D(
        shape=(21, 21, 100), origin_mm=(-10.5, -10.5, 0.0), spacing_mm=(1.0, 1.0, 2.0)
    )
    src = PencilBeamSource(150.0)
    lat = _lat(200.0)
    cpu = eng.run_scattering(
        src, lat, 1, seed=4, path="warp", device="cpu", dose_grid=dose, score_let=True
    )
    cuda = eng.run_scattering(
        src,
        lat,
        1,
        seed=4,
        path="warp",
        device="cuda:0",
        dose_grid=dose,
        score_let=True,
    )
    letd_c = dose.let_d_kev_um(cpu.let3d_num_mev_per_mm, cpu.dose3d_mev)
    letd_g = dose.let_d_kev_um(cuda.let3d_num_mev_per_mm, cuda.dose3d_mev)
    mask = letd_c > 0.0
    assert np.max(np.abs(letd_g[mask] - letd_c[mask]) / letd_c[mask]) < 5e-3
