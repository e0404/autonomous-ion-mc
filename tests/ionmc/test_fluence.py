"""Tests of energy-resolved fluence-spectrum scoring and lookup accumulation
(decision 0026).

A ``FluenceSpectrum`` histograms the track length ``Sum w*s`` by step-mean energy;
an on-the-fly lookup accumulator ``A_gate = Sum (w*s)*w_tab[bin]`` reproduces the
offline post-processing ``Sum_k counts*w_tab`` to round-off, and with ``w = S_lin``
recovers the deposited (step) energy. The tests cover the spectrum geometry, the
grid-3D requirement, the exact lookup identity, the stopping-power energy
reproduction, spectrum sanity, and reference/Warp CPU/CUDA agreement.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import (
    DepthDoseGrid,
    DoseGrid3D,
    FluenceSpectrum,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
)
from ionmc.transport.depth_dose import DepthLateralGrid


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


def _lat() -> DepthLateralGrid:
    return DepthLateralGrid(200.0, 400, 20.0, 40)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (11, 11, 400), (1.0, 1.0, 0.5), 1.0, origin_mm=(-5.5, -5.5, 0.0)
    )


def _dose() -> DoseGrid3D:
    return DoseGrid3D(
        shape=(11, 11, 100), origin_mm=(-5.5, -5.5, 0.0), spacing_mm=(1.0, 1.0, 2.0)
    )


def test_fluence_spectrum_geometry() -> None:
    fl = FluenceSpectrum(n_bins=160, e_lo_mev=0.0, e_hi_mev=160.0)
    assert fl.bin_width_mev == pytest.approx(1.0)
    assert fl.edges_mev.shape == (161,)
    assert fl.centers_mev.shape == (160,)
    assert fl.centers_mev[0] == pytest.approx(0.5)
    assert fl.empty().shape == (160,)
    # differential fluence scaling: counts * 100 / (V * dE * N)
    counts = np.ones(160)
    phi = fl.differential_fluence(counts, volume_mm3=1000.0, n_histories=2)
    assert phi[0] == pytest.approx(100.0 / (1000.0 * 1.0 * 2))
    # offline lookup is a plain dot product
    tab = np.arange(160, dtype=np.float64)
    assert fl.postprocess_lookup(counts, tab) == pytest.approx(float(tab.sum()))


def test_fluence_requires_grid3d(table) -> None:
    """A FluenceSpectrum is only supported on the VoxelGrid3D path."""
    from ionmc.transport.geometry import WaterSlab

    eng = TransportEngine(table, WaterSlab(200.0), DepthDoseGrid(200.0, 10))
    with pytest.raises(ValueError, match="FluenceSpectrum scorer is supported only"):
        eng.run_scattering(
            PencilBeamSource(150.0),
            _lat(),
            1,
            path="python",
            fluence=FluenceSpectrum(),
        )


def test_lookup_reproduces_offline(table) -> None:
    """The on-the-fly lookup accumulator equals the offline post-processing of the
    scored spectrum to round-off -- the V4 gate (decision 0026)."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(200.0, 10))  # scattering on
    fl = FluenceSpectrum()
    r = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(),
        30,
        seed=7,
        path="python",
        dose_grid=_dose(),
        fluence=fl,
    )
    a_gate = r.fluence_lookup_sum
    a_off = fl.postprocess_lookup(r.fluence_counts_mm, r.fluence_lookup_table)
    assert abs(a_gate - a_off) / abs(a_off) < 1e-12


def test_stopping_lookup_reproduces_energy(table) -> None:
    """With w = S_lin, the lookup accumulator equals the total step energy deposited
    (deterministic mode); it is the total dose minus the terminal energy-cut residual
    and the small binning error, so it is just below the deposited energy."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    dose = _dose()
    r = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(),
        50,
        seed=7,
        path="python",
        dose_grid=dose,
        fluence=FluenceSpectrum(),
    )
    dose_total = float(r.dose3d_mev.sum())
    ratio = r.fluence_lookup_sum / dose_total
    assert 0.99 <= ratio <= 1.0  # step energy: total minus terminal residual


def test_fluence_sanity(table) -> None:
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    fl = FluenceSpectrum()
    r = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(),
        20,
        seed=3,
        path="python",
        dose_grid=_dose(),
        fluence=fl,
    )
    counts = r.fluence_counts_mm
    assert np.all(counts >= 0.0)
    # support lies below the 150 MeV beam energy and above the ~0 cut
    occupied = np.nonzero(counts > 0.0)[0]
    assert fl.centers_mev[occupied].max() <= 150.0
    assert fl.centers_mev[occupied].min() >= 0.0
    # total track length is positive and finite
    assert counts.sum() > 0.0


@pytest.mark.warp
def test_warp_cpu_fluence_matches_reference_deterministic(warp_module, table) -> None:
    """Deterministic +z: the Warp CPU fluence spectrum and lookup accumulator match
    the reference per bin (float32 budget)."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    fl = FluenceSpectrum()
    ref = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(),
        30,
        seed=4,
        path="python",
        dose_grid=_dose(),
        fluence=fl,
    )
    cpu = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(),
        30,
        seed=4,
        path="warp",
        device="cpu",
        dose_grid=_dose(),
        fluence=fl,
    )
    peak = float(ref.fluence_counts_mm.max())
    assert np.max(np.abs(cpu.fluence_counts_mm - ref.fluence_counts_mm)) / peak < 5e-3
    assert (
        abs(cpu.fluence_lookup_sum - ref.fluence_lookup_sum)
        / abs(ref.fluence_lookup_sum)
        < 5e-3
    )


@pytest.mark.cuda
def test_warp_cuda_fluence_matches_cpu(warp_module, cuda_available, table) -> None:
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    fl = FluenceSpectrum()
    cpu = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(),
        30,
        seed=4,
        path="warp",
        device="cpu",
        dose_grid=_dose(),
        fluence=fl,
    )
    cuda = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(),
        30,
        seed=4,
        path="warp",
        device="cuda:0",
        dose_grid=_dose(),
        fluence=fl,
    )
    peak = float(cpu.fluence_counts_mm.max())
    assert np.max(np.abs(cuda.fluence_counts_mm - cpu.fluence_counts_mm)) / peak < 5e-3
