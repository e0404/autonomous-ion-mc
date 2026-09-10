"""Warp execution of the shared stopping-power source (skipped without warp).

Tolerances are those of decision 0005 (reference float64 Python versus the
float32 Warp instantiation) and decision 0001 (Warp CPU versus CUDA, by
expectation class), fixed before the first execution.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc import materials, particles
from ionmc.stopping_power import AnalyticStoppingPower

pytestmark = pytest.mark.warp

ENERGIES = np.array([2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 150.0, 200.0, 250.0, 400.0])

# decision 0005: float64 reference vs float32 Warp kernels
REF_VS_WARP_S = {"rtol": 1.0e-5, "atol": 0.0}
REF_VS_WARP_RANGE = {"rtol": 2.0e-5, "atol": 0.0}
# decision 0001 classes for these kernels: transcendental / iterative
CPU_VS_CUDA_S = {"rtol": 4.0e-6, "atol": 1.0e-6}
CPU_VS_CUDA_RANGE = {"rtol": 1.0e-5, "atol": 1.0e-6}


@pytest.fixture(scope="module")
def reference_python():
    return AnalyticStoppingPower(materials.WATER, particles.PROTON, path="python")


def _warp(device: str) -> AnalyticStoppingPower:
    return AnalyticStoppingPower(
        materials.WATER, particles.PROTON, path="warp", device=device
    )


def test_warp_cpu_matches_reference(warp_module, reference_python) -> None:
    cpu = _warp("cpu")
    np.testing.assert_allclose(
        cpu.mass_stopping_power(ENERGIES),
        reference_python.mass_stopping_power(ENERGIES),
        **REF_VS_WARP_S,
    )
    e = np.array([50.0, 100.0, 200.0])
    np.testing.assert_allclose(
        cpu.csda_range(e), reference_python.csda_range(e), **REF_VS_WARP_RANGE
    )


@pytest.mark.cuda
def test_warp_cuda_matches_cpu_and_reference(
    warp_module, cuda_available, reference_python
):
    cpu = _warp("cpu")
    cuda = _warp("cuda:0")
    s_cpu = cpu.mass_stopping_power(ENERGIES)
    s_cuda = cuda.mass_stopping_power(ENERGIES)
    np.testing.assert_allclose(s_cuda, s_cpu, **CPU_VS_CUDA_S)
    np.testing.assert_allclose(
        s_cuda, reference_python.mass_stopping_power(ENERGIES), **REF_VS_WARP_S
    )
    e = np.array([50.0, 100.0, 200.0])
    np.testing.assert_allclose(
        cuda.csda_range(e), cpu.csda_range(e), **CPU_VS_CUDA_RANGE
    )


def test_float32_error_of_stopping_power_stays_small(warp_module, reference_python):
    """Engineering guard behind the decision 0005 criterion (rtol 1e-5).

    The float32 cancellation in ``gamma^2 - 1`` (decision 0006) produced a
    relative error of 8.7e-6; the cancellation-free source gives about 4e-7.
    This bound fails for the former and passes for the latter, so the fix
    cannot be undone unnoticed. It is not an acceptance criterion.
    """
    grid = np.linspace(2.0, 400.0, 400)
    s_cpu = _warp("cpu").mass_stopping_power(grid)
    s_ref = reference_python.mass_stopping_power(grid)
    assert float(np.max(np.abs(s_cpu / s_ref - 1.0))) < 1.0e-6
