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


def test_python_scope_call_of_warp_function_matches_reference(
    warp_module, reference_python
):
    """Warp's Python-scope fallback executes the same source (float32 builtins)."""
    from ionmc.physics import stopping as warp_bound

    if warp_bound.m.name != "warp":
        pytest.skip("stopping module not bound to warp in this environment")
    args = reference_python.parameters.as_call_args()
    for e in (10.0, 100.0):
        value = warp_bound.mass_stopping_power(float(e), *args)
        ref = float(reference_python.mass_stopping_power(e)[0])
        assert abs(value / ref - 1.0) < 1e-5


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
