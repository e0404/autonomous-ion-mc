"""T12 (CI, reduced as frozen): python 200 histories versus warp-cpu 2e4 at 70 MeV, all physics
on, distinct seeds, with the frozen statistics (chi-square p > 0.001, Bonferroni max|z|,
|z| < 3.5 for R80, total deposit and lateral sigma at 0.5 R). Lateral bins are 0.5 mm here (the LV
and HR runs use 0.2 mm; python at 0.2 mm steps is too slow for CI)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ionmc.transport.parity import (
    bonferroni_z,
    t12_compare,
    t12_observables,
    t12_sample,
    wilson_hilferty_p,
)

KW = {"energy_mev": 70.0, "lateral_bin_mm": 0.5, "half_width_mm": 8.0}


def test_statistics_helpers() -> None:
    # Wilson-Hilferty against exact chi-square values: median of chi2(100) ~ 99.33
    assert wilson_hilferty_p(99.33, 100) == pytest.approx(0.5, abs=0.01)
    assert wilson_hilferty_p(149.45, 100) == pytest.approx(0.001, abs=3e-4)
    assert bonferroni_z(1) == pytest.approx(3.2905, abs=1e-3)
    assert bonferroni_z(100) > bonferroni_z(10)


def test_t12_statistics_have_power(make_config) -> None:  # type: ignore[no-untyped-def]
    """A 2 % scale error of one sample is detected (the test can fail)."""
    res, layout = t12_sample(
        backend="warp-cpu", precision="float32", seed=101, n_histories=20000, n_batches=20, **KW
    )
    ob = t12_observables(res, layout)
    scaled = type(ob)(
        {k: v * 1.02 for k, v in ob.arrays.items()},
        {k: v * 1.02 for k, v in ob.scalars.items()},
        ob.n_batches,
    )
    assert not t12_compare(ob, scaled)["pass"]
    assert t12_compare(ob, ob)["pass"]


def test_t12_ci_python_vs_warp_cpu() -> None:
    py, layout = t12_sample(
        backend="python", precision="float64", seed=2026100401, n_histories=200, n_batches=10,
        workers=4, timeout_s=900.0, **KW,
    )  # fmt: skip
    wp, _ = t12_sample(
        backend="warp-cpu", precision="float32", seed=2026100402, n_histories=20000,
        n_batches=20, **KW,
    )  # fmt: skip
    verdict = t12_compare(t12_observables(py, layout), t12_observables(wp, layout))
    assert verdict["pass"], verdict
    for v in verdict["arrays"].values():
        assert v["n_bins"] >= 2 and math.isfinite(v["chi2"])
    assert np.isfinite(verdict["scalars"]["r80_mm"]["z"])
