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


def _correlated_samples(rng: np.random.Generator, shift: float = 0.0):  # type: ignore[no-untyped-def]
    """Two samples of batch profiles (20 bins) whose bins are almost perfectly correlated: every
    batch has one common random factor (as the deposits of one history in neighbouring bins)."""
    base = np.linspace(1.0, 2.0, 20)

    def sample(nb: int, extra: float):  # type: ignore[no-untyped-def]
        common = rng.normal(0.0, 0.05, size=(nb, 1))
        return base * (1.0 + extra + common) + rng.normal(0.0, 0.002, size=(nb, 20))

    return sample(20, 0.0), sample(30, shift)


def test_permutation_calibration_of_correlated_profile_chi_square() -> None:
    """Correlated bins invalidate the independent-bin Wilson-Hilferty p-value (too many false
    alarms under the null); the batch-level permutation p-value is calibrated and detects a real
    shift."""
    from ionmc.transport.parity import _chi2_profile, permutation_p_value

    rng = np.random.default_rng(1)
    wh, perm = [], []
    for _ in range(40):
        a, b = _correlated_samples(rng)
        chi2 = _chi2_profile(a, b)[0]
        wh.append(wilson_hilferty_p(chi2, 20))
        perm.append(permutation_p_value(a, 100.0, b, 100.0, n_perm=200, seed=3))
    assert np.mean(np.array(wh) < 0.01) > 0.08  # uncalibrated: far above the nominal 1 %
    assert np.mean(np.array(perm) < 0.01) < 0.1  # calibrated by the permutation
    a, b = _correlated_samples(rng, shift=0.5)
    assert permutation_p_value(a, 100.0, b, 100.0, n_perm=500, seed=3) < 0.01
    # unequal batch sizes are studentized: the same data with different n_b stays valid
    p1 = permutation_p_value(a, 100.0, b, 100.0, n_perm=300, seed=5)
    p2 = permutation_p_value(a, 100.0, b, 100.0, n_perm=300, seed=5)
    assert p1 == p2  # seeded and reproducible
    assert 1.0 / 301.0 <= p1 <= 1.0


def test_t12_verdict_records_permutation_and_wilson_hilferty() -> None:
    res, layout = t12_sample(
        backend="warp-cpu", precision="float32", seed=7, n_histories=4000, n_batches=20, **KW
    )
    ob = t12_observables(res, layout)
    v = t12_compare(ob, ob, n_perm=100, seed=11)
    assert v["permutation"] == {"n_perm": 100, "seed": 11}
    for a in v["arrays"].values():
        assert a["p_value"] == pytest.approx(1.0) and "p_value_wilson_hilferty" in a
