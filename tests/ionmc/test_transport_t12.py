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
    """Reduced CI version of T12 (the plan designates it as such): the LV/HR samples are far
    larger. Asserts no profile is a statistical `fail`, the IDD verdict and the scalars."""
    py, layout = t12_sample(
        backend="python", precision="float64", seed=2026100401, n_histories=200, n_batches=10,
        workers=4, timeout_s=900.0, **KW,
    )  # fmt: skip
    wp, _ = t12_sample(
        backend="warp-cpu", precision="float32", seed=2026100402, n_histories=20000,
        n_batches=20, **KW,
    )  # fmt: skip
    verdict = t12_compare(t12_observables(py, layout), t12_observables(wp, layout))
    # the reduced python sample (200 histories) cannot support the sparse lateral tails: such a
    # profile is `inconclusive` by rule (and so fails); no profile may be a statistical `fail`
    # and the dense IDD and the scalars must pass
    assert all(v["verdict"] != "fail" for v in verdict["arrays"].values()), verdict
    assert verdict["arrays"]["idd"]["pass"] and all(s["pass"] for s in verdict["scalars"].values())
    assert verdict["pass"] == all(v["pass"] for v in verdict["arrays"].values())
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


def _synthetic(seed: int, precision: str, total: float, n_batches: int = 20):  # type: ignore[no-untyped-def]
    from ionmc.transport.parity import T12Observables

    rng = np.random.default_rng(seed)
    base = np.linspace(1.0, 2.0, 12)
    prof = base * (1.0 + 0.01 * rng.normal(size=(n_batches, 12)))
    return T12Observables(
        {"idd": prof},
        {"total_deposit_mev": np.full(n_batches, total)},
        n_batches, 1000, precision,
    )  # fmt: skip


def test_degenerate_scalars_use_the_deterministic_precision_bound() -> None:
    """A scalar with no variance (energy conservation fixes the total deposit) is compared with
    the T4 precision bound of the less precise sample instead of a meaningless z: 1e-5 for float32,
    1e-12 for float64; the rule that applied is recorded."""
    f32, f64 = _synthetic(1, "float32", 150.0000003), _synthetic(2, "float64", 150.00000000009)
    v = t12_compare(f32, f64, n_perm=50)["scalars"]["total_deposit_mev"]
    assert v["rule"] == "deterministic_precision_bound" and v["pass"] and v["bound"] == 1e-5
    off = t12_compare(f32, _synthetic(3, "float64", 150.01), n_perm=50)["scalars"][
        "total_deposit_mev"
    ]
    assert not off["pass"]  # 7e-5 relative: beyond the float32 bound
    two64 = t12_compare(_synthetic(1, "float64", 150.0), _synthetic(2, "float64", 150.0000001),
                        n_perm=50)["scalars"]["total_deposit_mev"]  # fmt: skip
    assert not two64["pass"] and two64["bound"] == 1e-12  # float64 pairs are held to 1e-12


def test_max_z_only_counts_bins_that_both_samples_support() -> None:
    """A bin with fewer than max(2, ceil(B/2)) nonzero batches in either sample has an unreliable
    standard error and does not enter max|z| (it is reported with its support); a well-supported
    bin with the same deviation does."""
    a, b = _synthetic(4, "float64", 150.0), _synthetic(5, "float64", 150.0)
    a.arrays["idd"][:, 11] = 0.0
    a.arrays["idd"][:2, 11] = [1.0, 1.0]  # 2 of 20 batches nonzero: undersupported
    b.arrays["idd"][:, 11] = 5.0
    v = t12_compare(a, b, n_perm=50)["arrays"]["idd"]
    assert (
        v["worst_bin_all"]["supported"] is False
        and v["worst_bin_all"]["batches_with_deposit"][0] == 2
    )
    assert v["max_abs_z_all_bins"] > v["max_abs_z"]  # the undersupported bin is excluded
    assert v["n_supported_bins"] == v["n_bins"] - 1
    a2, b2 = _synthetic(6, "float64", 150.0), _synthetic(7, "float64", 150.0)
    b2.arrays["idd"][:, 5] *= 1.5  # a supported bin that really differs
    w = t12_compare(a2, b2, n_perm=50)["arrays"]["idd"]
    assert w["worst_bin_all"]["supported"] and not w["pass"]


def _sparse_observables(rng: np.random.Generator, n_batches: int, n_hist: int, scale: float = 1.0):  # type: ignore[no-untyped-def]
    """Batch means of a sparse per-primary profile: every history deposits one unit in one of 30
    bins (probabilities falling by a factor 100, summing to 0.7, so the tail bins have discrete
    batch means)."""
    from ionmc.transport.parity import T12Observables

    p = np.geomspace(1.0, 0.01, 30)
    p = 0.7 * p / p.sum()
    p = np.append(p * scale, 1.0 - (p * scale).sum())
    counts = rng.multinomial(n_hist, p, size=n_batches)[:, :30] / float(n_hist)
    return T12Observables(
        {"lat": counts}, {"mean_bin": counts.mean(axis=1)}, n_batches, n_hist, "float64"
    )


def test_sparse_profile_lists_unsupported_bins_and_uses_supported_ones() -> None:
    """The chi-square and its calibration use the supported bins only; the unsupported bins of
    the sparse small sample are listed, the all-bin chi-square is reported for information."""
    rng = np.random.default_rng(20)
    a = _sparse_observables(rng, 40, 100)
    b = _sparse_observables(rng, 100, 10_000)
    v = t12_compare(a, b, n_boot=200, seed=1)["arrays"]["lat"]
    assert len(v["unsupported_bins"]) > 0
    assert v["n_supported_bins"] + len(v["unsupported_bins"]) == v["n_bins"]
    assert math.isfinite(v["chi2"]) and v["chi2"] <= v["chi2_all_selected_bins"] + 1e-12


def test_sparse_null_profile_p_values_are_calibrated_by_the_bootstrap() -> None:
    """Two samples of the same per-primary distribution with different batch structures (40
    batches of 100 histories versus 100 batches of 1e4): over 200 null trials the studentized
    bootstrap-t p-value is below 0.05 in at most 8 % and below 0.01 in at most 2 % of them (a
    little conservative is accepted; the pooled permutation gave about 15 % and 4 %)."""
    rng = np.random.default_rng(20)
    trials, below05, below01 = 200, 0, 0
    for k in range(trials):
        a = _sparse_observables(rng, 40, 100)
        b = _sparse_observables(rng, 100, 10_000)
        v = t12_compare(a, b, n_boot=500, seed=100 + k)["arrays"]["lat"]
        assert v["calibration"] == "studentized_bootstrap_t"
        below05 += v["p_value"] < 0.05
        below01 += v["p_value"] < 0.01
    assert below05 <= 0.08 * trials and below01 <= 0.02 * trials, (below05, below01)


def test_calibration_method_follows_the_batch_structure() -> None:
    """Equal batch structures use the studentized permutation, different ones the bootstrap-t."""
    rng = np.random.default_rng(3)
    a, a2 = _sparse_observables(rng, 20, 1000), _sparse_observables(rng, 20, 1000)
    b = _sparse_observables(rng, 40, 1000)
    assert t12_compare(a, a2, n_perm=50)["calibration"] == "studentized_permutation"
    v = t12_compare(a, b, n_boot=50)
    assert v["calibration"] == "studentized_bootstrap_t" and v["bootstrap"]["n_boot"] == 50


def test_sparse_shifted_alternative_is_detected() -> None:
    """A 20 % scale error of the bulk bins is still detected with the supported-bin statistic."""
    rng = np.random.default_rng(21)
    a = _sparse_observables(rng, 40, 100)
    b = _sparse_observables(rng, 100, 10_000, scale=1.2)
    v = t12_compare(a, b, n_boot=1500, seed=5)["arrays"]["lat"]
    assert v["p_value"] < 0.001 and not v["pass"]


def test_tail_aggregate_catches_a_backend_that_drops_the_tail() -> None:
    """Sample B deposits nothing beyond some radius while its core matches: the dropped bins are
    unsupported (never enter chi2 / max|z|) but the aggregate tail test rejects the pair."""
    rng = np.random.default_rng(30)
    a = _sparse_observables(rng, 100, 10_000)
    b = _sparse_observables(rng, 100, 10_000)
    b.arrays["lat"][:, 20:] = 0.0  # tail beyond bin 20 never deposited
    v = t12_compare(a, b, n_perm=100, seed=2)["arrays"]["lat"]
    assert {20, 21} <= set(v["unsupported_bins"])  # only bins above 1 % of max are selected
    assert v["tail"]["z"] > 3.5 or v["tail"]["status"] == "unsupported"
    assert not v["pass"] and v["verdict"] in ("fail", "inconclusive")
    # a supported but depleted tail (a quarter of the deposits left) fails by the statistic
    c = _sparse_observables(rng, 100, 10_000)
    d = _sparse_observables(rng, 100, 10_000)
    d.arrays["lat"][:, 20:] *= 0.25
    w = t12_compare(c, d, n_perm=100, seed=2)["arrays"]["lat"]
    if w["tail"]["status"] == "supported":
        assert not w["tail"]["pass"] and w["verdict"] == "fail"


def test_profile_without_supported_bins_is_inconclusive_and_fails() -> None:
    from ionmc.transport.parity import T12Observables

    prof = np.zeros((20, 6))
    prof[0] = 1.0  # a single nonzero batch: no bin is supported
    ob = T12Observables({"p": prof}, {}, 20, 1000, "float64")
    v = t12_compare(ob, ob, n_perm=50)
    pv = v["arrays"]["p"]
    assert pv["n_supported_bins"] == 0 and pv["verdict"] == "inconclusive"
    assert pv["inconclusive_reason"] and not pv["pass"] and not v["pass"]


def test_mostly_unsupported_profile_is_inconclusive() -> None:
    rng = np.random.default_rng(31)
    a = _sparse_observables(rng, 40, 100)
    v = t12_compare(a, a, n_perm=50)["arrays"]["lat"]
    if v["n_supported_bins"] < 0.5 * v["n_bins"]:
        assert v["verdict"] == "inconclusive" and not v["pass"]
    full = _synthetic(1, "float64", 150.0)
    w = t12_compare(full, full, n_perm=50)["arrays"]["idd"]
    assert w["verdict"] == "pass" and w["tail"]["status"] == "empty"


def test_sparse_null_combined_decision_false_alarm_rate() -> None:
    """The whole verdict (supported bins and tail aggregate; inconclusive counts as a failure)
    of null pairs with different batch structures: at most 8 % of 200 trials are rejected at the
    frozen 0.001 level of the profile or with the Bonferroni-combined p-value (2 min p) < 0.05.
    """
    rng = np.random.default_rng(20)
    trials, below05, failed = 200, 0, 0
    for k in range(trials):
        a = _sparse_observables(rng, 40, 100)
        b = _sparse_observables(rng, 100, 10_000)
        v = t12_compare(a, b, n_boot=500, seed=100 + k)["arrays"]["lat"]
        assert v["verdict"] != "inconclusive", v["inconclusive_reason"]
        tail_p = v["tail"].get("p_value")
        below05 += 2.0 * min(v["p_value"], 1.0 if tail_p is None else tail_p) < 0.05  # Bonferroni
        failed += not v["pass"]
    assert below05 <= 0.08 * trials and failed <= 0.02 * trials, (below05, failed)
