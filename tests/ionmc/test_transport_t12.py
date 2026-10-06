"""T12 (CI, reduced as frozen): python 2000 histories (20 batches) versus warp-cpu 2e4 at 70 MeV,
all physics on, distinct seeds, with the frozen statistics (chi-square p > 0.001, Bonferroni max|z|,
|z| < 3.5 for R80, total deposit and lateral sigma at 0.5 R), evaluated on the grouped profiles
of plan footnote 17: every profile must be a full `pass` (no `inconclusive` is allowed). Lateral
bins are 0.5 mm here (the LV and HR runs use 0.2 mm; python at 0.2 mm steps is too slow for CI).
The synthetic tests below measure the false-alarm rate of the null and the power against
tail-bin defects of the grouped rule."""

from __future__ import annotations

import math
import os

import numpy as np
import pytest
from numpy.typing import NDArray

from ionmc.transport.parity import (
    bonferroni_z,
    t12_compare,
    t12_observables,
    t12_sample,
    wilson_hilferty_p,
)


def _workers(n: int) -> int:
    """Worker processes used only for speed: 1 in the single-process diagnostic mode (histories
    and seeds unchanged)."""
    return 1 if os.environ.get("IONMC_SINGLE_PROCESS") == "1" else n


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
    larger. The full verdict must pass: every profile `pass` (grouped rule), every scalar. The
    python sample is 2000 histories in 20 batches (the pure-Python reference executor is about
    110 times faster; the earlier 200 histories gave 0.5*R inconclusive by the resolution rule,
    plan footnote 17(f)); seeds, bins and the warp-cpu sample are unchanged."""
    py, layout = t12_sample(
        backend="python", precision="float64", seed=2026100401, n_histories=2000, n_batches=20,
        workers=_workers(4), timeout_s=900.0, **KW,
    )  # fmt: skip
    wp, _ = t12_sample(
        backend="warp-cpu", precision="float32", seed=2026100402, n_histories=20000,
        n_batches=20, **KW,
    )  # fmt: skip
    verdict = t12_compare(t12_observables(py, layout), t12_observables(wp, layout))
    detail = {
        k: (v["verdict"], v["groups"], v["p_value"], v["p_value_max_z"], v["inconclusive_reason"])
        for k, v in verdict["arrays"].items()
    }
    assert verdict["pass"], detail
    assert all(v["verdict"] == "pass" for v in verdict["arrays"].values()), detail
    assert all(s["pass"] for s in verdict["scalars"].values())
    assert verdict["calibration"] == "studentized_bootstrap_t"
    assert verdict["max_z_calibration"] == "bootstrap_max_t"
    for v in verdict["arrays"].values():
        assert v["n_bins"] >= 2 and v["n_groups"] >= 1 and math.isfinite(v["chi2"])
    assert np.isfinite(verdict["scalars"]["r80_mm"]["z"])


def test_t12_compare_rejects_vacuous_resolution() -> None:
    ob = _synthetic(1, "float64", 150.0)
    for kw in ({"n_perm": 998}, {"n_boot": 998}):
        with pytest.raises(ValueError, match="999"):
            t12_compare(ob, ob, **kw)


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
    v = t12_compare(ob, ob, n_perm=999, seed=11)
    assert v["permutation"] == {"n_perm": 999, "seed": 11}
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
    v = t12_compare(f32, f64, n_perm=999)["scalars"]["total_deposit_mev"]
    assert v["rule"] == "deterministic_precision_bound" and v["pass"] and v["bound"] == 1e-5
    off = t12_compare(f32, _synthetic(3, "float64", 150.01), n_perm=999)["scalars"][
        "total_deposit_mev"
    ]
    assert not off["pass"]  # 7e-5 relative: beyond the float32 bound
    two64 = t12_compare(_synthetic(1, "float64", 150.0), _synthetic(2, "float64", 150.0000001),
                        n_perm=999)["scalars"]["total_deposit_mev"]  # fmt: skip
    assert not two64["pass"] and two64["bound"] == 1e-12  # float64 pairs are held to 1e-12


def test_deterministic_total_allows_the_fixed_point_rounding_noise() -> None:
    """The grid sum of a sample carries a zero-mean rounding error (~1e-8 MeV per history); two
    float64 samples whose totals differ within a few combined standard errors of that noise pass,
    a systematic offset (beyond 1e-12 relative plus the noise allowance) still fails."""
    from ionmc.transport.parity import T12Observables

    def sample(seed: int, n_batches: int, hpb: int, shift: float) -> T12Observables:
        rng = np.random.default_rng(seed)
        tot = 150.0 + shift + 1.14e-8 / np.sqrt(hpb) * rng.normal(size=n_batches)
        prof = np.linspace(1.0, 2.0, 12) * (1.0 + 0.01 * rng.normal(size=(n_batches, 12)))
        return T12Observables({"idd": prof}, {"total_deposit_mev": tot}, n_batches, hpb, "float64")

    a, b = sample(1, 40, 100, -1.8e-10), sample(2, 100, 10000, 0.0)
    v = t12_compare(a, b, n_perm=999)["scalars"]["total_deposit_mev"]
    assert v["rule"] == "deterministic_precision_bound" and v["pass"]
    bad = t12_compare(a, sample(3, 100, 10000, 5e-8), n_perm=999)["scalars"]["total_deposit_mev"]
    assert not bad["pass"]


def test_noise_allowance_is_restricted_to_the_total_deposit() -> None:
    """Any other scalar (R80, lateral sigma) keeps the frozen |z| < 3.5 rule even with a tiny
    standard error (the 1e-12 + 3.5 SE allowance would pass it); with no variance z is undefined
    and the comparison fails closed unless the values are exactly equal."""
    from ionmc.transport.parity import T12Observables

    def sample(seed: int, r80: float, spread: float, n: int = 20) -> T12Observables:
        rng = np.random.default_rng(seed)
        prof = np.linspace(1.0, 2.0, 12) * (1.0 + 0.01 * rng.normal(size=(n, 12)))
        tot = np.full(n, 150.0)
        return T12Observables(
            {"idd": prof},
            {"total_deposit_mev": tot, "r80_mm": r80 + spread * rng.normal(size=n) if spread
             else np.full(n, r80)},
            n, 1000, "float64",
        )  # fmt: skip

    # tiny standard error (~1e-12): the frozen z rule applies and no allowance is recorded
    v = t12_compare(sample(1, 100.0, 4e-12), sample(2, 100.0 + 1e-6, 4e-12), n_perm=999)
    r80 = v["scalars"]["r80_mm"]
    assert r80["rule"] == "z" and "noise_allowance" not in r80 and not r80["pass"]
    assert abs(r80["z"]) > 1e3
    assert v["scalars"]["total_deposit_mev"]["rule"] == "deterministic_precision_bound"
    # below the former 1e-9 relative "z = 0" shortcut: R80 ~ 100 mm, offset 5e-8 mm, se ~ 1e-12 mm
    low = t12_compare(sample(1, 100.0, 4e-12), sample(2, 100.0 + 5e-8, 4e-12), n_perm=999)
    r80 = low["scalars"]["r80_mm"]
    assert r80["rule"] == "z" and abs(r80["z"]) > 1e3 and not r80["pass"]
    same = t12_compare(sample(1, 100.0, 0.0), sample(2, 100.0, 0.0), n_perm=999)
    assert same["scalars"]["r80_mm"]["rule"] == "z_undefined" and same["scalars"]["r80_mm"]["pass"]
    off = t12_compare(sample(1, 100.0, 0.0), sample(2, 100.5, 0.0), n_perm=999)
    assert (
        off["scalars"]["r80_mm"]["rule"] == "z_undefined" and not off["scalars"]["r80_mm"]["pass"]
    )


def test_stochastic_total_deposit_keeps_the_z_rule() -> None:
    """The deterministic rule needs standard errors <= 1e-8 of the value (gate). A float32 total
    near 100 with se ~ 2e-6 each and a 9e-4 offset is within the 1e-5 relative allowance
    (1e-3) but |z| ~ 3e2: it must fail under the frozen |z| < 3.5 rule; a low-variance pair
    uses the deterministic rule."""
    from ionmc.transport.parity import T12Observables

    def sample(seed: int, shift: float, spread: float) -> T12Observables:
        rng = np.random.default_rng(seed)
        n = 20
        tot = 100.0 + shift + spread * rng.normal(size=n)
        prof = np.linspace(1.0, 2.0, 12) * (1.0 + 0.01 * rng.normal(size=(n, 12)))
        return T12Observables({"idd": prof}, {"total_deposit_mev": tot}, n, 1000, "float32")

    sd = 2e-6 * np.sqrt(20.0)  # standard error 2e-6
    v = t12_compare(sample(1, 0.0, sd), sample(2, 9e-4, sd), n_perm=999)["scalars"]
    r = v["total_deposit_mev"]
    assert r["rule"] == "z" and abs(r["z"]) > 100 and not r["pass"]
    low = t12_compare(sample(1, 0.0, 1e-12), sample(2, 3e-6, 1e-12), n_perm=999)["scalars"]
    assert low["total_deposit_mev"]["rule"] == "deterministic_precision_bound"
    assert low["total_deposit_mev"]["pass"]  # 3e-6 relative 3e-8 < 1e-5


def test_deterministic_scalar_verdict_both_precisions() -> None:
    """The rule of footnote 19 (the gated T12 total deposit; T10 keeps |z| < 3): bound of the less
    precise sample times the scale plus 3.5 combined standard errors."""
    from ionmc.transport.parity import deterministic_scalar_verdict

    rng = np.random.default_rng(5)
    a = 150.0 + 1e-6 * rng.normal(size=20)  # se ~ 2e-7
    for prec, bound in (("float32", 1e-5), ("float64", 1e-12)):
        same = deterministic_scalar_verdict(a, a + 1e-4 * bound * 150.0, prec, prec)
        assert same["pass"] and same["bound"] == bound and "noise_allowance" in same
        assert same["rule"] == "deterministic_precision_bound"
        # offset of the bound plus 20x the noise allowance fails, 0.5x passes
        allow = bound * 150.0 + same["noise_allowance"]
        assert not deterministic_scalar_verdict(a, a + 20 * allow, prec, prec)["pass"]
        assert deterministic_scalar_verdict(a, a + 0.5 * allow, prec, prec)["pass"]
    # the less precise sample sets the bound
    assert deterministic_scalar_verdict(a, a, "float32", "float64")["bound"] == 1e-5
    # zero variance and exact equality
    c = np.full(10, 150.0)
    assert deterministic_scalar_verdict(c, c, "float64", "float64")["pass"]
    assert not deterministic_scalar_verdict(c, c + 1e-6, "float64", "float64")["pass"]


def test_undersupported_bin_is_merged_not_excluded() -> None:
    """A bin with fewer than max(2, ceil(B/2)) nonzero batches in either sample has an unreliable
    standard error: it is not dropped but merged with its inner neighbour into one supported
    group (it is the mode, b's 5.0 dominating the average; the right walk forms no supported
    group, so its remainder joins the innermost left group). The group differs hugely: fail."""
    a, b = _synthetic(4, "float64", 150.0), _synthetic(5, "float64", 150.0)
    a.arrays["idd"][:, 11] = 0.0
    a.arrays["idd"][:2, 11] = [1.0, 1.0]  # 2 of 20 batches nonzero: undersupported
    b.arrays["idd"][:, 11] = 5.0
    v = t12_compare(a, b, n_perm=999)["arrays"]["idd"]
    assert [10, 12] in v["groups"]
    assert v["unsupported_bins"] == [11]
    assert v["worst_bin_all"]["supported"] is False
    assert v["worst_bin_all"]["batches_with_deposit"][0] == 2
    assert v["verdict"] == "fail" and not v["pass"]
    a2, b2 = _synthetic(6, "float64", 150.0), _synthetic(7, "float64", 150.0)
    b2.arrays["idd"][:, 5] *= 1.5  # a supported bin that really differs
    w = t12_compare(a2, b2, n_perm=999)["arrays"]["idd"]
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


def test_sparse_profile_lists_unsupported_bins_and_uses_supported_groups() -> None:
    """The groups are contiguous, ascending, cover the hull exactly and are supported by both
    samples; the unsupported bins of the sparse small sample are listed."""
    from ionmc.transport.parity import _support

    rng = np.random.default_rng(20)
    a = _sparse_observables(rng, 40, 100)
    b = _sparse_observables(rng, 100, 10_000)
    v = t12_compare(a, b, n_boot=999, seed=1)["arrays"]["lat"]
    assert len(v["unsupported_bins"]) > 0
    assert v["n_supported_bins"] + len(v["unsupported_bins"]) == v["n_bins"]
    lo, hi = v["hull"]
    groups = v["groups"]
    assert groups[0][0] == lo and groups[-1][1] == hi + 1
    assert all(g0 < g1 for g0, g1 in groups)
    assert all(groups[k][1] == groups[k + 1][0] for k in range(len(groups) - 1))
    assert len(groups) == v["n_groups"]
    for g0, g1 in groups:
        assert _support(a.arrays["lat"], b.arrays["lat"], g0, g1, 20, 50)
    assert math.isfinite(v["chi2"])


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
        v = t12_compare(a, b, n_boot=999, seed=100 + k)["arrays"]["lat"]
        assert v["calibration"] == "studentized_bootstrap_t"
        below05 += v["p_value"] < 0.05
        below01 += v["p_value"] < 0.01
    assert below05 <= 0.08 * trials and below01 <= 0.02 * trials, (below05, below01)


def test_calibration_method_follows_the_batch_structure() -> None:
    """Equal batch structures use the studentized permutation, different ones the bootstrap-t."""
    rng = np.random.default_rng(3)
    a, a2 = _sparse_observables(rng, 20, 1000), _sparse_observables(rng, 20, 1000)
    b = _sparse_observables(rng, 40, 1000)
    eq = t12_compare(a, a2, n_perm=999)
    assert eq["calibration"] == "studentized_permutation"
    assert eq["max_z_calibration"] == "permutation_max_t"
    assert eq["arrays"]["lat"]["max_z_rule"] == "permutation_max_t"
    v = t12_compare(a, b, n_boot=999)
    assert v["calibration"] == "studentized_bootstrap_t" and v["bootstrap"]["n_boot"] == 999
    assert v["max_z_calibration"] == "bootstrap_max_t"
    assert v["arrays"]["lat"]["max_z_rule"] == "bootstrap_max_t"
    assert v["arrays"]["lat"]["p_value_max_z"] is not None
    assert eq["arrays"]["lat"]["p_value_max_z"] is not None
    assert isinstance(eq["arrays"]["lat"]["bonferroni_informative_pass"], bool)


def test_sparse_shifted_alternative_is_detected() -> None:
    """A 20 % scale error of the bulk bins is still detected with the supported-bin statistic."""
    rng = np.random.default_rng(21)
    a = _sparse_observables(rng, 40, 100)
    b = _sparse_observables(rng, 100, 10_000, scale=1.2)
    v = t12_compare(a, b, n_boot=1500, seed=5)["arrays"]["lat"]
    assert v["p_value"] < 0.001 and not v["pass"]


def test_dropped_tail_is_in_the_hull_and_fails() -> None:
    """Sample B deposits nothing beyond bin 20 while its core matches (100 x 1e4 both): the bins
    above 1 % of a's maximum are selected (union) and the pair fails by the statistic."""
    rng = np.random.default_rng(30)
    a = _sparse_observables(rng, 100, 10_000)
    b = _sparse_observables(rng, 100, 10_000)
    b.arrays["lat"][:, 20:] = 0.0  # tail beyond bin 20 never deposited
    v = t12_compare(a, b, n_perm=999, seed=2)["arrays"]["lat"]
    ma = a.arrays["lat"].mean(axis=0)
    expected = np.nonzero(ma > 0.01 * ma.max())[0]
    assert set(range(20, int(expected.max()) + 1)) <= set(range(v["hull"][0], v["hull"][1] + 1))
    assert not v["pass"] and v["verdict"] == "fail"


def test_profile_without_supported_bins_is_inconclusive_and_fails() -> None:
    from ionmc.transport.parity import T12Observables

    prof = np.zeros((20, 6))
    prof[0] = 1.0  # a single nonzero batch: no group is supported
    ob = T12Observables({"p": prof}, {}, 20, 1000, "float64")
    v = t12_compare(ob, ob, n_perm=999)
    pv = v["arrays"]["p"]
    assert pv["n_supported_bins"] == 0 and pv["verdict"] == "inconclusive"
    assert "no contiguous group" in pv["inconclusive_reason"]
    assert not pv["pass"] and not v["pass"]


def test_mostly_unsupported_profile_is_inconclusive() -> None:
    """Identical sparse samples: the statistic passes (a sample against itself) but less than
    half of the hull's mass lies in single-bin groups: inconclusive. A dense profile passes with
    single-bin groups only."""
    from ionmc.transport.parity import T12Observables

    n_batches = 40
    base = np.geomspace(1.0, 0.012, 30)
    prof = np.zeros((n_batches, 30))
    # bins carry a deposit in exactly 5 of 40 batches (< 20 needed): no bin is supported alone,
    # so every group is a merged one
    for k in range(30):
        prof[(np.arange(5) + 3 * k) % n_batches, k] = base[k] * 8.0
    ob = T12Observables({"lat": prof}, {}, n_batches, 20, "float64")
    v = t12_compare(ob, ob, n_perm=999)["arrays"]["lat"]
    assert v["verdict"] == "inconclusive" and not v["pass"]
    assert v["singleton_mass_fraction"] < 0.5 and v["n_merged_groups"] >= 1
    assert (
        v["p_value"] == pytest.approx(1.0) and "individually supported" in v["inconclusive_reason"]
    )
    full = _synthetic(1, "float64", 150.0)
    w = t12_compare(full, full, n_perm=999)["arrays"]["idd"]
    assert w["verdict"] == "pass" and w["singleton_mass_fraction"] == 1.0
    assert w["n_groups"] == w["n_bins"] and w["n_merged_groups"] == 0


def test_negative_or_nonfinite_profile_values_raise() -> None:
    ob = _synthetic(1, "float64", 150.0)
    bad = _synthetic(2, "float64", 150.0)
    bad.arrays["idd"][3, 4] = -1e-9
    with pytest.raises(ValueError, match="negative or non-finite"):
        t12_compare(ob, bad, n_perm=999)
    bad.arrays["idd"][3, 4] = float("nan")
    with pytest.raises(ValueError, match="negative or non-finite"):
        t12_compare(bad, ob, n_perm=999)


def test_sparse_null_combined_decision_false_alarm_rate() -> None:
    """The whole verdict (grouped chi-square and the conditional max-T criterion; inconclusive
    counts as a failure) of null pairs with different batch structures: at most 2 % of 200 trials
    fail at the frozen 0.001 level of the profile and at most 8 % have the Bonferroni-combined
    p-value (2 min of the chi-square and max-T p-values) below 0.05."""
    rng = np.random.default_rng(20)
    trials, below05, failed = 200, 0, 0
    for k in range(trials):
        a = _sparse_observables(rng, 40, 100)
        b = _sparse_observables(rng, 100, 10_000)
        v = t12_compare(a, b, n_boot=1000, seed=100 + k)["arrays"]["lat"]
        assert v["verdict"] != "inconclusive", v["inconclusive_reason"]
        below05 += 2.0 * min(v["p_value"], v["p_value_max_z"]) < 0.05  # Bonferroni
        failed += not v["pass"]
    assert below05 <= 0.08 * trials and failed <= 0.02 * trials, (below05, failed)


# -- grouped sparse-bin rule (plan footnote 17): synthetic null and defect studies -------------
def _gaussian_bin_probs(n_bins: int, sigma_bins: float) -> NDArray[np.float64]:
    from statistics import NormalDist

    nd = NormalDist(0.0, sigma_bins)
    edges = np.arange(n_bins + 1) - n_bins / 2.0
    cdf = np.array([nd.cdf(float(e)) for e in edges])
    return np.diff(cdf)


def _gaussian_profile_observables(  # type: ignore[no-untyped-def]
    rng: np.random.Generator,
    n_batches: int,
    n_hist: int,
    *,
    n_bins: int = 61,
    sigma_bins: float = 6.0,
    weight_shape: float = 4.0,
    defect=None,
):
    """Batch means per primary of a Gaussian lateral profile: every history hits one bin with the
    normal-CDF probability and deposits a Gamma weight (mean 1, CV 1/sqrt(weight_shape) / sqrt(k)
    for k hits), so the batch values are sparse and skewed in the tails. ``defect`` maps the
    ``(B, n_bins)`` array to a defective one."""
    from ionmc.transport.parity import T12Observables

    p = _gaussian_bin_probs(n_bins, sigma_bins)
    p = np.append(p, 1.0 - p.sum())
    k = rng.multinomial(n_hist, p, size=n_batches)[:, :n_bins]
    x = np.where(k > 0, rng.gamma(np.maximum(k, 1) * weight_shape, 1.0 / weight_shape), 0.0)
    x = x / float(n_hist)
    if defect is not None:
        x = defect(x)
    return T12Observables({"lat": x}, {}, n_batches, n_hist, "float64")


def _outermost_right_bin(p: NDArray[np.float64], fraction: float) -> int:
    return int(np.nonzero(p >= fraction * p.max())[0].max())


def test_grouped_null_false_alarm_unequal() -> None:
    rng = np.random.default_rng(100)
    trials, not_pass, inconclusive, pmax01, informative_fail = 200, 0, 0, 0, 0
    for k in range(trials):
        a = _gaussian_profile_observables(rng, 40, 100)
        b = _gaussian_profile_observables(rng, 100, 2000)
        # 1000 replicates: the smallest attainable p is 1/1001 < 0.001, so the frozen level is
        # evaluable (with 500 replicates the bootstrap criteria could never fail)
        v = t12_compare(a, b, n_boot=1000, seed=100 + k)["arrays"]["lat"]
        not_pass += not v["pass"]
        inconclusive += v["verdict"] == "inconclusive"
        pmax01 += v["p_value_max_z"] < 0.01
        informative_fail += v["bonferroni_informative_pass"] is False
    msg = (not_pass, inconclusive, pmax01, informative_fail)
    assert not_pass <= 4 and inconclusive == 0 and pmax01 <= 6, msg


def test_grouped_null_false_alarm_equal() -> None:
    rng = np.random.default_rng(101)
    trials, not_pass = 200, 0
    for k in range(trials):
        a = _gaussian_profile_observables(rng, 100, 2000)
        b = _gaussian_profile_observables(rng, 100, 2000)
        not_pass += not t12_compare(a, b, n_perm=1000, seed=100 + k)["arrays"]["lat"]["pass"]
    assert not_pass <= 4, not_pass


def test_grouped_null_false_alarm_equal_few_batches() -> None:
    """10 batches of 2000 vs 10 batches of 2000 (the T9/T10 structure): the per-bin z is
    Student-t-like (about 18 degrees of freedom), so the normal Bonferroni bound would reject a
    null profile in several of 200 trials; the permutation max-T p-value does not."""
    rng = np.random.default_rng(102)
    trials, not_pass, normal_fail = 200, 0, 0
    for k in range(trials):
        a = _gaussian_profile_observables(rng, 10, 2000)
        b = _gaussian_profile_observables(rng, 10, 2000)
        v = t12_compare(a, b, n_perm=999, seed=100 + k)["arrays"]["lat"]
        assert v["max_z_rule"] == "permutation_max_t"
        not_pass += not v["pass"]
        normal_fail += v["bonferroni_informative_pass"] is False
    print("few-batch null: not pass", not_pass, "normal Bonferroni would fail", normal_fail)
    assert not_pass <= 4, (not_pass, normal_fail)
    assert normal_fail > not_pass, (not_pass, normal_fail)


def test_dropped_tail_bin_is_detected_equal() -> None:
    """A backend that empties the outermost bin above 1.5 % of the maximum (a 1 % bin could fall
    below the union threshold of the intact sample by sampling noise; the 1.5 % bin is selected
    reliably) is rejected, and that bin lies in the hull."""
    p = _gaussian_bin_probs(61, 6.0)
    k1 = _outermost_right_bin(p, 0.015)
    assert 2000 * 100 * p[k1] >= 50  # design guard: enough hits per sample

    def drop(x: NDArray[np.float64]) -> NDArray[np.float64]:
        x = x.copy()
        x[:, k1] = 0.0
        return x

    rng = np.random.default_rng(102)
    fails = 0
    for t in range(20):
        a = _gaussian_profile_observables(rng, 100, 2000)
        b = _gaussian_profile_observables(rng, 100, 2000, defect=drop)
        v = t12_compare(a, b, n_perm=999, seed=t)["arrays"]["lat"]
        assert v["hull"][0] <= k1 <= v["hull"][1]
        fails += v["verdict"] == "fail"
    assert fails >= 19, fails


def test_mass_shift_between_tail_bins_is_detected_equal() -> None:
    p = _gaussian_bin_probs(61, 6.0)
    k2 = _outermost_right_bin(p, 0.018)

    def shift_inward(x: NDArray[np.float64]) -> NDArray[np.float64]:
        x = x.copy()
        x[:, k2 - 3] += x[:, k2]
        x[:, k2] = 0.0
        return x

    rng = np.random.default_rng(103)
    fails = 0
    for t in range(20):
        a = _gaussian_profile_observables(rng, 100, 2000)
        b = _gaussian_profile_observables(rng, 100, 2000, defect=shift_inward)
        fails += t12_compare(a, b, n_perm=999, seed=t)["arrays"]["lat"]["verdict"] == "fail"
    assert fails >= 19, fails

    def shift_outward(x: NDArray[np.float64]) -> NDArray[np.float64]:
        x = x.copy()
        x[:, k2 + 1] += x[:, k2]
        x[:, k2] = 0.0
        return x

    # documented resolution limit: the emptied bin is unsupported in one sample and is merged
    # with its adjacent outward neighbour, so this defect is not resolved (verdict not asserted)
    a = _gaussian_profile_observables(rng, 100, 2000)
    b = _gaussian_profile_observables(rng, 100, 2000, defect=shift_outward)
    groups = t12_compare(a, b, n_perm=999, seed=1)["arrays"]["lat"]["groups"]
    assert any(g0 <= k2 and k2 + 1 < g1 for g0, g1 in groups), groups


def test_dropped_tail_region_is_detected_unequal() -> None:
    """Sample b (100 x 2000) has no deposit beyond 1.5 sigma on one side (6.7 % of the primaries);
    the small python-like sample (40 x 100) still resolves it."""

    def drop(x: NDArray[np.float64]) -> NDArray[np.float64]:
        x = x.copy()
        x[:, 30 + 9 + 1 :] = 0.0
        return x

    rng = np.random.default_rng(104)
    fails = 0
    for t in range(20):
        a = _gaussian_profile_observables(rng, 40, 100)
        b = _gaussian_profile_observables(rng, 100, 2000, defect=drop)
        fails += t12_compare(a, b, n_boot=1000, seed=t)["arrays"]["lat"]["verdict"] == "fail"
    assert fails >= 19, fails


def test_union_selection_keeps_a_bin_one_sample_dropped() -> None:
    from ionmc.transport.parity import T12Observables, _select_hull

    ma = np.array([0.015, 1.0, 0.5, 0.2])
    mb = np.array([0.0, 1.0, 0.5, 0.2])
    sel, lo, hi = _select_hull(ma, mb)  # type: ignore[misc]
    assert sel[0] and (lo, hi) == (0, 3)
    assert not (0.5 * (ma + mb) > 0.01 * (0.5 * (ma + mb)).max())[0]  # the old average rule
    assert _select_hull(np.zeros(4), np.zeros(4)) is None
    rng = np.random.default_rng(1)
    x = 1.0 + 0.01 * rng.normal(size=(20, 4))
    xa = x * ma
    xb = x * mb
    v = t12_compare(
        T12Observables({"p": xa}, {}, 20, 100, "float64"),
        T12Observables({"p": xb}, {}, 20, 100, "float64"),
        n_perm=999,
    )["arrays"]["p"]
    assert v["hull"] == [0, 3] and v["n_bins"] == 4


def test_profile_decision_truth_table() -> None:
    from ionmc.transport.parity import _profile_decision as d

    assert d(0.5, 0.2) is True  # both calibrated p-values above the level (any structure)
    assert d(0.5, 0.0005) is False  # the max-T p-value rules (the normal z is not an input)
    assert d(0.0005, 0.9) is False
    assert d(0.5, 0.001) is False  # strict: p must exceed 0.001
    assert d(0.5, None) is False
    assert d(float("nan"), 0.9) is False
    assert d(0.5, float("nan")) is False


def test_unequal_pair_normal_bonferroni_exceeded_but_max_t_passes() -> None:
    """The normal Bonferroni bound is anti-conservative for unequal sparse pairs: find a null pair
    (fixed generator) whose max|z| exceeds it although the bootstrap max-T p-value is large; the
    verdict follows the bootstrap."""
    rng = np.random.default_rng(7)
    found = None
    for _ in range(200):
        a = _gaussian_profile_observables(rng, 40, 100, n_bins=121, sigma_bins=12.5)
        b = _gaussian_profile_observables(rng, 100, 10_000, n_bins=121, sigma_bins=12.5)
        v = t12_compare(a, b, n_boot=999, seed=1)["arrays"]["lat"]
        if v["max_abs_z"] >= v["bonferroni_bound"]:
            found = (a, b)
            break
    assert found is not None, "no null pair exceeded the normal Bonferroni bound in 200 draws"
    v = t12_compare(*found, n_boot=1000, seed=1)["arrays"]["lat"]
    assert v["max_abs_z"] >= v["bonferroni_bound"]
    assert v["bonferroni_informative_pass"] is False
    assert v["p_value_max_z"] > 0.001 and v["verdict"] == "pass", v


def test_equal_pair_max_t_failure_is_not_rescued_by_chi2() -> None:
    """Many bins, one bin shifted to z = +5: the chi-square p-value stays above 0.001 (one bin of
    100) but the permutation max-T p-value is at its floor: fail."""
    from ionmc.transport.parity import T12Observables, _chi2_profile, _mean_se

    rng = np.random.default_rng(55)
    xa = 1.0 + 0.01 * rng.normal(size=(20, 100))
    xb = 1.0 + 0.01 * rng.normal(size=(20, 100))
    z5 = _chi2_profile(xa, xb)[2][50]
    _, sa = _mean_se(xa)
    _, sb = _mean_se(xb)
    xb[:, 50] -= (5.0 - z5) * float(np.sqrt(sa[50] ** 2 + sb[50] ** 2))  # z[50] = +5
    v = t12_compare(
        T12Observables({"p": xa}, {}, 20, 100, "float64"),
        T12Observables({"p": xb}, {}, 20, 100, "float64"),
        n_perm=999,
    )["arrays"]["p"]
    assert v["p_value"] > 0.001
    assert v["p_value_max_z"] <= 0.001 and v["verdict"] == "fail"
    assert v["max_abs_z"] >= v["bonferroni_bound"] and v["bonferroni_informative_pass"] is False


def test_unequal_pair_max_t_failure_not_rescued_by_chi2() -> None:
    p = _gaussian_bin_probs(61, 6.0)
    peak = int(np.argmax(p))

    def boost(x: NDArray[np.float64]) -> NDArray[np.float64]:
        x = x.copy()
        x[:, peak] *= 2.0  # +50 % gave max|z| of only 5-8 (not reliably rejected): doubled
        return x

    rng = np.random.default_rng(56)
    a = _gaussian_profile_observables(rng, 40, 100)
    b = _gaussian_profile_observables(rng, 100, 2000, defect=boost)
    v = t12_compare(a, b, n_boot=2000, seed=3)["arrays"]["lat"]
    assert v["p_value_max_z"] <= 0.001, v["p_value_max_z"]
    assert v["verdict"] == "fail"
