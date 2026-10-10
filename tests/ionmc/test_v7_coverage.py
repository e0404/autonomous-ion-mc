# ruff: noqa: E501
"""V7 replicate-coverage rule (plan Amendment 13, Codex REVIEW-be30e621 ... -22711e6f): the reference-free
paired design (3600 disjoint pairs of independent replicates, hit = |difference| <= sqrt(sum of squared SEMs),
nominal c0 = P(|t_38| <= 1), region c0 +- 0.03), the bin mask from the independent 1e6 reference only, exact
finite-sample Clopper-Pearson (scalar) / empirical Bernstein (profiles) bounds at alpha_tost = 0.04, combine of
shard/reference partials (fail closed), and (marker ``calibration``, minutes, not in CI) Monte Carlo calibration
of the whole procedure with Clopper-Pearson bounds of the false-acceptance rates: Gaussian profiles and scalar,
Gamma-block tallies with independent, correlated (shared-history) and zero-inflated bins."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"


def _load(name: str) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def v5b() -> ModuleType:
    return _load("steps_v5b")


def test_student_quantiles_and_nominal_coverage(v5b: ModuleType) -> None:
    for p, df, ref in ((0.95, 19, 1.7291), (0.95, 29, 1.6991), (0.95, 4499, 1.6453)):
        assert abs(v5b.student_t_quantile(p, df) - ref) < 5e-4
    assert abs(v5b.student_abs_prob(19, 1.0) - 0.670) < 0.002
    assert abs(v5b.student_abs_prob(1, 1.0) - 0.5) < 1e-9  # Cauchy: (2/pi) atan(1)


def test_region_c0_and_constants(v5b: ModuleType) -> None:
    c0, low, high = v5b.v7_region()
    assert (v5b.V7_PAIR_DF, v5b.V7_COV_MARGIN, v5b.V7_ALPHA_TOST) == (38, 0.03, 0.04)
    assert abs(c0 - v5b.student_abs_prob(38, 1.0)) < 1e-15 and abs(c0 - 0.6764) < 5e-4
    assert abs(low - (c0 - 0.03)) < 1e-15 and abs(high - (c0 + 0.03)) < 1e-15
    assert (round(c0, 4), round(low, 4), round(high, 4)) == (0.6764, 0.6464, 0.7064) or abs(
        c0 - 0.6764
    ) < 5e-4
    assert v5b.V7_REPLICATES == 7200 and v5b.V7_SHARDS == 8 and v5b.V7_REPS_PER_SHARD == 900


def test_equivalence_verdict_vectors(v5b: ModuleType) -> None:
    c0, low, high = v5b.v7_region()
    ok = v5b.v7_tost_verdict(low + 0.005, high - 0.005, 12, 3600, min_bins=10)
    assert (
        ok["pass"]
        and not ok["reasons"]
        and ok["c0"] == c0
        and ok["low"] == low
        and ok["high"] == high
    )
    lo = v5b.v7_tost_verdict(low - 1e-4, high - 0.005, 12, 3600, min_bins=10)
    assert not lo["pass"] and any("lower" in x for x in lo["reasons"])
    hi = v5b.v7_tost_verdict(low + 0.005, high + 1e-4, 12, 3600, min_bins=10)
    assert not hi["pass"] and any("upper" in x for x in hi["reasons"])
    assert v5b.v7_tost_verdict(low, high, 1, 3600)["pass"]  # closed region
    few = v5b.v7_tost_verdict(low + 0.005, high - 0.005, 12, 299, min_bins=10)
    assert not few["pass"] and any("pairs" in x for x in few["reasons"])
    bins = v5b.v7_tost_verdict(low + 0.005, high - 0.005, 9, 3600, min_bins=10)
    assert not bins["pass"] and any("bins_used" in x for x in bins["reasons"])


def test_clopper_pearson_known_values_and_binomial_consistency(v5b: ModuleType) -> None:
    """n = 10, k = 3, one-sided 95 %: [0.0873, 0.6066] (the two-sided 90 % Clopper-Pearson interval)."""
    assert (
        abs(v5b.cp_lower(3, 10, 0.05) - 0.08726) < 5e-5
        and abs(v5b.cp_upper(3, 10, 0.05) - 0.60662) < 5e-5
    )
    assert v5b.cp_lower(0, 10, 0.05) == 0.0 and v5b.cp_upper(10, 10, 0.05) == 1.0
    assert abs(v5b.cp_upper(0, 10, 0.05) - (1.0 - 0.05**0.1)) < 1e-9
    assert abs(v5b.cp_lower(10, 10, 0.05) - 0.05**0.1) < 1e-9
    for k, n, a in ((7, 40, 0.04), (2435, 3600, 0.04), (5, 3600, 0.04)):
        lo, up = v5b.cp_lower(k, n, a), v5b.cp_upper(k, n, a)
        assert lo < k / n < up

        def tail(p: float, js: range, n: int = n) -> float:
            return sum(
                math.exp(math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1)
                         + j * math.log(p) + (n - j) * math.log1p(-p))
                for j in js
            )  # fmt: skip

        assert (
            abs(tail(lo, range(k, n + 1)) - a) < 1e-8 and abs(tail(up, range(0, k + 1)) - a) < 1e-8
        )


def test_empirical_bernstein_hand_value(v5b: ModuleType) -> None:
    v, n, a = 0.0182, 3600, 0.04
    ln = math.log(50.0)
    want = math.sqrt(2 * v * ln / n) + 7 * ln / (3 * (n - 1))
    assert abs(v5b.empirical_bernstein_eps(v, n, a) - want) < 1e-15 and abs(want - 0.008826) < 1e-5
    assert v5b.empirical_bernstein_eps(0.0, n, a) == 7 * ln / (3 * (n - 1))  # variance-free term


def _scalar_acceptance(v5b: ModuleType, n: int = 3600) -> tuple[int, int]:
    """Hit counts ``kl <= k <= kc`` accepted by the scalar rule: ``CP_lower(k) >= low`` (monotone in k)
    and ``CP_upper(k) <= high`` (monotone in k), at ``alpha_tost`` = 0.04, ``n`` pairs."""
    _, low, high = v5b.v7_region()
    lo, hi = 0, n
    while lo < hi:  # smallest k with CP_lower(k) >= low
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if v5b.cp_lower(mid, n, 0.04) >= low else (mid + 1, hi)
    kl = lo
    lo, hi = 0, n
    while lo < hi:  # largest k with CP_upper(k) <= high
        mid = (lo + hi + 1) // 2
        lo, hi = (mid, hi) if v5b.cp_upper(mid, n, 0.04) <= high else (lo, mid - 1)
    return kl, lo


def _binom_interval_prob(n: int, p: float, a: int, b: int) -> float:
    return sum(
        math.exp(math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1)
                 + j * math.log(p) + (n - j) * math.log1p(-p))
        for j in range(a, b + 1)
    )  # fmt: skip


def test_scalar_exact_false_acceptance_at_both_margins(v5b: ModuleType) -> None:
    """The exact false-acceptance probability of the scalar rule (n = 3600 pairs): the probability that
    the hit count lies in the acceptance set {k : CP_lower(k) >= low and CP_upper(k) <= high} when the
    true paired coverage is p = low and p = high, evaluated with the binomial distribution. It is the
    gate of the scalar (<= alpha_tost = 0.04; about 0.0387); the Monte Carlo calibration is a
    consistency check of the implementation."""
    _, low, high = v5b.v7_region()
    kl, kc = _scalar_acceptance(v5b)
    assert 2300 < kl < kc < 2600 and v5b.cp_lower(kl, 3600, 0.04) >= low > v5b.cp_lower(
        kl - 1, 3600, 0.04
    )
    assert v5b.cp_upper(kc, 3600, 0.04) <= high < v5b.cp_upper(kc + 1, 3600, 0.04)
    fa_low = _binom_interval_prob(3600, low, kl, kc)
    fa_high = _binom_interval_prob(3600, high, kl, kc)
    print(
        "v7 scalar exact false acceptance",
        {"kl": kl, "kc": kc, "at_low": fa_low, "at_high": fa_high},
    )
    assert (
        fa_low <= 0.04
        and fa_high <= 0.04
        and abs(fa_low - 0.0387) < 5e-4
        and abs(fa_high - 0.0387) < 5e-4
    )


def test_pairing_is_fixed_and_disjoint(v5b: ModuleType) -> None:
    ia, ib = v5b.v7_pair_indices(7200)
    assert ia.size == ib.size == 3600
    allrows = np.concatenate([ia, ib])
    assert sorted(allrows.tolist()) == list(range(7200))  # every replicate in exactly one pair
    per = v5b.V7_REPS_PER_SHARD
    assert (ib // per == ia // per + 4).all() and (
        ia % per == ib % per
    ).all()  # (s, j) with (s + 4, j)
    assert (ia // per < 4).all() and (ia[:900] // per == 0).all()
    with pytest.raises(SystemExit, match="paired"):
        v5b.v7_pair_indices(7)


def test_paired_hit_rule_hand_case_scalar(v5b: ModuleType) -> None:
    """8 replicates -> 4 pairs (j, j + 4); scalar bin; hit iff |x_j - x_k| <= sqrt(s_j^2 + s_k^2)."""
    means = np.array([[1.00], [1.20], [0.80], [1.50], [1.10], [1.00], [0.20], [1.60]])
    sems = np.array([[0.10], [0.10], [0.20], [0.30], [0.10], [0.20], [0.20], [0.20]])
    e = v5b.v7_estimator_verdict("scalar", means, sems, np.array([1.0]), np.array([0.1]), None)
    # pairs (1.00,1.10): |.10| <= sqrt(.01+.01)=.1414 hit; (1.20,1.00): .20 <= sqrt(.01+.04)=.2236 hit;
    # (0.80,0.20): .60 <= sqrt(.04+.04)=.2828 miss; (1.50,1.60): .10 <= sqrt(.09+.04)=.3606 hit
    assert (
        e["per_pair_covered"] == [1, 1, 0, 1]
        and e["n_pairs"] == 4
        and e["covered"] == 3
        and e["m"] == 0.75
    )
    assert (
        e["bound_kind"] == "clopper-pearson" and e["bins_used"] == 1 and not e["pass"]
    )  # 4 pairs < 300
    assert (
        abs(e["ci"][0] - v5b.cp_lower(3, 4, 0.04)) < 1e-12
        and abs(e["ci"][1] - v5b.cp_upper(3, 4, 0.04)) < 1e-12
    )
    assert any("pairs" in r for r in e["reasons"]) and e["c0"] == v5b.v7_region()[0]
    # single replicates against the 1e6 reference 1.0: |x - 1| <= sem -> 1.0 yes, 1.2 no, .8 no(.2<=.2 yes), ...
    want = np.mean(np.abs(means[:, 0] - 1.0) <= sems[:, 0])
    assert abs(e["m_ref1e6"] - want) < 1e-15 and e["legacy_point_gate_pass"] is False


def test_degenerate_uncertainty_is_a_miss_and_is_counted(v5b: ModuleType) -> None:
    """A pair-bin with sem == 0 for either replicate (blocks identical, e.g. all zero) is a miss even if
    |difference| <= 0; counts are reported per bin, in total and as a fraction of pairs."""
    means = np.array(
        [
            [0.0, 5.0],
            [0.0, 5.0],
            [1.0, 5.0],
            [0.0, 5.0],
            [0.0, 5.0],
            [1.0, 5.0],
            [2.0, 5.0],
            [0.0, 5.0],
        ]
    )
    sems = np.array(
        [
            [0.0, 1.0],
            [0.0, 1.0],
            [0.5, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.5, 1.0],
            [0.5, 1.0],
            [0.3, 1.0],
        ]
    )
    ref = np.array([1.0, 5.0])
    e = v5b.v7_estimator_verdict("profile", means, sems, ref, np.array([0.01, 0.01]), None)
    # pairs (j, j + 4): (0,4) both sem 0 mean 0 -> degenerate miss; (1,5): sem 0 and .5 -> degenerate miss;
    # (2,6): |1-2| = 1 > sqrt(.25+.25) = .707 -> ordinary miss; (3,7): sem 0 and .3 -> degenerate miss
    assert e["degenerate_pair_bins"] == [3, 0] and e["degenerate_pairs_total"] == 3
    assert e["degenerate_pairs"] == 3 and e["degenerate_pair_fraction"] == 0.75
    assert e["per_pair_covered"] == [
        1,
        1,
        1,
        1,
    ]  # only bin 1 (identical means, sem 1) hits; bin 0 never
    assert e["per_bin_coverage"] == [0.0, 1.0]
    ok = np.zeros((8, 1))
    z = v5b.v7_estimator_verdict(
        "scalar", ok, np.zeros((8, 1)), np.array([1.0]), np.array([0.1]), None
    )
    assert (
        z["covered"] == 0
        and z["degenerate_pairs_total"] == 4
        and z["degenerate_pair_fraction"] == 1.0
    )
    assert (
        z["single"]["degenerate_total"] == 6 and z["single"]["m"] == 0.0
    )  # single-interval gate: all misses


def test_paired_region_maps_to_a_single_interval_coverage_band(v5b: ModuleType) -> None:
    """Gaussian batch means with error scale kappa: single-interval coverage P(|t_19| <= kappa), paired
    coverage P(|t_38| <= kappa), both strictly increasing; the paired region [c0 -/+ 0.03] is a kappa
    band and a single-interval coverage band (the justified link for Gaussian tallies)."""
    c0, low, high = v5b.v7_region()
    ks = np.linspace(0.5, 1.6, 23)
    single = np.array([v5b.student_abs_prob(19, k) for k in ks])
    paired = np.array([v5b.student_abs_prob(38, k) for k in ks])
    assert (np.diff(single) > 0).all() and (np.diff(paired) > 0).all()

    def kappa_for(df: int, cov: float) -> float:
        a, b = 0.1, 5.0
        for _ in range(80):
            m = 0.5 * (a + b)
            a, b = (m, b) if v5b.student_abs_prob(df, m) < cov else (a, m)
        return 0.5 * (a + b)

    k_lo, k_hi = kappa_for(38, low), kappa_for(38, high)
    s_lo, s_c, s_hi = (v5b.student_abs_prob(19, k) for k in (k_lo, 1.0, k_hi))
    print("v7 paired region -> kappa band and single-interval coverage band",
          {"kappa": (k_lo, k_hi), "single": (s_lo, s_c, s_hi), "paired": (low, c0, high)})  # fmt: skip
    assert k_lo < 1.0 < k_hi and abs(kappa_for(38, c0) - 1.0) < 1e-9
    assert abs(s_c - 0.670) < 0.002 and s_lo < s_c < s_hi
    assert abs(s_lo - 0.640) < 0.01 and abs(s_hi - 0.700) < 0.01  # close to 0.670 -/+ 0.03


def test_resolvability_criterion_excludes_exactly_the_sparse_bin(v5b: ModuleType) -> None:
    """Bin used iff ref mean > 0, spread > 0 and r_b = (SEM_ref / ref_mean) sqrt(n_ref / n_rep) < 0.5: a
    synthetic reference with one sparse bin excludes exactly that bin; the n_ref / n_rep scaling holds."""
    assert v5b.V7_MAX_REL_SE == 0.5
    rng = np.random.default_rng(8)
    means = 5.0 + rng.normal(0.0, 1.0, (16, 4))
    sems = np.full_like(means, 0.8)
    ref = np.array([5.0, 5.0, 5.0, 5.0])
    rsem = np.array([0.1, 0.1, 0.1, 0.1])
    rsem[2] = 0.3  # sparse bin: rel SEM_ref = 0.06 -> r = 0.6 at n_ref / n_rep = 100
    e = v5b.v7_estimator_verdict(
        "profile", means, sems, ref, rsem, None, n_ref=1_000_000, n_rep=10_000
    )
    assert e["bin_mask"] == [True, True, False, True] and e["bins_used"] == 3
    assert np.allclose(e["rel_se_per_replicate"], [0.2, 0.2, 0.6, 0.2])
    assert (
        list(e["excluded_bins"]) == [2]
        and "relative SE" in e["excluded_bins"][2]
        and e["max_rel_se"] == 0.5
    )
    # a smaller reference run (reduced rehearsal: n_ref / n_rep = 2) scales r_b by sqrt(2 / 100)
    e2 = v5b.v7_estimator_verdict(
        "profile", means, sems, ref, rsem, None, n_ref=20_000, n_rep=10_000
    )
    assert np.allclose(
        e2["rel_se_per_replicate"], np.array([0.2, 0.2, 0.6, 0.2]) * math.sqrt(2.0 / 100.0)
    )
    assert e2["bin_mask"] == [True] * 4
    # just below the cut a bin is used, exactly at the cut (strict inequality) it is excluded
    just = np.array([0.1, 0.1, 0.1, 0.1])
    just[0], just[1] = 0.249, 0.25  # r = 0.498 and 0.5
    e3 = v5b.v7_estimator_verdict(
        "profile", means, sems, ref, just, None, n_ref=1_000_000, n_rep=10_000
    )
    assert e3["bin_mask"] == [True, False, True, True]
    # ref mean 0 and zero spread keep their own reasons
    ref4 = np.array([5.0, 0.0, 5.0, 5.0])
    blocks = np.array([[5.0, 0.0, 4.0 + i % 3, 5.0 + (i % 3) * 0.1] for i in range(20)])
    e4 = v5b.v7_estimator_verdict(
        "profile", means, sems, ref4, rsem, None, n_ref=1_000_000, n_rep=10_000
    )
    assert e4["bin_mask"][1] is False and "reference mean" in e4["excluded_bins"][1]
    e5 = v5b.v7_estimator_verdict(
        "profile", means, sems, ref4, rsem, blocks, n_ref=1_000_000, n_rep=10_000
    )
    assert 1 in e5["excluded_bins"] and e5["bins_used"] <= 3
    # the scalar has no resolvability criterion
    sc = v5b.v7_estimator_verdict(
        "scalar", means[:, :1], sems[:, :1], np.array([5.0]), np.array([4.0]), None
    )
    assert sc["bin_mask"] == [True] and sc["excluded_bins"] == {}


def test_single_interval_region_and_hoeffding_radius(v5b: ModuleType) -> None:
    """One-sample target c1 = P(|t_19| <= 1) = 0.6701, region c1 -/+ 0.03; the fixed Hoeffding radius of
    the single gate for n = 5400 and alpha = 0.04 is sqrt(ln 25 / 10800) = 0.01726."""
    c1, low, high = v5b.v7_single_region()
    assert abs(c1 - v5b.student_abs_prob(19, 1.0)) < 1e-15 and abs(c1 - 0.6701) < 2e-4
    assert abs(low - (c1 - 0.03)) < 1e-15 and abs(high - (c1 + 0.03)) < 1e-15
    print("v7 single-interval nominal", {"c1": c1, "low": low, "high": high})
    eps = v5b.v7_hoeffding_radius(5400, 0.04)
    assert abs(eps - math.sqrt(math.log(25.0) / 10800)) < 1e-15 and abs(eps - 0.017256) < 1e-5
    assert not hasattr(
        v5b, "v7_welch"
    )  # the Welch width is gone: the interval is mean +- sem alone


def test_coverage_sweep_matches_brute_force(v5b: ModuleType) -> None:
    """The exact sweep equals a brute-force evaluation of C(delta) at all clipped breakpoints, at +-Z and
    at the midpoints between them (a fine grid cannot beat it); the returned shifts attain the extrema;
    left-out (degenerate) replicates count as misses through n_total."""
    rng = np.random.default_rng(11)
    for _ in range(200):
        n = int(rng.integers(1, 9))
        x = rng.normal(0.0, 2.0, n)
        half = rng.uniform(0.05, 2.0, n)
        z = float(rng.uniform(0.3, 4.0))
        lo, hi = x - half, x + half
        pts = np.unique(
            np.concatenate([[-z, z], lo[(lo > -z) & (lo < z)], hi[(hi > -z) & (hi < z)]])
        )
        cand = np.concatenate([pts, 0.5 * (pts[:-1] + pts[1:])])

        def cov(d: float, x=x, half=half) -> float:  # noqa: ANN001
            return float(np.mean(np.abs(x - d) <= half))

        brute = [cov(d) for d in cand]
        mn, mx, dmn, dmx = v5b.v7_coverage_extrema(x, half, z)
        assert mn == min(brute) and mx == max(brute)
        assert cov(dmn) == mn and cov(dmx) == mx and -z <= dmn <= z and -z <= dmx <= z
        fine = [cov(d) for d in np.linspace(-z, z, 4001)]
        assert mn <= min(fine) + 1e-12 and mx >= max(fine) - 1e-12
    mn, mx, _, _ = v5b.v7_coverage_extrema(np.array([0.0]), np.array([1.0]), 0.5, 4)  # 1 valid of 4
    assert (mn, mx) == (0.25, 0.25)
    assert v5b.v7_coverage_extrema(np.array([]), np.array([]), 2.0, 5) == (0.0, 0.0, 0.0, 0.0)


def test_box_quantiles_and_bootstrap(v5b: ModuleType) -> None:
    """Z_t = t_{1 - alpha_box / (2 B), 1799}: about 3.34 for 12 bins and 2.58 for one; the bootstrap
    quantile of Gaussian replicate means is close to it and reproducible for a fixed seed."""
    z12 = v5b.student_t_quantile(1.0 - 0.01 / 24.0, 1799)
    z1 = v5b.student_t_quantile(1.0 - 0.01 / 2.0, 1799)
    assert abs(z12 - 3.3415) < 0.03 and abs(z1 - 2.5758) < 0.02
    rng = np.random.default_rng(5)
    r, nh, b = 4000, 1000, 3
    means = 10.0 + rng.normal(0.0, 1.0, (r, b))
    sem_h = means[-nh:].std(axis=0, ddof=1) / math.sqrt(nh)
    z1_ = v5b.v7_boot_z(means, nh, means.mean(axis=0), sem_h, 7, 4000)
    z2_ = v5b.v7_boot_z(means, nh, means.mean(axis=0), sem_h, 7, 4000)
    zt = v5b.student_t_quantile(1.0 - 0.01 / (2 * b), nh - 1)
    assert z1_ == z2_ and abs(z1_ - zt) < 0.3, (z1_, zt)


def test_single_interval_gate_hand_case_and_heldout_independence(v5b: ModuleType) -> None:
    """Scalar, 8 replicates: E = first 6, H = last 2. Hit iff |x - mu_H| <= sem_j ALONE (no SEM_H in the
    width); the hit interval in the shift variable is [(x - sem) / SEM_H, (x + sem) / SEM_H]; nothing
    from E enters mu_H, SEM_H or the bin mask; row pass needs both gates."""
    means = np.array([[1.02], [1.21], [0.83], [0.98], [0.91], [1.13], [1.00], [1.20]])
    sems = np.full((8, 1), 0.10)
    ref, rsem = np.array([1.0]), np.array([0.01])
    e = v5b.v7_estimator_verdict("scalar", means, sems, ref, rsem, None, boot_seed=3, n_boot=300)
    sg = e["single"]
    mu_h = 1.1
    sem_h = float(np.std([1.0, 1.2], ddof=1)) / math.sqrt(2.0)  # 0.1
    assert abs(sg["mu_heldout"][0] - mu_h) < 1e-15 and abs(sg["sem_heldout"][0] - sem_h) < 1e-15
    want = (np.abs(means[:6, 0] - mu_h) <= 0.10).astype(
        int
    )  # 0.08 y, 0.11 n, 0.27 n, 0.12 n, 0.19 n, 0.03 y
    assert (
        sg["per_replicate_covered"] == want.tolist() == [1, 0, 0, 0, 0, 1]
        and abs(sg["m"] - 2 / 6) < 1e-15
    )
    c1, low, high = v5b.v7_single_region()
    assert sg["c1"] == c1 and sg["low"] == low and sg["high"] == high and "nu" not in sg
    zt = v5b.student_t_quantile(1.0 - 0.01 / 2.0, 1)
    assert (
        abs(sg["z_t"] - zt) < 1e-12
        and sg["z_box"] == max(sg["z_t"], sg["z_boot"])
        and sg["boot_seed"] == 3
    )
    x = (means[:6, 0] - mu_h) / sem_h
    mn, mx, _, _ = v5b.v7_coverage_extrema(
        x, np.full(6, 0.10 / sem_h), sg["z_box"]
    )  # half width sem / SEM_H
    assert (
        sg["per_bin_box_coverage_min_max"] == [[mn, mx]] and sg["m_lo"] == mn and sg["m_hi"] == mx
    )
    assert sg["bound_kind"] == "clopper-pearson" and sg["cp_count_min"] == round(mn * 6)
    assert abs(sg["ci"][0] - v5b.cp_lower(sg["cp_count_min"], 6, 0.04)) < 1e-12
    assert (
        e["pass"] == (e["pass_paired"] and e["pass_single"]) and not e["pass"]
    )  # 6 < 300 replicates
    assert any(r.startswith("single:") for r in e["reasons"]) and "replicates 6 < 300" in " ".join(
        e["reasons"]
    )
    m2 = means.copy()
    m2[2, 0] += 0.3  # an evaluation replicate does not move the reference
    e2 = v5b.v7_estimator_verdict("scalar", m2, sems, ref, rsem, None, boot_seed=3, n_boot=300)
    assert (
        e2["single"]["mu_heldout"] == sg["mu_heldout"]
        and e2["single"]["sem_heldout"] == sg["sem_heldout"]
    )
    assert e2["bin_mask"] == e["bin_mask"]
    m3 = means.copy()
    m3[6, 0] += 0.3  # a held-out replicate does
    e3 = v5b.v7_estimator_verdict("scalar", m3, sems, ref, rsem, None, boot_seed=3, n_boot=300)
    assert abs(e3["single"]["mu_heldout"][0] - (mu_h + 0.15)) < 1e-12
    s4 = sems.copy()
    s4[0, 0] = 0.0  # degenerate sem = 0 is a miss
    e4 = v5b.v7_estimator_verdict("scalar", means, s4, ref, rsem, None, boot_seed=3, n_boot=300)
    assert e4["single"]["per_replicate_covered"][0] == 0 and e4["single"]["degenerate_total"] == 1
    m5 = means.copy()
    m5[6:, 0] = 1.0  # identical H means: SEM_H = 0, no information: all misses
    e5 = v5b.v7_estimator_verdict("scalar", m5, sems, ref, rsem, None, boot_seed=3, n_boot=300)
    assert e5["single"]["m"] == 0.0 and e5["single"]["degenerate_total"] == 6


def test_single_gate_monotonicity_m_lo_m_hi_bracket_every_shift(v5b: ModuleType) -> None:
    """For any shift vector inside the box the mean coverage lies in [m_lo, m_hi] (the extrema over the
    box), so a bound on m_lo / m_hi (K_min / K_max) is conservative for the unknown true shift; the
    profile bound is the fixed Hoeffding radius (no data-selected variance)."""
    rng = np.random.default_rng(6)
    means = 5.0 + rng.normal(0.0, 1.0, (400, 4))
    sems = 0.9 * np.abs(1.0 + 0.2 * rng.normal(0.0, 1.0, (400, 4)))
    ref, rsem = np.full(4, 5.0), np.full(4, 0.01)
    e = v5b.v7_estimator_verdict("profile", means, sems, ref, rsem, None, boot_seed=5, n_boot=200)
    sg = e["single"]
    n_e = sg["replicates_eval"]
    me, se_ = means[:n_e], sems[:n_e]
    mu_h, sem_h = np.array(sg["mu_heldout"]), np.array(sg["sem_heldout"])
    for _ in range(300):
        delta = rng.uniform(-sg["z_box"], sg["z_box"], 4)
        m_delta = float(
            np.mean(
                [
                    (np.abs(me[:, b] - (mu_h[b] + delta[b] * sem_h[b])) <= se_[:, b]).mean()
                    for b in range(4)
                ]
            )
        )
        assert sg["m_lo"] - 1e-12 <= m_delta <= sg["m_hi"] + 1e-12
    assert (
        sg["bound_kind"] == "hoeffding"
        and abs(sg["radius"] - v5b.v7_hoeffding_radius(n_e, 0.04)) < 1e-15
    )
    assert (
        abs(sg["ci"][0] - (sg["m_lo"] - sg["radius"])) < 1e-15
        and abs(sg["ci"][1] - (sg["m_hi"] + sg["radius"])) < 1e-15
    )
    assert "variance_max" not in sg


def test_bin_mask_comes_from_the_reference_only(v5b: ModuleType) -> None:
    rng = np.random.default_rng(4)
    means = 5.0 + rng.normal(0.0, 1.0, (16, 3))
    means[:, 2] = -1.0  # the replicates are non-positive in bin 2: irrelevant for the mask
    sems = np.full_like(means, 0.8)
    ref = np.array([5.0, 0.0, 5.0])  # bin 1: reference mean 0 -> not used
    blocks = np.array(
        [[5.0, 0.0, 4.0 + i % 3] for i in range(20)]
    )  # bin 0: constant -> zero spread
    e = v5b.v7_estimator_verdict("profile", means, sems, ref, np.array([0.1, 0.1, 0.1]), blocks)
    assert (
        e["bin_mask"] == [False, False, True]
        and e["bins_used"] == 1
        and "reference" in e["bin_mask_source"]
    )
    means2 = means.copy()
    means2[:, 0] = -9.0  # changing the replicates does not change the mask
    e2 = v5b.v7_estimator_verdict("profile", means2, sems, ref, np.array([0.1, 0.1, 0.1]), blocks)
    assert e2["bin_mask"] == e["bin_mask"]
    blocks3 = np.array(
        [[5.0 + (i % 4) * 0.3, 1.0, 4.0 + (i % 3) * 0.3] for i in range(20)]
    )  # spreads positive, ref bin 1 is 0
    e3 = v5b.v7_estimator_verdict("profile", means, sems, ref, np.array([0.1, 0.1, 0.1]), blocks3)
    assert e3["bin_mask"] == [True, False, True] and e3["covered"] >= 0
    assert len(e3["per_bin_coverage"]) == 3 and e3["per_bin_coverage"][1] is None


def test_profile_empirical_bernstein_fields(v5b: ModuleType) -> None:
    rng = np.random.default_rng(1)
    big = 5.0 + rng.normal(0.0, 1.0, (400, 12))
    sems = np.full_like(big, 1.0)
    p = v5b.v7_estimator_verdict("profile", big, sems, np.full(12, 5.0), np.full(12, 0.1), None)
    assert p["bound_kind"] == "empirical-bernstein" and p["n_pairs"] == 200 and p["bins_used"] == 12
    f = np.array(p["per_pair_covered"]) / 12
    assert abs(p["m"] - f.mean()) < 1e-15 and abs(p["variance"] - f.var(ddof=1)) < 1e-15
    eps = v5b.empirical_bernstein_eps(p["variance"], 200, 0.04)
    assert abs(p["eb_eps"] - eps) < 1e-15 and abs(p["ci"][0] - (p["m"] - eps)) < 1e-15
    assert abs(p["ci"][1] - (p["m"] + eps)) < 1e-15 and len(p["replicate_mean_skewness"]) == 12


def _synthetic_partials(
    v5b: ModuleType, seed: int = 5, kappa: float = 1.0
) -> tuple[list[dict], dict]:
    """Shard and reference partials as the steps write them, from i.i.d. normal batch sampling."""
    rng = np.random.default_rng(seed)
    nb, per = v5b.V7_BATCHES_PER_REP, v5b.V7_REPS_PER_SHARD
    sig = {
        "sec_p": np.linspace(0.5, 2.0, 12),
        "nuclear_local": np.linspace(1.0, 3.0, 12),
        "escaped_neutral": np.array([2.5]),
    }
    mu = {
        "sec_p": 10.0 + np.arange(12.0),
        "nuclear_local": 20.0 + np.arange(12.0),
        "escaped_neutral": np.array([80.0]),
    }
    shards = []
    for k in range(v5b.V7_SHARDS):
        est = {}
        for nm in sig:
            x = mu[nm] + sig[nm] * np.sqrt(nb) * rng.standard_normal((per, nb, sig[nm].size))
            est[nm] = {
                "mean": x.mean(axis=1).tolist(),
                "sem": (kappa * x.std(axis=1, ddof=1) / math.sqrt(nb)).tolist(),
            }
        shards.append({"seed": v5b.seed_of("v7-rep", k), "shard": k, "valid": True,
                       "replicates": per,  # fmt: skip
                       "n": 9_000_000, "estimators": est})  # fmt: skip
    ref_est = {}
    for nm in sig:
        x = mu[nm] + sig[nm] * np.sqrt(nb) * 0.1 * rng.standard_normal((20, sig[nm].size))
        x = mu[nm] + (x - x.mean(axis=0)) / x.std(axis=0, ddof=1) * sig[nm] * np.sqrt(20) * 0.1
        ref_est[nm] = {"mean": mu[nm].tolist(), "sem": (0.1 * sig[nm]).tolist(),
                       "blocks": x.tolist()}  # fmt: skip
    ref = {
        "seed": v5b.seed_of("v7-rep", v5b.V7_SHARDS),
        "valid": True,
        "n": 1_000_000,
        "estimators": ref_est,
    }
    return shards, ref


def test_combine_document_and_fail_closed(v5b: ModuleType) -> None:
    shards, ref = _synthetic_partials(v5b)
    assert len(shards) == 8
    doc = v5b.v7_rep_combine(shards, ref)
    assert doc["replicates"] == 7200 and set(doc["estimators"]) == {
        "sec_p",
        "nuclear_local",
        "escaped_neutral",
    }
    rule = doc["rule"]
    c0, low, high = v5b.v7_region()
    assert (rule["c0"], rule["low"], rule["high"], rule["pairs"]) == (c0, low, high, 3600)
    assert rule["alpha_tost"] == 0.04 and rule["bins_from"] == "v7-rep-ref only"
    for name, kind, bins in (
        ("sec_p", "profile", 12),
        ("nuclear_local", "profile", 12),
        ("escaped_neutral", "scalar", 1),
    ):
        e = doc["estimators"][name]
        assert e["kind"] == kind and e["bins_used"] == bins and e["n_pairs"] == 3600
        assert len(e["per_pair_covered"]) == 3600 and all(
            isinstance(c, int) for c in e["per_pair_covered"]
        )
        for key in (
            *(
                "m",
                "c0",
                "low",
                "high",
                "bound_kind",
                "eb_eps",
                "variance",
                "ci",
                "reasons",
                "bin_mask",
            ),
            *(
                "per_bin_coverage",
                "replicate_mean_skewness",
                "legacy_point_gate_pass",
                "m_ref1e6",
                "pass",
            ),
            *("bin_mask_source", "alpha_tost", "margin", "pair_df"),
        ):
            assert key in e
        assert not {"z_box", "z_boot", "m_lo", "m_hi", "mu_heldout", "hoeffding_eps"} & set(e)
        assert e["bound_kind"] == ("clopper-pearson" if kind == "scalar" else "empirical-bernstein")
        assert (
            e["bin_mask"] == [True] * bins and abs(e["m"] - c0) < 0.03
        )  # correct SEMs: coverage near c0
        sg = e["single"]
        assert sg["pass"] and e["pass_paired"] and e["pass_single"] and e["pass"], (
            name,
            e["reasons"],
        )
        assert abs(sg["c1"] - 0.6701) < 2e-4 and sg["bound_kind"] == (
            "clopper-pearson" if kind == "scalar" else "hoeffding"
        )
        assert "nu" not in sg and abs(sg["low"] - (sg["c1"] - 0.03)) < 1e-15
        assert (
            sg["z_box"] == max(sg["z_t"], sg["z_boot"])
            and sg["boot_seed"] == ref["seed"]
            and sg["boot_n"] == 5000
        )
        assert (
            sg["replicates_eval"] == 5400
            and sg["replicates_heldout"] == 1800
            and abs(sg["m"] - sg["c1"]) < 0.05
        )
        assert (
            sg["m_lo"] <= sg["m"] <= sg["m_hi"] and len(sg["per_bin_box_coverage_min_max"]) == bins
        )
        assert 0.5 < e["m_ref1e6"] < 0.85
    assert doc["pass"] is True
    # a clearly overconfident estimator (SEM x 0.5: coverage ~ 0.38) fails the row
    bad_shards, bad_ref = _synthetic_partials(v5b, kappa=0.5)
    bad = v5b.v7_rep_combine(bad_shards, bad_ref)
    assert not bad["pass"] and not any(e["pass"] for e in bad["estimators"].values())
    # fail closed: missing estimator, wrong seed, invalid or missing shard, wrong shape
    import copy

    miss = copy.deepcopy(shards)
    del miss[1]["estimators"]["escaped_neutral"]
    with pytest.raises(SystemExit, match="escaped_neutral"):
        v5b.v7_rep_combine(miss, ref)
    with pytest.raises(SystemExit):
        v5b.v7_rep_combine([{**shards[0], "seed": 1}, *shards[1:]], ref)
    with pytest.raises(SystemExit):
        v5b.v7_rep_combine([{**shards[0], "valid": False}, *shards[1:]], ref)
    with pytest.raises(SystemExit):
        v5b.v7_rep_combine(shards[:7], ref)
    uneven = copy.deepcopy(shards)
    uneven[7]["replicates"] = 450
    with pytest.raises(SystemExit, match="equal replicate counts"):
        v5b.v7_rep_combine(uneven, ref)
    with pytest.raises(SystemExit):
        v5b.v7_rep_combine(shards, {**ref, "valid": False})
    short = copy.deepcopy(ref)
    short["estimators"]["sec_p"]["mean"] = short["estimators"]["sec_p"]["mean"][:5]
    with pytest.raises(SystemExit, match="inconsistent"):
        v5b.v7_rep_combine(shards, short)
    noblocks = copy.deepcopy(ref)
    del noblocks["estimators"]["sec_p"]["blocks"]
    with pytest.raises(SystemExit, match="sec_p"):
        v5b.v7_rep_combine(shards, noblocks)
    # the bin mask follows the reference partial alone: zero its mean in one bin
    masked = copy.deepcopy(ref)
    masked["estimators"]["sec_p"]["mean"][3] = 0.0
    md = v5b.v7_rep_combine(shards, masked)["estimators"]["sec_p"]
    assert md["bins_used"] == 11 and md["bin_mask"][3] is False


def test_replicate_grouping_and_estimator_verdict(v5b: ModuleType) -> None:
    b = np.arange(80.0).reshape(40, 2)  # 2 replicates of 20 batches, 2 bins
    m, s = v5b.v7_replicate_stats(b)
    assert m.shape == s.shape == (2, 2)
    assert np.allclose(m[0], b[:20].mean(axis=0)) and np.allclose(m[1], b[20:].mean(axis=0))
    assert np.allclose(s[1], b[20:].std(axis=0, ddof=1) / math.sqrt(20))
    m1, _ = v5b.v7_replicate_stats(np.arange(40.0))  # scalar estimator: [B] -> [R, 1]
    assert m1.shape == (2, 1)


# -- Monte Carlo calibration of the full procedure (marker ``calibration``: not in CI) --------
# The procedure is ``v7_estimator_verdict`` itself on simulated replicate means (R = 7200 -> 3600 pairs); no
# reference is simulated: the bin mask comes from the (independent) reference run and is all-true here.
def _kappa_for(v5b: ModuleType, coverage: float, df: int = 38) -> float:
    """kappa with P(|t_38| <= kappa) = coverage (paired one-sigma rule scaled by kappa, Gaussian batches)."""
    lo, hi = 0.1, 5.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if v5b.student_abs_prob(df, mid) < coverage else (lo, mid)
    return 0.5 * (lo + hi)


N_BOOT_CAL = 500  # bootstrap resamples inside the Monte Carlo (5000 in the implementation default)


def _verdicts(
    v5b: ModuleType, kind: str, mean: np.ndarray, sem0: np.ndarray, kappas: tuple[float, ...],
    blocks: np.ndarray | None = None,
) -> tuple[np.ndarray, list[bool]]:  # fmt: skip
    """``v7_estimator_verdict`` pass flags ``[row, paired gate, single gate]`` at each kappa (reported
    half width kappa * sem) ``[3, len(kappas)]`` and the bin mask of the first one. Reference: the simulated block means ``blocks`` (resolvability criterion
    applied) or, for the Gaussian cases, a dense dummy (ref 1, SEM 0.01: r_b = 0.1, all bins used)."""
    if blocks is not None:
        ref, rsem = blocks.mean(axis=0), blocks.std(axis=0, ddof=1) / math.sqrt(blocks.shape[0])
    else:
        ref, rsem = np.ones(mean.shape[1]), np.full(mean.shape[1], 0.01)
    res, zb = [], None
    for k in kappas:  # the bootstrap quantile does not depend on kappa: computed once per sample
        r = v5b.v7_estimator_verdict(kind, mean, k * sem0, ref, rsem, blocks, boot_seed=1234, n_boot=N_BOOT_CAL,
                                     z_boot_override=zb)  # fmt: skip
        zb = r["single"]["z_boot"]
        res.append(r)
    flags = np.array([[bool(r[k]) for r in res] for k in ("pass", "pass_paired", "pass_single")])
    return flags, list(res[0]["bin_mask"])


def _gauss_passes(v5b: ModuleType, kappas: tuple[float, ...], bins: int, trials: int,
                  rng: np.random.Generator, rho: float = 0.0, reps: int = 7200) -> np.ndarray:  # fmt: skip
    """Pass indicators ``[3 (row, paired, single), len(kappas), trials]``. Replicate means X_b = sigma_b (sqrt(rho) g + sqrt(1 -
    rho) e_b) with a replicate-level common factor g (bin correlation ``rho``, as for bins sharing
    histories), batch SEMs s = sigma_b sqrt(chi2_19 / 19 / 20-batch scaling); true mean 0 (shifted to 100)."""
    sigma = np.linspace(0.5, 2.0, bins) if bins > 1 else np.array([1.3])
    out = np.zeros((3, len(kappas), trials), dtype=bool)
    kind = "profile" if bins > 1 else "scalar"
    for i in range(trials):
        g = rng.standard_normal((reps, 1))
        x = sigma * (math.sqrt(rho) * g + math.sqrt(1.0 - rho) * rng.standard_normal((reps, bins)))
        s0 = sigma * np.sqrt(rng.chisquare(19, (reps, bins)) / 19.0)
        out[:, :, i] = _verdicts(v5b, kind, 100.0 + x, s0, kappas)[0]
    return out


def _rate_report(v5b: ModuleType, passes: np.ndarray) -> list[dict[str, float]]:
    """Rate and 95 % Clopper-Pearson upper bound per row of a ``[kappas, trials]`` pass array."""
    n = passes.shape[1]
    return [
        {"rate": float(r.mean()), "cp95_upper": v5b.cp_upper(int(r.sum()), n, 0.05), "n": n}
        for r in passes
    ]


GAUSS_CASES = (
    ("profile12_rho0.0", 12, 600, 0.0),
    ("profile12_rho0.8", 12, 600, 0.8),
    ("scalar", 1, 3000, 0.0),
)


@pytest.fixture(scope="module")
def gauss_rates(v5b: ModuleType) -> dict[str, np.ndarray]:
    """Pass indicators ``[3, trials]`` (true paired coverage low / c0 / high) per Gaussian case."""
    rng = np.random.default_rng(20471004)
    c0, low, high = v5b.v7_region()
    c1, s_low, s_high = v5b.v7_single_region()
    # kappas: paired true coverage at low / high, single-interval true coverage P(|t_19| <= kappa) at
    # c1 -/+ 0.03, and the nominal; the ROW (both gates) must reject at all four margins
    ks = (
        _kappa_for(v5b, low),
        _kappa_for(v5b, high),
        _kappa_for(v5b, s_low, 19),
        _kappa_for(v5b, s_high, 19),
        1.0,
    )
    passes = {lab: _gauss_passes(v5b, ks, b, n, rng, rho=rho) for lab, b, n, rho in GAUSS_CASES}
    print("v7 gaussian kappas [paired low, paired high, single low, single high, 1]", ks,
          "single true coverage at those kappas", [v5b.student_abs_prob(19, k) for k in ks],
          "paired true coverage", [v5b.student_abs_prob(38, k) for k in ks])  # fmt: skip
    for gate, gi in (("ROW", 0), ("PAIRED gate", 1), ("SINGLE gate", 2)):
        print(f"v7 gaussian {gate} rates at [paired low, paired high, single low, single high, nominal] with CP95 upper bounds",
              {k: _rate_report(v5b, v[gi]) for k, v in passes.items()}, "regions", (low, c0, high), (s_low, c1, s_high))  # fmt: skip
    return passes


@pytest.mark.calibration
@pytest.mark.parametrize("label", [c[0] for c in GAUSS_CASES])
def test_v7_pipeline_gaussian_false_acceptance(
    v5b: ModuleType, gauss_rates: dict, label: str
) -> None:
    """Gaussian pipeline (3600 pairs, 5400 + 1800 replicates), asserted separately: (a) the paired gate's
    false-acceptance rate, 95 % Clopper-Pearson upper bound <= 0.05 at the true paired coverage c0 -/+
    0.03 (profiles; the scalar paired gate is gated by its exact value, see below); (b) the single-interval gate's, <= 0.05 at the true single-interval coverage c1 -/+ 0.03
    (kappa calibrated on sem_j alone against the true mean); (c) the row (both gates) at all four
    margins. Scalar, paired margins: the paired-gate rate must not exceed the EXACT paired false
    acceptance (``test_scalar_exact_false_acceptance_at_both_margins``) by more than 3 standard errors."""
    row, paired, single = (_rate_report(v5b, gauss_rates[label][g]) for g in (0, 1, 2))
    assert all(r["cp95_upper"] <= 0.05 for r in row[:4]), ("row", label, row)
    assert all(r["cp95_upper"] <= 0.05 for r in (single[2], single[3])), ("single", label, single)
    if label != "scalar":
        assert all(r["cp95_upper"] <= 0.05 for r in (paired[0], paired[1])), (
            "paired",
            label,
            paired,
        )
        return
    # scalar paired gate alone: the exact value (<= 0.04, about 0.0387) is the gate; the simulated rate is
    # a consistency check, one-sided at 4 standard errors (it is compared at two margins of one stream)
    kl, kc = _scalar_acceptance(v5b)
    _, low, high = v5b.v7_region()
    for r, p in ((paired[0], low), (paired[1], high)):
        exact = _binom_interval_prob(3600, p, kl, kc)
        se = math.sqrt(exact * (1.0 - exact) / r["n"])
        assert r["rate"] <= exact + 4.0 * se, (r, exact, se)


@pytest.mark.calibration
def test_v7_pipeline_gaussian_power_targets(gauss_rates: dict) -> None:
    """Power at the nominal coverage (kappa = 1): profiles >= 0.95, scalar >= 0.90, joint row
    (independent draws of the two profiles and the scalar) >= 0.85."""
    power = {k: float(v[0][4].mean()) for k, v in gauss_rates.items()}
    single_power = {k: float(v[2][4].mean()) for k, v in gauss_rates.items()}
    n = 600
    joint = float((gauss_rates["profile12_rho0.0"][0][4][:n] & gauss_rates["profile12_rho0.8"][0][4][:n]
                   & gauss_rates["scalar"][0][4][:n]).mean())  # fmt: skip
    print(
        "v7 gaussian ROW power at nominal",
        power,
        "single gate alone",
        single_power,
        "joint row (independent draws)",
        joint,
    )
    assert power["profile12_rho0.0"] >= 0.95 and power["profile12_rho0.8"] >= 0.95, power
    assert power["scalar"] >= 0.90 and joint >= 0.85, (power, joint)


W_SHARED = 0.5  # shape of the shared per-block Gamma factor of the correlated sparse bins (mean 1)
ZERO_BINS = (3, 8)  # bins with 90 % exact zeros in the zero-inflated case
ZERO_KEEP = 0.1


def _gamma_blocks(shape: float, size: tuple[int, ...], rng: np.random.Generator, correlated: bool,
                  zero_bins: tuple[int, ...] = ()) -> np.ndarray:  # fmt: skip
    """Block values ``[..., bins]`` with mean ``shape``. Independent: Gamma(shape) per bin. Correlated
    sparse bins as from shared histories: the same Gamma draws times ONE factor W_b ~ Gamma(0.5) / 0.5
    (mean 1) per block shared by all bins, so that a block with few histories is sparse in every bin.
    Zero-inflated: in ``zero_bins`` a fraction 0.9 of the blocks are exactly zero (true mean 0.1 shape)."""
    g = rng.standard_gamma(shape, size, dtype=np.float32)
    if correlated:
        g = g * (rng.standard_gamma(W_SHARED, (*size[:-1], 1), dtype=np.float32) / W_SHARED)
    if zero_bins:
        keep = rng.random((*size[:-1], len(zero_bins))) < ZERO_KEEP
        g[..., list(zero_bins)] *= keep
    return g


def _ratios(
    shape: float, correlated: bool, zero: bool, rng: np.random.Generator, n: int = 300_000
) -> tuple[np.ndarray, np.ndarray]:
    """Paired ``|X_j - X_k| / sqrt(s_j^2 + s_k^2)`` of one bin over ``n`` pairs of independent replicates
    (20 Gamma blocks each) and the single-interval ``|X_j - true| / s_j`` (true mean ``shape`` or
    ``0.1 shape`` for a zero-inflated bin); a vanishing SEM (either replicate for the pair, the replicate
    for the single interval) gives infinity: a miss (undefined uncertainty never certifies coverage)."""
    g = _gamma_blocks(shape, (2, n, 20, 1), rng, correlated, (0,) if zero else ()).astype(float)[
        ..., 0
    ]
    mean, sem = g.mean(axis=2), g.std(axis=2, ddof=1) / math.sqrt(20)
    true = shape * (ZERO_KEEP if zero else 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        paired = np.abs(mean[0] - mean[1]) / np.sqrt(sem[0] ** 2 + sem[1] ** 2)
        single = np.abs(mean[0] - true) / sem[0]
    paired = np.where((sem[0] > 0) & (sem[1] > 0), paired, np.inf)
    single = np.where(sem[0] > 0, single, np.inf)
    return paired, single


def _sim_reference(
    shape: float, correlated: bool, zero: bool, rng: np.random.Generator, bins: int = 12
) -> np.ndarray:
    """Simulated reference run: 20 blocks of n_ref / 20 = 5e4 histories, each the mean of 100 replicate
    blocks (500 histories) of the same bin distributions (shared factor and zero inflation included);
    returns the block means ``[20, bins]``."""
    scale = np.linspace(0.5, 2.0, bins)
    g = _gamma_blocks(shape, (20, 100, bins), rng, correlated, ZERO_BINS if zero else ()) * scale
    return g.astype(float).mean(axis=1)


def _expected_mask(v5b: ModuleType, shape: float, correlated: bool, zero: bool, rng: np.random.Generator,
                   draws: int = 40) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """Expected bin mask and mean r_b of the resolvability criterion over simulated reference runs."""
    rel = np.zeros(12)
    for _ in range(draws):
        blk = _sim_reference(shape, correlated, zero, rng)
        rel += (
            blk.std(axis=0, ddof=1)
            / math.sqrt(20)
            / blk.mean(axis=0)
            * math.sqrt(1_000_000 / 10_000)
        )
    rel /= draws
    return rel < v5b.V7_MAX_REL_SE, rel


def _skew_kappas(
    v5b: ModuleType, shape: float, correlated: bool, zero: bool, rng: np.random.Generator,
    mask: np.ndarray | None = None,
) -> dict[str, float]:  # fmt: skip
    """kappa giving the TRUE bin-averaged paired coverage low and high (large-sample simulation of the
    paired statistic with the miss rule, bisection with common random numbers over the bin marginals: 10
    Gamma bins and, if ``zero``, 2 zero-inflated bins), the paired true coverage at kappa = 1, and the
    single-interval true coverage (against the true mean) at the same three kappas (report only)."""
    _, low, high = v5b.v7_region()
    mask = np.ones(12, dtype=bool) if mask is None else mask
    zb = np.zeros(12, dtype=bool)
    if zero:
        zb[list(ZERO_BINS)] = True
    n_g, n_z, n_all = int((mask & ~zb).sum()), int((mask & zb).sum()), int(mask.sum())
    marg = []
    if n_g:
        pa, si = _ratios(shape, correlated, False, rng)
        marg.append((n_g / n_all, np.sort(pa), np.sort(si)))
    if n_z:
        pz, sz = _ratios(shape, correlated, True, rng)
        marg.append((n_z / n_all, np.sort(pz), np.sort(sz)))

    def cov(k: float, which: int) -> float:
        return float(
            sum(m[0] * np.searchsorted(m[which], k, side="right") / m[which].size for m in marg)
        )

    def solve(target: float, which: int) -> float:
        lo, hi = 0.05, 30.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if cov(mid, which) < target else (lo, mid)
        return 0.5 * (lo + hi)

    c1, s_low, s_high = v5b.v7_single_region()
    k_lo, k_hi = solve(low, 1), solve(high, 1)  # paired margins
    ks_lo, ks_hi = solve(s_low, 2), solve(s_high, 2)  # single-interval margins c1 -/+ 0.03
    assert abs(cov(k_lo, 1) - low) < 2e-3 and abs(cov(k_hi, 1) - high) < 2e-3
    assert abs(cov(ks_lo, 2) - s_low) < 2e-3 and abs(cov(ks_hi, 2) - s_high) < 2e-3
    return {"kappa_low": k_lo, "kappa_high": k_hi, "kappa_single_low": ks_lo, "kappa_single_high": ks_hi,
            "paired_at_kappa1": cov(1.0, 1), "single_at_low": cov(k_lo, 2), "single_at_kappa1": cov(1.0, 2),
            "single_at_high": cov(k_hi, 2), "paired_at_single_low": cov(ks_lo, 1), "paired_at_single_high": cov(ks_hi, 1)}  # fmt: skip


def _skew_passes(v5b: ModuleType, shape: float, correlated: bool, zero: bool, kappas: tuple[float, ...],
                 trials: int, rng: np.random.Generator, bins: int = 12, reps: int = 7200) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """Pass indicators ``[3 (row, paired, single), len(kappas), trials]`` and the bin masks ``[trials, bins]`` on Gamma-block tallies (20 blocks per replicate and bin,
    replicate mean and SEM from the blocks, per-bin scales 0.5..2): the implemented paired procedure with
    a simulated reference run (20 blocks of 5e4 histories) for the bin mask (resolvability criterion)."""
    scale = np.linspace(0.5, 2.0, bins)
    out = np.zeros((3, len(kappas), trials), dtype=bool)
    masks: list[list[bool]] = []
    for i in range(trials):
        g = (
            _gamma_blocks(shape, (reps, 20, bins), rng, correlated, ZERO_BINS if zero else ())
            * scale
        )
        mean, sem = (
            g.mean(axis=1).astype(float),
            (g.std(axis=1, ddof=1) / math.sqrt(20)).astype(float),
        )
        blocks = _sim_reference(shape, correlated, zero, rng, bins)
        out[:, :, i], mask = _verdicts(v5b, "profile", mean, sem, kappas, blocks)
        masks.append(mask)
    return out, np.array(masks)


SKEW_CASES = [(sh, corr, False) for sh in (2.0, 0.5) for corr in (False, True)] + [
    (2.0, False, True),
    (0.5, True, True),
]


@pytest.fixture(scope="module", params=SKEW_CASES,
                ids=[f"gamma{sh}-{'correlated' if c else 'independent'}{'-zeroinflated' if z else ''}" for sh, c, z in SKEW_CASES])  # fmt: skip
def skew_result(request: pytest.FixtureRequest, v5b: ModuleType) -> dict:
    """Gamma shape 2 and 0.5, independent / correlated sparse bins and zero-inflated cases (two of the 12
    bins have 90 % exact zeros): kappa scales the paired half width so that the TRUE bin-averaged paired
    coverage is c0 -/+ 0.03; pass rates with CP bounds there, power at kappa = 1 and the true coverage there."""
    shape, corr, zero = request.param
    rng = np.random.default_rng(
        20471005 + int(2 * shape) + (7 if corr else 0) + (13 if zero else 0)
    )
    mask, rel = _expected_mask(v5b, shape, corr, zero, rng)
    if (
        mask.sum() < v5b.V7_MIN_PROFILE_BINS
    ):  # fewer than 10 resolvable bins: the profile fails closed
        passes, masks = _skew_passes(v5b, shape, corr, zero, (1.0,), 100, rng)
        passes = passes[0]
        out = {"shape": shape, "correlated": corr, "zero_inflated": zero, "expected_bins_used": int(mask.sum()),
               "excluded_bins": [int(b) for b in np.flatnonzero(~mask)], "r_b_mean": [round(float(x), 3) for x in rel],
               "fail_closed": True, "passes_at_kappa1": int(passes.sum()),
               "trials_with_fewer_than_10_bins": float(np.mean(masks.sum(axis=1) < v5b.V7_MIN_PROFILE_BINS))}  # fmt: skip
        print("v7 skewed-tally calibration", out)
        return out
    kap = _skew_kappas(v5b, shape, corr, zero, rng, mask)
    ks = (
        kap["kappa_low"],
        kap["kappa_high"],
        kap["kappa_single_low"],
        kap["kappa_single_high"],
        1.0,
    )
    passes, masks = _skew_passes(v5b, shape, corr, zero, ks, 300, rng)
    row, paired, single = (_rate_report(v5b, passes[g]) for g in (0, 1, 2))
    out = {"shape": shape, "correlated": corr, "zero_inflated": zero, "expected_bins_used": int(mask.sum()),
           "excluded_bins": [int(b) for b in np.flatnonzero(~mask)], "r_b_mean": [round(float(x), 3) for x in rel],
           "trials_with_expected_mask": float(np.mean((masks == mask).all(axis=1))), **kap,
           "PAIRED_FA_low": paired[0], "PAIRED_FA_high": paired[1], "SINGLE_FA_low": single[2],
           "SINGLE_FA_high": single[3], "ROW_FA_paired_low": row[0], "ROW_FA_paired_high": row[1],
           "ROW_FA_single_low": row[2], "ROW_FA_single_high": row[3], "ROW_power_kappa1": row[4],
           "PAIRED_power_kappa1": paired[4], "SINGLE_power_kappa1": single[4]}  # fmt: skip
    print("v7 skewed-tally calibration", out)
    return out


@pytest.mark.calibration
def test_v7_pipeline_skewed_tally_false_acceptance(skew_result: dict) -> None:
    """The 95 % Clopper-Pearson upper bound of the false-acceptance rate is <= 0.05 at the true paired
    coverage c0 - 0.03 and c0 + 0.03 for Gamma-block tallies: independent, correlated (shared-history)
    and zero-inflated bins."""
    if skew_result.get(
        "fail_closed"
    ):  # no resolvable profile: the gate never passes (FA = 0 by design)
        assert (
            skew_result["passes_at_kappa1"] == 0
            and skew_result["trials_with_fewer_than_10_bins"] > 0.9
        ), skew_result
        return
    for key in (
        "PAIRED_FA_low", "PAIRED_FA_high", "SINGLE_FA_low", "SINGLE_FA_high",
        "ROW_FA_paired_low", "ROW_FA_paired_high", "ROW_FA_single_low", "ROW_FA_single_high",
    ):  # fmt: skip
        assert skew_result[key]["cp95_upper"] <= 0.05, (key, skew_result)


def test_calibration_marker_registered_and_excluded_from_ci() -> None:
    """The Monte Carlo calibration is out of the lightweight CI selection and registered as a marker."""
    ci = (REPO / ".github" / "workflows" / "ci.yml").read_text()
    assert '-m "not cuda and not host and not calibration"' in ci
    assert '"calibration:' in (REPO / "pyproject.toml").read_text()


# -- full-scale block configurations and the fixed-point accumulator guard (C32) ------------
def _guard_env(v5b: ModuleType, fn):  # type: ignore[no-untyped-def]
    """Run ``fn`` (builds an effective configuration: needs the nuclear table of the cache); skip
    without it unless ``IONMC_REQUIRE_DATA=1`` (then a failure)."""
    import os

    try:
        return fn()
    except Exception as exc:  # noqa: BLE001
        if "table" in str(exc).lower() or type(exc).__name__.startswith("NuclearTable"):
            if os.environ.get("IONMC_REQUIRE_DATA") == "1":
                raise
            pytest.skip(f"no built nuclear table in the cache: {exc!r}")
        raise


def test_full_scale_block_configurations_pass_the_accumulator_guard(v5b: ModuleType) -> None:
    """The nominal two-batch configuration of a 9e6-history shard exceeds the fixed-point voxel
    accumulator guard of ``validate`` (the host failure of C32); ``batch_estimates`` now validates
    through ``block_effective`` (smallest batch count 2, 4, 8, ... that satisfies the guard), and the
    full-scale shard, reference and scan configurations all pass through the same entry path."""
    from ionmc.config import MAX_QUANTA, QUANTUM_MEV
    from ionmc.errors import UnsupportedCombinationError
    from ionmc.simulation import Simulation

    geo, grid = v5b.coarse_depth(v5b.V7_BINS)
    shard = v5b.v7_config(
        "warp-cpu", "float64", v5b.V7_SHARD_N, v5b.V7_SHARD_BATCHES, 1, grid=grid, geo=geo
    )
    assert (shard.run.n_histories, shard.run.n_batches) == (9_000_000, 18000)
    # pre-D2 anchor: the Warp backends have no elastic channel until C4
    shard = v5b.replace(shard, physics=v5b.replace(shard.physics, elastic=False))
    nominal2 = v5b.replace(shard, run=v5b.replace(shard.run, n_batches=2))
    with pytest.raises(UnsupportedCombinationError, match="fixed-point voxel accumulator"):
        _guard_env(v5b, lambda: Simulation(nominal2).effective)  # the old code path
    eff_nominal = _guard_env(
        v5b, lambda: Simulation(shard).effective
    )  # the shard's own batch structure
    e_cap = max(
        eff_nominal.nuclear.history_energy_bound_mev if eff_nominal.nuclear is not None else 0.0,
        0.0,
    )
    assert abs(e_cap - 3154.43) < 0.01
    eff = _guard_env(v5b, lambda: v5b.block_effective(shard))
    b = eff.requested.run.n_batches
    assert (
        b == 8
        and 9_000_000 // b * e_cap / QUANTUM_MEV
        < MAX_QUANTA
        <= 9_000_000 // (b // 2) * e_cap / QUANTUM_MEV
    )
    for n, nb in (
        (v5b.V7_REF_N, v5b.V7_REF_BATCHES),
        (1_000_000, v5b.V7_BATCHES),
    ):  # reference, scan
        cfg = v5b.v7_config("warp-cpu", "float64", n, nb, 1, grid=grid, geo=geo)
        cfg = v5b.replace(cfg, physics=v5b.replace(cfg.physics, elastic=False))  # pre-D2 anchor
        assert (
            _guard_env(v5b, lambda cfg=cfg: v5b.block_effective(cfg)).requested.run.n_batches == 2
        )
        assert _guard_env(v5b, lambda cfg=cfg: Simulation(cfg).effective) is not None


def test_guard_table_of_every_step_configuration(v5b: ModuleType) -> None:
    """Histories per batch x the history energy bound (3154.43 MeV, nuclear) against the capacity of the
    fixed-point voxel accumulator (4.29497e9 MeV) for every nominal configuration of the V3-005B steps:
    each must be below, or run through ``block_effective``."""
    from ionmc.config import MAX_QUANTA, QUANTUM_MEV

    cap = MAX_QUANTA * QUANTUM_MEV
    e_cap = 3154.43
    steps = {  # name -> (histories, batches) of the configuration validated at full scale
        "lv5b-throughput python": (v5b.THROUGHPUT_PY_N, 2),
        "lv5b-throughput warp / v2b / v5 geometry": (v5b.THROUGHPUT_WARP_N, 2),
        "v8-lv python vs warp-cpu (256, dense 96)": (v5b.V8_K, 2),
        "r1-nuc": (v5b.R1_NUC_N, 2),
        "v5-ionmc 150/200 on/off": (v5b.V5_N, v5b.V5_BATCHES),
        "v2b shard (per variant)": (v5b.V2B_N // v5b.V2B_SHARDS, v5b.V2B_BATCHES),
        "v7-scan N = 1e6 (block_effective, 2 batches)": (1_000_000, 2),
        "v7-scan nominal": (1_000_000, v5b.V7_BATCHES),
        "v7-shift": (v5b.V7_SHIFT_N, v5b.V7_SHIFT_BATCHES),
        "v7-rep-ref 1e6 (block_effective, 2 batches)": (v5b.V7_REF_N, 2),
        "v7-rep-s{k} 9e6 (block_effective, 8 batches)": (v5b.V7_SHARD_N, 8),
        "v7-rep-s{k} 9e6 nominal 18000": (v5b.V7_SHARD_N, v5b.V7_SHARD_BATCHES),
        "hr5 python sample": (v5b.HR5_PYTHON_N, v5b.HR5_PYTHON_BATCHES),
        "hr5 warp samples 1e6": (v5b.HR5_WARP_N, v5b.HR5_WARP_BATCHES),
        "hr5 v7-f32 1e6": (v5b.HR5_F32_N, v5b.V7_BATCHES),
    }
    table = {k: (n // b, n // b * e_cap) for k, (n, b) in steps.items()}
    print(f"v7 accumulator guard table (histories per batch, MeV; capacity {cap:.5g} MeV)", table)
    for k, (_per, mev) in table.items():
        assert mev < cap, (k, mev)
    assert 9_000_000 // 2 * e_cap >= cap  # the old nominal two-batch shard configuration fails
    assert abs(cap - 4.29497e9) < 1e4


def test_block_estimates_independent_of_the_effective_batch_count(v5b: ModuleType) -> None:
    """``batch_estimates`` over 4 blocks of a 2000-history run is bitwise identical for effective
    configurations of 2 (the old nominal two-batch configuration), 4 and the default (smallest valid)
    batch count: the trajectories do not depend on the batch structure and the accumulators are
    fixed-point integers."""
    geo, grid = v5b.coarse_depth(v5b.V7_BINS)
    cfg = v5b.v7_config("warp-cpu", "float64", 2000, 4, 20461004, grid=grid, geo=geo)
    # pre-D2 anchor: the Warp backends have no elastic channel until C4
    cfg = v5b.replace(cfg, physics=v5b.replace(cfg.physics, elastic=False))
    ref = _guard_env(v5b, lambda: v5b.batch_estimates(cfg, 2))
    for eb in (4, None):
        got = v5b.batch_estimates(cfg, eb)
        for key in ("sec_p", "nuc_local_dose", "idd", "nuclear_local", "escaped_neutral"):
            assert got[key].shape == ref[key].shape and np.array_equal(got[key], ref[key]), (
                key,
                eb,
            )
        assert got["counters_sum"] == ref["counters_sum"] == 0
