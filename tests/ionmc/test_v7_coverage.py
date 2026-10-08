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
        [[5.0 + i % 4, 1.0, 4.0 + i % 3] for i in range(20)]
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
        assert e["pass"], (name, e["reasons"], e["ci"])
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
def _kappa_for(v5b: ModuleType, coverage: float) -> float:
    """kappa with P(|t_38| <= kappa) = coverage (paired one-sigma rule scaled by kappa, Gaussian batches)."""
    lo, hi = 0.1, 5.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if v5b.student_abs_prob(38, mid) < coverage else (lo, mid)
    return 0.5 * (lo + hi)


def _verdicts(
    v5b: ModuleType, kind: str, mean: np.ndarray, sem0: np.ndarray, kappas: tuple[float, ...]
) -> list[bool]:
    """``v7_estimator_verdict`` pass flags at each kappa (half width kappa * sqrt(s_j^2 + s_k^2))."""
    ref, rsem = np.ones(mean.shape[1]), np.full(mean.shape[1], 0.1)
    return [
        bool(v5b.v7_estimator_verdict(kind, mean, k * sem0, ref, rsem, None)["pass"])
        for k in kappas
    ]


def _gauss_passes(v5b: ModuleType, kappas: tuple[float, ...], bins: int, trials: int,
                  rng: np.random.Generator, rho: float = 0.0, reps: int = 7200) -> np.ndarray:  # fmt: skip
    """Pass indicators ``[len(kappas), trials]``. Replicate means X_b = sigma_b (sqrt(rho) g + sqrt(1 -
    rho) e_b) with a replicate-level common factor g (bin correlation ``rho``, as for bins sharing
    histories), batch SEMs s = sigma_b sqrt(chi2_19 / 19 / 20-batch scaling); true mean 0 (shifted to 100)."""
    sigma = np.linspace(0.5, 2.0, bins) if bins > 1 else np.array([1.3])
    out = np.zeros((len(kappas), trials), dtype=bool)
    kind = "profile" if bins > 1 else "scalar"
    for i in range(trials):
        g = rng.standard_normal((reps, 1))
        x = sigma * (math.sqrt(rho) * g + math.sqrt(1.0 - rho) * rng.standard_normal((reps, bins)))
        s0 = sigma * np.sqrt(rng.chisquare(19, (reps, bins)) / 19.0)
        out[:, i] = _verdicts(v5b, kind, 100.0 + x, s0, kappas)
    return out


def _rate_report(v5b: ModuleType, passes: np.ndarray) -> list[dict[str, float]]:
    """Rate and 95 % Clopper-Pearson upper bound per row of a ``[kappas, trials]`` pass array."""
    n = passes.shape[1]
    return [
        {"rate": float(r.mean()), "cp95_upper": v5b.cp_upper(int(r.sum()), n, 0.05), "n": n}
        for r in passes
    ]


GAUSS_CASES = (
    ("profile12_rho0.0", 12, 1000, 0.0),
    ("profile12_rho0.8", 12, 1000, 0.8),
    ("scalar", 1, 10000, 0.0),
)


@pytest.fixture(scope="module")
def gauss_rates(v5b: ModuleType) -> dict[str, np.ndarray]:
    """Pass indicators ``[3, trials]`` (true paired coverage low / c0 / high) per Gaussian case."""
    rng = np.random.default_rng(20471004)
    c0, low, high = v5b.v7_region()
    ks = (_kappa_for(v5b, low), 1.0, _kappa_for(v5b, high))
    passes = {lab: _gauss_passes(v5b, ks, b, n, rng, rho=rho) for lab, b, n, rho in GAUSS_CASES}
    print("v7 gaussian rates [low, c0, high] with CP95 upper bounds",
          {k: _rate_report(v5b, v) for k, v in passes.items()}, "region", (low, c0, high))  # fmt: skip
    return passes


@pytest.mark.calibration
@pytest.mark.parametrize("label", [c[0] for c in GAUSS_CASES])
def test_v7_pipeline_gaussian_false_acceptance(
    v5b: ModuleType, gauss_rates: dict, label: str
) -> None:
    """Gaussian pipeline (3600 pairs). Profiles (12 bins, replicate-level bin correlation 0 and 0.8): the
    95 % Clopper-Pearson upper bound of the false-acceptance rate is <= 0.05 at the true paired coverage
    c0 - 0.03 and c0 + 0.03. Scalar: the EXACT false-acceptance probability (``test_scalar_exact_false_
    acceptance_at_both_margins``, 0.0387 <= 0.04) is the gate; here the simulated rate is only a
    consistency check of the implementation and must lie within 3 binomial standard errors of it."""
    rep = _rate_report(v5b, gauss_rates[label])
    if label != "scalar":
        assert rep[0]["cp95_upper"] <= 0.05 and rep[2]["cp95_upper"] <= 0.05, (label, rep)
        return
    kl, kc = _scalar_acceptance(v5b)
    _, low, high = v5b.v7_region()
    for r, p in ((rep[0], low), (rep[2], high)):
        exact = _binom_interval_prob(3600, p, kl, kc)
        se = math.sqrt(exact * (1.0 - exact) / r["n"])
        assert abs(r["rate"] - exact) <= 3.0 * se, (r, exact, se)


@pytest.mark.calibration
def test_v7_pipeline_gaussian_power_targets(gauss_rates: dict) -> None:
    """Power at the nominal coverage (kappa = 1): profiles >= 0.95, scalar >= 0.90, joint row
    (independent draws of the two profiles and the scalar) >= 0.85."""
    power = {k: float(v[1].mean()) for k, v in gauss_rates.items()}
    n = 1000
    joint = float((gauss_rates["profile12_rho0.0"][1][:n] & gauss_rates["profile12_rho0.8"][1][:n]
                   & gauss_rates["scalar"][1][:n]).mean())  # fmt: skip
    print("v7 gaussian power at nominal", power, "joint row (independent draws)", joint)
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
) -> np.ndarray:
    """Paired statistic ``|X_j - X_k| / sqrt(s_j^2 + s_k^2)`` of one bin over ``n`` pairs of independent
    replicates (20 Gamma blocks each; 0 where both SEMs vanish: an exact hit)."""
    g = _gamma_blocks(shape, (2, n, 20, 1), rng, correlated, (0,) if zero else ()).astype(float)[
        ..., 0
    ]
    mean, sem = g.mean(axis=2), g.std(axis=2, ddof=1) / math.sqrt(20)
    den = np.sqrt(sem[0] ** 2 + sem[1] ** 2)
    num = np.abs(mean[0] - mean[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(den > 0, num / den, np.where(num > 0, np.inf, 0.0))
    return r


def _skew_kappas(
    v5b: ModuleType, shape: float, correlated: bool, zero: bool, rng: np.random.Generator
) -> tuple[float, float, float]:
    """kappa giving the TRUE bin-averaged paired coverage low and high (large-sample simulation of the
    paired statistic, bisection with common random numbers over the bin marginals: 10 Gamma bins and,
    if ``zero``, 2 zero-inflated bins) and the true coverage at kappa = 1."""
    _, low, high = v5b.v7_region()
    marg = [(1.0, np.sort(_ratios(shape, correlated, False, rng)))]
    if zero:
        marg = [(10 / 12, marg[0][1]), (2 / 12, np.sort(_ratios(shape, correlated, True, rng)))]

    def cov(k: float) -> float:
        return float(sum(w * np.searchsorted(r, k, side="right") / r.size for w, r in marg))

    def solve(target: float) -> float:
        lo, hi = 0.05, 30.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if cov(mid) < target else (lo, mid)
        return 0.5 * (lo + hi)

    k_lo, k_hi = solve(low), solve(high)
    assert abs(cov(k_lo) - low) < 2e-3 and abs(cov(k_hi) - high) < 2e-3
    return k_lo, k_hi, cov(1.0)


def _skew_passes(v5b: ModuleType, shape: float, correlated: bool, zero: bool, kappas: tuple[float, ...],
                 trials: int, rng: np.random.Generator, bins: int = 12, reps: int = 7200) -> np.ndarray:  # fmt: skip
    """Pass indicators ``[len(kappas), trials]`` on Gamma-block tallies (20 blocks per replicate and bin,
    replicate mean and SEM from the blocks, per-bin scales 0.5..2): the implemented paired procedure."""
    scale = np.linspace(0.5, 2.0, bins)
    out = np.zeros((len(kappas), trials), dtype=bool)
    for i in range(trials):
        g = (
            _gamma_blocks(shape, (reps, 20, bins), rng, correlated, ZERO_BINS if zero else ())
            * scale
        )
        mean, sem = (
            g.mean(axis=1).astype(float),
            (g.std(axis=1, ddof=1) / math.sqrt(20)).astype(float),
        )
        out[:, i] = _verdicts(v5b, "profile", mean, sem, kappas)
    return out


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
    k_lo, k_hi, true_nominal = _skew_kappas(v5b, shape, corr, zero, rng)
    rep = _rate_report(v5b, _skew_passes(v5b, shape, corr, zero, (k_lo, 1.0, k_hi), 500, rng))
    out = {"shape": shape, "correlated": corr, "zero_inflated": zero, "kappa_low": k_lo, "kappa_high": k_hi,
           "true_coverage_at_kappa1": true_nominal, "FA_low": rep[0], "power_kappa1": rep[1], "FA_high": rep[2]}  # fmt: skip
    print("v7 skewed-tally calibration", out)
    return out


@pytest.mark.calibration
def test_v7_pipeline_skewed_tally_false_acceptance(skew_result: dict) -> None:
    """The 95 % Clopper-Pearson upper bound of the false-acceptance rate is <= 0.05 at the true paired
    coverage c0 - 0.03 and c0 + 0.03 for Gamma-block tallies: independent, correlated (shared-history)
    and zero-inflated bins."""
    assert (
        skew_result["FA_low"]["cp95_upper"] <= 0.05 and skew_result["FA_high"]["cp95_upper"] <= 0.05
    ), skew_result


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
    ref = _guard_env(v5b, lambda: v5b.batch_estimates(cfg, 2))
    for eb in (4, None):
        got = v5b.batch_estimates(cfg, eb)
        for key in ("sec_p", "nuc_local_dose", "idd", "nuclear_local", "escaped_neutral"):
            assert got[key].shape == ref[key].shape and np.array_equal(got[key], ref[key]), (
                key,
                eb,
            )
        assert got["counters_sum"] == ref["counters_sum"] == 0
