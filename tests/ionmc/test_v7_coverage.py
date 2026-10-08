# ruff: noqa: E501
"""V7 replicate-coverage rule (plan Amendment 13, Codex REVIEW-be30e621, -873ab9cd, -edb9970b, -f1b32614,
-852837a3): held-out reference (shards 6-7; evaluation replicates exactly i.i.d. given it), Bonferroni box with a t
quantile checked by a bootstrap, exact worst-case coverage by an endpoint sweep, Clopper-Pearson (scalar) and empirical
Bernstein (profiles) bounds at alpha_tost = 0.04, combine of shard/reference partials (fail closed), and (marker ``calibration``, minutes, not in CI) Monte Carlo
calibration of the whole pipeline with Clopper-Pearson bounds of the false-acceptance rates: Gaussian
profiles and scalar, Gamma-block tallies with independent and with correlated (shared-history) bins."""

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
PHI1 = 0.24197072451914337


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


def _synthetic_covered(mean: float, sd: float, bins: int = 12, r: int = 4500) -> np.ndarray:
    """Deterministic integer-free covered counts with exact mean and sd of f_j = covered / bins."""
    z = np.linspace(-1.0, 1.0, r)
    z = (z - z.mean()) / z.std(ddof=1)
    return (mean + sd * z) * bins


def test_student_quantiles_and_nominal_coverage(v5b: ModuleType) -> None:
    for p, df, ref in ((0.95, 19, 1.7291), (0.95, 29, 1.6991), (0.95, 4499, 1.6453)):
        assert abs(v5b.student_t_quantile(p, df) - ref) < 5e-4
    assert abs(v5b.student_abs_prob(19, 1.0) - 0.670) < 0.002
    assert abs(v5b.student_abs_prob(1, 1.0) - 0.5) < 1e-9  # Cauchy: (2/pi) atan(1)
    assert abs(v5b.PHI_1 - PHI1) < 1e-15


def test_reference_correction_matches_hand_values(v5b: ModuleType) -> None:
    sem = np.array([0.1, 0.2, 0.3, 0.1])
    sigma = np.array([1.0, 1.0, 1.0, 0.5])
    r = np.array([0.1, 0.2, 0.3, 0.2])
    c = v5b.v7_reference_correction(sem, sigma)
    assert np.allclose(c["r_b"], r)
    assert abs(c["b_ref"] - PHI1 * np.mean(r**2)) < 1e-15
    assert abs(c["b_ref"] - PHI1 * (0.01 + 0.04 + 0.09 + 0.04) / 4) < 1e-15
    assert abs(c["s_ref_bound"] - PHI1 * math.sqrt(2.0) * np.mean(r**2)) < 1e-15
    assert math.isnan(c["s_ref_cov"]) and math.isnan(c["ref_corr_mean_abs_offdiag"])
    one = v5b.v7_reference_correction(np.array([0.1]), np.array([1.0]))  # scalar: B = 1
    assert abs(one["b_ref"] - PHI1 * 0.01) < 1e-15
    assert abs(one["s_ref_bound"] - PHI1 * math.sqrt(2.0) * 1e-2) < 1e-15


def test_s_ref_cov_known_covariance_and_bound(v5b: ModuleType) -> None:
    k, bins = 20, 6
    rng = np.random.default_rng(3)
    sd = np.linspace(0.5, 2.0, bins)
    for rho in (0.0, 0.6, 0.999):
        corr = np.full((bins, bins), rho) + (1.0 - rho) * np.eye(bins)
        z = rng.standard_normal((k, bins)) @ np.linalg.cholesky(corr).T
        z = (z - z.mean(axis=0)) / z.std(axis=0, ddof=1)  # sample moments only fix the diagonal
        blocks = 7.0 + z * sd
        cov = np.cov(blocks, rowvar=False, ddof=1) / k  # covariance of the 20-block means
        sem = np.sqrt(np.diag(cov))
        assert np.allclose(sem, sd / math.sqrt(k))
        sigma = np.linspace(1.0, 3.0, bins)  # replicate sd
        c = v5b.v7_reference_correction(sem, sigma, blocks)
        expect = PHI1 / bins * math.sqrt(2.0 * np.sum((cov / np.outer(sigma, sigma)) ** 2))
        assert abs(c["s_ref_cov"] - expect) < 1e-14
        rc = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
        off = (np.abs(rc).sum() - bins) / (bins * (bins - 1))
        assert abs(c["ref_corr_mean_abs_offdiag"] - off) < 1e-14
        assert c["s_ref_cov"] <= c["s_ref_bound"] + 1e-15
    # exactly proportional blocks (all correlations 1): the bound is attained
    x = rng.standard_normal((k, 1))
    full = 1.0 + x * np.array([0.5, 1.0, 2.0])
    sem = full.std(axis=0, ddof=1) / math.sqrt(k)
    c = v5b.v7_reference_correction(sem, np.array([1.0, 1.0, 1.0]), full)
    assert abs(c["s_ref_cov"] - c["s_ref_bound"]) < 1e-12
    assert abs(c["ref_corr_mean_abs_offdiag"] - 1.0) < 1e-12
    # independent-bin term (the old s_ref) is a lower end of the bound
    r = sem / 1.0
    assert PHI1 * math.sqrt(2.0) * math.sqrt(np.sum(r**4)) / 3 <= c["s_ref_bound"] + 1e-15
    # single bin: covariance is the variance; cov value equals the bound
    one = v5b.v7_reference_correction(sem[:1], np.array([1.0]), full[:, :1])
    assert abs(one["s_ref_cov"] - one["s_ref_bound"]) < 1e-12 and math.isnan(
        one["ref_corr_mean_abs_offdiag"]
    )


def test_equivalence_verdict_vectors(v5b: ModuleType) -> None:
    assert (v5b.V7_COV_LOW, v5b.V7_COV_HIGH, v5b.V7_ALPHA_BOX, v5b.V7_ALPHA_TOST) == (
        0.640,
        0.700,
        0.01,
        0.04,
    )
    assert (v5b.V7_BOOT_N, v5b.V7_BOOT_Q, v5b.V7_REPLICATES) == (5000, 0.99, 7200)
    ok = v5b.v7_tost_verdict(0.650, 0.690, 12, 12 * 5400, min_bins=10)
    assert ok["pass"] and not ok["reasons"] and ok["ci"] == [0.650, 0.690]
    lo = v5b.v7_tost_verdict(0.6399, 0.690, 12, 12 * 5400, min_bins=10)
    assert not lo["pass"] and any("lower" in x for x in lo["reasons"])
    hi = v5b.v7_tost_verdict(0.650, 0.7001, 12, 12 * 5400, min_bins=10)
    assert not hi["pass"] and any("upper" in x for x in hi["reasons"])
    assert v5b.v7_tost_verdict(0.640, 0.700, 1, 5400)["pass"]  # closed region
    assert not v5b.v7_tost_verdict(0.650, 0.690, 12, 299, min_bins=10)["pass"]
    assert not v5b.v7_tost_verdict(0.650, 0.690, 9, 9 * 5400, min_bins=10)["pass"]


def test_clopper_pearson_known_values_and_binomial_consistency(v5b: ModuleType) -> None:
    """n = 10, k = 3, one-sided 95 %: [0.0873, 0.6065] (the two-sided 90 % Clopper-Pearson interval)."""
    assert (
        abs(v5b.cp_lower(3, 10, 0.05) - 0.08726) < 5e-5
        and abs(v5b.cp_upper(3, 10, 0.05) - 0.60662) < 5e-5
    )
    assert v5b.cp_lower(0, 10, 0.05) == 0.0 and v5b.cp_upper(10, 10, 0.05) == 1.0
    assert (
        abs(v5b.cp_upper(0, 10, 0.05) - (1.0 - 0.05**0.1)) < 1e-9
    )  # rule of three-like closed form
    assert abs(v5b.cp_lower(10, 10, 0.05) - 0.05**0.1) < 1e-9
    for k, n, a in ((7, 40, 0.04), (3000, 5400, 0.04), (5, 5400, 0.04)):
        lo, up = v5b.cp_lower(k, n, a), v5b.cp_upper(k, n, a)
        assert lo < k / n < up

        def tail_ge(p: float, k: int = k, n: int = n) -> float:  # P(X >= k | p) by direct summation
            return sum(
                math.exp(
                    math.lgamma(n + 1)
                    - math.lgamma(j + 1)
                    - math.lgamma(n - j + 1)
                    + j * math.log(p)
                    + (n - j) * math.log1p(-p)
                )
                for j in range(k, n + 1)
            )

        def tail_le(p: float, k: int = k, n: int = n) -> float:
            return sum(
                math.exp(
                    math.lgamma(n + 1)
                    - math.lgamma(j + 1)
                    - math.lgamma(n - j + 1)
                    + j * math.log(p)
                    + (n - j) * math.log1p(-p)
                )
                for j in range(0, k + 1)
            )

        assert abs(tail_ge(lo) - a) < 1e-8 and abs(tail_le(up) - a) < 1e-8


def test_empirical_bernstein_hand_value(v5b: ModuleType) -> None:
    v, n, a = 0.0182, 5400, 0.04
    ln = math.log(50.0)
    want = math.sqrt(2 * v * ln / n) + 7 * ln / (3 * (n - 1))
    assert abs(v5b.empirical_bernstein_eps(v, n, a) - want) < 1e-15 and abs(want - 0.00680) < 5e-5
    assert v5b.empirical_bernstein_eps(0.0, n, a) == 7 * ln / (3 * (n - 1))  # variance-free term


def test_box_quantile_t_and_normal(v5b: ModuleType) -> None:
    assert (
        abs(v5b.normal_quantile(0.975) - 1.959964) < 1e-5 and abs(v5b.normal_quantile(0.5)) < 1e-12
    )
    z12 = v5b.student_t_quantile(1.0 - 0.01 / 24.0, 1799)
    z1 = v5b.student_t_quantile(1.0 - 0.01 / 2.0, 1799)
    assert (
        abs(z12 - 3.3415) < 0.03 and abs(z1 - 2.5758) < 0.02
    )  # t_1799 is slightly above the normal


def test_coverage_sweep_matches_brute_force(v5b: ModuleType) -> None:
    """The exact sweep equals a brute-force evaluation of C(delta) at all clipped breakpoints, at
    +-Z and at the midpoints between them (a fine grid cannot beat it), on random small inputs; the
    returned shifts attain the extrema."""
    rng = np.random.default_rng(11)
    for _ in range(200):
        n = int(rng.integers(1, 9))
        x = rng.normal(0.0, 2.0, n)
        sem = rng.uniform(0.05, 2.0, n)
        z = float(rng.uniform(0.3, 4.0))
        lo, hi = x - sem, x + sem
        pts = np.unique(
            np.concatenate([[-z, z], lo[(lo > -z) & (lo < z)], hi[(hi > -z) & (hi < z)]])
        )
        cand = np.concatenate([pts, 0.5 * (pts[:-1] + pts[1:])])

        def cov(d: float, x=x, sem=sem) -> float:  # noqa: ANN001
            return float(np.mean(np.abs(x - d) <= sem))

        brute = [cov(d) for d in cand]
        mn, mx, dmn, dmx = v5b.v7_coverage_extrema(x, sem, z)
        assert mn == min(brute) and mx == max(brute)
        assert cov(dmn) == mn and cov(dmx) == mx and -z <= dmn <= z and -z <= dmx <= z
        fine = [cov(d) for d in np.linspace(-z, z, 4001)]  # a fine grid never goes beyond the sweep
        assert mn <= min(fine) + 1e-12 and mx >= max(fine) - 1e-12
    assert v5b.v7_coverage_extrema(np.array([0.5]), np.array([0.5]), 2.0)[:2] == (0.0, 1.0)
    assert v5b.v7_coverage_extrema(np.array([0.0]), np.array([5.0]), 2.0)[:2] == (
        1.0,
        1.0,
    )  # covers the box


def test_bootstrap_quantile_matches_t_on_gaussian_means(v5b: ModuleType) -> None:
    rng = np.random.default_rng(5)
    r, nh, b = 4000, 1000, 3
    means = 10.0 + rng.normal(0.0, 1.0, (r, b))
    sem_h = means[-nh:].std(axis=0, ddof=1) / math.sqrt(nh)
    z1 = v5b.v7_boot_z(means, nh, means.mean(axis=0), sem_h, 7, 4000)
    z2 = v5b.v7_boot_z(means, nh, means.mean(axis=0), sem_h, 7, 4000)
    z3 = v5b.v7_boot_z(means, nh, means.mean(axis=0), sem_h, 8, 4000)
    zt = v5b.student_t_quantile(1.0 - 0.01 / (2 * b), nh - 1)
    assert z1 == z2 and z1 != z3  # fixed seed reproducible
    assert abs(z1 - zt) < 0.3 and abs(z3 - zt) < 0.3, (z1, z3, zt)


def test_heldout_reference_never_contains_evaluation_replicates(v5b: ModuleType) -> None:
    """The reference mean, its SEM and the bin selection come from the held-out quarter alone."""
    rng = np.random.default_rng(2)
    x = 5.0 + rng.normal(0.0, 1.0, (16, 3))
    x[
        :12, 2
    ] = -1.0  # bin 2 is non-positive in E but positive in H: it is used (selection from H alone)
    x[12:, 1] = (
        np.abs(x[12:, 1]) * 0.0
    )  # bin 1 vanishes in H (zero spread): not used, although E is positive
    sem = np.full_like(x, 0.8)
    ref, rsem = np.array([5.0, 5.0, 5.0]), np.array([0.1, 0.1, 0.1])
    e0 = v5b.v7_estimator_verdict("profile", x, sem, ref, rsem, None, boot_seed=1, n_boot=200)
    assert (
        v5b.v7_heldout_split(16) == 12
        and e0["replicates_eval"] == 12
        and e0["replicates_heldout"] == 4
    )
    assert np.allclose(e0["mu_heldout"], x[12:].mean(axis=0))
    assert np.allclose(e0["sem_heldout"], x[12:].std(axis=0, ddof=1) / 2.0)
    assert e0["bins_used"] == 2 and e0["intervals"] == 12 * 2
    hits = (np.abs(x[:12] - x[12:].mean(axis=0)) <= 0.8)[:, [0, 2]]
    assert e0["per_replicate_covered"] == hits.sum(axis=1).tolist()
    x2 = x.copy()
    x2[3, 0] += 0.3  # an evaluation replicate: reference, SEM and bin selection are unchanged
    e1 = v5b.v7_estimator_verdict("profile", x2, sem, ref, rsem, None, boot_seed=1, n_boot=200)
    assert e1["mu_heldout"] == e0["mu_heldout"] and e1["sem_heldout"] == e0["sem_heldout"]
    assert e1["bins_used"] == e0["bins_used"]
    x3 = x.copy()
    x3[13, 0] += 0.3  # a held-out replicate moves the reference
    e2 = v5b.v7_estimator_verdict("profile", x3, sem, ref, rsem, None, boot_seed=1, n_boot=200)
    assert abs(e2["mu_heldout"][0] - (e0["mu_heldout"][0] + 0.075)) < 1e-12
    with pytest.raises(SystemExit, match="3:1"):
        v5b.v7_heldout_split(10)


def test_estimator_verdict_hand_case_and_fields(v5b: ModuleType) -> None:
    """Scalar, 8 replicates: E = first 6, H = last 2; hand-checked reference, hits, box and CP bounds."""
    means = np.array([[1.02], [1.21], [0.83], [0.98], [0.91], [1.13], [1.0], [1.2]])
    sems = np.full((8, 1), 0.14)
    ref, rsem = np.array([1.0]), np.array([0.1])
    e = v5b.v7_estimator_verdict("scalar", means, sems, ref, rsem, None, boot_seed=3, n_boot=300)
    mu_h, sem_h = 1.1, float(np.std([1.0, 1.2], ddof=1)) / math.sqrt(2.0)  # 0.1
    assert abs(e["mu_heldout"][0] - mu_h) < 1e-15 and abs(e["sem_heldout"][0] - sem_h) < 1e-15
    assert e["per_replicate_covered"] == [
        1,
        1,
        0,
        1,
        0,
        1,
    ]  # |x - 1.1| <= 0.15 for x = 1.02, 1.21, 0.83, 0.98, 0.91, 1.13
    zt = v5b.student_t_quantile(1.0 - 0.01 / 2.0, 1)
    assert (
        abs(e["z_t"] - zt) < 1e-12
        and e["z_box"] == max(e["z_t"], e["z_boot"])
        and e["boot_seed"] == 3
    )
    x = (means[:6, 0] - mu_h) / sem_h
    mn, mx, _, _ = v5b.v7_coverage_extrema(x, np.full(6, 0.14 / sem_h), e["z_box"])
    assert e["per_bin_box_coverage_min_max"] == [[mn, mx]] and e["m_lo"] == mn and e["m_hi"] == mx
    assert e["bound_kind"] == "clopper-pearson"
    assert e["cp_count_min"] == round(mn * 6) and e["cp_count_max"] == round(mx * 6)
    assert abs(e["ci"][0] - v5b.cp_lower(e["cp_count_min"], 6, 0.04)) < 1e-12
    assert abs(e["ci"][1] - v5b.cp_upper(e["cp_count_max"], 6, 0.04)) < 1e-12
    assert e["m_ref1e6"] == float(np.mean(np.abs(means[:6, 0] - 1.0) <= 0.14)) and not e["pass"]
    # profile: empirical Bernstein from the largest variance of f_j; ci = [m_lo - eps, m_hi + eps]
    rng = np.random.default_rng(1)
    big = 5.0 + rng.normal(0.0, 1.0, (40, 4))
    p = v5b.v7_estimator_verdict("profile", big, np.full_like(big, 0.9), np.full(4, 5.0), np.full(4, 0.1),
                                 None, boot_seed=2, n_boot=200)  # fmt: skip
    assert p["bound_kind"] == "empirical-bernstein" and p["replicates_eval"] == 30
    eps = v5b.empirical_bernstein_eps(p["variance_max"], 30, 0.04)
    assert abs(p["eb_eps"] - eps) < 1e-15 and abs(p["ci"][0] - (p["m_lo"] - eps)) < 1e-15
    assert abs(p["ci"][1] - (p["m_hi"] + eps)) < 1e-15 and p["m_lo"] <= p["m"] <= p["m_hi"]
    zt4 = v5b.student_t_quantile(1.0 - 0.01 / 8.0, 9)
    assert abs(p["z_t"] - zt4) < 1e-12 and p["z_box"] >= p["z_t"]
    assert len(p["per_bin_box_coverage_min_max"]) == 4 and len(p["z_b_eval_vs_heldout"]) == 4


def test_scalar_cp_rule_at_nominal_and_edges(v5b: ModuleType) -> None:
    """R_E = 5400: 3618 hits (0.670) has CP bounds at 0.04 per side inside [0.640, 0.700]; the lower
    bound crosses 0.640 between 3500 and 3530 hits."""
    assert v5b.cp_lower(3618, 5400, 0.04) > 0.640 and v5b.cp_upper(3618, 5400, 0.04) < 0.700
    assert v5b.cp_lower(3500, 5400, 0.04) < 0.640 < v5b.cp_lower(3530, 5400, 0.04)


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
    assert rule["alpha_box"] == 0.01 and rule["alpha_tost"] == 0.04 and rule["boot_n"] == 5000
    assert rule["evaluation_shards"] == "0-5" and rule["heldout_shards"] == "6-7"
    for name, kind, bins in (
        ("sec_p", "profile", 12),
        ("nuclear_local", "profile", 12),
        ("escaped_neutral", "scalar", 1),
    ):
        e = doc["estimators"][name]
        assert e["kind"] == kind and e["bins_used"] == bins and e["intervals"] == 5400 * bins
        assert e["replicates_eval"] == 5400 and e["replicates_heldout"] == 1800
        assert len(e["per_replicate_covered"]) == 5400 and all(
            isinstance(c, int) for c in e["per_replicate_covered"]
        )
        assert abs(e["nominal_coverage"] - 0.670) < 0.002
        for key in (
            *(
                "m",
                "m_lo",
                "m_hi",
                "z_t",
                "z_boot",
                "z_box",
                "alpha_box",
                "alpha_tost",
                "boot_seed",
            ),
            *(
                "bound_kind",
                "variance_max",
                "eb_eps",
                "cp_count_min",
                "cp_count_max",
                "ci",
                "reasons",
            ),
            *("mu_heldout", "sem_heldout", "per_bin_box_coverage_min_max"),
            *("z_b_eval_vs_heldout", "max_abs_z_eval_vs_heldout"),
            *("s_ref_cov", "ref_corr_mean_abs_offdiag", "replicate_mean_skewness"),
            *("legacy_point_gate_pass", "pass", "m_ref1e6", "z_b_ref1e6", "max_abs_z_ref1e6"),
        ):
            assert key in e
        assert not {"hoeffding_eps", "fold_means", "t", "r_b"} & set(e)
        assert "held-out" in e["reference_kind"] and e["boot_seed"] == ref["seed"]
        assert e["bound_kind"] == ("clopper-pearson" if kind == "scalar" else "empirical-bernstein")
        zt = v5b.student_t_quantile(1.0 - 0.01 / (2.0 * bins), 1799)
        assert abs(e["z_t"] - zt) < 1e-12 and e["z_box"] == max(e["z_t"], e["z_boot"])
        assert abs(e["z_t"] - (3.3415 if bins == 12 else 2.5758)) < 0.03
        assert abs(e["z_boot"] - e["z_t"]) < 0.5  # Gaussian replicates: bootstrap ~ t quantile
        assert len(e["per_bin_box_coverage_min_max"]) == bins and e["m_lo"] <= e["m"] <= e["m_hi"]
        assert e["max_abs_z_ref1e6"] < 6.0 and e["max_abs_z_eval_vs_heldout"] < 6.0
        assert len(e["replicate_mean_skewness"]) == bins
    assert doc["pass"] == all(e["pass"] for e in doc["estimators"].values())
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
    with pytest.raises(SystemExit, match="quarter"):
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
    badblocks = copy.deepcopy(ref)
    badblocks["estimators"]["sec_p"]["blocks"] = badblocks["estimators"]["sec_p"]["blocks"][:5]
    with pytest.raises(SystemExit, match="inconsistent"):
        v5b.v7_rep_combine(shards, badblocks)


def test_replicate_grouping_and_estimator_verdict(v5b: ModuleType) -> None:
    b = np.arange(80.0).reshape(40, 2)  # 2 replicates of 20 batches, 2 bins
    m, s = v5b.v7_replicate_stats(b)
    assert m.shape == s.shape == (2, 2)
    assert np.allclose(m[0], b[:20].mean(axis=0)) and np.allclose(m[1], b[20:].mean(axis=0))
    assert np.allclose(s[1], b[20:].std(axis=0, ddof=1) / math.sqrt(20))
    m1, _ = v5b.v7_replicate_stats(np.arange(40.0))  # scalar estimator: [B] -> [R, 1]
    assert m1.shape == (2, 1)


# -- Monte Carlo calibration of the full pipeline (marker ``calibration``: not in CI) ---------
# The procedure is ``v7_estimator_verdict`` itself on simulated replicate means (R = 7200: E = first 5400,
# H = last 1800), with Z_boot from 500 bootstrap resamples instead of 5000 (runtime) and computed once per
# trial for all kappas (the replicate means do not depend on kappa).
N_BOOT_CAL = 500


def _kappa_for(v5b: ModuleType, coverage: float) -> float:
    """kappa with P(|t_19| <= kappa) = coverage (intervals s = kappa sigma sqrt(chi2_19 / 19))."""
    lo, hi = 0.1, 5.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if v5b.student_abs_prob(19, mid) < coverage else (lo, mid)
    return 0.5 * (lo + hi)


def _verdicts(v5b: ModuleType, kind: str, mean: np.ndarray, sem0: np.ndarray, kappas: tuple[float, ...],
              seed: int) -> list[bool]:  # fmt: skip
    """``v7_estimator_verdict`` pass flags at each kappa (half width kappa * sem0) on one simulated
    sample; the 1e6 reference arguments are dummies (diagnostic only)."""
    ref, rsem = mean.mean(axis=0), np.full(mean.shape[1], 0.01)
    out, zb = [], None
    for kappa in kappas:
        e = v5b.v7_estimator_verdict(kind, mean, kappa * sem0, ref, rsem, None, boot_seed=seed,
                                     n_boot=N_BOOT_CAL, z_boot_override=zb)  # fmt: skip
        zb = e["z_boot"]
        out.append(bool(e["pass"]))
    return out


def _gauss_passes(v5b: ModuleType, kappas: tuple[float, ...], bins: int, trials: int,
                  rng: np.random.Generator, rho: float = 0.0, reps: int = 7200) -> np.ndarray:  # fmt: skip
    """Pass indicators ``[len(kappas), trials]``. Replicate means X_b = sigma_b (sqrt(rho) g + sqrt(1 -
    rho) e_b) with a replicate-level common factor g (bin correlation ``rho``, as for bins sharing
    histories), batch SEMs s = kappa sigma_b sqrt(chi2_19 / 19); true mean 0 (shifted to 100)."""
    sigma = np.linspace(0.5, 2.0, bins) if bins > 1 else np.array([1.3])
    out = np.zeros((len(kappas), trials), dtype=bool)
    kind = "profile" if bins > 1 else "scalar"
    for i in range(trials):
        g = rng.standard_normal((reps, 1))
        x = sigma * (math.sqrt(rho) * g + math.sqrt(1.0 - rho) * rng.standard_normal((reps, bins)))
        s0 = sigma * np.sqrt(rng.chisquare(19, (reps, bins)) / 19.0)
        out[:, i] = _verdicts(v5b, kind, 100.0 + x, s0, kappas, 1000 + i)
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
    ("scalar", 1, 2000, 0.0),
)


@pytest.fixture(scope="module")
def gauss_rates(v5b: ModuleType) -> dict[str, np.ndarray]:
    """Pass indicators ``[3, trials]`` (true coverage 0.640 / 0.670 / 0.700) per Gaussian case."""
    rng = np.random.default_rng(20471004)
    ks = (_kappa_for(v5b, 0.640), 1.0, _kappa_for(v5b, 0.700))
    assert abs(v5b.student_abs_prob(19, 1.0) - 0.670) < 0.002
    passes = {lab: _gauss_passes(v5b, ks, b, n, rng, rho=rho) for lab, b, n, rho in GAUSS_CASES}
    print("v7 gaussian rates [0.640, 0.670, 0.700] with CP95 upper bounds",
          {k: _rate_report(v5b, v) for k, v in passes.items()})  # fmt: skip
    return passes


@pytest.mark.calibration
@pytest.mark.parametrize("label", [c[0] for c in GAUSS_CASES])
def test_v7_pipeline_gaussian_false_acceptance(
    v5b: ModuleType, gauss_rates: dict, label: str
) -> None:
    """Gaussian pipeline (R_E = 5400, n_H = 1800), 12 bins with replicate-level bin correlation 0 and
    0.8 and the scalar: the 95 % Clopper-Pearson upper bound of the false-acceptance rate is <= 0.05
    at true coverage 0.640 and at 0.700."""
    rep = _rate_report(v5b, gauss_rates[label])
    assert rep[0]["cp95_upper"] <= 0.05 and rep[2]["cp95_upper"] <= 0.05, (label, rep)


@pytest.mark.calibration
def test_v7_pipeline_gaussian_power_targets(gauss_rates: dict) -> None:
    """Power at the nominal coverage (kappa = 1, true coverage 0.670): profiles >= 0.95, scalar >= 0.95,
    joint row (independent draws of the two profiles and the scalar) >= 0.90."""
    power = {k: float(v[1].mean()) for k, v in gauss_rates.items()}
    n = 600
    joint = float((gauss_rates["profile12_rho0.0"][1][:n] & gauss_rates["profile12_rho0.8"][1][:n]
                   & gauss_rates["scalar"][1][:n]).mean())  # fmt: skip
    print("v7 gaussian power at nominal", power, "joint row (independent draws)", joint)
    assert power["profile12_rho0.0"] >= 0.95 and power["profile12_rho0.8"] >= 0.95, power
    assert power["scalar"] >= 0.95 and joint >= 0.90, (power, joint)


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


def _ratios(shape: float, correlated: bool, zero: bool, rng: np.random.Generator) -> np.ndarray:
    """|replicate mean - true mean| / replicate SEM of one bin over 400000 replicates (inf for SEM 0)."""
    g = _gamma_blocks(shape, (400_000, 20, 1), rng, correlated, (0,) if zero else ()).astype(float)[
        ..., 0
    ]
    true = shape * (ZERO_KEEP if zero else 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.abs(g.mean(axis=1) - true) / (g.std(axis=1, ddof=1) / math.sqrt(20))
    return np.where(np.isnan(r), np.inf, r)


def _skew_kappas(
    shape: float, correlated: bool, zero: bool, rng: np.random.Generator
) -> tuple[float, float, float]:
    """kappa giving TRUE bin-averaged coverage 0.640 and 0.700 (bisection with common random numbers
    over the bin marginals: 10 Gamma bins and, if ``zero``, 2 zero-inflated bins) and the true coverage
    at kappa = 1."""
    marg = [(10 / 12, np.sort(_ratios(shape, correlated, False, rng)))]
    marg.append((2 / 12, np.sort(_ratios(shape, correlated, True, rng)))) if zero else None
    if not zero:
        marg = [(1.0, marg[0][1])]

    def cov(k: float) -> float:
        return float(sum(w * np.searchsorted(r, k, side="right") / r.size for w, r in marg))

    def solve(target: float) -> float:
        lo, hi = 0.05, 30.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if cov(mid) < target else (lo, mid)
        return 0.5 * (lo + hi)

    k_lo, k_hi = solve(0.640), solve(0.700)
    assert abs(cov(k_lo) - 0.640) < 2e-3 and abs(cov(k_hi) - 0.700) < 2e-3
    return k_lo, k_hi, cov(1.0)


def _skew_passes(v5b: ModuleType, shape: float, correlated: bool, zero: bool, kappas: tuple[float, ...],
                 trials: int, rng: np.random.Generator, bins: int = 12, reps: int = 7200) -> np.ndarray:  # fmt: skip
    """Pass indicators ``[len(kappas), trials]`` on Gamma-block tallies (20 blocks per replicate and bin,
    replicate mean and SEM from the blocks, half width kappa * SEM, per-bin scales 0.5..2)."""
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
        out[:, i] = _verdicts(v5b, "profile", mean, sem, kappas, 2000 + i)
    return out


SKEW_CASES = [(sh, corr, False) for sh in (2.0, 0.5) for corr in (False, True)] + [
    (2.0, False, True)
]


@pytest.fixture(scope="module", params=SKEW_CASES,
                ids=[f"gamma{sh}-{'correlated' if c else 'independent'}{'-zeroinflated' if z else ''}" for sh, c, z in SKEW_CASES])  # fmt: skip
def skew_result(request: pytest.FixtureRequest, v5b: ModuleType) -> dict:
    """Gamma shape 2 and 0.5, independent / correlated sparse bins and a zero-inflated case (two of the
    12 bins have 90 % exact zeros): kappa scales the SEM so that the TRUE bin-averaged coverage is 0.640
    and 0.700; pass rates with CP bounds there, power at kappa = 1 and the true coverage at kappa = 1."""
    shape, corr, zero = request.param
    rng = np.random.default_rng(
        20471005 + int(2 * shape) + (7 if corr else 0) + (13 if zero else 0)
    )
    k_lo, k_hi, true_nominal = _skew_kappas(shape, corr, zero, rng)
    rep = _rate_report(v5b, _skew_passes(v5b, shape, corr, zero, (k_lo, 1.0, k_hi), 400, rng))
    out = {"shape": shape, "correlated": corr, "zero_inflated": zero, "kappa_0.640": k_lo, "kappa_0.700": k_hi,
           "true_coverage_at_kappa1": true_nominal, "FA_0.640": rep[0], "power_kappa1": rep[1], "FA_0.700": rep[2]}  # fmt: skip
    print("v7 skewed-tally calibration", out)
    return out


@pytest.mark.calibration
def test_v7_pipeline_skewed_tally_false_acceptance(skew_result: dict) -> None:
    """The 95 % Clopper-Pearson upper bound of the false-acceptance rate is <= 0.05 at true coverage
    0.640 and 0.700 for Gamma-block tallies: independent, correlated (shared-history) and zero-inflated."""
    assert (
        skew_result["FA_0.640"]["cp95_upper"] <= 0.05
        and skew_result["FA_0.700"]["cp95_upper"] <= 0.05
    ), skew_result


def test_calibration_marker_registered_and_excluded_from_ci() -> None:
    """The Monte Carlo calibration is out of the lightweight CI selection and registered as a marker."""
    ci = (REPO / ".github" / "workflows" / "ci.yml").read_text()
    assert '-m "not cuda and not host and not calibration"' in ci
    assert '"calibration:' in (REPO / "pyproject.toml").read_text()
