# ruff: noqa: E501
"""V7 replicate-coverage rule (plan Amendment 13, Codex REVIEW-be30e621, REVIEW-873ab9cd, REVIEW-edb9970b):
leave-one-out pooled reference with a Bonferroni box over the bins, TOST at alpha_tost = 0.04, combine of
shard/reference partials (fail closed), and (marker ``calibration``, minutes, not in CI) Monte Carlo
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


def test_tost_rule_synthetic_vectors(v5b: ModuleType) -> None:
    assert (v5b.V7_COV_LOW, v5b.V7_COV_HIGH, v5b.V7_ALPHA_BOX, v5b.V7_ALPHA_TOST) == (
        0.640,
        0.700,
        0.01,
        0.04,
    )
    n = v5b.V7_REPLICATES * 12
    assert v5b.V7_REPLICATES == 5400

    def verdict(cov: np.ndarray, m_lo=None, m_hi=None, **kw: object) -> dict:  # type: ignore[no-untyped-def]
        return v5b.v7_tost_verdict(cov, 12, kw.pop("intervals", n), m_lo, m_hi, min_bins=10, **kw)

    ok = verdict(_synthetic_covered(0.670, 0.135, r=5400))
    assert ok["pass"] and not ok["reasons"] and abs(ok["m"] - 0.670) < 1e-9
    assert abs(ok["s"] - 0.135 / math.sqrt(5400)) < 1e-9 and ok["replicates"] == 5400
    assert abs(ok["t"] - v5b.student_t_quantile(0.98, 5399)) < 1e-12 and ok["alpha_tost"] == 0.04
    assert ok["m_lo"] == ok["m_hi"] == ok["m"]  # no box: collapses to the point
    assert abs(ok["ci"][1] - ok["ci"][0] - 2 * ok["t"] * ok["s"]) < 1e-12
    half = ok["t"] * ok["s"]
    low = verdict(_synthetic_covered(0.6425, 0.135, r=5400))  # CI crosses 0.640
    assert 0.640 < low["m"] < 0.700 and low["ci"][0] < 0.640 and not low["pass"]
    assert any("lower" in x for x in low["reasons"])
    up = verdict(_synthetic_covered(0.6985, 0.135, r=5400))
    assert up["ci"][1] > 0.700 and not up["pass"] and any("upper" in x for x in up["reasons"])
    # the box widens both sides: m_lo / m_hi enter the interval, s carries no reference term
    cov = _synthetic_covered(0.670, 0.135, r=5400)
    band = verdict(cov, 0.660, 0.680)
    assert (band["m_lo"], band["m_hi"]) == (0.660, 0.680) and band["pass"]
    assert (
        abs(band["ci"][0] - (0.660 - half)) < 1e-12 and abs(band["ci"][1] - (0.680 + half)) < 1e-12
    )
    assert not verdict(cov, 0.640 + 0.5 * half, 0.670)["pass"]
    lo_fail = verdict(cov, 0.640 + 0.5 * half, 0.670)
    assert any("lower" in x for x in lo_fail["reasons"])
    hi_fail = verdict(cov, 0.670, 0.700 - 0.5 * half)
    assert not hi_fail["pass"] and any("upper" in x for x in hi_fail["reasons"])
    br = verdict(cov, 0.69, 0.60)  # inconsistent inputs are clipped so that m_lo <= m <= m_hi
    assert br["m_lo"] == br["m_hi"] == br["m"]
    assert "s_tot" not in ok and "s_ref_bound" not in ok
    assert not verdict(cov, intervals=299)["pass"]
    assert not v5b.v7_tost_verdict(_synthetic_covered(0.670, 0.135, 9), 9, n, min_bins=10)["pass"]


def test_box_quantile_and_z_values(v5b: ModuleType) -> None:
    """Z_B = Phi^{-1}(1 - alpha_box / (2 B)): 3.34 for 12 bins, 2.58 for one bin."""
    assert abs(v5b.normal_quantile(0.975) - 1.959964) < 1e-5
    assert abs(v5b.normal_quantile(0.5)) < 1e-12
    z12 = v5b.normal_quantile(1.0 - 0.01 / 24.0)
    z1 = v5b.normal_quantile(1.0 - 0.01 / 2.0)
    assert abs(z12 - 3.3415) < 5e-4 and abs(z1 - 2.5758) < 5e-4


def test_leave_one_out_reference_is_independent_of_own_replicate(v5b: ModuleType) -> None:
    """mu_(-j) = (sum_k x_k - x_j) / (R - 1) equals the mean of the other replicates; replicate j's own
    value does not enter its reference (it enters the pooled mean only through 1 / R)."""
    x = np.array([[1.0, 5.0], [1.2, 4.0], [0.8, 6.0], [1.0, 5.5]])
    sem = np.full_like(x, 0.15)
    ref, rsem = np.array([1.0, 5.0]), np.array([0.1, 0.1])
    for j in range(4):
        direct = np.delete(x, j, axis=0).mean(axis=0)
        assert np.allclose(direct, (x.sum(axis=0) - x[j]) / 3)
    e0 = v5b.v7_estimator_verdict("profile", x, sem, ref, rsem, None)
    x2 = x.copy()
    x2[0, 0] += 0.05
    e1 = v5b.v7_estimator_verdict("profile", x2, sem, ref, rsem, None)
    assert abs(e1["pooled_mean"][0] - (e0["pooled_mean"][0] + 0.05 / 4)) < 1e-15
    hits0 = np.abs(x[:, 0] - (x[:, 0].sum() - x[:, 0]) / 3) <= 0.15  # loo = 1.0, .9333, 1.0667, 1.0
    assert hits0.tolist() == [True, False, False, True]
    assert abs(e0["pooled_sem"][0] - x[:, 0].std(ddof=1) / 2.0) < 1e-15


def test_estimator_verdict_box_and_cross_check_fields(v5b: ModuleType) -> None:
    """Hand-checkable scalar case, 4 replicates: Z_B = 2.5758, 9-point grid, per-bin min / max coverage."""
    means = np.array([[1.0], [1.2], [0.8], [1.0]])
    sems = np.array([[0.15]] * 4)
    ref, rsem = np.array([1.0]), np.array([0.1])
    e = v5b.v7_estimator_verdict("scalar", means, sems, ref, rsem, None)
    assert (
        e["per_replicate_covered"] == [1, 0, 0, 1] and e["m"] == 0.5
    )  # loo: 1.0, .9333, 1.0667, 1.0
    zb = v5b.normal_quantile(0.995)
    assert abs(e["z_box"] - zb) < 1e-12 and e["alpha_box"] == 0.01
    grid = np.linspace(-zb, zb, 9)
    assert np.allclose(e["box_grid"], grid) and grid[4] == 0.0
    sem_pool = float(means.std(ddof=1)) / 2.0
    loo = (means.sum() - means[:, 0]) / 3.0
    cov = [float((np.abs(means[:, 0] - (loo + d * sem_pool)) <= 0.15).mean()) for d in grid]
    assert e["per_bin_box_coverage_min_max"] == [[min(cov), max(cov)]]
    assert e["m_lo"] == min(cov) <= e["m"] <= e["m_hi"] == max(cov) and cov[4] == 0.5
    assert e["m_hi"] > e["m"] or e["m_lo"] < e["m"]
    assert not e["pass"] and len(e["replicate_mean_skewness"]) == 1
    # profile: m_lo / m_hi are the means over bins of the per-bin min / max
    x = np.array([[1.0, 5.0], [1.2, 4.0], [0.8, 6.0], [1.0, 5.5], [0.9, 5.2], [1.1, 4.7]])
    p = v5b.v7_estimator_verdict(
        "profile", x, np.full_like(x, 0.4), np.array([1.0, 5.0]), np.array([0.1, 0.1])
    )
    mm = np.array(p["per_bin_box_coverage_min_max"])
    assert abs(p["m_lo"] - mm[:, 0].mean()) < 1e-15 and abs(p["m_hi"] - mm[:, 1].mean()) < 1e-15
    assert abs(p["z_box"] - v5b.normal_quantile(1.0 - 0.01 / 4.0)) < 1e-12
    # cross-check against the 1e6 reference: hits |x - 1.0| <= 0.15 -> 1, 0, 0, 1
    assert e["m_ref1e6"] == 0.5
    z = (means.mean() - 1.0) / math.hypot(sem_pool, 0.1)
    assert abs(e["z_b_ref1e6"][0] - z) < 1e-15 and abs(e["max_abs_z_ref1e6"] - abs(z)) < 1e-15
    assert math.isnan(e["s_ref_cov"]) and e["b_ref"] > 0


def _binary_hits(frac: float, r: int = 5400) -> np.ndarray:
    h = np.zeros(r, dtype=np.int64)
    h[: int(round(frac * r))] = 1
    return h


def test_scalar_estimator_with_binary_hits(v5b: ModuleType) -> None:
    ok = v5b.v7_tost_verdict(_binary_hits(0.670), 1, 5400)
    assert ok["pass"] and abs(ok["m"] - 0.670) < 1e-3  # sd 0.47, s = 0.0064, CI +- 0.0105
    assert abs(ok["s"] - math.sqrt(0.67 * 0.33 * 5400 / 5399) / math.sqrt(5400)) < 1e-4
    edge = v5b.v7_tost_verdict(_binary_hits(0.650), 1, 5400)
    assert not edge["pass"] and any("lower" in x for x in edge["reasons"])


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
    assert len(shards) == 6
    doc = v5b.v7_rep_combine(shards, ref)
    assert doc["replicates"] == 5400 and set(doc["estimators"]) == {
        "sec_p",
        "nuclear_local",
        "escaped_neutral",
    }
    assert (
        doc["rule"]["alpha_box"] == 0.01
        and doc["rule"]["alpha_tost"] == 0.04
        and doc["rule"]["box_grid"] == 9
    )
    for name, kind, bins in (
        ("sec_p", "profile", 12),
        ("nuclear_local", "profile", 12),
        ("escaped_neutral", "scalar", 1),
    ):
        e = doc["estimators"][name]
        assert e["kind"] == kind and e["bins_used"] == bins and e["intervals"] == 5400 * bins
        assert len(e["per_replicate_covered"]) == 5400 and all(
            isinstance(c, int) for c in e["per_replicate_covered"]
        )
        assert len(e["r_b"]) == bins and abs(e["nominal_coverage"] - 0.670) < 0.002
        for key in (
            *(
                "m",
                "m_lo",
                "m_hi",
                "z_box",
                "alpha_box",
                "alpha_tost",
                "per_bin_box_coverage_min_max",
            ),
            *(
                "s",
                "t",
                "s_ref_cov",
                "b_ref",
                "ref_corr_mean_abs_offdiag",
                "replicate_mean_skewness",
            ),
            *("ci", "reasons", "legacy_point_gate_pass", "pass"),
            *("m_ref1e6", "z_b_ref1e6", "max_abs_z_ref1e6", "pooled_mean", "pooled_sem"),
        ):
            assert key in e
        assert not {"m_corr", "s_tot", "s_ref_bound", "means_ref_shift"} & set(e)
        assert "leave-one-out" in e["reference_kind"]
        zb = v5b.normal_quantile(1.0 - 0.01 / (2.0 * bins))
        assert abs(e["z_box"] - zb) < 1e-12 and (
            abs(zb - 3.3415) < 5e-4 if bins == 12 else abs(zb - 2.5758) < 5e-4
        )
        assert len(e["per_bin_box_coverage_min_max"]) == bins and e["m_lo"] <= e["m"] <= e["m_hi"]
        assert abs(e["ci"][0] - (e["m_lo"] - e["t"] * e["s"])) < 1e-12
        assert abs(e["ci"][1] - (e["m_hi"] + e["t"] * e["s"])) < 1e-12
        assert abs(e["t"] - v5b.student_t_quantile(0.98, 5399)) < 1e-12
        assert np.allclose(e["r_b"], 1.0 / math.sqrt(5400))
        assert e["max_abs_z_ref1e6"] < 6.0 and 0.5 < e["m_ref1e6"] < 0.85
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
        v5b.v7_rep_combine(shards[:5], ref)
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
    e = v5b.v7_estimator_verdict("scalar", np.array([[1.0], [3.0], [2.0]]), np.array([[0.5]] * 3),
                                 np.array([2.0]), np.array([0.1]))  # fmt: skip
    assert e["per_replicate_covered"] == [0, 0, 1] and e["intervals"] == 3 and not e["pass"]


# -- Monte Carlo calibration of the full pipeline (marker ``calibration``: not in CI) ---------
def _kappa_for(v5b: ModuleType, coverage: float) -> float:
    """kappa with P(|t_19| <= kappa) = coverage (intervals s = kappa sigma sqrt(chi2_19 / 19))."""
    lo, hi = 0.1, 5.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if v5b.student_abs_prob(19, mid) < coverage else (lo, mid)
    return 0.5 * (lo + hi)


def _cp_upper(v5b: ModuleType, k: int, n: int, level: float = 0.95) -> float:
    """One-sided Clopper-Pearson upper confidence bound of a binomial rate (k of n)."""
    if k >= n:
        return 1.0
    lo, hi = k / n, 1.0
    for _ in range(80):  # P(X <= k | p) = 1 - I_p(k + 1, n - k) = 1 - level
        mid = 0.5 * (lo + hi)
        if v5b.betainc_reg(k + 1.0, n - k, mid) < level:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _box_pass(v5b: ModuleType, tq: float, mean: np.ndarray, sem: np.ndarray) -> bool:
    """Vectorised form of ``v7_estimator_verdict`` (all bins used): leave-one-out pooled reference of
    each replicate, per-bin worst case over the 9-point grid of the Bonferroni box (Z_B SEM_pool),
    PASS iff m_lo - t s >= 0.640 and m_hi + t s <= 0.700 with t = t_{0.98, R - 1}."""
    reps, bins = mean.shape
    sig = mean.std(axis=0, ddof=1)
    loo = (mean.sum(axis=0) - mean) / (reps - 1)
    d0 = mean - loo
    grid = np.linspace(-1.0, 1.0, v5b.V7_BOX_GRID) * v5b.normal_quantile(
        1.0 - v5b.V7_ALPHA_BOX / (2.0 * bins)
    )
    sp = sig / math.sqrt(reps)
    cov = np.stack([(np.abs(d0 - d * sp) <= sem).mean(axis=0) for d in grid])
    s = (np.abs(d0) <= sem).mean(axis=1).std(ddof=1) / math.sqrt(reps)
    m_lo, m_hi = cov.min(axis=0).mean(), cov.max(axis=0).mean()
    return bool(m_lo - tq * s >= v5b.V7_COV_LOW and m_hi + tq * s <= v5b.V7_COV_HIGH)


def _t98(v5b: ModuleType, reps: int = 5400) -> float:
    return v5b.student_t_quantile(1.0 - v5b.V7_ALPHA_TOST / 2.0, reps - 1)


def _gauss_passes(v5b: ModuleType, kappas: tuple[float, ...], bins: int, trials: int,
                  rng: np.random.Generator, rho: float = 0.0, reps: int = 5400) -> np.ndarray:  # fmt: skip
    """Pass indicators ``[len(kappas), trials]``. Replicate means X_b = sigma_b (sqrt(rho) g + sqrt(1 -
    rho) e_b) with a replicate-level common factor g (the bins of one replicate share histories: bin
    correlation ``rho``), batch SEMs s = kappa sigma_b sqrt(chi2_19 / 19); true mean 0, reference = the
    pooled leave-one-out mean. The same draws serve all kappas. First trial checked against the pure
    ``v7_estimator_verdict``."""
    sigma = np.linspace(0.5, 2.0, bins) if bins > 1 else np.array([1.3])
    tq = _t98(v5b, reps)
    out = np.zeros((len(kappas), trials), dtype=bool)
    for i in range(trials):
        g = rng.standard_normal((reps, 1))
        x = sigma * (math.sqrt(rho) * g + math.sqrt(1.0 - rho) * rng.standard_normal((reps, bins)))
        s0 = sigma * np.sqrt(rng.chisquare(19, (reps, bins)) / 19.0)
        for j, kappa in enumerate(kappas):
            out[j, i] = _box_pass(v5b, tq, 100.0 + x, kappa * s0)
        if i == 0:
            e = v5b.v7_estimator_verdict("scalar", 100.0 + x, kappas[0] * s0, np.full(bins, 100.0),
                                         np.full(bins, 0.01))  # fmt: skip
            assert e["pass"] == out[0, 0]
    return out


def _rate_report(v5b: ModuleType, passes: np.ndarray) -> list[dict[str, float]]:
    """Rate and 95 % Clopper-Pearson upper bound per row of a ``[kappas, trials]`` pass array."""
    n = passes.shape[1]
    return [
        {"rate": float(row.mean()), "cp95_upper": _cp_upper(v5b, int(row.sum()), n), "n": n}
        for row in passes
    ]


GAUSS_CASES = (
    ("profile12_rho0.0", 12, 3000, 0.0),
    ("profile12_rho0.8", 12, 3000, 0.8),
    ("scalar", 1, 5000, 0.0),
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
    """Gaussian pipeline, R = 5400, pooled leave-one-out reference, 12 bins with replicate-level bin
    correlation 0 and 0.8 and the scalar: the 95 % Clopper-Pearson upper bound of the false-acceptance
    rate is <= 0.05 at true coverage 0.640 and at 0.700."""
    rep = _rate_report(v5b, gauss_rates[label])
    assert rep[0]["cp95_upper"] <= 0.05 and rep[2]["cp95_upper"] <= 0.05, (label, rep)


@pytest.mark.calibration
def test_v7_pipeline_gaussian_power_targets(gauss_rates: dict) -> None:
    """Power at the nominal coverage (kappa = 1, true coverage 0.670): profiles >= 0.95, scalar >= 0.90,
    joint row (independent draws of the two profiles and the scalar) >= 0.85."""
    power = {k: float(v[1].mean()) for k, v in gauss_rates.items()}
    n = 3000
    joint = float((gauss_rates["profile12_rho0.0"][1][:n] & gauss_rates["profile12_rho0.8"][1][:n]
                   & gauss_rates["scalar"][1][:n]).mean())  # fmt: skip
    print("v7 gaussian power at nominal", power, "joint row (independent draws)", joint)
    assert power["profile12_rho0.0"] >= 0.95 and power["profile12_rho0.8"] >= 0.95
    assert power["scalar"] >= 0.90 and joint >= 0.85, (power, joint)


W_SHARED = 0.5  # shape of the shared per-block Gamma factor of the correlated sparse bins (mean 1)


def _gamma_blocks(
    shape: float, size: tuple[int, ...], rng: np.random.Generator, correlated: bool
) -> np.ndarray:
    """Block values ``[..., bins]`` with mean ``shape``. Independent: Gamma(shape) per bin. Correlated
    sparse bins as from shared histories: the same Gamma draws times ONE factor W_b ~ Gamma(0.5) / 0.5
    (mean 1) per block shared by all bins, so that a block with few histories is sparse in every bin."""
    g = rng.standard_gamma(shape, size, dtype=np.float32)
    if correlated:
        g = g * (rng.standard_gamma(W_SHARED, (*size[:-1], 1), dtype=np.float32) / W_SHARED)
    return g


def _skew_kappas(
    shape: float, correlated: bool, rng: np.random.Generator
) -> tuple[float, float, float]:
    """kappa giving TRUE coverage 0.640 and 0.700 (large-sample quantiles of |mean - mu| / SEM of the
    marginal bin, the bisection with common random numbers) and the true coverage at kappa = 1."""
    g = _gamma_blocks(shape, (400_000, 20, 1), rng, correlated).astype(float)[..., 0]
    ratios = np.abs(g.mean(axis=1) - shape) / (g.std(axis=1, ddof=1) / math.sqrt(20))
    k_lo, k_hi = (float(np.quantile(ratios, c)) for c in (0.640, 0.700))
    assert (
        abs(np.mean(ratios <= k_lo) - 0.640) < 1e-4 and abs(np.mean(ratios <= k_hi) - 0.700) < 1e-4
    )
    return k_lo, k_hi, float(np.mean(ratios <= 1.0))


def _skew_passes(v5b: ModuleType, shape: float, correlated: bool, kappas: tuple[float, ...], trials: int,
                 rng: np.random.Generator, bins: int = 12, reps: int = 5400) -> np.ndarray:  # fmt: skip
    """Pass indicators ``[len(kappas), trials]`` on Gamma-block tallies (20 blocks per replicate and bin,
    replicate mean and SEM from the blocks, half width kappa * SEM, per-bin scales 0.5..2) against the
    pooled leave-one-out reference; ``correlated`` bins as in ``_gamma_blocks``."""
    scale = np.linspace(0.5, 2.0, bins)
    tq = _t98(v5b, reps)
    out = np.zeros((len(kappas), trials), dtype=bool)
    for i in range(trials):
        g = _gamma_blocks(shape, (reps, 20, bins), rng, correlated) * scale
        mean, sem = g.mean(axis=1), g.std(axis=1, ddof=1) / math.sqrt(20)
        for j, kappa in enumerate(kappas):
            out[j, i] = _box_pass(v5b, tq, mean.astype(float), (kappa * sem).astype(float))
        if i == 0:
            e = v5b.v7_estimator_verdict("profile", mean.astype(float), (kappas[0] * sem).astype(float),
                                         shape * scale, 0.01 * scale)  # fmt: skip
            assert e["pass"] == out[0, 0]
    return out


SKEW_CASES = [(sh, corr) for sh in (2.0, 0.5) for corr in (False, True)]


@pytest.fixture(
    scope="module",
    params=SKEW_CASES,
    ids=[f"gamma{sh}-{'correlated' if c else 'independent'}" for sh, c in SKEW_CASES],
)
def skew_result(request: pytest.FixtureRequest, v5b: ModuleType) -> dict:
    """Gamma shape 2 and 0.5, independent and correlated sparse bins, 12 bins, R = 5400: kappa scales
    the SEM so that the TRUE coverage is 0.640 and 0.700; pass rates with CP bounds there, power at
    kappa = 1 and the true coverage at kappa = 1."""
    shape, corr = request.param
    rng = np.random.default_rng(20471005 + int(2 * shape) + (7 if corr else 0))
    k_lo, k_hi, true_nominal = _skew_kappas(shape, corr, rng)
    rep = _rate_report(v5b, _skew_passes(v5b, shape, corr, (k_lo, 1.0, k_hi), 2000, rng))
    out = {"shape": shape, "correlated": corr, "kappa_0.640": k_lo, "kappa_0.700": k_hi,
           "true_coverage_at_kappa1": true_nominal, "FA_0.640": rep[0], "power_kappa1": rep[1], "FA_0.700": rep[2]}  # fmt: skip
    print("v7 skewed-tally calibration", out)
    return out


@pytest.mark.calibration
def test_v7_pipeline_skewed_tally_false_acceptance(skew_result: dict) -> None:
    """The 95 % Clopper-Pearson upper bound of the false-acceptance rate is <= 0.05 at true coverage
    0.640 and 0.700 for Gamma-block tallies, independent and correlated (shared-history) bins."""
    assert (
        skew_result["FA_0.640"]["cp95_upper"] <= 0.05
        and skew_result["FA_0.700"]["cp95_upper"] <= 0.05
    ), skew_result


def test_calibration_marker_registered_and_excluded_from_ci() -> None:
    """The Monte Carlo calibration is out of the lightweight CI selection and registered as a marker."""
    ci = (REPO / ".github" / "workflows" / "ci.yml").read_text()
    assert '-m "not cuda and not host and not calibration"' in ci
    assert '"calibration:' in (REPO / "pyproject.toml").read_text()
