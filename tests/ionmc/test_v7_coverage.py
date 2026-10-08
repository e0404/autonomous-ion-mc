# ruff: noqa: E501
"""V7 replicate-coverage rule (plan Amendment 13, Codex REVIEW-be30e621, -873ab9cd, -edb9970b, -f1b32614):
two-fold cross reference (conditionally independent replicates), Bonferroni box over both fold references with the
exact worst-case coverage by an endpoint sweep, Hoeffding deviation bound at alpha_tost = 0.04, combine of
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
    assert v5b.V7_REPLICATES == 7200

    def verdict(cov: np.ndarray, m_lo=None, m_hi=None, **kw: object) -> dict:  # type: ignore[no-untyped-def]
        return v5b.v7_tost_verdict(cov, 12, kw.pop("intervals", n), m_lo, m_hi, min_bins=10, **kw)

    ok = verdict(_synthetic_covered(0.670, 0.135, r=7200))
    eps = ok["hoeffding_eps"]
    assert ok["pass"] and not ok["reasons"] and abs(ok["m"] - 0.670) < 1e-9
    assert abs(ok["s"] - 0.135 / math.sqrt(7200)) < 1e-9 and ok["replicates"] == 7200
    assert abs(eps - math.sqrt(math.log(25.0) / (2 * 7200))) < 1e-15 and abs(eps - 0.014951) < 5e-5
    assert ok["m_lo"] == ok["m_hi"] == ok["m"] and "t" not in ok  # no box: collapses to the point
    assert abs(ok["ci"][0] - (0.670 - eps)) < 1e-9 and abs(ok["ci"][1] - (0.670 + eps)) < 1e-9
    low = verdict(_synthetic_covered(0.653, 0.135, r=7200))  # lower bound 0.638
    assert 0.640 < low["m"] < 0.700 and low["ci"][0] < 0.640 and not low["pass"]
    assert any("lower" in x for x in low["reasons"])
    up = verdict(_synthetic_covered(0.688, 0.135, r=7200))  # upper bound 0.703
    assert up["ci"][1] > 0.700 and not up["pass"] and any("upper" in x for x in up["reasons"])
    assert verdict(_synthetic_covered(0.656, 0.135, r=7200))["pass"]
    # the box widens both sides: m_lo / m_hi enter the interval
    cov = _synthetic_covered(0.670, 0.135, r=7200)
    band = verdict(cov, 0.660, 0.680)
    assert (band["m_lo"], band["m_hi"]) == (0.660, 0.680) and band["pass"]
    assert abs(band["ci"][0] - (0.660 - eps)) < 1e-12 and abs(band["ci"][1] - (0.680 + eps)) < 1e-12
    lo_fail = verdict(cov, 0.640 + 0.5 * eps, 0.670)
    assert not lo_fail["pass"] and any("lower" in x for x in lo_fail["reasons"])
    hi_fail = verdict(cov, 0.670, 0.700 - 0.5 * eps)
    assert not hi_fail["pass"] and any("upper" in x for x in hi_fail["reasons"])
    br = verdict(cov, 0.69, 0.60)  # inconsistent inputs are clipped so that m_lo <= m <= m_hi
    assert br["m_lo"] == br["m_hi"] == br["m"]
    assert not verdict(cov, intervals=299)["pass"]
    assert not v5b.v7_tost_verdict(_synthetic_covered(0.670, 0.135, 9), 9, n, min_bins=10)["pass"]


def test_box_quantile_hoeffding_and_z_values(v5b: ModuleType) -> None:
    """Z_B = Phi^{-1}(1 - alpha_box / (2 * 2B)): 3.53 for 12 bins, 2.81 for one bin; eps = 0.01495."""
    assert (
        abs(v5b.normal_quantile(0.975) - 1.959964) < 1e-5 and abs(v5b.normal_quantile(0.5)) < 1e-12
    )
    z12 = v5b.normal_quantile(1.0 - 0.01 / 48.0)
    z1 = v5b.normal_quantile(1.0 - 0.01 / 4.0)
    assert abs(z12 - 3.5293) < 2e-3 and abs(z1 - 2.8070) < 2e-3
    assert abs(v5b.v7_hoeffding_eps(7200) - 0.014951) < 1e-5
    assert abs(v5b.v7_hoeffding_eps(7200) - math.sqrt(math.log(1.0 / 0.04) / 14400)) < 1e-15


def test_coverage_sweep_matches_brute_force(v5b: ModuleType) -> None:
    """The exact sweep equals a brute-force evaluation of C(delta) at all clipped breakpoints, at
    +-Z and at the midpoints between them (a fine grid cannot beat it), on random small inputs."""
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
        mn, mx = v5b.v7_coverage_extrema(x, sem, z)
        assert mn == min(brute) and mx == max(brute)
        fine = [cov(d) for d in np.linspace(-z, z, 4001)]  # a fine grid never goes beyond the sweep
        assert mn <= min(fine) + 1e-12 and mx >= max(fine) - 1e-12
    # closed intervals: a single interval [0, 1], box [-2, 2]: max 1 at 0..1, min 0 outside
    assert v5b.v7_coverage_extrema(np.array([0.5]), np.array([0.5]), 2.0) == (0.0, 1.0)
    assert v5b.v7_coverage_extrema(np.array([0.5]), np.array([0.5]), 0.25) == (
        0.0,
        1.0,
    )  # box [-.25,.25]: lo = 0 inside
    assert v5b.v7_coverage_extrema(np.array([0.0]), np.array([5.0]), 2.0) == (
        1.0,
        1.0,
    )  # covers the box


def test_fold_reference_never_contains_own_replicate(v5b: ModuleType) -> None:
    """Replicates of fold A are judged against the mean of fold B and vice versa: changing a replicate's
    own value changes its own hit only, and never the reference of its own fold."""
    rng = np.random.default_rng(2)
    x = 5.0 + rng.normal(0.0, 1.0, (8, 2))
    sem = np.full_like(x, 0.8)
    ref, rsem = np.array([5.0, 5.0]), np.array([0.1, 0.1])
    e0 = v5b.v7_estimator_verdict("profile", x, sem, ref, rsem, None)
    mu_a, mu_b = x[:4].mean(axis=0), x[4:].mean(axis=0)
    assert np.allclose(e0["fold_means"]["A"], mu_a) and np.allclose(e0["fold_means"]["B"], mu_b)
    assert np.allclose(e0["fold_sems"]["A"], x[:4].std(axis=0, ddof=1) / 2.0)
    hits = np.concatenate([np.abs(x[:4] - mu_b) <= 0.8, np.abs(x[4:] - mu_a) <= 0.8])
    assert e0["per_replicate_covered"] == hits.sum(axis=1).tolist()
    x2 = x.copy()
    x2[1, 0] += 0.3  # replicate 1 (fold A): the reference of fold A is the fold-B mean, unchanged
    e1 = v5b.v7_estimator_verdict("profile", x2, sem, ref, rsem, None)
    assert e1["fold_means"]["B"] == e0["fold_means"]["B"]  # fold A's reference is untouched
    assert (
        abs(e1["fold_means"]["A"][0] - (mu_a[0] + 0.075)) < 1e-12
    )  # only the other fold's reference moves
    sd_a = x[:4].std(axis=0, ddof=1)
    z_f = (mu_a - mu_b) / np.sqrt(sd_a**2 / 4 + x[4:].std(axis=0, ddof=1) ** 2 / 4)
    assert (
        np.allclose(e0["z_b_folds"], z_f)
        and abs(e0["max_abs_z_folds"] - np.max(np.abs(z_f))) < 1e-12
    )


def test_estimator_verdict_box_and_cross_check_fields(v5b: ModuleType) -> None:
    """Hand-checkable scalar case, 4 replicates (folds of two): Z_B = 2.807, exact per-fold extrema."""
    means = np.array([[1.0], [1.2], [0.8], [1.0]])
    sems = np.array([[0.15]] * 4)
    ref, rsem = np.array([1.0]), np.array([0.1])
    e = v5b.v7_estimator_verdict("scalar", means, sems, ref, rsem, None)
    mu_a, mu_b = 1.1, 0.9  # fold A = rows 0-1, fold B = rows 2-3
    assert e["per_replicate_covered"] == [1, 0, 0, 1]  # |1.0-.9|, |1.2-.9|, |.8-1.1|, |1.0-1.1|
    zb = v5b.normal_quantile(1.0 - 0.01 / 4.0)
    assert abs(e["z_box"] - zb) < 1e-12 and e["alpha_box"] == 0.01
    sem_a = float(np.std([1.0, 1.2], ddof=1)) / math.sqrt(2.0)
    sem_b = float(np.std([0.8, 1.0], ddof=1)) / math.sqrt(2.0)
    lo_a, hi_a = (
        v5b.v7_coverage_extrema(np.array([1.0, 1.2]) - mu_b, np.array([0.15, 0.15]), 0.0)
        if False
        else (None, None)
    )
    ma = v5b.v7_coverage_extrema(
        (np.array([1.0, 1.2]) - mu_b) / sem_b, np.array([0.15, 0.15]) / sem_b, zb
    )
    mb = v5b.v7_coverage_extrema(
        (np.array([0.8, 1.0]) - mu_a) / sem_a, np.array([0.15, 0.15]) / sem_a, zb
    )
    assert e["per_fold_bin_box_coverage_min_max"] == {"A": [list(ma)], "B": [list(mb)]}
    assert (
        abs(e["m_lo"] - 0.5 * (ma[0] + mb[0])) < 1e-15
        and abs(e["m_hi"] - 0.5 * (ma[1] + mb[1])) < 1e-15
    )
    assert e["m_lo"] <= e["m"] <= e["m_hi"] and not e["pass"]
    assert abs(e["hoeffding_eps"] - v5b.v7_hoeffding_eps(4)) < 1e-15
    assert abs(e["ci"][0] - (e["m_lo"] - e["hoeffding_eps"])) < 1e-15
    # 1e6 cross-check: hits |x - 1.0| <= 0.15 -> 1, 0, 0, 1
    assert e["m_ref1e6"] == 0.5
    z = (means.mean() - 1.0) / math.hypot(float(means.std(ddof=1)) / 2.0, 0.1)
    assert abs(e["z_b_ref1e6"][0] - z) < 1e-15 and math.isnan(e["s_ref_cov"])
    # profile: m_lo / m_hi are the (fold-weighted) means over bins of the per-(fold, bin) minima / maxima
    x = np.array([[1.0, 5.0], [1.2, 4.0], [0.8, 6.0], [1.0, 5.5], [0.9, 5.2], [1.1, 4.7]])
    p = v5b.v7_estimator_verdict(
        "profile", x, np.full_like(x, 0.4), np.array([1.0, 5.0]), np.array([0.1, 0.1])
    )
    bx = p["per_fold_bin_box_coverage_min_max"]
    lo = 0.5 * (np.mean([b[0] for b in bx["A"]]) + np.mean([b[0] for b in bx["B"]]))
    assert (
        abs(p["m_lo"] - lo) < 1e-15
        and abs(p["z_box"] - v5b.normal_quantile(1.0 - 0.01 / 8.0)) < 1e-12
    )


def _binary_hits(frac: float, r: int = 7200) -> np.ndarray:
    h = np.zeros(r, dtype=np.int64)
    h[: int(round(frac * r))] = 1
    return h


def test_scalar_estimator_with_binary_hits(v5b: ModuleType) -> None:
    ok = v5b.v7_tost_verdict(_binary_hits(0.670), 1, 7200)
    assert ok["pass"] and abs(ok["m"] - 0.670) < 1e-3  # sd 0.47, s = 0.0064, CI +- 0.0105
    assert abs(ok["s"] - math.sqrt(0.67 * 0.33 * 7200 / 7199) / math.sqrt(7200)) < 1e-4
    edge = v5b.v7_tost_verdict(_binary_hits(0.650), 1, 7200)
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
    assert len(shards) == 8
    doc = v5b.v7_rep_combine(shards, ref)
    assert doc["replicates"] == 7200 and set(doc["estimators"]) == {
        "sec_p",
        "nuclear_local",
        "escaped_neutral",
    }
    assert (
        doc["rule"]["alpha_box"] == 0.01
        and doc["rule"]["alpha_tost"] == 0.04
        and abs(doc["rule"]["hoeffding_eps"] - 0.014951) < 1e-5
    )
    for name, kind, bins in (
        ("sec_p", "profile", 12),
        ("nuclear_local", "profile", 12),
        ("escaped_neutral", "scalar", 1),
    ):
        e = doc["estimators"][name]
        assert e["kind"] == kind and e["bins_used"] == bins and e["intervals"] == 7200 * bins
        assert len(e["per_replicate_covered"]) == 7200 and all(
            isinstance(c, int) for c in e["per_replicate_covered"]
        )
        assert abs(e["nominal_coverage"] - 0.670) < 0.002
        for key in (
            *("m", "m_lo", "m_hi", "z_box", "alpha_box", "alpha_tost", "hoeffding_eps", "s", "ci"),
            *(
                "per_fold_bin_box_coverage_min_max",
                "fold_means",
                "fold_sems",
                "z_b_folds",
                "max_abs_z_folds",
            ),
            *("s_ref_cov", "ref_corr_mean_abs_offdiag", "replicate_mean_skewness", "reasons"),
            *("legacy_point_gate_pass", "pass", "m_ref1e6", "z_b_ref1e6", "max_abs_z_ref1e6"),
        ):
            assert key in e
        assert not {"m_corr", "s_tot", "s_ref_bound", "t", "box_grid", "r_b"} & set(e)
        assert "two-fold" in e["reference_kind"]
        zb = v5b.normal_quantile(1.0 - 0.01 / (4.0 * bins))
        assert abs(e["z_box"] - zb) < 1e-12 and (
            abs(zb - 3.5293) < 2e-3 if bins == 12 else abs(zb - 2.8070) < 2e-3
        )
        assert (
            len(e["per_fold_bin_box_coverage_min_max"]["A"]) == bins
            and e["m_lo"] <= e["m"] <= e["m_hi"]
        )
        assert abs(e["ci"][0] - (e["m_lo"] - e["hoeffding_eps"])) < 1e-12
        assert abs(e["ci"][1] - (e["m_hi"] + e["hoeffding_eps"])) < 1e-12
        assert abs(e["hoeffding_eps"] - 0.014951) < 1e-5
        assert (
            e["max_abs_z_ref1e6"] < 6.0
            and e["max_abs_z_folds"] < 6.0
            and 0.5 < e["m_ref1e6"] < 0.85
        )
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
    e = v5b.v7_estimator_verdict("scalar", np.array([[1.0], [3.0], [2.0], [2.2]]), np.array([[0.5]] * 4),
                                 np.array([2.0]), np.array([0.1]))  # fmt: skip
    assert e["per_replicate_covered"] == [0, 0, 1, 1] and e["intervals"] == 4 and not e["pass"]
    with pytest.raises(SystemExit, match="two equal folds"):
        v5b.v7_estimator_verdict("scalar", np.array([[1.0], [3.0], [2.0]]), np.array([[0.5]] * 3),
                                 np.array([2.0]), np.array([0.1]))  # fmt: skip


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


def _box_pass(v5b: ModuleType, mean: np.ndarray, sem: np.ndarray) -> bool:
    """Vectorised form of ``v7_estimator_verdict`` (all bins used): the replicates of fold A (first half)
    are judged against the fold-B mean and vice versa, exact worst-case coverage per (fold, bin) over the
    box Z_B = Phi^{-1}(1 - alpha_box / (4 B)) of the reference (``v7_coverage_extrema``), PASS iff m_lo -
    eps >= 0.640 and m_hi + eps <= 0.700 with the Hoeffding eps of R replicates at alpha_tost."""
    reps, bins = mean.shape
    na = reps // 2
    z = v5b.normal_quantile(1.0 - v5b.V7_ALPHA_BOX / (4.0 * bins))
    folds = (slice(0, na), slice(na, reps))
    mu = [mean[sl].mean(axis=0) for sl in folds]
    se = [mean[sl].std(axis=0, ddof=1) / math.sqrt(mean[sl].shape[0]) for sl in folds]
    f = np.zeros(reps)
    lo = hi = 0.0
    for k in (0, 1):
        sl, o = folds[k], 1 - k
        x = mean[sl] - mu[o]
        f[sl] = (np.abs(x) <= sem[sl]).sum(axis=1)
        ext = [
            v5b.v7_coverage_extrema(x[:, b] / se[o][b], sem[sl][:, b] / se[o][b], z)
            for b in range(bins)
        ]
        lo += mean[sl].shape[0] / reps * float(np.mean([e[0] for e in ext]))
        hi += mean[sl].shape[0] / reps * float(np.mean([e[1] for e in ext]))
    f /= bins
    m = f.mean()
    eps = v5b.v7_hoeffding_eps(reps)
    return bool(min(lo, m) - eps >= v5b.V7_COV_LOW and max(hi, m) + eps <= v5b.V7_COV_HIGH)


def _gauss_passes(v5b: ModuleType, kappas: tuple[float, ...], bins: int, trials: int,
                  rng: np.random.Generator, rho: float = 0.0, reps: int = 7200) -> np.ndarray:  # fmt: skip
    """Pass indicators ``[len(kappas), trials]``. Replicate means X_b = sigma_b (sqrt(rho) g + sqrt(1 -
    rho) e_b) with a replicate-level common factor g (the bins of one replicate share histories: bin
    correlation ``rho``), batch SEMs s = kappa sigma_b sqrt(chi2_19 / 19); true mean 0, reference = the
    two-fold cross mean. The same draws serve all kappas. First trial checked against the pure
    ``v7_estimator_verdict``."""
    sigma = np.linspace(0.5, 2.0, bins) if bins > 1 else np.array([1.3])
    out = np.zeros((len(kappas), trials), dtype=bool)
    for i in range(trials):
        g = rng.standard_normal((reps, 1))
        x = sigma * (math.sqrt(rho) * g + math.sqrt(1.0 - rho) * rng.standard_normal((reps, bins)))
        s0 = sigma * np.sqrt(rng.chisquare(19, (reps, bins)) / 19.0)
        for j, kappa in enumerate(kappas):
            out[j, i] = _box_pass(v5b, 100.0 + x, kappa * s0)
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
    """Gaussian pipeline, R = 7200, two-fold cross reference, 12 bins with replicate-level bin
    correlation 0 and 0.8 and the scalar: the 95 % Clopper-Pearson upper bound of the false-acceptance
    rate is <= 0.05 at true coverage 0.640 and at 0.700."""
    rep = _rate_report(v5b, gauss_rates[label])
    assert rep[0]["cp95_upper"] <= 0.05 and rep[2]["cp95_upper"] <= 0.05, (label, rep)


@pytest.mark.calibration
def test_v7_pipeline_gaussian_power_profiles(gauss_rates: dict) -> None:
    """Power at the nominal coverage (kappa = 1, true coverage 0.670): profiles >= 0.95."""
    power = {k: float(v[1].mean()) for k, v in gauss_rates.items()}
    print("v7 gaussian power at nominal", power)
    assert power["profile12_rho0.0"] >= 0.95 and power["profile12_rho0.8"] >= 0.95, power


@pytest.mark.calibration
def test_v7_pipeline_gaussian_power_scalar_and_joint_targets(gauss_rates: dict) -> None:
    n = 1500
    scalar = float(gauss_rates["scalar"][1].mean())
    joint = float((gauss_rates["profile12_rho0.0"][1][:n] & gauss_rates["profile12_rho0.8"][1][:n]
                   & gauss_rates["scalar"][1][:n]).mean())  # fmt: skip
    print("v7 gaussian scalar power", scalar, "joint row (independent draws)", joint)
    assert scalar >= 0.90 and joint >= 0.85, (scalar, joint)


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
                 rng: np.random.Generator, bins: int = 12, reps: int = 7200) -> np.ndarray:  # fmt: skip
    """Pass indicators ``[len(kappas), trials]`` on Gamma-block tallies (20 blocks per replicate and bin,
    replicate mean and SEM from the blocks, half width kappa * SEM, per-bin scales 0.5..2) against the
    two-fold cross reference; ``correlated`` bins as in ``_gamma_blocks``."""
    scale = np.linspace(0.5, 2.0, bins)
    out = np.zeros((len(kappas), trials), dtype=bool)
    for i in range(trials):
        g = _gamma_blocks(shape, (reps, 20, bins), rng, correlated) * scale
        mean, sem = g.mean(axis=1), g.std(axis=1, ddof=1) / math.sqrt(20)
        for j, kappa in enumerate(kappas):
            out[j, i] = _box_pass(v5b, mean.astype(float), (kappa * sem).astype(float))
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
    """Gamma shape 2 and 0.5, independent and correlated sparse bins, 12 bins, R = 7200: kappa scales
    the SEM so that the TRUE coverage is 0.640 and 0.700; pass rates with CP bounds there, power at
    kappa = 1 and the true coverage at kappa = 1."""
    shape, corr = request.param
    rng = np.random.default_rng(20471005 + int(2 * shape) + (7 if corr else 0))
    k_lo, k_hi, true_nominal = _skew_kappas(shape, corr, rng)
    rep = _rate_report(v5b, _skew_passes(v5b, shape, corr, (k_lo, 1.0, k_hi), 800, rng))
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
