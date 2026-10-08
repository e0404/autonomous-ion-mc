"""V7 replicate-coverage rule (plan Amendment 13, Codex REVIEW-be30e621): reference-uncertainty
correction, TOST per estimator, combine of shard/reference partials (fail closed) and a Monte Carlo
calibration of the whole pipeline including the shared reference error."""

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
    assert abs(c["s_ref"] - PHI1 * math.sqrt(2.0) * math.sqrt(np.sum(r**4)) / 4) < 1e-15
    one = v5b.v7_reference_correction(np.array([0.1]), np.array([1.0]))  # scalar: B = 1
    assert (
        abs(one["b_ref"] - PHI1 * 0.01) < 1e-15
        and abs(one["s_ref"] - PHI1 * math.sqrt(2.0) * 1e-2) < 1e-15
    )


def test_tost_rule_synthetic_vectors(v5b: ModuleType) -> None:
    assert (v5b.V7_COV_LOW, v5b.V7_COV_HIGH, v5b.V7_TOST_ALPHA) == (0.640, 0.700, 0.05)
    n = 4500 * 12

    def verdict(cov: np.ndarray, b: float = 0.0, s: float = 0.0, **kw: object) -> dict:
        return v5b.v7_tost_verdict(cov, 12, kw.pop("intervals", n), b, s, min_bins=10, **kw)

    ok = verdict(_synthetic_covered(0.670, 0.135))
    assert ok["pass"] and not ok["reasons"] and abs(ok["m"] - 0.670) < 1e-9
    assert abs(ok["s"] - 0.135 / math.sqrt(4500)) < 1e-9 and ok["replicates"] == 4500
    assert abs(ok["t"] - v5b.student_t_quantile(0.95, 4499)) < 1e-12
    assert abs(ok["ci"][1] - ok["ci"][0] - 2 * ok["t"] * ok["s_tot"]) < 1e-12
    low = verdict(_synthetic_covered(0.6415, 0.135))  # point inside, CI crosses 0.640
    assert 0.640 < low["m"] < 0.700 and low["ci"][0] < 0.640 and not low["pass"]
    assert any("lower" in x for x in low["reasons"])
    up = verdict(_synthetic_covered(0.6985, 0.135))
    assert up["ci"][1] > 0.700 and not up["pass"] and any("upper" in x for x in up["reasons"])
    # bias correction moves the interval up: 0.6415 + 0.005 passes; s_ref widens it
    fixed = verdict(_synthetic_covered(0.6415, 0.135), b=0.005)
    assert fixed["pass"] and abs(fixed["m_corr"] - 0.6465) < 1e-12
    wide = verdict(_synthetic_covered(0.670, 0.135), s=0.02)
    assert abs(wide["s_tot"] - math.hypot(wide["s"], 0.02)) < 1e-15 and not wide["pass"]
    assert not verdict(_synthetic_covered(0.670, 0.135), intervals=299)["pass"]
    assert not v5b.v7_tost_verdict(
        _synthetic_covered(0.670, 0.135, 9), 9, n, 0.0, 0.0, min_bins=10
    )["pass"]


def _binary_hits(frac: float, r: int = 4500) -> np.ndarray:
    h = np.zeros(r, dtype=np.int64)
    h[: int(round(frac * r))] = 1
    return h


def test_scalar_estimator_with_binary_hits(v5b: ModuleType) -> None:
    ok = v5b.v7_tost_verdict(_binary_hits(0.670), 1, 4500, 0.0, 0.0)
    assert ok["pass"] and abs(ok["m"] - 0.670) < 1e-3  # sd 0.47, s = 0.0090, CI +- 0.0149
    assert abs(ok["s"] - math.sqrt(0.67 * 0.33 * 4500 / 4499) / math.sqrt(4500)) < 1e-4
    edge = v5b.v7_tost_verdict(_binary_hits(0.650), 1, 4500, 0.0, 0.0)
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
        ref_est[nm] = {"mean": mu[nm].tolist(), "sem": (0.1 * sig[nm]).tolist()}
    ref = {
        "seed": v5b.seed_of("v7-rep", v5b.V7_SHARDS),
        "valid": True,
        "n": 1_000_000,
        "estimators": ref_est,
    }
    return shards, ref


def test_combine_document_and_fail_closed(v5b: ModuleType) -> None:
    shards, ref = _synthetic_partials(v5b)
    doc = v5b.v7_rep_combine(shards, ref)
    assert doc["replicates"] == 4500 and set(doc["estimators"]) == {
        "sec_p",
        "nuclear_local",
        "escaped_neutral",
    }
    for name, kind, bins in (
        ("sec_p", "profile", 12),
        ("nuclear_local", "profile", 12),
        ("escaped_neutral", "scalar", 1),
    ):
        e = doc["estimators"][name]
        assert e["kind"] == kind and e["bins_used"] == bins and e["intervals"] == 4500 * bins
        assert len(e["per_replicate_covered"]) == 4500 and all(
            isinstance(c, int) for c in e["per_replicate_covered"]
        )
        assert len(e["r_b"]) == bins and abs(e["nominal_coverage"] - 0.670) < 0.002
        for key in ("m", "s", "b_ref", "s_ref", "ci", "reasons", "legacy_point_gate_pass", "pass"):
            assert key in e
        assert abs(e["m"] + e["b_ref"] - e["m_corr"]) < 1e-15
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
        v5b.v7_rep_combine(shards[:2], ref)
    with pytest.raises(SystemExit):
        v5b.v7_rep_combine(shards, {**ref, "valid": False})
    short = copy.deepcopy(ref)
    short["estimators"]["sec_p"]["mean"] = short["estimators"]["sec_p"]["mean"][:5]
    with pytest.raises(SystemExit, match="inconsistent"):
        v5b.v7_rep_combine(shards, short)


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


# -- Monte Carlo calibration of the full pipeline, including the shared reference ------------
def _kappa_for(v5b: ModuleType, coverage: float) -> float:
    """kappa with P(|t_19| <= kappa) = coverage (intervals s = kappa sigma sqrt(chi2_19 / 19))."""
    lo, hi = 0.1, 5.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if v5b.student_abs_prob(19, mid) < coverage else (lo, mid)
    return 0.5 * (lo + hi)


def _trial_passes(v5b: ModuleType, kappa: float, bins: int, trials: int, rng: np.random.Generator,
                 r_ref: float = 0.1, reps: int = 4500) -> np.ndarray:  # fmt: skip
    """Pass rate of the Amendment-13 rule. Per trial ONE reference error eps_b ~ N(0, r sigma_b) and
    one estimated reference SEM are drawn and shared by all replicates; means X ~ N(0, sigma_b),
    batch SEMs s = kappa sigma sqrt(chi2_19 / 19), hit = |X - eps| <= s. The rule is the vectorised
    form of ``v7_estimator_verdict`` (checked against it on the first trials)."""
    sigma = np.linspace(0.5, 2.0, bins) if bins > 1 else np.array([1.3])
    t = v5b.student_t_quantile(0.95, reps - 1)
    passes = np.zeros(trials, dtype=bool)
    for i in range(trials):
        eps = r_ref * sigma * rng.standard_normal(bins)
        sem_ref = r_ref * sigma * np.sqrt(rng.chisquare(19, bins) / 19.0)
        x = sigma * rng.standard_normal((reps, bins))
        s = kappa * sigma * np.sqrt(rng.chisquare(19, (reps, bins)) / 19.0)
        hit = np.abs(x - eps) <= s
        f = hit.mean(axis=1)
        sd = x.std(axis=0, ddof=1)
        rb = sem_ref / sd
        b_ref = PHI1 * np.mean(rb**2)
        s_ref = PHI1 * math.sqrt(2.0) * math.sqrt(np.sum(rb**4)) / bins
        m, se = f.mean() + b_ref, f.std(ddof=1) / math.sqrt(reps)
        half = t * math.hypot(se, s_ref)
        ok = bool(m - half >= v5b.V7_COV_LOW and m + half <= v5b.V7_COV_HIGH)
        if i < 3:  # the pure functions give the same verdict on the same inputs
            e = v5b.v7_estimator_verdict("scalar", 100.0 + x, s, 100.0 + eps, sem_ref)
            assert e["pass"] == ok and abs(e["m_corr"] - (f.mean() + b_ref)) < 1e-12
        passes[i] = ok
    return passes


def test_v7_pipeline_monte_carlo_calibration(v5b: ModuleType) -> None:
    """False acceptance <= 0.07 at both margins and power >= 0.95 at the nominal coverage, for the
    12-bin profile and the scalar estimator (R = 4500, shared reference error r = 0.1), and the
    joint row pass probability (two profiles, one scalar, independent draws) at the nominal."""
    rng = np.random.default_rng(20471004)
    k_lo, k_hi = _kappa_for(v5b, 0.640), _kappa_for(v5b, 0.700)
    assert abs(v5b.student_abs_prob(19, 1.0) - 0.670) < 0.002
    rates = {}
    for label, bins, trials in (("profile12", 12, 2000), ("scalar", 1, 20000)):
        rates[label] = {
            "0.640": float(_trial_passes(v5b, k_lo, bins, trials, rng).mean()),
            "0.670": float(_trial_passes(v5b, 1.0, bins, trials, rng).mean()),
            "0.700": float(_trial_passes(v5b, k_hi, bins, trials, rng).mean()),
        }
        assert rates[label]["0.640"] <= 0.07 and rates[label]["0.700"] <= 0.07, rates
        assert rates[label]["0.670"] >= 0.95, rates
    # joint pass of the row at the nominal coverage: sec_p profile, nuclear_local profile and
    # escaped_neutral scalar from independent draws (an estimate: the real estimators share
    # histories, and positive correlation raises the joint pass probability)
    n = 2000
    joint = (
        _trial_passes(v5b, 1.0, 12, n, rng)
        & _trial_passes(v5b, 1.0, 12, n, rng)
        & _trial_passes(v5b, 1.0, 1, n, rng)
    )
    rates["joint_row_0.670"] = float(joint.mean())
    assert rates["joint_row_0.670"] >= 0.90, rates
    print(f"v7 pipeline pass rates (kappa {k_lo:.4f}/{k_hi:.4f} for 0.640/0.700)", rates)
