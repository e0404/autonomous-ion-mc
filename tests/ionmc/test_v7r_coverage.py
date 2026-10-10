"""Calibration of the V7-R bootstrap-t gates of ``escaped_neutral`` (V3-005C C5; plan Amendment 14
(h), Amendment 17 (a)5-7). Marker ``calibration`` (host, exact SHA; excluded from CI). The step
``pytest-v7r-calibration`` runs this file and the unchanged Amendment 13 calibration
``tests/ionmc/test_v7_coverage.py``; its pass is the pytest exit status and ``summarize.py`` makes the
conjunction with ``v7r-combine`` for row V7-R.

Requirements (per surrogate, for the two bootstrap-t gates of ``escaped_neutral``)
-----------------------------------------------------------------------------------
* joint false acceptance of the two gates <= 0.02 at each margin (95 % Clopper-Pearson upper bound of
  the Monte Carlo rate); the four margins are the true single-interval coverage and the true pair
  coverage equal to the low and the high edge of [0.6527, 0.7127], reached by a scale ``kappa`` on the
  offsets of the intervals from the point estimate (as Amendment 13);
* each gate alone <= 0.04 (the exact alpha_tost of its construction): checked as "the 95 % CP LOWER
  bound of the Monte Carlo rate is <= 0.04" (the rate is not significantly above 0.04);
* power >= 0.9: the joint pass rate at kappa = 1 (the unscaled gates).

Surrogates (Amendment 14 (h))
-----------------------------
1. ``archive``: a zero-inflated Gamma block distribution (``x = B * Gamma(k)``, ``B ~ Bernoulli(1 -
   pi)``) fitted (grid, then local refinement, fixed fit seed) so that the studentised statistic
   ``(mean - mu) / sem`` of its simulated 20-block replicates reproduces the four empirical quantiles
   2.5 / 15.87 / 84.13 / 97.5 % of the archived 7200 lv5b replicates (7aae5bb; -2.228 / +2.005 at the
   outer pair, the inner pair computed deterministically from the archive). The fit record (family,
   parameters, achieved quantiles, seed) is printed and reported.
2. ``rehearsal``: 20 blocks drawn with replacement from the pooled rehearsal block sums (>= 2 shards
   of 4.5e6 histories, from the sidecars); truth = the pooled mean.

Fixtures (committed derived files; see ``FIXTURE_DIR``). Missing fixtures SKIP the test with a reason,
and FAIL when the environment variable ``IONMC_V7R_FIXTURES`` is ``required`` (the exact-SHA host run
of the evidence step sets it; CI and the sandbox do not):

* ``lv5b-7aae5bb-escaped-neutral-replicates.json`` (``v7r_fixture.py``): schema
  ``ionmc-v7r-lv5b-escaped-neutral-replicates-1``; keys ``mean``, ``sem`` (7200 floats each, shard
  order 0..7), ``sources`` (host run dir, file sha256, content_sha256, seed, git_sha, n of the 8 shard
  documents and the reference), ``git_sha_prefix`` ``7aae5bb1``.
* ``rehearsal/v7r-s{k}.json`` and ``rehearsal/v7r-s{k}-escaped_neutral.npz``, k = 0, 1, ...: the
  hash-sealed ``v7r-shard`` partials of a rehearsal run (seed base in 20505000-20509499) and their
  sidecars (schema of ``steps_v5c``: ``escaped_neutral_block_sums`` ``[450, 20]`` float64, MeV per
  500-history block); provenance fields in the partial: ``git_sha``, ``seed``, ``seed_base``,
  ``content_sha256``, ``sidecar.sha256`` (the orchestrator records the host run ids and the PARTIAL
  digests next to them, ``PROVENANCE.json`` in the same directory, and the hashed-source digest).

Reproducibility. Monte Carlo seeds are in the rehearsal family 20505000-20509499 (never 20509500 and
above, the E-shape family): repetition ``r`` uses data seed ``20507000 + r`` and bootstrap seed base
``20506000 + r``; the margin pre-simulations use 20508000-20508001 and the surrogate fit 20505000.
The number of repetitions is ``IONMC_V7R_CAL_REPS`` (default 200: with zero accepted repetitions the
95 % CP upper bound is 0.0149; with one it is 0.0236, above 0.02).

Recorded limitation (reported, not gated, under the chosen reading of Amendment 17 (a)5). For
``sec_p`` and ``nuclear_local`` the rule is the unchanged Amendment 13 rule and its calibration is
``tests/ionmc/test_v7_coverage.py``: the Gaussian 12-bin profile surrogates (the sampling distribution
of those profile means) have power 0.95-1.00 at the nominal; the skewed Gamma surrogates give a
single-gate power of 0.14, a recorded limitation of the distribution-free radius, not a gate. The
test ``test_v7r_readings_report_amendment13_skewed_power`` prints that number under the literal
reading ("each estimator and each surrogate", which cannot pass for these estimators irrespective of
the data) as the ``V7R-READINGS`` line and writes it to ``IONMC_V7R_REPORT_DIR/v7r-readings.json`` if
that variable is set. Under the literal reading the row cannot pass for sec_p/nuclear_local; the row
verdict uses the reading of Amendment 17 (a)5.
"""

# ruff: noqa: E501
from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
import time
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"
FIXTURE_DIR = REPO / "tests" / "ionmc" / "fixtures" / "v7r"
ARCHIVE_FIXTURE = FIXTURE_DIR / "lv5b-7aae5bb-escaped-neutral-replicates.json"
REHEARSAL_DIR = FIXTURE_DIR / "rehearsal"
SEED_FIT, SEED_KAPPA, SEED_BOOT0, SEED_DATA0 = 20505000, 20508000, 20506000, 20507000
QUANTILE_PROBS = (0.025, 0.1587, 0.8413, 0.975)
N_BOOT_Z = 1000  # resamples of the held-out box bootstrap inside the Monte Carlo (5000 by default)


def _load(name: str, path: Path | None = None) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, path or SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def v5b() -> ModuleType:
    return _load("steps_v5b")


@pytest.fixture(scope="module")
def v7r(v5b: ModuleType) -> ModuleType:
    return _load("v7r")


def _fixtures_required() -> bool:
    return os.environ.get("IONMC_V7R_FIXTURES", "") == "required"


def _need(path: Path, what: str) -> None:
    """Skip (or fail under ``IONMC_V7R_FIXTURES=required``) if a committed fixture is absent."""
    if not path.exists():
        msg = f"fixture missing: {what} ({path}); see the module docstring"
        if _fixtures_required():
            pytest.fail(msg)
        pytest.skip(msg)


# -- samplers ------------------------------------------------------------------------------------
class ZeroInflatedGamma:
    """Block values ``B * Gamma(k)``, ``B ~ Bernoulli(1 - pi)``; mean ``k (1 - pi)``."""

    def __init__(self, k: float, pi: float) -> None:
        self.k, self.pi = k, pi
        self.mu = k * (1.0 - pi)

    def __call__(self, rng: np.random.Generator, rows: int) -> np.ndarray:
        x = rng.standard_gamma(self.k, (rows, 20))
        return x * (rng.random((rows, 20)) >= self.pi)


class PooledBlocks:
    """20 blocks drawn with replacement from a pool; truth = pool mean."""

    def __init__(self, pool: np.ndarray) -> None:
        self.pool = np.asarray(pool, dtype=float).ravel()
        self.mu = float(self.pool.mean())

    def __call__(self, rng: np.random.Generator, rows: int) -> np.ndarray:
        return self.pool[rng.integers(0, self.pool.size, (rows, 20))]


def studentised_quantiles(x: np.ndarray, mu: float) -> np.ndarray:
    mean, se = x.mean(axis=1), x.std(axis=1, ddof=1) / math.sqrt(x.shape[1])
    ok = se > 0
    return np.quantile((mean[ok] - mu) / se[ok], QUANTILE_PROBS)


def fit_zero_inflated_gamma(targets: np.ndarray, seed: int = SEED_FIT, n: int = 40_000) -> dict:  # type: ignore[type-arg]
    """Least-squares fit of (k, pi) to the four studentised quantiles (common random numbers per
    candidate: the same seed; a coarse log grid and then two local refinements)."""

    def err(k: float, pi: float) -> tuple[float, np.ndarray]:
        s = ZeroInflatedGamma(k, pi)
        q = studentised_quantiles(s(np.random.default_rng(seed), n), s.mu)
        return float(np.sum((q - targets) ** 2)), q

    best = None
    for k in np.exp(np.linspace(math.log(0.1), math.log(8.0), 24)):
        for pi in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7):
            e, q = err(float(k), pi)
            if best is None or e < best[0]:
                best = (e, float(k), pi, q)
    assert best is not None
    for width in (0.15, 0.05):
        _, k0, p0, _ = best
        for k in k0 * np.exp(np.linspace(-width * 3, width * 3, 9)):
            for pi in np.clip(p0 + np.linspace(-width, width, 7), 0.0, 0.9):
                e, q = err(float(k), float(pi))
                if e < best[0]:
                    best = (e, float(k), float(pi), q)
    e, k, pi, q = best
    return {"family": "zero-inflated-gamma", "k": k, "pi": pi, "seed": seed, "n": n,
            "targets": targets.tolist(), "achieved": q.tolist(),
            "max_abs_error": float(np.max(np.abs(q - targets)))}  # fmt: skip


# -- calibration machinery ---------------------------------------------------------------------
def margin_kappas(
    v7r: ModuleType, sampler, kind: str, target: float, seed: int, n: int = 20_000
) -> float:  # type: ignore[no-untyped-def]
    """``kappa`` such that the true coverage of the bootstrap-t interval (``kind="single"``: of the
    mean ``sampler.mu``; ``"pair"``: of 0 for two independent replicates) equals ``target``, by
    bisection on the large-sample coverage (the coverage is monotone in kappa)."""
    rng = np.random.default_rng(seed)
    if kind == "single":
        x = sampler(rng, n)
        ing = v7r.single_ingredients(x, SEED_BOOT0 - 1, reps_per_shard=n, b=v7r.V7R_B)
        a = (ing["m_hat"] - sampler.mu) / ing["se"]
        lo, hi = ing["t_lo"], ing["t_hi"]
    else:
        x = sampler(rng, 2 * n)
        ing = v7r.pair_ingredients(x, SEED_BOOT0 - 1, reps_per_shard=n, b=v7r.V7R_B)
        a = ing["d_hat"] / ing["se_d"]
        lo, hi = ing["t_lo"], ing["t_hi"]
    ok = np.isfinite(a) & np.isfinite(lo) & np.isfinite(hi)

    def cov(k: float) -> float:
        return float(np.sum(ok & (k * hi >= a) & (k * lo <= a)) / a.size)

    lo_k, hi_k = 0.05, 20.0
    for _ in range(60):
        mid = 0.5 * (lo_k + hi_k)
        lo_k, hi_k = (mid, hi_k) if cov(mid) < target else (lo_k, mid)
    return 0.5 * (lo_k + hi_k)


def run_calibration(v7r: ModuleType, v5b: ModuleType, sampler, kappas: dict[str, float], reps: int, *,  # type: ignore[no-untyped-def]
                    replicates: int = 7200, reps_per_shard: int = 450, b: int = 1999,
                    n_boot: int = N_BOOT_Z) -> dict[str, np.ndarray]:  # fmt: skip
    """Monte Carlo of both gates at each named ``kappa``: ``{name: [reps, 2]}`` booleans (paired,
    single). One bootstrap per repetition serves every kappa (the resamples do not depend on it)."""
    out = {k: np.zeros((reps, 2), dtype=bool) for k in kappas}
    for r in range(reps):
        x = sampler(np.random.default_rng(SEED_DATA0 + r), replicates)
        seed_base = SEED_BOOT0 + r
        pi = v7r.pair_ingredients(x, seed_base, reps_per_shard, b)
        si = v7r.single_ingredients(x, seed_base, reps_per_shard, b)
        means = x.mean(axis=1)
        z_boot = None
        for name, kap in kappas.items():
            g = v7r.single_gate(si, means, kap, boot_seed=seed_base, n_boot=n_boot,
                                z_boot_override=z_boot)  # fmt: skip
            z_boot = g["z_boot"]
            out[name][r] = (v7r.pair_gate(pi, kap)["pass"], g["pass"])
    return out


def evaluate(v5b: ModuleType, rates: dict[str, np.ndarray], reps: int) -> dict:  # type: ignore[type-arg]
    """Verdict of the three requirements from ``run_calibration`` output (margins: every name but
    ``nominal``)."""
    rep = {}
    for name, p in rates.items():
        joint = p.all(axis=1)
        rep[name] = {
            "joint": int(joint.sum()), "paired": int(p[:, 0].sum()), "single": int(p[:, 1].sum()),
            "n": reps, "joint_cp95_upper": v5b.cp_upper(int(joint.sum()), reps, 0.05),
            "paired_cp95_lower": v5b.cp_lower(int(p[:, 0].sum()), reps, 0.05),
            "single_cp95_lower": v5b.cp_lower(int(p[:, 1].sum()), reps, 0.05),
            "joint_rate": float(joint.mean()),
        }  # fmt: skip
    failures = []
    for name, r in rep.items():
        if name == "nominal":
            if r["joint_rate"] < 0.9:
                failures.append(f"power {r['joint_rate']:.3f} < 0.9 at the nominal")
            continue
        if r["joint_cp95_upper"] > 0.02:
            failures.append(
                f"{name}: joint false acceptance CP95 upper {r['joint_cp95_upper']:.4f} > 0.02"
            )
        for gate in ("paired", "single"):
            if r[f"{gate}_cp95_lower"] > 0.04:
                failures.append(
                    f"{name}: {gate} gate false acceptance above 0.04 ({r[gate]}/{reps})"
                )
    return {"margins": rep, "failures": failures}


def _report(label: str, doc: dict) -> None:  # type: ignore[type-arg]
    print("V7R-CALIBRATION", json.dumps({label: doc}, sort_keys=True))
    d = os.environ.get("IONMC_V7R_REPORT_DIR")
    if d:
        Path(d).mkdir(parents=True, exist_ok=True)
        (Path(d) / f"v7r-calibration-{label}.json").write_text(
            json.dumps(doc, sort_keys=True, indent=1)
        )


def _calibrate(v7r: ModuleType, v5b: ModuleType, sampler, label: str, extra: dict) -> dict:  # type: ignore[no-untyped-def,type-arg]
    reps = int(os.environ.get("IONMC_V7R_CAL_REPS", "200"))
    lo, hi = v7r.V7R_REGION
    t0 = time.perf_counter()
    kappas = {
        "single_low": margin_kappas(v7r, sampler, "single", lo, SEED_KAPPA),
        "single_high": margin_kappas(v7r, sampler, "single", hi, SEED_KAPPA),
        "paired_low": margin_kappas(v7r, sampler, "pair", lo, SEED_KAPPA + 1),
        "paired_high": margin_kappas(v7r, sampler, "pair", hi, SEED_KAPPA + 1),
        "nominal": 1.0,
    }
    t1 = time.perf_counter()
    rates = run_calibration(v7r, v5b, sampler, kappas, reps)
    verdict = evaluate(v5b, rates, reps)
    doc = {**extra, "kappas": kappas, "reps": reps, **verdict, "wall_kappa_s": t1 - t0,
           "wall_mc_s": time.perf_counter() - t1, "seeds": {"fit": SEED_FIT, "kappa": SEED_KAPPA,
           "boot": SEED_BOOT0, "data": SEED_DATA0}, "rng": v7r.rng_record(), "B": v7r.V7R_B}  # fmt: skip
    _report(label, doc)
    return doc


# -- calibration tests (host) --------------------------------------------------------------------
@pytest.mark.calibration
def test_v7r_calibration_archive_fitted_skewed_surrogate(v7r: ModuleType, v5b: ModuleType) -> None:
    _need(ARCHIVE_FIXTURE, "archived lv5b replicates (mean, sem)")
    doc = json.loads(ARCHIVE_FIXTURE.read_text())
    assert doc["schema"] == "ionmc-v7r-lv5b-escaped-neutral-replicates-1"
    assert doc["git_sha_prefix"] == "7aae5bb1" and len(doc["mean"]) == len(doc["sem"]) == 7200
    m, s = np.array(doc["mean"]), np.array(doc["sem"])
    targets = np.quantile((m - m.mean()) / s, QUANTILE_PROBS)
    assert abs(targets[0] + 2.228) < 0.01 and abs(targets[3] - 2.005) < 0.01, targets  # plan values
    fit = fit_zero_inflated_gamma(targets)
    print("V7R-FIT", json.dumps(fit, sort_keys=True))
    assert fit["max_abs_error"] < 0.1, fit  # the surrogate reproduces the four quantiles
    sampler = ZeroInflatedGamma(fit["k"], fit["pi"])
    out = _calibrate(
        v7r,
        v5b,
        sampler,
        "archive",
        {"surrogate": "archive-fitted zero-inflated Gamma", "fit": fit},
    )
    assert not out["failures"], out["failures"]


@pytest.mark.calibration
def test_v7r_calibration_rehearsal_block_sums_surrogate(v7r: ModuleType, v5b: ModuleType) -> None:
    sidecars = (
        sorted(REHEARSAL_DIR.glob("v7r-s*-escaped_neutral.npz")) if REHEARSAL_DIR.exists() else []
    )
    if len(sidecars) < 2:
        _need(REHEARSAL_DIR / "v7r-s0-escaped_neutral.npz", "rehearsal sidecars (>= 2 shards)")
        msg = "fewer than two rehearsal sidecars"
        if _fixtures_required():
            pytest.fail(msg)
        pytest.skip(msg)
    steps = _load("steps_v5c")
    blocks = []
    for sc in sidecars:
        part = json.loads(
            sc.with_name(sc.name.replace("-escaped_neutral.npz", ".json")).read_text()
        )
        blocks.append(
            steps.load_v7r_sidecar(
                part, sc.with_name(sc.name.replace("-escaped_neutral.npz", ".json"))
            )
        )
        assert 20505000 <= part["seed_base"] <= 20509499, (
            "rehearsal seed base outside 20505000-20509499"
        )
    pool = np.concatenate(blocks).ravel()
    out = _calibrate(v7r, v5b, PooledBlocks(pool), "rehearsal",
                     {"surrogate": "pooled rehearsal block sums", "shards": len(sidecars), "blocks": int(pool.size)})  # fmt: skip
    assert not out["failures"], out["failures"]


@pytest.mark.calibration
def test_v7r_readings_report_amendment13_skewed_power(v5b: ModuleType) -> None:
    """Contrary reading (Amendment 17 (a)5): single-gate power of the Amendment 13 rule (sec_p /
    nuclear_local) on the skewed Gamma surrogate (shape 0.5, correlated bins), reported not gated."""
    cov = _load("test_v7_coverage_module", REPO / "tests" / "ionmc" / "test_v7_coverage.py")
    trials = int(os.environ.get("IONMC_V7R_READINGS_TRIALS", "40"))
    rng = np.random.default_rng(SEED_FIT + 1)
    passes, _ = cov._skew_passes(v5b, 0.5, True, False, (1.0,), trials, rng)
    power = {"row": float(passes[0, 0].mean()), "paired": float(passes[1, 0].mean()),
             "single": float(passes[2, 0].mean()), "trials": trials}  # fmt: skip
    doc = {"chosen_reading": "sec_p/nuclear_local power on the Amendment 13 Gaussian 12-bin profile "
                              "surrogates (test_v7_coverage.py, gated there)",
           "literal_reading": "each estimator and each surrogate: the skewed Gamma(0.5, correlated) "
                              "power at kappa = 1 below; cannot reach 0.9 irrespective of the data",
           "amendment13_skewed_gamma_power_kappa1": power, "gated": False}  # fmt: skip
    print("V7R-READINGS", json.dumps(doc, sort_keys=True))
    d = os.environ.get("IONMC_V7R_REPORT_DIR")
    if d:
        Path(d).mkdir(parents=True, exist_ok=True)
        (Path(d) / "v7r-readings.json").write_text(json.dumps(doc, sort_keys=True, indent=1))
    assert 0.0 <= power["single"] <= 1.0


# -- CI-scale tests of the calibration machinery (not marked: they run in CI) --------------------
def test_calibration_machinery_at_ci_scale(v7r: ModuleType, v5b: ModuleType) -> None:
    """The Monte Carlo, the margin solver and the evaluation run end to end on a tiny case."""
    sampler = ZeroInflatedGamma(6.0, 0.05)
    k = margin_kappas(v7r, sampler, "single", 0.6827, SEED_KAPPA, n=400)
    assert 0.5 < k < 2.0  # the nominal level is reached near kappa = 1
    k_lo = margin_kappas(v7r, sampler, "single", 0.55, SEED_KAPPA, n=400)
    k_hi = margin_kappas(v7r, sampler, "single", 0.80, SEED_KAPPA, n=400)
    assert k_lo < k < k_hi  # monotone in the target coverage
    rates = run_calibration(v7r, v5b, sampler, {"low": k_lo, "nominal": 1.0}, 3, replicates=120,
                            reps_per_shard=15, b=199, n_boot=50)  # fmt: skip
    assert all(v.shape == (3, 2) for v in rates.values())
    rep = evaluate(v5b, rates, 3)
    assert set(rep["margins"]) == {"low", "nominal"} and isinstance(rep["failures"], list)
    again = run_calibration(v7r, v5b, sampler, {"low": k_lo, "nominal": 1.0}, 3, replicates=120,
                            reps_per_shard=15, b=199, n_boot=50)  # fmt: skip
    assert all(np.array_equal(rates[n], again[n]) for n in rates)  # deterministic


def test_fit_recovers_a_known_surrogate_and_samplers_are_unbiased() -> None:
    truth = ZeroInflatedGamma(4.0, 0.1)
    targets = studentised_quantiles(truth(np.random.default_rng(5), 40_000), truth.mu)
    fit = fit_zero_inflated_gamma(targets, n=10_000)
    assert fit["max_abs_error"] < 0.08
    x = truth(np.random.default_rng(6), 5000)
    assert abs(x.mean() - truth.mu) < 0.05
    pool = PooledBlocks(x[:50])
    assert abs(pool(np.random.default_rng(7), 4000).mean() - pool.mu) < 0.2 * max(pool.mu, 1.0)


def test_fixture_absent_skips_and_required_fails(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    missing = tmp_path / "nope.json"
    monkeypatch.delenv("IONMC_V7R_FIXTURES", raising=False)
    with pytest.raises(pytest.skip.Exception):
        _need(missing, "x")
    monkeypatch.setenv("IONMC_V7R_FIXTURES", "required")
    with pytest.raises(pytest.fail.Exception):
        _need(missing, "x")
    present = tmp_path / "here.json"
    present.write_text("{}")
    _need(present, "x")  # no exception


def test_committed_archive_fixture_is_consistent_when_present() -> None:
    _need(ARCHIVE_FIXTURE, "archived lv5b replicates (mean, sem)")
    doc = json.loads(ARCHIVE_FIXTURE.read_text())
    assert len(doc["sources"]) == 9 and [s["seed"] for s in doc["sources"]] == list(
        range(20482004, 20482013)
    )
    assert all(s["git_sha"].startswith("7aae5bb1") for s in doc["sources"])
    assert len({s["content_sha256"] for s in doc["sources"]}) == 9
    assert all(math.isfinite(x) for x in doc["mean"] + doc["sem"]) and min(doc["sem"]) > 0
