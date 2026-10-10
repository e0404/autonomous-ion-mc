"""Calibration of the V7-R bootstrap-t gates of ``escaped_neutral`` (V3-005C C5; plan Amendment 14
(h), Amendment 17 (a)5-7). Marker ``calibration`` (host, exact SHA; excluded from CI). The step
``pytest-v7r-calibration`` runs this file and the unchanged Amendment 13 calibration
``tests/ionmc/test_v7_coverage.py``; its pass is the pytest exit status and ``summarize.py`` makes
the
conjunction with ``v7r-combine`` for row V7-R.

Requirements (per surrogate, for the two bootstrap-t gates of ``escaped_neutral``)
-----------------------------------------------------------------------------------
* joint false acceptance of the two gates <= 0.02 at each margin (95 % Clopper-Pearson upper bound
  of
  the Monte Carlo rate); the four margins are the true single-interval coverage and the true pair
  coverage equal to the low and the high edge of [0.6527, 0.7127], reached by a scale ``kappa`` on
  the
  offsets of the intervals from the point estimate (as Amendment 13);
* each gate alone <= 0.04 (the exact alpha_tost of its construction): the 95 % CP UPPER bound of the
  Monte Carlo acceptance rate of that gate is <= 0.04 at that gate's OWN margins (paired gate at the
  paired_low/high kappas, single gate at the single_low/high kappas); the other gate's values at a
  margin are reported, not gated;
* power >= 0.9: the 95 % CP LOWER bound of the joint pass rate at kappa = 1 (the unscaled gates) is
  >= 0.9.

Every "<=" requirement uses the one-sided 95 % Clopper-Pearson UPPER bound (``cp_upper(k, n,
0.05)``)
and the power requirement the LOWER bound. Maximum number of accepted repetitions ``k`` of ``n`` for
which the upper bound is still <= the limit, and the minimum number of accepted repetitions at the
nominal for a lower bound >= 0.9 (exact values, ``steps_v5b.cp_upper/cp_lower``):

====  =================  =================  ==================
n     joint <= 0.02      each gate <= 0.04  power >= 0.9
      (max accepted k)   (max accepted k)   (min accepted k)
====  =================  =================  ==================
200   0 (0.0149)         3 (0.0383)         188 (0.9046)
400   3 (0.0193)         9 (0.0389)         371 (0.9025)
800   9 (0.0195)         22 (0.0390)        735 (0.9011)
1600  22 (0.0196)        50 (0.0394)        1461 (0.9007)
====  =================  =================  ==================

One more acceptance than the maximum gives the bound 0.0235 / 0.0227 / 0.0211 / 0.0203 (joint) and
0.0452 / 0.0420 / 0.0405 / 0.0401 (each gate); one fewer than the minimum gives the lower bound
0.8986
/ 0.8996 / 0.8997 / 0.9000 (power).

(The tests ``test_requirement_bounds_*`` pin these numbers.) The repetition count is
``IONMC_V7R_CAL_REPS`` (default 200). A true false-acceptance rate at the limit (0.02) needs several
hundred repetitions to be certified; 200 repetitions can only certify a rate near 0 for the joint
requirement, which is why more repetitions may be chosen before the evidence run (Amendment 14 (h)
4).

Surrogates (Amendment 14 (h))
-----------------------------
1. ``archive``: a zero-inflated Gamma block distribution (``x = B * Gamma(k)``, ``B ~ Bernoulli(1 -
   pi)``) fitted (grid, then local refinement, fixed fit seed) so that the studentised statistic
   ``(mean - mu) / sem`` of its simulated 20-block replicates reproduces the four empirical
   quantiles
   2.5 / 15.87 / 84.13 / 97.5 % of the archived 7200 lv5b replicates (7aae5bb; -2.228 / +2.005 at
   the
   outer pair, the inner pair computed deterministically from the archive). The fit record (family,
   parameters, achieved quantiles, seed) is printed and reported.
2. ``rehearsal``: 20 blocks drawn with replacement from the pooled rehearsal block sums (>= 2 shards
   of 4.5e6 histories, from the sidecars); truth = the pooled mean.

Fixtures (committed derived files; see ``FIXTURE_DIR``). A missing or rejected fixture SKIPs the
test
with a reason, and FAILs when the environment variable ``IONMC_V7R_FIXTURES`` is ``required``. In
required mode NO calibration test skips for any reason (every skip path is ``_skip_or_fail``;
``test_no_skip_path_outside_skip_or_fail`` scans both calibration files); ``run_suite.py`` sets the
variable for every ``pytest-v7r-calibration*`` step and ``summarize.py`` rejects a step whose
archived pytest summary shows a skipped test. ``IONMC_V7R_REHEARSAL_DIR`` points the rehearsal
loader at another directory (tests and smoke runs only; it is REJECTED, a failure, in required
mode, and ``run_suite`` removes it from the environment of every calibration step and records the
digest of the committed fixture directory in the step header):

* ``lv5b-7aae5bb-escaped-neutral-replicates.json`` (``v7r_fixture.py``): schema
  ``ionmc-v7r-lv5b-escaped-neutral-replicates-1``; keys ``mean``, ``sem`` (7200 floats each, shard
  order 0..7), ``sources`` (host run dir, file sha256, content_sha256, seed, git_sha, n of the 8
  shard
  documents and the reference), ``git_sha_prefix`` ``7aae5bb1``.
* ``rehearsal/v7r-s{k}.json`` and ``rehearsal/v7r-s{k}-escaped_neutral.npz``, k = 0, 1, ...: the
  hash-sealed ``v7r-shard`` partials of a rehearsal run (seed base in 20505000-20509499) and their
  sidecars (schema of ``steps_v5c``: ``escaped_neutral_block_sums`` ``[450, 20]`` float64, MeV per
  500-history block), and ``PROVENANCE.json`` (required) with ``host_run_ids`` (non-empty list),
  ``partial_digests`` ({partial file name: content_sha256}, exactly the partials present),
  ``hashed_source_digest``, ``seed_base`` and ``code_sha``. ``load_rehearsal`` enforces (every
  violation
  is reported): each partial's ``content_sha256`` reproduces (``steps_v5.content_digest``),
  ``valid``
  true and ``reduced`` false, ``n`` 4.5e6, ``replicates`` 450, ``batch_histories`` 500,
  ``blocks_per_replicate`` 20, at least two shards with unique shard indices and unique seeds, the
  seed of shard k equal to ``seed_base + 1000 * 16 + k`` (the shard seeds themselves lie above the
  2050xxxx family; the seed BASE is the 20505000-20509499 quantity), ``seed_base`` equal to the
  provenance, ``git_sha`` equal to ``code_sha``, the sidecar file and array sha256 as bound in the
  partial (``steps_v5c.load_v7r_sidecar``), the digests of ``partial_digests``, and
  ``hashed_source_digest`` equal to the CURRENT ``run_suite.hashed_source_digest("lv5c")`` (the
  sha256
  over the sorted ``path sha256`` lines of ``source_file_list("lv5c")`` plus its prefixes, without
  the
  rehearsal fixtures themselves, so that committing them does not change the digest; Amendment 17
  (a)7: the rehearsal is valid only for identical hashed source). ``python
  validation/scripts/transport/run_suite.py --suite lv5c --out x --print-hashed-source-digest``
  prints
  it for the orchestrator.

Reproducibility. Monte Carlo seeds are in the rehearsal family 20505000-20509499 (never 20509500 and
above, the E-shape family): repetition ``r`` uses data seed ``20507000 + r`` and bootstrap seed base
``20506000 + r``; the margin pre-simulations use 20508000-20508001 and the surrogate fit 20505000.
The number of bootstrap resamples of the held-out box inside the gate is the frozen qualification
value ``steps_v5b.V7_BOOT_N`` (5000), the same constant ``v7r.single_gate`` uses for the evidence
(no
local override; ``test_calibration_uses_the_qualification_boot_count``).

Recorded limitation (reported, not gated, under the chosen reading of Amendment 17 (a)5). For
``sec_p`` and ``nuclear_local`` the rule is the unchanged Amendment 13 rule and its calibration is
``tests/ionmc/test_v7_coverage.py``: the Gaussian 12-bin profile surrogates (the sampling
distribution
of those profile means) have power 0.95-1.00 at the nominal; the skewed Gamma surrogates give a
single-gate power of 0.14, a recorded limitation of the distribution-free radius, not a gate. The
test ``test_v7r_readings_report_amendment13_skewed_power`` prints that number under the literal
reading ("each estimator and each surrogate", which cannot pass for these estimators irrespective of
the data) as the ``V7R-READINGS`` line and writes it to ``IONMC_V7R_REPORT_DIR/v7r-readings.json``
if
that variable is set. Under the literal reading the row cannot pass for sec_p/nuclear_local; the row
verdict uses the reading of Amendment 17 (a)5.
"""

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
DEFAULT_REPS = 200
EXPECTED_SHARD = {
    "n": 4_500_000,
    "replicates": 450,
    "batch_histories": 500,
    "blocks_per_replicate": 20,
}
SEED_BASE_RANGE = (20505000, 20509499)


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


def _skip_or_fail(msg: str) -> None:
    """The only skip path of the calibration tests: a failure in required mode."""
    if _fixtures_required():
        pytest.fail(msg)
    pytest.skip(msg)


def _need(path: Path, what: str) -> None:
    """Skip (or fail under ``IONMC_V7R_FIXTURES=required``) if a committed fixture is absent."""
    if not path.exists():
        _skip_or_fail(f"fixture missing: {what} ({path}); see the module docstring")


def rehearsal_dir() -> Path:
    """The rehearsal fixture directory (``IONMC_V7R_REHEARSAL_DIR`` overrides the committed one)."""
    env = os.environ.get("IONMC_V7R_REHEARSAL_DIR")
    if env and _fixtures_required():
        pytest.fail(
            "IONMC_V7R_REHEARSAL_DIR is set in required mode: the qualification step must use the "
            "committed rehearsal fixture directory (the override is for tests and smoke runs only)"
        )
    return Path(env) if env else REHEARSAL_DIR


def cal_reps() -> int:
    return int(os.environ.get("IONMC_V7R_CAL_REPS", str(DEFAULT_REPS)))


class RehearsalError(Exception):
    """The rehearsal fixture violates at least one requirement; ``str`` lists every violation."""


def load_rehearsal(
    directory: Path, steps: ModuleType, run_suite: ModuleType, digest: str | None = None
) -> tuple[np.ndarray, dict]:  # type: ignore[type-arg]
    """Verified pooled block sums ``[replicates, 20]`` of the rehearsal in ``directory`` and a
    record of what was verified, or ``RehearsalError`` listing every violation (module docstring).
    Seed reading (Amendment 14 (j)): the rehearsal BASE lies in 20505000-20509499 and shard seeds
    are ``base + 16000 + k`` (20521000+ for base 20505000). These coincide with the declared future
    hr5c base 20521004 only as a number, never as a used seed (hr5c seeds are
    ``base + 1000 r + k``, r >= 1).
    ``digest`` overrides the current hashed-source digest (tests only)."""
    bad: list[str] = []
    v7r = steps.v7r
    prov_path = directory / "PROVENANCE.json"
    prov: dict = {}  # type: ignore[type-arg]
    if not prov_path.is_file():
        bad.append(f"PROVENANCE.json missing in {directory}")
    else:
        try:
            prov = json.loads(prov_path.read_text())
        except ValueError as exc:
            bad.append(f"PROVENANCE.json unreadable: {exc}")
    if prov:
        for key in (
            "host_run_ids",
            "partial_digests",
            "hashed_source_digest",
            "seed_base",
            "code_sha",
        ):
            if key not in prov:
                bad.append(f"PROVENANCE.json lacks {key}")
        ids = prov.get("host_run_ids")
        if not (isinstance(ids, list) and ids and all(isinstance(i, str) and i for i in ids)):
            bad.append("PROVENANCE.json host_run_ids must be a non-empty list of run ids")
        if not (
            isinstance(prov.get("seed_base"), int)
            and SEED_BASE_RANGE[0] <= prov["seed_base"] <= SEED_BASE_RANGE[1]
        ):
            bad.append(
                f"PROVENANCE.json seed_base {prov.get('seed_base')!r} outside {SEED_BASE_RANGE}"
            )
        have = digest if digest is not None else run_suite.hashed_source_digest("lv5c")
        if prov.get("hashed_source_digest") != have:
            bad.append(
                "hashed_source_digest differs from the current hashed source set of lv5c "
                "(the rehearsal is not at the evidence code, Amendment 17 (a)7)"
            )
    parts = sorted(directory.glob("v7r-s*.json"))
    if len(parts) < 2:
        bad.append(f"{len(parts)} rehearsal partial(s) in {directory}, need >= 2 shards")
    recorded = prov.get("partial_digests") if isinstance(prov.get("partial_digests"), dict) else {}
    if prov and {p.name for p in parts} != set(recorded):
        bad.append(f"partial_digests names {sorted(recorded)} differ from the partials present")
    blocks, shards, seeds = [], [], []
    for path in parts:
        nm = path.name
        try:
            doc = json.loads(path.read_text())
        except ValueError as exc:
            bad.append(f"{nm}: unreadable: {exc}")
            continue
        dig = steps.v5.content_digest(doc)
        if doc.get("content_sha256") != dig:
            bad.append(f"{nm}: content_sha256 does not reproduce")
        if nm in recorded and recorded[nm] != dig:
            bad.append(f"{nm}: differs from the digest recorded in PROVENANCE.json")
        if doc.get("valid") is not True:
            bad.append(f"{nm}: valid is not true")
        if doc.get("reduced") is not False:
            bad.append(f"{nm}: reduced is not false")
        for key, want in EXPECTED_SHARD.items():
            if doc.get(key) != want:
                bad.append(f"{nm}: {key} is {doc.get(key)!r}, expected {want}")
        shards.append(doc.get("shard"))
        seeds.append(doc.get("seed"))
        sb = doc.get("seed_base")
        if not (isinstance(sb, int) and SEED_BASE_RANGE[0] <= sb <= SEED_BASE_RANGE[1]):
            bad.append(f"{nm}: seed_base {sb!r} outside {SEED_BASE_RANGE}")
        if prov and sb != prov.get("seed_base"):
            bad.append(f"{nm}: seed_base differs from PROVENANCE.json")
        if (
            isinstance(sb, int)
            and isinstance(doc.get("shard"), int)
            and (doc.get("seed") != sb + 1000 * v7r.V7R_R_INDEX + doc["shard"])
        ):
            bad.append(f"{nm}: seed is not seed_base + 1000 * 16 + shard")
        if prov and doc.get("git_sha") != prov.get("code_sha"):
            bad.append(f"{nm}: git_sha differs from code_sha of PROVENANCE.json")
        try:
            blocks.append(steps.load_v7r_sidecar(doc, path))
        except SystemExit as exc:
            bad.append(f"{nm}: sidecar rejected: {exc}")
    if len(set(shards)) != len(shards):
        bad.append(f"shard indices not unique: {shards}")
    if len(set(seeds)) != len(seeds):
        bad.append(f"seeds not unique: {seeds}")
    if bad:
        raise RehearsalError("; ".join(bad))
    info = {"shards": len(parts), **{k: prov[k] for k in
            ("seed_base", "code_sha", "host_run_ids", "hashed_source_digest")}}  # fmt: skip
    return np.concatenate(blocks), info


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


def run_calibration(  # type: ignore[no-untyped-def]
    v7r: ModuleType, v5b: ModuleType, sampler, kappas: dict[str, float], reps: int, *,
    replicates: int = 7200, reps_per_shard: int = 450, b: int = 1999, n_boot: int | None = None,
) -> dict[str, np.ndarray]:  # fmt: skip
    """Monte Carlo of both gates at each named ``kappa``: ``{name: [reps, 2]}`` booleans (paired,
    single). One bootstrap per repetition serves every kappa (the resamples do not depend on it).
    ``n_boot=None`` is the qualification value of the single gate (``steps_v5b.V7_BOOT_N``, read by
    ``v7r.single_gate``); only the CI-scale test passes a number."""
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


def _own_gates(margin: str) -> tuple[str, ...]:
    """Gates whose "alone <= 0.04" requirement applies at ``margin``: the paired gate at the
    ``paired_*`` margins and the single gate at the ``single_*`` margins (its alpha_tost
    construction). The other gate's values at that margin are reported (``*_cp95_upper``) but not
    gated; a margin name with neither prefix (synthetic tests) gates both."""
    for gate in ("paired", "single"):
        if margin.startswith(gate):
            return (gate,)
    return ("paired", "single")


def evaluate(v5b: ModuleType, rates: dict[str, np.ndarray], reps: int) -> dict:  # type: ignore[type-arg]
    """Verdict of the three requirements from ``run_calibration`` output (margins: every name but
    ``nominal``). The gate-alone requirement is checked at each gate's own margins only
    (``_own_gates``; cross-margin values are report-only). Every "<=" requirement uses the
    one-sided 95 % CP UPPER bound, the power the LOWER."""
    rep = {}
    for name, p in rates.items():
        joint = p.all(axis=1)
        kj, kp, ks = int(joint.sum()), int(p[:, 0].sum()), int(p[:, 1].sum())
        rep[name] = {
            "joint": kj, "paired": kp, "single": ks, "n": reps,
            "joint_cp95_upper": v5b.cp_upper(kj, reps, 0.05),
            "paired_cp95_upper": v5b.cp_upper(kp, reps, 0.05),
            "single_cp95_upper": v5b.cp_upper(ks, reps, 0.05),
            "joint_cp95_lower": v5b.cp_lower(kj, reps, 0.05),
            "joint_rate": float(joint.mean()),
        }  # fmt: skip
    failures = []
    for name, r in rep.items():
        if name == "nominal":
            if r["joint_cp95_lower"] < 0.9:
                failures.append(
                    f"power: joint CP95 lower {r['joint_cp95_lower']:.4f} < 0.9 at the nominal "
                    f"({r['joint']}/{reps})"
                )
            continue
        if r["joint_cp95_upper"] > 0.02:
            failures.append(
                f"{name}: joint false acceptance CP95 upper {r['joint_cp95_upper']:.4f} > 0.02 "
                f"({r['joint']}/{reps})"
            )
        for gate in _own_gates(name):
            if r[f"{gate}_cp95_upper"] > 0.04:
                failures.append(
                    f"{name}: {gate} gate false acceptance CP95 upper "
                    f"{r[f'{gate}_cp95_upper']:.4f} > 0.04 ({r[gate]}/{reps})"
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
    reps = cal_reps()
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
    seeds = {"fit": SEED_FIT, "kappa": SEED_KAPPA, "boot": SEED_BOOT0, "data": SEED_DATA0}
    doc = {**extra, "kappas": kappas, "reps": reps, **verdict, "wall_kappa_s": t1 - t0,
           "wall_mc_s": time.perf_counter() - t1, "seeds": seeds, "rng": v7r.rng_record(),
           "B": v7r.V7R_B, "boot_n": v5b.V7_BOOT_N}  # fmt: skip
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
    try:
        pool_blocks, info = load_rehearsal(rehearsal_dir(), _load("steps_v5c"), _load("run_suite"))
    except RehearsalError as exc:
        _skip_or_fail(f"rehearsal fixture rejected: {exc}")
        return
    pool = pool_blocks.ravel()
    extra = {"surrogate": "pooled rehearsal block sums", "blocks": int(pool.size), **info}
    out = _calibrate(v7r, v5b, PooledBlocks(pool), "rehearsal", extra)
    assert not out["failures"], out["failures"]


@pytest.mark.calibration
def test_v7r_readings_report_amendment13_skewed_power(v5b: ModuleType) -> None:
    """Contrary reading (Amendment 17 (a)5): single-gate power of the Amendment 13 rule (sec_p /
    nuclear_local) on the skewed Gamma surrogate (shape 0.5, correlated bins), reported not
    gated."""
    cov = _load("test_v7_coverage_module", REPO / "tests" / "ionmc" / "test_v7_coverage.py")
    trials = int(os.environ.get("IONMC_V7R_READINGS_TRIALS", "40"))
    rng = np.random.default_rng(SEED_FIT + 1)
    passes, _ = cov._skew_passes(v5b, 0.5, True, False, (1.0,), trials, rng)
    power = {"row": float(passes[0, 0].mean()), "paired": float(passes[1, 0].mean()),
             "single": float(passes[2, 0].mean()), "trials": trials}  # fmt: skip
    chosen = (
        "sec_p/nuclear_local power on the Amendment 13 Gaussian 12-bin profile surrogates "
        "(test_v7_coverage.py, gated there)"
    )
    literal = (
        "each estimator and each surrogate: the skewed Gamma(0.5, correlated) power at "
        "kappa = 1 below; cannot reach 0.9 irrespective of the data"
    )
    doc = {"chosen_reading": chosen, "literal_reading": literal,
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


# -- requirement bounds (Clopper-Pearson upper for "<=", lower for the power) ---------------------
def _rates(joint: int, paired: int, single: int, n: int) -> np.ndarray:
    """``[n, 2]`` booleans: the paired gate accepts ``paired`` rows, the single gate ``single``
    rows, and ``joint`` rows are accepted by both."""
    assert joint <= min(paired, single) and paired - joint + single <= n
    p = np.zeros((n, 2), dtype=bool)
    p[:paired, 0] = True
    p[paired - joint : paired - joint + single, 1] = True
    assert int(p.all(axis=1).sum()) == joint and int(p[:, 1].sum()) == single
    return p


def _verdict(v5b: ModuleType, margin: tuple[int, int, int], nominal: int, n: int) -> list[str]:
    rates = {"m": _rates(*margin, n), "nominal": _rates(nominal, nominal, nominal, n)}
    return evaluate(v5b, rates, n)["failures"]


@pytest.mark.parametrize(
    ("n", "joint_max", "gate_max", "power_min"),
    [(200, 0, 3, 188), (400, 3, 9, 371), (800, 9, 22, 735), (1600, 22, 50, 1461)],
)
def test_requirement_bounds_boundaries_match_the_docstring_table(
    v5b: ModuleType, n: int, joint_max: int, gate_max: int, power_min: int
) -> None:
    ok = (joint_max, gate_max, gate_max)
    assert _verdict(v5b, ok, power_min, n) == []
    assert _verdict(v5b, (joint_max + 1, gate_max, gate_max), power_min, n)  # joint: one more
    assert _verdict(v5b, (0, gate_max + 1, gate_max), power_min, n)  # paired gate alone
    assert _verdict(v5b, (0, gate_max, gate_max + 1), power_min, n)  # single gate alone
    assert _verdict(v5b, ok, power_min - 1, n)  # power one short
    assert v5b.cp_upper(joint_max, n, 0.05) <= 0.02 < v5b.cp_upper(joint_max + 1, n, 0.05)
    assert v5b.cp_upper(gate_max, n, 0.05) <= 0.04 < v5b.cp_upper(gate_max + 1, n, 0.05)
    assert v5b.cp_lower(power_min, n, 0.05) >= 0.9 > v5b.cp_lower(power_min - 1, n, 0.05)


def test_requirement_bounds_use_the_upper_bound_not_the_lower(v5b: ModuleType) -> None:
    """10/200 gate acceptances: the 95 % lower bound (0.027) is below 0.04, so the former
    lower-bound check would pass, the upper bound (0.083) fails; 1/200 (joint) likewise."""
    assert v5b.cp_lower(10, 200, 0.05) < 0.04 < v5b.cp_upper(10, 200, 0.05)
    assert (
        _verdict(v5b, (0, 10, 0), 200, 200)
        and "paired gate" in _verdict(v5b, (0, 10, 0), 200, 200)[0]
    )
    assert "single gate" in _verdict(v5b, (0, 0, 10), 200, 200)[0]
    assert v5b.cp_lower(1, 200, 0.05) < 0.02 < v5b.cp_upper(1, 200, 0.05)
    assert "joint" in _verdict(v5b, (1, 1, 1), 200, 200)[0]
    assert _verdict(v5b, (0, 0, 0), 200, 200) == []
    # power: 180/200 (rate 0.90, lower bound 0.85) fails, 190/200 passes
    assert "power" in _verdict(v5b, (0, 0, 0), 180, 200)[0]
    assert _verdict(v5b, (0, 0, 0), 190, 200) == []


def test_gate_alone_requirement_applies_at_each_gates_own_margins_only(v5b: ModuleType) -> None:
    n = 200
    nominal = _rates(190, 190, 190, n)
    # the paired gate accepts 10/200 at a SINGLE margin: reported, not gated
    rates = {"single_low": _rates(0, 10, 0, n), "nominal": nominal}
    out = evaluate(v5b, rates, n)
    assert out["failures"] == [] and out["margins"]["single_low"]["paired_cp95_upper"] > 0.04
    # the same acceptances at the PAIRED margin are gated
    out = evaluate(v5b, {"paired_high": _rates(0, 10, 0, n), "nominal": nominal}, n)
    assert len(out["failures"]) == 1 and "paired gate" in out["failures"][0]
    # the single gate at a single margin is gated; at a paired margin it is not
    out = evaluate(v5b, {"single_high": _rates(0, 0, 10, n), "nominal": nominal}, n)
    assert len(out["failures"]) == 1 and "single gate" in out["failures"][0]
    assert (
        evaluate(v5b, {"paired_low": _rates(0, 0, 10, n), "nominal": nominal}, n)["failures"] == []
    )
    # the joint requirement is gated at every margin
    out = evaluate(v5b, {"paired_low": _rates(1, 1, 1, n), "nominal": nominal}, n)
    assert "joint" in out["failures"][0]


def test_rehearsal_override_is_rejected_in_required_mode(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setenv("IONMC_V7R_REHEARSAL_DIR", str(tmp_path))
    monkeypatch.setenv("IONMC_V7R_FIXTURES", "required")
    with pytest.raises(pytest.fail.Exception, match="IONMC_V7R_REHEARSAL_DIR is set in required"):
        rehearsal_dir()
    monkeypatch.setenv("IONMC_V7R_FIXTURES", "")
    assert rehearsal_dir() == tmp_path  # tests and smoke runs only
    monkeypatch.delenv("IONMC_V7R_REHEARSAL_DIR")
    monkeypatch.setenv("IONMC_V7R_FIXTURES", "required")
    assert rehearsal_dir() == REHEARSAL_DIR  # the committed directory


def test_cal_reps_default_and_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("IONMC_V7R_CAL_REPS", raising=False)
    assert cal_reps() == DEFAULT_REPS == 200
    monkeypatch.setenv("IONMC_V7R_CAL_REPS", "400")
    assert cal_reps() == 400


# -- bootstrap count -----------------------------------------------------------------------------
def test_calibration_uses_the_qualification_boot_count(
    monkeypatch: pytest.MonkeyPatch, v7r: ModuleType
) -> None:
    """The calibration passes no resample count of its own: it and the evidence ``single_gate`` both
    read ``steps_v5b.V7_BOOT_N`` (5000), also when that constant is patched."""
    assert v7r.v5b.V7_BOOT_N == 5000
    assert not hasattr(sys.modules[__name__], "N_BOOT_Z")
    seen: list[int] = []
    real = v7r.v5b.v7_boot_z

    def spy(means, n_h, mu, sem, seed, n_boot):  # type: ignore[no-untyped-def]
        seen.append(int(n_boot))
        return real(means, n_h, mu, sem, seed, 7)  # small and fast; the requested count is recorded

    monkeypatch.setattr(v7r.v5b, "v7_boot_z", spy)
    sampler = ZeroInflatedGamma(6.0, 0.05)

    def counts() -> tuple[int, int, int]:
        seen.clear()
        x = sampler(np.random.default_rng(3), 120)
        si = v7r.single_ingredients(x, 5, reps_per_shard=15, b=99)
        gate = v7r.single_gate(si, x.mean(axis=1), 1.0, boot_seed=5)
        direct = seen[-1]
        run_calibration(v7r, v7r.v5b, sampler, {"nominal": 1.0}, 1, replicates=120,
                        reps_per_shard=15, b=99)  # fmt: skip
        return direct, seen[-1], gate["boot_n"]

    assert counts() == (5000, 5000, 5000)
    monkeypatch.setattr(v7r.v5b, "V7_BOOT_N", 37)
    assert counts() == (37, 37, 37)


# -- rehearsal fixture loader --------------------------------------------------------------------
def _build_rehearsal(  # type: ignore[no-untyped-def]
    d: Path,
    steps: ModuleType,
    digest: str,
    *,
    n_shards: int = 2,
    doc_hook=None,
    prov_hook=None,
    post=None,
):
    """Write a synthetic, valid rehearsal into ``d`` (450 x 20 block sums per shard) and return the
    provenance; ``doc_hook(k, doc)`` mutates a partial before it is sealed, ``prov_hook(prov)`` the
    provenance before it is written, ``post(d)`` tampers with the files afterwards."""
    d.mkdir(parents=True, exist_ok=True)
    base_seed, rng = 20505004, np.random.default_rng(11)
    prov = {"host_run_ids": ["RUN-SYNTH-1"], "partial_digests": {}, "hashed_source_digest": digest,
            "seed_base": base_seed, "code_sha": "c0de" * 10}  # fmt: skip
    for k in range(n_shards):
        sums = rng.gamma(2.0, 100.0, (450, 20))
        side = steps.write_v7r_sidecar(d / f"v7r-s{k}-escaped_neutral.npz", sums)
        per = sums / 500.0
        sem = (per.std(axis=1, ddof=1) / 20**0.5).tolist()
        doc = {"shard": k, "seed": base_seed + 16000 + k, "seed_base": base_seed,
               "git_sha": prov["code_sha"], "n": 4_500_000, "replicates": 450,
               "batch_histories": 500, "blocks_per_replicate": 20, "valid": True, "reduced": False,
               "estimators": {"escaped_neutral": {"mean": per.mean(axis=1).tolist(), "sem": sem}},
               "sidecar": side}  # fmt: skip
        if doc_hook:
            doc_hook(k, doc)
        doc["content_sha256"] = steps.v5.content_digest(doc)
        (d / f"v7r-s{k}.json").write_text(json.dumps(doc))
        prov["partial_digests"][f"v7r-s{k}.json"] = doc["content_sha256"]
    if prov_hook:
        prov_hook(prov)
    (d / "PROVENANCE.json").write_text(json.dumps(prov))
    if post:
        post(d)
    return prov


def _tamper_json(name: str, fn):  # type: ignore[no-untyped-def]
    def post(d: Path) -> None:
        p = d / name
        doc = json.loads(p.read_text())
        fn(doc)
        p.write_text(json.dumps(doc))

    return post


def _set(key: str, value: object):  # type: ignore[no-untyped-def]
    return lambda k, doc: doc.__setitem__(key, value) if k == 0 else None


def _drop_prov(key: str):  # type: ignore[no-untyped-def]
    return lambda prov: prov.pop(key)


def _sidecar_tamper(d: Path) -> None:
    f = d / "v7r-s1-escaped_neutral.npz"
    f.write_bytes(f.read_bytes() + b"x")


def _rm_prov(d: Path) -> None:
    (d / "PROVENANCE.json").unlink()


def _doc0(**kw: object):  # type: ignore[no-untyped-def]
    """Hook: update the fields of the partial of shard 0 before sealing."""
    return lambda k, doc: doc.update(kw) if k == 0 else None


def _prov(**kw: object):  # type: ignore[no-untyped-def]
    return lambda prov: prov.update(kw)


def _bad_digest(prov: dict) -> None:  # type: ignore[type-arg]
    prov["partial_digests"]["v7r-s0.json"] = "0" * 64


def _lose_partial(prov: dict) -> None:  # type: ignore[type-arg]
    prov["partial_digests"].pop("v7r-s1.json")


def _bad_sidecar_hash(k: int, doc: dict) -> None:  # type: ignore[type-arg]
    if k == 1:
        doc["sidecar"]["array_sha256"] = "2" * 64


def _wall(doc: dict) -> None:  # type: ignore[type-arg]
    doc["wall_s"] = 1.0


SEED0 = 20505004 + 16000
VIOLATIONS = {
    "partial tampered after sealing": (
        {"post": _tamper_json("v7r-s0.json", _wall)}, "content_sha256 does not reproduce"),
    "valid false": ({"doc_hook": _doc0(valid=False)}, "valid is not true"),
    "reduced true": ({"doc_hook": _doc0(reduced=True)}, "reduced is not false"),
    "n": ({"doc_hook": _doc0(n=4_000_000)}, "n is 4000000"),
    "replicates": ({"doc_hook": _doc0(replicates=449)}, "replicates is 449"),
    "batch_histories": ({"doc_hook": _doc0(batch_histories=250)}, "batch_histories is 250"),
    "blocks_per_replicate": (
        {"doc_hook": _doc0(blocks_per_replicate=10)}, "blocks_per_replicate is 10"),
    "one shard only": ({"n_shards": 1}, "need >= 2 shards"),
    "duplicate shard index": (
        {"doc_hook": lambda k, d: d.update(shard=0, seed=SEED0)}, "shard indices not unique"),
    "duplicate seed": ({"doc_hook": lambda k, d: d.update(seed=SEED0)}, "seeds not unique"),
    "seed base out of range": (
        {"doc_hook": lambda k, d: d.update(seed_base=20509500, seed=20509500 + 16000 + k),
         "prov_hook": _prov(seed_base=20509500)}, "outside (20505000, 20509499)"),
    "seed base differs from provenance": (
        {"prov_hook": _prov(seed_base=20505005)}, "seed_base differs from PROVENANCE"),
    "provenance missing": ({"post": _rm_prov}, "PROVENANCE.json missing"),
    "provenance lacks field": ({"prov_hook": _drop_prov("host_run_ids")}, "lacks host_run_ids"),
    "provenance empty run ids": ({"prov_hook": _prov(host_run_ids=[])}, "non-empty list"),
    "provenance lacks code_sha": ({"prov_hook": _drop_prov("code_sha")}, "lacks code_sha"),
    "code_sha differs": ({"prov_hook": _prov(code_sha="deadbeef")}, "differs from code_sha"),
    "provenance partial digest wrong": (
        {"prov_hook": _bad_digest}, "differs from the digest recorded"),
    "provenance misses a partial": ({"prov_hook": _lose_partial}, "partial_digests names"),
    "hashed source digest differs": (
        {"prov_hook": _prov(hashed_source_digest="1" * 64)}, "hashed_source_digest differs"),
    "sidecar file tampered": ({"post": _sidecar_tamper}, "sidecar rejected"),
    "sidecar array hash in the partial": ({"doc_hook": _bad_sidecar_hash}, "sidecar rejected"),
}  # fmt: skip


@pytest.fixture(scope="module")
def steps_c(v7r: ModuleType) -> ModuleType:
    return _load("steps_v5c")


@pytest.fixture(scope="module")
def run_suite_mod() -> ModuleType:
    return _load("run_suite")


def test_rehearsal_loader_accepts_a_valid_fixture(
    tmp_path: Path, steps_c: ModuleType, run_suite_mod: ModuleType
) -> None:
    digest = run_suite_mod.hashed_source_digest("lv5c")
    _build_rehearsal(tmp_path, steps_c, digest)
    blocks, info = load_rehearsal(tmp_path, steps_c, run_suite_mod)
    assert blocks.shape == (900, 20) and info["shards"] == 2 and info["seed_base"] == 20505004
    assert info["host_run_ids"] == ["RUN-SYNTH-1"] and info["hashed_source_digest"] == digest


@pytest.mark.parametrize("case", sorted(VIOLATIONS))
def test_rehearsal_loader_rejects_each_violation(
    case: str, tmp_path: Path, steps_c: ModuleType, run_suite_mod: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
) -> None:  # fmt: skip
    kwargs, expect = VIOLATIONS[case]
    _build_rehearsal(tmp_path, steps_c, run_suite_mod.hashed_source_digest("lv5c"), **kwargs)
    with pytest.raises(RehearsalError) as err:
        load_rehearsal(tmp_path, steps_c, run_suite_mod)
    assert expect in str(err.value), (case, str(err.value))
    # the calibration test turns it into a skip, or a failure in required mode
    monkeypatch.setenv("IONMC_V7R_REHEARSAL_DIR", str(tmp_path))
    monkeypatch.delenv("IONMC_V7R_FIXTURES", raising=False)
    with pytest.raises(pytest.skip.Exception, match="rehearsal fixture rejected"):
        _skip_or_fail(f"rehearsal fixture rejected: {err.value}")
    monkeypatch.setenv("IONMC_V7R_FIXTURES", "required")
    with pytest.raises(pytest.fail.Exception, match="rehearsal fixture rejected"):
        _skip_or_fail(f"rehearsal fixture rejected: {err.value}")


def test_hashed_source_digest_is_deterministic_and_excludes_the_rehearsal_fixtures(
    run_suite_mod: ModuleType, tmp_path: Path
) -> None:
    d = run_suite_mod.hashed_source_digest("lv5c")
    assert d == run_suite_mod.hashed_source_digest("lv5c") and len(d) == 64
    assert d != run_suite_mod.hashed_source_digest("lv5b")
    probe = REHEARSAL_DIR / "PROVENANCE-probe.json"
    assert not probe.exists()
    try:
        REHEARSAL_DIR.mkdir(parents=True, exist_ok=True)
        probe.write_text("{}")
        assert run_suite_mod.hashed_source_digest("lv5c") == d  # fixtures do not enter the digest
    finally:
        probe.unlink(missing_ok=True)
        if REHEARSAL_DIR.exists() and not any(REHEARSAL_DIR.iterdir()):
            REHEARSAL_DIR.rmdir()


def test_rehearsal_loader_missing_directory_is_a_violation(
    tmp_path: Path, steps_c: ModuleType, run_suite_mod: ModuleType
) -> None:
    with pytest.raises(RehearsalError, match="PROVENANCE.json missing"):
        load_rehearsal(tmp_path / "absent", steps_c, run_suite_mod)


def test_no_skip_path_outside_skip_or_fail() -> None:
    """Both calibration files can skip only through ``_skip_or_fail`` (a failure when required)."""
    import re

    pat = re.compile(r"pytest\.skip\(|importorskip|skipif|mark\.skip\b|xfail")
    for f in ("test_v7r_coverage.py", "test_v7_coverage.py"):
        lines = (REPO / "tests" / "ionmc" / f).read_text().splitlines()
        hits = [
            i for i, ln in enumerate(lines, 1) if pat.search(ln) and "pat = re.compile" not in ln
        ]
        if f == "test_v7r_coverage.py":
            hits = [i for i in hits if "pytest.skip(msg)" not in lines[i - 1]]
        else:  # the one guarded skip re-raises in required mode (IONMC_V7R_FIXTURES above it)
            hits = [i for i in hits if "IONMC_V7R_FIXTURES" not in "\n".join(lines[i - 8 : i])]
        assert hits == [], (f, hits)
