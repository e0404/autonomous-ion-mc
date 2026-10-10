"""V7-R: bootstrap-t coverage gates of the escaped neutral energy (V3-005C C5; plan Amendment 14 (h),
Amendment 17 (a)). Pure numpy; no transport, no file I/O (the sidecar and the steps are in steps_v5c.py).

Design (the plan's text with the Amendment 17 (a) readings)
-----------------------------------------------------------
* A replicate is 20 blocks of 500 histories; its data are the 20 block SUMS of the escaped neutral
  energy [MeV] (neutron + gamma). Every gate below is scale invariant (units declaration only).
* Per replicate there are ``B = 1999`` resamples of the 20 block sums (``rng.integers(0, 20, (B, 20))``
  with the default int64 dtype, drawn from a numpy ``PCG64`` generator), ``t*_b = (m*_b - m_hat) /
  se*_b`` and the interval ``[m_hat - t*_(1-a/2) se, m_hat - t*_(a/2) se]`` at ``a = 0.3173``.
  Order statistics (A17 (a)2): the sorted ``t*`` (ascending, 1-based) at ranks
  ``floor((B + 1) a / 2) = 317`` and ``ceil((B + 1)(1 - a / 2)) = 1683``.
  A resample with ``se*_b = 0`` (its 20 drawn blocks equal) gets ``t*_b = +inf`` if ``m*_b > m_hat``,
  ``-inf`` if ``m*_b < m_hat`` and 0 if equal; the count is recorded (A17 (a)3). A replicate with
  ``se = 0`` has no interval and is a MISS in both gates (counted).
* Streams (A17 (a)4): ``default_rng(SeedSequence([seed_base, 16, shard, replicate]))`` for the single
  interval of replicate ``replicate`` of shard ``shard`` and ``SeedSequence([seed_base, 16, 1000 + s,
  j])`` for the pair (replicate j of shard s with replicate j of shard s + 8); in a pair the indices
  of the first replicate are drawn first, then those of the second, from the one generator.
* Paired gate: a pair is a hit iff 0 lies in the closed bootstrap-t interval of ``mu_j - mu_k``
  (``t* = ((m*_j - m*_k) - (m_j - m_k)) / sqrt(se*_j^2 + se*_k^2)``). 16 shards of 450 replicates give
  the 3600 pairs of ``steps_v5b.v7_pair_indices(7200)`` (rows ``i`` and ``i + 3600``). Clopper-Pearson
  at ``alpha_tost = 0.04`` per side; pass iff ``lower >= 0.6527``, ``upper <= 0.7127``, >= 300 pairs.
* Single-interval gate: the Amendment 13 held-out reference. H = shards 12..15 (1800 replicates, the
  last quarter), E = shards 0..11 (5400). ``mu_H``, ``SEM_H`` from the H replicate means; the box
  ``Z_box = max(Z_t, Z_boot)`` (``Z_t`` at ``1 - 0.01/2`` with 1799 df; ``Z_boot`` = ``v7_boot_z`` over
  all 7200 replicate means, 5000 resamples, seed = reference seed). The bootstrap-t interval of an E
  replicate is asymmetric: with ``x_j = ((lo_j + hi_j) / 2 - mu_H) / SEM_H`` and ``half_j = (hi_j - lo_j)
  / (2 SEM_H)`` the closed hit set ``|x_j - delta| <= half_j`` is identical to ``lo_j <= mu_H + delta
  SEM_H <= hi_j``, and ``v7_coverage_extrema`` gives the exact minimum and maximum of the coverage over
  ``delta`` in the box; Clopper-Pearson on ``K_min`` / ``K_max`` at 0.04 per side, same region.

The frozen Amendment 13 functions are imported from ``steps_v5b`` and are not re-implemented. A scale
``kappa`` on the offsets of the interval from ``m_hat`` (calibration margins, as Amendment 13) is
applied by ``interval_from_quantiles`` so that the expensive resampling is done once per replicate.
"""

from __future__ import annotations

import math
import platform
from typing import Any

import numpy as np
from numpy.typing import NDArray

import steps_v5b as v5b

V7R_R_INDEX = 16
V7R_DIAG_R_INDEX = 17
V7R_SHARDS = 16
V7R_SHARD_N = 4_500_000
V7R_BATCH_N = 500
V7R_BLOCKS_PER_REP = 20
V7R_REPS_PER_SHARD = 450
V7R_REPLICATES = V7R_SHARDS * V7R_REPS_PER_SHARD  # 7200
V7R_PAIRS = V7R_REPLICATES // 2  # 3600
V7R_REF_INDEX = 16  # the reference is "shard" 16 of r_index 16
V7R_REF_N = 1_000_000
V7R_REF_BLOCKS = 20
V7R_SEED_FIRST = 20481004 + 1000 * V7R_R_INDEX  # 20497004: seed of shard k is this + k
V7R_REF_SEED = V7R_SEED_FIRST + V7R_REF_INDEX  # 20497020
V7R_DIAG_SEED_FIRST = 20481004 + 1000 * V7R_DIAG_R_INDEX  # 20498004
V7R_HELDOUT_SHARDS = (12, 13, 14, 15)
V7R_B = 1999
V7R_ALPHA = 0.3173
V7R_REGION = (0.6527, 0.7127)  # 0.6827 +- 0.03, both gates
V7R_RANK_LO = math.floor((V7R_B + 1) * V7R_ALPHA / 2.0)  # 317 (1-based)
V7R_RANK_HI = math.ceil((V7R_B + 1) * (1.0 - V7R_ALPHA / 2.0))  # 1683 (1-based)
V7R_AGG_B = (5, 10, 20)  # blocks per replicate of the exploratory diagnostic (exact aggregation)
V7R_UNITS = "MeV per 500-history block"
V7R_QUANTILE_RULE = (
    "sorted t*, ascending, 1-based: lower = floor((B+1)*alpha/2)-th, upper = ceil((B+1)*(1-alpha/2))-th "
    "(Amendment 17 (a)2); se*=0 resample: t* = +inf/-inf by sign(m*-m_hat), 0 if equal (A17 (a)3)"
)
SEED_RULE = (
    "single: SeedSequence([seed_base, 16, shard, replicate]); pair: SeedSequence([seed_base, 16, "
    "1000 + s, j]); numpy default_rng (PCG64), rng.integers(0, 20, (B, 20)) int64, pair: replicate j "
    "indices first, then replicate k (Amendment 17 (a)4)"
)
_BIG = 1e100  # finite stand-in for an infinite interval end in the exact sweep (cannot be hit-neutral)


def shard_seed(k: int) -> int:
    """Seed of shard ``k`` (0..15) of r_index 16 at the evidence base, ``base + 16000 + k``."""
    if not 0 <= k <= V7R_REF_INDEX:
        raise ValueError("shard index outside 0..16")
    return V7R_SEED_FIRST + k


def single_stream(seed_base: int, shard: int, replicate: int) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence([seed_base, V7R_R_INDEX, shard, replicate])
    )


def pair_stream(seed_base: int, s: int, j: int) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence([seed_base, V7R_R_INDEX, 1000 + s, j])
    )


def rng_record() -> dict[str, str]:
    return {"numpy": np.__version__, "bit_generator": "PCG64", "python": platform.python_version()}


def _ranks(b: int, alpha: float) -> tuple[int, int]:
    return math.floor((b + 1) * alpha / 2.0), math.ceil((b + 1) * (1.0 - alpha / 2.0))


# -- single bootstrap-t ------------------------------------------------------------------------
def _tstar(xs: NDArray[np.float64], m_hat: float) -> tuple[NDArray[np.float64], int]:
    """Resampled studentised statistic ``t*`` of the rows of ``xs`` ``[B, n]`` about ``m_hat``."""
    n = xs.shape[1]
    degenerate = xs.max(axis=1) == xs.min(axis=1)
    m = xs.mean(axis=1)
    m = np.where(degenerate, xs[:, 0], m)
    se = xs.std(axis=1, ddof=1) / math.sqrt(n)
    num = m - m_hat
    with np.errstate(divide="ignore", invalid="ignore"):
        t = num / se
    with np.errstate(invalid="ignore"):
        t = np.where(degenerate, np.sign(num) * np.inf, t)  # sign(0) * inf = nan, fixed next
    t = np.where(degenerate & (num == 0.0), 0.0, t)
    return t, int(degenerate.sum())


def boot_t_quantiles(
    sums: NDArray[np.float64], rng: np.random.Generator, b: int = V7R_B, alpha: float = V7R_ALPHA
) -> tuple[float, float, float, float, int]:
    """``(m_hat, se, t_lo, t_hi, n_degenerate)`` of one replicate: sample mean, standard error of
    its blocks and the two order statistics of the ``b`` resampled ``t*`` (``se = 0``: ``nan``
    quantiles, the replicate has no interval). The generator is advanced by one ``(b, n)`` draw only
    for a non-degenerate replicate (a degenerate replicate is a miss in either case)."""
    x = np.asarray(sums, dtype=np.float64)
    n = x.size
    m_hat = float(x.mean())
    if x.max() == x.min():
        return m_hat, 0.0, math.nan, math.nan, -1
    se = float(x.std(ddof=1) / math.sqrt(n))
    idx = rng.integers(0, n, size=(b, n))
    t, ndeg = _tstar(x[idx], m_hat)
    r_lo, r_hi = _ranks(b, alpha)
    part = np.partition(t, [r_lo - 1, r_hi - 1])
    return m_hat, se, float(part[r_lo - 1]), float(part[r_hi - 1]), ndeg


def interval_from_quantiles(
    m_hat: float, se: float, t_lo: float, t_hi: float, kappa: float = 1.0
) -> tuple[float, float]:
    """``[m_hat - kappa t_hi se, m_hat - kappa t_lo se]`` (``kappa`` scales the offsets from the point
    estimate; 1 is the gate). Infinite ends are returned as +-1e100 (finite, for the sweep)."""
    lo, hi = m_hat - kappa * t_hi * se, m_hat - kappa * t_lo * se
    return float(np.clip(lo, -_BIG, _BIG)), float(np.clip(hi, -_BIG, _BIG))


def boot_t_interval(
    sums: NDArray[np.float64], rng: np.random.Generator, b: int = V7R_B, alpha: float = V7R_ALPHA
) -> tuple[float, float, float, float, dict[str, Any]]:
    """Bootstrap-t interval ``(lo, hi, m_hat, se, diag)`` of the mean of the block values ``sums``
    (module docstring); ``se = 0`` gives ``lo = hi = nan`` (a miss) and ``diag["degenerate"]``."""
    m_hat, se, t_lo, t_hi, ndeg = boot_t_quantiles(sums, rng, b, alpha)
    diag = {"n_resamples_degenerate": max(ndeg, 0), "degenerate": se == 0.0, "t_lo": t_lo, "t_hi": t_hi}
    if se == 0.0:
        return math.nan, math.nan, m_hat, se, diag
    lo, hi = interval_from_quantiles(m_hat, se, t_lo, t_hi)
    return lo, hi, m_hat, se, diag


# -- pair bootstrap-t --------------------------------------------------------------------------
def boot_t_pair_quantiles(
    sums_j: NDArray[np.float64],
    sums_k: NDArray[np.float64],
    rng: np.random.Generator,
    b: int = V7R_B,
    alpha: float = V7R_ALPHA,
) -> tuple[float, float, float, float, int]:
    """``(d_hat, se_d, t_lo, t_hi, n_degenerate)`` of the pair: ``d_hat = m_j - m_k``, ``se_d =
    sqrt(se_j^2 + se_k^2)`` and the order statistics of ``t* = ((m*_j - m*_k) - d_hat) /
    sqrt(se*_j^2 + se*_k^2)``; indices of ``j`` are drawn first, then ``k``. A resample with both
    ``se* = 0`` follows the ``t*`` rule of the module docstring by the sign of its numerator."""
    xj, xk = (np.asarray(a, dtype=np.float64) for a in (sums_j, sums_k))
    n = xj.size
    mj, mk = float(xj.mean()), float(xk.mean())
    sej = float(xj.std(ddof=1) / math.sqrt(n))
    sek = float(xk.std(ddof=1) / math.sqrt(xk.size))
    d_hat, se_d = mj - mk, math.hypot(sej, sek)
    if sej == 0.0 or sek == 0.0:
        return d_hat, se_d, math.nan, math.nan, -1  # replicate-level degenerate se: a miss (counted)
    ij = rng.integers(0, n, size=(b, n))
    ik = rng.integers(0, xk.size, size=(b, xk.size))
    aj, ak = xj[ij], xk[ik]
    dj, dk = aj.max(axis=1) == aj.min(axis=1), ak.max(axis=1) == ak.min(axis=1)
    mjs = np.where(dj, aj[:, 0], aj.mean(axis=1))
    mks = np.where(dk, ak[:, 0], ak.mean(axis=1))
    sj = np.where(dj, 0.0, aj.std(axis=1, ddof=1) / math.sqrt(n))
    sk = np.where(dk, 0.0, ak.std(axis=1, ddof=1) / math.sqrt(xk.size))
    num = (mjs - mks) - d_hat
    den = np.sqrt(sj * sj + sk * sk)
    zero = den == 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        t = num / den
    with np.errstate(invalid="ignore"):
        t = np.where(zero, np.sign(num) * np.inf, t)  # sign(0) * inf = nan, fixed next
    t = np.where(zero & (num == 0.0), 0.0, t)
    r_lo, r_hi = _ranks(b, alpha)
    part = np.partition(t, [r_lo - 1, r_hi - 1])
    return d_hat, se_d, float(part[r_lo - 1]), float(part[r_hi - 1]), int(zero.sum())


def boot_t_pair_interval(
    sums_j: NDArray[np.float64],
    sums_k: NDArray[np.float64],
    rng: np.random.Generator,
    b: int = V7R_B,
    alpha: float = V7R_ALPHA,
) -> tuple[float, float, float, float, dict[str, Any]]:
    """Bootstrap-t interval ``(lo, hi, d_hat, se_d, diag)`` for ``mu_j - mu_k``; a hit iff ``0`` is in
    it. ``se_j = 0`` or ``se_k = 0`` gives ``lo = hi = nan`` (a miss, ``diag["degenerate"]``)."""
    d_hat, se_d, t_lo, t_hi, ndeg = boot_t_pair_quantiles(sums_j, sums_k, rng, b, alpha)
    diag = {"n_resamples_degenerate": max(ndeg, 0), "degenerate": ndeg < 0, "t_lo": t_lo, "t_hi": t_hi}
    if ndeg < 0:
        return math.nan, math.nan, d_hat, se_d, diag
    lo, hi = interval_from_quantiles(d_hat, se_d, t_lo, t_hi)
    return lo, hi, d_hat, se_d, diag


# -- ingredients and gates ---------------------------------------------------------------------
def pair_ingredients(
    sums: NDArray[np.float64], seed_base: int, reps_per_shard: int = V7R_REPS_PER_SHARD,
    b: int = V7R_B, alpha: float = V7R_ALPHA,
) -> dict[str, NDArray[np.float64]]:
    """Resampling of every pair ``p`` (rows ``p`` and ``p + R/2`` of ``sums`` ``[R, 20]``; shard ``s =
    p // reps_per_shard``, replicate ``j = p % reps_per_shard``): ``d_hat``, ``se_d``, ``t_lo``,
    ``t_hi``, ``n_degenerate`` (resamples, -1 for a degenerate-replicate pair)."""
    ia, ib = v5b.v7_pair_indices(sums.shape[0])
    out = np.empty((ia.size, 5))
    for p in range(ia.size):
        rng = pair_stream(seed_base, p // reps_per_shard, p % reps_per_shard)
        out[p] = boot_t_pair_quantiles(sums[ia[p]], sums[ib[p]], rng, b, alpha)
    return {k: out[:, i] for i, k in enumerate(("d_hat", "se_d", "t_lo", "t_hi", "n_degenerate"))}


def single_ingredients(
    sums: NDArray[np.float64], seed_base: int, reps_per_shard: int = V7R_REPS_PER_SHARD,
    b: int = V7R_B, alpha: float = V7R_ALPHA,
) -> dict[str, NDArray[np.float64]]:
    """Resampling of every evaluation replicate (rows ``0 .. n_E - 1``; shard ``r // reps_per_shard``,
    replicate ``r % reps_per_shard``): ``m_hat``, ``se``, ``t_lo``, ``t_hi``, ``n_degenerate`` (-1: se = 0)."""
    n_e = v5b.v7_heldout_split(sums.shape[0])
    out = np.empty((n_e, 5))
    for r in range(n_e):
        rng = single_stream(seed_base, r // reps_per_shard, r % reps_per_shard)
        out[r] = boot_t_quantiles(sums[r], rng, b, alpha)
    return {k: out[:, i] for i, k in enumerate(("m_hat", "se", "t_lo", "t_hi", "n_degenerate"))}


def pair_gate(ing: dict[str, NDArray[np.float64]], kappa: float = 1.0, *,
              region: tuple[float, float] = V7R_REGION) -> dict[str, Any]:  # fmt: skip
    """Paired gate from ``pair_ingredients`` at interval scale ``kappa``."""
    d, se, tl, th = ing["d_hat"], ing["se_d"], ing["t_lo"], ing["t_hi"]
    degenerate = ing["n_degenerate"] < 0
    with np.errstate(invalid="ignore"):
        lo, hi = d - kappa * th * se, d - kappa * tl * se
    hit = (lo <= 0.0) & (0.0 <= hi) & ~degenerate  # nan compares False: degenerate pairs are misses
    pairs, k = int(hit.size), int(hit.sum())
    lower = v5b.cp_lower(k, pairs, v5b.V7_ALPHA_TOST)
    upper = v5b.cp_upper(k, pairs, v5b.V7_ALPHA_TOST)
    reasons = []
    if lower < region[0]:
        reasons.append(f"CI lower bound {lower:.4f} below {region[0]:.4f}")
    if upper > region[1]:
        reasons.append(f"CI upper bound {upper:.4f} above {region[1]:.4f}")
    if pairs < v5b.V7_MIN_INTERVALS:
        reasons.append(f"pairs {pairs} < {v5b.V7_MIN_INTERVALS}")
    return {"n_pairs": pairs, "covered": k, "m": k / pairs, "ci": [lower, upper], "low": region[0],
            "high": region[1], "alpha_tost": v5b.V7_ALPHA_TOST, "bound_kind": "clopper-pearson",
            "degenerate_pairs": int(degenerate.sum()),
            "n_resamples_degenerate": int(np.maximum(ing["n_degenerate"], 0).sum()),
            "per_pair_hit": [int(h) for h in hit], "pass": not reasons, "reasons": reasons}  # fmt: skip


def single_gate(ing: dict[str, NDArray[np.float64]], means: NDArray[np.float64], kappa: float = 1.0, *,
                boot_seed: int = 0, n_boot: int | None = None, z_boot_override: float | None = None,
                region: tuple[float, float] = V7R_REGION) -> dict[str, Any]:  # fmt: skip
    """Single-interval gate from ``single_ingredients`` and ALL replicate means ``means`` ``[R]`` (E
    first, then H; module docstring). ``z_boot_override`` is for the calibration only. ``n_boot=None``
    is the frozen qualification value ``steps_v5b.V7_BOOT_N``, read at call time (the calibration uses
    the same default, so the two cannot differ)."""
    n_boot = v5b.V7_BOOT_N if n_boot is None else n_boot
    r = means.size
    n_e = v5b.v7_heldout_split(r)
    n_h = r - n_e
    if ing["m_hat"].size != n_e:
        raise ValueError("single_gate: ingredients do not cover the evaluation replicates")
    mh = means[n_e:]
    mu_h, sem_h = float(mh.mean()), float(mh.std(ddof=1) / math.sqrt(n_h))
    z_t = v5b.student_t_quantile(1.0 - v5b.V7_ALPHA_BOX / 2.0, n_h - 1)
    if z_boot_override is None:
        z_boot = v5b.v7_boot_z(means[:, None], n_h, np.array([means.mean()]), np.array([max(sem_h, 1e-300)]),
                               boot_seed, n_boot)  # fmt: skip
    else:
        z_boot = z_boot_override
    z_box = max(z_t, float(z_boot))
    degenerate = ing["n_degenerate"] < 0
    ok = ~degenerate & (sem_h > 0.0)
    nan = math.nan
    k_min = k_max = -1
    lower = upper = nan
    if ok.any():
        lo, hi = interval_from_quantiles_array(ing, kappa)
        x = ((lo + hi) / 2.0 - mu_h) / sem_h
        half = (hi - lo) / (2.0 * sem_h)
        c_min, c_max, _, _ = v5b.v7_coverage_extrema(x[ok], half[ok], z_box, n_e)
        k_min, k_max = int(round(c_min * n_e)), int(round(c_max * n_e))
        lower = v5b.cp_lower(k_min, n_e, v5b.V7_ALPHA_TOST)
        upper = v5b.cp_upper(k_max, n_e, v5b.V7_ALPHA_TOST)
    reasons = []
    if not lower >= region[0]:
        reasons.append(f"CI lower bound {lower:.4f} below {region[0]:.4f}")
    if not upper <= region[1]:
        reasons.append(f"CI upper bound {upper:.4f} above {region[1]:.4f}")
    if n_e < v5b.V7_MIN_INTERVALS:
        reasons.append(f"replicates {n_e} < {v5b.V7_MIN_INTERVALS}")
    return {"replicates_eval": n_e, "replicates_heldout": n_h, "mu_heldout": mu_h, "sem_heldout": sem_h,
            "alpha_box": v5b.V7_ALPHA_BOX, "boot_seed": boot_seed, "boot_n": n_boot, "z_t": z_t,
            "z_boot": float(z_boot), "z_box": z_box, "cp_count_min": k_min, "cp_count_max": k_max,
            "ci": [lower, upper], "low": region[0], "high": region[1], "alpha_tost": v5b.V7_ALPHA_TOST,
            "bound_kind": "clopper-pearson", "degenerate_replicates": int(degenerate.sum()),
            "n_resamples_degenerate": int(np.maximum(ing["n_degenerate"], 0).sum()),
            "pass": not reasons, "reasons": reasons}  # fmt: skip


def interval_from_quantiles_array(ing: dict[str, NDArray[np.float64]], kappa: float
                                  ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:  # fmt: skip
    """Vectorised ``interval_from_quantiles`` (non-degenerate rows meaningful, others nan)."""
    with np.errstate(invalid="ignore"):
        lo = ing["m_hat"] - kappa * ing["t_hi"] * ing["se"]
        hi = ing["m_hat"] - kappa * ing["t_lo"] * ing["se"]
    return np.clip(lo, -_BIG, _BIG), np.clip(hi, -_BIG, _BIG)


def escaped_neutral_verdict(
    sums: NDArray[np.float64], seed_base: int, *, boot_seed: int, ref_used: bool = True,
    reps_per_shard: int = V7R_REPS_PER_SHARD, b: int = V7R_B, alpha: float = V7R_ALPHA,
    n_boot: int | None = None,
) -> dict[str, Any]:  # fmt: skip
    """Both bootstrap-t gates of ``escaped_neutral`` on the block sums ``sums`` ``[R, 20]`` (rows in
    shard order); the document fields of ``step_v7r_combine``. ``ref_used`` is the bin mask of the
    independent reference (reference mean > 0 and block spread > 0), as for the frozen gates."""
    sums = np.asarray(sums, dtype=np.float64)
    if sums.ndim != 2 or sums.shape[1] != V7R_BLOCKS_PER_REP or not np.all(np.isfinite(sums)):
        raise SystemExit("v7r: block sums must be a finite [replicates, 20] array")
    pi = pair_ingredients(sums, seed_base, reps_per_shard, b, alpha)
    si = single_ingredients(sums, seed_base, reps_per_shard, b, alpha)
    paired = pair_gate(pi)
    single = single_gate(si, sums.mean(axis=1), boot_seed=boot_seed, n_boot=n_boot)
    if not ref_used:
        paired["pass"], single["pass"] = False, False
        paired["reasons"].append("reference bin unusable")
        single["reasons"].append("reference bin unusable")
    return {
        "kind": "scalar", "interval": "bootstrap-t", "B": b, "alpha": alpha,
        "quantile_rule": V7R_QUANTILE_RULE, "ranks": list(_ranks(b, alpha)),
        "seed_rule": SEED_RULE, "seed_base": seed_base, "rng": rng_record(), "units": V7R_UNITS,
        "n_pairs": paired["n_pairs"], "covered": paired["covered"], "m": paired["m"],
        "ci": paired["ci"], "low": paired["low"], "high": paired["high"],
        "n_resamples_degenerate": paired["n_resamples_degenerate"]
        + single["n_resamples_degenerate"],
        "n_resamples_degenerate_pairs": paired["n_resamples_degenerate"],
        "n_resamples_degenerate_single": single["n_resamples_degenerate"],
        "degenerate_pairs": paired["degenerate_pairs"],
        "per_pair_hit": paired.pop("per_pair_hit"), "paired": paired, "single": single,
        "pass_paired": paired["pass"], "pass_single": single["pass"],
        "reasons": paired["reasons"] + [f"single: {x}" for x in single["reasons"]],
        "pass": bool(paired["pass"] and single["pass"]),
    }  # fmt: skip


# -- exploratory diagnostic (r_index 17): exact aggregation of the recorded blocks ----------------
def aggregate_blocks(sums: NDArray[np.float64], b: int) -> NDArray[np.float64]:
    """Replicates of ``b`` blocks (of ``20 / b`` recorded blocks each) from the recorded block sums
    ``sums`` ``[R, 20]`` by exact aggregation: ``b = 20`` is the identity, ``b = 10`` sums adjacent
    pairs and ``b = 5`` adjacent groups of four (``math.fsum``, the correctly rounded exact sum).
    The result is ``[R, b]``; every replicate still covers the same 1e4 histories."""
    if b not in V7R_AGG_B:
        raise ValueError(f"b must be one of {V7R_AGG_B}")
    x = np.asarray(sums, dtype=np.float64)
    g = V7R_BLOCKS_PER_REP // b
    grouped = x.reshape(x.shape[0], b, g)
    return np.array([[math.fsum(c) for c in row] for row in grouped]).reshape(x.shape[0], b)


def excess_kurtosis(x: NDArray[np.float64]) -> float:
    """Excess kurtosis (population moments) of a sample; ``nan`` for zero variance."""
    x = np.asarray(x, dtype=np.float64).ravel()
    v = float(np.mean((x - x.mean()) ** 2))
    return math.nan if v == 0.0 else float(np.mean((x - x.mean()) ** 4) / v**2 - 3.0)
