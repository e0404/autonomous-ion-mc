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

import hashlib
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


def _alpha_tost(alpha_tost: float | None) -> float:
    """Per-side TOST level of a bootstrap-t gate. ``None`` (every archived suite) is the frozen
    ``steps_v5b.V7_ALPHA_TOST`` (0.04), read at call time exactly as before the parameter existed;
    V7-R2 passes ``V7R2_ALPHA_TOST`` (0.005) explicitly (repair plan section 6.5)."""
    a = v5b.V7_ALPHA_TOST if alpha_tost is None else float(alpha_tost)
    if not 0.0 < a < 0.5:
        raise ValueError(f"alpha_tost {a!r} outside (0, 0.5)")
    return a


def pair_gate(ing: dict[str, NDArray[np.float64]], kappa: float = 1.0, *,
              region: tuple[float, float] = V7R_REGION,
              alpha_tost: float | None = None) -> dict[str, Any]:  # fmt: skip
    """Paired gate from ``pair_ingredients`` at interval scale ``kappa`` (``alpha_tost``: per-side level,
    default the frozen 0.04)."""
    a_tost = _alpha_tost(alpha_tost)
    d, se, tl, th = ing["d_hat"], ing["se_d"], ing["t_lo"], ing["t_hi"]
    degenerate = ing["n_degenerate"] < 0
    with np.errstate(invalid="ignore"):
        lo, hi = d - kappa * th * se, d - kappa * tl * se
    hit = (lo <= 0.0) & (0.0 <= hi) & ~degenerate  # nan compares False: degenerate pairs are misses
    pairs, k = int(hit.size), int(hit.sum())
    lower = v5b.cp_lower(k, pairs, a_tost)
    upper = v5b.cp_upper(k, pairs, a_tost)
    reasons = []
    if lower < region[0]:
        reasons.append(f"CI lower bound {lower:.4f} below {region[0]:.4f}")
    if upper > region[1]:
        reasons.append(f"CI upper bound {upper:.4f} above {region[1]:.4f}")
    if pairs < v5b.V7_MIN_INTERVALS:
        reasons.append(f"pairs {pairs} < {v5b.V7_MIN_INTERVALS}")
    return {"n_pairs": pairs, "covered": k, "m": k / pairs, "ci": [lower, upper], "low": region[0],
            "high": region[1], "alpha_tost": a_tost, "bound_kind": "clopper-pearson",
            "degenerate_pairs": int(degenerate.sum()),
            "n_resamples_degenerate": int(np.maximum(ing["n_degenerate"], 0).sum()),
            "per_pair_hit": [int(h) for h in hit], "pass": not reasons, "reasons": reasons}  # fmt: skip


def single_gate(ing: dict[str, NDArray[np.float64]], means: NDArray[np.float64], kappa: float = 1.0, *,
                boot_seed: int = 0, n_boot: int | None = None, z_boot_override: float | None = None,
                region: tuple[float, float] = V7R_REGION,
                alpha_tost: float | None = None) -> dict[str, Any]:  # fmt: skip
    """Single-interval gate from ``single_ingredients`` and ALL replicate means ``means`` ``[R]`` (E
    first, then H; module docstring). ``z_boot_override`` is for the calibration only. ``n_boot=None``
    is the frozen qualification value ``steps_v5b.V7_BOOT_N``, read at call time (the calibration uses
    the same default, so the two cannot differ). ``alpha_tost``: per-side level of the Clopper-Pearson
    bounds on the count extrema, default the frozen 0.04 (the box level ``V7_ALPHA_BOX`` is unchanged)."""
    n_boot = v5b.V7_BOOT_N if n_boot is None else n_boot
    a_tost = _alpha_tost(alpha_tost)
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
        lower = v5b.cp_lower(k_min, n_e, a_tost)
        upper = v5b.cp_upper(k_max, n_e, a_tost)
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
            "ci": [lower, upper], "low": region[0], "high": region[1], "alpha_tost": a_tost,
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
    n_boot: int | None = None, alpha_tost: float | None = None,
) -> dict[str, Any]:  # fmt: skip
    """Both bootstrap-t gates of ``escaped_neutral`` on the block sums ``sums`` ``[R, 20]`` (rows in
    shard order); the document fields of ``step_v7r_combine``. ``ref_used`` is the bin mask of the
    independent reference (reference mean > 0 and block spread > 0), as for the frozen gates."""
    sums = np.asarray(sums, dtype=np.float64)
    if sums.ndim != 2 or sums.shape[1] != V7R_BLOCKS_PER_REP or not np.all(np.isfinite(sums)):
        raise SystemExit("v7r: block sums must be a finite [replicates, 20] array")
    pi = pair_ingredients(sums, seed_base, reps_per_shard, b, alpha)
    si = single_ingredients(sums, seed_base, reps_per_shard, b, alpha)
    paired = pair_gate(pi, alpha_tost=alpha_tost)
    single = single_gate(si, sums.mean(axis=1), boot_seed=boot_seed, n_boot=n_boot,
                         alpha_tost=alpha_tost)  # fmt: skip
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


# -- calibration sharding and verdict (plan Amendment 14a) -----------------------------------------
V7R_CAL_REPS = 3200
"""Monte Carlo replicates per surrogate and margin of the ``escaped_neutral`` calibration (Amendment 14a;
``run_suite.V7R_CAL_REPS`` is asserted equal by a test and exported to the steps as ``IONMC_V7R_CAL_REPS``)."""
CAL_SURROGATES = ("archive_fitted", "rehearsal_block_sums")
CAL_MARGINS = ("nominal", "paired_low", "paired_high", "single_low", "single_high")
CAL_PARTS = 15
CAL_SEED = 20507000
CAL_SCHEMA = "ionmc-v7r-calibration-record-1"
CAL_GATE_LIMIT = 0.04  # each gate alone: CP95 UPPER of its Monte Carlo rate at its own margin <= 0.04
CAL_JOINT_LIMIT = 0.02
CAL_POWER_MIN = 0.9


def cal_part_sizes(reps: int, parts: int = CAL_PARTS) -> tuple[int, ...]:
    """Replicates of each part: ``reps // parts`` each, the remainder one more to the first parts
    (3200 -> 214 for parts 0..4 and 213 for parts 5..14)."""
    q, r = divmod(reps, parts)
    return tuple(q + (1 if i < r else 0) for i in range(parts))


def cal_part_seeds(surrogate_idx: int, part: int) -> dict[str, int]:
    """Seeds of the data and of the bootstrap of one calibration part: the two 32-bit words generated by
    ``SeedSequence([20507000, surrogate_idx, part])``; replicate ``r`` of the part uses ``data + r`` (data
    stream ``default_rng``) and ``boot + r`` (the ``seed_base`` of the bootstrap streams). Distinct parts and
    surrogates draw independent replicates; the five margins of a part share each replicate's data and
    bootstrap (one bootstrap per replicate serves every kappa)."""
    w = np.random.SeedSequence([CAL_SEED, surrogate_idx, part]).generate_state(2, dtype=np.uint32)
    return {"data": int(w[0]), "boot": int(w[1])}


def _own_gates(margin: str) -> tuple[str, ...]:
    """Gates whose construction requirement applies at ``margin``: the paired gate at ``paired_*``, the
    single gate at ``single_*`` (the other gate's values there are reported only); a name with neither
    prefix (synthetic tests) names both."""
    for gate in ("paired", "single"):
        if margin.startswith(gate):
            return (gate,)
    return ("paired", "single")


def calibration_evaluate(counts: dict[str, dict[str, int]]) -> dict[str, Any]:
    """Verdict of one surrogate from the accepted counts per margin (``{margin: {n, joint, paired, single}}``).

    * joint false acceptance: the one-sided 95 % Clopper-Pearson UPPER bound <= 0.02 at every margin but
      ``nominal``;
    * power: the 95 % CP LOWER bound of the joint rate at ``nominal`` (kappa = 1) >= 0.9;
    * gate alone (Amendment 14a (b)4, as implemented before the rehearsal): the CP95 UPPER bound of the
      Monte Carlo rate of each gate at its OWN margins <= 0.04 (n = 3200: at most 109 accepted). The
      "construction reading" (the gate is refuted iff the CP95 LOWER bound exceeds 0.04) is reported as
      ``gate_alone_construction_reading`` with its verdict and never gates.

    Bounds at n = 3200 (``steps_v5b.cp_upper/cp_lower``): joint max 50 accepted (observed 0.0156, CP95
    upper 0.0197; 51 gives 0.0201), gate alone max 109 (upper 0.0398; 110 gives 0.0402), power min 2909
    (lower 0.9003; 2908 gives 0.8999). Table of n = 200..3200 in ``tests/ionmc/test_v7r_coverage.py``."""
    margins: dict[str, Any] = {}
    failures: list[str] = []
    for name, c in counts.items():
        n = int(c["n"])
        rep: dict[str, Any] = {"n": n}
        for key in ("joint", "paired", "single"):
            k = int(c[key])
            rep[key] = k
            rep[f"{key}_rate"] = k / n
            rep[f"{key}_cp95_lower"] = v5b.cp_lower(k, n, 0.05)
            rep[f"{key}_cp95_upper"] = v5b.cp_upper(k, n, 0.05)
        if name == "nominal":
            if rep["joint_cp95_lower"] < CAL_POWER_MIN:
                failures.append(
                    f"power: joint CP95 lower {rep['joint_cp95_lower']:.4f} < {CAL_POWER_MIN} "
                    f"at the nominal ({rep['joint']}/{n})"
                )
        else:
            if rep["joint_cp95_upper"] > CAL_JOINT_LIMIT:
                failures.append(
                    f"{name}: joint false acceptance CP95 upper {rep['joint_cp95_upper']:.4f} > "
                    f"{CAL_JOINT_LIMIT} ({rep['joint']}/{n})"
                )
            alone = {}
            reading = {}
            for gate in _own_gates(name):
                upper = rep[f"{gate}_cp95_upper"]
                lower = rep[f"{gate}_cp95_lower"]
                ok = upper <= CAL_GATE_LIMIT
                refuted = lower > CAL_GATE_LIMIT
                alone[gate] = {
                    "rate": rep[f"{gate}_rate"],
                    "cp95_lower": lower,
                    "cp95_upper": upper,
                    "cp95_upper_le_0.04": ok,
                }
                reading[gate] = {
                    "construction_refuted": refuted,
                    "verdict": "refuted (CP95 lower > 0.04)" if refuted else "not refuted",
                }
                if not ok:
                    failures.append(
                        f"{name}: {gate} gate alone: CP95 upper {upper:.4f} > {CAL_GATE_LIMIT} "
                        f"({rep[gate]}/{n})"
                    )
            rep["gate_alone"] = alone
            rep["gate_alone_construction_reading"] = reading
        margins[name] = rep
    return {"margins": margins, "failures": failures, "pass": not failures}


def calibration_combine(
    records: list[dict[str, Any]], reps: int = V7R_CAL_REPS, *, parts: int = CAL_PARTS,
    seed_fn: Any = None, kappa_record_sha256: str | None = None, kappa_fits: dict[str, Any] | None = None,
    layout: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Pooled verdict of the calibration shards. ``records`` are the sealed shard records (one per
    (surrogate, part), exactly once; fail closed with ``SystemExit`` otherwise), each holding the counts and the
    kappa of the five margins (``margins[m] = {kappa, accepted_joint, accepted_paired, accepted_single}``):
    identical fit across the parts of a surrogate and identical kappa per margin, the derived seeds, ``n`` per
    part and ``reps_env`` equal to the expected ones, ``n`` summing to ``reps``. Counts are summed per
    (surrogate, margin) and judged by :func:`calibration_evaluate`.

    Defaults are the V7-R layout (``parts`` 15, :func:`cal_part_seeds`, no kappa record), unchanged. V7-R2 passes
    ``parts=V7R2_CAL_PARTS``, ``seed_fn=r2_cal_part_seeds`` and the ``content_sha256`` of the sealed combined kappa
    record: every part must carry exactly that ``kappa_record_sha256`` (a part without it, or with another, fails
    the calibration closed; plan 6.7 item 6). ``kappa_fits`` ({surrogate: fit block} of the combined kappa record) makes
    every part's surrogate ``fit`` block (family / parameters or fixture digest and sidecar hashes) equal the one the
    kappa step used, fail closed. ``layout`` (V7-R2: ``r2_part_layout()``) must equal every part's recorded ``layout``
    (replicates, replicates per shard, parts), so a part produced for another design is rejected."""
    seed_fn = cal_part_seeds if seed_fn is None else seed_fn
    sizes = cal_part_sizes(reps, parts)
    by: dict[tuple[Any, Any], dict[str, Any]] = {}
    for r in records:
        key = (r.get("surrogate"), r.get("part"))
        if key in by:
            raise SystemExit(f"calibration records: duplicate {key}")
        by[key] = r
    want = {(s, p) for s in CAL_SURROGATES for p in range(parts)}
    if set(by) != want:
        raise SystemExit(f"calibration records: missing {sorted(want - set(by))[:3]} / "
                         f"unexpected {sorted(map(str, set(by) - want))[:3]}")  # fmt: skip
    out: dict[str, Any] = {}
    for si, s in enumerate(CAL_SURROGATES):
        fit0 = by[(s, 0)]["fit"]
        kappas = {m: by[(s, 0)]["margins"][m]["kappa"] for m in CAL_MARGINS}
        tot = {m: {"n": 0, "joint": 0, "paired": 0, "single": 0} for m in CAL_MARGINS}
        for p in range(parts):
            r = by[(s, p)]
            tag = f"{s}/part{p}"
            if r.get("reps_env") != reps:
                raise SystemExit(
                    f"calibration record {tag}: reps_env {r.get('reps_env')!r} is not {reps}"
                )
            if r.get("fit") != fit0:
                raise SystemExit(f"calibration record {tag}: fit differs within the surrogate")
            if r.get("n") != sizes[p]:
                raise SystemExit(f"calibration record {tag}: n {r.get('n')!r} is not {sizes[p]}")
            if r.get("seeds") != seed_fn(si, p):
                raise SystemExit(f"calibration record {tag}: seeds differ from the derived ones")
            if layout is not None and r.get("layout") != layout:
                raise SystemExit(f"calibration record {tag}: layout {r.get('layout')!r} is not {layout!r}")
            if kappa_fits is not None and r.get("fit") != kappa_fits.get(s):
                raise SystemExit(f"calibration record {tag}: fit block differs from the kappa record's")
            if kappa_record_sha256 is not None and r.get("kappa_record_sha256") != kappa_record_sha256:
                raise SystemExit(f"calibration record {tag}: kappa_record_sha256 differs from the kappa record")
            ms = r.get("margins")
            if not isinstance(ms, dict) or set(ms) != set(CAL_MARGINS):
                raise SystemExit(f"calibration record {tag}: margins are not exactly {CAL_MARGINS}")
            for m in CAL_MARGINS:
                if ms[m].get("kappa") != kappas[m]:
                    raise SystemExit(
                        f"calibration record {tag}: kappa differs within the margin {m}"
                    )
                for key in ("joint", "paired", "single"):
                    k = ms[m].get(f"accepted_{key}")
                    if not isinstance(k, int) or isinstance(k, bool) or not 0 <= k <= r["n"]:
                        raise SystemExit(
                            f"calibration record {tag}/{m}: accepted_{key} {k!r} invalid"
                        )
                    tot[m][key] += k
                tot[m]["n"] += r["n"]
        for m in CAL_MARGINS:
            if tot[m]["n"] != reps:
                raise SystemExit(f"calibration {s}/{m}: n sums to {tot[m]['n']}, not {reps}")
        out[s] = {**calibration_evaluate(tot), "kappas": kappas, "fit": fit0}
    return {"surrogates": out, "reps": reps, "part_sizes": list(sizes),
            "pass": all(v["pass"] for v in out.values())}  # fmt: skip


# -- V7-R2 (repair plan revision 4, validation/plans/v3-005d-v7r-calibration-repair.md) -------------
# A NEW row for the same estimator claim; V7-R above stays as recorded (archived suites keep every constant
# and code path above byte-identical). Everything below is additive. Nothing here is evidence by itself.
V7R2_ALPHA_TOST = 0.005  # per side, the two escaped_neutral bootstrap-t gates only (plan 3.4)
V7R2_SHARDS = 48
V7R2_SHARD_N = 2_250_000
V7R2_REPS_PER_SHARD = 225  # 225 replicates of 20 x 500 histories per shard
V7R2_REPLICATES = V7R2_SHARDS * V7R2_REPS_PER_SHARD  # 10 800
V7R2_PAIRS = V7R2_REPLICATES // 2  # 5400: replicate j of shard s with j of shard s + 24
V7R2_PAIR_SHIFT = V7R2_SHARDS // 2
V7R2_HELDOUT_SHARDS = tuple(range(36, 48))  # H = last quarter (2700 replicates); E = shards 0..35 (8100)
V7R2_REF_N = V7R_REF_N
V7R2_REF_INDEX = V7R2_SHARDS  # 48: distinct from every shard index 0..47 (V7R_REF_INDEX = 16 is shard 16 there)
V7R2_FALLBACK_SHARDS = 64  # 7200/0.005 fallback of the 3.8 selection rule: 64 x 2.25e6, reference k = 64
V7R2_FALLBACK_REF_INDEX = V7R2_FALLBACK_SHARDS
V7R2_SEED_BASE = 20531004  # lv5d base (integration route 4.4); r_index 16 as V7-R
V7R2_SEED_FIRST = V7R2_SEED_BASE + 1000 * V7R_R_INDEX  # 20547004: seed of shard k is this + k
V7R2_KAPPA_N = 1_000_000  # per kind and surrogate (full sample, no 3:1 split)
V7R2_KAPPA_SEED = 20980000  # data root; the bootstrap root is +1
V7R2_CAL_SEED = 20982000  # data root; the bootstrap root is +1
V7R2_A13_SEED = 20984000  # Amendment 13 Gaussian shard (sec_p, nuclear_local), if re-executed
V7R2_SEED_FAMILY = range(20980000, 20985000)  # guard: calibration / kappa roots (plan 6.7 item 5)
V7R2_CAL_PARTS = 39  # per surrogate: 78 parts of 83/82 repetitions (3200 = 2 * 83 + 37 * 82)
V7R2_KAPPA_KINDS = ("single", "pair")  # kind index 0 / 1 in the SeedSequence entropy
V7R2_KAPPA_PART_N = {"single": (1_000_000,), "pair": (500_000, 500_000)}
"""Rows of each kappa step (single replicates / pairs; a pair is two replicates, so a pair part draws 1e6 replicates). The split gives six steps of
about 645 s at the slow-host rates of plan 6.1 (1.29 ms per pair, 0.65 ms per single replicate); the plan states
"6 shards ~ 645 s" and the totals, and this split is the one consistent with both (flagged in the report)."""
CAL2_SCHEMA = "ionmc-v7r2-calibration-record-1"
V7R2_KAPPA_SCHEMA = "ionmc-v7r-kappa-record-1"
V7R2_KAPPA_COMBINED_SCHEMA = "ionmc-v7r-kappa-1"
V7R2_KAPPA_SE_MAX = 0.0005  # plan 6.7 item 4: coverage SE of the solved kappa, every margin
KAPPA_MARGINS = ("paired_low", "paired_high", "single_low", "single_high")
KAPPA_KIND_OF = {"paired_low": "pair", "paired_high": "pair", "single_low": "single", "single_high": "single"}
V7R_LEGACY_KAPPA_SEEDS = {"single": (20508000, 20505999), "pair": (20508001, 20505999)}
"""Legacy mode of the kappa step: (data seed, bootstrap seed base) of ``margin_kappas`` in
tests/ionmc/test_v7r_coverage.py (SEED_KAPPA, SEED_KAPPA + 1; SEED_BOOT0 - 1), n = 2e4 and the 3:1 split."""


def r2_shard_seed(k: int, shards: int = V7R2_SHARDS) -> int:
    """Seed of V7-R2 transport shard ``k`` (0 .. shards - 1) at the lv5d base: 20547004 + k."""
    if not 0 <= k < shards:
        raise ValueError(f"shard index outside 0..{shards - 1}")
    return V7R2_SEED_FIRST + k


def r2_ref_seed(shards: int = V7R2_SHARDS) -> int:
    """Seed of the reference of a ``shards``-shard layout: reference index = ``shards`` (48 -> 20547052,
    the 7200 fallback 64 -> 20547068); never equal to a shard seed."""
    return V7R2_SEED_FIRST + shards


def r2_seed_guard(root: int) -> int:
    """A calibration / kappa seed root must lie in 20980000-20984999 (plan 6.7 item 5); returns it."""
    if root not in V7R2_SEED_FAMILY:
        raise ValueError(f"seed root {root} outside {V7R2_SEED_FAMILY.start}-{V7R2_SEED_FAMILY.stop - 1}")
    return root


def _seed_word(root: int, *index: int) -> int:
    """One 32-bit word of ``SeedSequence([root, *index])`` (the integer seed of numpy streams below)."""
    return int(np.random.SeedSequence([r2_seed_guard(root), *index]).generate_state(1, dtype=np.uint32)[0])


def r2_kappa_part_seeds(surrogate_idx: int, kind: str, part: int) -> dict[str, int]:
    """Seeds of one kappa part: data = word of ``SeedSequence([20980000, surrogate_idx, kind_idx, part])``,
    bootstrap base = word of ``SeedSequence([20980001, ...])`` (kind_idx 0 single, 1 pair)."""
    ki = V7R2_KAPPA_KINDS.index(kind)
    return {"data": _seed_word(V7R2_KAPPA_SEED, surrogate_idx, ki, part),
            "boot": _seed_word(V7R2_KAPPA_SEED + 1, surrogate_idx, ki, part)}  # fmt: skip


def r2_cal_part_seeds(surrogate_idx: int, part: int) -> dict[str, int]:
    """Seeds of one V7-R2 calibration part: words of ``SeedSequence([20982000, surrogate_idx, part])`` (data) and
    ``SeedSequence([20982001, ...])`` (bootstrap base); replicate ``r`` of the part uses ``data + r`` and
    ``boot + r`` as in :func:`cal_part_seeds`."""
    return {"data": _seed_word(V7R2_CAL_SEED, surrogate_idx, part),
            "boot": _seed_word(V7R2_CAL_SEED + 1, surrogate_idx, part)}  # fmt: skip


def single_ingredients_all(
    sums: NDArray[np.float64], seed_base: int, reps_per_shard: int, b: int = V7R_B,
    alpha: float = V7R_ALPHA,
) -> dict[str, NDArray[np.float64]]:
    """Like :func:`single_ingredients` but over EVERY row of ``sums`` (no 3:1 evaluation split): the full-sample
    kappa solve of V7-R2 (plan 3.4). Row ``r`` uses ``single_stream(seed_base, r // reps_per_shard, r %
    reps_per_shard)``, so the first ``n_E`` rows equal ``single_ingredients`` row for row."""
    n = sums.shape[0]
    out = np.empty((n, 5))
    for r in range(n):
        rng = single_stream(seed_base, r // reps_per_shard, r % reps_per_shard)
        out[r] = boot_t_quantiles(sums[r], rng, b, alpha)
    return {k: out[:, i] for i, k in enumerate(("m_hat", "se", "t_lo", "t_hi", "n_degenerate"))}


def kappa_rows(
    sampler: Any, kind: str, n: int, seeds: dict[str, int], *, legacy: bool = False, b: int = V7R_B
) -> dict[str, NDArray[np.float64]]:
    """Per-row statistics ``a`` (studentised error of the truth), ``t_lo``, ``t_hi`` of ``n`` replicates
    (``kind="single"``: of the mean ``sampler.mu``; ``"pair"``: ``n`` pairs, i.e. ``2 n`` replicates, of 0).
    ``legacy=True`` is the code path of ``margin_kappas`` (single: the 3:1 evaluation split, rows 0..n_E-1,
    ``reps_per_shard = n``); V7-R2 mode uses every replicate. ``sampler(rng, rows)`` draws ``[rows, 20]`` block
    values and ``sampler.mu`` is the truth. Only the part's own data are read, so parts are independent."""
    rng = np.random.default_rng(seeds["data"])
    if kind == "single":
        x = sampler(rng, n)
        ing = (single_ingredients(x, seeds["boot"], reps_per_shard=n, b=b) if legacy
               else single_ingredients_all(x, seeds["boot"], n, b))  # fmt: skip
        a = (ing["m_hat"] - sampler.mu) / ing["se"]
    elif kind == "pair":
        x = sampler(rng, 2 * n)
        ing = pair_ingredients(x, seeds["boot"], reps_per_shard=n, b=b)
        a = ing["d_hat"] / ing["se_d"]
    else:
        raise ValueError(f"kind {kind!r} is not single/pair")
    return {"a": a, "t_lo": ing["t_lo"], "t_hi": ing["t_hi"]}


def kappa_coverage(rows: dict[str, NDArray[np.float64]], kappa: float) -> float:
    """Empirical coverage at interval scale ``kappa`` (the bisection function of ``margin_kappas``)."""
    a, lo, hi = rows["a"], rows["t_lo"], rows["t_hi"]
    ok = np.isfinite(a) & np.isfinite(lo) & np.isfinite(hi)
    return float(np.sum(ok & (kappa * hi >= a) & (kappa * lo <= a)) / a.size)


def kappa_solve(rows: dict[str, NDArray[np.float64]], target: float) -> float:
    """``kappa`` whose empirical coverage equals ``target``: the 60-step bisection on [0.05, 20] of
    ``margin_kappas`` (coverage is monotone in kappa), on the pooled rows."""
    lo_k, hi_k = 0.05, 20.0
    for _ in range(60):
        mid = 0.5 * (lo_k + hi_k)
        lo_k, hi_k = (mid, hi_k) if kappa_coverage(rows, mid) < target else (lo_k, mid)
    return 0.5 * (lo_k + hi_k)


def kappa_part_names(surrogate: str, kind: str, part: int) -> tuple[str, str]:
    """``(json, npz)`` partial file names of one kappa part (plan 6.5)."""
    stem = f"v7r-kappa-{surrogate}-{kind}-p{part}"
    return f"{stem}.json", f"{stem}.npz"


def kappa_combine(
    records: dict[tuple[str, str, int], dict[str, Any]],
    rows: dict[tuple[str, str, int], dict[str, NDArray[np.float64]]],
    *, part_n: dict[str, tuple[int, ...]] | None = None, legacy: bool = False,
    region: tuple[float, float] = V7R_REGION,
) -> dict[str, Any]:
    """Combine the kappa parts into the kappa record (``v7r-kappa.json``). ``records[(surrogate, kind, part)]`` is
    the part's sealed record (schema ``ionmc-v7r-kappa-record-1``: ``surrogate``, ``kind``, ``part``, ``n``,
    ``seeds``, ``arrays`` = {name: sha256 of the C-order float64 bytes}, ``fit``, ``mode``) and ``rows[...]`` its
    arrays. Fails closed (``SystemExit``) on a missing / extra part, a wrong ``n`` or mode, seeds that differ from
    the derived ones, an array whose sha256 differs from the record, or a fit that differs within a surrogate.
    The paired margins are solved on the pooled pair rows of a surrogate and the single margins on the pooled
    single rows; for each margin the coverage at the solved kappa and its binomial SE are reported (the
    requirement ``se <= 0.0005`` of the V7-R2 mode is reported per margin as ``se_ok``)."""
    part_n = part_n or ({"single": (V7R_LEGACY_KAPPA_N,), "pair": (V7R_LEGACY_KAPPA_N,)} if legacy
                        else V7R2_KAPPA_PART_N)  # fmt: skip
    mode = "legacy" if legacy else "v7r2"
    want = {(s, k, p) for s in CAL_SURROGATES for k in V7R2_KAPPA_KINDS for p in range(len(part_n[k]))}
    if set(records) != want or set(rows) != want:
        raise SystemExit(f"kappa records: missing {sorted(want - set(records))[:3]} / "
                         f"unexpected {sorted(map(str, set(records) - want))[:3]}")  # fmt: skip
    out: dict[str, Any] = {}
    names: list[dict[str, Any]] = []
    for si, s in enumerate(CAL_SURROGATES):
        fit0 = records[(s, "single", 0)].get("fit")
        pooled: dict[str, dict[str, NDArray[np.float64]]] = {}
        for kind in V7R2_KAPPA_KINDS:
            parts = []
            for p, n in enumerate(part_n[kind]):
                tag = f"{s}/{kind}/p{p}"
                r = records[(s, kind, p)]
                if r.get("schema") != V7R2_KAPPA_SCHEMA or r.get("mode") != mode:
                    raise SystemExit(f"kappa record {tag}: wrong schema or mode")
                if (r.get("surrogate"), r.get("kind"), r.get("part")) != (s, kind, p):
                    raise SystemExit(f"kappa record {tag}: identity differs from its key")
                if r.get("n") != n:
                    raise SystemExit(f"kappa record {tag}: n {r.get('n')!r} is not {n}")
                want_seeds = (r2_kappa_part_seeds(si, kind, p) if not legacy else
                              dict(zip(("data", "boot"), V7R_LEGACY_KAPPA_SEEDS[kind], strict=True)))  # fmt: skip
                if r.get("seeds") != want_seeds:
                    raise SystemExit(f"kappa record {tag}: seeds differ from the derived ones")
                if r.get("fit") != fit0:
                    raise SystemExit(f"kappa record {tag}: fit differs within the surrogate")
                arr = rows[(s, kind, p)]
                if set(arr) != {"a", "t_lo", "t_hi"} or set(r.get("arrays", {})) != set(arr):
                    raise SystemExit(f"kappa record {tag}: arrays are not exactly a / t_lo / t_hi")
                for name, v in arr.items():
                    digest = hashlib.sha256(np.ascontiguousarray(v, dtype="<f8").tobytes()).hexdigest()
                    if r["arrays"][name] != digest:
                        raise SystemExit(f"kappa record {tag}: array {name} does not match its sha256")
                names.append({"surrogate": s, "kind": kind, "part": p,
                              "content_sha256": r.get("content_sha256"), "arrays": r["arrays"]})  # fmt: skip
                parts.append(arr)
            pooled[kind] = {k: np.concatenate([q[k] for q in parts]) for k in ("a", "t_lo", "t_hi")}
        kappas: dict[str, float] = {}
        coverage: dict[str, Any] = {}
        for m, target in zip(KAPPA_MARGINS, (region[0], region[1], region[0], region[1]), strict=True):
            rr = pooled[KAPPA_KIND_OF[m]]
            kappas[m] = kappa_solve(rr, target)
            c = kappa_coverage(rr, kappas[m])
            n_rows = int(rr["a"].size)
            se = math.sqrt(c * (1.0 - c) / n_rows)
            coverage[m] = {"target": target, "coverage": c, "n": n_rows, "se": se,
                           "se_ok": bool(se <= V7R2_KAPPA_SE_MAX)}  # fmt: skip
        out[s] = {"kappas": kappas, "coverage": coverage, "fit": fit0}
    return {"schema": V7R2_KAPPA_COMBINED_SCHEMA, "mode": mode, "surrogates": out, "parts": names,
            "se_max": V7R2_KAPPA_SE_MAX,
            "se_ok": all(c["se_ok"] for v in out.values() for c in v["coverage"].values())}  # fmt: skip


V7R_LEGACY_KAPPA_N = 20_000


def r2_calibration_kappas(kappa_doc: dict[str, Any], surrogate: str) -> dict[str, float]:
    """The five margin kappas of one surrogate from a combined kappa record (``nominal`` = 1)."""
    if kappa_doc.get("schema") != V7R2_KAPPA_COMBINED_SCHEMA or kappa_doc.get("mode") != "v7r2":
        raise SystemExit("kappa record: not a V7-R2 combined kappa record")
    if not kappa_doc.get("se_ok"):
        raise SystemExit(f"kappa record: coverage SE exceeds {V7R2_KAPPA_SE_MAX}")
    k = kappa_doc["surrogates"][surrogate]["kappas"]
    return {"nominal": 1.0, **{m: float(k[m]) for m in KAPPA_MARGINS}}


# -- rep-based timing rule (route 4.1 (iii); plan 6.1, 6.7 item 7) ---------------------------------
def r2_part_layout() -> dict[str, Any]:
    """Layout every V7-R2 calibration part records (the 5400-pair design: 10 800 replicates of 225, 39 parts)."""
    return {"replicates": V7R2_REPLICATES, "reps_per_shard": V7R2_REPS_PER_SHARD, "parts": V7R2_CAL_PARTS}


def r2_calibration_layout(shards: int) -> dict[str, Any]:
    """Design a V7-R2 calibration status must be bound to for a ``shards``-shard transport layout: shards,
    replicates, pairs and calibration parts. Only the 48-shard design (5400 pairs, 39 parts per surrogate) has a
    calibration producer; for any other shard count ``parts`` is ``None``. A fabricated status could carry that
    same dictionary, so the rejection of unsupported shard counts rests on ``require_r2_calibration_producer``
    (``V7R2_CALIBRATED_SHARDS``), which the combine and both combine steps call before accepting any calibration
    layout (the combine inspects the status ``step`` and ``pass`` fields first, then runs the guard). TODO (repair plan 6.2, 7200/0.005 column: 102 parts of 63 repetitions, 14 400 replicates):
    the fallback calibration, to be implemented only if the selection rule of 3.8 picks 64 x 2.25e6."""
    reps = shards * V7R2_REPS_PER_SHARD
    return {"shards": shards, "replicates": reps, "pairs": reps // 2,
            "parts": V7R2_CAL_PARTS if shards == V7R2_SHARDS else None}  # fmt: skip


V7R2_CALIBRATED_SHARDS = (V7R2_SHARDS,)
"""Shard counts with a calibration producer. Only a test process may widen it (monkeypatch) to run the combine at CI
scale; production code never does."""


def require_r2_calibration_producer(shards: int) -> None:
    """Fail closed for a shard count without a calibration producer (everything except 48): the combine and the
    calibration-combine step call this FIRST, so that no status record, genuine or fabricated, can serve it. The
    7200-pair fallback (64 x 2.25e6) needs its own calibration (repair plan 6.2, 7200/0.005 column: 102 parts of 63
    repetitions, 14 400 replicates), which is not implemented."""
    if shards not in V7R2_CALIBRATED_SHARDS:
        raise SystemExit(
            f"v7r2: no calibration producer exists for {shards} shards (only {V7R2_CALIBRATED_SHARDS}); the fallback "
            "calibration (repair plan 6.2: 102 parts x 63 repetitions, 14 400 replicates) is not implemented "
            "(fail closed)"
        )


R2_STEP_TARGET_S = 1100.0  # planning target per step (27 % headroom under the 1500 s limit)
R2_STEP_LIMIT_S = 1500.0  # hard limit (run_suite.HOST_STEP_WALL_LIMIT_S)
R2_CAL_S_PER_REP = 12.2  # slow host, 10 800 replicates (plan 6.1; EXTRAPOLATED until one part is measured)
R2_CAL_FIXED_S = 70.0  # per-part fixed cost (surrogate fit or pool load, record load)
R2_KAPPA_S_PER_ROW = {"pair": 1.29e-3, "single": 0.65e-3}  # slow host, per pair / per single replicate


def rep_plan(name: str, items: int, s_per_item: float, fixed_s: float, declared_parts: int, *,
             target_s: float = R2_STEP_TARGET_S, limit_s: float = R2_STEP_LIMIT_S) -> dict[str, Any]:
    """Shard-table row of a repetition-based step family (Monte Carlo repetitions or kappa rows; the history-based
    ``steps_v5d.shard_plan_c`` does not cover them). ``items`` are split over ``declared_parts`` steps by
    ``cal_part_sizes`` (the largest part has ``ceil(items / parts)``); a part takes ``largest * s_per_item +
    fixed_s``. ``parts_needed`` is the smallest part count whose largest part fits the target, ``ceil(items /
    floor((target_s - fixed_s) / s_per_item))`` (3200 repetitions at 12.2 s: at most 84 per part, so 39 parts of 83 /
    82). ``ok`` iff ``declared_parts >= parts_needed`` and the planned time is <= ``target_s`` (< ``limit_s``)."""
    if target_s > limit_s or target_s <= fixed_s or s_per_item <= 0.0 or declared_parts < 1:
        raise ValueError("rep_plan: need fixed_s < target_s <= limit_s, s_per_item > 0 and parts >= 1")
    per_part_max = math.floor((target_s - fixed_s) / s_per_item)
    if per_part_max < 1:
        raise ValueError("rep_plan: not even one item fits the target")
    needed = max(1, math.ceil(items / per_part_max))
    largest = max(cal_part_sizes(items, declared_parts))
    planned = largest * s_per_item + fixed_s
    return {"name": name, "items": items, "s_per_item": s_per_item, "fixed_s": fixed_s,
            "declared_parts": declared_parts, "parts_needed": needed, "largest_part_items": largest,
            "planned_s_per_part": planned, "target_s": target_s, "limit_s": limit_s,
            "headroom_to_limit": 1.0 - planned / limit_s,
            "ok": bool(declared_parts >= needed and planned <= target_s)}  # fmt: skip


def r2_rep_plans(*, cal_s_per_rep: float = R2_CAL_S_PER_REP, fixed_s: float = R2_CAL_FIXED_S) -> dict[str, Any]:
    """The V7-R2 rep-based shard table: the 78 calibration parts and the six kappa steps. ``cal_s_per_rep`` is the
    slow-host rate (replace the extrapolated default by the measurement of one short calibration part, plan 6.1)."""
    plans = {"v7r2 calibration (3200 reps x 39 parts per surrogate)":
             rep_plan("v7r2-calibration", V7R_CAL_REPS, cal_s_per_rep, fixed_s, V7R2_CAL_PARTS)}  # fmt: skip
    for kind in V7R2_KAPPA_KINDS:
        for p, n in enumerate(V7R2_KAPPA_PART_N[kind]):
            plans[f"v7r2 kappa {kind} p{p}"] = rep_plan(  # n pairs (kind pair) / single replicates
                f"v7r2-kappa-{kind}-p{p}", n, R2_KAPPA_S_PER_ROW[kind], fixed_s, 1
            )
    return {"plans": plans, "ok": all(v["ok"] for v in plans.values())}
