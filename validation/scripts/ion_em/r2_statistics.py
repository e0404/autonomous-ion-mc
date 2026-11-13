"""V3-008 revision 2: estimators and TOST statistics of the TOPAS-agreement rows (numpy + stdlib).

Implements ``validation/plans/v3-008/r2-topas-target.md`` section 4 (observables and estimators,
identical algorithm for ionmc and TOPAS on identical bin edges) and section 5 (criteria; all TOST:
"PASS iff |x - r| + 1.645 sigma_comb <= tol; inconclusive = FAIL"). Input arrays are already on the
common bin edges (the loaders in ``r2_io`` verify the edges); nothing here reads a file.

Statistics (section 5, verbatim): "sigma_comb = sqrt(sigma_T^2 + sigma_i^2). sigma_T: sample SD over
the K = 4 TOPAS parts / sqrt(4), multiplied by t_{0.95,3}/1.645 = 1.430 (small-sample inflation).
sigma_i: ionmc batch SE (40 batches)."

Implementation decisions (the plan is silent on these details; the statistically standard option is
chosen and recorded here and in every report):

* D1 (point values): the graded value of each engine is the estimator applied to the POOLED data
  (ionmc: the mean of the equal-size batches; TOPAS: the run totals of the K parts summed and
  divided by the summed histories). The same estimator on the pooled data is the least biased
  choice for the non-linear estimators (parabola peak, interpolated R80, Gaussian ML fit).
* D2 (uncertainties, batch means): sigma_i = SD(ddof=1) over the per-batch values of the statistic /
  sqrt(B); sigma_T = 1.430 * SD(ddof=1) over the per-part values / sqrt(K). For a statistic that
  combines both engines (a ratio r, the mean ratio <r>), the per-batch (per-part) value is the
  statistic with that one batch (part) substituted for the pooled ionmc (TOPAS) data and the other
  engine held at its pooled value (first-order batch-means / delta method; correlations between
  depths of the same histories are kept because the substitution is made per batch / part).
  K must be 4 (the 1.430 factor is stated for K = 4 only; any other K fails closed).
* D3 (TOST form of the relative tolerances): every TD criterion is graded literally as
  ``|x - r| + 1.645 sigma_comb <= tol`` with x = ionmc, r = TOPAS (pooled) and the tolerance taken
  from the TOPAS value ("10 % of W_T" -> 0.10 W_T; "PE ratio 3 %" -> 0.03 PE_T; "plateau groups,
  E_i/E_T per group, 2 %" -> 0.02 E_T,group).
* D4 (bin selection): the entrance window, the plateau groups and the +-20 mm moment window are
  sets of whole bins whose centres lie inside the stated depth/lateral interval (all intervals are
  bin-aligned on the r2 grids); the +-2 sigma fit window is the set of bins whose centres satisfy
  ``|x| <= 2 sigma``, and the truncated likelihood is normalised to exactly the span of those bins.
* D5 (range class): the "R >= 30 mm" / "30-120 mm" / ">= 120 mm" classes of section 4.1 use the
  case's planning range R_plan (section 3, the same rule that fixed h_z); for all twelve cases
  R_plan and any plausible R80 fall in the same class (closest: 22.5 mm vs 30 mm).
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

Z95 = 1.645  # one-sided 95 % normal quantile of the TOST rule (section 5)
T_INFLATION_K4 = 1.430  # t_{0.95,3}/1.645 as stated in section 5 (K = 4 TOPAS parts)
K_TOPAS = 4
PEAK_GROUP_BINS = 5  # w_p = 5 h_z (section 4.1)
FIT_MAX_ITER = 50  # section 4.2
FIT_REL_TOL = 1e-4  # section 4.2
MOMENT_HALF_WINDOW_MM = 20.0  # section 4.2 "Sheppard-corrected second moment in +-20 mm"
FIT_WINDOW_SIGMAS = 2.0  # section 4.2 "truncated to |x| <= 2 sigma"
HALO_SIGMAS = 3.0  # section 4.2 "dose outside +-3 sigma_core over the full field"
_EDGE_EPS = 1e-9
FIT_CYCLE_MAX_AMPLITUDE = 0.02  # section 4.2, amendment R2-STOP3


class StatisticsError(Exception):
    """An estimator is undefined on the given data (fails closed: never a pass)."""


# ---------------------------------------------------------------------------------------------
# uncertainty and TOST
# ---------------------------------------------------------------------------------------------


def batch_se(values: Sequence[float] | np.ndarray) -> float:
    """ionmc batch standard error: SD(ddof=1) / sqrt(B) (B >= 2, all finite)."""
    v = np.asarray(values, dtype=np.float64)
    if v.ndim != 1 or v.size < 2 or not np.all(np.isfinite(v)):
        raise StatisticsError(f"batch values must be >= 2 finite numbers, got {v.size}")
    return float(np.std(v, ddof=1) / math.sqrt(v.size))


def topas_se(values: Sequence[float] | np.ndarray) -> float:
    """TOPAS part standard error: 1.430 * SD(ddof=1) / sqrt(4) (section 5; K = 4 only)."""
    v = np.asarray(values, dtype=np.float64)
    if v.ndim != 1 or v.size != K_TOPAS or not np.all(np.isfinite(v)):
        raise StatisticsError(
            f"sigma_T needs exactly K = {K_TOPAS} finite part values (the 1.430 inflation is "
            f"stated for K = 4 only), got {v.size}"
        )
    return float(T_INFLATION_K4 * np.std(v, ddof=1) / math.sqrt(v.size))


def tost(x: float, r: float, sigma_i: float, sigma_t: float, tol: float) -> dict[str, Any]:
    """Section 5: PASS iff ``|x - r| + 1.645 sigma_comb <= tol``, sigma_comb = sqrt(sT^2 + si^2).
    A non-finite input is a FAIL (invalid); ``|x - r| <= tol`` without the margin is reported as
    ``inconclusive`` and is a FAIL ("inconclusive = FAIL")."""
    vals = (x, r, sigma_i, sigma_t, tol)
    out: dict[str, Any] = {
        "x": x,
        "r": r,
        "sigma_i": sigma_i,
        "sigma_t": sigma_t,
        "tol": tol,
    }
    if not all(isinstance(v, int | float) and math.isfinite(v) for v in vals):
        return out | {"sigma_comb": None, "statistic": None, "verdict": "FAIL", "note": "invalid"}
    s = math.sqrt(sigma_t**2 + sigma_i**2)
    stat = abs(x - r) + Z95 * s
    passed = stat <= tol
    note = "pass" if passed else ("inconclusive" if abs(x - r) <= tol else "outside tolerance")
    return out | {
        "sigma_comb": s,
        "statistic": stat,
        "verdict": "PASS" if passed else "FAIL",
        "note": note,
    }


def band(r: float, sigma: float, lower: float, upper: float) -> dict[str, Any]:
    """TL (a): PASS iff ``(r - 1) - 1.645 sigma >= lower`` and ``(r - 1) + 1.645 sigma <= upper``
    (section 5.3: lower -5 %, upper +40 %). Non-finite input: FAIL (invalid)."""
    out: dict[str, Any] = {"r": r, "sigma": sigma, "lower": lower, "upper": upper}
    if not all(isinstance(v, int | float) and math.isfinite(v) for v in (r, sigma)):
        return out | {"verdict": "FAIL", "note": "invalid"}
    lo, hi = (r - 1.0) - Z95 * sigma, (r - 1.0) + Z95 * sigma
    passed = lo >= lower and hi <= upper
    inside = lower <= r - 1.0 <= upper
    note = "pass" if passed else ("inconclusive" if inside else "outside band")
    return out | {
        "low_edge": lo,
        "high_edge": hi,
        "verdict": "PASS" if passed else "FAIL",
        "note": note,
    }


def sampled(
    fn: Callable[[Any, Any], float],
    ion_pooled: Any,
    ion_batches: Sequence[Any],
    top_pooled: Any,
    top_parts: Sequence[Any],
) -> tuple[float, float, float]:
    """Decision D2 for a two-engine statistic ``fn(ionmc, topas)``: returns (value, sigma_i,
    sigma_T) by substituting one batch (part) at a time."""
    value = fn(ion_pooled, top_pooled)
    si = batch_se([fn(b, top_pooled) for b in ion_batches])
    st = topas_se([fn(ion_pooled, p) for p in top_parts])
    return value, si, st


# ---------------------------------------------------------------------------------------------
# 4.1 depth dose
# ---------------------------------------------------------------------------------------------


def range_class(r_plan_mm: float) -> dict[str, float]:
    """Section 4.1 classes (decision D5): entrance window and plateau group size g."""
    if r_plan_mm >= 120.0:
        return {"ent_lo": 1.0, "ent_hi": 3.0, "group_mm": 5.0}
    if r_plan_mm >= 30.0:
        return {"ent_lo": 1.0, "ent_hi": 3.0, "group_mm": 2.0}
    return {"ent_lo": 0.5, "ent_hi": 1.5, "group_mm": 1.0}


def _centres(n: int, hz: float) -> np.ndarray:
    return (np.arange(n, dtype=np.float64) + 0.5) * hz


def peak(e: np.ndarray, hz: float) -> tuple[float, float]:
    """Section 4.1 "Peak: E rebinned to w_p = 5 h_z; D_peak and z_peak from a parabola through
    the maximal group and its two neighbours." Groups start at z = 0 (a trailing partial group is
    dropped); a maximal group at either end or a non-concave triple is undefined (fails closed).
    Returns (z_peak, D_peak)."""
    e = np.asarray(e, dtype=np.float64)
    n = e.size // PEAK_GROUP_BINS
    if n < 3 or not np.all(np.isfinite(e)):
        raise StatisticsError("depth curve too short or non-finite for the peak estimator")
    g = e[: n * PEAK_GROUP_BINS].reshape(n, PEAK_GROUP_BINS).mean(axis=1)
    m = int(np.argmax(g))
    if m == 0 or m == n - 1:
        raise StatisticsError("maximal peak group at the end of the grid")
    y0, y1, y2 = g[m - 1], g[m], g[m + 1]
    den = y0 - 2.0 * y1 + y2
    if not den < 0.0:
        raise StatisticsError("peak triple is not concave")
    w = PEAK_GROUP_BINS * hz
    zc = (m + 0.5) * w
    z_peak = zc + 0.5 * w * (y0 - y2) / den
    d_peak = y1 - (y0 - y2) ** 2 / (8.0 * den)
    return float(z_peak), float(d_peak)


def distal_level(e: np.ndarray, hz: float, z_peak: float, level: float) -> float:
    """Section 4.1 "the first depth distal of z_peak at which E falls to <level> of D_peak, by
    linear interpolation on the h_z grid": the first bin centre beyond z_peak with E <= level,
    interpolated linearly against the preceding bin centre (which must lie above the level)."""
    e = np.asarray(e, dtype=np.float64)
    zc = _centres(e.size, hz)
    idx = np.nonzero((zc > z_peak) & (e <= level))[0]
    if idx.size == 0:
        raise StatisticsError("the curve does not fall to the level distal of the peak")
    j = int(idx[0])
    if j == 0 or not e[j - 1] > level:
        raise StatisticsError("no bracketing bin above the level before the crossing")
    return float(zc[j - 1] + (e[j - 1] - level) / (e[j - 1] - e[j]) * hz)


def _bins_in(n: int, hz: float, lo: float, hi: float) -> np.ndarray:
    zc = _centres(n, hz)
    return np.nonzero((zc > lo - _EDGE_EPS) & (zc < hi + _EDGE_EPS))[0]


def entrance_dose(e: np.ndarray, hz: float, r_plan_mm: float) -> float:
    """Section 4.1 "D_ent = mean of E over [1, 3] mm (R >= 30 mm) or [0.5, 1.5] mm (R < 30 mm)"."""
    c = range_class(r_plan_mm)
    idx = _bins_in(len(e), hz, c["ent_lo"], c["ent_hi"])
    if idx.size == 0:
        raise StatisticsError("empty entrance window")
    return float(np.mean(np.asarray(e)[idx]))


def plateau_groups(e: np.ndarray, hz: float, r_plan_mm: float, r80_ref: float) -> np.ndarray:
    """Section 4.1 "E averaged over groups of g = 5 mm (R >= 120 mm), 2 mm (30-120 mm), 1 mm
    (< 30 mm) from the entrance-window start to 0.8 R80_TOPAS (last partial group dropped)".
    ``r80_ref`` is the pooled TOPAS R80 of the row's configuration."""
    c = range_class(r_plan_mm)
    start, g, end = c["ent_lo"], c["group_mm"], 0.8 * r80_ref
    out = []
    j = 0
    while start + (j + 1) * g <= end + _EDGE_EPS:
        idx = _bins_in(len(e), hz, start + j * g, start + (j + 1) * g)
        if idx.size == 0:
            raise StatisticsError("empty plateau group")
        out.append(float(np.mean(np.asarray(e)[idx])))
        j += 1
    if not out:
        raise StatisticsError("no complete plateau group before 0.8 R80_T")
    return np.asarray(out)


def depth_quantities(e: np.ndarray, hz: float, r_plan_mm: float) -> dict[str, float]:
    """All single-curve quantities of section 4.1: z_peak, D_peak, R80, R20, W = R20 - R80, R90
    (report-only), D_ent and PE = D_peak / D_ent."""
    z_peak, d_peak = peak(e, hz)
    r80 = distal_level(e, hz, z_peak, 0.8 * d_peak)
    r20 = distal_level(e, hz, z_peak, 0.2 * d_peak)
    r90 = distal_level(e, hz, z_peak, 0.9 * d_peak)
    d_ent = entrance_dose(e, hz, r_plan_mm)
    if not d_ent > 0.0:
        raise StatisticsError("non-positive entrance dose")
    return {
        "z_peak": z_peak,
        "d_peak": d_peak,
        "r80": r80,
        "r20": r20,
        "w": r20 - r80,
        "r90": r90,
        "d_ent": d_ent,
        "pe": d_peak / d_ent,
    }


# ---------------------------------------------------------------------------------------------
# 4.2 lateral profiles
# ---------------------------------------------------------------------------------------------

_erfc = np.vectorize(math.erfc, otypes=[np.float64])
_SQRT2 = math.sqrt(2.0)


def _norm_interval(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Phi(hi) - Phi(lo) evaluated with erfc on the tail side (no cancellation)."""
    lo, hi = np.broadcast_arrays(np.asarray(lo, np.float64), np.asarray(hi, np.float64))
    out = np.empty(lo.shape)
    right = lo >= 0.0
    left = hi <= 0.0
    mid = ~(right | left)
    out[right] = 0.5 * (_erfc(lo[right] / _SQRT2) - _erfc(hi[right] / _SQRT2))
    out[left] = 0.5 * (_erfc(-hi[left] / _SQRT2) - _erfc(-lo[left] / _SQRT2))
    out[mid] = 1.0 - 0.5 * _erfc(hi[mid] / _SQRT2) - 0.5 * _erfc(-lo[mid] / _SQRT2)
    return out


def second_moment_sigma(w: np.ndarray, edges: np.ndarray) -> float:
    """Sheppard-corrected centred second moment over the bins inside +-20 mm (section 4.2):
    sqrt(sum w x^2 / sum w - h^2/12) with x the bin centre and h the bin width."""
    w = np.asarray(w, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.float64)
    c = 0.5 * (edges[:-1] + edges[1:])
    h = float(np.mean(np.diff(edges)))
    sel = (edges[:-1] >= -MOMENT_HALF_WINDOW_MM - _EDGE_EPS) & (
        edges[1:] <= MOMENT_HALF_WINDOW_MM + _EDGE_EPS
    )
    tot = float(np.sum(w[sel]))
    if not tot > 0.0:
        raise StatisticsError("no signal inside +-20 mm")
    m2 = float(np.sum(w[sel] * c[sel] ** 2)) / tot - h * h / 12.0
    if not m2 > 0.0:
        raise StatisticsError("non-positive Sheppard-corrected second moment")
    return math.sqrt(m2)


def _window_ml(w: np.ndarray, lo: np.ndarray, hi: np.ndarray, s0: float) -> float:
    """Binned ML scale of a centred normal truncated to [lo[0], hi[-1]] (golden section on
    ln sigma in [s0/4, 4 s0]); a maximum at the bracket edge is a non-convergence."""
    a, b = float(lo[0]), float(hi[-1])
    wt = float(np.sum(w))

    def loglik(log_s: float) -> float:
        s = math.exp(log_s)
        p = _norm_interval(lo / s, hi / s)
        z = float(_norm_interval(np.array([a / s]), np.array([b / s]))[0])
        if z <= 0.0:
            return -math.inf
        pos = w > 0.0
        if np.any(p[pos] <= 0.0):
            return -math.inf
        return float(np.sum(w[pos] * np.log(p[pos]))) - wt * math.log(z)

    x0, x1 = math.log(s0 / 4.0), math.log(s0 * 4.0)
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = x1 - gr * (x1 - x0), x0 + gr * (x1 - x0)
    fc, fd = loglik(c), loglik(d)
    for _ in range(200):
        if x1 - x0 < 1e-12:
            break
        if fc >= fd:
            x1, d, fd = d, c, fc
            c = x1 - gr * (x1 - x0)
            fc = loglik(c)
        else:
            x0, c, fc = c, d, fd
            d = x0 + gr * (x1 - x0)
            fd = loglik(d)
    xm = 0.5 * (x0 + x1)
    lo_b, hi_b = math.log(s0 / 4.0), math.log(s0 * 4.0)
    if xm - lo_b < 1e-6 or hi_b - xm < 1e-6:
        raise StatisticsError("ML scale at the edge of its bracket")
    return math.exp(xm)


@dataclass(frozen=True)
class CoreFit:
    sigma: float
    sigma0: float
    iterations: int
    cycle_period: int = 0  # > 0: sigma is the mean of a window cycle (amendment R2-STOP3)
    cycle_amplitude: float = 0.0  # (max - min) / mean of the cycle's ML scales


def core_sigma(w: np.ndarray, edges: np.ndarray) -> CoreFit:
    """Section 4.2 "binned maximum-likelihood fit of a centred Gaussian, integrated over each bin,
    to the profile truncated to |x| <= 2 sigma (window-truncated normal likelihood), iterated from
    the Sheppard-corrected second moment in +-20 mm until sigma changes by < 1e-4 relative
    (<= 50 iterations; non-convergence = invalid point)". Raises :class:`StatisticsError` for an
    invalid point.

    Section 4.2, amendment R2-STOP3 (window cycles). The window {bins with |centre| <= 2 s} is a
    step function of s, so the iteration can have no self-consistent window and 2-cycle (or
    longer-cycle) between windows. The convergence test above is applied first and unchanged, so
    every fit that converged before is bit-identical. If it fails and the window just used is the
    window of an earlier iteration j (not the immediately preceding one: that is a fixed window,
    handled by the ordinary iteration), the ML scales sigma_j ... sigma_{k-1} of the windows of the
    cycle are collected; sigma_core is their arithmetic mean, accepted iff the relative amplitude
    A = (max - min) / mean <= FIT_CYCLE_MAX_AMPLITUDE (0.02), else the point is invalid. The cycle
    length and A are recorded on the returned :class:`CoreFit`. Every other non-termination within
    FIT_MAX_ITER remains invalid."""
    w = np.asarray(w, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.float64)
    if w.ndim != 1 or edges.size != w.size + 1 or not np.all(np.isfinite(w)) or np.any(w < 0):
        raise StatisticsError("profile must be finite, non-negative and match its edges")
    c = 0.5 * (edges[:-1] + edges[1:])
    s0 = second_moment_sigma(w, edges)
    s = s0
    keys: list[tuple[int, int]] = []  # window (first, last bin) of each earlier iteration
    scales: list[float] = []  # ML scale produced on that window
    for it in range(1, FIT_MAX_ITER + 1):
        sel = np.nonzero(np.abs(c) <= FIT_WINDOW_SIGMAS * s + _EDGE_EPS)[0]
        if sel.size < 3 or not np.sum(w[sel]) > 0.0 or np.any(np.diff(sel) != 1):
            raise StatisticsError("fewer than 3 bins (or no signal) inside +-2 sigma")
        key = (int(sel[0]), int(sel[-1]))
        s_new = _window_ml(w[sel], edges[sel], edges[sel + 1], s)
        if abs(s_new - s) < FIT_REL_TOL * s:
            return CoreFit(s_new, s0, it)
        if key in keys[:-1]:
            cyc = scales[keys.index(key) :]
            mean = sum(cyc) / len(cyc)
            amp = (max(cyc) - min(cyc)) / mean
            if amp > FIT_CYCLE_MAX_AMPLITUDE:
                raise StatisticsError(f"window cycle amplitude {amp:.4g} > 0.02")
            return CoreFit(mean, s0, it, cycle_period=len(cyc), cycle_amplitude=amp)
        keys.append(key)
        scales.append(s_new)
        s = s_new
    raise StatisticsError(f"core fit did not converge in {FIT_MAX_ITER} iterations")


def halo_fraction(w: np.ndarray, edges: np.ndarray, sigma_core: float) -> float:
    """Report-only (section 4.2): fraction of the profile with |x_centre| > 3 sigma_core over the
    full field."""
    w = np.asarray(w, dtype=np.float64)
    c = 0.5 * (np.asarray(edges)[:-1] + np.asarray(edges)[1:])
    tot = float(np.sum(w))
    if not tot > 0.0:
        raise StatisticsError("empty profile")
    return float(np.sum(w[np.abs(c) > HALO_SIGMAS * sigma_core])) / tot


def core_xy_detail(
    px: np.ndarray, py: np.ndarray, ex: np.ndarray, ey: np.ndarray
) -> tuple[float, int]:
    """sigma_core = mean of the x and y fits (section 4.2) and the number (0-2) of those two fits
    that used the window-cycle rule (amendment R2-STOP3)."""
    fx, fy = core_sigma(px, ex), core_sigma(py, ey)
    return 0.5 * (fx.sigma + fy.sigma), int(fx.cycle_period > 0) + int(fy.cycle_period > 0)


def core_xy(px: np.ndarray, py: np.ndarray, ex: np.ndarray, ey: np.ndarray) -> float:
    """sigma_core = mean of the x and y fits (section 4.2)."""
    return core_xy_detail(px, py, ex, ey)[0]
