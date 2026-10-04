"""Depth-dose metrics for reference-engine comparisons (numpy only).

Depths are in mm and refer to bin centres: ``depth[i]`` is the centre of bin ``i``. Distal
levels are found by linear interpolation between the two samples that bracket the first
downward crossing after the peak bin.
"""

from __future__ import annotations

import numpy as np


def integral_depth_dose(dose3d: np.ndarray, axis: int, spacing: float | None = None) -> np.ndarray:
    """Sum a 3D dose array over the two axes other than ``axis``.

    ``spacing`` is the lateral voxel area factor (mm^2 per column cell); it is unused unless
    given, in which case the sum is multiplied by it (integral over area instead of a plain sum).
    """
    arr = np.asarray(dose3d, dtype=np.float64)
    if arr.ndim != 3 or not 0 <= axis < 3:
        raise ValueError("dose3d must be 3D and axis in 0..2")
    other = tuple(a for a in range(3) if a != axis)
    out = arr.sum(axis=other)
    return out if spacing is None else out * float(spacing)


def normalize_to_peak(curve: np.ndarray) -> np.ndarray:
    """Return ``curve / max(curve)``; fails on a non-positive maximum."""
    c = np.asarray(curve, dtype=np.float64)
    peak = float(np.max(c))
    if not peak > 0.0:
        raise ValueError("curve maximum must be positive")
    return c / peak


def _check(depth: np.ndarray, curve: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    d = np.asarray(depth, dtype=np.float64)
    c = np.asarray(curve, dtype=np.float64)
    if d.ndim != 1 or d.shape != c.shape or d.size < 3:
        raise ValueError("depth and curve must be 1D arrays of equal length >= 3")
    if np.any(np.diff(d) <= 0):
        raise ValueError("depth must be strictly increasing")
    if not np.all(np.isfinite(c)):
        raise ValueError("curve must be finite")
    return d, c


def peak_depth(depth: np.ndarray, curve: np.ndarray) -> float:
    """Depth (mm) of the maximum bin (resolution is the bin width)."""
    d, c = _check(depth, curve)
    return float(d[int(np.argmax(c))])


def distal_depth(depth: np.ndarray, curve: np.ndarray, fraction: float) -> float:
    """Depth (mm) where the curve falls to ``fraction`` of its maximum, distal of the peak."""
    d, c = _check(depth, curve)
    if not 0.0 < fraction < 1.0:
        raise ValueError("fraction must be in (0, 1)")
    ip = int(np.argmax(c))
    level = fraction * c[ip]
    for i in range(ip, d.size - 1):
        if c[i] >= level > c[i + 1]:
            t = (c[i] - level) / (c[i] - c[i + 1])
            return float(d[i] + t * (d[i + 1] - d[i]))
    raise ValueError(f"curve does not fall below {fraction:g} of the peak distally")


def r80(depth: np.ndarray, curve: np.ndarray) -> float:
    """Distal depth (mm) of 80 % of the peak."""
    return distal_depth(depth, curve, 0.8)


def r90(depth: np.ndarray, curve: np.ndarray) -> float:
    """Distal depth (mm) of 90 % of the peak."""
    return distal_depth(depth, curve, 0.9)


def distal_falloff_80_20(depth: np.ndarray, curve: np.ndarray) -> float:
    """Distal fall-off width (mm): depth(20 %) minus depth(80 %)."""
    return distal_depth(depth, curve, 0.2) - distal_depth(depth, curve, 0.8)


# --- lateral profile analysis ------------------------------------------------------------------
#
# Estimator (all lengths in mm). For a 1-D lateral profile P_i >= 0 (dose summed over the
# orthogonal lateral axis) on bins of width h with centres x_i measured from the grid centre
# (the beam axis):
#   1. window: keep |x_i| <= W (default: the whole grid, W = n h / 2);
#   2. second central moment of the dose profile m2 = sum P (x - mu)^2 / sum P,
#      mu = sum P x / sum P;
#   3. Sheppard's correction for binning: sigma^2 = m2 - h^2 / 12.
# There is NO background/pedestal subtraction: Monte Carlo dose has no additive background, the
# tails (nuclear halo, large-angle scatter) are physical, and the weights P are never negative
# (negative dose is rejected). The slab estimate is the mean of the x- and y-projection results.
#
# Window dependence. The estimator is the second moment of the dose profile inside the window:
# the full-field value (W = 60 mm) includes all tails by definition and is NOT the Gaussian core
# width; the W = 20 mm value excludes dose beyond +-20 mm. Both are reported; neither is the
# "core sigma". Bias bounds. (a) Binning: for bin-integrated dose Sheppard's correction is exact
# for a Gaussian (binned second moment = sigma^2 + h^2/12). (b) Truncation of a Gaussian core at
# |x| = W: m2 = sigma^2 [1 - 2 a phi(a) / (2 Phi(a) - 1)], a = W / sigma, with the dose fraction
# beyond +-W equal to erfc(a / sqrt 2). For W = 20 mm: sigma = 2 / 3 / 4 / 5 / 6 mm gives a = 10 /
# 6.67 / 5 / 4 / 3.33, fraction beyond 1.5e-23 / 2.6e-11 / 5.7e-7 / 6.3e-5 / 8.6e-4 and relative
# bias of sigma^2 of -1.5e-21 / -1.2e-9 / -1.5e-5 / -1.1e-3 / -1.0e-2. (c) Tails add
# (mass_tail / mass_total) * <x^2>_tail to m2 by construction (tested with a power-law tail
# against the analytic value). (d) Noise: the variance of m2 grows with the weight x^2 of the wing
# bins, so the statistical uncertainty is taken from the spread across independent seeds, never
# from this routine; for counting statistics the estimator is unbiased (tested with Poisson/
# multinomial histograms).


def bin_centres(n: int, h: float) -> np.ndarray:
    """Bin-centre coordinates (mm) of ``n`` bins of width ``h`` centred on the grid centre."""
    return (np.arange(n) + 0.5) * h - 0.5 * n * h


def lateral_variance_1d(
    x_mm: np.ndarray,
    profile: np.ndarray,
    bin_mm: float,
    *,
    window_mm: float | None = None,
) -> float:
    """Sheppard-corrected second moment (mm^2) of a non-negative 1-D lateral dose profile.

    No pedestal subtraction; ``window_mm`` keeps ``|x| <= window_mm`` (default whole grid).
    See the module notes above for the estimator, window dependence and bias bounds.
    """
    x = np.asarray(x_mm, dtype=np.float64)
    p = np.asarray(profile, dtype=np.float64)
    if x.ndim != 1 or x.shape != p.shape or x.size < 8:
        raise ValueError("x_mm and profile must be 1D arrays of equal length >= 8")
    if not np.all(np.isfinite(p)) or not bin_mm > 0:
        raise ValueError("invalid profile or bin width")
    if np.any(p < 0.0):
        raise ValueError("profile has negative dose")
    keep = np.ones(x.shape, dtype=bool) if window_mm is None else np.abs(x) <= float(window_mm)
    xs, ps = x[keep], p[keep]
    if xs.size < 8:
        raise ValueError("window too small")
    m0 = float(ps.sum())
    if not m0 > 0.0:
        raise ValueError("profile has non-positive integral in the window")
    mu = float((ps * xs).sum() / m0)
    m2 = float((ps * (xs - mu) ** 2).sum() / m0)
    return m2 - bin_mm**2 / 12.0


def lateral_sigma2_slab(
    slab: np.ndarray,
    bin_mm: tuple[float, float],
    *,
    window_mm: float | None = None,
) -> dict[str, float]:
    """Lateral variance (mm^2) of a 2-D slab ``slab[ix, iy]``: mean of the x and y projections."""
    a = np.asarray(slab, dtype=np.float64)
    if a.ndim != 2:
        raise ValueError("slab must be 2D")
    vx = lateral_variance_1d(
        bin_centres(a.shape[0], bin_mm[0]),
        a.sum(axis=1),
        bin_mm[0],
        window_mm=window_mm,
    )
    vy = lateral_variance_1d(
        bin_centres(a.shape[1], bin_mm[1]),
        a.sum(axis=0),
        bin_mm[1],
        window_mm=window_mm,
    )
    return {"sigma2_mm2": 0.5 * (vx + vy), "sigma2_x_mm2": vx, "sigma2_y_mm2": vy}


def lateral_sigma2_at_depth_fractions(
    dose3d: np.ndarray,
    bin_mm: tuple[float, float, float],
    r80_mm: float,
    fractions: tuple[float, ...] = (0.5, 0.9),
    *,
    window_mm: float | None = None,
) -> dict[float, dict[str, float]]:
    """Slab lateral variance at ``z = f * R80`` (the bin containing that depth, one bin thick).

    ``dose3d`` is ``(nx, ny, nz)``; ``r80_mm`` comes from the IDD of the same run. The slab
    centre depth and its ratio to R80 are reported with the variance.
    """
    d = np.asarray(dose3d, dtype=np.float64)
    if d.ndim != 3 or not r80_mm > 0:
        raise ValueError("dose3d must be 3D and r80_mm positive")
    out: dict[float, dict[str, float]] = {}
    for f in fractions:
        iz = int(np.floor(f * r80_mm / bin_mm[2]))
        if not 0 <= iz < d.shape[2]:
            raise ValueError("slab outside the dose grid")
        res = lateral_sigma2_slab(d[:, :, iz], (bin_mm[0], bin_mm[1]), window_mm=window_mm)
        centre = (iz + 0.5) * bin_mm[2]
        out[f] = {**res, "slab_centre_mm": centre, "slab_centre_over_r80": centre / r80_mm}
    return out
