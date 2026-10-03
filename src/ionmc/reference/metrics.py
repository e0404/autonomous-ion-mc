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
