"""Backend-parity measurements (T1, T12, T13 of ``validation/plans/v3-003-acceptance.md``).

Pure NumPy/``statistics`` functions (no SciPy: the host runner has none) that compare two
:class:`~ionmc.simulation.Result` objects and return JSON-serialisable verdicts. They are used
by the tests and by ``validation/scripts/transport/run_suite.py``, so the numbers a test
asserts and the numbers a validation archive records come from the same code.

* T1 (:func:`compare_traces`): discrete trace columns exactly equal, continuous columns within
  ``rtol = atol = 1e-10``.
* T13 (:func:`compare_partition`): counters, tallies and deposit grids bit-identical (the grids
  are int64 fixed-point accumulators); the old relative bounds (1e-5 float32, 1e-12 float64, see
  :func:`deposit_agreement`) are reported as secondary numbers.
* T12 (:func:`t12_compare`): per-bin ``z = (a - b) / sqrt(se_a^2 + se_b^2)`` of two independent
  samples (distinct seeds) on the IDD and the lateral profiles for bins above 1 % of the maximum,
  the chi-square of the z profile with a batch-level permutation p-value (> 0.001; correlated
  bins invalidate the independent-bin Wilson-Hilferty approximation, which is reported only for
  information), the Bonferroni bound on ``max |z|`` and
  ``|z| < 3.5`` for the scalars (R80, total deposit, lateral sigma at 0.5 R). Standard errors
  come from the batch method of each run, so ``z`` is Student-t rather than normal for few
  batches; the number of batches of every run is recorded with the verdict.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from statistics import NormalDist
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.scoring import ScoringGrid
from ionmc.simulation import Result
from ionmc.transport.tally import COUNTER_NAMES, TRACE_COLUMNS, TRACE_N_DISCRETE

TRACE_RTOL = 1e-10
TRACE_ATOL = 1e-10
DEPOSIT_RTOL = {"float32": 1e-5, "float64": 1e-12}
DEPOSIT_FLOOR_FRACTION = 0.01
P_VALUE_MIN = 0.001
SCALAR_Z_MAX = 3.5
DEGENERATE_RTOL = 1e-9
DEGENERATE_SE_RTOL = 1e-8
DETERMINISTIC_RTOL = {"float32": 1e-5, "float64": 1e-12}
DOSE_FRACTION = 0.01
T12_DEPTHS = (0.25, 0.5, 0.9)


# -- T1: trajectory parity ---------------------------------------------------------------------
def compare_traces(
    reference: dict[str, Any],
    other: dict[str, Any],
    *,
    rtol: float = TRACE_RTOL,
    atol: float = TRACE_ATOL,
) -> dict[str, Any]:
    """T1 verdict for two ``diagnostics`` dictionaries that contain a trace.

    Discrete columns (``history, step, ix, iy, iz, reason, blocks, attempts``) must be equal on
    every step of every traced history, the end codes equal, and every continuous column within
    ``|a - b| <= atol + rtol |b|``. A different number of steps fails.
    """
    tr_a, tr_b = reference["trace"], other["trace"]
    n_a, n_b = len(tr_a["step"]), len(tr_b["step"])
    out: dict[str, Any] = {
        "steps_reference": n_a,
        "steps_other": n_b,
        "rtol": rtol,
        "atol": atol,
        "discrete_mismatches": {},
        "continuous_max_abs_diff": {},
        "continuous_violations": {},
    }
    if n_a != n_b:
        out["pass"] = False
        out["reason"] = "different number of trace steps"
        return out
    ok = True
    for i, name in enumerate(TRACE_COLUMNS):
        a, b = np.asarray(tr_a[name]), np.asarray(tr_b[name])
        if i < TRACE_N_DISCRETE:
            bad = int(np.count_nonzero(a != b))
            out["discrete_mismatches"][name] = bad
            ok &= bad == 0
        else:
            diff = np.abs(a - b)
            out["continuous_max_abs_diff"][name] = float(diff.max()) if n_a else 0.0
            viol = int(np.count_nonzero(~(diff <= atol + rtol * np.abs(b))))
            out["continuous_violations"][name] = viol
            ok &= viol == 0
    codes_equal = bool(
        np.array_equal(reference["trace_end_code"], other["trace_end_code"])
        and np.array_equal(reference["trace_end_history"], other["trace_end_history"])
    )
    out["end_codes_equal"] = codes_equal
    out["pass"] = bool(ok and codes_equal and n_a > 0)
    return out


# -- repeatability -----------------------------------------------------------------------------
def compare_runs_bitwise(a: Result, b: Result) -> dict[str, Any]:
    """Bit-for-bit comparison of two runs of one configuration (tallies, counters, deposit grids
    and the per-history end state: position, direction, energy, end code). Both runs must have
    been made with ``track_end_positions``. ``first_difference`` names the first differing
    history (or voxel) so that a nondeterminism can be located."""
    out: dict[str, Any] = {"identical": True, "first_difference": None}

    def note(what: str, **where: Any) -> None:
        if out["identical"]:
            out["first_difference"] = {"what": what, **where}
        out["identical"] = False

    ea, eb = a.energy_balance, b.energy_balance
    for f in ("initial_mev", "step_deposit_mev", "cutoff_mev", "escaped_mev", "truncated_mev"):
        if getattr(ea, f) != getattr(eb, f):
            note(f"tally {f}", a=getattr(ea, f), b=getattr(eb, f))
    if a.counters.as_dict() != b.counters.as_dict():
        note("counters", a=a.counters.as_dict(), b=b.counters.as_dict())
    for ga, gb in zip(a.grids, b.grids, strict=True):
        x, y = np.asarray(ga.batch_energy_mev), np.asarray(gb.batch_energy_mev)
        bad = np.argwhere(x != y)
        if bad.size:
            note(
                f"deposit grid {ga.name}",
                voxel=bad[0].tolist(),
                a=float(x[tuple(bad[0])]),
                b=float(y[tuple(bad[0])]),
                n_differing=int(len(bad)),
            )
    keys = ("end_position_mm", "end_direction", "end_energy_mev", "end_code")
    for k in keys:
        x, y = np.asarray(a.diagnostics[k]), np.asarray(b.diagnostics[k])
        differing = np.nonzero((x != y).reshape(len(x), -1).any(axis=1))[0]
        if len(differing):
            h = int(differing[0])
            note(
                f"end state {k}",
                history=h,
                a=np.atleast_1d(x[h]).tolist(),
                b=np.atleast_1d(y[h]).tolist(),
                n_histories_differing=int(len(differing)),
            )
    return out


# -- T13: partition invariance -----------------------------------------------------------------
def deposit_agreement(
    a: NDArray[np.float64], b: NDArray[np.float64], rtol: float
) -> dict[str, Any]:
    """Elementwise agreement of two deposit arrays: ``|a - b| <= rtol * max(|a|, f max|a|)``
    with ``f = 1 %`` (voxels below 1 % of the peak are compared with an absolute tolerance of
    ``rtol`` times 1 % of the peak, because float32 position/direction rounding perturbs
    the deposit relative to the accumulated magnitude, not to a near-empty voxel).

    The verdict carries the diagnostics of the worst voxel: its index (into the array, e.g.
    ``(batch, ix, iy, iz)``), both values, the number of violating voxels and the largest
    difference relative to the peak.
    """
    peak = float(np.max(np.abs(a))) if a.size else 0.0
    floor = DEPOSIT_FLOOR_FRACTION * peak
    scale = np.maximum(np.abs(a), floor)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(scale > 0.0, np.abs(a - b) / scale, 0.0)
    worst = np.unravel_index(int(np.argmax(rel)), rel.shape) if a.size else ()
    dmax = np.unravel_index(int(np.argmax(np.abs(a - b))), a.shape) if a.size else ()
    return {
        "max_relative_difference": float(rel.max()) if a.size else 0.0,
        "worst_voxel_index": [int(i) for i in worst],
        "worst_voxel_value_a": float(a[worst]) if a.size else 0.0,
        "worst_voxel_value_b": float(b[worst]) if a.size else 0.0,
        "n_violating_voxels": int(np.count_nonzero(rel > rtol)),
        "n_compared_voxels": int(a.size),
        "max_difference_over_peak": float(np.abs(a - b).max() / peak) if peak > 0 else 0.0,
        "max_difference_voxel_index": [int(i) for i in dmax],
        "peak": peak,
        "total_relative_difference": float(abs(a.sum() - b.sum()) / abs(a.sum()))
        if a.sum() != 0.0
        else 0.0,
        "rtol": rtol,
        "pass": bool(np.all(rel <= rtol)),
    }


def format_partition_verdict(v: dict[str, Any]) -> str:
    """One readable block of the numbers of a :func:`compare_partition` verdict (for assertion
    messages)."""
    lines = [
        f"precision={v['precision']} tallies_identical={v['tallies_identical']} "
        f"counters_identical={v['counters_identical']} edep_identical={v['edep_identical']} "
        f"pass={v['pass']}"
    ]
    for name, g in v["deposit"].items():
        lines.append(
            f"  grid {name}: identical={g['edep_identical']} "
            f"max_rel={g['max_relative_difference']:.3e} (rtol {g['rtol']:.0e}) at "
            f"{g['worst_voxel_index']} a={g['worst_voxel_value_a']:.9g} "
            f"b={g['worst_voxel_value_b']:.9g}; violating {g['n_violating_voxels']}/"
            f"{g['n_compared_voxels']}; max|d|/peak={g['max_difference_over_peak']:.3e} at "
            f"{g['max_difference_voxel_index']}; total_rel={g['total_relative_difference']:.3e}"
        )
    return "\n".join(lines)


def compare_partition(a: Result, b: Result) -> dict[str, Any]:
    """T13 verdict for two runs of the same configuration that differ only in how the histories
    were split (workers, chunk size): counters and tallies identical, deposit grids within
    ``DEPOSIT_RTOL[precision]``."""
    precision = a.precision
    ea, eb = a.energy_balance, b.energy_balance
    tallies_a = {
        "initial": ea.initial_mev,
        "step_deposit": ea.step_deposit_mev,
        "cutoff": ea.cutoff_mev,
        "escaped": ea.escaped_mev,
        "truncated": ea.truncated_mev,
        "unaccounted": ea.unaccounted_mev,
        **{f"outside_{i}": v for i, v in enumerate(ea.outside_mev)},
    }
    tallies_b = {
        "initial": eb.initial_mev,
        "step_deposit": eb.step_deposit_mev,
        "cutoff": eb.cutoff_mev,
        "escaped": eb.escaped_mev,
        "truncated": eb.truncated_mev,
        "unaccounted": eb.unaccounted_mev,
        **{f"outside_{i}": v for i, v in enumerate(eb.outside_mev)},
    }
    tallies_identical = tallies_a == tallies_b
    counters_identical = a.counters.as_dict() == b.counters.as_dict()
    grids = {}
    for ga, gb in zip(a.grids, b.grids, strict=True):
        arr_a, arr_b = np.asarray(ga.batch_energy_mev), np.asarray(gb.batch_energy_mev)
        g = deposit_agreement(arr_a, arr_b, DEPOSIT_RTOL[precision])
        g["edep_identical"] = bool(np.array_equal(arr_a, arr_b))
        grids[ga.name] = g
    return {
        "precision": precision,
        "tallies_identical": tallies_identical,
        "tallies_a": tallies_a,
        "tallies_b": tallies_b,
        "counters_identical": counters_identical,
        "counters": a.counters.as_dict(),
        "counter_names": list(COUNTER_NAMES),
        "deposit": grids,
        "edep_identical": all(g["edep_identical"] for g in grids.values()),
        "pass": bool(
            tallies_identical
            and counters_identical
            and all(g["edep_identical"] for g in grids.values())
        ),
    }


# -- T12: statistical parity -------------------------------------------------------------------
def wilson_hilferty_p(chi2: float, dof: int) -> float:
    """Upper-tail probability ``P(X >= chi2)`` of a chi-square variable (Wilson-Hilferty)."""
    if dof < 1:
        raise ValueError("dof must be >= 1")
    c = 2.0 / (9.0 * dof)
    z = ((chi2 / dof) ** (1.0 / 3.0) - (1.0 - c)) / math.sqrt(c)
    return 1.0 - NormalDist().cdf(z)


def bonferroni_z(n: int, alpha: float = P_VALUE_MIN) -> float:
    """Two-sided Bonferroni bound ``Phi^-1(1 - alpha / (2 n))`` on ``max |z|`` of ``n`` bins."""
    return NormalDist().inv_cdf(1.0 - alpha / (2.0 * n))


@dataclass(frozen=True)
class T12Layout:
    """Scoring layout of the statistical-parity runs: the grid names and slab depths."""

    idd: str
    slabs: tuple[tuple[str, float], ...]  # (grid name, z/R)
    range_mm: float
    bin_mm: float
    lateral_bin_mm: float


def t12_scoring_grids(
    range_mm: float,
    *,
    half_width_mm: float,
    lateral_bin_mm: float = 0.2,
    slab_mm: float = 1.0,
    depth_mm: float | None = None,
) -> tuple[tuple[ScoringGrid, ...], T12Layout]:
    """Scoring grids of T12: an IDD ``(1, 1, nz)`` of ``slab_mm`` bins and lateral slabs of
    ``slab_mm`` thickness at ``z/R`` in (0.25, 0.5, 0.9) with ``lateral_bin_mm`` bins (at most
    four grids, the engine's limit). ``max_step_mm`` of the run must not exceed the smallest
    spacing (``lateral_bin_mm``)."""
    depth = depth_mm if depth_mm is not None else 1.3 * range_mm
    nz = int(math.ceil(depth / slab_mm))
    n_lat = int(round(2.0 * half_width_mm / lateral_bin_mm))
    grids = [
        ScoringGrid(
            (-half_width_mm, -half_width_mm, 0.0),
            (2.0 * half_width_mm, 2.0 * half_width_mm, slab_mm),
            (1, 1, nz),
            name="idd",
        )
    ]
    slabs = []
    for frac in T12_DEPTHS:
        z0 = math.floor(frac * range_mm / slab_mm) * slab_mm  # slab aligned with the IDD bins
        name = f"lat{int(round(100 * frac)):02d}"
        grids.append(
            ScoringGrid(
                (-half_width_mm, -half_width_mm, z0),
                (lateral_bin_mm, lateral_bin_mm, slab_mm),
                (n_lat, n_lat, 1),
                name=name,
            )
        )
        slabs.append((name, frac))
    return tuple(grids), T12Layout("idd", tuple(slabs), range_mm, slab_mm, lateral_bin_mm)


@dataclass(frozen=True)
class T12Observables:
    """Per-batch observables of one run: arrays ``(B, m)`` and scalars ``(B,)``."""

    arrays: dict[str, NDArray[np.float64]]
    scalars: dict[str, NDArray[np.float64]]
    n_batches: int
    histories_per_batch: int = 1
    precision: str = "float64"


def _r80(profile: NDArray[np.float64], dz: float) -> float:
    """Depth [mm] of the distal 80 % point of an IDD (linear interpolation between bin centres)."""
    i_max = int(np.argmax(profile))
    level = 0.8 * profile[i_max]
    for i in range(i_max, len(profile) - 1):
        if profile[i + 1] < level <= profile[i]:
            f = (profile[i] - level) / (profile[i] - profile[i + 1])
            return float((i + 0.5 + f) * dz)
    return float("nan")


def projected_idd(
    batch_energy: NDArray[np.float64],
    grid: ScoringGrid,
    start: tuple[float, float, float],
    direction: tuple[float, float, float],
    *,
    bin_mm: float = 0.1,
    n_bins: int,
    fine_mm: float = 0.01,
) -> NDArray[np.float64]:
    """Integrated depth-dose of a beam along ``direction`` from a 3-D deposit grid: the energy of
    every voxel is distributed over the depth ``t = (r - start) . u`` it covers. A cubic-voxel cell
    of sizes ``d_i`` projects onto the convolution of boxes of widths ``|u_i| d_i`` centred on the
    projected voxel centre (exact for a uniform deposit inside the voxel); the result is binned in
    ``bin_mm`` bins. ``batch_energy`` has shape ``(B, nx, ny, nz)`` (per-batch values); returns
    ``(B, n_bins)``. Needs no rotated grid, so oblique beams use the same estimator as axis-aligned
    ones."""
    u = np.asarray(direction, dtype=np.float64)
    u = u / np.linalg.norm(u)
    nx, ny, nz = grid.shape
    c = [
        grid.origin_mm[a] + (np.arange(grid.shape[a]) + 0.5) * grid.spacing_mm[a] - start[a]
        for a in range(3)
    ]
    t = (
        c[0][:, None, None] * u[0] + c[1][None, :, None] * u[1] + c[2][None, None, :] * u[2]
    ).reshape(-1)
    n_fine = int(round(n_bins * bin_mm / fine_mm))
    idx = np.floor(t / fine_mm).astype(np.int64)
    keep = (idx >= 0) & (idx < n_fine)
    kernel = np.ones(1)
    for a in range(3):
        width = abs(u[a]) * grid.spacing_mm[a]
        m = int(round(width / fine_mm))
        if m >= 2:
            kernel = np.convolve(kernel, np.full(m, 1.0 / m))
    per = int(round(bin_mm / fine_mm))
    out = np.zeros((batch_energy.shape[0], n_bins))
    for b in range(batch_energy.shape[0]):
        w = batch_energy[b].reshape(-1)
        fine = np.bincount(idx[keep], weights=w[keep], minlength=n_fine)[:n_fine]
        smooth = np.convolve(fine, kernel, mode="same") if kernel.size > 1 else fine
        out[b] = smooth[: n_bins * per].reshape(n_bins, per).sum(axis=1)
    return out


def r80_of(profile: NDArray[np.float64], bin_mm: float) -> float:
    """Distal 80 % depth [mm] of an IDD (the estimator of the T9/T12 observables)."""
    return _r80(profile, bin_mm)


def t12_observables(result: Result, layout: T12Layout) -> T12Observables:
    """Batch-wise IDD, lateral profiles and scalars (R80, total deposit, lateral sigma at 0.5 R)."""
    grids = {g.name: np.asarray(g.batch_energy_mev) for g in result.grids}
    return t12_observables_from_grids(
        grids,
        result.requested_config.scoring,
        layout,
        result.n_histories // result.n_batches,
        result.precision,
    )


def t12_observables_from_grids(
    batch_energy: dict[str, NDArray[np.float64]],
    scoring: tuple[ScoringGrid, ...],
    layout: T12Layout,
    histories_per_batch: int,
    precision: str = "float64",
) -> T12Observables:
    """:func:`t12_observables` from per-batch per-primary grid arrays ``(B, nx, ny, nz)`` by name
    (also used on samples merged from saved parts)."""
    idd_b = np.asarray(batch_energy[layout.idd])  # (B, 1, 1, nz)
    b = idd_b.shape[0]
    idd = idd_b.reshape(b, -1)
    arrays = {"idd": idd}
    scalars = {
        "total_deposit_mev": idd.sum(axis=1),
        "r80_mm": np.array([_r80(idd[k], layout.bin_mm) for k in range(b)]),
    }
    for name, frac in layout.slabs:
        slab = np.asarray(batch_energy[name])  # (B, nx, ny, 1)
        prof = slab.sum(axis=(2, 3))  # (B, nx): lateral x profile summed over y
        arrays[name] = prof
        if abs(frac - 0.5) < 1e-9:
            nx = prof.shape[1]
            grid = next(x for x in scoring if x.name == name)
            xc = grid.origin_mm[0] + (np.arange(nx) + 0.5) * grid.spacing_mm[0]
            tot = prof.sum(axis=1)
            mean = (prof * xc).sum(axis=1) / tot
            var = (prof * (xc[None, :] - mean[:, None]) ** 2).sum(axis=1) / tot
            scalars["sigma_lat_05R_mm"] = np.sqrt(var)
    return T12Observables(arrays, scalars, b, histories_per_batch, precision)


def _mean_se(x: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    b = x.shape[0]
    return x.mean(axis=0), x.std(axis=0, ddof=1) / math.sqrt(b)


def scalar_z(xa: NDArray[np.float64], xb: NDArray[np.float64]) -> dict[str, float]:
    """``z = (mean_a - mean_b) / sqrt(se_a^2 + se_b^2)`` of two samples of batch values (standard
    errors of the batch method). Values equal within ``DEGENERATE_RTOL`` (a scalar fixed by
    energy conservation has no variance) give ``z = 0``; a difference without variance gives
    ``inf``."""
    ma, sa = float(xa.mean()), float(xa.std(ddof=1)) / math.sqrt(len(xa))
    mb, sb = float(xb.mean()), float(xb.std(ddof=1)) / math.sqrt(len(xb))
    se = math.sqrt(sa**2 + sb**2)
    diff = ma - mb
    if abs(diff) <= DEGENERATE_RTOL * max(abs(ma), abs(mb)):
        z = 0.0
    else:
        z = diff / se if se > 0.0 else math.inf
    return {"a": ma, "b": mb, "se_a": sa, "se_b": sb, "z": z}


PERMUTATIONS = 2000
PERMUTATION_SEED = 20261004


def _chi2_profile(
    xa: NDArray[np.float64], xb: NDArray[np.float64]
) -> tuple[float, float, NDArray[np.float64]]:
    """``(chi2, max|z|, z)`` of the per-bin difference of two samples of batch means ``(B, m)``,
    with the standard errors of the batch method."""
    ma, sa = _mean_se(xa)
    mb, sb = _mean_se(xb)
    se = np.sqrt(sa**2 + sb**2)
    diff = ma - mb
    with np.errstate(divide="ignore", invalid="ignore"):
        z: NDArray[np.float64] = np.where(se > 0.0, diff / se, np.where(diff == 0.0, 0.0, np.inf))
    return float((z**2).sum()), float(np.abs(z).max()) if z.size else 0.0, z


def permutation_p_value(
    xa: NDArray[np.float64],
    na: float,
    xb: NDArray[np.float64],
    nb: float,
    *,
    n_perm: int = PERMUTATIONS,
    seed: int = PERMUTATION_SEED,
) -> float:
    """Batch-level studentized permutation p-value of the profile chi-square.

    ``xa`` (``Ba, m``) and ``xb`` (``Bb, m``) are per-batch per-primary means of ``m`` bins,
    ``na`` and ``nb`` the histories per batch of each sample. With the pooled mean
    ``mu`` (weighted by the batch sizes) the residuals ``r_b = (m_b - mu) sqrt(n_b)`` of all
    ``Ba + Bb`` batches are permuted across the two samples, ``m*_b = mu + r*_b / sqrt(n_b)``
    (the original size ``n_b`` of each batch) are rebuilt, and the same chi-square of the two
    sample means is recomputed. The p-value is ``(1 + #{chi2* >= chi2}) / (1 + n_perm)``.

    Assumption: under the null hypothesis the per-primary variance of a bin is the same in both
    samples, so that ``r_b`` are exchangeable; unlike the Wilson-Hilferty approximation this
    makes no assumption that the bins are independent (neighbouring bins of a profile are
    correlated because every history deposits in many bins). The generator is
    ``numpy.random.default_rng(seed)``; ``seed`` and ``n_perm`` are recorded by the caller.
    """
    ba, bb = xa.shape[0], xb.shape[0]
    sizes = np.concatenate([np.full(ba, float(na)), np.full(bb, float(nb))])
    m = np.concatenate([xa, xb], axis=0)
    mu = (sizes[:, None] * m).sum(axis=0) / sizes.sum()
    resid = (m - mu) * np.sqrt(sizes)[:, None]
    observed = _chi2_profile(xa, xb)[0]
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(ba + bb)
        star = mu + resid[perm] / np.sqrt(sizes)[:, None]
        if _chi2_profile(star[:ba], star[ba:])[0] >= observed:
            count += 1
    return (1.0 + count) / (1.0 + n_perm)


def t12_compare(
    a: T12Observables,
    b: T12Observables,
    *,
    n_perm: int = PERMUTATIONS,
    seed: int = PERMUTATION_SEED,
) -> dict[str, Any]:
    """T12 verdict for two independent samples (see the module docstring).

    For every profile the frozen statistic (chi-square of the per-bin ``z``, bins above 1 % of
    the maximum) is calibrated with :func:`permutation_p_value` (``p_value``, used for the
    verdict); the Wilson-Hilferty value is reported as ``p_value_wilson_hilferty`` only. The
    bound on ``max |z|`` is the frozen Bonferroni bound."""
    out: dict[str, Any] = {
        "n_batches": [a.n_batches, b.n_batches],
        "permutation": {"n_perm": n_perm, "seed": seed},
        "arrays": {},
        "scalars": {},
    }
    ok = True
    for i, name in enumerate(a.arrays):
        ma, _ = _mean_se(a.arrays[name])
        mb, _ = _mean_se(b.arrays[name])
        ref = 0.5 * (ma + mb)
        sel = ref > DOSE_FRACTION * ref.max()
        n = int(sel.sum())
        chi2, _zmax_all, z_sel = _chi2_profile(a.arrays[name][:, sel], b.arrays[name][:, sel])
        # defined-value rule: only bins that both samples support (at least max(2, ceil(B/2))
        # batches with a nonzero deposit) enter max|z|; the others have unreliable standard errors
        need_a, need_b = max(2, -(-a.n_batches // 2)), max(2, -(-b.n_batches // 2))
        sup_a = (a.arrays[name][:, sel] > 0.0).sum(axis=0)
        sup_b = (b.arrays[name][:, sel] > 0.0).sum(axis=0)
        supported = (sup_a >= need_a) & (sup_b >= need_b)
        n_supported = int(supported.sum())
        zmax = float(np.abs(z_sel[supported]).max()) if n_supported else 0.0
        worst_all = int(np.argmax(np.abs(z_sel))) if z_sel.size else 0
        sel_idx = np.nonzero(sel)[0]
        p_wh = wilson_hilferty_p(chi2, n) if n >= 1 else 1.0
        p = (
            permutation_p_value(
                a.arrays[name][:, sel],
                a.histories_per_batch,
                b.arrays[name][:, sel],
                b.histories_per_batch,
                n_perm=n_perm,
                seed=seed + i,
            )
            if n >= 1
            else 1.0
        )
        bound = bonferroni_z(n_supported) if n_supported else float("inf")
        passed = bool(p > P_VALUE_MIN and zmax < bound)
        out["arrays"][name] = {
            "n_bins": n,
            "n_supported_bins": n_supported,
            "chi2": chi2,
            "p_value": p,
            "p_value_wilson_hilferty": p_wh,
            "max_abs_z": zmax,
            "max_abs_z_all_bins": float(np.abs(z_sel[worst_all])) if z_sel.size else 0.0,
            "worst_bin_all": {
                "index": int(sel_idx[worst_all]) if z_sel.size else None,
                "batches_with_deposit": [int(sup_a[worst_all]), int(sup_b[worst_all])]
                if z_sel.size
                else None,
                "required": [need_a, need_b],
                "supported": bool(supported[worst_all]) if z_sel.size else None,
            },
            "bonferroni_bound": bound,
            "pass": passed,
        }
        ok &= passed
    for name in a.scalars:
        sc = scalar_z(a.scalars[name], b.scalars[name])
        scale = max(abs(sc["a"]), abs(sc["b"]))
        if scale > 0.0 and max(sc["se_a"], sc["se_b"]) <= DEGENERATE_SE_RTOL * scale:
            # a deterministic scalar (energy conservation makes the total deposit exact): z is
            # meaningless; compare with the precision bound of the less precise sample (the
            # deterministic T4 bound: 1e-5 relative for float32, 1e-12 for float64 / python)
            tol = max(DETERMINISTIC_RTOL[a.precision], DETERMINISTIC_RTOL[b.precision])
            rel = abs(sc["a"] - sc["b"]) / scale
            passed = bool(rel <= tol)
            out["scalars"][name] = {
                **sc,
                "rule": "deterministic_precision_bound",
                "relative_difference": rel,
                "bound": tol,
                "pass": passed,
            }
        else:
            passed = bool(abs(sc["z"]) < SCALAR_Z_MAX)
            out["scalars"][name] = {**sc, "rule": "z", "bound": SCALAR_Z_MAX, "pass": passed}
        ok &= passed
    out["pass"] = bool(ok)
    return out


def t12_config(
    *,
    energy_mev: float,
    backend: str,
    precision: str,
    seed: int,
    n_histories: int,
    n_batches: int,
    workers: int = 1,
    lateral_bin_mm: float = 0.2,
    half_width_mm: float = 20.0,
    max_step_mm: float | None = None,
    timeout_s: float | None = None,
    chunk_histories: int | None = None,
    memory_budget_bytes: int | None = None,
) -> tuple[Any, T12Layout]:
    """Configuration of one T12 sample: a pencil beam of protons in a water box, all physics on,
    offline analytic (Bethe, I = 78 eV) stopping, the T12 scoring grids. Each sample must use
    its own ``seed``. ``max_step_mm`` defaults to ``lateral_bin_mm``."""
    from ionmc.config import (
        DEFAULT_CHUNK_HISTORIES,
        PhysicsOptions,
        RunOptions,
        SimulationConfig,
    )
    from ionmc.geometry import BoxPhantom
    from ionmc.materials import WATER
    from ionmc.physics.projectiles import PROTON
    from ionmc.physics.stopping import BetheStoppingSource
    from ionmc.sources import PencilBeamSource
    from ionmc.transport.tables import TransportTables

    stopping = BetheStoppingSource()
    tab = TransportTables.from_stopping_tables([stopping.table(WATER, PROTON)])
    range_mm = tab.range_g_cm2(0, energy_mev) * 10.0 / WATER.density_g_cm3
    grids, layout = t12_scoring_grids(
        range_mm, half_width_mm=half_width_mm, lateral_bin_mm=lateral_bin_mm
    )
    depth = 1.3 * range_mm
    margin = half_width_mm + 10.0
    cfg = SimulationConfig(
        source=PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), energy_mev),
        geometry=BoxPhantom((-margin, -margin, 0.0), (2 * margin, 2 * margin, depth), WATER),
        scoring=grids,
        physics=PhysicsOptions(
            nuclear=False,
            stopping=stopping,
            max_step_mm=max_step_mm if max_step_mm is not None else lateral_bin_mm,
        ),
        run=RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision=precision,  # type: ignore[arg-type]
            seed=seed,
            n_histories=n_histories,
            n_batches=n_batches,
            cpu_workers=workers,
            worker_timeout_s=timeout_s,
            chunk_histories=chunk_histories or DEFAULT_CHUNK_HISTORIES,
            **({"memory_budget_bytes": memory_budget_bytes} if memory_budget_bytes else {}),  # type: ignore[arg-type]
        ),
    )
    return cfg, layout


def t12_sample(**kwargs: Any) -> tuple[Result, T12Layout]:
    """Run one T12 sample (see :func:`t12_config`)."""
    from ionmc.simulation import Simulation

    cfg, layout = t12_config(**kwargs)
    return Simulation(cfg).run(), layout
