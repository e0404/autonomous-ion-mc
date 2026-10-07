"""Scoring grids, voxel masses and the batch estimator (float64 numpy).

A scoring grid is an independent regular grid (same conventions as
:class:`ionmc.geometry.VoxelGeometry`; it need not align with the transport grid). The energy of
a step is deposited by track-length apportioning with a linear stopping-power ramp: along the
hinge path (both legs) the deposit density is taken to vary linearly from S(E_old) to S(E_new), and
each voxel receives the integral of that density over the path piece inside it (a per-grid
incremental voxel walk, see ``ionmc.transport.reference``). The result has a residual step
dependence of second order (the curvature of the stopping power along the step), instead of the
aliasing of a point deposit with the voxel edges or the first-order error of a uniform deposit;
the energy left at the cutoff is a point deposit at the end point. Deposits are quantized into int64
fixed-point accumulators (``ionmc.transport.tally``). The mass of a scoring
voxel is the exact overlap integral of the piecewise-constant density of the geometry over
the voxel (the geometry is vacuum outside, so partially covered voxels have a correspondingly
smaller mass). Dose uses 1 MeV/g = 1.602176634e-10 Gy.

Batch estimator: histories are assigned to batches by ``history_index mod n_batches``; with
per-batch energy sums ``E_b`` [MeV] and ``N/B`` histories per batch, ``x_b = E_b / (N/B)``
is the per-primary estimate of batch ``b``; the mean is ``mean(x_b)`` and the variance of
the mean ``sum (x_b - mean)^2 / (B (B - 1))``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from ionmc._validate import fail, int_triple, real, triple
from ionmc.geometry import VoxelGeometry

MEV_PER_G_TO_GY = 1.602176634e-10
MAX_SCORING_GRIDS = 4
TALLY_QUANTITIES = (
    "edep",
    "dose",
    "fluence",
    "let_t",
    "let_d",
    "let_d_eps",
    "lookup_sum",
    "lookup_dose_avg",
    "fluence_spectrum",
)
GENERATION_CHOICES = ("all", "primary", "secondary")


@dataclass(frozen=True)
class ScoringGrid:
    """A regular scoring grid: ``origin_mm`` (corner of voxel (0,0,0)), ``spacing_mm`` (> 0),
    ``shape`` ``(nx, ny, nz)``, C-ordered ``[ix, iy, iz]``; ``name`` labels the result."""

    origin_mm: tuple[float, float, float]
    spacing_mm: tuple[float, float, float]
    shape: tuple[int, int, int]
    name: str = "dose"

    def __post_init__(self) -> None:
        object.__setattr__(self, "origin_mm", triple("origin_mm", self.origin_mm))
        object.__setattr__(self, "spacing_mm", triple("spacing_mm", self.spacing_mm, positive=True))
        object.__setattr__(self, "shape", int_triple("shape", self.shape))
        if not isinstance(self.name, str) or not self.name:
            raise fail("name must be a non-empty string")

    @property
    def n_voxels(self) -> int:
        """Number of voxels."""
        return self.shape[0] * self.shape[1] * self.shape[2]

    def edges_mm(self, axis: int) -> NDArray[np.float64]:
        """Voxel boundary coordinates along ``axis`` [mm], ``shape[axis] + 1`` values."""
        return self.origin_mm[axis] + self.spacing_mm[axis] * np.arange(
            self.shape[axis] + 1, dtype=np.float64
        )


@dataclass(frozen=True)
class TallyRequest:
    """A request for one scored quantity on one scoring grid (decision 0040).

    ``quantity``: ``edep`` (MeV), ``dose`` (dose to medium), ``fluence`` (sum of path lengths per
    voxel volume, mm^-2), ``let_t`` / ``let_d`` (track- and dose-averaged LET in water, keV/um),
    ``let_d_eps`` (diagnostic, deposit-weighted), ``lookup_sum`` (sum of ``eps f``, MeV times the
    unit of ``f``), ``lookup_dose_avg`` (``lookup_sum / E_step``, the deposit-weighted average of
    the lookup ``f``) and ``fluence_spectrum`` (path length per energy bin per nucleon, with an
    underflow and an overflow bin). ``species`` are registry names (None: all transported charged
    species; for ``edep``/``dose`` also the non-transported pseudo-species so the total closes);
    ``generation`` selects primaries, secondaries or all. ``lookup`` names a
    :class:`~ionmc.lookup.LookupTable` of ``SimulationConfig.lookups`` (required for the
    ``lookup_*`` quantities only); ``energy_edges_mev_per_u`` are the uniform (linear or log) bin
    edges of ``fluence_spectrum`` (required for it only).

    ``let_medium`` and ``dose_reference`` exist so that unsupported choices fail closed instead of
    being unavailable to ask for: only ``"water"`` and ``"medium"`` are supported (the LET medium
    is part of the definition; dose-to-water needs a later channel kind).
    """

    name: str
    grid: str
    quantity: Literal[
        "edep",
        "dose",
        "fluence",
        "let_t",
        "let_d",
        "let_d_eps",
        "lookup_sum",
        "lookup_dose_avg",
        "fluence_spectrum",
    ]
    species: tuple[str, ...] | None = None
    generation: Literal["all", "primary", "secondary"] = "all"
    lookup: str | None = None
    energy_edges_mev_per_u: tuple[float, ...] | None = None
    let_medium: str = "water"
    dose_reference: str = "medium"

    def __post_init__(self) -> None:
        for f in ("name", "grid", "let_medium", "dose_reference"):
            v = getattr(self, f)
            if not isinstance(v, str) or not v:
                raise fail(f"tally request field {f!r} must be a non-empty string")
        if self.quantity not in TALLY_QUANTITIES:
            raise fail(f"tally quantity must be one of {TALLY_QUANTITIES}, got {self.quantity!r}")
        if self.generation not in GENERATION_CHOICES:
            raise fail(f"generation must be one of {GENERATION_CHOICES}, got {self.generation!r}")
        if self.species is not None:
            sp = self.species
            if (
                not isinstance(sp, tuple | list)
                or len(sp) == 0
                or not all(isinstance(x, str) for x in sp)
                or len(set(sp)) != len(sp)
            ):
                raise fail("species must be None or a non-empty tuple of distinct names")
            object.__setattr__(self, "species", tuple(sp))
        needs_lookup = self.quantity in ("lookup_sum", "lookup_dose_avg")
        if needs_lookup != (self.lookup is not None):
            raise fail(
                f"tally {self.name!r}: 'lookup' is required for lookup_* quantities and "
                "forbidden for all others"
            )
        if self.lookup is not None and (not isinstance(self.lookup, str) or not self.lookup):
            raise fail("lookup must be a non-empty table name")
        has_edges = self.energy_edges_mev_per_u is not None
        if (self.quantity == "fluence_spectrum") != has_edges:
            raise fail(
                f"tally {self.name!r}: 'energy_edges_mev_per_u' is required for fluence_spectrum "
                "and forbidden for all other quantities"
            )
        if self.energy_edges_mev_per_u is not None:
            e = tuple(
                real(f"energy_edges_mev_per_u[{i}]", v, positive=True)
                for i, v in enumerate(self.energy_edges_mev_per_u)
            )
            if len(e) < 2 or any(b <= a for a, b in zip(e, e[1:], strict=False)):
                raise fail("energy_edges_mev_per_u must be >= 2 strictly increasing values")
            if not all(math.isfinite(v) for v in e):  # pragma: no cover - real() checks
                raise fail("energy edges must be finite")
            object.__setattr__(self, "energy_edges_mev_per_u", e)


def overlap_matrix(
    edges_a: NDArray[np.float64], edges_b: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Overlap length [mm] ``O[i, j]`` of interval ``i`` of ``edges_a`` and ``j`` of ``edges_b``."""
    lo = np.maximum(edges_a[:-1, None], edges_b[None, :-1])
    hi = np.minimum(edges_a[1:, None], edges_b[None, 1:])
    out: NDArray[np.float64] = np.clip(hi - lo, 0.0, None)
    return out


def voxel_mass_g(grid: ScoringGrid, geometry: VoxelGeometry) -> NDArray[np.float64]:
    """Mass [g] of every scoring voxel, shape ``grid.shape``: the exact overlap integral of the
    geometry density [g/cm3] over the voxel (volumes converted from mm3 to cm3)."""
    ox = overlap_matrix(grid.edges_mm(0), _geometry_edges(geometry, 0))
    oy = overlap_matrix(grid.edges_mm(1), _geometry_edges(geometry, 1))
    oz = overlap_matrix(grid.edges_mm(2), _geometry_edges(geometry, 2))
    rho = geometry.densities_g_cm3()
    mass = np.einsum("ia,jb,kc,abc->ijk", ox, oy, oz, rho, optimize=True)
    out: NDArray[np.float64] = mass * 1.0e-3
    return out


def _geometry_edges(geometry: VoxelGeometry, axis: int) -> NDArray[np.float64]:
    edges: NDArray[np.float64] = geometry.origin_mm[axis] + geometry.spacing_mm[axis] * np.arange(
        geometry.shape[axis] + 1, dtype=np.float64
    )
    if axis == 2 and geometry.z_exit_mm is not None:  # voxels beyond the exit plane have no mass
        edges = np.minimum(edges, geometry.z_exit_mm)
    return edges


@dataclass(frozen=True, eq=False)
class BatchStatistics:
    """Batch-method estimate per voxel: ``mean`` and ``variance_of_mean`` (per primary,
    same units as the input divided by histories per batch), ``n_nonzero`` batches."""

    mean: NDArray[np.float64]
    variance_of_mean: NDArray[np.float64]
    n_nonzero: NDArray[np.int32]


def reduce_batches(batch_sums: NDArray[np.float64], histories_per_batch: int) -> BatchStatistics:
    """Reduce per-batch sums ``[B, n_voxels]`` (float64) to per-primary mean and variance."""
    s = np.asarray(batch_sums, dtype=np.float64)
    if s.ndim != 2 or s.shape[0] < 2:
        raise ValueError("batch sums must have shape [B >= 2, n_voxels]")
    if histories_per_batch < 1:
        raise ValueError("histories_per_batch must be >= 1")
    if not np.all(np.isfinite(s)):
        raise ValueError("batch sums must be finite")
    b = s.shape[0]
    x = s / float(histories_per_batch)
    mean = x.mean(axis=0)
    var = ((x - mean) ** 2).sum(axis=0) / (b * (b - 1))
    return BatchStatistics(mean, var, (s > 0.0).sum(axis=0).astype(np.int32))


@dataclass(frozen=True, eq=False)
class RatioStatistics:
    """Delta-method estimate of ``Xbar / Ybar`` per bin (see :func:`reduce_ratio`): ``mean`` and
    ``variance_of_mean`` are NaN where ``defined_mask`` is False (never 0); ``n_nonzero`` counts the
    batches with ``y_b > 0``."""

    mean: NDArray[np.float64]
    variance_of_mean: NDArray[np.float64]
    defined_mask: NDArray[np.bool_]
    n_nonzero: NDArray[np.int32]


def min_defined_batches(n_batches: int) -> int:
    """Smallest ``n_nonzero(y)`` of a defined ratio bin: ``max(2, ceil(B / 2))`` (decision 0040,
    section 5)."""
    return max(2, math.ceil(n_batches / 2))


def reduce_ratio(
    num_batch_sums: NDArray[np.float64], den_batch_sums: NDArray[np.float64]
) -> RatioStatistics:
    """Ratio of channel means with the delta-method variance (decision 0040, section 5).

    ``num_batch_sums`` and ``den_batch_sums`` are per-batch sums ``[B, n_bins]`` of the numerator
    and denominator channels (``B >= 2``; the common factor 1 / histories-per-batch cancels).
    With the batch means ``xbar``, ``ybar`` and ``R = xbar / ybar``::

        Var R = [Var xbar - 2 R Cov(xbar, ybar) + R^2 Var ybar] / ybar^2
        Cov   = sum_b (x_b - xbar)(y_b - ybar) / (B (B - 1))

    The biased mean of per-batch ratios is not used. A bin is defined iff ``ybar > 0`` and
    ``n_nonzero(y) >= max(2, ceil(B / 2))``; undefined bins are NaN. The variance is clipped at 0
    against rounding (the quadratic form is non-negative in exact arithmetic).
    """
    x = np.asarray(num_batch_sums, dtype=np.float64)
    y = np.asarray(den_batch_sums, dtype=np.float64)
    if x.ndim != 2 or x.shape != y.shape or x.shape[0] < 2:
        raise ValueError("batch sums must have the same shape [B >= 2, n_bins]")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("batch sums must be finite")
    b = x.shape[0]
    xbar = x.mean(axis=0)
    ybar = y.mean(axis=0)
    dx, dy = x - xbar, y - ybar
    denom = float(b * (b - 1))
    var_x = (dx * dx).sum(axis=0) / denom
    var_y = (dy * dy).sum(axis=0) / denom
    cov = (dx * dy).sum(axis=0) / denom
    n_nonzero = (y > 0.0).sum(axis=0).astype(np.int32)
    defined = (ybar > 0.0) & (n_nonzero >= min_defined_batches(b))
    safe_y = np.where(defined, ybar, 1.0)
    ratio = xbar / safe_y
    var = (var_x - 2.0 * ratio * cov + ratio * ratio * var_y) / (safe_y * safe_y)
    var = np.maximum(var, 0.0)
    return RatioStatistics(
        mean=np.where(defined, ratio, np.nan),
        variance_of_mean=np.where(defined, var, np.nan),
        defined_mask=defined,
        n_nonzero=n_nonzero,
    )


@dataclass(frozen=True, eq=False)
class QuantityResult:
    """One requested tally on its scoring grid (``GridResult.quantities[name]``).

    ``mean`` and ``std`` (standard error of the mean) have the grid shape, with a trailing bin axis
    of ``n_bins + 2`` for a ``fluence_spectrum`` (underflow, bins, overflow). Linear quantities
    are per primary; ratios are ratios of per-primary channel means. ``defined_mask`` is True where
    the value is defined (a ratio: ``ybar > 0`` and ``n_nonzero >= max(2, ceil(B/2))``; a linear
    quantity: positive mean, dose also positive mass); undefined ratio bins are NaN, never 0.
    ``n_nonzero`` counts the batches with a nonzero denominator (numerator for a linear quantity).
    ``rounding_bound`` is the deterministic fixed-point bound of ``mean``: every quantized piece
    adds at most ``q / 2``, so a channel's bound is ``(n_step [+ n_local]) q / 2`` per voxel and
    primary with the exact piece counts of the automatic channels N (class "step") and N_local
    (class "local"), taken according to the channel's class mask (both for edep and dose); ratios
    combine numerator and denominator bounds by first-order propagation. It bounds the
    quantization error of the stored sums, not the statistical error; ``quantum_exponents`` the
    ``k_c`` of the numerator and denominator channels; ``lookup`` the provenance of the lookup
    table (or None).
    """

    name: str
    quantity: str
    kind: str
    units: str
    definition: str
    mean: NDArray[np.float64]
    std: NDArray[np.float64]
    defined_mask: NDArray[np.bool_]
    n_nonzero: NDArray[np.int32]
    rounding_bound: NDArray[np.float64]
    quantum_exponents: tuple[int, ...]
    lookup: dict[str, object] | None = None
