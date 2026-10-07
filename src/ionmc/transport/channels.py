"""Channel compiler: scoring requests to deduplicated linear channels (decision 0040).

Every :class:`~ionmc.scoring.TallyRequest` is a combination of *linear channels*, each the sum over
scoring pieces of one per-piece increment ``x`` (kinds below), restricted to a scoring grid, a set
of species, a generation range and a class mask ("step" pieces of steps with ``s_act > 0``;
"local" point deposits: cutoff energy and ``eps > 0`` at ``s_act = 0``):

=====  ==========================  ===========================================================
kind   increment ``x`` per piece   used by
=====  ==========================  ===========================================================
E      ``eps_p``                   edep, dose, ``E_step`` (denominator), ``edep_excluded_from_let``
L      ``l_p``                     fluence, ``LET_t`` (denominator)
LS     ``l_p S_p``                 ``LET_t`` (numerator), ``LET_d`` (denominator)
LS2    ``l_p (S_p^2 + k^2 l_p^2/12)``  ``LET_d`` (numerator)
ES     ``eps_p S_p``               ``LET_d^eps`` (numerator)
FE     ``eps_p f(x_p)``            ``lookup_sum``, ``lookup_dose_avg`` (numerator)
FL     ``l_p`` per energy bin      ``fluence_spectrum`` (f = bin indicator, with under/overflow)
N      1                           piece count (automatic, bounds the fixed-point error);
                                   one channel counts class "step" pieces, one class "local"
=====  ==========================  ===========================================================

Ratios are ``LET_t = LS/L``, ``LET_d = LS2/LS``, ``LET_d^eps = ES/E_step`` and
``lookup_dose_avg = FE/E_step``; their reduction (delta method) is ``reduce_ratio`` (later step).
The automatic channels per used grid are ``N`` (class "step"), ``N_local`` (class "local", all
species including the pseudo-species: it counts the local point deposits that the E channels
quantize) and ``edep_excluded_from_let`` (E, class "local", all species including the
pseudo-species). Together the two counts give the deterministic rounding bound of every
channel according to its class mask.

Channels live in one int64 array ``acc[B, sum_c size_c]``; channel ``c`` owns a contiguous block
of ``size_c`` columns (``n_voxels``, times ``n_bins + 2`` for a spectrum) at ``offset_c``. An
increment is quantized as ``floor(x 2^k_c + 1/2)``. ``k_c`` follows the capacity rule of
decision 0040 section 4: energy channels use the fixed ``k = 30`` (as the qualified ``edep``
array), the others ``k = min(40, floor(62 - log2(hpb B_c)))`` with the per-history bound ``B_c``
of their kind; a quantum above the precision floor ``2^-16 u_c`` fails closed.

Everything here is host-side planning in float64/numpy and raises ``UnsupportedCombinationError``
before any transport.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc._frozen import freeze_array
from ionmc._validate import fail
from ionmc.geometry import VoxelGeometry
from ionmc.lookup import LookupTable
from ionmc.physics.projectiles import Projectile
from ionmc.physics.stopping import StoppingTable
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.species import (
    N_SPECIES_SLOTS,
    NUCLEAR_LOCAL,
    species_by_name,
    species_of_projectile,
    transported_ids,
)
from ionmc.transport.tables import TransportTables
from ionmc.transport.tally import MAX_QUANTA, QUANTUM_MEV

KINDS = ("E", "L", "LS", "LS2", "ES", "FE", "FL", "N")
CLASS_STEP = 1
CLASS_LOCAL = 2
GEN_MAX = 2**31 - 1
MAX_CHANNELS = 256
MAX_SPECTRUM_BINS = 4096
K_MAX = 40
CAPACITY_EXPONENT = 62
FLOOR_EXPONENT = -16
STRAGGLING_MARGIN = 1.25
STEP_LENGTH_MARGIN = 1.0 + 1.0e-6
"""Relative margin on the informative truncation bound ``max_steps * max_step_mm`` (float32
rounding of the step limit)."""
E_QUANTUM_EXPONENT = 30  # QUANTUM_MEV = 2**-30
FE_F_MIN_EXPONENT = -60  # admissible largest lookup value of a request: [2^-60, 2^60]
FE_F_MAX_EXPONENT = 60
EXCLUDED_CHANNEL_NAME = "edep_excluded_from_let"
PIECE_COUNT_NAME = "scoring_pieces"
LOCAL_PIECE_COUNT_NAME = "scoring_pieces_local"
# power of the stopping power S in the 1 mm entrance-piece scale u_c = 1 mm S_ref^j
_FLOOR_POWER = {"E": 1, "FE": 1, "L": 0, "FL": 0, "LS": 1, "LS2": 2, "ES": 2}


@dataclass(frozen=True)
class SpectrumSpec:
    """Uniform bin edges of a fluence spectrum: ``n_bins`` bins plus underflow (bin 0) and
    overflow (bin ``n_bins + 1``); ``log`` selects uniform in ``ln E``."""

    edges: tuple[float, ...]
    log: bool

    @property
    def n_bins(self) -> int:
        return len(self.edges) - 1


@dataclass(frozen=True)
class Channel:
    """One linear channel (see the module docstring). ``species`` is the int8 row of the match
    matrix; ``bound_per_history`` is ``B_c`` in the channel unit; ``residual_column`` is the index
    among the per-history residual columns (-1 for the exact N channel)."""

    kind: str
    grid: int
    class_mask: int
    species: tuple[int, ...]
    gen_lo: int
    gen_hi: int
    lookup: int
    spectrum: SpectrumSpec | None
    size: int
    offset: int
    k: int
    bound_per_history: float
    residual_column: int

    @property
    def quantum(self) -> float:
        """``2^-k`` in the channel unit."""
        return math.ldexp(1.0, -self.k)


@dataclass(frozen=True)
class QuantityDef:
    """How a request is assembled from channels: ``numerator`` (and ``denominator`` for a ratio)
    are channel indices; ``kind`` is ``"linear"`` or ``"ratio"``."""

    name: str
    request: TallyRequest
    kind: str
    numerator: int
    denominator: int | None
    units: str
    definition: str


@dataclass(frozen=True, eq=False)
class ChannelPlan:
    """The compiled channels: ``channels`` ordered by grid (``ch_begin[g]:ch_end[g]``),
    ``species_match[n_ch, N_SPECIES_SLOTS]`` (int8), ``total_size`` the columns of ``acc``,
    ``n_residual`` the per-history residual columns, ``bounds`` the quantities used for the
    quanta (recorded), ``memory_bytes`` the accumulator memory of the whole run (edep included)."""

    channels: tuple[Channel, ...]
    quantities: tuple[QuantityDef, ...]
    species_match: NDArray[np.int8]
    ch_begin: tuple[int, ...]
    ch_end: tuple[int, ...]
    total_size: int
    n_residual: int
    lookups: tuple[LookupTable, ...]
    bounds: dict[str, float]
    memory_bytes: int
    histories_per_batch: int
    path_bound_mm: float = math.inf
    """``B_L``: the per-history path length every history must not exceed (checked at runtime)."""

    def channel_index(self, kind: str, grid: int) -> list[int]:
        """Indices of the channels of ``kind`` on grid ``grid``."""
        return [i for i, c in enumerate(self.channels) if c.kind == kind and c.grid == grid]

    def count_channel(self, grid: int, class_mask: int) -> int:
        """Index of the automatic piece-count channel (kind N) of ``grid`` for the class
        ``class_mask`` (``CLASS_STEP`` or ``CLASS_LOCAL``)."""
        return next(
            i for i, c in enumerate(self.channels)
            if c.kind == "N" and c.grid == grid and c.class_mask == class_mask
        )  # fmt: skip

    def piece_count_indices(self, ci: int) -> tuple[int, ...]:
        """Indices of the automatic piece-count channels (kind N) whose counts bound channel
        ``ci``: N (class step) if its class mask contains "step", N_local if it contains "local"
        (both for edep and dose). The single source of every rounding bound."""
        c = self.channels[ci]
        return tuple(
            self.count_channel(c.grid, m) for m in (CLASS_STEP, CLASS_LOCAL) if c.class_mask & m
        )

    def summary(self) -> dict[str, Any]:
        """JSON-serialisable record for the effective configuration (quanta included)."""
        return {
            "n_channels": len(self.channels),
            "total_size": self.total_size,
            "n_residual_columns": self.n_residual,
            "histories_per_batch": self.histories_per_batch,
            "memory_bytes": self.memory_bytes,
            "bounds": dict(self.bounds),
            "lookups": [lk.provenance() for lk in self.lookups],
            "channels": [
                {
                    "kind": c.kind,
                    "grid": c.grid,
                    "class_mask": c.class_mask,
                    "species": [i for i, v in enumerate(c.species) if v],
                    "generation": [c.gen_lo, c.gen_hi],
                    "lookup": c.lookup,
                    "spectrum": None
                    if c.spectrum is None
                    else {"edges": list(c.spectrum.edges), "log": c.spectrum.log},
                    "size": c.size,
                    "offset": c.offset,
                    "k": c.k,
                    "bound_per_history": c.bound_per_history,
                }
                for c in self.channels
            ],
            "quantities": [
                {
                    "name": q.name,
                    "kind": q.kind,
                    "numerator": q.numerator,
                    "denominator": q.denominator,
                    "units": q.units,
                    "definition": q.definition,
                }
                for q in self.quantities
            ],
        }


# ---------------------------------------------------------------------------------------------
# axis parameters shared by every backend (the same floating-point numbers everywhere)


def lookup_axis_params(lk: LookupTable) -> tuple[float, float, int, int]:
    """``(a0, inv_da, n, log)`` of a lookup axis for ``lookup_bin``: ``a0`` is the first axis point
    (its logarithm for a log axis) and ``inv_da`` the inverse spacing in the axis variable."""
    n = int(lk.axis_values.size)
    log = lk.axis_spacing == "log"
    a0 = float(np.log(lk.axis_values[0]) if log else lk.axis_values[0])
    a1 = float(np.log(lk.axis_values[-1]) if log else lk.axis_values[-1])
    return a0, (n - 1) / (a1 - a0), n, 1 if log else 0


def spectrum_axis_params(spec: SpectrumSpec) -> tuple[float, float, int, int]:
    """``(a0, inv_da, n_bins, log)`` of a spectrum for ``spectrum_bin``."""
    e = spec.edges
    a0 = math.log(e[0]) if spec.log else e[0]
    a1 = math.log(e[-1]) if spec.log else e[-1]
    return a0, spec.n_bins / (a1 - a0), spec.n_bins, int(spec.log)


# ---------------------------------------------------------------------------------------------
# quantum rule


def adaptive_exponent(histories_per_batch: int, bound_per_history: float) -> int:
    """``k = min(40, floor(62 - log2(hpb B_c)))`` such that ``hpb B_c 2^k < 2^62``.

    ``k`` may be negative for an enormous bound (the precision floor then fails closed)."""
    if not (bound_per_history > 0.0 and math.isfinite(bound_per_history)):
        raise fail(
            f"channel bound per history must be positive and finite, got {bound_per_history}"
        )
    total = float(histories_per_batch) * bound_per_history
    k = min(K_MAX, math.floor(CAPACITY_EXPONENT - math.log2(total)))
    while math.ldexp(total, k) >= math.ldexp(1.0, CAPACITY_EXPONENT):  # float edge cases
        k -= 1
    return int(k)


def fe_exponent(f_max: float) -> int:
    """Exponent of the FE quantum ``2^-30 2^ceil(log2 f_max)`` (tied to the energy quantum);
    ``f_max`` must be within ``[2^-60, 2^60]`` (checked by the compiler), so ``k`` is in
    ``[-30, 90]`` and the scale ``2^k`` is a finite float64."""
    return E_QUANTUM_EXPONENT - math.ceil(math.log2(f_max))


# ---------------------------------------------------------------------------------------------
# water stopping power helpers (MeV/mm, E in MeV total kinetic energy of the ion)


def _water_s_linear(water: StoppingTable, e_mev: float, a: int) -> float:
    return float(water.stopping_at(e_mev / a)) * water.material.density_g_cm3 / 10.0


def water_s_extrema(
    water: StoppingTable, e_lo_mev: float, e_hi_mev: float, a: int
) -> tuple[float, float, float]:
    """``(S_min, S_max, S_ref)`` of the water linear stopping power [MeV/mm] over
    ``[e_lo, e_hi]`` (nodes inside plus both ends; the log-log interpolant is monotone between
    nodes); ``S_ref = S_w(e_hi)``. Raises if the water table does not cover the interval."""
    lo_t, hi_t = float(water.energy_per_u[0]), float(water.energy_per_u[-1])
    if e_lo_mev / a < lo_t * (1 - 1e-12) or e_hi_mev / a > hi_t * (1 + 1e-12):
        raise fail(
            f"the water stopping table covers [{lo_t}, {hi_t}] MeV/u, not the required "
            f"[{e_lo_mev / a}, {e_hi_mev / a}] MeV/u"
        )
    nodes = water.energy_per_u
    inside = nodes[(nodes > e_lo_mev / a) & (nodes < e_hi_mev / a)]
    s_nodes = water.s_el_linear[(nodes > e_lo_mev / a) & (nodes < e_hi_mev / a)]
    s_lo = _water_s_linear(water, e_lo_mev, a)
    s_ref = _water_s_linear(water, e_hi_mev, a)
    cand = np.concatenate([[s_lo, s_ref], s_nodes if inside.size else []])
    return float(cand.min()), float(cand.max()), s_ref


RAMP_LOSS_RATIO_MAX = 2.0
"""Largest ``f_E = dE_mean / E_mid`` of a step: ``dE_mean <= E`` (the CSDA mean loss is a
difference ``E - E_1`` of the range inversion with ``E_1 >= 0``, at most the whole energy in the
final range step) and ``E_mid = E - dE_mean / 2 >= E / 2``, so ``f_E <= 1 / (1/2) = 2``."""


def water_ramp_envelope(
    water: StoppingTable,
    tables: TransportTables,
    e_lo_mev: float,
    e_hi_mev: float,
    a: int,
) -> dict[str, float]:
    """Envelope of the water stopping power *actually scored* by the channels (decision 0040).

    A step starting at ``E in (e_lo, e_hi]`` has the hook midpoint energy ``E_mid = E - dE/2`` with
    ``0 < dE <= E`` (so ``E_mid >= E/2 > e_lo/2``: a cutoff-crossing step reaches below the
    cutoff) and the ramp ``S(tau) = S_mid + k tau``, ``k = -gamma S_mid dE / (E_mid s_act)``,
    ``|tau| <= s_act/2``; hence ``S in S_mid (1 -+ amp(E_mid))`` with
    ``amp = |gamma(E_mid)| f_E / 2`` (evaluated per table bin, the extrema taken over bins),
    ``f_E <= RAMP_LOSS_RATIO_MAX = 2``. ``S_mid``, ``gamma`` and ``S_ref``
    are evaluated on the *runtime* water row (``tables.water_ln_s_mass``, the interpolant the
    scorer uses; the source table only enters through the coverage check and, without a runtime
    row, a conservative fallback), per bin over ``[e_lo/2, e_hi]`` (clamped as the runtime clamps).
    Returns ``S_mid`` extrema, ``gamma_max``, ``amp`` and the ramp extrema ``S_bar_min/max``;
    fails closed if ``amp >= 1`` (the ramp could reach zero or negative: never clamped)."""
    floor_mev = float(water.energy_per_u[0]) * a
    e_lo_mid = max(0.5 * e_lo_mev, floor_mev)
    # coverage of the source table (fails closed); every value below comes from the runtime row
    s_lo_mid, s_hi_mid, s_ref = water_s_extrema(water, min(e_lo_mid, e_hi_mev), e_hi_mev, a)
    row = tables.water_ln_s_mass
    if row is not None and row.size >= 2:
        inv = tables.water_inv_dln_e
        n_bins = row.size - 1
        gam = np.diff(row) * inv  # d ln S / d ln E of every bin
        edges = tables.water_ln_e0 + np.arange(row.size) / inv  # ln E of the grid points
        lo, hi = math.log(0.5 * e_lo_mev), math.log(e_hi_mev)
        pts = np.concatenate([[lo], edges[(edges > lo) & (edges < hi)], [hi]])
        rho_w = tables.water_density_g_cm3

        def s_at(ln_e: float) -> float:  # as the runtime: clamped, log-log linear in the bin
            t = min(max((ln_e - tables.water_ln_e0) * inv, 0.0), float(n_bins))
            i = min(int(math.floor(t)), n_bins - 1)
            return math.exp(row[i] + (t - i) * (row[i + 1] - row[i])) * rho_w / 10.0

        s_bar_lo, s_bar_hi, gamma, amp = math.inf, 0.0, 0.0, 0.0
        s_lo_mid, s_hi_mid = math.inf, 0.0
        s_ref = s_at(hi)  # runtime row at E_hi
        for x0, x1 in zip(pts[:-1], pts[1:], strict=True):
            ib = min(max(int(math.floor(((0.5 * (x0 + x1)) - tables.water_ln_e0) * inv)), 0),
                     n_bins - 1)  # fmt: skip
            g_i = abs(float(gam[ib]))
            a_i = g_i * RAMP_LOSS_RATIO_MAX / 2.0
            ends = (s_at(float(x0)), s_at(float(x1)))
            s_lo_mid, s_hi_mid = min(s_lo_mid, min(ends)), max(s_hi_mid, max(ends))
            s_bar_lo = min(s_bar_lo, min(ends) * (1.0 - a_i))
            s_bar_hi = max(s_bar_hi, max(ends) * (1.0 + a_i))
            gamma, amp = max(gamma, g_i), max(amp, a_i)
    else:  # no runtime row: the nodes of the water table (one global slope, conservative)
        ln_s = np.log(water.s_el_linear)
        ln_e = np.log(water.energy_per_u)
        gamma = float(np.max(np.abs(np.diff(ln_s) / np.diff(ln_e))))
        amp = gamma * RAMP_LOSS_RATIO_MAX / 2.0
        s_bar_lo, s_bar_hi = s_lo_mid * (1.0 - amp), s_hi_mid * (1.0 + amp)
    if not math.isfinite(amp) or amp >= 1.0 or not s_bar_lo > 0.0:
        raise fail(
            f"the water stopping-power ramp is not guaranteed positive: the largest log-log slope "
            f"|gamma| = {gamma:g} of the water table over [{e_lo_mid:g}, {e_hi_mev:g}] MeV gives "
            f"the ramp amplification {amp:g} >= 1 (S_bar = S_mid (1 -+ |gamma| f_E / 2), "
            f"f_E <= {RAMP_LOSS_RATIO_MAX:g}), so S_bar could be zero or negative; LET channels "
            "cannot be scored with this stopping-power table"
        )
    return {
        "s_mid_min": s_lo_mid,
        "s_mid_max": s_hi_mid,
        "s_ref": s_ref,
        "gamma_max": gamma,
        "amp": amp,
        "s_bar_min": s_bar_lo,
        "s_bar_max": s_bar_hi,
    }


def _stopping_ratio_max(
    tables: TransportTables,
    material: int,
    water: StoppingTable,
    rho_min: float,
    e_lo_mev: float,
    e_hi_mev: float,
    a: int,
    amp_fallback: float,
) -> float:
    """``max_E (1 + amp(E)) S_w / (rho_min S_m)`` (mass stopping powers, the water density factors
    in) over ``[e_lo, e_hi]``: the bound of ``int S_bar dl`` per MeV deposited in the least dense
    voxel. The compiler passes ``e_lo = E_cut/2``, the lowest reachable *midpoint* energy of a
    cutoff-crossing step (``E_mid >= E/2 > E_cut/2``), not ``E_cut``: the scorer evaluates ``S_w``
    there, and a custom material table may drop sharply below the cutoff while the ratio stays
    ordinary above it. ``amp(E)`` is the ramp amplification ``|gamma_w(E)|`` of the runtime water
    row bin (the larger of the two bins at a node; ``amp_fallback`` without a row), applied to
    ``S_w`` as in :func:`water_ramp_envelope`. Both stopping powers are taken from the *runtime*
    rows the scorer uses (``tables.water_ln_s_mass`` and the material row; never the
    source-table interpolant, which differs between runtime nodes after resampling) and are
    log-log piecewise linear
    between the union of the runtime material and water nodes, so the ratio is monotone between
    them and the maximum is attained at a node or an end of the interval; the ratio candidate of
    the LS bound is therefore *proven* over the complete reachable domain (it is not dropped)."""
    n = int(tables.n_e)
    ln_e = tables.ln_e0[material] + np.arange(n) / tables.inv_dln_e[material]
    row = tables.water_ln_s_mass
    have_row = row is not None and row.size >= 2
    nodes = [np.exp(ln_e)]
    if have_row:
        assert row is not None
        nodes.append(np.exp(tables.water_ln_e0 + np.arange(row.size) / tables.water_inv_dln_e))
    e_all = np.concatenate(nodes)
    pts = np.concatenate([[e_lo_mev, e_hi_mev], e_all[(e_all > e_lo_mev) & (e_all < e_hi_mev)]])
    ln_s_m = np.interp(np.log(pts), ln_e, tables.ln_s_mass[material])
    s_m_mass = np.exp(ln_s_m)
    amp = np.full(pts.shape, amp_fallback)
    if have_row:
        assert row is not None
        # the runtime water interpolant (log-log linear on the runtime grid), as the scorer
        ln_ew = tables.water_ln_e0 + np.arange(row.size) / tables.water_inv_dln_e
        s_w_mass = np.exp(np.interp(np.log(pts), ln_ew, row))
        rho_w = tables.water_density_g_cm3
        gam = np.abs(np.diff(row)) * tables.water_inv_dln_e
        t = (np.log(pts) - tables.water_ln_e0) * tables.water_inv_dln_e
        i_hi = np.clip(np.floor(t + 1e-9).astype(int), 0, gam.size - 1)
        i_lo = np.clip(np.ceil(t - 1e-9).astype(int) - 1, 0, gam.size - 1)
        amp = np.maximum(gam[i_hi], gam[i_lo]) * RAMP_LOSS_RATIO_MAX / 2.0
    else:  # no runtime row: the source table (conservative global amplification)
        s_w_mass = water.stopping_at(pts / a)
        rho_w = water.material.density_g_cm3
    # S_w,lin = s_w_mass * rho_w / 10 ; S_m,lin(least dense) = s_m_mass * rho_min / 10
    ratio = (1.0 + amp) * s_w_mass * rho_w / (s_m_mass * rho_min)
    return float(ratio.max())


def mixed_path_bound_mm(
    tables: TransportTables, rho_min_by_material: dict[int, float], e_hi_mev: float
) -> float:
    """``int_0^{E_hi} dE / min_m S_lin,m(E)`` [mm]: the CSDA path of a particle that follows, at
    every energy, the lowest linear stopping power ``rho_min,m S_mass,m / 10`` of any material
    present (``rho_min,m`` the lowest density among the voxels of material ``m``).

    A homogeneous range ``max_m R_m(E_hi) / rho_min,m`` does *not* bound a path that changes
    material as the energy falls: ``int dE / min_m S_m`` can exceed ``max_m int dE / S_m`` by an
    arbitrary factor. The runtime rows are log-log piecewise linear on the union of the material
    nodes. On an interval between adjacent union nodes every ``ln S_m`` is linear in ``ln E``, so
    ``ln min_m S_m`` is concave and lies above its chord through the node values; replacing the
    minimum by that power-law chord can only *raise* ``1/S`` and so the integral, which is then
    exact: ``E_a/S_a (r^(1-g) - 1)/(1 - g)`` (``ln r`` for ``g = 1``), ``r = E_b/E_a``. Below the
    lowest node every row is clamped (constant), as the runtime clamps. The result is a rigorous
    upper bound on the CSDA path (the straggling margin is applied by the caller)."""
    mats = sorted(rho_min_by_material)
    n = int(tables.n_e)
    nodes = [
        np.exp(tables.ln_e0[m] + np.arange(n) / tables.inv_dln_e[m]) for m in mats
    ]  # fmt: skip
    e_all = np.unique(np.concatenate(nodes))
    pts = np.concatenate([e_all[e_all < e_hi_mev], [e_hi_mev]])
    ln_pts = np.log(pts)
    s_lin = np.min(
        [
            np.exp(
                np.interp(
                    ln_pts,
                    tables.ln_e0[m] + np.arange(n) / tables.inv_dln_e[m],
                    tables.ln_s_mass[m],
                )
            )
            * rho_min_by_material[m]
            / 10.0
            for m in mats
        ],
        axis=0,
    )
    total = float(pts[0] / s_lin[0])  # [0, lowest node]: every row clamped, S constant
    for e_a, e_b, s_a, s_b in zip(pts[:-1], pts[1:], s_lin[:-1], s_lin[1:], strict=True):
        g = math.log(s_b / s_a) / math.log(e_b / e_a)
        ln_r = math.log(e_b / e_a)
        if abs(1.0 - g) < 1e-9:
            total += e_a / s_a * ln_r
        else:
            total += e_a / s_a * math.expm1((1.0 - g) * ln_r) / (1.0 - g)
    return total


# ---------------------------------------------------------------------------------------------
# request validation helpers


def _spectrum_spec(name: str, edges: tuple[float, ...]) -> SpectrumSpec:
    e = np.asarray(edges, dtype=np.float64)
    if e.size - 1 > MAX_SPECTRUM_BINS:
        raise fail(f"tally {name!r}: more than {MAX_SPECTRUM_BINS} spectrum bins")
    if e.size == 2:
        return SpectrumSpec(tuple(edges), False)
    d_lin = np.diff(e)
    d_log = np.diff(np.log(e))
    if np.all(np.abs(d_lin - d_lin.mean()) <= 1e-9 * d_lin.mean()):
        return SpectrumSpec(tuple(edges), False)
    if np.all(np.abs(d_log - d_log.mean()) <= 1e-9 * d_log.mean()):
        return SpectrumSpec(tuple(edges), True)
    raise fail(
        f"tally {name!r}: spectrum edges must be uniform in energy or in ln energy "
        "(non-uniform edges are rejected)"
    )


def _species_row(species: tuple[str, ...] | None, *, include_pseudo: bool) -> tuple[int, ...]:
    row = [0] * N_SPECIES_SLOTS
    if species is None:
        for i in transported_ids():
            row[i] = 1
        if include_pseudo:
            row[species_by_name(NUCLEAR_LOCAL).id] = 1
    else:
        for name in species:
            row[species_by_name(name).id] = 1
    return tuple(row)


def _generation_range(generation: str) -> tuple[int, int]:
    return {"all": (0, GEN_MAX), "primary": (0, 0), "secondary": (1, GEN_MAX)}[generation]


def _check_species_generation(
    req: TallyRequest, producible: frozenset[tuple[str, str]]
) -> list[str]:
    """Fail closed for unknown or unproducible species and generations; return the species names
    the request selects (``None`` means every producible species)."""
    if req.species is not None:
        for name in req.species:
            try:
                sp = species_by_name(name)
            except ValueError as exc:
                raise fail(f"tally {req.name!r}: {exc}") from exc
            if not sp.transported and req.quantity not in ("edep", "dose"):
                raise fail(
                    f"tally {req.name!r}: pseudo-species {name!r} is not transported and "
                    "carries no LET, fluence or lookup contribution"
                )
    prod_names = sorted({s for s, _ in producible})
    selected = list(req.species) if req.species is not None else prod_names
    for name in selected:
        gens = {g for s, g in producible if s == name}
        wanted = {"all": gens, "primary": {"primary"}, "secondary": {"secondary"}}[req.generation]
        if not (gens & wanted):
            # also the pseudo-species (nuclear_local): an explicit request is accepted only once
            # the engine advertises it with its generation, never as an all-zero channel
            raise fail(
                f"tally {req.name!r}: species {name!r} (generation {req.generation!r}) is not "
                f"producible by this engine (producible: {sorted(producible)}; there is no "
                "secondary transport until V3-005A); nothing silently scores zero"
            )
    return selected


def _check_lookup(
    req: TallyRequest,
    table: LookupTable,
    selected: list[str],
    producible: frozenset[tuple[str, str]],
    e_lo_per_u: float,
    e_hi_per_u: float,
    s_lo: float,
    s_hi: float,
) -> None:
    need = [s for s in selected if s in {p for p, _ in producible}]
    missing = [s for s in need if s not in table.values]
    if missing:
        raise fail(
            f"tally {req.name!r}: lookup table {table.name!r} has no values for the producible "
            f"species {missing}"
        )
    lo, hi = (e_lo_per_u, e_hi_per_u) if table.axis == "energy_per_nucleon_mev" else (s_lo, s_hi)
    if not table.covers(lo, hi):
        a, b = table.axis_range
        raise fail(
            f"tally {req.name!r}: lookup table {table.name!r} covers [{a}, {b}] "
            f"({table.axis}) but the transport reaches [{lo}, {hi}]"
        )


# ---------------------------------------------------------------------------------------------
# the compiler


@dataclass
class _Spec:
    kind: str
    grid: int
    class_mask: int
    species: tuple[int, ...]
    gen: tuple[int, int]
    lookup: int
    spectrum: SpectrumSpec | None
    f_max: float = 0.0  # FE only: largest table value over the selected species (not in the key)

    def key(self) -> tuple[Any, ...]:
        sp = None if self.spectrum is None else (self.spectrum.edges, self.spectrum.log)
        return (self.kind, self.grid, self.class_mask, self.species, self.gen, self.lookup, sp)


def compile_channels(
    requests: Sequence[TallyRequest],
    lookups: Sequence[LookupTable],
    grids: Sequence[ScoringGrid],
    *,
    geometry: VoxelGeometry,
    tables: TransportTables,
    water: StoppingTable,
    projectile: Projectile,
    e_cut_mev: float,
    e_hi_mev: float,
    n_histories: int,
    n_batches: int,
    cpu_workers: int,
    memory_budget_bytes: int,
    max_steps: int,
    scoring_pieces: int,
    producible: frozenset[tuple[str, str]],
    max_step_mm: float | None = None,
) -> ChannelPlan:
    """Compile ``requests`` into a :class:`ChannelPlan` (deduplicated, ordered by grid, with
    offsets, quanta and memory guard) or raise ``UnsupportedCombinationError``.

    ``e_hi_mev`` is the largest source energy (mean plus 6 sigma, capped at the table maximum).
    ``producible`` is the engine capability (``ionmc.species.producible``). The bounds use the
    water stopping power over ``[E_cut, E_hi]`` as decision 0040 section 4 prescribes.
    """
    a = projectile.a
    species_of_projectile(projectile)  # registry membership
    grid_index = {g.name: i for i, g in enumerate(grids)}
    lookup_index = {lk.name: i for i, lk in enumerate(lookups)}
    if len(lookup_index) != len(lookups):
        raise fail("lookup table names must be unique")
    names = [r.name for r in requests]
    reserved = {EXCLUDED_CHANNEL_NAME, PIECE_COUNT_NAME, LOCAL_PIECE_COUNT_NAME} & set(names)
    if reserved:
        raise fail(f"tally request names {sorted(reserved)} are reserved for automatic channels")
    if len(set(names)) != len(names):
        raise fail(f"tally request names must be unique, got {names}")

    hpb = n_histories // n_batches
    env = water_ramp_envelope(water, tables, e_cut_mev, e_hi_mev, a)
    s_min, s_max, s_ref = env["s_bar_min"], env["s_bar_max"], env["s_ref"]
    e_floor_per_u = float(tables.e_min_mev.max()) / a
    env_floor = water_ramp_envelope(water, tables, float(tables.e_min_mev.max()), e_hi_mev, a)
    s_floor_min, s_floor_max = env_floor["s_bar_min"], env_floor["s_bar_max"]

    specs: dict[tuple[Any, ...], _Spec] = {}
    qdefs: list[tuple[TallyRequest, str, tuple[Any, ...], tuple[Any, ...] | None, str, str]] = []
    used_grids: set[int] = set()
    f_max = 0.0  # largest over all FE channels (recorded only)

    def get(spec: _Spec) -> tuple[Any, ...]:
        specs.setdefault(spec.key(), spec)
        return spec.key()

    for req in requests:
        if req.grid not in grid_index:
            raise fail(
                f"tally {req.name!r}: unknown scoring grid {req.grid!r} ({list(grid_index)})"
            )
        if req.let_medium != "water":
            raise fail(
                f"tally {req.name!r}: let_medium {req.let_medium!r} is not supported; the LET "
                "medium is water (decision 0040), a local-medium LET is a later channel kind"
            )
        if req.dose_reference != "medium":
            raise fail(
                f"tally {req.name!r}: dose_reference {req.dose_reference!r} is not supported; "
                "dose is dose-to-medium, dose-to-water needs a dedicated channel kind"
            )
        selected = _check_species_generation(req, producible)
        gi = grid_index[req.grid]
        used_grids.add(gi)
        gen = _generation_range(req.generation)
        li = -1
        f_req = 0.0
        spectrum = None
        if req.lookup is not None:
            if req.lookup not in lookup_index:
                raise fail(f"tally {req.name!r}: unknown lookup table {req.lookup!r}")
            li = lookup_index[req.lookup]
            table = lookups[li]
            _check_lookup(
                req, table, selected, producible,
                e_floor_per_u, e_hi_mev / a, s_floor_min, s_floor_max,
            )  # fmt: skip
            # per-channel domain (review 67f03e06 b): only the species this request selects
            need_sp = [s for s in selected if s in {p for p, _ in producible}]
            f_req = max(float(table.values[s].max()) for s in need_sp)
            if f_req > 0.0 and not (
                math.ldexp(1.0, FE_F_MIN_EXPONENT) <= f_req <= math.ldexp(1.0, FE_F_MAX_EXPONENT)
            ):
                raise fail(
                    f"tally {req.name!r}: the largest value {f_req:g} of lookup table "
                    f"{table.name!r} is outside the admissible magnitude range "
                    f"[2^{FE_F_MIN_EXPONENT}, 2^{FE_F_MAX_EXPONENT}] (the fixed-point quantum and "
                    "scale of the lookup channel must be representable); rescale the table "
                    "and its units"
                )
            f_max = max(f_max, f_req)
        if req.energy_edges_mev_per_u is not None:
            spectrum = _spectrum_spec(req.name, req.energy_edges_mev_per_u)

        def spec(
            kind: str,
            cls: int,
            *,
            li: int = -1,
            sp: SpectrumSpec | None = None,
            _f: float = 0.0,
            _req: TallyRequest = req,
            _gi: int = gi,
            _gen: tuple[int, int] = gen,
        ) -> tuple[Any, ...]:
            pseudo = kind == "E" and bool(cls & CLASS_LOCAL) and _req.species is None
            return get(
                _Spec(
                    kind,
                    _gi,
                    cls,
                    _species_row(_req.species, include_pseudo=pseudo),
                    _gen,
                    li,
                    sp,
                    _f,
                )
            )

        q = req.quantity
        both = CLASS_STEP | CLASS_LOCAL
        if q in ("edep", "dose"):
            qdefs.append((req, "linear", spec("E", both), None, "MeV", f"sum of deposits ({q})"))
        elif q == "fluence":
            qdefs.append((req, "linear", spec("L", CLASS_STEP), None, "mm", "sum of path lengths"))
        elif q == "let_t":
            qdefs.append(
                (req, "ratio", spec("LS", CLASS_STEP), spec("L", CLASS_STEP), "MeV/mm",
                 "sum(l S)/sum(l), S = unrestricted electronic stopping power in water")
            )  # fmt: skip
        elif q == "let_d":
            qdefs.append(
                (req, "ratio", spec("LS2", CLASS_STEP), spec("LS", CLASS_STEP), "MeV/mm",
                 "sum(l S^2)/sum(l S), S in water")
            )  # fmt: skip
        elif q == "let_d_eps":
            qdefs.append(
                (req, "ratio", spec("ES", CLASS_STEP), spec("E", CLASS_STEP), "MeV/mm",
                 "sum(eps S)/sum(eps) over step deposits, S in water")
            )  # fmt: skip
        elif q == "lookup_sum":
            qdefs.append(
                (req, "linear", spec("FE", CLASS_STEP, li=li, _f=f_req), None,
                 f"MeV*[{lookups[li].units}]", "sum(eps f)")
            )  # fmt: skip
        elif q == "lookup_dose_avg":
            qdefs.append(
                (req, "ratio", spec("FE", CLASS_STEP, li=li, _f=f_req), spec("E", CLASS_STEP),
                 lookups[li].units, "sum(eps f)/sum(eps) over step deposits")
            )  # fmt: skip
        else:  # fluence_spectrum
            qdefs.append(
                (req, "linear", spec("FL", CLASS_STEP, sp=spectrum), None, "mm",
                 "sum of path lengths per energy bin per nucleon (under/overflow included)")
            )  # fmt: skip

    # automatic channels of every used grid
    auto_excluded: dict[int, tuple[Any, ...]] = {}
    for gi in sorted(used_grids):
        get(
            _Spec(
                "N",
                gi,
                CLASS_STEP,
                _species_row(None, include_pseudo=False),
                (0, GEN_MAX),
                -1,
                None,
            )
        )
        get(
            _Spec(
                "N",
                gi,
                CLASS_LOCAL,
                _species_row(None, include_pseudo=True),
                (0, GEN_MAX),
                -1,
                None,
            )
        )
        auto_excluded[gi] = get(
            _Spec(
                "E",
                gi,
                CLASS_LOCAL,
                _species_row(None, include_pseudo=True),
                (0, GEN_MAX),
                -1,
                None,
            )
        )
    if any(s.kind == "FE" and s.f_max <= 0.0 for s in specs.values()):
        raise fail("a lookup table with no positive value cannot be used (all values are zero)")

    ordered = sorted(specs.values(), key=lambda s: s.grid)  # stable: creation order inside a grid
    if len(ordered) > MAX_CHANNELS:
        raise fail(f"{len(ordered)} channels exceed the limit of {MAX_CHANNELS}")

    # --- bounds ---------------------------------------------------------------------------
    dens = geometry.densities_g_cm3()
    r_max = 0.0
    rho_min_of: dict[int, float] = {}
    for m in range(len(geometry.materials)):
        mask = geometry.material_index == m
        if not mask.any():
            continue
        rho_min = float(dens[mask].min())
        rho_min_of[m] = rho_min
        r_max = max(
            r_max,
            _stopping_ratio_max(
                tables, m, water, rho_min, 0.5 * e_cut_mev, e_hi_mev, a, env["amp"]
            ),
        )
    # Per-history path bound B_L (a runtime-CHECKED assumption, decision 0040 section 4): the
    # heterogeneous CSDA path with a 25 % straggling margin. Neither straggling sampler
    # guarantees that a history stops within any multiple of its CSDA path (the legacy sampler
    # can clamp losses to zero, Gamma losses can be very small), so the bound is not provable
    # from the tables. It is enforced instead: the engine accumulates every history's scored path
    # (sum of s_act, the quantity of the L channel) and counts a history whose path exceeds B_L in
    # the per-history tally column ``path_bound_exceeded``; a nonzero count invalidates the
    # result. Hence either every history satisfied L_h <= B_L (then L, FL, LS = sum l S_bar <=
    # S_bar_max B_L and LS2 <= S_bar_max^2 B_L hold for every history and the int64 capacity proof
    # holds exactly as compiled) or the result is invalid and no accumulator is ever used, so a
    # wrap-around cannot produce a silently wrong result. ES <= S_bar_max E_hi needs no
    # assumption (sum of eps <= E_hi by energy conservation). The truncation bound
    # max_steps * max_step_mm (a step is at most max_step_mm long, a history at most max_steps
    # steps) is rigorous but about 20 times coarser; it is recorded only.
    b_l = STRAGGLING_MARGIN * mixed_path_bound_mm(tables, rho_min_of, e_hi_mev)
    b_l_trunc = (
        None if max_step_mm is None else float(max_steps) * float(max_step_mm) * STEP_LENGTH_MARGIN
    )
    b_ls = s_max * b_l  # S_bar <= S_bar_max on every piece
    bound = {
        "E": e_hi_mev,
        "L": b_l,
        "FL": b_l,
        "LS": b_ls,
        "LS2": s_max * b_ls,  # S_bar^2 <= S_bar_max^2 on every piece
        "ES": s_max * e_hi_mev,  # sum(eps S) <= S_bar_max sum(eps), sum(eps) <= E_hi (conservation)
        "N": float(max_steps) * 2.0 * scoring_pieces,
    }
    if hpb * e_hi_mev / QUANTUM_MEV >= MAX_QUANTA:
        raise fail(
            f"{hpb} histories per batch of up to {e_hi_mev:g} MeV exceed the capacity of a "
            "fixed-point accumulator; increase n_batches or reduce n_histories"
        )
    if hpb * bound["N"] >= 2.0**CAPACITY_EXPONENT:
        raise fail("a piece-count channel (N, N_local) could overflow its accumulator; reduce hpb")

    k_of: dict[str, int] = {"E": E_QUANTUM_EXPONENT, "N": 0}
    fe_of: dict[tuple[Any, ...], tuple[int, float]] = {}
    present = {s.kind for s in ordered}
    for kind in ("L", "LS", "LS2", "ES"):
        if kind in present or (kind == "L" and "FL" in present):
            k_of[kind] = adaptive_exponent(hpb, bound[kind])
    if "FL" in present:  # spectra share the exponent of L (sum over bins = L bitwise)
        k_of["FL"] = k_of["L"]

    # FE: bound, exponent and floor per channel from the table values of its own species
    for s in ordered:
        if s.kind == "FE":
            s_bound = s.f_max * e_hi_mev
            s_k = fe_exponent(s.f_max)
            fe_of[s.key()] = (s_k, s_bound)

    # --- precision floor ------------------------------------------------------------------
    checks: set[tuple[str, int, float]] = set()
    for s in ordered:
        if s.kind == "N":
            continue
        u = 1.0 * s_ref ** _FLOOR_POWER[s.kind]
        if s.kind == "FE":
            u *= s.f_max
            checks.add(("FE", fe_of[s.key()][0], u))
        else:
            checks.add((s.kind, k_of[s.kind], u))
    for kind, k, u in sorted(checks):
        if math.ldexp(1.0, -k) > math.ldexp(u, FLOOR_EXPONENT):
            raise fail(
                f"channel kind {kind}: the quantum 2^-{k} is above the precision floor "
                f"2^-16 u (u = {u:g}); too many histories per batch ({hpb}) for this "
                "configuration: increase n_batches or reduce n_histories"
            )

    # --- offsets, sizes, residual columns ---------------------------------------------------
    channels: list[Channel] = []
    key_to_index: dict[tuple[Any, ...], int] = {}
    offset = 0
    res_col = 0
    for s in ordered:
        n_vox = grids[s.grid].n_voxels
        size = n_vox * (s.spectrum.n_bins + 2 if s.spectrum is not None else 1)
        if s.kind == "N":
            rc = -1
        else:
            rc = res_col
            res_col += 1
        key_to_index[s.key()] = len(channels)
        if s.kind == "FE":
            k_c, b_c = fe_of[s.key()]
        else:
            k_c, b_c = k_of[s.kind], bound[s.kind]
        channels.append(
            Channel(
                s.kind, s.grid, s.class_mask, s.species, s.gen[0], s.gen[1], s.lookup, s.spectrum,
                size, offset, k_c, b_c, rc,
            )
        )  # fmt: skip
        offset += size

    total_voxels = sum(g.n_voxels for g in grids)
    memory = n_batches * (total_voxels + offset) * 8 * cpu_workers
    if memory > memory_budget_bytes:
        raise fail(
            f"per-batch accumulators including {len(channels)} scoring channels need {memory} "
            f"bytes (one private copy per worker process), above the memory budget of "
            f"{memory_budget_bytes} bytes (reduce n_batches, the grids, the spectra or workers)"
        )

    ch_begin: list[int] = []
    ch_end: list[int] = []
    for gi in range(len(grids)):
        idx = [i for i, c in enumerate(channels) if c.grid == gi]
        ch_begin.append(idx[0] if idx else 0)
        ch_end.append(idx[-1] + 1 if idx else 0)
    match = np.array([c.species for c in channels], dtype=np.int8).reshape(len(channels), -1)
    quantities = tuple(
        QuantityDef(
            r.name,
            r,
            kind,
            key_to_index[num],
            None if den is None else key_to_index[den],
            units,
            definition,
        )
        for r, kind, num, den, units, definition in qdefs
    )
    bounds = {
        "gamma_max": env["gamma_max"],
        "ramp_amplification": env["amp"],
        "S_w_min_mev_per_mm": s_min,
        "S_w_max_mev_per_mm": s_max,
        "S_w_ref_mev_per_mm": s_ref,
        "B_L_mm": b_l,
        "B_L_truncation_mm": math.nan if b_l_trunc is None else b_l_trunc,  # informative only
        "B_LS_mev": b_ls,
        "B_LS2": bound["LS2"],
        "B_ES": bound["ES"],
        "r_max": r_max,
        "f_max": f_max,
        "E_hi_mev": e_hi_mev,
    }
    return ChannelPlan(
        channels=tuple(channels),
        quantities=quantities,
        species_match=freeze_array(match, np.int8, "species_match"),
        ch_begin=tuple(ch_begin),
        ch_end=tuple(ch_end),
        total_size=offset,
        n_residual=res_col,
        lookups=tuple(lookups),
        bounds=bounds,
        memory_bytes=memory,
        histories_per_batch=hpb,
        path_bound_mm=b_l,
    )
