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
N      1                           piece count (automatic, bounds the fixed-point error)
=====  ==========================  ===========================================================

Ratios are ``LET_t = LS/L``, ``LET_d = LS2/LS``, ``LET_d^eps = ES/E_step`` and
``lookup_dose_avg = FE/E_step``; their reduction (delta method) is ``reduce_ratio`` (later step).
The automatic channels per used grid are ``N`` and ``edep_excluded_from_let`` (E, class "local",
all species including the pseudo-species).

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
E_QUANTUM_EXPONENT = 30  # QUANTUM_MEV = 2**-30
EXCLUDED_CHANNEL_NAME = "edep_excluded_from_let"
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

    def channel_index(self, kind: str, grid: int) -> list[int]:
        """Indices of the channels of ``kind`` on grid ``grid``."""
        return [i for i, c in enumerate(self.channels) if c.kind == kind and c.grid == grid]

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
    """Exponent of the FE quantum ``2^-30 2^ceil(log2 f_max)`` (tied to the energy quantum)."""
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


def _stopping_ratio_max(
    tables: TransportTables,
    material: int,
    water: StoppingTable,
    rho_min: float,
    e_lo_mev: float,
    e_hi_mev: float,
    a: int,
) -> float:
    """``max_E S_w / (rho_min S_m)`` (mass stopping powers, water density 1 g/cm3 factors in) over
    ``[e_lo, e_hi]``: the bound of ``int S_w dl`` per MeV deposited in the least dense voxel."""
    n = int(tables.n_e)
    ln_e = tables.ln_e0[material] + np.arange(n) / tables.inv_dln_e[material]
    e = np.exp(ln_e)
    pts = np.concatenate([[e_lo_mev, e_hi_mev], e[(e > e_lo_mev) & (e < e_hi_mev)]])
    ln_s_m = np.interp(np.log(pts), ln_e, tables.ln_s_mass[material])
    s_m_mass = np.exp(ln_s_m)
    s_w_mass = water.stopping_at(pts / a)
    # S_w,lin = s_w_mass * rho_w / 10 ; S_m,lin(least dense) = s_m_mass * rho_min / 10
    ratio = s_w_mass * water.material.density_g_cm3 / (s_m_mass * rho_min)
    return float(ratio.max())


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
            sp_obj = species_by_name(name)
            if sp_obj.transported or req.quantity not in ("edep", "dose"):
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
    if len(set(names)) != len(names):
        raise fail(f"tally request names must be unique, got {names}")

    hpb = n_histories // n_batches
    s_min, s_max, s_ref = water_s_extrema(water, e_cut_mev, e_hi_mev, a)
    e_floor_per_u = float(tables.e_min_mev.max()) / a
    s_floor_min, s_floor_max, _ = water_s_extrema(water, float(tables.e_min_mev.max()), e_hi_mev, a)

    specs: dict[tuple[Any, ...], _Spec] = {}
    qdefs: list[tuple[TallyRequest, str, tuple[Any, ...], tuple[Any, ...] | None, str, str]] = []
    used_grids: set[int] = set()
    f_max = 0.0

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
            f_max = max(f_max, *(float(v.max()) for v in table.values.values()))
        if req.energy_edges_mev_per_u is not None:
            spectrum = _spectrum_spec(req.name, req.energy_edges_mev_per_u)

        def spec(
            kind: str,
            cls: int,
            *,
            li: int = -1,
            sp: SpectrumSpec | None = None,
            _req: TallyRequest = req,
            _gi: int = gi,
            _gen: tuple[int, int] = gen,
        ) -> tuple[Any, ...]:
            pseudo = kind == "E" and bool(cls & CLASS_LOCAL) and _req.species is None
            return get(
                _Spec(
                    kind, _gi, cls, _species_row(_req.species, include_pseudo=pseudo), _gen, li, sp
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
                (req, "linear", spec("FE", CLASS_STEP, li=li), None,
                 f"MeV*[{lookups[li].units}]", "sum(eps f)")
            )  # fmt: skip
        elif q == "lookup_dose_avg":
            qdefs.append(
                (req, "ratio", spec("FE", CLASS_STEP, li=li), spec("E", CLASS_STEP),
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
    if f_max == 0.0 and any(s.kind == "FE" for s in specs.values()):
        raise fail("a lookup table with no positive value cannot be used (all values are zero)")

    ordered = sorted(specs.values(), key=lambda s: s.grid)  # stable: creation order inside a grid
    if len(ordered) > MAX_CHANNELS:
        raise fail(f"{len(ordered)} channels exceed the limit of {MAX_CHANNELS}")

    # --- bounds ---------------------------------------------------------------------------
    dens = geometry.densities_g_cm3()
    b_l = 0.0
    r_max = 0.0
    for m in range(len(geometry.materials)):
        mask = geometry.material_index == m
        if not mask.any():
            continue
        rho_min = float(dens[mask].min())
        b_l = max(b_l, STRAGGLING_MARGIN * tables.range_g_cm2(m, e_hi_mev) * 10.0 / rho_min)
        r_max = max(r_max, _stopping_ratio_max(tables, m, water, rho_min, e_cut_mev, e_hi_mev, a))
    b_ls = min(s_max * b_l, STRAGGLING_MARGIN * r_max * e_hi_mev)
    bound = {
        "E": e_hi_mev,
        "FE": f_max * e_hi_mev,
        "L": b_l,
        "FL": b_l,
        "LS": b_ls,
        "LS2": s_max * b_ls,
        "ES": s_max * e_hi_mev,
        "N": float(max_steps) * 2.0 * scoring_pieces,
    }
    if hpb * e_hi_mev / QUANTUM_MEV >= MAX_QUANTA:
        raise fail(
            f"{hpb} histories per batch of up to {e_hi_mev:g} MeV exceed the capacity of a "
            "fixed-point accumulator; increase n_batches or reduce n_histories"
        )
    if hpb * bound["N"] >= 2.0**CAPACITY_EXPONENT:
        raise fail("the piece-count channel N could overflow its accumulator; reduce hpb")

    k_of: dict[str, int] = {"E": E_QUANTUM_EXPONENT, "N": 0}
    if f_max > 0.0:
        k_of["FE"] = fe_exponent(f_max)
    present = {s.kind for s in ordered}
    for kind in ("L", "LS", "LS2", "ES"):
        if kind in present or (kind == "L" and "FL" in present):
            k_of[kind] = adaptive_exponent(hpb, bound[kind])
    if "FL" in present:  # spectra share the exponent of L (sum over bins = L bitwise)
        k_of["FL"] = k_of["L"]

    # --- precision floor ------------------------------------------------------------------
    for kind in sorted({s.kind for s in ordered}):
        if kind == "N":
            continue
        u = 1.0 * s_ref ** _FLOOR_POWER[kind]
        if kind == "FE":
            u *= f_max
        k = k_of[kind]
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
        channels.append(
            Channel(
                s.kind, s.grid, s.class_mask, s.species, s.gen[0], s.gen[1], s.lookup, s.spectrum,
                size, offset, k_of[s.kind], bound[s.kind], rc,
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
        "S_w_min_mev_per_mm": s_min,
        "S_w_max_mev_per_mm": s_max,
        "S_w_ref_mev_per_mm": s_ref,
        "B_L_mm": b_l,
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
    )
