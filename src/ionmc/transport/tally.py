"""Partition-independent reduction of per-history tallies and diagnostics (V3-003B).

Every backend writes, for each history, one row of independently accumulated float64 tallies
and one row of int32 transport-limit counters (``COUNTER_NAMES``). The tally row has
``N_FIXED_TALLIES + 2 n_grids + n_channel_columns`` columns, in this order:

* the fixed tallies (``TALLY_NAMES``);
* the energy deposited outside each scoring grid (``n_grids`` columns);
* the quantization residual of each grid (``n_grids`` columns);
* with scoring channels (``n_channel_columns = C_res + 2``, 0 without channels; decision 0040,
  ``ionmc.transport.channels``): ``C_res`` quantization-residual columns, one per channel but N;
  then the column ``C_res`` with the number of lookup arguments outside their table axis
  (``lookup_out_of_domain``); then the column ``C_res + 1`` with ``path_bound_exceeded``, 1.0 for a
  history whose scored path length exceeded the per-history bound ``B_L`` on which the channel
  quanta rest, else 0.0 (summed over histories it is ``ChannelRaw.path_bound_exceeded``; a nonzero
  value invalidates the result).

Column numbers of the channel block count from its first column. ``rows_to_partial`` stores the
columns unchanged and ``merge_partials`` indexes the last two from the end
(``n_channel_columns - 2`` and ``n_channel_columns - 1``). A history's row depends only on that
history, so it does not depend on how histories are split into chunks, worker processes or
threads. The host reduces the rows *exactly*: the sum of the float64 values of a column is
held as a Shewchuk-style expansion (a list of floats whose exact sum is the exact sum of the
column), expansions of different chunks and workers are concatenated, and the result is the
correctly rounded exact sum (``math.fsum``). The reduced tallies are therefore bit-identical
for any partition of the histories (test T13); counters are integer sums.

The per-voxel energy-deposit grids are int64 fixed-point accumulators (quantum ``QUANTUM_MEV`` =
2**-30 MeV): each deposit piece (see the track-length scoring in ``ionmc.transport.reference``)
is rounded to the nearest quantum by a deterministic function of the piece and added with an
integer (associative) addition, so the grids are bit-identical for any partition of the
histories and across chunk sizes and worker processes on one device kind (CPU against CUDA
is compared statistically, T12, never bit for bit). The rounding
error is at most q/2 per piece (random walk q/2 sqrt(N) over N pieces in a voxel) and is
tallied per history and grid (``quantization`` columns) so that the energy balance closes.
Conversion to float64 and the sum over workers happen at the reduction.

The scoring channels (decision 0040) are int64 accumulators ``acc[B, sum_c size_c]`` with the same
associativity argument; they are summed over chunks and workers as exact integers.

Diagnostic arrays are per history (end state) and per step (trace); they are concatenated in
history order.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray

END_CUTOFF = 0
END_ESCAPED = 1
END_TRUNCATED = 2
END_SOURCE_REJECTED = 3
END_MISSED_WORLD = 4

TRACE_COLUMNS = (
    "history",
    "step",
    "ix",
    "iy",
    "iz",
    "reason",
    "blocks",
    "attempts",
    "x_mm",
    "y_mm",
    "z_mm",
    "ux",
    "uy",
    "uz",
    "energy_mev",
    "deposit_mev",
    "step_mm",
)
"""Columns of the per-step trace (state after the step; ``blocks`` is the number of Philox
blocks drawn by the history so far, ``reason`` the step-limit reason of the shared
``select_step`` (0 geometry, 1 energy loss, 2 range, 3 maximum step)). The first eight columns
are discrete (``TRACE_N_DISCRETE``), the last nine continuous."""
TRACE_N_DISCRETE = 8
TRACE_N_CONTINUOUS = len(TRACE_COLUMNS) - TRACE_N_DISCRETE

COUNTER_NAMES = (
    "step_truncation",
    "stall",
    "straggling_rejection",
    "genealogy_overflow",
    "queue_overflow",
    "source_energy_out_of_range",
    "energy_inversion",
    "accumulator_overflow",
    "scoring_pieces_overflow",
)
TALLY_NAMES = ("initial", "cutoff", "step_deposit", "escaped", "truncated", "unaccounted")
N_FIXED_TALLIES = len(TALLY_NAMES)

QUANTUM_MEV = 2.0**-30
"""Fixed-point quantum of the deposit grids [MeV]: every deposit piece is rounded to the nearest
multiple (``floor(x / q + 1/2)``, a deterministic function of the piece) and accumulated in int64,
so the grids are bit-identical for any partition of the histories."""
QUANTUM_SCALE = 2.0**30
MAX_QUANTA = 2**62
"""Capacity bound of one voxel accumulator in quanta (validated before, checked after a run)."""


@dataclass
class ChannelRaw:
    """Reduced scoring channels of a run (decision 0040): ``acc`` the int64 accumulators
    ``[B, sum_c size_c]`` in the quanta of the channels, ``residual[r]`` the exact total of the
    quantization residual column ``r`` (``sum(x - n q_c)`` over all pieces, float64) and
    ``lookup_out_of_domain`` the number of lookup arguments outside their table axis (a nonzero
    value invalidates the result) and ``path_bound_exceeded`` the number of histories whose scored
    path length exceeded the per-history bound ``B_L`` on which the quanta rest (a nonzero value
    invalidates the result: the capacity proof is conditional on this check, decision 0040)."""

    acc: NDArray[np.int64]
    residual: list[float]
    lookup_out_of_domain: int
    path_bound_exceeded: int = 0


@dataclass
class RawTransport:
    """Raw reduced output of a transport run (float64).

    ``meta`` holds backend facts (device, compile time, chunking, worker reports); it is not
    part of the physics result. ``channels`` is None without scoring channels.
    """

    edep_mev: list[NDArray[np.float64]]
    tallies: dict[str, float]
    outside_mev: list[float]
    counters: dict[str, int]
    quantization_mev: list[float] = field(default_factory=list)
    diagnostics: dict[str, Any] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)
    channels: ChannelRaw | None = None


@dataclass
class HistoryDiagnostics:
    """End state of every history of a range and the per-step trace rows of the traced ones.

    ``trace_int`` (rows, 8) and ``trace_float`` (rows, 9) hold the trace rows of this range in
    history/step order (columns of ``TRACE_COLUMNS``). ``control_residual`` is the maximum over a
    history's steps of the snap displacement relative to the step and ``control_displacement`` the
    vector sum of its snap displacements [mm] (both zero unless ``truncated_hinge_diagnostic``).
    """

    end_position_mm: NDArray[np.float64]
    end_direction: NDArray[np.float64]
    end_energy_mev: NDArray[np.float64]
    end_code: NDArray[np.int8]
    control_residual: NDArray[np.float64]
    control_displacement: NDArray[np.float64]
    trace_int: NDArray[np.int32]
    trace_float: NDArray[np.float64]


@dataclass
class PartialTransport:
    """Result of the histories ``[h0, h1)`` (one chunk sequence of one worker).

    ``tally_components[c]`` is the exact-sum expansion of tally column ``c`` (``N_FIXED_TALLIES``
    fixed columns, then one per scoring grid for the outside deposit, then one per grid for the
    quantization residual ``sum(piece - quanta * q)`` of the deposits inside the grid);
    ``counter_sums`` the int64 sums of the counters; ``edep`` one ``(B, n_voxels)`` int64 array
    (quanta of ``QUANTUM_MEV``) per grid.
    """

    h0: int
    h1: int
    tally_components: list[list[float]]
    counter_sums: list[int]
    edep: list[NDArray[np.int64]]
    diagnostics: HistoryDiagnostics | None
    meta: dict[str, Any] = field(default_factory=dict)
    channel_acc: NDArray[np.int64] | None = None


def exact_components(values: NDArray[np.float64]) -> list[float]:
    """Floats whose exact sum equals the exact sum of ``values`` (all must be finite).

    Repeated correctly rounded sums of the residual (``math.fsum``) give a non-overlapping
    expansion; it ends when the residual rounds to zero, which for a sum of doubles means that
    it is exactly zero.
    """
    v = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(v)):
        raise ValueError("tally values must be finite")
    xs = v.reshape(-1).tolist()
    comps: list[float] = []
    for _ in range(64):
        s = math.fsum(xs + [-c for c in comps])
        if s == 0.0:
            return comps
        comps.append(s)
    raise ArithmeticError("exact summation did not terminate")  # pragma: no cover


def rows_to_partial(
    h0: int,
    h1: int,
    tally_rows: NDArray[np.float64],
    counter_rows: NDArray[np.int32],
    edep: list[NDArray[np.int64]],
    diagnostics: HistoryDiagnostics | None,
    meta: dict[str, Any] | None = None,
    channel_acc: NDArray[np.int64] | None = None,
) -> PartialTransport:
    """Reduce per-history rows ``(n, 6 + 2 G [+ C + 2])`` (fixed tallies, outside deposits,
    quantization residuals, and with scoring channels one residual column per channel but the
    count channels, the lookup out-of-domain count and the path-bound-exceeded flag) and
    ``(n, 9)`` (counters) of one range to a partial result."""
    return rows_to_partial_many(
        h0,
        h1,
        [exact_components(tally_rows[:, c]) for c in range(tally_rows.shape[1])],
        [int(x) for x in counter_rows.astype(np.int64).sum(axis=0)],
        edep,
        diagnostics,
        meta,
        channel_acc,
    )


def rows_to_partial_many(
    h0: int,
    h1: int,
    tally_components: list[list[float]],
    counter_sums: list[int],
    edep: list[NDArray[np.int64]],
    diagnostics: HistoryDiagnostics | None,
    meta: dict[str, Any] | None = None,
    channel_acc: NDArray[np.int64] | None = None,
) -> PartialTransport:
    """Assemble a :class:`PartialTransport` from already reduced pieces (chunked backends). The
    ``accumulator_overflow`` counter is set here: a voxel of ``edep`` or of a scoring channel at
    or above ``MAX_QUANTA`` (or negative, i.e. wrapped) invalidates the result."""
    counter_sums = list(counter_sums)
    grids = list(edep) + ([] if channel_acc is None else [channel_acc])
    bad = any(int(a.max()) >= MAX_QUANTA or int(a.min()) < 0 for a in grids if a.size)
    counter_sums[COUNTER_NAMES.index("accumulator_overflow")] += int(bad)
    return PartialTransport(
        h0, h1, tally_components, counter_sums, edep, diagnostics, dict(meta or {}), channel_acc
    )


def concat_partials(partials: list[PartialTransport]) -> PartialTransport:
    """One partial result for the contiguous range covered by ``partials`` (exact: expansions are
    concatenated, counters and integer grids are summed; diagnostics are dropped). Fails closed on
    gaps or overlaps."""
    parts = sorted(partials, key=lambda p: p.h0)
    if not parts:
        raise ValueError("no partial results")
    for a, b in zip(parts, parts[1:], strict=False):
        if a.h1 != b.h0:
            raise ValueError(f"partial results are not contiguous: {a.h1} then {b.h0}")
    n_cols = len(parts[0].tally_components)
    comps = [[c for p in parts for c in p.tally_components[col]] for col in range(n_cols)]
    counters = [sum(p.counter_sums[i] for p in parts) for i in range(len(COUNTER_NAMES))]
    edep = []
    for g in range(len(parts[0].edep)):
        acc = np.zeros_like(parts[0].edep[g], dtype=np.int64)
        for p in parts:
            acc += p.edep[g]
        edep.append(acc)
    chan = None
    if parts[0].channel_acc is not None:
        chan = np.zeros_like(parts[0].channel_acc, dtype=np.int64)
        for p in parts:
            assert p.channel_acc is not None
            chan += p.channel_acc
    return PartialTransport(parts[0].h0, parts[-1].h1, comps, counters, edep, None, {}, chan)


def merge_partials(
    partials: list[PartialTransport],
    n_histories: int,
    n_grids: int,
    n_channel_columns: int = 0,
) -> RawTransport:
    """Reduce the partial results of a complete, contiguous partition of ``[0, n_histories)``.

    ``n_channel_columns`` is 0 without scoring channels, else ``C + 2``: the residual columns of
    the ``C`` channels with a residual, the lookup out-of-domain count and the number of histories
    whose path exceeded the bound ``B_L``.

    Fails closed (``ValueError``) for gaps, overlaps, a different number of columns or grids.
    """
    parts = sorted(partials, key=lambda p: p.h0)
    pos = 0
    for p in parts:
        if p.h0 != pos or p.h1 < p.h0:
            raise ValueError(
                f"partial results do not tile [0, {n_histories}): gap or overlap at {pos}"
            )
        pos = p.h1
    if pos != n_histories:
        raise ValueError(f"partial results cover [0, {pos}), expected [0, {n_histories})")
    n_cols = N_FIXED_TALLIES + 2 * n_grids + n_channel_columns
    for p in parts:
        if (p.channel_acc is not None) != (n_channel_columns > 0):
            raise ValueError("partial result and configuration disagree about scoring channels")
    for p in parts:
        if (
            len(p.tally_components) != n_cols
            or len(p.counter_sums) != len(COUNTER_NAMES)
            or len(p.edep) != n_grids
        ):
            raise ValueError("partial result has the wrong number of tally columns or grids")
    totals = [math.fsum(c for p in parts for c in p.tally_components[col]) for col in range(n_cols)]
    counters = {
        name: int(sum(p.counter_sums[i] for p in parts)) for i, name in enumerate(COUNTER_NAMES)
    }
    edep = []
    for g in range(n_grids):
        acc = np.zeros_like(parts[0].edep[g], dtype=np.int64)
        for p in parts:  # integer sums are associative: any order gives the same grid
            acc += p.edep[g]
        edep.append(acc.astype(np.float64) * QUANTUM_MEV)
    tallies = dict(zip(TALLY_NAMES, totals[:N_FIXED_TALLIES], strict=True))
    base = N_FIXED_TALLIES + 2 * n_grids
    channels = None
    if n_channel_columns > 0:
        chan = np.zeros_like(parts[0].channel_acc, dtype=np.int64)
        for p in parts:
            assert p.channel_acc is not None
            chan += p.channel_acc
        channels = ChannelRaw(
            acc=chan,
            residual=totals[base : base + n_channel_columns - 2],
            lookup_out_of_domain=int(round(totals[base + n_channel_columns - 2])),
            path_bound_exceeded=int(round(totals[base + n_channel_columns - 1])),
        )
    return RawTransport(
        edep_mev=edep,
        tallies=tallies,
        outside_mev=totals[N_FIXED_TALLIES : N_FIXED_TALLIES + n_grids],
        quantization_mev=totals[N_FIXED_TALLIES + n_grids : base],
        counters=counters,
        channels=channels,
    )


def build_diagnostics(
    partials: list[PartialTransport],
    track_end_positions: bool,
    escape_records: bool,
    trace_histories: int,
) -> dict[str, Any]:
    """Diagnostics dictionary (the keys of the reference backend) from the partial results."""
    parts = sorted(partials, key=lambda p: p.h0)
    out: dict[str, Any] = {}
    if not (track_end_positions or escape_records or trace_histories > 0):
        return out
    diags = []
    for p in parts:
        if p.diagnostics is None:
            raise ValueError("diagnostics were requested but a partial result has none")
        diags.append(p.diagnostics)
    pos = np.concatenate([d.end_position_mm for d in diags], axis=0)
    direc = np.concatenate([d.end_direction for d in diags], axis=0)
    energy = np.concatenate([d.end_energy_mev for d in diags])
    code = np.concatenate([d.end_code for d in diags])
    if track_end_positions:
        out["end_position_mm"] = pos
        out["end_code"] = code
        out["end_direction"] = direc
        out["control_residual"] = np.concatenate([d.control_residual for d in diags])
        out["control_displacement"] = np.concatenate([d.control_displacement for d in diags])
        out["end_energy_mev"] = energy
    if escape_records:
        sel = np.nonzero(code == END_ESCAPED)[0]
        out["escape_history"] = sel.astype(np.int64)
        out["escape_position_mm"] = pos[sel]
        out["escape_direction"] = direc[sel]
        out["escape_energy_mev"] = energy[sel]
    if trace_histories > 0:
        ti = np.concatenate([d.trace_int for d in diags], axis=0)
        tf = np.concatenate([d.trace_float for d in diags], axis=0)
        table = np.concatenate([ti.astype(np.float64), tf], axis=1)
        out["trace"] = {name: table[:, i] for i, name in enumerate(TRACE_COLUMNS)}
        out["trace_columns"] = TRACE_COLUMNS
        k = min(trace_histories, len(code))
        out["trace_end_history"] = np.arange(k, dtype=np.int64)
        out["trace_end_code"] = code[:k].astype(np.int64)
        out["trace_end_energy_mev"] = energy[:k]
    return out
