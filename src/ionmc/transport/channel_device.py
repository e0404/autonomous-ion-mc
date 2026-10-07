"""Device-side packing of the scoring-channel plan for the Warp kernels (decision 0040).

:func:`make_channel_data` turns a :class:`~ionmc.transport.channels.ChannelPlan` and the water
row of the transport tables into the ``ChannelData`` struct of
:func:`ionmc.transport.kernels.make_kernel_support` (int64 accumulators ``acc[B, sum_c size_c]``,
per-channel integer and float parameter rows, the species match matrix, the lookup tables as a
flat value array and the water row); :func:`empty_channel_data` is the struct of a run without
tallies (``n_ch = 0``: the kernels run exactly the code of the qualified path). All numbers that
define a bin come from the same host functions as in the Python reference
(:func:`ionmc.transport.channels.lookup_axis_params` and ``spectrum_axis_params``).

Column layout of ``ch_i``: kind, class mask, generation lo, generation hi, offset, lookup index,
residual column (-1: none), spectrum bins (0: none), spectrum log flag; of ``ch_f``: ``2^k``,
``2^-k``, spectrum ``a0`` and inverse bin width.
"""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np
import warp as wp

from ionmc.species import N_SPECIES_SLOTS, species_by_name
from ionmc.transport.channels import ChannelPlan, lookup_axis_params, spectrum_axis_params
from ionmc.transport.scoring_funcs import KIND_CODES
from ionmc.transport.tables import TransportTables
from ionmc.transport.tally import N_FIXED_TALLIES

CH_I_COLUMNS = 9
CH_F_COLUMNS = 4


def _wp(a: np.ndarray, dtype: Any, device: str) -> Any:
    return wp.array(np.ascontiguousarray(a), dtype=dtype, device=device)


def empty_channel_data(support: SimpleNamespace, device: str) -> Any:
    """Channel struct of a run without scoring channels (``n_ch = 0``, one-element arrays)."""
    c = support.chan()
    c.n_ch = 0
    c.acc = wp.zeros((1, 1), dtype=wp.int64, device=device)
    c.ch_i = wp.zeros((1, CH_I_COLUMNS), dtype=wp.int32, device=device)
    c.ch_f = wp.zeros((1, CH_F_COLUMNS), dtype=wp.float64, device=device)
    c.ch_sp = wp.zeros((1, N_SPECIES_SLOTS), dtype=wp.int32, device=device)
    c.ch_begin = wp.zeros(1, dtype=int, device=device)
    c.ch_end = wp.zeros(1, dtype=int, device=device)
    c.ln_s_w = wp.zeros(2, dtype=wp.float64, device=device)
    c.lk_i = wp.zeros((1, 3), dtype=int, device=device)
    c.lk_f = wp.zeros((1, 2), dtype=wp.float64, device=device)
    c.lk_row = wp.zeros((1, N_SPECIES_SLOTS), dtype=int, device=device)
    c.lk_vals = wp.zeros(1, dtype=wp.float64, device=device)
    return c


def make_channel_data(
    support: SimpleNamespace,
    plan: ChannelPlan,
    tables: TransportTables,
    *,
    n_batches: int,
    n_grids: int,
    a_nucleon: int,
    species_id: int,
    generation: int,
    device: str,
) -> Any:
    """Channel struct of ``plan`` on ``device`` (the accumulator ``acc`` is zeroed)."""
    if tables.water_ln_s_mass is None:
        raise ValueError("scoring channels need transport tables with a water row")
    n_ch = len(plan.channels)
    ch_i = np.zeros((n_ch, CH_I_COLUMNS), dtype=np.int32)
    ch_f = np.zeros((n_ch, CH_F_COLUMNS))
    for i, c in enumerate(plan.channels):
        ch_i[i, :7] = (
            KIND_CODES[c.kind], c.class_mask, c.gen_lo, min(c.gen_hi, 2**31 - 1), c.offset,
            c.lookup, c.residual_column,
        )  # fmt: skip
        ch_f[i, 0] = np.ldexp(1.0, c.k)
        ch_f[i, 1] = np.ldexp(1.0, -c.k)
        if c.spectrum is not None:
            a0, inv, nb, lg = spectrum_axis_params(c.spectrum)
            ch_i[i, 7], ch_i[i, 8] = nb, lg
            ch_f[i, 2], ch_f[i, 3] = a0, inv
    n_lk = max(len(plan.lookups), 1)
    lk_i = np.zeros((n_lk, 3), dtype=np.int64)
    lk_f = np.zeros((n_lk, 2))
    lk_row = np.full((n_lk, N_SPECIES_SLOTS), -1, dtype=np.int64)
    vals: list[float] = [0.0]
    for li, lk in enumerate(plan.lookups):
        a0, inv, n, lg = lookup_axis_params(lk)
        lk_i[li] = (n, lg, 0 if lk.axis == "energy_per_nucleon_mev" else 1)
        lk_f[li] = (a0, inv)
        for name, row in lk.values.items():
            lk_row[li, species_by_name(name).id] = len(vals)
            vals.extend(float(x) for x in row)
    c = support.chan()
    c.n_ch = n_ch
    c.n_res = plan.n_residual
    c.path_bound = wp.float64(plan.path_bound_mm)
    c.res_base = N_FIXED_TALLIES + 2 * n_grids  # first channel column of a history's tally row
    c.n_water = tables.n_water
    c.a_nuc = a_nucleon
    c.species = species_id
    c.gen = generation
    c.ln_e0_w = wp.float64(tables.water_ln_e0)
    c.inv_dln_e_w = wp.float64(tables.water_inv_dln_e)
    c.rho_w = wp.float64(tables.water_density_g_cm3)
    c.acc = wp.zeros((n_batches, plan.total_size), dtype=wp.int64, device=device)
    c.ch_i = _wp(ch_i, wp.int32, device)
    c.ch_f = _wp(ch_f, wp.float64, device)
    c.ch_sp = _wp(plan.species_match.astype(np.int32), wp.int32, device)
    c.ch_begin = _wp(np.array(plan.ch_begin, dtype=np.int64), int, device)
    c.ch_end = _wp(np.array(plan.ch_end, dtype=np.int64), int, device)
    c.ln_s_w = _wp(np.asarray(tables.water_ln_s_mass, dtype=np.float64), wp.float64, device)
    c.lk_i = _wp(lk_i, int, device)
    c.lk_f = _wp(lk_f, wp.float64, device)
    c.lk_row = _wp(lk_row, int, device)
    c.lk_vals = _wp(np.array(vals), wp.float64, device)
    return c


__all__ = ["empty_channel_data", "make_channel_data"]
