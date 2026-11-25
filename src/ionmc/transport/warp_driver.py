"""Warp transport driver: tables and geometry on a device, chunked launches, exact reduction.

``run_warp_range(eff, h0, h1, device)`` transports the histories ``[h0, h1)`` on one Warp
device ("cpu" or "cuda:0") and returns a :class:`~ionmc.transport.tally.PartialTransport`.
The range is launched in chunks of at most ``RunOptions.chunk_histories`` histories (one
thread per history). After every chunk the per-history tally rows are reduced exactly
(:func:`ionmc.transport.tally.exact_components`) and the counters summed in int64; the
energy-deposit grids (int64 fixed-point quanta, integer atomic adds: associative) stay on the
device across chunks and are copied once at the end. The result is therefore bit-identical for
any chunk size within a precision.

:class:`WarpContext` holds everything that does not change between runs (kernel, tables,
materials, geometry, grid descriptors and, for a beamlet set, the per-beamlet parameter and
offset arrays); it is built once per run, and a ``Simulation`` keeps it across runs
(``run_transport(eff, cache)``). Beamlet sets (V3-007): the kernel derives the beamlet of
history ``h`` by bisection on the offsets ``s_j`` and reads that beamlet's source parameters
through the same expressions as a pencil source; with an influence request the deposits of the
influence grid are also added to int64 group accumulators ``[G B, V_inf]`` (row
``(j - j0) B + batch``). A group is launched as chunks that never cross the group boundary, the
accumulators are copied to the host once per group, reduced to one
:class:`~ionmc.transport.tally.BeamletAccumulation` per beamlet and zeroed.
"""

from __future__ import annotations

import bisect
import ctypes
import math
import time
from typing import Any

import numpy as np
import warp as wp

from ionmc.config import RESIDUAL_SUBSTEPS, EffectiveConfig
from ionmc.sources import BeamletSet, PencilBeamSource
from ionmc.species import species_of_projectile
from ionmc.transport.alloc import alloc_device_zeros, note_alloc, release
from ionmc.transport.channel_device import empty_channel_data, make_channel_data
from ionmc.transport.funcs import BIG_LENGTH_MM
from ionmc.transport.kernels import make_kernel_support, make_transport_kernel
from ionmc.transport.run import channel_columns
from ionmc.transport.tally import (
    COUNTER_NAMES,
    ION_COUNTER_NAMES,
    N_FIXED_TALLIES,
    QUANTUM_MEV,
    TALLY_NAMES,
    TRACE_N_CONTINUOUS,
    TRACE_N_DISCRETE,
    BeamletAccumulation,
    BeamletBalance,
    BeamletRaw,
    HistoryDiagnostics,
    PartialTransport,
    exact_components,
    rows_to_partial_many,
)

_CU_FUNC_ATTRIBUTE_NUM_REGS = 4


def device_for_backend(backend: str) -> str:
    """Warp device string of a Warp backend name."""
    if backend == "warp-cpu":
        return "cpu"
    if backend == "warp-cuda":
        return "cuda:0"
    raise ValueError(f"not a Warp backend: {backend!r}")


def wants_diagnostics(eff: EffectiveConfig) -> bool:
    """True if any diagnostic output is requested (selects the ``diag`` kernel variant)."""
    d = eff.requested.diagnostics
    return bool(d.track_end_positions or d.escape_records or d.trace_histories > 0)


def _register_count(kernel: Any, device: Any) -> int | None:
    """Registers per thread of the compiled CUDA kernel (``None`` if not obtainable)."""
    if not device.is_cuda:
        return None
    try:
        hooks = kernel.module.get_kernel_hooks(kernel, device)
        lib = ctypes.CDLL("libcuda.so.1")
        value = ctypes.c_int(0)
        res = lib.cuFuncGetAttribute(
            ctypes.byref(value), _CU_FUNC_ATTRIBUTE_NUM_REGS, ctypes.c_void_p(hooks.forward)
        )
        return int(value.value) if res == 0 else None
    except Exception:  # pragma: no cover - needs a CUDA host
        return None


def load_kernel(eff: EffectiveConfig, device: str) -> tuple[Any, float, int | None]:
    """Compile and load the kernel for ``eff`` on ``device``: ``(kernel, seconds, registers)``."""
    real = wp.float32 if eff.precision == "float32" else wp.float64
    ion_eval = eff.ion_mapping.evaluation if eff.ion_mapping is not None else "table"
    kernel = make_transport_kernel(real, wants_diagnostics(eff), eff.is_ion, ion_eval)
    dev = wp.get_device(device)
    t0 = time.perf_counter()
    wp.load_module(kernel.module, device=dev)
    seconds = time.perf_counter() - t0
    return kernel, seconds, _register_count(kernel, dev)


CONTEXT_STATS = {"builds": 0, "table_uploads": 0, "geometry_uploads": 0}
"""Instrumented counters of the device context (V3-007 REUSE-1): a second run on the same
``Simulation`` leaves them unchanged."""


class WarpContext:
    """Device context of one validated configuration: built once, reused by every run."""

    def __init__(self, eff: EffectiveConfig, device: str) -> None:
        cfg = eff.requested
        self.eff, self.device = eff, device
        self._runs = 0
        self.real = real = wp.float32 if eff.precision == "float32" else wp.float64
        np_real = np.float32 if eff.precision == "float32" else np.float64
        self.kernel, self.compile_s, self.registers = load_kernel(eff, device)
        self.support = make_kernel_support(real)
        geo, tab = eff.geometry, eff.tables
        CONTEXT_STATS["builds"] += 1
        CONTEXT_STATS["table_uploads"] += 1
        CONTEXT_STATS["geometry_uploads"] += 1
        self.n_g = len(cfg.scoring)
        self.sizes = [g.n_voxels for g in cfg.scoring]
        self.offsets = np.concatenate(([0], np.cumsum(self.sizes)[:-1])).astype(np.int32)
        self.total_vox = int(sum(self.sizes))
        g_origin = np.array([list(g.origin_mm) for g in cfg.scoring])
        g_inv = np.array([[1.0 / x for x in g.spacing_mm] for g in cfg.scoring])
        g_spacing = np.array([list(g.spacing_mm) for g in cfg.scoring])
        g_shape = np.array([list(g.shape) for g in cfg.scoring], dtype=np.int32)
        self.tab_dev = tab.to_warp(device, wp.float64)  # tables in double in every variant
        dummy = np.zeros((1, 1))
        if eff.ion_mapping is not None:  # exact-Omega mapping data (double precision)
            im = eff.ion_mapping
            self.map_arrays: list[Any] = [
                wp.array(a, dtype=wp.float64, device=device) for a in (im.ion, im.proton, im.comp)
            ] + [wp.array(im.quad, dtype=wp.float64, device=device)]
            # the pv_eq/pv table rows (evaluation "table"); a dummy row for "root" (never read)
            self.map_arrays.append(
                wp.array(dummy if im.table is None else im.table, dtype=wp.float64, device=device)
            )
        else:  # dummy arrays: the proton kernel never reads them
            self.map_arrays = [
                wp.array(dummy, dtype=wp.float64, device=device) for _ in range(3)
            ] + [
                wp.array(dummy[0], dtype=wp.float64, device=device),
                wp.array(dummy, dtype=wp.float64, device=device),
            ]
        self.arr_mat: Any = wp.array(
            np.ascontiguousarray(geo.material_index), dtype=wp.int32, device=device
        )
        self.arr_dens: Any = wp.array(
            np.ascontiguousarray(geo.densities_g_cm3(), dtype=np.float64),
            dtype=wp.float64,
            device=device,
        )
        self.arr_gorigin: Any = wp.array(g_origin.astype(np_real), dtype=real, device=device)
        self.arr_ginv: Any = wp.array(g_inv.astype(np_real), dtype=real, device=device)
        self.arr_gspacing: Any = wp.array(g_spacing.astype(np_real), dtype=real, device=device)
        self.arr_gshape: Any = wp.array(g_shape, dtype=int, device=device)
        self.arr_goff: Any = wp.array(self.offsets, dtype=int, device=device)
        # per-beamlet parameter and offset arrays (dummies for a pencil source)
        src = cfg.source
        self.beamlets: BeamletSet | None = src if isinstance(src, BeamletSet) else None
        self.beamlet_offsets: list[int] = []
        if self.beamlets is not None:
            bl = self.beamlets.beamlets
            self.beamlet_offsets = [int(x) for x in self.beamlets.offsets()]
            self.b_off: Any = wp.array(
                np.array(self.beamlet_offsets, dtype=np.int64), dtype=wp.int64, device=device
            )
            self.b_e0: Any = wp.array(
                np.array([b.kinetic_energy_mev for b in bl]), dtype=wp.float64, device=device
            )
            self.b_sigma_e: Any = wp.array(
                np.array([b.energy_sigma_mev for b in bl]), dtype=wp.float64, device=device
            )
            self.b_sigma_lat: Any = wp.array(
                np.array([b.lateral_sigma_mm for b in bl]), dtype=wp.float64, device=device
            )
            self.b_pos: Any = wp.array(
                np.array([list(b.position_mm) for b in bl]), dtype=wp.float64, device=device
            )
            self.b_dir: Any = wp.array(
                np.array([list(b.unit_direction()) for b in bl]).astype(np_real),
                dtype=real,
                device=device,
            )
        else:
            self.b_off = wp.zeros(1, dtype=wp.int64, device=device)
            self.b_e0 = wp.zeros(1, dtype=wp.float64, device=device)
            self.b_sigma_e = wp.zeros(1, dtype=wp.float64, device=device)
            self.b_sigma_lat = wp.zeros(1, dtype=wp.float64, device=device)
            self.b_pos = wp.zeros((1, 3), dtype=wp.float64, device=device)
            self.b_dir = wp.zeros((1, 3), dtype=real, device=device)
        self.chan0: Any = empty_channel_data(self.support, device)  # n_ch = 0: never written
        self.n_cols = N_FIXED_TALLIES + 2 * self.n_g + channel_columns(eff)
        self.counter_names = COUNTER_NAMES + (ION_COUNTER_NAMES if eff.is_ion else ())
        self.infl_grid = -1
        if cfg.influence is not None:
            self.infl_grid = [g.name for g in cfg.scoring].index(cfg.influence.grid)

    # -- control struct -----------------------------------------------------------------------
    def _control(self, h0: int, h1: int) -> Any:
        eff = self.eff
        cfg = eff.requested
        run, src, ph, diag = cfg.run, cfg.source, cfg.physics, cfg.diagnostics
        real, geo, tab = self.real, eff.geometry, eff.tables

        def vec(x: tuple[float, float, float]) -> Any:
            return self.support_v3(real(x[0]), real(x[1]), real(x[2]))

        ctl = self.support.control()
        ctl.seed0 = wp.uint32(run.seed & 0xFFFFFFFF)
        ctl.seed1 = wp.uint32(run.seed >> 32)
        ctl.n_batches = run.n_batches
        ctl.max_steps = eff.max_steps
        ctl.mcs = 1 if ph.multiple_scattering else 0
        ctl.straggling = 1 if ph.straggling else 0
        ctl.trunc_diag = 1 if ph.truncated_hinge_diagnostic else 0
        ctl.straggle_gamma = 1 if ph.straggling_model == "bohr_gamma_v1" else 0
        ctl.max_pieces = eff.scoring_pieces
        ctl.nx, ctl.ny, ctl.nz = geo.shape
        ctl.n_grids = self.n_g
        ctl.n_e = tab.n_e
        ctl.n_r = tab.n_r
        ctl.e_cut = wp.float64(eff.e_cut_mev)
        ctl.e_table_max = wp.float64(float(tab.e_max_mev.min()))
        ctl.c_alpha = wp.float64(ph.range_alpha)
        ctl.c_rho_f = wp.float64(ph.range_rho_f_mm)
        ctl.c_frac = wp.float64(ph.max_energy_loss_fraction)
        ctl.c_smax = wp.float64(ph.max_step_mm)
        ctl.c_fshort = wp.float64(ph.short_step_fraction)
        ctl.mass = wp.float64(src.projectile.mass_mev)
        ctl.z = wp.float64(float(src.projectile.z))
        ctl.ion = 1 if eff.is_ion else 0
        ctl.n_sub = RESIDUAL_SUBSTEPS
        ctl.eor = 1 if eff.end_of_range_model == "continuous_residual_v1" else 0
        if isinstance(src, PencilBeamSource):
            ctl.n_beamlets = 0
            ctl.e0 = wp.float64(src.kinetic_energy_mev)
            ctl.sigma_e = wp.float64(src.energy_sigma_mev)
            ctl.sigma_lat = real(src.lateral_sigma_mm)
            ctl.pos0 = vec(src.position_mm)
            ctl.dir = vec(eff.single_unit_direction())
        else:  # the source parameters are read from the per-beamlet arrays
            ctl.n_beamlets = src.n_beamlets
            ctl.e0 = wp.float64(0.0)
            ctl.sigma_e = wp.float64(0.0)
            ctl.sigma_lat = real(0.0)
            ctl.pos0 = vec((0.0, 0.0, 0.0))
            ctl.dir = vec((0.0, 0.0, 1.0))
        ctl.infl_grid = self.infl_grid
        ctl.j0 = 0
        ctl.origin = vec(geo.origin_mm)
        ctl.spacing = vec(geo.spacing_mm)
        ctl.lo = vec(geo.lower_mm)
        ctl.hi = vec(geo.world_upper_mm)
        ctl.z_clip = real(BIG_LENGTH_MM if geo.z_exit_mm is None else float(geo.z_exit_mm))
        ctl.trace_k = wp.uint32(min(diag.trace_histories, h1) if wants_diagnostics(eff) else 0)
        ctl.trace_h0 = wp.uint32(h0 if wants_diagnostics(eff) else 0)
        return ctl

    @property
    def support_v3(self) -> Any:
        from ionmc.transport.funcs import make_transport_funcs

        return make_transport_funcs(self.real).vec3

    # -- groups ---------------------------------------------------------------------------------
    def _groups(self, h0: int, h1: int) -> list[tuple[int, int, int]]:
        """``(first history, end history, first beamlet j0)`` of every launch group: the whole
        range for a pencil / a composite run, the consecutive beamlet ranges of size G for an
        influence run."""
        if self.beamlets is None:
            return [(h0, h1, 0)]
        offs = self.beamlet_offsets
        if not 0 <= h0 <= h1 <= offs[-1]:
            raise ValueError(f"histories [{h0}, {h1}) outside the beamlet offsets {offs[-1]}")
        if self.infl_grid < 0:
            return [(h0, h1, 0)]
        if (h0, h1) != (0, offs[-1]):
            raise ValueError("an influence run transports the whole history range")
        g = self.eff.influence_group_size or 1
        j_n = len(offs) - 1
        return [(offs[j0], offs[min(j0 + g, j_n)], j0) for j0 in range(0, j_n, g)]

    def run_range(self, h0: int, h1: int) -> PartialTransport:
        eff = self.eff
        cfg = eff.requested
        run, diag = cfg.run, cfg.diagnostics
        real, device = self.real, self.device
        t_start = time.perf_counter()
        use_diag = wants_diagnostics(eff)
        ctl = self._control(h0, h1)
        n_g, n_cols, sizes, offsets = self.n_g, self.n_cols, self.sizes, self.offsets
        bs, boffs = self.beamlets, self.beamlet_offsets
        n_b = run.n_batches
        groups = self._groups(h0, h1)
        group_size = eff.influence_group_size or 1

        edep = alloc_device_zeros("edep", (n_b, self.total_vox), wp.int64, device)
        if self.infl_grid >= 0:
            v_inf = sizes[self.infl_grid]
            infl = alloc_device_zeros("influence", (group_size * n_b, v_inf), wp.int64, device)
        else:
            infl = wp.zeros((1, 1), dtype=wp.int64, device=device)  # dummy (F), never written
        if eff.channels is None:
            chan = self.chan0
        else:
            chan = make_channel_data(
                self.support, eff.channels, eff.tables, n_batches=n_b, n_grids=n_g,
                a_nucleon=cfg.source.projectile.a,
                species_id=species_of_projectile(cfg.source.projectile).id,
                generation=0, device=device,
            )  # fmt: skip

        n_local = h1 - h0
        longest = max((b - a for a, b, _ in groups), default=0)
        chunk = min(run.chunk_histories, max(longest, 1))
        tally_rows = alloc_device_zeros("tally_rows", (chunk, n_cols), wp.float64, device)
        counter_names = self.counter_names
        counter_rows = alloc_device_zeros(
            "counter_rows", (chunk, len(counter_names)), wp.int32, device
        )
        if use_diag:
            end_state = wp.zeros((chunk, 11), dtype=real, device=device)
            end_code = wp.zeros(chunk, dtype=wp.int32, device=device)
            k_loc = max(0, min(h1, diag.trace_histories) - h0)
            t_shape = (max(k_loc, 1), eff.max_steps if k_loc else 1)
            trace_i = wp.zeros((*t_shape, TRACE_N_DISCRETE), dtype=wp.int32, device=device)
            trace_f = wp.zeros((*t_shape, TRACE_N_CONTINUOUS), dtype=wp.float64, device=device)
            trace_n = wp.zeros(max(k_loc, 1), dtype=wp.int32, device=device)
            host_pos = np.zeros((n_local, 3))
            host_dir = np.zeros((n_local, 3))
            host_energy = np.zeros(n_local)
            host_ctrl = np.zeros(n_local)
            host_ctrl_sum = np.zeros((n_local, 3))
            host_code = np.zeros(n_local, dtype=np.int8)
        else:
            end_state = wp.zeros((1, 11), dtype=real, device=device)
            end_code = wp.zeros(1, dtype=wp.int32, device=device)
            trace_i = wp.zeros((1, 1, TRACE_N_DISCRETE), dtype=wp.int32, device=device)
            trace_f = wp.zeros((1, 1, TRACE_N_CONTINUOUS), dtype=wp.float64, device=device)
            trace_n = wp.zeros(1, dtype=wp.int32, device=device)

        t = self.tab_dev
        comps: list[list[float]] = [[] for _ in range(n_cols)]
        counter_sums = [0] * len(counter_names)
        # per-beamlet exact reductions (beamlet runs): expansions of every column and counters
        b_comps: dict[int, list[list[float]]] = {}
        b_counters: dict[int, list[int]] = {}
        accs: list[BeamletAccumulation] = []
        chunk_seconds: list[float] = []
        is_cuda = wp.get_device(device).is_cuda
        t_loop = time.perf_counter()
        for g_lo, g_hi, j0 in groups:
            ctl.j0 = j0
            for c0 in range(g_lo, g_hi, chunk):
                n = min(chunk, g_hi - c0)
                t_c = time.perf_counter()
                tally_rows.zero_()
                counter_rows.zero_()
                ctl.h0 = wp.uint32(c0)
                wp.launch(
                    self.kernel,
                    dim=n,
                    inputs=[
                        ctl, self.arr_mat, self.arr_dens, t.ln_s_mass, t.r_mass, t.f_mass, t.d_f,
                        t.ln_e_of_r, t.ln_e0, t.inv_dln_e, t.ln_r0, t.inv_dln_r, t.z_over_a,
                        t.inv_rho_xs, *self.map_arrays,
                        self.arr_gorigin, self.arr_gspacing, self.arr_ginv, self.arr_gshape,
                        self.arr_goff,
                        self.b_off, self.b_e0, self.b_sigma_e, self.b_sigma_lat, self.b_pos,
                        self.b_dir,
                        edep, infl, tally_rows, counter_rows, end_state, end_code, trace_i,
                        trace_f, trace_n, chan,
                    ],
                    device=device,
                )  # fmt: skip
                trows_full = tally_rows.numpy()  # synchronises; a host copy on CUDA only
                crows_full = counter_rows.numpy()
                if is_cuda:  # payload class (iv): the chunk's host rows (released below)
                    note_alloc("tally_rows_host", trows_full)
                    note_alloc("counter_rows_host", crows_full)
                trows = trows_full[:n]
                crows = crows_full[:n]
                chunk_seconds.append(time.perf_counter() - t_c)
                if bs is None:
                    for c in range(n_cols):
                        comps[c].extend(exact_components(trows[:, c]))
                    for i, x in enumerate(crows.astype(np.int64).sum(axis=0)):
                        counter_sums[i] += int(x)
                else:  # reduce every beamlet segment of the chunk exactly, then concatenate
                    j = bisect.bisect_right(boffs, c0) - 1
                    while j < len(boffs) - 1 and boffs[j] < c0 + n:
                        lo, hi = max(boffs[j], c0) - c0, min(boffs[j + 1], c0 + n) - c0
                        bc = b_comps.setdefault(j, [[] for _ in range(n_cols)])
                        bk = b_counters.setdefault(j, [0] * len(counter_names))
                        for c in range(n_cols):
                            e = exact_components(trows[lo:hi, c])
                            bc[c].extend(e)
                            comps[c].extend(e)
                        for i, x in enumerate(crows[lo:hi].astype(np.int64).sum(axis=0)):
                            bk[i] += int(x)
                            counter_sums[i] += int(x)
                        j += 1
                if use_diag:
                    es = end_state.numpy()[:n].astype(np.float64)
                    sl = slice(c0 - h0, c0 - h0 + n)
                    host_pos[sl] = es[:, 0:3]
                    host_dir[sl] = es[:, 3:6]
                    host_energy[sl] = es[:, 6]
                    host_ctrl[sl] = es[:, 7]
                    host_ctrl_sum[sl] = es[:, 8:11]
                    host_code[sl] = end_code.numpy()[:n].astype(np.int8)
                del trows, crows  # A-8: released right after the exact reduction, before the
                if is_cuda:  # next readback and before the group copy
                    release("tally_rows_host", trows_full)
                    release("counter_rows_host", crows_full)
                del trows_full, crows_full
            if self.infl_grid >= 0:  # one device->host copy per group, then the host reduction
                host = infl.numpy()
                if is_cuda:
                    note_alloc("influence_host", host)
                j_end = min(j0 + group_size, len(boffs) - 1)
                for j in range(j0, j_end):
                    loc = j - j0
                    accs.append(
                        BeamletAccumulation(j, host[loc * n_b : (loc + 1) * n_b].astype(np.int64))
                    )
                if is_cuda:
                    release("influence_host", host)
                del host
                infl.zero_()
        loop_s = time.perf_counter() - t_loop

        beamlets_raw = None
        if bs is not None:
            beamlets_raw = self._beamlet_raw(b_comps, b_counters, accs, n_cols)
        release("tally_rows", tally_rows)
        release("counter_rows", counter_rows)
        if self.infl_grid >= 0:
            release("influence", infl)

        edep_all = edep.numpy().astype(np.int64)
        edep_grids = [edep_all[:, o : o + s].copy() for o, s in zip(offsets, sizes, strict=True)]
        diagnostics = None
        if use_diag:
            if k_loc:
                n_rows = trace_n.numpy()
                ti, tf = trace_i.numpy(), trace_f.numpy()
                trace_int = np.concatenate([ti[k, : n_rows[k]] for k in range(k_loc)], axis=0)
                trace_float = np.concatenate([tf[k, : n_rows[k]] for k in range(k_loc)], axis=0)
            else:
                trace_int = np.zeros((0, TRACE_N_DISCRETE), dtype=np.int32)
                trace_float = np.zeros((0, TRACE_N_CONTINUOUS))
            diagnostics = HistoryDiagnostics(
                end_position_mm=host_pos,
                end_direction=host_dir,
                end_energy_mev=host_energy,
                end_code=host_code,
                control_residual=host_ctrl,
                control_displacement=host_ctrl_sum,
                trace_int=trace_int,
                trace_float=trace_float,
            )
        meta: dict[str, Any] = {
            "device": str(wp.get_device(device)),
            "compile_s": self.compile_s if self._runs == 0 else 0.0,
            "register_count": self.registers,
            "chunk_histories": chunk,
            "n_chunks": len(chunk_seconds),
            "chunk_seconds": chunk_seconds,
            "transport_s": loop_s,
            "total_s": time.perf_counter() - t_start,
            "h0": h0,
            "h1": h1,
            "context_reused": self._runs > 0,
        }
        self._runs += 1
        channel_acc = None if eff.channels is None else chan.acc.numpy().astype(np.int64)
        part = rows_to_partial_many(
            h0, h1, comps, counter_sums, edep_grids, diagnostics, meta, channel_acc
        )
        part.beamlets = beamlets_raw
        return part

    def _beamlet_raw(
        self,
        b_comps: dict[int, list[list[float]]],
        b_counters: dict[int, list[int]],
        accs: list[BeamletAccumulation],
        n_cols: int,
    ) -> BeamletRaw:
        """Per-beamlet balances from the exact per-beamlet expansions (the same correctly rounded
        sums as the reference backend's) and the accumulations of the groups."""
        n_g = self.n_g
        in_grid = {
            a.beamlet_id: sum(int(x) for x in a.batch_quanta.sum(axis=1)) * QUANTUM_MEV
            for a in accs
        }
        balances = []
        for j in sorted(b_comps):
            tot = [math.fsum(b_comps[j][c]) for c in range(n_cols)]
            boffs = self.beamlet_offsets
            balances.append(
                BeamletBalance(
                    beamlet_id=j,
                    n_histories=boffs[j + 1] - boffs[j],
                    tallies=dict(zip(TALLY_NAMES, tot[:N_FIXED_TALLIES], strict=True)),
                    outside_mev=tuple(tot[N_FIXED_TALLIES : N_FIXED_TALLIES + n_g]),
                    quantization_mev=tuple(tot[N_FIXED_TALLIES + n_g : N_FIXED_TALLIES + 2 * n_g]),
                    counters=dict(zip(self.counter_names, b_counters[j], strict=True)),
                    in_grid_mev=in_grid.get(j),
                )
            )
        inf = self.eff.requested.influence
        return BeamletRaw(
            offsets=np.array(self.beamlet_offsets, dtype=np.int64),
            group_size=self.eff.influence_group_size,
            balances=tuple(balances),
            accumulations=tuple(accs),
            influence_grid=None if inf is None else inf.grid,
        )


def get_context(eff: EffectiveConfig, device: str, cache: dict[str, Any] | None) -> WarpContext:
    """The device context of ``eff``: taken from ``cache`` (a ``Simulation``'s, keyed by device)
    or built; a run without a cache builds its own."""
    if cache is not None:
        ctx = cache.get(device)
        if ctx is not None and ctx.eff is eff:
            return ctx  # type: ignore[no-any-return]
    ctx = WarpContext(eff, device)
    if cache is not None:
        cache[device] = ctx
    return ctx


def run_warp_range(
    eff: EffectiveConfig, h0: int, h1: int, device: str, cache: dict[str, Any] | None = None
) -> PartialTransport:
    """Transport the histories ``[h0, h1)`` on ``device`` (see the module docstring)."""
    return get_context(eff, device, cache).run_range(h0, h1)
