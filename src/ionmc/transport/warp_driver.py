"""Warp transport driver: tables and geometry on a device, chunked launches, exact reduction.

``run_warp_range(eff, h0, h1, device)`` transports the histories ``[h0, h1)`` on one Warp
device ("cpu" or "cuda:0") and returns a :class:`~ionmc.transport.tally.PartialTransport`.
The range is launched in chunks of at most ``RunOptions.chunk_histories`` histories (one
thread per history). After every chunk the per-history tally rows are reduced exactly
(:func:`ionmc.transport.tally.exact_components`) and the counters summed in int64; the
energy-deposit grids (int64 fixed-point quanta, integer atomic adds: associative) stay on the
device across chunks and are copied once at the end. The result is therefore bit-identical for
any chunk size within a precision.
"""

from __future__ import annotations

import ctypes
import time
from typing import Any

import numpy as np
import warp as wp

from ionmc.config import EffectiveConfig
from ionmc.species import species_of_projectile
from ionmc.transport.channel_device import empty_channel_data, make_channel_data
from ionmc.transport.funcs import BIG_LENGTH_MM
from ionmc.transport.kernels import make_kernel_support, make_transport_kernel
from ionmc.transport.run import channel_columns
from ionmc.transport.tally import (
    COUNTER_NAMES,
    N_FIXED_TALLIES,
    NUCLEAR_COUNTER_NAMES,
    NUCLEAR_TALLY_NAMES,
    TRACE_N_CONTINUOUS,
    TRACE_N_DISCRETE,
    HistoryDiagnostics,
    PartialTransport,
    exact_components,
    rows_to_partial_many,
)

_CU_FUNC_ATTRIBUTE_NUM_REGS = 4
NUCLEAR_CUDA_CHUNK_CAP = 1 << 14
"""Largest chunk of the nuclear variant on CUDA (per-history stack scratch, decision 0041)."""
_TABLE_FIELDS = (
    "ln_s_mass", "r_mass", "f_mass", "d_f", "ln_e_of_r", "ln_e0", "inv_dln_e", "ln_r0",
    "inv_dln_r", "z_over_a", "inv_rho_xs",
)  # fmt: skip


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
    if eff.nuclear is not None:
        kernel = make_transport_kernel(real, wants_diagnostics(eff), True)
    else:
        kernel = make_transport_kernel(real, wants_diagnostics(eff))
    dev = wp.get_device(device)
    t0 = time.perf_counter()
    wp.load_module(kernel.module, device=dev)
    seconds = time.perf_counter() - t0
    return kernel, seconds, _register_count(kernel, dev)


def _nuclear_inputs(
    eff: EffectiveConfig, real: Any, device: str, chunk: int, k_hist: int, n_cols: int
) -> tuple[Any, Any, Any]:
    """Nuclear kernel inputs: ``(tables, nd, nuc_device)``. ``tables`` are the proton and deuteron
    transport tables concatenated along the material axis (row ``m + species * n_materials``),
    ``nd`` the ``NucData`` struct (nuclear arrays, constants, per-history scratch and trace
    buffers for ``k_hist`` histories) and ``nuc_device`` the uploaded :class:`NuclearDevice`."""
    from types import SimpleNamespace

    from ionmc.config import NUCLEAR_MAX_ENERGY_MEV
    from ionmc.transport import kernels_nuclear as kn
    from ionmc.transport.kernels_nuclear import (
        MAX_EVENT_ROWS,
        STACK_COLUMNS,
        TRACE_WIDTH,
        make_nuclear_support,
    )
    from ionmc.transport.nuclear_device import NuclearDevice

    tab, nuc = eff.tables, eff.nuclear
    td = tab.deuteron
    if nuc is None or td is None or td.n_e != tab.n_e or td.n_r != tab.n_r:
        raise ValueError("nuclear kernels need deuteron tables on the proton grid sizes")
    tp_w, td_w = tab.to_warp(device, wp.float64), td.to_warp(device, wp.float64)
    cat = SimpleNamespace(n_e=tab.n_e, n_r=tab.n_r)
    for f in _TABLE_FIELDS:
        a, b = getattr(tp_w, f).numpy(), getattr(td_w, f).numpy()
        setattr(cat, f, wp.array(np.concatenate([a, b]), dtype=wp.float64, device=device))
    dev = NuclearDevice.from_table(
        nuc.table, eff.geometry.materials, real=real, device=device
    )  # fmt: skip
    ns = make_nuclear_support(real)
    nd = ns.nuc()
    nd.n_grid, nd.kmax, nd.n_mat = dev.n_grid, dev.kmax, len(eff.geometry.materials)
    nd.nt_base = n_cols - len(NUCLEAR_TALLY_NAMES)
    nd.mass_d = wp.float64(td.projectile.mass_mev)
    nd.e_cut_d = wp.float64(nuc.e_cut_deuteron_mev)
    nd.e_source_max = wp.float64(NUCLEAR_MAX_ENERGY_MEV)
    nd.stack_cap = kn.STACK_CAPACITY
    nd.child_limit = kn.CHILD_LIMIT
    nd.secondary_nuclear = int(kn.SECONDARY_NUCLEAR)
    nd.max_parent_gen = kn.MAX_PARENT_GENERATION
    nd.ledger_tol = wp.float64(kn.LEDGER_TOL)
    for f in (
        "grid", "lam", "edges", "rpre", "recoil", "tconst", "m_res", "sigma", "sigma_win",
        "sigma_end", "cum_sigma", "mat_ntargets", "mat_target",
    ):  # fmt: skip
        setattr(nd, f, getattr(dev, f))
    nd.stack = wp.zeros((chunk, 32, STACK_COLUMNS), dtype=wp.float64, device=device)
    nd.evi = wp.zeros(chunk * 16, dtype=int, device=device)
    nd.evf = wp.zeros(chunk * 8, dtype=real, device=device)
    nd.prod = wp.zeros(chunk * 32 * 8, dtype=real, device=device)
    k = max(k_hist, 1)
    nd.ev_tr = wp.zeros((k, MAX_EVENT_ROWS, TRACE_WIDTH), dtype=wp.float64, device=device)
    nd.ev_n = wp.zeros(k, dtype=wp.int32, device=device)
    nd.sec_tr = wp.zeros((k, 32, TRACE_WIDTH), dtype=wp.float64, device=device)
    nd.sec_n = wp.zeros(k, dtype=wp.int32, device=device)
    return cat, nd, dev


def _nuclear_trace(nd: Any, k_hist: int) -> dict[str, Any]:
    """The nuclear trace of the first ``k_hist`` histories (rows in history, then event order)."""
    out: dict[str, Any] = {}
    for key, tr, cnt in (("events", nd.ev_tr, nd.ev_n), ("secondaries", nd.sec_tr, nd.sec_n)):
        a, n = tr.numpy(), cnt.numpy()
        rows = [a[k, : n[k]] for k in range(k_hist)]
        out[key] = np.concatenate(rows, axis=0) if rows else np.zeros((0, a.shape[2]))
    return out


def run_warp_range(eff: EffectiveConfig, h0: int, h1: int, device: str) -> PartialTransport:
    """Transport the histories ``[h0, h1)`` on ``device`` (see the module docstring)."""
    cfg = eff.requested
    run, src, ph, diag = cfg.run, cfg.source, cfg.physics, cfg.diagnostics
    real = wp.float32 if eff.precision == "float32" else wp.float64
    np_real = np.float32 if eff.precision == "float32" else np.float64
    t_start = time.perf_counter()
    kernel, compile_s, registers = load_kernel(eff, device)
    support = make_kernel_support(real)
    from ionmc.transport.funcs import make_transport_funcs

    v3 = make_transport_funcs(real).vec3
    use_diag = wants_diagnostics(eff)

    def vec(x: tuple[float, float, float]) -> Any:
        return v3(real(x[0]), real(x[1]), real(x[2]))

    geo, tab = eff.geometry, eff.tables
    ctl = support.control()
    key_lo, key_hi = run.seed & 0xFFFFFFFF, run.seed >> 32
    ctl.seed0 = wp.uint32(key_lo)
    ctl.seed1 = wp.uint32(key_hi)
    ctl.n_batches = run.n_batches
    ctl.max_steps = eff.max_steps
    ctl.mcs = 1 if ph.multiple_scattering else 0
    ctl.straggling = 1 if ph.straggling else 0
    ctl.trunc_diag = 1 if ph.truncated_hinge_diagnostic else 0
    ctl.straggle_gamma = 1 if ph.straggling_model == "bohr_gamma_v1" else 0
    ctl.max_pieces = eff.scoring_pieces
    ctl.nx, ctl.ny, ctl.nz = geo.shape
    ctl.n_grids = len(cfg.scoring)
    ctl.n_e = tab.n_e
    ctl.n_r = tab.n_r
    ctl.e_cut = wp.float64(ph.e_cut_mev)
    ctl.e_table_max = wp.float64(float(tab.e_max_mev.min()))
    ctl.e0 = wp.float64(src.kinetic_energy_mev)
    ctl.sigma_e = wp.float64(src.energy_sigma_mev)
    ctl.sigma_lat = real(src.lateral_sigma_mm)
    ctl.c_alpha = wp.float64(ph.range_alpha)
    ctl.c_rho_f = wp.float64(ph.range_rho_f_mm)
    ctl.c_frac = wp.float64(ph.max_energy_loss_fraction)
    ctl.c_smax = wp.float64(ph.max_step_mm)
    ctl.c_fshort = wp.float64(ph.short_step_fraction)
    ctl.mass = wp.float64(src.projectile.mass_mev)
    ctl.pos0 = vec(src.position_mm)
    ctl.dir = vec(eff.unit_direction)
    ctl.origin = vec(geo.origin_mm)
    ctl.spacing = vec(geo.spacing_mm)
    ctl.lo = vec(geo.lower_mm)
    ctl.hi = vec(geo.world_upper_mm)
    ctl.z_clip = real(BIG_LENGTH_MM if geo.z_exit_mm is None else float(geo.z_exit_mm))

    n_g = len(cfg.scoring)
    sizes = [g.n_voxels for g in cfg.scoring]
    offsets = np.concatenate(([0], np.cumsum(sizes)[:-1])).astype(np.int32)
    total_vox = int(sum(sizes))
    g_origin = np.array([list(g.origin_mm) for g in cfg.scoring])
    g_inv = np.array([[1.0 / x for x in g.spacing_mm] for g in cfg.scoring])
    g_spacing = np.array([list(g.spacing_mm) for g in cfg.scoring])
    g_shape = np.array([list(g.shape) for g in cfg.scoring], dtype=np.int32)

    nuclear = eff.nuclear is not None
    if nuclear:
        n_nuc_t = len(NUCLEAR_TALLY_NAMES)
        counter_names: tuple[str, ...] = COUNTER_NAMES + NUCLEAR_COUNTER_NAMES
    else:
        n_nuc_t = 0
        counter_names = COUNTER_NAMES
    t = tab.to_warp(device, wp.float64)  # energy and range tables in double in every variant
    arr_mat: Any = wp.array(np.ascontiguousarray(geo.material_index), dtype=wp.int32, device=device)
    arr_dens: Any = wp.array(
        np.ascontiguousarray(geo.densities_g_cm3(), dtype=np.float64),
        dtype=wp.float64,
        device=device,
    )
    arr_gorigin: Any = wp.array(g_origin.astype(np_real), dtype=real, device=device)
    arr_ginv: Any = wp.array(g_inv.astype(np_real), dtype=real, device=device)
    arr_gspacing: Any = wp.array(g_spacing.astype(np_real), dtype=real, device=device)
    arr_gshape: Any = wp.array(g_shape, dtype=int, device=device)
    arr_goff: Any = wp.array(offsets, dtype=int, device=device)
    edep = wp.zeros((run.n_batches, total_vox), dtype=wp.int64, device=device)

    n_chan_cols = channel_columns(eff)
    n_cols = N_FIXED_TALLIES + 2 * n_g + n_chan_cols + n_nuc_t
    if eff.channels is None:
        chan = empty_channel_data(support, device)
    else:
        chan = make_channel_data(
            support, eff.channels, tab, n_batches=run.n_batches, n_grids=n_g,
            a_nucleon=src.projectile.a, species_id=species_of_projectile(src.projectile).id,
            generation=0, device=device,
        )  # fmt: skip

    n_local = h1 - h0
    chunk = min(run.chunk_histories, max(n_local, 1))
    if nuclear and wp.get_device(device).is_cuda:
        chunk = min(chunk, NUCLEAR_CUDA_CHUNK_CAP)
    extra: list[Any] = []
    if nuclear:
        k_hist = max(0, min(h1, diag.trace_histories) - h0) if use_diag else 0
        t, nd, nuc_dev = _nuclear_inputs(eff, real, device, chunk, k_hist, n_cols)
        assert tab.deuteron is not None
        if eff.channels is None:
            chan_d = empty_channel_data(support, device)
        else:
            chan_d = make_channel_data(
                support, eff.channels, tab.deuteron, n_batches=run.n_batches, n_grids=n_g,
                a_nucleon=tab.deuteron.projectile.a, species_id=1, generation=0, device=device,
            )  # fmt: skip
            chan_d.acc = chan.acc  # one accumulator, both species
        extra = [nd, chan_d]
    tally_rows = wp.zeros((chunk, n_cols), dtype=wp.float64, device=device)
    counter_rows = wp.zeros((chunk, len(counter_names)), dtype=wp.int32, device=device)
    if use_diag:
        end_state = wp.zeros((chunk, 11), dtype=real, device=device)
        end_code = wp.zeros(chunk, dtype=wp.int32, device=device)
        k_loc = max(0, min(h1, diag.trace_histories) - h0)
        t_shape = (max(k_loc, 1), eff.max_steps if k_loc else 1)
        trace_i = wp.zeros((*t_shape, TRACE_N_DISCRETE), dtype=wp.int32, device=device)
        trace_f = wp.zeros((*t_shape, TRACE_N_CONTINUOUS), dtype=wp.float64, device=device)
        trace_n = wp.zeros(max(k_loc, 1), dtype=wp.int32, device=device)
        ctl.trace_k = wp.uint32(min(diag.trace_histories, h1))
        ctl.trace_h0 = wp.uint32(h0)
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
        ctl.trace_k = wp.uint32(0)
        ctl.trace_h0 = wp.uint32(0)

    comps: list[list[float]] = [[] for _ in range(n_cols)]
    counter_sums = [0] * len(counter_names)
    chunk_seconds: list[float] = []
    t_loop = time.perf_counter()
    for c0 in range(h0, h1, chunk):
        n = min(chunk, h1 - c0)
        t_c = time.perf_counter()
        tally_rows.zero_()
        counter_rows.zero_()
        ctl.h0 = wp.uint32(c0)
        wp.launch(
            kernel,
            dim=n,
            inputs=[
                ctl, arr_mat, arr_dens, t.ln_s_mass, t.r_mass, t.f_mass, t.d_f, t.ln_e_of_r,
                t.ln_e0,
                t.inv_dln_e, t.ln_r0, t.inv_dln_r, t.z_over_a, t.inv_rho_xs,
                arr_gorigin, arr_gspacing, arr_ginv, arr_gshape, arr_goff,
                edep, tally_rows, counter_rows, end_state, end_code, trace_i, trace_f, trace_n,
                chan, *extra,
            ],
            device=device,
        )  # fmt: skip
        trows = tally_rows.numpy()[:n]  # synchronises
        crows = counter_rows.numpy()[:n]
        chunk_seconds.append(time.perf_counter() - t_c)
        for c in range(n_cols):
            comps[c].extend(exact_components(trows[:, c]))
        for i, x in enumerate(crows.astype(np.int64).sum(axis=0)):
            counter_sums[i] += int(x)
        if use_diag:
            es = end_state.numpy()[:n].astype(np.float64)
            sl = slice(c0 - h0, c0 - h0 + n)
            host_pos[sl] = es[:, 0:3]
            host_dir[sl] = es[:, 3:6]
            host_energy[sl] = es[:, 6]
            host_ctrl[sl] = es[:, 7]
            host_ctrl_sum[sl] = es[:, 8:11]
            host_code[sl] = end_code.numpy()[:n].astype(np.int8)
    loop_s = time.perf_counter() - t_loop

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
        "compile_s": compile_s,
        "register_count": registers,
        "chunk_histories": chunk,
        "n_chunks": len(chunk_seconds),
        "chunk_seconds": chunk_seconds,
        "transport_s": loop_s,
        "total_s": time.perf_counter() - t_start,
        "h0": h0,
        "h1": h1,
    }
    if nuclear:
        meta["nuclear_device_sha256"] = nuc_dev.sha256
        if use_diag:
            meta["nuclear_trace"] = _nuclear_trace(nd, k_hist)
    channel_acc = None if eff.channels is None else chan.acc.numpy().astype(np.int64)
    return rows_to_partial_many(
        h0, h1, comps, counter_sums, edep_grids, diagnostics, meta, channel_acc
    )
