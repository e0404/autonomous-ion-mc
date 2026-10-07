# mypy: ignore-errors
# (Warp kernel sources use runtime precision types in annotations; see _wpfunc.py.)
"""The Warp transport kernel: one thread transports one history (Warp CPU and CUDA).

``make_transport_kernel(real, diag)`` returns a cached kernel (``module="unique"``, backward
off) for precision ``real`` and diagnostics flag ``diag`` (end states and, in float64, the
per-step trace). Variants: float32, float32 with diagnostics, float64 with diagnostics (and
float64 without, for float64 production runs); only float64 is validated for traces.

The step loop below is a statement-by-statement copy of ``ionmc.transport.reference._history``:
the same branches in the same order and the same expression order, so that the float64 kernel
and the Python reference follow the same trajectory (test T1). All physical decisions are the
shared ``@wp.func`` functions of ``ionmc.physics`` and ``ionmc.transport.funcs``; the kernel
only holds the glue (table reads, state bookkeeping, random-number draws, scoring, tallies).
State lives in registers.

Outputs (no atomics except the energy-deposit grid): ``edep[B, sum(n_voxels)]`` of int64
fixed-point quanta (``QUANTUM_MEV`` = 2**-30 MeV), accumulated with integer ``atomic_add`` into
batch ``h % B``; per-history float64 ``tally_rows[chunk, 6 + 2 G]``
(``ionmc.transport.tally.TALLY_NAMES``, then the outside deposit of every grid, then the
quantization residual of every grid) and int32 ``counter_rows[chunk, 9]`` (the nine
``COUNTER_NAMES``; ``accumulator_overflow`` is set at the host reduction), each row written only
by its own thread; with ``diag``: ``end_state[chunk, 11]`` (position, direction, energy, control
residual, control displacement x, y, z), ``end_code[chunk]``, and for the first ``trace_k``
histories the trace ``trace_i[K, max_steps, 8]`` (int32), ``trace_f[K, max_steps, 9]``
(float64) and ``trace_n[K]`` (rows written). A history stops at ``max_steps`` so the trace buffer
cannot overflow.

Scoring channels (decision 0040): the last kernel argument is the ``ChannelData`` struct
(``ionmc.transport.channel_device``). Without tallies (``n_ch = 0``) no scoring code runs. With
channels, every in-grid piece is also added to the int64 channel accumulators ``chan.acc[B, sum
size]`` by ``score_piece`` (atomic adds), and the tally rows get ``n_res + 1`` extra columns after
the ``6 + 2 G`` ones: the quantization residual of every channel but N (per history) and the
lookup out-of-domain count.
"""

import functools
import math
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func
from ionmc.config import MAX_REJECTION_ATTEMPTS
from ionmc.physics.em import make_em
from ionmc.physics.kinematics import make_kinematics
from ionmc.rng.philox import PURPOSE_SOURCE, PURPOSE_TRANSPORT, make_philox
from ionmc.transport.funcs import make_transport_funcs
from ionmc.transport.scoring_funcs import make_scoring_funcs
from ionmc.transport.tally import (
    END_CUTOFF,
    END_ESCAPED,
    END_MISSED_WORLD,
    END_SOURCE_REJECTED,
    END_TRUNCATED,
    N_FIXED_TALLIES,
    QUANTUM_MEV,
    QUANTUM_SCALE,
)

wp.set_module_options({"enable_backward": False})

_TWO_PI = 2.0 * math.pi
_MAX_ATTEMPTS = MAX_REJECTION_ATTEMPTS


@functools.cache
def make_kernel_support(real: type) -> SimpleNamespace:
    """Control struct and helper functions of precision ``real`` (shared by the variants)."""
    name = check_real(real)
    R = real
    D = wp.float64  # energy, range and loss bookkeeping is double precision in every variant
    F = make_transport_funcs(real)
    FD = make_transport_funcs(D)
    v3 = F.vec3
    n_fixed = wp.constant(N_FIXED_TALLIES)

    class Control:
        pass

    Control.__annotations__ = {
        "seed0": wp.uint32,
        "seed1": wp.uint32,
        "h0": wp.uint32,
        "trace_k": wp.uint32,
        "trace_h0": wp.uint32,
        "n_batches": int,
        "max_steps": int,
        "mcs": int,
        "straggling": int,
        "trunc_diag": int,
        "straggle_gamma": int,
        "max_pieces": int,
        "z_clip": R,
        "nx": int,
        "ny": int,
        "nz": int,
        "n_grids": int,
        "n_e": int,
        "n_r": int,
        "e_cut": D,
        "e_table_max": D,
        "e0": D,
        "sigma_e": D,
        "sigma_lat": R,
        "c_alpha": D,
        "c_rho_f": D,
        "c_frac": D,
        "c_smax": D,
        "c_fshort": D,
        "mass": D,
        "pos0": v3,
        "dir": v3,
        "origin": v3,
        "spacing": v3,
        "lo": v3,
        "hi": v3,
    }
    Control.__name__ = f"TransportControl_{name}"
    Control.__qualname__ = Control.__name__
    control = wp.struct(Control)

    class Chan:
        pass

    Chan.__annotations__ = {
        "n_ch": int,  # 0 without tallies: the qualified path (no scoring code runs)
        "n_res": int,  # residual columns (the lookup out-of-domain count follows them)
        "res_base": int,  # first channel column of the per-history tally rows
        "n_water": int,
        "a_nuc": int,
        "species": int,  # species id and generation of the transported particle
        "gen": int,
        "ln_e0_w": D,
        "inv_dln_e_w": D,
        "rho_w": D,
        "acc": wp.array2d(dtype=wp.int64),
        "ch_i": wp.array2d(
            dtype=wp.int32
        ),  # kind, class, gen lo/hi, offset, lookup, res, bins, log
        "ch_f": wp.array2d(dtype=D),  # 2^k, 2^-k, spectrum a0, inv step
        "ch_sp": wp.array2d(dtype=wp.int32),  # species match
        "ch_begin": wp.array(dtype=int),
        "ch_end": wp.array(dtype=int),
        "ln_s_w": wp.array(dtype=D),
        "lk_i": wp.array2d(dtype=int),  # n, log, axis (0 energy per nucleon, 1 LET)
        "lk_f": wp.array2d(dtype=D),  # a0, inv step
        "lk_row": wp.array2d(dtype=int),  # offset of the species row in lk_vals (-1: none)
        "lk_vals": wp.array(dtype=D),
    }
    Chan.__name__ = f"ChannelData_{name}"
    Chan.__qualname__ = Chan.__name__
    chan_t = wp.struct(Chan)
    SC = make_scoring_funcs(D)
    half_d = wp.constant(wp.float64(0.5))
    ten_d = wp.constant(wp.float64(10.0))

    q_scale = wp.constant(wp.float64(QUANTUM_SCALE))
    q_mev = wp.constant(wp.float64(QUANTUM_MEV))

    @named_func(name)
    def step_state(chan: chan_t, e_mid: D, de_mean: D, s_act: D) -> tuple[D, D, D]:
        """Step quantities of the channel hook (as ``ReferenceChannelScorer.begin_step``): the
        water stopping power ``S_mid``, the ramp slope ``k`` and the energy rate ``Edot`` at the
        midpoint energy ``e_mid`` of a step of path length ``s_act`` and mean loss ``de_mean``."""
        iw, fw = FD.log_bin_index(e_mid, chan.ln_e0_w, chan.inv_dln_e_w, chan.n_water)
        ly0 = chan.ln_s_w[iw]
        ly1 = chan.ln_s_w[iw + 1]
        s_mass = FD.interp_exp(ly0, ly1, fw)
        gamma = SC.loglog_slope(ly0, ly1, chan.inv_dln_e_w)
        s_mid = s_mass * chan.rho_w / ten_d
        k = SC.let_ramp_slope(s_mid, gamma, de_mean, e_mid, s_act)
        return s_mid, k, de_mean / s_act

    @named_func(name)
    def score_piece(
        chan: chan_t,
        tally_rows: wp.array2d(dtype=wp.float64),
        tid: int,
        batch: int,
        g: int,
        vox: int,
        length: D,
        tau: D,
        eps: D,
        species: int,
        gen: int,
        cls: int,
        s_mid: D,
        k: D,
        e_mid: D,
        e_dot: D,
    ):
        """Add one piece (in grid ``g``, flat voxel ``vox``) to every channel of the grid that
        selects its class (1 step, 2 local), species and generation (as
        ``ReferenceChannelScorer.score_piece``: same functions, same operation order). The
        increment is computed in double precision in every kernel variant, quantized
        ``floor(x 2^k + 1/2)`` and added with an int64 atomic; the rounding residual and the
        lookup out-of-domain count go to the history's own tally row."""
        s_bar = D(0.0)
        e_bar = D(0.0)
        m1 = D(0.0)
        m2 = D(0.0)
        if cls == 1:
            s_bar, e_bar = SC.piece_state(s_mid, k, e_mid, e_dot, tau)
            m1, m2 = SC.piece_moments(s_bar, k, length)
            if s_bar <= D(0.0):  # compile-time envelope violated: invalidate the result
                c_neg = chan.res_base + chan.n_res
                tally_rows[tid, c_neg] = tally_rows[tid, c_neg] + D(1.0)
        for ci in range(chan.ch_begin[g], chan.ch_end[g]):
            sel = int(0)
            if (chan.ch_i[ci, 1] & cls) != 0 and chan.ch_sp[ci, species] != 0:
                if gen >= chan.ch_i[ci, 2] and gen <= chan.ch_i[ci, 3]:
                    sel = 1
            if sel == 1:
                kind = chan.ch_i[ci, 0]
                f = D(1.0)
                col = chan.ch_i[ci, 4] + vox
                if kind == 5:  # FE: lookup value at the piece argument
                    li = chan.ch_i[ci, 5]
                    xa = s_bar
                    if chan.lk_i[li, 2] == 0:
                        xa = e_bar / D(chan.a_nuc)
                    i, fr, inside = SC.lookup_bin(
                        xa, chan.lk_f[li, 0], chan.lk_f[li, 1], chan.lk_i[li, 0], chan.lk_i[li, 1]
                    )
                    if inside == 0:
                        c_ood = chan.res_base + chan.n_res
                        tally_rows[tid, c_ood] = tally_rows[tid, c_ood] + D(1.0)
                    row = chan.lk_row[li, species]
                    f = FD.lerp(chan.lk_vals[row + i], chan.lk_vals[row + i + 1], fr)
                if kind == 6:  # FL: energy bin of the piece (underflow 0, overflow n + 1)
                    nb = chan.ch_i[ci, 7]
                    b = SC.spectrum_bin(
                        e_bar / D(chan.a_nuc), chan.ch_f[ci, 2], chan.ch_f[ci, 3], nb,
                        chan.ch_i[ci, 8],
                    )  # fmt: skip
                    col = chan.ch_i[ci, 4] + vox * (nb + 2) + b
                x = SC.channel_value(kind, eps, length, m1, m2, s_bar, f)
                n = wp.int64(wp.floor(x * chan.ch_f[ci, 0] + half_d))
                wp.atomic_add(chan.acc, batch, col, n)
                rc = chan.ch_i[ci, 6]
                if rc >= 0:
                    c_res = chan.res_base + rc
                    tally_rows[tid, c_res] = tally_rows[tid, c_res] + (
                        x - wp.float64(n) * chan.ch_f[ci, 1]
                    )

    @named_func(name)
    def score_local(
        chan: chan_t,
        tally_rows: wp.array2d(dtype=wp.float64),
        tid: int,
        batch: int,
        g: int,
        ix: int,
        iy: int,
        iz: int,
        nxg: int,
        nyg: int,
        nzg: int,
        de: D,
    ):
        """Class "local" point deposit (cutoff energy, deposit at ``s_act = 0``; ``l = 0``) into
        the channels, if the voxel is inside the grid and channels are present."""
        if chan.n_ch > 0:
            if ix >= 0 and ix < nxg and iy >= 0 and iy < nyg and iz >= 0 and iz < nzg:
                score_piece(
                    chan, tally_rows, tid, batch, g, (ix * nyg + iy) * nzg + iz, D(0.0), D(0.0),
                    de, chan.species, chan.gen, 2, D(0.0), D(0.0), D(0.0), D(0.0),
                )  # fmt: skip

    @named_func(name)
    def deposit_voxel(
        edep: wp.array2d(dtype=wp.int64),
        tally_rows: wp.array2d(dtype=wp.float64),
        tid: int,
        batch: int,
        g: int,
        ix: int,
        iy: int,
        iz: int,
        nxg: int,
        nyg: int,
        nzg: int,
        off: int,
        n_grids: int,
        de: D,
    ):
        """Quantize and add ``de`` to voxel ``(ix, iy, iz)`` of grid ``g`` (int64 atomic add), or
        tally it as outside; the rounding residual is tallied per grid."""
        d64 = wp.float64(de)
        if ix >= 0 and ix < nxg and iy >= 0 and iy < nyg and iz >= 0 and iz < nzg:
            n = wp.int64(wp.floor(d64 * q_scale + wp.float64(0.5)))
            wp.atomic_add(edep, batch, off + (ix * nyg + iy) * nzg + iz, n)
            c = n_fixed + n_grids + g
            tally_rows[tid, c] = tally_rows[tid, c] + (d64 - wp.float64(n) * q_mev)
        else:
            tally_rows[tid, n_fixed + g] = tally_rows[tid, n_fixed + g] + d64

    @named_func(name)
    def deposit_point(
        edep: wp.array2d(dtype=wp.int64),
        tally_rows: wp.array2d(dtype=wp.float64),
        tid: int,
        batch: int,
        px: R,
        py: R,
        pz: R,
        de: D,
        g_origin: wp.array2d(dtype=R),
        g_inv: wp.array2d(dtype=R),
        g_shape: wp.array2d(dtype=int),
        g_off: wp.array(dtype=int),
        n_grids: int,
        chan: chan_t,
    ):
        """Point deposit (the energy left at the cutoff) in every grid."""
        p = v3(px, py, pz)
        for g in range(n_grids):
            go = v3(g_origin[g, 0], g_origin[g, 1], g_origin[g, 2])
            gi = v3(g_inv[g, 0], g_inv[g, 1], g_inv[g, 2])
            nxg = g_shape[g, 0]
            nyg = g_shape[g, 1]
            nzg = g_shape[g, 2]
            ix, iy, iz, inside = F.grid_index(p, go, gi, nxg, nyg, nzg)
            deposit_voxel(
                edep, tally_rows, tid, batch, g, ix, iy, iz, nxg, nyg, nzg, g_off[g], n_grids, de
            )
            score_local(chan, tally_rows, tid, batch, g, ix, iy, iz, nxg, nyg, nzg, de)

    @named_func(name)
    def ramp_weight(ta: D, tb: D, s_start: D, s_end: D, s_act: D) -> D:
        """Fraction of a step's energy between the path coordinates ``ta`` and ``tb`` for the linear
        stopping-power ramp ``s_start -> s_end`` (as ``reference.ramp_weight``)."""
        return (
            s_start * (tb - ta) + (s_end - s_start) * (tb * tb - ta * ta) / (D(2.0) * s_act)
        ) / (D(0.5) * (s_start + s_end) * s_act)

    @named_func(name)
    def deposit_leg(
        edep: wp.array2d(dtype=wp.int64),
        tally_rows: wp.array2d(dtype=wp.float64),
        tid: int,
        batch: int,
        g: int,
        px0: R,
        py0: R,
        pz0: R,
        ux: R,
        uy: R,
        uz: R,
        length: R,
        deposit: D,
        s_act: R,
        t_off: D,
        s_start: D,
        s_end: D,
        g_origin: wp.array2d(dtype=R),
        g_spacing: wp.array2d(dtype=R),
        g_inv: wp.array2d(dtype=R),
        g_shape: wp.array2d(dtype=int),
        g_off: wp.array(dtype=int),
        n_grids: int,
        max_pieces: int,
        chan: chan_t,
        s_mid: D,
        k_ramp: D,
        e_mid: D,
        e_dot: D,
    ) -> int:
        """Deposit the ramp-weighted share of ``deposit`` in every voxel of grid ``g`` crossed by
        the straight segment from ``p0`` along ``u`` (per-grid incremental DDA, ``seg_piece``);
        returns 1 if the walk exceeded ``max_pieces`` (the remainder goes to the last voxel)."""
        go = v3(g_origin[g, 0], g_origin[g, 1], g_origin[g, 2])
        gs = v3(g_spacing[g, 0], g_spacing[g, 1], g_spacing[g, 2])
        gi = v3(g_inv[g, 0], g_inv[g, 1], g_inv[g, 2])
        nxg = g_shape[g, 0]
        nyg = g_shape[g, 1]
        nzg = g_shape[g, 2]
        off = g_off[g]
        px = R(px0)
        py = R(py0)
        pz = R(pz0)
        uvec = v3(ux, uy, uz)
        ix, iy, iz, inside0 = F.grid_index(v3(px, py, pz), go, gi, nxg, nyg, nzg)
        remaining = R(length)
        tcur = D(t_off)
        for _it in range(max_pieces):
            if remaining > R(0.0):
                piece, axis = F.seg_piece(v3(px, py, pz), uvec, ix, iy, iz, go, gs, remaining)
                tb = tcur + D(piece)
                eps_p = deposit * ramp_weight(tcur, tb, s_start, s_end, D(s_act))
                deposit_voxel(
                    edep, tally_rows, tid, batch, g, ix, iy, iz, nxg, nyg, nzg, off, n_grids, eps_p
                )  # fmt: skip
                if chan.n_ch > 0:
                    if ix >= 0 and ix < nxg and iy >= 0 and iy < nyg and iz >= 0 and iz < nzg:
                        score_piece(
                            chan, tally_rows, tid, batch, g, (ix * nyg + iy) * nzg + iz, D(piece),
                            tcur + half_d * D(piece) - half_d * D(s_act), eps_p, chan.species,
                            chan.gen, 1, s_mid, k_ramp, e_mid, e_dot,
                        )  # fmt: skip
                tcur = tb
                remaining = remaining - piece
                if axis >= 0:
                    px = px + ux * piece
                    py = py + uy * piece
                    pz = pz + uz * piece
                    ua = R(ux)
                    if axis == 1:
                        ua = uy
                    if axis == 2:
                        ua = uz
                    upward = int(0)
                    step_i = int(-1)
                    if ua > R(0.0):
                        upward = 1
                        step_i = 1
                    if axis == 0:
                        px = F.plane_position(ix, upward, go[0], gs[0])
                        ix = ix + step_i
                    if axis == 1:
                        py = F.plane_position(iy, upward, go[1], gs[1])
                        iy = iy + step_i
                    if axis == 2:
                        pz = F.plane_position(iz, upward, go[2], gs[2])
                        iz = iz + step_i
        overflow = int(0)
        if remaining > R(0.0):
            overflow = 1
            eps_r = deposit * ramp_weight(tcur, tcur + D(remaining), s_start, s_end, D(s_act))
            deposit_voxel(
                edep, tally_rows, tid, batch, g, ix, iy, iz, nxg, nyg, nzg, off, n_grids, eps_r
            )  # fmt: skip
            if chan.n_ch > 0:
                if ix >= 0 and ix < nxg and iy >= 0 and iy < nyg and iz >= 0 and iz < nzg:
                    score_piece(
                        chan, tally_rows, tid, batch, g, (ix * nyg + iy) * nzg + iz,
                        D(remaining), tcur + half_d * D(remaining) - half_d * D(s_act), eps_r,
                        chan.species, chan.gen, 1, s_mid, k_ramp, e_mid, e_dot,
                    )  # fmt: skip
        return overflow

    @named_func(name)
    def deposit_step(
        edep: wp.array2d(dtype=wp.int64),
        tally_rows: wp.array2d(dtype=wp.float64),
        tid: int,
        batch: int,
        p0: v3,
        d0: v3,
        leg1: R,
        hinge: v3,
        d1: v3,
        leg2: R,
        deposit: D,
        s_act: R,
        s_start: D,
        s_end: D,
        g_origin: wp.array2d(dtype=R),
        g_spacing: wp.array2d(dtype=R),
        g_inv: wp.array2d(dtype=R),
        g_shape: wp.array2d(dtype=int),
        g_off: wp.array(dtype=int),
        n_grids: int,
        max_pieces: int,
        chan: chan_t,
        s_mid: D,
        k_ramp: D,
        e_mid: D,
        e_dot: D,
    ) -> int:
        """Track-length apportioning of a step deposit along both hinge legs in every grid;
        returns the number of legs that exceeded ``max_pieces``."""
        ovf = int(0)
        for g in range(n_grids):
            if s_act > R(0.0):
                ovf = ovf + deposit_leg(
                    edep, tally_rows, tid, batch, g, p0[0], p0[1], p0[2], d0[0], d0[1], d0[2],
                    leg1, deposit, s_act, D(0.0), s_start, s_end, g_origin, g_spacing, g_inv,
                    g_shape, g_off, n_grids, max_pieces, chan, s_mid, k_ramp, e_mid, e_dot,
                )  # fmt: skip
                ovf = ovf + deposit_leg(
                    edep, tally_rows, tid, batch, g, hinge[0], hinge[1], hinge[2], d1[0], d1[1],
                    d1[2], leg2, deposit, s_act, D(leg1), s_start, s_end, g_origin, g_spacing,
                    g_inv, g_shape, g_off, n_grids, max_pieces, chan, s_mid, k_ramp, e_mid, e_dot,
                )  # fmt: skip
            else:
                go = v3(g_origin[g, 0], g_origin[g, 1], g_origin[g, 2])
                gi = v3(g_inv[g, 0], g_inv[g, 1], g_inv[g, 2])
                nxg = g_shape[g, 0]
                nyg = g_shape[g, 1]
                nzg = g_shape[g, 2]
                ix, iy, iz, inside = F.grid_index(p0, go, gi, nxg, nyg, nzg)
                deposit_voxel(
                    edep, tally_rows, tid, batch, g, ix, iy, iz, nxg, nyg, nzg, g_off[g], n_grids,
                    deposit,
                )  # fmt: skip
                score_local(chan, tally_rows, tid, batch, g, ix, iy, iz, nxg, nyg, nzg, deposit)
        return ovf

    @named_func(name)
    def energy_from_range(
        ln_er: wp.array2d(dtype=D),
        ln_r0: wp.array(dtype=D),
        inv_dln_r: wp.array(dtype=D),
        n_r: int,
        m: int,
        r_g_cm2: D,
    ) -> D:
        """Energy [MeV] with CSDA range ``r_g_cm2`` (shared bin location, array reads; double)."""
        i, f = FD.log_bin_index(r_g_cm2, ln_r0[m], inv_dln_r[m], n_r)
        return FD.interp_exp(ln_er[m, i], ln_er[m, i + 1], f)

    return SimpleNamespace(
        control=control,
        chan=chan_t,
        score_piece=score_piece,
        step_state=step_state,
        deposit_point=deposit_point,
        deposit_step=deposit_step,
        energy_from_range=energy_from_range,
        real=name,
    )


@functools.cache
def make_transport_kernel(real: type, diag: bool):
    """Return the cached transport kernel for precision ``real`` and diagnostics flag ``diag``.

    ``enable_backward`` is off and ``module="unique"`` gives every variant its own Warp module
    (its own compile unit and cache entry).
    """
    name = check_real(real)
    R = real
    D = wp.float64  # energy / range bookkeeping in double precision (positions and angles in R)
    S = make_kernel_support(real)
    F = make_transport_funcs(real)
    FD = make_transport_funcs(D)
    EM = make_em(real)
    EMD = make_em(D)
    K = make_kinematics(D)
    PH = make_philox(real)
    v3 = F.vec3
    control = S.control
    deposit_point = S.deposit_point
    deposit_step = S.deposit_step
    step_state = S.step_state
    chan_t = S.chan
    energy_from_range = S.energy_from_range

    u_zero = wp.constant(wp.uint32(0))
    p_transport = wp.constant(wp.uint32(PURPOSE_TRANSPORT))
    p_source = wp.constant(wp.uint32(PURPOSE_SOURCE))
    two_pi = wp.constant(R(_TWO_PI))
    max_attempts = wp.constant(_MAX_ATTEMPTS)
    code_cutoff = wp.constant(END_CUTOFF)
    code_escaped = wp.constant(END_ESCAPED)
    code_truncated = wp.constant(END_TRUNCATED)
    code_rejected = wp.constant(END_SOURCE_REJECTED)
    code_missed = wp.constant(END_MISSED_WORLD)
    with_diag = wp.constant(diag)

    def transport(
        ctl: control,
        mat: wp.array3d(dtype=wp.int32),
        dens: wp.array3d(dtype=D),
        ln_s: wp.array2d(dtype=D),
        ln_r: wp.array2d(dtype=D),
        ln_er: wp.array2d(dtype=D),
        ln_e0: wp.array(dtype=D),
        inv_dln_e: wp.array(dtype=D),
        ln_r0: wp.array(dtype=D),
        inv_dln_r: wp.array(dtype=D),
        z_over_a: wp.array(dtype=D),
        inv_xs: wp.array(dtype=D),
        g_origin: wp.array2d(dtype=R),
        g_spacing: wp.array2d(dtype=R),
        g_inv: wp.array2d(dtype=R),
        g_shape: wp.array2d(dtype=int),
        g_off: wp.array(dtype=int),
        edep: wp.array2d(dtype=wp.int64),
        tally_rows: wp.array2d(dtype=wp.float64),
        counter_rows: wp.array2d(dtype=wp.int32),
        end_state: wp.array2d(dtype=R),
        end_code: wp.array(dtype=wp.int32),
        trace_i: wp.array3d(dtype=wp.int32),
        trace_f: wp.array3d(dtype=wp.float64),
        trace_n: wp.array(dtype=wp.int32),
        chan: chan_t,
    ):
        tid = wp.tid()
        h = ctl.h0 + wp.uint32(tid)
        batch = int(h % wp.uint32(ctl.n_batches))
        key = wp.vec2ui(ctl.seed0, ctl.seed1)
        nx = ctl.nx
        ny = ctl.ny
        nz = ctl.nz
        zero = R(0.0)
        zd = D(0.0)
        one_d = D(1.0)

        # per-history accumulators (a row depends on this history alone)
        t_initial = wp.float64(0.0)
        t_cutoff = wp.float64(0.0)
        t_step = wp.float64(0.0)
        t_escaped = wp.float64(0.0)
        t_truncated = wp.float64(0.0)
        c_trunc = int(0)
        c_stall = int(0)
        c_strag = int(0)
        c_src = int(0)
        c_inv = int(0)
        c_pieces = int(0)
        c_res = R(0.0)
        c_sx = R(0.0)
        c_sy = R(0.0)
        c_sz = R(0.0)
        code = int(-1)
        alive = int(1)

        # source sampling (purpose 1)
        w = PH.philox_block(h, u_zero, u_zero, p_source, key)
        z1, z2 = F.gauss_pair(PH.u01(w[0]), PH.u01(w[1]))
        z3 = F.gauss_one(PH.u01(w[2]), PH.u01(w[3]))
        d = ctl.dir
        e1v, e2v = F.orthonormal_basis(d)
        energy = D(ctl.e0)
        if ctl.sigma_e > zd:
            energy = energy + ctl.sigma_e * D(z3)
        px = R(ctl.pos0[0])
        py = R(ctl.pos0[1])
        pz = R(ctl.pos0[2])
        if ctl.sigma_lat > zero:
            sg = ctl.sigma_lat
            px = px + sg * (z1 * e1v[0] + z2 * e2v[0])
            py = py + sg * (z1 * e1v[1] + z2 * e2v[1])
            pz = pz + sg * (z1 * e1v[2] + z2 * e2v[2])
        ux = R(d[0])
        uy = R(d[1])
        uz = R(d[2])
        ix = int(0)
        iy = int(0)
        iz = int(0)
        p1v1 = D(0.0)
        t0 = R(0.0)
        src_ok = int(0)
        if energy >= ctl.e_cut and energy <= ctl.e_table_max:
            src_ok = 1
        if src_ok == 0:
            c_src = 1
            code = code_rejected
            alive = 0
        if alive == 1:
            t_initial = wp.float64(energy)
            p1v1 = K.pv_mev(energy, ctl.mass)
            # entry into the world (vacuum outside)
            t_in, t_out = F.ray_box(v3(px, py, pz), v3(ux, uy, uz), ctl.lo, ctl.hi)
            if t_in > t_out or t_out < zero:
                t_escaped = wp.float64(energy)
                code = code_missed
                alive = 0
            else:
                t0 = wp.max(t_in, zero)
                px = wp.min(wp.max(px + ux * t0, ctl.lo[0]), ctl.hi[0])
                py = wp.min(wp.max(py + uy * t0, ctl.lo[1]), ctl.hi[1])
                pz = wp.min(wp.max(pz + uz * t0, ctl.lo[2]), ctl.hi[2])
                ix = wp.min(wp.max(int(wp.floor((px - ctl.origin[0]) / ctl.spacing[0])), 0), nx - 1)
                iy = wp.min(wp.max(int(wp.floor((py - ctl.origin[1]) / ctl.spacing[1])), 0), ny - 1)
                iz = wp.min(wp.max(int(wp.floor((pz - ctl.origin[2]) / ctl.spacing[2])), 0), nz - 1)

        steps = int(0)
        birth = int(1)  # first step of the life: linearized analytic log-average of f_dM
        blocks = int(0)
        zero_run = int(0)
        while alive == 1:
            if energy <= ctl.e_cut:
                t_cutoff = t_cutoff + wp.float64(energy)
                deposit_point(
                    edep, tally_rows, tid, batch, px, py, pz, energy, g_origin, g_inv, g_shape,
                    g_off, ctl.n_grids, chan,
                )  # fmt: skip
                code = code_cutoff
                alive = 0
            if alive == 1 and steps >= ctl.max_steps:
                t_truncated = t_truncated + wp.float64(energy)
                c_trunc = c_trunc + 1
                code = code_truncated
                alive = 0
            if alive == 1:
                m = int(mat[ix, iy, iz])
                rho = dens[ix, iy, iz]
                ib, fb = FD.log_bin_index(energy, ln_e0[m], inv_dln_e[m], ctl.n_e)
                s0 = FD.interp_exp(ln_s[m, ib], ln_s[m, ib + 1], fb)
                r0 = FD.interp_exp(ln_r[m, ib], ln_r[m, ib + 1], fb)
                s_lin = s0 * rho / D(10.0)
                r_mm = r0 * D(10.0) / rho
                pvec = v3(px, py, pz)
                dvec = v3(ux, uy, uz)
                d_geo, axis_pre = F.dda_next_clip(
                    pvec, dvec, ix, iy, iz, ctl.origin, ctl.spacing, ctl.z_clip
                )
                s_el = FD.eloss_step_limit(energy, s_lin, ctl.c_frac)
                s_rg = FD.range_step_limit(r_mm, ctl.c_alpha, ctl.c_rho_f)
                s_d, reason = FD.select_step(D(d_geo), s_el, s_rg, ctl.c_smax)
                s = R(s_d)  # the step length in the geometry precision

                # block A
                wa = PH.philox_block(h, u_zero, wp.uint32(blocks), p_transport, key)
                blocks = blocks + 1
                ua0 = PH.u01(wa[0])
                ua1 = PH.u01(wa[1])
                ua2 = PH.u01(wa[2])
                leg1 = ua0 * s

                # multiple scattering
                d1x = R(ux)
                d1y = R(uy)
                d1z = R(uz)
                if ctl.mcs == 1:
                    e_mid = energy_from_range(
                        ln_er, ln_r0, inv_dln_r, ctl.n_r, m, r0 - rho * s_d / D(20.0)
                    )
                    inv_x = inv_xs[m]
                    pv_mid = K.pv_mev(e_mid, ctl.mass)
                    var = D(0.0)
                    if birth == 1:
                        e_end = wp.min(
                            energy_from_range(
                                ln_er, ln_r0, inv_dln_r, ctl.n_r, m, r0 - rho * s_d / D(10.0)
                            ),
                            energy,
                        )
                        var = EMD.scattering_variance_birth(
                            pv_mid, K.pv_mev(e_end, ctl.mass), p1v1, one_d, inv_x, rho, s_d
                        )
                    else:
                        t_pow = EMD.scattering_power_dm(pv_mid, p1v1, one_d, inv_x, rho)
                        var = t_pow * s_d
                    theta = EMD.polar_deflection(var, D(ua1))
                    nd = EM.rotate_dir(dvec, R(theta), two_pi * ua2)
                    d1x = nd[0]
                    d1y = nd[1]
                    d1z = nd[2]

                # hinge and second leg
                hx = px + ux * leg1
                hy = py + uy * leg1
                hz = pz + uz * leg1
                leg2 = R(0.0)
                axis2 = int(-1)
                if ctl.trunc_diag == 1:
                    # DIAGNOSTIC (T14 negative control, default off), "truncate-first": see
                    # reference._history; the planned step already ends at the straight-line plane
                    leg2 = s - leg1
                    if reason == 0:
                        axis2 = axis_pre
                else:
                    leg2, axis2 = F.leg2_limit_clip(
                        v3(hx, hy, hz), v3(d1x, d1y, d1z), ix, iy, iz, ctl.origin, ctl.spacing,
                        s - leg1, ctl.z_clip,
                    )  # fmt: skip
                nxp = hx + d1x * leg2
                nyp = hy + d1y * leg2
                nzp = hz + d1z * leg2
                exited = int(0)
                if axis2 == 3:  # the clip plane z = z_exit: the particle leaves the world there
                    if ctl.trunc_diag == 1 and s > zero:
                        c_res = wp.max(c_res, wp.abs(nzp - ctl.z_clip) / s)
                        c_sz = c_sz + (ctl.z_clip - nzp)
                    nzp = ctl.z_clip
                    exited = 1
                if axis2 >= 0 and axis2 < 3:
                    d1a = R(d1x)
                    if axis2 == 1:
                        d1a = d1y
                    if axis2 == 2:
                        d1a = d1z
                    if ctl.trunc_diag == 1:  # the plane crossed by the pre-hinge direction
                        d1a = R(ux)
                        if axis2 == 1:
                            d1a = uy
                        if axis2 == 2:
                            d1a = uz
                    upward = int(0)
                    if d1a > zero:
                        upward = 1
                    step_i = int(-1)
                    if upward == 1:
                        step_i = 1
                    if axis2 == 0:
                        plane_x = F.plane_position(ix, upward, ctl.origin[0], ctl.spacing[0])
                        if ctl.trunc_diag == 1 and s > zero:
                            c_res = wp.max(c_res, wp.abs(nxp - plane_x) / s)
                            c_sx = c_sx + (plane_x - nxp)
                        nxp = plane_x
                        ix = ix + step_i
                    if axis2 == 1:
                        plane_y = F.plane_position(iy, upward, ctl.origin[1], ctl.spacing[1])
                        if ctl.trunc_diag == 1 and s > zero:
                            c_res = wp.max(c_res, wp.abs(nyp - plane_y) / s)
                            c_sy = c_sy + (plane_y - nyp)
                        nyp = plane_y
                        iy = iy + step_i
                    if axis2 == 2:
                        plane_z = F.plane_position(iz, upward, ctl.origin[2], ctl.spacing[2])
                        if ctl.trunc_diag == 1 and s > zero:
                            c_res = wp.max(c_res, wp.abs(nzp - plane_z) / s)
                            c_sz = c_sz + (plane_z - nzp)
                        nzp = plane_z
                        iz = iz + step_i
                    inside = int(0)
                    if ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz:
                        inside = 1
                    if inside == 0:
                        exited = 1
                s_act = leg1 + leg2

                # energy loss
                s_act_d = D(s_act)
                tt = rho * s_act_d / D(10.0)
                e_r1 = energy_from_range(ln_er, ln_r0, inv_dln_r, ctl.n_r, m, r0 - tt)
                if ctl.c_fshort * r0 <= tt and tt < r0 and e_r1 > energy:
                    c_inv = c_inv + 1  # E1 > E0 from the inverse round trip
                    e_r1 = D(energy)
                s_branch = D(s0)  # stopping power of the linear branch: at the midpoint energy
                if tt < ctl.c_fshort * r0:
                    e_half = energy - D(0.5) * s0 * tt
                    ih, fh = FD.log_bin_index(e_half, ln_e0[m], inv_dln_e[m], ctl.n_e)
                    s_branch = FD.interp_exp(ln_s[m, ih], ln_s[m, ih + 1], fh)
                mean = EMD.csda_mean_loss(energy, e_r1, s_branch, tt, r0, ctl.c_fshort)
                attempts = int(1)
                loss = D(mean)
                if ctl.straggling == 1:
                    var_e = EMD.bohr_variance(
                        energy - D(0.5) * mean, ctl.mass, one_d, z_over_a[m], rho, s_act_d
                    )
                    accepted = int(0)
                    k = int(0)
                    while k < max_attempts and accepted == 0:
                        wb = PH.philox_block(h, u_zero, wp.uint32(blocks), p_transport, key)
                        blocks = blocks + 1
                        lw = D(0.0)
                        ok = int(0)
                        if ctl.straggle_gamma == 1:  # model bohr_gamma_v1
                            lw, ok = EMD.straggle_attempt_gamma(
                                mean,
                                var_e,
                                D(PH.u01(wb[0])),
                                D(PH.u01(wb[1])),
                                D(PH.u01(wb[2])),
                                D(PH.u01(wb[3])),
                            )
                        else:
                            lw, ok = EMD.straggle_attempt(
                                mean,
                                var_e,
                                D(PH.u01(wb[0])),
                                D(PH.u01(wb[1])),
                                D(PH.u01(wb[2])),
                                D(PH.u01(wb[3])),
                            )
                        attempts = k + 1
                        if ok == 1:
                            loss = lw
                            accepted = 1
                        k = k + 1
                    if accepted == 0:
                        c_strag = c_strag + 1
                        loss = D(mean)
                else:
                    blocks = blocks + 1  # one block is drawn and ignored
                loss = wp.min(loss, energy)
                e_new = energy - loss
                deposit = energy - e_new

                # scoring: the legs of a step with s_act > 0 are walked even when deposit = 0 if
                # scoring channels are present (fluence and LET do not depend on eps > 0)
                walk = int(0)
                if deposit > zd:
                    walk = 1
                if chan.n_ch > 0 and s_act > zero:
                    walk = 1
                if walk == 1:
                    ie, fe = FD.log_bin_index(e_new, ln_e0[m], inv_dln_e[m], ctl.n_e)
                    s_end_pw = FD.interp_exp(ln_s[m, ie], ln_s[m, ie + 1], fe)
                    ss_mid = zd
                    ss_k = zd
                    ss_e = zd
                    ss_dot = zd
                    if chan.n_ch > 0 and s_act > zero:
                        ss_e = energy - D(0.5) * mean
                        ss_mid, ss_k, ss_dot = step_state(chan, ss_e, mean, s_act_d)
                    c_pieces = c_pieces + deposit_step(
                        edep, tally_rows, tid, batch, v3(px, py, pz), v3(ux, uy, uz), leg1,
                        v3(hx, hy, hz), v3(d1x, d1y, d1z), leg2, deposit, s_act, s0, s_end_pw,
                        g_origin, g_spacing, g_inv, g_shape, g_off, ctl.n_grids, ctl.max_pieces,
                        chan, ss_mid, ss_k, ss_e, ss_dot,
                    )  # fmt: skip
                if deposit > zd:
                    t_step = t_step + wp.float64(deposit)

                px = R(nxp)
                py = R(nyp)
                pz = R(nzp)
                ux = R(d1x)
                uy = R(d1y)
                uz = R(d1z)
                energy = D(e_new)
                steps = steps + 1
                if s_act > zero:
                    birth = 0
                if with_diag:
                    if h < ctl.trace_k:
                        ti = int(h - ctl.trace_h0)
                        row = steps - 1
                        trace_i[ti, row, 0] = int(h)
                        trace_i[ti, row, 1] = steps
                        trace_i[ti, row, 2] = ix
                        trace_i[ti, row, 3] = iy
                        trace_i[ti, row, 4] = iz
                        trace_i[ti, row, 5] = reason
                        trace_i[ti, row, 6] = blocks
                        trace_i[ti, row, 7] = attempts
                        trace_f[ti, row, 0] = wp.float64(px)
                        trace_f[ti, row, 1] = wp.float64(py)
                        trace_f[ti, row, 2] = wp.float64(pz)
                        trace_f[ti, row, 3] = wp.float64(ux)
                        trace_f[ti, row, 4] = wp.float64(uy)
                        trace_f[ti, row, 5] = wp.float64(uz)
                        trace_f[ti, row, 6] = wp.float64(energy)
                        trace_f[ti, row, 7] = wp.float64(deposit)
                        trace_f[ti, row, 8] = wp.float64(s_act)
                        trace_n[ti] = steps
                if exited == 1:
                    t_escaped = t_escaped + wp.float64(energy)
                    code = code_escaped
                    alive = 0
                if alive == 1:
                    if s_act <= zero:
                        zero_run = zero_run + 1
                    else:
                        zero_run = 0
                    if zero_run > 3:
                        t_truncated = t_truncated + wp.float64(energy)
                        c_stall = c_stall + 1
                        code = code_truncated
                        alive = 0

        tally_rows[tid, 0] = t_initial
        tally_rows[tid, 1] = t_cutoff
        tally_rows[tid, 2] = t_step
        tally_rows[tid, 3] = t_escaped
        tally_rows[tid, 4] = t_truncated
        counter_rows[tid, 0] = c_trunc
        counter_rows[tid, 1] = c_stall
        counter_rows[tid, 2] = c_strag
        counter_rows[tid, 5] = c_src
        counter_rows[tid, 6] = c_inv
        counter_rows[tid, 8] = c_pieces
        if with_diag:
            end_state[tid, 0] = px
            end_state[tid, 1] = py
            end_state[tid, 2] = pz
            end_state[tid, 3] = ux
            end_state[tid, 4] = uy
            end_state[tid, 5] = uz
            end_state[tid, 6] = R(energy)
            end_state[tid, 7] = c_res
            end_state[tid, 8] = c_sx
            end_state[tid, 9] = c_sy
            end_state[tid, 10] = c_sz
            end_code[tid] = code

    transport.__name__ = f"transport_{name}_{'diag' if diag else 'plain'}"
    transport.__qualname__ = transport.__name__
    return wp.kernel(enable_backward=False, module="unique")(transport)
