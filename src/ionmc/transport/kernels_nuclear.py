# mypy: ignore-errors
# (Warp kernel sources use runtime precision types in annotations; see _wpfunc.py.)
"""The nuclear variant of the Warp transport kernel (decision 0041, V3-005B work package C11).

``make_nuclear_transport_kernel(real, diag)`` is the ``nuclear=True`` compile unit of
``ionmc.transport.kernels.make_transport_kernel(real, diag, nuclear=True)``. The electromagnetic
step (selection of the step, hinge legs, multiple scattering, energy loss, straggling, scoring)
is a statement-by-statement copy of the plain kernel; the nuclear code is added around it, in the
order of the Python reference ``ionmc.transport.reference._Reference`` (the semantic reference):

* a per-history LIFO particle stack ``nd.stack[chunk, 32, 16]`` (global scratch, float64 in every
  variant so that the energy ledger stays double precision; columns x, y, z, dx, dy, dz, T,
  species, genealogy id, generation, parent voxel ix, iy, iz, p1v1) with capacity 32
  (``queue_overflow`` as in Python) and the genealogy rule of ``child_genealogy_id``;
* species-aware stepping: the proton and deuteron tables are concatenated along the material
  axis by the driver (material row ``m + species * n_materials``), the deuteron cutoff and mass
  come from ``nd``; the scoring hook of the deuteron is the second channel struct ``chan_d``;
* majorant thinning for generation-0 protons only (``n_lambda`` bookkeeping, window or
  end-of-range majorant of the grid cell of ``E0`` by ``grid_locate``, candidate at the end of a
  nuclear-limited step, ``Sigma(E1) <= S^(E0)`` checked after EVERY primary step);
* the event of ``make_nuclear(real).sample_event`` on the PURPOSE_NUCLEAR stream with the slot
  layout of the Python reference (secondaries pushed with generation parent + 1, alpha and
  residual recoil deposited locally with species 64 / class local, neutron / gamma / binding /
  imbalance to the conditional tally columns, the per-event ledger).

All ledger and tally accumulation is float64 in both precisions. In the float32 variant the
event ledger ``binding`` and ``imbalance`` are recomputed in float64 from the sampler's float32
lab energies (the sampler's own float32 ledger entries carry ~1e-5 MeV of rounding that would
trip the 1e-9 conservation check); in float64 the sampler's entries are used as they are.
An event with more than 32 products (sampler status 2, never produced by the Python reference
below 80 products) fails closed like an exhausted rejection loop (``nuclear_rejection_limit``).

Counter columns after the nine of ``COUNTER_NAMES``: ``majorant_violation`` (9),
``nuclear_rejection_limit`` (10), ``nuclear_conservation`` (11). Tally columns: the unaccounted
column 5 is written, the six ``NUCLEAR_TALLY_NAMES`` are the last columns (``nd.nt_base``).
With ``diag`` and ``h < trace_k`` the nuclear trace ``nd.ev_tr`` (per accepted event: h, gid,
target, counts[5], Z_r, A_r, attempts, T1) and ``nd.sec_tr`` (per pushed secondary: h, gid,
parent gid, species, generation, T, direction, position) is written.
"""

import functools
import math
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func
from ionmc.config import MAX_REJECTION_ATTEMPTS
from ionmc.physics.em import make_em
from ionmc.physics.kinematics import make_kinematics
from ionmc.physics.nuclear import _nuc_u01_32, _nuc_u01_64, make_nuclear
from ionmc.rng.philox import PURPOSE_SOURCE, PURPOSE_TRANSPORT, make_philox
from ionmc.species import PSEUDO_BASE
from ionmc.transport.channels import MAX_PARTICLES
from ionmc.transport.funcs import make_transport_funcs
from ionmc.transport.kernels import make_kernel_support
from ionmc.transport.tally import (
    END_CUTOFF,
    END_ESCAPED,
    END_MISSED_WORLD,
    END_SOURCE_REJECTED,
    END_TRUNCATED,
)

_TWO_PI = 2.0 * math.pi
_MAX_ATTEMPTS = MAX_REJECTION_ATTEMPTS
END_NUCLEAR = 5
STACK_COLUMNS = 16
TRACE_WIDTH = 12
MAX_EVENT_ROWS = 2
STACK_CAPACITY = MAX_PARTICLES
CHILD_LIMIT = 31
LEDGER_TOL = 1.0e-9
"""Runtime limits written to ``NucData`` by the driver (tests lower them to force the counters)."""


@functools.cache
def make_nuclear_support(real: type) -> SimpleNamespace:
    """The ``NucData`` struct (nuclear arrays of precision ``real``, constants, scratch)."""
    name = check_real(real)
    R = real
    D = wp.float64

    class NucData:
        pass

    NucData.__annotations__ = {
        "n_grid": int,
        "kmax": int,
        "nt_base": int,
        "mass_d": D,
        "e_cut_d": D,
        "e_source_max": D,
        "n_mat": int,
        "stack_cap": int,
        "child_limit": int,
        "ledger_tol": D,
        "grid": wp.array(dtype=R),
        "lam": wp.array(dtype=R),
        "edges": wp.array(dtype=R),
        "rpre": wp.array(dtype=R),
        "recoil": wp.array(dtype=R),
        "tconst": wp.array(dtype=R),
        "m_res": wp.array(dtype=R),
        "sigma": wp.array(dtype=R),
        "sigma_win": wp.array(dtype=R),
        "sigma_end": wp.array(dtype=R),
        "cum_sigma": wp.array(dtype=R),
        "mat_ntargets": wp.array(dtype=wp.int32),
        "mat_target": wp.array(dtype=wp.int32),
        "stack": wp.array3d(dtype=D),
        "evi": wp.array(dtype=int),
        "evf": wp.array(dtype=R),
        "prod": wp.array(dtype=R),
        "ev_tr": wp.array3d(dtype=D),
        "ev_n": wp.array(dtype=wp.int32),
        "sec_tr": wp.array3d(dtype=D),
        "sec_n": wp.array(dtype=wp.int32),
    }
    NucData.__name__ = f"NuclearData_{name}"
    NucData.__qualname__ = NucData.__name__
    return SimpleNamespace(nuc=wp.struct(NucData), real=name)


@functools.cache
def make_nuclear_transport_kernel(real: type, diag: bool):
    """Return the cached nuclear transport kernel for precision ``real`` and diagnostics flag."""
    name = check_real(real)
    R = real
    D = wp.float64
    S = make_kernel_support(real)
    NS = make_nuclear_support(real)
    F = make_transport_funcs(real)
    FD = make_transport_funcs(D)
    EM = make_em(real)
    EMD = make_em(D)
    K = make_kinematics(D)
    PH = make_philox(real)
    NU = make_nuclear(real)
    ND = make_nuclear(D)
    v3 = F.vec3
    control = S.control
    deposit_point = S.deposit_point
    deposit_step = S.deposit_step
    step_state = S.step_state
    chan_t = S.chan
    nuc_t = NS.nuc
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
    code_nuclear = wp.constant(END_NUCLEAR)
    pseudo_local = wp.constant(PSEUDO_BASE)
    is_f32 = wp.constant(1 if name == "float32" else 0)
    with_diag = wp.constant(diag)

    nuc_u01 = _nuc_u01_64 if name == "float64" else _nuc_u01_32

    @named_func(name)
    def nuc_u(h: wp.uint32, gid: wp.uint32, blk: int, word: int, key: wp.vec2ui) -> R:
        """Uniform ``word`` of block ``blk`` of the PURPOSE_NUCLEAR stream of ``(h, gid)``."""
        return nuc_u01(h, gid, blk, word, key)

    @named_func(name)
    def sigma_at(e: D, grid: wp.array(dtype=R), n: int, row: wp.array(dtype=R), base: int) -> D:
        """``Sigma_mass(E)`` of a material row: step lookup of the grid cell by ``grid_locate``,
        lin-lin in E, clamped (as ``MaterialNuclear.sigma_at``), evaluated in double precision."""
        k = NU.grid_locate(R(e), grid, n)
        g0 = D(grid[k])
        g1 = D(grid[k + 1])
        t = wp.min(wp.max((e - g0) / (g1 - g0), D(0.0)), D(1.0))
        return (D(1.0) - t) * D(row[base + k]) + t * D(row[base + k + 1])

    def transport(
        ctl: control,
        mat: wp.array3d(dtype=wp.int32),
        dens: wp.array3d(dtype=D),
        ln_s: wp.array2d(dtype=D),
        r_mass: wp.array2d(dtype=D),
        f_mass: wp.array2d(dtype=D),
        d_f: wp.array2d(dtype=D),
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
        nd: nuc_t,
        chan_d: chan_t,
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
        t_unacc = wp.float64(0.0)
        n_local = wp.float64(0.0)
        n_alpha = wp.float64(0.0)
        n_neutron = wp.float64(0.0)
        n_gamma = wp.float64(0.0)
        n_binding = wp.float64(0.0)
        n_imbal = wp.float64(0.0)
        c_trunc = int(0)
        c_stall = int(0)
        c_strag = int(0)
        c_gen = int(0)
        c_queue = int(0)
        c_src = int(0)
        c_inv = int(0)
        c_pieces = int(0)
        c_major = int(0)
        c_rejlim = int(0)
        c_cons = int(0)
        c_res = R(0.0)
        c_sx = R(0.0)
        c_sy = R(0.0)
        c_sz = R(0.0)
        code = int(-1)
        alive = int(1)
        prim_done = int(0)
        n_st = int(0)
        n_ev_rows = int(0)
        n_sec_rows = int(0)
        gid = int(0)
        gen = int(0)
        species = int(0)

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
        if energy > nd.e_source_max:
            # nuclear runs: the cross sections end at 250 MeV (decision 0041 section 5): the
            # energy is booked as initial and unaccounted so that the ledger still closes
            c_src = 1
            t_initial = wp.float64(energy)
            t_unacc = wp.float64(energy)
            code = code_rejected
            alive = 0
        if alive == 1:
            if energy < ctl.e_cut or energy > ctl.e_table_max:
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

        path_mm = D(0.0)
        path_flag = int(0)
        cur_valid = int(alive)
        while cur_valid == 1:
            # ---- one particle of the stack (the step loop of decision 0039) ----
            alive = int(1)
            cc = chan
            mass = D(ctl.mass)
            ecut = D(ctl.e_cut)
            m_off = int(0)
            if species == 1:
                cc = chan_d
                mass = nd.mass_d
                ecut = nd.e_cut_d
                m_off = nd.n_mat
            cc.species = species
            cc.gen = gen
            gid_u = wp.uint32(gid)
            nuc_on = int(0)
            n_lam = D(0.0)
            nc = int(0)
            if gen == 0:
                nuc_on = 1
                u_birth = nuc_u(h, gid_u, 0, 0, key)
                n_lam = -wp.log(D(u_birth))
                nc = 1
            steps = int(0)
            birth = int(1)  # first step of the life: linearized analytic log-average of f_dM
            blocks = int(0)
            zero_run = int(0)
            while alive == 1:
                if energy <= ecut:
                    t_cutoff = t_cutoff + wp.float64(energy)
                    deposit_point(
                        edep, tally_rows, tid, batch, px, py, pz, energy, g_origin, g_inv,
                        g_shape, g_off, ctl.n_grids, cc,
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
                    mt = m + m_off
                    rho = dens[ix, iy, iz]
                    ib, fb = FD.log_bin_index(energy, ln_e0[mt], inv_dln_e[mt], ctl.n_e)
                    s0 = FD.interp_exp(ln_s[mt, ib], ln_s[mt, ib + 1], fb)
                    r0 = FD.range_in_bin(
                        r_mass[mt, ib], f_mass[mt, ib], d_f[mt, ib], D(1.0) / inv_dln_e[mt], fb
                    )
                    s_lin = s0 * rho / D(10.0)
                    r_mm = r0 * D(10.0) / rho
                    pvec = v3(px, py, pz)
                    dvec = v3(ux, uy, uz)
                    d_geo, axis_pre = F.dda_next_clip(
                        pvec, dvec, ix, iy, iz, ctl.origin, ctl.spacing, ctl.z_clip
                    )
                    s_el = FD.eloss_step_limit(energy, s_lin, ctl.c_frac)
                    s_rg = FD.range_step_limit(r_mm, ctl.c_alpha, ctl.c_rho_f)
                    s_hat = D(0.0)
                    s_d = D(0.0)
                    reason = int(0)
                    if nuc_on == 1:
                        # majorant of the step: window value at the grid cell of E0 (step lookup)
                        # or the end-of-range value for a range-limited step
                        k0 = NU.grid_locate(R(energy), nd.grid, nd.n_grid)
                        if s_rg < s_el:
                            s_hat = D(nd.sigma_end[m * nd.n_grid + k0])
                        else:
                            s_hat = D(nd.sigma_win[m * nd.n_grid + k0])
                        d_nuc = ND.nuclear_step_limit(n_lam, rho, s_hat)
                        s_d, reason = FD.select_step_nuclear(
                            D(d_geo), s_el, s_rg, ctl.c_smax, d_nuc
                        )
                    else:
                        s_d, reason = FD.select_step(D(d_geo), s_el, s_rg, ctl.c_smax)
                    s = R(s_d)  # the step length in the geometry precision

                    # block A
                    wa = PH.philox_block(h, gid_u, wp.uint32(blocks), p_transport, key)
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
                            ln_er, ln_r0, inv_dln_r, ctl.n_r, mt, r0 - rho * s_d / D(20.0)
                        )
                        inv_x = inv_xs[mt]
                        pv_mid = K.pv_mev(e_mid, mass)
                        var = D(0.0)
                        if birth == 1:
                            e_end = wp.min(
                                energy_from_range(
                                    ln_er, ln_r0, inv_dln_r, ctl.n_r, mt, r0 - rho * s_d / D(10.0)
                                ),
                                energy,
                            )
                            var = EMD.scattering_variance_birth(
                                pv_mid, K.pv_mev(e_end, mass), p1v1, one_d, inv_x, rho, s_d
                            )
                        else:
                            t_pow = EMD.scattering_power_dm(pv_mid, p1v1, one_d, inv_x, rho)
                            var = t_pow * s_d
                        theta = EMD.polar_deflection(var, D(ua1))
                        nd_dir = EM.rotate_dir(dvec, R(theta), two_pi * ua2)
                        d1x = nd_dir[0]
                        d1y = nd_dir[1]
                        d1z = nd_dir[2]

                    # hinge and second leg
                    hx = px + ux * leg1
                    hy = py + uy * leg1
                    hz = pz + uz * leg1
                    leg2 = R(0.0)
                    axis2 = int(-1)
                    if ctl.trunc_diag == 1:
                        # DIAGNOSTIC (T14 negative control, default off), "truncate-first": see
                        # reference._history; the planned step already ends at the straight-line
                        # plane
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
                    e_r1 = energy_from_range(ln_er, ln_r0, inv_dln_r, ctl.n_r, mt, r0 - tt)
                    if ctl.c_fshort * r0 <= tt and tt < r0 and e_r1 > energy:
                        c_inv = c_inv + 1  # E1 > E0 from the inverse round trip
                        e_r1 = D(energy)
                    s_branch = D(s0)  # stopping power of the linear branch: at the midpoint energy
                    if tt < ctl.c_fshort * r0:
                        e_half = energy - D(0.5) * s0 * tt
                        ih, fh = FD.log_bin_index(e_half, ln_e0[mt], inv_dln_e[mt], ctl.n_e)
                        s_branch = FD.interp_exp(ln_s[mt, ih], ln_s[mt, ih + 1], fh)
                    mean = EMD.csda_mean_loss(energy, e_r1, s_branch, tt, r0, ctl.c_fshort)
                    attempts = int(1)
                    loss = D(mean)
                    if ctl.straggling == 1:
                        var_e = EMD.bohr_variance(
                            energy - D(0.5) * mean, mass, one_d, z_over_a[mt], rho, s_act_d
                        )
                        accepted = int(0)
                        k = int(0)
                        while k < max_attempts and accepted == 0:
                            wb = PH.philox_block(h, gid_u, wp.uint32(blocks), p_transport, key)
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
                        ie, fe = FD.log_bin_index(e_new, ln_e0[mt], inv_dln_e[mt], ctl.n_e)
                        s_end_pw = FD.interp_exp(ln_s[mt, ie], ln_s[mt, ie + 1], fe)
                        ss_mid = zd
                        ss_k = zd
                        ss_e = zd
                        ss_dot = zd
                        if chan.n_ch > 0 and s_act > zero:
                            ss_e = energy - D(0.5) * mean
                            ss_mid, ss_k, ss_dot = step_state(cc, ss_e, mean, s_act_d)
                        c_pieces = c_pieces + deposit_step(
                            edep, tally_rows, tid, batch, v3(px, py, pz), v3(ux, uy, uz), leg1,
                            v3(hx, hy, hz), v3(d1x, d1y, d1z), leg2, deposit, s_act, s0, s_end_pw,
                            g_origin, g_spacing, g_inv, g_shape, g_off, ctl.n_grids, ctl.max_pieces,
                            cc, ss_mid, ss_k, ss_e, ss_dot,
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
                    if chan.n_ch > 0:
                        path_mm = path_mm + s_act_d  # the scored path of this history (L channel)
                        if path_mm > chan.path_bound:
                            path_flag = 1
                    if s_act > zero:
                        birth = 0
                    if with_diag:
                        if h < ctl.trace_k and gen == 0:
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
                    if nuc_on == 1:
                        n_lam = n_lam - rho * s_hat * s_act_d / D(10.0)
                        # majorant check after EVERY step (decision 0041 section 2): the Gamma
                        # straggling tail is unbounded, so Sigma(E1) <= S^(E0) is not guaranteed
                        sig_e = sigma_at(energy, nd.grid, nd.n_grid, nd.sigma, m * nd.n_grid)
                        acc_x, viol = ND.thinning_accept(D(0.5), sig_e, s_hat)
                        if viol == 1:
                            c_major = c_major + 1
                            t_unacc = t_unacc + wp.float64(energy)
                            code = code_nuclear
                            alive = 0
                        if alive == 1 and reason == 4 and axis2 < 0 and exited == 0:
                            if energy > ecut:
                                # candidate at the post-step point (leg 2 not truncated)
                                ub0 = nuc_u(h, gid_u, nc, 0, key)
                                ub1 = nuc_u(h, gid_u, nc, 1, key)
                                ub2 = nuc_u(h, gid_u, nc, 2, key)
                                nc = nc + 1
                                accepted_c, viol_c = ND.thinning_accept(D(ub0), sig_e, s_hat)
                                if viol_c == 1:  # fail closed: the majorant was violated
                                    c_major = c_major + 1
                                    t_unacc = t_unacc + wp.float64(energy)
                                    code = code_nuclear
                                    alive = 0
                                elif accepted_c == 0:  # fictitious: resample the optical depth
                                    n_lam = -wp.log(D(ub1))
                                else:
                                    alive = 0
                                    code = code_nuclear
                                    t1 = R(energy)
                                    j_t = NU.choose_target(
                                        ub2, t1, nd.grid, nd.n_grid, nd.sigma, nd.cum_sigma, m,
                                        nd.kmax, nd.mat_ntargets[m],
                                    )  # fmt: skip
                                    tgt = int(nd.mat_target[m * nd.kmax + j_t])
                                    status = NU.sample_event(
                                        h, gid_u, nc, key, tgt, t1, nd.grid, nd.n_grid, nd.lam,
                                        nd.edges, nd.rpre, nd.recoil, nd.tconst, nd.m_res,
                                        nd.evi, nd.evf, nd.prod, tid,
                                    )  # fmt: skip
                                    ib_e = tid * 16
                                    fb_e = tid * 8
                                    if status != 1:  # 64 attempts exhausted (or > 32 products)
                                        c_rejlim = c_rejlim + 1
                                        t_unacc = t_unacc + wp.float64(energy)
                                    else:
                                        tb_e = tgt * 16
                                        m_p_e = D(nd.tconst[tb_e])
                                        m_t_e = D(nd.tconst[tb_e + 1])
                                        n_prod = nd.evi[ib_e + 8]
                                        if with_diag:
                                            if h < ctl.trace_k and n_ev_rows < 2:
                                                ti = int(h - ctl.trace_h0)
                                                nd.ev_tr[ti, n_ev_rows, 0] = D(h)
                                                nd.ev_tr[ti, n_ev_rows, 1] = D(gid)
                                                nd.ev_tr[ti, n_ev_rows, 2] = D(tgt)
                                                for sp_c in range(5):
                                                    nd.ev_tr[ti, n_ev_rows, 3 + sp_c] = D(
                                                        nd.evi[ib_e + 1 + sp_c]
                                                    )
                                                nd.ev_tr[ti, n_ev_rows, 8] = D(nd.evi[ib_e + 6])
                                                nd.ev_tr[ti, n_ev_rows, 9] = D(nd.evi[ib_e + 7])
                                                nd.ev_tr[ti, n_ev_rows, 10] = D(nd.evi[ib_e])
                                                nd.ev_tr[ti, n_ev_rows, 11] = D(energy)
                                                n_ev_rows = n_ev_rows + 1
                                                nd.ev_n[ti] = n_ev_rows
                                        t_sum = D(0.0)
                                        alpha_t = D(0.0)
                                        sum_elab = D(0.0)
                                        m_out = D(0.0)
                                        n_children = int(0)
                                        # frame of the parent direction: e1, e2 perpendicular,
                                        # the event frame has z along the beam
                                        e1f, e2f = F.orthonormal_basis(v3(ux, uy, uz))
                                        gen_c = gen + 1
                                        for jp in range(n_prod):
                                            pb = (tid * 32 + jp) * 8
                                            si = int(nd.prod[pb])
                                            e_lab = D(nd.prod[pb + 4])
                                            mass_p = D(0.0)
                                            if si < 4:
                                                mass_p = D(nd.tconst[tb_e + 4 + si])
                                            t_lab = e_lab - mass_p
                                            t_sum = t_sum + t_lab
                                            sum_elab = sum_elab + e_lab
                                            if si == 0:
                                                n_neutron = n_neutron + t_lab
                                            elif si == 4:
                                                n_gamma = n_gamma + t_lab
                                            elif si == 3:
                                                alpha_t = alpha_t + t_lab
                                            else:  # p (1) or d (2): a secondary of generation + 1
                                                sp_new = si - 1
                                                cut_c = D(ctl.e_cut)
                                                mass_c = D(ctl.mass)
                                                if sp_new == 1:
                                                    cut_c = nd.e_cut_d
                                                    mass_c = nd.mass_d
                                                if t_lab < cut_c:
                                                    # below the species cutoff: deposited locally,
                                                    # tallied as cutoff
                                                    t_cutoff = t_cutoff + t_lab
                                                    cl = chan
                                                    cl.species = sp_new
                                                    cl.gen = gen_c
                                                    deposit_point(
                                                        edep, tally_rows, tid, batch, px, py, pz,
                                                        t_lab, g_origin, g_inv, g_shape, g_off,
                                                        ctl.n_grids, cl,
                                                    )  # fmt: skip
                                                else:
                                                    n_children = n_children + 1
                                                    gen_ok = int(1)
                                                    mult = int(1)
                                                    for _g in range(gen):
                                                        mult = mult * 32
                                                    cid = gid + n_children * mult
                                                    if n_children > nd.child_limit or gen > 6:
                                                        gen_ok = 0
                                                    if cid >= 1073741824:
                                                        gen_ok = 0
                                                    if gen_ok == 0:
                                                        c_gen = c_gen + 1
                                                        t_unacc = t_unacc + t_lab
                                                    elif n_st >= nd.stack_cap:
                                                        c_queue = c_queue + 1
                                                        t_unacc = t_unacc + t_lab
                                                    else:
                                                        qx = D(nd.prod[pb + 5])
                                                        qy = D(nd.prod[pb + 6])
                                                        qz = D(nd.prod[pb + 7])
                                                        pm = wp.sqrt(qx * qx + qy * qy + qz * qz)
                                                        dxs = ux
                                                        dys = uy
                                                        dzs = uz
                                                        if pm > D(0.0):
                                                            c0 = R(qx / pm)
                                                            c1 = R(qy / pm)
                                                            c2 = R(qz / pm)
                                                            dxs = (
                                                                c0 * e1f[0] + c1 * e2f[0] + c2 * ux
                                                            )
                                                            dys = (
                                                                c0 * e1f[1] + c1 * e2f[1] + c2 * uy
                                                            )
                                                            dzs = (
                                                                c0 * e1f[2] + c1 * e2f[2] + c2 * uz
                                                            )
                                                        pv_s = K.pv_mev(t_lab, mass_c)
                                                        nd.stack[tid, n_st, 0] = D(px)
                                                        nd.stack[tid, n_st, 1] = D(py)
                                                        nd.stack[tid, n_st, 2] = D(pz)
                                                        nd.stack[tid, n_st, 3] = D(dxs)
                                                        nd.stack[tid, n_st, 4] = D(dys)
                                                        nd.stack[tid, n_st, 5] = D(dzs)
                                                        nd.stack[tid, n_st, 6] = t_lab
                                                        nd.stack[tid, n_st, 7] = D(sp_new)
                                                        nd.stack[tid, n_st, 8] = D(cid)
                                                        nd.stack[tid, n_st, 9] = D(gen_c)
                                                        nd.stack[tid, n_st, 10] = D(ix)
                                                        nd.stack[tid, n_st, 11] = D(iy)
                                                        nd.stack[tid, n_st, 12] = D(iz)
                                                        nd.stack[tid, n_st, 13] = pv_s
                                                        n_st = n_st + 1
                                                        if with_diag:
                                                            if h < ctl.trace_k and n_sec_rows < 32:
                                                                ti2 = int(h - ctl.trace_h0)
                                                                nd.sec_tr[ti2, n_sec_rows, 0] = D(h)
                                                                nd.sec_tr[ti2, n_sec_rows, 1] = D(
                                                                    cid
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 2] = D(
                                                                    gid
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 3] = D(
                                                                    sp_new
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 4] = D(
                                                                    gen_c
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 5] = (
                                                                    t_lab
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 6] = D(
                                                                    dxs
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 7] = D(
                                                                    dys
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 8] = D(
                                                                    dzs
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 9] = D(
                                                                    px
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 10] = D(
                                                                    py
                                                                )
                                                                nd.sec_tr[ti2, n_sec_rows, 11] = D(
                                                                    pz
                                                                )
                                                                n_sec_rows = n_sec_rows + 1
                                                                nd.sec_n[ti2] = n_sec_rows
                                        for sp_m in range(5):
                                            mass_m = D(0.0)
                                            if sp_m < 4:
                                                mass_m = D(nd.tconst[tb_e + 4 + sp_m])
                                            m_out = m_out + D(nd.evi[ib_e + 1 + sp_m]) * mass_m
                                        t_r = D(nd.evf[fb_e])
                                        big_m = D(nd.evf[fb_e + 5])
                                        binding = D(nd.evf[fb_e + 2])
                                        imbal = D(nd.evf[fb_e + 1])
                                        if is_f32 == 1:
                                            # float64 ledger of the float32 sampler's products
                                            binding = m_out + big_m - m_p_e - m_t_e
                                            imbal = (
                                                D(energy) + m_p_e + m_t_e - sum_elab - big_m - t_r
                                            )
                                        local = alpha_t + t_r
                                        n_local = n_local + local
                                        n_alpha = n_alpha + alpha_t
                                        n_binding = n_binding + binding
                                        n_imbal = n_imbal + imbal
                                        cl2 = chan
                                        cl2.species = pseudo_local
                                        cl2.gen = gen + 1
                                        deposit_point(
                                            edep, tally_rows, tid, batch, px, py, pz, local,
                                            g_origin, g_inv, g_shape, g_off, ctl.n_grids, cl2,
                                        )  # fmt: skip
                                        # per-event ledger: T1 = sum T_lab + T_r + binding + imbal
                                        miss = energy - (t_sum + t_r + binding + imbal)
                                        if wp.abs(miss) > nd.ledger_tol * energy:
                                            c_cons = c_cons + 1
                                            t_unacc = t_unacc + miss
                    if alive == 1 and exited == 1:
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

            # ---- end of the particle: the primary's end record, then pop the stack ----
            if gen == 0:
                prim_done = 1
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
            cur_valid = 0
            if n_st > 0:
                n_st = n_st - 1
                px = R(nd.stack[tid, n_st, 0])
                py = R(nd.stack[tid, n_st, 1])
                pz = R(nd.stack[tid, n_st, 2])
                ux = R(nd.stack[tid, n_st, 3])
                uy = R(nd.stack[tid, n_st, 4])
                uz = R(nd.stack[tid, n_st, 5])
                energy = nd.stack[tid, n_st, 6]
                species = int(nd.stack[tid, n_st, 7])
                gid = int(nd.stack[tid, n_st, 8])
                gen = int(nd.stack[tid, n_st, 9])
                ix = int(nd.stack[tid, n_st, 10])
                iy = int(nd.stack[tid, n_st, 11])
                iz = int(nd.stack[tid, n_st, 12])
                p1v1 = nd.stack[tid, n_st, 13]
                cur_valid = 1

        if with_diag:
            if prim_done == 0:  # the primary never entered the step loop
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
        tally_rows[tid, 0] = t_initial
        tally_rows[tid, 1] = t_cutoff
        tally_rows[tid, 2] = t_step
        tally_rows[tid, 3] = t_escaped
        tally_rows[tid, 4] = t_truncated
        tally_rows[tid, 5] = t_unacc
        nb = nd.nt_base
        tally_rows[tid, nb] = n_local
        tally_rows[tid, nb + 1] = n_alpha
        tally_rows[tid, nb + 2] = n_neutron
        tally_rows[tid, nb + 3] = n_gamma
        tally_rows[tid, nb + 4] = n_binding
        tally_rows[tid, nb + 5] = n_imbal
        if chan.n_ch > 0 and path_flag == 1:
            tally_rows[tid, chan.res_base + chan.n_res + 1] = D(1.0)
        counter_rows[tid, 0] = c_trunc
        counter_rows[tid, 1] = c_stall
        counter_rows[tid, 2] = c_strag
        counter_rows[tid, 3] = c_gen
        counter_rows[tid, 4] = c_queue
        counter_rows[tid, 5] = c_src
        counter_rows[tid, 6] = c_inv
        counter_rows[tid, 8] = c_pieces
        counter_rows[tid, 9] = c_major
        counter_rows[tid, 10] = c_rejlim
        counter_rows[tid, 11] = c_cons

    transport.__name__ = f"transport_nuclear_{name}_{'diag' if diag else 'plain'}"
    transport.__qualname__ = transport.__name__
    return wp.kernel(enable_backward=False, module="unique")(transport)
