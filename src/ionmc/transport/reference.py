"""Reference backend: the step algorithm of decision 0039 in a Python history loop.

Every physical decision is made by the shared Warp functions called from Python scope with
``wp.float64`` arguments (``ionmc.physics`` and ``ionmc.transport.funcs``); this module only
contains glue: table reads from numpy arrays, state bookkeeping, random-number draws from the
Python Philox, scoring and tallies. The Warp backends (V3-003B) repeat this glue in a kernel,
and the per-step trace produced here is the reference for trajectory-level parity.

Draw order (counter ``(history, genealogy 0, block, purpose)``, key from the seed):

* source sampling, purpose 1, block 0: normals ``(u0, u1)`` give the lateral offsets along
  the deterministic orthonormal basis of the direction, ``(u2, u3)`` the energy offset (the
  block is always drawn; the values are ignored when the corresponding sigma is zero);
* transport, purpose 0: the block counter starts at 0 and increases by one per Philox call;
  each step draws block A (hinge fraction, polar-angle uniform, azimuth, reserved), then
  block B (straggling attempts, at most 64; exactly one ignored block when straggling is
  off). A step that finds ``E <= E_cut`` draws nothing.

Tallies (MeV, accumulated independently of the grids): ``initial``, ``cutoff`` (local
deposition below ``E_cut``; scored into the grids like any deposit), ``step_deposit``,
``escaped``, ``truncated`` (energy of histories stopped by the step limit or a stall, never
scored), ``unaccounted`` (genealogy overflow; always 0 here) and ``outside[g]`` (deposits
outside grid ``g``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import warp as wp
from numpy.typing import NDArray

from ionmc.config import MAX_REJECTION_ATTEMPTS, EffectiveConfig
from ionmc.physics.em import make_em
from ionmc.physics.kinematics import make_kinematics
from ionmc.rng.philox import (
    PURPOSE_SOURCE,
    PURPOSE_TRANSPORT,
    draw_block,
    key_from_seed,
    u01_py,
)
from ionmc.transport.funcs import make_transport_funcs

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
``select_step`` (0 geometry, 1 energy loss, 2 range, 3 maximum step))."""

COUNTER_NAMES = (
    "step_truncation",
    "stall",
    "straggling_rejection",
    "genealogy_overflow",
    "queue_overflow",
    "source_energy_out_of_range",
    "energy_inversion",
)


@dataclass
class RawTransport:
    """Raw output of the reference transport (float64)."""

    edep_mev: list[NDArray[np.float64]]
    tallies: dict[str, float]
    outside_mev: list[float]
    counters: dict[str, int]
    diagnostics: dict[str, Any] = field(default_factory=dict)


def run_reference(eff: EffectiveConfig) -> RawTransport:
    """Transport ``n_histories`` primaries with the Python reference loop (float64)."""
    return _Reference(eff).run()


def f_short_t(tt: float, r0: float, f_short: object) -> bool:
    """True for a step in the inverse-range branch (not the short linear one, not end of range)."""
    return float(f_short) * r0 <= tt < r0  # type: ignore[arg-type]


class _Reference:
    def __init__(self, eff: EffectiveConfig) -> None:
        cfg = eff.requested
        self.eff = eff
        self.cfg = cfg
        self.tab = eff.tables
        self.geo = eff.geometry
        self.F = make_transport_funcs(wp.float64)
        self.EM = make_em(wp.float64)
        self.K = make_kinematics(wp.float64)
        self.R = wp.float64
        self.V = self.F.vec3
        self.key = key_from_seed(cfg.run.seed)
        self.n_batches = cfg.run.n_batches
        ph = cfg.physics
        r = self.R
        self.e_cut = ph.e_cut_mev
        self.c_alpha = r(ph.range_alpha)
        self.c_rho_f = r(ph.range_rho_f_mm)
        self.c_frac = r(ph.max_energy_loss_fraction)
        self.c_smax = r(ph.max_step_mm)
        self.c_fshort = r(ph.short_step_fraction)
        self.mass = r(cfg.source.projectile.mass_mev)
        self.one = r(1.0)
        self.mat = self.geo.material_index
        self.dens = self.geo.densities_g_cm3()
        self.shape = self.geo.shape
        self.origin = self.V(*(r(x) for x in self.geo.origin_mm))
        self.spacing = self.V(*(r(x) for x in self.geo.spacing_mm))
        self.lo = self.V(*(r(x) for x in self.geo.lower_mm))
        self.hi = self.V(*(r(x) for x in self.geo.upper_mm))
        self.grids = cfg.scoring
        self.g_origin = [self.V(*(r(x) for x in g.origin_mm)) for g in self.grids]
        self.g_inv = [self.V(*(r(1.0 / x) for x in g.spacing_mm)) for g in self.grids]
        self.edep = [np.zeros((self.n_batches, g.n_voxels), dtype=np.float64) for g in self.grids]
        self.outside = [0.0] * len(self.grids)
        self.tallies = {
            k: 0.0
            for k in ("initial", "cutoff", "step_deposit", "escaped", "truncated", "unaccounted")
        }
        self.counters = dict.fromkeys(COUNTER_NAMES, 0)
        self.e_table_max = float(self.tab.e_max_mev.min())

    # -- table reads (shared bin location and interpolation, numpy memory access) -------------
    def _stopping_range(self, m: int, e: float) -> tuple[float, float]:
        t, r = self.tab, self.R
        i, f = self.F.log_bin_index(r(e), r(t.ln_e0[m]), r(t.inv_dln_e[m]), t.n_e)
        ls, lr = t.ln_s_mass[m], t.ln_r_mass[m]
        s = self.F.interp_exp(r(ls[i]), r(ls[i + 1]), f)
        rng = self.F.interp_exp(r(lr[i]), r(lr[i + 1]), f)
        return float(s), float(rng)

    def _energy_from_range(self, m: int, r_g_cm2: float) -> float:
        t, r = self.tab, self.R
        i, f = self.F.log_bin_index(r(r_g_cm2), r(t.ln_r0[m]), r(t.inv_dln_r[m]), t.n_r)
        row = t.ln_e_of_r[m]
        return float(self.F.interp_exp(r(row[i]), r(row[i + 1]), f))

    # -- scoring ------------------------------------------------------------------------------
    def _score(self, batch: int, pos: tuple[float, float, float], deposit: float) -> None:
        r = self.R
        p = self.V(r(pos[0]), r(pos[1]), r(pos[2]))
        for g, grid in enumerate(self.grids):
            nx, ny, nz = grid.shape
            ix, iy, iz, inside = self.F.grid_index(p, self.g_origin[g], self.g_inv[g], nx, ny, nz)
            if inside:
                self.edep[g][batch, (ix * ny + iy) * nz + iz] += deposit
            else:
                self.outside[g] += deposit

    # -- driver -------------------------------------------------------------------------------
    def run(self) -> RawTransport:
        cfg = self.cfg
        n = cfg.run.n_histories
        diag = cfg.diagnostics
        self.end_pos = np.full((n, 3), np.nan) if diag.track_end_positions else None
        self.end_code = np.full(n, -1, dtype=np.int8) if diag.track_end_positions else None
        self.end_energy = np.full(n, np.nan) if diag.track_end_positions else None
        self.escapes: list[tuple[int, float, float, float, float, float, float, float]] = []
        self.trace: list[list[float]] = []
        self.trace_end: list[tuple[int, int, float]] = []
        for h in range(n):
            self._history(h)
        diagnostics: dict[str, Any] = {}
        if diag.track_end_positions:
            diagnostics["end_position_mm"] = self.end_pos
            diagnostics["end_code"] = self.end_code
            diagnostics["end_energy_mev"] = self.end_energy
        if diag.escape_records:
            arr = np.array(self.escapes, dtype=np.float64).reshape(-1, 8)
            diagnostics["escape_history"] = arr[:, 0].astype(np.int64)
            diagnostics["escape_position_mm"] = arr[:, 1:4]
            diagnostics["escape_direction"] = arr[:, 4:7]
            diagnostics["escape_energy_mev"] = arr[:, 7]
        if diag.trace_histories > 0:
            tr = np.array(self.trace, dtype=np.float64).reshape(-1, len(TRACE_COLUMNS))
            diagnostics["trace"] = {name: tr[:, i] for i, name in enumerate(TRACE_COLUMNS)}
            diagnostics["trace_columns"] = TRACE_COLUMNS
            te = np.array(self.trace_end, dtype=np.float64).reshape(-1, 3)
            diagnostics["trace_end_history"] = te[:, 0].astype(np.int64)
            diagnostics["trace_end_code"] = te[:, 1].astype(np.int64)
            diagnostics["trace_end_energy_mev"] = te[:, 2]
        return RawTransport(self.edep, self.tallies, self.outside, self.counters, diagnostics)

    def _end(self, h: int, code: int, pos: tuple[float, float, float], energy: float) -> None:
        if self.end_pos is not None and self.end_code is not None and self.end_energy is not None:
            self.end_pos[h] = pos
            self.end_code[h] = code
            self.end_energy[h] = energy
        if h < self.cfg.diagnostics.trace_histories:
            self.trace_end.append((h, code, energy))

    # -- one history --------------------------------------------------------------------------
    def _history(self, h: int) -> None:
        cfg, F, EM, K, V, r = self.cfg, self.F, self.EM, self.K, self.V, self.R
        src, ph = cfg.source, cfg.physics
        key = self.key
        batch = h % self.n_batches
        nx, ny, nz = self.shape
        trace_this = h < cfg.diagnostics.trace_histories

        # source sampling (purpose 1)
        w = draw_block(key, h, 0, 0, PURPOSE_SOURCE)
        u = [u01_py(x, "float64") for x in w]
        z1, z2 = F.gauss_pair(r(u[0]), r(u[1]))
        z3 = F.gauss_one(r(u[2]), r(u[3]))
        d = self.eff.unit_direction
        e1v, e2v = F.orthonormal_basis(V(r(d[0]), r(d[1]), r(d[2])))
        energy = src.kinetic_energy_mev
        if src.energy_sigma_mev > 0.0:
            energy = energy + src.energy_sigma_mev * float(z3)
        px, py, pz = src.position_mm
        if src.lateral_sigma_mm > 0.0:
            sg = src.lateral_sigma_mm
            zz1, zz2 = float(z1), float(z2)
            px = px + sg * (zz1 * float(e1v[0]) + zz2 * float(e2v[0]))
            py = py + sg * (zz1 * float(e1v[1]) + zz2 * float(e2v[1]))
            pz = pz + sg * (zz1 * float(e1v[2]) + zz2 * float(e2v[2]))
        if not (self.e_cut <= energy <= self.e_table_max):
            self.counters["source_energy_out_of_range"] += 1
            self._end(h, END_SOURCE_REJECTED, (px, py, pz), energy)
            return
        self.tallies["initial"] += energy
        p1v1 = float(K.pv_mev(r(energy), self.mass))
        ux, uy, uz = d

        # entry into the world (vacuum outside)
        t_in, t_out = F.ray_box(V(r(px), r(py), r(pz)), V(r(ux), r(uy), r(uz)), self.lo, self.hi)
        t_in, t_out = float(t_in), float(t_out)
        if t_in > t_out or t_out < 0.0:
            self.tallies["escaped"] += energy
            self._end(h, END_MISSED_WORLD, (px, py, pz), energy)
            return
        t0 = max(t_in, 0.0)
        lo, hi = self.geo.lower_mm, self.geo.upper_mm
        px = min(max(px + ux * t0, lo[0]), hi[0])
        py = min(max(py + uy * t0, lo[1]), hi[1])
        pz = min(max(pz + uz * t0, lo[2]), hi[2])
        o, sp = self.geo.origin_mm, self.geo.spacing_mm
        ix = min(max(int(math.floor((px - o[0]) / sp[0])), 0), nx - 1)
        iy = min(max(int(math.floor((py - o[1]) / sp[1])), 0), ny - 1)
        iz = min(max(int(math.floor((pz - o[2]) / sp[2])), 0), nz - 1)

        steps = 0
        birth = True  # first step of the particle's life (exact log average of f_dM)
        blocks = 0
        zero_run = 0
        max_steps = self.eff.max_steps
        while True:
            if energy <= self.e_cut:
                self.tallies["cutoff"] += energy
                self._score(batch, (px, py, pz), energy)
                self._end(h, END_CUTOFF, (px, py, pz), energy)
                return
            if steps >= max_steps:
                self.tallies["truncated"] += energy
                self.counters["step_truncation"] += 1
                self._end(h, END_TRUNCATED, (px, py, pz), energy)
                return

            m = int(self.mat[ix, iy, iz])
            rho = float(self.dens[ix, iy, iz])
            s0, r0 = self._stopping_range(m, energy)
            s_lin = s0 * rho / 10.0
            r_mm = r0 * 10.0 / rho
            pvec = V(r(px), r(py), r(pz))
            dvec = V(r(ux), r(uy), r(uz))
            d_geo, _axis = F.dda_next(pvec, dvec, ix, iy, iz, self.origin, self.spacing)
            s_el = F.eloss_step_limit(r(energy), r(s_lin), self.c_frac)
            s_rg = F.range_step_limit(r(r_mm), self.c_alpha, self.c_rho_f)
            s_w, reason = F.select_step(d_geo, s_el, s_rg, self.c_smax)
            s = float(s_w)

            # block A
            wa = draw_block(key, h, 0, blocks, PURPOSE_TRANSPORT)
            blocks += 1
            ua = [u01_py(x, "float64") for x in wa]
            leg1 = ua[0] * s

            # multiple scattering
            d1x, d1y, d1z = ux, uy, uz
            if ph.multiple_scattering:
                e_mid = self._energy_from_range(m, r0 - rho * s / 20.0)
                inv_xs = r(self.tab.inv_rho_xs_cm2_g[m])
                pv_mid = K.pv_mev(r(e_mid), self.mass)
                if birth:
                    e_end = min(self._energy_from_range(m, r0 - rho * s / 10.0), energy)
                    var = float(
                        EM.scattering_variance_birth(
                            pv_mid,
                            K.pv_mev(r(e_end), self.mass),
                            r(p1v1),
                            self.one,
                            inv_xs,
                            r(rho),
                            r(s),
                        )
                    )
                else:
                    t_pow = EM.scattering_power_dm(pv_mid, r(p1v1), self.one, inv_xs, r(rho))
                    var = float(t_pow) * s
                theta = EM.polar_deflection(r(var), r(ua[1]))
                nd = EM.rotate_dir(dvec, theta, r(2.0 * math.pi * ua[2]))
                d1x, d1y, d1z = float(nd[0]), float(nd[1]), float(nd[2])

            # hinge and second leg
            hx, hy, hz = px + ux * leg1, py + uy * leg1, pz + uz * leg1
            leg2_w, axis2 = F.leg2_limit(
                V(r(hx), r(hy), r(hz)),
                V(r(d1x), r(d1y), r(d1z)),
                ix,
                iy,
                iz,
                self.origin,
                self.spacing,
                r(s - leg1),
            )
            leg2 = float(leg2_w)
            nxp, nyp, nzp = hx + d1x * leg2, hy + d1y * leg2, hz + d1z * leg2
            exited = False
            if axis2 >= 0:
                d1a = (d1x, d1y, d1z)[axis2]
                upward = 1 if d1a > 0.0 else 0
                idx = [ix, iy, iz]
                plane = float(F.plane_position(idx[axis2], upward, r(o[axis2]), r(sp[axis2])))
                idx[axis2] += 1 if upward else -1
                if axis2 == 0:
                    nxp = plane
                elif axis2 == 1:
                    nyp = plane
                else:
                    nzp = plane
                ix, iy, iz = idx
                exited = not (0 <= ix < nx and 0 <= iy < ny and 0 <= iz < nz)
            s_act = leg1 + leg2

            # energy loss
            tt = rho * s_act / 10.0
            e_r1 = self._energy_from_range(m, r0 - tt)
            if f_short_t(tt, r0, self.c_fshort) and e_r1 > energy:
                self.counters["energy_inversion"] += 1  # E1 > E0 from the inverse round trip
                e_r1 = energy
            mean = EM.csda_mean_loss(r(energy), r(e_r1), r(s0), r(tt), r(r0), self.c_fshort)
            mean_f = float(mean)
            attempts = 1
            if ph.straggling:
                var_e = EM.bohr_variance(
                    r(energy - 0.5 * mean_f),
                    self.mass,
                    self.one,
                    r(self.tab.z_over_a[m]),
                    r(rho),
                    r(s_act),
                )
                loss = mean_f
                accepted = False
                for k in range(MAX_REJECTION_ATTEMPTS):
                    wb = draw_block(key, h, 0, blocks, PURPOSE_TRANSPORT)
                    blocks += 1
                    ub = [u01_py(x, "float64") for x in wb]
                    lw, ok = EM.straggle_attempt(
                        mean, var_e, r(ub[0]), r(ub[1]), r(ub[2]), r(ub[3])
                    )
                    attempts = k + 1
                    if ok:
                        loss = float(lw)
                        accepted = True
                        break
                if not accepted:
                    self.counters["straggling_rejection"] += 1
                    loss = mean_f
            else:
                draw_block(key, h, 0, blocks, PURPOSE_TRANSPORT)
                blocks += 1
                loss = mean_f
            loss = min(loss, energy)
            e_new = energy - loss
            deposit = energy - e_new

            # scoring at the hinge-path midpoint
            if deposit > 0.0:
                mid = F.point_on_hinge(
                    pvec, dvec, r(leg1), V(r(d1x), r(d1y), r(d1z)), r(0.5 * s_act)
                )
                self._score(batch, (float(mid[0]), float(mid[1]), float(mid[2])), deposit)
                self.tallies["step_deposit"] += deposit

            px, py, pz = nxp, nyp, nzp
            ux, uy, uz = d1x, d1y, d1z
            energy = e_new
            steps += 1
            if s_act > 0.0:
                birth = False
            if trace_this:
                self.trace.append(
                    [
                        h,
                        steps,
                        ix,
                        iy,
                        iz,
                        int(reason),
                        blocks,
                        attempts,
                        px,
                        py,
                        pz,
                        ux,
                        uy,
                        uz,
                        energy,
                        deposit,
                        s_act,
                    ]
                )
            if exited:
                self.tallies["escaped"] += energy
                if self.cfg.diagnostics.escape_records:
                    self.escapes.append((h, px, py, pz, ux, uy, uz, energy))
                self._end(h, END_ESCAPED, (px, py, pz), energy)
                return
            zero_run = zero_run + 1 if s_act <= 0.0 else 0
            if zero_run > 3:
                self.tallies["truncated"] += energy
                self.counters["stall"] += 1
                self._end(h, END_TRUNCATED, (px, py, pz), energy)
                return
