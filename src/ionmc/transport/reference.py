"""Reference backend: the step algorithm of decision 0039 in a Python history loop.

Every physical decision is made by the shared functions of ``ionmc.physics`` and
``ionmc.transport.funcs`` in their pure-Python float64 twins (``ionmc._wpfunc.python_twin``:
the same source text re-executed without Warp; no Warp call of any kind is made at Python scope,
because Warp 1.17 Python-scope dispatch segfaulted intermittently in worker processes);
this module only
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

Scoring (decision amended in V3-003B): the deposit of a step (after straggling) is distributed
along both hinge legs for every scoring grid by an incremental DDA over the leg (``seg_piece``):
a piece of the path inside one voxel receives the integral over the piece of a linear
stopping-power ramp from S(E_old) to S(E_new), normalised over the step (``ramp_weight``), so the
pieces sum to the deposit and the residual step dependence is second order (the curvature of S
along the step); the energy left at the cutoff is a point deposit at the end point. Pieces are
quantized to int64 quanta (``tally``).

Tallies (MeV, accumulated independently of the grids): ``initial``, ``cutoff`` (local
deposition below ``E_cut``; scored into the grids like any deposit), ``step_deposit``,
``escaped``, ``truncated`` (energy of histories stopped by the step limit or a stall, never
scored), ``unaccounted`` (genealogy overflow; always 0 here) and ``outside[g]`` (deposits
outside grid ``g``).
"""

from __future__ import annotations

import math

import numpy as np

from ionmc._wpfunc import python_twin
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
from ionmc.species import species_of_projectile
from ionmc.transport.channels import CLASS_LOCAL, CLASS_STEP
from ionmc.transport.funcs import BIG_LENGTH_MM, make_transport_funcs
from ionmc.transport.run import channel_columns
from ionmc.transport.scoring_ref import ReferenceChannelScorer
from ionmc.transport.tally import (
    COUNTER_NAMES,
    END_CUTOFF,
    END_ESCAPED,
    END_MISSED_WORLD,
    END_SOURCE_REJECTED,
    END_TRUNCATED,
    N_FIXED_TALLIES,
    QUANTUM_MEV,
    QUANTUM_SCALE,
    TALLY_NAMES,
    TRACE_COLUMNS,
    HistoryDiagnostics,
    PartialTransport,
    RawTransport,
    build_diagnostics,
    merge_partials,
    rows_to_partial,
)

__all__ = [
    "COUNTER_NAMES",
    "END_CUTOFF",
    "END_ESCAPED",
    "END_MISSED_WORLD",
    "END_SOURCE_REJECTED",
    "END_TRUNCATED",
    "TRACE_COLUMNS",
    "RawTransport",
    "run_reference",
    "run_reference_range",
]


def ramp_weight(ta: float, tb: float, s_start: float, s_end: float, s_act: float) -> float:
    """Fraction of a step's energy deposited between the path coordinates ``ta`` and ``tb`` when the
    deposit density follows the linear stopping-power ramp ``s_start -> s_end`` along the step.
    (The same expression, in the same order, as ``ionmc.transport.kernels``.)"""
    return (s_start * (tb - ta) + (s_end - s_start) * (tb * tb - ta * ta) / (2.0 * s_act)) / (
        0.5 * (s_start + s_end) * s_act
    )


def run_reference(eff: EffectiveConfig) -> RawTransport:
    """Transport ``n_histories`` primaries with the Python reference loop (float64)."""
    n = eff.requested.run.n_histories
    part = run_reference_range(eff, 0, n)
    diag = eff.requested.diagnostics
    raw = merge_partials([part], n, len(eff.requested.scoring), channel_columns(eff))
    raw.diagnostics = build_diagnostics(
        [part], diag.track_end_positions, diag.escape_records, diag.trace_histories
    )
    raw.meta = {"workers": 1}
    return raw


def run_reference_range(eff: EffectiveConfig, h0: int, h1: int) -> PartialTransport:
    """Transport the histories ``[h0, h1)``; the result depends on no other history."""
    return _Reference(eff).run_range(h0, h1)


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
        # pure-Python twins of the shared Warp functions (same source text, no Warp call)
        self.F = python_twin(make_transport_funcs)
        self.EM = python_twin(make_em)
        self.K = python_twin(make_kinematics)
        self.R = float
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
        self.path_bound_mm = math.inf if eff.channels is None else eff.channels.path_bound_mm
        self.path_exceeded = 0
        self.c_fshort = r(ph.short_step_fraction)
        self.mass = r(cfg.source.projectile.mass_mev)
        self.trunc_diag = ph.truncated_hinge_diagnostic
        self.straggle_attempt = (
            self.EM.straggle_attempt_gamma
            if ph.straggling_model == "bohr_gamma_v1"
            else self.EM.straggle_attempt
        )
        self.ctrl_res = 0.0
        self.ctrl_sum = [0.0, 0.0, 0.0]
        self.max_pieces = eff.scoring_pieces
        self.one = r(1.0)
        self.mat = self.geo.material_index
        self.dens = self.geo.densities_g_cm3()
        self.shape = self.geo.shape
        self.origin = self.V(*(r(x) for x in self.geo.origin_mm))
        self.spacing = self.V(*(r(x) for x in self.geo.spacing_mm))
        self.lo = self.V(*(r(x) for x in self.geo.lower_mm))
        self.hi = self.V(*(r(x) for x in self.geo.world_upper_mm))
        zc = self.geo.z_exit_mm
        self.z_clip = r(BIG_LENGTH_MM if zc is None else float(zc))
        self.grids = cfg.scoring
        self.g_origin = [self.V(*(r(x) for x in g.origin_mm)) for g in self.grids]
        self.g_inv = [self.V(*(r(1.0 / x) for x in g.spacing_mm)) for g in self.grids]
        self.g_spacing = [self.V(*(r(x) for x in g.spacing_mm)) for g in self.grids]
        self.edep = [np.zeros((self.n_batches, g.n_voxels), dtype=np.int64) for g in self.grids]
        self.outside = [0.0] * len(self.grids)
        self.quant = [0.0] * len(self.grids)
        self.tallies = dict.fromkeys(TALLY_NAMES, 0.0)
        self.counters = dict.fromkeys(COUNTER_NAMES, 0)
        self.e_table_max = float(self.tab.e_max_mev.min())
        # scoring channels (decision 0040): None without tallies, then the qualified path is
        # untouched (no step precompute, no extra leg walks, no extra tally columns)
        self.scorer: ReferenceChannelScorer | None = None
        if eff.channels is not None:
            self.scorer = ReferenceChannelScorer(
                eff.channels, self.grids, self.n_batches, tables=self.tab,
                a_nucleon=cfg.source.projectile.a,
            )  # fmt: skip
            self.species_id = species_of_projectile(cfg.source.projectile).id
            self.generation = 0  # primaries only until V3-005A

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
    # Deposits are quantized (nearest multiple of QUANTUM_MEV, floor(x / q + 1/2)) and added to
    # int64 grids; the rounding residual of every in-grid piece is tallied per grid, deposits
    # outside a grid are tallied in float64 (see ionmc.transport.tally).
    def _deposit_voxel(self, batch: int, g: int, ix: int, iy: int, iz: int, de: float) -> int:
        """Deposit ``de`` into voxel ``(ix, iy, iz)`` of grid ``g``; returns its flat index, or -1
        if the voxel is outside the grid (the deposit is then tallied as outside)."""
        nx, ny, nz = self.grids[g].shape
        if 0 <= ix < nx and 0 <= iy < ny and 0 <= iz < nz:
            n = math.floor(de * QUANTUM_SCALE + 0.5)
            vox = (ix * ny + iy) * nz + iz
            self.edep[g][batch, vox] += n
            self.quant[g] += de - n * QUANTUM_MEV
            return vox
        self.outside[g] += de
        return -1

    def _score_local(self, batch: int, g: int, vox: int, de: float) -> None:
        """Class "local" deposit (cutoff energy, deposit at ``s_act = 0``) into the channels."""
        if self.scorer is not None and vox >= 0:
            self.scorer.score_piece(
                batch, g, vox, 0.0, 0.0, de, self.species_id, self.generation, CLASS_LOCAL
            )

    def _deposit_point(self, batch: int, pos: tuple[float, float, float], de: float) -> None:
        """Point deposit (the energy left at the cutoff) in every grid."""
        r = self.R
        p = self.V(r(pos[0]), r(pos[1]), r(pos[2]))
        for g, grid in enumerate(self.grids):
            nx, ny, nz = grid.shape
            ix, iy, iz, _inside = self.F.grid_index(p, self.g_origin[g], self.g_inv[g], nx, ny, nz)
            self._score_local(batch, g, self._deposit_voxel(batch, g, ix, iy, iz, de), de)

    def _deposit_leg(
        self,
        batch: int,
        g: int,
        p: tuple[float, float, float],
        u: tuple[float, float, float],
        length: float,
        deposit: float,
        s_act: float,
        t_off: float,
        s_start: float,
        s_end: float,
    ) -> None:
        """Deposit the share of ``deposit`` that belongs to each piece in every voxel of grid ``g``
        crossed by the straight segment from ``p`` along ``u`` of the given length (per-grid
        incremental DDA). The leg starts at the path coordinate ``t_off`` of the step; the share
        of a piece ``[ta, tb]`` is the integral of the linear stopping-power ramp
        ``w(t) = s_start + (s_end - s_start) t / s_act`` over it, divided by the integral over the
        whole step (the pieces of both legs sum to the deposit)."""
        r, F, V = self.R, self.F, self.V
        nx, ny, nz = self.grids[g].shape
        go, gs = self.g_origin[g], self.g_spacing[g]
        og, sg = self.grids[g].origin_mm, self.grids[g].spacing_mm
        px, py, pz = p
        ux, uy, uz = u
        uvec = V(r(ux), r(uy), r(uz))
        ix, iy, iz, _inside = F.grid_index(V(r(px), r(py), r(pz)), go, self.g_inv[g], nx, ny, nz)
        remaining = length
        tcur = t_off
        for _ in range(self.max_pieces):
            if remaining > 0.0:
                piece_w, axis = F.seg_piece(V(r(px), r(py), r(pz)), uvec, ix, iy, iz, go, gs,
                                            r(remaining))  # fmt: skip
                piece = float(piece_w)
                tb = tcur + piece
                eps_p = deposit * ramp_weight(tcur, tb, s_start, s_end, s_act)
                vox = self._deposit_voxel(batch, g, ix, iy, iz, eps_p)
                if self.scorer is not None and vox >= 0:
                    self.scorer.score_piece(
                        batch, g, vox, piece, tcur + 0.5 * piece - 0.5 * s_act, eps_p,
                        self.species_id, self.generation, CLASS_STEP,
                    )  # fmt: skip
                tcur = tb
                remaining = remaining - piece
                if axis >= 0:
                    px = px + ux * piece
                    py = py + uy * piece
                    pz = pz + uz * piece
                    ua = (ux, uy, uz)[axis]
                    upward = 1 if ua > 0.0 else 0
                    idx = (ix, iy, iz)[axis]
                    plane = float(F.plane_position(idx, upward, r(og[axis]), r(sg[axis])))
                    step = 1 if upward else -1
                    if axis == 0:
                        px = plane
                        ix = ix + step
                    elif axis == 1:
                        py = plane
                        iy = iy + step
                    else:
                        pz = plane
                        iz = iz + step
        if remaining > 0.0:  # exceeds the validated piece bound: conserve energy, invalidate
            self.counters["scoring_pieces_overflow"] += 1
            eps_p = deposit * ramp_weight(tcur, tcur + remaining, s_start, s_end, s_act)
            vox = self._deposit_voxel(batch, g, ix, iy, iz, eps_p)
            if self.scorer is not None and vox >= 0:
                self.scorer.score_piece(
                    batch, g, vox, remaining, tcur + 0.5 * remaining - 0.5 * s_act, eps_p,
                    self.species_id, self.generation, CLASS_STEP,
                )  # fmt: skip

    def _deposit_step(
        self,
        batch: int,
        p0: tuple[float, float, float],
        d0: tuple[float, float, float],
        leg1: float,
        hinge: tuple[float, float, float],
        d1: tuple[float, float, float],
        leg2: float,
        deposit: float,
        s_act: float,
        s_start: float,
        s_end: float,
    ) -> None:
        """Track-length apportioning of a step deposit along both hinge legs in every grid, with
        the linear stopping-power ramp from ``s_start`` (energy at the step start) to ``s_end``
        (energy after the step)."""
        for g in range(len(self.grids)):
            if s_act > 0.0:
                self._deposit_leg(batch, g, p0, d0, leg1, deposit, s_act, 0.0, s_start, s_end)
                self._deposit_leg(batch, g, hinge, d1, leg2, deposit, s_act, leg1, s_start, s_end)
            else:
                self._deposit_point_grid(batch, g, p0, deposit)

    def _deposit_point_grid(
        self, batch: int, g: int, pos: tuple[float, float, float], de: float
    ) -> None:
        r = self.R
        nx, ny, nz = self.grids[g].shape
        ix, iy, iz, _inside = self.F.grid_index(
            self.V(r(pos[0]), r(pos[1]), r(pos[2])), self.g_origin[g], self.g_inv[g], nx, ny, nz
        )
        self._score_local(batch, g, self._deposit_voxel(batch, g, ix, iy, iz, de), de)

    # -- driver -------------------------------------------------------------------------------
    def run_range(self, h0: int, h1: int) -> PartialTransport:
        n = h1 - h0
        diag = self.cfg.diagnostics
        n_g = len(self.grids)
        self.want_diag = diag.track_end_positions or diag.escape_records or diag.trace_histories > 0
        self.end_pos = np.full((n, 3), np.nan)
        self.end_dir = np.full((n, 3), np.nan)
        self.end_code = np.full(n, -1, dtype=np.int8)
        self.end_energy = np.full(n, np.nan)
        self.end_ctrl = np.zeros(n)
        self.end_ctrl_sum = np.zeros((n, 3))
        self.trace: list[list[float]] = []
        self.h_base = h0
        n_chan_cols = channel_columns(self.eff)
        tally_rows = np.zeros((n, N_FIXED_TALLIES + 2 * n_g + n_chan_cols))
        counter_rows = np.zeros((n, len(COUNTER_NAMES)), dtype=np.int32)
        for h in range(h0, h1):
            # per-history accumulators: a row depends on this history alone
            self.tallies = dict.fromkeys(TALLY_NAMES, 0.0)
            self.outside = [0.0] * n_g
            self.quant = [0.0] * n_g
            self.counters = dict.fromkeys(COUNTER_NAMES, 0)
            self.ctrl_res = 0.0
            self.ctrl_sum = [0.0, 0.0, 0.0]
            self.path_exceeded = 0
            if self.scorer is not None:
                self.scorer.begin_history()
            self._history(h)
            row = h - h0
            tally_rows[row, :N_FIXED_TALLIES] = [self.tallies[k] for k in TALLY_NAMES]
            tally_rows[row, N_FIXED_TALLIES : N_FIXED_TALLIES + n_g] = self.outside
            tally_rows[row, N_FIXED_TALLIES + n_g : N_FIXED_TALLIES + 2 * n_g] = self.quant
            if self.scorer is not None:
                tally_rows[row, N_FIXED_TALLIES + 2 * n_g : -2] = self.scorer.residual
                tally_rows[row, -2] = self.scorer.lookup_ood
                tally_rows[row, -1] = self.path_exceeded
            counter_rows[row] = [self.counters[k] for k in COUNTER_NAMES]
        diagnostics = None
        if self.want_diag:
            tr = np.array(self.trace, dtype=np.float64).reshape(-1, len(TRACE_COLUMNS))
            diagnostics = HistoryDiagnostics(
                end_position_mm=self.end_pos,
                end_direction=self.end_dir,
                end_energy_mev=self.end_energy,
                end_code=self.end_code,
                control_residual=self.end_ctrl,
                control_displacement=self.end_ctrl_sum,
                trace_int=tr[:, :8].astype(np.int32),
                trace_float=tr[:, 8:],
            )
        return rows_to_partial(
            h0, h1, tally_rows, counter_rows, self.edep, diagnostics,
            channel_acc=None if self.scorer is None else self.scorer.acc,
        )  # fmt: skip

    def _end(
        self,
        h: int,
        code: int,
        pos: tuple[float, float, float],
        direction: tuple[float, float, float],
        energy: float,
    ) -> None:
        i = h - self.h_base
        self.end_pos[i] = pos
        self.end_dir[i] = direction
        self.end_code[i] = code
        self.end_energy[i] = energy
        self.end_ctrl[i] = self.ctrl_res
        self.end_ctrl_sum[i] = self.ctrl_sum

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
            self._end(h, END_SOURCE_REJECTED, (px, py, pz), d, energy)
            return
        self.tallies["initial"] += energy
        p1v1 = float(K.pv_mev(r(energy), self.mass))
        ux, uy, uz = d

        # entry into the world (vacuum outside)
        t_in, t_out = F.ray_box(V(r(px), r(py), r(pz)), V(r(ux), r(uy), r(uz)), self.lo, self.hi)
        t_in, t_out = float(t_in), float(t_out)
        if t_in > t_out or t_out < 0.0:
            self.tallies["escaped"] += energy
            self._end(h, END_MISSED_WORLD, (px, py, pz), (ux, uy, uz), energy)
            return
        t0 = max(t_in, 0.0)
        lo, hi = self.geo.lower_mm, self.geo.world_upper_mm
        px = min(max(px + ux * t0, lo[0]), hi[0])
        py = min(max(py + uy * t0, lo[1]), hi[1])
        pz = min(max(pz + uz * t0, lo[2]), hi[2])
        o, sp = self.geo.origin_mm, self.geo.spacing_mm
        ix = min(max(int(math.floor((px - o[0]) / sp[0])), 0), nx - 1)
        iy = min(max(int(math.floor((py - o[1]) / sp[1])), 0), ny - 1)
        iz = min(max(int(math.floor((pz - o[2]) / sp[2])), 0), nz - 1)

        steps = 0
        path_mm = 0.0
        # first step of the particle's life: linearized analytic log-average of f_dM (~1e-3)
        birth = True
        blocks = 0
        zero_run = 0
        max_steps = self.eff.max_steps
        while True:
            if energy <= self.e_cut:
                self.tallies["cutoff"] += energy
                self._deposit_point(batch, (px, py, pz), energy)
                self._end(h, END_CUTOFF, (px, py, pz), (ux, uy, uz), energy)
                return
            if steps >= max_steps:
                self.tallies["truncated"] += energy
                self.counters["step_truncation"] += 1
                self._end(h, END_TRUNCATED, (px, py, pz), (ux, uy, uz), energy)
                return

            m = int(self.mat[ix, iy, iz])
            rho = float(self.dens[ix, iy, iz])
            s0, r0 = self._stopping_range(m, energy)
            s_lin = s0 * rho / 10.0
            r_mm = r0 * 10.0 / rho
            pvec = V(r(px), r(py), r(pz))
            dvec = V(r(ux), r(uy), r(uz))
            d_geo, axis_pre = F.dda_next_clip(
                pvec, dvec, ix, iy, iz, self.origin, self.spacing, self.z_clip
            )
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
            if self.trunc_diag:
                # DIAGNOSTIC (T14 negative control, default off), "truncate-first": the planned
                # step already ends at the first plane the straight line (pre-hinge direction)
                # reaches (reason 0, axis_pre); the angle was sampled for that length, and the two
                # legs are travelled without cutting leg 2 again; the end point is then snapped
                # onto that plane (below) keeping the lateral displacement of the hinge path.
                leg2 = s - leg1
                axis2 = int(axis_pre) if int(reason) == 0 else -1
            else:
                leg2_w, axis2 = F.leg2_limit_clip(
                    V(r(hx), r(hy), r(hz)),
                    V(r(d1x), r(d1y), r(d1z)),
                    ix,
                    iy,
                    iz,
                    self.origin,
                    self.spacing,
                    r(s - leg1),
                    self.z_clip,
                )
                leg2 = float(leg2_w)
            nxp, nyp, nzp = hx + d1x * leg2, hy + d1y * leg2, hz + d1z * leg2
            exited = False
            if axis2 == 3:  # the clip plane z = z_exit: the particle leaves the world there
                if self.trunc_diag and s > 0.0:
                    self.ctrl_res = max(self.ctrl_res, abs(nzp - float(self.z_clip)) / s)
                    self.ctrl_sum[2] += float(self.z_clip) - nzp
                nzp = float(self.z_clip)
                exited = True
            elif axis2 >= 0:
                # the plane crossed: by the bent leg 2 normally, by the pre-hinge direction in the
                # truncate-first control
                d1a = (ux, uy, uz)[axis2] if self.trunc_diag else (d1x, d1y, d1z)[axis2]
                upward = 1 if d1a > 0.0 else 0
                idx = [ix, iy, iz]
                plane = float(F.plane_position(idx[axis2], upward, r(o[axis2]), r(sp[axis2])))
                if self.trunc_diag and s > 0.0:  # snap displacement relative to the step
                    got = (nxp, nyp, nzp)[axis2]
                    self.ctrl_res = max(self.ctrl_res, abs(got - plane) / s)
                    self.ctrl_sum[axis2] += plane - got
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
            s_branch = s0  # stopping power of the linear branch: at the midpoint energy
            if tt < float(self.c_fshort) * r0:
                e_half = energy - 0.5 * s0 * tt
                s_branch = self._stopping_range(m, e_half)[0]
            mean = EM.csda_mean_loss(r(energy), r(e_r1), r(s_branch), r(tt), r(r0), self.c_fshort)
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
                    lw, ok = self.straggle_attempt(
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

            # scoring: the deposit is apportioned along both legs by the linear stopping-power ramp
            # integrated over the path inside each voxel
            # (with scoring channels the legs of a step with s_act > 0 are walked even when
            # deposit = 0: fluence and LET do not depend on eps > 0)
            walk = deposit > 0.0 or (self.scorer is not None and s_act > 0.0)
            if walk:
                if self.scorer is not None and s_act > 0.0:
                    self.scorer.begin_step(energy - 0.5 * mean_f, mean_f, s_act)
                self._deposit_step(
                    batch, (px, py, pz), (ux, uy, uz), leg1, (hx, hy, hz), (d1x, d1y, d1z), leg2,
                    deposit, s_act, s0, self._stopping_range(m, e_new)[0],
                )  # fmt: skip
            if deposit > 0.0:
                self.tallies["step_deposit"] += deposit

            px, py, pz = nxp, nyp, nzp
            ux, uy, uz = d1x, d1y, d1z
            energy = e_new
            steps += 1
            if self.scorer is not None:
                path_mm += s_act  # the scored path of this history (the L channel's quantity)
                if path_mm > self.path_bound_mm:  # the capacity proof assumes L_h <= B_L
                    self.path_exceeded = 1
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
                self._end(h, END_ESCAPED, (px, py, pz), (ux, uy, uz), energy)
                return
            zero_run = zero_run + 1 if s_act <= 0.0 else 0
            if zero_run > 3:
                self.tallies["truncated"] += energy
                self.counters["stall"] += 1
                self._end(h, END_TRUNCATED, (px, py, pz), (ux, uy, uz), energy)
                return
