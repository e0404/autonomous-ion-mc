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
from typing import Any

import numpy as np

from ionmc._wpfunc import python_twin
from ionmc.config import MAX_REJECTION_ATTEMPTS, NUCLEAR_MAX_ENERGY_MEV, EffectiveConfig
from ionmc.errors import CounterOverflowError
from ionmc.physics.em import make_em
from ionmc.physics.kinematics import make_kinematics
from ionmc.physics.nuclear import make_nuclear
from ionmc.rng.philox import (
    PURPOSE_NUCLEAR,
    PURPOSE_SOURCE,
    PURPOSE_TRANSPORT,
    child_genealogy_id,
    draw_block,
    key_from_seed,
    u01_py,
)
from ionmc.species import PSEUDO_BASE, species_of_projectile
from ionmc.transport.channels import CLASS_LOCAL, CLASS_STEP, MAX_PARTICLES
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
    NUCLEAR_COUNTER_NAMES,
    NUCLEAR_TALLY_NAMES,
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


STACK_CAPACITY = MAX_PARTICLES
"""Capacity of the per-history LIFO particle stack (decision 0041 section 3)."""
NUC_SPECIES_KEYS = ("n", "p", "d", "a", "g")
"""Light-product keys of the nuclear event diagnostics (``ionmc.nuclear.events.SPECIES``)."""
END_NUCLEAR = 5
"""History end code (diagnostics only, not part of ``tally``): the primary ended in a nuclear
event."""
EVENT_BLOCKS_PER_ATTEMPT = 164
"""Philox blocks reserved per event attempt: the largest slot ``particle_slot(79, 4) = 644`` of 80
products lies in block 161 (4 uniforms per block)."""
_BIG = 1.0e30


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
    raw = merge_partials(
        [part], n, len(eff.requested.scoring), channel_columns(eff), eff.nuclear is not None
    )
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
        self.species_id = 0
        self.generation = 0
        self.nuc = eff.nuclear
        self.counter_names: tuple[str, ...] = COUNTER_NAMES
        self.ntallies: dict[str, float] = {}
        self.nuc_diag: dict[str, dict[str, Any]] = {}  # nuclear event diagnostics (nuclear runs)
        # diagnostics-only nuclear trace of the first trace_histories histories: accepted events
        # (h, gid, target, counts[5], Z_r, A_r, attempts, T1) and pushed secondaries (h, gid,
        # parent gid, species, generation, T, direction, position); no effect on the transport
        self.want_diag = False  # set by run_range
        self.nuc_trace: dict[str, list[list[float]]] = {"events": [], "secondaries": []}
        if self.nuc is not None:  # conditional nuclear blocks (decision 0041 sections 2-5)
            self.counter_names = COUNTER_NAMES + NUCLEAR_COUNTER_NAMES
            self.NU = python_twin(make_nuclear)
            self.tab_p, self.tab_d = self.tab, self.nuc.deuteron_tables
            self.mass_p, self.mass_d = self.mass, r(self.tab_d.projectile.mass_mev)
            self.e_cut_p, self.e_cut_d = self.e_cut, float(self.nuc.e_cut_deuteron_mev)
            self._models: dict[int, tuple[object, dict[str, np.ndarray]]] = {}
            self.counters = dict.fromkeys(self.counter_names, 0)
            self.ntallies = dict.fromkeys(NUCLEAR_TALLY_NAMES, 0.0)
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

    # -- table reads (shared bin location and interpolation, numpy memory access) -------------
    def _stopping_range(self, m: int, e: float) -> tuple[float, float]:
        t, r = self.tab, self.R
        i, f = self.F.log_bin_index(r(e), r(t.ln_e0[m]), r(t.inv_dln_e[m]), t.n_e)
        ls = t.ln_s_mass[m]
        s = self.F.interp_exp(r(ls[i]), r(ls[i + 1]), f)
        rng = self.F.range_in_bin(
            r(t.r_mass[m, i]), r(t.f_mass[m, i]), r(t.d_f[m, i]), r(1.0 / t.inv_dln_e[m]), f
        )
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

    def _score_local(
        self, batch: int, g: int, vox: int, de: float, species: int = -1, gen: int = -1
    ) -> None:
        """Class "local" deposit (cutoff energy, deposit at ``s_act = 0``) into the channels; the
        species and generation default to those of the particle being transported."""
        if self.scorer is not None and vox >= 0:
            self.scorer.score_piece(
                batch, g, vox, 0.0, 0.0, de,
                self.species_id if species < 0 else species,
                self.generation if gen < 0 else gen,
                CLASS_LOCAL,
            )  # fmt: skip

    def _deposit_point(
        self,
        batch: int,
        pos: tuple[float, float, float],
        de: float,
        species: int = -1,
        gen: int = -1,
    ) -> None:
        """Point deposit (the energy left at the cutoff) in every grid."""
        r = self.R
        p = self.V(r(pos[0]), r(pos[1]), r(pos[2]))
        for g, grid in enumerate(self.grids):
            nx, ny, nz = grid.shape
            ix, iy, iz, _inside = self.F.grid_index(p, self.g_origin[g], self.g_inv[g], nx, ny, nz)
            self._score_local(
                batch, g, self._deposit_voxel(batch, g, ix, iy, iz, de), de, species, gen
            )

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
        counter_rows = np.zeros((n, len(self.counter_names)), dtype=np.int32)
        n_nuc_t = len(NUCLEAR_TALLY_NAMES) if self.nuc is not None else 0
        tally_rows = np.zeros((n, N_FIXED_TALLIES + 2 * n_g + n_chan_cols + n_nuc_t))
        c_chan = N_FIXED_TALLIES + 2 * n_g
        for h in range(h0, h1):
            # per-history accumulators: a row depends on this history alone
            self.tallies = dict.fromkeys(TALLY_NAMES, 0.0)
            self.outside = [0.0] * n_g
            self.quant = [0.0] * n_g
            self.counters = dict.fromkeys(self.counter_names, 0)
            self.ntallies = dict.fromkeys(NUCLEAR_TALLY_NAMES, 0.0) if self.nuc is not None else {}
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
                tally_rows[row, c_chan : c_chan + n_chan_cols - 2] = self.scorer.residual
                tally_rows[row, c_chan + n_chan_cols - 2] = self.scorer.lookup_ood
                tally_rows[row, c_chan + n_chan_cols - 1] = self.path_exceeded
            if self.nuc is not None:
                tally_rows[row, -n_nuc_t:] = [self.ntallies[k] for k in NUCLEAR_TALLY_NAMES]
            counter_rows[row] = [self.counters[k] for k in self.counter_names]
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
        part = rows_to_partial(
            h0, h1, tally_rows, counter_rows, self.edep, diagnostics,
            channel_acc=None if self.scorer is None else self.scorer.acc,
        )  # fmt: skip
        if self.nuc is not None:
            part.meta["nuclear_diagnostics"] = self.nuc_diag
            if self.want_diag:
                part.meta["nuclear_trace"] = {
                    k: np.array(v, dtype=np.float64).reshape(-1, 12)
                    for k, v in self.nuc_trace.items()
                }
        return part

    def _end(
        self,
        h: int,
        code: int,
        pos: tuple[float, float, float],
        direction: tuple[float, float, float],
        energy: float,
    ) -> None:
        if self.generation != 0:  # the end record is that of the primary
            return
        i = h - self.h_base
        self.end_pos[i] = pos
        self.end_dir[i] = direction
        self.end_code[i] = code
        self.end_energy[i] = energy
        self.end_ctrl[i] = self.ctrl_res
        self.end_ctrl_sum[i] = self.ctrl_sum

    # -- nuclear events (decision 0041 sections 2-4) ------------------------------------------
    def _nuc_u(self, h: int, gid: int, block: int) -> list[float]:
        """The four uniforms of block ``block`` of the PURPOSE_NUCLEAR stream of the particle."""
        w = draw_block(self.key, h, gid, block, PURPOSE_NUCLEAR)
        return [u01_py(x, "float64") for x in w]

    def _sigma(self, rows: object, e: float) -> float:
        """``Sigma_mass(E)`` [cm2/g] of a material's rows: step lookup of the cell by the shared
        ``grid_locate``, lin-lin in E (clamped), as ``MaterialNuclear.sigma_at``."""
        g = rows.grid_e_mev  # type: ignore[attr-defined]
        k = int(self.NU.grid_locate(e, g, g.size))
        t = min(max((e - g[k]) / (g[k + 1] - g[k]), 0.0), 1.0)
        row = rows.sigma_mass_cm2_g  # type: ignore[attr-defined]
        return float((1.0 - t) * row[k] + t * row[k + 1])

    def _target_model(self, tgt: int) -> tuple[object, dict[str, np.ndarray]]:
        """Event model (masses, separation energies, residual masses) and product rows of a table
        target, built on first use from the cached AME2020 file."""
        if tgt not in self._models:
            from ionmc.data import cache
            from ionmc.data.ame import load_ame2020
            from ionmc.nuclear.events import build_event_model

            if not hasattr(self, "_ame"):
                path = cache.verify("ame2020-mass", cache.resolve_cache_dir(None))
                self._ame = load_ame2020(path.read_text(encoding="ascii"))
            t = self.nuc.table  # type: ignore[union-attr]
            info = t.info["targets"][tgt]
            self._models[tgt] = (
                build_event_model(self._ame, int(info["z"]), int(info["a"])),
                t.product_rows(tgt),
            )
        return self._models[tgt]

    def _candidate(
        self,
        h: int,
        batch: int,
        gid: int,
        nc: int,
        rows: object,
        s_hat: float,
        energy: float,
        pos: tuple[float, float, float],
        direction: tuple[float, float, float],
        vox: tuple[int, int, int],
        stack: list[tuple[float, ...]],
    ) -> tuple[int, float, bool]:
        """A nuclear candidate at the post-step point with energy ``energy`` = E1: accept with
        probability ``Sigma(E1) / S^(E0)`` (``thinning_accept``). Returns the next nuclear block
        counter, the (resampled) optical depth and whether the primary has ended."""
        u = self._nuc_u(h, gid, nc)  # (u_accept, u_n_lambda, u_target, .)
        nc += 1
        r = self.R
        accepted, violation = self.NU.thinning_accept(
            r(u[0]), r(self._sigma(rows, energy)), r(s_hat)
        )
        if violation:  # fail closed: the majorant was violated
            self.counters["majorant_violation"] += 1
            self.tallies["unaccounted"] += energy
            self._end(h, END_NUCLEAR, pos, direction, energy)
            return nc, 0.0, True
        if not accepted:  # fictitious: resample the optical depth
            return nc, -math.log(u[1]), False
        self._event(h, batch, gid, nc, rows, u[2], energy, pos, direction, vox, stack)
        self._end(h, END_NUCLEAR, pos, direction, energy)
        return nc, 0.0, True

    def _event(
        self,
        h: int,
        batch: int,
        gid: int,
        nc: int,
        rows: object,
        u_target: float,
        t1: float,
        pos: tuple[float, float, float],
        direction: tuple[float, float, float],
        vox: tuple[int, int, int],
        stack: list[tuple[float, ...]],
    ) -> None:
        """An accepted nuclear event of the primary (energy ``t1`` = T1): the parent ends, the
        products are scored and the secondary p and d pushed (module docstring, decision 0041
        section 3 with the amendment of 2026-10-07)."""
        from ionmc.nuclear.events import interp_rows, sample_event_scalar

        r, nt, ctr = self.R, self.ntallies, self.counters
        cum = rows.cum_fraction_at(t1)  # type: ignore[attr-defined]
        tgt_list = rows.target_index  # type: ignore[attr-defined]
        chosen = -1
        for k in range(len(tgt_list)):
            chosen = int(
                self.NU.select_target(
                    r(u_target), r(float(cum[k])), k, chosen, int(k == len(tgt_list) - 1)
                )
            )  # fmt: skip
        model, prow = self._target_model(tgt_list[chosen])
        grid = self.nuc.table.arrays["grid_e_mev"]  # type: ignore[union-attr]
        erows = interp_rows(
            grid, prow["lam"], prow["edges_mev"], prow["r_pre"], prow["recoil_t_cm_mev"], t1
        )
        cache_blk: dict[int, list[float]] = {}

        def uni(attempt: int, slot: int) -> float:
            blk = nc + (attempt - 1) * EVENT_BLOCKS_PER_ATTEMPT + slot // 4
            if blk not in cache_blk:
                cache_blk[blk] = self._nuc_u(h, gid, blk)
            return cache_blk[blk][slot % 4]

        ev = sample_event_scalar(model, erows, t1, uni, self.NU)  # type: ignore[arg-type]
        if not ev.accepted:  # 64 attempts exhausted: fail closed
            ctr["nuclear_rejection_limit"] += 1
            self.tallies["unaccounted"] += t1
            return
        # independent event record (per table target): events, light products per species and
        # residual (Z_r, A_r), from which the binding of the run is recomputed from AME2020 masses
        rec = self.nuc_diag.setdefault(
            str(int(tgt_list[chosen])),
            {"events": 0, "light": dict.fromkeys(NUC_SPECIES_KEYS, 0), "residual": {}},
        )
        if self.want_diag and h < self.cfg.diagnostics.trace_histories:
            self.nuc_trace["events"].append(
                [h, gid, int(tgt_list[chosen]), *(int(c) for c in ev.counts), int(ev.z_r),
                 int(ev.a_r), ev.attempts, t1]
            )  # fmt: skip
        rec["events"] += 1
        for key, cnt in zip(NUC_SPECIES_KEYS, ev.counts, strict=True):
            rec["light"][key] += int(cnt)
        res_key = f"{int(ev.z_r)},{int(ev.a_r)}"
        rec["residual"][res_key] = rec["residual"].get(res_key, 0) + 1
        masses = [float(x) for x in model.species_mass_mev]  # type: ignore[attr-defined]
        gen = self.generation + 1
        # e_lab is the total energy of the lab product, T_lab = e_lab - m
        t_sum = alpha_t = 0.0
        # frame of the parent direction: e1, e2 perpendicular, the event frame has z along the beam
        d3 = self.V(r(direction[0]), r(direction[1]), r(direction[2]))
        e1, e2 = self.F.orthonormal_basis(d3)
        n_children = 0
        for sp_i, _ecm, _mu, _phi, e_lab, qx, qy, qz in ev.particles:
            si = int(sp_i)
            t_lab = e_lab - masses[si] if si < 4 else e_lab
            t_sum += t_lab
            if si == 0:
                nt["nuclear_escaped_neutron"] += t_lab
            elif si == 4:
                nt["nuclear_escaped_gamma"] += t_lab
            elif si == 3:
                alpha_t += t_lab
            else:  # p (1) or d (2): a secondary of generation + 1
                species = si - 1
                cut = self.e_cut_p if species == 0 else self.e_cut_d
                if t_lab < cut:  # below the species cutoff: deposited locally, tallied as cutoff
                    self.tallies["cutoff"] += t_lab
                    self._deposit_point(batch, pos, t_lab, species, gen)
                    continue
                n_children += 1
                try:
                    cid = child_genealogy_id(gid, self.generation, n_children)
                except CounterOverflowError:
                    ctr["genealogy_overflow"] += 1
                    self.tallies["unaccounted"] += t_lab
                    continue
                if len(stack) >= STACK_CAPACITY:
                    ctr["queue_overflow"] += 1
                    self.tallies["unaccounted"] += t_lab
                    continue
                pm = math.sqrt(qx * qx + qy * qy + qz * qz)
                if pm > 0.0:
                    c = (qx / pm, qy / pm, qz / pm)
                    dx = c[0] * float(e1[0]) + c[1] * float(e2[0]) + c[2] * direction[0]
                    dy = c[0] * float(e1[1]) + c[1] * float(e2[1]) + c[2] * direction[1]
                    dz = c[0] * float(e1[2]) + c[1] * float(e2[2]) + c[2] * direction[2]
                else:
                    dx, dy, dz = direction
                mass = self.mass_p if species == 0 else self.mass_d
                pv = float(self.K.pv_mev(r(t_lab), mass))
                if self.want_diag and h < self.cfg.diagnostics.trace_histories:
                    self.nuc_trace["secondaries"].append(
                        [h, cid, gid, species, gen, t_lab, dx, dy, dz, pos[0], pos[1], pos[2]]
                    )
                stack.append(
                    (pos[0], pos[1], pos[2], dx, dy, dz, t_lab, float(species), float(cid),
                     float(gen), float(vox[0]), float(vox[1]), float(vox[2]), pv)
                )  # fmt: skip
        local = alpha_t + ev.recoil_t_mev
        nt["nuclear_local"] += local
        nt["nuclear_alpha_local"] += alpha_t
        nt["nuclear_binding"] += ev.binding_mev
        nt["nuclear_imbalance"] += ev.imbalance_mev
        self._deposit_point(batch, pos, local, PSEUDO_BASE, gen)
        # per-event ledger: T1 = sum T_lab + T_r + binding + imbalance
        miss = t1 - (t_sum + ev.recoil_t_mev + ev.binding_mev + ev.imbalance_mev)
        if abs(miss) > 1.0e-9 * t1:
            ctr["nuclear_conservation"] += 1
            self.tallies["unaccounted"] += miss

    # -- one history --------------------------------------------------------------------------
    def _history(self, h: int) -> None:
        cfg, F, K, V, r = self.cfg, self.F, self.K, self.V, self.R
        src = cfg.source
        key = self.key
        nx, ny, nz = self.shape

        self.generation = 0
        self._set_species(0)
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
        if self.nuc is not None and energy > NUCLEAR_MAX_ENERGY_MEV:
            # nuclear runs: the cross sections end at 250 MeV and the Gaussian source is unbounded
            # (decision 0041 section 5, amended 2026-10-08), whatever the stopping tables cover;
            # the energy is booked as initial and unaccounted so that the ledger still closes
            self.counters["source_energy_out_of_range"] += 1
            self.tallies["initial"] += energy
            self.tallies["unaccounted"] += energy
            self._end(h, END_SOURCE_REJECTED, (px, py, pz), d, energy)
            return
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

        self._path_mm = 0.0
        stack: list[tuple[float, ...]] = [
            (px, py, pz, ux, uy, uz, energy, 0.0, 0.0, 0.0, ix, iy, iz, p1v1)
        ]
        # per-history LIFO stack (decision 0041 section 3): entries are (position, direction, T,
        # species id, genealogy id, generation, the parent's geometry voxel indices, p1v1)
        while stack:
            self._particle(h, stack.pop(), stack)

    def _set_species(self, species: int) -> None:
        """Switch the per-particle table, mass, cutoff and species id (nuclear runs: p or d)."""
        if self.nuc is not None:
            if species == 0:
                self.tab, self.mass, self.e_cut = self.tab_p, self.mass_p, self.e_cut_p
            else:
                self.tab, self.mass, self.e_cut = self.tab_d, self.mass_d, self.e_cut_d
            self.species_id = species  # without nuclear the species id is the source's

    def _particle(self, h: int, entry: tuple[float, ...], stack: list[tuple[float, ...]]) -> None:
        """Transport one particle of the stack (the step loop of decision 0039)."""
        cfg, F, EM, K, V, r = self.cfg, self.F, self.EM, self.K, self.V, self.R
        ph = cfg.physics
        key = self.key
        batch = h % self.n_batches
        nx, ny, nz = self.shape
        o, sp = self.geo.origin_mm, self.geo.spacing_mm
        px, py, pz, ux, uy, uz, energy = entry[:7]
        species, gid, gen = int(entry[7]), int(entry[8]), int(entry[9])
        ix, iy, iz = int(entry[10]), int(entry[11]), int(entry[12])
        p1v1 = entry[13]
        self.generation = gen
        self._set_species(species)
        trace_this = h < cfg.diagnostics.trace_histories and gen == 0
        # nuclear interactions of the PRIMARY proton only (decision 0041 section 3: secondaries
        # and deuterons have none in slice A)
        nuc_on = self.nuc is not None and gen == 0
        nu_rows = None
        n_lam = 0.0
        nc = 0
        if nuc_on:  # birth: the optical depth to the first candidate (nuclear block 0, slot 0)
            n_lam = -math.log(self._nuc_u(h, gid, 0)[0])
            nc = 1
        steps = 0
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
            s_hat = 0.0
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
            if nuc_on:
                # majorant of the step: window value at the grid cell of E0 (step lookup) or the
                # end-of-range value for a range-limited step; d_nuc in mm (rate rho * S^ / cm)
                assert self.nuc is not None
                nu_rows = self.nuc.rows[m]
                k0 = int(self.NU.grid_locate(energy, nu_rows.grid_e_mev, nu_rows.grid_e_mev.size))
                majorant_row = (
                    nu_rows.sigma_hat_end if float(s_rg) < float(s_el) else nu_rows.sigma_hat_window
                )
                s_hat = float(majorant_row[k0])
                d_nuc = float(self.NU.nuclear_step_limit(r(n_lam), r(rho), r(s_hat)))
                s_w, reason = F.select_step_nuclear(d_geo, s_el, s_rg, self.c_smax, r(d_nuc))
            else:
                s_w, reason = F.select_step(d_geo, s_el, s_rg, self.c_smax)
            s = float(s_w)

            # block A
            wa = draw_block(key, h, gid, blocks, PURPOSE_TRANSPORT)
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
                    wb = draw_block(key, h, gid, blocks, PURPOSE_TRANSPORT)
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
                draw_block(key, h, gid, blocks, PURPOSE_TRANSPORT)
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
                    if self.nuc is None:
                        self.scorer.begin_step(energy - 0.5 * mean_f, mean_f, s_act)
                    else:  # the species' water row and mass number
                        self.scorer.begin_step(energy - 0.5 * mean_f, mean_f, s_act, species)
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
                self._path_mm += s_act  # scored path of this history (the L channel's quantity)
                if self._path_mm > self.path_bound_mm:  # the capacity proof assumes L_h <= B_L
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
            if nuc_on:
                n_lam -= rho * s_hat * s_act / 10.0
                # majorant check after EVERY step (decision 0041 section 2): the Gamma straggling
                # tail is unbounded, so Sigma(E1) <= S^(E0) is not guaranteed by the window
                _acc, viol = self.NU.thinning_accept(
                    r(0.5), r(self._sigma(nu_rows, energy)), r(s_hat)
                )
                if viol:
                    self.counters["majorant_violation"] += 1
                    self.tallies["unaccounted"] += energy
                    self._end(h, END_NUCLEAR, (px, py, pz), (ux, uy, uz), energy)
                    return
                if int(reason) == 4 and axis2 < 0 and not exited and energy > self.e_cut:
                    # candidate at the post-step point (leg 2 not truncated)
                    nc, n_lam, done = self._candidate(
                        h, batch, gid, nc, nu_rows, s_hat, energy, (px, py, pz), (ux, uy, uz),
                        (ix, iy, iz), stack,
                    )  # fmt: skip
                    if done:
                        return
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
