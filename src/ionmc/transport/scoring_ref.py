"""Reference (Python, float64) implementation of the scoring-channel hook (decision 0040).

:class:`ReferenceChannelScorer` owns the int64 channel accumulators ``acc[B, sum_c size_c]``, the
per-history float64 residual columns and the lookup out-of-domain count. The transport calls

* :meth:`begin_step` once per step with ``s_act > 0`` (one water-row lookup at ``E_mid``) and
* :meth:`score_piece` for every scoring piece that lies inside a grid, right after the deposit of
  the same piece ``eps_p`` into the qualified ``edep`` grid (class "step"), and for the point
  deposits (cutoff energy, deposits at ``s_act = 0``; class "local", ``l = 0``).

The scientifically relevant arithmetic is in the shared functions of
:mod:`ionmc.transport.scoring_funcs` (twin of the Warp functions); this class is glue: channel
selection by grid, class mask, species and generation, table reads from numpy arrays,
quantization ``n = floor(x 2^k + 1/2)`` in float64, the integer add and the residual ``x - n q``.
V3-005A supplies other ``species`` and ``gen`` values; nothing here assumes primaries.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

from ionmc._wpfunc import python_twin
from ionmc.scoring import ScoringGrid
from ionmc.species import species_by_id
from ionmc.transport.channels import CLASS_STEP, KINDS, ChannelPlan
from ionmc.transport.funcs import make_transport_funcs
from ionmc.transport.scoring_funcs import KIND_CODES, make_scoring_funcs
from ionmc.transport.tables import TransportTables

_KIND_ASSERT = tuple(KIND_CODES[k] for k in KINDS)
assert _KIND_ASSERT == tuple(range(len(KINDS))), "channel kind codes follow channels.KINDS"


class ReferenceChannelScorer:
    """The scoring hook of the Python reference backend (see the module docstring)."""

    def __init__(
        self,
        plan: ChannelPlan,
        grids: Sequence[ScoringGrid],
        n_batches: int,
        *,
        tables: TransportTables | None = None,
        a_nucleon: int = 1,
    ) -> None:
        self.plan = plan
        self.grids = tuple(grids)
        self.a = a_nucleon
        self.tables = tables
        self.SC = python_twin(make_scoring_funcs)
        self.F = python_twin(make_transport_funcs)
        self.acc = np.zeros((n_batches, plan.total_size), dtype=np.int64)
        self.residual = np.zeros(plan.n_residual)
        self.lookup_ood = 0
        self._scale = [math.ldexp(1.0, c.k) for c in plan.channels]
        self._quantum = [math.ldexp(1.0, -c.k) for c in plan.channels]
        self._kind = [KIND_CODES[c.kind] for c in plan.channels]
        # per lookup table: axis start (ln for a log axis), inverse step, points, per-species rows
        self._lookups = []
        for lk in plan.lookups:
            n = int(lk.axis_values.size)
            log = lk.axis_spacing == "log"
            a0 = float(np.log(lk.axis_values[0]) if log else lk.axis_values[0])
            a1 = float(np.log(lk.axis_values[-1]) if log else lk.axis_values[-1])
            rows = {species_of: v for species_of, v in lk.values.items()}
            self._lookups.append((a0, (n - 1) / (a1 - a0), n, 1 if log else 0, lk.axis, rows))
        self._spectra: list[tuple[float, float, int, int] | None] = []
        for c in plan.channels:
            if c.spectrum is None:
                self._spectra.append(None)
                continue
            e = c.spectrum.edges
            lg = c.spectrum.log
            a0 = math.log(e[0]) if lg else e[0]
            a1 = math.log(e[-1]) if lg else e[-1]
            self._spectra.append((a0, c.spectrum.n_bins / (a1 - a0), c.spectrum.n_bins, int(lg)))
        self._ny_nz = [(g.shape[1], g.shape[2]) for g in grids]
        self.begin_history()
        self.set_step_state(0.0, 0.0, 0.0, 0.0)

    # -- per history / per step ------------------------------------------------------------
    def begin_history(self) -> None:
        """Zero the per-history residual columns and the out-of-domain count."""
        self.residual = np.zeros(self.plan.n_residual)
        self.lookup_ood = 0

    def water_state(self, e_mid: float) -> tuple[float, float]:
        """``(S_mid [MeV/mm], gamma)`` of the water row at ``e_mid`` (one lookup per step)."""
        t = self.tables
        if t is None or t.water_ln_s_mass is None:
            raise ValueError("the scorer needs transport tables with a water row")
        row = t.water_ln_s_mass
        i, f = self.F.log_bin_index(e_mid, t.water_ln_e0, t.water_inv_dln_e, int(row.size))
        ly0, ly1 = float(row[i]), float(row[i + 1])
        s_mass = self.F.interp_exp(ly0, ly1, f)
        gamma = self.SC.loglog_slope(ly0, ly1, t.water_inv_dln_e)
        return float(s_mass) * t.water_density_g_cm3 / 10.0, float(gamma)

    def begin_step(self, e_mid: float, de_mean: float, s_act: float) -> None:
        """Step-level quantities from the step's midpoint energy, CSDA mean loss and path length
        (``s_act > 0``): ``S_mid``, the ramp slope ``k`` and the energy rate ``Edot``."""
        s_mid, gamma = self.water_state(e_mid)
        k = float(self.SC.let_ramp_slope(s_mid, gamma, de_mean, e_mid, s_act))
        e_dot = de_mean / s_act
        self.set_step_state(s_mid, k, e_mid, e_dot)

    def set_step_state(self, s_mid: float, k: float, e_mid: float, e_dot: float) -> None:
        """Set the step quantities directly (hook-level tests with a synthetic stream)."""
        self.s_mid, self.k, self.e_mid, self.e_dot = s_mid, k, e_mid, e_dot

    # -- per piece ---------------------------------------------------------------------------
    def score_piece(
        self,
        batch: int,
        g: int,
        vox: int,
        length: float,
        tau: float,
        eps: float,
        species: int,
        gen: int,
        cls: int,
    ) -> None:
        """Add one piece to every channel of grid ``g`` that selects its class, species and
        generation. ``vox`` is the flat voxel index (the piece must lie inside the grid);
        ``length``/``tau`` are 0 for a local point deposit."""
        plan, SC = self.plan, self.SC
        if cls == CLASS_STEP:
            s_bar, e_bar = SC.piece_state(self.s_mid, self.k, self.e_mid, self.e_dot, tau)
            m1, m2 = SC.piece_moments(s_bar, self.k, length)
        else:
            s_bar = e_bar = m1 = m2 = 0.0
        for ci in range(plan.ch_begin[g], plan.ch_end[g]):
            c = plan.channels[ci]
            if not (c.class_mask & cls) or not plan.species_match[ci, species]:
                continue
            if not (c.gen_lo <= gen <= c.gen_hi):
                continue
            f = 1.0
            col = c.offset + vox
            if c.kind == "FE":
                f = self._lookup_value(c.lookup, species, s_bar, e_bar)
            elif c.kind == "FL":
                spec = self._spectra[ci]
                assert spec is not None
                a0, inv, nb, lg = spec
                b = int(SC.spectrum_bin(e_bar / self.a, a0, inv, nb, lg))
                col = c.offset + vox * (nb + 2) + b
            x = float(SC.channel_value(self._kind[ci], eps, length, m1, m2, s_bar, f))
            n = math.floor(x * self._scale[ci] + 0.5)
            self.acc[batch, col] += n
            if c.residual_column >= 0:
                self.residual[c.residual_column] += x - n * self._quantum[ci]

    def _lookup_value(self, li: int, species: int, s_bar: float, e_bar: float) -> float:
        a0, inv, n, lg, axis, rows = self._lookups[li]
        x = e_bar / self.a if axis == "energy_per_nucleon_mev" else s_bar
        i, fr, inside = self.SC.lookup_bin(x, a0, inv, n, lg)
        if not inside:
            self.lookup_ood += 1
        v = rows[species_by_id(species).name]
        return float(self.F.lerp(float(v[i]), float(v[i + 1]), fr))
