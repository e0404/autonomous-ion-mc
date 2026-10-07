"""V3-004 acceptance rows A1 to A6, A10, A14 and the ratio reduction on the python backend
(plan validation/plans/v3-004-acceptance.md; criteria are frozen, CI tier, fixed small seeds).

Notation of the plan: ``X_v`` is the sum over batches and pieces of channel X in voxel v, ``n_v``
the piece count of channel N and ``delta_X = n_v q_X / 2`` the fixed-point bound.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace
from fractions import Fraction
from typing import Any

import numpy as np
import pytest

from ionmc.config import (
    DiagnosticsOptions,
    EffectiveConfig,
    SimulationConfig,
    validate,
)
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.lookup import LookupTable
from ionmc.materials import ALUMINIUM, WATER, Material
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource, StoppingTable, build_table
from ionmc.scoring import ScoringGrid, TallyRequest, min_defined_batches, reduce_ratio
from ionmc.simulation import Simulation
from ionmc.species import species_by_name
from ionmc.transport.channels import CLASS_LOCAL, CLASS_STEP, ChannelPlan, compile_channels
from ionmc.transport.reference import run_reference_range
from ionmc.transport.run import channel_columns, run_transport
from ionmc.transport.scoring_ref import ReferenceChannelScorer
from ionmc.transport.tally import RawTransport, build_diagnostics, merge_partials

MakeConfig = Callable[..., SimulationConfig]
RHO = 1.0e-12  # float rounding term of python / float64


class ConstSource:
    """Synthetic stopping source with a constant mass stopping power per material (water
    ``c_water``, anything else ``c_other``, MeV cm2/g): constant S_w, exact linear range."""

    name = "constant-s"

    def __init__(self, c_water: float = 10.0, c_other: float = 5.0) -> None:
        self.c_water, self.c_other = c_water, c_other

    def table(self, material: Material, projectile: Any) -> StoppingTable:
        e = np.geomspace(1.0, 500.0, 400)
        c = self.c_water if material.name == WATER.name else self.c_other
        return build_table(projectile, material, e, np.full_like(e, c), e[0] / c, {"source": "c"})


def _grid(nz: int = 30, dz: float = 2.0, oz: float = 0.0, name: str = "dose") -> ScoringGrid:
    return ScoringGrid((-10.0, -10.0, oz), (20.0, 20.0, dz), (1, 1, nz), name=name)


def _reqs(*quantities: str, grid: str = "dose") -> tuple[TallyRequest, ...]:
    return tuple(TallyRequest(q, grid, q) for q in quantities)  # type: ignore[arg-type]


def _cfg(make_config: MakeConfig, tallies: tuple[TallyRequest, ...], **kw: Any) -> SimulationConfig:
    kw.setdefault("position", (1.0, 1.0, 0.0))
    lookups = kw.pop("lookups", ())
    stopping = kw.pop("stopping", None)
    cfg = make_config(**kw)
    cfg = replace(cfg, tallies=tallies, lookups=lookups)
    if stopping is not None:
        cfg = replace(cfg, physics=replace(cfg.physics, stopping=stopping))
    return cfg


def _slab(grid_dz: float, nz: int, material: Material = WATER) -> VoxelGeometry:
    return VoxelGeometry(
        origin_mm=(-30.0, -30.0, 0.0),
        spacing_mm=(60.0, 60.0, grid_dz),
        shape=(1, 1, nz),
        materials=(material,),
        material_index=np.zeros((1, 1, nz), dtype=np.int32),
    )


class Run:
    """A validated configuration and its raw channels (sums over batches in channel units)."""

    def __init__(self, cfg: SimulationConfig) -> None:
        self.eff: EffectiveConfig = validate(cfg)
        n = cfg.run.n_histories
        self.part = run_reference_range(self.eff, 0, n)  # per-batch integer arrays
        self.raw: RawTransport = merge_partials(
            [self.part], n, len(cfg.scoring), channel_columns(self.eff)
        )
        self.raw.diagnostics = build_diagnostics(
            [self.part],
            cfg.diagnostics.track_end_positions,
            cfg.diagnostics.escape_records,
            cfg.diagnostics.trace_histories,
        )
        assert self.eff.channels is not None and self.raw.channels is not None
        self.plan: ChannelPlan = self.eff.channels
        self.acc = self.raw.channels.acc

    def q_of(self, name: str) -> Any:
        return next(d for d in self.plan.quantities if d.name == name)

    def total(self, ci: int) -> np.ndarray:
        """Sum over batches of channel ``ci`` [channel unit], per voxel (and bin)."""
        c = self.plan.channels[ci]
        return self.acc[:, c.offset : c.offset + c.size].sum(axis=0) * c.quantum

    def n_pieces(self, grid: int = 0, local: bool = False) -> np.ndarray:
        """Exact piece counts per voxel of the automatic N channels (step or local class)."""
        return self.total(self.plan.count_channel(grid, CLASS_LOCAL if local else CLASS_STEP))

    def delta(self, ci: int) -> np.ndarray:
        """Rounding bound ``n_v q / 2`` of channel ``ci`` with the piece counts of its classes."""
        c = self.plan.channels[ci]
        n_v = (self.n_pieces(c.grid) if c.class_mask & CLASS_STEP else 0.0) + (
            self.n_pieces(c.grid, local=True) if c.class_mask & CLASS_LOCAL else 0.0
        )
        return n_v * c.quantum / 2.0

    def residual(self, ci: int) -> float:
        assert self.raw.channels is not None
        return self.raw.channels.residual[self.plan.channels[ci].residual_column]


# --------------------------------------------------------------------------- A1


@pytest.mark.parametrize("material", [WATER, ALUMINIUM], ids=["water", "aluminium"])
def test_a1_constant_s_closed_forms(make_config: MakeConfig, material: Material) -> None:
    cfg = _cfg(
        make_config, _reqs("let_t", "let_d", "let_d_eps"), energy=50.0, n=8, n_batches=4,
        straggling=False, mcs=False, max_step=1.0, stopping=ConstSource(),
        geometry=BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, 60.0), material),
        scoring=(_grid(),),
    )  # fmt: skip
    r = Run(cfg)
    s_w = r.eff.tables.s_water(20.0)
    assert s_w == pytest.approx(1.0, rel=1e-12)
    n_v = r.n_pieces()
    traversed = n_v > 0
    assert traversed.sum() >= 15
    for qname in ("let_t", "let_d", "let_d_eps"):
        qd = r.q_of(qname)
        x, y = r.total(qd.numerator), r.total(qd.denominator)
        dx, dy = r.delta(qd.numerator), r.delta(qd.denominator)
        ok = y > 0
        assert ok.any()
        err = np.abs(x[ok] / y[ok] - s_w)
        bound = RHO * s_w + (dx[ok] + s_w * dy[ok]) / y[ok]
        assert np.all(err <= bound), (qname, err.max(), bound[err > bound][:3])
    # global closure: sum_v LS + res_LS = S (sum_v L + res_L)
    ls = r.q_of("let_t").numerator
    ll = r.q_of("let_t").denominator
    lhs = r.total(ls).sum() + r.residual(ls)
    rhs = s_w * (r.total(ll).sum() + r.residual(ll))
    assert abs(lhs - rhs) <= 1e-12 * abs(rhs)


# --------------------------------------------------------------------------- A3


def test_a3_pencil_fluence_aligned_grid(make_config: MakeConfig) -> None:
    n, nb = 8, 4
    cfg = _cfg(
        make_config, _reqs("fluence"), energy=60.0, n=n, n_batches=nb, straggling=False,
        mcs=False, max_step=2.0, scoring=(_grid(),),
        diagnostics=DiagnosticsOptions(track_end_positions=True),
    )  # fmt: skip
    r = Run(cfg)
    cfg2 = r.eff.requested
    grid = cfg2.scoring[0]
    dz = grid.spacing_mm[2]
    vol = float(np.prod(grid.spacing_mm))
    area = grid.spacing_mm[0] * grid.spacing_mm[1]
    ci = r.q_of("fluence").numerator
    l_v = r.total(ci)
    z_end = float(r.raw.diagnostics["end_position_mm"][0, 2])
    full = np.arange(grid.n_voxels)[(np.arange(grid.n_voxels) + 1) * dz <= z_end - 1e-9]
    assert full.size >= 10
    bound = RHO / area + r.delta(ci)[full] / (n * vol)
    err = np.abs(l_v[full] / (n * vol) - 1.0 / area)
    assert np.all(err <= bound), err.max()
    # batch standard deviation of the per-primary fluence is at most the same bound
    q = Simulation(cfg).run().grids[0].quantities["fluence"]
    # (std of the mean over batches; hinge splits are identical here, so it is rounding level)
    assert np.all(q.std.reshape(-1)[full] <= bound)
    chord = n * z_end
    tot = l_v.sum() + r.residual(ci)
    assert abs(tot - chord) <= 1e-12 * chord


# --------------------------------------------------------------------------- A4 / A4b


def _energy_at_depth(eff: EffectiveConfig, e0: float, z_mm: np.ndarray) -> np.ndarray:
    """Table CSDA energy ``Rinv(R(E0) - rho z)`` at depth ``z`` in water (single-shot inversion of
    the run's own range table). Informative only (plan footnote 1): it carries the trapezoid bias
    of the range table (decision 0039 follow-up V3-003D) and is not the reference of A4/A4b."""
    t = eff.tables
    r0 = t.range_g_cm2(0, e0)
    return np.array([t.energy_from_range(0, r0 - z / 10.0) for z in z_mm])


def _energy_transported(r: Run, z_mm: np.ndarray) -> np.ndarray:
    """Energy of the transported particle at depth ``z`` (plan footnote 1 of A4/A4b):
    ``E(z) = Rinv(R(E_b) - rho (z - z_b))`` anchored at the last engine step boundary
    ``(z_b, E_b)`` at or before ``z``, from the trace of the deterministic history (the scored
    ``LET_t`` is the path average of ``S`` of this particle, so ``LET_t dz`` must equal its own
    energy loss). Water, ``rho`` = 1 g/cm3 (the transport slab is water along z)."""
    t = r.eff.tables
    tr = r.raw.diagnostics["trace"]
    sel = tr["history"] == 0
    z_b = np.concatenate([[0.0], tr["z_mm"][sel]])
    e_b = np.concatenate([[r.eff.requested.source.kinetic_energy_mev], tr["energy_mev"][sel]])
    out = np.empty(len(z_mm))
    for j, z in enumerate(np.asarray(z_mm, dtype=np.float64)):
        b = int(np.searchsorted(z_b, z + 1e-12, side="right")) - 1
        out[j] = t.energy_from_range(0, t.range_g_cm2(0, e_b[b]) - (z - z_b[b]) / 10.0)
    return out


def _csda_run(make_config: MakeConfig, tallies: tuple[TallyRequest, ...], e0: float, dz_t: float,
              grid: ScoringGrid) -> Run:  # fmt: skip
    nz = int(round(170.0 / dz_t))
    return Run(
        _cfg(
            make_config,
            tallies,
            energy=e0,
            n=2,
            n_batches=2,
            straggling=False,
            mcs=False,
            max_step=dz_t,
            geometry=_slab(dz_t, nz),
            scoring=(grid,),
            diagnostics=DiagnosticsOptions(trace_histories=1),
        )  # fmt: skip
    )


def _z_cut(r: Run) -> float:
    return float(r.eff.requested.source.kinetic_energy_mev * 0 + 1) and _cut_depth(r)


def _cut_depth(r: Run) -> float:
    t = r.eff.tables
    e0 = r.eff.requested.source.kinetic_energy_mev
    ec = r.eff.requested.physics.e_cut_mev
    return (t.range_g_cm2(0, e0) - t.range_g_cm2(0, ec)) * 10.0


def test_a4_csda_let_t_and_let_d(make_config: MakeConfig) -> None:
    grid = _grid(nz=85, dz=2.0)
    r = _csda_run(make_config, _reqs("let_t", "let_d"), 150.0, 2.0, grid)
    z_cut = _cut_depth(r)
    iz = np.arange(85)
    sel = iz[((iz + 1) * 2.0 <= z_cut) & (iz * 2.0 >= 100.0)]
    assert sel.size >= 10
    n_prim = 2
    e_lo = _energy_transported(r, (sel + 1) * 2.0)
    e_hi = _energy_transported(r, sel * 2.0)
    de = e_hi - e_lo
    de_tab = _energy_at_depth(r.eff, 150.0, sel * 2.0) - _energy_at_depth(
        r.eff, 150.0, (sel + 1) * 2.0
    )
    lt, ld = r.q_of("let_t"), r.q_of("let_d")
    l_v, ls_v, ls2_v = r.total(lt.denominator), r.total(lt.numerator), r.total(ld.numerator)
    let_t = ls_v / l_v
    assert np.all(np.abs(let_t[sel] * 2.0 - de) / de <= 1e-4), np.abs(
        let_t[sel] * 2.0 / de - 1
    ).max()
    # informative, non-gating (plan footnote 1): deviation from the table CSDA energy
    print("A4 worst vs transported E(z):", float(np.max(np.abs(let_t[sel] * 2.0 - de) / de)))
    print("A4 worst vs table CSDA E(z) (non-gating):",
          float(np.max(np.abs(let_t[sel] * 2.0 - de_tab) / de_tab)))  # fmt: skip
    # LET_d = int S dE / dE on the same water table (numpy quadrature)
    for j, v in enumerate(sel):
        e = np.linspace(e_lo[j], e_hi[j], 2001)
        s = np.array([r.eff.tables.s_water(x) for x in e])
        quad = float(np.trapezoid(s, e) / (e_hi[j] - e_lo[j]))
        assert abs(ls2_v[v] / ls_v[v] - quad) / quad <= 1e-3
    assert n_prim == 2


def _a4b_worst(
    r: Run, z_cut: float, z_min: float, z_max: float, *, table_csda: bool = False
) -> float:
    """Worst relative error of LS_v / N against E(z1) - E(z2) of the transported particle (plan
    footnote 1; ``table_csda`` selects the single-shot table inversion, informative) over the
    1.5 mm bins (offset 0.25 mm) fully inside ``[z_min, min(z_max, z_cut)]``."""
    iz = np.arange(113)
    z1, z2 = 0.25 + iz * 1.5, 0.25 + (iz + 1) * 1.5
    sel = iz[(z2 <= min(z_cut, z_max)) & (z1 >= z_min)]
    assert sel.size >= 5
    if table_csda:
        de = _energy_at_depth(r.eff, 150.0, z1[sel]) - _energy_at_depth(r.eff, 150.0, z2[sel])
    else:
        de = _energy_transported(r, z1[sel]) - _energy_transported(r, z2[sel])
    ls = r.total(r.q_of("let_t").numerator)[sel] / 2.0  # per primary
    return float(np.max(np.abs(ls - de) / de))


A4B_GRID = ScoringGrid((-10.0, -10.0, 0.25), (20.0, 20.0, 1.5), (1, 1, 113), name="dose")


def test_a4b_offset_grid_all_bins_proximal_to_cutoff(make_config: MakeConfig) -> None:
    r = _csda_run(make_config, _reqs("let_t"), 150.0, 1.0, A4B_GRID)
    z_cut = _cut_depth(r)
    assert _a4b_worst(r, z_cut, 100.0, 1e9) <= 1e-4
    # informative, non-gating record of the table-CSDA comparison (range-table trapezoid bias,
    # follow-up V3-003D; plan footnote 1)
    print("A4b worst vs transported E(z):", _a4b_worst(r, z_cut, 100.0, 1e9))
    print("A4b worst vs table CSDA E(z) (non-gating):",
          _a4b_worst(r, z_cut, 100.0, 1e9, table_csda=True))  # fmt: skip


@pytest.mark.parametrize("max_step", [0.5, 0.25])
def test_a4b_smaller_transport_steps(make_config: MakeConfig, max_step: float) -> None:
    """Supplementary (the plan footnote reports 1 / 0.5 / 0.25 mm): the same criterion with
    smaller transport steps."""
    r = _csda_run(make_config, _reqs("let_t"), 150.0, max_step, A4B_GRID)
    assert _a4b_worst(r, _cut_depth(r), 100.0, 1e9) <= 1e-4


def test_a4b_offset_grid_plateau_bins_and_negative_control(
    make_config: MakeConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Supplementary (not the frozen row): bins between 100 and 140 mm meet 1e-4 with the ramp;
    the negative control (k forced to 0) does not in the bins up to the cutoff."""
    r = _csda_run(make_config, _reqs("let_t"), 150.0, 1.0, A4B_GRID)
    z_cut = _cut_depth(r)
    assert _a4b_worst(r, z_cut, 100.0, 140.0) <= 1e-4
    ramp_all = _a4b_worst(r, z_cut, 100.0, 1e9)
    from ionmc._wpfunc import python_twin
    from ionmc.transport.scoring_funcs import make_scoring_funcs

    monkeypatch.setattr(python_twin(make_scoring_funcs), "let_ramp_slope", lambda *a: 0.0)
    control = _csda_run(make_config, _reqs("let_t"), 150.0, 1.0, A4B_GRID)
    assert _a4b_worst(control, z_cut, 100.0, 1e9) > 1e-4
    print("A4b worst: ramp", ramp_all, "control", _a4b_worst(control, z_cut, 100.0, 1e9))


# --------------------------------------------------------------------------- A5 / A6


def _all_physics_run(make_config: MakeConfig, tallies: tuple[TallyRequest, ...], **kw: Any) -> Run:
    kw.setdefault("energy", 60.0)
    return Run(
        _cfg(make_config, tallies, n=kw.pop("n", 40), n_batches=kw.pop("nb", 8), seed=7,
             scoring=(_grid(nz=20, dz=2.0),), **kw)
    )  # fmt: skip


def test_rounding_bound_counts_step_and_local_pieces(make_config: MakeConfig) -> None:
    """``QuantityResult.rounding_bound`` uses the exact piece counts of the channel's classes:
    edep/dose and the excluded-energy channel include the local point deposits (N_local)."""
    tallies = (TallyRequest("edep", "dose", "edep"), TallyRequest("fl", "dose", "fluence"))
    r = _all_physics_run(make_config, tallies)
    res = Simulation(r.eff.requested).run()
    qs = res.grids[0].quantities
    hpb = res.n_histories // res.n_batches
    n_step = r.n_pieces() / res.n_batches / hpb
    n_loc = r.n_pieces(local=True) / res.n_batches / hpb
    assert n_loc.sum() > 0
    e_ci = r.q_of("edep").numerator
    q = r.plan.channels[e_ci].quantum
    np.testing.assert_allclose(
        qs["edep"].rounding_bound.reshape(-1), (n_step + n_loc) * q / 2.0, rtol=1e-12
    )
    f_ci = r.q_of("fl").numerator
    np.testing.assert_allclose(
        qs["fl"].rounding_bound.reshape(-1) * np.prod(r.eff.requested.scoring[0].spacing_mm),
        n_step * r.plan.channels[f_ci].quantum / 2.0, rtol=1e-12,
    )  # fmt: skip
    ex = qs["edep_excluded_from_let"].rounding_bound.reshape(-1)
    np.testing.assert_allclose(ex, n_loc * q / 2.0, rtol=1e-12)
    assert qs["scoring_pieces_local"].mean.sum() == pytest.approx(n_loc.sum(), rel=1e-12)


def test_a5_integer_conservation_all_physics(make_config: MakeConfig) -> None:
    tallies = (
        TallyRequest("edep", "dose", "edep"),
        TallyRequest("edep_p", "dose", "edep", species=("proton",)),
        TallyRequest("edep_prim", "dose", "edep", generation="primary"),
        TallyRequest("e_step", "dose", "let_d_eps"),
    )
    r = _all_physics_run(make_config, tallies)
    edep = r.part.edep[0]  # the qualified int64 array [B, nvox]
    # per batch: the qualified edep array is bit-identical to the channel (E, both classes)
    ch = r.plan

    def block(ci: int) -> np.ndarray:
        c = ch.channels[ci]
        return r.acc[:, c.offset : c.offset + c.size]

    e_all = block(r.q_of("edep").numerator)
    assert np.array_equal(e_all, edep)
    assert np.array_equal(block(r.q_of("edep_p").numerator), e_all)
    assert np.array_equal(block(r.q_of("edep_prim").numerator), e_all)
    e_step = block(r.q_of("e_step").denominator)
    local = [i for i, c in enumerate(ch.channels) if c.class_mask == CLASS_LOCAL and c.kind == "E"][
        0
    ]
    assert np.array_equal(e_step + block(local), e_all)
    assert block(local).sum() > 0  # the cutoff energy is excluded from the step channel


def test_a6_cauchy_schwarz_and_a14_undefined_bins(make_config: MakeConfig) -> None:
    r = _all_physics_run(make_config, _reqs("let_t", "let_d"), n=48, nb=12)
    lt, ld = r.q_of("let_t"), r.q_of("let_d")
    l_v, ls_v, ls2_v = r.total(lt.denominator), r.total(lt.numerator), r.total(ld.numerator)
    dl, dls, dls2 = (r.delta(i) for i in (lt.denominator, lt.numerator, ld.numerator))
    res = Simulation(r.eff.requested).run()
    q_ld = res.grids[0].quantities["let_d"]
    q_lt = res.grids[0].quantities["let_t"]
    defined = q_ld.defined_mask.reshape(-1) & q_lt.defined_mask.reshape(-1)
    assert defined.sum() >= 5
    with np.errstate(divide="ignore", invalid="ignore"):
        let_t, let_d = ls_v / l_v, ls2_v / ls_v
        delta_v = let_t * (dls2 / ls2_v + 2.0 * dls / ls_v + dl / l_v)
    assert np.all(let_d[defined] >= let_t[defined] * (1.0 - RHO) - delta_v[defined])
    # A14: voxels without steps are NaN with defined_mask False (never 0)
    empty = r.n_pieces() == 0
    assert empty.any()
    for q in (q_ld, q_lt):
        assert np.all(np.isnan(q.mean.reshape(-1)[empty]))
        assert not q.defined_mask.reshape(-1)[empty].any()
        assert np.all(np.isnan(q.std.reshape(-1)[empty]))


def test_a14_reduce_ratio_threshold_and_never_zero() -> None:
    b = 10
    x = np.zeros((b, 4))
    y = np.zeros((b, 4))
    y[:, 0] = 2.0  # defined
    x[:, 0] = 1.0
    y[: min_defined_batches(b), 1] = 1.0  # exactly at the threshold: defined
    x[: min_defined_batches(b), 1] = 3.0
    y[: min_defined_batches(b) - 1, 2] = 1.0  # one below: undefined
    x[: min_defined_batches(b) - 1, 2] = 3.0
    # column 3: no denominator at all (numerator nonzero is impossible in the engine)
    st = reduce_ratio(x, y)
    assert st.defined_mask.tolist() == [True, True, False, False]
    assert st.mean[0] == pytest.approx(0.5)
    assert np.isnan(st.mean[2:]).all() and np.isnan(st.variance_of_mean[2:]).all()
    assert st.n_nonzero.tolist() == [10, 5, 4, 0]
    assert min_defined_batches(3) == 2 and min_defined_batches(20) == 10


def test_delta_method_matches_jackknife_on_synthetic_batches() -> None:
    rng = np.random.default_rng(20351004)
    b = 200
    base = rng.normal(0.0, 1.0, b)
    y = 10.0 * (1.0 + 0.10 * base + 0.05 * rng.normal(size=b))
    x = 4.0 * (1.0 + 0.12 * base + 0.08 * rng.normal(size=b))  # correlated with y
    st = reduce_ratio(x[:, None], y[:, None])
    r_all = x.sum() / y.sum()
    assert st.mean[0] == pytest.approx(r_all, rel=1e-14)
    r_i = (x.sum() - x) / (y.sum() - y)
    jack = (b - 1) / b * np.sum((r_i - r_i.mean()) ** 2)
    assert st.variance_of_mean[0] == pytest.approx(jack, rel=0.03)
    # the covariance term matters: ignoring it is clearly different
    vx = x.var(ddof=1) / b
    vy = y.var(ddof=1) / b
    naive = (vx + r_all**2 * vy) / y.mean() ** 2
    assert abs(naive - jack) / jack > 0.3


# --------------------------------------------------------------------------- A10


def test_a10_lookup_linear_in_s_and_constant_energy_table(
    make_config: MakeConfig, tmp_path: Any
) -> None:
    import hashlib
    from pathlib import Path

    fixture = Path(__file__).resolve().parents[1] / "data" / "synthetic_lookup.json"
    lk = LookupTable.from_file(fixture)  # f = 1 + 0.1 L on 0..100 keV/um
    a, b = 1.0, 0.1
    const = LookupTable(
        "const2", "q", "1", "energy_per_nucleon_mev", "log", np.geomspace(0.5, 200.0, 50),
        {"proton": np.full(50, 2.0)}, "c", "l", "s", True,
    )  # fmt: skip
    tallies = (
        TallyRequest("fe", "dose", "lookup_sum", lookup=lk.name),
        TallyRequest("avg", "dose", "lookup_dose_avg", lookup=lk.name),
        TallyRequest("es", "dose", "let_d_eps"),
        TallyRequest("fe2", "dose", "lookup_sum", lookup="const2"),
    )
    r = _all_physics_run(make_config, tallies, lookups=(lk, const))
    fe, es = r.q_of("fe").numerator, r.q_of("es").numerator
    e_step = r.q_of("es").denominator
    assert r.q_of("avg").denominator == e_step
    lhs = r.total(fe)
    rhs = a * r.total(e_step) + b * r.total(es)
    bound = RHO * lhs + r.delta(fe) + a * r.delta(e_step) + b * r.delta(es)
    assert np.all(np.abs(lhs - rhs) <= bound), np.abs(lhs - rhs).max()
    assert (lhs > 0).sum() >= 5
    # constant table f = 2: FE = 2 E_step bitwise (integer accumulators in quanta of 2^-29 / 2^-30)
    fe2 = r.q_of("fe2").numerator
    c2, ce = r.plan.channels[fe2], r.plan.channels[e_step]
    assert c2.quantum == 2.0 * ce.quantum
    a2 = r.acc[:, c2.offset : c2.offset + c2.size]
    ae = r.acc[:, ce.offset : ce.offset + ce.size]
    assert np.array_equal(a2 * c2.quantum, 2.0 * ae * ce.quantum)
    # provenance: the recorded sha256 is the file hash
    res = Simulation(r.eff.requested).run()
    prov = {p["name"]: p for p in res.lookups}
    assert prov[lk.name]["file_sha256"] == hashlib.sha256(fixture.read_bytes()).hexdigest()
    assert prov[lk.name]["synthetic"] is True and res.valid
    assert res.channel_raw is not None and res.channel_raw.lookup_out_of_domain == 0


# --------------------------------------------------------------------------- A2 (python half)


def _hook_plan(
    make_config: MakeConfig, requests: tuple[TallyRequest, ...], lookups: tuple[LookupTable, ...]
) -> tuple[ChannelPlan, EffectiveConfig]:
    cfg = make_config(energy=60.0, n=40, n_batches=4, scoring=(_grid(nz=4),))
    eff = validate(cfg)
    producible = frozenset((s, g) for s in ("proton", "deuteron") for g in ("primary", "secondary"))
    water = BetheStoppingSource().table(WATER, PROTON)
    plan = compile_channels(
        requests, lookups, cfg.scoring, geometry=eff.geometry, tables=eff.tables, water=water,
        projectile=PROTON, e_cut_mev=2.0, e_hi_mev=60.0, n_histories=40, n_batches=4,
        cpu_workers=1, memory_budget_bytes=2**31, max_steps=1000, scoring_pieces=8,
        producible=producible,
    )  # fmt: skip
    return plan, eff


def test_a2_hook_two_species_two_generations_exact(make_config: MakeConfig) -> None:
    x = np.arange(65.0)
    lk = LookupTable(
        "dy", "q", "1", "let_water_kev_um", "linear", x,
        {"proton": 1.0 + 0.125 * x, "deuteron": 0.5 + 0.25 * x}, "c", "l", "s", True,
    )  # fmt: skip
    reqs: list[TallyRequest] = []
    for tag, sp in (("all", None), ("p", ("proton",)), ("d", ("deuteron",))):
        for gen in ("all", "primary", "secondary"):
            kw: dict[str, Any] = {"species": sp, "generation": gen}
            reqs += [
                TallyRequest(f"lt_{tag}_{gen}", "dose", "let_t", **kw),
                TallyRequest(f"ld_{tag}_{gen}", "dose", "let_d", **kw),
                TallyRequest(f"le_{tag}_{gen}", "dose", "let_d_eps", **kw),
                TallyRequest(f"fe_{tag}_{gen}", "dose", "lookup_sum", lookup="dy", **kw),
                TallyRequest(f"ed_{tag}_{gen}", "dose", "edep", **kw),
            ]
    plan, eff = _hook_plan(make_config, tuple(reqs), (lk,))
    sc = ReferenceChannelScorer(plan, eff.requested.scoring, 4, tables=eff.tables)
    rng = np.random.default_rng(20351004)
    s_mid, k, e_mid, e_dot = 4.0, 0.75, 32.0, 0.5
    sc.set_step_state(s_mid, k, e_mid, e_dot)
    pieces = []
    for _ in range(300):
        pc = {
            "batch": int(rng.integers(0, 4)), "vox": int(rng.integers(0, 4)),
            "l": float(rng.choice([1.0, 0.5, 0.25])), "tau": float(rng.integers(-8, 9) / 8.0),
            "eps": float(rng.integers(1, 64) / 64.0),
            "species": int(rng.integers(0, 2)), "gen": int(rng.integers(0, 2)),
        }  # fmt: skip
        pieces.append(pc)
        sc.score_piece(pc["batch"], 0, pc["vox"], pc["l"], pc["tau"], pc["eps"],
                       pc["species"], pc["gen"], CLASS_STEP)  # fmt: skip
    local = [(1, 0, 0.375, 0, 0), (2, 3, 0.125, 1, 1)]
    for b, v, e, spc, gn in local:
        sc.score_piece(b, 0, v, 0.0, 0.0, e, spc, gn, CLASS_LOCAL)
    assert sc.lookup_ood == 0

    def frac(v: float) -> Fraction:
        return Fraction(v)

    def expected(ci: int) -> np.ndarray:
        c = plan.channels[ci]
        out = [[Fraction(0)] * c.size for _ in range(4)]
        sel = [
            pc for pc in pieces
            if plan.species_match[ci, pc["species"]] and c.gen_lo <= pc["gen"] <= c.gen_hi
        ]  # fmt: skip
        for pc in sel:
            if not c.class_mask & CLASS_STEP:
                continue
            sb = frac(s_mid) + frac(k) * frac(pc["tau"])
            ln = frac(pc["l"])
            val = {
                "E": frac(pc["eps"]), "L": ln, "LS": ln * sb,
                "LS2": ln * (sb * sb + frac(k) ** 2 * ln * ln / 12), "ES": frac(pc["eps"]) * sb,
                "N": Fraction(1),
            }  # fmt: skip
            if c.kind == "FE":
                fv = {0: 1 + frac(0.125) * sb, 1: Fraction(1, 2) + frac(0.25) * sb}[pc["species"]]
                val["FE"] = frac(pc["eps"]) * fv
            out[pc["batch"]][pc["vox"]] += val[c.kind]
        if c.class_mask & CLASS_LOCAL and c.kind in ("E", "N"):
            for b, v, e, spc, gn in local:
                if plan.species_match[ci, spc] and c.gen_lo <= gn <= c.gen_hi:
                    out[b][v] += frac(e) if c.kind == "E" else Fraction(1)
        return np.array(out, dtype=object)

    for ci, c in enumerate(plan.channels):
        want = expected(ci)
        got = sc.acc[:, c.offset : c.offset + c.size]
        for b in range(4):
            for v in range(c.size):
                exp_q = want[b, v] * Fraction(2) ** c.k
                assert exp_q.denominator == 1, (c.kind, "not an exact multiple of the quantum")
                assert int(got[b, v]) == int(exp_q), (c.kind, b, v)

    # closed forms of the ratios (0 tolerance) and recombination of the species partials
    def tot(name: str, den: bool = False) -> np.ndarray:
        qd = next(d for d in plan.quantities if d.name == name)
        c = plan.channels[qd.denominator if den else qd.numerator]
        return sc.acc[:, c.offset : c.offset + c.size].sum(axis=0) * c.quantum

    for gen in ("all", "primary", "secondary"):
        for kind in ("lt", "ld", "le"):
            ln = f"{kind}_all_{gen}"
            num, den = tot(ln), tot(ln, True)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = num / den
            for v in range(4):
                if den[v] > 0:
                    assert ratio[v] == float(Fraction(float(num[v])) / Fraction(float(den[v])))
        for nm in ("lt", "ld", "le"):
            for den in (False, True):
                a_all = tot(f"{nm}_all_{gen}", den)
                a_sum = tot(f"{nm}_p_{gen}", den) + tot(f"{nm}_d_{gen}", den)
                assert np.array_equal(a_all, a_sum), (nm, gen, den)
    assert np.array_equal(tot("ed_all_all"), tot("ed_p_all") + tot("ed_d_all"))
    assert np.array_equal(tot("fe_all_all"), tot("fe_p_all") + tot("fe_d_all"))
    # unproducible requests outside the test producible set raise
    with pytest.raises(Exception, match="not producible"):
        validate(
            replace(
                make_config(), tallies=(TallyRequest("x", "dose", "let_t", species=("deuteron",)),)
            )
        )
    assert species_by_name("deuteron").id == 1


# --------------------------------------------------------------------------- A16 (python) / misc


def test_a16_tallies_do_not_change_the_qualified_outputs(make_config: MakeConfig) -> None:
    """Edep, tallies, counters and traces are bit-identical with and without scoring channels
    (the no-tally path is the code of a524f209: no step precompute, no extra leg walks)."""
    diag = DiagnosticsOptions(track_end_positions=True, trace_histories=2)
    base = _cfg(make_config, (), energy=60.0, n=12, n_batches=4, seed=3, diagnostics=diag)
    with_t = replace(base, tallies=_reqs("let_t", "let_d", "let_d_eps", "fluence"))
    a = validate(base)
    b = validate(with_t)
    ra, rb = run_transport(a), run_transport(b)
    assert ra.channels is None and rb.channels is not None
    for ea, eb in zip(ra.edep_mev, rb.edep_mev, strict=True):
        assert np.array_equal(ea, eb)
    assert ra.tallies == rb.tallies and ra.counters == rb.counters
    assert ra.outside_mev == rb.outside_mev and ra.quantization_mev == rb.quantization_mev
    for k, va in ra.diagnostics["trace"].items():
        assert np.array_equal(va, rb.diagnostics["trace"][k]), k
    assert np.array_equal(ra.diagnostics["end_energy_mev"], rb.diagnostics["end_energy_mev"])


def test_water_row_identity_and_species_interface(make_config: MakeConfig) -> None:
    eff = validate(_cfg(make_config, _reqs("let_t")))
    t = eff.tables
    assert t.has_water and t.water_identity["source"] == "bethe"
    assert eff.summary()["tables"]["water_row"]["material"] == WATER.name
    assert t.s_water(100.0, "proton") == pytest.approx(
        BetheStoppingSource().table(WATER, PROTON).stopping_at(100.0) * 0.1, rel=1e-4
    )
    with pytest.raises(ValueError, match="species"):
        t.s_water(100.0, "alpha")
    assert not validate(make_config()).tables.has_water


def test_reserved_automatic_channel_names(make_config: MakeConfig) -> None:
    from ionmc.errors import UnsupportedCombinationError

    with pytest.raises(UnsupportedCombinationError, match="reserved"):
        validate(_cfg(make_config, (TallyRequest("edep_excluded_from_let", "dose", "edep"),)))
