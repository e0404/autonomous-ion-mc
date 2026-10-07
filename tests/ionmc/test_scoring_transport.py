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
from ionmc.materials import ALUMINIUM, PMMA, WATER, Material
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource, StoppingTable, build_table
from ionmc.scoring import ScoringGrid, TallyRequest, min_defined_batches, reduce_ratio
from ionmc.simulation import Simulation, ratio_rounding_bound
from ionmc.species import species_by_name
from ionmc.transport.channels import (
    CLASS_LOCAL,
    CLASS_STEP,
    ChannelPlan,
    compile_channels,
    mixed_path_bound_mm,
)
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
        n_v = sum((self.total(ni) for ni in self.plan.piece_count_indices(ci)), np.zeros(1))
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
    the run's own range table). Gating since V3-003D (plan amendment 6, row D3 of the V3-003D
    plan): the range table is the exact integral of the interpolated S, so the transported energy
    and this energy agree to 1e-4. Before V3-003D the trapezoid bias of the range table made it
    informative only."""
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
    # deviation from the table CSDA energy (informative here; gated by the D3 rows)
    print("A4 worst vs transported E(z):", float(np.max(np.abs(let_t[sel] * 2.0 - de) / de)))
    print("A4 worst vs table CSDA E(z):",
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
    footnote 1; ``table_csda`` selects the single-shot table inversion, gating since V3-003D) over
    the
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
    # gating since V3-003D (row D3 of the V3-003D plan; the trapezoid control is in
    # test_range_quadrature.py): the table CSDA energy agrees as well
    print("A4b worst vs transported E(z):", _a4b_worst(r, z_cut, 100.0, 1e9))
    worst_table = _a4b_worst(r, z_cut, 100.0, 1e9, table_csda=True)
    print("A4b worst vs table CSDA E(z):", worst_table)
    assert worst_table <= 1e-4


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
        producible=producible, max_step_mm=2.0,
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


def test_ratio_rounding_bound_is_interval_arithmetic() -> None:
    """The ratio bound covers the corner error that the first-order formula underestimates for a
    small denominator, is +inf when ``Y - b <= 0`` and agrees with first order for large Y."""
    a = np.array([0.5])
    b = np.array([0.5])
    xhat, yhat = np.array([2.0]), np.array([1.0])  # R_hat = 2; true corner (X + a)/(Y - b) = 5
    rhat = xhat / yhat
    first_order = (a + rhat * b) / yhat  # 1.5
    corner = np.abs((xhat + a) / (yhat - b) - rhat)  # 3.0
    got = ratio_rounding_bound(xhat, a, yhat, b, rhat)
    assert corner[0] > first_order[0]  # discriminating: first order is not conservative
    assert got[0] == pytest.approx(corner[0], rel=1e-15) and got[0] >= corner[0]
    # Y_hat - b <= 0: unbounded, never 0 (also for a vanishing numerator)
    for y in (0.5, 0.25, 0.0):
        assert np.isinf(ratio_rounding_bound(xhat, a, np.array([y]), b, np.array([0.0]))[0])
    # large denominator: interval and first-order bounds agree to O(b / Y)
    y = np.array([1.0e4])
    r = xhat / y
    ratio = ratio_rounding_bound(xhat, a, y, b, r)[0] / ((a + r * b) / y)[0]
    assert abs(ratio - 1.0) <= 4.0 * b[0] / y[0]
    # random exhaustive check: the bound covers every point of the interval corners
    rng = np.random.default_rng(1)
    for _ in range(200):
        x, yy = rng.uniform(0, 5), rng.uniform(0.2, 5)
        aa, bb = rng.uniform(0, 0.15), rng.uniform(0, 0.15)
        bnd = ratio_rounding_bound(np.array([x]), np.array([aa]), np.array([yy]),
                                   np.array([bb]), np.array([x / yy]))[0]  # fmt: skip
        for xt in (x - aa, x + aa):
            for yt in (yy - bb, yy + bb):
                assert abs(xt / yt - x / yy) <= bnd * (1 + 1e-12)


def test_piece_counts_follow_the_class_mask_of_each_channel(make_config: MakeConfig) -> None:
    """One helper (``ChannelPlan.piece_count_indices`` / ``piece_counts_for``) gives every rounding
    bound its counts: step-only channels use N, local-only channels N_local, edep both."""
    from ionmc.transport.parity_channels import piece_counts_for

    tallies = (TallyRequest("edep", "dose", "edep"), TallyRequest("lt", "dose", "let_t"))
    r = _all_physics_run(make_config, tallies)
    res = Simulation(r.eff.requested).run()
    plan = r.plan
    n_s, n_l = plan.count_channel(0, CLASS_STEP), plan.count_channel(0, CLASS_LOCAL)
    s_cnt, l_cnt = res.channel_batches(n_s) * 1.0, res.channel_batches(n_l) * 1.0
    assert s_cnt.sum() != l_cnt.sum() and l_cnt.sum() > 0  # the two counts differ
    q = r.q_of
    e_both = q("edep").numerator
    e_loc = next(
        i for i, c in enumerate(plan.channels) if c.kind == "E" and c.class_mask == CLASS_LOCAL
    )
    l_step = q("lt").denominator
    assert plan.piece_count_indices(e_both) == (n_s, n_l)
    assert plan.piece_count_indices(e_loc) == (n_l,)
    assert plan.piece_count_indices(l_step) == (n_s,)
    assert plan.piece_count_indices(n_s) == (n_s,) and plan.piece_count_indices(n_l) == (n_l,)
    assert np.array_equal(piece_counts_for(res, e_both), s_cnt + l_cnt)
    assert np.array_equal(piece_counts_for(res, e_loc), l_cnt)
    assert np.array_equal(piece_counts_for(res, l_step), s_cnt)
    assert not np.array_equal(piece_counts_for(res, e_loc), piece_counts_for(res, l_step))
    # the test helper and the public bound agree with it
    np.testing.assert_allclose(
        r.delta(e_both), (r.n_pieces() + r.n_pieces(local=True)) * plan.channels[e_both].quantum / 2
    )
    np.testing.assert_allclose(
        r.delta(e_loc), r.n_pieces(local=True) * plan.channels[e_loc].quantum / 2
    )


# --------------------------------------------------------------------------- ramp envelope


class SteepSource:
    """Synthetic stopping source ``S_mass = c E^-1.5`` (log-log slope -1.5 everywhere): the LET
    ramp ``S_mid (1 -+ |gamma| f_E / 2)`` of a long final step can reach zero or go negative."""

    name = "steep-s"

    def table(self, material: Material, projectile: Any) -> StoppingTable:
        e = np.geomspace(1.0, 500.0, 400)
        c = 10.0
        return build_table(
            projectile, material, e, c * e**-1.5, e[0] ** 2.5 / (2.5 * c), {"source": "s"}
        )


def test_ramp_envelope_rejects_a_steep_table_at_validate(make_config: MakeConfig) -> None:
    from ionmc.errors import UnsupportedCombinationError

    cfg = _cfg(make_config, _reqs("let_t"), energy=60.0, scoring=(_grid(),), stopping=SteepSource())
    with pytest.raises(UnsupportedCombinationError, match="not guaranteed positive"):
        validate(cfg)
    # without LET channels the qualified path does not look at the ramp
    validate(replace(cfg, tallies=()))


def test_ramp_envelope_covers_cutoff_crossing_steps(
    make_config: MakeConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A large cutoff makes steps cross it: observed E_mid reaches below E_cut (but stays above
    E_cut/2), every observed ramp value lies inside the recorded envelope, the adaptive-channel
    accumulators stay non-negative and no ramp guard fires."""
    seen: list[tuple[float, float, float]] = []
    orig = ReferenceChannelScorer.begin_step

    def spy(self: ReferenceChannelScorer, e_mid: float, de_mean: float, s_act: float) -> None:
        orig(self, e_mid, de_mean, s_act)
        half = 0.5 * self.k * s_act
        seen.append((e_mid, self.s_mid - abs(half), self.s_mid + abs(half)))

    monkeypatch.setattr(ReferenceChannelScorer, "begin_step", spy)
    e_cut = 4.0
    cfg = _cfg(
        make_config, _reqs("let_t", "let_d", "let_d_eps"), energy=12.0, n=300, n_batches=6,
        seed=5, e_cut=e_cut, scoring=(_grid(nz=30, dz=0.5),),
    )  # fmt: skip
    r = Run(cfg)
    assert len(seen) > 100
    e_mid = np.array([x[0] for x in seen])
    assert e_mid.min() < e_cut and e_mid.min() >= 0.5 * e_cut  # cutoff-crossing steps exist
    lo = np.array([x[1] for x in seen])
    hi = np.array([x[2] for x in seen])
    b = r.plan.bounds
    assert lo.min() >= b["S_w_min_mev_per_mm"] > 0.0
    assert hi.max() <= b["S_w_max_mev_per_mm"]
    assert (r.acc >= 0).all()
    assert r.raw.channels is not None and r.raw.channels.lookup_out_of_domain == 0


def test_nonpositive_ramp_is_flagged_at_runtime(make_config: MakeConfig) -> None:
    """The hook never scores a nonpositive ramp silently: a piece with ``S_bar <= 0`` increments
    the out-of-domain counter (python hook; the kernel has the same guard on the tally row), which
    invalidates the result. A state inside the envelope does not."""
    plan, eff = _hook_plan(make_config, _reqs("let_t"), ())
    sc = ReferenceChannelScorer(plan, eff.requested.scoring, 4, tables=eff.tables)
    sc.set_step_state(4.0, -3.0, 30.0, 0.5)  # S_bar = 4 - 3 tau
    sc.score_piece(0, 0, 0, 1.0, 0.5, 0.1, 0, 0, CLASS_STEP)  # tau = 0.5: S_bar = 2.5 > 0
    assert sc.lookup_ood == 0
    sc.score_piece(0, 0, 0, 1.0, 2.0, 0.1, 0, 0, CLASS_STEP)  # tau = 2: S_bar = -2 <= 0
    assert sc.lookup_ood == 1


class DropSource:
    """Water: constant mass stopping power 10 MeV cm2/g (gamma = 0). Any other material: 5 MeV cm2/g
    above ``e_drop`` MeV/u and falling as ``(E/e_drop)^3`` below it (log-log slope 3, ordinary above
    the cutoff, sharply lower stopping below it, so ``S_w / S_m`` grows by 8 between E_cut and
    E_cut/2)."""

    name = "drop-s"

    def __init__(self, e_drop: float) -> None:
        self.e_drop = e_drop

    def table(self, material: Material, projectile: Any) -> StoppingTable:
        e = np.geomspace(1.0, 500.0, 400)
        if material.name == WATER.name:
            return build_table(projectile, material, e, np.full_like(e, 10.0), e[0] / 10.0, {})
        s = 5.0 * np.minimum(1.0, (e / self.e_drop) ** 3)
        # range integral dE / S on the grid (exact log-log quadrature in build_table), start
        # range from
        # the first node
        r0 = e[0] / (4.0 * s[0])
        return build_table(projectile, material, e, s, r0, {"source": "drop"})


def test_ratio_bound_covers_the_cutoff_crossing_midpoint_domain(make_config: MakeConfig) -> None:
    """Review f6bd2d53: the stopping-ratio candidate of the LS bound must be evaluated down to
    E_cut/2. Heterogeneous phantom (water, then a material whose stopping power drops sharply below
    the cutoff): the recorded r_max exceeds the value of the old domain [E_cut, E_hi], and the
    recorded LS/LS2 bounds cover every observed per-history LS/LS2 (one history per batch)."""
    from ionmc.transport.channels import _stopping_ratio_max

    e_cut, e_hi = 4.0, 12.0
    nz = 30
    mat_index = np.concatenate([np.zeros(10), np.ones(nz - 10)]).astype(np.int32)
    mat_index = mat_index.reshape(1, 1, nz)
    geo = VoxelGeometry(
        origin_mm=(-30.0, -30.0, 0.0), spacing_mm=(60.0, 60.0, 1.0), shape=(1, 1, nz),
        materials=(WATER, ALUMINIUM),
        material_index=mat_index,
    )  # fmt: skip
    cfg = _cfg(
        make_config, _reqs("let_t", "let_d"), energy=e_hi, n=240, n_batches=240, seed=11,
        e_cut=e_cut, geometry=geo, stopping=DropSource(e_drop=e_cut),
        scoring=(_grid(nz=nz, dz=1.0),), straggling=False, mcs=False,
    )  # fmt: skip
    r = Run(cfg)
    b = r.plan.bounds
    water = DropSource(e_cut).table(WATER, PROTON)
    amp = b["ramp_amplification"]
    rho_min = float(geo.densities_g_cm3()[geo.material_index == 1].min())
    old = _stopping_ratio_max(r.eff.tables, 1, water, rho_min, e_cut, e_hi, 1, amp)
    new = _stopping_ratio_max(r.eff.tables, 1, water, rho_min, 0.5 * e_cut, e_hi, 1, amp)
    assert new > 1.5 * old  # discriminating: the old domain underestimates the ratio
    assert b["r_max"] >= new * (1 - 1e-12)
    # observed per-history sums (one history per batch) stay below the recorded bounds
    ls_ci = r.q_of("let_t").numerator
    ls2_ci = r.q_of("let_d").numerator
    for ci, key in ((ls_ci, "B_LS_mev"), (ls2_ci, "B_LS2")):
        c = r.plan.channels[ci]
        per_hist = r.acc[:, c.offset : c.offset + c.size].sum(axis=1) * c.quantum
        assert per_hist.max() <= b[key], (key, per_hist.max(), b[key])
        assert per_hist.max() > 0.0
    assert (r.acc >= 0).all()


class KinkSource:
    """Water: constant 10 MeV cm2/g except a narrow V-shaped dip of ``ln S`` (``depth``) around
    ``e_kink`` that lies entirely *inside one runtime-grid bin* (source nodes at ``e_kink exp(+-w)``
    with ``w`` a fraction of the bin). The runtime water row samples the source at its own nodes,
    so it is flat 10 there (zero slope in every bin: no ramp amplification anywhere), while the
    source interpolant dips below it. Other materials: 5 MeV cm2/g above ``e_drop`` and falling as
    ``(E/e_drop)^power`` below it (the ratio ``S_w/S_m`` grows steeply towards low energies)."""

    name = "kink-s"

    def __init__(self, e_kink: float, w: float, e_drop: float, depth: float = 0.5,
                 power: float = 8.0, skip: tuple[float, float] = (0.0, 0.0)) -> None:  # fmt: skip
        self.e_kink, self.w, self.e_drop = e_kink, w, e_drop
        self.depth, self.power, self.skip = depth, power, skip

    def table(self, material: Material, projectile: Any) -> StoppingTable:
        e = np.geomspace(1.0, 500.0, 40)
        if material.name != WATER.name:
            s = 5.0 * np.minimum(1.0, (e / self.e_drop) ** self.power)
            return build_table(projectile, material, e, s, e[0] / (4.0 * s[0]), {"source": "drop"})
        e = e[(e <= self.skip[0]) | (e >= self.skip[1])]  # keep the runtime bin of the dip clean
        v = self.e_kink * np.exp([-self.w, 0.0, self.w])
        e = np.unique(np.concatenate([e, v]))
        s = np.full_like(e, 10.0)
        s[e == v[1]] = 10.0 * np.exp(-self.depth)
        return build_table(projectile, material, e, s, e[0] / (4.0 * s[0]), {"source": "kink"})


def _ratio_max_on_runtime_rows(tb: Any, material: int, rho_min: float, e_lo: float, e_hi: float,
                               water: StoppingTable | None = None) -> float:  # fmt: skip
    """Independent re-computation of the stopping-ratio bound. ``water=None``: S_w from the runtime
    row (the fixed computation); otherwise S_w from the source table ``water.stopping_at`` with the
    same candidate set and the same ``(1 + amp)`` as the pre-fix code at 952f248 (the legacy
    bound)."""
    row = tb.water_ln_s_mass
    ln_ew = tb.water_ln_e0 + np.arange(row.size) / tb.water_inv_dln_e
    ln_em = tb.ln_e0[material] + np.arange(tb.ln_s_mass[material].size) / tb.inv_dln_e[material]
    cand = np.exp(np.concatenate([ln_ew, ln_em]))
    pts = np.concatenate([[e_lo, e_hi], cand[(cand > e_lo) & (cand < e_hi)]])
    s_m = np.exp(np.interp(np.log(pts), ln_em, tb.ln_s_mass[material]))
    if water is None:
        s_w = np.exp(np.interp(np.log(pts), ln_ew, row))
        rho_w = tb.water_density_g_cm3
    else:
        s_w = water.stopping_at(pts)
        rho_w = water.material.density_g_cm3
    gam = np.abs(np.diff(row)) * tb.water_inv_dln_e  # amplification |gamma| f_E / 2, f_E = 2
    t = (np.log(pts) - tb.water_ln_e0) * tb.water_inv_dln_e
    i_hi = np.clip(np.floor(t + 1e-9).astype(int), 0, gam.size - 1)
    i_lo = np.clip(np.ceil(t - 1e-9).astype(int) - 1, 0, gam.size - 1)
    amp = np.maximum(gam[i_hi], gam[i_lo])
    return float(((1.0 + amp) * s_w * rho_w / (s_m * rho_min)).max())


def test_ratio_bound_uses_the_runtime_water_row_not_the_source_interpolant(
    make_config: MakeConfig,
) -> None:
    """Review fb43008e / f1baadf0: the stopping-ratio bound must take S_w from the runtime water
    row. A narrow dip of the source table lies inside one runtime bin at E_cut/2 (all runtime bins
    flat, so amp = 0 everywhere and no amplified neighbour can compensate): the legacy bound
    (source interpolant, same candidates and amplification as before the fix) is lower than the
    fixed one by a stated margin, and the compiled r_max equals the fixed computation.

    The observed per-history LS stays far below even the legacy-derived bound in a CI-sized run
    (the bound is conservative by design), so the gating discriminator is the legacy-vs-fixed
    comparison (b); the recorded bounds are still checked to cover the observed LS/LS2."""
    nz, e_hi = 30, 12.0
    mat_index = np.concatenate([np.zeros(10), np.ones(nz - 10)]).astype(np.int32)
    geo = VoxelGeometry(
        origin_mm=(-30.0, -30.0, 0.0), spacing_mm=(60.0, 60.0, 1.0), shape=(1, 1, nz),
        materials=(WATER, ALUMINIUM), material_index=mat_index.reshape(1, 1, nz),
    )  # fmt: skip

    def config(src: Any, e_cut: float) -> SimulationConfig:
        return _cfg(
            make_config, _reqs("let_t", "let_d"), energy=e_hi, n=240, n_batches=240, seed=11,
            e_cut=e_cut, geometry=geo, stopping=src, scoring=(_grid(nz=nz, dz=1.0),),
            straggling=False, mcs=False,
        )  # fmt: skip

    t0 = validate(
        config(KinkSource(2.0, 0.0, 4.0, depth=0.0), 4.0)
    ).tables  # the runtime grid (energy range)
    assert t0 is not None and t0.water_ln_s_mass is not None
    nodes = np.exp(t0.water_ln_e0 + np.arange(t0.water_ln_s_mass.size) / t0.water_inv_dln_e)
    k = int(np.searchsorted(nodes, 2.0)) - 1
    dln = 1.0 / t0.water_inv_dln_e
    e_kink = float(np.sqrt(nodes[k] * nodes[k + 1]))  # middle of the bin; E_cut/2 sits on the dip
    src = KinkSource(e_kink, 0.2 * dln, e_drop=2.0 * e_kink, skip=(nodes[k], nodes[k + 1]))
    r = Run(config(src, 2.0 * e_kink))
    tb, b = r.eff.tables, r.plan.bounds
    assert np.allclose(tb.water_ln_s_mass, np.log(10.0), atol=1e-9)  # runtime row is flat
    assert b["ramp_amplification"] == pytest.approx(0.0, abs=1e-9)
    water = src.table(WATER, PROTON)
    e_lo = 0.5 * (2.0 * e_kink)
    s_run = float(np.exp(np.interp(np.log(e_lo), np.log(nodes), tb.water_ln_s_mass)))
    assert s_run > 1.5 * float(water.stopping_at(e_lo))  # the source interpolant dips at E_cut/2
    rho_al = float(geo.densities_g_cm3()[geo.material_index == 1].min())
    fixed = _ratio_max_on_runtime_rows(tb, 1, rho_al, e_lo, e_hi)
    legacy = _ratio_max_on_runtime_rows(tb, 1, rho_al, e_lo, e_hi, water=water)
    assert b["r_max"] == pytest.approx(max(fixed, 1.0), rel=1e-9)  # (a) compiled = fixed (runtime)
    assert fixed > 1.04 * legacy  # (b) the legacy bound is lower by more than 4 %
    for ci, key in ((r.q_of("let_t").numerator, "B_LS_mev"), (r.q_of("let_d").numerator, "B_LS2")):
        c = r.plan.channels[ci]
        per_hist = r.acc[:, c.offset : c.offset + c.size].sum(axis=1) * c.quantum
        assert 0.0 < per_hist.max() <= b[key]
    assert (r.acc >= 0).all()


class SwapSource:
    """Water: constant 10 MeV cm2/g (the LET row). Aluminium and PMMA have *opposite* energy
    dependence of the linear stopping power (MeV/mm): PMMA 1 above ``e0`` and 10 below, aluminium
    10 above and 1 below. Each homogeneous CSDA range is short (about 0.55 of the mixed path), but
    a path through PMMA first and aluminium after follows the lowest linear stopping power at every
    energy."""

    name = "swap-s"

    def __init__(self, e0: float) -> None:
        self.e0 = e0

    def table(self, material: Material, projectile: Any) -> StoppingTable:
        e = np.geomspace(1.0, 500.0, 400)
        if material.name == WATER.name:
            return build_table(projectile, material, e, np.full_like(e, 10.0), e[0] / 10.0, {})
        high = e > self.e0
        lin = np.where(high == (material.name == PMMA.name), 1.0, 10.0)  # MeV/mm
        s = lin * 10.0 / material.density_g_cm3  # MeV cm2/g
        return build_table(projectile, material, e, s, e[0] / (4.0 * s[0]), {"source": "swap"})


def test_path_bound_covers_a_path_that_changes_material_with_energy(
    make_config: MakeConfig,
) -> None:
    """Review 0352ebfd: B_L = 1.25 max_m R_m(E_hi)/rho_min,m (legacy) is not a bound for a path that
    changes material as the energy falls: PMMA (1 MeV/mm above 6 MeV) then aluminium (1 MeV/mm
    below 6 MeV) is followed at the lowest linear stopping power at every energy. The observed
    per-history path exceeds the legacy bound, never the compiled B_L (the integral of
    1/min_m S_lin,m(E), capped by max_steps * max_step); accumulators stay non-negative and no
    overflow counter fires."""
    nz, e_hi, e0 = 30, 12.0, 6.0
    mat_index = np.concatenate([np.zeros(6), np.ones(nz - 6)]).astype(np.int32)  # 6 mm PMMA, Al
    geo = VoxelGeometry(
        origin_mm=(-30.0, -30.0, 0.0), spacing_mm=(60.0, 60.0, 1.0), shape=(1, 1, nz),
        materials=(PMMA, ALUMINIUM), material_index=mat_index.reshape(1, 1, nz),
    )  # fmt: skip
    cfg = _cfg(
        make_config, _reqs("fluence", "let_t"), energy=e_hi, n=4, n_batches=4, seed=3, e_cut=2.0,
        geometry=geo, stopping=SwapSource(e0), scoring=(_grid(nz=nz, dz=1.0),),
        straggling=False, mcs=False, max_step=0.5,
    )  # fmt: skip
    r = Run(cfg)
    b, tb = r.plan.bounds, r.eff.tables
    dens = geo.densities_g_cm3()
    legacy = 1.25 * max(
        tb.range_g_cm2(m, e_hi) * 10.0 / float(dens[geo.material_index == m].min()) for m in (0, 1)
    )
    c = r.plan.channels[r.q_of("fluence").numerator]
    per_hist = r.acc[:, c.offset : c.offset + c.size].sum(axis=1) * c.quantum
    observed = float(per_hist.max())
    print(f"B_L legacy {legacy:.3f} mm, new {b['B_L_mm']:.3f} mm, observed {observed:.3f} mm")
    assert observed > legacy  # (b) the legacy bound is exceeded at runtime
    assert observed <= b["B_L_mm"]  # (a) the compiled bound covers it
    assert (r.acc >= 0).all() and r.raw.counters["accumulator_overflow"] == 0  # (c)
    assert r.raw.channels is not None and r.raw.channels.path_bound_exceeded == 0
    # the recorded bounds are the CSDA-envelope formulas (the truncation bound is informative)
    s_max = b["S_w_max_mev_per_mm"]
    assert b["B_LS_mev"] == pytest.approx(s_max * b["B_L_mm"], rel=1e-12)
    assert b["B_LS2"] == pytest.approx(s_max**2 * b["B_L_mm"], rel=1e-12)
    assert b["B_ES"] == pytest.approx(s_max * 12.0, rel=1e-12)
    assert b["B_L_truncation_mm"] == pytest.approx(r.eff.max_steps * 0.5 * (1.0 + 1e-6), rel=1e-12)
    assert r.plan.path_bound_mm == b["B_L_mm"]


@pytest.mark.parametrize("straggling", [True, False])
def test_bounds_are_the_checked_formulas_with_and_without_straggling(
    make_config: MakeConfig, straggling: bool
) -> None:
    """B_L = 1.25 int dE / min S_lin (heterogeneous envelope), B_LS = S_bar_max B_L,
    B_LS2 = S_bar_max^2 B_L, B_ES = S_bar_max E_hi, whatever the straggling setting."""
    cfg = _cfg(
        make_config, _reqs("fluence", "let_t", "let_d", "let_d_eps"), energy=12.0, n=8,
        n_batches=4, straggling=straggling, scoring=(_grid(nz=10, dz=1.0),),
    )  # fmt: skip
    plan = validate(cfg).channels
    assert plan is not None
    b, tb = plan.bounds, validate(cfg).tables
    rho = {0: float(WATER.density_g_cm3)}
    expected = 1.25 * mixed_path_bound_mm(tb, rho, 12.0)
    assert b["B_L_mm"] == pytest.approx(expected, rel=1e-12) and plan.path_bound_mm == b["B_L_mm"]
    s_max = b["S_w_max_mev_per_mm"]
    assert b["B_LS_mev"] == pytest.approx(s_max * expected, rel=1e-12)
    assert b["B_LS2"] == pytest.approx(s_max**2 * expected, rel=1e-12)
    assert b["B_ES"] == pytest.approx(s_max * 12.0, rel=1e-12)
    assert b["B_L_truncation_mm"] > 5.0 * expected  # rigorous but coarse: informative only


def test_low_loss_sampler_violates_the_path_bound_and_invalidates(
    make_config: MakeConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Review 36271bf4 failure mode: a straggling sampler that returns near-zero losses lets a
    history travel far beyond 1.25 x its CSDA path. The bound B_L is therefore checked at
    runtime: such histories are counted in ``path_bound_exceeded`` and the result is invalid
    (fail closed). The ordinary sampler never trips it (1e3 histories)."""
    from ionmc.simulation import _assemble
    from ionmc.transport import reference as ref_mod

    def config(n: int) -> SimulationConfig:
        return _cfg(
            make_config, _reqs("fluence", "let_t"), energy=12.0, n=n, n_batches=4, seed=9,
            scoring=(_grid(nz=30, dz=1.0),), straggling=True, mcs=False,
        )  # fmt: skip

    normal = Run(config(1000))
    assert normal.raw.channels is not None and normal.raw.channels.path_bound_exceeded == 0
    assert _assemble(normal.eff, normal.raw).valid

    orig_init = ref_mod._Reference.__init__

    def patched(self: Any, eff: Any) -> None:
        orig_init(self, eff)
        self.straggle_attempt = lambda mean, var, *u: (0.01 * mean, True)  # near-zero loss

    monkeypatch.setattr(ref_mod._Reference, "__init__", patched)
    r = Run(config(12))
    b = r.plan.bounds
    c = r.plan.channels[r.q_of("fluence").numerator]
    per_hist_budget = b["B_L_mm"]
    # observed path of the near-zero-loss histories: the truncation bound, far above B_L
    assert r.raw.counters["step_truncation"] > 0
    assert b["B_L_truncation_mm"] > per_hist_budget
    assert r.raw.channels is not None and r.raw.channels.path_bound_exceeded > 0
    path = r.acc[:, c.offset : c.offset + c.size].sum() * c.quantum / 12.0  # mean per history
    assert path > per_hist_budget  # the CSDA-derived bound is exceeded: it is not a hard bound
    assert not _assemble(r.eff, r.raw).valid  # fail closed: no accumulator value is used
