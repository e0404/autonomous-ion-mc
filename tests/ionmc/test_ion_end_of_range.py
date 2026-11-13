"""End-of-range model of ion histories (V3-008 R2-STOP3): ``continuous_residual_v1`` deposits the
residual energy below the cutoff along the straight continuation over the table CSDA range with
the table's energy density, instead of the point deposit of ``point_v1``.

Low-energy ion starts are outside the validated source domain (>= 50 MeV/u), so the tests widen
``ION_DOMAIN_MEV_PER_U`` by monkeypatch and run with MCS off (the exact-Omega mapping is undefined
there) and, for the profile tests, straggling off: a deterministic straight track from 4.2 MeV/u
down to the cutoff at 2 MeV/u and then the residual segment.
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest
import warp as wp

from ionmc import config as cfgmod
from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    validate,
)
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.lookup import LookupTable
from ionmc.materials import ALUMINIUM, BONE_COMPACT_ICRU, WATER
from ionmc.physics.projectiles import ALPHA, CARBON12, PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport.tally import END_CUTOFF, END_ESCAPED

wp.config.log_level = wp.LOG_WARNING
CELL = 0.01
N_CELLS = 60
E0_U = 4.2


def _config(
    proj,
    model: str = "default",
    backend: str = "warp-cpu",
    e0_u: float = E0_U,
    n: int = 8,
    sigma_u: float = 0.0,
    physics_off: bool = True,
    straggling: bool | None = None,
    n_cells: int = N_CELLS,
    depth: float | None = None,
) -> SimulationConfig:
    depth = CELL * n_cells if depth is None else depth
    e0 = e0_u * proj.a
    return SimulationConfig(
        source=PencilBeamSource(proj, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), e0, sigma_u * proj.a, 0.0),
        geometry=BoxPhantom((-1.0, -1.0, 0.0), (2.0, 2.0, depth), WATER),
        scoring=(ScoringGrid((-1.0, -1.0, 0.0), (2.0, 2.0, depth / n_cells), (1, 1, n_cells)),),
        physics=PhysicsOptions(
            nuclear=False,
            stopping=BetheStoppingSource(),
            straggling=(not physics_off) if straggling is None else straggling,
            multiple_scattering=False,
            end_of_range_model=model,
        ),
        run=RunOptions(backend=backend, precision="float64", seed=5, n_histories=n, n_batches=2),  # type: ignore[arg-type]
    )


@pytest.fixture
def low_domain(monkeypatch):
    for name in ("alpha", "carbon12"):
        monkeypatch.setitem(cfgmod.ION_DOMAIN_MEV_PER_U, name, (1.0, 430.0))


def _profile(cfg: SimulationConfig):
    res = Simulation(cfg).run()
    return res, np.asarray(res.grids[0].energy_mev).reshape(-1)  # per primary


def _table(cfg: SimulationConfig):
    eff = validate(cfg)
    return eff, eff.tables


@pytest.mark.parametrize("proj", [ALPHA, CARBON12], ids=lambda p: p.name)
def test_profile_matches_table_beyond_cutoff(low_domain, proj) -> None:
    """Deterministic track: cumulative deposit beyond depth z equals the table prediction
    ``E_tab(R_cut - (z - z_cut))`` (the last cell with step deposit is excluded), the residual
    sums to E_cut, and the last non-zero cell is at z_cut + R_tab(E_cut) within one cell."""
    cfg = _config(proj)
    eff, tab = _table(cfg)
    res, prof = _profile(cfg)
    assert res.valid
    e_cut = eff.e_cut_mev
    assert e_cut == pytest.approx(2.0 * proj.a)
    # the residual energy of every (identical) history: the last step ends at or below E_cut
    e_res = res.energy_balance.cutoff_mev / cfg.run.n_histories
    assert 0.9 * e_cut < e_res <= e_cut
    r_cut_mm = tab.range_g_cm2(0, e_res) * 10.0 / WATER.density_g_cm3
    r0_mm = tab.range_g_cm2(0, E0_U * proj.a) * 10.0 / WATER.density_g_cm3
    z_cut = r0_mm - r_cut_mm
    # the grid holds the whole residual: the track total is E0 - 0 (no escape)
    assert prof.sum() == pytest.approx(E0_U * proj.a, rel=1e-9)
    assert res.energy_balance.escaped_mev == 0.0
    assert abs(res.energy_balance.closure_residual_mev) < 1e-9 * res.energy_balance.initial_mev
    # last non-zero cell
    last = int(np.nonzero(prof > 0)[0].max())
    assert abs(last - int(math.floor((z_cut + r_cut_mm) / CELL))) <= 1
    # cumulative beyond cell boundaries (skip the cell holding z_cut)
    first = int(math.floor(z_cut / CELL)) + 1
    r_lo = float(tab.r_mass[0, 0])
    e_lo = math.exp(float(tab.ln_e0[0]))
    for j in range(first, last + 1):
        z_b = j * CELL
        r_g = (z_cut + r_cut_mm - z_b) * WATER.density_g_cm3 / 10.0
        e_pred = tab.energy_from_range(0, r_g) if r_g >= r_lo else e_lo * r_g / r_lo
        assert prof[j:].sum() == pytest.approx(e_pred, abs=2e-3 * e_cut)
    assert prof[first:].sum() < e_res


@pytest.mark.parametrize("proj", [ALPHA, CARBON12], ids=lambda p: p.name)
def test_point_v1_puts_residual_in_one_cell(low_domain, proj) -> None:
    """``point_v1`` deposits all of E_cut in the cell of the cutoff point; both models agree on
    everything before the cutoff cell and on the total."""
    cfg_c = _config(proj, "continuous_residual_v1")
    cfg_p = _config(proj, "point_v1")
    eff, tab = _table(cfg_c)
    res_p, prof_p = _profile(cfg_p)
    res_c, prof_c = _profile(cfg_c)
    e_cut = eff.e_cut_mev
    e_res = res_p.energy_balance.cutoff_mev / cfg_p.run.n_histories
    assert 0.9 * e_cut < e_res <= e_cut
    assert res_c.energy_balance.cutoff_mev == pytest.approx(
        res_p.energy_balance.cutoff_mev, rel=1e-9
    )
    z_cut = (tab.range_g_cm2(0, E0_U * proj.a) - tab.range_g_cm2(0, e_res)) * 10.0
    j = int(math.floor(z_cut / CELL))
    # beyond the cutoff cell point_v1 deposits nothing
    assert np.all(prof_p[j + 1 :] == 0.0)
    assert prof_p[j] >= e_res
    np.testing.assert_array_equal(prof_p[:j], prof_c[:j])
    assert prof_p.sum() == pytest.approx(prof_c.sum(), rel=1e-9)
    assert prof_c[j + 1 :].sum() > 0.5 * e_cut


def test_energy_conservation_default_beam(low_domain) -> None:
    """Straggling and MCS on, energy spread, small box: total deposited plus escapes equals the
    source energy with the new default (cutoff and escaped tallies split the residual)."""
    cfg = _config(CARBON12, n=64, physics_off=False, sigma_u=0.05, n_cells=40)
    assert validate(cfg).end_of_range_model == "continuous_residual_v1"
    res, prof = _profile(cfg)
    b = res.energy_balance
    assert res.valid
    assert abs(b.closure_residual_mev) < 1e-9 * b.initial_mev
    grid_total = b.in_grid_mev[0] + b.outside_mev[0] + b.quantization_mev[0]
    assert grid_total == pytest.approx(b.step_deposit_mev + b.cutoff_mev, rel=1e-12)
    assert grid_total + b.escaped_mev + b.truncated_mev == pytest.approx(b.initial_mev, rel=1e-9)


def test_residual_beyond_boundary_is_counted_as_escaped(low_domain) -> None:
    """A phantom ending inside the residual segment: the energy beyond the boundary is tallied as
    escaped, the rest in the grid, total conserved."""
    cfg = _config(ALPHA, depth=0.2, n_cells=20)
    eff, tab = _table(cfg)
    r_cut = tab.range_g_cm2(0, eff.e_cut_mev) * 10.0
    r0 = tab.range_g_cm2(0, E0_U * ALPHA.a) * 10.0
    assert r0 - r_cut < 0.2 < r0  # the boundary cuts the residual segment
    res, prof = _profile(cfg)
    b = res.energy_balance
    assert b.escaped_mev > 0.0
    assert abs(b.closure_residual_mev) < 1e-9 * b.initial_mev
    assert prof.sum() * cfg.run.n_histories == pytest.approx(
        b.initial_mev - b.escaped_mev, rel=1e-9
    )


def test_python_reference_matches_warp(low_domain) -> None:
    """Float64 Warp kernel and Python reference give the same profile (energy per cell to the
    quantum) with the continuous model."""
    _r, a = _profile(_config(ALPHA, backend="warp-cpu", n=4))
    _r, b = _profile(_config(ALPHA, backend="python", n=4))
    np.testing.assert_allclose(a, b, rtol=0, atol=1e-8)


def test_option_recorded_and_validated(low_domain) -> None:
    cfg = _config(ALPHA)
    s = validate(cfg).summary()["physics"]
    assert s["end_of_range_model"] == "default"
    assert s["end_of_range_model_effective"] == "continuous_residual_v1"
    s = validate(_config(ALPHA, "point_v1")).summary()["physics"]
    assert (s["end_of_range_model"], s["end_of_range_model_effective"]) == (
        "point_v1",
        "point_v1",
    )
    with pytest.raises(Exception, match="end_of_range_model"):
        PhysicsOptions(nuclear=False, stopping=BetheStoppingSource(), end_of_range_model="x")
    pcfg = SimulationConfig(
        source=PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 100.0, 0.0, 0.0),
        geometry=BoxPhantom((-1.0, -1.0, 0.0), (2.0, 2.0, 100.0), WATER),
        scoring=(ScoringGrid((-1.0, -1.0, 0.0), (2.0, 2.0, 10.0), (1, 1, 10)),),
        physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource()),
        run=RunOptions(backend="python", precision="float64", seed=1, n_histories=20),
    )
    assert validate(pcfg).end_of_range_model == "point_v1"
    bad = replace(pcfg, physics=replace(pcfg.physics, end_of_range_model="continuous_residual_v1"))
    with pytest.raises(UnsupportedCombinationError, match="ion sources only"):
        validate(bad)


PRE_CHANGE_ALPHA_POINT_PROFILE = [
    0.3618932340759784, 0.3682054771343246, 0.37487980106379837, 0.3819557068636641,
    0.3894737270893529, 0.39748195244465023, 0.40605085832066834, 0.415230609010905,
    0.42510996374767274, 0.4357754747616127, 0.44735017977654934, 0.4599569550482556,
    0.47376923984847963, 0.48898997041396797, 0.5058791041374207, 0.5247872635954991,
    0.5461636778200045, 0.5705741306301206, 0.5988018296193331, 8.227670843945816,
]  # fmt: skip


def test_point_v1_is_bit_identical_to_the_pre_change_code(low_domain) -> None:
    """The depth-deposit array of the alpha run below, taken with the code BEFORE the end-of-range
    option existed (HEAD 76f5bda, warp-cpu float64), is reproduced bit for bit by ``point_v1``
    (checked once the same way for the carbon run, and for a proton run on both backends)."""
    _res, prof = _profile(_config(ALPHA, "point_v1"))
    np.testing.assert_array_equal(prof[:20], np.array(PRE_CHANGE_ALPHA_POINT_PROFILE))
    assert np.all(prof[20:] == 0.0)


# --- review REVIEW-00cd55d4: F2 material stepping, F3 channel semantics, F4 escape record -------

# sha256 of the depth-deposit arrays (float64 bytes) taken with the code of e6f17bc2, BEFORE the
# material stepping, the residual channel semantics and the escape record existed
PRE_STEPPING_SHA256 = {
    ("alpha", "warp-cpu"): "03312af5d317bd228fae04fc51b74a58a955780967c6e6574f5f90a0c7783efe",
    ("alpha", "python"): "c040c4f1ab63cc4ac645122f5f4c3f900c59c88f4f1481c31c6ef2e1adafe0dd",
    ("carbon12", "warp-cpu"): "d9ea216b1064f3641806506f695b6c3748b754a073627ac6d5de86b17deda29e",
    ("carbon12", "python"): "0ca74783ff3581e303df216fb344d92f8bb256a929b287a05b09a7faaa96d9f7",
}
PRE_STEPPING_STRAGGLED_SHA256 = "edad86e9caaa23f72d7f634b49e79f017b20dada0d50d9ae9c509663d978e806"
# the clipped (escape) run changed with the exact voxel walk of REVIEW-b6484e56: the leg cut by
# the world boundary is refined over its own length and the escaped energy is the table energy at
# the range left at the boundary (before: the 32-sub-step grid of the full range, the crossing
# sub-step's energy split linearly by length); pinned with the code of the exact walk
EXACT_WALK_ESCAPE_SHA256 = "edd0c2667039017177a526ee2555082b6784c291dcdad5876090923518ae6712"


def _sha(a: np.ndarray) -> str:
    import hashlib

    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float64).tobytes()).hexdigest()


@pytest.mark.parametrize("backend", ["warp-cpu", "python"])
@pytest.mark.parametrize("proj", [ALPHA, CARBON12], ids=lambda p: p.name)
def test_homogeneous_water_is_bit_identical_to_the_pre_stepping_code(
    low_domain, proj, backend
) -> None:
    _r, prof = _profile(_config(proj, backend=backend, n=8 if backend == "warp-cpu" else 4))
    assert _sha(prof) == PRE_STEPPING_SHA256[(proj.name, backend)]


def test_homogeneous_straggled_and_escape_runs_are_bit_identical(low_domain) -> None:
    _r, prof = _profile(_config(CARBON12, n=64, physics_off=False, sigma_u=0.05, n_cells=40))
    assert _sha(prof) == PRE_STEPPING_STRAGGLED_SHA256
    for backend in ("warp-cpu", "python"):
        res, prof = _profile(_config(ALPHA, backend=backend, depth=0.2, n_cells=20))
        assert _sha(prof) == EXACT_WALK_ESCAPE_SHA256
        # The profile hash above is the bit-identity check; the escaped-energy sum differs by one
        # ulp between platforms (CI run 42941591557 gave 60.76560568262537), so it is compared
        # at 1e-12 relative.
        assert res.energy_balance.escaped_mev == pytest.approx(60.76560568262535, rel=1e-12)


def _residual_e(tab, m: int, r_g: float) -> float:
    r_lo = float(tab.r_mass[m, 0])
    e_lo = math.exp(float(tab.ln_e0[m]))
    return tab.energy_from_range(m, r_g) if r_g >= r_lo else e_lo * r_g / r_lo


def _residual_r(tab, m: int, e: float) -> float:
    e_lo = math.exp(float(tab.ln_e0[m]))
    return tab.range_g_cm2(m, e) if e >= e_lo else float(tab.r_mass[m, 0]) * e / e_lo


def _two_material_config(backend: str, second, frac: float = 0.4, n: int = 4):
    """Alpha straight track in water up to the cutoff; the water/second-material boundary lies
    ``frac`` of the residual range beyond the cutoff point."""
    base = _config(ALPHA, backend=backend, n=n)
    tab = validate(base).tables
    e_cut = validate(base).e_cut_mev
    r0 = tab.range_g_cm2(0, E0_U * ALPHA.a) * 10.0
    r_cut = tab.range_g_cm2(0, e_cut) * 10.0
    z_b = (r0 - r_cut) + frac * r_cut
    geo = VoxelGeometry(
        (-1.0, -1.0, 0.0), (2.0, 2.0, z_b), (1, 1, 2), (WATER, second),
        np.array([[[0, 1]]], dtype=np.int32),
    )  # fmt: skip
    n_cells = int(math.ceil(2.0 * z_b / CELL))
    grid = ScoringGrid((-1.0, -1.0, 0.0), (2.0, 2.0, CELL), (1, 1, n_cells))
    return replace(base, geometry=geo, scoring=(grid,)), z_b


@pytest.mark.parametrize("backend", ["warp-cpu", "python"])
@pytest.mark.parametrize("frac", [0.1, 0.5, 0.9])
@pytest.mark.parametrize("second", [BONE_COMPACT_ICRU, ALUMINIUM], ids=lambda m: m.name)
def test_residual_steps_through_materials(low_domain, second, frac, backend) -> None:
    """F2 / REVIEW-b6484e56: the boundary lies at 10 %, 50 % or 90 % of the residual range. The
    exact voxel walk switches material at the boundary itself: the energy there is the water
    table energy at the range left, the end point is that of the second material's table from
    it (within one scoring cell), the cumulative deposit beyond the boundary follows the second
    material's table, and energy is conserved."""
    cfg, z_b = _two_material_config(backend, second, frac=frac)
    eff = validate(cfg)
    tab = eff.tables
    res, prof = _profile(cfg)
    b = res.energy_balance
    assert res.valid
    assert b.escaped_mev == 0.0
    assert abs(b.closure_residual_mev) < 1e-9 * b.initial_mev
    deposited = prof.sum() * cfg.run.n_histories + b.quantization_mev[0]  # pieces round
    assert deposited == pytest.approx(b.initial_mev, rel=1e-12)
    # the water leg from the cutoff point up to the boundary
    e_res = b.cutoff_mev / cfg.run.n_histories
    r_w = _residual_r(tab, 0, e_res)
    z_cut = (tab.range_g_cm2(0, E0_U * ALPHA.a) - r_w) * 10.0 / WATER.density_g_cm3
    assert z_cut < z_b < z_cut + r_w * 10.0 / WATER.density_g_cm3
    e_b = _residual_e(tab, 0, r_w - (z_b - z_cut) * WATER.density_g_cm3 / 10.0)
    rho2 = second.density_g_cm3
    r2 = _residual_r(tab, 1, e_b)
    z_end = z_b + r2 * 10.0 / rho2
    last = int(np.nonzero(prof > 0)[0].max())
    assert abs(last - int(math.floor(z_end / CELL))) <= 1
    first = int(math.ceil(z_b / CELL))
    for c in range(first, last + 1):
        r_rem = r2 - (c * CELL - z_b) * rho2 / 10.0
        assert prof[c:].sum() == pytest.approx(_residual_e(tab, 1, r_rem), abs=3e-3 * e_res)
    assert prof[first:].sum() <= e_b + 1e-9
    if backend == "warp-cpu":  # python == warp-cpu
        cfg_py, _ = _two_material_config("python", second, frac=frac)
        _r, prof_py = _profile(cfg_py)
        np.testing.assert_allclose(prof, prof_py, rtol=0, atol=1e-8)


def _thin_slab_config(backend: str, second, n: int = 2):
    """Water with one transport layer of ``second`` inside the residual range, half a residual
    sub-step thick (layers and scoring cells of that thickness, so the slab is one cell)."""
    base = _config(ALPHA, backend=backend, n=n)
    eff = validate(base)
    tab = eff.tables
    r0 = tab.range_g_cm2(0, E0_U * ALPHA.a) * 10.0
    r_cut = tab.range_g_cm2(0, eff.e_cut_mev) * 10.0
    t = r_cut / cfgmod.RESIDUAL_SUBSTEPS / 2.0
    k = round((r0 - r_cut + 0.3 * r_cut) / t)
    n_z = int(math.ceil((r0 + r_cut) / t))
    ids = np.zeros((1, 1, n_z), dtype=np.int32)
    ids[0, 0, k] = 1
    geo = VoxelGeometry((-1.0, -1.0, 0.0), (2.0, 2.0, t), (1, 1, n_z), (WATER, second), ids)
    grid = ScoringGrid((-1.0, -1.0, 0.0), (2.0, 2.0, t), (1, 1, n_z))
    return replace(base, geometry=geo, scoring=(grid,)), k, t


@pytest.mark.parametrize("backend", ["warp-cpu", "python"])
def test_thin_voxel_inside_the_residual_range_is_not_skipped(low_domain, backend) -> None:
    """REVIEW-b6484e56: a transport voxel of aluminium thinner than one residual sub-step (half
    of one) inside the residual range receives the table energy loss over its thickness in
    aluminium from the energy at its entry, and the end point is that of the water table from
    the energy left behind it; python == warp-cpu."""
    cfg, k, t = _thin_slab_config(backend, ALUMINIUM)
    tab = validate(cfg).tables
    res, prof = _profile(cfg)
    b = res.energy_balance
    assert res.valid
    assert abs(b.closure_residual_mev) < 1e-9 * b.initial_mev
    deposited = prof.sum() * cfg.run.n_histories + b.quantization_mev[0]  # pieces round
    assert deposited == pytest.approx(b.initial_mev, rel=1e-12)
    e_res = b.cutoff_mev / cfg.run.n_histories
    r_w = _residual_r(tab, 0, e_res)
    z_cut = (tab.range_g_cm2(0, E0_U * ALPHA.a) - r_w) * 10.0
    z_b = k * t
    assert z_cut < z_b and z_b + t < z_cut + r_w * 10.0
    e_b = _residual_e(tab, 0, r_w - (z_b - z_cut) / 10.0)
    rho2 = ALUMINIUM.density_g_cm3
    e_a = _residual_e(tab, 1, _residual_r(tab, 1, e_b) - t * rho2 / 10.0)
    assert prof[k] == pytest.approx(e_b - e_a, rel=1e-5)
    z_end = z_b + t + _residual_r(tab, 0, e_a) * 10.0
    last = int(np.nonzero(prof > 0)[0].max())
    assert abs(last - int(math.floor(z_end / t))) <= 1
    if backend == "warp-cpu":
        cfg_py, _k, _t = _thin_slab_config("python", ALUMINIUM)
        _r, prof_py = _profile(cfg_py)
        np.testing.assert_allclose(prof, prof_py, rtol=0, atol=1e-8)


def test_clipped_residual_records_an_escape(low_domain) -> None:
    """F4: a residual segment cut by the world boundary ends the history as an escape at the
    crossing point with the escaped energy; the records agree with ``escaped_mev`` and between
    the backends."""
    out = {}
    for backend in ("warp-cpu", "python"):
        cfg = _config(ALPHA, backend=backend, depth=0.2, n_cells=20)
        cfg = replace(
            cfg, diagnostics=DiagnosticsOptions(track_end_positions=True, escape_records=True)
        )
        res = Simulation(cfg).run()
        d = res.diagnostics
        b = res.energy_balance
        n = cfg.run.n_histories
        assert len(d["escape_history"]) == n  # every (identical) history leaves
        assert np.all(d["end_code"] == END_ESCAPED)
        assert d["escape_energy_mev"].sum() == pytest.approx(b.escaped_mev, rel=1e-12)
        np.testing.assert_allclose(d["escape_position_mm"][:, 2], 0.2, atol=1e-9)
        np.testing.assert_array_equal(d["end_energy_mev"], d["escape_energy_mev"])
        out[backend] = d
    np.testing.assert_allclose(
        out["warp-cpu"]["escape_energy_mev"], out["python"]["escape_energy_mev"], rtol=1e-12
    )
    # a residual segment that stays inside still ends as a cutoff
    cfg = replace(_config(ALPHA), diagnostics=DiagnosticsOptions(track_end_positions=True))
    d = Simulation(cfg).run().diagnostics
    assert np.all(d["end_code"] == END_CUTOFF)


def _let_config(backend: str, tallies, lookups=()):
    cfg = _config(ALPHA, backend=backend, n=4)
    grid = replace(cfg.scoring[0], name="g")
    return replace(cfg, scoring=(grid,), tallies=tuple(tallies), lookups=tuple(lookups))


@pytest.mark.parametrize("backend", ["warp-cpu", "python"])
def test_let_on_residual_substeps_below_the_floor_is_the_floor_stopping_power(
    low_domain, backend
) -> None:
    """F3: cells crossed only by residual sub-steps below the table floor have
    ``LET_t = S_w(E_lo)`` (constant stopping power, ``k = 0``)."""
    cfg = _let_config(backend, [TallyRequest("lt", "g", "let_t")])
    res = Simulation(cfg).run()
    assert res.valid
    tab = res.effective_config.tables
    e_res = res.energy_balance.cutoff_mev / cfg.run.n_histories
    r_w = _residual_r(tab, 0, e_res)
    n_sub = cfgmod.RESIDUAL_SUBSTEPS
    r_lo = float(tab.r_mass[0, 0])
    k_f = next(k for k in range(n_sub) if r_w * (n_sub - 1 - k) / n_sub < r_lo)
    z_cut = (tab.range_g_cm2(0, E0_U * ALPHA.a) - r_w) * 10.0
    l_w = r_w * 10.0 / n_sub
    z_f, z_end = z_cut + k_f * l_w, z_cut + n_sub * l_w
    cells = [c for c in range(N_CELLS) if c * CELL >= z_f and c * CELL < z_end - 1e-9]
    assert cells
    s_floor = tab.s_water(math.exp(float(tab.ln_e0[0])), "alpha")
    let = res.grids[0].quantities["lt"].mean.reshape(-1)
    for c in cells:
        assert let[c] == pytest.approx(s_floor, rel=1e-6)


@pytest.mark.parametrize("backend", ["warp-cpu", "python"])
def test_lookup_domain_starting_at_the_floor_has_no_out_of_domain(low_domain, backend) -> None:
    """F3: an energy-axis lookup whose domain starts exactly at the table floor is accepted and
    produces zero out-of-domain counts with the continuous model (the sub-steps below the floor
    are evaluated at the floor, clamped by definition)."""
    eff0 = validate(_config(ALPHA))
    e_floor_u = float(eff0.tables.e_min_mev.max()) / ALPHA.a
    x = np.linspace(e_floor_u, 10.0, 11)
    lk = LookupTable(
        "lk", "q", "1", "energy_per_nucleon_mev", "linear", x, {"alpha": 1.0 + x},
        "c", "l", "s", True,
    )  # fmt: skip
    cfg = _let_config(backend, [TallyRequest("ls", "g", "lookup_sum", lookup="lk")], [lk])
    res = Simulation(cfg).run()
    assert res.valid
    assert res.channel_raw is not None and res.channel_raw.lookup_out_of_domain == 0
