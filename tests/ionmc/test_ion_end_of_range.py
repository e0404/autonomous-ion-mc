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
from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig, validate
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import ALPHA, CARBON12, PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource

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
