"""V3-005A work packages C10 to C12: the per-history particle stack, thinning of the primary's
nuclear candidates and the nuclear events of the Python reference (decision 0041 sections 2 to 5,
rows P4, V3-CI, V2 smoke and the fail-closed counters of the acceptance plan).

The data-backed tests use a built nuclear table of the cache (``IONMC_CACHE_DIR``) and skip
without one (``IONMC_REQUIRE_DATA=1`` makes that a failure)."""

from __future__ import annotations

import math
import os
from collections.abc import Callable
from dataclasses import replace
from typing import Any

import numpy as np
import pytest

import ionmc.nuclear.events as events
import ionmc.transport.reference as reference
from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
)
from ionmc.data import cache
from ionmc.errors import CounterOverflowError, TransportLimitError
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.nuclear.events import ScalarEvent
from ionmc.nuclear.tables import MaterialNuclear, NuclearTable, NuclearTableError, majorant
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import NuclearEnergyBalance, Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport.reference import END_NUCLEAR, _Reference

MakeConfig = Callable[..., SimulationConfig]
DEEP_WATER = BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, 200.0), WATER)
GRID = (ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 100)),)


def _table_id() -> str:
    try:
        cdir = cache.resolve_cache_dir(None)
    except Exception:
        cdir = None
    found: list[tuple[bool, str]] = []
    if cdir is not None:
        for p in sorted((cdir / "derived").glob("nuclear-proton-*.json")):
            tid = p.stem.removeprefix("nuclear-proton-")
            try:
                info = NuclearTable.load(None, tid).info
            except NuclearTableError:
                continue
            if "transport_path_bound_terms" in info:
                found.append((info["options"]["diagnostic_nodes_mev"] is None, tid))
    if found:
        return max(found)[1]
    if os.environ.get("IONMC_REQUIRE_DATA") == "1":
        pytest.fail("no loadable nuclear table in the cache")
    pytest.skip("no built nuclear table in the cache (IONMC_CACHE_DIR)")


@pytest.fixture(scope="module")
def tid() -> str:
    return _table_id()


def _nuc(cfg: SimulationConfig, tid: str, **run: Any) -> SimulationConfig:
    return replace(
        cfg,
        physics=replace(cfg.physics, nuclear=True, nuclear_table_id=tid),
        run=replace(cfg.run, **run),
    )


def _scaled_rows(monkeypatch: pytest.MonkeyPatch, factor: float) -> None:
    """Multiply every Sigma (and majorant) of ``material_rows`` by ``factor``: many events."""
    orig = NuclearTable.material_rows

    def scaled(self: NuclearTable, material: Any, f_e: float = 0.02) -> MaterialNuclear:
        m = orig(self, material, f_e)
        return replace(
            m,
            sigma_mass_cm2_g=m.sigma_mass_cm2_g * factor,
            sigma_hat_window=m.sigma_hat_window * factor,
            sigma_hat_end=m.sigma_hat_end * factor,
            cum_sigma_mass_cm2_g=m.cum_sigma_mass_cm2_g * factor,
        )

    monkeypatch.setattr(NuclearTable, "material_rows", scaled)


def _expected_event_probability(sim: Simulation, e0: float) -> float:
    """``1 - exp(-int Sigma(E) / S_mass(E) dE)`` along the deterministic CSDA path (the tables of
    the run), from the cutoff to ``e0`` (trapezoid on 4000 log-spaced nodes)."""
    eff = sim.effective
    assert eff.nuclear is not None
    rows = eff.nuclear.rows[0]
    e = np.geomspace(eff.requested.physics.e_cut_mev, e0, 4000)
    f = np.array([rows.sigma_at(x) / eff.tables.stopping_mass(0, x) for x in e])
    return 1.0 - math.exp(-float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(e))))


def _config(
    tid: str,
    *,
    energy: float,
    n: int,
    seed: int = 20351004,
    mcs: bool = True,
    straggling: bool = True,
    max_step: float = 2.0,
    e_cut: float = 2.0,
    tallies: tuple[TallyRequest, ...] = (),
    allow_invalid: bool = False,
    track_end: bool = True,
) -> SimulationConfig:
    return SimulationConfig(
        source=PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), energy, 0.0, 0.0),
        geometry=DEEP_WATER,
        scoring=GRID,
        tallies=tallies,
        physics=PhysicsOptions(
            nuclear=True,
            stopping=BetheStoppingSource(),
            straggling=straggling,
            multiple_scattering=mcs,
            e_cut_mev=e_cut,
            max_step_mm=max_step,
            nuclear_table_id=tid,
        ),
        run=RunOptions(
            backend="python",
            precision="float64",
            seed=seed,
            n_histories=n,
            n_batches=2,
            max_steps=None,
            allow_invalid_result=allow_invalid,
        ),
        diagnostics=DiagnosticsOptions(track_end_positions=track_end),
    )


TALLIES = (
    TallyRequest("fl_sec", "dose", "fluence", generation="secondary"),
    TallyRequest("fl_prim", "dose", "fluence", species=("proton",), generation="primary"),
    TallyRequest("fl_d", "dose", "fluence", species=("deuteron",), generation="secondary"),
)


# ---- V3-CI, V2 smoke --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def ci_run(tid: str) -> tuple[Any, Simulation]:
    sim = Simulation(_config(tid, energy=150.0, n=2000, tallies=TALLIES))
    return sim.run(), sim


def test_v3_ci_energy_balance_and_counters(ci_run: tuple[Any, Simulation]) -> None:
    """Row V3-CI: the amended balance closes to 1e-12 relative, the per-grid identity with the
    nuclear-local destination closes, every fail-closed counter is 0 (2e3 histories, 150 MeV)."""
    res, _ = ci_run
    b = res.energy_balance
    print(f"V3-CI closure {b.relative_residual:.3e}, grid {b.grid_relative_residual(0):.3e}")
    print("counters", res.counters.as_dict())
    assert isinstance(b, NuclearEnergyBalance) and res.valid
    assert b.relative_residual <= 1e-12 and b.grid_relative_residual(0) <= 1e-12
    assert not any(res.counters.as_dict().values())
    n = b.nuclear_mev
    assert n["nuclear_local"] > 0 and n["nuclear_escaped_neutron"] > 0
    assert 0 < n["nuclear_alpha_local"] <= n["nuclear_local"]
    assert b.unaccounted_mev == 0.0 and b.truncated_mev == 0.0


def test_v2_smoke_primary_attenuation(ci_run: tuple[Any, Simulation]) -> None:
    """CI smoke of V2: the fraction of primaries ending in a nuclear event equals
    1 - exp(-int Sigma dR) along the CSDA path within 4 sigma (binomial, 2e3 histories)."""
    res, sim = ci_run
    code = res.diagnostics["end_code"]
    frac = float(np.mean(code == END_NUCLEAR))
    expected = _expected_event_probability(sim, 150.0)
    sigma = math.sqrt(expected * (1.0 - expected) / code.size)
    print(f"V2 smoke: measured {frac:.5f} expected {expected:.5f}, ratio {frac / expected:.4f}"
          f" +- {sigma / expected:.4f}")  # fmt: skip
    assert abs(frac - expected) <= 4.0 * sigma


def test_secondaries_visible_in_generation_resolved_fluence(ci_run: tuple[Any, Simulation]) -> None:
    res, sim = ci_run
    plan = sim.effective.channels
    assert plan is not None
    idx = {q.name: q.numerator for q in plan.quantities}
    sec = res.channel_batches(idx["fl_sec"]).sum()
    prim = res.channel_batches(idx["fl_prim"]).sum()
    assert sec > 0 and prim > 0
    assert res.channel_batches(idx["fl_d"]).sum() > 0  # deuterons are produced and scored


# ---- stack, cutoffs and event geometry ---------------------------------------------------------
def test_secondary_cutoffs_and_generation_labels(tid: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """Many events (Sigma x 40): proton secondaries stop at E_cut = 2 MeV, deuterons at
    E_cut,d = 4 MeV; the nuclear-local deposit has species 64 and generation 1."""
    _scaled_rows(monkeypatch, 40.0)
    seen: list[tuple[int, int, float]] = []
    local: list[tuple[int, int]] = []
    orig = _Reference._deposit_point

    def spy(self: _Reference, batch: int, pos: Any, de: float, species: int = -1, gen: int = -1):
        if species == 64:
            local.append((species, gen))
        elif de > 0.0:
            seen.append(
                (self.species_id if species < 0 else species,
                 self.generation if gen < 0 else gen, de)
            )  # fmt: skip
        return orig(self, batch, pos, de, species, gen)

    monkeypatch.setattr(_Reference, "_deposit_point", spy)
    res = Simulation(_config(tid, energy=60.0, n=40, e_cut=2.0, track_end=False)).run()
    assert res.valid and res.energy_balance.relative_residual <= 1e-12
    d = [x for x in seen if x[0] == 1]
    p_sec = [x for x in seen if x[0] == 0 and x[1] >= 1]
    assert d and p_sec and all(g == 1 for _, g, _ in d)
    assert all(0.0 < e <= 4.0 for _, _, e in d)  # deuteron cutoff deposits (end of transport)
    bad = [e for _, _, e in p_sec if not 0.0 < e <= 2.0]
    assert not bad, f"secondary-proton point deposits above E_cut = 2 MeV: {bad[:8]}"
    assert local and set(local) == {(64, 1)}


def _fake_event(particles: tuple[tuple[float, ...], ...], local: float = 0.0) -> Any:
    def sample(model: Any, rows: Any, t_lab: float, uni: Any, nu: Any = None) -> ScalarEvent:
        return ScalarEvent(
            True, 1, (0,) * 5, 0, 0, 0.0, t_lab - sum(p[4] - 0.0 for p in particles) * 0.0, 0.0,
            local, particles,
        )  # fmt: skip

    return sample


def _reference(tid: str) -> _Reference:
    return _Reference(Simulation(_config(tid, energy=100.0, n=2, track_end=False)).effective)


def test_event_frame_rotation_is_bitwise_for_a_parent_along_z(
    tid: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A parent along z reproduces the sampler's (normalised) lab directions bitwise; for a tilted
    parent the angle to the parent direction is preserved and the direction stays a unit vector."""
    ref = _reference(tid)
    t_lab = 30.0  # kinetic energy of the product
    calls = {}

    def sample(model: Any, rows: Any, t1: float, uni: Any, nu: Any = None) -> ScalarEvent:
        # ledger: T1 = 30 + 0 + 0 + imbalance
        prods = ((1.0, 0.0, 0.0, 0.0, float(model.species_mass_mev[1]) + 30.0, 11.0, -7.0, 53.0),)
        return ScalarEvent(True, 1, (0, 1, 0, 0, 0), 0, 0, 0.0, t1 - t_lab, 0.0, 0.0, prods)

    monkeypatch.setattr(events, "sample_event_scalar", sample)
    stack: list[tuple[float, ...]] = []
    rows = ref.nuc.rows[0]
    ref._event(0, 0, 0, 1, rows, 0.5, 100.0, (0.0, 0.0, 10.0), (0.0, 0.0, 1.0), (15, 15, 5), stack)
    pm = math.sqrt(11.0**2 + 7.0**2 + 53.0**2)
    assert stack[0][3:6] == (11.0 / pm, -7.0 / pm, 53.0 / pm)
    assert stack[0][6] == pytest.approx(t_lab) and stack[0][8] == 1.0 and stack[0][9] == 1.0
    assert stack[0][10:13] == (15.0, 15.0, 5.0)  # the parent's voxel indices, not recomputed
    assert ref.counters["nuclear_conservation"] == 0
    d = np.array([0.6, -0.48, 0.64])
    d /= np.linalg.norm(d)
    stack2: list[tuple[float, ...]] = []
    ref._event(0, 0, 0, 1, rows, 0.5, 100.0, (0.0, 0.0, 10.0), tuple(d), (1, 2, 3), stack2)
    out = np.array(stack2[0][3:6])
    assert np.linalg.norm(out) == pytest.approx(1.0, abs=1e-14)
    assert float(out @ d) == pytest.approx(53.0 / pm, abs=1e-14)
    del calls


# ---- forced fail-closed counters ----------------------------------------------------------------
def _forced(tid: str, monkeypatch: pytest.MonkeyPatch, n: int = 12, energy: float = 100.0) -> Any:
    _scaled_rows(monkeypatch, 40.0)
    cfg = _config(tid, energy=energy, n=n, allow_invalid=True, track_end=False)
    return Simulation(cfg).run()


@pytest.mark.parametrize("counter", ["queue_overflow", "genealogy_overflow"])
def test_stack_and_genealogy_overflow_fail_closed(
    tid: str, monkeypatch: pytest.MonkeyPatch, counter: str
) -> None:
    if counter == "queue_overflow":
        monkeypatch.setattr(reference, "STACK_CAPACITY", 0)
    else:

        def boom(*a: Any) -> int:
            raise CounterOverflowError("forced")

        monkeypatch.setattr(reference, "child_genealogy_id", boom)
    res = _forced(tid, monkeypatch)
    assert not res.valid and res.counters.as_dict()[counter] > 0
    b = res.energy_balance
    assert b.unaccounted_mev > 0 and b.relative_residual <= 1e-12  # routed, never dropped
    with pytest.raises(TransportLimitError):
        Simulation(replace(res.requested_config, run=replace(res.requested_config.run,
                   allow_invalid_result=False))).run()  # fmt: skip


def test_rejection_limit_fails_closed(tid: str, monkeypatch: pytest.MonkeyPatch) -> None:
    def sample(model: Any, rows: Any, t1: float, uni: Any, nu: Any = None) -> ScalarEvent:
        return ScalarEvent(False, 64, (0,) * 5, -1, -1, math.nan, math.nan, math.nan, math.nan, ())

    monkeypatch.setattr(events, "sample_event_scalar", sample)
    res = _forced(tid, monkeypatch)
    c = res.counters.as_dict()
    assert not res.valid and c["nuclear_rejection_limit"] > 0
    assert res.energy_balance.unaccounted_mev > 0 and res.energy_balance.relative_residual <= 1e-12


def test_majorant_violation_fails_closed(tid: str, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(_Reference, "_sigma", lambda self, rows, e: 1.0e9)
    res = _forced(tid, monkeypatch)
    assert not res.valid and res.counters.as_dict()["majorant_violation"] > 0
    assert res.energy_balance.unaccounted_mev > 0 and res.energy_balance.relative_residual <= 1e-12


def test_majorant_checked_on_every_step_not_only_candidates(
    tid: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A step that is not a candidate (the stub never reaches one) with Sigma(E1) above the
    majorant (an arbitrarily low straggling draw) is counted and ends the history invalid."""
    calls: list[int] = []

    def never(self: Any, h: int, batch: int, gid: int, nc: int, *args: Any) -> Any:
        calls.append(h)
        return nc, 1.0, False

    monkeypatch.setattr(_Reference, "_sigma", lambda self, rows, e: 1.0e9)
    monkeypatch.setattr(_Reference, "_candidate", never)
    res = _forced(tid, monkeypatch)
    assert not res.valid and res.counters.as_dict()["majorant_violation"] > 0
    assert calls == []  # the violation was found before any candidate


def test_nuclear_diagnostics_block_counts_events(ci_run: tuple[Any, Simulation]) -> None:
    res, _ = ci_run
    nuc = res.diagnostics["nuclear"]
    n_ev = sum(r["events"] for r in nuc.values())
    assert n_ev > 0
    for r in nuc.values():
        assert sum(r["residual"].values()) == r["events"]


def test_nuclear_conservation_counter_fires(tid: str, monkeypatch: pytest.MonkeyPatch) -> None:
    real = events.sample_event_scalar

    def broken(model: Any, rows: Any, t1: float, uni: Any, nu: Any = None) -> ScalarEvent:
        ev = real(model, rows, t1, uni, nu)
        return replace(ev, imbalance_mev=ev.imbalance_mev + 0.5) if ev.accepted else ev

    monkeypatch.setattr(events, "sample_event_scalar", broken)
    res = _forced(tid, monkeypatch)
    assert not res.valid and res.counters.as_dict()["nuclear_conservation"] > 0
    assert res.energy_balance.relative_residual <= 1e-12  # the miss is routed to unaccounted


# ---- P4: thinning harness through the reference glue -------------------------------------------
@pytest.mark.parametrize("shape", ["constant", "inv_sqrt"])
@pytest.mark.parametrize("s_max", [0.1, 1.0])
def test_p4_thinning_reproduces_the_survival(
    tid: str, monkeypatch: pytest.MonkeyPatch, shape: str, s_max: float
) -> None:
    """Synthetic Sigma (constant, or proportional to E^-1/2, whose majorant window is exercised
    on a rising Sigma) in water, MCS and straggling off; the events are empty (no products):
    the fraction of primaries ending in an event before depth z equals 1 - exp(-int rho Sigma dz)
    along the deterministic path within 4 sigma at 10 depths."""
    sigma0 = 1.0  # cm2/g at 25 MeV
    orig = NuclearTable.material_rows

    def synthetic(self: NuclearTable, material: Any, f_e: float = 0.02) -> MaterialNuclear:
        m = orig(self, material, f_e)
        g = m.grid_e_mev
        sig = np.full_like(g, sigma0) if shape == "constant" else sigma0 * np.sqrt(25.0 / g)
        win, end = majorant(g, sig, f_e)
        return replace(
            m, sigma_mass_cm2_g=sig, sigma_hat_window=win, sigma_hat_end=end,
            cum_sigma_mass_cm2_g=m.cum_fraction * sig,
        )  # fmt: skip

    def empty(model: Any, rows: Any, t1: float, uni: Any, nu: Any = None) -> ScalarEvent:
        return ScalarEvent(True, 1, (0,) * 5, 0, 0, 0.0, t1, 0.0, 0.0, ())

    monkeypatch.setattr(NuclearTable, "material_rows", synthetic)
    monkeypatch.setattr(events, "sample_event_scalar", empty)
    n = 1500
    sim = Simulation(
        _config(
            tid, energy=25.0, n=n, mcs=False, straggling=False, max_step=s_max, allow_invalid=True
        )  # fmt: skip
    )
    res = sim.run()
    assert res.counters.as_dict()["majorant_violation"] == 0
    code, pos = res.diagnostics["end_code"], res.diagnostics["end_position_mm"]
    eff = sim.effective
    assert eff.nuclear is not None
    rows = eff.nuclear.rows[0]
    r0 = eff.tables.range_g_cm2(0, 25.0)
    depths = np.linspace(0.5, 0.9, 10) * r0 * 10.0  # mm
    zz = np.linspace(0.0, depths[-1], 4001)
    e_z = np.array([eff.tables.energy_from_range(0, r0 - z / 10.0) for z in zz])
    sig_z = np.array([rows.sigma_at(e) for e in e_z])  # cm2/g, rho = 1 g/cm3
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (sig_z[1:] + sig_z[:-1]) * np.diff(zz) / 10.0)))
    worst = 0.0
    for dz in depths:
        p = 1.0 - math.exp(-float(np.interp(dz, zz, cum)))
        frac = float(np.mean((code == END_NUCLEAR) & (pos[:, 2] <= dz)))
        z = (frac - p) / math.sqrt(p * (1.0 - p) / n)
        worst = max(worst, abs(z))
    print(f"P4 {shape} s_max={s_max}: max |z| = {worst:.2f}")
    assert worst <= 4.0


def test_binding_tally_recomputed_from_event_counts_and_ame(
    ci_run: tuple[Any, Simulation], tid: str
) -> None:
    """The ``nuclear_binding`` tally equals the AME2020 recomputation from the per-target event,
    light-product and residual counts (1e-9 relative + 1e-9 MeV); a shifted count fails it."""
    import importlib.util
    from pathlib import Path

    from ionmc.data.ame import load_ame2020

    path = Path(__file__).resolve().parents[2] / "validation/scripts/transport/nuclear_checks.py"
    spec = importlib.util.spec_from_file_location("nuclear_checks_loop", path)
    assert spec and spec.loader
    chk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(chk)
    res, _ = ci_run
    ame = load_ame2020(cache.verify("ame2020-mass", cache.resolve_cache_dir(None)).read_text())
    targets = NuclearTable.load(None, tid).info["targets"]
    tally = float(res.energy_balance.nuclear_mev["nuclear_binding"])
    nuc = res.diagnostics["nuclear"]
    got = chk.recompute_binding_total(nuc, targets, ame)
    assert abs(got - tally) <= 1e-9 * abs(tally) + 1e-9
    k = next(iter(nuc))
    nuc[k]["light"]["a"] += 1
    assert abs(chk.recompute_binding_total(nuc, targets, ame) - tally) > 1.0
    nuc[k]["light"]["a"] -= 1


class _FixedNormal:
    """Proxy of the transport function namespace whose ``gauss_one`` returns ``z``."""

    def __init__(self, inner: Any, z: float) -> None:
        self._inner, self._z = inner, z

    def gauss_one(self, u0: Any, u1: Any) -> float:
        return self._z

    def __getattr__(self, name: str) -> Any:
        return getattr(self._inner, name)


def _source_draw_run(tid: str, monkeypatch: pytest.MonkeyPatch, z: float) -> Any:
    """Two histories of a 249 MeV source with sigma 0.1 MeV (E0 + 6 sigma passes validate) whose
    energy normal draw is forced to ``z``: 249 + 0.1 z MeV."""
    cfg = _config(tid, energy=249.0, n=2, allow_invalid=True)
    cfg = replace(cfg, source=replace(cfg.source, energy_sigma_mev=0.1))
    orig = _Reference.__init__

    def init(self: Any, *args: Any, **kw: Any) -> None:
        orig(self, *args, **kw)
        self.F = _FixedNormal(self.F, z)

    monkeypatch.setattr(_Reference, "__init__", init)
    return Simulation(cfg).run()


def test_nuclear_source_energy_above_the_limit_invalidates(
    tid: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Decision 0041 section 5 (amended): with nuclear=True a sampled source energy above 250 MeV
    (here 250.1 MeV from a forced +11 sigma draw, which the 6 sigma validation cannot exclude for
    an unbounded Gaussian) increments source_energy_out_of_range, books the energy as unaccounted
    and invalidates the result; +4 sigma (249.4 MeV) stays valid (Codex finding 1)."""
    bad = _source_draw_run(tid, monkeypatch, 11.0)
    assert bad.counters.source_energy_out_of_range == 2 and not bad.valid
    assert bad.energy_balance.unaccounted_mev == pytest.approx(2 * 250.1)
    assert bad.energy_balance.relative_residual <= 1e-12
    with monkeypatch.context() as m:
        ok = _source_draw_run(tid, m, 4.0)
    assert ok.counters.source_energy_out_of_range == 0 and ok.valid
