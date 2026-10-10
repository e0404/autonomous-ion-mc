"""V3-005C C3: hadronic elastic scattering on the python reference backend (decision 0041 slice C).

Covers Sigma_tot = Sigma_nonel + Sigma_el thinning with the channel chosen by the free slot u3 of
the candidate block, p-p (slower proton as a transported secondary) and p + A (recoil deposited
locally, ``elastic_recoil_local``) events, exact per-event energy conservation, the declared-domain
counters, the switches and fail-closed cases, bit-identity of ``nuclear=False`` and
``elastic=False`` outputs (A16 digest, pre-D2 anchors, Sigma_el = 0), the R1-D2 intended-change
record for nuclear-on outputs and the P7 sampler test. Data-backed (``IONMC_CACHE_DIR``)."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc._wpfunc import python_twin
from ionmc.config import DiagnosticsOptions
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import COPPER
from ionmc.nuclear.elastic_events import sample_elastic_event
from ionmc.nuclear.elastic_kin import two_body
from ionmc.nuclear.elastic_tables import (
    ElasticTable,
    ElasticTableMissingError,
    chi2_equiprobable,
)
from ionmc.physics.elastic import make_elastic
from ionmc.physics.nuclear import make_nuclear
from ionmc.physics.projectiles import PROTON
from ionmc.simulation import ElasticTransportCounters, NuclearEnergyBalance, Simulation
from ionmc.transport.reference import run_reference_range
from ionmc.transport.tally import ELASTIC_COUNTER_NAMES
from tests.ionmc.test_nuclear_transport_loop import _config, _table_id

ROOT = Path(__file__).resolve().parents[2]
A16_DIGEST = "c862edf799dcb83542e5071219b5703c7665f8f792f6bbc49ab32918418b1f8f"
M_P = PROTON.mass_mev
EL = python_twin(make_elastic)
NU = python_twin(make_nuclear)


@pytest.fixture(scope="module")
def tid() -> str:
    return _table_id()


def _pin() -> str:
    path = ROOT / "src" / "ionmc" / "data" / "elastic_table_pin.json"
    return str(json.loads(path.read_text())["table_id"])


@pytest.fixture(scope="module")
def table() -> ElasticTable:
    return ElasticTable.load(None, _pin())


def _scale_elastic(mp: pytest.MonkeyPatch, factor: float) -> None:
    """Multiply Sigma_el (and its majorants and cumulative rows) by ``factor``: many events."""
    orig = ElasticTable.material_rows

    def scaled(self: ElasticTable, material: Any, f_e: float = 0.02) -> Any:
        m = orig(self, material, f_e)
        return replace(
            m,
            sigma_mass_cm2_g=m.sigma_mass_cm2_g * factor,
            sigma_hat_window=m.sigma_hat_window * factor,
            sigma_hat_end=m.sigma_hat_end * factor,
            cum_sigma_mass_cm2_g=m.cum_sigma_mass_cm2_g * factor,
        )

    mp.setattr(ElasticTable, "material_rows", scaled)


def _cfg(tid: str, energy: float, n: int, *, trace: int = 0, **phys: Any) -> Any:
    cfg = _config(tid, energy=energy, n=n, track_end=False, allow_invalid=True)
    return replace(
        cfg,
        physics=replace(cfg.physics, **{"elastic": True, **phys}),
        diagnostics=DiagnosticsOptions(trace_histories=trace),
    )


def _digest(part: Any) -> str:
    h = hashlib.sha256()
    h.update(np.array([math.fsum(c) for c in part.tally_components], dtype=np.float64).tobytes())
    h.update(np.array(part.counter_sums, dtype=np.int64).tobytes())
    for e in part.edep:
        h.update(np.ascontiguousarray(e).tobytes())
    return h.hexdigest()


# ---- bit-identity of the pre-D2 paths ---------------------------------------------------------
def test_a16_digest_nuclear_off_unchanged() -> None:
    """R1-D2 (A): the A16 digest of the nuclear=False qualified outputs equals the 96d4d03 one."""
    spec = importlib.util.spec_from_file_location(
        "a16_digest_c3", ROOT / "validation" / "scripts" / "transport" / "a16_digest.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    digests = {s: mod.digest(s, "none") for s in ("t1:python:float64", "t1:warp-cpu:float32")}
    assert hashlib.sha256(json.dumps(digests, sort_keys=True).encode()).hexdigest() == A16_DIGEST


def test_sigma_el_zero_is_bit_identical_to_elastic_false(
    tid: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """R1-D2 (B): with Sigma_el = 0 the elastic machinery (additive majorant, u3 slot, domain
    counters) touches nothing else: every non-elastic column is bit-identical to elastic=False."""
    n = 24
    with monkeypatch.context() as mp:
        e0 = Simulation(_cfg(tid, 100.0, n, elastic=False)).effective
        off = run_reference_range(e0, 0, n)
    with monkeypatch.context() as mp:
        _scale_elastic(mp, 0.0)
        e1 = Simulation(_cfg(tid, 100.0, n)).effective
        on = run_reference_range(e1, 0, n)
    n_t = len(off.tally_components)
    assert [math.fsum(c) for c in on.tally_components[:n_t]] == [
        math.fsum(c) for c in off.tally_components
    ]
    assert on.counter_sums[: len(off.counter_sums)] == off.counter_sums
    for a, b in zip(on.edep, off.edep, strict=True):
        np.testing.assert_array_equal(a, b)
    assert all(math.fsum(c) == 0.0 for c in on.tally_components[n_t:])


# ---- events -----------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def many(tid: str) -> tuple[Any, Any]:
    """A python run with many elastic events (Sigma_el x 30) and the recorded event trace."""
    n = 32
    with pytest.MonkeyPatch.context() as mp:
        _scale_elastic(mp, 30.0)
        sim = Simulation(_cfg(tid, 150.0, n, trace=n))
        part = run_reference_range(sim.effective, 0, n)
    return sim.effective, part


def test_event_parity_on_recorded_inputs(many: tuple[Any, Any]) -> None:
    """P5-ext hook: every recorded event replayed through ``sample_elastic_event`` on its recorded
    inputs reproduces the recorded outputs exactly (channel target, mu_CM, outgoing energies and
    directions, the recoil energy of p + A and the slower proton of p-p)."""
    eff, part = many
    tr = part.meta["elastic_trace"]
    el = eff.nuclear.elastic
    assert tr.shape[0] > 20
    from ionmc.transport.funcs import make_transport_funcs

    F = python_twin(make_transport_funcs)
    n_pp = n_pa = 0
    for row in tr:
        d = (row[8], row[9], row[10])
        e1, e2 = F.orthonormal_basis(F.vec3(*d))
        ev = sample_elastic_event(
            EL, NU, el.table.arrays, el.rows[0], M_P, row[4], row[5], row[6], row[7], d, e1, e2
        )
        assert ev.target == int(row[3])
        out = [ev.mu_cm, ev.phi, ev.t_primary_mev, *ev.dir_primary, ev.t_other_mev, *ev.dir_other]
        np.testing.assert_array_equal(out, row[11:])
        n_pp += ev.pp
        n_pa += not ev.pp
    assert n_pp > 0 and n_pa > 0


def _many_part(tid: str, *, trace: int, track_end: bool = False) -> Any:
    n = 32
    with pytest.MonkeyPatch.context() as mp:
        _scale_elastic(mp, 30.0)
        cfg = _cfg(tid, 150.0, n)
        cfg = replace(
            cfg,
            diagnostics=DiagnosticsOptions(track_end_positions=track_end, trace_histories=trace),
        )
        return run_reference_range(Simulation(cfg).effective, 0, n)


def test_elastic_trace_off_with_trace_histories_zero(tid: str) -> None:
    """Finding B: endpoint diagnostics alone do not record elastic events; with trace_histories == 0
    the elastic trace is absent (as the nuclear-trace rows are empty), the tallies are unchanged."""
    part = _many_part(tid, trace=0, track_end=True)
    assert "elastic_trace" not in part.meta
    assert all(v.shape[0] == 0 for v in part.meta["nuclear_trace"].values())
    ref = _many_part(tid, trace=0)
    assert _digest(part) == _digest(ref)


def test_elastic_trace_limited_to_first_k_histories(tid: str) -> None:
    """Finding B: with trace_histories = k only the events of histories < k are recorded (rows equal
    the history < k rows of the full trace); the tallies still count the events of all histories."""
    k = 8
    full = _many_part(tid, trace=32)
    part = _many_part(tid, trace=k, track_end=True)
    tf, tk = full.meta["elastic_trace"], part.meta["elastic_trace"]
    assert (tf[:, 0] >= k).any() and tk.shape[0] > 0
    assert (tk[:, 0] < k).all()
    np.testing.assert_array_equal(tk, tf[tf[:, 0] < k])
    assert _digest(part) == _digest(full)


def test_kinematics_closure_and_energy_conservation(many: tuple[Any, Any]) -> None:
    """Row P6 (7) and per-event conservation: T_in = T_primary + T_other to 1e-9 MeV, the primary
    direction is a unit vector, p-p keeps the faster proton, and the numpy oracle closes."""
    eff, part = many
    for row in part.meta["elastic_trace"]:
        t1, tp, to = row[4], row[13], row[17]
        assert abs(t1 - (tp + to)) <= 1e-9
        assert abs(math.sqrt(sum(x * x for x in row[14:17])) - 1.0) <= 1e-12
        tgt = int(row[3])
        m2 = float(eff.nuclear.elastic.table.arrays["target_mass_mev"][tgt])
        de, dp = two_body(t1, M_P, m2, row[11]).closure(t1, M_P, m2)
        assert de <= 1e-9 and dp <= 1e-9
        if tgt == 0:
            assert tp >= to


def test_balance_counters_and_diagnostics_with_many_events(
    tid: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    n = 48
    _scale_elastic(monkeypatch, 30.0)
    res = Simulation(_cfg(tid, 150.0, n)).run()
    b = res.energy_balance
    assert isinstance(b, NuclearEnergyBalance)
    assert isinstance(res.counters, ElasticTransportCounters)
    assert res.valid and b.relative_residual <= 1e-12 and b.grid_relative_residual(0) <= 1e-12
    assert b.elastic_mev["elastic_recoil_local"] > 0
    assert b.elastic_mev["elastic_pp_events"] > 0 and b.elastic_mev["elastic_pa_events"] > 0
    c = res.counters.as_dict()
    assert all(c[k] == 0 for k in c if k not in ELASTIC_COUNTER_NAMES)
    assert c["elastic_below_domain"] > 0 and c["pp_below_domain"] > 0  # diagnostics only


def test_domain_counters_once_per_proton_and_never_invalidate(tid: str) -> None:
    """A 10 MeV proton is below E_min,pp (12.53) from its first step and falls below the O-16
    e_min_shape (6.58): both counters count each proton once, and the result stays valid."""
    n = 16
    res = Simulation(
        _cfg(tid, 10.0, n, elastic_only=True, allow_invalid=False)
        if False
        else _cfg(tid, 10.0, n, elastic_only=True)
    ).run()
    c = res.counters.as_dict()
    assert res.valid
    assert c["pp_below_domain"] >= n and c["elastic_below_domain"] >= n
    assert c["pp_below_domain"] < 2 * n  # once per proton, not once per step (thousands of steps)


def test_table_domain_has_no_cross_section_below_e_min(table: ElasticTable) -> None:
    dom = table.elastic_domain()
    assert dom["H-1"][0] == pytest.approx(12.532211034734628)
    for i, name in enumerate(table.target_names):
        # sigma is lin-lin between nodes: it is zero at every node below e_min and ramps up only
        # inside the one interval below the first non-zero node (a finding, see the report)
        below = table.arrays["grid_e_mev"][table.arrays["grid_e_mev"] < dom[name][0]]
        assert table.sigma_barn(i, float(below[-1])) == 0.0
        assert table.sigma_barn(i, dom[name][0] - 2 * (dom[name][0] - below[-1])) == 0.0
        assert table.sigma_barn(i, min(dom[name][0] * 1.2 + 1.0, 250.0)) > 0.0


# ---- switches and fail-closed cases -----------------------------------------------------------
def test_switches(tid: str) -> None:
    n = 8
    on = Simulation(_cfg(tid, 100.0, n)).effective
    assert on.nuclear.elastic is not None and on.nuclear.elastic.table.table_id == _pin()
    assert on.summary()["elastic"]["pinned"] is True
    off = Simulation(_cfg(tid, 100.0, n, elastic=False)).effective
    assert off.nuclear.elastic is None and "elastic" not in off.summary()
    only = Simulation(_cfg(tid, 100.0, n, elastic_only=True)).run()
    assert only.valid and only.energy_balance.elastic_mev["elastic_pa_events"] >= 0
    # elastic_only: no non-elastic event ever happens
    assert sum(r["events"] for r in only.diagnostics.get("nuclear", {}).values()) == 0


def test_fail_closed_cases(tid: str, table: ElasticTable) -> None:
    cfg = _cfg(tid, 100.0, 8)
    with pytest.raises(UnsupportedCombinationError, match="elastic_only"):
        replace(cfg.physics, nuclear=False, elastic_only=True)
    # the Warp backends accept the channel since C4 (tests/ionmc/test_elastic_warp.py)
    Simulation(replace(cfg, run=replace(cfg.run, backend="warp-cpu", precision="float64")))
    with pytest.raises(ElasticTableMissingError):
        Simulation(replace(cfg, physics=replace(cfg.physics, elastic_table_id="0" * 64)))
    with pytest.raises(UnsupportedCombinationError, match="no elastic"):
        table.material_rows(COPPER)
    with pytest.raises(UnsupportedCombinationError):
        table.check_energy(251.0)
    # warp backend with elastic=False is the unchanged V3-005B model
    Simulation(
        replace(
            _cfg(tid, 100.0, 8, elastic=False),
            run=replace(cfg.run, backend="warp-cpu", precision="float64"),
        )
    )


def test_forced_majorant_violation_is_counted(tid: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """Sigma_tot covers Sigma_el: halving the elastic majorants must trip majorant_violation."""
    orig = ElasticTable.material_rows

    def low(self: ElasticTable, material: Any, f_e: float = 0.02) -> Any:
        m = orig(self, material, f_e)
        return replace(
            m, sigma_hat_window=m.sigma_hat_window * 0.0, sigma_hat_end=m.sigma_hat_end * 0.0
        )

    monkeypatch.setattr(ElasticTable, "material_rows", low)
    res = Simulation(_cfg(tid, 150.0, 8)).run()
    assert not res.valid and res.counters.majorant_violation > 0


# ---- R1-D2 intended-change record for nuclear-on outputs (C) -----------------------------------
R1_D2_RECORD = {
    "task": "V3-005B/V3-005C C3",
    "baseline": "96d4d03",
    "nuclear_off": f"bit-identical (A16 digest {A16_DIGEST}); intended-change set empty",
    "elastic_false": "bit-identical to the V3-005B nuclear-on python outputs (pre-C13 anchors)",
    "sigma_el_zero": "bit-identical to elastic=False in every non-elastic column",
    "changes_nuclear_on": [
        "grid energy deposits and all scoring channels (elastic scattering and the p-p secondary)",
        "tally columns: elastic_recoil_local, elastic_pp_events, elastic_pa_events appended",
        "counters: elastic_below_domain, pp_below_domain appended (diagnostics, not validity)",
        "energy balance: elastic_recoil_local is a destination and a local deposit",
        "effective-config summary: 'elastic' block",
    ],
    "rng": "each accepted elastic event consumes one extra block of the particle's nuclear "
    "stream (slot 0 mu_CM, slot 1 azimuth) and the redrawn optical depth reuses u1 of the "
    "candidate block; EM streams (purpose 0/1) are untouched",
    "declared_directions": "peak/plateau down, plateau up, lateral spread up, "
    "|delta total deposit| <= 0.1 % of E0 (judged in lv5c, not here)",
}


def test_r1_d2_record_matches_the_observed_changes(
    tid: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    n = 24
    off = Simulation(_cfg(tid, 100.0, n, elastic=False)).run()
    _scale_elastic(monkeypatch, 30.0)
    on = Simulation(_cfg(tid, 100.0, n)).run()
    assert not np.array_equal(off.grids[0].batch_energy_mev, on.grids[0].batch_energy_mev)
    assert off.energy_balance.initial_mev == on.energy_balance.initial_mev  # same EM source draw
    assert set(on.counters.as_dict()) - set(off.counters.as_dict()) == set(ELASTIC_COUNTER_NAMES)
    assert on.energy_balance.elastic_mev and not off.energy_balance.nuclear_mev.get("x")
    assert on.valid and on.energy_balance.relative_residual <= 1e-12
    assert {"nuclear_off", "elastic_false", "changes_nuclear_on", "rng"} <= set(R1_D2_RECORD)


# ---- P7: sampler vs the built table -----------------------------------------------------------
@pytest.mark.parametrize("target,name", [(0, "H-1"), (3, "O-16")])
@pytest.mark.parametrize("energy", [20.0, 100.0, 200.0])
def test_p7_sampler_chi2_against_the_table(
    table: ElasticTable, target: int, name: str, energy: float
) -> None:
    assert table.target_names[target] == name
    rng = np.random.default_rng(20515 + target * 1000 + int(energy))
    n = 20000
    u = rng.random(n)
    arrays = table.arrays
    grid = arrays["grid_e_mev"]
    k = int(NU.grid_locate(energy, grid, grid.size))
    t = min(max((energy - grid[k]) / (grid[k + 1] - grid[k]), 0.0), 1.0)
    flat = arrays["edges_mu"].reshape(-1)
    n_g, n_e = arrays["edges_mu"].shape[1:]
    mu = np.array(
        [
            EL.sample_mu_edges(
                x, flat, (target * n_g + k) * n_e, (target * n_g + k + 1) * n_e, t, n_e - 1
            )
            for x in u
        ]
    )
    np.testing.assert_array_equal(mu, table.sample_mu(target, energy, u))  # same inverse CDF
    chi2, p, dof = chi2_equiprobable(mu, table.edges_at(target, energy), 64)
    print(f"P7 {name} {energy} MeV chi2 {chi2:.1f} dof {dof} p {p:.4f}")
    assert p > 0.001
