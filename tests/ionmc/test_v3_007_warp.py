"""Rows INF-1, INF-2, INF-3, INF-5, INF-7, MEM-1 (warp-cpu) and REUSE-1 (CI part) of plan v3-007
(step S3): the Warp driver transports beamlet sets (device context built once, per-beamlet
parameter and offset arrays, per-history beamlet id in the kernel, grouped int64 accumulators
with host reduction, instrumented allocation). warp-cpu only here; warp-cuda runs the same code
with the device parameter (HR tier). REG-1 (A16 digests, T1 trace fixture) is the unchanged
existing tests. Fixed small seeds, not evidence. Not a clinical claim."""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))

from test_v3_007_reference import (  # noqa: E402
    BETHE,
    BOX,
    GRID,
    POS,
    beamlets_of,
    closure_errors,
    quanta,
    rinf,
    sum_columns,
)

import ionmc.transport.warp_driver as wd  # noqa: E402
from ionmc.config import (  # noqa: E402
    DiagnosticsOptions,
    InfluenceOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    allocation_estimate,
    validate,
)
from ionmc.physics.projectiles import ALPHA, PROTON  # noqa: E402
from ionmc.simulation import Simulation  # noqa: E402
from ionmc.sources import Beamlet, BeamletSet, PencilBeamSource  # noqa: E402
from ionmc.transport.alloc import record_allocations  # noqa: E402
from ionmc.transport.parity import compare_traces  # noqa: E402
from ionmc.transport.run import run_transport  # noqa: E402
from ionmc.transport.tables import TransportTables  # noqa: E402
from ionmc.transport.tally import (  # noqa: E402
    QUANTUM_MEV,
    BeamletAccumulation,
    BeamletRaw,
    RawTransport,
)


def on_warp(
    c: SimulationConfig, precision: str = "float64", chunk: int | None = None
) -> SimulationConfig:
    """``c`` on the warp-cpu backend (optionally another precision / chunk size)."""
    kw: dict[str, Any] = {"backend": "warp-cpu", "precision": precision}
    if chunk is not None:
        kw["chunk_histories"] = chunk
    return dataclasses.replace(c, run=dataclasses.replace(c.run, **kw))


def lowe(n: int, b: int = 4, g: int | None = None, influence: bool = True, **kw: Any) -> Any:
    """R-INF geometry, low energies (cheap): 4 beamlets of ``n`` histories."""
    inf = InfluenceOptions("dose", group_size=g) if influence else None
    return rinf(n=n, b=b, influence=inf, energies=(20.0, 24.0, 28.0, 32.0), **kw)


def tr(c: SimulationConfig) -> RawTransport:
    return run_transport(validate(c))


@pytest.fixture(scope="module")
def inf_runs() -> dict[str, RawTransport]:
    """R-INF CI variant (n_j = 200, B = 4) on warp-cpu: with influence, composite, python."""
    return {
        "warp": tr(on_warp(rinf())),
        "composite": tr(on_warp(rinf(influence=None))),
        "python": tr(rinf()),
    }


# ---- INF-1 --------------------------------------------------------------------------------------


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_inf1_one_beamlet_equals_the_pencil_run_bitwise_on_warp_cpu(precision: str) -> None:
    b = Beamlet((0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 70.0, 0.5, 3.0, 40)
    ph = PhysicsOptions(nuclear=False, stopping_model=BETHE, multiple_scattering=True)

    def run(n: int | None) -> RunOptions:
        return RunOptions(
            backend="warp-cpu",
            precision=precision,
            seed=7,
            n_histories=n,
            n_batches=4,  # type: ignore[arg-type]
        )

    pencil = SimulationConfig(
        PencilBeamSource(PROTON, b.position_mm, b.direction, 70.0, 0.5, 3.0), BOX, (GRID,), ph,
        run(40),
    )  # fmt: skip
    one = SimulationConfig(
        BeamletSet(PROTON, (b,)), BOX, (GRID,), ph, run(None), influence=InfluenceOptions("dose")
    )
    rp, rb = tr(pencil), tr(one)
    assert np.array_equal(rp.edep_mev[0], rb.edep_mev[0])
    assert rp.tallies == rb.tallies and rp.counters == rb.counters
    assert rp.outside_mev == rb.outside_mev and rp.quantization_mev == rb.quantization_mev
    assert rp.quantization_mev == rb.quantization_mev
    bl = beamlets_of(rb)
    col = bl.accumulations[0].batch_quanta
    assert np.array_equal(col, quanta(rp))  # the column at tau = 0 is the composite grid
    assert np.array_equal(
        np.flatnonzero(col.sum(axis=0)), np.flatnonzero(rp.edep_mev[0].sum(axis=0))
    )
    assert bl.balances[0].tallies == rp.tallies and bl.group_size == 1


# ---- INF-2 / INF-7 ------------------------------------------------------------------------------


def test_inf2_composite_identity_on_warp_cpu(inf_runs: dict[str, RawTransport]) -> None:
    raw, comp_run = inf_runs["warp"], inf_runs["composite"]
    bl = beamlets_of(raw)
    assert [a.beamlet_id for a in bl.accumulations] == [0, 1, 2, 3]
    assert np.array_equal(sum_columns(bl), quanta(comp_run))  # (a) integer sum identity
    n_total = 800
    comp = raw.edep_mev[0].sum(axis=0) / n_total  # (b) per primary, relative 1e-12
    cols = sum(
        (200 / n_total) * (a.batch_quanta.sum(axis=0) * QUANTUM_MEV / 200) for a in bl.accumulations
    )
    nz = comp > 0
    assert np.all(np.abs(cols[nz] - comp[nz]) <= 1e-12 * comp[nz])
    assert np.array_equal(raw.edep_mev[0], comp_run.edep_mev[0])  # (c) composite bitwise
    assert raw.tallies == comp_run.tallies and raw.counters == comp_run.counters
    assert comp_run.beamlets is not None and comp_run.beamlets.accumulations == ()


def test_warp_cpu_equals_the_python_reference_on_the_r_inf_run(
    inf_runs: dict[str, RawTransport],
) -> None:
    w, p = beamlets_of(inf_runs["warp"]), beamlets_of(inf_runs["python"])
    assert all(
        np.array_equal(a.batch_quanta, b.batch_quanta)
        for a, b in zip(w.accumulations, p.accumulations, strict=True)
    )
    assert [x.n_histories for x in w.balances] == [x.n_histories for x in p.balances]
    for x, y in zip(w.balances, p.balances, strict=True):
        assert x.tallies == y.tallies and x.outside_mev == y.outside_mev


def test_inf7_per_beamlet_closure_and_negative_controls_on_warp_cpu(
    inf_runs: dict[str, RawTransport], monkeypatch: pytest.MonkeyPatch
) -> None:
    bl = beamlets_of(inf_runs["warp"])
    assert max(closure_errors(bl)) < 1e-12
    assert [x.n_histories for x in bl.balances] == [200] * 4
    for x in bl.balances:
        assert abs(x.closure_residual_mev) < 1e-9 * x.tallies["initial"]
    a = bl.accumulations  # control 1: columns permuted
    perm = tuple(
        BeamletAccumulation(a[i].beamlet_id, a[(i + 1) % 4].batch_quanta) for i in range(4)
    )
    assert min(closure_errors(dataclasses.replace(bl, accumulations=perm))) > 1e-6
    # control 2: the influence row of the kernel off by one (group base j0 + 1 in the launch only):
    # the columns hold other beamlets' deposits, so the closure of the group fails
    orig = wd.wp.launch

    def shifted(*a: Any, **k: Any) -> Any:
        k["inputs"][0].j0 = k["inputs"][0].j0 + 1
        return orig(*a, **k)

    with monkeypatch.context() as mp:
        mp.setattr(wd.wp, "launch", shifted)
        bad = beamlets_of(tr(on_warp(lowe(8))))
    assert max(closure_errors(bad)) > 1e-6
    assert max(closure_errors(beamlets_of(tr(on_warp(lowe(8)))))) < 1e-12  # the unshifted run


# ---- INF-3 --------------------------------------------------------------------------------------


def test_inf3_chunk_and_group_invariance_on_warp_cpu() -> None:
    """J = 4, n_j = 600, B = 4: chunk_histories in {2^10, 2^11} x G in {1, 2, J} (six variants;
    chunks span beamlets inside a group) give bitwise equal per-beamlet accumulators and
    balances and composite grids."""
    ref: RawTransport | None = None
    for chunk in (1024, 2048):
        for g in (1, 2, 4):
            raw = tr(on_warp(lowe(600, g=g), chunk=chunk))
            bl = beamlets_of(raw)
            assert bl.group_size == g
            if ref is None:
                ref = raw
                continue
            rb = beamlets_of(ref)
            for x, y in zip(bl.accumulations, rb.accumulations, strict=True):
                assert x.beamlet_id == y.beamlet_id and np.array_equal(
                    x.batch_quanta, y.batch_quanta
                )
            assert bl.balances == rb.balances
            assert np.array_equal(raw.edep_mev[0], ref.edep_mev[0])
            assert raw.tallies == ref.tallies and raw.counters == ref.counters


# ---- INF-5 --------------------------------------------------------------------------------------


def test_inf5_trajectory_parity_python_vs_warp_cpu_float64_per_beamlet() -> None:
    """The first 3 histories of beamlets 0, 1 and J - 1 (J = 4, n_j = 4): python and warp-cpu
    float64 traces agree with the T1 tolerances (discrete columns exactly, continuous columns
    rtol = atol = 1e-10)."""
    bs = tuple(
        Beamlet((x, y, -5.0), (0.0, 0.0, 1.0), e, 0.3, 1.0, 4)
        for (x, y), e in zip(POS, (20.0, 22.0, 24.0, 26.0), strict=True)
    )
    ph = PhysicsOptions(nuclear=False, stopping_model=BETHE, multiple_scattering=True)

    def run(backend: str) -> Any:
        r = RunOptions(backend=backend, precision="float64", seed=5, n_histories=None, n_batches=2)  # type: ignore[arg-type]
        c = SimulationConfig(
            BeamletSet(PROTON, bs), BOX, (GRID,), ph, r,
            diagnostics=DiagnosticsOptions(trace_histories=16, track_end_positions=True),
        )  # fmt: skip
        return Simulation(c).run().diagnostics

    py, wp_ = run("python"), run("warp-cpu")
    want = {h for s in (0, 4, 12) for h in (s, s + 1, s + 2)}

    def pick(d: dict[str, Any]) -> dict[str, Any]:
        t = d["trace"]
        m = np.isin(np.asarray(t["history"]), sorted(want))
        me = np.isin(np.asarray(d["trace_end_history"]), sorted(want))
        out: dict[str, Any] = {"trace": {k: np.asarray(v)[m] for k, v in t.items()}}
        out.update({k: np.asarray(d[k])[me] for k in d if k.startswith("trace_end_")})
        return out

    a, b = pick(py), pick(wp_)
    assert set(np.unique(a["trace"]["history"])) == want
    v = compare_traces(a, b)
    assert v["pass"], v
    assert compare_traces(py, wp_)["pass"]  # all 16 traced histories as well


# ---- MEM-1 (warp-cpu) ---------------------------------------------------------------------------


def _estimate(eff: Any, g: int, backend: str = "warp-cpu") -> int:
    cfg = eff.requested
    return allocation_estimate(
        backend, cfg.run.n_batches, GRID.n_voxels, 1, eff, eff.n_histories,
        cfg.run.chunk_histories, GRID.n_voxels, g, [int(x) for x in cfg.source.offsets()],
    )  # fmt: skip


@pytest.mark.parametrize("precision", ["float64", "float32"])
@pytest.mark.parametrize("g", [1, 2, 4])
def test_mem1_warp_cpu_measured_peak_matches_the_estimate(precision: str, g: int) -> None:
    eff = validate(on_warp(lowe(8, g=g), precision))
    with record_allocations() as rec:
        run_transport(eff)
    est = _estimate(eff, g)
    assert rec.peak_bytes <= est <= 1.5 * rec.peak_bytes
    kinds = {(k, t) for k, t, _ in rec.events}
    assert {("alloc", "tally_rows"), ("release", "tally_rows"), ("release", "influence")} <= kinds
    assert rec.live_bytes == 8 * 4 * GRID.n_voxels  # only the grid accumulators stay live


def test_mem1_warp_cpu_ion_run() -> None:
    bs = tuple(
        Beamlet((x, y, -5.0), (0.0, 0.0, 1.0), e, 0.5, 1.0, 4)
        for (x, y), e in zip(POS[:2], (900.0, 950.0), strict=True)
    )
    ph = PhysicsOptions(nuclear=False, stopping_model=BETHE, multiple_scattering=True)
    r = RunOptions(backend="warp-cpu", precision="float64", seed=3, n_histories=None, n_batches=2)
    c = SimulationConfig(
        BeamletSet(ALPHA, bs), BOX, (GRID,), ph, r, influence=InfluenceOptions("dose", group_size=2)
    )
    eff = validate(c)
    with record_allocations() as rec:
        raw = run_transport(eff)
    cfg = eff.requested
    est = allocation_estimate(
        "warp-cpu", 2, GRID.n_voxels, 1, eff, 8, cfg.run.chunk_histories, GRID.n_voxels, 2,
        [int(x) for x in cfg.source.offsets()],
    )  # fmt: skip
    assert rec.peak_bytes <= est <= 1.5 * rec.peak_bytes
    assert max(closure_errors(beamlets_of(raw))) < 1e-12


def test_mem1_warp_cpu_a_lingering_array_fails_the_audit() -> None:
    from ionmc.transport.alloc import alloc_zeros

    eff = validate(on_warp(lowe(8, g=4)))
    with record_allocations() as rec:
        alloc_zeros("rogue", (100000,), np.int64)
        run_transport(eff)
    assert rec.peak_bytes > _estimate(eff, 4)


# ---- REUSE-1 (CI part) --------------------------------------------------------------------------


def test_reuse1_second_run_builds_and_uploads_nothing(monkeypatch: pytest.MonkeyPatch) -> None:
    sim = Simulation(on_warp(lowe(8, g=2)))
    first = sim.run()
    calls = {"to_warp": 0, "load_kernel": 0}
    o_to, o_load = TransportTables.to_warp, wd.load_kernel

    def to_warp(self: Any, *a: Any, **k: Any) -> Any:
        calls["to_warp"] += 1
        return o_to(self, *a, **k)

    def load_kernel(*a: Any, **k: Any) -> Any:
        calls["load_kernel"] += 1
        return o_load(*a, **k)

    monkeypatch.setattr(TransportTables, "to_warp", to_warp)
    monkeypatch.setattr(wd, "load_kernel", load_kernel)
    stats = dict(wd.CONTEXT_STATS)
    second = sim.run()
    assert dict(wd.CONTEXT_STATS) == stats and calls == {"to_warp": 0, "load_kernel": 0}
    fresh = Simulation(on_warp(lowe(8, g=2))).run()  # a fresh object builds its own context
    assert wd.CONTEXT_STATS["builds"] == stats["builds"] + 1 and calls["to_warp"] == 1
    for res in (second, fresh):
        assert np.array_equal(res.grids[0].batch_energy_mev, first.grids[0].batch_energy_mev)
        assert res.beamlet_raw is not None and first.beamlet_raw is not None
        a, b = res.beamlet_raw, first.beamlet_raw
        assert a.balances == b.balances
        for x, y in zip(a.accumulations, b.accumulations, strict=True):
            assert np.array_equal(x.batch_quanta, y.batch_quanta)
        assert res.counters == first.counters
    assert isinstance(first.beamlet_raw, BeamletRaw)
