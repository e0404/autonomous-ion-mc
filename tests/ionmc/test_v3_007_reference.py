"""Rows BEAM-4, BEAM-5, BEAM-6, INF-1, INF-2, INF-6, INF-7 and MEM-1 (python) of plan v3-007 (step
S2): the python reference backend transports beamlet sets (per-history beamlet id, per-beamlet
integer accumulators and exact per-beamlet balances). Fixed small seeds, not evidence.
Not a clinical claim."""

from __future__ import annotations

import dataclasses
import os
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))

import ionmc.transport.reference as reference  # noqa: E402
from ionmc.config import (  # noqa: E402
    InfluenceOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    allocation_estimate,
    validate,
)
from ionmc.geometry import BoxPhantom  # noqa: E402
from ionmc.materials import WATER  # noqa: E402
from ionmc.physics.projectiles import PROTON  # noqa: E402
from ionmc.rng.philox import PURPOSE_SOURCE, key_from_seed  # noqa: E402
from ionmc.scoring import ScoringGrid  # noqa: E402
from ionmc.simulation import Simulation  # noqa: E402
from ionmc.sources import Beamlet, BeamletSet, PencilBeamSource  # noqa: E402
from ionmc.transport.alloc import alloc_zeros, record_allocations  # noqa: E402
from ionmc.transport.run import run_transport  # noqa: E402
from ionmc.transport.tally import (  # noqa: E402
    QUANTUM_MEV,
    BeamletAccumulation,
    BeamletRaw,
    RawTransport,
)

BETHE = "bethe-ionmc-v1"
GRID = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 100))
BOX = BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, 200.0), WATER)
_DEFAULT_INFLUENCE = InfluenceOptions("dose")
POS = ((-10.0, 0.0), (10.0, 0.0), (0.0, -10.0), (0.0, 10.0))


def rinf(
    n: int = 200,
    b: int = 4,
    influence: InfluenceOptions | None = _DEFAULT_INFLUENCE,
    energies: tuple[float, ...] = (70.0, 90.0, 110.0, 130.0),
    seed: int = 11,
    nuclear: bool = False,
    tid: str | None = None,
    sigma_e: float = 0.5,
) -> SimulationConfig:
    """R-INF CI variant: 4 beamlets (lateral +-10 mm), sigma_E 0.5 MeV, sigma_lat 3 mm, +z."""
    bs = tuple(
        Beamlet((x, y, -5.0), (0.0, 0.0, 1.0), e, sigma_e, 3.0, n)
        for (x, y), e in zip(POS, energies, strict=False)
    )
    ph = PhysicsOptions(nuclear=nuclear, stopping_model=BETHE, multiple_scattering=True)
    if nuclear:
        ph = dataclasses.replace(ph, nuclear_table_id=tid)
    run = RunOptions(
        backend="python", precision="float64", seed=seed, n_histories=None, n_batches=b
    )
    return SimulationConfig(BeamletSet(PROTON, bs), BOX, (GRID,), ph, run, influence=influence)


def transport(c: SimulationConfig) -> RawTransport:
    return run_transport(validate(c))


def quanta(raw: RawTransport) -> np.ndarray:
    """The composite integer accumulators [B, V] recovered exactly from the grid energies."""
    return np.rint(raw.edep_mev[0] / QUANTUM_MEV).astype(np.int64)


def beamlets_of(raw: RawTransport) -> BeamletRaw:
    assert raw.beamlets is not None
    return raw.beamlets


def sum_columns(bl: BeamletRaw) -> np.ndarray:
    return sum(a.batch_quanta for a in bl.accumulations)  # type: ignore[no-any-return]


def closure_errors(bl: BeamletRaw, grid: int = 0) -> list[float]:
    """INF-7: relative mismatch of (sum_v column quanta Q + outside + quantization) - deposited."""
    bal = {x.beamlet_id: x for x in bl.balances}
    out = []
    for a in bl.accumulations:
        x = bal[a.beamlet_id]
        tot = sum(int(v) for v in a.batch_quanta.sum(axis=1)) * QUANTUM_MEV
        lhs = tot + x.outside_mev[grid] + x.quantization_mev[grid]
        out.append(abs(lhs - x.deposited_mev) / abs(x.deposited_mev))
    return out


@pytest.fixture(scope="module")
def inst() -> dict[str, Any]:
    """The R-INF CI run with influence, instrumented: source-draw counters, sampled source values
    and the payload allocation events (BEAM-4, BEAM-5, MEM-1 with G = J)."""
    draws: list[tuple[int, int]] = []
    samples: list[tuple[int, tuple[Any, ...], tuple[Any, ...]]] = []
    orig_draw, orig_sample = reference.draw_block, reference.sample_source

    def draw(key: Any, c0: int, c1: int, c2: int, c3: int) -> Any:
        if c3 == PURPOSE_SOURCE:
            draws.append((c0, c3))
        return orig_draw(key, c0, c1, c2, c3)

    def sample(F: Any, key: Any, h: int, *args: Any) -> Any:
        out = orig_sample(F, key, h, *args)
        samples.append((h, args, out))
        return out

    c = rinf()
    eff = validate(c)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(reference, "draw_block", draw)
        mp.setattr(reference, "sample_source", sample)
        with record_allocations() as rec:
            raw = run_transport(eff)
    return {"cfg": c, "eff": eff, "raw": raw, "draws": draws, "samples": samples, "rec": rec}


@pytest.fixture(scope="module")
def composite() -> RawTransport:
    return transport(rinf(influence=None))


# ---- BEAM-4 / BEAM-5 ----------------------------------------------------------------------------


def test_beam4_source_counters_are_unique_and_equal_range_of_total(inst: dict[str, Any]) -> None:
    n_total = 800
    assert len(inst["draws"]) == n_total
    assert sorted(inst["draws"]) == [(h, PURPOSE_SOURCE) for h in range(n_total)]


def test_beam5_first_history_of_every_beamlet_samples_its_own_parameters(
    inst: dict[str, Any],
) -> None:
    from ionmc._wpfunc import python_twin
    from ionmc.transport.funcs import make_transport_funcs

    F = python_twin(make_transport_funcs)
    key = key_from_seed(inst["cfg"].run.seed)
    bs = inst["cfg"].source
    offs = [int(x) for x in bs.offsets()]
    by_h = {h: (args, out) for h, args, out in inst["samples"]}
    assert len(by_h) == 800
    energies = []
    for j, b in enumerate(bs.beamlets):
        args, out = by_h[offs[j]]
        want_args = (
            b.position_mm,
            b.unit_direction(),
            b.kinetic_energy_mev,
            b.energy_sigma_mev,
            b.lateral_sigma_mm,
        )
        assert args == want_args
        assert out == reference.sample_source(F, key, offs[j], *want_args)  # bitwise
        energies.append(out[0])
    assert len(set(energies)) == 4


# ---- INF-1 --------------------------------------------------------------------------------------


def test_inf1_one_beamlet_equals_the_pencil_run_bitwise() -> None:
    b = Beamlet((0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 70.0, 0.5, 3.0, 40)
    ph = PhysicsOptions(nuclear=False, stopping_model=BETHE, multiple_scattering=True)

    def run(n: int | None) -> RunOptions:
        return RunOptions(backend="python", precision="float64", seed=7, n_histories=n, n_batches=4)

    pencil = SimulationConfig(
        PencilBeamSource(PROTON, b.position_mm, b.direction, 70.0, 0.5, 3.0),
        BOX,
        (GRID,),
        ph,
        run(40),
    )
    one = SimulationConfig(
        BeamletSet(PROTON, (b,)), BOX, (GRID,), ph, run(None), influence=InfluenceOptions("dose")
    )
    rp, rb = transport(pencil), transport(one)
    assert np.array_equal(rp.edep_mev[0], rb.edep_mev[0])
    assert rp.tallies == rb.tallies and rp.counters == rb.counters
    assert rp.outside_mev == rb.outside_mev and rp.quantization_mev == rb.quantization_mev
    bl = beamlets_of(rb)
    col = bl.accumulations[0].batch_quanta
    assert np.array_equal(col, quanta(rp))
    assert np.array_equal(
        np.flatnonzero(col.sum(axis=0)), np.flatnonzero(rp.edep_mev[0].sum(axis=0))
    )
    assert bl.group_size == 1


# ---- INF-2 / INF-7 ------------------------------------------------------------------------------


def test_inf2_composite_identity(inst: dict[str, Any], composite: RawTransport) -> None:
    raw = inst["raw"]
    bl = beamlets_of(raw)
    assert [a.beamlet_id for a in bl.accumulations] == [0, 1, 2, 3]
    assert all(a.batch_quanta.shape == (4, 90000) for a in bl.accumulations)
    assert np.array_equal(sum_columns(bl), quanta(composite))  # (a)
    n_total = 800
    comp = raw.edep_mev[0].sum(axis=0) / n_total  # (b) per primary
    cols = sum(
        (200 / n_total) * (a.batch_quanta.sum(axis=0) * QUANTUM_MEV / 200) for a in bl.accumulations
    )
    nz = comp > 0
    assert np.all(np.abs(cols[nz] - comp[nz]) <= 1e-12 * comp[nz])
    assert np.array_equal(raw.edep_mev[0], composite.edep_mev[0])  # (c)
    assert raw.tallies == composite.tallies and raw.counters == composite.counters
    assert composite.beamlets is not None and composite.beamlets.accumulations == ()


def test_inf7_per_beamlet_closure_and_negative_controls(inst: dict[str, Any]) -> None:
    bl = beamlets_of(inst["raw"])
    assert max(closure_errors(bl)) < 1e-12
    assert [x.n_histories for x in bl.balances] == [200] * 4
    for x in bl.balances:
        assert abs(x.closure_residual_mev) < 1e-9 * x.tallies["initial"]
    # columns permuted: the identity must fail
    a = bl.accumulations
    perm = tuple(
        BeamletAccumulation(a[i].beamlet_id, a[(i + 1) % 4].batch_quanta) for i in range(4)
    )
    bad = dataclasses.replace(bl, accumulations=perm)
    assert min(closure_errors(bad)) > 1e-6


def _small(**kw: Any) -> SimulationConfig:
    return rinf(n=8, **kw)


def test_inf6_negative_controls_and_inf7(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg_i = _small()
    cfg_c = _small(influence=None)
    base = transport(cfg_i)
    ref = quanta(transport(cfg_c))
    assert np.array_equal(sum_columns(beamlets_of(base)), ref)
    assert max(closure_errors(beamlets_of(base))) < 1e-12
    cls = reference._Reference

    # (i) c0 reset per beamlet: BEAM-4 fails
    offs = [int(x) for x in cfg_i.source.offsets()]
    pairs: list[tuple[int, int]] = []
    orig = reference.draw_block

    def reset(key: Any, c0: int, c1: int, c2: int, c3: int) -> Any:
        if c3 == PURPOSE_SOURCE:
            c0 = c0 - max(s for s in offs if s <= c0)
            pairs.append((c0, c3))
        return orig(key, c0, c1, c2, c3)

    with monkeypatch.context() as mp:
        mp.setattr(reference, "draw_block", reset)
        transport(cfg_i)
    assert sorted(pairs) != [(h, PURPOSE_SOURCE) for h in range(32)]
    assert len(set(pairs)) < 32

    # (ii) accumulator rows rotated with correct source: INF-2(a) passes, INF-7 fails
    def rotated(self: Any, jl: int, batch: int) -> int:
        g = self.eff.influence_group_size
        return ((jl + 1) % g) * self.n_batches + batch

    with monkeypatch.context() as mp:
        mp.setattr(cls, "_influence_row", rotated)
        bad = beamlets_of(transport(cfg_i))
    assert np.array_equal(sum_columns(bad), ref)
    assert max(closure_errors(bad)) > 1e-6

    # accumulator row off by one (spills into the neighbour's block): INF-7 fails
    def off_by_one(self: Any, jl: int, batch: int) -> int:
        g = self.eff.influence_group_size
        return (jl * self.n_batches + batch + 1) % (g * self.n_batches)

    with monkeypatch.context() as mp:
        mp.setattr(cls, "_influence_row", off_by_one)
        bad = beamlets_of(transport(cfg_i))
    assert max(closure_errors(bad)) > 1e-6

    # (ii) beamlet id off by one at a boundary: INF-7 fails
    real_id = cls._beamlet_id

    def boundary(self: Any, h: int) -> int:
        j = real_id(self, h)
        return j - 1 if j > 0 and h == self.offsets_list[j] else j

    with monkeypatch.context() as mp:
        mp.setattr(cls, "_beamlet_id", boundary)
        bad = beamlets_of(transport(cfg_i))
    assert max(closure_errors(bad)) > 1e-6

    # (iii) accumulator stride B - 1: INF-2(a) fails; stride B + 1 overruns the rows
    def stride(delta: int) -> Any:
        def row(self: Any, jl: int, batch: int) -> int:
            return jl * (self.n_batches + delta) + batch

        return row

    with monkeypatch.context() as mp:
        mp.setattr(cls, "_influence_row", stride(-1))
        bad = beamlets_of(transport(cfg_i))
    assert not np.array_equal(sum_columns(bad), ref)
    with monkeypatch.context() as mp:
        mp.setattr(cls, "_influence_row", stride(+1))
        with pytest.raises(IndexError):
            transport(cfg_i)


def test_inf_group_size_does_not_change_columns() -> None:
    a = beamlets_of(transport(_small(influence=InfluenceOptions("dose", group_size=1))))
    b = beamlets_of(transport(_small(influence=InfluenceOptions("dose", group_size=3))))
    assert (a.group_size, b.group_size) == (1, 3)
    for x, y in zip(a.accumulations, b.accumulations, strict=True):
        assert x.beamlet_id == y.beamlet_id and np.array_equal(x.batch_quanta, y.batch_quanta)
    assert a.balances == b.balances


# ---- BEAM-6 (nuclear, small) --------------------------------------------------------------------


def test_beam6_nuclear_per_beamlet_balance_and_sum_identity() -> None:
    from test_nuclear_transport_loop import _table_id

    tid = _table_id()
    kw: dict[str, Any] = {
        "n": 20,
        "energies": (150.0, 200.0),
        "nuclear": True,
        "tid": tid,
        "seed": 70001,
    }
    c = rinf(**kw)
    bs = BeamletSet(PROTON, c.source.beamlets[:2])  # type: ignore[union-attr]
    c = dataclasses.replace(c, source=bs)
    raw = transport(c)
    comp = transport(dataclasses.replace(c, influence=None))
    bl = beamlets_of(raw)
    assert np.array_equal(sum_columns(bl), quanta(comp))
    assert np.array_equal(raw.edep_mev[0], comp.edep_mev[0])
    assert raw.tallies == comp.tallies
    assert max(closure_errors(bl)) < 1e-9
    for x in bl.balances:
        assert abs(x.closure_residual_mev) <= 1e-9 * x.tallies["initial"]
        assert "nuclear_binding" in x.tallies
    assert sum(x.tallies["nuclear_local"] for x in bl.balances) > 0.0
    for k in raw.tallies:
        tot = sum(x.tallies[k] for x in bl.balances)
        assert abs(tot - raw.tallies[k]) <= 1e-9 * max(abs(raw.tallies[k]), 1.0)


# ---- MEM-1 (python) -----------------------------------------------------------------------------


def _estimate(eff: Any, g: int) -> int:
    cfg = eff.requested
    return allocation_estimate(
        "python", cfg.run.n_batches, GRID.n_voxels, 1, eff, eff.n_histories,
        cfg.run.chunk_histories, GRID.n_voxels, g, [int(x) for x in cfg.source.offsets()],
    )  # fmt: skip


@pytest.mark.parametrize("nuclear", [False, True])
@pytest.mark.parametrize("g", [1, 2, 4])
def test_mem1_python_measured_peak_matches_the_estimate(g: int, nuclear: bool) -> None:
    tid = None
    if nuclear:
        from test_nuclear_transport_loop import _table_id

        tid = _table_id()
    c = rinf(n=4, influence=InfluenceOptions("dose", group_size=g), nuclear=nuclear, tid=tid,
             energies=(60.0, 70.0, 80.0, 90.0))  # fmt: skip
    eff = validate(c)
    with record_allocations() as rec:
        run_transport(eff)
    est = _estimate(eff, g)
    assert rec.peak_bytes <= est <= 1.5 * rec.peak_bytes
    kinds = {(k, t) for k, t, _ in rec.events}
    assert {("alloc", "tally_rows"), ("release", "tally_rows"), ("release", "influence")} <= kinds
    assert rec.live_bytes == 8 * 4 * GRID.n_voxels  # only the grid accumulators stay live


def test_mem1_group_j_on_the_r_inf_run_and_a_lingering_array_fails_the_audit(
    inst: dict[str, Any],
) -> None:
    est = _estimate(inst["eff"], 4)
    assert inst["rec"].peak_bytes <= est <= 1.5 * inst["rec"].peak_bytes
    with record_allocations() as rec:
        alloc_zeros("rogue", (1000,), np.int64)  # a lingering payload copy
        run_transport(validate(rinf(n=4)))
    assert rec.peak_bytes > _estimate(validate(rinf(n=4)), 4)


# ---- Result level (S1 fix 5 note) ---------------------------------------------------------------


def test_result_of_a_toml_without_n_histories_carries_one_derived_row(tmp_path: Path) -> None:
    from ionmc import configio

    base = (Path(__file__).parent / "data" / "v3_006" / "proton_box.toml").read_text()
    head = base.split("[geometry]")[0].split("[source]")[0]
    rest = "[geometry]" + base.split("[geometry]")[1]
    rest = re.sub(r"n_histories = \d+\n", "", rest)
    toml = (
        head
        + '[source]\ntype = "beamlets"\nprojectile = "proton"\n'
        + "[[source.beamlets]]\nposition_mm = [0.0, 0.0, -5.0]\ndirection = [0.0, 0.0, 1.0]\n"
        + "kinetic_energy_mev = 20.0\nn_histories = 4\n"
        + "[[source.beamlets]]\nposition_mm = [1.0, 0.0, -5.0]\ndirection = [0.0, 0.0, 1.0]\n"
        + "kinetic_energy_mev = 21.0\nn_histories = 8\n\n"
        + '[influence]\ngrid = "dose"\n\n'
        + rest
    )
    p = tmp_path / "b.toml"
    p.write_text(toml)
    cfgd = configio.load_toml(p)[0]
    assert cfgd.run.n_histories is None
    sim = Simulation(cfgd)
    res = sim.run()
    rows = [r for r in sim.resolutions if "run.n_histories" in r.requested]
    assert len(rows) == 1 and (rows[0].cls, rows[0].rule) == ("derived", "sum-of-beamlets")
    assert res.n_histories == res.effective_config.n_histories == 12
    assert res.beamlet_raw is not None and len(res.beamlet_raw.accumulations) == 2
    assert os.path.exists(p)
