"""V3-004 Warp channel scoring: A2 (Warp half), A11 (CI), A15, A16 on warp-cpu and the A11-HR
comparison function (plan validation/plans/v3-004-acceptance.md). Deliberately without
``from __future__ import annotations``: the Warp test kernel needs real annotations."""

import math
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import warp as wp

import ionmc.config as config_module
from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    validate,
)
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import VoxelGeometry
from ionmc.lookup import LookupTable
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport.channel_device import make_channel_data
from ionmc.transport.channels import CLASS_LOCAL, CLASS_STEP
from ionmc.transport.kernels import make_kernel_support
from ionmc.transport.parity_channels import (
    channel_t12_observables,
    compare_channel_partition,
    compare_channel_runs,
)
from ionmc.transport.run import run_transport
from ionmc.transport.scoring_ref import ReferenceChannelScorer
from ionmc.transport.tally import N_FIXED_TALLIES

SYNTHETIC = Path(__file__).resolve().parents[1] / "data" / "synthetic_lookup.json"
LK_LET = LookupTable.from_file(SYNTHETIC)
LK_E = LookupTable(
    "e_lin", "q", "1", "energy_per_nucleon_mev", "log", np.geomspace(0.5, 200.0, 50),
    {"proton": np.linspace(1.0, 3.0, 50)}, "c", "l", "s", True,
)  # fmt: skip
EDGES = tuple(np.geomspace(2.0, 110.0, 9))


def _tallies(grid: str = "dose") -> tuple[TallyRequest, ...]:
    g = grid
    return (
        TallyRequest("edep", g, "edep"),
        TallyRequest("lt", g, "let_t"),
        TallyRequest("ld", g, "let_d"),
        TallyRequest("le", g, "let_d_eps"),
        TallyRequest("fl", g, "fluence"),
        TallyRequest("fe", g, "lookup_sum", lookup=LK_LET.name),
        TallyRequest("fa", g, "lookup_dose_avg", lookup="e_lin"),
        TallyRequest("sp", g, "fluence_spectrum", energy_edges_mev_per_u=EDGES),
        TallyRequest("lt_p", g, "let_t", species=("proton",), generation="primary"),
    )


def _cfg(
    backend: str,
    precision: str = "float64",
    *,
    n: int = 16,
    nb: int = 2,
    seed: int = 20261004,
    energy: float = 100.0,
    tallies: tuple[TallyRequest, ...] | None = None,
    chunk: int | None = None,
    workers: int = 1,
    diag: DiagnosticsOptions | None = None,
    grid_nz: int = 80,
    sigma: tuple[float, float] = (0.5, 1.0),
) -> SimulationConfig:
    shape = (12, 12, 32)
    geo = VoxelGeometry(
        origin_mm=(-30.0, -30.0, 0.0), spacing_mm=(5.0,) * 3, shape=shape, materials=(WATER,),
        material_index=np.zeros(shape, dtype=np.int32),
    )  # fmt: skip
    extra: dict[str, Any] = {} if chunk is None else {"chunk_histories": chunk}
    return SimulationConfig(
        source=PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), energy, *sigma),
        geometry=geo,
        scoring=(
            ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, grid_nz), name="dose"),
        ),
        physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource(), max_step_mm=2.0),
        run=RunOptions(
            backend=backend, precision=precision, n_histories=n, n_batches=nb, seed=seed,
            cpu_workers=workers, **extra,
        ),  # type: ignore[arg-type]
        tallies=_tallies() if tallies is None else tallies,
        lookups=() if tallies == () else (LK_LET, LK_E),
        diagnostics=diag or DiagnosticsOptions(),
    )  # fmt: skip


# --------------------------------------------------------------------------- A11 (CI)


def test_a11_ci_python_vs_warp_cpu_float64_all_channels() -> None:
    a = Simulation(_cfg("python")).run()
    b = Simulation(_cfg("warp-cpu")).run()
    plan = a.effective_config.channels
    assert plan is not None
    n_ci = {c.grid: i for i, c in enumerate(plan.channels) if c.kind == "N"}
    n_per = a.channel_batches(n_ci[0]) * plan.channels[n_ci[0]].quantum  # [B, nvox]
    assert np.array_equal(a.channel_batches(n_ci[0]), b.channel_batches(n_ci[0]))
    assert n_per.sum() > 1000
    for ci, c in enumerate(plan.channels):
        x = a.channel_batches(ci) * c.quantum
        y = b.channel_batches(ci) * c.quantum
        nv = np.repeat(n_per, c.size // n_per.shape[1], axis=1)
        bound = 1e-10 * np.maximum(np.abs(x), np.abs(y)) + nv * c.quantum
        assert np.all(np.abs(x - y) <= bound), (ci, c.kind)
    assert a.channel_raw is not None and b.channel_raw is not None
    assert a.channel_raw.lookup_out_of_domain == b.channel_raw.lookup_out_of_domain == 0
    np.testing.assert_allclose(
        a.channel_raw.residual, b.channel_raw.residual, rtol=1e-9, atol=1e-18
    )
    qa, qb = a.grids[0].quantities, b.grids[0].quantities
    assert set(qa) == set(qb)
    for name in qa:
        assert np.array_equal(qa[name].n_nonzero, qb[name].n_nonzero), name
        assert np.array_equal(qa[name].defined_mask, qb[name].defined_mask), name


# --------------------------------------------------------------------------- A15 / A16


@pytest.mark.parametrize("precision", ["float32", "float64"])
def test_a15_chunk_size_invariance_warp_cpu(precision: str) -> None:
    kw: dict[str, Any] = {"n": 2500, "nb": 5, "energy": 60.0, "seed": 20351004, "grid_nz": 40}
    small = Simulation(_cfg("warp-cpu", precision, chunk=2**10, **kw)).run()
    large = Simulation(_cfg("warp-cpu", precision, chunk=2**18, **kw)).run()
    assert small.transport_report["partials"][0]["n_chunks"] == 3
    v = compare_channel_partition(small, large)
    assert v["pass"], v
    assert small.channel_raw is not None and int(small.channel_raw.acc.sum()) > 0


@pytest.mark.multiprocess
def test_a15_one_vs_three_workers() -> None:
    kw: dict[str, Any] = {"n": 1200, "nb": 4, "energy": 60.0, "seed": 20351004, "grid_nz": 40}
    one = Simulation(_cfg("warp-cpu", "float64", workers=1, **kw)).run()
    three = Simulation(_cfg("warp-cpu", "float64", workers=3, **kw)).run()
    assert compare_channel_partition(one, three)["pass"]


@pytest.mark.parametrize("precision", ["float32", "float64"])
@pytest.mark.parametrize("n,nb,chunk", [(16, 2, None), (2500, 5, 2**10)])
def test_a16_tallies_do_not_change_qualified_outputs_warp_cpu(
    precision: str, n: int, nb: int, chunk: int | None
) -> None:
    diag = DiagnosticsOptions(
        track_end_positions=True, trace_histories=2 if precision == "float64" else 0
    )
    kw: dict[str, Any] = {"n": n, "nb": nb, "seed": 3, "chunk": chunk, "diag": diag,
                          "energy": 60.0, "grid_nz": 40}  # fmt: skip
    ra = run_transport(validate(_cfg("warp-cpu", precision, tallies=(), **kw)))
    rb = run_transport(validate(_cfg("warp-cpu", precision, **kw)))
    assert ra.channels is None and rb.channels is not None
    for ea, eb in zip(ra.edep_mev, rb.edep_mev, strict=True):
        assert np.array_equal(ea, eb)
    assert ra.tallies == rb.tallies and ra.counters == rb.counters
    assert ra.outside_mev == rb.outside_mev and ra.quantization_mev == rb.quantization_mev
    for k, va in ra.diagnostics.get("trace", {}).items():
        assert np.array_equal(va, rb.diagnostics["trace"][k]), k
    assert np.array_equal(ra.diagnostics["end_energy_mev"], rb.diagnostics["end_energy_mev"])


# --------------------------------------------------------------------------- A12


def test_a12_backend_without_channels_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    assert set(config_module.CHANNEL_BACKENDS) == {"python", "warp-cpu", "warp-cuda"}
    monkeypatch.setattr(config_module, "CHANNEL_BACKENDS", ("python",), raising=True)
    for prec in ("float32", "float64"):
        with pytest.raises(UnsupportedCombinationError, match="does not implement scoring"):
            validate(_cfg("warp-cpu", prec))
    validate(_cfg("python"))  # still accepted


# --------------------------------------------------------------------------- A2 (Warp half)

chan_t = make_kernel_support(wp.float64).chan
SUP = make_kernel_support(wp.float64)


@wp.kernel
def _a2_kernel(
    chan: chan_t,
    rows: wp.array2d(dtype=wp.float64),
    batch: wp.array(dtype=int),
    vox: wp.array(dtype=int),
    length: wp.array(dtype=wp.float64),
    tau: wp.array(dtype=wp.float64),
    eps: wp.array(dtype=wp.float64),
    species: wp.array(dtype=int),
    gen: wp.array(dtype=int),
    cls: wp.array(dtype=int),
    st: wp.array(dtype=wp.float64),
):
    i = wp.tid()
    SUP.score_piece(
        chan, rows, i, batch[i], 0, vox[i], length[i], tau[i], eps[i], species[i], gen[i], cls[i],
        st[0], st[1], st[2], st[3],
    )  # fmt: skip


def test_a2_warp_half_bitwise_vs_python_hook(make_config: Callable[..., SimulationConfig]) -> None:
    x = np.arange(65.0)
    lk = LookupTable(
        "dy", "q", "1", "let_water_kev_um", "linear", x,
        {"proton": 1.0 + 0.125 * x, "deuteron": 0.5 + 0.25 * x}, "c", "l", "s", True,
    )  # fmt: skip
    from fractions import Fraction  # noqa: F401  (documentation of exactness: dyadic values)

    from ionmc.materials import WATER as W
    from ionmc.physics.stopping import BetheStoppingSource as B
    from ionmc.transport.channels import compile_channels

    reqs: list[TallyRequest] = []
    for tag, sp in (("all", None), ("p", ("proton",)), ("d", ("deuteron",))):
        for g in ("all", "primary", "secondary"):
            kw: dict[str, Any] = {"species": sp, "generation": g}
            reqs += [
                TallyRequest(f"lt_{tag}_{g}", "dose", "let_t", **kw),
                TallyRequest(f"ld_{tag}_{g}", "dose", "let_d", **kw),
                TallyRequest(f"le_{tag}_{g}", "dose", "let_d_eps", **kw),
                TallyRequest(f"fe_{tag}_{g}", "dose", "lookup_sum", lookup="dy", **kw),
                TallyRequest(f"ed_{tag}_{g}", "dose", "edep", **kw),
            ]
    cfg = make_config(
        energy=60.0,
        n=40,
        n_batches=4,
        scoring=(ScoringGrid((-10.0, -10.0, 0.0), (20.0, 20.0, 2.0), (1, 1, 4), name="dose"),),
    )
    eff = validate(cfg)
    producible = frozenset((s, g) for s in ("proton", "deuteron") for g in ("primary", "secondary"))
    plan = compile_channels(
        tuple(reqs), (lk,), cfg.scoring, geometry=eff.geometry, tables=eff.tables,
        water=B().table(W, PROTON), projectile=PROTON, e_cut_mev=2.0, e_hi_mev=60.0,
        n_histories=40, n_batches=4, cpu_workers=1, memory_budget_bytes=2**31, max_steps=1000,
        scoring_pieces=8, producible=producible,
    )  # fmt: skip
    # stream (dyadic values, exact multiples of every quantum)
    rng = np.random.default_rng(20351004)
    npc = 300
    pb = rng.integers(0, 4, npc)
    pv = rng.integers(0, 4, npc)
    pl = rng.choice([1.0, 0.5, 0.25], npc)
    pt = rng.integers(-8, 9, npc) / 8.0
    pe = rng.integers(1, 64, npc) / 64.0
    ps = rng.integers(0, 2, npc)
    pg = rng.integers(0, 2, npc)
    local = [(1, 0, 0.375, 0, 0), (2, 3, 0.125, 1, 1)]
    state = (4.0, 0.75, 32.0, 0.5)
    ref = ReferenceChannelScorer(plan, cfg.scoring, 4, tables=eff.tables)
    ref.set_step_state(*state)
    for i in range(npc):
        ref.score_piece(int(pb[i]), 0, int(pv[i]), float(pl[i]), float(pt[i]), float(pe[i]),
                        int(ps[i]), int(pg[i]), CLASS_STEP)  # fmt: skip
    for b, v, e, sp, g in local:
        ref.score_piece(b, 0, v, 0.0, 0.0, e, sp, g, CLASS_LOCAL)

    # the same stream through the kernel-support score_piece on warp-cpu
    n_all = npc + len(local)
    cat = np.concatenate
    arrs = {
        "batch": cat([pb, [t[0] for t in local]]).astype(np.int64),
        "vox": cat([pv, [t[1] for t in local]]).astype(np.int64),
        "length": cat([pl, np.zeros(2)]),
        "tau": cat([pt, np.zeros(2)]),
        "eps": cat([pe, [t[2] for t in local]]),
        "species": cat([ps, [t[3] for t in local]]).astype(np.int64),
        "gen": cat([pg, [t[4] for t in local]]).astype(np.int64),
        "cls": cat([np.full(npc, CLASS_STEP), np.full(2, CLASS_LOCAL)]).astype(np.int64),
    }
    tables_w = validate(replace(cfg, tallies=(TallyRequest("x", "dose", "let_t"),))).tables
    chan = make_channel_data(
        SUP, plan, tables_w, n_batches=4, n_grids=1, a_nucleon=1, species_id=0, generation=0,
        device="cpu",
    )  # fmt: skip
    chan.res_base = 0
    rows = wp.zeros((n_all, plan.n_residual + 1), dtype=wp.float64, device="cpu")

    def dev(a: np.ndarray, t: Any) -> Any:
        return wp.array(a, dtype=t, device="cpu")

    wp.launch(
        _a2_kernel, dim=n_all,
        inputs=[chan, rows, dev(arrs["batch"], int), dev(arrs["vox"], int),
                dev(arrs["length"], wp.float64), dev(arrs["tau"], wp.float64),
                dev(arrs["eps"], wp.float64), dev(arrs["species"], int), dev(arrs["gen"], int),
                dev(arrs["cls"], int), dev(np.array(state), wp.float64)],
        device="cpu",
    )  # fmt: skip
    got = chan.acc.numpy()
    assert np.array_equal(got, ref.acc)  # bitwise equal int64 accumulators
    res = rows.numpy().sum(axis=0)
    assert np.array_equal(res[:-1], ref.residual)
    assert res[-1] == ref.lookup_ood == 0
    assert got.sum() > 0


# --------------------------------------------------------------------------- A11-HR function


def _idd_cfg(seed: int) -> SimulationConfig:
    base = _cfg("warp-cpu", "float32", n=1200, nb=12, seed=seed, energy=60.0,
                sigma=(0.0, 0.0), tallies=())  # fmt: skip
    grid = ScoringGrid((-30.0, -30.0, 0.0), (60.0, 60.0, 1.0), (1, 1, 40), name="idd")
    tallies = tuple(
        replace(t, name=t.name) for t in _tallies("idd") if t.name in ("lt", "ld", "le", "fe", "sp")
    )
    return replace(base, scoring=(grid,), tallies=tallies, lookups=(LK_LET,))


def test_a11_hr_comparison_function_cpu_only() -> None:
    a = Simulation(_idd_cfg(20351004)).run()
    b = Simulation(_idd_cfg(20351005)).run()
    obs = channel_t12_observables(a)
    assert set(obs.arrays) == {"L", "LS", "LS2", "ES", "E_step", "FE", "spectrum"}
    assert set(obs.scalars) == {"L_total", "LS_total", "LS2_total", "ES_total"}
    v = compare_channel_runs(a, b, n_perm=999, n_boot=999)
    assert set(v["t12"]["arrays"]) == set(obs.arrays)
    assert math.isfinite(v["ratio_z"]["let_d"]["max_abs_z"])
    assert all(not r["gating"] for r in v["ratio_z"].values())
    assert v["pass"] == v["t12"]["pass"]
    # negative control: a 5 % bias of LS in one sample must fail the pair
    from ionmc.transport.parity import t12_compare

    oa, ob = channel_t12_observables(a), channel_t12_observables(b)
    biased = replace(ob, arrays={**ob.arrays, "LS": ob.arrays["LS"] * 1.05},
                     scalars={**ob.scalars, "LS_total": ob.scalars["LS_total"] * 1.05})  # fmt: skip
    assert not t12_compare(oa, biased, n_perm=999, n_boot=999)["pass"]
    # one-sample partition comparison sanity: identical runs are bit-identical
    assert compare_channel_partition(a, a)["pass"]
    assert N_FIXED_TALLIES == 6


# --------------------------------------------------------------------------- CUDA (HR host)


def _cuda_kw() -> dict[str, Any]:
    return {"n": 20000, "nb": 4, "energy": 60.0, "seed": 9, "grid_nz": 40}


@pytest.mark.cuda
@pytest.mark.parametrize("precision", ["float32", "float64"])
def test_cuda_a15_chunk_size_invariance_of_channels(precision: str) -> None:
    """A15 (HR, strict): chunk sizes 2^10 and 2^18 on warp-cuda give bit-identical channels,
    residuals, counters and deposits."""
    s = Simulation(_cfg("warp-cuda", precision, chunk=2**10, **_cuda_kw())).run()
    t = Simulation(_cfg("warp-cuda", precision, chunk=2**18, **_cuda_kw())).run()
    v = compare_channel_partition(s, t)
    assert v["pass"], v


@pytest.mark.cuda
@pytest.mark.parametrize("precision", ["float32", "float64"])
def test_cuda_a5_integer_identities_and_internal_consistency(precision: str) -> None:
    """A5 on CUDA: the species and generation E channels equal edep per batch exactly (integer
    identity), ``E_step + edep_excluded_from_let = edep``; channel totals are consistent with the
    run's own tallies (N > 0 where edep > 0, residual columns finite)."""
    tallies = (
        TallyRequest("edep", "dose", "edep"),
        TallyRequest("edep_p", "dose", "edep", species=("proton",)),
        TallyRequest("edep_prim", "dose", "edep", generation="primary"),
        TallyRequest("e_step", "dose", "let_d_eps"),
        TallyRequest("lt", "dose", "let_t"),
    )
    r = Simulation(_cfg("warp-cuda", precision, tallies=tallies, **_cuda_kw())).run()
    plan = r.effective_config.channels
    assert plan is not None and r.channel_raw is not None
    q = {d.name: d for d in plan.quantities}
    e_all = r.channel_batches(q["edep"].numerator)
    assert np.array_equal(r.channel_batches(q["edep_p"].numerator), e_all)
    assert np.array_equal(r.channel_batches(q["edep_prim"].numerator), e_all)
    local = next(i for i, c in enumerate(plan.channels) if c.class_mask == CLASS_LOCAL)
    e_step = r.channel_batches(q["e_step"].denominator)
    assert np.array_equal(e_step + r.channel_batches(local), e_all)
    n_ci = next(i for i, c in enumerate(plan.channels) if c.kind == "N")
    assert np.all(r.channel_batches(n_ci)[e_step > 0] > 0)
    assert np.all(np.isfinite(r.channel_raw.residual))
    assert r.channel_raw.lookup_out_of_domain == 0


@pytest.mark.cuda
def test_cuda_vs_python_float64_informative_record() -> None:
    """NON-gating: per-voxel python vs warp-cuda float64 differences are printed for the HR
    archive. CPU and CUDA float64 trajectories are not bitwise identical (ulp-level differences
    in transcendental functions amplify through sampling and voxel crossings), so channel parity
    across devices is statistical (A11 HR, ``compare_channel_runs``)."""
    a = Simulation(_cfg("python")).run()
    b = Simulation(_cfg("warp-cuda", "float64")).run()
    plan = a.effective_config.channels
    assert plan is not None
    for ci, c in enumerate(plan.channels):
        d = np.abs(a.channel_batches(ci) - b.channel_batches(ci))
        print(
            f"INFO channel {ci} {c.kind}: max |python - cuda| = {int(d.max())} quanta, "
            f"{int((d > 0).sum())} of {d.size} differing"
        )
