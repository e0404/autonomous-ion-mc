"""V3-005B C10: NuclearDevice upload and the event sampler as a Warp-scope function (P5 extension).

* the packed table uploads to ``wp.array`` objects of the kernel precision on an explicit device
  and its ``sha256`` identity is reproduced from the bytes read back from the device;
* event-level parity (V8 precondition): 1e4 recorded inputs (table target, T1, Philox counter
  ``(h, gid, base block)``; fixed CI seed ``INPUT_SEED``) are sampled by the pure-Python twin of
  ``sample_event`` and by a Warp CPU kernel (float64: identical integer outputs and species order,
  continuous outputs within 1e-12 relative; float32 kernel against the float64 twin fed with the
  float32 uniforms and the float32-rounded table: discrete mismatches counted, continuous maxima
  reported);
* ``choose_target`` against the reference backend's selection loop;
* the numpy batch sampler of the table builder against the shared sampler (one algorithm).

Needs a qualified nuclear table in the cache (``IONMC_CACHE_DIR``; ``IONMC_REQUIRE_DATA=1`` turns
its absence into a failure).
"""

import os
from typing import Any

import numpy as np
import pytest
import warp as wp

from ionmc._wpfunc import python_twin
from ionmc.data import cache
from ionmc.materials import PMMA, WATER
from ionmc.nuclear.tables import NuclearTable, NuclearTableError
from ionmc.physics.nuclear import (
    EVENT_ACCEPTED,
    EVENT_BLOCKS_PER_ATTEMPT,
    EVF_STRIDE,
    EVI_STRIDE,
    MAX_PRODUCTS,
    PROD_STRIDE,
    make_nuclear,
)
from ionmc.rng.philox import key_from_seed, philox4x32_10_py, u01_py
from ionmc.transport.nuclear_device import NuclearDevice, NuclearHost, pack_nuclear
from ionmc.transport.reference import EVENT_BLOCKS_PER_ATTEMPT as REFERENCE_BLOCKS

wp.config.log_level = wp.LOG_WARNING
DEV = "cpu"
INPUT_SEED = 20261008  # recorded inputs of the event-level parity test (CI)
KEY_SEED = 987654321  # Philox key of the uniform stream
N_EVENTS = int(os.environ.get("IONMC_C10_N", "10000"))
NU64 = python_twin(make_nuclear)


def _table() -> NuclearTable:
    cdir = cache.resolve_cache_dir(None)
    best: tuple[bool, str] | None = None
    for p in sorted((cdir / "derived").glob("nuclear-proton-*.json")):
        tid = p.stem.removeprefix("nuclear-proton-")
        try:
            info = NuclearTable.load(None, tid).info
        except NuclearTableError:
            continue
        if "transport_path_bound_terms" in info:
            cand = (info["options"]["diagnostic_nodes_mev"] is None, tid)
            best = cand if best is None else max(best, cand)
    if best is not None:
        return NuclearTable.load(None, best[1])
    if os.environ.get("IONMC_REQUIRE_DATA") == "1":
        pytest.fail("no loadable nuclear table in the cache")
    pytest.skip("no built nuclear table in the cache (IONMC_CACHE_DIR)")


@pytest.fixture(scope="module")
def host() -> NuclearHost:
    return pack_nuclear(_table(), (WATER, PMMA))


def test_c10_namespace_has_the_event_functions() -> None:
    assert {"sample_event", "choose_target"} <= set(NU64.__dict__)
    assert {"sample_event", "choose_target"} <= set(make_nuclear(wp.float64).__dict__)
    assert EVENT_BLOCKS_PER_ATTEMPT == REFERENCE_BLOCKS == 164


@pytest.mark.parametrize("real", [wp.float64, wp.float32])
def test_nuclear_device_upload_identity(host: NuclearHost, real: Any) -> None:
    dev = NuclearDevice(host, real, DEV)
    assert dev.device == DEV and dev.n_grid == host.n_grid and dev.n_targets == host.n_targets
    assert dev.grid.dtype == real and dev.m_res.dtype == real and dev.mat_target.dtype == wp.int32
    assert dev.grid.device.alias == DEV
    assert dev.sha256 == dev.device_sha256()  # bytes read back from the device
    back = dev.readback()
    np_real = np.float64 if real is wp.float64 else np.float32
    for k, v in back.items():
        if k.startswith("mat_"):
            assert np.array_equal(v, host.arrays[k].astype(np.int32)), k
        else:
            assert np.array_equal(v, host.arrays[k].astype(np_real)), k
    again = NuclearDevice(pack_nuclear(_table(), (WATER, PMMA)), real, DEV)
    assert again.sha256 == dev.sha256
    other = NuclearDevice(host, wp.float32 if real is wp.float64 else wp.float64, DEV)
    assert other.sha256 != dev.sha256
    only_water = NuclearDevice(pack_nuclear(_table(), (WATER,)), real, DEV)
    assert only_water.sha256 != dev.sha256
    assert dev.bounds["history_energy_bound_mev"] > 0.0


# ---------------------------------------------------------------------------------------------
# recorded inputs
# ---------------------------------------------------------------------------------------------
def _inputs(host: NuclearHost) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(INPUT_SEED)
    n = N_EVENTS
    return {
        "h": rng.integers(0, 2**32, n, dtype=np.uint64).astype(np.uint32),
        "gid": rng.integers(0, 2**30, n, dtype=np.uint64).astype(np.uint32),
        "base": rng.integers(0, 2000, n).astype(np.int32),
        "tgt": rng.integers(0, host.n_targets, n).astype(np.int32),
        "t1": np.exp(rng.uniform(np.log(1.5), np.log(250.0), n)),
    }


def _run_twin(
    host: NuclearHost, inp: dict[str, np.ndarray], *, precision: str
) -> dict[str, np.ndarray]:
    """The pure-Python twin on every input. ``precision='float32'``: the float32 uniform stream
    and the float32-rounded table and T1 (so that only the arithmetic differs from the kernel)."""
    n = inp["h"].size
    k0, k1 = key_from_seed(KEY_SEED)
    arr = {k: np.asarray(v) for k, v in host.arrays.items()}
    t1 = inp["t1"]
    if precision == "float32":
        arr = {k: v.astype(np.float32).astype(np.float64) if v.dtype.kind == "f" else v
               for k, v in arr.items()}  # fmt: skip
        t1 = t1.astype(np.float32).astype(np.float64)

        def key(h: int, gid: int, blk: int, word: int) -> float:
            return u01_py(philox4x32_10_py((h, gid, blk, 2), (k0, k1))[word], "float32")

        key_arg: Any = key
    else:
        key_arg = (k0, k1)
    evi = np.zeros(n * EVI_STRIDE, dtype=np.int64)
    evf = np.zeros(n * EVF_STRIDE)
    prod = np.zeros(n * MAX_PRODUCTS * PROD_STRIDE)
    status = np.zeros(n, dtype=np.int64)
    for i in range(n):
        status[i] = NU64.sample_event(
            int(inp["h"][i]), int(inp["gid"][i]), int(inp["base"][i]), key_arg,
            int(inp["tgt"][i]), float(t1[i]), arr["grid"], host.n_grid, arr["lam"],
            arr["edges"], arr["rpre"], arr["recoil"], arr["tconst"], arr["m_res"],
            evi, evf, prod, i,
        )  # fmt: skip
    return {"status": status, "evi": evi, "evf": evf, "prod": prod}


def _kernel(real: Any) -> Any:
    nu = make_nuclear(real)

    @wp.kernel(module="unique")
    def kernel(
        h: wp.array(dtype=wp.uint32),  # type: ignore[valid-type]
        gid: wp.array(dtype=wp.uint32),  # type: ignore[valid-type]
        base: wp.array(dtype=int),  # type: ignore[valid-type]
        tgt: wp.array(dtype=int),  # type: ignore[valid-type]
        t1: wp.array(dtype=real),  # type: ignore[valid-type]
        key: wp.vec2ui,
        n_grid: int,
        grid: wp.array(dtype=real),  # type: ignore[valid-type]
        lam: wp.array(dtype=real),  # type: ignore[valid-type]
        edges: wp.array(dtype=real),  # type: ignore[valid-type]
        rpre: wp.array(dtype=real),  # type: ignore[valid-type]
        recoil: wp.array(dtype=real),  # type: ignore[valid-type]
        tconst: wp.array(dtype=real),  # type: ignore[valid-type]
        m_res: wp.array(dtype=real),  # type: ignore[valid-type]
        evi: wp.array(dtype=int),  # type: ignore[valid-type]
        evf: wp.array(dtype=real),  # type: ignore[valid-type]
        prod: wp.array(dtype=real),  # type: ignore[valid-type]
        status: wp.array(dtype=int),  # type: ignore[valid-type]
    ) -> None:
        i = wp.tid()
        status[i] = nu.sample_event(
            h[i], gid[i], base[i], key, tgt[i], t1[i], grid, n_grid, lam, edges, rpre, recoil,
            tconst, m_res, evi, evf, prod, i,
        )  # fmt: skip

    return kernel


def _run_kernel(host: NuclearHost, inp: dict[str, np.ndarray], real: Any) -> dict[str, np.ndarray]:
    n = inp["h"].size
    d = NuclearDevice(host, real, DEV)
    np_real = np.float64 if real is wp.float64 else np.float32
    evi = wp.zeros(n * EVI_STRIDE, dtype=int, device=DEV)
    evf = wp.zeros(n * EVF_STRIDE, dtype=real, device=DEV)
    prod = wp.zeros(n * MAX_PRODUCTS * PROD_STRIDE, dtype=real, device=DEV)
    status = wp.zeros(n, dtype=int, device=DEV)
    k0, k1 = key_from_seed(KEY_SEED)
    wp.launch(
        _kernel(real),
        dim=n,
        inputs=[
            wp.array(inp["h"], dtype=wp.uint32, device=DEV),
            wp.array(inp["gid"], dtype=wp.uint32, device=DEV),
            wp.array(inp["base"], dtype=int, device=DEV),
            wp.array(inp["tgt"], dtype=int, device=DEV),
            wp.array(inp["t1"].astype(np_real), dtype=real, device=DEV),
            wp.vec2ui(k0, k1), d.n_grid, d.grid, d.lam, d.edges, d.rpre, d.recoil, d.tconst,
            d.m_res, evi, evf, prod, status,
        ],
        device=DEV,
    )  # fmt: skip
    return {
        "status": status.numpy(),
        "evi": evi.numpy().astype(np.int64),
        "evf": evf.numpy().astype(np.float64),
        "prod": prod.numpy().astype(np.float64),
    }


def _compare(
    inp: dict[str, np.ndarray],
    host: NuclearHost,
    a: dict[str, np.ndarray],
    b: dict[str, np.ndarray],
    *,
    floor_rel: float,
) -> dict[str, Any]:
    """Discrete mismatches (events) and the maximum continuous differences (a = twin, b = kernel)
    over the events whose discrete outputs agree. Difference measure ``|a - b| / max(|b|, ref)``
    with ``ref`` = ``floor_rel`` x the event energy scale (T1 + m_p + m_t) for energies and the
    ledger columns (their values are differences of ~1e4 MeV masses), ``floor_rel`` x the lab
    momentum magnitude for momentum components, ``floor_rel`` for ``mu`` and ``phi``. float64:
    ``floor_rel = 1e-12`` (the 1e-12 relative class); float32: ``floor_rel = 1e-6`` (decision 0039
    float32 class)."""
    n = inp["h"].size
    ia, ib = a["evi"].reshape(n, EVI_STRIDE), b["evi"].reshape(n, EVI_STRIDE)
    pa = a["prod"].reshape(n, MAX_PRODUCTS, PROD_STRIDE)
    pb = b["prod"].reshape(n, MAX_PRODUCTS, PROD_STRIDE)
    acc = a["status"] == EVENT_ACCEPTED
    n_prod = np.where(acc, ia[:, 8], 0)
    jj = np.arange(MAX_PRODUCTS)[None, :] < n_prod[:, None]
    same = (a["status"] == b["status"]) & np.all(ia[:, :10] == ib[:, :10], axis=1)
    same &= ~acc | np.all(~jj | (pa[:, :, 0] == pb[:, :, 0]), axis=1)
    ok = same & acc
    tc = host.arrays["tconst"].reshape(host.n_targets, -1)
    scale = inp["t1"] + tc[inp["tgt"], 0] + tc[inp["tgt"], 1]
    fa, fb = a["evf"].reshape(n, EVF_STRIDE), b["evf"].reshape(n, EVF_STRIDE)
    worst: dict[str, float] = {}
    for c, name in enumerate(("recoil", "imbalance", "binding", "local", "alpha_t", "m_res")):
        d = np.abs(fa[ok, c] - fb[ok, c]) / np.maximum(np.abs(fb[ok, c]), floor_rel * scale[ok])
        worst[name] = float(d.max(initial=0.0))
    mask = jj & ok[:, None]
    sc = np.broadcast_to(scale[:, None], mask.shape)
    p_lab = np.sqrt(sum(pb[:, :, c] ** 2 for c in (5, 6, 7)))
    np.seterr(invalid="ignore", divide="ignore")  # padded (unused) product rows are zero
    refs = {
        "e_prime": floor_rel * sc, "mu": floor_rel, "phi": floor_rel, "e_lab": floor_rel * sc,
        "px": floor_rel * p_lab, "py": floor_rel * p_lab, "pz": floor_rel * p_lab,
    }  # fmt: skip
    for c, name in enumerate(("species", "e_prime", "mu", "phi", "e_lab", "px", "py", "pz")):
        if c == 0:
            continue
        if name in ("px", "py", "pz"):  # relative to the momentum magnitude of the product
            d = np.abs(pa[:, :, c] - pb[:, :, c]) / np.maximum(p_lab, refs[name])
        else:
            d = np.abs(pa[:, :, c] - pb[:, :, c]) / np.maximum(np.abs(pb[:, :, c]), refs[name])
        worst[name] = float(d[mask].max(initial=0.0))
    return {
        "events": n,
        "accepted": int(acc.sum()),
        "discrete_mismatches": int((~same).sum()),
        "products": int(n_prod[ok].sum()),
        "max_rel": worst,
        "max_rel_all": max(worst.values()),
    }


@pytest.fixture(scope="module")
def recorded(host: NuclearHost) -> dict[str, np.ndarray]:
    return _inputs(host)


def test_p5_event_sampler_twin_equals_warp_cpu_float64(
    host: NuclearHost, recorded: dict[str, np.ndarray]
) -> None:
    a = _run_twin(host, recorded, precision="float64")
    b = _run_kernel(host, recorded, wp.float64)
    res = _compare(recorded, host, a, b, floor_rel=1e-12)
    print("P5 event-level float64:", res)
    assert res["accepted"] > 0.99 * res["events"] and res["products"] > res["events"]
    assert res["discrete_mismatches"] == 0
    assert res["max_rel_all"] <= 1e-12, res["max_rel"]


def test_p5_event_sampler_float32_kernel_against_float64_twin(
    host: NuclearHost, recorded: dict[str, np.ndarray]
) -> None:
    a = _run_twin(host, recorded, precision="float32")
    b = _run_kernel(host, recorded, wp.float32)
    res = _compare(recorded, host, a, b, floor_rel=1e-6)
    print("P5 event-level float32:", res)
    assert res["accepted"] > 0.99 * res["events"]
    # float32 rounding of the interpolated multiplicity means can move an event across a
    # Bernoulli threshold (measured and reported; a handful at most in 1e4 events)
    assert res["discrete_mismatches"] <= max(5, res["events"] // 1000)
    # measured float32 differences (decision 0039 float32 class, 1e-6, plus the documented
    # exceptions): energies and phi at a few 1e-7; momentum components relative to |p| at 1e-5
    # (mu carries the 2e-6 Newton step tolerance of float32 ``kalbach_mu``); the ledger columns
    # are differences of ~1e4 MeV masses whose float32 rounding is ~1e-3 MeV absolute
    worst = res["max_rel"]
    for name in ("recoil", "e_prime", "phi", "e_lab", "m_res"):
        assert worst[name] <= 1e-6, (name, worst[name])
    for name in ("px", "py", "pz"):
        assert worst[name] <= 5e-5, (name, worst[name])
    assert worst["mu"] <= 5e-3, worst["mu"]
    for name in ("imbalance", "binding", "local", "alpha_t"):
        assert worst[name] <= 0.1, (name, worst[name])


def test_event_ledger_is_summed_like_python_sum(
    host: NuclearHost, recorded: dict[str, np.ndarray]
) -> None:
    """The ledger of the shared sampler equals the ``sum()`` recomputation bit for bit (CPython
    >= 3.12 compensated float sum: the arithmetic of the reference backend since slice A)."""
    sub = {k: v[:300] for k, v in recorded.items()}
    a = _run_twin(host, sub, precision="float64")
    n = 300
    evi = a["evi"].reshape(n, EVI_STRIDE)
    evf = a["evf"].reshape(n, EVF_STRIDE)
    prod = a["prod"].reshape(n, MAX_PRODUCTS, PROD_STRIDE)
    tc = host.arrays["tconst"].reshape(host.n_targets, -1)
    checked = 0
    for i in range(n):
        if a["status"][i] != EVENT_ACCEPTED:
            continue
        c = tc[sub["tgt"][i]]
        m_p, m_t, masses = float(c[0]), float(c[1]), [float(x) for x in c[4:9]]
        rows = prod[i, : evi[i, 8]]
        e_lab = [float(x) for x in rows[:, 4]]
        big_m, recoil = float(evf[i, 5]), float(evf[i, 0])
        t1 = float(sub["t1"][i])
        assert float(evf[i, 1]) == t1 + m_p + m_t - sum(e_lab) - big_m - recoil
        m_out = sum(int(evi[i, 1 + s]) * masses[s] for s in range(5))
        assert float(evf[i, 2]) == m_out + big_m - m_p - m_t
        alpha = sum(float(r[4]) - masses[3] for r in rows if r[0] == 3.0)
        assert float(evf[i, 3]) == alpha + recoil
        checked += 1
    assert checked > 250


def test_event_sampler_status_and_overflow_contract(host: NuclearHost) -> None:
    """A recorded batch has only accepted events (status 1) and no product overflow."""
    inp = _inputs(host)
    sub = {k: v[:400] for k, v in inp.items()}
    a = _run_twin(host, sub, precision="float64")
    assert set(np.unique(a["status"])) <= {0, 1}
    assert (a["evi"].reshape(-1, EVI_STRIDE)[:, 8] <= MAX_PRODUCTS).all()


# ---------------------------------------------------------------------------------------------
# choose_target
# ---------------------------------------------------------------------------------------------
def test_choose_target_matches_the_reference_selection(host: NuclearHost) -> None:
    table = _table()
    rng = np.random.default_rng(INPUT_SEED + 1)
    n = 3000
    e = np.exp(rng.uniform(np.log(1.5), np.log(250.0), n))
    u = rng.uniform(0.0, 1.0, n)
    arr = host.arrays
    grid = arr["grid"]
    mats = (WATER, PMMA)
    for mi, mat in enumerate(mats):
        rows = table.material_rows(mat)
        kc = int(arr["mat_ntargets"][mi])
        tg = arr["mat_target"].reshape(host.n_materials, host.kmax)[mi]
        for i in range(n):
            cum = rows.cum_fraction_at(float(e[i]))
            chosen = -1
            for k in range(kc):
                chosen = int(
                    NU64.select_target(float(u[i]), float(cum[k]), k, chosen, int(k == kc - 1))
                )
            got = NU64.choose_target(
                float(u[i]), float(e[i]), grid, host.n_grid, arr["sigma"], arr["cum_sigma"],
                mi, host.kmax, kc,
            )  # fmt: skip
            assert got == chosen
            if kc:
                assert int(tg[got]) == rows.target_index[chosen]
    assert int(arr["mat_ntargets"][0]) >= 1
