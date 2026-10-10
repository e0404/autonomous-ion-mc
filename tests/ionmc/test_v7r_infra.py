"""CI-scale tests of the V7-R re-test infrastructure (V3-005C C5): the ``batch_estimates`` column
fix, and (with ``validation/scripts/transport/v7r.py`` and ``steps_v5c.py``) the block-sum
sidecar, the bootstrap-t gates and the suite registration. No transport beyond a few hundred
histories."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"


def _load(name: str) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def v5b() -> ModuleType:
    return _load("steps_v5b")


def _cfg(v5b: ModuleType, elastic: bool, n: int = 600, nb: int = 3):  # type: ignore[no-untyped-def]
    geo, grid = v5b.coarse_depth(v5b.V7_BINS)
    cfg = v5b.v7_config("warp-cpu", "float64", n, nb, 20505001, grid=grid, geo=geo)
    return v5b.replace(cfg, physics=v5b.replace(cfg.physics, elastic=elastic))


def _old_batch_estimates(v5b: ModuleType, cfg):  # type: ignore[no-untyped-def]
    """The pre-fix column addressing (the LAST six columns), valid only with the elastic block
    off."""
    from ionmc.transport.run import run_range
    from ionmc.transport.tally import NUCLEAR_TALLY_NAMES

    eff = v5b.block_effective(cfg, None)
    hpb = cfg.run.n_histories // cfg.run.n_batches
    k0 = len(NUCLEAR_TALLY_NAMES)
    ix = {name: i for i, name in enumerate(NUCLEAR_TALLY_NAMES)}
    loc, esc = [], []
    for b in range(cfg.run.n_batches):
        part = run_range(eff, b * hpb, (b + 1) * hpb)
        comps = part.tally_components
        col = {name: math.fsum(comps[len(comps) - k0 + i]) for name, i in ix.items()}
        loc.append(col["nuclear_local"] / hpb)
        esc.append((col["nuclear_escaped_neutron"] + col["nuclear_escaped_gamma"]) / hpb)
    return np.array(loc), np.array(esc)


def test_batch_estimates_elastic_off_is_bit_identical_to_the_last_six_columns(
    v5b: ModuleType,
) -> None:
    """Frozen lv5b/hr5 path (elastic off): the by-name addressing equals the old tail addressing."""
    cfg = _cfg(v5b, elastic=False)
    new = v5b.batch_estimates(cfg)
    loc, esc = _old_batch_estimates(v5b, cfg)
    assert np.array_equal(new["nuclear_local"], loc)
    assert np.array_equal(new["escaped_neutral"], esc)
    # existing fields unchanged by the additive ``block_sums`` output; the sum is the per-primary
    # value times the block histories to 1e-14 relative
    hpb = cfg.run.n_histories // cfg.run.n_batches
    assert set(new) >= {
        "sec_p",
        "nuc_local_dose",
        "idd",
        "nuclear_local",
        "escaped_neutral",
        "block_sums",
        "counters_sum",
        "batch_assignment_ok",
    }
    np.testing.assert_allclose(new["block_sums"] / hpb, new["escaped_neutral"], rtol=1e-14, atol=0)


def test_batch_estimates_elastic_on_reads_the_nuclear_block_not_the_elastic_tail(
    v5b: ModuleType,
) -> None:
    """With elastic on the three elastic columns follow the nuclear block: the per-block values
    equal the ``merge_partials`` tallies of the same block (1e-12 relative), and differ from the
    old tail addressing (the defect)."""
    from ionmc.transport.run import channel_columns, run_range
    from ionmc.transport.tally import merge_partials

    cfg = _cfg(v5b, elastic=True)
    eff = v5b.block_effective(cfg, None)
    assert eff.nuclear is not None and eff.nuclear.elastic is not None
    new = v5b.batch_estimates(cfg)
    hpb = cfg.run.n_histories // cfg.run.n_batches
    n_grids = len(eff.requested.scoring)
    for b in range(cfg.run.n_batches):
        part = run_range(eff, b * hpb, (b + 1) * hpb)
        # a block is a complete partition of [0, hpb) only after a shift of its range: reuse the
        # partial with h0 = 0 for the reduction of the exact column sums
        part.h1 -= part.h0
        part.h0 = 0
        raw = merge_partials([part], hpb, n_grids, channel_columns(eff), True, True)
        t = raw.tallies
        want_loc = t["nuclear_local"] / hpb
        want_esc = (t["nuclear_escaped_neutron"] + t["nuclear_escaped_gamma"]) / hpb
        assert new["nuclear_local"][b] == pytest.approx(want_loc, rel=1e-12, abs=0.0)
        assert new["escaped_neutral"][b] == pytest.approx(want_esc, rel=1e-12, abs=0.0)
    loc, esc = _old_batch_estimates(v5b, cfg)
    assert not np.allclose(esc, new["escaped_neutral"], rtol=1e-6)  # the defect was real


def test_tally_column_index_matches_the_merge_partials_layout() -> None:
    from ionmc.transport.tally import ELASTIC_TALLY_NAMES, NUCLEAR_TALLY_NAMES, tally_column_index

    n = 20
    for name in NUCLEAR_TALLY_NAMES:
        i = NUCLEAR_TALLY_NAMES.index(name)
        assert tally_column_index(n, False, name) == n - 6 + i
        assert tally_column_index(n, True, name) == n - 3 - 6 + i
    for i, name in enumerate(ELASTIC_TALLY_NAMES):
        assert tally_column_index(n, True, name) == n - 3 + i
    with pytest.raises(ValueError):
        tally_column_index(n, False, "elastic_recoil_local")
    with pytest.raises(ValueError):
        tally_column_index(n, True, "initial")


# -- v7r.py: constants, streams, bootstrap-t, gates -------------------------------------------
@pytest.fixture(scope="module")
def v7r(v5b: ModuleType) -> ModuleType:  # v5b first: v7r imports steps_v5b
    return _load("v7r")


def test_v7r_constants_seeds_ranks_pairing_and_heldout_split(
    v5b: ModuleType, v7r: ModuleType
) -> None:
    assert (v7r.V7R_SHARDS, v7r.V7R_REPS_PER_SHARD, v7r.V7R_REPLICATES, v7r.V7R_PAIRS) == (
        16,
        450,
        7200,
        3600,
    )
    assert v7r.V7R_SHARD_N == v7r.V7R_REPS_PER_SHARD * v7r.V7R_BLOCKS_PER_REP * v7r.V7R_BATCH_N
    assert v7r.V7R_SEED_FIRST == 20497004 and v7r.V7R_REF_SEED == 20497020
    assert [v7r.shard_seed(k) for k in (0, 15, 16)] == [20497004, 20497019, 20497020]
    assert v7r.V7R_DIAG_SEED_FIRST == 20498004
    assert (v7r.V7R_RANK_LO, v7r.V7R_RANK_HI) == (317, 1683)
    assert (v7r.V7R_B, v7r.V7R_ALPHA, v7r.V7R_REGION) == (1999, 0.3173, (0.6527, 0.7127))
    # pairing: replicate j of shard s with replicate j of shard s + 8, i.e. rows i and i + 3600
    ia, ib = v5b.v7_pair_indices(v7r.V7R_REPLICATES)
    for s in range(8):
        for j in (0, 449):
            p = s * 450 + j
            assert (ia[p], ib[p]) == (s * 450 + j, (s + 8) * 450 + j)
    # held-out H = shards 12-15, evaluation E = shards 0-11
    n_e = v5b.v7_heldout_split(v7r.V7R_REPLICATES)
    assert n_e == 12 * 450 and v7r.V7R_REPLICATES - n_e == 4 * 450
    assert v7r.V7R_HELDOUT_SHARDS == (12, 13, 14, 15)


def test_v7r_streams_are_the_literal_seed_sequences_and_deterministic(v7r: ModuleType) -> None:
    a = v7r.single_stream(20481004, 3, 7).integers(0, 20, size=(4, 20))
    b = np.random.default_rng(np.random.SeedSequence([20481004, 16, 3, 7])).integers(
        0, 20, size=(4, 20)
    )
    assert np.array_equal(a, b)
    c = v7r.pair_stream(20481004, 2, 9).integers(0, 20, size=5)
    d = np.random.default_rng(np.random.SeedSequence([20481004, 16, 1002, 9])).integers(
        0, 20, size=5
    )
    assert np.array_equal(c, d)
    # a rehearsal base gives different streams
    e = v7r.single_stream(20505000, 3, 7).integers(0, 20, size=(4, 20))
    assert not np.array_equal(a, e)
    rec = v7r.rng_record()
    assert rec["bit_generator"] == "PCG64" and rec["numpy"] == np.__version__


def _gauss_blocks(
    rows: int, rng: np.random.Generator, mu: float = 100.0, sd: float = 20.0
) -> np.ndarray:
    return rng.normal(mu, sd, size=(rows, 20))


def test_boot_t_interval_is_deterministic_and_about_the_t_interval_on_gaussian_blocks(
    v7r: ModuleType,
) -> None:
    rng = np.random.default_rng(5)
    x = _gauss_blocks(1, rng)[0]
    i1 = v7r.boot_t_interval(x, v7r.single_stream(1, 0, 0))
    i2 = v7r.boot_t_interval(x, v7r.single_stream(1, 0, 0))
    assert i1[:4] == i2[:4]
    lo, hi, m, se, diag = i1
    assert m == pytest.approx(x.mean()) and se == pytest.approx(x.std(ddof=1) / math.sqrt(20))
    assert lo < m < hi and diag["n_resamples_degenerate"] == 0
    # the coverage of the TRUE mean over 400 synthetic replicates is the nominal 0.6827 within
    # 4 sigma
    hits = 0
    for r in range(400):
        xr = _gauss_blocks(1, rng)[0]
        lo, hi, *_ = v7r.boot_t_interval(xr, v7r.single_stream(2, 0, r))
        hits += int(lo <= 100.0 <= hi)
    cov = hits / 400
    assert abs(cov - 0.6827) < 4.0 * math.sqrt(0.6827 * 0.3173 / 400), cov


def test_boot_t_pair_interval_gaussian_coverage_and_draw_order(v7r: ModuleType) -> None:
    rng = np.random.default_rng(6)
    hits = 0
    for r in range(400):
        xj, xk = _gauss_blocks(2, rng)
        lo, hi, d, se, _ = v7r.boot_t_pair_interval(xj, xk, v7r.pair_stream(3, 0, r))
        hits += int(lo <= 0.0 <= hi)
    assert abs(hits / 400 - 0.6827) < 4.0 * math.sqrt(0.6827 * 0.3173 / 400)
    # j's indices are drawn first, then k's, from the one generator
    xj, xk = _gauss_blocks(2, rng)
    g = v7r.pair_stream(3, 1, 1)
    ij, ik = g.integers(0, 20, size=(1999, 20)), g.integers(0, 20, size=(1999, 20))
    mj, mk = xj[ij], xk[ik]
    t = ((mj.mean(1) - mk.mean(1)) - (xj.mean() - xk.mean())) / np.sqrt(
        mj.std(1, ddof=1) ** 2 / 20 + mk.std(1, ddof=1) ** 2 / 20
    )
    t.sort()
    d_hat = xj.mean() - xk.mean()
    se_d = math.sqrt(xj.var(ddof=1) / 20 + xk.var(ddof=1) / 20)
    lo, hi, *_ = v7r.boot_t_pair_interval(xj, xk, v7r.pair_stream(3, 1, 1))
    assert lo == pytest.approx(d_hat - t[1682] * se_d, rel=1e-12)
    assert hi == pytest.approx(d_hat - t[316] * se_d, rel=1e-12)


def test_degenerate_se_is_a_miss_and_se_star_zero_follows_the_amendment_rule(
    v7r: ModuleType,
) -> None:
    flat = np.full(20, 7.0)
    lo, hi, m, se, diag = v7r.boot_t_interval(flat, v7r.single_stream(1, 0, 0))
    assert se == 0.0 and math.isnan(lo) and math.isnan(hi) and diag["degenerate"]
    lo, hi, d, se, diag = v7r.boot_t_pair_interval(flat, np.arange(20.0), v7r.pair_stream(1, 0, 0))
    assert math.isnan(lo) and diag["degenerate"]  # miss: no interval can contain 0
    # t* of resamples with se* = 0: +inf above, -inf below, 0 on m_hat
    xs = np.array([[5.0] * 4, [1.0] * 4, [3.0] * 4])
    t, n = v7r._tstar(xs, 3.0)
    assert n == 3 and t[0] == math.inf and t[1] == -math.inf and t[2] == 0.0
    # one odd block among 19 equal: 36 % of the resamples are all-equal; counted, interval finite
    x = np.full(20, 10.0)
    x[0] = 30.0
    lo, hi, m, se, diag = v7r.boot_t_interval(x, v7r.single_stream(4, 0, 0))
    assert 200 < diag["n_resamples_degenerate"] < 1000 and se > 0.0
    # gates: a degenerate replicate or pair is counted as a miss
    sums = np.vstack([_gauss_blocks(7196, np.random.default_rng(1)), np.tile(flat, (4, 1))])
    pi = v7r.pair_ingredients(sums, 11, reps_per_shard=450)
    assert (pi["n_degenerate"] < 0).sum() == 4
    assert v7r.pair_gate(pi)["degenerate_pairs"] == 4


def test_kappa_monotone_and_gates_on_synthetic_gaussian_blocks(v7r: ModuleType) -> None:
    rng = np.random.default_rng(7)
    sums = _gauss_blocks(1200, rng)  # 1200 replicates (150 per "shard"): pairs 600, E 900, H 300
    pi = v7r.pair_ingredients(sums, 21, reps_per_shard=75)
    si = v7r.single_ingredients(sums, 21, reps_per_shard=75)
    ms = [v7r.pair_gate(pi, k)["m"] for k in (0.8, 0.9, 1.0, 1.1, 1.2)]
    assert ms == sorted(ms) and ms[0] < ms[-1]
    assert 0.55 < ms[2] < 0.80
    means = sums.mean(axis=1)
    cs = [
        v7r.single_gate(si, means, k, boot_seed=1, n_boot=200)["cp_count_min"]
        for k in (0.8, 1.0, 1.2)
    ]
    assert cs == sorted(cs) and cs[0] < cs[-1]
    g = v7r.single_gate(si, means, boot_seed=1, n_boot=200)
    assert g["replicates_eval"] == 900 and g["replicates_heldout"] == 300 and g["z_box"] >= g["z_t"]
    # determinism of the whole verdict at CI scale
    v1 = v7r.escaped_neutral_verdict(sums, 21, boot_seed=1, reps_per_shard=75, n_boot=200)
    v2 = v7r.escaped_neutral_verdict(sums, 21, boot_seed=1, reps_per_shard=75, n_boot=200)
    assert v1 == v2 and v1["interval"] == "bootstrap-t" and v1["ranks"] == [317, 1683]
    assert v1["n_pairs"] == 600


def test_aggregate_blocks_is_exact_grouping(v7r: ModuleType) -> None:
    rng = np.random.default_rng(8)
    s = rng.gamma(2.0, 3.0, size=(5, 20))
    assert np.array_equal(v7r.aggregate_blocks(s, 20), s)
    a10, a5 = v7r.aggregate_blocks(s, 10), v7r.aggregate_blocks(s, 5)
    assert a10.shape == (5, 10) and a5.shape == (5, 5)
    assert a10[2, 3] == math.fsum(s[2, 6:8]) and a5[1, 4] == math.fsum(s[1, 16:20])
    for a in (a10, a5):
        assert np.allclose(a.sum(axis=1), s.sum(axis=1), rtol=1e-14)
    with pytest.raises(ValueError):
        v7r.aggregate_blocks(s, 40)
    assert math.isnan(v7r.excess_kurtosis(np.ones(5)))


# -- steps_v5c.py: configuration builders, sidecar, combine -------------------------------------
@pytest.fixture(scope="module")
def v5c(v5b: ModuleType, v7r: ModuleType) -> ModuleType:
    return _load("steps_v5c")


def test_elastic_table_id_equals_the_pin_and_builders_set_elastic(
    v5b: ModuleType, v5c: ModuleType
) -> None:
    import json

    pin = json.loads((REPO / "src" / "ionmc" / "data" / "elastic_table_pin.json").read_text())
    assert v5c.ELASTIC_TABLE_ID == pin["table_id"]
    geo, grid = v5b.coarse_depth(v5b.V7_BINS)
    kw = dict(energy=150.0, n=100, seed=1, geometry=geo, grid=grid)
    frozen = v5b.wcfg("warp-cpu", "float64", **kw)
    assert frozen.physics.elastic is False  # frozen lv5/lv5b/hr5 path unchanged
    c = v5c.nuc_config_c(**kw)
    assert c.physics.elastic is True and c.physics.elastic_only is False
    assert c.physics.elastic_table_id == pin["table_id"]
    eo = v5c.nuc_config_c_elastic_only(**kw)
    assert eo.physics.elastic_only is True and eo.physics.elastic
    w = v5c.wcfg_c("warp-cpu", "float64", **kw)
    assert (w.run.backend, w.run.precision, w.physics.elastic) == ("warp-cpu", "float64", True)
    cfg = v5c.v7r_config(10_000, 20, 5)
    assert (
        cfg.physics.elastic
        and cfg.run.backend == "warp-cpu"
        and cfg.source.kinetic_energy_mev == 150.0
    )
    from ionmc.simulation import Simulation

    summary = Simulation(c).effective.summary()
    assert "elastic" in summary and summary["elastic"]["table_id"] == pin["table_id"]


def _doc_with_sidecar(
    v5c: ModuleType,
    v7r: ModuleType,
    tmp_path: Path,
    k: int,
    sums: np.ndarray,
    seed_base: int = 20505000,
    name: str | None = None,
) -> tuple[dict, Path]:  # type: ignore[type-arg]
    import json

    side = v5c.write_v7r_sidecar(tmp_path / f"v7r-s{k}-escaped_neutral.npz", sums)
    per = sums / 500.0
    doc = {
        "shard": k,
        "seed": seed_base + 16000 + k,
        "seed_base": seed_base,
        "valid": True,
        "n": sums.shape[0] * 10_000,
        "reduced": True,
        "replicates": sums.shape[0],
        "batch_histories": 500,
        "sidecar": side,
        "estimators": {
            "escaped_neutral": {
                "mean": per.mean(axis=1)[:, None].tolist(),
                "sem": (per.std(axis=1, ddof=1) / math.sqrt(20))[:, None].tolist(),
            }
        },
    }
    p = tmp_path / f"v7r-s{k}.json"
    p.write_text(json.dumps(doc))
    return doc, p


def test_sidecar_roundtrip_determinism_and_fail_closed(
    v5c: ModuleType, v7r: ModuleType, tmp_path: Path
) -> None:
    import json

    sums = np.random.default_rng(1).gamma(2.0, 50.0, size=(6, 20))
    doc, p = _doc_with_sidecar(v5c, v7r, tmp_path, 0, sums)
    rec = doc["sidecar"]
    assert (
        rec["shape"] == [6, 20]
        and rec["dtype"] == "float64"
        and rec["units"] == "MeV per 500-history block"
    )
    assert rec["array"] == "escaped_neutral_block_sums" and rec["batch_histories"] == 500
    assert rec["array_sha256"] == hashlib_sha(np.ascontiguousarray(sums, dtype="<f8").tobytes())
    assert np.array_equal(v5c.load_v7r_sidecar(doc, p), sums)
    # byte-determinism: writing twice gives the same file bytes
    again = v5c.write_v7r_sidecar(tmp_path / "again.npz", sums)
    assert again["sha256"] == rec["sha256"]
    # a tampered npz fails closed
    raw = bytearray((tmp_path / rec["file"]).read_bytes())
    raw[len(raw) // 2] ^= 0xFF
    (tmp_path / rec["file"]).write_bytes(bytes(raw))
    with pytest.raises(SystemExit):
        v5c.load_v7r_sidecar(doc, p)
    # a missing file, a wrong shape, a wrong dtype, wrong units, a wrong JSON mean
    doc, p = _doc_with_sidecar(v5c, v7r, tmp_path, 1, sums)
    bad = json.loads(json.dumps(doc))
    bad["sidecar"]["shape"] = [5, 20]
    with pytest.raises(SystemExit):
        v5c.load_v7r_sidecar(bad, p)
    bad = json.loads(json.dumps(doc))
    bad["sidecar"]["units"] = "MeV per primary"
    with pytest.raises(SystemExit):
        v5c.load_v7r_sidecar(bad, p)
    bad = json.loads(json.dumps(doc))
    bad["estimators"]["escaped_neutral"]["mean"][0][0] *= 1.0 + 1e-9
    with pytest.raises(SystemExit):
        v5c.load_v7r_sidecar(bad, p)
    bad = json.loads(json.dumps(doc))
    bad["sidecar"]["file"] = "missing.npz"
    with pytest.raises(SystemExit):
        v5c.load_v7r_sidecar(bad, p)
    with pytest.raises(SystemExit):  # wrong dtype written by hand (float32 array, matching hashes)
        arr32 = tmp_path / "f32.npz"
        np.savez(arr32, escaped_neutral_block_sums=sums.astype(np.float32))
        b = json.loads(json.dumps(doc))
        b["sidecar"].update(file="f32.npz", sha256=hashlib_sha(arr32.read_bytes()))
        v5c.load_v7r_sidecar(b, p)
    with pytest.raises(SystemExit):  # non-finite
        v5c.write_v7r_sidecar(tmp_path / "nan.npz", np.full((2, 20), np.nan))


def hashlib_sha(b: bytes) -> str:
    import hashlib

    return hashlib.sha256(b).hexdigest()


def _synthetic_v7r(
    v5c: ModuleType,
    v7r: ModuleType,
    tmp_path: Path,
    n_shards: int = 4,
    rps: int = 75,
    seed_base: int = 20505000,
):  # type: ignore[no-untyped-def]
    rng = np.random.default_rng(3)
    shards = []
    allsums = rng.gamma(6.0, 20.0, size=(n_shards * rps, 20))
    for k in range(n_shards):
        shards.append(
            _doc_with_sidecar(v5c, v7r, tmp_path, k, allsums[k * rps : (k + 1) * rps], seed_base)
        )

    def est(rows: int) -> dict:  # type: ignore[type-arg]
        return {"mean": [[1.0] * 12], "sem": [[0.1] * 12]}

    # shard documents need the profile estimators too (12 bins, replicates x 12)
    for doc, _ in shards:
        r = doc["replicates"]
        g = rng.gamma(20.0, 1.0, size=(r, 12))
        doc["estimators"]["sec_p"] = {
            "mean": g.tolist(),
            "sem": (0.2 * g.mean(axis=0) + 0 * g).tolist(),
        }
        doc["estimators"]["nuclear_local"] = {
            "mean": g.tolist(),
            "sem": (0.2 * g.mean(axis=0) + 0 * g).tolist(),
        }
    blocks = rng.gamma(20.0, 1.0, size=(20, 12))
    refest = {
        nm: {
            "mean": blocks.mean(axis=0).tolist(),
            "sem": (blocks.std(axis=0, ddof=1) / math.sqrt(20)).tolist(),
            "blocks": blocks.tolist(),
        }
        for nm in ("sec_p", "nuclear_local")
    }
    eb = rng.gamma(6.0, 20.0, size=(20, 1))
    refest["escaped_neutral"] = {
        "mean": [float(eb.mean())],
        "sem": [float(eb.std(ddof=1) / math.sqrt(20))],
        "blocks": eb.tolist(),
    }
    ref = {"seed": seed_base + 16000 + 16, "valid": True, "n": 100_000, "estimators": refest}
    # the escaped_neutral estimator of the shards is 1-bin (scalar): [r, 1] already
    return shards, ref, allsums


def test_v7r_combine_dry_run_document_and_fail_closed(
    v5c: ModuleType, v7r: ModuleType, tmp_path: Path
) -> None:
    shards, ref, allsums = _synthetic_v7r(v5c, v7r, tmp_path)
    doc = v5c.v7r_combine(shards, ref, 20505000, n_shards=4, n_boot=100)
    en = doc["estimators"]["escaped_neutral"]
    assert en["interval"] == "bootstrap-t" and en["B"] == 1999 and en["alpha"] == 0.3173
    assert "amendment13_rule_report_only" in en and "pass_paired" in en and "pass_single" in en
    assert en["pass"] == (
        en["pass_paired"] and en["pass_single"]
    )  # the report-only rule is not in pass
    assert doc["estimators"]["sec_p"]["kind"] == "profile" and doc["replicates"] == 300
    assert doc["seeds"]["shards"] == [20521000 + k for k in range(4)] and doc["seeds"][
        "reference"
    ] == [20521016]
    assert doc["rng"]["bit_generator"] == "PCG64" and doc["units"] == "MeV per 500-history block"
    assert doc["pass"] == doc["gates_pass"]
    # fail closed: seed mismatch, invalid shard, unequal replicates, wrong shard count
    with pytest.raises(SystemExit):
        v5c.v7r_combine(shards, {**ref, "seed": 1}, 20505000, n_shards=4)
    with pytest.raises(SystemExit):
        v5c.v7r_combine(shards[:3], ref, 20505000, n_shards=4)
    bad = [dict(shards[0][0], valid=False), *[s[0] for s in shards[1:]]]
    with pytest.raises(SystemExit):
        v5c.v7r_combine(
            [(d, p) for d, (_, p) in zip(bad, shards, strict=True)], ref, 20505000, n_shards=4
        )
    # a seed base whose lists intersect a consumed evidence seed set is refused
    with pytest.raises(SystemExit):
        shards2, ref2, _ = _synthetic_v7r(v5c, v7r, tmp_path, seed_base=20471004 - 16000 + 5000)
        v5c.v7r_combine(shards2, ref2, 20471004 - 16000 + 5000, n_shards=4)
    # consumed sets contain the lv5b V7-replicate seeds 20482004..20482012
    cs = v5c.consumed_seed_set()
    assert set(range(20482004, 20482013)) <= cs and not (cs & set(range(20497004, 20497021)))


def test_invalid_counters_excludes_exactly_the_elastic_diagnostics_by_name(v5c: ModuleType) -> None:
    from ionmc.simulation import ElasticTransportCounters
    from ionmc.transport.tally import ELASTIC_COUNTER_NAMES

    names = v5c.counter_names()
    assert set(ElasticTransportCounters.DIAGNOSTIC) == set(ELASTIC_COUNTER_NAMES)
    zero = {"counter_vector": np.zeros(len(names), dtype=np.int64)}
    assert v5c.invalid_counters(zero) == 0
    for nm in names:
        vec = np.zeros(len(names), dtype=np.int64)
        vec[names.index(nm)] = 7
        est = {"counter_vector": vec}
        assert v5c.counters_by_name(est)[nm] == 7
        want = 0 if nm in ELASTIC_COUNTER_NAMES else 7
        assert v5c.invalid_counters(est) == want, nm
    for nm in ("majorant_violation", "queue_overflow"):
        vec = np.zeros(len(names), dtype=np.int64)
        vec[names.index(nm)] = 1
        assert v5c.invalid_counters({"counter_vector": vec}) == 1
    with pytest.raises(SystemExit):
        v5c.counters_by_name({"counter_vector": np.zeros(3, dtype=np.int64)})


# -- suite registration (lv5c, V7-R steps only) --------------------------------------------------
def test_lv5c_suite_registration_names_seeds_timeouts_tags(
    v7r: ModuleType, v5c: ModuleType
) -> None:
    rs, sm = _load("run_suite"), _load("summarize")
    assert rs.V7R_SHARDS == v5c.V7R_SHARDS == v7r.V7R_SHARDS == 16
    assert rs.DEFAULT_SEED_BASES["lv5c"] == 20481004 == sm.V5C_QUALIFICATION_SEED_BASE
    assert "lv5c" in rs.SUITES and rs.DEFERRED_STEPS["lv5c"] == ()
    names = [n.split("-", 1)[1] for n in rs.full_step_names("lv5c", 2)]
    assert names == [*(f"v7r-s{k}" for k in range(16)), "v7r-ref", "pytest-v7r-calibration",
                     "v7r-combine", "v7r-diag", *v5c.V11_STEP_NAMES, "v11-compare"]  # fmt: skip
    steps = rs.suite_steps("lv5c", 1, 1.0)
    for name, cmd, env in steps:
        assert env["IONMC_REQUIRE_DATA"] == "1"
        if "pytest" in name:
            assert "tests/ionmc/test_v7r_coverage.py" in cmd and "--seed-base" not in cmd
            assert cmd[cmd.index("-m", 3) + 1] == "calibration"
            continue
        assert "steps_v5c.py" in " ".join(cmd) and cmd[cmd.index("--seed-base") + 1] == "20481004"
        assert "--timeout" in cmd
    by = {n.split("-", 1)[1]: c for n, c, _ in steps}
    for k in (0, 15):
        c = by[f"v7r-s{k}"]
        assert c[c.index("v7r-shard") + 1 :][:2] == ["--shard", str(k)] and "--out-dir" in c
        # the seed of shard k derives from the base: 20481004 + 1000 * 16 + k
        assert 20481004 + 1000 * v7r.V7R_R_INDEX + k == 20497004 + k
    assert "--dirs" in by["v7r-combine"] and "--dirs" in by["v7r-diag"]
    full = rs.full_step_names("lv5c", 2)
    for n in full:
        base_name = n.split("-", 1)[1]
        want = 1800 if base_name in ("v7r-combine", "v7r-diag", "v11-compare") else 3300
        assert rs.step_timeout_s("lv5c", n, 1500) == want, n
    tags = {nm: sm.expected_tag(nm) for nm in full}
    assert tags["01-v7r-s0"] == tags["16-v7r-s15"] == "v7r-shard"
    assert tags["17-v7r-ref"] == "v7r-ref" and tags["18-pytest-v7r-calibration"] is None
    assert tags["19-v7r-combine"] == "v7r-combine" and tags["20-v7r-diag"] == "v7r-diag"
    got = {str(f.relative_to(rs.REPO)) for f in rs.source_files("lv5c")}
    assert {"validation/scripts/transport/steps_v5c.py", "validation/scripts/transport/v7r.py",
            "validation/scripts/transport/steps_v5b.py"} <= got  # fmt: skip
    assert "validation/scripts/transport/steps_v5c.py" not in rs.source_file_list("lv5b")
    assert sm.seed_blockers(20481004, "lv5c") == []
    assert sm.seed_blockers(20505000, "lv5c") and sm.seed_blockers(20471004, "lv5c")
    with pytest.raises(SystemExit):
        rs.main(
            ["--suite", "lv5c", "--out", "x", "--expected-sha", "0" * 40, "--seed-base", "20471004"]
        )


def _args(shard: int | None = None) -> str:
    rs = _load("run_suite")
    targets, k_expr, _ = rs.calibration_selection(shard)
    return "python -m pytest " + " ".join(rs.calibration_args(targets, k_expr))


def _cal(  # type: ignore[no-untyped-def]
    ok: bool = True, counts: dict | None = None, cmd: str | None = None, shard: int | None = None
) -> dict:  # type: ignore[type-arg]
    rs = _load("run_suite")
    want = rs.calibration_selection(shard)[2]
    counts = {"passed": want, "deselected": 37} if counts is None else counts
    return {"pass": ok, "status": "executed", "pytest_counts": counts,
            "pytest_command": _args(shard) if cmd is None else cmd}  # fmt: skip


CAL_CMD = _args(None)


def test_v7r_row_is_the_conjunction_of_calibration_and_gates() -> None:
    sm = _load("summarize")

    def st(ok: bool) -> dict:  # type: ignore[type-arg]
        return {"pass": ok, "status": "executed"}

    full = {
        "17-v7r-ref": st(True),
        "18-pytest-v7r-calibration": _cal(),
        "19-v7r-combine": st(True),
    }
    assert sm.slice_c_rows("lv5c", full)["V7-R"]["pass"] is True
    for key, bad_step in (
        ("18-pytest-v7r-calibration", _cal(False)),
        ("19-v7r-combine", st(False)),
    ):
        row = sm.slice_c_rows("lv5c", {**full, key: bad_step})["V7-R"]
        assert row["pass"] is False
        assert row["verdict"] == "uncertainty-coverage evidence not established"
    no_cal = {k: v for k, v in full.items() if "pytest" not in k}
    assert (
        sm.slice_c_rows("lv5c", no_cal)["V7-R"]["pass"] is False
    )  # gates alone never pass the row
    shards = {
        **full,
        "18-pytest-v7r-calibration-s0": _cal(),
        "19-pytest-v7r-calibration-s1": _cal(False),
    }
    assert (
        sm.slice_c_rows("lv5c", shards)["V7-R"]["pass"] is False
    )  # every calibration shard must pass
    assert sm.slice_c_rows("lv5b", full) == {}


ONLY_V7R = CAL_CMD.replace(" tests/ionmc/test_v7_coverage.py", "")
SHARD0 = _args(0)


@pytest.mark.parametrize(
    "step",
    [
        _cal(counts={"passed": 12, "skipped": 1, "deselected": 37}),  # a skip exits 0, no pass
        _cal(counts={"passed": 13, "deselected": 37, "xfailed": 1}),
        _cal(counts={"passed": 12, "failed": 1}),
        _cal(counts={"passed": 3, "deselected": 37}),  # the old 3-case count of one file
        _cal(counts={"passed": 12}),  # count mismatch: not the complete expected case set
        _cal(counts={"passed": 14}),
        _cal(counts={"skipped": 13}),  # nothing passed at all
        _cal(counts={}),
        {**_cal(), "pytest_counts": None},  # unreadable output: not established
        {k: v for k, v in _cal().items() if k != "pytest_counts"},
        _cal(counts={"passed": 3}, cmd=ONLY_V7R),  # only the V7-R target: no Amendment 13 file
        _cal(cmd=ONLY_V7R),  # right count for the full set but the wrong command
        _cal(cmd=CAL_CMD + " -k archive"),  # a narrowed selection deselects calibration tests
        _cal(cmd=CAL_CMD + " --deselect tests/ionmc/test_v7r_coverage.py::test_x"),
        _cal(cmd=CAL_CMD.replace("-m calibration", "-m 'not calibration'")),
        _cal(cmd=CAL_CMD.replace(" -m calibration", "")),
        _cal(cmd=CAL_CMD + " tests/ionmc/test_other.py"),
        _cal(shard=0, cmd=SHARD0.replace("archive_fitted", "archive")),  # not the fixed -k
    ],
)
def test_v7r_calibration_row_rejects_skips_deselection_and_unreadable_output(
    step: dict,  # type: ignore[type-arg]
) -> None:
    sm = _load("summarize")
    steps = {"17-v7r-ref": _cal(), "19-v7r-combine": {"pass": True, "status": "executed"},
             "18-pytest-v7r-calibration": step}  # fmt: skip
    row = sm.slice_c_rows("lv5c", steps)["V7-R"]
    assert row["gates_pass"] is True and row["calibration_pass"] is False and row["pass"] is False
    assert row["verdict"] == "uncertainty-coverage evidence not established"


def test_pytest_summary_parser_and_parse_step_record_the_counts(tmp_path: Path) -> None:
    sm = _load("summarize")
    pc = sm.pytest_counts
    assert pc("...\n4 passed, 37 deselected, 1 warning in 12.30s (0:00:12)") == {
        "passed": 4, "deselected": 37, "warning": 1}  # fmt: skip
    assert pc("=== 3 passed, 2 skipped in 1.0s ===") == {"passed": 3, "skipped": 2}
    assert pc("1 failed, 2 passed, 1 error in 3.00s") == {"failed": 1, "passed": 2, "error": 1}
    assert pc("no summary here") is None and pc("") is None
    sha = "a" * 40
    p = tmp_path / "18-pytest-v7r-calibration.txt"
    body = ["# command: " + CAL_CMD, f"# git_sha: {sha}", "# started_utc: t", "# step_timeout_s: 1",
            "", "....", "13 passed, 37 deselected in 12.30s", "", "# exit=0"]  # fmt: skip
    p.write_text("\n".join(body))
    st = sm.parse_step(p, sha, "lv5c", "20481004")
    assert st["pass"] and st["pytest_counts"] == {"passed": 13, "deselected": 37}
    assert st["pytest_command"] == CAL_CMD
    assert sm.calibration_step_established(st)
    p.write_text(
        "\n".join(body).replace("13 passed, 37 deselected", "12 passed, 1 skipped, 37 deselected")
    )
    st = sm.parse_step(p, sha, "lv5c", "20481004")
    assert st["pass"] is True and not sm.calibration_step_established(st)  # exit 0 with a skip
    p.write_text("\n".join(body).replace("13 passed, 37 deselected in 12.30s", "garbled"))
    st = sm.parse_step(p, sha, "lv5c", "20481004")
    assert st["pytest_counts"] is None and not sm.calibration_step_established(st)


def test_lv5c_calibration_steps_run_with_the_fixtures_required() -> None:
    rs = _load("run_suite")
    steps = rs.suite_steps("lv5c", 1, 1.0)
    cal = [(n, e) for n, _, e in steps if "pytest-v7r-calibration" in n]
    assert cal and all(e["IONMC_V7R_FIXTURES"] == "required" for _, e in cal)
    assert all("IONMC_V7R_REHEARSAL_DIR" not in e for _, e in cal)
    assert all("IONMC_V7R_FIXTURES" not in e for n, _, e in steps if "pytest" not in n)
    shared = rs.NUCLEAR_ENV
    assert "IONMC_V7R_FIXTURES" not in shared  # the shared dict is not mutated


def test_sub_block_run_range_250_plus_250_equals_the_500_history_block_bit_for_bit(
    v5c: ModuleType, v5b: ModuleType
) -> None:
    """Exploratory sub-blocks (b = 40, 80) need no change of the RNG keying: the counter-based RNG
    is per history, and the exact component expansions of the sub-ranges, summed with
    ``math.fsum``, give the correctly rounded total of the block: every tally column (elastic on,
    warp-cpu float64), the escaped-neutral block sum and the integer grids / counters are
    identical."""
    from ionmc.transport.run import run_range
    from ionmc.transport.tally import NUCLEAR_TALLY_NAMES, tally_column_index

    cfg = v5c.v7r_config(1000, 2, 20505003)
    eff = v5b.block_effective(cfg, None)
    whole = run_range(eff, 500, 1000)
    parts = [run_range(eff, 500, 750), run_range(eff, 750, 1000)]
    assert len(whole.tally_components) == len(parts[0].tally_components)
    for col, comp in enumerate(whole.tally_components):
        joined = math.fsum(c for p in parts for c in p.tally_components[col])
        assert math.fsum(comp) == joined, col
    assert np.array_equal(
        np.asarray(whole.counter_sums),
        np.asarray(parts[0].counter_sums) + np.asarray(parts[1].counter_sums),
    )
    assert np.array_equal(whole.edep[0], parts[0].edep[0] + parts[1].edep[0])
    n = len(whole.tally_components)
    esc = [
        tally_column_index(n, True, nm)
        for nm in ("nuclear_escaped_neutron", "nuclear_escaped_gamma")
    ]
    assert all(i in range(n) for i in esc) and set(NUCLEAR_TALLY_NAMES) >= {
        "nuclear_escaped_neutron", "nuclear_escaped_gamma"}  # fmt: skip
    s_whole = math.fsum(math.fsum(whole.tally_components[i]) for i in esc)
    s_parts = math.fsum(math.fsum(c for p in parts for c in p.tally_components[i]) for i in esc)
    assert s_whole == s_parts


# -- calibration step environment, expected case counts and sharding ------------------------------
def test_calibration_step_environment_never_inherits_the_rehearsal_override() -> None:
    rs = _load("run_suite")
    host = {"IONMC_V7R_REHEARSAL_DIR": "/tmp/evil", "PATH": "/usr/bin"}
    for name in ("18-pytest-v7r-calibration", "18-pytest-v7r-calibration-s1"):
        env = rs.step_environment(host, name, {"IONMC_V7R_FIXTURES": "required"})
        assert "IONMC_V7R_REHEARSAL_DIR" not in env and env["IONMC_V7R_FIXTURES"] == "required"
        assert env["PATH"] == "/usr/bin"
    other = rs.step_environment(host, "17-v7r-ref", {})
    assert other["IONMC_V7R_REHEARSAL_DIR"] == "/tmp/evil"  # only the calibration steps are pinned
    assert host["IONMC_V7R_REHEARSAL_DIR"] == "/tmp/evil"  # the input is not mutated


def test_rehearsal_fixture_digest_is_recorded_for_the_committed_directory(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    rs = _load("run_suite")
    monkeypatch.setattr(rs, "REPO", tmp_path)
    assert rs.rehearsal_fixture_digest() == "absent"
    d = tmp_path / rs.REHEARSAL_FIXTURE_DIR
    d.mkdir(parents=True)
    (d / "a.json").write_text("1")
    first = rs.rehearsal_fixture_digest()
    assert len(first) == 64 and first == rs.rehearsal_fixture_digest()
    (d / "a.json").write_text("2")
    assert rs.rehearsal_fixture_digest() != first
    assert rs.REHEARSAL_FIXTURE_DIR == "tests/ionmc/fixtures/v7r/rehearsal"


def _collect(*args: str) -> list[str]:
    import subprocess

    rs = _load("run_suite")
    out = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider",
         "-m", "calibration", *args],
        cwd=rs.REPO, capture_output=True, text=True, check=True,
    ).stdout  # fmt: skip
    return [ln for ln in out.splitlines() if "::" in ln]


def test_expected_calibration_cases_match_the_collection() -> None:
    rs = _load("run_suite")
    for target, n in rs.EXPECTED_CALIBRATION_CASES.items():
        assert len(_collect(target)) == n, target
    full = _collect(*rs.calibration_selection(None)[0])
    assert len(full) == sum(rs.EXPECTED_CALIBRATION_CASES.values()) == 13


def test_calibration_partition_is_disjoint_and_complete() -> None:
    rs = _load("run_suite")
    full = _collect(*rs.calibration_selection(None)[0])
    seen: list[str] = []
    for i, (targets, k_expr, want) in enumerate(rs.V7R_CALIBRATION_PARTITION):
        ids = _collect(*targets, "-k", k_expr)
        assert len(ids) == want, i
        seen += ids
    assert sorted(seen) == sorted(full) and len(set(seen)) == len(seen)


def test_sharded_registration_only_when_the_constant_is_set(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rs = _load("run_suite")
    names = [n.split("-", 1)[1] for n, _, _ in rs.suite_steps("lv5c", 1, 1.0)]
    assert "pytest-v7r-calibration" in names and not any("calibration-s" in n for n in names)
    monkeypatch.setattr(rs, "V7R_CALIBRATION_SHARDS", 3)
    steps = rs.suite_steps("lv5c", 1, 1.0)
    cal = [(n.split("-", 1)[1], c, e) for n, c, e in steps if "pytest-v7r-calibration" in n]
    assert [n for n, _, _ in cal] == [f"pytest-v7r-calibration-s{k}" for k in range(3)]
    for k, (_, cmd, env) in enumerate(cal):
        assert " ".join(cmd[3:]).startswith(
            " ".join(rs.calibration_args(*rs.calibration_selection(k)[:2]))
        )
        assert env["IONMC_V7R_FIXTURES"] == "required"
    for n in rs.full_step_names("lv5c", 2):
        assert rs.step_timeout_s("lv5c", n, 1500) >= 1500
    monkeypatch.setattr(rs, "V7R_CALIBRATION_SHARDS", 2)
    with pytest.raises(SystemExit):
        rs.suite_steps("lv5c", 1, 1.0)


def _sharded() -> dict:  # type: ignore[type-arg]
    steps = {"17-v7r-ref": _cal(), "19-v7r-combine": {"pass": True, "status": "executed"}}
    for k in range(3):
        steps[f"{18 + k}-pytest-v7r-calibration-s{k}"] = _cal(shard=k)
    return steps


def test_sharded_calibration_set_must_be_complete_exact_and_disjoint() -> None:
    sm = _load("summarize")
    ok = sm.slice_c_rows("lv5c", _sharded())["V7-R"]
    assert ok["pass"] is True and ok["calibration_steps"] == [
        f"pytest-v7r-calibration-s{k}" for k in range(3)]  # fmt: skip

    def verdict(steps: dict) -> bool:  # type: ignore[type-arg]
        return sm.slice_c_rows("lv5c", steps)["V7-R"]["calibration_pass"]

    missing = _sharded()
    del missing["20-pytest-v7r-calibration-s2"]
    assert verdict(missing) is False
    double = _sharded()
    double["23-pytest-v7r-calibration-s1"] = _cal(shard=1)  # shard 1 twice: double cover
    assert verdict(double) is False
    swapped = _sharded()
    swapped["19-pytest-v7r-calibration-s1"] = _cal(shard=0)  # s1 name with the s0 selection
    assert verdict(swapped) is False
    failed = _sharded()
    failed["20-pytest-v7r-calibration-s2"] = _cal(ok=False, shard=2)
    assert verdict(failed) is False
    short = _sharded()
    short["20-pytest-v7r-calibration-s2"] = _cal(shard=2, counts={"passed": 10})
    assert verdict(short) is False
    mixed = _sharded()
    mixed["24-pytest-v7r-calibration"] = _cal()  # un-sharded step next to the shards
    assert verdict(mixed) is False
    assert verdict({**_sharded(), "18-pytest-v7r-calibration-s3": _cal(shard=2)}) is False
    only = {"17-v7r-ref": _cal(), "19-v7r-combine": {"pass": True, "status": "executed"},
            "18-pytest-v7r-calibration": _cal()}  # fmt: skip
    assert verdict(only) is True  # the un-sharded default
