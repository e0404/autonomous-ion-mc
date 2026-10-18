"""CI-scale tests of the V7-R2 code (repair plan revision 4, ``validation/plans/v3-005d-v7r-calibration-repair.md``):
the per-gate ``alpha_tost`` parameter and its bit-identity at the default, the 48 x 2.25e6 layout and the seeds,
the exact operating characteristic of the 5400-pair / alpha = 0.005 gate, the kappa machinery and the fail-closed
calibration combine, and the rep-based timing rule. Pure numpy, no transport. Sandbox results here are never
evidence."""

# ruff: noqa: E501

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"
DD976C5 = REPO / ".ionmc-cache" / "transport-samples" / "dd976c5"


def _load(name: str, path: Path | None = None) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, path or SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def v5b() -> ModuleType:
    return _load("steps_v5b")


@pytest.fixture(scope="module")
def v7r(v5b: ModuleType) -> ModuleType:
    return _load("v7r")


@pytest.fixture(scope="module")
def cov() -> ModuleType:
    return _load("test_v7r_coverage_module", REPO / "tests" / "ionmc" / "test_v7r_coverage.py")


def _canon(d: Any) -> str:
    return hashlib.sha256(json.dumps(d, sort_keys=True, default=float).encode()).hexdigest()


def _binom_pmf(n: int, p: float) -> np.ndarray:
    k = np.arange(n + 1)
    logc = np.array([math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) for i in k])
    return np.exp(logc + k * math.log(p) + (n - k) * math.log1p(-p))


def _band(
    v5b: ModuleType, n: int, alpha: float, lo: float = 0.6527, hi: float = 0.7127
) -> tuple[int, int]:
    """Acceptance set [k_lo, k_hi] of the Clopper-Pearson TOST (paired gate form)."""
    k0 = _first(lambda k: v5b.cp_lower(k, n, alpha) >= lo, n)  # cp_lower increases with k
    k1 = _first(lambda k: v5b.cp_upper(k, n, alpha) > hi, n) - 1  # cp_upper increases with k
    return k0, k1


def _first(pred: Any, n: int) -> int:
    """Smallest k in 0..n+1 with ``pred(k)`` (monotone predicate; n + 1 if none), by bisection."""
    lo, hi = 0, n + 1
    while lo < hi:
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if mid <= n and pred(mid) else (mid + 1, hi)
    return lo


def _oc(v5b: ModuleType, n: int, alpha: float, c: float) -> float:
    k0, k1 = _band(v5b, n, alpha)
    return float(_binom_pmf(n, c)[k0 : k1 + 1].sum())


# -- 1. alpha_tost parameter: bit identity at the default ------------------------------------------
GOLDEN_SYNTHETIC = [  # sha256 of the canonical JSON of the gate documents, computed with the pre-parameter v7r.py (HEAD b7b269f3)
    "a54e73767bdcaa93b32464b6f67e067b0d99e889ecbdd1a3894af2911cfa39d2",
    "28080ea05157d2ae4dd8bbfe402c28bff96e33954b128c690695d7b59d268fda",
    "5f14e61c26ff7eddb320bc15f0f4367bd6c57ca6467a4cfc58c629f859ad0bd5",
    "d4054f9ce353d8f69380be0003315be28661d315182d75717f1108bff1f88680",
    "a54e73767bdcaa93b32464b6f67e067b0d99e889ecbdd1a3894af2911cfa39d2",
    "c70f6f7951277c08aff6016745d709bc9457c851d59ca73edd4d2241d3ea2d04",
    "0021a3f5f908c1bbc8b37220fc5780ed79d162604eb5c39927990434e7df44e7",
]


def _synthetic() -> np.ndarray:
    rng = np.random.default_rng(123)
    return rng.gamma(5.0, 100.0, (120, 20)) * (rng.random((120, 20)) > 0.2)


def test_default_alpha_reproduces_the_pre_parameter_gates_bit_for_bit(v7r: ModuleType) -> None:
    x = _synthetic()
    got = []
    for k in (1.0, 0.93, 1.07):
        pi = v7r.pair_ingredients(x, 777, 15, 199)
        si = v7r.single_ingredients(x, 777, 15, 199)
        got += [
            _canon(v7r.pair_gate(pi, k)),
            _canon(v7r.single_gate(si, x.mean(axis=1), k, boot_seed=5, n_boot=50)),
        ]
    got.append(
        _canon(
            v7r.escaped_neutral_verdict(x, 777, boot_seed=5, reps_per_shard=15, b=199, n_boot=50)
        )
    )
    assert got == GOLDEN_SYNTHETIC


def test_default_equals_explicit_0_04_and_the_patched_constant(
    v5b: ModuleType, v7r: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    x = _synthetic()
    pi = v7r.pair_ingredients(x, 777, 15, 199)
    si = v7r.single_ingredients(x, 777, 15, 199)
    base = (v7r.pair_gate(pi, 1.0), v7r.single_gate(si, x.mean(axis=1), boot_seed=5, n_boot=50))
    assert base[0]["alpha_tost"] == base[1]["alpha_tost"] == 0.04
    explicit = (v7r.pair_gate(pi, 1.0, alpha_tost=0.04),
                v7r.single_gate(si, x.mean(axis=1), boot_seed=5, n_boot=50, alpha_tost=0.04))  # fmt: skip
    assert _canon(explicit) == _canon(base)
    # the sandbox joint-power script patched steps_v5b.V7_ALPHA_TOST; None still reads it at call time
    monkeypatch.setattr(v7r.v5b, "V7_ALPHA_TOST", 0.005)
    patched = (v7r.pair_gate(pi, 1.0), v7r.single_gate(si, x.mean(axis=1), boot_seed=5, n_boot=50))
    param = (v7r.pair_gate(pi, 1.0, alpha_tost=0.005),
             v7r.single_gate(si, x.mean(axis=1), boot_seed=5, n_boot=50, alpha_tost=0.005))  # fmt: skip
    assert _canon(patched) == _canon(param) and patched[0]["alpha_tost"] == 0.005
    # a smaller alpha is strictly wider-bound: the paired interval of the same hits contains the 0.04 one
    assert param[0]["ci"][0] <= base[0]["ci"][0] and param[0]["ci"][1] >= base[0]["ci"][1]
    with pytest.raises(ValueError):
        v7r.pair_gate(pi, 1.0, alpha_tost=0.0)
    with pytest.raises(ValueError):
        v7r.single_gate(si, x.mean(axis=1), alpha_tost=0.7)


def _archived_sums() -> np.ndarray:
    files = [DD976C5 / str(30 + k) / "samples" / f"v7r-s{k}-escaped_neutral.npz" for k in range(16)]
    if not all(f.is_file() for f in files):
        pytest.skip("dd976c5a lv5c sidecars are not in the local cache")
    out = []
    for f in files:
        with np.load(f, allow_pickle=False) as z:
            out.append(z["escaped_neutral_block_sums"])
    return np.concatenate(out)


def test_archived_dd976c5a_escaped_neutral_verdict_is_reproduced(v7r: ModuleType) -> None:
    """The 16 archived dd976c5a sidecars (7200 replicates) through ``escaped_neutral_verdict`` at the default alpha:
    paired 2467/3600 = 0.6853 with CP bounds [0.6715, 0.6988] (repair plan 1.5), digest pinned from the pre-parameter code."""
    sums = _archived_sums()
    v = v7r.escaped_neutral_verdict(sums, 20481004, boot_seed=20497020)
    assert (v["covered"], v["n_pairs"]) == (2467, 3600)
    assert [round(c, 4) for c in v["ci"]] == [0.6715, 0.6988] and v["pass"] is True
    assert v["paired"]["alpha_tost"] == v["single"]["alpha_tost"] == 0.04
    assert _canon(v) == ARCHIVED_DIGEST


ARCHIVED_DIGEST = "b22acab60944247d80292255c33d0d35e009daa082424b09178a0d1dd4b068a6"


# -- 3. exact OC pins ---------------------------------------------------------------------------
def test_exact_oc_pins_of_the_5400_pair_gate(v5b: ModuleType) -> None:
    assert _band(v5b, 5400, 0.005) == (3615, 3762)
    assert _oc(v5b, 5400, 0.005, 0.6527) == pytest.approx(0.0050, abs=5e-5)
    assert _oc(v5b, 5400, 0.005, 0.7127) == pytest.approx(0.0050, abs=5e-5)
    assert _oc(v5b, 5400, 0.005, 0.6827) == pytest.approx(0.969, abs=5e-4)
    assert _band(v5b, 3600, 0.04) == (2401, 2517)  # the frozen V7-R band (plan 1.2)
    assert _oc(v5b, 3600, 0.04, 0.6527) == pytest.approx(0.0374, abs=5e-5)


def test_the_calibration_clause_table_at_r3200_is_unchanged(
    v5b: ModuleType, v7r: ModuleType
) -> None:
    assert v7r.V7R_CAL_REPS == 3200 and v7r.CAL_GATE_LIMIT == 0.04 and v7r.CAL_JOINT_LIMIT == 0.02
    assert v7r.CAL_POWER_MIN == 0.9
    assert _first(lambda k: v5b.cp_upper(k, 3200, 0.05) > 0.04, 3200) - 1 == 109
    assert _first(lambda k: v5b.cp_upper(k, 3200, 0.05) > 0.02, 3200) - 1 == 50
    assert _first(lambda k: v5b.cp_lower(k, 3200, 0.05) >= 0.9, 3200) == 2909


def test_planted_kappa_error_changes_the_oc_as_expected(v5b: ModuleType) -> None:
    """A kappa error of 2 sigma of the coverage (SE 0.0005) at the margin moves the false acceptance by about +-0.002."""
    margin = _oc(v5b, 5400, 0.005, 0.6527)
    inside = _oc(v5b, 5400, 0.005, 0.6527 + 2 * 0.0005)
    outside = _oc(v5b, 5400, 0.005, 0.6527 - 2 * 0.0005)
    assert (
        outside < margin < inside < 0.01 and inside > 1.2 * margin
    )  # still far under the 0.04 clause


# -- 2./4. layout, reference index, seeds ----------------------------------------------------------
def test_layout_constants(v5b: ModuleType, v7r: ModuleType) -> None:
    assert (v7r.V7R2_SHARDS, v7r.V7R2_SHARD_N, v7r.V7R2_REPS_PER_SHARD) == (48, 2_250_000, 225)
    assert (
        v7r.V7R2_SHARDS * v7r.V7R2_SHARD_N == 24 * 4_500_000
    )  # the same total histories as the 24 x 4.5e6 layout of revision 1
    assert v7r.V7R2_SHARD_N // v7r.V7R_BATCH_N // v7r.V7R_BLOCKS_PER_REP == v7r.V7R2_REPS_PER_SHARD
    assert (v7r.V7R2_REPLICATES, v7r.V7R2_PAIRS) == (10_800, 5400)
    ia, ib = v5b.v7_pair_indices(v7r.V7R2_REPLICATES)
    assert ia.size == 5400 and (ib - ia == 5400).all()
    shard = lambda i: i // v7r.V7R2_REPS_PER_SHARD  # noqa: E731
    assert (shard(ib) - shard(ia) == 24).all() and shard(ia).max() == 23  # s <-> s + 24, s = 0..23
    n_e = v5b.v7_heldout_split(v7r.V7R2_REPLICATES)
    assert n_e == 8100 and shard(np.arange(n_e, v7r.V7R2_REPLICATES)).tolist()[0] == 36
    assert set(shard(np.arange(n_e, 10_800)).tolist()) == set(v7r.V7R2_HELDOUT_SHARDS)
    assert (v7r.V7R2_CAL_PARTS, sum(v7r.cal_part_sizes(3200, v7r.V7R2_CAL_PARTS))) == (39, 3200)
    assert sorted(set(v7r.cal_part_sizes(3200, 39))) == [82, 83]
    assert v7r.V7R2_ALPHA_TOST == 0.005 and v7r.V7R2_KAPPA_N == 1_000_000
    # the V7-R constants are untouched
    assert (v7r.V7R_SHARDS, v7r.V7R_SHARD_N, v7r.V7R_REF_INDEX, v7r.CAL_PARTS, v7r.CAL_SEED) == (
        16,
        4_500_000,
        16,
        15,
        20507000,
    )


def test_reference_index_is_distinct_from_every_shard_index(v7r: ModuleType) -> None:
    assert v7r.V7R2_REF_INDEX == 48 and v7r.V7R2_REF_INDEX not in range(v7r.V7R2_SHARDS)
    assert v7r.V7R_REF_INDEX in range(
        v7r.V7R2_SHARDS
    )  # the collision the constant V7R_REF_INDEX would cause
    shards = [v7r.r2_shard_seed(k) for k in range(48)]
    assert shards == list(range(20547004, 20547052)) and v7r.r2_ref_seed() == 20547052
    assert v7r.r2_ref_seed() not in shards and len(set(shards)) == 48
    assert (
        v7r.r2_ref_seed(v7r.V7R2_FALLBACK_SHARDS) == 20547068
    )  # 7200 fallback: k = 0..63, reference k = 64
    assert v7r.V7R2_FALLBACK_REF_INDEX == 64 and v7r.r2_shard_seed(63, 64) == 20547067
    with pytest.raises(ValueError):
        v7r.r2_shard_seed(48)
    with pytest.raises(ValueError):
        v7r.r2_shard_seed(-1)


def test_calibration_and_kappa_seed_roots_and_guard(v7r: ModuleType) -> None:
    for root in (
        v7r.V7R2_KAPPA_SEED,
        v7r.V7R2_KAPPA_SEED + 1,
        v7r.V7R2_CAL_SEED,
        v7r.V7R2_CAL_SEED + 1,
        v7r.V7R2_A13_SEED,
    ):
        assert v7r.r2_seed_guard(root) == root
    assert (v7r.V7R2_KAPPA_SEED, v7r.V7R2_CAL_SEED, v7r.V7R2_A13_SEED) == (
        20980000,
        20982000,
        20984000,
    )
    for bad in (20979999, 20985000, 20507000, 20989901, 20547004):
        with pytest.raises(ValueError):
            v7r.r2_seed_guard(bad)
    seeds = {(s, k, p): v7r.r2_kappa_part_seeds(s, k, p)
             for s in (0, 1) for k in v7r.V7R2_KAPPA_KINDS for p in range(2)}  # fmt: skip
    words = [w for d in seeds.values() for w in d.values()]
    cal = [w for s in (0, 1) for p in range(39) for w in v7r.r2_cal_part_seeds(s, p).values()]
    assert len(set(words)) == len(words) and len(set(cal)) == len(cal) and not set(words) & set(cal)
    # deterministic and entropy-distinct from the V7-R calibration seeds
    assert v7r.r2_cal_part_seeds(1, 7) == v7r.r2_cal_part_seeds(1, 7)
    assert v7r.r2_cal_part_seeds(0, 0)["data"] != v7r.cal_part_seeds(0, 0)["data"]
    assert not {v7r.r2_cal_part_seeds(0, p)["data"] for p in range(15)} & {
        v7r.cal_part_seeds(0, p)["data"] for p in range(15)
    }


# -- kappa machinery -----------------------------------------------------------------------------
def _legacy_records(v7r: ModuleType, cov: ModuleType, sampler: Any, fit: dict[str, Any], n: int):  # type: ignore[no-untyped-def]
    recs, rows = {}, {}
    for s in v7r.CAL_SURROGATES:
        for kind in v7r.V7R2_KAPPA_KINDS:
            seeds = dict(zip(("data", "boot"), v7r.V7R_LEGACY_KAPPA_SEEDS[kind], strict=True))
            r = v7r.kappa_rows(sampler, kind, n, seeds, legacy=True)
            rows[(s, kind, 0)] = r
            recs[(s, kind, 0)] = {
                "schema": v7r.V7R2_KAPPA_SCHEMA, "mode": "legacy", "surrogate": s, "kind": kind, "part": 0,
                "n": n, "seeds": seeds, "fit": fit, "content_sha256": "x" * 64,
                "arrays": {k: hashlib.sha256(np.ascontiguousarray(v, dtype="<f8").tobytes()).hexdigest() for k, v in r.items()},
            }  # fmt: skip
    return recs, rows


def test_legacy_kappa_mode_equals_margin_kappas_bit_for_bit(
    v7r: ModuleType, cov: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Acceptance test 2 at CI scale: the per-row statistics + pooled bisection reproduce ``margin_kappas`` exactly
    (single: the 3:1 split; pair: ``2 n`` replicates). The five kappa of record of both surrogates are reproduced by the
    marked test ``test_legacy_kappa_of_record`` of ``test_v7r2_calibration.py`` (n = 2e4, minutes)."""
    sampler = cov.ZeroInflatedGamma(6.0, 0.05)
    n = 400
    monkeypatch.setattr(v7r, "V7R_LEGACY_KAPPA_N", n)
    recs, rows = _legacy_records(v7r, cov, sampler, {"k": 6.0}, n)
    out = v7r.kappa_combine(recs, rows, legacy=True)
    for s in v7r.CAL_SURROGATES:
        k = out["surrogates"][s]["kappas"]
        lo, hi = v7r.V7R_REGION
        assert k["paired_low"] == cov.margin_kappas(
            v7r, sampler, "pair", lo, cov.SEED_KAPPA + 1, n=n
        )
        assert k["paired_high"] == cov.margin_kappas(
            v7r, sampler, "pair", hi, cov.SEED_KAPPA + 1, n=n
        )
        assert k["single_low"] == cov.margin_kappas(v7r, sampler, "single", lo, cov.SEED_KAPPA, n=n)
        assert k["single_high"] == cov.margin_kappas(
            v7r, sampler, "single", hi, cov.SEED_KAPPA, n=n
        )
    assert out["surrogates"]["archive_fitted"]["coverage"]["single_low"]["n"] == v5b_split(n)


def v5b_split(n: int) -> int:
    return n - n // 4  # the 3:1 evaluation split of the legacy single mode


def test_full_sample_single_rows_extend_the_split_rows(v7r: ModuleType, cov: ModuleType) -> None:
    sampler = cov.ZeroInflatedGamma(6.0, 0.05)
    x = sampler(np.random.default_rng(5), 80)
    full = v7r.single_ingredients_all(x, 11, 80, b=99)
    split = v7r.single_ingredients(x, 11, reps_per_shard=80, b=99)
    assert full["m_hat"].size == 80 and split["m_hat"].size == 60
    for key in split:
        assert np.array_equal(full[key][:60], split[key], equal_nan=True)


def _r2_records(v7r: ModuleType, cov: ModuleType, part_n: dict[str, tuple[int, ...]]):  # type: ignore[no-untyped-def]
    sampler = cov.ZeroInflatedGamma(6.0, 0.05)
    recs, rows = {}, {}
    for si, s in enumerate(v7r.CAL_SURROGATES):
        for kind in v7r.V7R2_KAPPA_KINDS:
            for p, n in enumerate(part_n[kind]):
                seeds = v7r.r2_kappa_part_seeds(si, kind, p)
                r = v7r.kappa_rows(sampler, kind, n, seeds, b=199)
                rows[(s, kind, p)] = r
                recs[(s, kind, p)] = {
                    "schema": v7r.V7R2_KAPPA_SCHEMA, "mode": "v7r2", "surrogate": s, "kind": kind, "part": p,
                    "n": n, "seeds": seeds, "fit": {"k": 6.0}, "content_sha256": "y" * 64,
                    "arrays": {k: hashlib.sha256(np.ascontiguousarray(v, dtype="<f8").tobytes()).hexdigest() for k, v in r.items()},
                }  # fmt: skip
    return recs, rows


def test_kappa_combine_v7r2_mode_pools_full_samples_and_reports_the_se(
    v7r: ModuleType, cov: ModuleType
) -> None:
    part_n = {"single": (300,), "pair": (120, 120)}
    recs, rows = _r2_records(v7r, cov, part_n)
    out = v7r.kappa_combine(recs, rows, part_n=part_n)
    assert out["schema"] == v7r.V7R2_KAPPA_COMBINED_SCHEMA and out["mode"] == "v7r2"
    for s in v7r.CAL_SURROGATES:
        c = out["surrogates"][s]["coverage"]
        assert (
            c["single_low"]["n"] == 300 and c["paired_low"]["n"] == 240
        )  # full sample, pairs pooled over the parts
        for m in v7r.KAPPA_MARGINS:
            assert (
                abs(c[m]["coverage"] - c[m]["target"]) < 0.01
            )  # the empirical coverage hits the target (grid 1/n)
            assert c[m]["se"] == pytest.approx(
                math.sqrt(c[m]["coverage"] * (1 - c[m]["coverage"]) / c[m]["n"])
            )
            assert c[m]["se_ok"] is False  # n is tiny: SE > 0.0005
    assert out["se_ok"] is False
    with pytest.raises(SystemExit, match="coverage SE"):
        v7r.r2_calibration_kappas(out, "archive_fitted")
    # the production layout: SE of every margin <= 0.0005 at n = 1e6
    for target in (0.6527, 0.7127):
        assert math.sqrt(target * (1 - target) / v7r.V7R2_KAPPA_N) <= v7r.V7R2_KAPPA_SE_MAX
    assert (
        sum(v7r.V7R2_KAPPA_PART_N["pair"])
        == sum(v7r.V7R2_KAPPA_PART_N["single"])
        == v7r.V7R2_KAPPA_N
    )


def test_kappa_combine_planted_error_in_the_coverage_is_about_two_se(
    v7r: ModuleType, cov: ModuleType
) -> None:
    sampler = cov.ZeroInflatedGamma(6.0, 0.05)
    r = v7r.kappa_rows(sampler, "single", 4000, v7r.r2_kappa_part_seeds(0, "single", 0), b=199)
    k = v7r.kappa_solve(r, 0.6527)
    c0 = v7r.kappa_coverage(r, k)
    se = math.sqrt(c0 * (1 - c0) / 4000)
    k2 = v7r.kappa_solve(r, 0.6527 + 2 * se)
    assert k2 > k and v7r.kappa_coverage(r, k2) - c0 == pytest.approx(2 * se, abs=2e-3)


def test_kappa_combine_fails_closed(v7r: ModuleType, cov: ModuleType) -> None:
    part_n = {"single": (200,), "pair": (60, 60)}
    recs, rows = _r2_records(v7r, cov, part_n)
    ok = v7r.kappa_combine(recs, rows, part_n=part_n)
    assert ok["mode"] == "v7r2"
    k0 = ("archive_fitted", "pair", 1)

    def fails(match: str, mutate: Any) -> None:
        r2 = {k: dict(v) for k, v in recs.items()}
        rw = {k: {a: b.copy() for a, b in v.items()} for k, v in rows.items()}
        mutate(r2, rw)
        with pytest.raises(SystemExit, match=match):
            v7r.kappa_combine(r2, rw, part_n=part_n)

    fails("missing", lambda r, w: (r.pop(k0), w.pop(k0)))
    fails("wrong schema or mode", lambda r, w: r[k0].update(mode="legacy"))
    fails("is not 60", lambda r, w: r[k0].update(n=61))
    fails("seeds differ", lambda r, w: r[k0].update(seeds={"data": 1, "boot": 2}))
    fails("fit differs", lambda r, w: r[k0].update(fit={"k": 7.0}))
    fails("sha256", lambda r, w: w[k0]["a"].__setitem__(0, w[k0]["a"][0] + 1.0))
    fails("exactly a / t_lo / t_hi", lambda r, w: w[k0].pop("t_hi"))
    with pytest.raises(
        SystemExit, match="missing"
    ):  # a part layout different from the declared one
        v7r.kappa_combine(recs, rows, part_n={"single": (200,), "pair": (60,)})


# -- calibration combine: kappa record hash, part counts, R ----------------------------------------
def _cal_records(v7r: ModuleType, sha: str | None) -> list[dict[str, Any]]:
    sizes = v7r.cal_part_sizes(3200, v7r.V7R2_CAL_PARTS)
    out = []
    for si, s in enumerate(v7r.CAL_SURROGATES):
        for p in range(v7r.V7R2_CAL_PARTS):
            n = sizes[p]
            margins = {m: {"kappa": 1.0 if m == "nominal" else 0.9, "accepted_joint": 0,
                           "accepted_paired": 0, "accepted_single": 0} for m in v7r.CAL_MARGINS}  # fmt: skip
            margins["nominal"].update(accepted_joint=n, accepted_paired=n, accepted_single=n)
            rec: dict[str, Any] = {"surrogate": s, "part": p, "n": n, "margins": margins, "fit": {"k": 1.0},
                                   "reps_env": 3200, "seeds": v7r.r2_cal_part_seeds(si, p)}  # fmt: skip
            if sha is not None:
                rec["kappa_record_sha256"] = sha
            out.append(rec)
    return out


def test_r2_calibration_combine_is_fail_closed(v7r: ModuleType) -> None:
    sha = "ab" * 32
    kw: dict[str, Any] = {
        "parts": v7r.V7R2_CAL_PARTS,
        "seed_fn": v7r.r2_cal_part_seeds,
        "kappa_record_sha256": sha,
    }
    good = v7r.calibration_combine(_cal_records(v7r, sha), **kw)
    assert good["pass"] is True and good["part_sizes"] == list(v7r.cal_part_sizes(3200, 39))
    assert good["surrogates"]["archive_fitted"]["margins"]["nominal"]["n"] == 3200
    recs = _cal_records(v7r, sha)
    del recs[5]["kappa_record_sha256"]  # missing kappa hash
    with pytest.raises(SystemExit, match="kappa_record_sha256"):
        v7r.calibration_combine(recs, **kw)
    recs = _cal_records(v7r, sha)
    recs[7]["kappa_record_sha256"] = "cd" * 32  # a part with a different kappa record
    with pytest.raises(SystemExit, match="kappa_record_sha256"):
        v7r.calibration_combine(recs, **kw)
    recs = _cal_records(v7r, sha)
    recs[3]["n"] -= 1  # unequal part counts
    with pytest.raises(SystemExit, match="is not"):
        v7r.calibration_combine(recs, **kw)
    recs = _cal_records(v7r, sha)
    recs[3]["reps_env"] = 800  # R != 3200
    with pytest.raises(SystemExit, match="reps_env"):
        v7r.calibration_combine(recs, **kw)
    with pytest.raises(SystemExit, match="missing"):
        v7r.calibration_combine(_cal_records(v7r, sha)[:-1], **kw)
    recs = _cal_records(v7r, sha)
    recs[2]["seeds"] = v7r.cal_part_seeds(0, 2)  # the V7-R seeds are not the V7-R2 ones
    with pytest.raises(SystemExit, match="seeds differ"):
        v7r.calibration_combine(recs, **kw)


def test_legacy_calibration_combine_defaults_are_unchanged(v7r: ModuleType) -> None:
    """Called with the V7-R defaults it still demands 15 parts per surrogate with the V7-R seeds and ignores any hash."""
    sizes = v7r.cal_part_sizes(3200)
    assert len(sizes) == 15 and sizes[0] == 214 and sizes[14] == 213
    recs = _cal_records(v7r, None)
    with pytest.raises(SystemExit, match="unexpected|missing"):
        v7r.calibration_combine(recs)  # 39-part records are not a V7-R layout


# -- 5. rep-based timing rule ---------------------------------------------------------------------
def test_rep_plan_calibration_and_kappa_parts_fit_the_1100_s_target(v7r: ModuleType) -> None:
    plans = v7r.r2_rep_plans()
    assert plans["ok"] is True
    cal = plans["plans"]["v7r2 calibration (3200 reps x 39 parts per surrogate)"]
    assert (cal["largest_part_items"], cal["parts_needed"], cal["declared_parts"]) == (83, 39, 39)
    assert (
        cal["planned_s_per_part"] == pytest.approx(83 * 12.2 + 70.0)
        and cal["planned_s_per_part"] <= 1100.0
    )
    assert all(
        p["planned_s_per_part"] <= v7r.R2_STEP_TARGET_S < v7r.R2_STEP_LIMIT_S
        for p in plans["plans"].values()
    )
    assert (
        len(plans["plans"]) == 4
    )  # one calibration family + the three kappa step shapes (six steps)
    assert (
        sum(1 for _ in v7r.V7R2_KAPPA_PART_N["single"])
        + sum(1 for _ in v7r.V7R2_KAPPA_PART_N["pair"])
        == 3
    )  # per surrogate


def test_rep_plan_rejects_a_part_that_is_too_slow_or_too_few(v7r: ModuleType) -> None:
    assert v7r.r2_rep_plans(cal_s_per_rep=14.0)["ok"] is False  # 83 * 14 + 70 = 1232 s > 1100 s
    assert (
        v7r.r2_rep_plans(cal_s_per_rep=14.0)["plans"][
            "v7r2 calibration (3200 reps x 39 parts per surrogate)"
        ]["parts_needed"]
        == 44
    )
    slow = v7r.rep_plan("x", 3200, 12.2, 70.0, 30)
    assert (
        slow["ok"] is False and slow["parts_needed"] == 39 and slow["planned_s_per_part"] > 1100.0
    )
    assert v7r.rep_plan("x", 3200, 12.2, 70.0, 38)["ok"] is False  # largest part 85: 1107 s
    assert v7r.rep_plan("x", 3200, 12.2, 70.0, 39)["ok"] is True
    assert v7r.rep_plan("x", 3200, 12.2, 70.0, 40)["ok"] is True
    with pytest.raises(ValueError):
        v7r.rep_plan("x", 10, 1.0, 2000.0, 1)  # fixed cost above the target
    with pytest.raises(ValueError):
        v7r.rep_plan("x", 10, 1.0, 10.0, 0)


# -- the V7-R2 combine (design="v7r2") ------------------------------------------------------------
def test_v7r2_combine_layout_alpha_and_fail_closed(tmp_path: Path) -> None:
    from tests.ionmc.test_v7r_infra import _load as load_infra
    from tests.ionmc.test_v7r_infra import _synthetic_v7r

    v5c = load_infra("steps_v5c")
    v7r = v5c.v7r
    shards, ref, _ = _synthetic_v7r(v5c, v7r, tmp_path, n_shards=4, rps=225, seed_base=20531004)
    ref = {**ref, "seed": 20531004 + 16000 + 4}  # reference index = the shard count
    cal = {
        "step": "v7r2-calibration-combine",
        "pass": True,
        "reps": 3200,
        "content_sha256": "0" * 64,
    }
    doc = v5c.v7r_combine(
        shards, ref, 20531004, n_shards=4, n_boot=100, calibration=cal, design="v7r2"
    )
    en = doc["estimators"]["escaped_neutral"]
    assert doc["row"] == "V7-R2" and doc["replicates"] == 900 and doc["rule"]["pairs"] == 450
    assert doc["rule"]["alpha_tost"] == 0.005 and doc["rule"]["alpha_tost_profiles"] == 0.04
    assert en["paired"]["alpha_tost"] == en["single"]["alpha_tost"] == 0.005
    assert doc["rule"]["pairing"] == "(shard s, replicate j) with (shard s + 2, replicate j)"
    assert doc["rule"]["heldout"] == "shards 3-3" and doc["seeds"]["reference"] == [
        20531004 + 16004
    ]
    assert doc["calibration_status"]["partial"] == "v7r2-calibration-status.json"
    assert "v7r2-calibration-status" in doc["row_pass_condition"]
    with pytest.raises(
        SystemExit, match="calibration status"
    ):  # a V7-R status partial is not accepted
        v5c.v7r_combine(
            shards,
            ref,
            20531004,
            n_shards=4,
            calibration={**cal, "step": "v7r-calibration-combine"},
            design="v7r2",
        )
    with pytest.raises(
        SystemExit, match="seed/mismatch|reference"
    ):  # the V7-R reference index (16) is not the V7-R2 one
        v5c.v7r_combine(
            shards,
            {**ref, "seed": 20531004 + 16016},
            20531004,
            n_shards=4,
            calibration=cal,
            design="v7r2",
        )
    s2, r2, _ = _synthetic_v7r(v5c, v7r, tmp_path / "x", n_shards=4, rps=75, seed_base=20531004)
    with pytest.raises(SystemExit, match="225"):
        v5c.v7r_combine(
            s2, {**r2, "seed": ref["seed"]}, 20531004, n_shards=4, calibration=cal, design="v7r2"
        )
    with pytest.raises(SystemExit, match="design"):
        v5c.v7r_combine(shards, ref, 20531004, n_shards=4, calibration=cal, design="v7r3")
