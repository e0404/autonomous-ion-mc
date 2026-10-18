"""V7-R2 kappa and calibration parts (repair plan revision 4, 3.4, 6.5; planning notes PN-V7R2-1..5).

Host parts (marker ``calibration``, excluded from CI; run by the lv5d pytest steps once the suite is registered):

* ``test_v7r2_kappa_part``: 6 cases (per surrogate: ``single`` 1 part of 1e6 replicates, ``pair`` 2 parts of 5e5 pairs);
  each writes ``v7r-kappa-<surrogate>-<kind>-p<k>.json`` (sealed) and the sidecar ``.npz`` of the per-row statistics
  into ``IONMC_V7R_CAL_OUT`` and prints ``PARTIAL <name> <sha256>``;
* ``test_v7r2_calibration_part``: 78 cases (surrogate x 39 parts of 83/82 repetitions of R = 3200, 10 800 replicates of
  225 per shard, alpha_tost 0.005); each REFUSES to run without the sealed combined kappa record ``v7r-kappa.json`` in
  ``IONMC_V7R_KAPPA_DIR``, uses its kappa, checks that its surrogate ``fit`` block equals the kappa record's, records the
  record's ``content_sha256`` as ``kappa_record_sha256`` and prints ``PARTIAL`` for its sealed record;
* ``test_v7r2_legacy_kappa_of_record``: the five kappa of record of the archive-fitted surrogate (n = 2e4, old seeds,
  3:1 split) are reproduced bit for bit through the V7-R2 kappa machinery.

CI-scale tests below (not marked) run the same helpers on tiny layouts. Sandbox results are never evidence. The verdict
is made by ``v7r2-calibration-combine`` (``steps_v5c``)."""

# ruff: noqa: E501

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import platform
import sys
import time
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"


def _load(name: str, path: Path) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def cov() -> ModuleType:
    return _load("test_v7r_coverage_module_r2", REPO / "tests" / "ionmc" / "test_v7r_coverage.py")


@pytest.fixture(scope="module")
def s5c() -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    return _load("steps_v5c_r2", SCRIPTS / "steps_v5c.py")


@pytest.fixture(scope="module")
def v7r(s5c: ModuleType) -> ModuleType:
    return s5c.v7r  # the instance the steps use


def _env_run() -> dict[str, str]:
    return {
        "run_sha": os.environ.get("IONMC_RUN_SHA", "unknown"),
        "suite": os.environ.get("IONMC_RUN_SUITE", "unknown"),
    }


def _write_sealed(cov: ModuleType, out_dir: Path, name: str, rec: dict[str, Any]) -> dict[str, Any]:
    sealed = cov._seal(rec)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / name).write_text(json.dumps(sealed, sort_keys=True))
    print(f"PARTIAL {name} {sealed['content_sha256']}", flush=True)
    return sealed


def r2_surrogate(cov: ModuleType, v7r: ModuleType, name: str) -> tuple[Any, dict[str, Any]]:
    """``(sampler, fit block)`` of a surrogate; the fit block also binds the fixture: the sha256 of the archived
    replicate file (archive-fitted) or the V7-R2 rehearsal pool with its partial digests and sidecar sha256."""
    if name == "archive_fitted":
        sampler, fit = cov._surrogate(v7r, name)
        return sampler, {
            **fit,
            "fixture_sha256": hashlib.sha256(cov.ARCHIVE_FIXTURE.read_bytes()).hexdigest(),
        }
    s5c = cov._load("steps_v5c")
    try:
        blocks, info = cov.load_rehearsal(
            cov.rehearsal_dir(), s5c, cov._load("run_suite"), layout="v7r2"
        )
    except cov.RehearsalError as exc:
        cov._skip_or_fail(f"V7-R2 rehearsal fixture rejected: {exc}")
        raise
    pool = blocks.ravel()
    return cov.PooledBlocks(pool), {
        "surrogate": "pooled rehearsal block sums (V7-R2, 4 x 225 x 20)",
        "blocks": int(pool.size),
        **info,
    }


def run_kappa_part(cov: ModuleType, v7r: ModuleType, sampler: Any, fit: dict[str, Any], surrogate: str, kind: str, part: int,
                   out_dir: Path, *, n: int | None = None, b: int | None = None) -> dict[str, Any]:  # fmt: skip
    """One kappa part: per-row statistics of ``n`` replicates / pairs (layout value by default), the sidecar npz and
    the sealed record; prints ``PARTIAL``."""
    si = v7r.CAL_SURROGATES.index(surrogate)
    n = v7r.V7R2_KAPPA_PART_N[kind][part] if n is None else n
    seeds = v7r.r2_kappa_part_seeds(si, kind, part)
    t0 = time.perf_counter()
    rows = v7r.kappa_rows(sampler, kind, n, seeds, b=v7r.V7R_B if b is None else b)
    wall = time.perf_counter() - t0
    nm, npz_name = v7r.kappa_part_names(surrogate, kind, part)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_dir / npz_name,
        **{k: np.ascontiguousarray(v, dtype=np.float64) for k, v in rows.items()},
    )
    rec = {"schema": v7r.V7R2_KAPPA_SCHEMA, "mode": "v7r2", "surrogate": surrogate, "kind": kind, "part": part, "n": n,
           "n_rows": int(rows["a"].size), "seeds": seeds, "seed_roots": {"data": v7r.V7R2_KAPPA_SEED, "boot": v7r.V7R2_KAPPA_SEED + 1},
           "fit": fit, "arrays": {k: hashlib.sha256(np.ascontiguousarray(v, dtype="<f8").tobytes()).hexdigest() for k, v in rows.items()},
           "sidecar": {"file": npz_name, "sha256": hashlib.sha256((out_dir / npz_name).read_bytes()).hexdigest()},
           **_env_run(), "wall_s": wall, "wall_per_row_s": wall / n, "rng": v7r.rng_record(), "B": v7r.V7R_B,
           "python": platform.python_version()}  # fmt: skip
    return _write_sealed(cov, out_dir, nm, rec)


def load_kappa_doc(cov: ModuleType, directory: Path | None) -> dict[str, Any]:
    """The sealed combined kappa record; a part REFUSES to run without it (``pytest.fail``, never a skip)."""
    hit = (
        [
            p
            for p in (directory / "v7r-kappa.json", directory / "samples" / "v7r-kappa.json")
            if p.is_file()
        ]
        if directory
        else []
    )
    if len(hit) != 1:
        pytest.fail(
            "the sealed kappa record v7r-kappa.json is missing in IONMC_V7R_KAPPA_DIR: the part refuses to run"
        )
    doc = json.loads(hit[0].read_text())
    again = cov._seal({k: v for k, v in doc.items() if k != "content_sha256"})["content_sha256"]
    if (
        doc.get("content_sha256") != again
        or doc.get("schema") != "ionmc-v7r-kappa-1"
        or not doc.get("pass")
    ):
        pytest.fail(
            "the kappa record is not a sealed passing combined record: the part refuses to run"
        )
    if any(doc.get(k) != v for k, v in _env_run().items()):  # exact-SHA evidence rule
        pytest.fail("the kappa record belongs to another run_sha / suite: the part refuses to run")
    return doc


def run_cal_part(cov: ModuleType, v7r: ModuleType, v5b: ModuleType, sampler: Any, fit: dict[str, Any], surrogate: str,
                 part: int, kappa_doc: dict[str, Any], out_dir: Path, *, reps: int | None = None, parts: int | None = None,
                 replicates: int | None = None, reps_per_shard: int | None = None, b: int = 1999,
                 n_boot: int | None = None) -> dict[str, Any]:  # fmt: skip
    """One calibration part: the slice of the Monte Carlo at the kappas of the kappa record (alpha_tost 0.005)."""
    reps = v7r.V7R_CAL_REPS if reps is None else reps
    parts = v7r.V7R2_CAL_PARTS if parts is None else parts
    si = v7r.CAL_SURROGATES.index(surrogate)
    kfit = json.loads(json.dumps(kappa_doc["surrogates"][surrogate]["fit"]))
    if json.loads(json.dumps(fit)) != kfit:
        pytest.fail(
            "the surrogate fit block differs from the one of the kappa record: the part refuses to run"
        )
    kappas = v7r.r2_calibration_kappas(kappa_doc, surrogate)
    n = v7r.cal_part_sizes(reps, parts)[part]
    seeds = v7r.r2_cal_part_seeds(si, part)
    t0 = time.perf_counter()
    rates = cov.run_calibration(
        v7r, v5b, sampler, kappas, n, replicates=v7r.V7R2_REPLICATES if replicates is None else replicates,
        reps_per_shard=v7r.V7R2_REPS_PER_SHARD if reps_per_shard is None else reps_per_shard, b=b, n_boot=n_boot,
        data_seed0=seeds["data"], boot_seed0=seeds["boot"], alpha_tost=v7r.V7R2_ALPHA_TOST)  # fmt: skip
    wall = time.perf_counter() - t0
    margins = {}
    for m, kap in kappas.items():
        p = rates[m]
        assert p.shape == (n, 2) and np.isfinite(kap) and kap > 0.0
        margins[m] = {"kappa": kap, "accepted_joint": int(p.all(axis=1).sum()), "accepted_paired": int(p[:, 0].sum()),
                      "accepted_single": int(p[:, 1].sum())}  # fmt: skip
        assert margins[m]["accepted_joint"] <= min(
            margins[m]["accepted_paired"], margins[m]["accepted_single"]
        )
    rec = {"schema": v7r.CAL2_SCHEMA, "surrogate": surrogate, "part": part, "n": n, "margins": margins, "fit": fit,
           "reps_env": reps, "seeds": seeds, "seed_roots": {"data": v7r.V7R2_CAL_SEED, "boot": v7r.V7R2_CAL_SEED + 1},
           "kappa_record_sha256": kappa_doc["content_sha256"], "alpha_tost": v7r.V7R2_ALPHA_TOST,
           "layout": {"replicates": v7r.V7R2_REPLICATES if replicates is None else replicates,
                      "reps_per_shard": v7r.V7R2_REPS_PER_SHARD if reps_per_shard is None else reps_per_shard, "parts": parts},
           **_env_run(), "wall_mc_s": wall, "wall_per_replicate_s": wall / n, "rng": v7r.rng_record(), "B": b,
           "boot_n": v5b.V7_BOOT_N if n_boot is None else n_boot, "python": platform.python_version()}  # fmt: skip
    return _write_sealed(cov, out_dir, f"v7r2-calibration-{surrogate}-p{part:02d}.json", rec)


# -- host parts ----------------------------------------------------------------------------------
KAPPA_CASES = [
    (s, k, p)
    for s in ("archive_fitted", "rehearsal_block_sums")
    for k in ("single", "pair")
    for p in range({"single": 1, "pair": 2}[k])
]
CAL2_CASES = [(s, p) for s in ("archive_fitted", "rehearsal_block_sums") for p in range(39)]


def _out_dir(cov: ModuleType) -> Path:
    d = os.environ.get("IONMC_V7R_CAL_OUT")
    if cov._fixtures_required() and not d:
        pytest.fail("IONMC_V7R_CAL_OUT is not set in required mode: the record has nowhere to go")
    return Path(d or ".")


@pytest.mark.calibration
@pytest.mark.parametrize(
    ("surrogate", "kind", "part"), KAPPA_CASES, ids=[f"{s}-{k}-p{p}" for s, k, p in KAPPA_CASES]
)
def test_v7r2_kappa_part(
    cov: ModuleType, v7r: ModuleType, surrogate: str, kind: str, part: int
) -> None:
    sampler, fit = r2_surrogate(cov, v7r, surrogate)
    rec = run_kappa_part(cov, v7r, sampler, fit, surrogate, kind, part, _out_dir(cov))
    assert rec["n_rows"] == v7r.V7R2_KAPPA_PART_N[kind][part]


@pytest.mark.calibration
@pytest.mark.parametrize(
    ("surrogate", "part"), CAL2_CASES, ids=[f"{s}-part{p:02d}" for s, p in CAL2_CASES]
)
def test_v7r2_calibration_part(cov: ModuleType, v7r: ModuleType, surrogate: str, part: int) -> None:
    kdir = os.environ.get("IONMC_V7R_KAPPA_DIR")
    kdoc = load_kappa_doc(cov, Path(kdir) if kdir else None)
    sampler, fit = r2_surrogate(cov, v7r, surrogate)
    run_cal_part(cov, v7r, v7r.v5b, sampler, fit, surrogate, part, kdoc, _out_dir(cov))


@pytest.mark.calibration
def test_v7r2_legacy_kappa_of_record(cov: ModuleType, v7r: ModuleType) -> None:
    """The five kappa of record of the archive-fitted surrogate (repair plan 2.1: 0.938606 / 1.066375 paired,
    0.929174 / 1.054671 single), reproduced bit for bit by the legacy mode of the kappa machinery."""
    sampler, fit = cov._surrogate(v7r, "archive_fitted")
    lo, hi = v7r.V7R_REGION
    out = {}
    for kind in ("single", "pair"):
        seeds = dict(zip(("data", "boot"), v7r.V7R_LEGACY_KAPPA_SEEDS[kind], strict=True))
        rows = v7r.kappa_rows(sampler, kind, v7r.V7R_LEGACY_KAPPA_N, seeds, legacy=True)
        out[kind] = (v7r.kappa_solve(rows, lo), v7r.kappa_solve(rows, hi))
    assert [round(x, 6) for x in (*out["pair"], *out["single"])] == [
        0.938606,
        1.066375,
        0.929174,
        1.054671,
    ]
    assert out["pair"][0] == cov.margin_kappas(
        v7r, sampler, "pair", lo, cov.SEED_KAPPA + 1
    )  # same code path value


# -- CI-scale tests of the helpers -----------------------------------------------------------------
@pytest.fixture
def tiny(
    monkeypatch: pytest.MonkeyPatch, cov: ModuleType, v7r: ModuleType
) -> tuple[Any, dict[str, Any]]:
    monkeypatch.setattr(v7r, "V7R2_KAPPA_PART_N", {"single": (80,), "pair": (40, 40)})
    monkeypatch.setattr(
        v7r, "V7R2_KAPPA_SE_MAX", 1.0
    )  # tiny samples: the SE requirement is not under test here
    return cov.ZeroInflatedGamma(6.0, 0.05), {
        "surrogate": "tiny",
        "k": 6.0,
        "fixture_sha256": "ab" * 32,
    }


def _make_kappa_dir(cov: ModuleType, v7r: ModuleType, s5c: ModuleType, tiny: Any, d: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:  # fmt: skip
    sampler, fit = tiny
    for s in v7r.CAL_SURROGATES:
        for kind in v7r.V7R2_KAPPA_KINDS:
            for p in range(len(v7r.V7R2_KAPPA_PART_N[kind])):
                run_kappa_part(cov, v7r, sampler, fit, s, kind, p, d, b=99)
    sealed: dict[str, Any] = {}

    def fake_write(a: Any, name: str, doc: dict[str, Any]) -> str:
        sealed.update(cov._seal({**_env_run(), **doc}))  # write_partial_c binds run_sha / suite
        path = Path(a.out_dir) / f"{name}.json"
        path.write_text(json.dumps(sealed, sort_keys=True))
        print(f"PARTIAL {path.name} {sealed['content_sha256']}")
        return str(path)

    monkeypatch.setattr(s5c, "write_partial_c", fake_write)
    monkeypatch.setattr(s5c.v5, "attestation_block", lambda a: {})
    monkeypatch.setattr(s5c.v5b, "finish5b", lambda doc, *a: 0 if doc.get("pass") else 1)
    a = SimpleNamespace(dirs=[str(d)], out_dir=str(d))
    assert s5c.step_v7r_kappa_combine(a) == 0
    return json.loads((d / "v7r-kappa.json").read_text())


def test_kappa_parts_print_partial_lines_and_combine_end_to_end(
    cov: ModuleType, v7r: ModuleType, s5c: ModuleType, tiny: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:  # fmt: skip
    kdoc = _make_kappa_dir(cov, v7r, s5c, tiny, tmp_path, monkeypatch)
    out = capsys.readouterr().out
    assert out.count("PARTIAL v7r-kappa-") == 6 + 0 and "PARTIAL v7r-kappa.json" in out
    assert (
        kdoc["schema"] == "ionmc-v7r-kappa-1" and kdoc["pass"] is True and len(kdoc["parts"]) == 6
    )
    rec = json.loads((tmp_path / "v7r-kappa-archive_fitted-pair-p1.json").read_text())
    assert rec["seeds"] == v7r.r2_kappa_part_seeds(0, "pair", 1) and rec["seed_roots"] == {
        "data": 20980000,
        "boot": 20980001,
    }
    assert (
        rec["fit"]["fixture_sha256"] == "ab" * 32
        and rec["n_rows"] == 40
        and rec["sidecar"]["file"].endswith("-p1.npz")
    )
    # tampering with a sidecar is caught by the combine
    npz = tmp_path / "v7r-kappa-archive_fitted-single-p0.npz"
    with np.load(npz) as z:
        arr = {k: z[k].copy() for k in z.files}
    arr["a"][0] += 1.0
    np.savez(npz, **arr)
    with pytest.raises(SystemExit, match="sha256"):
        s5c.step_v7r_kappa_combine(SimpleNamespace(dirs=[str(tmp_path)], out_dir=str(tmp_path)))


def test_calibration_parts_consume_the_kappa_record_and_refuse_otherwise(
    cov: ModuleType, v7r: ModuleType, s5c: ModuleType, tiny: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:  # fmt: skip
    sampler, fit = tiny
    with pytest.raises(pytest.fail.Exception, match="refuses to run"):
        load_kappa_doc(cov, tmp_path)  # no kappa record
    with pytest.raises(pytest.fail.Exception, match="refuses to run"):
        load_kappa_doc(cov, None)
    kdoc = _make_kappa_dir(cov, v7r, s5c, tiny, tmp_path, monkeypatch)
    assert load_kappa_doc(cov, tmp_path)["content_sha256"] == kdoc["content_sha256"]
    capsys.readouterr()
    kw: dict[str, Any] = {
        "reps": 78,
        "parts": 39,
        "replicates": 40,
        "reps_per_shard": 5,
        "b": 49,
        "n_boot": 20,
    }
    recs = {}
    for s in v7r.CAL_SURROGATES:
        for p in range(39):
            recs[(s, p)] = run_cal_part(cov, v7r, v7r.v5b, sampler, fit, s, p, kdoc, tmp_path, **kw)
    out = capsys.readouterr().out
    assert out.count("PARTIAL v7r2-calibration-") == 78
    r = recs[("rehearsal_block_sums", 7)]
    assert (
        r["kappa_record_sha256"] == kdoc["content_sha256"]
        and r["alpha_tost"] == 0.005
        and r["n"] == 2
    )
    assert set(r["margins"]) == set(v7r.CAL_MARGINS) and r["margins"]["nominal"]["kappa"] == 1.0
    assert (
        r["margins"]["paired_low"]["kappa"]
        == kdoc["surrogates"]["rehearsal_block_sums"]["kappas"]["paired_low"]
    )
    # a part with another fit block refuses
    with pytest.raises(pytest.fail.Exception, match="fit block"):
        run_cal_part(
            cov, v7r, v7r.v5b, sampler, {**fit, "k": 7.0}, "archive_fitted", 0, kdoc, tmp_path, **kw
        )
    # the consumer: 78 records + kappa record -> status partial; fit mismatch against the kappa record is fail closed
    monkeypatch.setattr(v7r, "V7R_CAL_REPS", 78)
    monkeypatch.setattr(
        v7r, "r2_part_layout", lambda: {"replicates": 40, "reps_per_shard": 5, "parts": 39}
    )  # tiny test layout
    a = SimpleNamespace(dirs=[str(tmp_path)], out_dir=str(tmp_path))
    assert s5c.step_v7r2_calibration_combine(a) in (0, 1)
    status = json.loads((tmp_path / "v7r2-calibration-status.json").read_text())
    assert (
        status["step"] == "v7r2-calibration-combine"
        and status["kappa_record"]["content_sha256"] == kdoc["content_sha256"]
    )
    assert status["records"] == 78 and status["reps"] == 78
    bad = json.loads((tmp_path / "v7r2-calibration-archive_fitted-p03.json").read_text())
    bad["fit"] = {**bad["fit"], "k": 9.0}
    (tmp_path / "v7r2-calibration-archive_fitted-p03.json").write_text(
        json.dumps(cov._seal({k: v for k, v in bad.items() if k != "content_sha256"}))
    )
    with pytest.raises(SystemExit, match="fit"):
        s5c.step_v7r2_calibration_combine(a)
    ok_records = [
        json.loads((tmp_path / f"v7r2-calibration-{s}-p{p:02d}.json").read_text())
        for s in v7r.CAL_SURROGATES
        for p in range(39)
    ]
    with pytest.raises(SystemExit, match="fit block"):
        v7r.calibration_combine(ok_records, reps=78, parts=39, seed_fn=v7r.r2_cal_part_seeds,
                                kappa_record_sha256=kdoc["content_sha256"],
                                kappa_fits={s: {"other": 1} for s in v7r.CAL_SURROGATES})  # fmt: skip


def test_v7r2_rehearsal_loader_layout(cov: ModuleType, v7r: ModuleType, tmp_path: Path) -> None:
    """The V7-R2 rehearsal layout (PN-V7R2-3): 4 shards x 225 x 20 blocks, shards only; the V7-R layout still loads
    through the unchanged default (``test_v7r_coverage.py`` covers it)."""
    s5c = cov._load("steps_v5c")
    rs = cov._load("run_suite")
    digest = "d" * 64
    rng = np.random.default_rng(4)
    prov: dict[str, Any] = {
        "host_run_ids": ["RUN-X"],
        "hashed_source_digest": digest,
        "seed_base": 20505000,
        "code_sha": "c0de" * 10,
    }
    digests = {}
    for k in range(4):
        sums = rng.gamma(5.0, 100.0, (225, 20))
        side = s5c.write_v7r_sidecar(tmp_path / f"v7r2-s{k}-escaped_neutral.npz", sums)
        per = sums / 500.0
        doc = {"row": "v7-r2", "kind": "shard", "shard": k, "seed": 20505000 + 16000 + k, "seed_base": 20505000, "n": 2_250_000,
               "replicates": 225, "batch_histories": 500, "blocks_per_replicate": 20, "valid": True, "reduced": False,
               "git_sha": prov["code_sha"], "sidecar": side,
               "estimators": {"escaped_neutral": {"mean": per.mean(axis=1).tolist(), "sem": (per.std(axis=1, ddof=1) / 20**0.5).tolist()}}}  # fmt: skip
        sealed = cov._seal(doc)
        (tmp_path / f"v7r2-s{k}.json").write_text(json.dumps(sealed))
        digests[f"v7r2-s{k}.json"] = sealed["content_sha256"]
    prov["partial_digests"] = digests
    (tmp_path / "PROVENANCE.json").write_text(json.dumps(prov))
    blocks, info = cov.load_rehearsal(tmp_path, s5c, rs, digest=digest, layout="v7r2")
    assert (
        blocks.shape == (900, 20)
        and blocks.size == 18_000
        and info["layout"] == "v7r2"
        and info["shards"] == 4
    )
    assert set(info["sidecar_sha256"]) == set(digests) and info["partial_digests"] == digests
    with pytest.raises(cov.RehearsalError, match="hashed_source_digest"):
        cov.load_rehearsal(tmp_path, s5c, rs, digest="e" * 64, layout="v7r2")
    with pytest.raises(cov.RehearsalError, match="need >= 2|need exactly|no v7r-s"):
        cov.load_rehearsal(
            tmp_path, s5c, rs, digest=digest
        )  # the V7-R layout sees no v7r-s partials
    (tmp_path / "v7r2-ref.json").write_text("{}")
    with pytest.raises(cov.RehearsalError, match="shards only"):
        cov.load_rehearsal(tmp_path, s5c, rs, digest=digest, layout="v7r2")


def test_foreign_sha_kappa_record_is_rejected_by_parts_and_consumers(
    cov: ModuleType,
    v7r: ModuleType,
    s5c: ModuleType,
    tiny: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _make_kappa_dir(cov, v7r, s5c, tiny, tmp_path, monkeypatch)
    assert load_kappa_doc(cov, tmp_path)["run_sha"] == "unknown"
    a = SimpleNamespace(dirs=[str(tmp_path)], out_dir=str(tmp_path))
    assert s5c.load_kappa_record(a)["suite"] == "unknown"
    monkeypatch.setenv(
        "IONMC_RUN_SHA", "f" * 40
    )  # the current run is another SHA than the record's
    with pytest.raises(pytest.fail.Exception, match="another run_sha"):
        load_kappa_doc(cov, tmp_path)
    with pytest.raises(SystemExit, match="run_sha"):
        s5c.load_kappa_record(a)
    with pytest.raises(SystemExit, match="run_sha"):
        s5c.step_v7r2_calibration_combine(a)
    monkeypatch.delenv("IONMC_RUN_SHA")
    monkeypatch.setenv("IONMC_RUN_SUITE", "lv5x")
    with pytest.raises(SystemExit, match="suite"):
        s5c.load_kappa_record(a)


def test_calibration_layout_is_bound_to_the_design(v7r: ModuleType) -> None:
    assert v7r.r2_calibration_layout(48) == {
        "shards": 48,
        "replicates": 10800,
        "pairs": 5400,
        "parts": 39,
    }
    fb = v7r.r2_calibration_layout(
        64
    )  # the 7200-pair fallback has no calibration producer: parts is None
    assert fb["replicates"] == 14400 and fb["pairs"] == 7200 and fb["parts"] is None
    assert v7r.r2_part_layout() == {"replicates": 10800, "reps_per_shard": 225, "parts": 39}


def test_calibration_combine_rejects_parts_of_another_layout(v7r: ModuleType) -> None:
    kw: dict[str, Any] = {
        "parts": 39,
        "seed_fn": v7r.r2_cal_part_seeds,
        "layout": v7r.r2_part_layout(),
    }
    recs = []
    for si, s in enumerate(v7r.CAL_SURROGATES):
        for p in range(39):
            n = v7r.cal_part_sizes(3200, 39)[p]
            m = {
                k: {"kappa": 1.0, "accepted_joint": 0, "accepted_paired": 0, "accepted_single": 0}
                for k in v7r.CAL_MARGINS
            }
            m["nominal"].update(accepted_joint=n, accepted_paired=n, accepted_single=n)
            recs.append({"surrogate": s, "part": p, "n": n, "margins": m, "fit": {}, "reps_env": 3200,
                         "seeds": v7r.r2_cal_part_seeds(si, p), "layout": v7r.r2_part_layout()})  # fmt: skip
    assert v7r.calibration_combine(recs, **kw)["pass"] is True
    recs[4]["layout"] = {**v7r.r2_part_layout(), "replicates": 14400}
    with pytest.raises(SystemExit, match="layout"):
        v7r.calibration_combine(recs, **kw)
    recs[4].pop("layout")
    with pytest.raises(SystemExit, match="layout"):
        v7r.calibration_combine(recs, **kw)
