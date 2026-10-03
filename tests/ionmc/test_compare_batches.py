"""Batch protocol of compare_batches.py on synthetic reference runs."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from ionmc.reference import RunError, dose_3d, load_run
from tests.ionmc.test_reference_lateral import gauss_binned

SCRIPT = Path(__file__).resolve().parents[2] / "validation/scripts/reference/compare_batches.py"
NX, NZ, H = 80, 40, 0.5


def _load():
    spec = importlib.util.spec_from_file_location("compare_batches", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _amp(z: np.ndarray, r: float) -> np.ndarray:
    return np.where(z <= r - 2, 1 + 0.01 * z, np.clip((r + 2 - z) / 4.0, 0, None) * (1 + 0.01 * r))


def _write(tmp: Path, engine: str, files: dict[str, bytes], case: dict, run_id: str) -> Path:
    files = {**files, "inputs/case.json": json.dumps(case).encode()}
    for rel, blob in files.items():
        p = tmp / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(blob)
    manifest = {
        "engine": engine,
        "run_id": run_id,
        "source_sha": "abc",
        "files": {
            rel: {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
            for rel, b in files.items()
        },
    }
    (tmp / "transfer-manifest.json").write_text(json.dumps(manifest))
    return tmp


def topas_run(
    tmp: Path,
    seed: int,
    *,
    range_mm: float = 30.0,
    s_scale: float = 1.0,
    histories: int = 100_000,
    native_seed: int | None = None,
    lateral: bool = True,
    idd_wrong: bool = False,
    primary: bool | str = False,
) -> Path:
    z1 = (np.arange(NZ) + 0.5) * 1.0
    amp = _amp(z1, range_mm)
    slabs = [
        np.outer(*(2 * [gauss_binned(NX, H, (1.5 + 0.05 * z) * s_scale)])) * a
        for z, a in zip(z1, amp, strict=True)
    ]
    d3 = np.stack(slabs, axis=2)
    lat = d3.mean(axis=(0, 1))  # IDD voxel spans the field: dose = lateral mean
    if idd_wrong:
        lat = lat[::-1]
    idd = np.repeat(lat, 2)
    rows = "".join(f"0,0,{i},{float(v)!r},0.0\n" for i, v in enumerate(idd))
    csv = (
        f"# X in 1 bin of 120 mm\n# Y in 1 bin of 120 mm\n# Z in {idd.size} bins of 0.5 mm\n"
        f"# DoseToMedium ( Gy ) : Sum Standard_Deviation\n{rows}"
    )
    hdr = (
        f"# X in {NX} bins of 0.5 mm\n# Y in {NX} bins of 0.5 mm\n# Z in {NZ} bins of 1 mm\n"
        "# DoseToMedium ( Gy ) : Sum\n"
    )
    files = {
        "work/idd_dose.csv": csv.encode(),
        "inputs/input.txt": f"i:Ts/Seed = {native_seed or seed}\n".encode(),
    }
    if lateral:
        files["work/dose3d.bin"] = d3.astype("<f8").ravel(order="F").tobytes()
        files["work/dose3d.binheader"] = hdr.encode()
    if lateral and primary:
        amp_p = 0.9 if primary is True else 1.5  # "exceeds" -> more dose than all particles
        d3p = d3 * amp_p * 0.8
        files["work/dose3d_primary.bin"] = d3p.astype("<f8").ravel(order="F").tobytes()
        files["work/dose3d_primary.binheader"] = hdr.encode()
    case = {"histories": histories, "input": "input.txt", "seeds": [seed]}
    return _write(tmp, "topas", files, case, f"REF-topas-{seed}")


def fred_run(tmp: Path, seed: int) -> Path:
    z = (np.arange(100) + 0.5) * 0.5
    arr = _amp(z, 25.0).reshape(100, 1, 1).astype("<f4")
    head = (
        b"NDims = 3\nDimSize = 1 1 100\nOffset = 0 0 0.25\nElementSpacing = 1 1 0.5\n"
        b"ElementType = MET_FLOAT\nElementByteOrderMSB = False\nElementDataFile = LOCAL\n"
    )
    case = {
        "histories": 100_000,
        "input": "fred.inp",
        "seeds": [seed],
        "arguments": ["-rseed", str(seed)],
    }
    return _write(
        tmp,
        "fred",
        {"work/out/score/Phantom.Dose.mhd": head + arr.tobytes(), "inputs/fred.inp": b"x\n"},
        case,
        f"REF-fred-{seed}",
    )


def test_batch_statistics_and_lateral(tmp_path: Path) -> None:
    mod = _load()
    dirs = [
        topas_run(tmp_path / f"t{i}", 100 + i, range_mm=30.0 + 0.2 * i, s_scale=1 + 0.02 * i)
        for i in range(3)
    ]
    doc = mod.batch_analysis(dirs)
    assert doc["evidence_status"] == "batched"
    st = doc["engines"]["topas"]
    assert st["r80_mm"]["n"] == 3 and st["r80_mm"]["sd"] > 0
    lo, hi = st["r80_mm"]["ci95"]
    assert lo < st["r80_mm"]["mean"] < hi
    assert st["r80_mm"]["se"] == pytest.approx(st["r80_mm"]["sd"] / 3**0.5)
    assert hi - st["r80_mm"]["mean"] == pytest.approx(4.3027 * st["r80_mm"]["se"])
    assert set(st["lateral"]) == {"full@0.5", "full@0.9", "w20@0.5", "w20@0.9"}
    s = st["lateral"]["w20@0.5"]["sigma_mm"]
    # sigma(z) = (1.5 + 0.05 z) * scale, z about 0.5 R80 (~15 mm): ~2.25 mm * (1..1.04)
    assert 2.2 < s["mean"] < 2.45 and s["se"] > 0
    assert mod.t95(1) == 12.7062 and mod.t95(30) == 2.0423 and mod.t95(99) == 1.96


def test_primary_scorer_reported_and_checked(tmp_path: Path) -> None:
    mod = _load()
    dirs = [
        topas_run(tmp_path / f"p{i}", 200 + i, s_scale=1 + 0.02 * i, primary=True) for i in range(2)
    ]
    doc = mod.batch_analysis(dirs)
    st = doc["engines"]["topas"]
    assert set(st["lateral_primary"]) == set(st["lateral"])
    # the synthetic primary dose has the same shape as the all-particle dose
    a, b = st["lateral"]["w20@0.9"], st["lateral_primary"]["w20@0.9"]
    assert b["sigma_mm"]["mean"] == pytest.approx(a["sigma_mm"]["mean"], rel=1e-6)
    d3p = dose_3d(load_run(dirs[0]), "dose3d_primary")
    assert d3p.dose.shape == (NX, NX, NZ)
    with pytest.raises(RunError, match="exceeds"):
        dose_3d(load_run(topas_run(tmp_path / "x", 5, primary="exceeds")), "dose3d_primary")
    with pytest.raises(mod.BatchError, match="3-D"):  # primary scored in only one run
        mod.batch_analysis([dirs[0], topas_run(tmp_path / "y", 6)])
    with pytest.raises(RunError, match="unknown"):
        dose_3d(load_run(dirs[0]), "other")


def test_single_run_fails(tmp_path: Path) -> None:
    mod = _load()
    with pytest.raises(mod.BatchError, match="at least 2"):
        mod.batch_analysis([topas_run(tmp_path / "a", 1)])
    with pytest.raises(SystemExit):
        mod.main(
            ["--runs", str(topas_run(tmp_path / "b", 2)), "--output", str(tmp_path / "o.json")]
        )
    assert not (tmp_path / "o.json").exists()


def test_duplicate_seeds_fail(tmp_path: Path) -> None:
    mod = _load()
    with pytest.raises(mod.BatchError, match="duplicate seeds"):
        mod.batch_analysis([topas_run(tmp_path / "a", 7), topas_run(tmp_path / "b", 7)])


def test_native_seed_must_match_case_json(tmp_path: Path) -> None:
    mod = _load()
    bad = topas_run(tmp_path / "a", 8, native_seed=9)
    with pytest.raises(mod.BatchError, match="native input seed"):
        mod.batch_analysis([bad, topas_run(tmp_path / "b", 10)])


def test_mismatched_runs_fail(tmp_path: Path) -> None:
    mod = _load()
    a = topas_run(tmp_path / "a", 1)
    # engine with a single run alongside a valid batch
    with pytest.raises(mod.BatchError, match="engine fred"):
        mod.batch_analysis([a, topas_run(tmp_path / "b", 2), fred_run(tmp_path / "f", 3)])
    # mixed histories
    with pytest.raises(mod.BatchError, match="histories"):
        mod.batch_analysis([a, topas_run(tmp_path / "c", 4, histories=200_000)])
    # too few histories
    with pytest.raises(mod.BatchError, match="histories"):
        mod.batch_analysis(
            [
                topas_run(tmp_path / "d", 5, histories=200),
                topas_run(tmp_path / "e", 6, histories=200),
            ]
        )
    # 3-D dose present in only some runs
    with pytest.raises(mod.BatchError, match="3-D"):
        mod.batch_analysis([a, topas_run(tmp_path / "g", 11, lateral=False)])
    # same run directory twice (distinct seeds impossible)
    with pytest.raises(mod.BatchError):
        mod.batch_analysis([a, a])


def test_dose3d_must_reproduce_idd(tmp_path: Path) -> None:
    ok = load_run(topas_run(tmp_path / "a", 1))
    assert dose_3d(ok).dose.shape == (NX, NX, NZ)
    with pytest.raises(RunError, match="IDD"):
        dose_3d(load_run(topas_run(tmp_path / "b", 2, idd_wrong=True)))
    with pytest.raises(RunError, match="engine"):
        dose_3d(load_run(fred_run(tmp_path / "f", 3)))
