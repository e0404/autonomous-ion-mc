"""Tests of depth-dose metrics and run loading on synthetic data."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from ionmc.reference import RunError, depth_dose, load_run
from ionmc.reference.metrics import (
    distal_falloff_80_20,
    integral_depth_dose,
    normalize_to_peak,
    peak_depth,
    r80,
    r90,
)
from ionmc.reference.runs import file_hashes


def _curve() -> tuple[np.ndarray, np.ndarray]:
    d = np.arange(0.5, 200.0, 1.0)
    c = np.where(d <= 150.5, 0.3 + 0.7 * d / 150.5, np.clip(1 - (d - 150.5) / 10.0, 0, None))
    return d, c


def test_metrics_analytic() -> None:
    d, c = _curve()
    assert peak_depth(d, c) == 150.5
    # distal ramp 1 - (d-150.5)/10: 90 % at 151.5, 80 % at 152.5, 20 % at 158.5
    assert r90(d, c) == pytest.approx(151.5, abs=1e-9)
    assert r80(d, c) == pytest.approx(152.5, abs=1e-9)
    assert distal_falloff_80_20(d, c) == pytest.approx(6.0, abs=1e-9)
    n = normalize_to_peak(c * 7)
    assert n.max() == 1.0


def test_metrics_errors() -> None:
    d, c = _curve()
    with pytest.raises(ValueError):
        r80(d, np.ones_like(d))  # never falls
    with pytest.raises(ValueError):
        peak_depth(d[::-1], c)
    with pytest.raises(ValueError):
        normalize_to_peak(np.zeros(4))


def test_integral_depth_dose() -> None:
    a = np.arange(24, dtype=float).reshape(2, 3, 4)
    np.testing.assert_allclose(integral_depth_dose(a, 2), a.sum(axis=(0, 1)))
    np.testing.assert_allclose(integral_depth_dose(a, 0, 2.0), 2 * a.sum(axis=(1, 2)))


def _run(tmp: Path, engine: str, files: dict[str, bytes], skip: tuple[str, ...] = ()) -> Path:
    """Write fixture files plus a manifest with their real size and sha256."""
    case = json.dumps({"histories": 10, "input": "x"}).encode()
    files = {**files, "inputs/case.json": case}
    for rel, blob in files.items():
        p = tmp / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(blob)
    manifest = {
        "engine": engine,
        "run_id": "REF-test",
        "source_sha": "abc",
        "files": {
            rel: {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
            for rel, b in files.items()
            if rel not in skip
        },
    }
    (tmp / "transfer-manifest.json").write_text(json.dumps(manifest))
    return tmp


def test_topas_run(tmp_path: Path) -> None:
    csv = (
        "# X in 1 bin of 12 cm\n# Y in 1 bin of 12 cm\n# Z in 3 bins of 0.1 cm\n"
        "# DoseToMedium ( Gy ) : Sum Standard_Deviation\n"
        "0,0,0,1,0.1\n0,0,1,3,0.1\n0,0,2,2,0.1\n"
    )
    run = load_run(_run(tmp_path, "topas", {"work/idd_dose.csv": csv.encode()}))
    dd = depth_dose(run)
    np.testing.assert_allclose(dd.depth_mm, [0.5, 1.5, 2.5])
    np.testing.assert_allclose(dd.dose, [1, 3, 2])
    assert dd.histories == 10


def _mhd(
    dims: tuple[int, int, int], sp: float, offset: str, local: bool, arr: np.ndarray
) -> tuple[bytes, bytes]:
    head = (
        f"ObjectType = Image\nNDims = 3\nDimSize = {dims[0]} {dims[1]} {dims[2]}\n"
        f"Offset = {offset}\nElementSpacing = {sp} {sp} {sp}\nElementType = MET_FLOAT\n"
        f"ElementByteOrderMSB = False\nElementDataFile = {'LOCAL' if local else 'Dose.raw'}\n"
    )
    return head.encode(), arr.astype("<f4").tobytes()


def test_fred_and_mcsquare_runs(tmp_path: Path) -> None:
    arr = np.ones((4, 2, 2))
    arr[:, :, :] *= np.arange(4)[:, None, None]  # z-dependent
    h, b = _mhd((2, 2, 4), 1.0, "-0.5 -0.5 0.5", True, arr)
    fred = load_run(_run(tmp_path / "f", "fred", {"work/out/score/Phantom.Dose.mhd": h + b}))
    dd = depth_dose(fred)
    np.testing.assert_allclose(dd.depth_mm, [0.5, 1.5, 2.5, 3.5])
    np.testing.assert_allclose(dd.dose, [0, 4, 8, 12])

    arr2 = np.ones((2, 4, 2)) * np.arange(4)[None, :, None]  # (nz, ny, nx), y-dependent
    h, b = _mhd((2, 4, 2), 2.0, "0 0 0", False, arr2)
    mc = load_run(
        _run(
            tmp_path / "m",
            "mcsquare",
            {"work/Outputs/Dose.mhd": h, "work/Outputs/Dose.raw": b},
        )
    )
    dd = depth_dose(mc)
    np.testing.assert_allclose(dd.depth_mm, [1, 3, 5, 7])
    np.testing.assert_allclose(dd.dose, [12, 8, 4, 0])  # reversed: depth = y_max - y


def test_run_fail_closed(tmp_path: Path) -> None:
    with pytest.raises(RunError):
        load_run(tmp_path)  # no manifest
    run = load_run(_run(tmp_path, "topas", {"work/other.csv": b"x"}))
    with pytest.raises(RunError):
        depth_dose(run)


_CSV = (
    "# X in 1 bin of 12 cm\n# Y in 1 bin of 12 cm\n# Z in 3 bins of 0.1 cm\n"
    "# DoseToMedium ( Gy ) : Sum Standard_Deviation\n0,0,0,1,0.1\n0,0,1,3,0.1\n0,0,2,2,0.1\n"
)


def test_hashes_reported_for_bytes_read(tmp_path: Path) -> None:
    run = load_run(_run(tmp_path, "topas", {"work/idd_dose.csv": _CSV.encode()}))
    dd = depth_dose(run)
    hashes = file_hashes(run, [*dd.files, "inputs/case.json"])
    assert hashes["work/idd_dose.csv"] == hashlib.sha256(_CSV.encode()).hexdigest()
    with pytest.raises(RunError):
        file_hashes(run, ["work/never-read.csv"])


def test_corrupted_file_rejected(tmp_path: Path) -> None:
    d = _run(tmp_path, "topas", {"work/idd_dose.csv": _CSV.encode()})
    run = load_run(d)
    blob = bytearray((d / "work/idd_dose.csv").read_bytes())
    blob[-3] = ord("9")  # same size, modified byte
    (d / "work/idd_dose.csv").write_bytes(bytes(blob))
    with pytest.raises(RunError, match="sha256"):
        depth_dose(run)
    (d / "work/idd_dose.csv").write_bytes(_CSV.encode() + b"\n")
    with pytest.raises(RunError, match="size"):
        depth_dose(run)


def test_corrupted_case_json_rejected(tmp_path: Path) -> None:
    d = _run(tmp_path, "topas", {"work/idd_dose.csv": _CSV.encode()})
    (d / "inputs/case.json").write_text('{"histories": 11, "input": "x"}')
    with pytest.raises(RunError):
        load_run(d)


def test_missing_manifest_entry_rejected(tmp_path: Path) -> None:
    d = _run(tmp_path, "topas", {"work/idd_dose.csv": _CSV.encode()}, skip=("work/idd_dose.csv",))
    with pytest.raises(RunError, match="manifest"):
        depth_dose(load_run(d))
    d2 = _run(tmp_path / "x", "topas", {}, skip=("inputs/case.json",))
    with pytest.raises(RunError, match="manifest"):
        load_run(d2)


def test_mcsquare_raw_must_be_manifested(tmp_path: Path) -> None:
    arr = np.ones((2, 4, 2))
    h, b = _mhd((2, 4, 2), 2.0, "0 0 0", False, arr)
    d = _run(
        tmp_path,
        "mcsquare",
        {"work/Outputs/Dose.mhd": h, "work/Outputs/Dose.raw": b},
        skip=("work/Outputs/Dose.raw",),
    )
    with pytest.raises(RunError, match="manifest"):
        depth_dose(load_run(d))
    h2 = h.replace(b"Dose.raw", b"Other.raw")
    d2 = _run(
        tmp_path / "o",
        "mcsquare",
        {"work/Outputs/Dose.mhd": h2, "work/Outputs/Other.raw": b},
    )
    with pytest.raises(RunError, match="Dose.raw"):
        depth_dose(load_run(d2))
