"""Tests of depth-dose metrics and run loading on synthetic data."""

from __future__ import annotations

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


def _run(tmp: Path, engine: str, files: dict[str, bytes]) -> Path:
    for rel, blob in files.items():
        p = tmp / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(blob)
    manifest = {
        "engine": engine,
        "run_id": "REF-test",
        "source_sha": "abc",
        "files": {rel: {"bytes": len(b), "sha256": "0" * 64} for rel, b in files.items()},
    }
    (tmp / "transfer-manifest.json").write_text(json.dumps(manifest))
    (tmp / "inputs").mkdir(exist_ok=True)
    (tmp / "inputs" / "case.json").write_text(json.dumps({"histories": 10, "input": "x"}))
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
        load_run(tmp_path)
    run = load_run(_run(tmp_path, "topas", {"work/other.csv": b"x"}))
    with pytest.raises(RunError):
        depth_dose(run)
