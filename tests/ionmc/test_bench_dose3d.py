"""Reference-path regression tests for the 3-D dose benchmark driver (decision
0031). These guard the driver's own logic — per-history normalisation, digest
determinism, and dose-grid independence between runs — on the reference path, so
they run in CI without a GPU. The harness itself is tested in test_benchmarking.py.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_DRIVER = Path(__file__).resolve().parents[2] / "benchmarks" / "bench_dose3d.py"


def _load_driver():
    spec = importlib.util.spec_from_file_location("bench_dose3d", _DRIVER)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def driver():
    return _load_driver()


@pytest.fixture(scope="module")
def built(driver, pstar_cache_root):
    return driver._build(pstar_cache_root)


def test_per_history_dose_is_deterministic_and_conserving(driver, built) -> None:
    """The deterministic (scattering/straggling off) 3-D dose is per-history exact:
    the per-history integral is the beam energy (150 MeV), and two runs agree bit for
    bit (identical digest) — the property the cross-backend gate relies on."""
    eng, src, lat = built
    d1 = driver._run(eng, src, lat, 1, "python", "cpu")
    d2 = driver._run(eng, src, lat, 1, "python", "cpu")
    assert float(d1.sum()) == pytest.approx(150.0, rel=1e-9)
    assert driver.array_digest(d1.ravel()) == driver.array_digest(d2.ravel())


def test_dose_grids_do_not_leak_between_runs(driver, built) -> None:
    """Each run allocates a fresh dose grid, so per-history dose is independent of
    the history count (no accumulation leak): dose(n)/n == dose(1)/1."""
    eng, src, lat = built
    d1 = driver._run(eng, src, lat, 1, "python", "cpu")
    d5 = driver._run(eng, src, lat, 5, "python", "cpu")
    assert np.allclose(d1, d5, rtol=1e-12, atol=0.0)


def test_digest_structure(driver, built) -> None:
    eng, src, lat = built
    dig = driver._digest(driver._run(eng, src, lat, 1, "python", "cpu"))
    assert dig["integral_per_history_mev"] == pytest.approx(150.0, rel=1e-9)
    assert dig["peak_voxel_mev"] > 0.0
    assert dig["n_nonzero_voxels"] > 0
    assert isinstance(dig["shape_digest"], str)


def test_main_reference_only_passes_and_is_well_formed(
    driver, pstar_cache_root, monkeypatch, capsys
) -> None:
    """End-to-end ``main()`` on the reference path (Warp forced absent, no
    --require-cuda): the physics gate passes, the return code is 0, and the emitted
    report is well-formed JSON with provenance, config and the reference digest. This
    covers the pass/fail decision logic that the unit helpers do not."""
    import json

    monkeypatch.setattr(driver.mathlib, "HAVE_WARP", False)
    argv = ["bench_dose3d"]
    if pstar_cache_root is not None:
        argv += ["--cache-dir", str(pstar_cache_root)]
    monkeypatch.setattr("sys.argv", argv)

    rc = driver.main()
    assert rc == 0

    report = json.loads(capsys.readouterr().out)
    assert report["physics_gate_passed"] is True
    assert report["benchmark"] == "dose3d"
    assert report["warp"] == {"available": False}
    assert report["config"]["reference_precision"] == "float64"
    assert report["digest"]["reference"]["integral_per_history_mev"] == pytest.approx(
        150.0, rel=1e-9
    )
    # no Warp backend -> no cross-backend comparisons and an empty scaling sweep
    assert report["cross_backend"] == {}
    assert report["scaling"] == []
