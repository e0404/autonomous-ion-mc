"""V5-SS ionmc thin-slab sampler ``steps_v8.py ss-sample``: seeds, argument checks and a tiny
warp-cpu run (800 histories, CI-sized; ``reduced``, never evidence)."""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "validation" / "scripts" / "transport"


def _load(name: str) -> ModuleType:
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


V8 = _load("steps_v8")


def test_constants_and_seeds() -> None:
    assert V8.SS_SEED_BASE == 20631004
    assert V8.SS_CASES == (22, 25, 26, 27, 28, 29, 30, 31)
    assert V8.SS_R_INDEX == {c: i for i, c in enumerate(V8.SS_CASES)}
    assert V8.SS_HISTORIES == 400_000
    assert [V8.ss_seed(c) for c in V8.SS_CASES] == [20631004 + 100 * i for i in range(8)]
    assert V8.ss_seed(31, 3) == 20631707
    assert "ss-sample" in V8.STEPS
    assert not set(V8.ROW_R_INDEX) & {"v5ss", "ss-sample"}  # different base, not an lv8 row
    with pytest.raises(ValueError):
        V8.ss_seed(23)
    with pytest.raises(ValueError):
        V8.ss_seed(22, 100)


def test_configuration_is_the_td_tl_default_physics() -> None:
    cfg = V8.ss_config(29)
    ref = V8.sample_config(6)  # C 430 TD/TL sample
    assert cfg.physics == ref.physics
    assert cfg.run.backend == "warp-cpu" and cfg.run.precision == "float64"
    assert cfg.run.n_batches == 40 and cfg.run.n_histories == 400_000
    assert cfg.run.seed == 20631004 + 500
    assert cfg.diagnostics.escape_records is True
    assert tuple(cfg.source.position_mm) == (0.0, 0.0, -0.001)


@pytest.mark.parametrize("argv", [
    ["ss-sample", "--case", "23"], ["ss-sample", "--case", "1"],
    ["ss-sample", "--case", "22", "--k", "100"], ["ss-sample", "--case", "22", "--scale", "0"],
])  # fmt: skip
def test_argument_checks(argv: list[str]) -> None:
    with pytest.raises((SystemExit, ValueError)):
        V8.main(argv)


def test_pooled_rms() -> None:
    tx, ty = np.array([1.0, -3.0]), np.array([2.0, 10.0])
    assert V8.pooled_rms(tx, ty) == (pytest.approx(np.sqrt((1 + 9 + 4 + 100) / 4)), 4)
    assert V8.pooled_rms(tx, ty, 3.0) == (pytest.approx(np.sqrt(14 / 3)), 3)


def test_tiny_run_writes_npz_and_document(tmp_path: Path) -> None:
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = V8.main(["ss-sample", "--case", "28", "--scale", "0.002", "--timeout", "100",
                      "--out-dir", str(tmp_path)])  # fmt: skip
    doc = json.loads(buf.getvalue().split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])
    assert rc == 0 and doc["pass"] is True and doc["valid"] is True, doc.get("problems")
    assert doc["row"] == "v5ss" and doc["case"] == 28 and doc["r_index"] == 4 and doc["k"] == 0
    assert doc["seed"] == 20631404 and doc["seed_base"] == 20631004
    assert doc["histories"] == 800 and doc["frozen_histories"] == 400_000 and doc["reduced"]
    assert doc["n_escaped"] == 800 and not any(doc["counters"].values())
    assert doc["effective_config"]["cpu_workers"] == 1
    assert doc["backend"] == "warp-cpu" and doc["precision"] == "float64"
    npz = tmp_path / "ionmc-ss-case28-k0.npz"
    assert doc["npz"] == npz.name
    assert doc["npz_sha256"] == hashlib.sha256(npz.read_bytes()).hexdigest()
    with np.load(npz) as z:
        tx, ty = z["theta_x"], z["theta_y"]
        assert tx.dtype == np.float64 and tx.shape == (800,)
        assert np.all(np.abs(z["z_exit_mm"] - 1.0) <= 1e-9)
        assert np.array_equal(z["escape_history"], np.arange(800))
        anchor = float(z["theta_rms_untruncated"])
    assert anchor == pytest.approx(doc["theta_anchor_rad"])
    assert anchor == pytest.approx(V8.pooled_rms(tx, ty)[0])
    assert 1e-4 < anchor < 1e-2  # He 250 MeV/u in 1 mm water: ~1 mrad
    assert doc["theta_w_rad"] == pytest.approx(V8.pooled_rms(tx, ty, 2.5 * anchor)[0])
    assert doc["theta_w_rad"] < anchor
