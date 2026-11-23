"""V3-013 R1-ion local-validation script (``validation/scripts/stopping/lv13_r1ion.py``)."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "validation/scripts/stopping/lv13_r1ion.py"


def _load():
    spec = importlib.util.spec_from_file_location("lv13_r1ion", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_r1ion_comparison_python_tiny() -> None:
    mod = _load()
    doc = mod.compare_selections("python", 20)
    assert doc["pass"], doc["comparisons"]
    assert doc["seed"] == 20612004 + 1100
    assert set(doc["comparisons"]) == {"A1-vs-A2", "B1-vs-B2", "A1-vs-B1", "A2-vs-B2"}
    assert doc["identity"]["A1"]["selected_by"] != doc["identity"]["B1"]["selected_by"]
    assert all(v["model_id_effective"] == "bethe-ionmc-v1" for v in doc["identity"].values())


def test_r1ion_seeds() -> None:
    mod = _load()
    assert mod.seed_for("python") == 20613104
    assert mod.seed_for("warp-cpu") == 20613204


@pytest.mark.host
@pytest.mark.parametrize("backend", ["python", "warp-cpu"])
def test_r1ion_full_script_host(backend: str, tmp_path: Path) -> None:
    out = subprocess.run(
        [sys.executable, str(SCRIPT), "--backend", backend, "--out-dir", str(tmp_path)],
        capture_output=True,
        text=True,
        env=dict(os.environ),
        check=False,
    )
    doc = json.loads(out.stdout.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])
    assert out.returncode == 0 and doc["pass"], doc["comparisons"]
