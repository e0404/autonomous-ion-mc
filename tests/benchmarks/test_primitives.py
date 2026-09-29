import importlib.util
import json
from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[2] / "benchmarks/calibration/primitives.py"
spec = importlib.util.spec_from_file_location("calibration_primitives", MODULE)
primitives = importlib.util.module_from_spec(spec)
spec.loader.exec_module(primitives)

pytest.importorskip("warp")


def test_cpu_step_kernel_report_has_expected_shape(tmp_path):
    out = tmp_path / "primitives.json"
    assert (
        primitives.main(
            [
                "--devices",
                "cpu",
                "--particles",
                "256",
                "--steps",
                "8",
                "--dtypes",
                "float32",
                "float64",
                "--repeats",
                "2",
                "--output",
                str(out),
            ]
        )
        == 0
    )
    report = json.loads(out.read_text())
    assert report["schema_version"] == 1
    assert "cpu" in report["devices"]
    steps = report["measurements"]["step_kernel"]
    assert len(steps) == 4  # 2 dtypes x deposit on/off
    for entry in steps:
        assert entry["best_seconds"] > 0
        assert entry["particle_steps_per_second_upper"] > 0
        assert entry["mean_final_energy"] < 150.0
        if entry["deposit"]:
            assert entry["grid_total_deposit"] > 0
        else:
            assert entry["grid_total_deposit"] is None
    assert report["measurements"]["host_peak_rss_kib"] > 0


def test_float32_and_float64_agree_on_energy_loss(tmp_path):
    out = tmp_path / "p.json"
    primitives.main(
        ["--devices", "cpu", "--particles", "512", "--steps", "20", "--repeats", "1",
         "--output", str(out)]
    )
    report = json.loads(out.read_text())
    by_dtype = {
        (e["dtype"], e["deposit"]): e for e in report["measurements"]["step_kernel"]
    }
    f32 = by_dtype[("float32", True)]["mean_final_energy"]
    f64 = by_dtype[("float64", True)]["mean_final_energy"]
    assert abs(f32 - f64) / f64 < 0.05
