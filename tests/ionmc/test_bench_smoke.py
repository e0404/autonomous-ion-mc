"""Smoke test of the EM-only benchmark harness: 1e3 histories, one repeat, JSON written."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def test_bench_em_smoke() -> None:
    out = REPO / "benchmarks" / "generated" / "transport" / "test-smoke.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.unlink(missing_ok=True)
    try:
        cmd = [sys.executable, str(REPO / "benchmarks/transport/bench_em.py"), "--out", str(out),
               "--histories", "1000", "--workers", "1", "--repeats", "1"]  # fmt: skip
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=REPO)
        assert r.returncode == 0, r.stderr[-2000:]
        doc = json.loads(out.read_text())
        assert doc["targets"] is None and "EM-only" in doc["label"]
        kinds = [x["kind"] for x in doc["results"]]
        assert kinds == ["cold_start", "repeat_1"]
        assert all(x["valid"] and x["label"] == "EM-only" for x in doc["results"])
        again = subprocess.run(cmd, capture_output=True, text=True, timeout=60, cwd=REPO)
        assert again.returncode != 0  # refuses to overwrite
    finally:
        out.unlink(missing_ok=True)
