"""D-3 range attribution: reproduces the committed aggregate JSON; no raw stopping values."""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "validation/scripts/ion_em/d3_range_attribution.py"
COMMITTED = ROOT / "validation/results/stopping/v3-008-d3-range-attribution.json"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("d3_range_attribution", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    sys.modules["d3_range_attribution"] = m
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def result(mod):
    try:
        mod.find_icru(ROOT / ".ionmc-cache")
    except FileNotFoundError:
        pytest.skip("registered ICRU 90 object not cached")
    return mod.compute(ROOT, ROOT / ".ionmc-cache")


def _same(a, b, path=""):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a:
            _same(a[k], b[k], f"{path}/{k}")
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _same(x, y, f"{path}[{i}]")
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9), path
    else:
        assert a == b, path


def test_reproduces_committed_json(result):
    _same(result, json.loads(COMMITTED.read_text()))


def test_no_raw_stopping_values(result, mod):
    text = json.dumps(result, indent=2, sort_keys=True)
    icru = mod.load_icru90_water(mod.find_icru(ROOT / ".ionmc-cache"))
    for arr in (icru.proton_stopping, icru.alpha_stopping, icru.proton_energy_mev):
        for v in arr:
            assert f"{v:.6g}" not in text.replace(" ", "") or len(f"{v:.6g}") < 5
    for row in result["cases"]:
        assert all(not isinstance(v, list | dict) for v in row.values())
        assert not any("stopping" in k and not k.endswith("_pct") for k in row)


def test_sanity(result):
    rows = {(r["ion"], r["mev_u"]): r for r in result["cases"]}
    assert len(rows) == 12
    he = rows[("he", 250.0)]
    assert 0.0 < he["icru90_alpha_direct_vs_ionmc_pct"] < 0.15
    for r in rows.values():
        assert 0.0 < r["prestep_euler_0p1mm_pct"] <= r["prestep_pure_euler_0p1mm_pct"] + 1e-12
    assert rows[("oxygen", 430.0)]["topas_r80_over_ionmc_r80_pct"] > 0.25
