"""Test of the energy-window helper of validation/scripts/stopping/compare_nist.py."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[2] / "validation/scripts/stopping/compare_nist.py"


def _load():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location("compare_nist", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_window_selects_2_to_500_mev_per_u() -> None:
    mod = _load()
    grid = np.concatenate(([0.5, 1.0, 1.99], np.geomspace(2.0, 10000.0, 60), [500.0]))
    selected = grid[mod.window_mask(grid)]
    assert selected.size > 0
    assert selected.min() >= 2.0 and selected.max() <= 500.0
    assert 2.0 in selected and 500.0 in selected
    assert np.all(~mod.window_mask(np.array([1.99, 500.1, 10000.0])))


def test_aggregate_reports_actual_window() -> None:
    mod = _load()
    e = np.array([2.0, 10.0, 100.0, 500.0])
    out = mod._aggregate(e, np.array([0.01, -0.02, 0.0, 0.005]))
    assert out["energy_per_u_range_mev"] == [2.0, 500.0]
    assert out["n_points"] == 4 and out["n_points_ge_10_MeV_u"] == 3
