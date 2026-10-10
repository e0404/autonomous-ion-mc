"""Unit tests of validation/scripts/transport/nuclear_d9.py on synthetic data (no cache)."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

import numpy as np


def _mod() -> Any:
    path = Path(__file__).resolve().parents[2] / "validation/scripts/transport/nuclear_d9.py"
    spec = importlib.util.spec_from_file_location("nuclear_d9", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_window_stats_known_values() -> None:
    m = _mod()
    s = m.window_stats(np.array([1.0, 1.2]), np.array([0.1, 0.1]))
    assert abs(s["ratio"] - 1.1) < 1e-12
    assert s["ndf"] == 1 and abs(s["chi2_scatter"] - 2.0) < 1e-12
    assert abs(s["pdg_scale"] - np.sqrt(2.0)) < 1e-12
    assert abs(s["chi2_vs_unity"] - 4.0) < 1e-9
    one = m.window_stats(np.array([0.9]), np.array([0.05]))
    assert one["chi2_per_ndf"] is None and one["pdg_scale"] == 1.0


def test_windows_half_open_and_last_closed() -> None:
    m = _mod()
    e = np.array([19.9, 20.0, 39.99, 40.0, 250.0, 250.1])
    got = [m.in_window(e, lo, hi).tolist() for lo, hi in m.WINDOWS]
    assert got[0] == [False, True, True, False, False, False]
    assert got[1][3] and not got[0][3]
    assert got[4] == [False, False, False, False, True, False]


def test_run_d9_recovers_injected_ratios_and_rule() -> None:
    m = _mod()
    e = np.array([25.0, 50.0, 80.0, 120.0, 200.0])
    truth = 100.0 + e
    data = {"exfor-d0356": {"C-12": np.column_stack((e, truth, 0.02 * truth))}}
    models = {
        "la150": lambda t: lambda x: 0.85 * (100.0 + x),
        "tripathi": lambda t: lambda x: 1.00 * (100.0 + x),
        "bad": lambda t: lambda x: 1.30 * (100.0 + x),
    }
    res = m.run_d9(data, models)
    rec = res["sets"]["exfor-d0356"]["targets"]["C-12"]
    for k in ("20-40", "40-70", "70-110", "110-160", "160-250"):
        assert abs(rec["models"]["la150"][k]["ratio"] - 0.85) < 1e-12
        assert abs(rec["models"]["tripathi"][k]["ratio"] - 1.0) < 1e-12
    assert res["sets"]["exfor-d0356"]["targets"]["O-16"] == {"data": "none"}
    dec = res["decision"]["exfor-d0356"]["C-12"]
    assert dec["tripathi"]["preferred"] is True and dec["tripathi"]["n_windows"] == 5
    assert res["decision"]["exfor-d0356"]["N-14"] == {"status": "no data"}
    assert res["bgg"]["evaluable"] is False
    assert "URL" in m.markdown(res) and "not evaluable" in m.markdown(res)


def test_preferred_needs_every_window_and_a_margin() -> None:
    m = _mod()

    def w(r: float, err: float = 0.02) -> dict[str, Any]:
        return {"ratio": r, "err": err}

    la = {"a": w(0.85), "b": w(0.90), "c": None}
    assert m.preferred(la, {"a": w(1.0), "b": w(1.0), "c": w(1.0)})["preferred"]
    assert not m.preferred(la, {"a": w(1.0), "b": w(0.89)})["preferred"]  # worse in one window
    assert not m.preferred(la, {"a": w(0.86), "b": w(0.91)})["preferred"]  # inside 1 sigma_comb
    assert not m.preferred({"a": None}, {"a": None})["preferred"]  # no data -> not preferred


def test_roles_are_kept_separate() -> None:
    m = _mod()
    assert m.SETS["exfor-d0356"].startswith("evaluation")
    assert m.SETS["exfor-c1862"].startswith("report-only")
    assert m.EVALUATION_SET == "exfor-d0356"
