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


# --- row I2 (He table reproduction) on synthetic ASTAR / ICRU 90 inputs ----------------------


def _synthetic_alpha(s_factor: float = 1.0, r_factor: float = 1.0, icru_factor: float = 1.0):  # type: ignore[no-untyped-def]
    """ASTAR-like table from Bethe 75 eV and ICRU-90-like alpha arrays from Bethe 78 eV."""
    from ionmc.data.nist_star import StarTable
    from ionmc.materials import water
    from ionmc.physics.projectiles import ALPHA
    from ionmc.physics.stopping import BetheStoppingSource

    t75 = BetheStoppingSource().table(water(75.0), ALPHA)
    t78 = BetheStoppingSource().table(water(78.0), ALPHA)
    e_u = np.geomspace(1.0, 250.0, 122)
    s = t75.stopping_at(e_u) * s_factor
    r = t75.range_at(e_u) * r_factor
    ones = np.ones_like(e_u)
    star = StarTable("ASTAR", "WATER, LIQUID", 4.0 * e_u, s, 1e-3 * s, 1.001 * s, r, r, ones)
    e_icru = np.geomspace(1.0, 250.0, 49)
    return star, 4.0 * e_icru, t78.stopping_at(e_icru) * icru_factor


def test_i2_synthetic_reference_passes() -> None:
    mod = _load()
    out = mod.i2_block(*_synthetic_alpha())
    assert out["pass"] is True and all(out["gates"].values()) and len(out["gates"]) == 6
    for ref in ("vs_astar_75eV", "vs_icru90_78eV"):
        low = out["stopping"][ref][mod.I2_BAND_LOW]
        high = out["stopping"][ref][mod.I2_BAND_HIGH]
        assert low["energy_per_u_range_mev"][0] >= 2.0 and low["energy_per_u_range_mev"][1] < 10
        assert high["energy_per_u_range_mev"][0] >= 10.0
        assert high["energy_per_u_range_mev"][1] <= 250.0 + 1e-9
        assert low["n_points"] > 0 and high["n_points"] > 0
    rng = out["csda_range"]
    assert rng["energies_per_u_mev"] == [50.0, 100.0, 150.0, 200.0, 250.0]
    assert rng["vs_astar_75eV"]["max_abs_relative_deviation"] < 1e-4
    assert rng["vs_icru90_integral_78eV"]["max_abs_relative_deviation"] < 1e-3
    assert rng["vs_astar_75eV"]["energy_per_u_mev_at_max"] in rng["energies_per_u_mev"]
    shift = rng["bethe_78eV_minus_75eV_range_shift_mm_by_energy_per_u_mev"]
    assert sorted(shift, key=float) == ["50", "100", "150", "200", "250"]
    assert all(v > 0.0 for v in shift.values())  # higher I: lower S, longer range
    assert "75 eV" in out["i_value_statement"] and "78 eV" in out["i_value_statement"]


def test_i2_gates_respond_to_the_right_perturbation() -> None:
    mod = _load()
    # ASTAR S scaled by 1/1.005: |dev| ~ 0.5 % -> >= 10 MeV/u gate (0.30 %) fails, 2-10 passes;
    # the ASTAR range column is unchanged, so its range gate passes.
    g = mod.i2_block(*_synthetic_alpha(s_factor=1.0 / 1.005))["gates"]
    assert g["stopping_vs_astar_75eV_ge_10_MeV_u_max_abs_le_0.003"] is False
    assert g["stopping_vs_astar_75eV_2_to_10_MeV_u_max_abs_le_0.01"] is True
    assert g["csda_range_vs_astar_75eV_max_abs_le_0.002"] is True
    assert g["stopping_vs_icru90_78eV_ge_10_MeV_u_max_abs_le_0.003"] is True
    # ASTAR range scaled by 1.003 fails only the ASTAR range gate (0.20 %)
    out = mod.i2_block(*_synthetic_alpha(r_factor=1.003))
    failed = [k for k, v in out["gates"].items() if not v]
    assert failed == ["csda_range_vs_astar_75eV_max_abs_le_0.002"] and out["pass"] is False
    # ICRU 90 S scaled by 1.004 fails its >= 10 MeV/u S gate and the integral range gate (0.25 %)
    g = mod.i2_block(*_synthetic_alpha(icru_factor=1.004))["gates"]
    assert g["stopping_vs_icru90_78eV_ge_10_MeV_u_max_abs_le_0.003"] is False
    assert g["csda_range_vs_icru90_integral_78eV_max_abs_le_0.0025"] is False
    assert g["stopping_vs_astar_75eV_ge_10_MeV_u_max_abs_le_0.003"] is True


def test_i2_gates_fail_closed_on_empty_band() -> None:
    mod = _load()
    empty = mod.band_aggregate(np.array([20.0]), np.array([0.0]), 2.0, 10.0, include_hi=False)
    assert empty["n_points"] == 0 and empty["max_abs_relative_deviation"] is None
    ok = {"max_abs_relative_deviation": 0.0}
    stopping = {
        ref: {mod.I2_BAND_LOW: empty, mod.I2_BAND_HIGH: ok}
        for ref in ("vs_astar_75eV", "vs_icru90_78eV")
    }
    gates = mod.i2_gates(stopping, {ref: ok for ref in mod.I2_RANGE_MAX_ABS})
    assert gates["pass"] is False
    assert gates["criteria"]["stopping_vs_astar_75eV_2_to_10_MeV_u_max_abs_le_0.01"] is False
    assert mod.i2_gates({}, {})["pass"] is False


def test_i2_band_boundaries() -> None:
    mod = _load()
    e = np.array([2.0, 9.99, 10.0, 500.0])
    dev = np.array([0.001, 0.002, 0.003, 0.004])
    low = mod.band_aggregate(e, dev, 2.0, 10.0, include_hi=False)
    high = mod.band_aggregate(e, dev, 10.0, 500.0, include_hi=True)
    assert low["n_points"] == 2 and low["max_abs_relative_deviation"] == 0.002
    assert high["n_points"] == 2 and high["energy_per_u_range_mev"] == [10.0, 500.0]
