"""Pencil-beam angular spread and x-theta correlation (V3-005C C8, Amendment 17 (c)9)."""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import replace
from typing import Any

import numpy as np
import pytest

from ionmc.config import DiagnosticsOptions, SimulationConfig
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport.reference import run_reference_range
from ionmc.transport.warp_driver import run_warp_range

MakeConfig = Callable[..., SimulationConfig]
SX, ST, RHO = 5.948, 31.55, -0.521


def _miss_config(make_config: MakeConfig, n: int, **src: Any) -> SimulationConfig:
    """A beam that misses the world: the end record is the sampled source state."""
    geo = BoxPhantom((900.0, 900.0, 0.0), (20.0, 20.0, 20.0), WATER)
    cfg = make_config(
        energy=100.0, n=n, n_batches=2, geometry=geo, position=(0.0, 0.0, -10.0),
        diagnostics=DiagnosticsOptions(track_end_positions=True), allow_invalid=True,
    )  # fmt: skip
    return replace(cfg, source=replace(cfg.source, **src))


def _ends(cfg: SimulationConfig, backend: str) -> tuple[np.ndarray, np.ndarray]:
    eff = Simulation(cfg).effective
    n = cfg.run.n_histories
    p = run_reference_range(eff, 0, n) if backend == "python" else run_warp_range(eff, 0, n, "cpu")
    d = p.diagnostics
    return np.asarray(d.end_position_mm), np.asarray(d.end_direction)


def test_validation_fails_closed() -> None:
    base: dict[str, Any] = {"projectile": PROTON, "position_mm": (0, 0, 0),
                            "direction": (0, 0, 1), "kinetic_energy_mev": 100.0}  # fmt: skip
    PencilBeamSource(**base, lateral_sigma_mm=1.0, angular_sigma_mrad=5.0, x_theta_correlation=-0.5)
    for bad in (
        {"angular_sigma_mrad": -1.0},
        {"angular_sigma_mrad": 101.0},
        {"angular_sigma_mrad": math.nan},
        {"angular_sigma_mrad": 5.0, "lateral_sigma_mm": 1.0, "x_theta_correlation": 1.0},
        {"angular_sigma_mrad": 5.0, "lateral_sigma_mm": 1.0, "x_theta_correlation": -1.0},
        {"angular_sigma_mrad": 0.0, "lateral_sigma_mm": 1.0, "x_theta_correlation": 0.3},
        {"angular_sigma_mrad": 5.0, "lateral_sigma_mm": 0.0, "x_theta_correlation": 0.3},
    ):
        with pytest.raises(UnsupportedCombinationError):
            PencilBeamSource(**{**base, **bad})


def test_defaults_are_unchanged_and_summary_has_no_new_keys(make_config: MakeConfig) -> None:
    cfg = make_config(lateral_sigma=0.5, energy_sigma=0.2)
    s = Simulation(cfg).effective.summary()["source"]
    assert "angular_sigma_mrad" not in s and "x_theta_correlation" not in s
    assert cfg.source.angular_sigma_mrad == 0.0 and cfg.source.x_theta_correlation == 0.0
    cfg2 = replace(cfg, source=replace(cfg.source, angular_sigma_mrad=3.0))
    assert Simulation(cfg2).effective.summary()["source"]["angular_sigma_mrad"] == 3.0


def test_default_draws_no_extra_number_bit_identical(make_config: MakeConfig) -> None:
    """With angular_sigma_mrad = 0 the sample equals a source built without the new fields."""
    a = _ends(_miss_config(make_config, 64, lateral_sigma_mm=SX, energy_sigma_mev=1.0), "python")
    b = _ends(
        _miss_config(make_config, 64, lateral_sigma_mm=SX, energy_sigma_mev=1.0,
                     angular_sigma_mrad=0.0, x_theta_correlation=0.0), "python")  # fmt: skip
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
    assert np.all(a[1] == np.array([0.0, 0.0, 1.0]))


def test_direction_convention_and_statistics(make_config: MakeConfig) -> None:
    n = 20000
    cfg = _miss_config(make_config, n, lateral_sigma_mm=SX, angular_sigma_mrad=ST,
                       x_theta_correlation=RHO)  # fmt: skip
    pos, d = _ends(cfg, "python")
    assert np.allclose(np.linalg.norm(d, axis=1), 1.0, atol=1e-14)
    tx, ty = d[:, 0] / d[:, 2], d[:, 1] / d[:, 2]  # (tan thx, tan thy, 1) normalised
    x, y = pos[:, 0], pos[:, 1]
    for xs, ts in ((x, tx), (y, ty)):
        assert abs(xs.std() / SX - 1) < 0.03
        assert abs(ts.std() / (ST * 1e-3) - 1) < 0.03
        assert abs(np.corrcoef(xs, ts)[0, 1] - RHO) < 0.02
    assert abs(np.corrcoef(x, ty)[0, 1]) < 0.03


def test_crn_variants_scale_the_same_deviates(make_config: MakeConfig) -> None:
    n = 256
    nom = _ends(_miss_config(make_config, n, lateral_sigma_mm=SX, angular_sigma_mrad=ST,
                             x_theta_correlation=RHO), "python")  # fmt: skip
    sx = _ends(_miss_config(make_config, n, lateral_sigma_mm=1.2 * SX, angular_sigma_mrad=ST,
                            x_theta_correlation=RHO), "python")  # fmt: skip
    assert np.allclose(sx[0][:, :2], 1.2 * nom[0][:, :2], rtol=1e-13, atol=1e-13)
    assert np.array_equal(sx[1], nom[1])
    st = _ends(_miss_config(make_config, n, lateral_sigma_mm=SX, angular_sigma_mrad=1.2 * ST,
                            x_theta_correlation=RHO), "python")  # fmt: skip
    assert np.array_equal(st[0], nom[0])
    # theta scales linearly: tan(theta) ratio 1.2 to first order
    t_nom, t_st = np.arctan(nom[1][:, 0] / nom[1][:, 2]), np.arctan(st[1][:, 0] / st[1][:, 2])
    assert np.allclose(t_st, 1.2 * t_nom, rtol=1e-12, atol=1e-15)
    # rho variant: theta = s (rho g1 + sqrt(1-rho^2) g2); recover g2 from the nominal and
    # check the variant equals the same-deviate formula
    g1 = nom[0][:, 0] / SX
    g2 = (t_nom / (ST * 1e-3) - RHO * g1) / math.sqrt(1 - RHO**2)
    for r2 in (-0.417, -0.625):
        v = _ends(_miss_config(make_config, n, lateral_sigma_mm=SX, angular_sigma_mrad=ST,
                               x_theta_correlation=r2), "python")  # fmt: skip
        assert np.array_equal(v[0], nom[0])
        t_v = np.arctan(v[1][:, 0] / v[1][:, 2])
        assert np.allclose(t_v, ST * 1e-3 * (r2 * g1 + math.sqrt(1 - r2**2) * g2), atol=1e-14)


def test_python_vs_warp_cpu_f64_source_parity(make_config: MakeConfig) -> None:
    n = 2000
    cfg = _miss_config(make_config, n, lateral_sigma_mm=SX, angular_sigma_mrad=ST,
                       x_theta_correlation=RHO, energy_sigma_mev=0.8)  # fmt: skip
    py = _ends(cfg, "python")
    wcfg = replace(cfg, run=replace(cfg.run, backend="warp-cpu", precision="float64"))
    wc = _ends(wcfg, "warp-cpu")
    dpos = float(np.max(np.abs(py[0] - wc[0])))
    ddir = float(np.max(np.abs(py[1] - wc[1])))
    print(f"source parity python vs warp-cpu f64: max|dpos|={dpos:.3e} mm, max|ddir|={ddir:.3e}")
    assert dpos <= 1e-12 and ddir <= 1e-14
    f32 = _ends(replace(wcfg, run=replace(wcfg.run, precision="float32")), "warp-cpu")
    assert float(np.max(np.abs(py[1] - f32[1]))) < 1e-5
