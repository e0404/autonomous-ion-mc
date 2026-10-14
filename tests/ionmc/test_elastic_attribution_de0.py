"""Exploratory D-E0 kernels (Amendment 19, report-only): analytic cases."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validation.scripts.elastic_attribution import de0_kernels as de0  # noqa: E402
from validation.scripts.elastic_attribution import de1_shapes as de1  # noqa: E402


def _lin_range(e: object) -> object:
    return np.asarray(e, dtype=float)  # toy range R(E) = E mm


def _ev(mu: np.ndarray, m_t: float = 1.0e12, e: float = 100.0) -> dict[str, np.ndarray]:
    from ionmc.nuclear.elastic_kin import two_body

    tb = two_body(e, 938.272, m_t, mu)
    return {"mu": mu, "t1": tb.t1, "cos": tb.cos_lab1, "t_other": tb.t2}


def test_isotropic_infinite_target_exceedance() -> None:
    """Isotropic CM, infinite target: Delta = R (1 - mu), P(Delta > w) = 1 - w / (2 R)."""
    u = (np.arange(4096) + 0.5) / 4096
    ev = _ev(2.0 * u - 1.0)
    d = de0.deficit_mm(100.0, ev, _lin_range)
    for w in (1.0, 2.0, 5.0):
        assert np.mean(d > w) == pytest.approx(1.0 - w / 200.0, abs=1e-3)


def test_forward_delta_has_no_deficit() -> None:
    ev = _ev(np.full(64, 1.0))
    assert np.all(de0.deficit_mm(100.0, ev, _lin_range) < 1e-6)
    assert de0.f_pred(0.0, 0.0) == 1.0


def test_threshold_angle() -> None:
    th = de0.theta_w_deg(150.0, 2.0, _lin_range)
    assert 150.0 * (1.0 - np.cos(np.radians(th))) == pytest.approx(2.0)


def test_disk_family_normalised_and_quantiles() -> None:
    """Disk form: 1 at q = 0, first zero at x = 3.8317; the families are 1 at q = 0."""
    assert float(de1.disk_a(np.array([0.0]))[0]) == 1.0
    assert float(de1.disk_a(np.array([de1.J1_FIRST_ZERO]))[0]) < 1e-9
    assert float(de1.family_weight("s0.5", np.array([0.0]), 2.6)[0]) == 1.0
    assert float(de1.family_weight("gauss", np.array([0.0]), 2.6)[0]) == 1.0


@pytest.mark.requires_tables("elastic")
def test_kernels_on_pinned_table() -> None:
    from validation.scripts.elastic_attribution import common as c

    tab = c.load_table()
    sh, water = de1.Shapes(tab), c.Water(tab)
    k = de0.kernels(sh, water, 150, "table", ("O-16",), False, n_z=40)["O-16"]
    assert k["P_el"] > k["K_1mm"] >= k["K_2mm"] >= k["K_3mm"] >= k["K_5mm"] > 0.0
    soft = de0.kernels(sh, water, 150, "s0.9", ("O-16",), False, n_z=40)["O-16"]
    assert soft["K_2mm"] < k["K_2mm"]


def test_kstar_sign_pinned_on_asymmetric_idd() -> None:
    """Synthetic IDD (slow rise, sharp fall): shortened range gives d(z_pk) = IDD(z_pk + Delta)
    (shifted, physical); the as-frozen z_pk - Delta samples the rise."""

    def idd(z: object) -> object:
        z = np.asarray(z, dtype=float)
        return np.where(z <= 100.0, z / 100.0, np.clip(1.0 - (z - 100.0) / 10.0, 0.0, None))

    lit, shifted = de0.kstar_node(idd, 100.0, np.array([5.0, 5.0]))
    assert shifted == pytest.approx(1.0 - 0.5)  # IDD(105) = 0.5
    assert lit == pytest.approx(1.0 - 0.95)  # IDD(95) = 0.95
    assert shifted > lit
