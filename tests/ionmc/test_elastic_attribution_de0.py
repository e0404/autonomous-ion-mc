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


# --- Git projection of the exploratory records (decision 0041: no nuclear data values in Git) ---
from validation.scripts.elastic_attribution import common as de_common  # noqa: E402

RESULTS = Path(__file__).resolve().parents[2] / "validation" / "results" / "elastic_attribution"


def test_git_projection_removes_values_and_per_point_arrays() -> None:
    doc = {
        "banner": "x",
        "datasets": [
            {
                "dataset": "A",
                "n_points": 12,
                "theta_cm_deg": list(range(12)),
                "data_mb_sr": [1.0] * 12,
                "rel_err": [0.05] * 12,
                "model_total_mb_sr": [2.0] * 12,
                "model_ratio_minus_1": [1.0] * 12,
                "theta_range_deg": [0.0, 11.0],
                "minima": {
                    "measured": [
                        {"theta_cm_deg": 3.0, "data_mb_sr": 1.0, "neighbour_ratio": [1.2, 1.1]}
                    ]
                },
                "partial": {"10-30": {"model": {"sigma_ratio_to_data": 1.4}}},
            }
        ],
        "pdg_lowest_points": [{"t_mev": 325.0, "sigma_mb": 23.2}],
        "pdg_lowest_energies_mev": [325.0, 418.0],
    }
    proj = de_common.git_projection(doc)
    assert de_common.git_safety_violations(proj) == []
    ds = proj["datasets"][0]
    assert set(ds) == {"dataset", "n_points", "theta_range_deg", "minima", "partial"}
    assert ds["minima"]["measured"][0] == {"theta_cm_deg": 3.0, "neighbour_ratio": [1.2, 1.1]}
    assert "pdg_lowest_points" not in proj and proj["pdg_lowest_energies_mev"] == [325.0, 418.0]
    assert "note_git" not in proj  # the writer adds the note, the projection stays pure
    # the original is flagged, and the projection is idempotent
    assert len(de_common.git_safety_violations(doc)) == 7
    assert de_common.git_projection(proj) == proj


@pytest.mark.parametrize(
    "path", sorted(RESULTS.glob("**/*.json")), ids=lambda p: str(p.relative_to(RESULTS))
)
def test_committed_records_are_git_safe(path: Path) -> None:
    import json

    doc = json.loads(path.read_text())
    assert de_common.git_safety_violations(doc) == []
    assert de_common.git_projection(doc) == doc
