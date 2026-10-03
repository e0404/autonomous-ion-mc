"""Formula checks of the shared functions against independent NumPy implementations.

``test_shared_funcs.py`` shows that Python scope and kernels agree; these tests show that the
formulas themselves are the documented ones (independent of the Warp code path).
"""

import math

import numpy as np
import pytest
import warp as wp

from ionmc.physics.em import make_em
from ionmc.physics.kinematics import make_kinematics
from ionmc.physics.projectiles import ELECTRON_MASS_MEV, PROJECTILES, PROTON
from ionmc.physics.stopping import kinematics as numpy_kinematics
from ionmc.transport.funcs import make_transport_funcs

wp.config.log_level = wp.LOG_WARNING
R = wp.float64
KIN = make_kinematics(R)
EM = make_em(R)
TF = make_transport_funcs(R)
V3 = TF.vec3
RE_CM = 2.8179403262e-13
N_A = 6.02214076e23


def _vec(a: np.ndarray | tuple[float, float, float]) -> object:
    return V3(R(float(a[0])), R(float(a[1])), R(float(a[2])))


@pytest.mark.parametrize("projectile", ["proton", "alpha", "carbon12"])
def test_kinematics_match_the_stopping_module(projectile: str) -> None:
    p = PROJECTILES[projectile]
    for e_u in (1.0, 10.0, 100.0, 400.0):
        beta2, bg2, tmax = numpy_kinematics(e_u, p)
        t = e_u * p.a
        assert float(KIN.beta2(R(t), R(p.mass_mev))) == pytest.approx(float(beta2), rel=1e-13)
        assert float(KIN.tmax_mev(R(t), R(p.mass_mev))) == pytest.approx(float(tmax), rel=1e-13)
        assert float(KIN.gamma(R(t), R(p.mass_mev))) ** 2 - 1.0 == pytest.approx(
            float(bg2), rel=1e-12
        )
        pv = t * (t + 2 * p.mass_mev) / (t + p.mass_mev)
        assert float(KIN.pv_mev(R(t), R(p.mass_mev))) == pytest.approx(pv, rel=1e-14)
        # p v = p c beta with p c = sqrt(T (T + 2 M))
        pc = math.sqrt(t * (t + 2 * p.mass_mev))
        assert pv == pytest.approx(pc * math.sqrt(float(beta2)), rel=1e-12)


def test_bohr_variance_is_the_geant4_dispersion() -> None:
    """2 pi r_e^2 m_e c^2 n_el z^2 x T_max (1/beta^2 - 1/2) with n_el = N_A (Z/A) rho."""
    for e, rho, s in ((10.0, 1.0, 0.5), (100.0, 1.19, 1.0), (200.0, 0.0012, 5.0)):
        z_over_a = 0.5551
        beta2, _, tmax = numpy_kinematics(e, PROTON)
        n_el = N_A * z_over_a * rho
        expected = (
            2.0
            * math.pi
            * RE_CM**2
            * ELECTRON_MASS_MEV
            * n_el
            * (s / 10.0)
            * float(tmax)
            * (1.0 / float(beta2) - 0.5)
        )
        got = float(EM.bohr_variance(R(e), R(PROTON.mass_mev), R(1.0), R(z_over_a), R(rho), R(s)))
        assert got == pytest.approx(expected, rel=2e-6)  # K = 0.307075 is rounded to 7 digits


def test_mean_loss_branches() -> None:
    f = EM.csda_mean_loss
    # telescoping branch: difference of the two energies
    assert float(f(R(100.0), R(99.0), R(7.0), R(0.2), R(7.7), R(1e-3))) == 1.0
    # short branch: S t
    assert float(f(R(100.0), R(100.0), R(7.0), R(0.001), R(7.7), R(1e-3))) == pytest.approx(0.007)
    # whole residual range travelled: everything is lost
    assert float(f(R(5.0), R(0.0), R(30.0), R(0.9), R(0.9), R(1e-3))) == 5.0
    assert float(f(R(5.0), R(0.0), R(30.0), R(2.0), R(0.9), R(1e-3))) == 5.0
    # never negative and never above the energy
    assert float(f(R(100.0), R(100.5), R(7.0), R(0.2), R(7.7), R(1e-3))) == 0.0
    assert float(f(R(100.0), R(-5.0), R(7.0), R(0.2), R(7.7), R(1e-3))) == 100.0


def test_fdm_and_scattering_power_follow_gottschalk() -> None:
    pv, p1v1 = 400.0, 540.0
    l1 = math.log10(1.0 - (pv / p1v1) ** 2)
    l2 = math.log10(pv)
    expected = 0.5244 + 0.1975 * l1 + 0.2320 * l2 - 0.0098 * l1 * l2
    assert float(EM.fdm(R(pv), R(p1v1))) == pytest.approx(expected, rel=1e-14)
    # clamped to zero where the fit diverges (pv -> p1v1) and for pv >= p1v1
    assert float(EM.fdm(R(p1v1 * (1 - 1e-12)), R(p1v1))) == 0.0
    assert float(EM.fdm(R(p1v1), R(p1v1))) == 0.0
    assert float(EM.fdm(R(1.01 * p1v1), R(p1v1))) == 0.0
    inv_xs, rho = 1.0 / 46.88, 1.0
    t = float(EM.scattering_power_dm(R(pv), R(p1v1), R(1.0), R(inv_xs), R(rho)))
    assert t == pytest.approx(expected * (15.0 / pv) ** 2 * rho * inv_xs / 10.0, rel=1e-14)
    # the charge number enters squared, the density linearly
    t2 = float(EM.scattering_power_dm(R(pv), R(p1v1), R(2.0), R(inv_xs), R(2.0 * rho)))
    assert t2 == pytest.approx(8.0 * t, rel=1e-14)


def test_step_limit_formulas() -> None:
    a, rho_f = 0.2, 0.1
    for r in (0.5, 10.0, 77.0):
        expected = a * r + rho_f * (1 - a) * (2 - rho_f / r)
        assert float(TF.range_step_limit(R(r), R(a), R(rho_f))) == pytest.approx(expected)
    assert float(TF.range_step_limit(R(0.05), R(a), R(rho_f))) == 0.05  # final range step
    assert float(TF.range_step_limit(R(rho_f), R(a), R(rho_f))) == rho_f
    assert float(TF.eloss_step_limit(R(100.0), R(0.73), R(0.02))) == pytest.approx(2.0 / 0.73)
    # geometry wins ties; otherwise the smallest limit with its reason
    assert TF.select_step(R(1.0), R(1.0), R(2.0), R(3.0))[1] == 0
    assert TF.select_step(R(9.0), R(4.0), R(2.0), R(3.0))[1] == 2
    assert TF.select_step(R(9.0), R(4.0), R(5.0), R(3.0))[1] == 3
    assert TF.select_step(R(9.0), R(2.0), R(5.0), R(3.0))[1] == 1
    s, reason = TF.select_step(R(0.4), R(2.0), R(5.0), R(3.0))
    assert float(s) == 0.4 and reason == 0


def _dda_reference(
    p: np.ndarray, u: np.ndarray, idx: np.ndarray, origin: np.ndarray, spacing: np.ndarray
) -> tuple[float, int]:
    best, axis = 1e30, 0
    for a in range(3):
        if abs(u[a]) < 1e-12:
            continue
        plane = origin[a] + spacing[a] * (idx[a] + (1 if u[a] > 0 else 0))
        d = max((plane - p[a]) / u[a], 0.0)
        if d < best:
            best, axis = d, a
    return best, axis


def test_dda_matches_a_numpy_reference_and_walks_a_ray_exactly() -> None:
    rng = np.random.default_rng(3)
    origin, spacing = np.array([0.5, -1.0, 2.0]), np.array([1.0, 2.0, 0.5])
    for _ in range(300):
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        idx = rng.integers(-3, 5, 3)
        p = origin + spacing * (idx + rng.uniform(0, 1, 3))
        ix, iy, iz = (int(i) for i in idx)
        d, ax = TF.dda_next(_vec(p), _vec(u), ix, iy, iz, _vec(origin), _vec(spacing))
        d_ref, ax_ref = _dda_reference(p, u, idx, origin, spacing)
        assert float(d) == pytest.approx(d_ref, rel=1e-13) and ax == ax_ref
    # walking a diagonal through corners: every crossing is found, positions stay on the ray
    u = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
    p, idx = np.array([0.5, -1.0, 2.0]), np.array([0, 0, 0])
    for _ in range(12):
        d, ax = TF.dda_next(_vec(p), _vec(u), *(int(i) for i in idx), _vec(origin), _vec(spacing))
        p = p + u * float(d)
        idx[ax] += 1
        assert np.allclose(np.cross(p - np.array([0.5, -1.0, 2.0]), u), 0.0, atol=1e-12)


def test_ray_box_and_grid_index() -> None:
    lo, hi = np.array([0.0, 0.0, 0.0]), np.array([2.0, 3.0, 4.0])
    cases = [
        ((-1.0, 1.0, 1.0), (1.0, 0.0, 0.0), (1.0, 3.0)),
        ((1.0, 1.0, 1.0), (0.0, 0.0, 1.0), (-1.0, 3.0)),
        ((1.0, 5.0, 1.0), (0.0, 0.0, 1.0), None),  # parallel and outside the slab
        ((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0), (1.0, 3.0)),  # through the lower corner
    ]
    for p, u, expected in cases:
        norm = np.linalg.norm(u)
        t_in, t_out = TF.ray_box(_vec(p), _vec(np.array(u) / norm), _vec(lo), _vec(hi))
        if expected is None:
            assert float(t_in) > float(t_out)
        else:
            assert float(t_in) == pytest.approx(expected[0] * norm)
            assert float(t_out) == pytest.approx(expected[1] * norm)
    # grid index: floor, far edge outside, negative side
    origin, inv = _vec((0.0, 0.0, 0.0)), _vec((1.0, 1.0, 1.0))
    assert TF.grid_index(_vec((0.5, 1.5, 2.5)), origin, inv, 3, 3, 3) == (0, 1, 2, 1)
    assert TF.grid_index(_vec((3.0, 1.5, 2.5)), origin, inv, 3, 3, 3)[3] == 0
    assert TF.grid_index(_vec((-0.1, 1.5, 2.5)), origin, inv, 3, 3, 3) == (-1, 1, 2, 0)
    assert TF.grid_index(_vec((0.0, 0.0, 0.0)), origin, inv, 3, 3, 3) == (0, 0, 0, 1)


def test_hinge_point_and_leg_cut() -> None:
    p0, d0, d1 = _vec((0.0, 0.0, 0.0)), _vec((0.0, 0.0, 1.0)), _vec((1.0, 0.0, 0.0))
    mid = TF.point_on_hinge(p0, d0, R(2.0), d1, R(1.0))
    assert [float(c) for c in mid] == [0.0, 0.0, 1.0]
    end = TF.point_on_hinge(p0, d0, R(2.0), d1, R(5.0))
    assert [float(c) for c in end] == [3.0, 0.0, 2.0]
    origin, spacing = _vec((0.0, 0.0, 0.0)), _vec((1.0, 1.0, 1.0))
    ph = _vec((0.5, 0.5, 0.5))
    leg, axis = TF.leg2_limit(ph, d1, 0, 0, 0, origin, spacing, R(2.0))
    assert float(leg) == 0.5 and axis == 0  # cut at x = 1
    leg, axis = TF.leg2_limit(ph, d1, 0, 0, 0, origin, spacing, R(0.25))
    assert float(leg) == 0.25 and axis == -1  # does not reach the plane
    leg, axis = TF.leg2_limit(ph, d1, 0, 0, 0, origin, spacing, R(0.5 * (1.0 - 1e-12)))
    assert axis == 0  # planned boundary step within the tolerance still crosses


def test_log_bin_index_and_interpolation() -> None:
    i, f = TF.log_bin_index(R(math.exp(2.5)), R(0.0), R(2.0), 10)
    assert i + float(f) == pytest.approx(5.0, abs=1e-12)  # ln x * inv_dl = 5 on the grid point
    i, f = TF.log_bin_index(R(math.exp(3.3)), R(0.0), R(2.0), 10)
    assert i == 6 and float(f) == pytest.approx(0.6, abs=1e-12)
    i, f = TF.log_bin_index(R(1e-9), R(0.0), R(2.0), 10)
    assert (i, float(f)) == (0, 0.0)  # clamped below the grid
    i, f = TF.log_bin_index(R(1e9), R(0.0), R(2.0), 10)
    assert (i, float(f)) == (8, 1.0)  # clamped above the grid
    assert float(TF.interp_exp(R(math.log(2.0)), R(math.log(8.0)), R(0.5))) == pytest.approx(4.0)
    assert float(TF.lerp(R(1.0), R(3.0), R(0.25))) == 1.5


def test_box_muller_formula() -> None:
    for u0, u1 in ((0.5, 0.0), (0.1, 0.3), (0.99, 0.77), (1e-6, 0.25)):
        g0, g1 = TF.gauss_pair(R(u0), R(u1))
        r = math.sqrt(-2.0 * math.log(u0))
        assert float(g0) == pytest.approx(r * math.cos(2 * math.pi * u1), abs=1e-13)
        assert float(g1) == pytest.approx(r * math.sin(2 * math.pi * u1), abs=1e-13)
        assert float(TF.gauss_one(R(u0), R(u1))) == pytest.approx(float(g0), abs=1e-14)
