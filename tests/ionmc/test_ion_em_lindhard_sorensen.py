"""Tests of lindhard_sorensen.py (V3-008 check P3).

The module ``validation/scripts/ion_em/lindhard_sorensen.py`` is loaded by file path (like
``test_ion_em_moliere_hanson.py``). The Kats et al. Table 1 values are published theory values
(arXiv:2008.04819v2); the golden table is script-generated (regression only).
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

_PATH = Path(__file__).resolve().parents[2] / "validation/scripts/ion_em/lindhard_sorensen.py"
_SPEC = importlib.util.spec_from_file_location("lindhard_sorensen", _PATH)
assert _SPEC is not None and _SPEC.loader is not None
ls = importlib.util.module_from_spec(_SPEC)
sys.modules["lindhard_sorensen"] = ls
_SPEC.loader.exec_module(ls)

A = ls.ALPHA_FS
EULER_GAMMA = 0.57721566490153286061


# ---------------------------------------------------------------- complex log-Gamma
@pytest.mark.parametrize("y", [0.01, 0.3, 1.0, 2.5, 10.0])
def test_clgamma_modulus_on_line_re1(y: float) -> None:
    """|Gamma(1+iy)|^2 = pi y / sinh(pi y)."""
    got = float(ls.clgamma(np.array([1 + 1j * y]))[0].real)
    ref = 0.5 * math.log(math.pi * y / math.sinh(math.pi * y))
    assert abs(got - ref) < 1e-15 * max(1.0, abs(ref))


@pytest.mark.parametrize("x", [0.5, 1.0, 1.5, 3.7, 15.9, 20.2, 150.0])
def test_clgamma_real_axis(x: float) -> None:
    got = ls.clgamma(np.array([x + 0j]))[0]
    assert abs(float(got.real) - math.lgamma(x)) < 2e-15 * max(1.0, abs(math.lgamma(x)))
    assert abs(float(got.imag)) < 1e-18


@pytest.mark.parametrize("y", [0.05, 0.4, 1.3])
def test_clgamma_phase_series(y: float) -> None:
    """arg Gamma(1+iy) = -gamma_E y + sum_n (y/n - atan(y/n))  (mod 2 pi)."""
    n = np.arange(1, 2_000_001, dtype=np.float64)
    s = -EULER_GAMMA * y + float(np.sum(y / n - np.arctan(y / n))) + y**3 / (6 * 2_000_000.0**2)
    got = float(ls.clgamma(np.array([1 + 1j * y]))[0].imag)
    d = (got - s + math.pi) % (2 * math.pi) - math.pi
    assert abs(d) < 1e-12


def test_clgamma_recurrence() -> None:
    z = np.array([0.3 + 0.7j, 2.2 + 0.05j, 7.5 + 3.0j, 1.0 + 0.9j], dtype=np.clongdouble)
    d = ls.clgamma(z + 1) - ls.clgamma(z) - np.log(z)
    re = np.abs(d.real.astype(float))
    im = np.abs(((d.imag.astype(float) + math.pi) % (2 * math.pi)) - math.pi)
    assert np.all(re < 1e-15) and np.all(im < 1e-15)


# ---------------------------------------------------------------- (a) small-eta limit, k-sum
@pytest.mark.parametrize("beta", [0.3, 0.6, 0.9])
def test_small_eta_first_order_mott_emerges(beta: float) -> None:
    """LS - Bloch -> (pi/2) alpha z beta with relative error O(eta) that halves with z."""
    errs = []
    for z in (1.0, 0.5, 0.25):
        az = A * z
        eta = az / beta
        d = ls.delta_l_ls(az, beta).delta_l - ls.bloch_exact(eta)
        rel = d / ls.mott1_ionmc(beta, z) - 1.0
        assert 0.0 < rel < 2.0 * eta  # tolerance shrinks with eta
        errs.append(rel)
    assert errs[1] / errs[0] == pytest.approx(0.5, rel=0.05)
    assert errs[2] / errs[1] == pytest.approx(0.5, rel=0.05)


@pytest.mark.parametrize("beta", [0.2, 0.5, 0.8])
def test_small_eta_ls_equals_bloch_plus_mott1_to_second_order(beta: float) -> None:
    """For z = 1: |LS - (Bloch + Mott1)| <= 3 (alpha z)^2, and the residual scales as z^2."""
    r = {}
    for z in (1.0, 0.5):
        az = A * z
        r[z] = ls.delta_l_ls(az, beta).delta_l - ls.bloch_exact(az / beta) - ls.mott1_ionmc(beta, z)
        assert abs(r[z]) <= 3.0 * az**2
    assert r[1.0] / r[0.5] == pytest.approx(4.0, rel=0.03)


@pytest.mark.parametrize(("z", "beta"), [(1, 0.5), (6, 0.15), (8, 0.0654), (8, 0.72)])
def test_k_sum_converges(z: int, beta: float) -> None:
    """Tail-corrected sums agree for K = 1000 ... 16000; the raw partial sums approach them as
    c2/K with c2 = first-order Mott term."""
    ref = ls.delta_l_ls(A * z, beta, k_start=16000, tol=1e-11)
    lo = ls.delta_l_ls(A * z, beta, k_start=1000, tol=1e-11)
    assert abs(ref.delta_l - lo.delta_l) < 1e-10
    assert ref.tail_coeff_c2 == pytest.approx(ls.mott1_ionmc(beta, z), rel=1e-6)
    t, b2, _ = ls.ls_summand(A * z, beta, 16000)
    raw16 = float(np.sum(t, dtype=np.longdouble)) + b2
    raw8 = float(np.sum(t[:8000], dtype=np.longdouble)) + b2
    # untreated truncation error falls as 1/K (summand ~ c2/k^2): halving K doubles it
    assert (ref.delta_l - raw8) / (ref.delta_l - raw16) == pytest.approx(2.0, rel=0.01)
    assert ref.delta_l - raw16 == pytest.approx(ref.tail_coeff_c2 / 16000, rel=0.01)


def test_convergence_failure_is_reported() -> None:
    with pytest.raises(RuntimeError):
        ls.delta_l_ls(A * 8, 0.5, k_start=8, k_limit=16, tol=1e-30)


# ---------------------------------------------------------------- (b) non-relativistic limit
@pytest.mark.parametrize("eta", [0.1, 0.5, 1.0, 2.0])
def test_nonrelativistic_limit_is_bloch(eta: float) -> None:
    """At fixed eta and beta -> 0 (alpha z = eta beta -> 0, gamma -> 1) LS -> Bloch, with a
    deviation that falls as beta^2; the first-order Mott term vanishes as beta^2 at fixed eta."""
    dev = []
    for beta in (1e-2, 1e-3):
        az = eta * beta
        d = ls.delta_l_ls(az, beta).delta_l - ls.bloch_exact(eta)
        mott = ls.mott1_ionmc(beta, az / A)
        assert mott == pytest.approx(0.5 * math.pi * eta * beta**2, rel=1e-12)
        assert abs(d) < 5.0 * eta * beta**2 * (1 + eta)
        dev.append(d)
    assert abs(dev[1]) < 0.02 * abs(dev[0])


def test_nonrelativistic_summand_is_bloch_summand() -> None:
    eta = 0.7
    beta = 1e-5
    t, b2, _ = ls.ls_summand(eta * beta * 1.0, beta, 50)
    k = np.arange(1, 51, dtype=np.float64)
    bloch_t = k / (k * k + eta * eta) - 1 / k
    assert np.max(np.abs(t.astype(float) - bloch_t)) < 1e-8
    assert b2 < 1e-9


# ---------------------------------------------------------------- (c) residual behaviour
E_GRID = np.geomspace(2.0, 430.0, 61)


@pytest.mark.parametrize("z", [6, 8])
def test_residual_finite_positive_smooth(z: int) -> None:
    r = np.array([ls.residual(z, float(e)) for e in E_GRID])
    assert np.all(np.isfinite(r))
    assert np.all(r > 0) and np.all(r < 0.01)
    d2 = np.diff(r, 2)  # second difference on a uniform ln E grid
    assert np.max(np.abs(d2)) < 0.02 * np.max(r)


def test_residual_scales_as_z_squared_at_high_velocity() -> None:
    """At equal velocity (E/u fixed, beta >~ 0.3) r(8)/r(6) ~ (8/6)^2."""
    for e in (100.0, 290.0, 430.0):
        b6, _ = ls.kinematics(6, e)
        b8, _ = ls.kinematics(8, e)
        r6 = ls.delta_l_ls(A * 6, b6).delta_l - ls.bloch_ionmc(b6**2, 6) - ls.mott1_ionmc(b6, 6)
        r8 = ls.delta_l_ls(A * 8, b6).delta_l - ls.bloch_ionmc(b6**2, 8) - ls.mott1_ionmc(b6, 8)
        assert r8 / r6 == pytest.approx((8 / 6) ** 2, rel=0.15)
        assert abs(b6 - b8) < 1e-3
        assert 1.5 < r6 / (A * 6) ** 2 < 2.5


# ---------------------------------------------------------------- literature check (theory)
# Kats, Halenka, Voskresenskaya, arXiv:2008.04819v2, Table 1, rows "Delta L_LS" (point
# nucleus, eq. 10, Mathematica, summed to k = 5000, 3 decimals). Published theory values.
KATS_BETA = (0.85, 0.90, 0.95, 0.97, 0.99)
KATS_LS = {
    6: (0.059, 0.063, 0.067, 0.068, 0.070),
    12: (0.120, 0.128, 0.136, 0.139, 0.142),
    26: (0.267, 0.288, 0.309, 0.317, 0.325),
    36: (0.377, 0.411, 0.443, 0.455, 0.467),
    52: (0.562, 0.621, 0.676, 0.698, 0.718),
    60: (0.659, 0.733, 0.802, 0.829, 0.855),
    80: (0.903, 1.024, 1.140, 1.184, 1.228),
    92: (1.040, 1.196, 1.345, 1.404, 1.461),
    104: (1.145, 1.338, 1.527, 1.601, 1.675),
    114: (1.170, 1.392, 1.614, 1.702, 1.789),
}


@pytest.mark.parametrize("z", sorted(KATS_LS))
def test_against_kats_table1(z: int) -> None:
    for beta, ref in zip(KATS_BETA, KATS_LS[z], strict=True):
        got = ls.delta_l_ls(A * z, beta).delta_l
        # rounding (5e-4) + their truncation at k = 5000 without tail (c2/5000 = Mott1/5000)
        tol = 5e-4 + ls.mott1_ionmc(beta, z) / 5000 + 1e-6
        assert abs(got - ref) <= tol, (z, beta, got, ref)


# ---------------------------------------------------------------- (d) golden regression values
# GOLDEN values SCRIPT-GENERATED by the module itself (lindhard_sorensen.residual, bare charge,
# 2026-11-01). Regression only, not an independent reference.
GOLDEN = {
    (6, 2.0): 0.0019509452112704152,
    (6, 10.0): 0.0033215801140141204,
    (6, 50.0): 0.0038368927989918707,
    (6, 100.0): 0.0039118771851832865,
    (6, 290.0): 0.00392887432382872,
    (6, 430.0): 0.003914109439582515,
    (8, 2.0): 0.0023198118008166593,
    (8, 10.0): 0.005240583448276102,
    (8, 50.0): 0.006610538553542766,
    (8, 100.0): 0.006835359071730121,
    (8, 290.0): 0.006942949073304422,
    (8, 430.0): 0.006933679345884944,
}


@pytest.mark.parametrize(("z", "e"), sorted(GOLDEN))
def test_golden_residual(z: int, e: float) -> None:
    assert ls.residual(z, e) == pytest.approx(GOLDEN[(z, e)], rel=1e-6)


# ---------------------------------------------------------------- consistency with ionmc
def _ionmc():
    from ionmc.physics import projectiles, stopping

    return stopping, projectiles


def test_forms_match_ionmc() -> None:
    st, pj = _ionmc()
    assert ls.ALPHA_FS == st.ALPHA_FS
    proj = {1: pj.PROTON, 2: pj.ALPHA, 6: pj.CARBON12, 8: pj.OXYGEN16}
    for z, p in proj.items():
        assert ls.ION_MASS[z] == (p.a, p.mass_mev)
        e = np.array([2.0, 37.0, 430.0])
        beta2, _, _ = st.kinematics(e, p)
        q = st.pierce_blann_charge(np.sqrt(beta2), p) if z >= 2 else np.ones(3)
        rep = st.bloch_term(beta2, q)
        for i in range(3):
            assert ls.kinematics(z, float(e[i]))[0] ** 2 == pytest.approx(beta2[i], rel=1e-14)
            assert ls.bloch_ionmc(float(beta2[i]), float(q[i])) == pytest.approx(rep[i], rel=1e-14)
            if z >= 2:
                b = math.sqrt(beta2[i])
                assert ls.effective_charge_ionmc(z, b) == pytest.approx(q[i], rel=1e-14)
    # Mott term: ionmc adds 0.5*pi*alpha*beta*q inside the bracket (z >= 2)
    assert ls.mott1_ionmc(0.5, 6.0) == 0.5 * math.pi * st.ALPHA_FS * 0.5 * 6.0
