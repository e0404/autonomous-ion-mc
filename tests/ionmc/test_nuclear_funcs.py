"""V3-005A shared nuclear functions (decision 0041; acceptance rows P2, P3, P5).

* P2 - samplers: Kalbach mu (a x r grid), Poisson multiplicities, inverse-CDF E' on a synthetic
  64-bin table; 10^6 draws each from the fixed Philox seed ``SEED`` (purpose 2),
  chi-square goodness of fit (Wilson-Hilferty p-value, ``ionmc.transport.parity``) and the Kalbach
  first moment within 4 sigma. The draws are made by a compiled Warp CPU float64 kernel that calls
  the SHARED functions (the same source text as the twins; kernel and twin agree to 1e-12, P5), so
  that 10^6 draws take seconds in single-threaded CI instead of minutes with scalar twins.
* P3 - ledger closure of the event sampler: ``test_nuclear_table.py`` (synthetic tables).
* P5 - python twin against a Warp CPU kernel for EVERY function of ``make_nuclear`` (float64:
  discrete outputs equal, continuous within 1e-12 relative; float32 kernel against the float64
  twin within analytic budgets), and the U1 completeness check of ``test_shared_funcs``.
"""

import math
import os
from typing import Any

import numpy as np
import pytest
import warp as wp

from ionmc._wpfunc import python_twin
from ionmc.physics.nuclear import (
    KALBACH_I,
    KALBACH_M_B,
    kalbach_separation_energy,
    make_nuclear,
)
from ionmc.rng.philox import PURPOSE_NUCLEAR, PURPOSE_RESERVED, key_from_seed, make_philox
from ionmc.transport.parity import wilson_hilferty_p

wp.config.log_level = wp.LOG_WARNING

SEED = 20261007
N_DRAWS = 1_000_000
P_MIN = 0.001
U = 2.0**-24
N_ARGS = int(os.environ.get("IONMC_U1_N", "10000"))

# masses [MeV] (AME2020 atomic masses x u - Z m_e; hard-coded, the cache is not needed in CI)
U_MEV = 931.49410242
M_E = 0.51099895
M_P = 938.27208816
M_N = 939.56542052
M_D = 1875.61294257
M_A = 3727.3794066
M_C12 = 12.0 * U_MEV - 6 * M_E
M_O16 = 15.99491461957 * U_MEV - 8 * M_E

NU64 = python_twin(make_nuclear)


def test_purpose_nuclear_alias() -> None:
    assert PURPOSE_NUCLEAR == 2 and PURPOSE_RESERVED == PURPOSE_NUCLEAR


# ---------------------------------------------------------------------------------------------
# P2: samplers (Warp CPU float64 kernels over the shared functions)
# ---------------------------------------------------------------------------------------------
_PH = make_philox(wp.float64)
_NU = make_nuclear(wp.float64)
_NU32 = make_nuclear(wp.float32)
_KEY = key_from_seed(SEED)


@wp.kernel(module="unique")
def _k_kalbach(
    out: wp.array(dtype=wp.float64),  # type: ignore[valid-type]
    a: wp.float64,
    r: wp.float64,
    stream: wp.uint32,
    s0: wp.uint32,
    s1: wp.uint32,
) -> None:
    i = wp.tid()
    key = wp.vec2ui(s0, s1)
    w = _PH.philox_block(wp.uint32(i), stream, wp.uint32(0), wp.uint32(2), key)
    out[i] = _NU.kalbach_mu(_PH.u01(w[0]), a, r)


@wp.kernel(module="unique")
def _k_poisson(
    out: wp.array(dtype=int),  # type: ignore[valid-type]
    lam: wp.float64,
    stream: wp.uint32,
    s0: wp.uint32,
    s1: wp.uint32,
) -> None:
    i = wp.tid()
    key = wp.vec2ui(s0, s1)
    w = _PH.philox_block(wp.uint32(i), stream, wp.uint32(0), wp.uint32(2), key)
    out[i] = _NU.poisson_inverse(_PH.u01(w[0]), lam, 16)


@wp.kernel(module="unique")
def _k_mult(
    out: wp.array(dtype=int),  # type: ignore[valid-type]
    lam: wp.float64,
    stream: wp.uint32,
    s0: wp.uint32,
    s1: wp.uint32,
) -> None:
    i = wp.tid()
    key = wp.vec2ui(s0, s1)
    w = _PH.philox_block(wp.uint32(i), stream, wp.uint32(0), wp.uint32(2), key)
    out[i] = _NU.multiplicity_round(_PH.u01(w[0]), lam, 16)


@wp.kernel(module="unique")
def _k_invcdf(
    out: wp.array(dtype=wp.float64),  # type: ignore[valid-type]
    edges: wp.array(dtype=wp.float64),  # type: ignore[valid-type]
    stream: wp.uint32,
    s0: wp.uint32,
    s1: wp.uint32,
) -> None:
    i = wp.tid()
    key = wp.vec2ui(s0, s1)
    w = _PH.philox_block(wp.uint32(i), stream, wp.uint32(0), wp.uint32(2), key)
    k = _NU.inv_cdf_bin(_PH.u01(w[0]), 64)
    out[i] = _NU.inv_cdf_sample(_PH.u01(w[1]), edges[k], edges[k + 1])


def _chi2_uniform(x: np.ndarray, n_bins: int) -> tuple[float, float]:
    """``(chi2, p)`` of ``x`` in [0, 1) against the uniform distribution on ``n_bins`` equiprobable
    bins (expected count N / n_bins >= 5 required)."""
    counts = np.bincount(np.minimum((x * n_bins).astype(np.int64), n_bins - 1), minlength=n_bins)
    exp = len(x) / n_bins
    assert exp >= 5.0
    chi2 = float(np.sum((counts - exp) ** 2 / exp))
    return chi2, wilson_hilferty_p(chi2, n_bins - 1)


def _kalbach_cdf_np(mu: np.ndarray, a: float, r: float) -> np.ndarray:
    """Independent numpy CDF ``[sinh(a mu) + sinh a + r (cosh(a mu) - cosh a)] / (2 sinh a)``."""
    return (np.sinh(a * mu) + np.sinh(a) + r * (np.cosh(a * mu) - np.cosh(a))) / (2.0 * np.sinh(a))


def _kalbach_moments(a: float, r: float) -> tuple[float, float]:
    """Analytic mean and variance of mu.

    ``<mu> = int mu f dmu`` with ``f = a [cosh(a mu) + r sinh(a mu)] / (2 sinh a)``: the cosh part
    is odd in mu and integrates to zero. ``int_{-1}^{1} mu sinh(a mu) dmu =
    [mu cosh(a mu) / a - sinh(a mu) / a^2]_{-1}^{1} = 2 cosh a / a - 2 sinh a / a^2``, so
    ``<mu> = r a / (2 sinh a) (2 cosh a / a - 2 sinh a / a^2) = r (coth a - 1 / a)``.
    Likewise the sinh part of ``<mu^2>`` is odd and
    ``int mu^2 cosh(a mu) dmu = [mu^2 sinh(a mu) / a - 2 mu cosh(a mu) / a^2 + 2 sinh(a mu) / a^3]``
    ``= 2 sinh a / a - 4 cosh a / a^2 + 4 sinh a / a^3``, so
    ``<mu^2> = 1 - 2 coth a / a + 2 / a^2`` (independent of r) and ``var = <mu^2> - <mu>^2``."""
    coth = 1.0 / math.tanh(a)
    mean = r * (coth - 1.0 / a)
    second = 1.0 - 2.0 * coth / a + 2.0 / (a * a)
    return mean, second - mean * mean


KALBACH_CASES = [(a, r) for a in (0.01, 1.0, 5.0) for r in (0.0, 0.5, 1.0)]
OBSERVED: dict[str, float] = {}


@pytest.mark.parametrize(("a", "r"), KALBACH_CASES)
def test_p2_kalbach_mu_chi2_and_first_moment(a: float, r: float) -> None:
    case = KALBACH_CASES.index((a, r))
    out = wp.zeros(N_DRAWS, dtype=wp.float64, device="cpu")
    wp.launch(
        _k_kalbach,
        dim=N_DRAWS,
        inputs=[out, a, r, case, _KEY[0], _KEY[1]],
        device="cpu",
    )
    mu = out.numpy()
    assert mu.min() >= -1.0 and mu.max() <= 1.0
    chi2, p = _chi2_uniform(_kalbach_cdf_np(mu, a, r), 100)
    mean, var = _kalbach_moments(a, r)
    z = (mu.mean() - mean) / math.sqrt(var / N_DRAWS)
    OBSERVED[f"kalbach a={a} r={r}"] = p
    print(f"P2 kalbach a={a} r={r}: chi2={chi2:.1f} (99 dof) p={p:.4f} moment z={z:+.2f}")
    assert p > P_MIN
    assert abs(z) < 4.0


def test_p2_kalbach_moment_formulas_against_quadrature() -> None:
    """The analytic mean and variance (derivation in ``_kalbach_moments``) against a midpoint
    quadrature of the density."""
    m = 200_000
    mu = -1.0 + (np.arange(m) + 0.5) * 2.0 / m
    for a, r in KALBACH_CASES:
        f = a * (np.cosh(a * mu) + r * np.sinh(a * mu)) / (2.0 * np.sinh(a))
        w = 2.0 / m
        mean_q = float(np.sum(mu * f) * w)
        second_q = float(np.sum(mu * mu * f) * w)
        mean, var = _kalbach_moments(a, r)
        assert float(np.sum(f) * w) == pytest.approx(1.0, abs=1e-9)
        assert mean_q == pytest.approx(mean, abs=1e-8)
        assert second_q - mean_q**2 == pytest.approx(var, abs=1e-8)


POISSON_LAMBDAS = (0.1, 1.0, 4.0)


@pytest.mark.parametrize("lam", POISSON_LAMBDAS)
def test_p2_poisson_chi2(lam: float) -> None:
    """Categories 0..15 and the cap 16 (tail mass assigned to 16), grouped from the right until
    every group expects >= 5 counts (a discrete variable has no 50 equiprobable bins)."""
    case = 10 + POISSON_LAMBDAS.index(lam)
    out = wp.zeros(N_DRAWS, dtype=int, device="cpu")
    wp.launch(_k_poisson, dim=N_DRAWS, inputs=[out, lam, case, _KEY[0], _KEY[1]], device="cpu")
    n = out.numpy()
    assert n.min() >= 0 and n.max() <= 16
    counts = np.bincount(n, minlength=17).astype(float)
    pmf = np.array([math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1)) for k in range(17)])
    pmf[16] = 1.0 - pmf[:16].sum()  # the cap: tail mass above 15 is assigned to 16
    exp = pmf * N_DRAWS
    # group the right tail until the expected count is at least 5
    last = 16
    while exp[: last + 1][last:].sum() < 5.0 and last > 0:
        exp[last - 1] += exp[last]
        counts[last - 1] += counts[last]
        exp[last], counts[last] = 0.0, 0.0
        last -= 1
    k = np.flatnonzero(exp > 0)
    assert np.all(exp[k] >= 5.0)
    chi2 = float(np.sum((counts[k] - exp[k]) ** 2 / exp[k]))
    dof = len(k) - 1
    p = wilson_hilferty_p(chi2, dof)
    mean = n.mean()
    z = (mean - lam) / math.sqrt(lam / N_DRAWS)
    OBSERVED[f"poisson lam={lam}"] = p
    print(
        f"P2 poisson lam={lam}: groups={len(k)} chi2={chi2:.2f} ({dof} dof) p={p:.4f} "
        f"mean z={z:+.2f}"
    )
    assert p > P_MIN
    assert abs(z) < 4.0


@pytest.mark.parametrize("lam", [0.1, 1.0, 2.4])
def test_p2_multiplicity_round_two_point_frequencies(lam: float) -> None:
    """floor + Bernoulli: exactly two values, floor(lam) and floor(lam) + 1, with probabilities
    1 - frac and frac (chi-square on the two cells, p > 0.001)."""
    out = wp.zeros(N_DRAWS, dtype=int, device="cpu")
    case = 20 + [0.1, 1.0, 2.4].index(lam)
    wp.launch(_k_mult, dim=N_DRAWS, inputs=[out, lam, case, _KEY[0], _KEY[1]], device="cpu")
    n = out.numpy()
    lo = math.floor(lam)
    frac = lam - lo
    assert set(np.unique(n)) <= {lo, lo + 1}
    n_hi = int(np.sum(n == lo + 1))
    if frac == 0.0:  # integer mean: deterministic
        assert n_hi == 0 and np.all(n == lo)
        OBSERVED[f"multiplicity_round lam={lam}"] = 1.0
        return
    exp_hi, exp_lo = N_DRAWS * frac, N_DRAWS * (1.0 - frac)
    chi2 = (n_hi - exp_hi) ** 2 / exp_hi + (N_DRAWS - n_hi - exp_lo) ** 2 / exp_lo
    p = wilson_hilferty_p(chi2, 1)
    OBSERVED[f"multiplicity_round lam={lam}"] = p
    print(f"P2 multiplicity_round lam={lam}: chi2={chi2:.3f} (1 dof) p={p:.4f}")
    assert p > P_MIN
    assert abs(float(n.mean()) - lam) < 5.0 * math.sqrt(frac * (1.0 - frac) / N_DRAWS) + 1e-12


def test_multiplicity_round_exact_cases_and_cap() -> None:
    f = NU64.multiplicity_round
    assert f(0.5, 2.0, 16) == 2 and f(0.999, 2.0, 16) == 2 and f(0.001, 0.0, 16) == 0
    assert f(0.29, 2.3, 16) == 3 and f(0.31, 2.3, 16) == 2  # u < frac chooses the upper value
    assert f(0.01, 2.3, 2) == 2 and f(0.01, 40.0, 99) == 16  # n_max and the cap 16


def _synthetic_edges() -> np.ndarray:
    """65 edges of 64 equiprobable bins of a density ``exp(-E / 8)`` truncated to [0, 40] MeV (the
    quantile function is analytic)."""
    q = np.arange(65) / 64.0
    tail = 1.0 - math.exp(-40.0 / 8.0)
    return -8.0 * np.log(1.0 - q * tail)


def test_p2_inverse_cdf_chi2_and_mean() -> None:
    edges = _synthetic_edges()
    out = wp.zeros(N_DRAWS, dtype=wp.float64, device="cpu")
    wp.launch(
        _k_invcdf,
        dim=N_DRAWS,
        inputs=[out, wp.array(edges, dtype=wp.float64, device="cpu"), 20, _KEY[0], _KEY[1]],
        device="cpu",
    )
    e = out.numpy()
    assert e.min() >= edges[0] and e.max() <= edges[-1]
    # the sampling distribution is piecewise uniform with 1/64 per bin: its CDF is the linear
    # interpolation of (edge_k, k / 64), which maps the sample to the uniform distribution
    chi2, p = _chi2_uniform(np.interp(e, edges, np.arange(65) / 64.0), 128)
    # mean = average of the bin mid points; variance = mean of (w^2/12 + (mid - mean)^2)
    mid = 0.5 * (edges[:-1] + edges[1:])
    w = np.diff(edges)
    mean = float(mid.mean())
    var = float(np.mean(w * w / 12.0 + (mid - mean) ** 2))
    z = (e.mean() - mean) / math.sqrt(var / N_DRAWS)
    OBSERVED["inverse cdf"] = p
    print(f"P2 inverse cdf: chi2={chi2:.1f} (127 dof) p={p:.4f} mean z={z:+.2f}")
    assert p > P_MIN
    assert abs(z) < 4.0


# ---------------------------------------------------------------------------------------------
# Function-level checks (twin)
# ---------------------------------------------------------------------------------------------
def test_kalbach_cdf_derivative_is_the_density_and_cdf_matches_closed_form() -> None:
    h = 1e-6
    for a, r in KALBACH_CASES + [(2.5, 0.3)]:
        mus = np.linspace(-0.99, 0.99, 25)
        for mu in mus:
            d = (NU64.kalbach_cdf(mu + h, a, r) - NU64.kalbach_cdf(mu - h, a, r)) / (2 * h)
            assert d == pytest.approx(NU64.kalbach_pdf(mu, a, r), rel=1e-6, abs=1e-9)
            direct = _kalbach_cdf_np(np.array([mu]), a, r)[0]
            assert NU64.kalbach_cdf(mu, a, r) == pytest.approx(direct, rel=1e-9, abs=1e-12)
        assert NU64.kalbach_cdf(-1.0, a, r) == 0.0
        assert NU64.kalbach_cdf(1.0, a, r) == pytest.approx(1.0, abs=1e-15)


def test_kalbach_mu_equals_closed_form_for_r0_and_inverts_the_cdf() -> None:
    us = (np.arange(2000) + 0.5) / 2000.0
    for a in (0.01, 0.3, 1.0, 2.0, 5.0, 12.0):
        closed = np.arcsinh((2 * us - 1) * np.sinh(a)) / a
        mu0 = np.array([NU64.kalbach_mu(u, a, 0.0) for u in us])
        assert np.max(np.abs(mu0 - closed)) < 1e-12
        for r in (0.0, 0.5, 1.0, -0.5):
            mu = np.array([NU64.kalbach_mu(u, a, r) for u in us])
            assert np.max(np.abs(_kalbach_cdf_np(mu, a, r) - us)) < 2e-13
            assert np.all(np.diff(mu) >= 0.0)  # the inverse CDF is monotone
    # extreme uniforms of the Philox map (float32 limit 2^-24) and the uniform branch a < 1e-6
    for a, r in ((5.0, 1.0), (5.0, 0.0), (0.01, 1.0)):
        for u in (2.0**-25, 1.0 - 2.0**-25):
            mu = NU64.kalbach_mu(u, a, r)
            assert abs(_kalbach_cdf_np(np.array([mu]), a, r)[0] - u) < 1e-12
    assert NU64.kalbach_mu(0.25, 0.0, 0.7) == -0.5 and NU64.kalbach_mu(0.25, -1.0, 0.7) == -0.5


def test_kalbach_a_formula_constants_and_range() -> None:
    # hand evaluation: e_a = 60, e_b = 20 MeV, m_b = 1: X1 = min(60, 130) 20 / 60 = 20,
    # X3 = min(60, 41) 20 / 60 = 41/3: a = 0.04 20 + 1.8e-6 20^3 + 6.7e-7 (41/3)^4
    x3 = 41.0 / 3.0
    assert NU64.kalbach_a(60.0, 20.0, 1.0) == pytest.approx(
        0.8 + 1.8e-6 * 8000.0 + 6.7e-7 * x3**4, rel=1e-14
    )
    # e_a above both thresholds: X1 = 130 e_b / e_a, X3 = 41 e_b / e_a; alpha has m_b = 2
    x1, x3 = 130.0 * 90.0 / 200.0, 41.0 * 90.0 / 200.0
    assert NU64.kalbach_a(200.0, 90.0, 2.0) == pytest.approx(
        0.04 * x1 + 1.8e-6 * x1**3 + 6.7e-7 * 2.0 * x3**4, rel=1e-14
    )
    # e_a <= E_t3: X1 = X3 = e_b, so a is independent of e_a
    assert NU64.kalbach_a(30.0, 10.0, 0.5) == NU64.kalbach_a(20.0, 10.0, 0.5)
    # for every kinematically possible pair e_b <= e_a <= 400 MeV a stays far below A_MAX = 20
    ea, eb = np.meshgrid(np.linspace(1.0, 400.0, 120), np.linspace(0.0, 1.0, 40))
    amax = max(
        NU64.kalbach_a(float(x), float(y * x), mb)
        for x, y in zip(ea.ravel(), eb.ravel(), strict=True)
        for mb in (0.5, 1.0, 2.0)
    )
    assert 0.0 < amax < 14.0
    assert KALBACH_M_B == {"n": 0.5, "p": 1.0, "d": 1.0, "t": 1.0, "he3": 1.0, "a": 2.0}


def test_kalbach_separation_energy_formula() -> None:
    """Hand evaluation of the liquid-drop form for p + 12C (compound 13N = (7, 13), target
    (6, 12)); the six terms are evaluated separately."""
    zc, ac, zr, ar = 7, 13, 6, 12
    nc, nr = ac - zc, ar - zr
    terms = [
        15.68 * (ac - ar),
        -28.07 * ((nc - zc) ** 2 / ac - (nr - zr) ** 2 / ar),
        -18.56 * (ac ** (2 / 3) - ar ** (2 / 3)),
        33.22 * ((nc - zc) ** 2 / ac ** (4 / 3) - (nr - zr) ** 2 / ar ** (4 / 3)),
        -0.717 * (zc * zc / ac ** (1 / 3) - zr * zr / ar ** (1 / 3)),
        1.211 * (zc * zc / ac - zr * zr / ar),
    ]
    assert kalbach_separation_energy(zc, ac, zr, ar, KALBACH_I["p"]) == pytest.approx(
        sum(terms), rel=1e-13
    )
    # the binding correction is subtracted
    assert kalbach_separation_energy(zc, ac, zr, ar, 2.22) == pytest.approx(sum(terms) - 2.22)
    # liquid-drop accuracy only: the proton separation energy of 13N is 1.94 MeV (measured); the
    # formula gives a value within a few MeV of it (documented limitation, see the report)
    assert abs(sum(terms) - 1.943) < 6.0


def test_thinning_select_target_and_step_limit() -> None:
    assert NU64.thinning_accept(0.3, 0.5, 1.0) == (1, 0)
    assert NU64.thinning_accept(0.5, 0.5, 1.0) == (0, 0)  # u S^ = S(E1): not strictly smaller
    assert NU64.thinning_accept(0.99, 1.2, 1.0) == (1, 1)  # violation: always accepted, flagged
    assert NU64.thinning_accept(0.5, 0.0, 0.0) == (0, 0)
    assert NU64.thinning_accept(0.5, 0.1, 0.0) == (0, 1)
    cum = [0.2, 0.5, 0.8, 1.0 - 1e-17]
    cum = [0.2, 0.5, 0.8, 0.9999999999999999]
    for u, want in ((0.1, 0), (0.2, 1), (0.49, 1), (0.5, 2), (0.79, 2), (0.8, 3), (0.99999999, 3)):
        cur = -1
        for k, c in enumerate(cum):
            cur = NU64.select_target(u, c, k, cur, int(k == len(cum) - 1))
        assert cur == want
    # a row summing to less than u still selects the last element
    cur = -1
    for k, c in enumerate([0.1, 0.2]):
        cur = NU64.select_target(0.9999, c, k, cur, int(k == 1))
    assert cur == 1
    # d = 10 n / (rho S^) [mm]: n = 1, rho = 1 g/cm3, S^ = 0.005 cm2/g -> 200 cm = 2000 mm
    assert NU64.nuclear_step_limit(1.0, 1.0, 0.005) == pytest.approx(2000.0)
    assert NU64.nuclear_step_limit(2.0, 0.0, 0.005) == 1.0e30
    assert NU64.nuclear_step_limit(2.0, 1.0, 0.0) == 1.0e30


def test_poisson_cap_and_exact_small_cases() -> None:
    lam = 1.0
    cdf = np.cumsum([math.exp(-lam) * lam**k / math.factorial(k) for k in range(8)])
    for k in range(7):
        assert NU64.poisson_inverse(cdf[k] * (1 - 1e-9), lam, 16) == k
        assert NU64.poisson_inverse(cdf[k] * (1 + 1e-9), lam, 16) == k + 1
    assert NU64.poisson_inverse(1.0 - 1e-15, 4.0, 5) == 5  # n_max caps the multiplicity
    assert NU64.poisson_inverse(1.0 - 1e-15, 4.0, 99) <= 16
    # probability mass above the cap of 16 at the largest intended lambda
    tail = 1.0 - sum(math.exp(-4.0) * 4.0**k / math.factorial(k) for k in range(16))
    assert tail < 1e-5


def test_residual_and_boost_functions() -> None:
    assert NU64.residual_invariant_mass(10.0, 0.0, 0.0, 0.0) == 10.0
    assert NU64.residual_invariant_mass(5.0, 3.0, 0.0, 0.0) == pytest.approx(4.0, rel=1e-15)
    assert NU64.residual_invariant_mass(3.0, 0.0, 4.0, 0.0) == pytest.approx(-math.sqrt(7.0))
    assert NU64.residual_mass_ok(10.0, 10.0) == 1 and NU64.residual_mass_ok(9.999, 10.0) == 0
    beta, gamma, sqrt_s = NU64.cm_boost(100.0, M_P, M_C12)
    e_tot = 100.0 + M_P + M_C12
    p = math.sqrt(100.0 * (100.0 + 2 * M_P))
    assert sqrt_s == pytest.approx(math.sqrt(e_tot**2 - p**2), rel=1e-14)
    assert beta == pytest.approx(p / e_tot, rel=1e-14)
    assert gamma == pytest.approx(1.0 / math.sqrt(1.0 - beta**2), rel=1e-13)
    # the boost of the CM four-momentum (sqrt s, 0) gives the lab total (E, p)
    e_b, pz_b = NU64.boost_z(sqrt_s, 0.0, beta, gamma)
    assert e_b == pytest.approx(e_tot, rel=1e-14) and pz_b == pytest.approx(p, rel=1e-14)
    # cm_to_lab: a particle at rest in the CM moves with the CM; direction/azimuth conventions
    e, px, py, pz = NU64.cm_to_lab(0.0, 0.3, 1.0, M_P, beta, gamma)
    assert e == pytest.approx(gamma * M_P) and px == 0.0 and py == 0.0
    assert pz == pytest.approx(gamma * beta * M_P)
    e, px, py, pz = NU64.cm_to_lab(10.0, 0.0, 0.0, M_N, 0.0, 1.0)  # no boost, mu = 0, phi = 0
    pcm = math.sqrt(10.0 * (10.0 + 2 * M_N))
    assert (px, py, pz) == pytest.approx((pcm, 0.0, 0.0), abs=1e-12)
    assert e == pytest.approx(10.0 + M_N)


# ---------------------------------------------------------------------------------------------
# select_step_nuclear
# ---------------------------------------------------------------------------------------------
def test_select_step_nuclear_reasons_ties_and_equality_with_select_step() -> None:
    from ionmc.transport.funcs import make_transport_funcs

    tf = python_twin(make_transport_funcs)
    big = 1.0e30
    # nuclear strictly smallest -> reason 4
    assert tf.select_step_nuclear(5.0, 4.0, 3.0, 2.0, 1.0) == (1.0, 4)
    # geometry wins ties with the nuclear limit and with the others
    assert tf.select_step_nuclear(1.0, 4.0, 3.0, 2.0, 1.0) == (1.0, 0)
    assert tf.select_step_nuclear(2.0, 4.0, 3.0, 2.0, 5.0) == (2.0, 0)
    # a tie between the nuclear limit and an EM limit goes to the EM limit (strict <)
    assert tf.select_step_nuclear(9.0, 4.0, 3.0, 2.0, 2.0) == (2.0, 3)
    assert tf.select_step_nuclear(9.0, 4.0, 3.0, 5.0, 3.0) == (3.0, 2)
    assert tf.select_step_nuclear(9.0, 3.0, 4.0, 5.0, 3.0) == (3.0, 1)
    # larger nuclear limit: the EM reasons 1, 2, 3 are unchanged
    assert tf.select_step_nuclear(9.0, 1.0, 3.0, 5.0, 7.0) == (1.0, 1)
    # +inf and BIG_LENGTH nuclear limits reproduce select_step exactly (random inputs)
    rng = np.random.default_rng(5)
    for _ in range(2000):
        d, s1, s2, s3 = (float(v) for v in rng.choice([0.5, 1.0, 2.0, 3.0], 4))
        for d_nuc in (math.inf, big):
            assert tf.select_step_nuclear(d, s1, s2, s3, d_nuc) == tf.select_step(d, s1, s2, s3)
    # the unchanged select_step still returns reasons 0-3 only
    assert tf.select_step(9.0, 4.0, 3.0, 2.0) == (2.0, 3)


def test_select_step_source_is_unchanged() -> None:
    """``select_step`` keeps its exact body (decision 0041 section 8): reasons 0-3."""
    import inspect

    from ionmc.transport import funcs

    src = inspect.getsource(funcs.make_transport_funcs)
    body = src.split("def select_step(")[1].split("@named_func(name)")[0]
    assert "reason = 4" not in body
    assert "d_nuc" not in body


# ---------------------------------------------------------------------------------------------
# P5: twin against Warp CPU kernel for every function
# ---------------------------------------------------------------------------------------------
COVERED = {
    "nuclear_step_limit", "thinning_accept", "select_target", "poisson_inverse", "inv_cdf_bin",
    "multiplicity_round", "grid_locate", "inv_cdf_sample", "kalbach_a", "kalbach_cdf",
    "kalbach_pdf", "kalbach_mu",
    "residual_invariant_mass", "residual_mass_ok", "cm_boost", "boost_z", "cm_to_lab",
}  # fmt: skip

NXC = 31  # real argument columns
NIC = 5  # int argument columns
NOR = 16  # real outputs
NOIC = 7  # int outputs


def _p5_args(precision: str) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(20261007)
    n = N_ARGS
    dt = np.float64 if precision == "float64" else np.float32
    x = np.zeros((n, NXC))
    ii = np.zeros((n, NIC), dtype=np.int32)

    def logu(lo: float, hi: float) -> np.ndarray:
        return np.exp(rng.uniform(math.log(lo), math.log(hi), n))

    def uni() -> np.ndarray:
        w = rng.integers(0, 2**32, n, dtype=np.uint64)
        if precision == "float64":
            return ((w >> np.uint64(8)).astype(np.float64) + 0.5) * 2.0**-24
        return ((w >> np.uint64(9)).astype(np.float64) + 0.5) * 2.0**-23

    x[:, 0] = uni()  # u
    x[:, 1] = logu(1e-3, 8.0)  # a (the first rows are the edge values below)
    x[:, 2] = np.where(
        rng.random(n) < 0.5, rng.uniform(0.0, 1.0, n), rng.choice([0.0, 0.5, 1.0, -0.4, 0.97], n)
    )
    x[:8, 1] = [0.0, 5e-7, 7e-7, 1.5e-6, 20.0, 25.0, 0.01, 5.0]
    x[:8, 2] = [0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0]
    x[:, 3] = rng.exponential(1.0, n) + 1e-3  # n_lambda
    x[:, 4] = rng.uniform(0.0, 3.0, n)  # rho
    x[0, 4] = 0.0
    x[:, 5] = logu(1e-4, 1.0)  # sigma_hat
    x[1, 5] = 0.0
    x[:, 6] = x[:, 5] * rng.uniform(0.0, 1.15, n)  # sigma_e1 (some above the majorant)
    x[:, 7] = rng.uniform(0.0, 1.0, n)  # cum_k
    x[:, 8] = logu(0.05, 8.0)  # lambda
    x[:, 9] = uni()  # u_frac
    x[:, 10] = rng.uniform(0.0, 40.0, n)  # e_lo
    x[:, 11] = x[:, 10] + logu(1e-3, 5.0)  # e_hi
    x[:, 12] = rng.uniform(5.0, 300.0, n)  # e_a
    x[:, 13] = x[:, 12] * rng.uniform(0.0, 1.0, n)  # e_b
    x[:, 14] = rng.choice([0.5, 1.0, 2.0], n)  # m_b
    x[:, 15] = rng.uniform(-1.0, 1.0, n)  # mu for cdf / pdf / cm_to_lab
    x[:6, 15] = [-1.0, 1.0, 0.0, -0.999999, 0.999999, 1e-9]
    # residual: m, p log-uniform, e = sqrt(m^2 + p^2); a few rows below the light cone
    mm = logu(1.0, 3000.0)
    pp = logu(0.1, 1000.0)
    x[:, 16] = np.sqrt(mm * mm + pp * pp)
    x[:3, 16] = [1.0, 5.0, 3.0]
    x[:, 17] = pp * rng.uniform(-1.0, 1.0, n) / math.sqrt(3.0)
    x[:, 18] = pp * rng.uniform(-1.0, 1.0, n) / math.sqrt(3.0)
    x[:, 19] = pp * rng.uniform(-1.0, 1.0, n) / math.sqrt(3.0)
    x[:3, 17:20] = [[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]]
    x[:, 20] = logu(1.0, 250.0)  # T
    x[:, 21] = M_P
    x[:, 22] = rng.choice([M_C12, M_O16, 2.5e4, 3.7e4], n)  # m_t
    x[:, 23] = rng.uniform(0.0, 2.0 * math.pi, n)  # phi
    x[:, 24] = rng.choice([M_N, M_P, M_D, M_A], n)  # m
    x[:, 25] = logu(1e-2, 60.0)  # t_cm
    b = rng.uniform(0.0, 0.9, n)
    b[0] = 0.0
    x[:, 27] = b
    x[:, 28] = 1.0 / np.sqrt(1.0 - b * b)
    # residual_mass_ok: m_r (column 29) against M_r (column 30), with exact ties
    x[:, 29] = logu(1.0, 5000.0)
    x[:, 30] = x[:, 29] * rng.uniform(0.9, 1.1, n)
    x[:4, 30] = x[:4, 29]
    ii[:, 0] = rng.integers(0, 6, n)  # k
    ii[:, 1] = rng.choice([-1, -1, 0, 3], n)  # current
    ii[:, 2] = rng.integers(0, 2, n)  # is_last
    ii[:, 3] = rng.choice([0, 1, 3, 8, 16, 40], n)  # n_max
    ii[:, 4] = 64 if precision == "float32" else rng.choice([64, 50, 128, 7], n)
    return x.astype(dt), ii


def _p5_kernel(real: Any) -> Any:
    nu = make_nuclear(real)

    @wp.kernel(module="unique")
    def kernel(
        x: wp.array2d(dtype=real),  # type: ignore[valid-type]
        ii: wp.array2d(dtype=int),  # type: ignore[valid-type]
        o: wp.array2d(dtype=real),  # type: ignore[valid-type]
        oi: wp.array2d(dtype=int),  # type: ignore[valid-type]
    ) -> None:
        i = wp.tid()
        o[i, 0] = nu.nuclear_step_limit(x[i, 3], x[i, 4], x[i, 5])
        acc, vio = nu.thinning_accept(x[i, 0], x[i, 6], x[i, 5])
        oi[i, 0] = acc
        oi[i, 1] = vio
        oi[i, 2] = nu.select_target(x[i, 0], x[i, 7], ii[i, 0], ii[i, 1], ii[i, 2])
        oi[i, 3] = nu.poisson_inverse(x[i, 0], x[i, 8], ii[i, 3])
        oi[i, 4] = nu.inv_cdf_bin(x[i, 0], ii[i, 4])
        o[i, 1] = nu.inv_cdf_sample(x[i, 9], x[i, 10], x[i, 11])
        o[i, 2] = nu.kalbach_a(x[i, 12], x[i, 13], x[i, 14])
        o[i, 3] = nu.kalbach_cdf(x[i, 15], wp.max(x[i, 1], real(1e-3)), x[i, 2])
        o[i, 4] = nu.kalbach_pdf(x[i, 15], wp.max(x[i, 1], real(1e-3)), x[i, 2])
        o[i, 5] = nu.kalbach_mu(x[i, 0], x[i, 1], x[i, 2])
        o[i, 6] = nu.residual_invariant_mass(x[i, 16], x[i, 17], x[i, 18], x[i, 19])
        oi[i, 5] = nu.residual_mass_ok(x[i, 29], x[i, 30])
        oi[i, 6] = nu.multiplicity_round(x[i, 0], x[i, 8], ii[i, 3])
        bt, gm, ss = nu.cm_boost(x[i, 20], x[i, 21], x[i, 22])
        o[i, 7] = bt
        o[i, 8] = gm
        o[i, 9] = ss
        be, bp = nu.boost_z(x[i, 16], x[i, 19], x[i, 27], x[i, 28])
        o[i, 10] = be
        o[i, 11] = bp
        el, lx, ly, lz = nu.cm_to_lab(x[i, 25], x[i, 15], x[i, 23], x[i, 24], x[i, 27], x[i, 28])
        o[i, 12] = el
        o[i, 13] = lx
        o[i, 14] = ly
        o[i, 15] = lz

    return kernel


def _p5_twin(x: np.ndarray, ii: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nu = NU64
    n = x.shape[0]
    o = np.zeros((n, NOR))
    oi = np.zeros((n, NOIC), dtype=np.int64)
    for i in range(n):
        r = [float(c) for c in x[i]]
        k = [int(c) for c in ii[i]]
        o[i, 0] = nu.nuclear_step_limit(r[3], r[4], r[5])
        acc, vio = nu.thinning_accept(r[0], r[6], r[5])
        oi[i, 0], oi[i, 1] = acc, vio
        oi[i, 2] = nu.select_target(r[0], r[7], k[0], k[1], k[2])
        oi[i, 3] = nu.poisson_inverse(r[0], r[8], k[3])
        oi[i, 4] = nu.inv_cdf_bin(r[0], k[4])
        o[i, 1] = nu.inv_cdf_sample(r[9], r[10], r[11])
        o[i, 2] = nu.kalbach_a(r[12], r[13], r[14])
        o[i, 3] = nu.kalbach_cdf(r[15], max(r[1], 1e-3), r[2])
        o[i, 4] = nu.kalbach_pdf(r[15], max(r[1], 1e-3), r[2])
        o[i, 5] = nu.kalbach_mu(r[0], r[1], r[2])
        o[i, 6] = nu.residual_invariant_mass(r[16], r[17], r[18], r[19])
        oi[i, 5] = nu.residual_mass_ok(r[29], r[30])
        oi[i, 6] = nu.multiplicity_round(r[0], r[8], k[3])
        o[i, 7], o[i, 8], o[i, 9] = nu.cm_boost(r[20], r[21], r[22])
        o[i, 10], o[i, 11] = nu.boost_z(r[16], r[19], r[27], r[28])
        o[i, 12:16] = nu.cm_to_lab(r[25], r[15], r[23], r[24], r[27], r[28])
    return o, oi


NAMES_REAL = {
    0: "nuclear_step_limit", 1: "inv_cdf_sample", 2: "kalbach_a", 3: "kalbach_cdf",
    4: "kalbach_pdf", 5: "kalbach_mu", 6: "residual_invariant_mass", 7: "cm_boost_beta",
    8: "cm_boost_gamma", 9: "cm_boost_sqrt_s", 10: "boost_z_e", 11: "boost_z_pz",
    12: "cm_to_lab_e", 13: "cm_to_lab_px", 14: "cm_to_lab_py", 15: "cm_to_lab_pz",
}  # fmt: skip
NAMES_INT = {
    0: "thinning_accepted", 1: "thinning_violation", 2: "select_target", 3: "poisson_inverse",
    4: "inv_cdf_bin", 5: "residual_mass_ok", 6: "multiplicity_round",
}  # fmt: skip


def _float32_budgets(x: np.ndarray) -> dict[int, np.ndarray]:
    """Absolute float32 error budgets (kernel float32 against twin float64 on the float32
    arguments), per real output column; all derived from the unit roundoff ``U = 2^-24`` and the
    number of roundings, never from observed differences. Columns not listed: relative 16.8 U
    (the frozen U1 default 1e-6) plus 1e-6 absolute."""
    out: dict[int, np.ndarray] = {}
    a = np.maximum(x[:, 1], 1e-3)
    mu = x[:, 15]
    r = x[:, 2]
    # kalbach_a = c1 x1 + c2 x1^3 + c3 m x3^4 with float32 constants (each of relative error U
    # against the double constants of the twin): x1, x3 two roundings, the powers 3 and 4 more, the
    # sum two: every term within 8 U, constants add U: 9 U + the sum 2 U < 12 U relative
    # CDF F = sinh p (cosh q - r sinh q) / sinh a, p, q = a (1 +- mu) / 2: p and q carry 3 roundings
    # (relative 2 U, plus the rounding of 1 +- mu, absolute U); sinh/cosh of an argument of
    # relative error eps have relative error eps |arg| coth/tanh (<= eps (1 + |arg|)) plus 2 U for
    # the library function; the difference cosh q - r sinh q has absolute error
    # (4 + 2 q) U cosh q, the product and quotient 3 U: |dF| <= (12 + 4 a) U Fs,
    # Fs = sinh p (cosh q + |r| sinh q) / sinh a (magnitudes of the terms)
    p_ = 0.5 * a * (1.0 + mu)
    q_ = 0.5 * a * (1.0 - mu)
    fs = np.sinh(p_) * (np.cosh(q_) + np.abs(r) * np.sinh(q_)) / np.sinh(a)
    out[3] = (12.0 + 4.0 * a) * U * fs + 1e-7
    ps = a * (np.cosh(a * mu) + np.abs(r) * np.abs(np.sinh(a * mu))) / (2.0 * np.sinh(a))
    out[4] = (12.0 + 4.0 * a) * U * ps + 1e-7
    # kalbach_mu: the solution of F(mu) = u from a float32 evaluation of F with |dF| as above, so
    # |d mu| <= 2 |dF| / f(mu) (factor 2: the last Newton iterate may sit one step off the root),
    # plus the step tolerance 2e-6 of the float32 iteration and the float32 resolution of mu (4 U)
    mu_s = np.array(
        [
            NU64.kalbach_mu(float(u), float(aa), float(rr))
            for u, aa, rr in zip(x[:, 0], x[:, 1], x[:, 2], strict=True)
        ]
    )
    aa_ = np.clip(x[:, 1], 1e-6, 20.0)
    p2, q2 = 0.5 * aa_ * (1.0 + mu_s), 0.5 * aa_ * (1.0 - mu_s)
    fs2 = np.sinh(p2) * (np.cosh(q2) + np.abs(r) * np.sinh(q2)) / np.sinh(aa_)
    pdf2 = aa_ * (np.cosh(aa_ * mu_s) + r * np.sinh(aa_ * mu_s)) / (2.0 * np.sinh(aa_))
    big_a = x[:, 1] >= 1e-6
    out[5] = np.where(
        big_a,
        2.0 * (12.0 + 4.0 * aa_) * U * fs2 / np.maximum(pdf2, 1e-30) + 2e-6 + 4.0 * U,
        4.0 * U,
    )
    # residual_invariant_mass = sqrt((e - p)(e + p)): p has 4 roundings (3 U p), e - p one (U e),
    # so |d m2| <= (e + p)(U e + 4 U p) + 2 U |m2| <= 8 U (e + p)^2; |d m| from sqrt
    e, p3 = x[:, 16], np.sqrt(x[:, 17] ** 2 + x[:, 18] ** 2 + x[:, 19] ** 2)
    m2 = (e - p3) * (e + p3)
    dm2 = 8.0 * U * (e + p3) ** 2 + 2.0 * U * np.abs(m2)
    out[6] = 2.0 * (np.sqrt(np.abs(m2) + dm2) - np.sqrt(np.maximum(np.abs(m2) - dm2, 0.0)))
    # inv_cdf_sample = lo + u (hi - lo): three roundings of terms of at most max(|lo|, |hi|)
    out[1] = 4.0 * U * (np.abs(x[:, 10]) + np.abs(x[:, 11])) + 1e-7
    # boost_z: gamma (e + beta pz) and gamma (pz + beta e): products, sum, product: 4 U of the
    # sum of the magnitudes (cancellation of the terms is the only ill-conditioning)
    g, bt = x[:, 28], x[:, 27]
    out[10] = 6.0 * U * g * (np.abs(x[:, 16]) + bt * np.abs(x[:, 19])) + 1e-7
    out[11] = 6.0 * U * g * (np.abs(x[:, 19]) + bt * np.abs(x[:, 16])) + 1e-7
    # cm_to_lab: p (3 U), st = sqrt(1 - mu^2) (abs error U / (2 st) + U st, at most sqrt(U) near 0),
    # cos / sin phi (abs 3 U): px, py within p (dst + 8 U); e_lab, pz_lab like boost_z with p |mu|
    pcm = np.sqrt(x[:, 25] * (x[:, 25] + 2.0 * x[:, 24]))
    st = np.sqrt(np.maximum(1.0 - mu * mu, 0.0))
    dst = np.where(st < math.sqrt(U), math.sqrt(U), U / (2.0 * np.maximum(st, 1e-30)) + U)
    out[13] = pcm * (dst + 8.0 * U) + 1e-7
    out[14] = out[13]
    ecm = x[:, 25] + x[:, 24]
    out[12] = 10.0 * U * g * (ecm + bt * pcm * np.abs(mu)) + 1e-7
    out[15] = 10.0 * U * g * (pcm * np.abs(mu) + bt * ecm) + 1e-7
    return out


def _poisson_near(x: np.ndarray, n_twin: np.ndarray) -> np.ndarray:
    """Rows whose uniform lies within the float32 error of a Poisson CDF value (where the discrete
    output may legitimately differ between precisions): relative ``8 (k + 2) U`` of the CDF."""
    lam, u = x[:, 8], x[:, 0]
    p = np.exp(-lam)
    cdf = p.copy()
    near = np.zeros(len(u), dtype=bool)
    near |= np.abs(u - cdf) <= 8.0 * 2.0 * U * cdf
    for k in range(1, 16):
        p = p * lam / k
        cdf = cdf + p
        near |= np.abs(u - cdf) <= 8.0 * (k + 2) * U * cdf
    return near


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_p5_nuclear_twin_equals_warp_cpu_kernel(precision: str) -> None:
    real = wp.float64 if precision == "float64" else wp.float32
    x, ii = _p5_args(precision)
    n = x.shape[0]
    kernel = _p5_kernel(real)
    ox = wp.zeros((n, NOR), dtype=real, device="cpu")
    oi = wp.zeros((n, NOIC), dtype=int, device="cpu")
    wp.launch(
        kernel,
        dim=n,
        inputs=[
            wp.array(x, dtype=real, device="cpu"),
            wp.array(ii, dtype=int, device="cpu"),
            ox,
            oi,
        ],
        device="cpu",
    )
    k_o, k_oi = ox.numpy().astype(np.float64), oi.numpy()
    p_o, p_oi = _p5_twin(x, ii)
    xx = x.astype(np.float64)
    assert COVERED == set(NU64.__dict__) - {"real"}, "every make_nuclear function is exercised"
    worst: dict[str, float] = {}
    budgets = _float32_budgets(xx) if precision == "float32" else {}
    for col, name in NAMES_REAL.items():
        a, b = p_o[:, col], k_o[:, col]
        assert np.all(np.isfinite(a)) and np.all(np.isfinite(b)), f"{name}: NaN or inf"
        if precision == "float64":
            tol_abs = 1e-13 if name in ("kalbach_mu",) else 1e-14
            bad = ~(np.abs(a - b) <= tol_abs + 1e-12 * np.abs(b))
            worst[name] = float(np.max(np.abs(a - b) / (tol_abs + np.abs(b))))
        else:
            tol = budgets.get(col, 16.8 * U * np.abs(b) + 1e-6)
            if col in budgets:
                tol = budgets[col]
            bad = ~(np.abs(a - b) <= tol)
            worst[name] = float(np.max(np.abs(a - b) / np.where(tol > 0, tol, 1.0)))
        assert not bad.any(), (
            f"{name}: {bad.sum()} of {n} mismatches; first twin={a[bad][:3]} kernel={b[bad][:3]}"
        )
    near_pois = _poisson_near(xx, p_oi[:, 3]) if precision == "float32" else np.zeros(n, bool)
    sig = np.abs(xx[:, 0] * xx[:, 5] - xx[:, 6])
    near_thin = (
        sig <= 4.0 * U * np.maximum(xx[:, 6], 1e-30)
        if precision == "float32"
        else np.zeros(n, bool)
    )
    for col, name in NAMES_INT.items():
        skip = np.zeros(n, bool)
        if name == "poisson_inverse":
            skip = near_pois
        if name == "thinning_accepted":
            skip = near_thin
        assert np.array_equal(p_oi[~skip, col], k_oi[~skip, col]), f"{name}: integer outputs differ"
    if precision == "float32":
        assert (~near_pois).mean() > 0.99
    print(f"P5 {precision}: max |twin - kernel| / budget per function: {worst}")
    if precision == "float64":
        assert max(worst.values()) <= 1.0


@wp.kernel(module="unique")
def _k_locate64(
    out: wp.array(dtype=int),  # type: ignore[valid-type]
    e: wp.array(dtype=wp.float64),  # type: ignore[valid-type]
    grid: wp.array(dtype=wp.float64),  # type: ignore[valid-type]
    n: int,
) -> None:
    i = wp.tid()
    out[i] = _NU.grid_locate(e[i], grid, n)


@wp.kernel(module="unique")
def _k_locate32(
    out: wp.array(dtype=int),  # type: ignore[valid-type]
    e: wp.array(dtype=wp.float32),  # type: ignore[valid-type]
    grid: wp.array(dtype=wp.float32),  # type: ignore[valid-type]
    n: int,
) -> None:
    i = wp.tid()
    out[i] = _NU32.grid_locate(e[i], grid, n)


@pytest.mark.parametrize("precision", ["float64", "float32"])
@pytest.mark.parametrize("n", [2, 3, 17, 64, 604])
def test_grid_locate_twin_kernel_and_searchsorted(precision: str, n: int) -> None:
    """P5 for ``grid_locate``: the twin, the Warp CPU kernel and ``searchsorted(side='right') - 1``
    (clamped to [0, n - 2]) agree exactly on random energies, every node and the clamps."""
    dt = np.float64 if precision == "float64" else np.float32
    rng = np.random.default_rng(20261007 + n)
    grid = np.unique(np.exp(rng.uniform(0.0, math.log(250.0), n)).astype(dt))
    grid[0], grid[-1] = dt(1.0), dt(250.0)
    n_g = grid.size
    e = np.concatenate(
        [
            np.exp(rng.uniform(math.log(0.5), math.log(300.0), 3000)).astype(dt),
            grid,
            np.nextafter(grid, dt(0.0)),
            np.nextafter(grid, dt(1e9)),
            np.array([0.0, 0.5, 1.0, 249.999, 250.0, 251.0, 1e9], dtype=dt),
        ]
    )
    ref = np.clip(np.searchsorted(grid, e, side="right") - 1, 0, n_g - 2)
    twin = np.array([NU64.grid_locate(float(v), grid.astype(np.float64), n_g) for v in e])
    assert np.array_equal(
        twin,
        np.clip(
            np.searchsorted(grid.astype(np.float64), e.astype(np.float64), side="right") - 1,
            0,
            n_g - 2,
        ),
    )
    out = wp.zeros(e.size, dtype=int, device="cpu")
    kernel = _k_locate64 if precision == "float64" else _k_locate32
    real = wp.float64 if precision == "float64" else wp.float32
    wp.launch(
        kernel,
        dim=e.size,
        inputs=[
            out,
            wp.array(e, dtype=real, device="cpu"),
            wp.array(grid, dtype=real, device="cpu"),
            n_g,
        ],
        device="cpu",
    )
    assert np.array_equal(out.numpy(), ref)
