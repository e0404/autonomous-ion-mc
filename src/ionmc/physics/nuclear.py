# mypy: ignore-errors
# (Warp function sources use runtime precision types in annotations; see _wpfunc.py.)
"""Proton non-elastic nuclear interactions: shared sampling, kinematics and bookkeeping functions.

Decision 0041 (sections 2-4). ``make_nuclear(R)`` returns precision-generic Warp functions for
kernels (``R`` = ``wp.float32`` or ``wp.float64``); the Python reference backend uses their
pure-Python twins (``ionmc._wpfunc.python_twin(make_nuclear)``: the same source text, no Warp
call at Python scope). The functions are pure: they never index arrays and never draw random
numbers (uniforms and table values are arguments; the backend reads the table rows and the
Philox blocks of ``PURPOSE_NUCLEAR`` in the layout of decision 0041 section 4). Units: energies,
masses and momenta in MeV (c = 1), lengths in mm, mass cross sections in cm2/g, densities in
g/cm3.

Functions
---------

``nuclear_step_limit``      distance [mm] to the next candidate interaction, ``10 n / (rho S^)``
``thinning_accept``         majorant thinning: accept with probability ``S(E1) / S^(E0)``
``select_target``           inverse-CDF choice of the target element from cumulative fractions
``poisson_inverse``         Poisson multiplicity by CDF inversion (``n <= 16``; test reference)
``multiplicity_round``      floor + Bernoulli multiplicity ``floor(lam) + [u < lam - floor(lam)]``
                            (the sampler's multiplicity, decision 0041 section 3 amended)
``grid_locate``             fixed-step binary search ``k`` with ``G_k <= e < G_{k+1}`` on the
                            non-uniform union grid of the nuclear table
``inv_cdf_bin`` / ``inv_cdf_sample``  E' by inverse CDF of a table of equiprobable bins
``kalbach_a``               Kalbach (1988) slope ``a(e_a, e_b)`` of the angular distribution
``kalbach_cdf`` / ``kalbach_pdf`` / ``kalbach_mu``  Kalbach-Mann angular distribution and its
                            exact inverse sampling
``residual_invariant_mass`` / ``residual_mass_ok``  residual four-momentum mass check (no longer
                            used by the sampler since the 2026-10-07 amendment; kept, with
                            their twin and kernel tests, as shared kinematics helpers)
``cm_boost`` / ``boost_z`` / ``cm_to_lab``  p + target centre-of-mass boost

Kalbach-Mann angular distribution (ENDF-6 LAW=1, LANG=2, NA=1)
--------------------------------------------------------------

With ``mu`` the cosine of the emission angle of the ejectile in the centre of mass relative to
the incident direction::

    f(mu | a, r) = a / (2 sinh a) * [cosh(a mu) + r sinh(a mu)],     mu in [-1, 1],

``r`` the pre-compound fraction tabulated in MF6 (the LIST row of NA=1 is ``E', f_0, r``: ``f_0``
is the energy distribution, ``r = b_1``) and ``a = a(e_a, e_b)`` from the Kalbach (1988)
systematics (S. Kalbach, Phys. Rev. C 37 (1988) 2350; ENDF-6 Formats Manual BNL-203218-2018-INRE,
section 6.2.3.2 -- the transcription of the constants and of the separation-energy formula below
is from the published form and is TO BE CROSS-CHECKED AGAINST THE MANUAL PDF)::

    a = C1 X1 + C2 X1^3 + C3 M_a m_b X3^4
    X1 = min(e_a, E_t1) e_b / e_a,     X3 = min(e_a, E_t3) e_b / e_a
    e_a = E_cm,a + S_a,     e_b = E_cm,b + S_b

    C1 = 0.04 MeV^-1,  C2 = 1.8e-6 MeV^-3,  C3 = 6.7e-7 MeV^-4,  E_t1 = 130 MeV,  E_t3 = 41 MeV

``E_cm,a`` is the incident-channel energy in the centre of mass (``E A_t / (A_t + 1)`` for a
proton, non-relativistic), ``E_cm,b`` the ejectile energy ``E'`` of the table (LCT=3: taken as
the p + target CM, decision 0041 section 3), ``S_a`` and ``S_b`` the separation energies of the
incident particle and of the ejectile from the compound nucleus (:func:`kalbach_separation_energy`).
``M_a`` (incident proton) = 1; ``m_b`` = 1/2 (n), 1 (p, d, t, 3He), 2 (alpha) (``KALBACH_M_B``).
The CDF is ``F(mu) = [sinh(a mu) + sinh a + r (cosh(a mu) - cosh a)] / (2 sinh a)``; with
``p = a (1 + mu) / 2`` and ``q = a (1 - mu) / 2`` the identical cancellation-free form
``F = sinh p (cosh q - r sinh q) / sinh a`` is evaluated.

Sampling ``F(mu) = u`` (:func:`kalbach_mu`): the initial value is the exact inverse of the
``r = 0`` case, ``mu_0 = asinh((2u - 1) sinh a) / a``; then a guarded Newton iteration on
``F(mu) - u`` inside the bracket ``[lo, hi]`` (initially ``[-1, 1]``, tightened by the sign of
``F - u`` at every iterate; a Newton step outside the bracket is replaced by the bisection step),
a FIXED maximum of 40 iterations and a step tolerance ``tol`` (``1e-14`` float64, ``2e-6``
float32) after which no further iteration changes ``mu``. The algorithm, the constants and the
order of operations are identical in the twin and in the kernels. For ``r = 0`` the first
Newton step is zero up to rounding, so the sample equals the closed form to rounding. For
``a < A_MIN = 1e-6`` the distribution is taken as uniform (error of the CDF at most
``a r / 4 < 2.5e-7``); ``a`` is limited to ``A_MAX = 20`` so that ``cosh a`` stays finite (a is
about 0.05-3 for 1-250 MeV protons; tested).
"""

import functools
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func

wp.set_module_options({"enable_backward": False})

BIG_LENGTH_MM = 1.0e30
"""Nuclear step limit returned for a vanishing interaction rate (as ``transport.funcs``)."""
POISSON_N_MAX = 16
"""Largest Poisson multiplicity (decision 0041 section 3): probability mass above the cap is
assigned to the cap (negligible for lambda <= 4, see ``tests/ionmc/test_nuclear_funcs.py``)."""
KALBACH_NEWTON_ITERATIONS = 40
KALBACH_TOLERANCE = {"float64": 1.0e-14, "float32": 2.0e-6}
KALBACH_A_MIN = 1.0e-6
KALBACH_A_MAX = 20.0

KALBACH_C1 = 0.04  # MeV^-1
KALBACH_C2 = 1.8e-6  # MeV^-3
KALBACH_C3 = 6.7e-7  # MeV^-4
KALBACH_E_T1 = 130.0  # MeV
KALBACH_E_T3 = 41.0  # MeV
KALBACH_M_A_PROTON = 1.0
KALBACH_M_B = {"n": 0.5, "p": 1.0, "d": 1.0, "t": 1.0, "he3": 1.0, "a": 2.0}
"""Outgoing-particle factor ``m_b`` of the Kalbach (1988) systematics."""
KALBACH_I = {"n": 0.0, "p": 0.0, "d": 2.22, "t": 8.48, "he3": 7.72, "a": 28.3}
"""Binding-energy corrections ``I`` [MeV] of the separation-energy formula."""

# Liquid-drop coefficients of the separation-energy formula [MeV] (Kalbach 1988; ENDF-6 manual):
_SEP = (15.68, 28.07, 18.56, 33.22, 0.717, 1.211)


def kalbach_separation_energy(z_c: int, a_c: int, z_r: int, a_r: int, binding_i: float) -> float:
    """Separation energy [MeV] of a particle taking the compound nucleus ``(z_c, a_c)`` to the
    residual ``(z_r, a_r)`` in the liquid-drop form used by Kalbach (1988) and the ENDF-6 manual
    (to be cross-checked against the manual PDF)::

        S = 15.68 (A_C - A_R)
            - 28.07 [(N_C - Z_C)^2 / A_C - (N_R - Z_R)^2 / A_R]
            - 18.56 [A_C^(2/3) - A_R^(2/3)]
            + 33.22 [(N_C - Z_C)^2 / A_C^(4/3) - (N_R - Z_R)^2 / A_R^(4/3)]
            - 0.717 [Z_C^2 / A_C^(1/3) - Z_R^2 / A_R^(1/3)]
            + 1.211 [Z_C^2 / A_C - Z_R^2 / A_R]
            - I

    ``binding_i`` is ``KALBACH_I[particle]``. For the incident proton ``(z_r, a_r)`` is the target
    and the compound nucleus ``(z_t + 1, a_t + 1)``; for the ejectile it is the residual left
    after the emission. Build-time helper (plain Python floats)."""
    c1, c2, c3, c4, c5, c6 = _SEP
    n_c, n_r = a_c - z_c, a_r - z_r

    def t(zz: int, nn: int, aa: int) -> tuple[float, float, float, float, float]:
        return (
            (nn - zz) ** 2 / aa,
            aa ** (2.0 / 3.0),
            (nn - zz) ** 2 / aa ** (4.0 / 3.0),
            zz * zz / aa ** (1.0 / 3.0),
            zz * zz / aa,
        )

    tc, tr = t(z_c, n_c, a_c), t(z_r, n_r, a_r)
    return (
        c1 * (a_c - a_r)
        - c2 * (tc[0] - tr[0])
        - c3 * (tc[1] - tr[1])
        + c4 * (tc[2] - tr[2])
        - c5 * (tc[3] - tr[3])
        + c6 * (tc[4] - tr[4])
        - binding_i
    )


@functools.cache
def make_nuclear(real: type) -> SimpleNamespace:
    """Return the shared nuclear functions for precision ``real``."""
    name = check_real(real)
    big = wp.constant(real(BIG_LENGTH_MM))
    c1 = wp.constant(real(KALBACH_C1))
    c2 = wp.constant(real(KALBACH_C2))
    c3 = wp.constant(real(KALBACH_C3))
    et1 = wp.constant(real(KALBACH_E_T1))
    et3 = wp.constant(real(KALBACH_E_T3))
    m_a = wp.constant(real(KALBACH_M_A_PROTON))
    a_min = wp.constant(real(KALBACH_A_MIN))
    a_max = wp.constant(real(KALBACH_A_MAX))
    tol = wp.constant(real(KALBACH_TOLERANCE[name]))

    @named_func(name)
    def nuclear_step_limit(n_lambda: real, rho_g_cm3: real, sigma_hat_cm2_g: real) -> real:
        """Distance [mm] after which ``n_lambda`` remaining mean free paths are used up at the
        majorant ``sigma_hat_cm2_g`` [cm2/g] in a medium of density ``rho_g_cm3``:
        ``10 n_lambda / (rho S^)`` (the 10 converts cm to mm); ``BIG_LENGTH_MM`` for a vanishing
        rate."""
        rate = rho_g_cm3 * sigma_hat_cm2_g
        d = big
        if rate > real(0.0):
            d = real(10.0) * n_lambda / rate
        return d

    @named_func(name)
    def thinning_accept(u: real, sigma_e1: real, sigma_hat_e0: real) -> tuple[int, int]:
        """Majorant thinning: ``(accepted, violation)``. A candidate is accepted with probability
        ``S(E1) / S^(E0)``, i.e. iff ``u S^ < S(E1)`` (``u`` in (0, 1)); ``violation`` = 1 iff the
        majorant was violated, ``S(E1) > S^(E0)`` (the caller counts it and invalidates the
        result; the candidate is then accepted, ``u S^ < S^ < S(E1)``). A non-positive majorant
        accepts nothing."""
        accepted = int(0)
        violation = int(0)
        if sigma_hat_e0 > real(0.0):
            if u * sigma_hat_e0 < sigma_e1:
                accepted = 1
        if sigma_e1 > sigma_hat_e0:
            violation = 1
        return accepted, violation

    @named_func(name)
    def select_target(u: real, cum_k: real, k: int, current: int, is_last: int) -> int:
        """One step of the inverse-CDF choice of the target element from a row of cumulative
        fractions: the caller loops ``k = 0 .. n-1`` and passes the running result ``current``
        (-1 at the start). Returns ``k`` if nothing was chosen yet and ``u < cum_k`` (or, for the
        last element, always: guards a row that sums to 1 - rounding), else ``current``."""
        chosen = current
        if current < 0:
            if u < cum_k or is_last == 1:
                chosen = k
        return chosen

    @named_func(name)
    def poisson_inverse(u: real, lam: real, n_max: int) -> int:
        """Poisson(``lam``) variate by inversion of the CDF with the recurrence
        ``p_k = p_{k-1} lam / k`` in the real type, ``n = min(#{k : u > CDF(k-1)}, n_max,
        POISSON_N_MAX = 16)``. Float budget (analytic): ``p_k`` carries relative error
        ``(3 k + 2) u_r`` and the CDF at most ``(4 k + 4) u_r`` (``u_r`` the unit roundoff), so
        ``n`` can differ between precisions only when ``u`` lies within that relative distance
        of a CDF value. ``lam`` <= 20 (exp(-lam) must not underflow)."""
        p = wp.exp(-lam)
        cdf = p
        n = int(0)
        for _k in range(16):
            if u > cdf and n < n_max:
                n = n + 1
                p = p * lam / real(n)
                cdf = cdf + p
        return n

    @named_func(name)
    def multiplicity_round(u: real, lam: real, n_max: int) -> int:
        """Floor + Bernoulli multiplicity ``floor(lam) + [u < lam - floor(lam)]`` capped at
        ``min(n_max, POISSON_N_MAX = 16)`` (``u`` in (0, 1), ``lam`` >= 0): the minimum-variance
        integer variate of mean ``lam`` (decision 0041 section 3, amendment 2026-10-07).
        ``lam - floor(lam)`` is exact in either precision, so kernels and twin agree exactly."""
        base = wp.floor(lam)
        n = int(base)
        if u < lam - base:
            n = n + 1
        return wp.min(n, wp.min(n_max, 16))

    @named_func(name)
    def grid_locate(e: real, grid: wp.array(dtype=real), n: int) -> int:
        """Index ``k`` in ``[0, n - 2]`` with ``grid[k] <= e < grid[k + 1]`` of a strictly
        increasing ``grid`` of ``n >= 2`` nodes, clamped (``e < grid[0]`` gives 0,
        ``e >= grid[n-1]`` gives ``n - 2``). Bisection of ``[lo, hi] = [0, n - 1]`` with the
        invariant ``grid[lo] <= e`` or ``lo = 0`` and ``e < grid[hi]`` or ``hi = n - 1``:
        exactly ``ceil(log2(n - 1))`` halving steps are active whatever ``e`` is (the loop bound
        31 only fixes the compiled trip count), so the work is data independent."""
        lo = int(0)
        hi = n - 1
        for _step in range(31):
            if hi - lo > 1:
                mid = (lo + hi) // 2
                if grid[mid] <= e:
                    lo = mid
                else:
                    hi = mid
        return lo

    @named_func(name)
    def inv_cdf_bin(u: real, n_bins: int) -> int:
        """Equiprobable bin ``floor(u n_bins)`` (clamped to ``n_bins - 1``) of an inverse-CDF
        table; exact in any precision for ``n_bins = 64`` (a power of two)."""
        k = int(wp.floor(u * real(n_bins)))
        k = wp.max(wp.min(k, n_bins - 1), 0)
        return k

    @named_func(name)
    def inv_cdf_sample(u_frac: real, e_lo: real, e_hi: real) -> real:
        """Value inside the bin ``[e_lo, e_hi]`` chosen by :func:`inv_cdf_bin`: linear in the bin
        (a piecewise-constant density, decision 0041 section 3), ``e_lo + u_frac (e_hi - e_lo)``."""
        return e_lo + u_frac * (e_hi - e_lo)

    @named_func(name)
    def kalbach_a(e_a: real, e_b: real, m_b: real) -> real:
        """Kalbach (1988) slope ``a`` for an incident proton (module docstring). ``e_a`` and
        ``e_b`` are the channel energies INCLUDING the separation energies, ``e_a = E_cm,a + S_a``
        and ``e_b = E_cm,b + S_b`` [MeV]; ``m_b`` is ``KALBACH_M_B`` of the ejectile."""
        x1 = wp.min(e_a, et1) * e_b / e_a
        x3 = wp.min(e_a, et3) * e_b / e_a
        return c1 * x1 + c2 * x1 * x1 * x1 + c3 * m_a * m_b * x3 * x3 * x3 * x3

    @named_func(name)
    def kalbach_cdf(mu: real, a: real, r: real) -> real:
        """CDF of the Kalbach-Mann distribution, ``sinh p (cosh q - r sinh q) / sinh a`` with
        ``p = a (1 + mu) / 2``, ``q = a (1 - mu) / 2`` (module docstring); ``a > 0``."""
        p = real(0.5) * a * (real(1.0) + mu)
        q = real(0.5) * a * (real(1.0) - mu)
        return wp.sinh(p) * (wp.cosh(q) - r * wp.sinh(q)) / wp.sinh(a)

    @named_func(name)
    def kalbach_pdf(mu: real, a: real, r: real) -> real:
        """Density ``a [cosh(a mu) + r sinh(a mu)] / (2 sinh a)``; ``a > 0``."""
        return a * (wp.cosh(a * mu) + r * wp.sinh(a * mu)) / (real(2.0) * wp.sinh(a))

    @named_func(name)
    def kalbach_mu(u: real, a_in: real, r: real) -> real:
        """Cosine ``mu`` with ``F(mu) = u``, ``u`` in (0, 1): guarded Newton iteration with a
        bisection fallback (module docstring; 40 iterations at most)."""
        a = wp.min(wp.max(a_in, real(0.0)), a_max)
        mu = real(2.0) * u - real(1.0)
        lo = real(-1.0)
        hi = real(1.0)
        done = int(0)
        if a >= a_min:
            s = (real(2.0) * u - real(1.0)) * wp.sinh(a)
            t = wp.abs(s)
            g = wp.log(t + wp.sqrt(t * t + real(1.0))) / a
            if s < real(0.0):
                g = -g
            mu = wp.min(wp.max(g, real(-1.0)), real(1.0))
            for _it in range(40):
                if done == 0:
                    diff = kalbach_cdf(mu, a, r) - u
                    if diff > real(0.0):
                        hi = mu
                    else:
                        lo = mu
                    pdf = kalbach_pdf(mu, a, r)
                    nxt = real(0.5) * (lo + hi)
                    if pdf > real(0.0):
                        cand = mu - diff / pdf
                        if cand >= lo and cand <= hi:
                            nxt = cand
                    if wp.abs(nxt - mu) <= tol:
                        done = 1
                    mu = nxt
        return wp.min(wp.max(mu, real(-1.0)), real(1.0))

    @named_func(name)
    def residual_invariant_mass(e_r: real, px: real, py: real, pz: real) -> real:
        """Signed invariant mass [MeV] of the residual four-momentum ``(e_r, px, py, pz)``:
        ``sqrt((E - p)(E + p))``, negative (``-sqrt(p^2 - E^2)``) if ``E < p`` (unphysical: fails
        :func:`residual_mass_ok`)."""
        p = wp.sqrt(px * px + py * py + pz * pz)
        m2 = (e_r - p) * (e_r + p)
        m = wp.sqrt(wp.abs(m2))
        if m2 < real(0.0):
            m = -m
        return m

    @named_func(name)
    def residual_mass_ok(m_r: real, big_m_r: real) -> int:
        """1 iff the residual can exist: invariant mass ``m_r`` >= the ground-state mass
        ``big_m_r`` (excitation ``E* = m_r - M_r >= 0``; decision 0041 section 3)."""
        ok = int(0)
        if m_r >= big_m_r:
            ok = 1
        return ok

    @named_func(name)
    def cm_boost(t_mev: real, m_p: real, m_t: real) -> tuple[real, real, real]:
        """Velocity ``beta`` and Lorentz factor ``gamma`` of the p + target centre of mass along
        the incident direction and ``sqrt(s)`` [MeV] for a projectile of kinetic energy ``t_mev``
        and mass ``m_p`` on a target of mass ``m_t`` at rest."""
        e_lab = t_mev + m_p
        p_lab = wp.sqrt(t_mev * (t_mev + real(2.0) * m_p))
        e_tot = e_lab + m_t
        s = m_p * m_p + m_t * m_t + real(2.0) * e_lab * m_t
        sqrt_s = wp.sqrt(s)
        return p_lab / e_tot, e_tot / sqrt_s, sqrt_s

    @named_func(name)
    def boost_z(e: real, pz: real, beta: real, gamma: real) -> tuple[real, real]:
        """Boost ``(e, pz)`` along z with velocity ``beta`` (the transverse momentum is
        unchanged): ``E' = gamma (E + beta p_z)``, ``p_z' = gamma (p_z + beta E)``."""
        return gamma * (e + beta * pz), gamma * (pz + beta * e)

    @named_func(name)
    def cm_to_lab(
        t_cm: real, mu: real, phi: real, m: real, beta: real, gamma: real
    ) -> tuple[real, real, real, real]:
        """Lab four-momentum ``(E, px, py, pz)`` [MeV], z along the incident direction, of an
        ejectile of mass ``m`` with kinetic energy ``t_cm`` and direction ``(mu, phi)`` (polar
        cosine and azimuth about z) in the p + target centre of mass moving with ``beta``,
        ``gamma`` (:func:`cm_boost`). Relativistic: ``p = sqrt(T (T + 2 m))``."""
        p = wp.sqrt(t_cm * (t_cm + real(2.0) * m))
        st = wp.sqrt(wp.max(real(1.0) - mu * mu, real(0.0)))
        px = p * st * wp.cos(phi)
        py = p * st * wp.sin(phi)
        e_cm = t_cm + m
        e_lab = gamma * (e_cm + beta * p * mu)
        pz_lab = gamma * (p * mu + beta * e_cm)
        return e_lab, px, py, pz_lab

    return SimpleNamespace(
        nuclear_step_limit=nuclear_step_limit,
        thinning_accept=thinning_accept,
        select_target=select_target,
        poisson_inverse=poisson_inverse,
        multiplicity_round=multiplicity_round,
        grid_locate=grid_locate,
        inv_cdf_bin=inv_cdf_bin,
        inv_cdf_sample=inv_cdf_sample,
        kalbach_a=kalbach_a,
        kalbach_cdf=kalbach_cdf,
        kalbach_pdf=kalbach_pdf,
        kalbach_mu=kalbach_mu,
        residual_invariant_mass=residual_invariant_mass,
        residual_mass_ok=residual_mass_ok,
        cm_boost=cm_boost,
        boost_z=boost_z,
        cm_to_lab=cm_to_lab,
        real=name,
    )
