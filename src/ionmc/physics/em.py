# mypy: ignore-errors
# (Warp function sources use runtime precision types in annotations; see _wpfunc.py.)
"""Electromagnetic condensed-history physics as precision-generic Warp functions.

Models and sources (decision 0039; derivations in ``docs/physics/em-transport.md``):

* mean energy loss over a path of mass thickness ``t`` [g/cm2]: the CSDA range inversion
  ``E1 = Rinv(R(E0) - t)``, with the linear branch ``S(E0) t`` when ``t < f_short R(E0)``;
* straggling: Bohr variance ``(K/2)(Z/A) rho x z^2 Tmax (1/beta^2 - 1/2)``; Gamma with exactly the
  Bohr mean and variance for mean/sigma < 3, Gaussian clamped to ``[0, 2 mean]`` for
  mean/sigma >= 3 (mean preserved, variance reduced by at most about 0.5 %);
* multiple Coulomb scattering: the differential Moliere scattering power ``T_dM`` of
  B. Gottschalk, Med. Phys. 37 (2010) 352 (arXiv:0908.1413), ``E_s = 15.0 MeV``, applied
  as a two-dimensional Gaussian polar angle, with the rotation of ``G4ThreeVector::rotateUz``.

Units: energies MeV, lengths mm, mass thickness and ranges g/cm2, densities g/cm3,
scattering power rad^2/mm (projected angle), angles rad. Functions never index arrays or
draw random numbers; uniforms and table values are arguments.
"""

import functools
import math
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func
from ionmc.physics.kinematics import make_kinematics
from ionmc.physics.stopping import K_MEV_CM2_MOL

wp.set_module_options({"enable_backward": False})

E_S_MEV = 15.0
"""Scattering constant E_s of Gottschalk's differential Moliere power [MeV]."""
FDM_COEFFICIENTS = (0.5244, 0.1975, 0.2320, -0.0098)
"""Coefficients of f_dM (Gottschalk 2010, eq. fdM)."""


@functools.cache
def make_em(real: type) -> SimpleNamespace:
    """Return the electromagnetic physics functions for precision ``real``."""
    name = check_real(real)
    kin = make_kinematics(real)
    v3 = wp.types.vector(length=3, dtype=real)
    k_const = wp.constant(real(K_MEV_CM2_MOL))
    es_const = wp.constant(real(E_S_MEV))
    two_pi = wp.constant(real(2.0 * math.pi))
    pi_const = wp.constant(real(math.pi))
    c0 = wp.constant(real(FDM_COEFFICIENTS[0]))
    c1 = wp.constant(real(FDM_COEFFICIENTS[1]))
    c2 = wp.constant(real(FDM_COEFFICIENTS[2]))
    c3 = wp.constant(real(FDM_COEFFICIENTS[3]))

    @named_func(name)
    def csda_mean_loss(
        e0_mev: real,
        e_r1_mev: real,
        s_mass0: real,
        mass_thickness: real,
        r0_g_cm2: real,
        f_short: real,
    ) -> real:
        """Mean energy loss [MeV] over mass thickness ``t`` = ``mass_thickness`` [g/cm2].

        ``e_r1_mev`` is Rinv(R0 - t) (clamped to the table), ``s_mass0`` the mass stopping
        power at ``e0_mev`` [MeV cm2/g], ``r0_g_cm2`` = R(e0). Result in [0, e0].
        """
        loss = e0_mev - e_r1_mev
        if mass_thickness >= r0_g_cm2:
            loss = e0_mev
        elif mass_thickness < f_short * r0_g_cm2:
            loss = s_mass0 * mass_thickness
        loss = wp.max(loss, real(0.0))
        loss = wp.min(loss, e0_mev)
        return loss

    @named_func(name)
    def bohr_variance(
        e_mid_mev: real, m_mev: real, z: real, z_over_a: real, rho_g_cm3: real, s_mm: real
    ) -> real:
        """Bohr energy-loss variance [MeV^2] of a step of length ``s_mm`` at ``e_mid_mev``."""
        b2 = kin.beta2(e_mid_mev, m_mev)
        tm = kin.tmax_mev(e_mid_mev, m_mev)
        x_cm = s_mm / real(10.0)
        return (
            real(0.5)
            * k_const
            * z_over_a
            * rho_g_cm3
            * x_cm
            * z
            * z
            * tm
            * (real(1.0) / b2 - real(0.5))
        )

    @named_func(name)
    def straggle_attempt(
        mean_mev: real, var_mev2: real, u0: real, u1: real, u2: real, u3: real
    ) -> tuple[real, int]:
        """One sampling attempt of the energy loss [MeV] and an accepted flag (1 or 0).

        With ratio = mean / sigma, ``ratio >= 3`` is Gaussian clamped to ``[0, 2 mean]`` (always
        accepted; mean preserved, P(x < 0) = 0.13 %, variance reduced by at most about 0.5 %);
        otherwise Gamma (exact mean and Bohr variance) with shape ``k = ratio^2`` and scale
        ``sigma^2 / mean`` by Marsaglia-Tsang (shape ``k + 1`` and a ``u3^(1/k)`` factor for
        ``k < 1``), accepted with the Marsaglia-Tsang test. All four uniforms are arguments (one
        Philox block); the caller draws a new block after a rejection (at most 64 attempts).
        """
        loss = real(0.0)
        ok = int(0)
        if mean_mev <= real(0.0):
            ok = 1
        elif var_mev2 <= real(0.0):
            loss = mean_mev
            ok = 1
        else:
            sigma = wp.sqrt(var_mev2)
            ratio = mean_mev / sigma
            x = wp.sqrt(real(-2.0) * wp.log(u0)) * wp.cos(two_pi * u1)
            if ratio >= real(3.0):
                loss = wp.min(wp.max(mean_mev + sigma * x, real(0.0)), real(2.0) * mean_mev)
                ok = 1
            else:
                k = ratio * ratio
                a = k
                if k < real(1.0):
                    a = k + real(1.0)
                d = a - real(1.0) / real(3.0)
                c = real(1.0) / wp.sqrt(real(9.0) * d)
                v1 = real(1.0) + c * x
                if v1 > real(0.0):
                    v = v1 * v1 * v1
                    lhs = wp.log(u2)
                    rhs = real(0.5) * x * x + d - d * v + d * wp.log(v)
                    if lhs < rhs:
                        g = d * v
                        if k < real(1.0):
                            g = g * wp.pow(u3, real(1.0) / k)
                        loss = mean_mev / k * g
                        ok = 1
        return loss, ok

    @named_func(name)
    def straggle_attempt_gamma(
        mean_mev: real, var_mev2: real, u0: real, u1: real, u2: real, u3: real
    ) -> tuple[real, int]:
        """One sampling attempt of the energy loss [MeV] for the model ``bohr_gamma_v1``: a Gamma
        distribution with exactly the Bohr mean and variance for EVERY ratio mean / sigma (shape
        ``k = ratio^2``, scale ``sigma^2 / mean``; Marsaglia-Tsang, shape ``k + 1`` and a
        ``u3^(1/k)`` factor for ``k < 1``): positive, no clamp, no Gaussian branch. With a common
        scale ``theta = sigma^2 / mean`` along a path the sum of Gamma steps is exactly Gamma with
        the summed shape, so the whole loss distribution (not only its first two moments) is
        independent of the step length wherever ``theta`` varies slowly. Same inputs, outputs and
        rejection protocol as :func:`straggle_attempt`."""
        loss = real(0.0)
        ok = int(0)
        if mean_mev <= real(0.0):
            ok = 1
        elif var_mev2 <= real(0.0):
            loss = mean_mev
            ok = 1
        else:
            sigma = wp.sqrt(var_mev2)
            ratio = mean_mev / sigma
            x = wp.sqrt(real(-2.0) * wp.log(u0)) * wp.cos(two_pi * u1)
            k = ratio * ratio
            a = k
            if k < real(1.0):
                a = k + real(1.0)
            d = a - real(1.0) / real(3.0)
            c = real(1.0) / wp.sqrt(real(9.0) * d)
            v1 = real(1.0) + c * x
            if v1 > real(0.0):
                v = v1 * v1 * v1
                lhs = wp.log(u2)
                rhs = real(0.5) * x * x + d - d * v + d * wp.log(v)
                if lhs < rhs:
                    g = d * v
                    if k < real(1.0):
                        g = g * wp.pow(u3, real(1.0) / k)
                    loss = mean_mev / k * g
                    ok = 1
        return loss, ok

    @named_func(name)
    def fdm(pv_mev: real, p1v1_mev: real) -> real:
        """Gottschalk's f_dM(pv, p1v1) >= 0 (zero where the fit diverges, pv -> p1v1)."""
        ratio = pv_mev / p1v1_mev
        om = real(1.0) - ratio * ratio
        f = real(0.0)
        if om > real(0.0):
            l1 = wp.log10(om)
            l2 = wp.log10(pv_mev)
            f = c0 + c1 * l1 + c2 * l2 + c3 * l1 * l2
        return wp.max(f, real(0.0))

    @named_func(name)
    def scattering_power_dm(
        pv_mev: real, p1v1_mev: real, z: real, inv_rho_xs_cm2_g: real, rho_g_cm3: real
    ) -> real:
        """Projected differential Moliere scattering power T_dM [rad^2/mm].

        ``f_dM (z E_s / pv)^2 / X_S`` with ``1/X_S = rho * inv_rho_xs_cm2_g`` per cm,
        divided by 10 for per mm.
        """
        q = z * es_const / pv_mev
        return fdm(pv_mev, p1v1_mev) * q * q * rho_g_cm3 * inv_rho_xs_cm2_g / real(10.0)

    @named_func(name)
    def scattering_variance_birth(
        pv_mid_mev: real,
        pv_end_mev: real,
        p1v1_mev: real,
        z: real,
        inv_rho_xs_cm2_g: real,
        rho_g_cm3: real,
        s_mm: real,
    ) -> real:
        """Projected variance [rad^2] of the first step of a particle's life (birth step).

        On this step ``1 - (pv/p1v1)^2`` grows linearly from 0, so the logarithmic term of
        f_dM is replaced by the analytic average ``L1 = lg(1 - (pv_end/p1v1)^2) - 1/ln(10)``,
        which is exact only if ``1 - (pv/p1v1)^2`` grew linearly along the step (linearized
        birth-step approximation; under CSDA slowing the growth is slightly nonlinear, so a
        residual of order 1e-3 remains for steps up to 1 mm, bounded by the quadrature test). The
        smooth terms (``lg pv``) and the prefactor are taken at the midpoint energy; f_dM
        is clamped at 0 after averaging. ``pv_end_mev`` belongs to the end of the step.
        """
        om = real(1.0) - (pv_end_mev / p1v1_mev) * (pv_end_mev / p1v1_mev)
        var = real(0.0)
        if om > real(0.0):
            l1 = wp.log10(om) - real(0.4342944819032518)
            l2 = wp.log10(pv_mid_mev)
            f = wp.max(c0 + c1 * l1 + c2 * l2 + c3 * l1 * l2, real(0.0))
            q = z * es_const / pv_mid_mev
            var = f * q * q * rho_g_cm3 * inv_rho_xs_cm2_g / real(10.0) * s_mm
        return var

    @named_func(name)
    def polar_deflection(var_rad2: real, u_r: real) -> real:
        """Polar angle [rad] of a 2-D Gaussian with projected variance ``var_rad2``."""
        theta = wp.sqrt(real(-2.0) * var_rad2 * wp.log(u_r))
        return wp.min(theta, pi_const)

    @named_func(name)
    def rotate_dir(u: v3, theta: real, phi: real) -> v3:
        """Rotate unit vector ``u`` by polar angle ``theta`` and azimuth ``phi`` (rotateUz)."""
        ct = wp.cos(theta)
        st = wp.sin(theta)
        px = st * wp.cos(phi)
        py = st * wp.sin(phi)
        up2 = u[0] * u[0] + u[1] * u[1]
        out = v3(real(0.0), real(0.0), real(1.0))
        if up2 < real(1.0e-14):
            sgn = real(1.0)
            if u[2] < real(0.0):
                sgn = real(-1.0)
            out = v3(px, sgn * py, sgn * ct)
        else:
            up = wp.sqrt(up2)
            out = v3(
                (u[0] * u[2] * px - u[1] * py) / up + u[0] * ct,
                (u[1] * u[2] * px + u[0] * py) / up + u[1] * ct,
                -up * px + u[2] * ct,
            )
        return wp.normalize(out)

    return SimpleNamespace(
        csda_mean_loss=csda_mean_loss,
        bohr_variance=bohr_variance,
        straggle_attempt=straggle_attempt,
        straggle_attempt_gamma=straggle_attempt_gamma,
        fdm=fdm,
        scattering_power_dm=scattering_power_dm,
        scattering_variance_birth=scattering_variance_birth,
        polar_deflection=polar_deflection,
        rotate_dir=rotate_dir,
        vec3=v3,
        real=name,
    )
