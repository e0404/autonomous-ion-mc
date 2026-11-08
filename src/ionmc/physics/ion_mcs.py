# mypy: ignore-errors
# (Warp function sources use runtime precision types in annotations; see _wpfunc.py.)
"""Per-step multiple-scattering variance of light ions: the exact-Omega mapping of f_dM
(decision 0043 section 3) as shared Warp functions (python twin and Warp kernel).

For an ion (z >= 2) the projected scattering power of a step at the midpoint ``pv`` is::

    T_ion = f_dM(pv_eq, p1v1_eq) (z E_s / pv)^2 / X_S,   p1v1_eq / pv_eq = p1v1 / pv

where ``pv_eq`` is the root of ``b_p(pv_eq, p1v1_eq) = b_ion(pv, p1v1)``
(:func:`ionmc.physics.exact_omega.make_exact_omega`). The birth step (first step of a life)
keeps the analytic log-average of the memory term of ``scattering_variance_birth``
(:func:`ionmc.physics.em.make_em`): the memory ratio is the same for the ion and the equivalent
proton, so only the smooth ``lg pv`` term of f_dM and the prefactor change (``lg pv_eq`` at the
mapped midpoint state).

Thin-slab guard (P11 note: ``t = R(E1) - R(E)`` cancels for a state close to the initial
energy): the mapping is evaluated at ``pv_map = min(pv_mid, p1v1 / (1 + THIN_SLAB_R))`` with
``THIN_SLAB_R = 1e-6``; ``p1v1_eq = pv_eq p1v1 / pv_map`` keeps the memory ratio of the
evaluated state. The guard only acts on steps whose slab from the initial state is thinner than
1e-6 of the momentum, where the variance (proportional to the step length) is negligible. The
mapping runs in the precision passed to the factory (the kernel passes float64, as for every
energy quantity).

Two evaluations of ``pv_eq`` (``PhysicsOptions.mcs_mapping_eval``, decision 0043 section 3,
"Evaluation"), both the exact root on every step: :func:`mapped_state` / ``mcs_variance`` solve
it from the full proton-path bracket (``"root"``); ``mapped_state_table`` /
``mcs_variance_table`` (``"table"``, the default) solve the same equation to the same tolerance
from a tight bracket around the bilinear ``pv_eq/pv`` table value
(:func:`ionmc.physics.exact_omega.pv_equivalent_seeded`; the table is an accelerator only) and
return a bracket-miss flag (the full bracket was needed; a diagnostic, not a failure).
``interpolated_state`` is the bare table value (diagnostics and tests). All apply the same
thin-slab guard and the same variance formula (``variance_of_state``).

A nonzero status (see ``exact_omega.STATUS_NAMES``) returns variance 0; the caller counts it
(``mcs_mapping_failure``) and the result is invalid (fail closed).
"""

import functools
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func
from ionmc.physics.em import E_S_MEV, FDM_COEFFICIENTS, make_em
from ionmc.physics.exact_omega import make_exact_omega

wp.set_module_options({"enable_backward": False})

THIN_SLAB_R = 1.0e-6
"""Smallest relative momentum excess ``p1v1/pv - 1`` at which the mapping is evaluated."""


@functools.cache
def make_ion_mcs(real: type) -> SimpleNamespace:
    """Return the ion MCS functions for precision ``real``."""
    name = check_real(real)
    em = make_em(real)
    xo = make_exact_omega(real)
    es_const = wp.constant(real(E_S_MEV))
    c0 = wp.constant(real(FDM_COEFFICIENTS[0]))
    c1 = wp.constant(real(FDM_COEFFICIENTS[1]))
    c2 = wp.constant(real(FDM_COEFFICIENTS[2]))
    c3 = wp.constant(real(FDM_COEFFICIENTS[3]))
    thin = wp.constant(real(THIN_SLAB_R))

    @named_func(name)
    def mapped_state(
        ion_pk: wp.array2d(dtype=real),
        pro_pk: wp.array2d(dtype=real),
        comp_pk: wp.array2d(dtype=real),
        quad: wp.array(dtype=real),
        m: int,
        z: real,
        mass_mev: real,
        pv_mid: real,
        p1v1: real,
    ) -> tuple[real, real, int]:
        """``(pv_eq, p1v1_eq, status)`` of the ion state (``pv_mid``, ``p1v1``) in material ``m``
        (rows of the packed arrays); thin-slab guard in the module docstring."""
        pv_map = wp.min(pv_mid, p1v1 / (real(1.0) + thin))
        peq, status, iters = xo.pv_equivalent(
            ion_pk[m], pro_pk[m], comp_pk[m], quad, z, mass_mev, pv_map, p1v1, 0
        )
        return peq, peq * (p1v1 / pv_map), status

    @named_func(name)
    def interpolated_state(
        tab: wp.array2d(dtype=real), m: int, pv_mid: real, p1v1: real
    ) -> tuple[real, real, int]:
        """``(pv_eq, p1v1_eq, status)`` interpolated from the ``pv_eq/pv`` table row ``m`` (status
        5 outside the table); the same thin-slab guard as :func:`mapped_state`. The starting
        guess of :func:`mapped_state_table`, not a transport quantity."""
        pv_map = wp.min(pv_mid, p1v1 / (real(1.0) + thin))
        peq, status = xo.pv_equivalent_table(tab[m], pv_map, p1v1)
        return peq, peq * (p1v1 / pv_map), status

    @named_func(name)
    def mapped_state_table(
        ion_pk: wp.array2d(dtype=real),
        pro_pk: wp.array2d(dtype=real),
        comp_pk: wp.array2d(dtype=real),
        quad: wp.array(dtype=real),
        tab: wp.array2d(dtype=real),
        m: int,
        z: real,
        mass_mev: real,
        pv_mid: real,
        p1v1: real,
    ) -> tuple[real, real, int, int]:
        """``(pv_eq, p1v1_eq, status, bracket_miss)``: the exact root of :func:`mapped_state`
        (same equation, tolerance and status codes) started from a tight bracket around the
        table value of row ``m`` (half-width ``tab[m, 8]``); the same thin-slab guard."""
        pv_map = wp.min(pv_mid, p1v1 / (real(1.0) + thin))
        row = tab[m]
        q0 = real(0.0)
        g, st_t = xo.pv_equivalent_table(row, pv_map, p1v1)
        if st_t == 0:
            q0 = g / pv_map
        peq, status, iters, miss = xo.pv_equivalent_seeded(
            ion_pk[m], pro_pk[m], comp_pk[m], quad, z, mass_mev, pv_map, p1v1, q0, row[8]
        )
        return peq, peq * (p1v1 / pv_map), status, miss

    @named_func(name)
    def variance_of_state(
        peq: real,
        p1eq: real,
        status: int,
        z: real,
        pv_mid: real,
        pv_end: real,
        p1v1: real,
        birth: int,
        inv_rho_xs_cm2_g: real,
        rho_g_cm3: real,
        s_mm: real,
    ) -> real:
        """Projected variance [rad^2] of the step from the mapped state (0 for a nonzero status);
        ``birth`` = 1 selects the analytic memory-term average of the first step."""
        var = real(0.0)
        if status == 0:
            q = z * es_const / pv_mid
            pref = q * q * rho_g_cm3 * inv_rho_xs_cm2_g / real(10.0)
            if birth == 1:
                om = real(1.0) - (pv_end / p1v1) * (pv_end / p1v1)
                if om > real(0.0):
                    l1 = wp.log10(om) - real(0.4342944819032518)
                    l2 = wp.log10(peq)
                    f = wp.max(c0 + c1 * l1 + c2 * l2 + c3 * l1 * l2, real(0.0))
                    var = f * pref * s_mm
            else:
                f = em.fdm(peq, p1eq)
                var = f * pref * s_mm
        return var

    @named_func(name)
    def mcs_variance(
        ion_pk: wp.array2d(dtype=real),
        pro_pk: wp.array2d(dtype=real),
        comp_pk: wp.array2d(dtype=real),
        quad: wp.array(dtype=real),
        m: int,
        z: real,
        mass_mev: real,
        pv_mid: real,
        pv_end: real,
        p1v1: real,
        birth: int,
        inv_rho_xs_cm2_g: real,
        rho_g_cm3: real,
        s_mm: real,
    ) -> tuple[real, int]:
        """``(projected variance [rad^2] of the step, mapping status)`` for an ion, root solve;
        ``birth`` = 1 selects the analytic memory-term average of the first step (``pv_end`` its
        end)."""
        peq, p1eq, status = mapped_state(
            ion_pk, pro_pk, comp_pk, quad, m, z, mass_mev, pv_mid, p1v1
        )
        var = variance_of_state(
            peq, p1eq, status, z, pv_mid, pv_end, p1v1, birth, inv_rho_xs_cm2_g, rho_g_cm3, s_mm
        )
        return var, status

    @named_func(name)
    def mcs_variance_table(
        ion_pk: wp.array2d(dtype=real),
        pro_pk: wp.array2d(dtype=real),
        comp_pk: wp.array2d(dtype=real),
        quad: wp.array(dtype=real),
        tab: wp.array2d(dtype=real),
        m: int,
        z: real,
        mass_mev: real,
        pv_mid: real,
        pv_end: real,
        p1v1: real,
        birth: int,
        inv_rho_xs_cm2_g: real,
        rho_g_cm3: real,
        s_mm: real,
    ) -> tuple[real, int, int]:
        """``(variance, status, bracket_miss)``: as :func:`mcs_variance` with the table-accelerated
        exact root (:func:`mapped_state_table`, ``mcs_mapping_eval = "table"``)."""
        peq, p1eq, status, miss = mapped_state_table(
            ion_pk, pro_pk, comp_pk, quad, tab, m, z, mass_mev, pv_mid, p1v1
        )
        var = variance_of_state(
            peq, p1eq, status, z, pv_mid, pv_end, p1v1, birth, inv_rho_xs_cm2_g, rho_g_cm3, s_mm
        )
        return var, status, miss

    return SimpleNamespace(
        mapped_state=mapped_state,
        interpolated_state=interpolated_state,
        mapped_state_table=mapped_state_table,
        variance_of_state=variance_of_state,
        mcs_variance=mcs_variance,
        mcs_variance_table=mcs_variance_table,
        real=name,
    )
