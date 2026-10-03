"""Scattering length ``X_S`` for the differential Moliere scattering power.

Gottschalk, Med. Phys. 37 (2010) 352 (arXiv:0908.1413), eq. (XS): for an element of atomic
number ``Z`` and molar mass ``A`` [g/mol]::

    1 / X_S = alpha N_A r_e^2 (Z^2 / A) { 2 ln( 33219 (A Z)^(-1/3) ) - 1 }    [cm^2/g]

with the natural logarithm (the test ``test_scattering_length_matches_gottschalk`` confirms
this against the tabulated values). Compounds follow Bragg additivity by weight,
``1/X_S = sum_i w_i / X_S,i``. ``X_S`` is a mass thickness [g/cm^2]; the scattering length in
a medium of density rho is ``X_S / rho`` [cm].
"""

from __future__ import annotations

import math

from ionmc.materials import N_A, Material

ALPHA_FINE_STRUCTURE = 1.0 / 137.035999084
CLASSICAL_ELECTRON_RADIUS_CM = 2.8179403262e-13


def inverse_scattering_length_cm2_per_g(material: Material) -> float:
    """``1/X_S`` of ``material`` [cm^2/g] (Bragg-additive Gottschalk formula)."""
    prefactor = ALPHA_FINE_STRUCTURE * N_A * CLASSICAL_ELECTRON_RADIUS_CM**2
    total = 0.0
    for element, weight in material._fractions:
        a = element.A_g_mol
        z = float(element.Z)
        bracket = 2.0 * math.log(33219.0 * (a * z) ** (-1.0 / 3.0)) - 1.0
        total += weight * prefactor * z * z / a * bracket
    return total


def scattering_length_g_cm2(material: Material) -> float:
    """Scattering length ``X_S`` of ``material`` as a mass thickness [g/cm^2]."""
    return 1.0 / inverse_scattering_length_cm2_per_g(material)
