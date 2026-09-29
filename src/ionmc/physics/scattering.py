"""Multiple Coulomb scattering: Gottschalk's differential Molière scattering power.

Following Gottschalk (arXiv:0908.1413, 2010) the mean squared projected
scattering angle grows per unit path length as

    T_dM = f_dM(pv, p₁v₁) · (E_s / pv)² / X_S,      E_s = 15.0 MeV,

    f_dM = 0.5244 + 0.1975 lg(1 − (pv/p₁v₁)²) + 0.2320 lg(pv)
           − 0.0098 lg(pv) lg(1 − (pv/p₁v₁)²),

with pv the product of momentum and velocity (MeV) at the point of interest,
p₁v₁ at the particle's birth, lg the decimal logarithm and X_S the material's
scattering length (cm),

    1/(ρ X_S) = α N_A r_e² Σ_i w_i (Z_i²/A_i) {2 ln(33219 (A_i Z_i)^{-1/3}) − 1}.

The fit reproduces Molière/Fano/Hanson theory for protons of 3–300 MeV in
slabs with 0.001 ≤ 1 − (pv/p₁v₁)² ≤ 0.97 to a few per cent and is non-local
(depends on the initial energy), which removes the step-size dependence of
per-step Highland sampling. For ions the scattering power scales with z²
(Rutherford); this Gaussian core carries no single-scattering tail.
"""

from __future__ import annotations

import math

import numpy as np

from ionmc.materials import ELEMENTS, Material
from ionmc.species import Species

E_S_MEV = 15.0
ALPHA = 1.0 / 137.035999084
N_A = 6.02214076e23
R_E_CM = 2.8179403262e-13


def scattering_length_g_cm2(material: Material) -> float:
    """ρ X_S in g/cm² for a compound (Bragg rule over elements)."""
    total = 0.0
    for symbol, w in material.composition.items():
        z, a, _ = ELEMENTS[symbol]
        total += (
            w * (z * z / a) * (2.0 * math.log(33219.0 * (a * z) ** (-1.0 / 3.0)) - 1.0)
        )
    return 1.0 / (ALPHA * N_A * R_E_CM**2 * total)


def pv_mev(species: Species, t_mev_per_u: np.ndarray | float) -> np.ndarray:
    """Momentum times velocity, pv = β² γ m c² = T (T + 2 m c²)/(T + m c²), in MeV."""
    t = np.asarray(t_mev_per_u, dtype=np.float64) * species.a
    m = species.mass_mev
    return t * (t + 2.0 * m) / (t + m)


def f_dm(pv: np.ndarray | float, p1v1: float) -> np.ndarray:
    ratio = np.clip(1.0 - (np.asarray(pv, dtype=np.float64) / p1v1) ** 2, 1e-3, 0.97)
    lg_r = np.log10(ratio)
    lg_pv = np.log10(np.asarray(pv, dtype=np.float64))
    return 0.5244 + 0.1975 * lg_r + 0.2320 * lg_pv - 0.0098 * lg_pv * lg_r


def scattering_power_per_cm(
    species: Species,
    pv: np.ndarray | float,
    p1v1: float,
    rho_x_s_g_cm2: float,
    density_g_cm3: float,
) -> np.ndarray:
    """T_dM in rad²/cm for a projectile of charge z: z² f_dM (E_s/pv)² ρ / (ρX_S)."""
    x_s_cm = rho_x_s_g_cm2 / density_g_cm3
    return (
        species.z**2
        * f_dm(pv, p1v1)
        * (E_S_MEV / np.asarray(pv, dtype=np.float64)) ** 2
        / x_s_cm
    )


def fermi_eyges_lateral_sigma_mm(
    species: Species,
    energies_mev_per_u: np.ndarray,
    path_cm: np.ndarray,
    rho_x_s_g_cm2: float,
    density_g_cm3: float,
) -> float:
    """Lateral (one-axis) sigma at the end of a path from Fermi–Eyges A₂ = ∫ (x−x')² T dx'.

    ``energies_mev_per_u`` gives the kinetic energy along the path at the
    positions ``path_cm`` (both increasing in depth); used for validation.
    """
    pv = pv_mev(species, energies_mev_per_u)
    t = scattering_power_per_cm(species, pv, float(pv[0]), rho_x_s_g_cm2, density_g_cm3)
    x_end = path_cm[-1]
    integrand = (x_end - path_cm) ** 2 * t
    a2 = np.trapezoid(integrand, path_cm)
    return float(np.sqrt(a2)) * 10.0
