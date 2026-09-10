"""Fermi-Eyges lateral-spread oracle for multiple-scattering validation.

The projected (one transverse axis) spatial variance of a zero-size,
zero-divergence proton pencil beam at depth ``z`` in a homogeneous medium is
the Fermi-Eyges A2 moment (decision 0011):

    sigma_x^2(z) = integral_0^z (z - u)^2 T(u) du

with the projected scattering power ``T(u) = (13.6 / pv(u))^2 / X0`` [rad^2/cm]
(the Highland scattering-power form, no log term), ``pv`` the momentum-velocity
product and ``X0`` the radiation length. The proton energy ``E(u)`` at depth
``u`` is obtained by inverting the CSDA range table. This is a **numpy
reference oracle** (not transport shared source): it is the analytic quantity
the Monte Carlo lateral spread is validated against.

Units: depths in mm, energies in MeV, radiation length in g/cm^2, density in
g/cm^3, result in mm.
"""

from __future__ import annotations

import numpy as np

from ionmc.constants import PROTON_MASS_MEV
from ionmc.physics.transport import HIGHLAND_CONSTANT_MEV
from ionmc.tabulated_stopping_power import TabulatedStoppingPower

MM_PER_CM = 10.0


def energy_at_depth(
    model: TabulatedStoppingPower,
    energy0_mev: float,
    depth_mm: np.ndarray,
    density_g_per_cm3: float = 1.0,
    energy_floor_mev: float = 1.0,
) -> np.ndarray:
    """Proton kinetic energy [MeV] at each depth, from the CSDA range table.

    ``E(u)`` such that the residual mass range ``R(E0) - rho*u`` (g/cm^2) equals
    the CSDA range of ``E(u)``.
    """
    e_grid = np.linspace(energy_floor_mev, energy0_mev, 6000)
    r_grid = model.csda_range(e_grid)  # g/cm^2, increasing in E
    r0 = float(model.csda_range(energy0_mev)[0])
    residual = r0 - density_g_per_cm3 * (depth_mm / MM_PER_CM)
    residual = np.clip(residual, r_grid[0], r_grid[-1])
    return np.interp(residual, r_grid, e_grid)


def lateral_sigma_x_mm(
    model: TabulatedStoppingPower,
    energy0_mev: float,
    depths_mm: np.ndarray,
    radiation_length_g_per_cm2: float,
    density_g_per_cm3: float = 1.0,
    n_integration: int = 4000,
) -> np.ndarray:
    """Projected Fermi-Eyges lateral sigma_x [mm] at each requested depth."""
    out = np.empty_like(np.asarray(depths_mm, dtype=np.float64))
    for i, z_mm in enumerate(np.atleast_1d(depths_mm)):
        u_cm = np.linspace(1.0e-4, z_mm / MM_PER_CM, n_integration)
        e_u = energy_at_depth(model, energy0_mev, u_cm * MM_PER_CM, density_g_per_cm3)
        pv = e_u * (e_u + 2.0 * PROTON_MASS_MEV) / (e_u + PROTON_MASS_MEV)
        # projected scattering power T = (13.6/pv)^2 * rho / X0 [rad^2/cm]
        scattering_power = (
            (HIGHLAND_CONSTANT_MEV / pv) ** 2
            * density_g_per_cm3
            / radiation_length_g_per_cm2
        )
        integrand = (z_mm / MM_PER_CM - u_cm) ** 2 * scattering_power  # cm^2/cm
        out[i] = np.sqrt(np.trapezoid(integrand, u_cm)) * MM_PER_CM
    return out


def _density_at_depth(
    depth_cm: np.ndarray, z_boundaries_mm: np.ndarray, densities: np.ndarray
) -> np.ndarray:
    """Piecewise-constant density [g/cm^3] at each depth [cm] for a VoxelSlab
    profile (decision 0016)."""
    z_cm = np.asarray(z_boundaries_mm, dtype=np.float64) / MM_PER_CM
    idx = np.clip(
        np.searchsorted(z_cm, depth_cm, side="right") - 1, 0, densities.size - 1
    )
    return np.asarray(densities, dtype=np.float64)[idx]


def lateral_sigma_x_heterogeneous_mm(
    model: TabulatedStoppingPower,
    energy0_mev: float,
    depths_mm: np.ndarray,
    radiation_length_g_per_cm2: float,
    z_boundaries_mm: np.ndarray,
    densities: np.ndarray,
    n_integration: int = 6000,
) -> np.ndarray:
    """Projected Fermi-Eyges lateral sigma_x [mm] for a piecewise-density profile.

    The energy vs depth follows the integrated water-equivalent thickness
    ``integral rho dl`` and the scattering power uses the *local* density
    (decision 0016), so a density interface produces the correct sigma_x kink.
    """
    z_max = float(np.max(depths_mm))
    u_cm = np.linspace(0.0, z_max / MM_PER_CM, n_integration)
    rho_u = _density_at_depth(u_cm, z_boundaries_mm, densities)
    # integrated water-equivalent thickness [g/cm^2] from 0 to each u
    step_wet = 0.5 * (rho_u[1:] + rho_u[:-1]) * np.diff(u_cm)
    wet = np.concatenate([[0.0], np.cumsum(step_wet)])
    e_grid = np.linspace(1.0, energy0_mev, 6000)
    r_grid = model.csda_range(e_grid)  # g/cm^2
    r0 = float(model.csda_range(energy0_mev)[0])
    residual = np.clip(r0 - wet, r_grid[0], r_grid[-1])
    e_u = np.interp(residual, r_grid, e_grid)
    pv = e_u * (e_u + 2.0 * PROTON_MASS_MEV) / (e_u + PROTON_MASS_MEV)
    scattering_power = (
        (HIGHLAND_CONSTANT_MEV / pv) ** 2 * rho_u / radiation_length_g_per_cm2
    )
    out = np.empty_like(np.asarray(depths_mm, dtype=np.float64))
    for i, z_mm in enumerate(np.atleast_1d(depths_mm)):
        z_c = z_mm / MM_PER_CM
        mask = u_cm <= z_c
        integrand = (z_c - u_cm[mask]) ** 2 * scattering_power[mask]
        out[i] = np.sqrt(np.trapezoid(integrand, u_cm[mask])) * MM_PER_CM
    return out
