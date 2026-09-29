"""Energy-loss straggling: Bohr variance with the relativistic (Bethe–Livingston) factor.

Per unit mass path length the variance of the electronic energy loss is

    dΩ²/d(ρs) = 0.1569 MeV² cm²/g · z_eff² · (Z/A) · (1 − β²/2) / (1 − β²)

(4π N_A r_e² (m_e c²)² = 0.1569 MeV² cm²/g). The distribution of the loss in
a step is sampled as a Gamma distribution with the CSDA mean and this
variance (decision 0041): it is strictly positive, tends to the Gaussian
limit when κ = ξ/T_max is large, and avoids the bias that truncating a
Gaussian at zero would introduce for thin steps. Per-step spectra are not
Landau/Vavilov-shaped; only the sum over many steps (ranges, depth doses) is
represented, which is the documented domain.
"""

from __future__ import annotations

import numpy as np

BOHR_CONSTANT_MEV2_CM2_G = 0.1569


def bohr_variance_per_g_cm2(
    z_eff: np.ndarray | float, z_over_a: float, beta: np.ndarray | float
) -> np.ndarray:
    beta2 = np.asarray(beta, dtype=np.float64) ** 2
    return (
        BOHR_CONSTANT_MEV2_CM2_G
        * np.asarray(z_eff) ** 2
        * z_over_a
        * (1.0 - 0.5 * beta2)
        / (1.0 - beta2)
    )
