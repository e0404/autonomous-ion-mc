"""Exact relativistic two-body elastic kinematics (units MeV, c = 1; decision 0041 slice C).

A projectile of mass ``m1`` and kinetic energy ``t_mev`` hits a target of mass ``m2`` at rest and
scatters to the CM cosine ``mu_cm`` (azimuth 0, scattering plane x-z, beam along z). The projectile
is
boosted from the CM; the recoil is boosted from its own CM momentum, independently, so that the sum
of the two lab four-momenta can be compared with the incident one (closure, row P6 (7)). For equal
masses ``T_recoil = T (1 - mu_cm) / 2`` exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class TwoBody:
    """Lab four-momenta ``(E, p_x, p_z)`` [MeV] of both particles and the kinetic energies."""

    e1: NDArray[np.float64]
    px1: NDArray[np.float64]
    pz1: NDArray[np.float64]
    e2: NDArray[np.float64]
    px2: NDArray[np.float64]
    pz2: NDArray[np.float64]
    t1: NDArray[np.float64]
    t2: NDArray[np.float64]
    cos_lab1: NDArray[np.float64]

    def closure(self, t_mev: float, m1: float, m2: float) -> tuple[float, float]:
        """``(max |sum E - (T + m1 + m2)|, max |sum p - p_beam|)`` over the events (MeV)."""
        p0 = float(np.sqrt((t_mev + m1) ** 2 - m1**2))
        de = np.abs(self.e1 + self.e2 - (t_mev + m1 + m2))
        dp = np.maximum(np.abs(self.px1 + self.px2), np.abs(self.pz1 + self.pz2 - p0))
        return float(de.max()), float(dp.max())


def two_body(t_mev: float, m1: float, m2: float, mu_cm: Any) -> TwoBody:
    """Scatter to ``mu_cm`` (array-like in [-1, 1]); see the module docstring."""
    mu = np.asarray(mu_cm, dtype=np.float64)
    e1 = t_mev + m1
    p = np.sqrt(e1 * e1 - m1 * m1)
    rs = np.sqrt(m1 * m1 + m2 * m2 + 2.0 * e1 * m2)
    gam, bgam = (e1 + m2) / rs, p / rs  # gamma and beta*gamma of the CM
    pcm = p * m2 / rs
    sin = np.sqrt(np.clip(1.0 - mu * mu, 0.0, None))
    ecm1 = np.sqrt(pcm * pcm + m1 * m1)
    ecm2 = np.sqrt(pcm * pcm + m2 * m2)
    e1p = gam * ecm1 + bgam * pcm * mu
    pz1 = gam * pcm * mu + bgam * ecm1
    px1 = pcm * sin
    e2p = gam * ecm2 - bgam * pcm * mu
    pz2 = -gam * pcm * mu + bgam * ecm2
    px2 = -pcm * sin
    pl = np.sqrt(px1 * px1 + pz1 * pz1)
    return TwoBody(e1p, px1, pz1, e2p, px2, pz2, e1p - m1, e2p - m2, pz1 / pl)
