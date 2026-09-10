"""Shared-source proton nonelastic nuclear physics (decision 0012).

Task DEV-007 adds catastrophic proton nonelastic nuclear reactions: along the
track a primary is removed at the macroscopic nonelastic rate, a local fraction
of its energy is deposited (short-range recoils and fragments) and the rest
(secondary protons, neutrons, gammas) is booked to an escaping/deferred channel
for later secondary transport. The nonelastic cross section on oxygen-16 is an
analytic parameterization of the ICRU-63 shape (threshold 7 MeV, ~550 mb peak
at 20 MeV, ~340-400 mb plateau over 100-250 MeV); it reproduces the published
primary survival to the Bragg peak (~0.80 at 150 MeV, ~0.73 at 200 MeV). The
hydrogen channel is proton-proton elastic scattering, which deflects rather
than removes the primary, so it is not counted here (it belongs with secondary
production). Written against the math namespace ``m`` so it runs as the Warp
kernel and the float64 reference identically.

Units: energy in MeV, cross section in barn, number density in cm^-3,
macroscopic cross section in 1/cm.
"""

from __future__ import annotations

from ionmc.backend import mathlib

m = mathlib.current()
func = mathlib.func

#: barn to cm^2.
BARN_TO_CM2: float = 1.0e-24
#: mm per cm.
MM_PER_CM: float = 10.0
#: Nonelastic threshold [MeV]; below it protons do not react.
NONELASTIC_THRESHOLD_MEV: float = 7.0


@func
def nonelastic_cross_section_oxygen(energy: float) -> float:
    """Proton nonelastic cross section on oxygen-16 [barn] (ICRU-63 shape).

    ``turnon(E) [plateau(E) + peak(E)]``: a linear threshold turn-on from 7 to
    20 MeV, a plateau declining mildly from 0.40 b at 100 MeV to ~0.34 b at
    250 MeV, and a Gaussian resonance bump of 0.15 b centred at 20 MeV.
    """
    turn = m.min(m.max((energy - NONELASTIC_THRESHOLD_MEV) / 13.0, 0.0), 1.0)
    plateau = m.max(0.40 - 0.0004 * m.max(energy - 100.0, 0.0), 0.30)
    reduced = (energy - 20.0) / 15.0
    bump = 0.15 * m.exp(-reduced * reduced)
    return turn * (plateau + bump)


@func
def macroscopic_nonelastic(energy: float, oxygen_density_per_cm3: float) -> float:
    """Macroscopic nonelastic cross section ``Sigma`` [1/cm] for the medium.

    ``Sigma = n_O sigma_nonel(E)``; only oxygen contributes to catastrophic
    primary removal (the hydrogen channel is elastic).
    """
    return (
        oxygen_density_per_cm3 * nonelastic_cross_section_oxygen(energy) * BARN_TO_CM2
    )


@func
def nonelastic_step_probability(
    energy: float, step_mm: float, oxygen_density_per_cm3: float
) -> float:
    """Probability of a nonelastic reaction over a step, ``Sigma(E) * dl``.

    The leading-order (thin-step) form; for the ~1 mm steps used here it is
    accurate to better than 0.1 % versus ``1 - exp(-Sigma dl)`` (decision 0012).
    ``energy`` is the step-entry energy. Zero below the nonelastic threshold.
    Clamped to at most 1: with ``Sigma ~ 0.012/cm`` and 1 mm steps ``P ~ 1e-3``,
    but the clamp keeps the returned value a valid probability if a caller ever
    raises the step cap far enough that the thin-step form would exceed 1.
    """
    sigma = macroscopic_nonelastic(energy, oxygen_density_per_cm3)
    dl_cm = step_mm / MM_PER_CM
    p = m.min(sigma * dl_cm, 1.0)
    return m.where(energy > NONELASTIC_THRESHOLD_MEV, p, 0.0)
