"""Tabulated Ashley-Ritchie-Brandt function for the Barkas (z^3) correction.

``BARKAS_W`` / ``BARKAS_F`` reproduce the 47-point table used by Geant4's
``G4EmCorrections`` (``fTable`` in ``G4EmCorrections.cc``), which tabulates
the function ``F_ARB(W)`` of J. C. Ashley, R. H. Ritchie and W. Brandt,
Phys. Rev. B 5 (1972) 2393, as adopted by ICRU Report 49. The correction is
evaluated as ``L1 = 1.29 * F(W) / (sqrt(Z X) X)`` with ``X = beta^2 / (alpha^2 Z)``
and ``W = b / sqrt(X)``; for ``W`` above the last tabulated point the value
is scaled by ``W_max / W`` as in Geant4. The element-dependent parameter
``b`` follows the same source (:func:`barkas_b_parameter`).

Provenance: Geant4 (Geant4 Software License) is the tabulation source; the
underlying function is published in the paper above. The numbers are data,
not code, and are reproduced here with attribution.
"""

from __future__ import annotations

BARKAS_W: tuple[float, ...] = (
    0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4,
    0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
    2.0, 2.1, 2.4, 3.0, 3.08, 3.1, 3.3, 3.5, 3.8, 4.0, 4.1, 4.8, 5.0, 5.1,
    6.0, 6.5, 7.0, 7.1, 8.0, 9.0, 10.0,
)  # fmt: skip

BARKAS_F: tuple[float, ...] = (
    21.5, 20.0, 18.0, 15.6, 15.0, 14.0, 13.5, 13.0, 12.2, 9.25, 7.0, 6.0,
    4.5, 3.5, 3.0, 2.5, 2.0, 1.7, 1.2, 1.0, 0.86, 0.7, 0.61, 0.52, 0.5,
    0.43, 0.42, 0.3, 0.2, 0.13, 0.1, 0.09, 0.08, 0.07, 0.06, 0.051, 0.04,
    0.03, 0.024, 0.02, 0.013, 0.01, 0.009, 0.008, 0.006, 0.0032, 0.0025,
)  # fmt: skip

assert len(BARKAS_W) == 47 and len(BARKAS_F) == 47

#: Overall normalisation applied to F_ARB (ICRU 49 / Geant4).
BARKAS_SCALE: float = 1.29


def barkas_b_parameter(atomic_number: int) -> float:
    """Parameter ``b`` of the ARB Barkas correction for a target element.

    Values follow Geant4's ``G4EmCorrections::BarkasCorrection`` (which cites
    ICRU 49): 1.8 for Z <= 10 (0.6 for liquid hydrogen and helium, not used
    here), 1.4 for 11 <= Z <= 17 and 19 <= Z <= 25, 1.8 for argon, 1.35 for
    26 <= Z <= 50, 1.3 otherwise.
    """
    z = atomic_number
    if z == 2:
        return 0.6
    if z <= 10:
        return 1.8
    if z <= 17:
        return 1.4
    if z == 18:
        return 1.8
    if z <= 25:
        return 1.4
    if z <= 50:
        return 1.35
    return 1.3
