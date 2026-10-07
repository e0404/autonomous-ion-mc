"""Tripathi light-system p + X cross section (NASA/TP-1999-209726; build-time numpy)."""

import math

import numpy as np
import pytest

from ionmc.physics.tripathi import extension_factor, sigma_tripathi_light_p

ENERGIES = (50.0, 100.0, 150.0, 200.0)


def _hand(e: float, z_t: int, a_t: int, rt: float) -> float:
    """Independent scalar evaluation of eqs. (3)-(15) with ``math`` only (hand evaluation: the
    terms are computed in the order of the paper: S, D, C_E, dE, X_m, R, B, sigma)."""
    mt = a_t * 931.49410242 - z_t * 0.51099895
    mp = 938.27208816
    ecm = math.sqrt(mp**2 + mt**2 + 2 * (e + mp) * mt) - mt - mp
    at = a_t ** (1 / 3)
    s = at / (1 + at)
    d = 1.85 + 0.16 / (1 + math.exp((500 - e) / 200))
    ce = d * (1 - math.exp(-e / 23)) - 0.292 * math.exp(-e / 792) * math.cos(0.229 * e**0.453)
    de = 1.85 * s + 0.16 * s / ecm ** (1 / 3) - ce + 0.91 * (a_t - 2 * z_t) / a_t
    x1 = 2.83 - 0.031 * a_t + 1.7e-4 * a_t**2
    xm = 1 - x1 * math.exp(-e / (x1 * (1.2 + 1.6 * (1 - math.exp(-e / 15)))))
    rad = 1.29 * (0.875 + rt) + 1.2 * (1 + at) / ecm ** (1 / 3)
    b = 1.44 * z_t / rad
    return math.pi * (1.1 * (1 + at + de)) ** 2 * (1 - b / ecm) * xm * 0.01


@pytest.mark.parametrize(("z", "a", "rt"), [(6, 12, 2.472), (8, 16, 2.730)])
def test_matches_hand_evaluation(z: int, a: int, rt: float) -> None:
    for e in ENERGIES:
        assert float(sigma_tripathi_light_p(e, z, a)) == pytest.approx(_hand(e, z, a, rt), rel=1e-6)


def test_sanity_and_extension() -> None:
    e = np.linspace(20.0, 250.0, 231)
    for z, a in ((6, 12), (8, 16), (13, 27)):
        s = sigma_tripathi_light_p(e, z, a)
        assert np.all(s > 0) and np.all(np.abs(np.diff(s, 2)) < 0.02)
    assert float(sigma_tripathi_light_p(200.0, 6, 12)) == pytest.approx(0.220, rel=0.2)
    assert extension_factor(150.0, 6, 12) == 1.0
    assert float(extension_factor(200.0, 8, 16)) == pytest.approx(
        float(sigma_tripathi_light_p(200.0, 8, 16) / sigma_tripathi_light_p(150.0, 8, 16))
    )
    # radius sensitivity: +5 % radii change sigma(100 MeV) by < 0.5 %
    s0 = float(sigma_tripathi_light_p(100.0, 6, 12))
    assert abs(float(sigma_tripathi_light_p(100.0, 6, 12, 1.05)) / s0 - 1.0) < 0.005
    with pytest.raises(ValueError):
        sigma_tripathi_light_p(100.0, 11, 23)
    with pytest.raises(ValueError):
        sigma_tripathi_light_p(0.0, 6, 12)


def test_print_values() -> None:
    for z, a in ((6, 12), (8, 16)):
        print(
            "TRIPATHI",
            a,
            [round(float(sigma_tripathi_light_p(e, z, a)) * 1000, 3) for e in ENERGIES],
        )
