"""P1 (V3-008): the numpy Moliere-Hanson reference ``validation/scripts/ion_em/moliere_hanson.py``.

Pure-function identities run in CI. The proton validation against theta_Hanson of Gottschalk,
Med. Phys. 37 (2010) 352 (arXiv:0908.1413, ``tbl:theta0Single``, Be/Al/Cu/Pb at 158.6 MeV) needs
the hash-pinned arXiv source; it is looked up in the repository cache and in
``$IONMC_CACHE_DIR/objects``. When absent the test is skipped, and FAILS under
``IONMC_REQUIRE_DATA=1`` (repository cached-data convention).
"""

import importlib.util
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from ionmc.materials import ALUMINIUM, BERYLLIUM, COPPER, LEAD, WATER, Material
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, PROTON
from ionmc.physics.scattering import inverse_scattering_length_cm2_per_g
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.transport.mcs_checks import (
    GOTTSCHALK_CACHE,
    GOTTSCHALK_SHA256,
    gottschalk_tex,
    parse_tables,
)

_PATH = Path(__file__).resolve().parents[2] / "validation/scripts/ion_em/moliere_hanson.py"
_SPEC = importlib.util.spec_from_file_location("moliere_hanson", _PATH)
assert _SPEC is not None and _SPEC.loader is not None
MH = importlib.util.module_from_spec(_SPEC)
sys.modules["moliere_hanson"] = MH
_SPEC.loader.exec_module(MH)

PROTON_TOLERANCE = 0.01
"""Plan P1: <= 1.0 % at every tabulated point (NOT loosened; see the P1 note for the outcome)."""


def _elements(m: Material) -> tuple[tuple[int, float, float], ...]:
    return tuple((e.Z, e.A_g_mol, w) for e, w in m._fractions)


def _path(bethe: BetheStoppingSource, m: Material, e0: float):  # type: ignore[no-untyped-def]
    table = bethe.table(m, PROTON)
    r1 = float(table.range_at(e0))

    def energy_at(t: np.ndarray) -> np.ndarray:
        rem = np.maximum(r1 - t, float(table.csda_range_g_cm2[0]))
        return np.asarray(table.energy_from_range(rem))

    return r1, energy_at


def _constant_energy(e: float):  # type: ignore[no-untyped-def]
    return lambda t: np.full_like(t, e)


# --------------------------------------------------------------------------- identities (CI)
def test_constants() -> None:
    assert MH.K_CHI_C == pytest.approx(0.156915, rel=2e-5)
    assert MH.OMEGA_FACTOR == pytest.approx(1.167, abs=5e-4)
    # chi_0^2 (1.13) Z^(2/3)/p^2 reproduces Lynch & Dahl's 2.007e-5 Z^(2/3)/p^2 (MeV/c)
    assert MH.CHI_0_MEV**2 * 1.13 == pytest.approx(2.007e-5, rel=2e-3)
    assert MH.FDM == pytest.approx((0.5244, 0.1975, 0.2320, -0.0098))


def test_water_composition_and_xs_match_ionmc() -> None:
    for ours, theirs in zip(_elements(WATER), MH.WATER, strict=True):
        assert ours[0] == theirs[0]
        assert ours[1:] == pytest.approx(theirs[1:], rel=1e-12)
    for m in (WATER, BERYLLIUM, COPPER, LEAD):
        assert MH.inverse_scattering_length(_elements(m)) == pytest.approx(
            inverse_scattering_length_cm2_per_g(m), rel=1e-12
        )


def test_solve_b() -> None:
    for b in (2.0, 5.0, 10.0, 20.0):
        big = MH.solve_b(b)
        assert big - math.log(big) == pytest.approx(b, rel=1e-13)
    with pytest.raises(ValueError):
        MH.solve_b(0.9)


def test_z_scaling_identity_at_fixed_velocity() -> None:
    """At fixed beta and t: chi_c^2 ~ z^2/(pv)^2 exactly; ln chi_a^2 changes only by the mass
    (through p) and by ln F_scr (through alpha')."""
    w = MH.WATER
    gamma = 1.3
    for proj in (ALPHA, CARBON12):
        e_p = PROTON.mass_mev * (gamma - 1.0)
        e_i = proj.mass_mev * (gamma - 1.0)
        sp = MH.moliere_slab(1.0, PROTON.mass_mev, w, 0.5, _constant_energy(e_p))
        si = MH.moliere_slab(float(proj.z), proj.mass_mev, w, 0.5, _constant_energy(e_i))
        mass_ratio = proj.mass_mev / PROTON.mass_mev
        assert si["chi_c2"] / sp["chi_c2"] == pytest.approx(proj.z**2 / mass_ratio**2, rel=1e-12)
        b2 = float(MH.beta2(e_p, PROTON.mass_mev))
        expected = -2.0 * math.log(mass_ratio) + math.log(
            float(MH.f_scr(proj.z, np.array([b2]), w)[0])
        )
        assert si["ln_chi_a2"] - sp["ln_chi_a2"] == pytest.approx(expected, abs=1e-12)
        # Omega_ion / Omega_p = (M/m_p) / F_scr at equal velocity and x/R (t_ion = (M/m_p)/z^2 t_p)
        si2 = MH.moliere_slab(
            float(proj.z), proj.mass_mev, w, 0.5 * mass_ratio / proj.z**2, _constant_energy(e_i)
        )
        f = float(MH.f_scr(proj.z, np.array([b2]), w)[0])
        assert si2["b"] - sp["b"] == pytest.approx(math.log(mass_ratio / f), abs=1e-12)


def test_f_scr_proton_is_one_and_grows_with_z() -> None:
    b2 = np.array([0.3, 0.5])
    assert MH.f_scr(1.0, b2, MH.WATER) == pytest.approx(np.ones(2), rel=1e-15)
    values = [MH.f_scr(z, b2, MH.WATER) for z in (2.0, 6.0, 8.0)]
    assert np.all(np.diff(np.array(values), axis=0) > 0.0)


def test_monotone_in_thickness(bethe: BetheStoppingSource) -> None:
    r1, energy_at = _path(bethe, WATER, 150.0)
    thetas = [
        MH.theta0_hanson(1.0, PROTON.mass_mev, MH.WATER, f * r1, energy_at)
        for f in (0.01, 0.1, 0.3, 0.6, 0.9)
    ]
    assert np.all(np.diff(thetas) > 0.0)


def test_compound_rule_reduces_to_element() -> None:
    """Water vs its components: a 'compound' of one element split in two equal parts equals the
    element; and water's chi_c^2 is the mass-weighted sum of the H and O contributions."""
    o = ((8, 15.9994, 1.0),)
    o_split = ((8, 15.9994, 0.5), (8, 15.9994, 0.5))
    e = _constant_energy(150.0)
    a = MH.moliere_slab(1.0, PROTON.mass_mev, o, 2.0, e)
    b = MH.moliere_slab(1.0, PROTON.mass_mev, o_split, 2.0, e)
    assert a["theta0"] == pytest.approx(b["theta0"], rel=1e-13)
    w = MH.moliere_slab(1.0, PROTON.mass_mev, MH.WATER, 2.0, e)
    h = MH.moliere_slab(1.0, PROTON.mass_mev, ((1, 1.00794, 1.0),), 2.0, e)
    wh, wo = MH.WATER[0][2], MH.WATER[1][2]
    assert w["chi_c2"] == pytest.approx(wh * h["chi_c2"] + wo * a["chi_c2"], rel=1e-12)


def test_fano_term_and_electron_switch() -> None:
    """P1b: for an element Fano's delta_b is D/Z with D = ln[1130 Z^-4/3 beta^2/(1-beta^2)]
    - u_in - beta^2/2 (Fano 1954 (b)); -u_in = 3.6 (H), 5.0 otherwise; the "bethe" variant has
    no delta_b and Z(Z+1); unknown switches are rejected."""
    b2 = np.array([0.2, 0.4])
    for zt, neg_u in ((1, 3.6), (4, 5.0), (82, 5.0)):
        d = np.log(1130.0 * zt ** (-4.0 / 3.0) * b2 / (1.0 - b2)) + neg_u - 0.5 * b2
        el = ((zt, 2.0 * zt, 1.0),)
        assert MH.fano_delta_b(el, b2) == pytest.approx(d / zt, rel=1e-13)
    e = _constant_energy(150.0)
    o = ((8, 15.9994, 1.0),)
    fano = MH.moliere_slab(1.0, PROTON.mass_mev, o, 2.0, e)
    bethe = MH.moliere_slab(1.0, PROTON.mass_mev, o, 2.0, e, electrons="bethe")
    assert bethe["delta_b"] == 0.0
    assert fano["chi_c2"] / bethe["chi_c2"] == pytest.approx(8.0 / 9.0, rel=1e-13)
    with pytest.raises(ValueError):
        MH.moliere_slab(1.0, PROTON.mass_mev, o, 2.0, e, electrons="none")


def test_candidate_a_is_velocity_scaling_and_protons_agree() -> None:
    """For z = 1, M = m_p the candidates a, b, d and e coincide (correction ln(1/1) = 0)."""
    e = _constant_energy(150.0)
    vals = {
        c: MH.theta0_candidate(
            c, 1.0, PROTON.mass_mev, MH.WATER, 2.0, 150.0, e, lambda x: x, lambda x: x
        )
        for c in ("a", "b", "d", "e")
    }
    assert max(vals.values()) == pytest.approx(min(vals.values()), rel=1e-12)


@pytest.mark.xfail(
    strict=True,
    reason="P1 finding (validation/plans/v3-008/p1-moliere-hanson-reference.md): the decision's "
    "estimate n_Omega ~ 0.77 / c_Omega ~ 0.13 was non-relativistic; the complete-Omega formula "
    "(screening and Fano terms included) gives n_Omega = 1.15-1.38 and c_Omega ~ 0.07-0.09. "
    "Band kept as declared, not loosened.",
)
def test_c_omega_close_to_planning_estimate(bethe: BetheStoppingSource) -> None:
    """c_Omega at the proton level (water, 50-300 MeV, L_1 = lg 0.5) lies near the decision's
    estimate 0.13; frozen formula, not a fitted number."""
    t = bethe.table(WATER, PROTON)
    e = np.array([50.0, 150.0, 300.0])
    n_om = MH.n_omega(e, t.range_at, t.stopping_at)
    pv = MH.pv_mev(e, PROTON.mass_mev)
    c = MH.c_omega(pv, 1e9, n_om)  # L_1 -> 0
    assert np.all((n_om > 0.6) & (n_om < 0.9)), n_om
    assert np.all((c > 0.11) & (c < 0.17)), c


def _ln_omega_proton(
    bethe: BetheStoppingSource, e: np.ndarray, electrons: str = "fano"
) -> np.ndarray:
    """Complete local proton ln Omega at fixed x/R (t ~ R): ln R - ln beta^2 - sum_i W_i ln(1.13 +
    3.76 (Z_i alpha)^2/beta^2) + delta_b, built from the module's own screening/Fano functions."""
    tab = bethe.table(WATER, PROTON)
    m = PROTON.mass_mev
    b2 = MH.beta2(e, m)
    out = np.log(tab.range_at(e)) - np.log(b2)
    for wgt, (zt, _a, _w) in zip(MH.screening_weights(MH.WATER, electrons), MH.WATER, strict=True):
        out = out - wgt * MH.ln_screening_factor(1.0, zt, b2)
    if electrons == "fano":
        out = out + MH.fano_delta_b(MH.WATER, b2)
    return np.asarray(out)


@pytest.mark.parametrize("electrons", ["fano", "bethe"])
def test_n_omega_matches_numerical_derivative_of_full_omega(
    bethe: BetheStoppingSource, electrons: str
) -> None:
    """n_Omega = d ln Omega_p / d ln(pv) of the COMPLETE proton Omega (Codex REVIEW-98315d4e):
    analytic result vs an independent central difference in ln E of the full ln Omega, at several
    proton energies and at the velocity-scaled energies of He/C/O at x/R 0.1-0.9, rel 1e-3."""
    tab = bethe.table(WATER, PROTON)
    m = PROTON.mass_mev
    energies = [30.0, 60.0, 100.0, 150.0, 250.0, 400.0]
    for name, e0_u in (("alpha", 150.0), ("carbon12", 290.0), ("oxygen16", 350.0)):
        proj = {"alpha": ALPHA, "carbon12": CARBON12, "oxygen16": OXYGEN16}[name]
        table = bethe.table(WATER, proj)
        r1 = float(table.range_at(e0_u))
        for frac in MH.I11_FRACTIONS:
            e_u = float(table.energy_from_range((1.0 - frac) * r1))  # MeV/u at depth x/R
            energies.append(e_u * m / (proj.mass_mev / proj.a))  # velocity-scaled proton E*
    e = np.array(energies)
    analytic = MH.n_omega(e, tab.range_at, tab.stopping_at, electrons=electrons)
    h = 1e-3
    lo, hi = e * (1.0 - h), e * (1.0 + h)
    num = (_ln_omega_proton(bethe, hi, electrons) - _ln_omega_proton(bethe, lo, electrons)) / (
        np.log(MH.pv_mev(hi, m)) - np.log(MH.pv_mev(lo, m))
    )
    assert analytic == pytest.approx(num, rel=1e-3), (e, analytic, num)
    # the omitted terms are not small: the range/kinematics-only value differs by > 10 %
    old = (e / (tab.stopping_at(e) * tab.range_at(e)) - (1 + e / (e + 2 * m) - 2 * e / (e + m))) / (
        1 + e / (e + 2 * m) - e / (e + m)
    )
    if electrons == "fano":
        assert np.max(np.abs(analytic / old - 1.0)) > 0.1


# --------------------------------------------------------------------------- I11 grid
_PLAN = Path(__file__).resolve().parents[2] / "validation/plans/v3-008-acceptance.md"


def _plan_i11_grid() -> tuple[dict[str, tuple[float, ...]], tuple[float, ...]]:
    """Energies and depth fractions parsed from the acceptance plan (Notation + row I11-abs)."""
    text = _PLAN.read_text(encoding="utf-8")
    std = re.search(r"\*\*Standard points\*\*: He (\d+), C (\d+), O (\d+) MeV/u", text)
    ends = re.search(r"\*\*Ends\*\*: He (\d+)/(\d+), C (\d+)/(\d+), O (\d+)/(\d+) MeV/u", text)
    row = re.search(
        r"^\| I11(?:-abs)? \|.*?x/R ∈ \{([^}]*)\}.*?standard points and ends", text, re.M
    )
    assert std and ends and row, "plan wording for the I11 grid changed; update the P1 driver grid"
    s, e = [float(g) for g in std.groups()], [float(g) for g in ends.groups()]
    grid = {
        "alpha": tuple(sorted({e[0], s[0], e[1]})),
        "carbon12": tuple(sorted({e[2], s[1], e[3]})),
        "oxygen16": tuple(sorted({e[4], s[2], e[5]})),
    }
    return grid, tuple(float(x) for x in row.group(1).split(","))


def test_driver_grid_equals_plan_i11_definition() -> None:
    """The P1 report driver evaluates exactly the I11 grid of the acceptance plan."""
    grid, fractions = _plan_i11_grid()
    assert {k: tuple(sorted(v)) for k, v in MH.I11_GRID_MEV_U.items()} == grid
    assert MH.I11_FRACTIONS == fractions == (0.1, 0.25, 0.5, 0.9)


# --------------------------------------------------------------------------- proton validation (LV)
_MATERIALS = {"Be": BERYLLIUM, "Al": ALUMINIUM, "Cu": COPPER, "Pb": LEAD}


def _tex() -> str:
    dirs = [GOTTSCHALK_CACHE]
    env = os.environ.get("IONMC_CACHE_DIR")
    if env:
        dirs.append(Path(env) / "objects")
    for d in dirs:
        if (d / GOTTSCHALK_SHA256).is_file():
            return gottschalk_tex(d)
    msg = (
        f"Gottschalk arXiv:0908.1413 source object sha256 {GOTTSCHALK_SHA256[:8]}... (registry "
        "reference data, arXiv e-print of 0908.1413) is not in the repository reference cache "
        f"({GOTTSCHALK_CACHE}) or $IONMC_CACHE_DIR/objects; materialize it (acquire/materialize "
        "the reference-data object with that sha256 into the cache, or set IONMC_CACHE_DIR to a "
        "cache holding it) before the P1 proton validation can run"
    )
    if os.environ.get("IONMC_REQUIRE_DATA") == "1":
        pytest.fail("IONMC_REQUIRE_DATA=1 but " + msg)
    pytest.skip(msg)


def proton_deviations(
    bethe: BetheStoppingSource, name: str, electrons: str = "fano"
) -> dict[float, float]:
    rows, _, _ = parse_tables(_tex())
    m = _MATERIALS[name]
    r1, energy_at = _path(bethe, m, 158.6)
    out = {}
    for frac, (theta_hanson_mrad, _dm) in sorted(rows[name].items()):
        th = MH.theta0_hanson(
            1.0, PROTON.mass_mev, _elements(m), frac * r1, energy_at, electrons=electrons
        )
        out[frac] = 1e3 * th / theta_hanson_mrad - 1.0
    return out


_P1_OUTCOME = {
    "Pb": "FAILS P1 with Fano (P1b; max +1.32 % at x/R1 0.97, +1.10 % at 0.9, <= 0.61 % for "
    "x/R1 <= 0.5): residual grows toward the end of range; ICRU 49 vs ionmc Bethe range table "
    "suspected (OPEN 2)",
}


@pytest.mark.parametrize(
    "name",
    [
        pytest.param(n, marks=pytest.mark.xfail(strict=True, reason=_P1_OUTCOME[n]))
        if n in _P1_OUTCOME
        else n
        for n in _MATERIALS
    ],
)
def test_p1_proton_reference_against_gottschalk_theta_hanson(
    bethe: BetheStoppingSource, name: str
) -> None:
    """Plan P1: |theta_MH / theta_Hanson - 1| <= 1.0 % at every tabulated x/R1 (0.001-0.97),
    reference with the Fano correction (P1b)."""
    dev = proton_deviations(bethe, name)
    assert max(abs(v) for v in dev.values()) <= PROTON_TOLERANCE, dev


def test_p1_bethe_variant_record(bethe: BetheStoppingSource) -> None:
    """Record of P1 (no Fano, Z(Z+1)): Be max +5.40 %, the low-Z excess removed by Fano (P1b)."""
    dev = proton_deviations(bethe, "Be", electrons="bethe")
    assert max(dev.values()) == pytest.approx(0.0540, abs=5e-4)
