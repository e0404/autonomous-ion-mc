"""Tests for the analytic and tabulated electronic stopping power (offline)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from ionmc.data.nist_star import StarTable, parse_star_text
from ionmc.materials import AIR, WATER, Material, fractions_from_atom_counts, water
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, PROTON
from ionmc.physics.stopping import (
    BetheOptions,
    BetheStoppingSource,
    NistStarStoppingSource,
    StoppingTable,
    bethe_mass_stopping,
    bloch_term,
    build_table,
    kinematics,
    pierce_blann_charge,
)

# Anchor values: ICRU Report 90 (2016) electronic stopping power of liquid water
# (I = 78 eV), as embedded in the Geant4 source file G4ICRU90StoppingData.cc, v11.4.2
# (Geant4 Software License; see src/ionmc/THIRD_PARTY_NOTICES.md). Arrays e1_proton and
# e1_alpha (index 1 = G4_WATER). Tuples: (kinetic energy MeV, S_el MeV cm2/g); the alpha
# energies are total kinetic energies of the alpha particle. All are >= 10 MeV/u.
# Tolerances are pre-chosen (0.3 % protons, 0.5 % alphas) and were not tuned to the model.
ICRU90_PROTON_ANCHORS = [
    (10.0, 45.32), (20.0, 25.89), (50.0, 12.38), (100.0, 7.250),
    (150.0, 5.417), (200.0, 4.470), (300.0, 3.504), (500.0, 2.731),
]  # fmt: skip
ICRU90_ALPHA_ANCHORS = [
    (40.0, 180.2), (80.0, 103.0), (200.0, 49.23), (400.0, 28.85), (800.0, 17.80), (1000.0, 15.50),
]  # fmt: skip
PROTON_TOLERANCE = 0.003
ALPHA_TOLERANCE = 0.005

W78 = water(78.0)
W75 = water(75.0)


@pytest.mark.parametrize(("e", "s"), ICRU90_PROTON_ANCHORS)
def test_proton_stopping_matches_icru90_anchors(e: float, s: float) -> None:
    assert abs(bethe_mass_stopping(e, PROTON, W78)[0] / s - 1.0) < PROTON_TOLERANCE


@pytest.mark.parametrize(("e", "s"), ICRU90_ALPHA_ANCHORS)
def test_alpha_stopping_matches_icru90_anchors(e: float, s: float) -> None:
    assert abs(bethe_mass_stopping(e / ALPHA.a, ALPHA, W78)[0] / s - 1.0) < ALPHA_TOLERANCE


def test_csda_range_equals_independent_integral() -> None:
    """R(200 MeV) - R(100 MeV) from the table equals an independent Simpson integral of 1/S."""
    table = BetheStoppingSource().table(W78, PROTON)
    e = np.linspace(100.0, 200.0, 2001)
    f = 1.0 / bethe_mass_stopping(e, PROTON, W78)
    h = e[1] - e[0]
    simpson = h / 3.0 * (f[0] + f[-1] + 4.0 * f[1:-1:2].sum() + 2.0 * f[2:-1:2].sum())
    assert abs((table.range_at(200.0) - table.range_at(100.0)) / simpson - 1.0) < 1e-4


def test_larger_I_lowers_stopping_and_lengthens_range() -> None:
    e = np.array([10.0, 100.0, 250.0])
    s75 = bethe_mass_stopping(e, PROTON, water(75.0))
    s78 = bethe_mass_stopping(e, PROTON, water(78.0))
    assert np.all(s78 < s75)
    r75 = BetheStoppingSource().table(water(75.0), PROTON).range_at(200.0)
    r78 = BetheStoppingSource().table(water(78.0), PROTON).range_at(200.0)
    assert r78 > r75


def test_domain_limit() -> None:
    with pytest.raises(ValueError):
        bethe_mass_stopping(0.99, PROTON, WATER)
    assert bethe_mass_stopping(1.0, PROTON, WATER)[0] > 0.0


def test_table_units_and_monotonicity() -> None:
    t = BetheStoppingSource().table(WATER, PROTON)
    assert t.energy_per_u[0] == 1.0 and t.energy_per_u[-1] == pytest.approx(500.0)
    assert np.all(np.diff(t.energy_per_u) > 0)
    assert np.all(np.diff(t.csda_range_g_cm2) > 0)
    assert np.all(np.diff(t.range_mm) > 0)
    np.testing.assert_allclose(t.s_el_linear, t.s_el_mass * WATER.density_g_cm3 / 10.0)
    np.testing.assert_allclose(t.range_mm, t.csda_range_g_cm2 / WATER.density_g_cm3 * 10.0)
    assert t.metadata["I_eV"] == 78.0 and t.metadata["source"] == "bethe"
    # stopping power falls from 1 MeV/u to the minimum ionisation region
    assert np.all(np.diff(t.s_el_mass[t.energy_per_u < 300.0]) < 0)


def test_range_integration_converges_with_grid() -> None:
    fine = BetheStoppingSource(points_per_decade=400).table(WATER, PROTON)
    coarse = BetheStoppingSource(points_per_decade=100).table(WATER, PROTON)
    mid = BetheStoppingSource(points_per_decade=200).table(WATER, PROTON)
    err_coarse = abs(coarse.range_at(200.0) / fine.range_at(200.0) - 1.0)
    err_mid = abs(mid.range_at(200.0) / fine.range_at(200.0) - 1.0)
    assert err_coarse < 3e-4 and err_mid < err_coarse


@pytest.mark.parametrize("projectile", [PROTON, ALPHA, CARBON12])
def test_inverse_round_trip(projectile) -> None:  # type: ignore[no-untyped-def]
    t = BetheStoppingSource().table(WATER, projectile)
    e = np.geomspace(1.5, 400.0, 57)
    back = t.energy_from_range(t.range_at(e))
    assert np.max(np.abs(back - e) / e) < 1e-6
    with pytest.raises(ValueError):
        t.energy_from_range(1e-9)
    with pytest.raises(ValueError):
        t.stopping_at(600.0)


def test_stopping_at_reproduces_grid_values() -> None:
    t = BetheStoppingSource().table(WATER, PROTON)
    np.testing.assert_allclose(t.stopping_at(t.energy_per_u[::50]), t.s_el_mass[::50], rtol=1e-12)


def test_carbon_differs_from_z2_scaled_protons_by_barkas_and_bloch() -> None:
    e = np.array([10.0, 20.0, 50.0, 100.0])
    s_c = bethe_mass_stopping(e, CARBON12, WATER)
    s_p = bethe_mass_stopping(e, PROTON, WATER)
    dev = s_c / (36.0 * s_p) - 1.0
    # Full carbon stopping is below z^2-scaled protons (Bloch dominates the z^3/z^4 terms).
    assert np.all(dev < 0.0)
    assert np.all(np.abs(dev) > 0.003) and np.all(np.abs(dev) < 0.03)


def test_barkas_and_bloch_signs() -> None:
    e = np.array([10.0])
    full = bethe_mass_stopping(e, CARBON12, WATER, effective_charge=False)
    no_barkas = bethe_mass_stopping(e, CARBON12, WATER, barkas=False, effective_charge=False)
    no_bloch = bethe_mass_stopping(e, CARBON12, WATER, bloch=False, effective_charge=False)
    assert full > no_barkas  # Barkas increases stopping
    assert full < no_bloch  # Bloch decreases stopping
    beta2, _, _ = kinematics(e, CARBON12)
    assert bloch_term(beta2, np.array([6.0]))[0] < 0.0


def test_effective_charge_tends_to_z() -> None:
    for proj in (CARBON12, OXYGEN16, ALPHA):
        beta = np.sqrt(kinematics(np.array([1.0, 10.0, 100.0, 400.0]), proj)[0])
        q = pierce_blann_charge(beta, proj) / proj.z
        assert q[0] < 0.995 and q[2] > 0.9999
        assert abs(q[-1] - 1.0) < 1e-6
        assert np.all(np.diff(q) > 0)
    e = np.array([2.0])
    low = bethe_mass_stopping(e, CARBON12, WATER)
    unscreened = bethe_mass_stopping(e, CARBON12, WATER, effective_charge=False)
    assert low < unscreened


def test_options_dataclass_and_nonwater_materials() -> None:
    opts = BetheOptions(shell=False)
    t = BetheStoppingSource(options=opts).table(WATER, PROTON)
    assert t.metadata["options"]["shell"] is False
    t_air = BetheStoppingSource().table(AIR, PROTON)
    assert (
        t_air.range_at(100.0)
        > 1000.0 * BetheStoppingSource().table(WATER, PROTON).range_at(100.0) / 1000.0
    )
    assert bethe_mass_stopping(100.0, PROTON, AIR)[0] < bethe_mass_stopping(100.0, PROTON, WATER)[0]


def _synthetic_star(program: str, energy_scale: float) -> StarTable:
    """Smooth synthetic STAR-like table from the analytic model (test of plumbing only)."""
    proj = PROTON if program == "PSTAR" else ALPHA
    e_u = np.geomspace(1.0, 600.0, 120)
    s = bethe_mass_stopping(e_u, proj, W75)
    text = [f"{program}: x", "WATER, LIQUID", "", "h", "h", "h", ""]
    r = np.cumsum(np.concatenate(([0.0], proj.a * np.diff(e_u) / s[1:])))
    for eu, si, ri in zip(e_u, s, r + 1e-3, strict=True):
        text.append(f"{eu * proj.a * energy_scale:.6E} {si:.6E} 0.0 {si:.6E} {ri:.6E} 0.0 1.0")
    return parse_star_text("\n".join(text))


def test_nist_star_source_builds_consistent_table() -> None:
    star = _synthetic_star("PSTAR", 1.0)
    src = NistStarStoppingSource(star)
    t = src.table(WATER, PROTON)
    assert t.metadata["source"] == "nist-star"
    ref = BetheStoppingSource().table(W75, PROTON)
    np.testing.assert_allclose(t.s_el_mass, ref.s_el_mass[: t.s_el_mass.size], rtol=1e-3)
    assert np.all(np.diff(t.csda_range_g_cm2) > 0)
    e = np.geomspace(1.5, 300.0, 20)
    assert np.max(np.abs(t.energy_from_range(t.range_at(e)) - e) / e) < 1e-6
    with pytest.raises(ValueError):
        src.table(WATER, ALPHA)
    with pytest.raises(ValueError):
        src.table(AIR, PROTON)
    astar = _synthetic_star("ASTAR", 1.0)
    ta = NistStarStoppingSource(astar, e_max_per_u=200.0).table(WATER, ALPHA)
    assert ta.energy_per_u[-1] == pytest.approx(200.0)


def test_counterfeit_water_rejected_and_renamed_water_accepted() -> None:
    star = _synthetic_star("PSTAR", 1.0)
    src = NistStarStoppingSource(star)
    fake = Material("water_fake", 1.0, {"C": 1.0}, 75.0)
    with pytest.raises(ValueError, match="liquid water only"):
        src.table(fake, PROTON)
    wrong_density = Material("water", 1.05, dict(WATER.mass_fractions), 75.0)
    with pytest.raises(ValueError):
        src.table(wrong_density, PROTON)
    renamed = Material(
        "my_h2o", 1.0, fractions_from_atom_counts({"H": 2, "O": 1}), 78.0, WATER.sternheimer
    )
    t = src.table(renamed, PROTON)
    assert t.metadata["requested_material"]["name"] == "my_h2o"
    assert t.metadata["requested_material"]["I_eV"] == 78.0
    assert t.metadata["effective_material"]["I_eV"] == 75.0


def test_table_metadata_carries_provenance() -> None:
    star = replace(
        _synthetic_star("PSTAR", 1.0),
        dataset_id="ds",
        version="v",
        sha256="a" * 64,
        retrieved_at="2026-01-01T00:00:00+00:00",
    )
    meta = NistStarStoppingSource(star).table(WATER, PROTON).metadata
    assert meta["dataset_id"] == "ds" and meta["sha256"] == "a" * 64
    assert meta["content_sha256"] == star.content_sha256 and len(meta["content_sha256"]) == 64
    assert meta["retrieved_at"] == "2026-01-01T00:00:00+00:00" and meta["version"] == "v"
    bmeta = BetheStoppingSource().table(WATER, PROTON).metadata
    assert bmeta["I_eV"] == 78.0 and bmeta["options"]["bloch"] is True
    assert "G4_WATER" in bmeta["material_source"]


def test_nonfinite_inputs_rejected() -> None:
    with pytest.raises(ValueError):
        bethe_mass_stopping(float("nan"), PROTON, WATER)
    with pytest.raises(ValueError):
        bethe_mass_stopping(np.array([10.0, np.inf]), PROTON, WATER)
    e = np.geomspace(1.0, 100.0, 10)
    s = bethe_mass_stopping(e, PROTON, WATER)
    bad = s.copy()
    bad[3] = np.nan
    with pytest.raises(ValueError):
        build_table(PROTON, WATER, e, bad, 1e-3, {})
    with pytest.raises(ValueError):
        build_table(PROTON, WATER, e, s, float("nan"), {})
    bad_e = e.copy()
    bad_e[2] = np.nan
    with pytest.raises(ValueError):
        build_table(PROTON, WATER, bad_e, s, 1e-3, {})
    t = build_table(PROTON, WATER, e, s, 1e-3, {})
    with pytest.raises(ValueError):
        t.stopping_at(float("nan"))
    with pytest.raises(ValueError):
        t.energy_from_range(float("nan"))
    with pytest.raises(ValueError):
        StoppingTable(PROTON, WATER, e, bad, bad * 0.1, t.csda_range_g_cm2, t.range_mm)


def _valid_arrays() -> dict[str, np.ndarray]:
    e = np.geomspace(1.0, 100.0, 10)
    s = bethe_mass_stopping(e, PROTON, WATER)
    t = build_table(PROTON, WATER, e, s, 1e-3, {})
    return {
        "energy_per_u": t.energy_per_u,
        "s_el_mass": t.s_el_mass,
        "s_el_linear": t.s_el_linear,
        "csda_range_g_cm2": t.csda_range_g_cm2,
        "range_mm": t.range_mm,
    }


@pytest.mark.parametrize(
    ("field", "mutate", "message"),
    [
        ("energy_per_u", lambda a: -a[::-1] * 1.0, "energy_per_u"),
        ("energy_per_u", lambda a: a[:-1], "length"),
        ("energy_per_u", lambda a: a * np.array([1.0] * 9 + [0.5]), "strictly increasing"),
        ("s_el_mass", lambda a: np.where(np.arange(a.size) == 2, -1.0, a), "s_el_mass"),
        ("s_el_linear", lambda a: a * 1.001, "s_el_linear"),
        ("range_mm", lambda a: a * 1.001, "range_mm"),
        ("csda_range_g_cm2", lambda a: a[::-1].copy(), "csda_range_g_cm2"),
        ("csda_range_g_cm2", lambda a: np.where(np.arange(a.size) == 0, -1e-3, a), "non-negative"),
        ("energy_per_u", lambda a: a.reshape(2, 5), "one-dimensional"),
        ("range_mm", lambda a: np.where(np.arange(a.size) == 4, np.nan, a), "NaN"),
    ],
)
def test_stopping_table_fail_closed(field: str, mutate, message: str) -> None:  # type: ignore[no-untyped-def]
    arrays = _valid_arrays()
    arrays[field] = mutate(arrays[field])
    with pytest.raises(ValueError, match=message):
        StoppingTable(PROTON, WATER, **arrays)


def test_stopping_table_minimum_length() -> None:
    arrays = {k: v[:1] for k, v in _valid_arrays().items()}
    with pytest.raises(ValueError, match="at least two"):
        StoppingTable(PROTON, WATER, **arrays)
    StoppingTable(PROTON, WATER, **_valid_arrays())


def test_negative_increasing_energies_rejected() -> None:
    e = np.array([-3.0, -2.0, -1.0])
    with pytest.raises(ValueError):
        build_table(PROTON, WATER, e, np.array([5.0, 4.0, 3.0]), 0.0, {})
