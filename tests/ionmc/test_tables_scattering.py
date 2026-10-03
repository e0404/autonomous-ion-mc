"""U2-U5 and sampling-function checks: transport tables, X_S, scattering power integrals."""

import math
from pathlib import Path

import numpy as np
import pytest
import warp as wp

from ionmc.materials import (
    ALUMINIUM,
    BERYLLIUM,
    COPPER,
    LEAD,
    WATER,
    Material,
)
from ionmc.physics.em import make_em
from ionmc.physics.kinematics import make_kinematics
from ionmc.physics.projectiles import PROTON
from ionmc.physics.scattering import (
    inverse_scattering_length_cm2_per_g,
    scattering_length_g_cm2,
)
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.rng.philox import make_philox
from ionmc.transport.funcs import make_transport_funcs
from ionmc.transport.mcs_checks import (
    GOTTSCHALK_SHA256,
    U4_FRACTIONS,
    U4_TOLERANCE,
    ProtonPath,
    gottschalk_tex,
    highland_cross_check,
    parse_tables,
    radiation_length_g_cm2,
    u4_deviations,
    u5_deviations,
    u5_negative_control,
    xs_deviations,
)
from ionmc.transport.tables import TransportTables

wp.config.log_level = wp.LOG_WARNING

F64 = wp.float64
EM = make_em(F64)
KIN = make_kinematics(F64)
TF = make_transport_funcs(F64)
PH = make_philox(F64)
M_P = PROTON.mass_mev


@pytest.fixture(scope="module")
def water_tables(bethe: BetheStoppingSource) -> TransportTables:
    return TransportTables.from_stopping_tables([bethe.table(WATER, PROTON)])


# ---------------------------------------------------------------------------------------------
# U2: table round trip and monotonicity
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("material", [WATER, COPPER], ids=lambda m: m.name)
def test_u2_round_trip_and_monotonicity(bethe: BetheStoppingSource, material: Material) -> None:
    source_table = bethe.table(material, PROTON)
    t = TransportTables.from_stopping_tables([source_table])
    assert t.n_e >= 200 * math.log10(t.e_max_mev[0] / t.e_min_mev[0])
    assert t.n_r >= 800 * math.log10(t.r_max_g_cm2[0] / t.r_min_g_cm2[0])  # default 800 / decade
    assert np.all(np.diff(t.ln_r_mass[0]) > 0.0) and np.all(np.diff(t.ln_e_of_r[0]) > 0.0)
    energies = np.geomspace(t.e_min_mev[0], t.e_max_mev[0], 20001)
    ranges = np.array([t.range_g_cm2(0, e) for e in energies])
    assert np.all(np.diff(ranges) > 0.0)  # R strictly monotone as evaluated
    back = np.array([t.energy_from_range(0, r) for r in ranges])
    assert np.all(np.diff(back) > 0.0)  # Rinv strictly monotone as evaluated
    assert np.max(np.abs(back / energies - 1.0)) <= 1e-5
    # ranges follow the source table (same stopping powers, integrated on the same grid)
    probe = np.geomspace(1.5, 400.0, 40)
    for e in probe:
        assert t.range_g_cm2(0, e) == pytest.approx(float(source_table.range_at(e)), rel=1e-4)
        assert t.stopping_mass(0, e) == pytest.approx(float(source_table.stopping_at(e)), rel=1e-4)
    # R is the integral of the same stopping power: dR/dE = 1/S
    for e in (3.0, 20.0, 100.0, 250.0):
        h = 1e-3 * e
        dr = (t.range_g_cm2(0, e + h) - t.range_g_cm2(0, e - h)) / (2 * h)
        assert dr == pytest.approx(1.0 / t.stopping_mass(0, e), rel=2e-3)


def test_u2_shared_funcs_reproduce_table_reads(water_tables: TransportTables) -> None:
    t = water_tables
    for e in np.geomspace(1.2, 450.0, 60):
        i, f = TF.log_bin_index(F64(e), F64(t.ln_e0[0]), F64(t.inv_dln_e[0]), t.n_e)
        row = t.ln_r_mass[0]
        r = float(TF.interp_exp(F64(row[i]), F64(row[i + 1]), f))
        assert r == pytest.approx(t.range_g_cm2(0, e), rel=1e-13)


def test_tables_reject_non_protons_and_clamp_outside(bethe: BetheStoppingSource) -> None:
    from ionmc.physics.projectiles import ALPHA

    with pytest.raises(ValueError, match="protons only"):
        TransportTables.from_stopping_tables([bethe.table(WATER, ALPHA)])
    with pytest.raises(ValueError):
        TransportTables.from_stopping_tables([])
    t = TransportTables.from_stopping_tables([bethe.table(WATER, PROTON)])
    assert t.energy_from_range(0, 1e-9) == pytest.approx(float(t.e_min_mev[0]), rel=1e-9)
    assert t.energy_from_range(0, 1e9) == pytest.approx(float(t.e_max_mev[0]), rel=1e-9)
    assert len(t.sha256) == 64
    w = t.to_warp("cpu", wp.float32)
    assert w.ln_s_mass.shape == (1, t.n_e) and w.ln_e_of_r.shape == (1, t.n_r)
    with pytest.raises(ValueError):
        t.to_warp("cpu", wp.int32)


# ---------------------------------------------------------------------------------------------
# U3: scattering length X_S against Gottschalk's tabulated values
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("material", "x_s"),
    [(WATER, 46.88), (BERYLLIUM, 92.60), (ALUMINIUM, 28.75), (COPPER, 14.62), (LEAD, 6.62)],
    ids=lambda v: v.name if isinstance(v, Material) else str(v),
)
def test_u3_scattering_length_matches_gottschalk(material: Material, x_s: float) -> None:
    """Natural logarithm in the formula (log10 would miss by far more than 0.3 %)."""
    assert scattering_length_g_cm2(material) == pytest.approx(x_s, rel=3e-3)
    assert scattering_length_g_cm2(material) * inverse_scattering_length_cm2_per_g(
        material
    ) == pytest.approx(1.0, rel=1e-14)


def test_radiation_length_helper_reproduces_pdg_values() -> None:
    """Sanity of the independent helper used below (PDG: water 36.08, Al 24.01, Cu 12.86,
    Pb 6.37 g/cm2)."""
    assert radiation_length_g_cm2(WATER) == pytest.approx(36.08, rel=3e-3)
    assert radiation_length_g_cm2(ALUMINIUM) == pytest.approx(24.01, rel=3e-3)
    assert radiation_length_g_cm2(COPPER) == pytest.approx(12.86, rel=3e-3)
    assert radiation_length_g_cm2(LEAD) == pytest.approx(6.37, rel=3e-3)


# ---------------------------------------------------------------------------------------------
# U4: theta_dM(x) against the tabulated theta_Hanson (1 + dM %) of Gottschalk (2010)
# (helpers in ionmc.transport.mcs_checks, shared with validation/scripts/transport/mcs_checks.py)
# ---------------------------------------------------------------------------------------------
def _gottschalk_tex() -> str:
    try:
        return gottschalk_tex()
    except FileNotFoundError as exc:
        pytest.skip(f"{exc}; the U4 table comparison runs locally (LV), not in CI")


_U4_MATERIALS = {"Be": BERYLLIUM, "Al": ALUMINIUM, "Cu": COPPER, "Pb": LEAD}


@pytest.mark.parametrize("name", list(_U4_MATERIALS))
def test_u4_theta_dm_against_gottschalk_table(bethe: BetheStoppingSource, name: str) -> None:
    """Frozen U4: theta_dM by quadrature of T_dM along the CSDA path (158.6 MeV) at
    x/R1 = 0.01, 0.1, 0.5, 0.9 within 1.5 % of theta_Hanson (1 + dM %/100). The reference is
    parsed at run time from the cached LaTeX source of arXiv:0908.1413; the test is skipped
    when that file is absent. The table has no water block (Be, Al, Cu, Pb only)."""
    rows, _, _ = parse_tables(_gottschalk_tex())
    material = _U4_MATERIALS[name]
    tables = TransportTables.from_stopping_tables([bethe.table(material, PROTON)])
    assert set(rows[name]) >= set(U4_FRACTIONS)
    deviations, _ = u4_deviations(tables, material, rows[name])
    assert max(abs(d) for d in deviations.values()) <= U4_TOLERANCE, deviations


def test_scattering_length_against_gottschalk_table() -> None:
    """Our X_S against ``tbl:LS`` (water, Be, Al, Cu, Pb) within 0.5 %."""
    _, _, x_s = parse_tables(_gottschalk_tex())
    ours = {"Be": BERYLLIUM, "H2O": WATER, "Al": ALUMINIUM, "Cu": COPPER, "Pb": LEAD}
    for name, dev in xs_deviations(x_s, ours).items():
        assert abs(dev) <= 0.005, (name, dev)


def test_gottschalk_loader_fails_closed(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        gottschalk_tex(tmp_path)
    (tmp_path / GOTTSCHALK_SHA256).write_bytes(b"not the pinned file")
    with pytest.raises(ValueError, match="sha256"):
        gottschalk_tex(tmp_path)


def test_u4_generalised_highland_cross_check(bethe: BetheStoppingSource) -> None:
    """Cache-independent cross-check (a related model, not U4): theta_dM within 6 % of the
    generalised Highland formula of Gottschalk et al. (1993) with a Tsai radiation length, for
    water at 158.6 MeV; the deviations are systematic in depth (about +3.7 % at x/R1 = 0.01,
    -3.4 % at 0.9)."""
    tables = TransportTables.from_stopping_tables([bethe.table(WATER, PROTON)])
    for frac, dev in highland_cross_check(tables, WATER).items():
        assert abs(dev) < 0.06, frac


# ---------------------------------------------------------------------------------------------
# U5: step independence of the MCS integrator
# ---------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def u5_deviation(water_tables: TransportTables) -> dict[tuple[float, float], float]:
    """Relative deviation of the stepped sum from the quadrature, keyed by (x/R1, s [mm])."""
    return u5_deviations(water_tables)


def test_u5_frozen_criterion(u5_deviation: dict[tuple[float, float], float]) -> None:
    """Frozen criterion U5: |d theta^2 / theta^2| <= 2e-3 for every s and x/R1 >= 0.05."""
    worst = max(u5_deviation.items(), key=lambda kv: abs(kv[1]))
    assert abs(worst[1]) <= 2e-3, f"worst (x/R1, s) = {worst[0]}: {worst[1]:+.2e}"


def test_birth_variance_matches_fine_quadrature(water_tables: TransportTables) -> None:
    """``scattering_variance_birth`` against the quadrature of the clamped T over the first
    step (150 MeV, water): relative error <= 1e-3 for 0.1 <= s <= 1 mm. For s = 0.01 mm the
    step changes pv by less than the round-trip accuracy of the tables (1e-5 relative), so the
    1 - (pv/p1v1)^2 of the quadrature is no longer linear from 0 and the error is about 2e-3."""
    path = ProtonPath(water_tables, WATER, 150.0)
    for s in (0.1, 0.25, 0.5, 1.0):
        exact = path.theta2_quadrature(s)
        assert abs(path.birth_variance(s) / exact - 1.0) <= 1e-3, s
    assert (
        float(
            EM.scattering_variance_birth(
                F64(500.0), F64(540.0), F64(540.0), F64(1.0), F64(0.02), F64(1.0), F64(1.0)
            )
        )
        == 0.0
    )  # zero-length step: nothing to scatter


def test_u5_negative_control_per_step_highland_depends_on_step(
    water_tables: TransportTables,
) -> None:
    """Per-step Highland (log term evaluated per step) varies by far more than 5 % with s."""
    nc = u5_negative_control(water_tables)
    assert nc["per_step_highland_relative_spread"] >= 0.05
    # the engine's stepped T_dM varies far less over the steps up to 1 mm (instrument power)
    assert nc["stepped_tdm_relative_spread_up_to_1mm"] < 1e-3


# ---------------------------------------------------------------------------------------------
# sampling functions: moments of the straggling regimes and the polar deflection
# ---------------------------------------------------------------------------------------------
@wp.kernel(module="unique")
def _straggle_kernel(
    mean: wp.float64,
    var: wp.float64,
    key: wp.vec2ui,
    loss: wp.array(dtype=wp.float64),  # type: ignore[valid-type]
    attempts: wp.array(dtype=wp.int32),  # type: ignore[valid-type]
) -> None:
    i = wp.tid()
    value = wp.float64(0.0)
    used = int(0)
    done = int(0)
    for k in range(64):
        if done == 0:
            r = PH.philox_block(wp.uint32(i), wp.uint32(0), wp.uint32(k), wp.uint32(0), key)
            lw, ok = EM.straggle_attempt(
                mean, var, PH.u01(r[0]), PH.u01(r[1]), PH.u01(r[2]), PH.u01(r[3])
            )
            used = k + 1
            if ok == 1:
                value = lw
                done = 1
    loss[i] = value
    attempts[i] = used


@pytest.mark.parametrize(
    ("ratio", "kind"),
    [
        (12.0, "sigma"),  # Gaussian, truncation negligible: variance sigma^2
        (3.0, "truncated"),  # Gaussian cut at 2 sigma: variance slightly reduced
        (1.5, "sigma"),  # Gamma, n = 2.25
        (0.8, "sigma"),  # Gamma, n = 0.64 < 1 (u3^(1/n) branch)
        (0.05, "uniform"),  # uniform on (0, 2 mean): variance mean^2 / 3
    ],
)
def test_straggling_regimes_preserve_mean_and_variance(ratio: float, kind: str) -> None:
    n = 400_000
    mean, sigma = 2.0, 2.0 / ratio
    loss = wp.zeros(n, dtype=wp.float64, device="cpu")
    attempts = wp.zeros(n, dtype=wp.int32, device="cpu")
    wp.launch(
        _straggle_kernel,
        dim=n,
        inputs=[F64(mean), F64(sigma**2), wp.vec2ui(11, 0), loss, attempts],
        device="cpu",
    )
    x = loss.numpy()
    assert np.all(np.isfinite(x)) and np.all(x >= 0.0) and attempts.numpy().max() < 64
    assert abs(x.mean() - mean) < 5 * x.std() / math.sqrt(n)  # statistical error of the mean
    if kind == "sigma":
        assert x.var() == pytest.approx(sigma**2, rel=0.03)
    elif kind == "uniform":
        assert x.var() == pytest.approx(mean**2 / 3.0, rel=0.01)
        assert x.max() < 2.0 * mean
    else:
        assert 0.8 * sigma**2 < x.var() < sigma**2
        assert x.max() < 2.0 * mean


def test_polar_deflection_and_rotation_statistics() -> None:
    """theta^2 of the 2-D Gaussian has mean 2 var; rotation preserves norm and angle."""
    rng = np.random.default_rng(5)
    var = 1.7e-4
    u = rng.uniform(0.0, 1.0, 20000)
    u = np.clip(u, 1e-12, 1 - 1e-12)
    th = np.array([float(EM.polar_deflection(F64(var), F64(v))) for v in u[:2000]])
    assert np.mean(th**2) == pytest.approx(2.0 * var, rel=0.1)
    assert float(EM.polar_deflection(F64(1.0e3), F64(0.5))) == pytest.approx(math.pi)
    v3 = EM.vec3
    for d in ([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [0.36, 0.48, 0.8], [1.0, 0.0, 0.0]):
        out = EM.rotate_dir(v3(F64(d[0]), F64(d[1]), F64(d[2])), F64(0.3), F64(1.1))
        o = np.array([float(c) for c in out])
        assert np.linalg.norm(o) == pytest.approx(1.0, abs=1e-14)
        assert float(np.dot(o, d)) == pytest.approx(math.cos(0.3), abs=1e-12)
