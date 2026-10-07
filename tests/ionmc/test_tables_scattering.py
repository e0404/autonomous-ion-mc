"""U2-U5 and sampling-function checks: transport tables, X_S, scattering power integrals."""

import math
from pathlib import Path

import numpy as np
import pytest
import warp as wp

from ionmc._wpfunc import python_twin
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
    U4B_TOLERANCE,
    ProtonPath,
    gottschalk_tex,
    highland_cross_check,
    parse_tables,
    radiation_length_g_cm2,
    u4_deviations,
    u4b_deviations,
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
# Python-side evaluation uses the pure-Python twins; Warp functions only run inside kernels
EMP = python_twin(make_em)
TFP = python_twin(make_transport_funcs)
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
        i, f = TFP.log_bin_index(float(e), float(t.ln_e0[0]), float(t.inv_dln_e[0]), t.n_e)
        # the range is the exact closed form in the bin (V3-003D), not a log-log interpolation
        r = float(
            TFP.range_in_bin(
                float(t.r_mass[0, i]),
                float(t.f_mass[0, i]),
                float(t.d_f[0, i]),
                1.0 / float(t.inv_dln_e[0]),
                f,
            )
        )
        assert r == pytest.approx(t.range_g_cm2(0, e), rel=1e-13)
        # the stopping power read is unchanged: log-log interpolation of ln S
        row = t.ln_s_mass[0]
        s = float(TFP.interp_exp(float(row[i]), float(row[i + 1]), f))
        assert s == pytest.approx(t.stopping_mass(0, e), rel=1e-13)


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


@pytest.mark.parametrize("name", list(_U4_MATERIALS))
def test_u4b_theta_dm_against_theta_hanson(bethe: BetheStoppingSource, name: str) -> None:
    """U4b: our theta_dM against the paper's theta_Hanson column itself, within 4.5 %.

    theta_Hanson (Moliere/Fano/Hanson theory) is independent of the dM fit. The paper's own
    T_dM is within 2.74 % of it over the frozen points (largest |dM %|) and U4 bounds our
    reproduction of T_dM at 1.5 %, so (1.0274)(1.015) - 1 = 4.3 % -> frozen tolerance 4.5 %."""
    rows, _, _ = parse_tables(_gottschalk_tex())
    material = _U4_MATERIALS[name]
    tables = TransportTables.from_stopping_tables([bethe.table(material, PROTON)])
    deviations, paper_dm = u4b_deviations(tables, material, rows[name])
    assert paper_dm <= 2.74 + 1e-9, paper_dm
    assert max(abs(d) for d in deviations.values()) <= U4B_TOLERANCE, deviations


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
            EMP.scattering_variance_birth(
                float(500.0),
                float(540.0),
                float(540.0),
                float(1.0),
                float(0.02),
                float(1.0),
                float(1.0),
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
    "ratio", [30.0, 8.0, 4.0, 3.0, 2.99, 2.0, 1.0, 0.5, 0.3], ids=lambda r: f"ratio{r}"
)
@pytest.mark.parametrize("mean", [0.05, 2.0])
def test_straggling_gamma_moments_exact_and_gaussian_variance_loss_bounded(
    ratio: float, mean: float
) -> None:
    """Gamma branch (ratio < 3): mean and variance of the step energy loss equal the Bohr mean
    and variance (sigma = mean / ratio) within three standard errors over 4e5 draws. Gaussian
    branch (ratio >= 3): the mean is preserved, but the clamp to [0, 2 mean] deliberately
    lowers the variance (about 0.5 % at ratio 3, less at larger ratios); the assertion allows a
    1 % variance deficit there, so this branch is NOT a two-moment check. Losses are positive and
    the branch boundary is 3. Warp's Philox is bit-identical to the Python-integer Philox
    (test_rng)."""
    n = 400_000
    sigma = mean / ratio
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
    se_mean = x.std() / math.sqrt(n)
    assert abs(x.mean() - mean) <= 3.0 * se_mean + 1e-3 * sigma * (ratio >= 3.0)
    centred = (x - x.mean()) ** 2
    se_var = centred.std() / math.sqrt(n)
    allowance = 0.01 * sigma**2 if ratio >= 3.0 else 0.0
    assert abs(x.var() - sigma**2) <= 3.0 * se_var + allowance, (x.var(), sigma**2)
    if ratio >= 3.0:
        assert x.max() <= 2.0 * mean and attempts.numpy().max() == 1  # always accepted
    else:
        assert attempts.numpy().max() >= 1


def test_straggling_branch_boundary() -> None:
    """Exactly at ratio 3 the Gaussian branch is used (always accepted, clamped), just below
    it the Gamma branch (a rejection is possible)."""
    f = EMP.straggle_attempt
    sigma = 1.0
    hi = f(float(3.0), float(sigma**2), float(0.5), float(0.0), float(0.5), float(0.5))
    assert hi[1] == 1 and float(hi[0]) == pytest.approx(
        3.0 + sigma * math.sqrt(2.0 * math.log(2.0))
    )
    clamped = f(float(3.0), float(sigma**2), float(1e-300), float(0.5), float(0.5), float(0.5))
    assert float(clamped[0]) == 0.0 and clamped[1] == 1  # z = -37: clamped at 0
    upper = f(float(3.0), float(sigma**2), float(1e-300), float(0.0), float(0.5), float(0.5))
    assert float(upper[0]) == 6.0  # clamped at 2 mean
    # below the boundary a Gamma sample with a failing acceptance test is rejected
    rejected = f(float(2.9), float(sigma**2), float(1e-300), float(0.5), float(0.5), float(0.5))
    assert rejected[1] == 0  # 1 + c x <= 0: Marsaglia-Tsang rejects
    assert f(float(0.0), float(1.0), float(0.5), float(0.5), float(0.5), float(0.5)) == (0.0, 1)
    assert f(float(2.0), float(0.0), float(0.5), float(0.5), float(0.5), float(0.5)) == (2.0, 1)


def test_polar_deflection_and_rotation_statistics() -> None:
    """theta^2 of the 2-D Gaussian has mean 2 var; rotation preserves norm and angle."""
    rng = np.random.default_rng(5)
    var = 1.7e-4
    u = rng.uniform(0.0, 1.0, 20000)
    u = np.clip(u, 1e-12, 1 - 1e-12)
    th = np.array([float(EMP.polar_deflection(float(var), float(v))) for v in u[:2000]])
    assert np.mean(th**2) == pytest.approx(2.0 * var, rel=0.1)
    assert float(EMP.polar_deflection(float(1.0e3), float(0.5))) == pytest.approx(math.pi)
    v3 = EMP.vec3
    for d in ([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [0.36, 0.48, 0.8], [1.0, 0.0, 0.0]):
        out = EMP.rotate_dir(v3(float(d[0]), float(d[1]), float(d[2])), float(0.3), float(1.1))
        o = np.array([float(c) for c in out])
        assert np.linalg.norm(o) == pytest.approx(1.0, abs=1e-14)
        assert float(np.dot(o, d)) == pytest.approx(math.cos(0.3), abs=1e-12)
