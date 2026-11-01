"""P11 (V3-008): the exact-Omega mapping ``validation/scripts/ion_em/exact_omega.py``.

Decision 0043 section 3. Rows I11-map (deterministic self-consistency: z = 1 identity, root
residual <= 1e-10, Omega identity, tabulation bound <= 0.1 %) are asserted. Rows I11-rel (gating
2 %, evaluated pre-freeze, SEEN, not blind) and I11-abs (informative) are only COMPUTED here: the
test checks that the evaluation runs and is finite and prints the maximum |ratio - 1|; the gate is a
plan row, never a unit-test assertion.
"""

import importlib.util
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc.physics.em import FDM_COEFFICIENTS

_DIR = Path(__file__).resolve().parents[2] / "validation/scripts/ion_em"
if "moliere_hanson" not in sys.modules:
    _s = importlib.util.spec_from_file_location("moliere_hanson", _DIR / "moliere_hanson.py")
    assert _s is not None and _s.loader is not None
    _m = importlib.util.module_from_spec(_s)
    sys.modules["moliere_hanson"] = _m
    _s.loader.exec_module(_m)
_SPEC = importlib.util.spec_from_file_location("exact_omega", _DIR / "exact_omega.py")
assert _SPEC is not None and _SPEC.loader is not None
X = importlib.util.module_from_spec(_SPEC)
sys.modules["exact_omega"] = X
_SPEC.loader.exec_module(X)
MH = X.MH

IONS = ("alpha", "carbon12", "oxygen16")


@pytest.fixture(scope="module")
def proton() -> Any:
    return X.bethe_path("proton", X.PROTON_EMAX_PER_U)


@pytest.fixture(scope="module")
def ions() -> dict[str, Any]:
    return {name: X.bethe_path(name) for name in IONS}


def _state(path: Any, e_u: float, frac: float) -> tuple[float, float]:
    """(pv, p1v1) of the ion at depth fraction ``frac`` of its range, initial energy ``e_u``."""
    r1 = float(path.range_of(np.float64(e_u)))
    e = float(path.energy_of(np.float64((1.0 - frac) * r1)))
    return X.pv_of_path(path, e), X.pv_of_path(path, e_u)


def test_proton_fdm_twin_coefficients() -> None:
    """The numpy f_dM used by the model is the project's proton f_dM (same coefficients)."""
    assert tuple(MH.FDM) == tuple(FDM_COEFFICIENTS)


def test_b_slab_matches_reference_moliere_slab(proton: Any) -> None:
    """b of the slab equals ``moliere_hanson.moliere_slab`` (n = 400, same path)."""
    e1, frac = 150.0, 0.5
    r1 = float(proton.range_of(np.float64(e1)))

    def energy_at(t: np.ndarray) -> np.ndarray:
        return np.asarray(proton.energy_of(r1 - t))

    ref = MH.moliere_slab(1.0, MH.M_P_MEV, MH.WATER, frac * r1, energy_at, n=400)["b"]
    pv, p1v1 = _state(proton, e1, frac)
    assert X.b_proton(proton, pv, p1v1, n=400) == pytest.approx(ref, abs=1e-9)
    assert X.b_proton(proton, pv, p1v1) == pytest.approx(ref, abs=1e-6)  # N_SLAB = 64


def test_z1_mp_reproduces_proton_fdm_exactly(proton: Any) -> None:
    """I11-map: z = 1, M = m_p gives pv_eq = pv and the proton f_dM, to <= 1e-12 (20 points)."""
    rng = np.random.default_rng(11)
    for _ in range(20):
        e1 = float(rng.uniform(30.0, 480.0))
        frac = float(rng.uniform(0.05, 0.95))
        pv, p1v1 = _state(proton, e1, frac)
        peq = X.pv_eq(proton, proton, pv, p1v1)
        assert abs(peq / pv - 1.0) <= 1e-12
        f_map = X.f_dm_equivalent(proton, proton, pv, p1v1)
        f_p = float(MH.fdm(np.array([pv]), p1v1)[0])
        assert abs(f_map / f_p - 1.0) <= 1e-12


def test_z1_theta_equals_proton_model(proton: Any) -> None:
    """theta0_ion with the proton as the ion equals the proton model's theta0 (z = 1 path)."""
    for e1, frac in ((100.0, 0.5), (200.0, 0.9)):
        a = X.theta0_ion(proton, proton, e1, frac, n=60)
        b = X.theta0_proton_model(proton, e1, frac, n=60)
        assert abs(a / b - 1.0) <= 1e-12


def test_root_residual_on_the_i11_grid(proton: Any, ions: dict[str, Any]) -> None:
    """I11-map: |b_p(pv_eq) - b_ion| <= 1e-10 at every I11 grid state and 7 more depths."""
    worst = 0.0
    for name, energies in MH.I11_GRID_MEV_U.items():
        for e_u in energies:
            for frac in (*MH.I11_FRACTIONS, 0.02, 0.05, 0.2, 0.7, 0.8, 0.95, 0.97):
                pv, p1v1 = _state(ions[name], e_u, frac)
                peq = X.pv_eq(ions[name], proton, pv, p1v1)
                res = X.residual(ions[name], proton, pv, p1v1, peq)
                worst = max(worst, abs(res))
                assert abs(res) <= X.RESIDUAL_TOL, (name, e_u, frac, res)
    print(f"max root residual on grid: {worst:.2e}")


def test_root_fails_closed_without_bracket(proton: Any, ions: dict[str, Any]) -> None:
    with pytest.raises(ValueError, match="no bracket"):
        X.brent(lambda x: x * x + 1.0, -1.0, 1.0)
    with pytest.raises(ValueError, match="pv below|below p1v1|need 0"):
        X.pv_eq(ions["alpha"], proton, 100.0, 90.0)


def test_omega_identity_at_equal_velocity(proton: Any) -> None:
    """Omega_ion/Omega_p = A*/F_scr at equal velocity and x/R, hence (Fano cancels)
    b_ion - b_p = ln A* - <ln F_scr>_(1/(pv)^2 slab average), exactly, on a path with the exact
    scaling R_ion(E_u) = (A*/z^2) R_p(E_u)."""
    for z, a_star in ((2.0, 3.97), (6.0, 11.9), (8.0, 15.85)):
        ion = X.scaled_path(proton, z, a_star)
        for electrons in ("fano", "bethe"):
            e1, frac = 150.0, 0.5
            pv_p, p1_p = _state(proton, e1, frac)
            pv_i, p1_i = _state(ion, e1, frac)
            assert pv_i == pytest.approx(a_star * pv_p, rel=1e-12)
            d_b = X.b_slab(ion, pv_i, p1_i, electrons=electrons) - X.b_slab(
                proton, pv_p, p1_p, electrons=electrons
            )
            # <ln F_scr> on the proton slab nodes (the ion nodes map onto the same fractions)
            r1 = float(proton.range_of(np.float64(e1)))
            r_e = float(proton.range_of(np.float64(_e_of(proton, pv_p))))
            t = r1 - r_e
            nodes, w = MH._gauss_nodes(t, X.N_SLAB)
            en = np.asarray(proton.energy_of(r_e + (t - nodes)))
            wt = w / MH.pv_mev(en, MH.M_P_MEV) ** 2
            ln_f = np.log(MH.f_scr(z, MH.beta2(en, MH.M_P_MEV), MH.WATER, electrons))
            expect = math.log(a_star) - float(np.sum(wt * ln_f) / np.sum(wt))
            assert d_b == pytest.approx(expect, abs=1e-10), (z, electrons)


def _e_of(path: Any, pv: float) -> float:
    return float(X.kinetic_from_pv(pv, path.mass_mev) / path.a)


def test_pv_eq_monotone_and_continuous_in_pv(proton: Any, ions: dict[str, Any]) -> None:
    """pv_eq(pv; p1v1 fixed) is strictly monotone and continuous (no jumps) on 120 points."""
    for name, e_u in (("alpha", 150.0), ("carbon12", 290.0), ("oxygen16", 350.0)):
        ion = ions[name]
        p1v1 = X.pv_of_path(ion, e_u)
        fr = np.linspace(0.05, 0.95, 120)
        vals = []
        pvs = []
        for f in fr:
            pv, _ = _state(ion, e_u, float(f))
            pvs.append(pv)
            vals.append(X.pv_eq(ion, proton, pv, p1v1))
        pvs_a, vals_a = np.array(pvs), np.array(vals)  # pv decreases along fr
        d = np.diff(vals_a)
        assert np.all(d > 0) or np.all(d < 0), name
        # continuity: adjacent relative steps vary smoothly (second difference of ln pv_eq small)
        lv = np.log(vals_a)
        lp = np.log(pvs_a)
        slope = np.diff(lv) / np.diff(lp)
        assert np.max(np.abs(np.diff(slope))) < 0.2, name


@pytest.fixture(scope="module")
def table_carbon(proton: Any, ions: dict[str, Any]) -> Any:
    ion = ions["carbon12"]
    pv_lo = X.pv_of_path(ion, 100.0)
    pv_hi = X.pv_of_path(ion, 250.0)
    return X.PvEqTable(
        ion,
        proton,
        (pv_lo, pv_hi),
        (math.log(1.02), math.log(1.3)),
        n_pv=30,
        n_r=30,
    )


def test_tabulation_theta_error_bound(proton: Any, ions: dict[str, Any], table_carbon: Any) -> None:
    """Allowed alternative: the interpolated pv_eq changes f_dM (hence T, and theta by half) by
    <= 0.1 % against the direct root on 40 random states of the table domain."""
    ion = ions["carbon12"]
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(40):
        pv = math.exp(rng.uniform(table_carbon.x[0], table_carbon.x[-1]))
        p1v1 = pv * math.exp(rng.uniform(table_carbon.y[0], table_carbon.y[-1]))
        direct = X.f_dm_equivalent(ion, proton, pv, p1v1)
        tab = X.f_dm_equivalent(ion, proton, pv, p1v1, mapping=table_carbon)
        worst = max(worst, abs(tab / direct - 1.0))
    print(f"tabulation max |f_dM error| = {worst:.2e}")
    assert worst <= 1e-3


def test_i11_rel_proton_is_at_equal_velocity() -> None:
    """I11-rel regression: beta_p(E_p) == beta_ion(E_u) to 1e-12 at every I11 grid point (the proton
    energy is E_u A m_p / M, not the ion's numeric MeV/u)."""
    for name, energies in MH.I11_GRID_MEV_U.items():
        ion = X.bethe_path(name)
        for e_u in energies:
            e_p = X.equal_velocity_proton_energy(ion, e_u)
            b_ion = X.beta_of_kinetic(e_u * ion.a, ion.mass_mev)
            b_p = X.beta_of_kinetic(e_p, MH.M_P_MEV)
            assert abs(b_p - b_ion) <= 1e-12, (name, e_u)
            assert abs(e_p / e_u - 1.0) > 1e-4  # genuinely not the numeric MeV/u


def test_i11_rel_and_abs_are_computed_and_finite() -> None:
    """The I11-rel / I11-abs evaluation runs on a sub-grid and is finite. NOT a gate: the 2 % row is
    not blind; a disclosed gate (plan row I11-rel). Prints max |rel| and |abs| per ion."""
    sub = {name: (MH.I11_GRID_MEV_U[name][1],) for name in IONS}
    rows = X.i11_tables(fractions=(0.1, 0.9), grid=sub, n=60)
    assert len(rows) == 6
    for r in rows:
        for key in ("rel", "abs", "map_residual", "pv_eq_over_pv"):
            assert math.isfinite(float(r[key])), (r["ion"], key)
        assert float(r["map_residual"]) <= X.RESIDUAL_TOL
    for name in IONS:
        sel = [r for r in rows if r["ion"] == name]
        print(
            f"{name}: max|rel-1| {max(abs(float(r['rel'])) for r in sel):.4f}, "
            f"max|abs| {max(abs(float(r['abs'])) for r in sel):.4f} (sub-grid, n=60)"
        )
