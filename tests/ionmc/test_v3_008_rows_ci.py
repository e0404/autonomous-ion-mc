"""V3-008 CI rows of ``validation/plans/v3-008-acceptance.md`` section A for the ion engine
(12C, 16O and 4He as transported projectiles; decision 0043): I0, I1, I4, I4b, I5, I6, I7, I10,
I11-map (transport level), I16, C-res and T1-ion (CI part).

Each test asserts the criterion of its row; no tolerance is changed. Deterministic or very short
runs (seconds). Row text is quoted in the test docstrings where it matters.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Any

import numpy as np
import pytest
import warp as wp

from ionmc._wpfunc import python_twin
from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    validate,
)
from ionmc.data.nist_star import parse_star_text
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.materials import MATERIALS, MUSCLE_SKELETAL_ICRP, WATER, Material, water
from ionmc.physics import exact_omega as XO
from ionmc.physics.em import make_em
from ionmc.physics.ion_mcs import THIN_SLAB_R, make_ion_mcs
from ionmc.physics.kinematics import make_kinematics
from ionmc.physics.projectiles import (
    ALPHA,
    CARBON12,
    DEUTERON,
    HELIUM3,
    OXYGEN16,
    PROTON,
    TRITON,
    Projectile,
)
from ionmc.physics.scattering import inverse_scattering_length_cm2_per_g
from ionmc.physics.stopping import (
    BetheOptions,
    BetheStoppingSource,
    NistStarStoppingSource,
    bethe_mass_stopping,
    build_table,
)
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport import ion_mcs_checks
from ionmc.transport.mcs_checks import radiation_length_g_cm2
from ionmc.transport.parity import compare_traces
from ionmc.transport.tables import TransportTables

wp.config.log_level = wp.LOG_WARNING
D = wp.float64
IONS = (ALPHA, CARBON12, OXYGEN16)
ION_IDS = [p.name for p in IONS]
STANDARD_E_U = {"alpha": 150.0, "carbon12": 290.0, "oxygen16": 350.0}
CI_E_U = {"alpha": 100.0, "carbon12": 150.0, "oxygen16": 150.0}
BOX_MM = {"alpha": 100.0, "carbon12": 100.0, "oxygen16": 100.0}  # >= range at CI_E_U (78, 53, 50)
SIDE = 60.0
INV_XS_WATER = inverse_scattering_length_cm2_per_g(WATER)
"""1/X_S of water [cm2/g], the value the transport loads into ``TransportTables.inv_rho_xs``."""
RHO_WATER = WATER.density_g_cm3
"""Density of water [g/cm3] (1.0), the voxel density of the water phantoms."""


# --------------------------------------------------------------------------------------------
# helpers
def _slab(dz: float, nz: int, material: Material = WATER) -> VoxelGeometry:
    return VoxelGeometry(
        origin_mm=(-SIDE / 2, -SIDE / 2, 0.0),
        spacing_mm=(SIDE, SIDE, dz),
        shape=(1, 1, nz),
        materials=(material,),
        material_index=np.zeros((1, 1, nz), dtype=np.int32),
    )


def _cfg(
    proj: Projectile,
    e_u: float | None = None,
    *,
    backend: str = "python",
    precision: str = "float64",
    n: int = 4,
    nb: int = 2,
    seed: int = 7,
    mcs: bool = True,
    straggling: bool = True,
    max_step: float = 1.0,
    stopping: Any = None,
    geometry: Any = None,
    scoring: tuple[ScoringGrid, ...] | None = None,
    tallies: tuple[TallyRequest, ...] = (),
    diagnostics: DiagnosticsOptions | None = None,
    sigma_mev: float = 0.0,
    nuclear: bool = False,
    allow_invalid: bool = False,
    physics_kwargs: dict[str, Any] | None = None,
    energy_mev: float | None = None,
) -> SimulationConfig:
    e0 = energy_mev if energy_mev is not None else (e_u or CI_E_U[proj.name]) * proj.a
    geo = geometry or BoxPhantom(
        (-SIDE / 2, -SIDE / 2, 0.0), (SIDE, SIDE, BOX_MM.get(proj.name, 100.0)), WATER
    )
    grids = scoring or (ScoringGrid((-SIDE / 2, -SIDE / 2, 0.0), (2.0, 2.0, 2.0), (30, 30, 50)),)
    return SimulationConfig(
        source=PencilBeamSource(proj, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), e0, sigma_mev, 0.0),
        geometry=geo,
        scoring=grids,
        physics=PhysicsOptions(
            nuclear=nuclear,
            stopping=stopping or BetheStoppingSource(),
            straggling=straggling,
            multiple_scattering=mcs,
            max_step_mm=max_step,
            **(physics_kwargs or {}),
        ),
        run=RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision=precision,  # type: ignore[arg-type]
            seed=seed,
            n_histories=n,
            n_batches=nb,
            allow_invalid_result=allow_invalid,
        ),
        diagnostics=diagnostics or DiagnosticsOptions(),
        tallies=tallies,
    )


def _star(program: str) -> Any:
    proj = PROTON if program == "PSTAR" else ALPHA
    e_u = np.geomspace(1.0, 600.0, 120)
    s = bethe_mass_stopping(e_u, proj, water(75.0))
    text = [f"{program}: x", "WATER, LIQUID", "", "h", "h", "h", ""]
    r = np.cumsum(np.concatenate(([0.0], proj.a * np.diff(e_u) / s[1:])))
    for eu, si, ri in zip(e_u, s, r + 1e-3, strict=True):
        text.append(f"{eu * proj.a:.6E} {si:.6E} 1.0E-3 {si:.6E} {ri:.6E} {ri * 0.9:.6E} 1.0")
    return parse_star_text("\n".join(text))


# --------------------------------------------------------------------------------------------
# I0: capability contract
@pytest.mark.parametrize(
    ("backend", "precision"),
    [("python", "float64"), ("warp-cpu", "float64"), ("warp-cpu", "float32")],
)
@pytest.mark.parametrize("proj", IONS, ids=ION_IDS)
def test_i0_sources_run_and_persist_fields(proj: Projectile, backend: str, precision: str) -> None:
    """He/C/O sources run on python and warp-cpu (cuda: ``test_i0_cuda_source_runs``); the
    effective configuration persists projectile, z, mass, ``e_cut_ion_mev_per_u``, the MCS model
    id and the f_dM domain status; a run is valid with every counter 0 and a closed balance."""
    cfg = _cfg(proj, backend=backend, precision=precision, n=4, nb=2)
    sim = Simulation(cfg)
    res = sim.run()
    assert res.valid and not res.counters.any_nonzero
    assert res.counters.as_dict()["mcs_mapping_failure"] == 0
    s = sim.effective.summary()
    assert s["source"]["projectile"] == proj.name
    assert (s["source"]["z"], s["source"]["a"], s["source"]["mass_mev"]) == (
        proj.z,
        proj.a,
        proj.mass_mev,
    )
    assert s["physics"]["e_cut_ion_mev_per_u"] == 2.0
    assert s["physics"]["e_cut_effective_mev"] == 2.0 * proj.a
    assert s["physics"]["mcs_model"] == "differential_moliere"
    assert (s["physics"]["ls_requested"], s["physics"]["ls_effective"]) == ("off", "off")
    m = s["mcs_mapping"]
    assert m["mapping"].startswith("exact-Omega")
    assert m["p1v1_eq_kinetic_energy_mev"] > 0.0
    assert m["fdm_extrapolated"] == (m["p1v1_eq_kinetic_energy_mev"] > m["fdm_domain_edge_mev"])
    assert sim.effective.requested.physics.e_cut_ion_mev_per_u == 2.0
    # tables: the ion's own, on the total kinetic energy axis
    assert sim.effective.tables.projectile == proj


@pytest.mark.cuda
@pytest.mark.parametrize("proj", IONS, ids=ION_IDS)
def test_i0_cuda_source_runs(proj: Projectile) -> None:
    res = Simulation(_cfg(proj, backend="warp-cuda", precision="float32", n=64, nb=2)).run()
    assert res.valid


@pytest.mark.parametrize("bad", [HELIUM3, DEUTERON, TRITON], ids=lambda p: p.name)
def test_i0_other_ions_fail_closed(bad: Projectile) -> None:
    with pytest.raises(UnsupportedCombinationError, match="not supported"):
        validate(_cfg(bad, energy_mev=bad.a * 100.0))
    forged = Projectile("carbon12", "C12", 6, 12, CARBON12.mass_mev * 1.001)
    with pytest.raises(UnsupportedCombinationError):
        validate(_cfg(forged, energy_mev=12 * 150.0))


@pytest.mark.parametrize(
    ("proj", "e_u", "sigma_u"),
    [
        (ALPHA, 49.9, 0.0), (ALPHA, 250.1, 0.0), (ALPHA, 245.0, 1.0),  # E0 + 6 sigma_E = 251
        (CARBON12, 79.0, 0.0), (CARBON12, 431.0, 0.0), (CARBON12, 425.0, 1.5),
        (OXYGEN16, 99.0, 0.0), (OXYGEN16, 431.0, 0.0), (OXYGEN16, 424.0, 1.2),
    ],
)  # fmt: skip
def test_i0_energy_outside_domain_fails_closed(
    proj: Projectile, e_u: float, sigma_u: float
) -> None:
    cfg = _cfg(proj, energy_mev=e_u * proj.a, sigma_mev=sigma_u * proj.a)
    with pytest.raises(UnsupportedCombinationError, match="supported domain"):
        validate(cfg)


@pytest.mark.parametrize(
    ("proj", "e_u"), [(ALPHA, 50.0), (ALPHA, 250.0), (CARBON12, 430.0), (OXYGEN16, 100.0)]
)
def test_i0_domain_edges_accepted(proj: Projectile, e_u: float) -> None:
    geo = BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, 500.0), WATER)
    eff = validate(_cfg(proj, e_u, geometry=geo))
    assert eff.fdm_domain.p1v1_eq_kinetic_energy_hi_mev is not None


def test_i0_nuclear_with_ions_fails_closed() -> None:
    for proj in IONS:
        with pytest.raises(UnsupportedCombinationError, match="nuclear"):
            validate(_cfg(proj, nuclear=True, physics_kwargs={"nuclear_table_id": "x" * 64}))


def test_i0_nist_star_for_c_and_o_and_astar_non_water_fail_closed() -> None:
    pstar, astar = _star("PSTAR"), _star("ASTAR")
    for proj in (CARBON12, OXYGEN16):
        with pytest.raises(UnsupportedCombinationError, match="no stopping table"):
            validate(_cfg(proj, stopping=NistStarStoppingSource(astar)))
    bone = MATERIALS["bone_compact_icru"]
    geo = BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, 100.0), bone)
    with pytest.raises(UnsupportedCombinationError, match="no stopping table"):
        validate(_cfg(ALPHA, stopping=NistStarStoppingSource(astar), geometry=geo))
    with pytest.raises(UnsupportedCombinationError):  # a PSTAR table cannot describe He either
        validate(_cfg(ALPHA, stopping=NistStarStoppingSource(pstar)))
    # control: ASTAR with water is accepted for helium
    eff = validate(_cfg(ALPHA, stopping=NistStarStoppingSource(astar, e_max_per_u=250.0)))
    assert eff.requested.physics.stopping.name == "nist-star"


def test_i0_let_medium_and_unproducible_tallies_fail_closed() -> None:
    def tally(**kw: Any) -> tuple[TallyRequest, ...]:
        return (TallyRequest("t", "dose", "let_t", **kw),)

    ok = validate(_cfg(CARBON12, tallies=tally(species=("carbon12",), generation="primary")))
    assert ok.channels is not None
    cases = [
        (tally(let_medium="local"), "let_medium"),
        (tally(species=("proton",)), "not producible"),
        (tally(species=("alpha",)), "not producible"),
        (tally(generation="secondary"), "not producible"),
        (tally(species=("oxygen16",)), "not producible"),
    ]
    for t, msg in cases:
        with pytest.raises(UnsupportedCombinationError, match=msg):
            validate(_cfg(CARBON12, tallies=t))
    with pytest.raises(UnsupportedCombinationError, match="not producible"):
        validate(_cfg(ALPHA, tallies=tally(species=("carbon12",))))


def test_i0_e_cut_rules() -> None:
    """``e_cut_mev`` is a proton option: an ion run refuses a non-default value; the ion cutoff
    must be at least twice the table floor (C-res)."""
    with pytest.raises(UnsupportedCombinationError, match="e_cut_ion_mev_per_u"):
        validate(_cfg(CARBON12, physics_kwargs={"e_cut_mev": 3.0}))
    with pytest.raises(UnsupportedCombinationError, match="twice the table floor"):
        validate(_cfg(CARBON12, physics_kwargs={"e_cut_ion_mev_per_u": 1.5}))
    eff = validate(_cfg(CARBON12, physics_kwargs={"e_cut_ion_mev_per_u": 3.0}))
    assert eff.e_cut_mev == 36.0 and eff.summary()["physics"]["e_cut_ion_mev_per_u"] == 3.0


@pytest.mark.parametrize("evaluation", ["table", "root"])
@pytest.mark.parametrize("backend", ["python", "warp-cpu"])
def test_i0_mapping_failure_is_a_named_counter_and_invalidates(
    backend: str, evaluation: str
) -> None:
    """A failed mapping fails closed: ``valid=False`` with the named counter
    ``mcs_mapping_failure``; the step scatters with zero variance (an impossible proton path,
    status 3, for both evaluations: the table-accelerated root solves the same equation). A table
    without any valid node is not a failure: every step has no starting guess and uses the full
    bracket from the start, as ``"root"``; no tight bracket is tried, so the diagnostic
    ``mcs_table_bracket_miss`` stays 0 (Codex REVIEW-16c63514: a step without a guess is not a
    miss; formerly counted), and the result stays valid."""
    kw = {"mcs_mapping_eval": evaluation}
    sim = Simulation(
        _cfg(CARBON12, backend=backend, n=2, nb=2, allow_invalid=True, physics_kwargs=kw)
    )
    im = sim.effective.ion_mapping
    assert im is not None and im.evaluation == evaluation
    pro = im.proton.copy()
    pro[:, 7] = 1.5  # the proton path ends below its own start: status 3 (path too short)
    sim.effective = replace(sim.effective, ion_mapping=replace(im, proton=pro))
    res = sim.run()
    assert not res.valid
    assert res.counters.mcs_mapping_failure > 0  # type: ignore[attr-defined]
    good = Simulation(_cfg(CARBON12, backend=backend, n=2, nb=2, physics_kwargs=kw)).run()
    assert good.valid and good.counters.as_dict()["mcs_mapping_failure"] == 0
    assert good.counters.as_dict()["mcs_table_bracket_miss"] == 0
    if evaluation == "table":
        assert im.table is not None
        tab = im.table.copy()
        tab[:, XO.TABLE_HEADER :] = XO.TABLE_INVALID_NODE  # no starting guess anywhere
        sim2 = Simulation(_cfg(CARBON12, backend=backend, n=2, nb=2, physics_kwargs=kw))
        sim2.effective = replace(sim2.effective, ion_mapping=replace(im, table=tab))
        miss = sim2.run()
        assert miss.valid  # a diagnostic, never a failure
        assert miss.counters.as_dict()["mcs_table_bracket_miss"] == 0  # no guess, no miss
        assert miss.counters.as_dict()["mcs_mapping_failure"] == 0


def test_i0_proton_summary_and_counters_unchanged() -> None:
    """z = 1: no ion block in the counters, no mapping, the proton path of the kernel unchanged."""
    cfg = _cfg(PROTON, energy_mev=20.0, geometry=BoxPhantom((-30, -30, 0), (60, 60, 20), WATER))
    sim = Simulation(cfg)
    assert sim.effective.ion_mapping is None and not sim.effective.is_ion
    res = sim.run()
    assert "mcs_mapping_failure" not in res.counters.as_dict()
    assert sim.effective.summary()["mcs_mapping"]["fdm_extrapolated"] is False


# --------------------------------------------------------------------------------------------
# I1: shared functions with a charge argument and the mapped MCS step
CHARGES = ((PROTON, 1.0), (ALPHA, 2.0), (CARBON12, 6.0), (OXYGEN16, 8.0))
EM_NAMES = ("bohr_variance", "scattering_power_dm", "scattering_variance_birth")


def _charge_kernel(real: Any) -> Any:
    EM = make_em(real)
    K = make_kinematics(real)
    inv_xs = wp.constant(real(INV_XS_WATER))
    rho_w = wp.constant(real(RHO_WATER))

    @wp.kernel(module="unique")
    def k(
        e_mid: wp.array(dtype=real),
        e_end: wp.array(dtype=real),
        e_0: wp.array(dtype=real),
        m: wp.array(dtype=real),
        z: wp.array(dtype=real),
        out: wp.array2d(dtype=real),
    ):  # type: ignore[no-untyped-def]
        i = wp.tid()
        pv = K.pv_mev(e_mid[i], m[i])
        pv_e = K.pv_mev(e_end[i], m[i])
        p1v1 = K.pv_mev(e_0[i], m[i])
        out[i, 0] = EM.bohr_variance(e_mid[i], m[i], z[i], real(0.55), real(1.0), real(1.0))
        out[i, 1] = EM.scattering_power_dm(pv, p1v1, z[i], inv_xs, rho_w)
        out[i, 2] = EM.scattering_variance_birth(pv, pv_e, p1v1, z[i], inv_xs, rho_w, real(0.7))
        out[i, 3] = K.beta2(e_mid[i], m[i])
        out[i, 4] = K.tmax_mev(e_mid[i], m[i])

    return k


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_i1_charge_argument_functions_twin_equals_warp_cpu(precision: str) -> None:
    """Shared functions with a charge argument, z in {1, 2, 6, 8}, at velocity-scaled arguments
    (energies per nucleon 2 to 430 MeV/u times A): python twin against a Warp CPU kernel. float64:
    decision 0001 class (rtol 1e-12, atol 1e-14); float32 against the float64 twin on float32
    rounded inputs, relative 5e-4 (the 2 MeV/u beta^2 carries the 1.2e-7 / (gamma^2 - 1) = 3e-5
    cancellation of float32 gamma^2 - 1, analytic)."""
    real = wp.float64 if precision == "float64" else wp.float32
    dt = np.float64 if precision == "float64" else np.float32
    rtol, atol = (1e-12, 1e-14) if precision == "float64" else (5e-4, 1e-9)
    rows = []
    for proj, z in CHARGES:
        for e_u in (2.0, 20.0, 150.0, 430.0):
            e = e_u * proj.a
            rows.append((0.98 * e, 0.9 * e, e, proj.mass_mev, z))
    arr = np.array(rows, dtype=dt)
    n = len(rows)
    out = wp.zeros((n, 5), dtype=real, device="cpu")
    wp.launch(
        _charge_kernel(real), dim=n,
        inputs=[wp.array(arr[:, j].copy(), dtype=real, device="cpu") for j in range(5)] + [out],
        device="cpu",
    )  # fmt: skip
    got = out.numpy().astype(np.float64)
    em, kin = python_twin(make_em), python_twin(make_kinematics)
    for i in range(n):
        e_mid, e_end, e0, m, z = (float(x) for x in arr[i])
        pv, pv_e, p1 = kin.pv_mev(e_mid, m), kin.pv_mev(e_end, m), kin.pv_mev(e0, m)
        ref = (
            em.bohr_variance(e_mid, m, z, 0.55, 1.0, 1.0),
            em.scattering_power_dm(pv, p1, z, float(dt(INV_XS_WATER)), float(dt(RHO_WATER))),
            em.scattering_variance_birth(
                pv, pv_e, p1, z, float(dt(INV_XS_WATER)), float(dt(RHO_WATER)), 0.7
            ),
            kin.beta2(e_mid, m),
            kin.tmax_mev(e_mid, m),
        )
        for j, r in enumerate(ref):
            assert abs(got[i, j] - r) <= atol + rtol * abs(r), (precision, i, j, got[i, j], r)
    # the charge enters as z^2: ratio of the z = 6 and z = 1 variances at equal arguments
    e = 150.0
    v1 = em.bohr_variance(e, CARBON12.mass_mev, 1.0, 0.55, 1.0, 1.0)
    v6 = em.bohr_variance(e, CARBON12.mass_mev, 6.0, 0.55, 1.0, 1.0)
    assert v6 / v1 == pytest.approx(36.0, rel=1e-14)


def _ion_kernel(real: Any) -> Any:
    IM = make_ion_mcs(real)
    inv_xs = wp.constant(real(INV_XS_WATER))
    rho_w = wp.constant(real(RHO_WATER))

    @wp.kernel(module="unique")
    def k(
        ion: wp.array2d(dtype=real),
        pro: wp.array2d(dtype=real),
        comp: wp.array2d(dtype=real),
        quad: wp.array(dtype=real),
        z: real,
        mass: real,
        pv: wp.array(dtype=real),
        pv_end: wp.array(dtype=real),
        p1: wp.array(dtype=real),
        birth: wp.array(dtype=int),
        out: wp.array2d(dtype=real),
        oi: wp.array(dtype=int),
    ):  # type: ignore[no-untyped-def]
        i = wp.tid()
        v, st = IM.mcs_variance(
            ion, pro, comp, quad, 0, z, mass, pv[i], pv_end[i], p1[i], birth[i],
            inv_xs, rho_w, real(0.7),
        )  # fmt: skip
        peq, p1eq, st2 = IM.mapped_state(ion, pro, comp, quad, 0, z, mass, pv[i], p1[i])
        out[i, 0] = v
        out[i, 1] = peq
        out[i, 2] = p1eq
        oi[i] = st * 10 + st2

    return k


def _ion_states(mp: XO.ExactOmegaMapping) -> list[tuple[float, float, float, int]]:
    out = []
    for fr in (0.01, 0.1, 0.25, 0.5, 0.75, 0.9):
        pv, p1 = mp.state_at_depth(150.0 * mp.projectile.a, fr)
        if p1 / pv > 1.05:
            out.append((pv, pv * 0.99999, p1, 0))
            out.append((pv, pv * 0.99, p1, 1))  # a birth-step state (analytic memory average)
    return out


def test_i1_new_callables_are_completely_covered() -> None:
    """Completeness: every callable of the new factory ``make_ion_mcs`` is exercised below (root
    path; ``variance_of_state`` through both) or in ``test_ion_mcs_table.py`` (table path)."""
    exercised = {  # the table path in tests/ionmc/test_ion_mcs_table.py
        "mapped_state", "mcs_variance", "mapped_state_table", "mcs_variance_table",
        "variance_of_state", "interpolated_state",
    }  # fmt: skip

    def names(ns: Any) -> set[str]:
        return {k for k, v in vars(ns).items() if callable(v)}

    assert names(python_twin(make_ion_mcs)) == exercised
    assert names(make_ion_mcs(wp.float64)) == exercised


@pytest.mark.parametrize("precision", ["float64", "float32"])
@pytest.mark.parametrize("proj", [ALPHA, CARBON12, OXYGEN16], ids=ION_IDS)
def test_i1_mapped_mcs_step_twin_equals_warp_cpu(proj: Projectile, precision: str) -> None:
    real = wp.float64 if precision == "float64" else wp.float32
    dt = np.float64 if precision == "float64" else np.float32
    rtol, atol = (1e-12, 1e-14) if precision == "float64" else (5e-4, 1e-12)
    mp = XO.ExactOmegaMapping.build(proj, WATER)
    st = _ion_states(mp)
    n = len(st)
    cols = [np.array([s[j] for s in st], dtype=dt) for j in range(3)]
    birth = np.array([s[3] for s in st], dtype=np.int32)

    def A(a: np.ndarray) -> Any:
        return wp.array(np.atleast_2d(a).astype(dt), dtype=real, device="cpu")

    out = wp.zeros((n, 3), dtype=real, device="cpu")
    oi = wp.zeros(n, dtype=int, device="cpu")
    wp.launch(
        _ion_kernel(real), dim=n,
        inputs=[A(mp.ion_path.pack), A(mp.proton.pack), A(mp.comp),
                wp.array(mp.quad.astype(dt), dtype=real, device="cpu"),
                real(proj.z), real(proj.mass_mev),
                *[wp.array(c, dtype=real, device="cpu") for c in cols],
                wp.array(birth, dtype=int, device="cpu"), out, oi],
        device="cpu",
    )  # fmt: skip
    got, status = out.numpy().astype(np.float64), oi.numpy()
    tw = python_twin(make_ion_mcs)

    def L(a: np.ndarray) -> list[list[float]]:
        return [a.astype(dt).astype(np.float64).tolist()]

    quad = mp.quad.astype(dt).astype(np.float64).tolist()
    for i in range(n):
        pv, pe, p1 = (float(c[i]) for c in cols)
        v, s1 = tw.mcs_variance(
            L(mp.ion_path.pack), L(mp.proton.pack), L(mp.comp), quad, 0, float(proj.z),
            float(dt(proj.mass_mev)), pv, pe, p1, int(birth[i]), float(dt(INV_XS_WATER)),
            float(dt(RHO_WATER)),
            float(dt(0.7)),
        )  # fmt: skip
        peq, p1eq, s2 = tw.mapped_state(
            L(mp.ion_path.pack), L(mp.proton.pack), L(mp.comp), quad, 0, float(proj.z),
            float(dt(proj.mass_mev)), pv, p1,
        )  # fmt: skip
        assert s1 == s2 == 0 and status[i] == 0
        assert v > 0.0
        for a, b in ((v, got[i, 0]), (peq, got[i, 1]), (p1eq, got[i, 2])):
            assert abs(a - b) <= atol + rtol * abs(b), (precision, proj.name, i, a, b)


def test_i1_thin_slab_guard_and_z1_flag() -> None:
    """A state at the initial momentum (pv = p1v1) is mapped at p1v1 / (1 + THIN_SLAB_R) instead of
    failing; the host flag ``ion_is_proton`` of ``exact_omega`` keeps the identity for z = 1."""
    mp = XO.ExactOmegaMapping.build(CARBON12, WATER)
    pv, p1 = mp.state_at_depth(150.0 * 12, 0.5)
    tw = python_twin(make_ion_mcs)
    args = (
        [mp.ion_path.pack.tolist()],
        [mp.proton.pack.tolist()],
        [mp.comp.tolist()],
        mp.quad.tolist(),
    )
    peq, p1eq, st = tw.mapped_state(*args, 0, 6.0, CARBON12.mass_mev, p1, p1)
    assert st == 0 and p1eq / peq == pytest.approx(1.0 + THIN_SLAB_R, rel=1e-12)
    v, st_v = tw.mcs_variance(
        *args, 0, 6.0, CARBON12.mass_mev, p1, p1 * 0.999, p1, 1, 0.01, 1.0, 1.0
    )
    assert st_v == 0 and v > 0.0
    pmap = XO.ExactOmegaMapping.build(PROTON, WATER)
    assert pmap.fdm_equivalent(pv, p1)[1] == pv  # z = 1: pv_eq = pv exactly


def test_t1_ion_trajectory_parity_python_vs_warp_cpu_f64() -> None:
    """T1-ion (CI part): K = 16 He 100 MeV/u histories, all EM on, python against warp-cpu float64:
    discrete trace columns identical for 100 % of steps, continuous within rtol = atol = 1e-10."""
    diag = DiagnosticsOptions(trace_histories=16, track_end_positions=True)
    a = Simulation(_cfg(ALPHA, n=16, nb=2, diagnostics=diag)).run()
    b = Simulation(_cfg(ALPHA, backend="warp-cpu", n=16, nb=2, diagnostics=diag)).run()
    v = compare_traces(a.diagnostics, b.diagnostics)
    assert v["pass"], v
    assert all(c == 0 for c in v["discrete_mismatches"].values())


# --------------------------------------------------------------------------------------------
# I4 / I4b: theory identities (all corrections off)
OFF = BetheOptions(shell=False, barkas=False, bloch=False, mott=False, effective_charge=False)
M_P = PROTON.mass_mev


def _surrogate(proj: Projectile) -> Projectile:
    """z = 1 particle with the ion's mass and mass number: the proton's formula with the Tmax mass
    term evaluated exactly for the ion."""
    return Projectile("surrogate", "X", 1, proj.a, proj.mass_mev)


@pytest.mark.parametrize("material", [WATER, MUSCLE_SKELETAL_ICRP], ids=lambda m: m.name)
@pytest.mark.parametrize("proj", IONS, ids=ION_IDS)
def test_i4_velocity_scaling_identity(proj: Projectile, material: Material) -> None:
    """S_ion(beta) = z^2 S_p(beta) and R_ion(beta) = (M_ion / m_p) / z^2 R_p(beta) at every grid
    node 2 to 500 MeV/u (corrections off). Exact version (Tmax mass term of the ion): rel <= 1e-6;
    against the real proton at the equal-velocity energy (Tmax with m_p): <= 2e-3 (reported)."""
    z = float(proj.z)
    src = BetheStoppingSource(options=OFF)
    tab = src.table(material, proj)
    sel = tab.energy_per_u >= 2.0
    e_u = tab.energy_per_u[sel]
    sur = _surrogate(proj)
    s_sur = bethe_mass_stopping(e_u, sur, material, **vars_of(OFF))
    assert np.max(np.abs(tab.s_el_mass[sel] / (z * z * s_sur) - 1.0)) <= 1e-6
    tab_sur = src.table(material, sur)
    r_ratio = tab.csda_range_g_cm2[sel] * z * z / tab_sur.csda_range_g_cm2[sel]
    assert np.max(np.abs(r_ratio - 1.0)) <= 1e-6
    # real proton at equal velocity: T_p = T_ion m_p / M_ion
    e_p = e_u * proj.a * M_P / proj.mass_mev
    s_p = bethe_mass_stopping(e_p, PROTON, material, **vars_of(OFF))
    dev_s = np.max(np.abs(tab.s_el_mass[sel] / (z * z * s_p) - 1.0))
    # proton table on the equal-velocity grid (start range from the same construction)
    s_p_all = bethe_mass_stopping(tab.energy_per_u * proj.a * M_P / proj.mass_mev, PROTON, material,
                                  **vars_of(OFF))  # fmt: skip
    r0 = tab.energy_per_u[0] * proj.a * M_P / proj.mass_mev / s_p_all[0]
    tab_p = build_table(
        PROTON, material, tab.energy_per_u * proj.a * M_P / proj.mass_mev, s_p_all, r0, {}
    )
    dev_r = np.max(
        np.abs(
            tab.csda_range_g_cm2[sel]
            / ((proj.mass_mev / M_P) / (z * z) * tab_p.csda_range_g_cm2[sel])
            - 1.0
        )
    )
    print(f"I4 {proj.name} {material.name}: S vs proton {dev_s:.2e}, R vs proton {dev_r:.2e}")
    assert dev_s <= 2e-3 and dev_r <= 2e-3


def vars_of(o: BetheOptions) -> dict[str, Any]:
    return {
        "shell": o.shell, "barkas": o.barkas, "bloch": o.bloch, "mott": o.mott,
        "effective_charge": o.effective_charge,
    }  # fmt: skip


def _sigma_r_over_r(
    proj: Projectile, e0_total: float, e_lo_total: float, material: Material
) -> float:
    """sigma_R / R from numpy quadrature sigma_R^2 = int (d sigma^2 / dt) / S^3 dE with the
    implemented Bohr variance (shared function, python twin), corrections off, total energy axis."""
    em = python_twin(make_em)
    n = 4001
    u = np.linspace(math.log(e_lo_total), math.log(e0_total), n)
    e = np.exp(u)
    s = bethe_mass_stopping(e / proj.a, proj, material, **vars_of(OFF))
    var_t = np.array(
        [
            em.bohr_variance(float(x), proj.mass_mev, float(proj.z), material.z_over_a, 1.0, 10.0)
            for x in e
        ]
    )  # MeV^2 per g/cm2 (rho = 1 g/cm3, 10 mm = 1 cm)
    h = u[1] - u[0]
    w = np.full(n, 2.0)
    w[1::2] = 4.0
    w[0] = w[-1] = 1.0

    def simpson(f: np.ndarray) -> float:
        return float(h / 3.0 * np.sum(w * f))

    rng = simpson(e / s)
    sig2 = simpson(e * var_t / s**3)
    return math.sqrt(sig2) / rng


@pytest.mark.parametrize("proj", IONS, ids=ION_IDS)
def test_i4b_range_straggling_scaling(proj: Projectile) -> None:
    """sigma_R/R(ion) / [sigma_R/R(p) (m_p/M)^(1/2)] - 1 within 1e-3 at the standard point (He 150,
    C 290, O 350 MeV/u; proton at the equal-velocity energy, lower limit 2 MeV/u scaled)."""
    e_u = STANDARD_E_U[proj.name]
    scale = M_P / proj.mass_mev
    ion = _sigma_r_over_r(proj, e_u * proj.a, 2.0 * proj.a, WATER)
    pro = _sigma_r_over_r(PROTON, e_u * proj.a * scale, 2.0 * proj.a * scale, WATER)
    dev = ion / (pro * math.sqrt(scale)) - 1.0
    print(f"I4b {proj.name}: {dev:.2e}")
    assert abs(dev) <= 1e-3


# --------------------------------------------------------------------------------------------
# I5 / I6 / C-res / I7 / I16: transport-level deterministic rows
def _tables(proj: Projectile, material: Material = WATER) -> TransportTables:
    return TransportTables.from_stopping_tables(
        [BetheStoppingSource().table(material, proj)], projectile=proj
    )


@pytest.mark.parametrize("proj", IONS, ids=ION_IDS)
def test_i5_range_table_consistency_for_ions(proj: Projectile) -> None:
    """V3-003D rows D1b and D1 re-pointed to the ion tables: strictly monotone R(E) and Rinv(r),
    round trips within 1e-5, the nodes equal the StoppingTable nodes (same grid) and the exact
    construction identity is recorded."""
    table = BetheStoppingSource().table(WATER, proj)
    t = TransportTables.from_stopping_tables([table], projectile=proj)
    energies = np.geomspace(t.e_min_mev[0], t.e_max_mev[0], 6001)
    ranges = np.array([t.range_g_cm2(0, e) for e in energies])
    assert np.all(np.diff(ranges) > 0.0)
    back = np.array([t.energy_from_range(0, r) for r in ranges])
    assert np.all(np.diff(back) > 0.0)
    assert np.max(np.abs(back / energies - 1.0)) <= 1e-5
    r_grid = np.geomspace(t.r_min_g_cm2[0], t.r_max_g_cm2[0], 6001)
    e_of_r = np.array([t.energy_from_range(0, r) for r in r_grid])
    assert np.all(np.diff(e_of_r) > 0.0)
    assert np.max(np.abs(np.array([t.range_g_cm2(0, e) for e in e_of_r]) / r_grid - 1.0)) <= 1e-5
    assert np.max(np.abs(table.csda_range_g_cm2 / t.r_mass[0] - 1.0)) <= 1e-12
    assert t.identity[0]["range_construction"] == "exact-loglog-quadrature-v1"
    assert t.identity[0]["energy_axis"] == "total_kinetic_mev"
    assert t.e_min_mev[0] == pytest.approx(proj.a * 1.0)


def _det_run(proj: Projectile, backend: str, s_max: float, e_u: float | None = None) -> Any:
    length = BOX_MM[proj.name]
    nz = int(round(length / s_max))
    cfg = _cfg(
        proj, e_u, backend=backend, n=2, nb=2, straggling=False, mcs=False, max_step=s_max,
        geometry=_slab(s_max, nz), scoring=(ScoringGrid((-30, -30, 0), (60, 60, 2.0), (1, 1, 50)),),
        diagnostics=DiagnosticsOptions(trace_histories=1, track_end_positions=True),
    )  # fmt: skip
    return Simulation(cfg).run()


@pytest.mark.parametrize(
    ("backend", "s_max"), [("python", 1.0), ("python", 0.25), ("warp-cpu", 1.0)]
)
@pytest.mark.parametrize("proj", [ALPHA, CARBON12], ids=lambda p: p.name)
def test_i5_branch_consistency_d2(proj: Projectile, backend: str, s_max: float) -> None:
    """D2 re-pointed: for every step boundary (z_b, E_b) with E_b > E_cut,
    |R(E_b) - (R(E0) - z_b)| <= 0.5 um (deterministic history, straggling and MCS off)."""
    res = _det_run(proj, backend, s_max)
    t = res.effective_config.tables
    e0 = res.effective_config.requested.source.kinetic_energy_mev
    ec = res.effective_config.e_cut_mev
    tr = res.diagnostics["trace"]
    sel = (tr["history"] == 0) & (tr["energy_mev"] > ec)
    z_b, e_b = tr["z_mm"][sel], tr["energy_mev"][sel]
    assert z_b.size >= 20
    r_b = np.array([t.range_g_cm2(0, e) for e in e_b]) * 10.0
    resid = np.abs(r_b - (t.range_g_cm2(0, e0) * 10.0 - z_b))
    assert float(np.max(resid)) * 1000.0 <= 0.5


@pytest.mark.parametrize("backend", ["python", "warp-cpu"])
def test_i6_deterministic_end_depth_he_100(backend: str) -> None:
    """I6 (CI: He 100 MeV/u only): R(E0) - R(E_cut) - 1e-4 R <= z_end <= R(E0) + 1e-4 R."""
    res = _det_run(ALPHA, backend, 1.0, 100.0)
    t = res.effective_config.tables
    e0 = res.effective_config.requested.source.kinetic_energy_mev
    ec = res.effective_config.e_cut_mev
    r0, rc = t.range_g_cm2(0, e0) * 10.0, t.range_g_cm2(0, ec) * 10.0
    z_end = float(res.diagnostics["end_position_mm"][0, 2])
    assert r0 - rc - 1e-4 * r0 <= z_end <= r0 + 1e-4 * r0
    assert res.diagnostics["end_code"][0] == 0  # END_CUTOFF


def test_i6_residual_range_at_e_cut_reported() -> None:
    """I6 (report): the residual CSDA range at E_cut per ion against the decision's estimates
    (He ~ 0.075, C ~ 0.025-0.03, O ~ 0.02-0.025 mm). The row only reports; the measured values
    (printed; He 0.092, C 0.043 mm) lie 20-70 % above the velocity-scaling estimates, a finding for
    the report. The assertion is a sanity band (factor 3), not a criterion."""
    est = {"alpha": 0.075, "carbon12": 0.0275, "oxygen16": 0.0225}
    for proj in IONS:
        rc = _tables(proj).range_g_cm2(0, 2.0 * proj.a) * 10.0
        print(f"I6 residual range at 2 MeV/u: {proj.name} {rc:.4f} mm (estimate {est[proj.name]})")
        assert est[proj.name] / 3.0 <= rc <= est[proj.name] * 3.0


@pytest.mark.parametrize(
    ("backend", "precision", "bound"),
    [("python", "float64", 1e-12), ("warp-cpu", "float64", 1e-12), ("warp-cpu", "float32", 1e-5)],
)
@pytest.mark.parametrize("proj", IONS, ids=ION_IDS)
def test_i7_energy_balance_ion_sources(
    proj: Projectile, backend: str, precision: str, bound: float
) -> None:
    """I7: energy balance, ion sources, all EM on: rel <= 1e-12 (f64), 1e-5 (f32); counters 0."""
    res = Simulation(
        _cfg(proj, backend=backend, precision=precision, n=8, nb=2, sigma_mev=0.2)
    ).run()
    eb = res.energy_balance
    assert eb.relative_residual <= bound
    for g in range(len(res.grids)):
        assert eb.grid_relative_residual(g) <= bound
    assert not res.counters.any_nonzero and res.valid
    assert eb.cutoff_mev > 0.0  # C-res: the residual energy below E_cut is tallied locally


def test_c_res_cutoff_energy_local_and_exact() -> None:
    """C-res: every history that stops deposits exactly E_cut-ish residual (<= 2 MeV/u x A) locally;
    the cutoff tally equals the in-grid part of the local deposits."""
    proj = CARBON12
    res = Simulation(
        _cfg(
            proj,
            straggling=False,
            mcs=False,
            n=4,
            nb=2,
            diagnostics=DiagnosticsOptions(track_end_positions=True),
        )
    ).run()
    eb = res.energy_balance
    ec = 2.0 * proj.a
    ends = res.diagnostics["end_energy_mev"]
    assert np.all(ends <= ec) and np.all(res.diagnostics["end_code"] == 0)
    assert eb.cutoff_mev == pytest.approx(float(ends.sum()), rel=1e-14)
    assert eb.initial_mev == pytest.approx(4 * 150.0 * 12.0)


@pytest.mark.parametrize("proj", IONS, ids=ION_IDS)
def test_i16_primary_ion_let_in_thin_slab(proj: Projectile) -> None:
    """I16: the track-averaged water LET of the primary ion in 1 mm slabs (deterministic, MCS and
    straggling off) equals S_w,ion(E_mid) of the ion's own water row within 1e-3."""
    grid = ScoringGrid((-30.0, -30.0, 0.0), (60.0, 60.0, 1.0), (1, 1, 10), name="dose")
    req = (TallyRequest("let", "dose", "let_t", species=(proj.name,), generation="primary"),)
    cfg = _cfg(
        proj, n=2, nb=2, straggling=False, mcs=False, max_step=0.5, scoring=(grid,), tallies=req
    )
    sim = Simulation(cfg)
    res = sim.run()
    q = res.grids[0].quantities["let"]
    t = sim.effective.tables
    e0 = cfg.source.kinetic_energy_mev
    r0 = t.range_g_cm2(0, e0)
    worst = 0.0
    n_ok = 0
    for k in range(10):
        if not q.defined_mask.reshape(-1)[k]:
            continue
        e_mid = t.energy_from_range(0, r0 - (k + 0.5) / 10.0)
        ref = t.s_water(e_mid, proj.name)
        worst = max(worst, abs(float(q.mean.reshape(-1)[k]) / ref - 1.0))
        n_ok += 1
    assert n_ok >= 8
    assert worst <= 1e-3
    assert sim.effective.summary()["tables"]["water_row"]["projectile"] == proj.name


# --------------------------------------------------------------------------------------------
# I10: MCS integrator step independence for T_ion (U5 analogue) and I11-map at transport level
# The CSDA path of the ion (I10, I11-map and the lv8 steps I9 / I12 / I15 through
# ``steps_v8.ion_path``) lives in ``ionmc.transport.ion_mcs_checks``: 1/X_S and rho are the
# material values the transport uses (``inverse_scattering_length_cm2_per_g`` and the density).
_mapping_arrays = ion_mcs_checks.mapping_arrays
_mcs_kernel = ion_mcs_checks.mcs_kernel
_IonPath = ion_mcs_checks.IonPath


def test_ion_path_uses_the_material_scattering_length() -> None:
    """The I9 / I10 / I12 / I15 reference path carries the transport's 1/X_S of water
    (``inverse_scattering_length_cm2_per_g(WATER)``, loaded into ``TransportTables.inv_rho_xs``),
    not the former placeholder 0.0123 cm2/g, and the water density; its T_ion scales as 1/X_S."""
    path = _IonPath(CARBON12, 290.0)
    assert path.inv_rho_xs_cm2_g == inverse_scattering_length_cm2_per_g(WATER)
    assert path.inv_rho_xs_cm2_g == pytest.approx(0.021331, rel=1e-4)
    assert abs(path.inv_rho_xs_cm2_g - 0.0123) > 1e-3
    assert path.inv_rho_xs_cm2_g == float(path.t.inv_rho_xs_cm2_g[0])
    assert path.rho_g_cm3 == WATER.density_g_cm3
    assert INV_XS_WATER == path.inv_rho_xs_cm2_g
    x = np.array([0.1, 0.5]) * path.r1_mm
    other = _IonPath(CARBON12, 290.0, inv_rho_xs_cm2_g=0.0123)
    np.testing.assert_allclose(
        path.power(x) / other.power(x), path.inv_rho_xs_cm2_g / 0.0123, rtol=1e-12
    )


STEPS_MM = (0.01, 0.1, 0.5, 1.0, 2.0, 5.0)
DEPTH_FRACTIONS = (0.05, 0.25, 0.5, 0.9)


@pytest.mark.parametrize("proj", IONS, ids=ION_IDS)
def test_i10_mcs_step_independence_and_highland_control(proj: Projectile) -> None:
    """I10: |stepped theta^2 / quadrature - 1| <= 2e-3 for s in {0.01, 0.1, 0.5, 1, 2, 5} mm at
    x/R in {0.05, 0.25, 0.5, 0.9}, standard point; negative control: the per-step Highland formula
    (with the charge) varies by >= 5 % over the same steps."""
    path = _IonPath(proj, STANDARD_E_U[proj.name])
    worst = 0.0
    for fr in DEPTH_FRACTIONS:
        x = fr * path.r1_mm
        exact = path.theta2_quadrature(x)
        assert exact > 0.0 and np.all(path.last_status == 0)
        for s in STEPS_MM:
            dev = abs(path.theta2_steps(x, s) / exact - 1.0)
            assert np.all(path.last_status == 0)
            worst = max(worst, dev)
            assert dev <= 2e-3, (proj.name, fr, s, dev)
    print(f"I10 {proj.name}: worst step dependence {worst:.2e}")
    # negative control
    x0 = radiation_length_g_cm2(WATER)
    x = 0.5 * path.r1_mm
    totals = []
    for s in STEPS_MM:
        n = int(round(x / s))
        step = x / n
        mids = (np.arange(n) + 0.5) * step
        pv = path._pv(path.energy(mids))
        step_g = step / 10.0
        th2 = (14.1 * proj.z / pv) ** 2 * (step_g / x0) * (1.0 + np.log10(step_g / x0) / 9.0) ** 2
        totals.append(float(np.sum(th2)))
    assert (max(totals) - min(totals)) / min(totals) >= 0.05


@pytest.mark.parametrize("proj", [ALPHA, CARBON12], ids=lambda p: p.name)
def test_i11_map_residual_at_transport_states(proj: Projectile) -> None:
    """I11-map at transport level: for the states of a traced transported history (python backend,
    the packed arrays of the effective configuration), the Warp CPU kernel path (``mapped_state``
    with the thin-slab guard) satisfies |b_p(pv_eq, r pv_eq) - b_ion(pv, p1v1)| <= 1e-10 with
    status 0 at every step."""
    cfg = _cfg(proj, n=2, nb=2, diagnostics=DiagnosticsOptions(trace_histories=2))
    sim = Simulation(cfg)
    res = sim.run()
    im = sim.effective.ion_mapping
    assert im is not None
    tr = res.diagnostics["trace"]
    e = tr["energy_mev"][(tr["history"] == 0)]
    e = e[e > 2.0 * proj.a * 1.5]
    assert e.size >= 20
    path = _IonPath(proj, CI_E_U[proj.name])
    path.arrs = (im.ion[:1], im.proton[:1], im.comp[:1], im.quad)
    pv = path._pv(e)
    out = path._run(pv, pv, np.ones_like(pv), 0)
    assert np.all(path.last_status == 0)
    assert float(np.max(np.abs(out[:, 1]))) <= 1e-10


@pytest.mark.parametrize("proj", [ALPHA, CARBON12], ids=lambda p: p.name)
def test_i11_map_table_against_root_at_transport_states(proj: Projectile) -> None:
    """I11-map, table counterpart (Codex REVIEW-ee53867c: the table only seeds the root): at the
    states of a traced transported history (python backend, the default ``"table"`` evaluation)
    the table-accelerated root of the effective configuration satisfies the same residual bound
    as the root path, |b_p(pv_eq, r pv_eq) - b_ion(pv, p1v1)| <= 1e-10, with status 0, agrees
    with the root path's pv_eq within 1e-8 and its step variance within 1e-8."""
    cfg = _cfg(proj, n=2, nb=2, diagnostics=DiagnosticsOptions(trace_histories=2))
    sim = Simulation(cfg)
    res = sim.run()
    im = sim.effective.ion_mapping
    assert im is not None and im.evaluation == "table" and im.table is not None
    tr = res.diagnostics["trace"]
    e = tr["energy_mev"][(tr["history"] == 0)]
    e = e[e > 2.0 * proj.a * 1.5]
    assert e.size >= 20
    path = _IonPath(proj, CI_E_U[proj.name])
    path.arrs = (im.ion[:1], im.proton[:1], im.comp[:1], im.quad)
    pv = path._pv(e)
    root = path._run(pv, pv, np.ones_like(pv), 0)
    assert np.all(path.last_status == 0)
    tw = python_twin(make_ion_mcs)
    xo = python_twin(XO.make_exact_omega)
    ion, pro, comp = ([a[0].tolist()] for a in (im.ion, im.proton, im.comp))
    quad = im.quad.tolist()
    tab = [im.table[0].tolist()]
    z, mass, p1 = float(proj.z), proj.mass_mev, float(path.p1v1)
    worst = worst_g = 0.0
    for k in range(pv.size):
        peq, _p1eq, st, _miss = tw.mapped_state_table(
            ion, pro, comp, quad, tab, 0, z, mass, float(pv[k]), p1
        )
        var, st_v, _m = tw.mcs_variance_table(
            ion, pro, comp, quad, tab, 0, z, mass, float(pv[k]), float(pv[k]), p1, 0,
            INV_XS_WATER, RHO_WATER, 1.0,
        )  # fmt: skip
        assert st == st_v == 0
        pv_map = min(float(pv[k]), p1 / (1.0 + THIN_SLAB_R))
        target, s_t = xo.b_slab(ion[0], comp[0], quad, z, mass, pv_map, p1)
        g, s_g = xo.b_residual(pro[0], comp[0], quad, math.log(peq), p1 / pv_map, target)
        assert s_t == s_g == 0
        worst_g = max(worst_g, abs(g))
        assert abs(peq / root[k, 2] - 1.0) <= 1e-8
        worst = max(worst, abs(var / root[k, 0] - 1.0))
    print(
        f"I11-map table {proj.name}: max residual {worst_g:.2e}, "
        f"max |var_table/var_root - 1| = {worst:.2e}"
    )
    assert worst_g <= 1e-10
    assert worst <= 1e-8


def test_i0_same_name_materials_are_evaluated_separately() -> None:
    """Distinct materials sharing a name keep separate ion tables in the f_dM domain scan (the
    tables are keyed by the full material identity, not by ``Material.name``) and the persisted
    mapping metadata is the maximum over both."""
    a = Material("dup", 1.0, {"H": 0.111894, "O": 0.888106})
    b = Material("dup", 1.19, {"H": 0.080538, "C": 0.599848, "O": 0.319614})
    stop = BetheStoppingSource(options=BetheOptions(ls="eff"))  # table differs from the default
    e_u = 290.0

    def status(mats: tuple[Material, ...]) -> Any:
        idx = np.arange(len(mats), dtype=np.int32).reshape(1, 1, len(mats))
        geo = VoxelGeometry(
            origin_mm=(-SIDE / 2, -SIDE / 2, 0.0),
            spacing_mm=(SIDE, SIDE, 50.0),
            shape=(1, 1, len(mats)),
            materials=mats,
            material_index=idx,
        )
        eff = validate(_cfg(CARBON12, e_u, geometry=geo, stopping=stop))
        return eff, eff.summary()["mcs_mapping"]

    _, sa = status((a,))
    _, sb = status((b,))
    eff, sab = status((a, b))
    assert sa["p1v1_eq_kinetic_energy_mev"] != sb["p1v1_eq_kinetic_energy_mev"]
    for key in ("p1v1_eq_kinetic_energy_mev", "p1v1_eq_kinetic_energy_hi_mev"):
        assert sab[key] == max(sa[key], sb[key])
    assert sab["fdm_extrapolated"] is True
    assert eff.ion_mapping is not None and eff.ion_mapping.ion.shape[0] == 2
    # the identity key separates what the name does not
    from ionmc.transport.tables import material_fingerprint  # noqa: PLC0415

    assert a.name == b.name and material_fingerprint(a) != material_fingerprint(b)
