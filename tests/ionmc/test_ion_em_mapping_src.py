"""V3-008 src exact-Omega mapping (``ionmc.physics.exact_omega``; rows I0, I1, I11-map).

* I11-map: b residual <= 1e-10 (float64 twin root) on a (pv, p1v1) set covering the ion domain
  (ion energy >= 2 MeV/u) for He/C/O in water and every ``materials.py`` material; src b and
  pv_eq against the committed reference script; z = 1 identity exact.
* I1: python twin against a Warp CPU kernel for EVERY callable of ``make_exact_omega`` (float64:
  decision 0001 class; float32 against the float64 twin on float32-rounded inputs, analytic
  budgets) and the completeness check.
* I0: proton table edge (fail closed), ``fdm_extrapolated`` and the proton-equivalent energy.
* LS metadata: requested/effective modes outside the table hash.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import warp as wp

from ionmc._wpfunc import python_twin
from ionmc.config import validate
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import MATERIALS, WATER
from ionmc.physics import exact_omega as XO
from ionmc.physics.em import make_em
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, PROTON
from ionmc.physics.stopping import BetheOptions, BetheStoppingSource
from ionmc.transport.tables import TransportTables

wp.config.log_level = wp.LOG_WARNING
ION_LIST = (ALPHA, CARBON12, OXYGEN16)
U = 2.0**-24
E_U = (20.0, 150.0, 430.0)
FRACTIONS = (0.001, 0.02, 0.1, 0.25, 0.5, 0.9, 0.97)
_MAPS: dict[tuple[str, str], XO.ExactOmegaMapping] = {}


def mapping(proj: Any, mat: Any) -> XO.ExactOmegaMapping:
    key = (proj.name, mat.name)
    if key not in _MAPS:
        _MAPS[key] = XO.ExactOmegaMapping.build(proj, mat)
    return _MAPS[key]


def domain_states(mp: XO.ExactOmegaMapping) -> list[tuple[float, float]]:
    """(pv, p1v1) at depth fractions of E0 in E_U, restricted to ion energy >= 2 MeV/u."""
    out = []
    a = mp.projectile.a
    for e_u in E_U:
        for fr in FRACTIONS:
            pv, p1 = mp.state_at_depth(e_u * a, fr)
            e = XO.proton_kinetic_from_pv  # noqa: F841 (kept local: proton only)
            m = mp.projectile.mass_mev
            q = 2 * m - pv
            t = 2 * pv * m / (q + math.sqrt(q * q + 4 * pv * m))
            if t / a >= 2.0 and pv < p1 * (1 - 1e-9):
                out.append((pv, p1))
    return out


@pytest.mark.parametrize("proj", ION_LIST, ids=lambda p: p.name)
@pytest.mark.parametrize("mat", list(MATERIALS.values()), ids=lambda m: m.name)
def test_i11_map_residual_all_materials(proj: Any, mat: Any) -> None:
    mp = mapping(proj, mat)
    states = domain_states(mp)
    assert len(states) >= 12
    for pv, p1 in states:
        peq, status, iters = mp.pv_equivalent(pv, p1)
        assert status == 0, (proj.name, mat.name, pv, p1, XO.STATUS_NAMES[status])
        assert iters <= XO.MAX_ITERATIONS
        assert abs(mp.residual(pv, p1, peq)) <= 1e-10
        assert 0.0 < peq < pv  # ions are the heavier-collision case: pv_eq < pv


def test_src_matches_reference_script() -> None:
    d = Path(__file__).resolve().parents[2] / "validation/scripts/ion_em"
    for nm in ("moliere_hanson", "exact_omega"):
        if nm not in sys.modules:
            s = importlib.util.spec_from_file_location(nm, d / f"{nm}.py")
            assert s is not None and s.loader is not None
            m = importlib.util.module_from_spec(s)
            sys.modules[nm] = m
            s.loader.exec_module(m)
    ref = sys.modules["exact_omega"]
    prot = ref.bethe_path("proton", 8000.0)
    for proj in ION_LIST:
        ion = ref.bethe_path(proj.name)
        mp = mapping(proj, WATER)
        for e_u, fr in ((150.0, 0.25), (430.0, 0.9)):
            pv, p1 = mp.state_at_depth(e_u * proj.a, fr)
            assert abs(mp.b_ion(pv, p1)[0] - ref.b_ion(ion, pv, p1)) <= 1e-4
            peq = mp.pv_equivalent(pv, p1)[0]
            assert peq == pytest.approx(ref.pv_eq(ion, prot, pv, p1), rel=1e-4)


def test_z1_identity_exact_and_root_consistent() -> None:
    mp = mapping(PROTON, WATER)
    assert mp.ion_is_proton == 1
    twin = python_twin(XO.make_exact_omega)
    em = python_twin(make_em)
    rng = np.random.default_rng(5)
    for _ in range(20):
        pv, p1 = mp.state_at_depth(float(rng.uniform(30.0, 480.0)), float(rng.uniform(0.05, 0.95)))
        peq, st, it = mp.pv_equivalent(pv, p1)
        assert (peq, st, it) == (pv, 0, 0)
        f, peq2, p1eq, st2 = mp.fdm_equivalent(pv, p1)
        assert peq2 == pv and p1eq == p1 and st2 == 0 and f == em.fdm(pv, p1)
        # without the flag the root solve of identical functions returns the same state
        ion, pro, comp, quad = mp._lists
        r, s, _ = twin.pv_equivalent(ion, pro, comp, quad, 1.0, PROTON.mass_mev, pv, p1, 0)
        assert s == 0 and r == pytest.approx(pv, rel=1e-9)


def test_packing_equals_transport_tables_row_and_default_hash_untouched() -> None:
    src = BetheStoppingSource()
    tab = src.table(WATER, PROTON)
    tt = TransportTables.from_stopping_tables([tab])
    assert np.array_equal(XO.pack_transport_row(tt), XO.pack_path(tab).pack)
    assert src.e_max_per_u == 500.0  # the default proton transport tables are unchanged
    p = XO.proton_path(WATER)
    assert p.e_max_mev == XO.PROTON_PATH_EMAX_MEV >= 3000.0 and p.e_min_mev == 1.0
    assert int(p.pack[4]) == len(
        BetheStoppingSource(e_max_per_u=XO.PROTON_PATH_EMAX_MEV).table(WATER, PROTON).energy_per_u
    )


def test_fail_closed_table_edge() -> None:
    with pytest.raises(UnsupportedCombinationError):
        XO.check_proton_table_edge(3001.0, 3000.0)
    XO.check_proton_table_edge(2999.0, 3000.0)
    with pytest.raises(UnsupportedCombinationError):  # C 430 MeV/u needs ~2.5 GeV in water
        XO.fdm_domain_status(CARBON12, [WATER], 430.0 * 12, 0.0, proton_edge_mev=1500.0)
    # the E0 + 6 sigma_E bound alone trips the check (nominal fits, upper does not)
    ok = XO.fdm_domain_status(OXYGEN16, [WATER], 300.0 * 16, 0.0)
    with pytest.raises(UnsupportedCombinationError):
        XO.fdm_domain_status(
            OXYGEN16,
            [WATER],
            300.0 * 16,
            20.0,
            proton_edge_mev=ok.p1v1_eq_kinetic_energy_mev * 1.001,
        )


@pytest.mark.parametrize(
    ("proj", "e_u", "approx_mev"),
    [
        (ALPHA, 150.0, 458.0),
        (ALPHA, 250.0, 807.0),
        (CARBON12, 290.0, 1486.0),
        (CARBON12, 430.0, 2471.0),
        (OXYGEN16, 350.0, 1876.0),
        (OXYGEN16, 430.0, 2482.0),
    ],
)
def test_fdm_extrapolated_flag(proj: Any, e_u: float, approx_mev: float) -> None:
    st = XO.fdm_domain_status(proj, [WATER], e_u * proj.a)
    assert st.fdm_extrapolated is True
    assert st.p1v1_eq_kinetic_energy_mev == pytest.approx(approx_mev, rel=0.02)  # P11b values
    assert st.p1v1_eq_kinetic_energy_mev <= st.p1v1_eq_kinetic_energy_hi_mev <= 4000.0
    d = st.to_dict()
    assert d["fdm_extrapolated"] is True and d["fdm_domain_edge_mev"] == 300.0


def test_fdm_flag_false_for_proton_and_persisted(make_config: Any) -> None:
    assert XO.fdm_domain_status(PROTON, [WATER], 150.0).fdm_extrapolated is False
    s = validate(make_config()).summary()["mcs_mapping"]
    assert s["fdm_extrapolated"] is False and s["p1v1_eq_kinetic_energy_mev"] is None


def test_ls_metadata_hash_rules(make_config: Any) -> None:
    cfg = make_config()

    def run(ls: str) -> Any:
        src = BetheStoppingSource(options=BetheOptions(ls=ls))
        return validate(replace(cfg, physics=replace(cfg.physics, stopping=src)))

    off, eff, mixed = run("off"), run("eff"), run("bare_mixed")
    assert off.tables.sha256 == eff.tables.sha256 == mixed.tables.sha256  # protons
    for e, req in ((off, "off"), (eff, "eff"), (mixed, "bare_mixed")):
        ph = e.summary()["physics"]
        assert (ph["ls_requested"], ph["ls_effective"]) == (req, "off")

    # ions: the mode is part of the table identity and effective = requested
    def ion_hash(ls: str) -> str:
        t = BetheStoppingSource(options=BetheOptions(ls=ls), e_max_per_u=20.0, points_per_decade=20)
        return TransportTables.__name__ + str(t.table(WATER, CARBON12).metadata["options"])

    assert ion_hash("eff") != ion_hash("off")
    ion_cfg = replace(
        off,
        requested=replace(off.requested, source=replace(off.requested.source, projectile=CARBON12)),
    )
    assert ion_cfg.ls_effective == "off" and ion_cfg.ls_requested == "off"
    src = BetheStoppingSource(options=BetheOptions(ls="eff"))
    ion_eff = replace(
        eff,
        requested=replace(
            eff.requested,
            source=replace(eff.requested.source, projectile=CARBON12),
            physics=replace(eff.requested.physics, stopping=src),
        ),
    )
    assert (ion_eff.ls_requested, ion_eff.ls_effective) == ("eff", "eff")


# --- I1: twin vs Warp CPU ----------------------------------------------------------------------
EXERCISED = {
    "kinetic_from_pv", "path_range", "path_energy", "ln_chi_a2", "fano_delta_b", "b_slab",
    "b_residual", "pv_equivalent", "fdm_equivalent",
    "pv_equivalent_table",  # the table evaluation: tests/ionmc/test_ion_mcs_table.py
}  # fmt: skip
NO = 20


def _kernel(real: Any) -> Any:
    F = XO.make_exact_omega(real)

    @wp.kernel(module="unique")
    def k(
        ion: wp.array(dtype=real),
        pro: wp.array(dtype=real),
        comp: wp.array(dtype=real),
        quad: wp.array(dtype=real),
        z: real,
        m: real,
        pv: wp.array(dtype=real),
        p1: wp.array(dtype=real),
        flag: int,
        out: wp.array2d(dtype=real),
        oi: wp.array2d(dtype=int),
    ):  # type: ignore[no-untyped-def]
        i = wp.tid()
        e = F.kinetic_from_pv(pv[i], m)
        out[i, 0] = e
        out[i, 1] = F.path_range(ion, e)
        out[i, 2] = F.path_energy(ion, out[i, 1])
        pc2 = e * (e + real(2.0) * m)
        b2 = pc2 / ((e + m) * (e + m))
        out[i, 3] = F.ln_chi_a2(comp, z, b2, pc2)
        out[i, 4] = F.fano_delta_b(comp, b2, pc2 / (m * m))
        b, sb = F.b_slab(ion, comp, quad, z, m, pv[i], p1[i])
        out[i, 5] = b
        peq, st, it = F.pv_equivalent(ion, pro, comp, quad, z, m, pv[i], p1[i], flag)
        out[i, 6] = peq
        g, sg = F.b_residual(pro, comp, quad, wp.log(peq), p1[i] / pv[i], b)
        out[i, 7] = g
        f, peq2, p1eq, st2 = F.fdm_equivalent(ion, pro, comp, quad, z, m, pv[i], p1[i], flag)
        out[i, 8] = f
        out[i, 9] = p1eq
        oi[i, 0] = sb
        oi[i, 1] = st
        oi[i, 2] = it
        oi[i, 3] = sg
        oi[i, 4] = st2

    return k


def _states() -> list[tuple[Any, Any, list[tuple[float, float]]]]:
    out = []
    for proj, mat in (
        (ALPHA, WATER),
        (CARBON12, WATER),
        (OXYGEN16, MATERIALS["bone_compact_icru"]),
    ):
        mp = mapping(proj, mat)
        out.append((proj, mp, [s for s in domain_states(mp) if s[1] / s[0] > 1.05]))
    return out


def test_i1_completeness() -> None:
    def names(ns: Any) -> set[str]:
        return {k for k, v in vars(ns).items() if callable(v)}

    assert names(python_twin(XO.make_exact_omega)) == EXERCISED
    assert names(XO.make_exact_omega(wp.float64)) == EXERCISED


@pytest.mark.parametrize("precision", ["float64", "float32"])
@pytest.mark.parametrize("device", ["cpu", pytest.param("cuda:0", marks=pytest.mark.cuda)])
def test_i1_twin_equals_warp_kernel(precision: str, device: str) -> None:
    real = wp.float64 if precision == "float64" else wp.float32
    dt = np.float64 if precision == "float64" else np.float32
    kern = _kernel(real)
    twin = python_twin(XO.make_exact_omega)
    rtol, atol = (1e-12, 1e-14) if precision == "float64" else (4e-5, 1e-6)
    for proj, mp, st in _states():
        pv = np.array([s[0] for s in st], dtype=dt).astype(np.float64)
        p1 = np.array([s[1] for s in st], dtype=dt).astype(np.float64)
        n = len(pv)

        def A(a: np.ndarray) -> Any:
            return wp.array(a.astype(dt), dtype=real, device=device)

        out = wp.zeros((n, NO), dtype=real, device=device)
        oi = wp.zeros((n, 5), dtype=int, device=device)
        wp.launch(
            kern, dim=n,
            inputs=[A(mp.ion_path.pack), A(mp.proton.pack), A(mp.comp), A(mp.quad),
                    real(proj.z), real(proj.mass_mev), A(pv), A(p1), 0, out, oi],
            device=device,
        )  # fmt: skip
        ko, ki = out.numpy().astype(np.float64), oi.numpy()
        r32 = (
            (lambda a: a)
            if precision == "float64"
            else (lambda a: a.astype(np.float32).astype(np.float64))
        )
        ion, pro, comp, quad = (
            r32(a).tolist() for a in (mp.ion_path.pack, mp.proton.pack, mp.comp, mp.quad)
        )
        z, m = float(real(proj.z)), float(np.float64(dt(proj.mass_mev)))
        for i in range(n):
            e = twin.kinetic_from_pv(pv[i], m)
            b, sb = twin.b_slab(ion, comp, quad, z, m, pv[i], p1[i])
            peq, st, _it = twin.pv_equivalent(ion, pro, comp, quad, z, m, pv[i], p1[i], 0)
            f, _pq, p1eq, st2 = twin.fdm_equivalent(ion, pro, comp, quad, z, m, pv[i], p1[i], 0)
            assert sb == ki[i, 0] == 0 and st == ki[i, 1] == 0 and st2 == ki[i, 4] == 0
            assert 0 < ki[i, 2] <= XO.MAX_ITERATIONS
            if precision == "float64":
                assert ko[i, 0] == pytest.approx(e, rel=rtol, abs=atol)
                tw = (b, peq, f, p1eq)
                kv = (ko[i, 5], ko[i, 6], ko[i, 8], ko[i, 9])
                for a, c in zip(tw, kv, strict=True):
                    assert abs(a - c) <= atol + rtol * abs(c)
                assert ko[i, 1] == pytest.approx(twin.path_range(ion, e), rel=rtol, abs=atol)
                assert abs(ko[i, 7]) <= 1e-10  # g at the kernel's own root vs its own b
            else:
                assert ko[i, 0] == pytest.approx(e, rel=rtol, abs=atol)
                assert ko[i, 1] == pytest.approx(twin.path_range(ion, e), rel=rtol, abs=atol)
                bud_b = 80.0 * U * (abs(b) + 60.0)  # 64-term sums, ~10 further roundings per term
                assert abs(ko[i, 5] - b) <= bud_b
                h = 1e-3
                u = math.log(peq)
                gp = [
                    twin.b_slab(
                        pro,
                        comp,
                        quad,
                        1.0,
                        PROTON.mass_mev,
                        math.exp(u + s * h),
                        p1[i] / pv[i] * math.exp(u + s * h),
                    )[0]
                    for s in (-1, 1)
                ]
                slope = (gp[1] - gp[0]) / (2 * h)
                assert slope > 0.5
                bud_u = 2.0 * bud_b / slope + 4 * U
                assert abs(math.log(ko[i, 6] / peq)) <= bud_u
                assert (
                    abs(ko[i, 8] - f) <= 1e-6 + 0.7 * bud_u * f + 4 * U * f
                )  # |df/dln pv| <= 0.7 f


def test_root_iteration_bound_and_status_codes() -> None:
    mp = mapping(CARBON12, WATER)
    its = []
    for pv, p1 in domain_states(mp):
        its.append(mp.pv_equivalent(pv, p1)[2])
    assert 1 <= max(its) <= 40  # cap is 100; observed <= 20
    # invalid ion state and a state below the ion table: status 1
    assert mp.pv_equivalent(100.0, 50.0)[1] == 1
    # an ion state whose b exceeds the proton b range (deep, below the domain): status 2, f = 0
    pv, p1 = mp.state_at_depth(5.0 * 12, 0.97)
    peq, st, _ = mp.pv_equivalent(pv, p1)
    assert st == 2 and mp.fdm_equivalent(pv, p1)[0] == 0.0
