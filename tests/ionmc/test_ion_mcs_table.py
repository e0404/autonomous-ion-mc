"""Tabulated exact-Omega mapping (decision 0043 section 3, "Evaluation", the allowed alternative):
the committed bound test of the ``pv_eq/pv`` table and its twin/kernel equality.

Decision 0043: the tabulation of ``pv_eq/pv`` read by bilinear interpolation is allowed only if a
committed test bounds the interpolation error of the resulting theta (or T) by <= 0.1 % over the
domain. :func:`test_table_variance_error_bound` is that test: the per-step variance (proportional
to T and theta^2) from the table against the root solve, relative error <= 1e-3, on random states
over the table domain of each ion (He, C, O) in water and in a tissue (skeletal muscle), including
the thin-slab edge, the lowest midpoint state and the band edges of the initial state.
In addition, every produced table is certified at build by a domain-wide analytic bound
(:func:`ionmc.physics.exact_omega.certify_pv_eq_table`: bilinear error bound from second
differences on a 2x-refined root-solve grid, safety factor 2, converted exactly into the variance
error of later and first steps over all physical end points, cross-checked by 5 x 5 interior
samples per cell); :func:`test_certification_every_builtin_and_user_material` runs it for every
built-in material and a user-defined one, and the forced-failure tests show it fails closed.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Any

import numpy as np
import pytest
import warp as wp

import ionmc
from ionmc._wpfunc import python_twin
from ionmc.config import (
    MCS_MAPPING_EVALS,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    validate,
)
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import MATERIALS, MUSCLE_SKELETAL_ICRP, WATER, Material
from ionmc.physics import exact_omega as XO
from ionmc.physics.ion_mcs import THIN_SLAB_R, make_ion_mcs
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, PROTON, Projectile
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.sources import PencilBeamSource

wp.config.log_level = wp.LOG_WARNING
D = wp.float64
IONS = (ALPHA, CARBON12, OXYGEN16)
DOMAIN_U = {"alpha": (50.0, 250.0), "carbon12": (80.0, 430.0), "oxygen16": (100.0, 430.0)}
BOUND = 1.0e-3
"""Decision 0043: the table's error of T (theta^2) is at most 0.1 %."""
USER_MATERIAL = Material(
    "user_mix_test",
    1.37,
    {"H": 0.08, "C": 0.40, "N": 0.05, "O": 0.40, "Ca": 0.07},
)
"""An arbitrary user-defined material (composition, density; I from Bragg additivity)."""
UNTABULABLE = {
    ("oxygen16", "copper", "_lower_case"),
    ("oxygen16", "lead", "_lower_case"),
}
"""(ion, material, case) whose table node build fails closed (a node of a band cell has no root:
oxygen at 100 MeV/u in copper and lead); the root option remains. Not a certification result."""
INV_XS = 0.0123
RHO = 1.0
S_MM = 0.7


def _config(
    proj: Projectile,
    material: Material,
    e0_mev: float,
    sigma_mev: float,
    evaluation: str = "table",
    backend: str = "python",
) -> SimulationConfig:
    return SimulationConfig(
        source=PencilBeamSource(proj, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), e0_mev, sigma_mev, 0.0),
        geometry=BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, 100.0), material),
        scoring=(ScoringGrid((-30.0, -30.0, 0.0), (60.0, 60.0, 10.0), (1, 1, 10)),),
        physics=PhysicsOptions(
            nuclear=False, stopping=BetheStoppingSource(), mcs_mapping_eval=evaluation
        ),
        run=RunOptions(backend=backend, precision="float64", seed=3, n_histories=4, n_batches=2),  # type: ignore[arg-type]
    )


def _upper_case(proj: Projectile, material: Material) -> SimulationConfig:
    """The domain's upper end: E0 + 6 sigma_E at the maximum energy per nucleon, sigma 1 %."""
    e0 = DOMAIN_U[proj.name][1] * proj.a / 1.0601
    return _config(proj, material, e0, 0.01 * e0)


def _lower_case(proj: Projectile, material: Material) -> SimulationConfig:
    """The domain's lower end, monoenergetic."""
    return _config(proj, material, DOMAIN_U[proj.name][0] * proj.a, 0.0)


def _compare_kernel() -> Any:
    IM = make_ion_mcs(D)

    @wp.kernel(module="unique")
    def k(
        ion: wp.array2d(dtype=D),
        pro: wp.array2d(dtype=D),
        comp: wp.array2d(dtype=D),
        quad: wp.array(dtype=D),
        tab: wp.array2d(dtype=D),
        z: D,
        mass: D,
        pv: wp.array(dtype=D),
        pv_end: wp.array(dtype=D),
        p1: wp.array(dtype=D),
        birth: wp.array(dtype=int),
        out: wp.array2d(dtype=D),
        st: wp.array2d(dtype=int),
    ):  # type: ignore[no-untyped-def]
        i = wp.tid()
        vr, sr = IM.mcs_variance(
            ion, pro, comp, quad, 0, z, mass, pv[i], pv_end[i], p1[i], birth[i], D(INV_XS),
            D(RHO), D(S_MM),
        )  # fmt: skip
        vt, stt = IM.mcs_variance_table(
            tab, 0, z, pv[i], pv_end[i], p1[i], birth[i], D(INV_XS), D(RHO), D(S_MM)
        )
        pr, p1r, s1 = IM.mapped_state(ion, pro, comp, quad, 0, z, mass, pv[i], p1[i])
        pt, p1t, s2 = IM.mapped_state_table(tab, 0, pv[i], p1[i])
        out[i, 0] = vr
        out[i, 1] = vt
        out[i, 2] = pr
        out[i, 3] = pt
        st[i, 0] = sr
        st[i, 1] = stt
        st[i, 2] = s1
        st[i, 3] = s2

    return k


def _evaluate(
    proj: Projectile, im: Any, pv: np.ndarray, pv_end: np.ndarray, p1: np.ndarray, birth: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    n = len(pv)
    out = wp.zeros((n, 4), dtype=D, device="cpu")
    st = wp.zeros((n, 4), dtype=int, device="cpu")

    def a2(x: np.ndarray) -> Any:
        return wp.array(np.ascontiguousarray(x[:1]), dtype=D, device="cpu")

    def a1(x: np.ndarray) -> Any:
        return wp.array(np.ascontiguousarray(x, dtype=np.float64), dtype=D, device="cpu")

    wp.launch(
        _compare_kernel(), dim=n,
        inputs=[a2(im.ion), a2(im.proton), a2(im.comp), a1(im.quad), a2(im.table),
                D(float(proj.z)), D(proj.mass_mev), a1(pv), a1(pv_end), a1(p1),
                wp.array(birth.astype(np.int32), dtype=int, device="cpu"), out, st],
        device="cpu",
    )  # fmt: skip
    return out.numpy(), st.numpy()


def _states(im: Any, rng: np.random.Generator, n: int) -> tuple[np.ndarray, ...]:
    """Random states over the table domain plus the edge states (``pv`` in [pv_lo, p1v1], ``p1v1``
    in the band of the initial states)."""
    pv_lo, p1_min, p1_max = im.table_domain
    m = 1.0 + XO.TABLE_PV_MARGIN
    lo1, hi1 = math.log(p1_min * m), math.log(p1_max / m)  # the run's initial states proper
    lp1 = rng.uniform(lo1, hi1, n)
    span = lp1 - math.log(pv_lo)
    kind = rng.random(n)
    frac = np.where(
        kind < 0.6,
        rng.uniform(0.0, 1.0, n),  # uniform in ln(p1v1/pv)
        np.exp(rng.uniform(math.log(1e-7), 0.0, n)),  # thin slabs (log-uniform)
    )
    lr = np.maximum(frac * span, math.log1p(THIN_SLAB_R))
    p1 = np.exp(lp1)
    pv = p1 / np.exp(lr)
    # edges: thin-slab guard (pv = p1v1 and r = 1 + 1e-6), lowest midpoint, band edges, grid nodes
    edge_p1 = np.array([p1_min * m, p1_max / m, math.sqrt(p1_min * p1_max)])
    e_pv: list[float] = []
    e_p1: list[float] = []
    for q in edge_p1:
        for pvv in (q, q / (1.0 + THIN_SLAB_R), q / (1.0 + 2e-6), q / 1.001, pv_lo, pv_lo * 1.0001):
            e_pv.append(pvv)
            e_p1.append(q)
    h = im.table[0]
    xs = h[0] + np.arange(int(h[2])) / h[1]
    for x in xs[::7]:
        pvx = math.exp(x)
        if pv_lo <= pvx <= p1_max / m:
            e_pv.append(min(pvx, p1_max / m))
            e_p1.append(p1_max / m)
    pv = np.concatenate((pv, e_pv))
    p1 = np.concatenate((p1, e_p1))
    return pv, p1


@pytest.mark.parametrize("material", [WATER, MUSCLE_SKELETAL_ICRP], ids=lambda m: m.name)
@pytest.mark.parametrize("proj", IONS, ids=lambda p: p.name)
def test_table_variance_error_bound(proj: Projectile, material: Material) -> None:
    """Decision 0043 bound: |var_table / var_root - 1| <= 1e-3 (T, theta^2) on >= 2000 random
    states (plus edge states) over the table domain of the upper and the lower domain end; every
    state maps with status 0 on both paths. The measured maximum is printed."""
    rng = np.random.default_rng(20261108)
    worst = 0.0
    worst_q = 0.0
    total = 0
    for cfg, n in ((_upper_case(proj, material), 1600), (_lower_case(proj, material), 600)):
        im = validate(cfg).ion_mapping
        assert im is not None and im.evaluation == "table" and im.table is not None
        pv, p1 = _states(im, rng, n)
        for birth in (0, 1):
            pv_end = pv * (1.0 - 0.01 * rng.random(pv.size)) if birth else pv
            b = np.full(pv.size, birth)
            out, st = _evaluate(proj, im, pv, pv_end, p1, b)
            assert np.all(st == 0), (proj.name, material.name, np.argwhere(st != 0)[:5])
            assert np.all(out[:, 0] > 0.0)
            rel = np.abs(out[:, 1] / out[:, 0] - 1.0)
            worst = max(worst, float(rel.max()))
            worst_q = max(worst_q, float(np.max(np.abs(out[:, 3] / out[:, 2] - 1.0))))
            total += pv.size
    print(
        f"table bound {proj.name} in {material.name}: max |var_table/var_root - 1| = {worst:.3e} "
        f"(pv_eq: {worst_q:.3e}) over {total} state evaluations"
    )
    assert total >= 2 * 2000
    assert worst <= BOUND


def test_table_outside_domain_fails_closed() -> None:
    """States outside the table (below the lowest midpoint, above the initial-state band, beyond
    the memory-ratio axis) return status 5 and variance 0 (counted as a mapping failure)."""
    cfg = _lower_case(CARBON12, WATER)
    im = validate(cfg).ion_mapping
    assert im is not None and im.table_domain is not None
    pv_lo, _p1_min, p1_max = im.table_domain
    p1 = np.array([p1_max / 1.01, p1_max * 1.5, p1_max / 1.01])
    pv = np.array([pv_lo / 1.5, p1_max / 1.6, p1_max / 1.01 / 1.0001])
    out, st = _evaluate(CARBON12, im, pv, pv, p1, np.zeros(3, dtype=int))
    assert list(st[:2, 1]) == [5, 5] and list(st[:2, 3]) == [5, 5]
    assert np.all(out[:2, 1] == 0.0)
    assert st[2, 1] == 0 and out[2, 1] > 0.0  # a state inside stays valid
    # a cell touching an invalid node: status 5 as well
    tab = im.table.copy()
    tab[0, XO.TABLE_HEADER :] = XO.TABLE_INVALID_NODE
    bad = replace(im, table=tab)
    out, st = _evaluate(CARBON12, bad, pv[2:], pv[2:], p1[2:], np.zeros(1, dtype=int))
    assert st[0, 1] == 5 and out[0, 1] == 0.0


def test_table_nodes_are_the_root_solve() -> None:
    """The table nodes are pv_eq/pv of the root solve: the python twin of ``pv_equivalent`` at
    the nodes reproduces them (float64, relative 1e-13); the table is cached per input."""
    cfg = _upper_case(ALPHA, WATER)
    im = validate(cfg).ion_mapping
    assert im is not None
    h = im.table[0]
    n_x, n_y = int(h[2]), int(h[5])
    q = im.table[0, XO.TABLE_HEADER : XO.TABLE_HEADER + n_x * n_y].reshape(n_x, n_y)
    xo = python_twin(XO.make_exact_omega)
    ion, pro, comp = (a[0].tolist() for a in (im.ion, im.proton, im.comp))
    quad = im.quad.tolist()
    valid = np.argwhere(q > 0.0)
    assert len(valid) > 100
    for i, j in valid[:: max(1, len(valid) // 12)]:
        x = h[0] + i * (h[6] - h[0]) / (n_x - 1)
        y = h[3] + j * (h[7] - h[3]) / (n_y - 1)
        pv = math.exp(x)
        peq, st, _it = xo.pv_equivalent(
            ion, pro, comp, quad, 2.0, ALPHA.mass_mev, pv, pv * math.exp(y), 0
        )
        assert st == 0
        assert peq / pv == pytest.approx(q[i, j], rel=1e-13)
    again = validate(cfg).ion_mapping
    assert again is not None and np.array_equal(again.table, im.table)


def _table_kernel(real: Any) -> Any:
    IM = make_ion_mcs(real)

    @wp.kernel(module="unique")
    def k(
        tab: wp.array2d(dtype=real),
        z: real,
        pv: wp.array(dtype=real),
        pv_end: wp.array(dtype=real),
        p1: wp.array(dtype=real),
        birth: wp.array(dtype=int),
        out: wp.array2d(dtype=real),
        oi: wp.array(dtype=int),
    ):  # type: ignore[no-untyped-def]
        i = wp.tid()
        v, st = IM.mcs_variance_table(
            tab, 0, z, pv[i], pv_end[i], p1[i], birth[i], real(INV_XS), real(RHO), real(S_MM)
        )
        peq, p1eq, st2 = IM.mapped_state_table(tab, 0, pv[i], p1[i])
        out[i, 0] = v
        out[i, 1] = peq
        out[i, 2] = p1eq
        oi[i] = st * 10 + st2

    return k


@pytest.mark.parametrize("precision", ["float64", "float32"])
@pytest.mark.parametrize("proj", IONS, ids=lambda p: p.name)
def test_table_twin_equals_warp_cpu(proj: Projectile, precision: str) -> None:
    """The table path (``mapped_state_table``, ``mcs_variance_table``): python twin against a Warp
    CPU kernel; float64 relative 1e-12 on all states including the thin-slab edge, float32
    relative 1e-5 (twin on the float32-rounded inputs and table) on the states with
    ``p1v1/pv >= 1.01`` off the lowest table edge. Below that the float32 memory term
    ``1 - (pv/p1v1)^2`` of f_dM loses ``~6e-8 / (2 (r - 1))`` relative accuracy by cancellation
    (3e-2 at r - 1 = 1e-6), a property of f_dM in float32, not of the table; the transport
    kernel evaluates the mapping in float64 for every precision."""
    real = wp.float64 if precision == "float64" else wp.float32
    dt = np.float64 if precision == "float64" else np.float32
    rtol = 1e-12 if precision == "float64" else 1e-5
    im = validate(_upper_case(proj, WATER)).ion_mapping
    assert im is not None
    pv, p1 = _states(im, np.random.default_rng(11), 60)
    if precision == "float32":
        # well-conditioned states, not on a table edge (float32 rounding of an edge state
        # crosses the edge by ~1e-6 grid units, above the float64 edge tolerance of the twin)
        assert im.table_domain is not None
        keep = (p1 / pv >= 1.01) & (pv >= im.table_domain[0] * 1.001)
        pv, p1 = pv[keep], p1[keep]
    n = pv.size
    assert n >= 30
    birth = (np.arange(n) % 3 == 0).astype(np.int32)
    pv_end = pv * np.where(birth == 1, 0.995, 1.0)
    cols = [np.asarray(c, dtype=dt) for c in (pv, pv_end, p1)]
    tab = im.table.astype(dt)
    out = wp.zeros((n, 3), dtype=real, device="cpu")
    oi = wp.zeros(n, dtype=int, device="cpu")
    wp.launch(
        _table_kernel(real), dim=n,
        inputs=[wp.array(tab, dtype=real, device="cpu"), real(proj.z),
                *[wp.array(c, dtype=real, device="cpu") for c in cols],
                wp.array(birth, dtype=int, device="cpu"), out, oi],
        device="cpu",
    )  # fmt: skip
    got, status = out.numpy().astype(np.float64), oi.numpy()
    tw = python_twin(make_ion_mcs)
    tl = [tab[0].astype(np.float64).tolist()]
    for i in range(n):
        a, b, c = (float(x[i]) for x in cols)
        v, s1 = tw.mcs_variance_table(
            tl, 0, float(proj.z), a, b, c, int(birth[i]), float(dt(INV_XS)), float(dt(RHO)),
            float(dt(S_MM)),
        )  # fmt: skip
        peq, p1eq, s2 = tw.mapped_state_table(tl, 0, a, c)
        assert s1 == s2 == 0 and status[i] == 0, (i, s1, s2, status[i])
        for x, y in ((v, got[i, 0]), (peq, got[i, 1]), (p1eq, got[i, 2])):
            assert abs(x - y) <= rtol * abs(x), (precision, proj.name, i, x, y)


def test_option_validation_and_summary() -> None:
    """``mcs_mapping_eval`` fails closed on unknown values; the summary records the evaluation,
    the table sizes, axes bounds, domain and sha256; ``"root"`` carries no table; protons ignore
    the option (no mapping, summary unchanged)."""
    assert MCS_MAPPING_EVALS == ("table", "root")
    for bad in ("Table", "interp", "", 1):
        with pytest.raises(UnsupportedCombinationError):
            PhysicsOptions(nuclear=False, stopping=BetheStoppingSource(), mcs_mapping_eval=bad)  # type: ignore[arg-type]
    eff = validate(_lower_case(CARBON12, WATER))
    s = eff.summary()["mcs_mapping"]
    assert s["evaluation"] == "table"
    t = s["table"]
    assert (t["n_pv"], t["n_r"]) == (XO.TABLE_N_PV, XO.TABLE_N_R)
    assert t["ln_r_range"][0] == pytest.approx(math.log1p(THIN_SLAB_R), rel=1e-12)
    assert t["pv_lo_mev"] < t["p1v1_min_mev"] < t["p1v1_max_mev"]
    assert len(t["sha256"]) == 64
    root = validate(_config(CARBON12, WATER, 80.0 * 12, 0.0, evaluation="root"))
    assert root.ion_mapping is not None and root.ion_mapping.table is None
    rs = root.summary()["mcs_mapping"]
    assert rs["evaluation"] == "root" and "table" not in rs
    for ev in MCS_MAPPING_EVALS:
        p = validate(_config(PROTON, WATER, 100.0, 0.0, evaluation=ev))
        assert p.ion_mapping is None
        assert "evaluation" not in p.summary()["mcs_mapping"]


def test_table_build_fails_closed_beyond_the_paths() -> None:
    """A band of initial states whose cells need nodes beyond the ion path cannot be tabulated:
    ``UnsupportedCombinationError`` (the root option remains)."""
    mp = XO.ExactOmegaMapping.build(CARBON12, WATER)
    e_top = mp.ion_path.e_max_mev
    dom = XO.table_domain(CARBON12, 24.0, e_top * 0.99, e_top * 0.9999)
    with pytest.raises(UnsupportedCombinationError, match="mcs_mapping_eval='root'"):
        XO.build_pv_eq_table(
            mp.ion_path.pack, mp.proton.pack, mp.comp, mp.quad, 6.0, CARBON12.mass_mev, *dom
        )
    with pytest.raises(UnsupportedCombinationError):
        XO.build_pv_eq_table(
            mp.ion_path.pack, mp.proton.pack, mp.comp, mp.quad, 6.0, CARBON12.mass_mev,
            100.0, 50.0, 200.0,
        )  # fmt: skip


def test_certificate_constants() -> None:
    """The threshold is decision 0043's 0.1 % on the variance (T, theta^2); safety factor 2 on the
    second-difference estimate; 5 x 5 samples per cell."""
    assert XO.TABLE_VARIANCE_BOUND == BOUND
    assert XO.TABLE_SAFETY_FACTOR == 2.0
    assert XO.TABLE_SAMPLE_N == 5
    assert XO.TABLE_MIN_LOSS == 1.0e-4
    assert XO.LOG10_E == pytest.approx(math.log10(math.e), rel=1e-15)


def test_fdm_sensitivity_maximum_is_at_a_corner() -> None:
    """``fdm_sensitivity_max`` (|c2 + c3 l1| / f over an (l1, l2) rectangle, exact at the corners)
    is not exceeded on a dense grid of each rectangle; a rectangle where f <= 0 gives inf."""
    rng = np.random.default_rng(5)
    for _ in range(200):
        l1a = rng.uniform(-7.0, -0.01)
        l1b = min(l1a + rng.uniform(0.0, 2.0), 0.0)
        l2a = rng.uniform(0.0, 3.5)
        l2b = l2a + rng.uniform(0.0, 0.5)
        got = float(XO.fdm_sensitivity_max(*(np.array([v]) for v in (l1a, l1b, l2a, l2b)))[0])
        g1, g2 = np.meshgrid(np.linspace(l1a, l1b, 41), np.linspace(l2a, l2b, 41))
        c0, c1, c2, c3 = XO.FDM_COEFFICIENTS
        f = c0 + c1 * g1 + c2 * g2 + c3 * g1 * g2
        if np.all(f > 0.0):
            brute = float(np.max(np.abs(c2 + c3 * g1) / f))
            assert brute <= got * (1.0 + 1e-12)
            assert got == pytest.approx(brute, rel=1e-12)  # a grid corner is the maximum
        elif np.any(f <= 0.0):
            assert got == math.inf


def _check_certificate(im: Any, label: str) -> Any:
    assert im is not None and im.certificates is not None and len(im.certificates) == 1
    assert im.table_domain is not None
    c = im.certificates[0]
    h = im.table[0]
    xs = np.linspace(h[0], h[6], int(h[2]))
    ys = np.linspace(h[3], h[7], int(h[5]))
    n_band = int(XO._band_cells(xs, ys, *im.table_domain[1:]).sum())
    assert c.n_cells == n_band > 100
    assert c.n_sampled_points == XO.TABLE_SAMPLE_N**2 * c.n_cells
    assert c.n_cells < c.n_refined_points <= 9 * c.n_cells
    assert c.bound == XO.TABLE_VARIANCE_BOUND and c.safety_factor == XO.TABLE_SAFETY_FACTOR
    assert c.certified_bound_analytic == max(c.bound_later_steps, c.bound_birth_step)
    assert 0.0 < c.observed_max_sampled <= c.certified_bound_analytic <= XO.TABLE_VARIANCE_BOUND
    print(
        f"certified {label}: analytic bound {c.certified_bound_analytic:.3e} (later steps "
        f"{c.bound_later_steps:.3e}, first step {c.bound_birth_step:.3e}), observed sampled max "
        f"{c.observed_max_sampled:.3e}; {c.n_cells} cells, {c.n_refined_points} stencil roots, "
        f"{c.n_sampled_points} samples"
    )
    return c


@pytest.mark.parametrize("material", [*MATERIALS.values(), USER_MATERIAL], ids=lambda m: m.name)
def test_certification_every_builtin_and_user_material(material: Material) -> None:
    """Every built-in material (the ``ionmc.materials.MATERIALS`` registry) and a user-defined one:
    the tables of He, C and O at the upper and the lower domain end are certified at build with
    the analytic bound <= ``TABLE_VARIANCE_BOUND`` and the sampled maximum <= the analytic bound
    (per cell, enforced by the build); the certificate covers every solved cell. The values are
    printed."""
    for proj in IONS:
        for case in (_upper_case, _lower_case):
            label = f"{proj.name} in {material.name} ({case.__name__[1:6]})"
            if (proj.name, material.name, case.__name__) in UNTABULABLE:
                # pre-existing node-build limit (before any certification): fails closed
                with pytest.raises(UnsupportedCombinationError, match="cannot be built"):
                    validate(case(proj, material))
                print(f"not tabulable {label}: the node build fails closed (root option)")
                continue
            im = validate(case(proj, material)).ion_mapping
            _check_certificate(im, label)


def test_certification_fails_closed_above_the_threshold(monkeypatch: pytest.MonkeyPatch) -> None:
    """A threshold below the analytic bound fails closed (``UnsupportedCombinationError``,
    suggesting the root option), for a fresh build and for a cached table; the root option still
    validates."""
    cfg = _lower_case(ALPHA, WATER)
    validate(cfg)  # cached and certified
    monkeypatch.setattr(XO, "TABLE_VARIANCE_BOUND", 1.0e-9)
    with pytest.raises(UnsupportedCombinationError, match="mcs_mapping_eval='root'"):
        validate(cfg)  # cached table, re-checked
    monkeypatch.setattr(XO, "_TABLES", {})
    with pytest.raises(UnsupportedCombinationError, match="interpolation bound"):
        validate(cfg)  # fresh build
    root = validate(_config(ALPHA, WATER, 50.0 * 4, 0.0, evaluation="root"))
    assert root.ion_mapping is not None and root.ion_mapping.table is None


def test_certification_fails_closed_when_samples_exceed_the_analytic_bound(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The sampled cross-check falsifies an analytic bound that is too small: with the safety
    factor forced to 0 (bound = the root tolerance only) the build fails closed."""
    monkeypatch.setattr(XO, "_TABLES", {})
    monkeypatch.setattr(XO, "TABLE_SAFETY_FACTOR", 0.0)
    with pytest.raises(UnsupportedCombinationError, match="smoothness assumption"):
        validate(_lower_case(OXYGEN16, WATER))


def test_certification_detects_a_corrupted_node() -> None:
    """The analytic bound is computed from fresh root solves, so a corrupted table node (+1 %) is
    caught by the sampled cross-check (observed error above the analytic bound of its cell); an
    invalid node makes the certification fail closed as well."""
    im = validate(_lower_case(CARBON12, WATER)).ion_mapping
    assert im is not None and im.table_domain is not None
    args = (im.ion[0], im.proton[0], im.comp[0], im.quad, 6.0, CARBON12.mass_mev)
    row = im.table[0].copy()
    good = XO.certify_pv_eq_table(row, *args, *im.table_domain[1:])
    assert good == im.certificates[0]  # deterministic
    XO.check_table_certificate(good)
    h = row[: XO.TABLE_HEADER]
    n_y = int(h[5])
    xs = np.linspace(h[0], h[6], int(h[2]))
    ys = np.linspace(h[3], h[7], n_y)
    i, j = np.argwhere(XO._band_cells(xs, ys, *im.table_domain[1:]))[len(xs) // 2]
    node = XO.TABLE_HEADER + int(i) * n_y + int(j)
    bad = row.copy()
    bad[node] *= 1.01
    with pytest.raises(UnsupportedCombinationError, match="exceeds the analytic bound"):
        XO.certify_pv_eq_table(bad, *args, *im.table_domain[1:])
    bad[node] = XO.TABLE_INVALID_NODE
    with pytest.raises(UnsupportedCombinationError, match="cannot be certified"):
        XO.certify_pv_eq_table(bad, *args, *im.table_domain[1:])


@pytest.mark.parametrize("proj", IONS, ids=lambda p: p.name)
def test_birth_step_endpoint_dependence_is_covered(proj: Projectile) -> None:
    """The first step's ratio var_table / var_root depends on the step's end point (through the
    memory term l1, c3 != 0): it differs between pv_end = pv_mid, pv_end at the table's lowest
    pv (cutoff side) and pv_end -> 0, and every one of them is within the certified first-step
    bound (the certificate covers pv_end in (0, min(pv_mid, (1 - 1e-4) p1v1)])."""
    im = validate(_lower_case(proj, WATER)).ion_mapping
    assert im is not None and im.table_domain is not None and im.certificates is not None
    cert = im.certificates[0]
    pv, p1 = _states(im, np.random.default_rng(17), 400)
    keep = p1 * (1.0 - XO.TABLE_MIN_LOSS) >= pv  # the contract scope (TABLE_MIN_LOSS)
    pv, p1 = pv[keep], p1[keep]
    ones = np.ones(pv.size, dtype=int)
    ratios = []
    for pv_end in (pv, np.minimum(pv, im.table_domain[0]), pv * 1e-6):
        out, st = _evaluate(proj, im, pv, pv_end, p1, ones)
        assert np.all(st == 0) and np.all(out[:, 0] > 0.0)
        r = out[:, 1] / out[:, 0]
        assert np.max(np.abs(r - 1.0)) <= cert.bound_birth_step
        ratios.append(r)
    assert np.max(np.abs(ratios[0] - ratios[2])) > 0.0  # not endpoint-independent


def test_summary_and_capabilities_carry_the_certificate() -> None:
    """The summary's ``mcs_mapping.table`` block records the certificate (analytic bound and
    sampled maximum, overall and per material, cells, stencil and sample points, safety factor,
    threshold, contract) for a two-material geometry; the capability report advertises the
    contract with its threshold."""
    from ionmc.geometry import VoxelGeometry

    geo = VoxelGeometry(
        origin_mm=(-30.0, -30.0, 0.0),
        spacing_mm=(60.0, 60.0, 50.0),
        shape=(1, 1, 2),
        materials=(WATER, MUSCLE_SKELETAL_ICRP),
        material_index=np.array([0, 1], dtype=np.int32).reshape(1, 1, 2),
    )
    cfg = replace(_lower_case(OXYGEN16, WATER), geometry=geo)
    eff = validate(cfg)
    im = eff.ion_mapping
    assert im is not None and im.certificates is not None and len(im.certificates) == 2
    cs = im.certificates
    t = eff.summary()["mcs_mapping"]["table"]
    per = t["certified_bound_analytic_per_material"]
    assert per == [c.certified_bound_analytic for c in cs]
    assert t["certified_bound_analytic"] == max(per) <= t["variance_bound"] == BOUND
    obs = t["observed_max_sampled_per_material"]
    assert obs == [c.observed_max_sampled for c in cs]
    assert t["observed_max_sampled"] == max(obs) <= t["certified_bound_analytic"]
    assert t["certified_cells"] == sum(c.n_cells for c in cs)
    assert t["refined_points"] == sum(c.n_refined_points for c in cs)
    assert t["sampled_points"] == 25 * t["certified_cells"]
    assert t["samples_per_cell"] == 25 and t["safety_factor"] == 2.0
    assert t["refinement"].startswith("2x")
    assert t["certification"] == XO.TABLE_CERT_CONTRACT
    for key in ("every produced table", "2x-refined", "safety factor 2", "5 x 5", "first step"):
        assert key in t["certification"]
    cap = ionmc.capabilities()["species_details"]
    for proj in IONS:
        tc = cap[proj.name]["mcs"]["table_certification"]
        assert tc["scope"] == "certified at build for every produced table"
        assert tc["contract"] == XO.TABLE_CERT_CONTRACT
        assert tc["variance_bound"] == BOUND and tc["safety_factor"] == 2.0
        assert tc["samples_per_cell"] == 25
        assert "certified at build for every produced table" in cap[proj.name]["mcs"]["evaluation"]
