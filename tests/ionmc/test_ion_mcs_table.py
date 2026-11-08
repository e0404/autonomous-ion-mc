"""Table-accelerated exact-Omega mapping (decision 0043 section 3, "Evaluation"; Codex
REVIEW-ee53867c): ``mcs_mapping_eval = "table"`` solves the same root as ``"root"`` on every step,
to the same tolerance and with the same status codes, started from a tight bracket around the
bilinear ``pv_eq/pv`` table value; the table is an accelerator, its interpolation error a
performance diagnostic. The decision-0043 bound test is therefore the per-step root guarantee:
:func:`test_accelerated_root_equals_root` checks, on >= 2000 random states per ion and material
(every built-in material and a user material, both domain ends, thin-slab and band-edge states,
later and first steps with both end-point extremes), that the accelerated result equals the root
path's within the solver tolerance (pv_eq relative <= 1e-8, variance ratio - 1 <= 1e-8 times
the f_dM sensitivity max(1, |c2 + c3 l1| / f), equal status). The twin/kernel equality,
fall-back, corrupted-table and summary tests follow.
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
from ionmc.physics.em import FDM_COEFFICIENTS
from ionmc.physics.ion_mcs import THIN_SLAB_R, make_ion_mcs
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, PROTON, Projectile
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.sources import PencilBeamSource

wp.config.log_level = wp.LOG_WARNING
D = wp.float64
IONS = (ALPHA, CARBON12, OXYGEN16)
DOMAIN_U = {"alpha": (50.0, 250.0), "carbon12": (80.0, 430.0), "oxygen16": (100.0, 430.0)}
PV_EQ_TOL = 1.0e-8
"""Agreement of the accelerated and the root path's pv_eq (both converge to the solver's residual
tolerance 1e-13; observed differences <= ~1e-9)."""
VAR_TOL = 1.0e-8
"""Agreement of the step variances, relative, divided by the f_dM sensitivity
max(1, |c2 + c3 l1| / f) (the exact propagation of a pv_eq difference; it exceeds 1 only near the
fit's zero, first steps ending close to p1v1)."""
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
oxygen at 100 MeV/u in copper and lead); the root option remains."""
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
        vt, stt, miss = IM.mcs_variance_table(
            ion, pro, comp, quad, tab, 0, z, mass, pv[i], pv_end[i], p1[i], birth[i], D(INV_XS),
            D(RHO), D(S_MM),
        )  # fmt: skip
        pr, p1r, s1 = IM.mapped_state(ion, pro, comp, quad, 0, z, mass, pv[i], p1[i])
        pt, p1t, s2, miss2 = IM.mapped_state_table(
            ion, pro, comp, quad, tab, 0, z, mass, pv[i], p1[i]
        )
        pi, p1i, s3 = IM.interpolated_state(tab, 0, pv[i], p1[i])
        out[i, 0] = vr
        out[i, 1] = vt
        out[i, 2] = pr
        out[i, 3] = pt
        out[i, 4] = pi
        st[i, 0] = sr
        st[i, 1] = stt
        st[i, 2] = s1
        st[i, 3] = s2
        st[i, 4] = s3
        st[i, 5] = miss
        st[i, 6] = miss2

    return k


def _evaluate(
    proj: Projectile, im: Any, pv: np.ndarray, pv_end: np.ndarray, p1: np.ndarray, birth: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Columns: var root, var table, pv_eq root, pv_eq table, pv_eq interpolated; status root,
    status table, mapped status root, mapped status table, interpolation status, miss (variance),
    miss (mapped state)."""
    n = len(pv)
    out = wp.zeros((n, 5), dtype=D, device="cpu")
    st = wp.zeros((n, 7), dtype=int, device="cpu")

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


def _compare(proj: Projectile, im: Any, pv: np.ndarray, p1: np.ndarray, pv_end: np.ndarray,
             birth: int) -> tuple[float, float, int]:  # fmt: skip
    """Accelerated against root: equal status everywhere, pv_eq and variance within the solver
    tolerance where converged. Returns (max pv_eq difference, max variance difference, misses)."""
    out, st = _evaluate(proj, im, pv, pv_end, p1, np.full(pv.size, birth))
    assert np.array_equal(st[:, 0], st[:, 1]) and np.array_equal(st[:, 2], st[:, 3])
    ok = st[:, 2] == 0
    assert np.all(ok), (proj.name, np.argwhere(~ok)[:5])  # the table domain maps everywhere
    # f_dM is clamped to 0 where its fit vanishes (first steps that end very close to p1v1):
    # both paths then give 0; elsewhere the variances agree relatively
    pos = out[:, 0] > 0.0
    assert np.all(out[~pos, 1] == 0.0) and np.all(out[~pos, 0] == 0.0)
    assert np.mean(pos) > 0.5
    d_q = float(np.max(np.abs(out[:, 3] / out[:, 2] - 1.0)))
    # the variance is f_dM = c0 + c1 l1 + c2 l2 + c3 l1 l2 (l2 = log10 pv_eq) times a common
    # prefactor, so a pv_eq difference d propagates exactly as (c2 + c3 l1) log10(1 + d) / f:
    # amplified by S = |c2 + c3 l1| / f near the fit's zero. Compare with S factored out.
    c0, c1, c2, c3 = FDM_COEFFICIENTS
    pv_map = np.minimum(pv, p1 / (1.0 + THIN_SLAB_R))
    with np.errstate(divide="ignore", invalid="ignore"):  # om = 0 only where var = 0 (masked)
        if birth:
            l1 = np.log10(1.0 - (pv_end / p1) ** 2) - math.log10(math.e)
        else:
            l1 = np.log10(1.0 - (pv_map / p1) ** 2)
        f = c0 + c1 * l1 + (c2 + c3 * l1) * np.log10(out[:, 2])
        sens = np.maximum(np.abs(c2 + c3 * l1) / np.where(f > 0.0, f, np.inf), 1.0)
    d_v = float(np.max(np.abs(out[pos, 1] / out[pos, 0] - 1.0) / sens[pos]))
    assert np.array_equal(st[:, 5], st[:, 6])
    return d_q, d_v, int(st[:, 6].sum())


@pytest.mark.parametrize("material", [*MATERIALS.values(), USER_MATERIAL], ids=lambda m: m.name)
def test_accelerated_root_equals_root(material: Material) -> None:
    """The decision-0043 test of the table evaluation (an exact root on every step): for He, C
    and O at the upper and the lower domain end of every built-in material and a user material,
    >= 2000 random states plus thin-slab, band-edge and node states, later steps and first steps
    with both end-point extremes (pv_end = pv_mid and pv_end -> 0): the table-accelerated result
    equals the root path's (same status; pv_eq within 1e-8, variance within 1e-8 x the f_dM
    sensitivity, ``VAR_TOL``). Bracket misses
    (a diagnostic) stay below 0.1 % of the states. The maxima are printed. Oxygen at 100 MeV/u
    in copper and lead cannot be tabulated (``UNTABULABLE``): the build fails closed."""
    rng = np.random.default_rng(20261108)
    for proj in IONS:
        for case in (_upper_case, _lower_case):
            label = f"{proj.name} in {material.name} ({case.__name__[1:6]})"
            if (proj.name, material.name, case.__name__) in UNTABULABLE:
                with pytest.raises(UnsupportedCombinationError, match="cannot be built"):
                    validate(case(proj, material))
                print(f"not tabulable {label}: the build fails closed (root option)")
                continue
            im = validate(case(proj, material)).ion_mapping
            assert im is not None and im.evaluation == "table" and im.table is not None
            pv, p1 = _states(im, rng, 2000)
            assert pv.size >= 2000
            worst_q = worst_v = 0.0
            misses = 0
            for birth, pv_end in ((0, pv), (1, pv), (1, pv * 1e-6)):
                d_q, d_v, n_miss = _compare(proj, im, pv, p1, pv_end, birth)
                worst_q, worst_v = max(worst_q, d_q), max(worst_v, d_v)
                misses += n_miss
            print(
                f"accelerated root {label}: max |pv_eq/pv_eq_root - 1| = {worst_q:.2e}, max "
                f"|var/var_root - 1| = {worst_v:.2e}, bracket misses {misses} of {3 * pv.size}"
            )
            assert worst_q <= PV_EQ_TOL and worst_v <= VAR_TOL
            assert misses <= 1e-3 * 3 * pv.size


def test_outside_table_states_fall_back_to_the_root() -> None:
    """States outside the table (below the lowest midpoint, above the initial-state band) have no
    starting guess (interpolation status 5): the accelerated path uses the full bracket (a
    counted bracket miss) and returns the root path's result; a cell with an invalid node does
    the same."""
    cfg = _lower_case(CARBON12, WATER)
    im = validate(cfg).ion_mapping
    assert im is not None and im.table_domain is not None
    pv_lo, _p1_min, p1_max = im.table_domain
    p1 = np.array([p1_max / 1.01, p1_max * 1.5, p1_max / 1.01])
    pv = np.array([pv_lo / 1.5, p1_max / 1.6, p1_max / 1.01 / 1.0001])
    out, st = _evaluate(CARBON12, im, pv, pv, p1, np.zeros(3, dtype=int))
    assert list(st[:, 4]) == [5, 5, 0]  # interpolation
    assert list(st[:, 6]) == [1, 1, 0]  # bracket misses
    assert np.array_equal(st[:, 0], st[:, 1])
    ok = st[:, 0] == 0
    assert np.all(np.abs(out[ok, 3] / out[ok, 2] - 1.0) <= PV_EQ_TOL)
    tab = im.table.copy()
    tab[0, XO.TABLE_HEADER :] = XO.TABLE_INVALID_NODE
    bad = replace(im, table=tab)
    out, st = _evaluate(CARBON12, bad, pv[2:], pv[2:], p1[2:], np.zeros(1, dtype=int))
    assert st[0, 4] == 5 and st[0, 6] == 1 and st[0, 1] == st[0, 0] == 0
    assert abs(out[0, 3] / out[0, 2] - 1.0) <= PV_EQ_TOL


def test_a_corrupted_table_cannot_change_the_result() -> None:
    """The table only seeds the root: corrupted nodes (+5 %) or a vanishing starting bracket
    change the number of iterations and bracket misses, never the result (equal to the root
    path within the solver tolerance)."""
    im = validate(_upper_case(OXYGEN16, WATER)).ion_mapping
    assert im is not None
    pv, p1 = _states(im, np.random.default_rng(3), 300)
    for mutate in ("nodes", "width"):
        tab = im.table.copy()
        if mutate == "nodes":
            tab[0, XO.TABLE_HEADER :] *= 1.05
        else:
            tab[0, 8] = 1e-14
        bad = replace(im, table=tab)
        d_q, d_v, misses = _compare(OXYGEN16, bad, pv, p1, pv, 0)
        print(f"corrupted table ({mutate}): {misses} bracket misses of {pv.size}")
        assert d_q <= PV_EQ_TOL and d_v <= VAR_TOL
        assert misses > 0


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


def _twin_kernel(real: Any) -> Any:
    IM = make_ion_mcs(real)

    @wp.kernel(module="unique")
    def k(
        ion: wp.array2d(dtype=real),
        pro: wp.array2d(dtype=real),
        comp: wp.array2d(dtype=real),
        quad: wp.array(dtype=real),
        tab: wp.array2d(dtype=real),
        z: real,
        mass: real,
        pv: wp.array(dtype=real),
        pv_end: wp.array(dtype=real),
        p1: wp.array(dtype=real),
        birth: wp.array(dtype=int),
        out: wp.array2d(dtype=real),
        oi: wp.array2d(dtype=int),
    ):  # type: ignore[no-untyped-def]
        i = wp.tid()
        v, st, miss = IM.mcs_variance_table(
            ion, pro, comp, quad, tab, 0, z, mass, pv[i], pv_end[i], p1[i], birth[i],
            real(INV_XS), real(RHO), real(S_MM),
        )  # fmt: skip
        peq, p1eq, st2, miss2 = IM.mapped_state_table(
            ion, pro, comp, quad, tab, 0, z, mass, pv[i], p1[i]
        )
        pi, p1i, st3 = IM.interpolated_state(tab, 0, pv[i], p1[i])
        out[i, 0] = v
        out[i, 1] = peq
        out[i, 2] = p1eq
        out[i, 3] = pi
        out[i, 4] = p1i
        oi[i, 0] = st
        oi[i, 1] = st2
        oi[i, 2] = st3
        oi[i, 3] = miss + miss2

    return k


@pytest.mark.parametrize("precision", ["float64", "float32"])
@pytest.mark.parametrize("proj", IONS, ids=lambda p: p.name)
def test_table_twin_equals_warp_cpu(proj: Projectile, precision: str) -> None:
    """The table path (``mapped_state_table``, ``mcs_variance_table``, ``interpolated_state``):
    python twin against a Warp CPU kernel. float64: relative 1e-12 on all states including the
    thin-slab edge. float32 (the transport kernel evaluates the mapping in float64 for every
    precision; float32 is a twin/kernel equality check only): relative 1e-5 for the interpolated
    state on the states with ``p1v1/pv >= 1.01`` off the lowest table edge, and the accelerated
    root of the twin and the kernel agree within 1e-4 with equal status (float32 slab sums)."""
    real = wp.float64 if precision == "float64" else wp.float32
    dt = np.float64 if precision == "float64" else np.float32
    rtol = 1e-12 if precision == "float64" else 1e-5
    im = validate(_upper_case(proj, WATER)).ion_mapping
    assert im is not None
    pv, p1 = _states(im, np.random.default_rng(11), 60)
    if precision == "float32":
        assert im.table_domain is not None
        keep = (p1 / pv >= 1.01) & (pv >= im.table_domain[0] * 1.001)
        pv, p1 = pv[keep], p1[keep]
    n = pv.size
    assert n >= 30
    birth = (np.arange(n) % 3 == 0).astype(np.int32)
    pv_end = pv * np.where(birth == 1, 0.995, 1.0)
    cols = [np.asarray(c, dtype=dt) for c in (pv, pv_end, p1)]
    arrs = [np.asarray(a[:1], dtype=dt) for a in (im.ion, im.proton, im.comp, im.table)]
    quad = np.asarray(im.quad, dtype=dt)
    out = wp.zeros((n, 5), dtype=real, device="cpu")
    oi = wp.zeros((n, 4), dtype=int, device="cpu")
    wp.launch(
        _twin_kernel(real), dim=n,
        inputs=[*[wp.array(a, dtype=real, device="cpu") for a in arrs[:3]],
                wp.array(quad, dtype=real, device="cpu"),
                wp.array(arrs[3], dtype=real, device="cpu"), real(proj.z),
                real(proj.mass_mev), *[wp.array(c, dtype=real, device="cpu") for c in cols],
                wp.array(birth, dtype=int, device="cpu"), out, oi],
        device="cpu",
    )  # fmt: skip
    got, status = out.numpy().astype(np.float64), oi.numpy()
    tw = python_twin(make_ion_mcs)
    ion, pro, comp, tab = ([a[0].astype(np.float64).tolist()] for a in arrs)
    ql = quad.astype(np.float64).tolist()
    z, mass = float(proj.z), float(dt(proj.mass_mev))
    for i in range(n):
        a, b, c = (float(x[i]) for x in cols)
        pi, p1i, s3 = tw.interpolated_state(tab, 0, a, c)
        assert s3 == status[i, 2] == 0
        for x, y in ((pi, got[i, 3]), (p1i, got[i, 4])):
            assert abs(x - y) <= rtol * abs(x), (precision, proj.name, i, x, y)
        if precision == "float32":
            continue
        v, s1, m1 = tw.mcs_variance_table(
            ion, pro, comp, ql, tab, 0, z, mass, a, b, c, int(birth[i]), float(dt(INV_XS)),
            float(dt(RHO)), float(dt(S_MM)),
        )  # fmt: skip
        peq, p1eq, s2, m2 = tw.mapped_state_table(ion, pro, comp, ql, tab, 0, z, mass, a, c)
        assert s1 == s2 == status[i, 0] == status[i, 1] == 0, (i, s1, s2, status[i])
        assert m1 + m2 == status[i, 3]
        for x, y in ((v, got[i, 0]), (peq, got[i, 1]), (p1eq, got[i, 2])):
            assert abs(x - y) <= rtol * abs(x), (precision, proj.name, i, x, y)


def test_option_validation_and_summary() -> None:
    """``mcs_mapping_eval`` fails closed on unknown values; the summary records the evaluation,
    the table sizes, axes bounds, domain, sha256 and the evaluation contract with the build-time
    interpolation estimate and the bracket widths; ``"root"`` carries no table; protons ignore
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


def test_estimate_sets_the_starting_bracket() -> None:
    """The build-time interpolation estimate (3 x 3 root-solved stencil per solved cell, safety
    factor 2) sets the row's starting bracket half-width max(10 x estimate, 1e-5) (header entry
    8); constants as documented."""
    assert XO.TABLE_HEADER == 9
    assert XO.TABLE_BRACKET_FACTOR == 10.0 and XO.TABLE_BRACKET_MIN == 1.0e-5
    assert XO.TABLE_SAFETY_FACTOR == 2.0
    for case in (_upper_case, _lower_case):
        im = validate(case(CARBON12, WATER)).ion_mapping
        assert im is not None and im.estimates is not None and im.table_domain is not None
        e = im.estimates[0]
        h = im.table[0]
        xs = np.linspace(h[0], h[6], int(h[2]))
        ys = np.linspace(h[3], h[7], int(h[5]))
        assert e.n_cells == int(XO._band_cells(xs, ys, *im.table_domain[1:]).sum()) > 100
        assert e.n_cells < e.n_refined_points <= 9 * e.n_cells
        assert 0.0 < e.interp_estimate < 1e-2
        assert e.bracket_half_width == max(10.0 * e.interp_estimate, 1e-5) == h[8]
        print(
            f"carbon12 water {case.__name__[1:6]}: interpolation estimate "
            f"{e.interp_estimate:.2e}, bracket half-width {e.bracket_half_width:.2e}"
        )


def test_summary_and_capabilities_carry_the_contract() -> None:
    """The summary's ``mcs_mapping.table`` block records the evaluation contract (exact root on
    every step, the table only seeds the bracket), the solver tolerances, the bracket constants
    and widths and the interpolation estimate (overall and per material) for a two-material
    geometry; the capability report states the same contract."""
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
    assert im is not None and im.estimates is not None and len(im.estimates) == 2
    es = im.estimates
    t = eff.summary()["mcs_mapping"]["table"]
    assert t["contract"] == XO.TABLE_EVAL_CONTRACT
    for key in ("exact root solve on every step", "only supplies the starting bracket",
                "performance diagnostic", "not a correctness claim"):  # fmt: skip
        assert key in t["contract"]
    assert t["root_residual_tol"] == XO.RESIDUAL_TOL["float64"]
    assert t["root_step_tol"] == XO.STEP_TOL["float64"]
    assert (t["bracket_factor"], t["bracket_min"]) == (10.0, 1e-5)
    assert t["bracket_half_width_per_material"] == [e.bracket_half_width for e in es]
    assert t["interp_estimate_per_material"] == [e.interp_estimate for e in es]
    assert t["interp_estimate"] == max(e.interp_estimate for e in es)
    assert t["estimated_cells"] == sum(e.n_cells for e in es)
    assert t["refined_points"] == sum(e.n_refined_points for e in es)
    assert "certified_bound_analytic" not in t and "variance_bound" not in t
    cap = ionmc.capabilities()["species_details"]
    for proj in IONS:
        mcs = cap[proj.name]["mcs"]
        assert mcs["table_evaluation"]["contract"] == XO.TABLE_EVAL_CONTRACT
        assert "exact root solve on every step" in mcs["evaluation"]
        assert "mcs_table_bracket_miss" in mcs["table_evaluation"]["diagnostic_counter"]
