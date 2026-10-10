"""Axisymmetric (r, z) annular scoring grid (V3-005C C8): geometry functions, python vs Warp CPU
float64 parity, and neutrality for the Cartesian grids."""

# ruff: noqa: E501

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import replace

import numpy as np
import pytest

from ionmc._wpfunc import python_twin
from ionmc.config import SimulationConfig
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.scoring import ScoringGrid, TallyRequest, annular_cell_volume_mm3, voxel_mass_g
from ionmc.simulation import Simulation
from ionmc.transport.funcs import make_transport_funcs
from ionmc.transport.reference import run_reference_range
from ionmc.transport.warp_driver import run_warp_range

MakeConfig = Callable[..., SimulationConfig]
ANN = ScoringGrid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (40, 1, 90), name="rz", kind="annular")
CART = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 45), name="xyz")


def test_grid_validation_and_mass() -> None:
    with pytest.raises(UnsupportedCombinationError):
        ScoringGrid((0, 0, 0), (1, 1, 1), (4, 2, 4), kind="annular")
    with pytest.raises(UnsupportedCombinationError):
        ScoringGrid((0, 0, 0), (1, 2, 1), (4, 1, 4), kind="annular")
    with pytest.raises(UnsupportedCombinationError):
        ScoringGrid((0, 0, 0), (1, 1, 1), (4, 1, 4), kind="polar")  # type: ignore[arg-type]
    geo = BoxPhantom((-50.0, -50.0, 0.0), (100.0, 100.0, 100.0), WATER)
    m = voxel_mass_g(ANN, geo.to_geometry())
    assert m.shape == ANN.shape
    assert math.isclose(float(m.sum()), math.pi * 40.0**2 * 90.0 * 1e-3, rel_tol=1e-12)
    assert math.isclose(
        float(annular_cell_volume_mm3(ANN).sum()) * 1e-3, float(m.sum()), rel_tol=1e-12
    )
    small = BoxPhantom((-20.0, -20.0, 0.0), (40.0, 40.0, 100.0), WATER)  # r = 40 mm leaves it
    with pytest.raises(UnsupportedCombinationError):
        voxel_mass_g(ANN, small.to_geometry())


def test_cylinder_walk_matches_dense_sampling() -> None:
    F = python_twin(make_transport_funcs)
    V = F.vec3
    rng = np.random.default_rng(7)
    origin, spacing = V(0.3, -0.2, 0.0), V(1.0, 1.0, 1.0)
    inv = V(1.0, 1.0, 1.0)
    nr, nz = 12, 20
    for _ in range(60):
        p = np.array([rng.uniform(-6, 6), rng.uniform(-6, 6), rng.uniform(0.5, 4)])
        u = rng.normal(size=3)
        if _ % 5 == 0:
            u[2] = 0.0  # lateral grazing chords
        if _ % 7 == 0:
            u[:2] *= 1e-4  # near axial
        u /= np.linalg.norm(u)
        length = rng.uniform(5, 25)
        # walk
        ir, iz, _ = F.cyl_index(V(*p), origin, inv, nr, nz)
        cur, rem, walked = p.copy(), length, {}
        for _it in range(400):
            if rem <= 0:
                break
            pc, code = F.seg_piece_cyl(V(*cur), V(*u), ir, iz, origin, spacing, rem)
            walked[(ir, iz)] = walked.get((ir, iz), 0.0) + float(pc)
            cur = cur + u * float(pc)
            rem -= float(pc)
            if code == 0:
                ir += 1
            elif code == 1:
                ir -= 1
            elif code == 2:
                iz += 1
            elif code == 3:
                iz -= 1
        assert rem <= 1e-9
        t = (np.arange(400000) + 0.5) * (length / 400000)
        pts = p[None, :] + t[:, None] * u[None, :]
        rr = np.hypot(pts[:, 0] - 0.3, pts[:, 1] + 0.2)
        cells = np.stack([np.floor(rr).astype(int), np.floor(pts[:, 2]).astype(int)], axis=1)
        dense: dict[tuple[int, int], float] = {}
        for c, n in zip(*np.unique(cells, axis=0, return_counts=True), strict=True):
            dense[(int(c[0]), int(c[1]))] = n * length / 400000
        assert math.isclose(sum(walked.values()), length, rel_tol=1e-12)
        for k in set(dense) | {k for k, v in walked.items() if v > 1e-6}:
            assert abs(dense.get(k, 0.0) - walked.get(k, 0.0)) < 2e-3 * max(length, 1.0), k


def _cfg(make_config: MakeConfig, grids: tuple[ScoringGrid, ...], n: int = 24) -> SimulationConfig:
    cfg = make_config(
        energy=100.0, n=n, n_batches=4, lateral_sigma=3.0, energy_sigma=0.5, scoring=grids,
        geometry=BoxPhantom((-50.0, -50.0, 0.0), (100.0, 100.0, 100.0), WATER), max_step=1.0,
    )  # fmt: skip
    return replace(
        cfg,
        source=replace(cfg.source, angular_sigma_mrad=20.0, x_theta_correlation=-0.5),
    )


def test_python_vs_warp_parity_and_cartesian_neutral(make_config: MakeConfig) -> None:
    both = _cfg(make_config, (CART, ANN))
    eff = Simulation(both).effective
    py = run_reference_range(eff, 0, both.run.n_histories)
    wcfg = replace(both, run=replace(both.run, backend="warp-cpu", precision="float64"))
    wc = run_warp_range(Simulation(wcfg).effective, 0, both.run.n_histories, "cpu")
    a_py, a_wc = np.asarray(py.edep[1]), np.asarray(wc.edep[1])
    dmax = float(np.abs(a_py - a_wc).max()) * 2.0**-30
    nonzero = a_py != 0
    print(
        f"annular parity python vs warp-cpu f64: max|d edep|={dmax:.3e} MeV, cells {nonzero.sum()}"
    )
    assert dmax <= 1e-6  # fixed-point quanta (2^-30 MeV); identical trajectories
    assert (
        np.array_equal(np.asarray(py.edep[0]), np.asarray(wc.edep[0]))
        or float(np.abs(np.asarray(py.edep[0]) - np.asarray(wc.edep[0])).max()) * 2.0**-30 <= 1e-6
    )
    # neutrality: the Cartesian grid is bit-identical with and without the annular grid
    only = _cfg(make_config, (CART,))
    p0 = run_reference_range(Simulation(only).effective, 0, only.run.n_histories)
    assert np.array_equal(np.asarray(p0.edep[0]), np.asarray(py.edep[0]))
    # energy: annular total (in grid + outside) equals the Cartesian one
    res = Simulation(both).run()
    eb = res.energy_balance
    tot = [i + o for i, o in zip(eb.in_grid_mev, eb.outside_mev, strict=True)]
    assert math.isclose(tot[0], tot[1], rel_tol=1e-9)
    # the radial profile of the beam lies inside r < 40 mm: annular in-grid fraction ~ 1
    assert eb.in_grid_mev[1] > 0.99 * tot[1]
    # annulus-cell total equals the Cartesian in-grid total when both contain the deposit
    assert math.isclose(eb.in_grid_mev[1], eb.in_grid_mev[0], rel_tol=2e-2)


def test_tally_requests_on_annular_grid_are_refused(make_config: MakeConfig) -> None:
    cfg = _cfg(make_config, (ANN,))
    cfg = replace(cfg, tallies=(TallyRequest("e", "rz", "edep"),))
    with pytest.raises(UnsupportedCombinationError):
        Simulation(cfg)


def test_float32_warp_runs_and_agrees(make_config: MakeConfig) -> None:
    cfg = _cfg(make_config, (ANN,), n=64)
    py = Simulation(cfg).run().grid("rz")
    w32 = (
        Simulation(replace(cfg, run=replace(cfg.run, backend="warp-cpu", precision="float32")))
        .run()
        .grid("rz")
    )
    assert math.isclose(float(py.energy_mev.sum()), float(w32.energy_mev.sum()), rel_tol=2e-3)


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_u1_new_funcs_twin_equals_warp_cpu(precision: str) -> None:
    """The python twins of ``source_direction``, ``cyl_index``, ``cyl_next`` and ``seg_piece_cyl`` equal the
    Warp CPU functions on random inputs (the harness of test_shared_funcs lists them as exercised here)."""
    import warp as wp

    real = wp.float64 if precision == "float64" else wp.float32
    npr = np.float64 if precision == "float64" else np.float32
    tf = make_transport_funcs(real)
    v3 = tf.vec3
    F = python_twin(make_transport_funcs)

    @wp.kernel(module="unique")
    def kernel(
        p: wp.array(dtype=v3),
        u: wp.array(dtype=v3),
        z: wp.array2d(dtype=real),
        idx: wp.array2d(dtype=int),
        o: wp.array2d(dtype=real),
        oi: wp.array2d(dtype=int),
        ov: wp.array(dtype=v3),
    ) -> None:  # type: ignore[valid-type]
        i = wp.tid()
        org = v3(real(0.3), real(-0.2), real(0.0))
        spc = v3(real(1.0), real(1.0), real(1.0))
        d = u[i]
        e1, e2 = tf.orthonormal_basis(d)
        ov[i] = tf.source_direction(d, e1, e2, z[i, 0], z[i, 1], z[i, 2], z[i, 3], z[i, 4], z[i, 5])
        ir, iz, ins = tf.cyl_index(p[i], org, spc, 12, 20)
        oi[i, 0] = ir
        oi[i, 1] = iz
        oi[i, 2] = ins
        dist, code = tf.cyl_next(p[i], u[i], idx[i, 0], idx[i, 1], org, spc)
        o[i, 0] = dist
        oi[i, 3] = code
        pc, c2 = tf.seg_piece_cyl(p[i], u[i], idx[i, 0], idx[i, 1], org, spc, z[i, 6])
        o[i, 1] = pc
        oi[i, 4] = c2

    rng = np.random.default_rng(11)
    n = 200
    pp = rng.uniform([-5, -5, 0.2], [5, 5, 15], size=(n, 3)).astype(npr)
    uu = rng.normal(size=(n, 3))
    uu[::7, 2] = 0.0
    uu /= np.linalg.norm(uu, axis=1, keepdims=True)
    uu = uu.astype(npr)
    zz = np.concatenate([rng.normal(size=(n, 4)), np.full((n, 1), 0.03), np.full((n, 1), -0.5),
                         rng.uniform(0.1, 6.0, size=(n, 1))], axis=1)  # fmt: skip
    zz[:, 4] = 0.03
    zz[:, 5] = -0.5
    zz = zz.astype(npr)
    ia = np.array([[int(np.hypot(a[0] - 0.3, a[1] + 0.2)), int(a[2])] for a in pp], dtype=np.int32)
    o = wp.zeros((n, 2), dtype=real)
    oi = wp.zeros((n, 5), dtype=int)
    ov = wp.zeros(n, dtype=v3)
    wp.launch(kernel, dim=n, inputs=[wp.array(pp, dtype=v3), wp.array(uu, dtype=v3), wp.array(zz, dtype=real),
                                     wp.array(ia, dtype=int), o, oi, ov], device="cpu")  # fmt: skip
    o_n, oi_n, ov_n = o.numpy(), oi.numpy(), ov.numpy()
    rtol = 1e-12 if precision == "float64" else 2e-5
    V = F.vec3
    org, spc = V(0.3, -0.2, 0.0), V(1.0, 1.0, 1.0)
    for i in range(n):
        d = V(*map(float, uu[i]))
        e1, e2 = F.orthonormal_basis(d)
        dv = F.source_direction(d, e1, e2, *map(float, zz[i, :6]))
        assert np.allclose([dv[0], dv[1], dv[2]], ov_n[i], rtol=rtol, atol=rtol)
        p = V(*map(float, pp[i]))
        ir, iz, ins = F.cyl_index(p, org, spc, 12, 20)
        assert (ir, iz, ins) == tuple(oi_n[i, :3])
        dist, code = F.cyl_next(p, d, int(ia[i, 0]), int(ia[i, 1]), org, spc)
        assert math.isclose(dist, o_n[i, 0], rel_tol=rtol, abs_tol=rtol) and code == oi_n[i, 3]
        pc, c2 = F.seg_piece_cyl(p, d, int(ia[i, 0]), int(ia[i, 1]), org, spc, float(zz[i, 6]))
        assert math.isclose(pc, o_n[i, 1], rel_tol=rtol, abs_tol=rtol) and c2 == oi_n[i, 4]
