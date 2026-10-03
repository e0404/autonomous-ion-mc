"""U1: every shared physics/transport function, Python scope versus a Warp CPU kernel.

Both paths evaluate the same ``@wp.func`` on identical arguments (including edge cases):
Python scope with ``wp.float64``/``wp.float32`` objects, and a compiled kernel. Tolerances are
the decision 0001 classes (float64 rtol 1e-12 / atol 1e-14, float32 rtol 1e-6 / atol 1e-6);
a NaN in either path fails; integer outputs (voxel indices, step reasons, accept flags) must
be identical.
"""

import math
import os
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
import warp as wp

from ionmc.physics.em import make_em
from ionmc.physics.kinematics import make_kinematics
from ionmc.transport.funcs import make_transport_funcs

wp.config.log_level = wp.LOG_WARNING

# Frozen criterion U1 asks for 1e4 arguments per function. Python-scope Warp calls cost about
# 3 ms per argument set, so the default is 2000 (about 12 s for both precisions); set the
# environment variable IONMC_U1_N=10000 for the full count (about 60 s).
N_ARGS = int(os.environ.get("IONMC_U1_N", "2000"))
TOL = {"float64": (1e-12, 1e-14), "float32": (1e-6, 1e-6)}
M_P = 938.27208816

# column layout of the argument pools -------------------------------------------------------
NX = 56  # real columns
NV = 16  # vec3 columns
NI = 9  # int columns
NO = 40  # real outputs
NOI = 16  # int outputs
NOV = 8  # vec3 outputs


def _unit(rng: np.random.Generator, n: int) -> np.ndarray:
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    # edge cases: axes, poles, zero components, diagonals
    edge = np.array(
        [
            [0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
            [1, 1, 0], [1, 1, 1], [-1, 1, -1], [0, 1, 1], [1e-8, 0, 1], [1e-4, 1e-4, 1],
        ],
        dtype=np.float64,
    )  # fmt: skip
    edge /= np.linalg.norm(edge, axis=1)[:, None]
    v[: len(edge)] = edge
    return v


def _args(precision: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Deterministic argument pools (identical for both paths), values cast to ``precision``."""
    rng = np.random.default_rng(20261003)
    n = N_ARGS
    dt = np.float64 if precision == "float64" else np.float32
    x = np.zeros((n, NX))
    v = np.zeros((n, NV, 3))
    ii = np.zeros((n, NI), dtype=np.int32)

    def logu(lo: float, hi: float) -> np.ndarray:
        return np.exp(rng.uniform(math.log(lo), math.log(hi), n))

    def uni() -> np.ndarray:  # a uniform of the Philox mapping of the precision
        w = rng.integers(0, 2**32, n, dtype=np.uint64)
        if precision == "float64":
            return ((w >> np.uint64(8)).astype(np.float64) + 0.5) * 2.0**-24
        return ((w >> np.uint64(9)).astype(np.float64) + 0.5) * 2.0**-23

    x[:, 0] = logu(1e-3, 1e3)  # kinetic energy
    x[:, 1] = rng.choice([M_P, 1875.61294257, 3727.3794066, 0.51099895], n)
    # csda_mean_loss: e0, e_r1, s0, tt, r0, f_short
    x[:, 2] = rng.uniform(2.0, 300.0, n)
    x[:, 3] = x[:, 2] * rng.uniform(0.0, 1.0, n)
    x[:, 4] = rng.uniform(1.0, 400.0, n)
    x[:, 6] = logu(0.01, 50.0)
    x[:, 5] = x[:, 6] * rng.uniform(0.0, 1.5, n)
    x[:, 7] = 1e-3
    x[0:4, 5] = [0.0, x[1, 6], 1e-3 * x[2, 6], x[3, 6] * 2.0]
    # bohr_variance: e_mid, m, z, z/a, rho, s
    x[:, 8] = logu(2.0, 300.0)
    x[:, 9] = M_P
    x[:, 10] = 1.0
    x[:, 11] = 0.555
    x[:, 12] = rng.uniform(0.001, 12.0, n)
    x[:, 13] = logu(1e-3, 5.0)
    # straggle_attempt: mean, var, u0..u3
    x[:, 14] = logu(1e-3, 50.0)
    ratio = logu(0.02, 40.0)
    ratio[:8] = [3.0, 2.999999, 3.000001, 1.0, 0.1, 1e3, 1e-3, 2.0]
    x[:, 15] = (x[:, 14] / ratio) ** 2
    x[0, 14] = 0.0
    x[1, 15] = 0.0
    for c in range(16, 20):
        x[:, c] = uni()
    # fdm / scattering power
    x[:, 21] = logu(50.0, 900.0)  # p1v1
    x[:, 20] = x[:, 21] * rng.uniform(0.05, 1.05, n)
    x[0:3, 20] = x[0:3, 21]
    x[:, 22] = rng.uniform(0.01, 0.2, n)  # inv_rho_xs
    x[:, 23] = rng.uniform(0.001, 12.0, n)  # rho
    x[:, 24] = 1.0  # z
    # polar_deflection, rotate_dir
    x[:, 25] = logu(1e-12, 1e-2)
    x[:, 26] = uni()
    x[:, 27] = rng.uniform(0.0, math.pi, n)
    x[:, 28] = rng.uniform(0.0, 2.0 * math.pi, n)
    # log_bin_index: x, lx0, inv_dl ; lerp ; interp_exp
    x[:, 29] = logu(1e-2, 1e4)
    k = rng.integers(0, 560, 40)
    x[:40, 29] = np.exp(k / 86.0)  # exact grid points
    x[:, 30] = 0.0
    x[:, 31] = 86.0
    x[:, 32] = rng.uniform(-5, 5, n)
    x[:, 33] = rng.uniform(-5, 5, n)
    x[:, 34] = rng.uniform(0.0, 1.0, n)
    x[:, 35] = rng.uniform(-5, 10, n)
    x[:, 36] = rng.uniform(-5, 10, n)
    # plane_position origin, spacing
    x[:, 37] = np.round(rng.uniform(-0.5, 0.5, n) * 64.0) / 64.0
    x[:, 38] = rng.integers(6, 17, n) / 64.0
    # leg2 target
    x[:, 39] = logu(1e-3, 0.5)
    # range / eloss / select / hinge
    x[:, 40] = logu(0.01, 300.0)
    x[:, 41] = rng.uniform(0.05, 1.0, n)
    x[:, 42] = rng.uniform(0.01, 1.0, n)
    x[0, 40] = x[0, 42]
    x[:, 43] = logu(2.0, 300.0)
    x[:, 44] = logu(0.1, 100.0)
    x[:, 45] = rng.uniform(0.001, 0.2, n)
    for c in range(46, 50):
        x[:, c] = logu(1e-3, 5.0)
    x[0:4, 46:50] = [[1, 1, 1, 1], [1, 2, 1, 3], [2, 1, 1, 1], [3, 2, 1, 1]]
    x[:, 50] = logu(1e-3, 5.0)
    x[:, 51] = x[:, 50] * rng.uniform(0.0, 2.0, n)
    x[0:3, 51] = [x[0, 50], 0.0, x[2, 50] * 2.0]
    x[:, 52] = uni()
    x[:, 53] = uni()
    # ints: ix, iy, iz in the DDA grid and indices
    ii[:, 0:3] = rng.integers(-2, 9, (n, 3))
    ii[:, 3] = rng.integers(-3, 50, n)
    ii[:, 4] = rng.integers(0, 2, n)
    ii[:, 5] = 541
    ii[:, 6:9] = rng.integers(1, 30, (n, 3))
    # vectors: origin, spacing, positions on planes
    # Coordinates of order one, so that the float32 absolute tolerance is meaningful, with
    # dyadic origin and spacing so that ``index * spacing + origin`` is exact in float32
    # (a compiled kernel may fuse it into one rounding, Python scope rounds twice).
    origin = np.round(rng.uniform(-0.5, 0.5, (n, 3)) * 64.0) / 64.0
    spacing = rng.integers(6, 17, (n, 3)) / 64.0
    origin[:6] = [[0, 0, 0], [0, 0, 0], [-0.5, -0.5, -0.5], [0, 0, 0], [0.25, 0.5, 0.75], [0, 0, 0]]
    spacing[:6] = [[0.25] * 3, [0.125] * 3, [0.25] * 3, [0.5] * 3, [0.125] * 3, [0.125, 0.25, 0.5]]
    v[:, 2] = origin
    v[:, 3] = spacing
    pos = origin + spacing * (ii[:, 0:3] + rng.uniform(0, 1, (n, 3)))
    on_plane = rng.random((n, 3)) < 0.3
    pos = np.where(on_plane, origin + spacing * np.round((pos - origin) / spacing), pos)
    v[:, 0] = pos
    v[:, 1] = _unit(rng, n)
    v[:, 4] = pos
    v[:, 5] = _unit(rng, n)
    v[:, 6] = rng.uniform(-1, 1, (n, 3))
    v[:, 7] = _unit(rng, n)
    v[:, 8] = _unit(rng, n)
    v[:, 9] = rng.uniform(-1, 3, (n, 3))
    v[:9, 9] = origin[:9] + spacing[:9] * ii[:9, 6:9]  # far edge
    v[:, 10] = 1.0 / spacing
    lo = rng.uniform(-1, 0, (n, 3))
    hi = lo + rng.uniform(0.2, 2, (n, 3))
    v[:, 13] = lo
    v[:, 14] = hi
    v[:, 11] = np.where(rng.random((n, 3)) < 0.3, lo, rng.uniform(-2, 2, (n, 3)))
    v[:, 12] = _unit(rng, n)
    v[:, 15] = _unit(rng, n)
    return x.astype(dt), v.astype(dt), ii


class _Fn(SimpleNamespace):
    kin: Any
    em: Any
    tf: Any


def _funcs(real: Any) -> _Fn:
    return _Fn(kin=make_kinematics(real), em=make_em(real), tf=make_transport_funcs(real))


def _make_kernel(real: Any) -> Any:
    f = _funcs(real)
    kin, em, tf = f.kin, f.em, f.tf
    v3 = tf.vec3
    two = real(2.0)
    del two

    @wp.kernel(module="unique")
    def kernel(
        x: wp.array2d(dtype=real),  # type: ignore[valid-type]
        v: wp.array2d(dtype=v3),  # type: ignore[valid-type]
        ii: wp.array2d(dtype=int),  # type: ignore[valid-type]
        o: wp.array2d(dtype=real),  # type: ignore[valid-type]
        oi: wp.array2d(dtype=int),  # type: ignore[valid-type]
        ov: wp.array2d(dtype=v3),  # type: ignore[valid-type]
    ) -> None:
        i = wp.tid()
        o[i, 0] = kin.pv_mev(x[i, 0], x[i, 1])
        o[i, 1] = kin.beta2(x[i, 0], x[i, 1])
        o[i, 2] = kin.gamma(x[i, 0], x[i, 1])
        o[i, 3] = kin.tmax_mev(x[i, 0], x[i, 1])
        o[i, 4] = em.csda_mean_loss(x[i, 2], x[i, 3], x[i, 4], x[i, 5], x[i, 6], x[i, 7])
        o[i, 5] = em.bohr_variance(x[i, 8], x[i, 9], x[i, 10], x[i, 11], x[i, 12], x[i, 13])
        sl, sok = em.straggle_attempt(x[i, 14], x[i, 15], x[i, 16], x[i, 17], x[i, 18], x[i, 19])
        o[i, 6] = sl
        oi[i, 0] = sok
        o[i, 7] = em.fdm(x[i, 20], x[i, 21])
        o[i, 8] = em.scattering_power_dm(x[i, 20], x[i, 21], x[i, 24], x[i, 22], x[i, 23])
        o[i, 9] = em.polar_deflection(x[i, 25], x[i, 26])
        ov[i, 0] = em.rotate_dir(v[i, 15], x[i, 27], x[i, 28])
        lb, lf = tf.log_bin_index(x[i, 29], x[i, 30], x[i, 31], ii[i, 5])
        oi[i, 1] = lb
        o[i, 10] = lf
        o[i, 11] = tf.lerp(x[i, 32], x[i, 33], x[i, 34])
        o[i, 12] = tf.interp_exp(x[i, 35], x[i, 36], x[i, 34])
        o[i, 13] = tf.plane_position(ii[i, 3], ii[i, 4], x[i, 37], x[i, 38])
        dd, dax = tf.dda_next(v[i, 0], v[i, 1], ii[i, 0], ii[i, 1], ii[i, 2], v[i, 2], v[i, 3])
        o[i, 14] = dd
        oi[i, 2] = dax
        l2, l2ax = tf.leg2_limit(
            v[i, 4], v[i, 5], ii[i, 0], ii[i, 1], ii[i, 2], v[i, 2], v[i, 3], x[i, 39]
        )
        o[i, 15] = l2
        oi[i, 3] = l2ax
        o[i, 16] = tf.range_step_limit(x[i, 40], x[i, 41], x[i, 42])
        o[i, 17] = tf.eloss_step_limit(x[i, 43], x[i, 44], x[i, 45])
        ss, sr = tf.select_step(x[i, 46], x[i, 47], x[i, 48], x[i, 49])
        o[i, 18] = ss
        oi[i, 4] = sr
        ov[i, 1] = tf.point_on_hinge(v[i, 6], v[i, 7], x[i, 50], v[i, 8], x[i, 51])
        gx, gy, gz, gin = tf.grid_index(v[i, 9], v[i, 2], v[i, 10], ii[i, 6], ii[i, 7], ii[i, 8])
        oi[i, 5] = gx
        oi[i, 6] = gy
        oi[i, 7] = gz
        oi[i, 8] = gin
        rt0, rt1 = tf.ray_box(v[i, 11], v[i, 12], v[i, 13], v[i, 14])
        o[i, 19] = rt0
        o[i, 20] = rt1
        g0, g1 = tf.gauss_pair(x[i, 52], x[i, 53])
        o[i, 21] = g0
        o[i, 22] = g1
        o[i, 23] = tf.gauss_one(x[i, 52], x[i, 53])
        b1, b2 = tf.orthonormal_basis(v[i, 15])
        ov[i, 2] = b1
        ov[i, 3] = b2

    return kernel


def _python_scope(real: Any, x: np.ndarray, v: np.ndarray, ii: np.ndarray) -> tuple[Any, ...]:
    f = _funcs(real)
    kin, em, tf = f.kin, f.em, f.tf
    v3 = tf.vec3
    n = x.shape[0]
    o = np.zeros((n, NO))
    oi = np.zeros((n, NOI), dtype=np.int64)
    ov = np.zeros((n, NOV, 3))

    def vec(a: np.ndarray) -> Any:
        return v3(real(float(a[0])), real(float(a[1])), real(float(a[2])))

    for i in range(n):
        r = [real(float(c)) for c in x[i]]
        vv = [vec(c) for c in v[i]]
        k = [int(c) for c in ii[i]]
        o[i, 0] = float(kin.pv_mev(r[0], r[1]))
        o[i, 1] = float(kin.beta2(r[0], r[1]))
        o[i, 2] = float(kin.gamma(r[0], r[1]))
        o[i, 3] = float(kin.tmax_mev(r[0], r[1]))
        o[i, 4] = float(em.csda_mean_loss(r[2], r[3], r[4], r[5], r[6], r[7]))
        o[i, 5] = float(em.bohr_variance(r[8], r[9], r[10], r[11], r[12], r[13]))
        sl, sok = em.straggle_attempt(r[14], r[15], r[16], r[17], r[18], r[19])
        o[i, 6] = float(sl)
        oi[i, 0] = int(sok)
        o[i, 7] = float(em.fdm(r[20], r[21]))
        o[i, 8] = float(em.scattering_power_dm(r[20], r[21], r[24], r[22], r[23]))
        o[i, 9] = float(em.polar_deflection(r[25], r[26]))
        ov[i, 0] = [float(c) for c in em.rotate_dir(vv[15], r[27], r[28])]
        lb, lf = tf.log_bin_index(r[29], r[30], r[31], k[5])
        oi[i, 1] = int(lb)
        o[i, 10] = float(lf)
        o[i, 11] = float(tf.lerp(r[32], r[33], r[34]))
        o[i, 12] = float(tf.interp_exp(r[35], r[36], r[34]))
        o[i, 13] = float(tf.plane_position(k[3], k[4], r[37], r[38]))
        dd, dax = tf.dda_next(vv[0], vv[1], k[0], k[1], k[2], vv[2], vv[3])
        o[i, 14] = float(dd)
        oi[i, 2] = int(dax)
        l2, l2ax = tf.leg2_limit(vv[4], vv[5], k[0], k[1], k[2], vv[2], vv[3], r[39])
        o[i, 15] = float(l2)
        oi[i, 3] = int(l2ax)
        o[i, 16] = float(tf.range_step_limit(r[40], r[41], r[42]))
        o[i, 17] = float(tf.eloss_step_limit(r[43], r[44], r[45]))
        ss, sr = tf.select_step(r[46], r[47], r[48], r[49])
        o[i, 18] = float(ss)
        oi[i, 4] = int(sr)
        ov[i, 1] = [float(c) for c in tf.point_on_hinge(vv[6], vv[7], r[50], vv[8], r[51])]
        gx, gy, gz, gin = tf.grid_index(vv[9], vv[2], vv[10], k[6], k[7], k[8])
        oi[i, 5:9] = [int(gx), int(gy), int(gz), int(gin)]
        rt0, rt1 = tf.ray_box(vv[11], vv[12], vv[13], vv[14])
        o[i, 19] = float(rt0)
        o[i, 20] = float(rt1)
        g0, g1 = tf.gauss_pair(r[52], r[53])
        o[i, 21] = float(g0)
        o[i, 22] = float(g1)
        o[i, 23] = float(tf.gauss_one(r[52], r[53]))
        b1, b2 = tf.orthonormal_basis(vv[15])
        ov[i, 2] = [float(c) for c in b1]
        ov[i, 3] = [float(c) for c in b2]
    return o, oi, ov


REAL_OUTPUTS = {
    0: "pv_mev", 1: "beta2", 2: "gamma", 3: "tmax_mev", 4: "csda_mean_loss",
    5: "bohr_variance", 6: "straggle_loss", 7: "fdm", 8: "scattering_power_dm",
    9: "polar_deflection", 10: "log_bin_frac", 11: "lerp", 12: "interp_exp",
    13: "plane_position", 14: "dda_distance", 15: "leg2_length", 16: "range_step_limit",
    17: "eloss_step_limit", 18: "select_step", 19: "ray_box_in", 20: "ray_box_out",
    21: "gauss_pair_0", 22: "gauss_pair_1", 23: "gauss_one",
}  # fmt: skip
INT_OUTPUTS = {
    0: "straggle_ok", 1: "log_bin_index", 2: "dda_axis", 3: "leg2_axis", 4: "select_reason",
    5: "grid_ix", 6: "grid_iy", 7: "grid_iz", 8: "grid_inside",
}  # fmt: skip
VEC_OUTPUTS = {0: "rotate_dir", 1: "point_on_hinge", 2: "basis_e1", 3: "basis_e2"}


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_u1_python_scope_equals_warp_cpu_kernel(precision: str) -> None:
    real = wp.float64 if precision == "float64" else wp.float32
    rtol, atol = TOL[precision]
    x, v, ii = _args(precision)
    n = x.shape[0]
    kernel = _make_kernel(real)
    v3 = _funcs(real).tf.vec3
    dev = "cpu"
    ox = wp.zeros((n, NO), dtype=real, device=dev)
    oi = wp.zeros((n, NOI), dtype=int, device=dev)
    ov = wp.zeros((n, NOV), dtype=v3, device=dev)
    wp.launch(
        kernel,
        dim=n,
        inputs=[
            wp.array(x, dtype=real, device=dev),
            wp.array(v, dtype=v3, device=dev),
            wp.array(ii, dtype=int, device=dev),
            ox,
            oi,
            ov,
        ],
        device=dev,
    )
    k_o, k_oi, k_ov = ox.numpy().astype(np.float64), oi.numpy(), ov.numpy().astype(np.float64)
    p_o, p_oi, p_ov = _python_scope(real, x, v, ii)
    for col, name in REAL_OUTPUTS.items():
        a, b = p_o[:, col], k_o[:, col]
        assert np.all(np.isfinite(a)) and np.all(np.isfinite(b)), f"{name}: NaN or inf"
        bad = ~np.isclose(a, b, rtol=rtol, atol=atol)
        assert not bad.any(), (
            f"{name}: {bad.sum()} of {n} mismatches; first python={a[bad][:3]} kernel={b[bad][:3]}"
        )
    for col, name in INT_OUTPUTS.items():
        assert np.array_equal(p_oi[:, col], k_oi[:, col]), f"{name}: integer outputs differ"
    for col, name in VEC_OUTPUTS.items():
        a, b = p_ov[:, col], k_ov[:, col]
        assert np.all(np.isfinite(a)) and np.all(np.isfinite(b)), f"{name}: NaN or inf"
        assert np.allclose(a, b, rtol=rtol, atol=atol), f"{name}: components differ"
