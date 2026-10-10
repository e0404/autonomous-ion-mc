"""U1: every shared physics/transport function, Python twin versus a Warp CPU kernel.

Warp functions are never called at Python scope; the Python side is the float64 pure-Python twin
(``python_twin``, the reference backend's functions: the same source text) evaluated on identical
arguments (including edge cases), the kernel side a compiled float64 or float32 kernel. See
``TOL`` for the tolerances (float64: decision 0001 class; float32: float32 rounding against the
float64 twin). A NaN in either path fails; integer outputs (voxel indices, step reasons, accept
flags) must be identical, except the discontinuous float32 log-bin index (see the test).
"""

import math
import os
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
import warp as wp

from ionmc._wpfunc import python_twin
from ionmc.physics.em import FDM_COEFFICIENTS, make_em
from ionmc.physics.kinematics import make_kinematics
from ionmc.physics.nuclear import make_nuclear
from ionmc.transport.funcs import make_transport_funcs
from ionmc.transport.scoring_funcs import make_scoring_funcs

wp.config.log_level = wp.LOG_WARNING

# Frozen criterion U1: 1e4 arguments per function (the environment variable only allows a
# smaller count for quick local runs; CI uses the default).
N_ARGS = int(os.environ.get("IONMC_U1_N", "10000"))
TOL = {"float64": (1e-12, 1e-14), "float32": (1e-6, 1e-6)}
"""(rtol, atol). float64: twin against kernel, the same double arithmetic (decision 0001 class).
float32: the float32 kernel against the FLOAT64 twin evaluated on the float32-rounded arguments
(plan amendment 18): the default is the frozen 1e-6 / 1e-6 (16.8 u, u = 2^-24 = 5.96e-8, so a
function of up to about 16 roundings without ill-conditioning meets it); the functions whose
conditioning is worse have the analytic budgets of ``_float32_budgets`` (derived from the
condition number times u times the number of roundings, not from observed differences)."""
U = 2.0**-24
LN10 = math.log(10.0)
M_P = 938.27208816

# column layout of the argument pools -------------------------------------------------------
NX = 87  # real columns
NV = 16  # vec3 columns
NI = 13  # int columns
NO = 40  # real outputs
NOI = 17  # int outputs
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
    # scattering_variance_birth: pv_mid, pv_end, p1v1, z, inv_rho_xs, rho, s
    x[:, 56] = logu(50.0, 900.0)  # p1v1
    r_end = rng.uniform(0.02, 0.97, n)
    x[:, 55] = x[:, 56] * r_end  # pv_end
    x[:, 54] = x[:, 55] + (x[:, 56] - x[:, 55]) * rng.uniform(0.0, 0.5, n)  # pv_mid
    x[:, 57] = 1.0
    x[:, 58] = rng.uniform(0.01, 0.2, n)
    x[:, 59] = rng.uniform(0.001, 12.0, n)
    x[:, 60] = logu(1e-3, 5.0)
    # edges: E1 == E0 (no loss), clamp active (pv_end within 1e-9 of p1v1, beyond it), tiny step,
    # near the end of the range (low pv, large p1v1)
    x[0, 54:57] = [400.0, 400.0, 400.0]
    x[1, 54:57] = [399.9999, 400.0 * (1 - 1e-9), 400.0]
    x[2, 54:57] = [400.0, 400.0 * 1.01, 400.0]
    x[3, 54:57] = [12.0, 10.0, 600.0]
    x[4, 54:57] = [399.99, 399.98, 400.0]
    x[4, 60] = 1e-6
    x[5, 54:57] = [30.0, 25.0, 30.0]
    x[6, 54:57] = [500.0, 480.0, 550.0]  # ordinary step with the clamp inactive
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
    # clipped DDA, clipped leg and scoring piece: the clip plane z_clip (column 61), the length of
    # the segment still to be walked (column 62). Random rows: the clip plane behind, ahead of
    # and far beyond the position (30 % at BIG_LENGTH), the remaining length from 1e-3 to 2.
    big = 1.0e30
    sz = spacing[:, 2]
    zc = pos[:, 2] + rng.uniform(-0.3, 1.5, n) * sz * rng.choice([-1.0, 1.0], n)
    x[:, 61] = np.where(rng.random(n) < 0.3, big, zc)
    x[:, 62] = logu(1e-3, 2.0)
    # Edge rows 100..499 (exact in float32 and float64, so that ties are the same in both): the
    # position on a dyadic point of the voxel, the direction an axis, so that the distance to the
    # next voxel plane d is exact. Remaining length 0, d (the plane does NOT end the piece), d / 2,
    # 1.5 d; leg target d (the plane at the leg end), d / 2, 1.5 d; clip plane absent, exactly on
    # the next voxel plane (tie: the clip plane wins), before it, after it, ignored (u_z = 0).
    for r in range(100, 500):
        a = r % 3
        sgn = 1.0 if (r // 3) % 2 == 0 else -1.0
        pe = origin[r] + spacing[r] * (ii[r, 0:3] + ((r // 6) % 9) / 8.0)
        ue = np.zeros(3)
        ue[a] = sgn
        v[r, 0] = v[r, 4] = pe
        v[r, 1] = v[r, 5] = ue
        up = 1 if sgn > 0 else 0
        plane = (ii[r, a] + up) * spacing[r, a] + origin[r, a]
        d = abs(plane - pe[a])
        x[r, 62] = d * (0.0, 1.0, 0.5, 1.5)[(r // 4) % 4]
        x[r, 39] = max(d, 1.0 / 64.0) * (1.0, 0.5, 1.5)[(r // 5) % 3]
        if a == 2:
            x[r, 61] = (
                big,
                plane,
                pe[2] + 0.5 * (plane - pe[2]),
                plane + 0.25 * sgn,
            )[(r // 7) % 4]
        else:
            x[r, 61] = pe[2] + 0.3  # u_z = 0: the clip plane is ignored
    # scoring functions (V3-004): columns 63..80 of x, 9..12 of ii
    x[:, 63] = rng.uniform(-5, 5, n)  # ly0
    x[:, 64] = rng.uniform(-5, 5, n)  # ly1
    x[:, 65] = 86.0  # inv_dl
    x[:, 66] = logu(0.05, 30.0)  # s_mid
    x[:, 67] = rng.uniform(-1.2, 0.3, n)  # gamma
    x[:, 68] = logu(1e-4, 5.0)  # de_mean
    x[:, 69] = logu(1.0, 300.0)  # e_mid
    x[:, 70] = logu(1e-3, 2.0)  # s_act
    x[0:3, 70] = [0.0, 0.0, 1e-3]
    x[1, 69] = 0.0  # e_mid <= 0 -> k = 0
    x[:, 71] = rng.uniform(-0.1, 0.1, n)  # k
    x[4:8, 71] = 0.0
    x[:, 72] = logu(0.05, 5.0)  # e_dot
    x[:, 73] = rng.uniform(-1.0, 1.0, n)  # tau
    x[0:6, 73] = [0.0, 0.5, -0.5, 1.0, -1.0, 0.25]
    x[:, 74] = logu(1e-3, 2.0)  # length
    x[0, 74] = 0.0
    x[:, 75] = logu(1e-4, 5.0)  # eps
    x[:, 76] = logu(0.05, 30.0)  # s_bar
    x[:, 77] = rng.uniform(0.0, 20.0, n)  # f
    ii[:, 9] = rng.integers(-1, 9, n)  # kind, including invalid codes -1 and 8
    ii[:, 9][:10] = np.arange(-1, 9)
    ii[:, 10] = rng.integers(2, 41, n)  # axis points n >= 2 (n = 2 included)
    ii[:, 11] = rng.integers(0, 2, n)  # log axis
    ii[:, 12] = rng.integers(1, 41, n)  # spectrum bins
    ii[0:4, 10] = 2
    a0 = rng.uniform(-1.0, 1.0, n)
    inv_da = rng.uniform(0.5, 2.0, n)
    t_target = rng.uniform(-0.6, ii[:, 10] - 0.4)
    arg = a0 + t_target / inv_da
    # edge rows: dyadic axis (a0 = 0, inv_da = 1), x exactly on, below and above the edges
    a0[:12], inv_da[:12] = 0.0, 1.0
    arg[:12] = [0.0, 4.0, -0.5, 4.5, 2.5, 1.0, 3.0, 0.25, 5.0, 4.999, 2.0, -1e-3]
    ii[:12, 10] = 5
    ii[:12, 12] = 4
    ii[:12, 11] = 0
    x[:, 78] = np.where(ii[:, 11] == 1, np.exp(arg), arg)  # lookup / spectrum argument
    x[:, 79] = a0
    x[:, 80] = inv_da
    # range_in_bin: r_i, f_i, d_i, h, phi (|d phi| spans 0 to 0.2 and crosses the series guard 1e-5)
    x[:, 81] = logu(1e-3, 200.0)
    x[:, 82] = logu(0.05, 300.0)
    x[:, 83] = logu(1e-9, 0.2) * rng.choice([-1.0, 1.0], n)
    x[:, 84] = logu(5e-3, 2e-2)
    x[:, 85] = rng.uniform(0.0, 1.0, n)
    x[0:8, 85] = [0.0, 1.0, 1.0, 0.5, 1.0, 1.0, 1e-3, 1.0]
    x[0:8, 83] = [0.1, 0.1, 0.0, -0.1, 1.0e-5, 0.99999e-5, 1.0e-5, -1.0e-5]
    x[8:12, 85] = 1.0
    x[8:12, 83] = [1.0e-4, -1.0e-4, 1.0e-6, -3.0e-5]
    # select_step_nuclear (V3-005A): column 86 is the nuclear limit d_nuc.
    # Rows 0..3 tie the nuclear limit with each of the EM limits (EM wins) and the
    # geometry (geometry wins); rows 4..5 are below every limit and BIG.
    x[:, 86] = logu(1e-3, 5.0)
    x[0:4, 86] = [x[0, 47], x[1, 48], x[2, 49], x[3, 46]]
    x[4, 86] = 0.5 * min(x[4, 46:50])
    x[5, 86] = 1.0e30
    x[6, 86] = x[6, 46]
    return x.astype(dt), v.astype(dt), ii


class _Fn(SimpleNamespace):
    kin: Any
    em: Any
    tf: Any
    sc: Any


def _funcs(real: Any) -> _Fn:
    return _Fn(
        kin=make_kinematics(real),
        em=make_em(real),
        tf=make_transport_funcs(real),
        sc=make_scoring_funcs(real),
    )


def _make_kernel(real: Any) -> Any:
    f = _funcs(real)
    kin, em, tf, sc = f.kin, f.em, f.tf, f.sc
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
        o[i, 24] = em.scattering_variance_birth(
            x[i, 54], x[i, 55], x[i, 56], x[i, 57], x[i, 58], x[i, 59], x[i, 60]
        )
        b1, b2 = tf.orthonormal_basis(v[i, 15])
        ov[i, 2] = b1
        ov[i, 3] = b2
        cd, cax = tf.dda_next_clip(
            v[i, 0], v[i, 1], ii[i, 0], ii[i, 1], ii[i, 2], v[i, 2], v[i, 3], x[i, 61]
        )
        o[i, 25] = cd
        oi[i, 9] = cax
        cl, clax = tf.leg2_limit_clip(
            v[i, 4], v[i, 5], ii[i, 0], ii[i, 1], ii[i, 2], v[i, 2], v[i, 3], x[i, 39], x[i, 61]
        )
        o[i, 26] = cl
        oi[i, 10] = clax
        sp, spax = tf.seg_piece(
            v[i, 0], v[i, 1], ii[i, 0], ii[i, 1], ii[i, 2], v[i, 2], v[i, 3], x[i, 62]
        )
        o[i, 27] = sp
        oi[i, 11] = spax
        gl, gok = em.straggle_attempt_gamma(
            x[i, 14], x[i, 15], x[i, 16], x[i, 17], x[i, 18], x[i, 19]
        )
        o[i, 28] = gl
        oi[i, 12] = gok
        o[i, 29] = sc.loglog_slope(x[i, 63], x[i, 64], x[i, 65])
        o[i, 30] = sc.let_ramp_slope(x[i, 66], x[i, 67], x[i, 68], x[i, 69], x[i, 70])
        ps, pe = sc.piece_state(x[i, 66], x[i, 71], x[i, 69], x[i, 72], x[i, 73])
        o[i, 31] = ps
        o[i, 32] = pe
        pm1, pm2 = sc.piece_moments(x[i, 76], x[i, 71], x[i, 74])
        o[i, 33] = pm1
        o[i, 34] = pm2
        o[i, 35] = sc.channel_value(ii[i, 9], x[i, 75], x[i, 74], pm1, pm2, x[i, 76], x[i, 77])
        lbi, lbf, lbin = sc.lookup_bin(x[i, 78], x[i, 79], x[i, 80], ii[i, 10], ii[i, 11])
        oi[i, 13] = lbi
        o[i, 36] = lbf
        o[i, 37] = tf.range_in_bin(x[i, 81], x[i, 82], x[i, 83], x[i, 84], x[i, 85])
        oi[i, 14] = lbin
        oi[i, 15] = sc.spectrum_bin(x[i, 78], x[i, 79], x[i, 80], ii[i, 12], ii[i, 11])
        sn, snr = tf.select_step_nuclear(x[i, 46], x[i, 47], x[i, 48], x[i, 49], x[i, 86])
        o[i, 38] = sn
        oi[i, 16] = snr

    return kernel


def _python_scope(real: Any, x: np.ndarray, v: np.ndarray, ii: np.ndarray) -> tuple[Any, ...]:
    """Evaluate every function with the pure-Python float64 twins used by the reference backend
    (same source text as the Warp functions, no Warp call). Warp functions run only in kernels."""
    f = _Fn(
        kin=python_twin(make_kinematics),
        em=python_twin(make_em),
        tf=python_twin(make_transport_funcs),
        sc=python_twin(make_scoring_funcs),
    )
    real = float
    kin, em, tf, sc = f.kin, f.em, f.tf, f.sc
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
        o[i, 24] = float(
            em.scattering_variance_birth(r[54], r[55], r[56], r[57], r[58], r[59], r[60])
        )
        b1, b2 = tf.orthonormal_basis(vv[15])
        ov[i, 2] = [float(c) for c in b1]
        ov[i, 3] = [float(c) for c in b2]
        cd, cax = tf.dda_next_clip(vv[0], vv[1], k[0], k[1], k[2], vv[2], vv[3], r[61])
        o[i, 25] = float(cd)
        oi[i, 9] = int(cax)
        cl, clax = tf.leg2_limit_clip(vv[4], vv[5], k[0], k[1], k[2], vv[2], vv[3], r[39], r[61])
        o[i, 26] = float(cl)
        oi[i, 10] = int(clax)
        sp, spax = tf.seg_piece(vv[0], vv[1], k[0], k[1], k[2], vv[2], vv[3], r[62])
        o[i, 27] = float(sp)
        oi[i, 11] = int(spax)
        gl, gok = em.straggle_attempt_gamma(r[14], r[15], r[16], r[17], r[18], r[19])
        o[i, 28] = float(gl)
        oi[i, 12] = int(gok)
        o[i, 29] = float(sc.loglog_slope(r[63], r[64], r[65]))
        o[i, 30] = float(sc.let_ramp_slope(r[66], r[67], r[68], r[69], r[70]))
        ps, pe = sc.piece_state(r[66], r[71], r[69], r[72], r[73])
        o[i, 31], o[i, 32] = float(ps), float(pe)
        pm1, pm2 = sc.piece_moments(r[76], r[71], r[74])
        o[i, 33], o[i, 34] = float(pm1), float(pm2)
        o[i, 35] = float(sc.channel_value(k[9], r[75], r[74], pm1, pm2, r[76], r[77]))
        lbi, lbf, lbin = sc.lookup_bin(r[78], r[79], r[80], k[10], k[11])
        oi[i, 13], o[i, 36], oi[i, 14] = int(lbi), float(lbf), int(lbin)
        oi[i, 15] = int(sc.spectrum_bin(r[78], r[79], r[80], k[12], k[11]))
        o[i, 37] = float(tf.range_in_bin(r[81], r[82], r[83], r[84], r[85]))
        sn, snr = tf.select_step_nuclear(r[46], r[47], r[48], r[49], r[86])
        o[i, 38] = float(sn)
        oi[i, 16] = int(snr)
    return o, oi, ov


REAL_OUTPUTS = {
    0: "pv_mev", 1: "beta2", 2: "gamma", 3: "tmax_mev", 4: "csda_mean_loss",
    5: "bohr_variance", 6: "straggle_loss", 7: "fdm", 8: "scattering_power_dm",
    9: "polar_deflection", 10: "log_bin_frac", 11: "lerp", 12: "interp_exp",
    13: "plane_position", 14: "dda_distance", 15: "leg2_length", 16: "range_step_limit",
    17: "eloss_step_limit", 18: "select_step", 19: "ray_box_in", 20: "ray_box_out",
    21: "gauss_pair_0", 22: "gauss_pair_1", 23: "gauss_one", 24: "scattering_variance_birth",
    25: "dda_clip_distance", 26: "leg2_clip_length", 27: "seg_piece_length",
    28: "straggle_gamma_loss", 29: "loglog_slope", 30: "let_ramp_slope", 31: "piece_s_bar",
    32: "piece_e_bar", 33: "piece_m1", 34: "piece_m2", 35: "channel_value", 36: "lookup_frac",
    37: "range_in_bin", 38: "select_step_nuclear",
}  # fmt: skip
INT_OUTPUTS = {
    0: "straggle_ok", 1: "log_bin_index", 2: "dda_axis", 3: "leg2_axis", 4: "select_reason",
    5: "grid_ix", 6: "grid_iy", 7: "grid_iz", 8: "grid_inside",
    9: "dda_clip_axis", 10: "leg2_clip_axis", 11: "seg_piece_axis", 12: "straggle_gamma_ok",
    13: "lookup_index", 14: "lookup_in_domain", 15: "spectrum_bin",
    16: "select_reason_nuclear",
}  # fmt: skip
VEC_OUTPUTS = {0: "rotate_dir", 1: "point_on_hinge", 2: "basis_e1", 3: "basis_e2"}


def _fdm_abs_budget(
    pv_om: np.ndarray, p1v1: np.ndarray, pv_log: np.ndarray, shift: float
) -> np.ndarray:
    """Absolute float32 error budget of f_dM = c0 + c1 l1 + c2 l2 + c3 l1 l2 (clamped at 0), with
    ``l1 = lg(om) - shift``, ``om = 1 - (pv / p1v1)^2`` and ``l2 = lg(pv)``.

    Derivation: the ratio (relative error u), its square (u) and the difference ``1 - r^2`` give
    an absolute error of ``om`` of at most ``3 u`` (r^2 <= 1), hence an absolute error of
    ``l1`` of ``3 u / (om ln 10)``: the cancellation near pv -> p1v1 is the ill-conditioned part
    and gives ``(|c1| + |c3 l2|) 3 u / (om ln 10)`` in f. The remaining roundings (two logarithms,
    two products, three sums, about 6) each contribute u times the magnitude of the terms."""
    c = FDM_COEFFICIENTS
    om = 1.0 - (pv_om / p1v1) ** 2
    pos = om > 0.0
    omp = np.where(pos, om, 1.0)
    l1, l2 = np.log10(omp) - shift, np.log10(pv_log)
    terms = abs(c[0]) + np.abs(c[1] * l1) + np.abs(c[2] * l2) + np.abs(c[3] * l1 * l2)
    cond = (abs(c[1]) + abs(c[3] * l2)) * 3.0 * U / (omp * LN10)
    return np.where(pos, cond + 6.0 * U * terms, 0.0)


def _box_muller_budget(u0: np.ndarray) -> np.ndarray:
    """Absolute float32 error of ``r cos(2 pi u1)`` / ``r sin(2 pi u1)``, ``r = sqrt(-2 ln u0)``.

    Derivation: the angle ``2 pi u1`` has an absolute error of at most ``2 pi (u + u)`` (rounding
    of the constant, rounding of the product; u1 < 1), the sine/cosine adds u; ``r`` has relative
    error 1.5 u (logarithm of an exact float32 argument, which stays relatively accurate near
    u0 -> 1, then the square root). The error is absolute (``r (4 pi u + u + 1.5 u)``): near a zero
    of the cosine the relative error diverges, whatever u0 is. ``r`` is largest for the smallest
    u0: ``r <= sqrt(2 |ln u0|_max)`` (5.2 in the pool)."""
    r = np.sqrt(-2.0 * np.log(u0))
    return r * (4.0 * math.pi + 2.5) * U


def _float32_budgets(x: np.ndarray) -> dict[int, tuple[float, np.ndarray | float]]:
    """Per-function float32 budgets ``{output column: (rtol, atol)}`` (atol may be per element);
    columns not listed use ``TOL['float32']``. All derived analytically (see each function)."""
    out: dict[int, tuple[float, np.ndarray | float]] = {}
    # interp_exp = exp(ly0 (1 - f) + ly1 f): the exponent is a sum of terms of at most 10 in
    # magnitude with 4 roundings (1 - f, two products, the sum), so its absolute error is
    # (3 * 10 + 1) u and the relative error of the exponential is that plus u for exp itself.
    out[12] = (32.0 * U, 1e-6)
    # range_in_bin = r_i + h f_i g, g = (exp(x) - 1) / d with x = d phi (exact float32 arguments).
    # exp branch: x has relative error u (one product), exp(x) a further u, the subtraction of 1 a
    # rounding u |exp(x) - 1|: absolute error of exp(x) - 1 is at most
    # u (e^x (|x| + 1) + |e^x - 1|);
    # divided by |d| plus u |g| for the division. Series branch (|x| < 1e-5): about 5 roundings of
    # terms of at most 1, 5 u |g|. The two further products and the final sum add 3 u of the
    # result, and the budget is doubled (conservative: the kernel may fuse operations, the
    # compiled and Python-scope evaluations of the exp may differ by an ulp).
    r_i, f_i, d_i, h_i, ph = (x[:, c].astype(np.float64) for c in range(81, 86))
    xx_ = d_i * ph
    gg = np.where(np.abs(xx_) >= 1e-5, np.expm1(xx_) / np.where(d_i == 0.0, 1.0, d_i), ph)
    exp_err = U * (np.exp(xx_) * (np.abs(xx_) + 1.0) + np.abs(np.expm1(xx_)))
    g_err = np.where(
        np.abs(xx_) >= 0.9e-5,
        exp_err / np.maximum(np.abs(d_i), 1e-300) + U * np.abs(gg),
        5.0 * U * np.abs(gg),
    )
    out[37] = (0.0, 2.0 * (h_i * f_i * g_err + 3.0 * U * np.abs(r_i + h_i * f_i * gg)) + 1e-6)
    # f_dM: cancellation of 1 - (pv/p1v1)^2 (see _fdm_abs_budget)
    fdm_abs = _fdm_abs_budget(x[:, 20], x[:, 21], x[:, 20], 0.0)
    out[7] = (1e-6, fdm_abs + 1e-6)
    # scattering power: f_dM (z E_s / pv)^2 rho inv_xs / 10; the error of f_dM is multiplied by the
    # weight W, the 8 further roundings (q, q^2, three products, a division) meet the default rtol
    q = 1.0 * 15.0 / x[:, 20]
    w = q * q * x[:, 23] * x[:, 22] / 10.0
    out[8] = (1e-6, w * fdm_abs + 1e-6)
    # birth variance: the same with om from pv_end, l2 from pv_mid, l1 shifted by 1/ln 10, times s
    shift = 0.4342944819032518
    fdm_b = _fdm_abs_budget(x[:, 55], x[:, 56], x[:, 54], shift)
    wb = (15.0 / x[:, 54]) ** 2 * x[:, 59] * x[:, 58] / 10.0 * x[:, 60]
    out[24] = (1e-6, wb * fdm_b + 1e-6)
    # Box-Muller normals (pair and single use the same u0, u1 columns 52/53)
    bm = _box_muller_budget(x[:, 52])
    out[21] = (1e-6, bm + 1e-6)
    out[22] = (1e-6, bm + 1e-6)
    out[23] = (1e-6, bm + 1e-6)
    out[6] = _straggle_budget(x)
    # straggle_attempt_gamma (the production default, bohr_gamma_v1): the Gamma branch of
    # straggle_attempt at EVERY ratio (no Gaussian branch); the same first-order relative budget.
    out[28] = _straggle_budget(x, gamma_everywhere=True)
    # loglog_slope = (ly1 - ly0) inv_dl: the operands are exact float32 arguments, so the difference
    # has one rounding (u |diff|); the budget 4 u max(|ly0|, |ly1|) inv_dl covers it and the
    # product's rounding with a margin for operands that were themselves computed (table reads)
    out[29] = (1e-6, 4.0 * U * np.maximum(np.abs(x[:, 63]), np.abs(x[:, 64])) * x[:, 65] + 1e-6)
    # dda_next_clip, leg2_limit_clip, seg_piece: dda_next adds a subtraction, a division and a
    # product-sum per axis (about 4 roundings, no cancellation of computed quantities: the
    # operands are exact float32 arguments); the clip distance (z_clip - p_z) / u_z is one rounded
    # difference and one division; comparisons and min/max are exact. All stay within the default.
    return out


def _straggle_budget(x: np.ndarray, *, gamma_everywhere: bool = False) -> tuple[float, np.ndarray]:
    """Absolute float32 error budget of the sampled energy loss (column 6).

    Gaussian branch (ratio >= 3): ``loss = clamp(mean + sigma x)``, x a Box-Muller normal of
    absolute error ``_box_muller_budget``: error ``sigma (eps_x + 3 u |x|) + 4 u |loss|``.

    Gamma branch: first-order relative error of ``loss = (mean / k) d v [u3^(1/k)]`` with
    ``ratio = mean / sqrt(var)`` (2 u), ``k = ratio^2`` (5 u), ``d = a - 1/3`` (at most 10 u, as
    ``k / d <= 1.5``), ``c = 1 / sqrt(9 d)`` (10 u), ``v1 = 1 + c x`` with
    ``dv1 = c eps_x + 10 u c |x| + u v1`` and ``v = v1^3`` (``3 dv1 / v1 + 2 u``), the product and
    quotient (about 9 u) and, for k < 1, the power ``u3^(1/k)`` whose exponent ``1/k`` (relative
    error 6 u) multiplies ``y = |ln u3| / k``: ``rel = 3 dv1 / v1 + 21 u + 6 u y``. The relative
    error diverges as v1 -> 0 (such draws are rare and tiny: v = v1^3), so it is an error budget
    proportional to the value, never a constant."""
    mean, var = x[:, 14], x[:, 15]
    u0, u1, u3 = x[:, 16], x[:, 17], x[:, 19]
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma = np.sqrt(var)
        ratio = np.where(sigma > 0.0, mean / np.where(sigma > 0.0, sigma, 1.0), np.inf)
        r = np.sqrt(-2.0 * np.log(u0))
        xn = r * np.cos(2.0 * math.pi * u1)
        eps_x = r * (4.0 * math.pi + 2.5) * U
        gauss = (ratio >= 3.0) & (not gamma_everywhere)
        k = ratio * ratio
        a = np.where(k < 1.0, k + 1.0, k)
        d = a - 1.0 / 3.0
        c = 1.0 / np.sqrt(9.0 * d)
        v1 = 1.0 + c * xn
        dv1 = c * eps_x + 10.0 * U * c * np.abs(xn) + U * np.abs(v1)
        y = np.where(k < 1.0, np.abs(np.log(u3)) / k, 0.0)
        rel_g = 3.0 * dv1 / np.where(v1 > 0.0, v1, 1.0) + 21.0 * U + 6.0 * U * y
        loss = np.where(gauss, 0.0, mean / k * d * np.where(v1 > 0, v1, 0.0) ** 3)
        abs_gauss = sigma * (eps_x + 3.0 * U * np.abs(xn)) + 4.0 * U * (mean + sigma * np.abs(xn))
        atol = np.where(gauss, abs_gauss, rel_g * np.abs(loss))
    atol = np.where(np.isfinite(atol), atol, 0.0)
    return 1e-6, atol + 1e-6


EXERCISED = {
    "kin": {"pv_mev", "beta2", "gamma", "tmax_mev"},
    "em": {
        "csda_mean_loss", "bohr_variance", "straggle_attempt", "straggle_attempt_gamma", "fdm",
        "scattering_power_dm", "scattering_variance_birth", "polar_deflection", "rotate_dir",
    },
    "sc": {
        "loglog_slope", "let_ramp_slope", "piece_state", "piece_moments", "channel_value",
        "lookup_bin", "spectrum_bin",
    },
    "tf": {
        "lerp", "interp_exp", "range_in_bin", "log_bin_index", "plane_position", "dda_next",
        "dda_next_clip",
        "leg2_limit", "leg2_limit_clip", "seg_piece", "range_step_limit", "eloss_step_limit",
        "select_step", "select_step_nuclear", "point_on_hinge", "grid_index", "ray_box",
        "gauss_pair", "gauss_one",
        "orthonormal_basis",
    },
    "nu": {
        "nuclear_step_limit", "thinning_accept", "select_target", "poisson_inverse",
        "multiplicity_round", "grid_locate", "inv_cdf_bin", "inv_cdf_sample", "kalbach_a",
        "kalbach_cdf", "kalbach_pdf", "kalbach_mu",
        "residual_invariant_mass", "residual_mass_ok", "cm_boost", "boost_z", "cm_to_lab",
        "choose_target", "sample_event",
    },
}  # fmt: skip
"""Functions exercised by the harness (kernel and twin sides, both in ``_make_kernel`` and
``_python_scope``), per namespace. The nuclear namespace ``nu`` (V3-005A) is
exercised, kernel against twin on recorded inputs, by ``tests/ionmc/test_nuclear_funcs.py``
(P5), which asserts it covers
exactly this set. ``test_u1_harness_covers_every_shared_function`` requires this to
equal the callables of the twin and of the Warp namespace, so that a new shared function cannot be
left out of U1 silently."""


def test_u1_harness_covers_every_shared_function() -> None:
    from ionmc.physics.em import make_em as _em

    def names(ns: Any) -> set[str]:
        return {k for k, val in vars(ns).items() if callable(val) and k != "vec3"}

    for key, factory in (
        ("kin", make_kinematics),
        ("em", _em),
        ("tf", make_transport_funcs),
        ("sc", make_scoring_funcs),
        ("nu", make_nuclear),
    ):
        assert names(python_twin(factory)) == EXERCISED[key], key
        assert names(factory(wp.float64)) == EXERCISED[key], key


@pytest.mark.parametrize("precision", ["float64", "float32"])
def test_u1_python_twin_equals_warp_cpu_kernel(precision: str) -> None:
    real = wp.float64 if precision == "float64" else wp.float32
    rtol, atol = TOL[precision]
    x, v, ii = _args(precision)
    budgets = _float32_budgets(x) if precision == "float32" else {}
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
    if precision == "float32":
        # the bin index and fraction are discontinuous at bin edges (a float32 rounding moves a
        # value to the neighbouring bin: index i + 1 / fraction 0 against index i / fraction 1);
        # the continuous position i + f is compared instead (t = ln x * inv_dl, relative error
        # <= 3 u, so the default rtol applies)
        t_py, t_k = p_oi[:, 1] + p_o[:, 10], k_oi[:, 1] + k_o[:, 10]
        assert np.allclose(t_py, t_k, rtol=rtol, atol=atol), "log_bin position i + f differs"
    if precision == "float32":
        # bin functions on a uniform axis: t = (v - a0) inv_da with v = x (linear) or ln x. The
        # float32 error of t is at most 4 u ((|v| + |a0|) inv_da + |t|) (one rounded difference,
        # one product, the logarithm); the comparisons below skip the rows whose t lies within
        # that distance of a bin edge (an integer), where the discontinuous outputs may differ
        xx = x.astype(np.float32).astype(np.float64)
        a0, inv = xx[:, 79], xx[:, 80]
        vv = np.where(ii[:, 11] == 1, np.log(np.maximum(xx[:, 78], 1e-30)), xx[:, 78])
        tt = (vv - a0) * inv
        tol_t = 4.0 * U * ((np.abs(vv) + np.abs(a0)) * inv + np.abs(tt))
        near = np.abs(tt - np.round(tt)) <= tol_t
        cont_py = p_oi[:, 13] + p_o[:, 36]
        cont_k = k_oi[:, 13] + k_o[:, 36]
        assert np.all(np.abs(cont_py - cont_k) <= tol_t + rtol * np.abs(cont_k) + atol), (
            "lookup_bin position i + f differs"
        )
        assert np.array_equal(p_oi[~near, 14], k_oi[~near, 14]), "lookup_in_domain differs"
        assert np.array_equal(p_oi[~near, 15], k_oi[~near, 15]), "spectrum_bin differs"
        assert (~near).sum() > 0.9 * n
    for col, name in REAL_OUTPUTS.items():
        if precision == "float32" and col in (10, 36):
            continue
        a, b = p_o[:, col], k_o[:, col]
        assert np.all(np.isfinite(a)) and np.all(np.isfinite(b)), f"{name}: NaN or inf"
        col_rtol, col_atol = budgets.get(col, (rtol, atol))
        bad = ~(np.abs(a - b) <= col_atol + col_rtol * np.abs(b))
        assert not bad.any(), (
            f"{name}: {bad.sum()} of {n} mismatches; first python={a[bad][:3]} kernel={b[bad][:3]}"
        )
    for col, name in INT_OUTPUTS.items():
        if precision == "float32" and col in (1, 13, 14, 15):
            continue
        assert np.array_equal(p_oi[:, col], k_oi[:, col]), f"{name}: integer outputs differ"
    for col, name in VEC_OUTPUTS.items():
        a, b = p_ov[:, col], k_ov[:, col]
        assert np.all(np.isfinite(a)) and np.all(np.isfinite(b)), f"{name}: NaN or inf"
        assert np.allclose(a, b, rtol=rtol, atol=atol), f"{name}: components differ"
