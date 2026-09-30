"""Shared condensed-history step physics (single source for all backends).

Every function here is written in the Warp kernel subset and decorated with
``@wp.func``. The module is loaded twice by :mod:`ionmc.transport.shared`:
once with ``MATH`` bound to Warp builtins (compiled into CPU/CUDA kernels in
float32 or float64) and once with ``MATH`` bound to Python's ``math`` module
(executed by the CPython interpreter in float64 by the reference backend).
The physics text is therefore identical across the three execution paths
(decision 0038).

Numeric literals are written as ``type(x)(value)`` so that generic kernels
instantiate them in the precision of their arguments (Warp types bare
literals as float32); in Python scope this is simply ``float(value)``.

Conventions: lengths mm, energies MeV (kinetic energy per nucleon ``t`` in
MeV/u), densities g/cm³, mass stopping powers MeV cm²/g, ranges g/cm².
Tables are 3-D arrays indexed ``[species, material, grid]`` on logarithmic
grids described by ``(log_min, inv_dlog, n)``.
"""

from typing import Any

import warp as wp

MATH: Any  # injected by ionmc.transport.shared before the module body executes

# Overshoot applied after reaching a voxel face so the next lookup lands in the
# neighbouring voxel in every precision: float32 positions up to ~500 mm have a
# resolution of ~6e-5 mm, so 2.5e-4 mm is at least four ulps. The overshoot is
# part of the step length (energy loss is computed for it), so it only
# reassigns a 0.25 µm slice per crossing to the neighbouring voxel.
BOUNDARY_OVERSHOOT_MM = 2.5e-4


@wp.func
def boundary_overshoot(x: Any) -> Any:
    return type(x)(2.5e-4)


@wp.func
def log_lookup(
    table: Any, s: int, m: int, x: Any, log_min: Any, inv_dlog: Any, n: int
) -> Any:
    """Linear interpolation on a logarithmic grid; clamps to the table ends."""
    u = (MATH.log(x) - log_min) * inv_dlog
    u = MATH.max(u, type(x)(0.0))
    u = MATH.min(u, type(x)(n - 1) - type(x)(1.0e-9))
    i = int(MATH.floor(u))
    f = u - type(x)(i)
    return table[s, m, i] * (type(x)(1.0) - f) + table[s, m, i + 1] * f


@wp.func
def range_lookup(
    log_range: Any, s: int, m: int, t: Any, log_min: Any, inv_dlog: Any, n: int
) -> Any:
    """CSDA range (g/cm²) by log-log interpolation of the stored ln R table."""
    u = (MATH.log(t) - log_min) * inv_dlog
    u = MATH.max(u, type(t)(0.0))
    u = MATH.min(u, type(t)(n - 1) - type(t)(1.0e-9))
    i = int(MATH.floor(u))
    f = u - type(t)(i)
    return MATH.exp(
        log_range[s, m, i] * (type(t)(1.0) - f) + log_range[s, m, i + 1] * f
    )


@wp.func
def energy_after_path(
    log_range: Any,
    s: int,
    m: int,
    t: Any,
    rho_path_g_cm2: Any,
    log_min: Any,
    inv_dlog: Any,
    n: int,
) -> Any:
    """Kinetic energy per nucleon after a mass path, by *exact* inversion of the
    piecewise log-log range table used by :func:`range_lookup` (no separate
    inverse table, so arbitrarily small steps have no interpolation bias).
    Returns 0 when the residual range falls below the table."""
    u = (MATH.log(t) - log_min) * inv_dlog
    u = MATH.max(u, type(t)(0.0))
    u = MATH.min(u, type(t)(n - 1) - type(t)(1.0e-9))
    i = int(MATH.floor(u))
    f = u - type(t)(i)
    r = (
        MATH.exp(log_range[s, m, i] * (type(t)(1.0) - f) + log_range[s, m, i + 1] * f)
        - rho_path_g_cm2
    )
    if r <= MATH.exp(log_range[s, m, 0]):
        return type(t)(0.0)
    ln_r = MATH.log(r)
    # walk down to the cell containing ln_r (ln R increases with index)
    for _ in range(n):
        if i <= 0 or log_range[s, m, i] <= ln_r:
            break
        i = i - 1
    span = log_range[s, m, i + 1] - log_range[s, m, i]
    g = (ln_r - log_range[s, m, i]) / span
    return MATH.exp(log_min + (type(t)(i) + g) / inv_dlog)


@wp.func
def beta_squared(t: Any, mass_per_nucleon_mev: Any) -> Any:
    gamma = type(t)(1.0) + t / mass_per_nucleon_mev
    return type(t)(1.0) - type(t)(1.0) / (gamma * gamma)


@wp.func
def pv_mev(t: Any, mass_nucleus_mev: Any, a: Any) -> Any:
    """Momentum times velocity for total kinetic energy T = a t: T (T + 2m)/(T + m)."""
    big_t = t * a
    return (
        big_t * (big_t + type(t)(2.0) * mass_nucleus_mev) / (big_t + mass_nucleus_mev)
    )


@wp.func
def effective_charge(z: Any, beta: Any) -> Any:
    """Pierce–Blann effective charge z (1 − exp(−125 β z^{−2/3}))."""
    return z * (
        type(beta)(1.0)
        - MATH.exp(
            -type(beta)(125.0)
            * beta
            / MATH.exp(MATH.log(z) * type(beta)(2.0) / type(beta)(3.0))
        )
    )


@wp.func
def straggling_variance_mev2(
    z_eff: Any, z_over_a: Any, beta2: Any, rho_path_g_cm2: Any
) -> Any:
    """Bohr variance (0.1569 MeV² cm²/g) with the relativistic factor over a mass path."""
    return (
        type(beta2)(0.1569)
        * z_eff
        * z_eff
        * z_over_a
        * (type(beta2)(1.0) - type(beta2)(0.5) * beta2)
        / (type(beta2)(1.0) - beta2)
        * rho_path_g_cm2
    )


@wp.func
def sample_gamma(state: Any, k: Any):
    """Marsaglia–Tsang sampler for Gamma(shape k, scale 1); k > 0 (boost for k < 1).

    Returns ``(sample, state)``: Warp passes the RNG state by value, so the
    advanced state must be handed back to the caller (a state passed into a
    ``wp.func`` and not returned would replay the same draws at every call).
    """
    boost = type(k)(1.0)
    kk = k
    if kk < type(k)(1.0):
        boost = MATH.exp(
            MATH.log(MATH.max(type(k)(MATH.randf(state)), type(k)(1.0e-12))) / kk
        )
        kk = kk + type(k)(1.0)
    d = kk - type(k)(1.0) / type(k)(3.0)
    c = type(k)(1.0) / MATH.sqrt(type(k)(9.0) * d)
    x = d  # fallback (the mode) if 64 rejections occur; probability negligible
    for _ in range(64):
        z = type(k)(MATH.randn(state))
        v = type(k)(1.0) + c * z
        if v > type(k)(0.0):
            v = v * v * v
            u = type(k)(MATH.randf(state))
            if MATH.log(MATH.max(u, type(k)(1.0e-12))) < type(k)(
                0.5
            ) * z * z + d - d * v + d * MATH.log(v):
                x = d * v
                break
    return x * boost, state


@wp.func
def sample_energy_loss(state: Any, mean_mev: Any, variance_mev2: Any):
    """Gamma-distributed loss with the given mean and variance (positive, mean-preserving).

    Returns ``(loss, state)`` with the advanced RNG state (see sample_gamma).
    """
    if variance_mev2 <= type(mean_mev)(0.0) or mean_mev <= type(mean_mev)(0.0):
        return mean_mev, state
    k = mean_mev * mean_mev / variance_mev2
    theta = variance_mev2 / mean_mev
    x, state = sample_gamma(state, k)
    return x * theta, state


@wp.func
def f_dm(pv: Any, p1v1: Any) -> Any:
    """Gottschalk's differential-Molière correction factor (arXiv:0908.1413 eq. 40)."""
    ratio = type(pv)(1.0) - (pv / p1v1) * (pv / p1v1)
    ratio = MATH.max(ratio, type(pv)(1.0e-3))
    ratio = MATH.min(ratio, type(pv)(0.97))
    ln10 = type(pv)(2.302585092994046)
    lg_r = MATH.log(ratio) / ln10
    lg_pv = MATH.log(pv) / ln10
    return (
        type(pv)(0.5244)
        + type(pv)(0.1975) * lg_r
        + type(pv)(0.2320) * lg_pv
        - type(pv)(0.0098) * lg_pv * lg_r
    )


@wp.func
def mcs_theta0_squared(
    z: Any, pv: Any, p1v1: Any, rho_x_s_g_cm2: Any, rho_path_g_cm2: Any
) -> Any:
    """Mean squared projected angle over a mass path: z² f_dM (E_s/pv)² ρs/(ρX_S), E_s = 15 MeV."""
    es = type(pv)(15.0)
    return (
        z * z * f_dm(pv, p1v1) * (es / pv) * (es / pv) * rho_path_g_cm2 / rho_x_s_g_cm2
    )


@wp.func
def rotate_x(ux: Any, uy: Any, uz: Any, cos_t: Any, sin_t: Any, phi: Any) -> Any:
    """x component of the direction rotated by polar angle (cos_t, sin_t) and azimuth phi."""
    cp = MATH.cos(phi)
    sp = MATH.sin(phi)
    if MATH.abs(uz) < type(uz)(0.999999):
        pn = MATH.sqrt(ux * ux + uy * uy)
        px = -uy / pn
        py = ux / pn
        qx = -uz * py
        return cos_t * ux + sin_t * (cp * px + sp * qx)
    return cos_t * ux + sin_t * cp


@wp.func
def rotate_y(ux: Any, uy: Any, uz: Any, cos_t: Any, sin_t: Any, phi: Any) -> Any:
    cp = MATH.cos(phi)
    sp = MATH.sin(phi)
    if MATH.abs(uz) < type(uz)(0.999999):
        pn = MATH.sqrt(ux * ux + uy * uy)
        px = -uy / pn
        py = ux / pn
        qy = uz * px
        return cos_t * uy + sin_t * (cp * py + sp * qy)
    return cos_t * uy + sin_t * sp


@wp.func
def rotate_z(ux: Any, uy: Any, uz: Any, cos_t: Any, sin_t: Any, phi: Any) -> Any:
    sp = MATH.sin(phi)
    if MATH.abs(uz) < type(uz)(0.999999):
        pn = MATH.sqrt(ux * ux + uy * uy)
        px = -uy / pn
        py = ux / pn
        qz = ux * py - uy * px
        return cos_t * uz + sin_t * sp * qz
    return cos_t * uz


@wp.func
def distance_to_voxel_boundary(
    x: Any,
    y: Any,
    z: Any,
    ux: Any,
    uy: Any,
    uz: Any,
    ox: Any,
    oy: Any,
    oz: Any,
    dx: Any,
    dy: Any,
    dz: Any,
) -> Any:
    """Distance (mm) along the direction to the next face of the voxel containing (x, y, z)."""
    ix = MATH.floor((x - ox) / dx)
    iy = MATH.floor((y - oy) / dy)
    iz = MATH.floor((z - oz) / dz)
    best = type(x)(1.0e30)
    if ux > type(x)(0.0):
        best = MATH.min(best, (ox + (ix + type(x)(1.0)) * dx - x) / ux)
    elif ux < type(x)(0.0):
        best = MATH.min(best, (ox + ix * dx - x) / ux)
    if uy > type(x)(0.0):
        best = MATH.min(best, (oy + (iy + type(x)(1.0)) * dy - y) / uy)
    elif uy < type(x)(0.0):
        best = MATH.min(best, (oy + iy * dy - y) / uy)
    if uz > type(x)(0.0):
        best = MATH.min(best, (oz + (iz + type(x)(1.0)) * dz - z) / uz)
    elif uz < type(x)(0.0):
        best = MATH.min(best, (oz + iz * dz - z) / uz)
    return MATH.max(best, type(x)(0.0))


@wp.func
def voxel_axis_index(x: Any, o: Any, d: Any, n: int) -> int:
    """Voxel index along one axis, or -1 when outside [o, o + n d)."""
    u = (x - o) / d
    if u < type(x)(0.0):
        return -1
    i = int(MATH.floor(u))
    if i >= n:
        return -1
    return i


@wp.func
def distance_to_box_entry(
    x: Any,
    y: Any,
    z: Any,
    ux: Any,
    uy: Any,
    uz: Any,
    x0: Any,
    y0: Any,
    z0: Any,
    x1: Any,
    y1: Any,
    z1: Any,
) -> Any:
    """Distance along the direction to enter the axis-aligned box [x0,x1]×[y0,y1]×[z0,z1].

    Returns 0 when already inside, and −1 when the ray misses the box. The
    entry point is nudged 1e-6 mm inside so the voxel lookup succeeds.
    """
    inside = x >= x0 and x < x1 and y >= y0 and y < y1 and z >= z0 and z < z1
    if inside:
        return type(x)(0.0)
    tmin = type(x)(0.0)
    tmax = type(x)(1.0e30)
    big = type(x)(1.0e30)
    # x slab
    if MATH.abs(ux) < type(x)(1.0e-12):
        if x < x0 or x >= x1:
            return type(x)(-1.0)
    else:
        ta = (x0 - x) / ux
        tb = (x1 - x) / ux
        tmin = MATH.max(tmin, MATH.min(ta, tb))
        tmax = MATH.min(tmax, MATH.max(ta, tb))
    if MATH.abs(uy) < type(x)(1.0e-12):
        if y < y0 or y >= y1:
            return type(x)(-1.0)
    else:
        ta = (y0 - y) / uy
        tb = (y1 - y) / uy
        tmin = MATH.max(tmin, MATH.min(ta, tb))
        tmax = MATH.min(tmax, MATH.max(ta, tb))
    if MATH.abs(uz) < type(x)(1.0e-12):
        if z < z0 or z >= z1:
            return type(x)(-1.0)
    else:
        ta = (z0 - z) / uz
        tb = (z1 - z) / uz
        tmin = MATH.max(tmin, MATH.min(ta, tb))
        tmax = MATH.min(tmax, MATH.max(ta, tb))
    if tmax < tmin or tmax >= big:
        return type(x)(-1.0)
    return tmin + type(x)(1.0e-6)


@wp.func
def deposit_segment(
    edep: Any,
    ox: Any,
    oy: Any,
    oz: Any,
    dx: Any,
    dy: Any,
    dz: Any,
    nx: int,
    ny: int,
    nz: int,
    x0: Any,
    y0: Any,
    z0: Any,
    x1: Any,
    y1: Any,
    z1: Any,
    energy: Any,
) -> Any:
    """Distribute ``energy`` over the scoring voxels crossed by the segment, proportionally
    to the path length in each (uniform loss rate along the step). Returns the energy
    actually scored (parts outside the grid are not scored)."""
    scored = type(x0)(0.0)
    if energy <= type(x0)(0.0):
        return scored
    length = MATH.sqrt(
        (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0)
    )
    if length <= type(x0)(0.0):
        i = voxel_axis_index(x0, ox, dx, nx)
        j = voxel_axis_index(y0, oy, dy, ny)
        k = voxel_axis_index(z0, oz, dz, nz)
        if i >= 0 and j >= 0 and k >= 0:
            MATH.add3(edep, i, j, k, energy)
            scored = energy
        return scored
    ux = (x1 - x0) / length
    uy = (y1 - y0) / length
    uz = (z1 - z0) / length
    s = type(x0)(0.0)
    eps = boundary_overshoot(x0)
    for _ in range(4096):
        if s >= length:
            break
        x = x0 + ux * s
        y = y0 + uy * s
        z = z0 + uz * s
        seg = (
            distance_to_voxel_boundary(x, y, z, ux, uy, uz, ox, oy, oz, dx, dy, dz)
            + eps
        )
        seg = MATH.min(seg, length - s)
        i = voxel_axis_index(x, ox, dx, nx)
        j = voxel_axis_index(y, oy, dy, ny)
        k = voxel_axis_index(z, oz, dz, nz)
        if i >= 0 and j >= 0 and k >= 0:
            part = energy * seg / length
            MATH.add3(edep, i, j, k, part)
            scored = scored + part
        s = s + seg
    return scored
