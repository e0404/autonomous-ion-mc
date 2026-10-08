# mypy: ignore-errors
# (Warp function sources use runtime precision types in annotations; see _wpfunc.py.)
"""Shared transport helper functions: table-bin location, DDA, step limits, hinge, scoring.

``make_transport_funcs(R)`` returns precision-generic Warp functions for kernels; the reference
backend uses their pure-Python twins (``ionmc._wpfunc.python_twin``: the same source text without
Warp, so no Warp call happens at Python scope). They are pure: table reads are split
into the shared bin location ``log_bin_index``, a backend memory read and the shared
interpolation ``interp_exp`` / ``lerp``. Units: mm for lengths, MeV for energies; a "plane
index" ``i`` bounds voxel ``i`` between ``origin + i*spacing`` and ``origin + (i+1)*spacing``.
Reasons returned by ``select_step``: 0 geometry, 1 energy-loss limit, 2 range limit, 3 maximum
step.
"""

import functools
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func

wp.set_module_options({"enable_backward": False})

REASON_GEOMETRY = 0
REASON_ELOSS = 1
REASON_RANGE = 2
REASON_MAX_STEP = 3
REASON_NUCLEAR = 4
BIG_LENGTH_MM = 1.0e30
DIRECTION_EPS = 1.0e-12
LEG_TOLERANCE = {"float64": 1.0e-9, "float32": 1.0e-5}
"""Relative tolerance with which a second hinge leg that ends at a voxel plane within
rounding is treated as reaching it (so planned boundary steps always cross)."""


@functools.cache
def make_transport_funcs(real: type) -> SimpleNamespace:
    """Return the shared transport functions for precision ``real``."""
    name = check_real(real)
    v3 = wp.types.vector(length=3, dtype=real)
    big = wp.constant(real(BIG_LENGTH_MM))
    eps = wp.constant(real(DIRECTION_EPS))
    leg_tol = wp.constant(real(LEG_TOLERANCE[name]))
    two_pi = wp.constant(real(6.283185307179586))

    @named_func(name)
    def lerp(y0: real, y1: real, f: real) -> real:
        """Linear interpolation ``y0 (1 - f) + y1 f``."""
        return y0 * (real(1.0) - f) + y1 * f

    @named_func(name)
    def interp_exp(ly0: real, ly1: real, f: real) -> real:
        """``exp`` of the linear interpolation of two logarithms (log-log interpolation)."""
        return wp.exp(ly0 * (real(1.0) - f) + ly1 * f)

    @named_func(name)
    def range_in_bin(r_i: real, f_i: real, d_i: real, h: real, phi: real) -> real:
        """CSDA range [g/cm2] at the fraction ``phi`` of a log-log interval: ``R_i + h f_i g`` with
        ``g = (exp(d_i phi) - 1) / d_i``, the exact integral of the log-log interpolated
        ``a E / S`` (``f_i`` its value at the node, ``d_i`` the change of its logarithm over the
        interval, ``h`` the interval width in ``ln E``; decision 0039, V3-003D). For
        ``|d_i phi| < 1e-5`` the series ``phi (1 + x/2 + x^2/6)``, ``x = d_i phi``, is used
        (no ``expm1`` is guaranteed in Warp)."""
        x = d_i * phi
        g = phi * (real(1.0) + x * real(0.5) + x * x / real(6.0))
        if wp.abs(x) >= real(1.0e-5):
            g = (wp.exp(x) - real(1.0)) / d_i
        return r_i + h * f_i * g

    @named_func(name)
    def log_bin_index(x: real, lx0: real, inv_dl: real, n: int) -> tuple[int, real]:
        """Bin ``i`` in [0, n-2] and fraction ``f`` in [0, 1] of ``ln x`` on a uniform grid.

        ``lx0`` is ``ln`` of the first grid point, ``inv_dl`` the inverse step in ``ln x``,
        ``n`` the number of grid points. Arguments outside the grid are clamped (no
        extrapolation).
        """
        lx = wp.log(wp.max(x, real(1.0e-30)))
        t = (lx - lx0) * inv_dl
        i = int(wp.floor(t))
        i = wp.max(wp.min(i, n - 2), 0)
        f = wp.min(wp.max(t - real(i), real(0.0)), real(1.0))
        return i, f

    @named_func(name)
    def plane_position(index: int, upward: int, origin_a: real, spacing_a: real) -> real:
        """Coordinate of the voxel plane bounding voxel ``index`` (upper if ``upward`` = 1)."""
        return real(index + upward) * spacing_a + origin_a

    @named_func(name)
    def dda_next(
        p: v3, u: v3, ix: int, iy: int, iz: int, origin: v3, spacing: v3
    ) -> tuple[real, int]:
        """Distance [mm] to the next voxel plane and its axis (0, 1, 2); ties go to the lowest
        axis. Directions with ``|u_a| < 1e-12`` never reach a plane of axis ``a``. The
        distance to the plane of voxel ``i`` in the travel sense is
        ``max(((i + (u > 0)) spacing + origin - p) / u, 0)``. (Written without nested function
        calls: calls are expensive in Python scope.)"""
        dx = big
        if wp.abs(u[0]) >= eps:
            upx = 0
            if u[0] > real(0.0):
                upx = 1
            dx = wp.max((real(ix + upx) * spacing[0] + origin[0] - p[0]) / u[0], real(0.0))
        dy = big
        if wp.abs(u[1]) >= eps:
            upy = 0
            if u[1] > real(0.0):
                upy = 1
            dy = wp.max((real(iy + upy) * spacing[1] + origin[1] - p[1]) / u[1], real(0.0))
        dz = big
        if wp.abs(u[2]) >= eps:
            upz = 0
            if u[2] > real(0.0):
                upz = 1
            dz = wp.max((real(iz + upz) * spacing[2] + origin[2] - p[2]) / u[2], real(0.0))
        best = dx
        axis = int(0)
        if dy < best:
            best = dy
            axis = 1
        if dz < best:
            best = dz
            axis = 2
        return best, axis

    @named_func(name)
    def leg2_limit(
        ph: v3,
        d1: v3,
        ix: int,
        iy: int,
        iz: int,
        origin: v3,
        spacing: v3,
        target_mm: real,
    ) -> tuple[real, int]:
        """Length of the second hinge leg [mm] and the crossed axis (-1 if it is not cut).

        The leg is cut at the next voxel plane if that plane lies within the target length
        (relative tolerance for rounding of planned boundary steps)."""
        dd, ax = dda_next(ph, d1, ix, iy, iz, origin, spacing)
        leg = target_mm
        axis = int(-1)
        if dd <= target_mm * (real(1.0) + leg_tol):
            leg = dd
            axis = ax
        return leg, axis

    @named_func(name)
    def dda_next_clip(
        p: v3, u: v3, ix: int, iy: int, iz: int, origin: v3, spacing: v3, z_clip: real
    ) -> tuple[real, int]:
        """:func:`dda_next` with an extra world face at ``z = z_clip`` (axis 3): the transport
        distance to the next voxel plane or to the clip plane, whichever is nearer (the clip plane
        wins ties). With ``z_clip`` at ``BIG_LENGTH_MM`` it is :func:`dda_next`."""
        d, ax = dda_next(p, u, ix, iy, iz, origin, spacing)
        if u[2] > eps:
            dclip = wp.max((z_clip - p[2]) / u[2], real(0.0))
            if dclip <= d:
                d = dclip
                ax = 3
        return d, ax

    @named_func(name)
    def leg2_limit_clip(
        ph: v3,
        d1: v3,
        ix: int,
        iy: int,
        iz: int,
        origin: v3,
        spacing: v3,
        target_mm: real,
        z_clip: real,
    ) -> tuple[real, int]:
        """:func:`leg2_limit` with the clip plane of :func:`dda_next_clip` (axis 3)."""
        dd, ax = dda_next_clip(ph, d1, ix, iy, iz, origin, spacing, z_clip)
        leg = target_mm
        axis = int(-1)
        if dd <= target_mm * (real(1.0) + leg_tol):
            leg = dd
            axis = ax
        return leg, axis

    @named_func(name)
    def seg_piece(
        p: v3, u: v3, ix: int, iy: int, iz: int, origin: v3, spacing: v3, remaining_mm: real
    ) -> tuple[real, int]:
        """Length [mm] of the next piece of a straight segment inside voxel ``(ix, iy, iz)`` of a
        grid and the axis of the plane that ends it (-1 if the segment ends inside the voxel):
        the shared per-grid DDA step of the track-length scoring. ``remaining_mm`` is the length
        of the segment still to be walked; a plane at exactly that distance does not end the
        piece (the segment ends there)."""
        d, ax = dda_next(p, u, ix, iy, iz, origin, spacing)
        piece = remaining_mm
        axis = int(-1)
        if d < remaining_mm:
            piece = d
            axis = ax
        return piece, axis

    @named_func(name)
    def range_step_limit(r_mm: real, alpha: real, rho_f_mm: real) -> real:
        """Geant4 range step function: ``alpha R + rho_f (1 - alpha)(2 - rho_f/R)`` for
        ``R > rho_f``, else ``R`` (final range step); ``r_mm`` is the residual CSDA range."""
        s = r_mm
        if r_mm > rho_f_mm:
            s = alpha * r_mm + rho_f_mm * (real(1.0) - alpha) * (real(2.0) - rho_f_mm / r_mm)
        return s

    @named_func(name)
    def eloss_step_limit(e_mev: real, s_lin_mev_mm: real, frac: real) -> real:
        """Step length [mm] at which the linear energy-loss estimate reaches ``frac`` of E."""
        return frac * e_mev / s_lin_mev_mm

    @named_func(name)
    def select_step(
        d_geo_mm: real, s_eloss_mm: real, s_range_mm: real, s_max_mm: real
    ) -> tuple[real, int]:
        """Step length [mm] and reason: the smallest limit; geometry wins ties."""
        s = s_eloss_mm
        reason = int(1)
        if s_range_mm < s:
            s = s_range_mm
            reason = 2
        if s_max_mm < s:
            s = s_max_mm
            reason = 3
        if d_geo_mm <= s:
            s = d_geo_mm
            reason = 0
        return s, reason

    @named_func(name)
    def select_step_nuclear(
        d_geo_mm: real,
        s_eloss_mm: real,
        s_range_mm: real,
        s_max_mm: real,
        d_nuc_mm: real,
    ) -> tuple[real, int]:
        """:func:`select_step` plus the nuclear limit ``d_nuc_mm`` (reason 4; decision 0041): the
        smallest limit, geometry wins ties exactly as in :func:`select_step`; the nuclear limit
        replaces the EM limits only when STRICTLY smaller (an EM limit wins a tie with it), so
        ``d_nuc_mm = BIG_LENGTH_MM`` or ``+inf`` returns exactly the result of
        :func:`select_step`."""
        s = s_eloss_mm
        reason = int(1)
        if s_range_mm < s:
            s = s_range_mm
            reason = 2
        if s_max_mm < s:
            s = s_max_mm
            reason = 3
        if d_nuc_mm < s:
            s = d_nuc_mm
            reason = 4
        if d_geo_mm <= s:
            s = d_geo_mm
            reason = 0
        return s, reason

    @named_func(name)
    def point_on_hinge(p0: v3, d0: v3, leg1_mm: real, d1: v3, l_mm: real) -> v3:
        """Position at path length ``l_mm`` along leg 1 (direction ``d0``, length ``leg1_mm``)
        followed by leg 2 (direction ``d1``)."""
        out = p0 + d0 * l_mm
        if l_mm > leg1_mm:
            out = p0 + d0 * leg1_mm + d1 * (l_mm - leg1_mm)
        return out

    @named_func(name)
    def grid_index(
        p: v3, origin: v3, inv_spacing: v3, nx: int, ny: int, nz: int
    ) -> tuple[int, int, int, int]:
        """Voxel indices of ``p`` (floor) and ``inside`` (1/0); the far edge is outside."""
        ix = int(wp.floor((p[0] - origin[0]) * inv_spacing[0]))
        iy = int(wp.floor((p[1] - origin[1]) * inv_spacing[1]))
        iz = int(wp.floor((p[2] - origin[2]) * inv_spacing[2]))
        inside = 0
        if ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz:
            inside = 1
        return ix, iy, iz, inside

    @named_func(name)
    def ray_box(p: v3, u: v3, lo: v3, hi: v3) -> tuple[real, real]:
        """Slab intersection of the line ``p + t u`` with the box [lo, hi]: ``(t_in, t_out)``.
        A miss returns ``t_in > t_out``."""
        t_in = -big
        t_out = big
        for a in range(3):
            if wp.abs(u[a]) < eps:
                if p[a] < lo[a] or p[a] > hi[a]:
                    t_in = big
                    t_out = -big
            else:
                t1 = (lo[a] - p[a]) / u[a]
                t2 = (hi[a] - p[a]) / u[a]
                t_in = wp.max(t_in, wp.min(t1, t2))
                t_out = wp.min(t_out, wp.max(t1, t2))
        return t_in, t_out

    @named_func(name)
    def gauss_pair(u0: real, u1: real) -> tuple[real, real]:
        """Two independent standard normals by the Box-Muller transform."""
        r = wp.sqrt(real(-2.0) * wp.log(u0))
        return r * wp.cos(two_pi * u1), r * wp.sin(two_pi * u1)

    @named_func(name)
    def gauss_one(u0: real, u1: real) -> real:
        """One standard normal (cosine branch of Box-Muller)."""
        return wp.sqrt(real(-2.0) * wp.log(u0)) * wp.cos(two_pi * u1)

    @named_func(name)
    def orthonormal_basis(n: v3) -> tuple[v3, v3]:
        """Deterministic orthonormal pair perpendicular to the unit vector ``n``
        (Duff, Burgess, Christensen, Hery, Kensler, Max, Tatarchuk, JCGT 6 (2017) 1)."""
        sgn = real(1.0)
        if n[2] < real(0.0):
            sgn = real(-1.0)
        a = real(-1.0) / (sgn + n[2])
        b = n[0] * n[1] * a
        e1 = v3(real(1.0) + sgn * n[0] * n[0] * a, sgn * b, -sgn * n[0])
        e2 = v3(b, sgn + n[1] * n[1] * a, -n[1])
        return e1, e2

    return SimpleNamespace(
        lerp=lerp,
        interp_exp=interp_exp,
        range_in_bin=range_in_bin,
        log_bin_index=log_bin_index,
        plane_position=plane_position,
        dda_next=dda_next,
        dda_next_clip=dda_next_clip,
        leg2_limit_clip=leg2_limit_clip,
        leg2_limit=leg2_limit,
        seg_piece=seg_piece,
        range_step_limit=range_step_limit,
        eloss_step_limit=eloss_step_limit,
        select_step=select_step,
        select_step_nuclear=select_step_nuclear,
        point_on_hinge=point_on_hinge,
        grid_index=grid_index,
        ray_box=ray_box,
        gauss_pair=gauss_pair,
        gauss_one=gauss_one,
        orthonormal_basis=orthonormal_basis,
        vec3=v3,
        real=name,
    )
