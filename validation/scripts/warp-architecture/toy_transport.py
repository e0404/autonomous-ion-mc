"""Cost-representative condensed-history proton kernel (NOT validated physics).
Per step: 3D DDA to next 1 mm voxel boundary, material/density fetch, 3 log-spaced table
lookups (range, inverse range, nuclear sigma), Philox RNG (8 uniforms), Gaussian straggling,
Highland-type MCS with direction rotation, nuclear interaction w/ per-thread secondary stack,
2 atomic adds (dose, LET numerator) into a batch-indexed 2 mm scoring grid."""
import os, sys, time, math, resource
import numpy as np
import warp as wp
wp.config.log_level = getattr(wp, "LOG_WARNING", None) or 30
from philox_lib import make_philox

def build(real, stack_n=8, tag=""):
    philox4x32_10, u01 = make_philox(real)
    vstack = wp.types.vector(length=stack_n, dtype=real)
    R = real
    V3 = wp.types.vector(length=3, dtype=real)
    @wp.func
    def table_lookup(tab: wp.array2d(dtype=R), m: int, x: R, lx0: R, inv_dl: R, nb: int) -> R:
        # log-spaced table, linear interpolation in log(x) index space
        t = (wp.log(x) - lx0) * inv_dl
        t = wp.clamp(t, R(0.0), R(nb - 1) - R(1.0e-3))
        i = int(t)
        f = t - R(i)
        return tab[m, i] * (R(1.0) - f) + tab[m, i + 1] * f

    @wp.func
    def dist_to_boundary(p: R, u: R, i: int, dx: R) -> R:
        if u > R(1.0e-12):
            return (R(i + 1) * dx - p) / u
        if u < R(-1.0e-12):
            return (R(i) * dx - p) / u
        return R(1.0e30)

    @wp.func
    def rotate(u: V3, ct: R, phi: R) -> V3:
        st = wp.sqrt(wp.max(R(0.0), R(1.0) - ct * ct)); cp = wp.cos(phi); sp = wp.sin(phi)
        ux = u[0]; uy = u[1]; uz = u[2]
        den = wp.sqrt(wp.max(R(1.0) - uz * uz, R(0.0)))
        if den < R(1.0e-6):
            return V3(st * cp, st * sp, wp.sign(uz) * ct)
        nx = ux * ct + st * (ux * uz * cp - uy * sp) / den
        ny = uy * ct + st * (uy * uz * cp + ux * sp) / den
        nz = uz * ct - st * cp * den
        n = V3(nx, ny, nz)
        return n / wp.length(n)

    @wp.kernel
    def transport(seed: wp.uint32, hist_offset: int, n_batches: int, e0: R,
                  mat: wp.array3d(dtype=wp.int32), rho: wp.array3d(dtype=R),
                  rng_tab: wp.array2d(dtype=R), irng_tab: wp.array2d(dtype=R), sig_tab: wp.array2d(dtype=R),
                  lE0: R, inv_dlE: R, lR0: R, inv_dlR: R, nb: int,
                  dose: wp.array2d(dtype=wp.float32), letn: wp.array2d(dtype=wp.float32),
                  counters: wp.array(dtype=wp.int64), max_steps: int):
        tid = wp.tid()
        h = tid + hist_offset
        b = h % n_batches
        nx = mat.shape[0]; ny = mat.shape[1]; nz = mat.shape[2]
        dx = R(1.0)
        key = wp.vec2ui(seed, wp.uint32(0))
        # secondary stack (SoA of local fixed-size vectors)
        sx = vstack(); sy = vstack(); sz = vstack(); sE = vstack(); sux = vstack(); suy = vstack(); suz = vstack()
        top = int(0); overflow = int(0)
        pid = int(0)
        p = V3(R(60.0), R(60.0), R(0.0))
        u = V3(R(0.0), R(0.0), R(1.0))
        E = e0
        alive = int(1)
        nstep = int(0); ctr = wp.uint32(0)
        ix = int(wp.floor(p[0])); iy = int(wp.floor(p[1])); iz = int(wp.floor(p[2]))
        while alive == 1:
            if nstep >= max_steps:
                wp.atomic_add(counters, 1, wp.int64(1))   # truncated particle (reported, never silent)
                alive = 0
            if alive == 1 and (ix < 0 or iy < 0 or iz < 0 or ix >= nx or iy >= ny or iz >= nz):
                wp.atomic_add(counters, 2, wp.int64(1))   # escaped
                alive = 0
            if alive == 1:
                m = mat[ix, iy, iz]; dens = rho[ix, iy, iz]
                r4 = philox4x32_10(wp.vec4ui(wp.uint32(h), wp.uint32(pid), ctr, wp.uint32(0)), key)
                r5 = philox4x32_10(wp.vec4ui(wp.uint32(h), wp.uint32(pid), ctr, wp.uint32(1)), key)
                ctr = ctr + wp.uint32(1)
                tx = wp.max(dist_to_boundary(p[0], u[0], ix, dx), R(0.0)); ty = wp.max(dist_to_boundary(p[1], u[1], iy, dx), R(0.0)); tz = wp.max(dist_to_boundary(p[2], u[2], iz, dx), R(0.0))
                tb = wp.min(tx, wp.min(ty, tz))
                rng = table_lookup(rng_tab, m, E, lE0, inv_dlE, nb) / dens
                sig = table_lookup(sig_tab, m, E, lE0, inv_dlE, nb) * dens
                s_nuc = -wp.log(u01(r4[0])) / sig
                s = wp.min(rng, R(0.25) * rng + R(0.05))
                cross = int(0)
                if tb <= s:
                    s = tb; cross = 1
                nuc = int(0)
                if s_nuc < s:
                    s = s_nuc; nuc = 1; cross = 0
                r_old = rng * dens
                r_new = r_old - s * dens
                e_old_t = table_lookup(irng_tab, 0, r_old, lR0, inv_dlR, nb)
                e_new = R(0.0)
                if r_new > R(1.0e-3):
                    e_new = E - (e_old_t - table_lookup(irng_tab, 0, r_new, lR0, inv_dlR, nb))
                e_new = wp.max(e_new, R(0.0))
                # Box-Muller: straggling + 2 MCS gaussians
                g_r = wp.sqrt(R(-2.0) * wp.log(u01(r4[1])))
                g1 = g_r * wp.cos(R(6.283185307) * u01(r4[2]))
                g2 = g_r * wp.sin(R(6.283185307) * u01(r4[2]))
                dE = E - e_new
                sigE = R(0.03) * wp.sqrt(dE * s)
                dE = wp.clamp(dE + sigE * g1 * R(0.1), R(0.0), E)
                e_new = E - dE
                th0 = R(13.6) / (E + R(1.0e-3)) * wp.sqrt(s * dens / R(361.0)) * R(1.0)
                theta = th0 * wp.abs(g2)
                pmid = p + u * (R(0.5) * s)
                p = p + u * s
                if cross == 1:
                    if tb == tx:
                        if u[0] > R(0.0):
                            ix += 1
                        else:
                            ix -= 1
                    elif tb == ty:
                        if u[1] > R(0.0):
                            iy += 1
                        else:
                            iy -= 1
                    else:
                        if u[2] > R(0.0):
                            iz += 1
                        else:
                            iz -= 1
                u = rotate(u, wp.cos(theta), R(6.283185307) * u01(r4[3]))
                # scoring on 2 mm grid at step midpoint
                jx = int(pmid[0] * R(0.5)); jy = int(pmid[1] * R(0.5)); jz = int(pmid[2] * R(0.5))
                if jx >= 0 and jy >= 0 and jz >= 0 and jx < nx / 2 and jy < ny / 2 and jz < nz / 2:
                    v = (jx * (ny / 2) + jy) * (nz / 2) + jz
                    let = dE / wp.max(s, R(1.0e-6))
                    wp.atomic_add(dose, b, v, wp.float32(dE))
                    wp.atomic_add(letn, b, v, wp.float32(dE * let))
                E = e_new
                if nuc == 1:
                    esec = R(0.4) * E * u01(r5[0]); eloc = R(0.1) * E
                    E = E - esec - eloc
                    if top < stack_n:
                        sx[top] = p[0]; sy[top] = p[1]; sz[top] = p[2]; sE[top] = esec
                        su = rotate(u, R(0.9), R(6.283185307) * u01(r5[1]))
                        sux[top] = su[0]; suy[top] = su[1]; suz[top] = su[2]
                        top += 1
                    else:
                        overflow += 1
                if E < R(0.5):
                    alive = 0
                nstep += 1
                wp.atomic_add(counters, 4, wp.int64(1))
            if alive == 0 and top > 0:
                top -= 1
                p = V3(sx[top], sy[top], sz[top]); E = sE[top]
                u = V3(sux[top], suy[top], suz[top])
                pid += 1; ctr = wp.uint32(0); alive = 1; nstep = 0
                ix = int(wp.floor(p[0])); iy = int(wp.floor(p[1])); iz = int(wp.floor(p[2]))
        wp.atomic_add(counters, 0, wp.int64(pid + 1))
        if overflow > 0:
            wp.atomic_add(counters, 3, wp.int64(overflow))
    return transport

def tables(real_np, nb=512, nmat=4):
    E = np.geomspace(0.1, 400.0, nb)
    rng = np.stack([0.0225 * E ** 1.77 * (1.0 + 0.05 * m) for m in range(nmat)])
    Rg = np.geomspace(1e-3, 0.0225 * 400 ** 1.77 * 1.2, nb)
    irng = np.stack([(Rg / 0.0225) ** (1 / 1.77) for m in range(nmat)])
    sig = np.stack([np.full(nb, 0.0012) * (1 + 0.1 * m) for m in range(nmat)])
    lE0 = math.log(0.1); inv_dlE = (nb - 1) / (math.log(400.0) - math.log(0.1))
    lR0 = math.log(1e-3); inv_dlR = (nb - 1) / (math.log(Rg[-1]) - math.log(1e-3))
    return [a.astype(real_np) for a in (rng, irng, sig)], (lE0, inv_dlE, lR0, inv_dlR, nb)
