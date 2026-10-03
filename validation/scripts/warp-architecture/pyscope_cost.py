import time, math, numpy as np, warp as wp
wp.init()
@wp.func
def lerp_log_index(x: float, lx0: float, inv_dl: float, nb: int) -> wp.vec2:
    t = wp.clamp((wp.log(x) - lx0) * inv_dl, 0.0, float(nb - 1) - 1.0e-3)
    i = wp.floor(t)
    return wp.vec2(i, t - i)       # (bin index as float, fraction) -> adapter fetches tab[m,i], tab[m,i+1]
@wp.func
def interp(y0: float, y1: float, f: float) -> float:
    return y0 * (1.0 - f) + y1 * f
@wp.func
def plane_dist(p: float, u: float, i: int, dx: float) -> float:
    if u > 1.0e-12:
        return (float(i + 1) * dx - p) / u
    if u < -1.0e-12:
        return (float(i) * dx - p) / u
    return 1.0e30
@wp.func
def highland_theta0(e_kin: float, mass: float, charge: float, s: float, x0: float) -> float:
    pc = wp.sqrt(e_kin * (e_kin + 2.0 * mass)); beta = pc / (e_kin + mass); t = s / x0
    return 13.6 / (beta * pc) * wp.abs(charge) * wp.sqrt(t) * (1.0 + 0.038 * wp.log(t * charge * charge / (beta * beta)))
@wp.func
def rotate(u: wp.vec3, ct: float, phi: float) -> wp.vec3:
    st = wp.sqrt(wp.max(0.0, 1.0 - ct * ct)); cp = wp.cos(phi); sp = wp.sin(phi)
    den = wp.sqrt(wp.max(1.0 - u[2] * u[2], 0.0))
    if den < 1.0e-6:
        return wp.vec3(st * cp, st * sp, wp.sign(u[2]) * ct)
    n = wp.vec3(u[0] * ct + st * (u[0] * u[2] * cp - u[1] * sp) / den,
                u[1] * ct + st * (u[1] * u[2] * cp + u[0] * sp) / den,
                u[2] * ct - st * cp * den)
    return wp.normalize(n)
tab = np.linspace(1, 2, 512).astype(np.float32).reshape(1, -1)
n = 5000; u = wp.vec3(0.0, 0.0, 1.0)
t0 = time.perf_counter()
for k in range(n):
    iv = lerp_log_index(150.0 - 0.01 * k, -2.3, 60.0, 512); i = int(iv[0])
    y = interp(float(tab[0, i]), float(tab[0, i + 1]), iv[1])
    d = min(plane_dist(60.3, 0.01, 60, 1.0), plane_dist(60.2, 0.02, 60, 1.0), plane_dist(float(k % 7) + 0.3, 0.999, k % 7, 1.0))
    th = highland_theta0(150.0 - 0.01 * k, 938.272, 1.0, 1.0, 360.8)
    u = rotate(u, math.cos(th), 1.0)
dt = (time.perf_counter() - t0) / n
print(f"python-scope shared @wp.func: {dt*1e6:.0f} us per step-equivalent (1 lookup+3 plane dists+highland+rotate)")
