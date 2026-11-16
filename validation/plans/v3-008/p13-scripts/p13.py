import math
import sys

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[4] / "src"))
from ionmc.materials import water
from ionmc.physics import projectiles as P
from ionmc.physics.stopping import bethe_mass_stopping, kinematics, pierce_blann_charge

W = water(78.0)
F1, F2 = 0.03, 0.048


def S(eu, pr, **kw):
    return float(bethe_mass_stopping(eu, pr, W, **kw)[0])


def q2(eu, pr):
    if pr.z < 2:
        return 1.0
    b2 = kinematics(np.atleast_1d(float(eu)), pr)[0]
    return float(pierce_blann_charge(np.sqrt(b2), pr)[0]) ** 2


def rk(
    e0, L, fn, n=400
):  # energy per u, path in g/cm2 per u-mass units: dE_u/dx = -S/A handled by caller
    h = L / n
    e = e0
    for _ in range(n):
        k1 = fn(e)
        k2 = fn(e - 0.5 * h * k1)
        k3 = fn(e - 0.5 * h * k2)
        k4 = fn(e - h * k3)
        e -= h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return e0 - e


def rich(g1, g2):
    return (F2 * g1 - F1 * g2) / (F2 - F1)


def analyse(pr, A, e0):
    q0 = q2(e0, pr)
    s_ex = lambda e: S(e, pr) / A
    s_nc = lambda e: S(e, pr, barkas=False, bloch=False, mott=False) / A  # q^2 * leading-order part
    out = {}
    for tag, alg in (("frozen_all", lambda e: s_ex(e) * q0 / q2(e, pr)),):
        g = []
        for f in (F1, F2, 0.005, 0.01):
            L = f * e0 / s_ex(e0)
            g.append(rk(e0, L, alg) / rk(e0, L, s_ex) - 1)
        out[tag] = g
    # G4-like split: table part frozen q^2 on leading-order part, correction part (C=S-S_nc) at midpoint, q^2 at mid
    g = []
    for f in (F1, F2, 0.005, 0.01):
        L = f * e0 / s_ex(e0)
        dtab = rk(e0, L, lambda e: s_nc(e) * q0 / q2(e, pr))
        em = e0 - 0.5 * dtab
        dalg = dtab + (s_ex(em) - s_nc(em)) * L
        g.append(dalg / rk(e0, L, s_ex) - 1)
    out["g4split"] = g
    # log slopes
    h = 1e-3
    dq = (math.log(q2(e0 * (1 + h), pr)) - math.log(q2(e0 * (1 - h), pr))) / (2 * h)
    dS = (math.log(s_ex(e0 * (1 + h))) - math.log(s_ex(e0 * (1 - h)))) / (2 * h)
    return q0, dq, dS, out


# ---- G4 table emulation (spline, not-a-knot as ComputeSecDerivative1; Interpolation as G4PhysicsVector.icc L125-147)
def secder(x, y):
    n = len(x) - 1
    u = np.zeros(n)
    s = np.zeros(n + 1)
    u[1] = (y[2] - y[1]) / (x[2] - x[1]) - (y[1] - y[0]) / (x[1] - x[0])
    u[1] = 6 * u[1] * (x[2] - x[1]) / ((x[2] - x[0]) ** 2)
    s[1] = (2 * x[1] - x[0] - x[2]) / (2 * x[2] - x[0] - x[1])
    for i in range(2, n - 1):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * s[i - 1] + 2
        s[i] = (sig - 1) / p
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])
        u[i] = 6 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1] / p
    sig = (x[n - 1] - x[n - 2]) / (x[n] - x[n - 2])
    p = sig * s[n - 3] + 2
    u[n - 1] = (y[n] - y[n - 1]) / (x[n] - x[n - 1]) - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])
    u[n - 1] = 6 * sig * u[n - 1] / (x[n] - x[n - 2]) - (2 * sig - 1) * u[n - 2] / p
    p = (1 + sig) + (2 * sig - 1) * s[n - 2]
    s[n - 1] = u[n - 1] / p
    for k in range(n - 2, 1, -1):
        s[k] *= s[k + 1] - u[k] * (x[k + 1] - x[k - 1]) / (x[k + 1] - x[k])
    s[n] = (s[n - 1] - (1 - sig) * s[n - 2]) / sig
    sig = 1 - ((x[2] - x[1]) / (x[2] - x[0]))
    s[1] *= s[2] - u[1] / (1 - sig)
    s[0] = (s[1] - sig * s[2]) / (1 - sig)
    return s


def interp(x, y, s, e):
    i = int(np.searchsorted(x, e) - 1)
    i = min(max(i, 0), len(x) - 2)
    dl = x[i + 1] - x[i]
    b = (e - x[i]) / dl
    return (
        y[i]
        + b * (y[i + 1] - y[i])
        + b * (b - 1) * ((2 - b) * s[i] + (1 + b) * s[i + 1]) * dl * dl / 6
    )


PR = P.PROTON


def sp(e):  # proxy base stopping (MeV cm2/g) of scaled energy e [MeV]
    if e >= 1.0:
        return S(e, PR)
    return S(1.0, PR) * (e**0.25 if e >= 0.1 else (0.1**0.25) * (e / 0.1) ** 0.45) / 1.0


xs = 1e-4 * 10 ** (np.arange(161) / 20.0)
ys = np.array([sp(e) for e in xs])
s_d = secder(xs, ys)
R = np.zeros(161)
R[0] = 2 * xs[0] / ys[0]
for j in range(1, 161):
    de = (xs[j] - xs[j - 1]) / 100
    e = xs[j] + de / 2
    sm = 0
    for k in range(100):
        e -= de
        sm += de / interp(xs, ys, s_d, e)
    R[j] = R[j - 1] + sm
s_r = secder(xs, R)
s_i = secder(R, xs)


def Rex(e1, e2):  # exact range difference of proxy
    n = 2000
    h = (e2 - e1) / n
    return sum(h / sp(e1 + (k + 0.5) * h) for k in range(n))


def eps_int(es, f):  # relative error of table-branch loss vs exact proxy integral, same path
    dE = f * es
    Lr = Rex(es - dE, es)  # path (scaled units) giving exact loss dE
    Rpre = interp(xs, R, s_r, es)
    ep = interp(R, xs, s_i, Rpre - Lr)
    return (es - ep) / dE - 1


if __name__ == "__main__":
    ions = {"C": (P.CARBON12, 12), "O": (P.OXYGEN16, 16), "He": (P.ALPHA, 4)}
    pts = {
        "C": [2.4, 3, 5, 7, 10, 15, 20, 30, 50],
        "O": [2.4, 3, 5, 7, 10, 15, 20, 30, 50],
        "He": [10, 50],
    }
    for k, (pr, A) in ions.items():
        for e0 in pts[k]:
            q0, dq, dS, o = analyse(pr, A, e0)
            for tag, g in o.items():
                a = g[2] / 0.005
                b = (g[3] / 0.01 - g[2] / 0.005) / 0.005
                print(
                    f"{k} {e0:5} {tag:10} dlnq2={dq:+.4f} dlnS={dS:+.4f} g1={g[0] * 100:+.4f}% g2={g[1] * 100:+.4f}% a={a:+.4f} b={b:+.3f} Rich={rich(g[0], g[1]) * 100:+.5f}% bf1f2={b * F1 * F2 * 100:+.5f}%"
                )
    m = 0
    for k, (pr, A) in ions.items():
        for e0 in pts[k]:
            es = e0 * 1.00728
            e1 = eps_int(es, F1)
            e2 = eps_int(es, F2)
            r = rich(e1, e2)
            m = max(m, abs(r))
            print(f"INT {k} {e0} eps1={e1 * 100:+.5f}% eps2={e2 * 100:+.5f}% Rich={r * 100:+.5f}%")
    sc = [abs(rich(eps_int(e, F1), eps_int(e, F2))) for e in np.geomspace(2.3, 60, 120)]
    print(
        "INT max over points %.5f%%, scan 2.3-60 MeV max %.5f%% p95 %.5f%%"
        % (m * 100, max(sc) * 100, np.percentile(sc, 95) * 100)
    )
    # pure-linear (no spline) comparison for scale
