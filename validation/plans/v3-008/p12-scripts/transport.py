import math

import numpy as np

from ionmc.materials import water
from ionmc.physics import projectiles as P
from ionmc.physics import stopping as S

W = water(78.0)
MP = P.PROTON.mass_mev


def ziegler_q(E, ion, vF):
    er = E * ion.a * MP / ion.mass_mev * 1000.0  # keV per proton mass
    eF = 25.0 * vF**2
    v1sq = er / eF
    vFsq = vF**2
    z = ion.z
    z13 = z ** (1 / 3)
    z23 = z13**2
    y = np.where(
        v1sq > 1,
        vF * np.sqrt(v1sq) * (1 + 0.2 / v1sq) / z23,
        0.692308 * vF * (1 + 0.666666 * v1sq + v1sq**2 / 15) / z23,
    )
    y3 = y**0.3
    q = 1 - np.exp(0.803 * y3 - 1.3167 * y3 * y3 - 0.38157 * y - 0.008983 * y * y)
    tq = 7.6 - np.log(er)
    sq = 1 + (0.18 + 0.0015 * z) * np.exp(-tq * tq) / z**2
    lam = 10 * vF * (1 - q) ** (2 / 3) / (z13 * (6 + q))
    xx = (0.5 / q - 0.5) * np.log(1 + lam**2) / vFsq
    return z * q, z * q * np.sqrt(
        sq * (1 + xx)
    )  # bare q; q with charge correction (G4 multiplies q^2 by chargeCorrection)


def sg_q(E, ion, zt=7.22):
    vv0 = np.sqrt(E * 1000 / 25.0 * MP * ion.a / ion.mass_mev)
    zp = ion.z
    x = (vv0 * zp**-0.52 * zt ** (-0.019 * zp**-0.52 * vv0) / 1.68) ** (1 + 1.8 / zp)
    return zp * (12 * x + x**4) / (0.07 / x + 6 + 0.3 * x**0.5 + 10.37 * x + x**4)


E = np.array([2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0])
for ion in (P.CARBON12, P.OXYGEN16):
    b = np.sqrt(S.kinematics(E, ion)[0])
    qb = S.pierce_blann_charge(b, ion) / ion.z
    print(
        ion.name,
        " E  q/z:Barkas125  Ziegler(vF=1) Ziegler+corr(vF=1)  Ziegler+corr(vF=.94)  SG2001(mean)",
    )
    for i, e in enumerate(E):
        z1, z1c = ziegler_q(np.array([e]), ion, 1.0)
        _, z2c = ziegler_q(np.array([e]), ion, 0.94)
        print(
            f"{e:5.0f} {qb[i]:.4f} {z1[0] / ion.z:.4f} {z1c[0] / ion.z:.4f} {z2c[0] / ion.z:.4f} {sg_q(e, ion) / ion.z:.4f}"
        )


# transport
def grid(e0):
    return np.exp(np.linspace(np.log(2.0), np.log(e0), 6000))


def residual(e0, ion, delta):
    e = grid(e0)
    s = S.bethe_mass_stopping(e, ion, W) * (
        1 + delta(e)
    )  # MeV cm2/g per nucleon-energy? check units
    se = s / ion.a  # dE_u/dx in MeV/u per g/cm2
    inv = 1 / se
    r = (
        np.concatenate([[0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(e))]) * 10
    )  # mm from 2 MeV/u
    return e, s, r


for ion, e0, d2 in (
    (P.CARBON12, 290, 0.10),
    (P.OXYGEN16, 350, 0.135),
    (P.CARBON12, 80, 0.10),
    (P.OXYGEN16, 100, 0.135),
):
    z = lambda e: 0 * e
    e, s, r = residual(e0, ion, z)
    R = r[-1]
    i20 = np.searchsorted(e, 20)
    i50 = np.searchsorted(e, 50)
    i10 = np.searchsorted(e, 10)
    print(
        f"{ion.name} {e0}: R(from 2)={R:.2f} mm; R(10)={r[i10]:.3f} R(20)={r[i20]:.3f} R(50)={r[i50]:.3f} mm; frac<20 {100 * r[i20] / R:.2f}% <50 {100 * r[i50] / R:.2f}%"
    )
    # observed band profile (band averages from the results note)
    bands = {
        6: [(2, 5, 0.086), (5, 10, 0.043), (10, 20, 0.018), (20, 50, 0.003)],
        8: [(2, 5, 0.122), (5, 10, 0.075), (10, 20, 0.039), (20, 50, 0.008)],
    }[ion.z]
    prof = lambda e, b=bands: sum(np.where((e >= lo) & (e < hi), d, 0) for lo, hi, d in b)
    for name, dl in (
        ("+5%<20", lambda e: np.where(e < 20, 0.05, 0)),
        ("+10%<20", lambda e: np.where(e < 20, 0.10, 0)),
        ("G4-smooth d2*2/E", lambda e, d2=d2: d2 * 2 / e),
        ("band-avg profile", prof),
    ):
        e2, s2, r2 = residual(e0, ion, dl)
        # depth dose: energy deposited vs depth with Gaussian straggling (Bohr), point deposit at 2 MeV/u
        out = []
        for ee, ss, rr in ((e, s, r), (e2, s2, r2)):
            depth = rr[-1] - rr  # depth where ion has energy ee
            g = np.sqrt(np.maximum(S.kinematics(ee, ion)[0], 1e-12))
            b2 = g**2
            om = (
                0.1569 * ion.z**2 * W.z_over_a * (1 - b2 / 2) / (1 - b2)
            )  # MeV^2 cm2/g (total energy)
            st = ss  # total-energy mass stopping MeV cm2/g
            # sigma_R^2(E) from E down to 2 MeV/u: int om/S^3 dE_total
            dEt = np.diff(ee) * ion.a
            f = om / st**3
            sig2 = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * dEt)]) * 100  # mm^2
            sR = math.sqrt(sig2[-1])
            x = np.arange(0, rr[-1] + 3, 0.02)
            # deposit profile: stopping at depth smeared with straggling at that depth
            dep = np.zeros_like(x)
            dz = np.diff(depth)  # negative
            mid = 0.5 * (depth[1:] + depth[:-1])
            smid = 0.5 * (st[1:] + st[:-1]) * 0.1
            sigm = np.sqrt(np.maximum(sig2[-1] - 0.5 * (sig2[1:] + sig2[:-1]), 1e-8))
            for k in range(0, len(mid), 3):
                w = -np.sum(dz[k : k + 3] * smid[k : k + 3])
                dep += w * np.exp(-0.5 * ((x - mid[k]) / sigm[k]) ** 2) / sigm[k]
            # point deposit 2 MeV/u * A at end
            dep += (
                2.0
                * ion.a
                * np.exp(-0.5 * ((x - rr[-1]) / sR) ** 2)
                / sR
                / math.sqrt(2 * math.pi)
                / math.sqrt(2 * math.pi)
                * math.sqrt(2 * math.pi)
            )
            dep /= math.sqrt(2 * math.pi)
            grp = dep[: len(dep) // 5 * 5].reshape(-1, 5).mean(1)
            ent = dep[(x > 0.5) & (x < 2.5)].mean()
            ip = np.argmax(grp)
            out.append((rr[-1], 0.1 * ip + 0.05, grp[ip] / ent, sR))
        (R0, p0, pe0, sR), (R1, p1, pe1, _) = out
        print(
            f"   {name:18s} dR={R1 - R0:+.3f} mm ({100 * (R1 / R0 - 1):+.3f}%), dz_peak={p1 - p0:+.2f} mm, PE {pe0:.2f}->{pe1:.2f} ({100 * (pe1 / pe0 - 1):+.2f}%), sigmaR~{sR:.3f} mm"
        )
