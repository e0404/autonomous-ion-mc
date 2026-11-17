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


def charge_table() -> None:
    """Print q/z of the committed effective charge (the Barkas 1963 coefficient-125 form, exposed
    by the callable `pierce_blann_charge`; see the P12 record) against literature parameterisations."""
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
                f"{e:5.0f} {qb[i]:.4f} {z1[0] / ion.z:.4f} {z1c[0] / ion.z:.4f} "
                f"{z2c[0] / ion.z:.4f} {sg_q(e, ion) / ion.z:.4f}"
            )


# transport
def grid(e0):
    return np.exp(np.linspace(np.log(2.0), np.log(e0), 6000))


def residual(e0, ion, delta):
    e = grid(e0)
    s = S.bethe_mass_stopping(e, ion, W) * (1 + delta(e))  # MeV cm2/g (total energy)
    se = s / ion.a  # dE_u/dx in MeV/u per g/cm2
    inv = 1 / se
    r = np.concatenate([[0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(e))]) * 10  # mm
    return e, s, r


E_CUT_PER_U = 2.0  # MeV/u: the energy deposited at a point by the committed transport
POINT_TOL = 2e-3  # the cutoff Gaussian must integrate to E_cut * A within 0.2 %
PROFILE_TOL = 1e-2  # the whole profile must integrate to the deposited energy within 1 %


def cutoff_point_deposit(x, r_end, sigma_r, a):
    """Normalised Gaussian carrying the E_cut * A point deposit, straggled by sigma_r (mm).

    Codex REVIEW-792cbd49 / REVIEW-c43a2e52: the first version normalised this term twice;
    its integral is now asserted separately (``check_point_deposit``)."""
    return (
        E_CUT_PER_U
        * a
        * np.exp(-0.5 * ((x - r_end) / sigma_r) ** 2)
        / (sigma_r * math.sqrt(2 * math.pi))
    )


def check_point_deposit(x, term, a, tol=POINT_TOL):
    """Fail closed unless the point term integrates to E_cut * A within ``tol``."""
    integ = float(np.sum(term) * (x[1] - x[0]))
    if abs(integ / (E_CUT_PER_U * a) - 1) > tol:
        raise RuntimeError(f"cutoff deposit integrates to {integ:.4f}, not {E_CUT_PER_U * a}")
    return integ


def depth_dose(ee, ss, rr, ion):
    """(range mm, peak depth mm, peak-to-entrance, sigma_R mm) of a 0.1 mm-grouped depth-dose
    profile: continuous deposits straggled with the Bohr sigma at depth, plus the cutoff point
    deposit at the end of range."""
    depth = rr[-1] - rr  # depth where the ion has energy ee
    b2 = np.maximum(S.kinematics(ee, ion)[0], 1e-12)
    om = 0.1569 * ion.z**2 * W.z_over_a * (1 - b2 / 2) / (1 - b2)  # MeV^2 cm2/g
    st = ss
    dEt = np.diff(ee) * ion.a
    f = om / st**3
    sig2 = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * dEt)]) * 100  # mm^2
    sR = math.sqrt(sig2[-1])
    x = np.arange(0, rr[-1] + 3, 0.02)
    dep = np.zeros_like(x)
    dz = np.diff(depth)  # negative
    mid = 0.5 * (depth[1:] + depth[:-1])
    smid = 0.5 * (st[1:] + st[:-1]) * 0.1
    sigm = np.sqrt(np.maximum(sig2[-1] - 0.5 * (sig2[1:] + sig2[:-1]), 1e-8))
    for k in range(0, len(mid), 3):
        w = -np.sum(dz[k : k + 3] * smid[k : k + 3])
        dep += w * np.exp(-0.5 * ((x - mid[k]) / sigm[k]) ** 2) / (sigm[k] * math.sqrt(2 * math.pi))
    point = cutoff_point_deposit(x, rr[-1], sR, ion.a)
    check_point_deposit(x, point, ion.a)
    dep += point
    e_expected = float(-np.sum(dz * smid) + E_CUT_PER_U * ion.a)
    integ = float(np.sum(dep) * (x[1] - x[0]))
    if abs(integ / e_expected - 1) > PROFILE_TOL:
        raise RuntimeError(f"energy conservation {integ:.2f} vs {e_expected:.2f} MeV")
    grp = dep[: len(dep) // 5 * 5].reshape(-1, 5).mean(1)
    ent = dep[(x > 0.5) & (x < 2.5)].mean()
    ip = int(np.argmax(grp))
    return rr[-1], 0.1 * ip + 0.05, float(grp[ip] / ent), sR


BANDS = {
    6: [(2, 5, 0.086), (5, 10, 0.043), (10, 20, 0.018), (20, 50, 0.003)],
    8: [(2, 5, 0.122), (5, 10, 0.075), (10, 20, 0.039), (20, 50, 0.008)],
}


def band_profile(ion):
    """Observed band-mean deficit profile of the plan results note (hypothesis), as a stopping
    increase delta(e)."""
    b = BANDS[ion.z]
    return lambda e: sum(np.where((e >= lo) & (e < hi), d, 0) for lo, hi, d in b)


def profiles(ion, d2):
    return (
        ("+5%<20", lambda e: np.where(e < 20, 0.05, 0)),
        ("+10%<20", lambda e: np.where(e < 20, 0.10, 0)),
        ("G4-smooth d2*2/E", lambda e, d2=d2: d2 * 2 / e),
        ("band-avg profile", band_profile(ion)),
    )


def run_case(ion, e0, d2, names=None):
    """Per profile: dR (mm), dR (%), peak shift (mm), PE before/after and change (%), sigma_R."""
    e, s, r = residual(e0, ion, lambda e: 0 * e)
    R = r[-1]
    i10, i20, i50 = (np.searchsorted(e, v) for v in (10, 20, 50))
    info = {"R": R, "R10": r[i10], "R20": r[i20], "R50": r[i50]}
    base = depth_dose(e, s, r, ion)
    out = {}
    for name, dl in profiles(ion, d2):
        if names is not None and name not in names:
            continue
        e2, s2, r2 = residual(e0, ion, dl)
        R1, p1, pe1, _ = depth_dose(e2, s2, r2, ion)
        R0, p0, pe0, sR = base
        out[name] = {
            "dR_mm": R1 - R0,
            "dR_pct": 100 * (R1 / R0 - 1),
            "dz_peak_mm": p1 - p0,
            "pe0": pe0,
            "pe1": pe1,
            "pe_change_pct": 100 * (pe1 / pe0 - 1),
            "sigma_r_mm": sR,
        }
    return info, out


CASES = (
    (P.CARBON12, 290, 0.10),
    (P.OXYGEN16, 350, 0.135),
    (P.CARBON12, 80, 0.10),
    (P.OXYGEN16, 100, 0.135),
)


def main() -> None:
    charge_table()
    for ion, e0, d2 in CASES:
        info, out = run_case(ion, e0, d2)
        R = info["R"]
        print(
            f"{ion.name} {e0}: R(from 2)={R:.2f} mm; R(10)={info['R10']:.3f} R(20)={info['R20']:.3f} "
            f"R(50)={info['R50']:.3f} mm; frac<20 {100 * info['R20'] / R:.2f}% <50 {100 * info['R50'] / R:.2f}%"
        )
        for name, o in out.items():
            print(
                f"   {name:18s} dR={o['dR_mm']:+.3f} mm ({o['dR_pct']:+.3f}%), dz_peak={o['dz_peak_mm']:+.2f} mm, "
                f"PE {o['pe0']:.2f}->{o['pe1']:.2f} ({o['pe_change_pct']:+.2f}%), sigmaR~{o['sigma_r_mm']:.3f} mm"
            )


if __name__ == "__main__":
    main()
