"""Planning estimate (revision 2): Fermi-Eyges core-only beam sensitivity and the expected V6-sens-F screen
outcome (coverage expectation of design 6). Uses no ionmc dose output.

Part 1 (beam parameters and analytic MCS only):
  * sigma^2(z) = sx^2 + 2 rho sx st z + st^2 z^2 + MCS(z).
  * MCS(z) = int T_dM(z') (z - z')^2 dz', with Gottschalk's differential Moliere scattering power
    T_dM = f_dM(pv, p1v1) (Es/pv)^2 / X_S, Es = 15.0 MeV, X_S = 46.88 cm (water).
  * pv comes from a Bragg-Kleeman residual range scaled to R(E0) = 21.321 cm.
  * The core is a 2-D Gaussian averaged over the 0.4 cm sphere (the V6 observable). The output is
    R_core - 1 per variant.
  * |R_total - 1| <= |R_core - 1| holds ONLY if the halo itself is beam-insensitive. This is an
    assumption; the halo is a convolution with the beam and is not exactly invariant.

Part 2 (adds the measured table and the precision column of se_table.json; no model mean):
  * Core dose at (r, z) = D_meas(0, z) x Gbar(r) / Gbar(0). This assumes the axis dose is all core,
    which over-states the core slightly.
  * Screen metric of design D3 = |dD_core| / D_meas(r, z) + 2 sigma' / sqrt(k), where sigma' is the
    paired-difference SE over D_meas per 1e6 (st-hi proxy for every variant).
  * Counts the expected H_g, radii and sign-test depths on H_F after the aura exclusion (D7).

Usage: python3 -I core_sens_estimate.py   (reads se_table.json in this directory)
"""

import json
import math
from pathlib import Path

import numpy as np

trapz = getattr(np, "trapezoid", None) or np.trapz  # NumPy 2 renamed trapz
MP = 938.272
E0, R0 = 177.84075737, 21.321
ALPHA_BK, P_BK = 0.0022, 1.77
SCALE = R0 / (ALPHA_BK * E0**P_BK)
ES, XS = 15.0, 46.88


def energy(z):
    rres = np.clip(R0 - z, 1e-6, None) / SCALE
    return (rres / ALPHA_BK) ** (1.0 / P_BK)


def pv(e):
    return e * (e + 2 * MP) / (e + MP)


PV1 = pv(E0)


def t_dm(z):
    x = pv(energy(z))
    l1 = np.log10(np.clip(1.0 - (x / PV1) ** 2, 1e-12, None))
    lp = np.log10(x)
    f = 0.5244 + 0.1975 * l1 + 0.2320 * lp - 0.0098 * lp * l1
    return f * (ES / x) ** 2 / XS


def mcs_var(z, n=4000):
    if z <= 0:
        return 0.0
    zz = np.linspace(0.0, min(z, R0 - 1e-3), n)
    return float(trapz(t_dm(zz) * (z - zz) ** 2, zz))


NOM = {"sx": 0.5948, "st": 0.03155, "rho": -0.5206011525150516}
VAR = {"sx-lo": ("sx", 0.8), "sx-hi": ("sx", 1.2), "st-lo": ("st", 0.8), "st-hi": ("st", 1.2),
       "rho-lo": ("rho", 0.8), "rho-hi": ("rho", 1.2)}


def sig2(z, b):
    return b["sx"] ** 2 + 2 * b["rho"] * b["sx"] * b["st"] * z + (b["st"] * z) ** 2 + mcs_var(z)


# sphere (radius 0.2 cm) sample points in the transverse plane, weighted by chord length
_g = (np.arange(24) + 0.5) / 24 * 0.4 - 0.2
_X, _Y = np.meshgrid(_g, _g)
_IN = _X**2 + _Y**2 <= 0.04
_W = np.sqrt(np.clip(0.04 - _X**2 - _Y**2, 0, None))[_IN]
_DX, _DY = _X[_IN], _Y[_IN]


def gbar(r, s2):
    rr2 = (r + _DX) ** 2 + _DY**2
    return float(np.sum(_W * np.exp(-rr2 / (2 * s2)) / (2 * np.pi * s2)) / np.sum(_W))


def variants(z):
    out = {}
    for k, (key, f) in VAR.items():
        b = dict(NOM)
        b[key] = NOM[key] * f
        out[k] = sig2(z, b)
    return out


def part1():
    print("Part 1: R_core - 1 (sphere-averaged), X_S = 46.88 cm")
    print("z_cm  r   sigma_cm " + " ".join(f"{k:>8s}" for k in VAR))
    for z in (1.533, 5.0, 10.0, 15.0, 18.0, 20.0, 21.0):
        s2n = sig2(z, NOM)
        sv = variants(z)
        for r in (0.0, 1.0, 2.0, 3.0):
            g0 = gbar(r, s2n)
            print(f"{z:5.3g} {r:3.1f} {math.sqrt(s2n):7.3f} " + " ".join(f"{gbar(r, sv[k]) / g0 - 1:+8.3f}" for k in VAR))


def part2(k_shards=8):
    d = json.loads((Path(__file__).resolve().parent / "se_table.json").read_text())
    pts = [p for p in d["points"] if not p["missing"]]
    dmax = d["data"]["r0_max_measured_MeV_per_g"]
    r0 = {p["z_cm"]: p["measured_fraction_of_r0_max"] * dmax for p in pts if p["group"] == "r0"}
    zs0 = sorted(r0)
    hf = [p for p in pts if p["group"] == "halo" and p["z_cm"] <= 21.6 + 1e-9
          and p["measured_fraction_of_r0_max"] >= 1e-4 and p["aura_share_of_measured"] <= 0.05]
    kept = []
    for p in hf:
        z, r = p["z_cm"], p["r_cm"]
        zr = min(zs0, key=lambda q: abs(q - z))
        s2n = sig2(z, NOM)
        sv = variants(z)
        core_nom = r0[zr] * gbar(r, s2n) / gbar(0.0, s2n)
        dmeas = p["measured_fraction_of_r0_max"] * dmax
        sig = p["sens_sigma_diff_over_meas_1e6_st_hi"] / math.sqrt(k_shards)
        worst = max(abs(core_nom * (gbar(r, sv[v]) / gbar(r, s2n) - 1)) / dmeas for v in VAR) + 2 * sig
        if worst < 0.15:
            kept.append(p)
    depths = {p["z_cm"] for p in kept if p["r_cm"] <= 3} & {p["z_cm"] for p in kept if p["r_cm"] >= 4}
    by_r = {int(rr): sum(1 for p in kept if p["r_cm"] == rr) for rr in (1, 2, 3, 4, 5)}
    print(f"Part 2: H_F after aura exclusion {len(hf)} points; expected H_g {len(kept)} points, by r {by_r}, "
          f"radii {sum(1 for v in by_r.values() if v)}, sign-test depths {len(depths)}")
    r1 = sorted(p["z_cm"] for p in kept if p["r_cm"] == 1)
    print(f"  r = 1 cm kept at z = {r1[0] if r1 else None} .. {r1[-1] if r1 else None} cm")


if __name__ == "__main__":
    part1()
    part2()
