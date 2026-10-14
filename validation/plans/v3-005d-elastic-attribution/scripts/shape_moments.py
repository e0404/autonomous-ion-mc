"""Exploratory (not evidence): <1-mu_CM> and theta_rms of black-disk |2J1(qR)/(qR)|^2 with optional Gaussian
surface smearing exp(-q^2 s^2), p+O at 150 MeV (p_CM 514 MeV/c, R 2.68 fm as in the built table)."""
import numpy as np, math
def j1(x):  # Bessel J1 by series/asymptotic via numpy integral representation
    t = np.linspace(0, math.pi, 801)
    return np.trapz(np.cos(t[None, :] - x[:, None] * np.sin(t[None, :])), t, axis=1) / math.pi
k = 514.1 / 197.327
mu = np.linspace(-1, 1, 8001)
q = 2 * k * np.sqrt((1 - mu) / 2)
x = np.maximum(q * 2.68, 1e-9)
bd = (2 * j1(x) / x) ** 2
for R in (2.68, 3.0, 3.3):
    xr = np.maximum(q * R, 1e-9); f = (2 * j1(xr) / xr) ** 2
    for s in (0.0, 0.5, 0.9):
        w = f * np.exp(-(q * s) ** 2)
        w /= np.trapz(w, mu)
        m1 = np.trapz(w * (1 - mu), mu)
        th = np.arccos(mu); rms = math.degrees(math.sqrt(np.trapz(w * th ** 2, mu)))
        print(f"R={R:4.2f} fm s={s:3.1f} fm  <1-mu_CM>={m1:.4f}  theta_rms_CM={rms:5.2f} deg  frac(theta>20deg)={np.trapz(w*(th>math.radians(20)),mu):.3f}")
g = np.exp(-(q * 2.68) ** 2 / 4); g /= np.trapz(g, mu)
print("Gaussian exp(-q^2R^2/4), R 2.68: <1-mu>=%.4f theta_rms=%.2f" % (np.trapz(g * (1 - mu), mu), math.degrees(math.sqrt(np.trapz(g * np.arccos(mu) ** 2, mu)))))
