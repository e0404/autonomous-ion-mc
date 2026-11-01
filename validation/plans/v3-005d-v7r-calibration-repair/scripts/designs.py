"""EXPLORATORY design evaluation for the repaired V7-R calibration (planner, not evidence).
Model: paired gate = exact binomial CP-TOST (n_pairs); single gate bounded by the binomial CP-TOST on n_E with an
effective box shift d_box fitted to the observed nominal single power (0.985/0.989 at alpha 0.04, n_E 5400);
kappa errors e ~ N(0, sd_k) per (kind), low/high correlated rho = 0.87 (one sample per kind);
joint false acceptance = r x own-gate rate (r = 0.32 = max observed conditional ratio; r = 1 = worst case bound).
Clause thresholds from exact CP95 bounds at R repetitions. Power = paired OC x single OC at the measured kappa=1
coverage (c_p0, c_s0)."""
import math, sys
import numpy as np
sys.path.insert(0, sys.argv[1] + "/scripts")
from binom_gate import accept_set, p_accept, p_le, cp_lower, LO, HI, NOM
rng = np.random.default_rng(20989990)

def cp_up(k, n, a=0.05): return 1.0 - cp_lower(n, n - k, a) if k < n else 1.0
def kmax(R, lim):
    k = 0
    while cp_up(k + 1, R) <= lim: k += 1
    return k
def kmin_pow(R, lim=0.9):
    k = R
    while cp_lower(R, k - 1, 0.05) >= lim: k -= 1
    return k
def single_oc(nE, a, dbox):
    hmin, hmax = accept_set(nE, a)
    pm = lambda c: np.exp(__import__('binom_gate').logpmf(nE)(c))
    def f(c):  # accept iff K - dK >= hmin and K + dK <= hmax, dK = dbox * nE
        dk = int(round(dbox * nE)); p = pm(c)
        return float(p[hmin + dk: hmax - dk + 1].sum())
    return f

c_p0, c_s0 = float(sys.argv[2]), float(sys.argv[3])
# fit d_box: single power at alpha .04 nE 5400 observed ~0.987 (mean of 0.9853, 0.98875)
best = min((abs(single_oc(5400, 0.04, d)(c_s0) - 0.987), d) for d in np.arange(0, 0.006, 0.0002))
dbox = best[1]
print(f"c_p0 {c_p0:.4f} c_s0 {c_s0:.4f} fitted d_box {dbox:.4f}")
rows = []
for npair in (3600, 7200):
    nE = npair * 3 // 2
    for a in (0.04, 0.025, 0.02, 0.015, 0.012, 0.01):
        hp = accept_set(npair, a); fp = lambda c, hp=hp: p_accept(npair, *hp, c); fs = single_oc(nE, a, dbox)
        fa_p = max(fp(LO), fp(HI)); fa_s = max(fs(LO), fs(HI)); pw_p = fp(c_p0); pw_s = fs(c_s0); pw = pw_p * pw_s
        for R in (3200, 6400):
            kg, kj, kp = kmax(R, 0.04), kmax(R, 0.02), kmin_pow(R)
            ppow = 1.0 - p_le(R, kp - 1, pw)
            for sdk in (0.0034, 0.00047, 0.0):
                for r in (0.32, 1.0):
                    tot = 0.0; draws = 200 if sdk > 0 else 1
                    for _ in range(draws):
                        z = rng.multivariate_normal([0, 0, 0, 0], [[1, .87, 0, 0], [.87, 1, 0, 0], [0, 0, 1, .87], [0, 0, .87, 1]]) * sdk
                        rates = [fp(LO + z[0]), fp(HI + z[1]), fs(LO + z[2]), fs(HI + z[3])]  # one kappa sample per kind: true - target has the same sign at both margins
                        ok = 1.0
                        for q in rates:
                            ok *= p_le(R, kg, q) * p_le(R, kj, min(1.0, r * q))
                        tot += ok
                    p_one = tot / draws * ppow
                    rows.append((npair, a, R, sdk, r, fa_p, fa_s, pw_p, pw_s, pw, p_one, p_one ** 2))
print("pairs alpha  R    sd_k    r   | FA_p   FA_s(bound) | power p  s  joint | P(clauses, one surrogate) P(both)")
for row in rows:
    print("%5d %5.3f %5d %.5f %.2f | %.4f %.4f | %.4f %.4f %.4f | %.4f %.4f" % row)
