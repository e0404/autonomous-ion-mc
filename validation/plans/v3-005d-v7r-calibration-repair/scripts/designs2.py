"""EXPLORATORY (planner): probability that ALL frozen V7-R calibration clauses pass (both surrogates), per design.
Paired gate: exact binomial CP-TOST on n_pairs. Single gate: binomial CP-TOST on n_E with an effective box shrink of
0.0019*n_E counts per side (fitted to the 10 archived single-gate rates; shrink 0 = the rigorous upper bound on FA).
kappa error: true - target ~ N(0, sqrt(c(1-c)/n_kappa)) per kind, the same draw sign at both margins (rho 0.87).
Joint false acceptance = r * own-gate rate (r = 0.32: largest observed conditional ratio 117/370; r = 1: bound).
Power = OC_paired(c_p0) * OC_single(c_s0) at the measured kappa=1 coverages of each surrogate."""
import sys, math
import numpy as np
sys.path.insert(0, sys.argv[1] + "/v7r-scripts")
from fastbin import oc, cdf, kmax, kmin_pow
LO, HI = 0.6527, 0.7127
SUR = {"archive_fitted": (0.68396, 0.68664), "rehearsal_block_sums": (0.68305, 0.68083)}
rng = np.random.default_rng(20989990)
Z = rng.multivariate_normal([0, 0], [[1, .87], [.87, 1]], size=(300, 2))  # [draw, kind, margin]
def run(npair, a, nk_p, nk_s, R, r, single_shrink=True, draws=300):
    nE = npair * 3 // 2; sh = int(round(0.0019 * nE)) if single_shrink else 0
    kg, kj, kp = kmax(R, 0.04), kmax(R, 0.02), kmin_pow(R)
    sdp = math.sqrt(LO * (1 - LO) / nk_p) if nk_p else 0.0; sds = math.sqrt(LO * (1 - LO) / nk_s) if nk_s else 0.0
    tot = 1.0; info = {}
    for name, (cp0, cs0) in SUR.items():
        pw = oc(npair, a, cp0) * oc(nE, a, cs0, sh)
        ppow = 1.0 - cdf(R, kp - 1, pw)
        acc = 0.0
        for d in range(draws if (sdp or sds) else 1):
            ep, es = Z[d, 0] * sdp, Z[d, 1] * sds
            rates = [oc(npair, a, LO + ep[0]), oc(npair, a, HI + ep[1]), oc(nE, a, LO + es[0], sh), oc(nE, a, HI + es[1], sh)]
            p = 1.0
            for q in rates:
                p *= cdf(R, kg, q) * cdf(R, kj, min(1.0, r * q))
            acc += p
        acc /= (draws if (sdp or sds) else 1)
        info[name] = (pw, ppow, acc * ppow); tot *= acc * ppow
    fa_p = max(oc(npair, a, LO), oc(npair, a, HI)); fa_s = max(oc(nE, a, LO, sh), oc(nE, a, HI, sh))
    return fa_p, fa_s, info, tot
designs = [("D0 frozen", 3600, .04, 20000, 15000, 3200), ("D0b frozen alpha, n_kappa 1e6", 3600, .04, 10**6, 10**6, 3200),
           ("D0c alpha .04, n_kappa 2e7, R 40000", 3600, .04, 2 * 10**7, 2 * 10**7, 40000),
           ("R1 alpha .025, n_kappa 1e6", 3600, .025, 10**6, 10**6, 3200), ("R1' alpha .025, n_kappa 2e5", 3600, .025, 2 * 10**5, 2 * 10**5, 3200),
           ("R2 alpha .02, n_kappa 1e6", 3600, .02, 10**6, 10**6, 3200), ("R3 alpha .015, n_kappa 1e6", 3600, .015, 10**6, 10**6, 3200),
           ("R4 alpha .03, n_kappa 1e6", 3600, .03, 10**6, 10**6, 3200),
           ("R5 2x transport, alpha .02", 7200, .02, 10**6, 10**6, 3200), ("R6 2x transport, alpha .01", 7200, .01, 10**6, 10**6, 3200)]
print("design | FA paired, FA single(model) | per surrogate: joint power, P(power clause), P(all clauses) | P(all, both surrogates)  [r=0.32 | r=1]")
for lab, npair, a, nkp, nks, R in designs:
    out = []
    for r in (0.32, 1.0):
        fa_p, fa_s, info, tot = run(npair, a, nkp, nks, R, r)
        out.append((info, tot))
    i0 = out[0][0]
    print(f"{lab:38s} | {fa_p:.4f}, {fa_s:.4f} | " + "; ".join(f"{k[:4]} pw {v[0]:.4f} Ppow {v[1]:.3f} Pall {v[2]:.3f}" for k, v in i0.items())
          + f" | {out[0][1]:.3f} | r=1: {out[1][1]:.3f}", flush=True)
# real-data paired pass probability for an estimator with true paired coverage c (escaped_neutral observed 0.6853)
for npair in (3600, 7200):
    for a in (.04, .03, .025, .02, .01):
        print(f"paired pass prob n={npair} alpha={a}: " + " ".join(f"c={c}:{oc(npair, a, c):.3f}" for c in (0.6764, 0.6827, 0.6853, 0.69)))
