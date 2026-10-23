"""Reviewer: re-derive P(all calibration clauses) for 5400 pairs / alpha 0.005 (n_E 8100), R in {3200, 1600, 800},
kappa error n=1e6 (and 2e4), single gate with box shrink 0.0019*nE (planner model) or none; r = 1 bound on joint.
Seeds: reviewer diagnostics 20989995 (inside the planner-diagnostic block; exploratory, not evidence)."""
import sys, math, numpy as np
sys.path.insert(0, sys.argv[1])
from oc import pmf, acc, binom_cdf, LO, HI
def kmax(R, lim):
    # largest k with CP95 upper <= lim: P(X<=k | lim) >= 0.05  <=> cp_upper(k) <= lim
    p = pmf(R, lim); c = np.cumsum(p); return int(np.nonzero(c <= 1 - 0.95)[0].max()) if False else int(np.nonzero(c >= 0.05)[0].min()) - 0
def kmax2(R, lim):
    # CP upper(k) <= lim  <=>  P(X <= k | lim) <= 0.05
    c = np.cumsum(pmf(R, lim)); return int(np.nonzero(c <= 0.05)[0].max())
def kmin2(R, lim):
    # CP lower(k) >= lim <=> P(X >= k | lim) <= 0.05
    p = pmf(R, lim); sf = np.cumsum(p[::-1])[::-1]; return int(np.nonzero(sf <= 0.05)[0].min())
def ocs(n, a, c, sh):
    lo, hi = acc(n, a); return float(pmf(n, c)[lo + sh: hi - sh + 1].sum())
rng = np.random.default_rng(20989995)
for R in (3200, 1600, 800):
    kg, kj, kp = kmax2(R, .04), kmax2(R, .02), kmin2(R, .9)
    for nk in (10**6, 2*10**4):
        sdp = math.sqrt(LO*(1-LO)/nk)
        tot = 1.0
        for name, (cp0, cs0) in {"arch": (.68396, .68664), "reh": (.68305, .68083)}.items():
            pw = ocs(5400, .005, cp0, 0) * ocs(8100, .005, cs0, 15)
            ppow = 1 - binom_cdf(R, kp - 1, pw)
            acc_ = 0.0; D = 400
            for _ in range(D):
                z = rng.multivariate_normal([0, 0], [[1, .87], [.87, 1]], size=2) * sdp
                rates = [ocs(5400, .005, LO + z[0, 0], 0), ocs(5400, .005, HI + z[0, 1], 0),
                         ocs(8100, .005, LO + z[1, 0], 0), ocs(8100, .005, HI + z[1, 1], 0)]  # single: no shrink (upper bound)
                p = 1.0
                for q in rates:
                    p *= binom_cdf(R, kg, q) * binom_cdf(R, kj, q)   # r = 1 : joint == own-gate rate (worst case)
                acc_ += p / D
            tot *= acc_ * ppow
            print(f"R={R} n_kappa={nk} {name}: power {pw:.4f} P(power clause) {ppow:.4f} P(margins, r=1) {acc_:.4f}")
        print(f"R={R} n_kappa={nk} limits gate<={kg} joint<={kj} power>={kp}: P(all, both surrogates, r=1, single unshrunk) = {tot:.4f}")
