"""Exact binomial operating characteristics of a Clopper-Pearson TOST gate (scalar paired gate)."""
import math
import numpy as np

LO, HI, NOM = 0.6527, 0.7127, 0.6827

def logpmf(n):
    k = np.arange(n + 1)
    lg = np.array([math.lgamma(i + 1) for i in range(n + 1)])
    return lambda p: lg[n] - lg[k] - lg[n - k] + k * math.log(p) + (n - k) * math.log1p(-p)

def betainc_cdf(n, k, p):  # P(X <= k)
    lp = logpmf(n)(p)
    return float(np.exp(lp[: k + 1]).sum())

def cp_lower(n, h, a):  # one-sided lower bound at level a
    if h == 0: return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        m = 0.5 * (lo + hi)
        # P(X >= h | m) = a
        if 1.0 - betainc_cdf(n, h - 1, m) < a: lo = m
        else: hi = m
    return 0.5 * (lo + hi)

def accept_set(n, a):
    # find hmin with cp_lower >= LO, hmax with cp_upper <= HI (cp_upper(h) = 1 - cp_lower(n, n-h) on complement)
    hs = np.arange(n + 1)
    # bisection on h for monotone bounds
    def first(pred):
        lo, hi = 0, n
        while lo < hi:
            mid = (lo + hi) // 2
            if pred(mid): hi = mid
            else: lo = mid + 1
        return lo
    hmin = first(lambda h: cp_lower(n, h, a) >= LO)
    hmax = first(lambda h: (1.0 - cp_lower(n, n - h, a)) > HI) - 1
    return hmin, hmax

def p_accept(n, hmin, hmax, p):
    pm = np.exp(logpmf(n)(p))
    return float(pm[hmin: hmax + 1].sum())

def p_le(n, k, p):
    return betainc_cdf(n, k, p)

if __name__ == "__main__":
    print("n_pairs alpha  [hmin,hmax]  FA@lo  FA@hi  power@nom  FA@lo+0.001 FA@lo+0.0034 FA@lo+0.0044  P(cal gate-alone clause passes | FA@margin)")
    # calibration clause: CP95 upper <= 0.04 at n_cal=3200 <=> accepted <= 109
    for n in (3600, 7200, 10800):
        for a in (0.04, 0.025, 0.02, 0.015, 0.01):
            hmin, hmax = accept_set(n, a)
            fl, fh = p_accept(n, hmin, hmax, LO), p_accept(n, hmin, hmax, HI)
            pw = p_accept(n, hmin, hmax, NOM)
            s = [p_accept(n, hmin, hmax, LO + e) for e in (0.001, 0.0034, 0.0044)]
            fa = max(fl, fh)
            cal = p_le(3200, 109, fa)
            print(f"{n:6d} {a:5.3f} [{hmin},{hmax}] {fl:.4f} {fh:.4f} {pw:.4f}  " + " ".join(f"{x:.4f}" for x in s) + f"  {cal:.3f}")
