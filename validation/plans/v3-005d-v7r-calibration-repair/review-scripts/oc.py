"""Reviewer: independent OC of the CP-TOST gate. Acceptance iff P(X>=k|LO)<=a and P(X<=k|HI)<=a
(equivalent to CP lower(k)>=LO and CP upper(k)<=HI by monotonicity of the binomial tails in p)."""
import math, numpy as np
LO, HI = 0.6527, 0.7127
def logpmf(n, p):
    k = np.arange(n + 1)
    lg = np.array([math.lgamma(i + 1) for i in range(n + 1)])
    return lg[n] - lg - lg[::-1] + k * math.log(p) + (n - k) * math.log1p(-p)
def pmf(n, p):
    return np.exp(logpmf(n, p))
def acc(n, a):
    pl, ph = pmf(n, LO), pmf(n, HI)
    sf_lo = np.cumsum(pl[::-1])[::-1]   # P(X>=k|LO)
    cdf_hi = np.cumsum(ph)              # P(X<=k|HI)
    ok = (sf_lo <= a) & (cdf_hi <= a)
    ks = np.nonzero(ok)[0]
    return int(ks.min()), int(ks.max())
def oc(n, a, c):
    lo, hi = acc(n, a); return float(pmf(n, c)[lo:hi + 1].sum())
def binom_cdf(n, k, p):
    return float(pmf(n, p)[:k + 1].sum())
if __name__ == "__main__":
    for n, a in [(3600, .04), (3600, .025), (3600, .02), (5400, .04), (5400, .005), (5400, .0075), (5400, .01), (8100, .005), (8100, .04), (7200, .005)]:
        lo, hi = acc(n, a)
        print(f"n={n} a={a}: accept {lo}-{hi} ({lo/n:.4f}-{hi/n:.4f}); OC LO {oc(n,a,LO):.4f} HI {oc(n,a,HI):.4f}; "
              f"c=.6764 {oc(n,a,.6764):.4f} .6827 {oc(n,a,.6827):.4f} .6840 {oc(n,a,.6840):.4f} .6853 {oc(n,a,.6853):.4f} "
              f".678 {oc(n,a,.678):.4f} .675 {oc(n,a,.675):.4f} .69 {oc(n,a,.69):.4f} .695 {oc(n,a,.695):.4f}")
    for q in (0.0374, 0.0383, 0.0390, 0.0400):
        print(f"P(<=109 of 3200 | {q}) = {binom_cdf(3200,109,q):.4f}; P(<=50|{q})={binom_cdf(3200,50,q):.2e}")
    # real dd976c5 data at frozen design: 2467/3600
    lo, hi = acc(3600, .04); print("2467 in frozen set:", lo <= 2467 <= hi)
    # P(all four paired margins pass clause) with exact kappa
    p1 = binom_cdf(3200,109,oc(3600,.04,LO)); p2 = binom_cdf(3200,109,oc(3600,.04,HI))
    print("P(both paired margins, one surrogate)", p1*p2, "both surrogates", (p1*p2)**2)
    # 5400/.005 clause probabilities
    f = oc(5400,.005,LO); print("5400/.005 FA", f, "P(<=109)", binom_cdf(3200,109,f), "P(<=50)", binom_cdf(3200,50,f))
    # power clause: need >=2909 of 3200
    for pw in (0.914, 0.92, 0.93, 0.954):
        print(f"P(>=2909 | {pw}) = {1-binom_cdf(3200,2908,pw):.4f}")
