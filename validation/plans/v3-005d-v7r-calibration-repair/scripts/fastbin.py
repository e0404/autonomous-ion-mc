import math, functools
import numpy as np
@functools.lru_cache(None)
def lf(n):
    return np.array([math.lgamma(i + 1) for i in range(n + 1)])
def pmf(n, p):
    k = np.arange(n + 1); L = lf(n)
    return np.exp(L[n] - L - L[::-1] + k * math.log(p) + (n - k) * math.log1p(-p))
def cdf(n, k, p):
    return float(pmf(n, p)[: k + 1].sum()) if k >= 0 else 0.0
def sf(n, k, p):  # P(X >= k)
    return float(pmf(n, p)[k:].sum())
def cp_lower(k, n, a):  # one-sided lower bound: P(X >= k | p) = a
    if k == 0: return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(50):
        m = 0.5 * (lo + hi)
        if sf(n, k, m) < a: lo = m
        else: hi = m
    return 0.5 * (lo + hi)
def cp_upper(k, n, a):
    return 1.0 - cp_lower(n - k, n, a)
@functools.lru_cache(None)
def accept_set(n, a, LO=0.6527, HI=0.7127):
    hmin = next(h for h in range(int(LO * n) - 5, n) if cp_lower(h, n, a) >= LO)
    hmax = next(h for h in range(int(HI * n) + 5, 0, -1) if cp_upper(h, n, a) <= HI)
    return hmin, hmax
def oc(n, a, c, shrink=0):
    hmin, hmax = accept_set(n, a)
    return float(pmf(n, c)[hmin + shrink: hmax - shrink + 1].sum())
@functools.lru_cache(None)
def kmax(R, lim):
    k = int(lim * R)
    while cp_upper(k, R, 0.05) > lim: k -= 1
    while cp_upper(k + 1, R, 0.05) <= lim: k += 1
    return k
@functools.lru_cache(None)
def kmin_pow(R, lim=0.9):
    k = int(0.9 * R)
    while cp_lower(k, R, 0.05) < lim: k += 1
    while cp_lower(k - 1, R, 0.05) >= lim: k -= 1
    return k
