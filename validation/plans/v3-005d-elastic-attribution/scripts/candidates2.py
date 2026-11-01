"""Expected probability that the frozen calibration requirements pass, integrating over the kappa-solve error
e ~ N(0, se_k) per margin (paired and single solved separately), for candidate gate designs. Joint false acceptance
taken as 0.32 x the smaller gate rate (largest conditional ratio observed at dd976c5a: 117/370); power at a true
nominal coverage 0.6827 +- d with the observed joint/predicted ratio 0.951/0.961 applied."""
import math, numpy as np
from binom_gate import accept_set, p_accept, p_le, LO, HI, NOM
rng = np.random.default_rng(20990001)
def run(npair, nE, a, se_k, d=0.0, draws=400):
    hp, hs = accept_set(npair, a), accept_set(nE, a)
    fp = lambda c: p_accept(npair, *hp, c); fs = lambda c: p_accept(nE, *hs, c)
    pw = fp(NOM + d) * fs(NOM + d) * (0.951 / 0.961)
    p_pow = 1 - p_le(3200, 2908, pw)
    tot = 0.0
    for _ in range(draws):
        e = rng.normal(0, se_k, 4)
        rates = [fp(LO + e[0]), fp(HI - e[1]), fs(LO + e[2]), fs(HI - e[3])]
        ok = np.prod([p_le(3200, 109, r) for r in rates]) * np.prod([p_le(3200, 50, 0.32 * r) for r in rates])
        tot += ok
    return pw, p_pow, tot / draws * p_pow
for name, npair, nE, a in (("D0 frozen a=0.04", 3600, 5400, 0.04), ("D1 a=0.025", 3600, 5400, 0.025),
                           ("D2 a=0.02", 3600, 5400, 0.02), ("D3 2x reps a=0.02", 7200, 10800, 0.02)):
    for se_k in (0.00337, 0.00048, 0.00034):
        for d in (0.0, 0.003):
            pw, pp, pt = run(npair, nE, a, se_k, d, draws=150)
            print(f"{name:20s} se_kappa={se_k:.5f} nominal offset {d:+.3f}: joint power {pw:.3f}, P(power clause) {pp:.3f}, P(all clauses) {pt:.3f}")
