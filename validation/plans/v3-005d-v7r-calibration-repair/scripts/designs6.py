"""EXPLORATORY (planner, Revision 2): (B) reduced R for 5400/0.005; power windows of the kappa=1 paired coverage;
P(<=109/3200) at the exact frozen FA values. Produces v7r-data/designs6.out."""
import sys; S=sys.argv[1]; sys.argv=[None,S]
exec(open(S+"/v7r-scripts/designs2.py").read().split("designs = [")[0])
from fastbin import accept_set, kmax, kmin_pow
from fastbin import cdf
for R in (800, 1600):
    _, _, _, t1 = run(5400, 0.005, 10**6, 10**6, R, 1.0, single_shrink=False, draws=200)
    _, _, _, t32 = run(5400, 0.005, 10**6, 10**6, R, 0.32, single_shrink=True, draws=200)
    print(f"(B) 5400/0.005 R {R}: P(all) r=1 bound {t1:.3f} r=.32 {t32:.3f}")
for p in (0.0374, 0.0383): print(f"P(<=109/3200 | {p}) = {cdf(3200, 109, p):.3f}")
kp = kmin_pow(3200)
for npair in (5400, 7200):
    nE = npair * 3 // 2; sh = int(round(0.0019 * nE))
    for c in (0.674, 0.676, 0.678, 0.679, 0.6827, 0.687, 0.689, 0.690, 0.692):
        pw = oc(npair, 0.005, c) * oc(nE, 0.005, c, sh)
        print(f"{npair}/0.005 kappa=1 coverage c={c} (paired=single): paired OC {oc(npair, 0.005, c):.3f} joint power {pw:.3f} P(power clause) {1 - cdf(3200, kp - 1, pw):.3f}")
