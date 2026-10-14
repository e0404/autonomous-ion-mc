"""EXPLORATORY (planner): (A) 5400-pair designs; (B) reduced R -- NOTE: as executed this block used 7200/0.005
(label in designs5.out is correct); the 5400/0.005 numbers cited in the plan come from designs6.py; (C) exact OCs."""
import sys; S=sys.argv[1]; sys.argv=[None,S]
exec(open(S+"/v7r-scripts/designs2.py").read().split("designs = [")[0])
from fastbin import accept_set, kmax, kmin_pow
print("== (A) MC clauses as written, R=3200, n_kappa 1e6: 5400 pairs")
for a in (0.02, 0.015, 0.01, 0.0075, 0.005):
    fa_p, fa_s, info, t1 = run(5400, a, 10**6, 10**6, 3200, 1.0, single_shrink=False, draws=200)
    _, _, info2, t32 = run(5400, a, 10**6, 10**6, 3200, 0.32, single_shrink=True, draws=200)
    print(f"5400 pairs alpha {a}: FA_p {fa_p:.4f} FA_s(bound) {fa_s:.4f} joint power(model) " + " ".join(f"{k[:4]} {v[0]:.4f}" for k,v in info2.items()) + f" | P(all) r=1 bound {t1:.3f} r=.32 {t32:.3f}")
print("== (B) reduced R with 7200/alpha .005")
for R in (400, 800, 1600):
    fa_p, fa_s, info, t1 = run(7200, 0.005, 10**6, 10**6, R, 1.0, single_shrink=False, draws=200)
    _, _, _, t32 = run(7200, 0.005, 10**6, 10**6, R, 0.32, single_shrink=True, draws=200)
    print(f"R {R}: kmax gate {kmax(R,.04)} joint {kmax(R,.02)} power kmin {kmin_pow(R)} | P(all) r=1 bound {t1:.3f} r=.32 {t32:.3f}")
print("== (C) deterministic paired certification")
for npair in (3600, 5400, 7200):
    nE = npair * 3 // 2
    for a in (0.04, 0.02, 0.01, 0.005):
        fa = max(oc(npair, a, LO), oc(npair, a, HI)); sb = max(oc(nE, a, LO), oc(nE, a, HI)) + 0.01
        pp = [oc(npair, a, c) for c in (0.6827, 0.68396, 0.68305)]; ps = [oc(nE, a, c, int(round(0.0019 * nE))) for c in (0.68664, 0.68083)]
        print(f"pairs {npair} alpha {a}: paired FA exact {fa:.4f}; single FA rigorous bound {sb:.4f}; paired power {pp[0]:.4f}/{pp[1]:.4f}/{pp[2]:.4f}; single power model {ps[0]:.4f}/{ps[1]:.4f}; joint ~ {min(pp[1:])*min(ps):.4f}")
