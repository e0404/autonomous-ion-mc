"""EXPLORATORY (planner): 7200/10800-pair variants at R 3200/4800 (produced v7r-data/designs3.out). Usage: python3 -I designs3.py <scratch dir>"""
import sys; S=sys.argv[1]; sys.argv=[None,S]
exec(open(S+"/v7r-scripts/designs2.py").read().split("designs = [")[0])
from fastbin import accept_set, kmax, kmin_pow
for R in (3200, 4800, 6400): print("R", R, "gate kmax", kmax(R,.04), "joint kmax", kmax(R,.02), "power kmin", kmin_pow(R))
for npair in (7200, 10800):
  for a in (0.02, 0.015, 0.0125, 0.01):
    for R in (3200, 4800):
      res=[]
      for r in (0.32, 1.0):
        for shr in (True, False):
          fa_p, fa_s, info, tot = run(npair, a, 10**6, 10**6, R, r, single_shrink=shr, draws=200)
          res.append(tot)
      pw = min(v[0] for v in info.values())
      print(f"pairs {npair} alpha {a} R {R}: FA_p {fa_p:.4f} FA_s(bound) {fa_s:.4f} min joint power {pw:.4f} | P(all) r=.32 model {res[0]:.3f} bound {res[1]:.3f} | r=1 model {res[2]:.3f} bound {res[3]:.3f}", flush=True)
for nk in (2*10**5, 10**6):
    fa_p, fa_s, info, tot = run(7200, 0.01, nk, nk, 3200, 1.0, single_shrink=False, draws=200)
    print("7200/0.01 n_kappa", nk, "P(all) r=1 bound", round(tot,3))
