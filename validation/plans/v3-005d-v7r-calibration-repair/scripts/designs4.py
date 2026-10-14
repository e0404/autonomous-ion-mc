"""EXPLORATORY (planner): 7200 pairs at alpha 0.01/0.0075/0.005 (produced v7r-data/designs4.out)."""
import sys; S=sys.argv[1]; sys.argv=[None,S]
exec(open(S+"/v7r-scripts/designs2.py").read().split("designs = [")[0])
from fastbin import accept_set, kmax, kmin_pow
for a in (0.01, 0.0075, 0.005):
    for nk in (2*10**5, 10**6):
        fa_p, fa_s, info, tot1 = run(7200, a, nk, nk, 3200, 1.0, single_shrink=False, draws=200)
        _, _, _, tot32 = run(7200, a, nk, nk, 3200, 0.32, single_shrink=True, draws=200)
        pw = {k[:4]: round(v[0], 4) for k, v in info.items()}
        print(f"7200 pairs alpha {a} n_kappa {nk}: accept set {accept_set(7200, a)} FA_p {fa_p:.4f} FA_s(bound) {fa_s:.4f} joint power {pw} | P(all) r=1 bound {tot1:.3f} | r=.32 model {tot32:.3f}")
    print("  real-data paired pass prob:", " ".join(f"c={c}:{oc(7200, a, c):.3f}" for c in (0.6764, 0.6827, 0.6853, 0.69, 0.695)))
