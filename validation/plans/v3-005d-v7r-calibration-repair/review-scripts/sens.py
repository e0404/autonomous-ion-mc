"""Reviewer: sensitivity of the recommended 5400/0.005 design to kappa error and to the surrogate's nominal coverage."""
import sys, math
sys.path.insert(0, sys.argv[1])
from oc import oc, binom_cdf, acc, LO, HI
se = math.sqrt(LO*(1-LO)/1e6)
for m in (0, 1, 2, 3, 5):
    d = m*se
    print(f"kappa err {m} SE ({d:.5f} inside): FA paired LO {oc(5400,.005,LO+d):.4f} HI {oc(5400,.005,HI-d):.4f};"
          f" P(<=109) {binom_cdf(3200,109,oc(5400,.005,LO+d)):.4f}")
# where would FA reach 109/3200-type trouble (rate ~0.03) and joint (rate ~0.0125)?
for target in (0.0125, 0.03):
    lo, hi = LO, 0.6827
    for _ in range(50):
        mid = (lo+hi)/2
        if oc(5400,.005,mid) < target: lo = mid
        else: hi = mid
    print(f"FA {target} reached at c = {lo:.5f} (= {lo-LO:.5f} inside, {(lo-LO)/se:.1f} SE_kappa(1e6), {(lo-LO)/math.sqrt(LO*(1-LO)/2e4):.2f} SE_kappa(2e4))")
# nominal-coverage window for power: paired OC x single(model) >= 0.92 (P(>=2909/3200) ~0.99)
print("paired OC vs nominal coverage c0 (power side)")
for c in (0.676, 0.678, 0.679, 0.680, 0.681, 0.6827, 0.685, 0.687, 0.688, 0.689, 0.690):
    print(f"  c0={c}: 5400/.005 {oc(5400,.005,c):.4f}  3600/.04 {oc(3600,.04,c):.4f}  7200/.005 {oc(7200,.005,c):.4f}")
# single gate model at nE=8100 with box shrink 0.0019*nE per side
lo, hi = acc(8100, .005); sh = round(0.0019*8100)
import numpy as np
from oc import pmf
for cs in (0.6808, 0.6866, 0.676, 0.690):
    p = pmf(8100, cs)
    print(f"single nE=8100 a=.005 c_s={cs}: OC no shrink {p[lo:hi+1].sum():.4f}, shrink {sh}: {p[lo+sh:hi-sh+1].sum():.4f}")
