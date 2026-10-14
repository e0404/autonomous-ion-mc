"""Reviewer: consistency of observed paired-gate counts (3200 reps) with exact OC at directly measured coverage.
Propagates the direct-coverage SE (delta method) and the binomial SE of the observed count."""
import sys, math
sys.path.insert(0, sys.argv[1])
from oc import oc
rows = [  # surrogate/margin, measured true paired coverage, its SE, observed count
 ("arch nominal", .68396, .00104, 3096), ("arch paired_high", .71464, .00101, 67), ("arch paired_low", .65329, .00106, 150),
 ("arch @single_high", .70937, .00102, 309), ("arch @single_low", .64831, .00107, 35),
 ("reh nominal", .68305, .00104, 3073), ("reh paired_high", .71700, .00101, 69), ("reh paired_low", .65814, .00106, 369),
 ("reh @single_high", .70871, .00102, 432), ("reh @single_low", .64729, .00107, 15)]
chi2 = 0
for lab, c, se, k in rows:
    p = oc(3600, .04, c); d = (oc(3600, .04, c + 1e-4) - oc(3600, .04, c - 1e-4)) / 2e-4
    obs = k / 3200; s = math.sqrt((d * se) ** 2 + p * (1 - p) / 3200)
    z = (obs - p) / s; chi2 += z * z
    # implied coverage
    lo, hi = (0.60, c) if (obs > p) == (c > 0.6827) else (c, 0.75)
    lo, hi = 0.60, 0.6827
    if c > 0.6827: lo, hi = 0.6827, 0.76
    for _ in range(60):
        m = (lo + hi) / 2; v = oc(3600, .04, m)
        if (c < 0.6827 and v < obs) or (c > 0.6827 and v > obs): lo = m
        else: hi = m
    print(f"{lab:20s} pred {p:.4f} obs {obs:.4f} z {z:+.2f}  implied c {(lo+hi)/2:.5f} vs direct {c:.5f}")
print("chi2", chi2, "dof", len(rows))
