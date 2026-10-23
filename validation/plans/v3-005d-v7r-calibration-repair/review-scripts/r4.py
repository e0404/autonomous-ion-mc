"""Reviewer round 4: recount joint_power_5400.jsonl, verify seeds, threshold arithmetic (independent CP via tail identity)."""
import json, sys, math, collections
import numpy as np
sys.path.insert(0, sys.argv[2])
from oc import pmf, binom_cdf
rows = [json.loads(l) for l in open(sys.argv[1])]
by = collections.defaultdict(list)
for x in rows: by[x["surrogate"]].append(x)
roots = {"archive_fitted": (20989920, 20989921), "rehearsal_block_sums": (20989922, 20989923)}
for s, r in by.items():
    rs = sorted(x["r"] for x in r)
    assert rs == list(range(400)), (s, len(rs))
    assert all(x["alpha_tost"] == 0.005 and x["design"] == 5400 for x in r)
    assert all(x["joint"] == (x["paired"] and x["single"]) for x in r)
    rd, rb = roots[s]
    ok = all(x["data_seed"] == int(np.random.SeedSequence([rd, x["r"]]).generate_state(1)[0]) and
             x["boot_seed_base"] == int(np.random.SeedSequence([rb, x["r"]]).generate_state(1)[0]) for x in r)
    print(s, "n", len(r), "joint", sum(x["joint"] for x in r), "paired", sum(x["paired"] for x in r),
          "single", sum(x["single"] for x in r), "seeds ok", ok, "wall h", round(sum(x["wall_s"] for x in r)/3600, 3))
def cpl(k, n, a=0.05):
    lo, hi = 0.0, k / n
    for _ in range(60):
        m = (lo + hi) / 2
        sf = float(pmf(n, m)[k:].sum())
        lo, hi = (m, hi) if sf < a else (lo, m)
    return (lo + hi) / 2
for k in (376, 377, 378, 383, 384):
    p = cpl(k, 400); print(f"k={k}: CP95 lower {p:.4f}, P(Bin(3200,p)>=2909) {1 - binom_cdf(3200, 2908, p):.4f}")
for t in (0.95, 0.955, 0.96):
    print(f"P(>=378/400 | {t}) = {1 - binom_cdf(400, 377, t):.3f}")
