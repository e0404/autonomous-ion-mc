"""Implied true paired coverage at the calibration kappas, from the archived acceptance counts (exact binomial OC)."""
from binom_gate import accept_set, p_accept, LO, HI
import math
n = 3600; hmin, hmax = accept_set(n, 0.04)
def solve(rate, side):
    # side 'lo': acceptance increases with c for c near LO; 'hi': acceptance decreases with c near HI
    a, b = (LO - 0.02, LO + 0.03) if side == "lo" else (HI - 0.03, HI + 0.02)
    f = lambda c: p_accept(n, hmin, hmax, c) - rate
    for _ in range(60):
        m = 0.5 * (a + b)
        if (f(m) < 0) == (side == "lo"): a = m
        else: b = m
    return 0.5 * (a + b)
for sur, pl, ph in (("rehearsal_block_sums", 369, 69), ("archive_fitted", 150, 67)):
    for name, k, side, marg in (("paired_low", pl, "lo", LO), ("paired_high", ph, "hi", HI)):
        r = k / 3200; se = math.sqrt(r * (1 - r) / 3200)
        c = solve(r, side); c1 = solve(max(r - 1.96 * se, 1e-4), side); c2 = solve(r + 1.96 * se, side)
        print(f"{sur:22s} {name:12s} accepted {k:4d}/3200 -> implied true paired coverage {c:.5f} (margin {marg}; shift {c - marg:+.5f}; 95% [{min(c1,c2) - marg:+.5f},{max(c1,c2) - marg:+.5f}])")
print("kappa-solve MC coverage SE at n=20000:", math.sqrt(LO*(1-LO)/20000), math.sqrt(HI*(1-HI)/20000))
