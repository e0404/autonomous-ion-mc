"""EXPLORATORY: true coverage at the recorded kappas (fresh 2e5-pair / 1.5e5-single samples) vs the margin targets;
predicted gate acceptance from the exact binomial OC vs the archived counts."""
import json, math, sys
from pathlib import Path
import numpy as np
S = Path(sys.argv[1]); sys.path.insert(0, str(S / "scripts"))
from binom_gate import accept_set, p_accept, LO, HI, NOM
st = json.loads(Path("/home/wahln/aiprojects/ion-mc-worktrees/v3-005b/.ionmc-cache/transport-samples/dd976c5/78b/samples/v7r-calibration-status.json").read_text())
TGT = {"nominal": NOM, "paired_low": LO, "paired_high": HI, "single_low": LO, "single_high": HI}
ocp = accept_set(3600, 0.04); ocs = accept_set(5400, 0.04)
for name in ("archive_fitted", "rehearsal_block_sums"):
    z = np.load(S / f"v7r-data/direct-{name}.npz"); K = st["surrogates"][name]["kappas"]; M = st["surrogates"][name]["margins"]
    print(f"== {name}")
    for kind, a, lo, hi, deg, oc, n in (("paired", z["pa"], z["plo"], z["phi"], z["pdeg"], ocp, 3600), ("single", z["sa"], z["slo"], z["shi"], z["sdeg"], ocs, 5400)):
        fin = np.isfinite(a) & np.isfinite(lo) & np.isfinite(hi)
        print(f"  {kind}: N={a.size}, non-finite quantile rows {int((~fin).sum())}, degenerate rows {int((deg<0).sum())}")
        def cov(k, strict):
            okk = (fin if strict else (deg >= 0))
            with np.errstate(invalid="ignore"):
                h = okk & (k * hi >= a) & (k * lo <= a)
            return h.mean()
        for m, k in K.items():
            c = cov(k, True); c2 = cov(k, False); se = math.sqrt(c * (1 - c) / a.size)
            pred = p_accept(n, *oc, c)
            obs = M[m][kind] / 3200
            own = m.startswith(kind) or m == "nominal"
            print(f"    {m:12s} kappa {k:.5f} target {TGT[m]:.4f} true {c:.5f} (gate-def {c2:.5f}) +- {se:.5f}  shift {c - TGT[m]:+.5f}"
                  f"  binomial OC {pred:.4f}  observed {kind} {obs:.4f} ({M[m][kind]}/3200){'  <- own margin' if own else ''}")
        # accurate kappa at the targets from the big sample
        for t in (LO, HI):
            l, h = 0.05, 20.0
            for _ in range(60):
                mid = 0.5 * (l + h); l, h = (mid, h) if cov(mid, True) < t else (l, mid)
            print(f"    big-sample kappa for {t}: {0.5*(l+h):.5f}")
