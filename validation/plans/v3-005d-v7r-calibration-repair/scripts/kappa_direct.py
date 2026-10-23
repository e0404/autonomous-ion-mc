"""EXPLORATORY planner diagnostic (not evidence). For the two V7-R calibration surrogates of the dd976c5a
calibration: (1) reproduce the recorded kappas with the repository's own margin_kappas (n = 20000, frozen seeds);
(2) estimate the TRUE paired / single coverage at the recorded kappas from a fresh, independent large sample
(seeds 2098990x, planner-diagnostic sub-block), using the repository's bootstrap-t functions unchanged.
Arrays are saved so that all further analysis is offline."""
import importlib.util, json, sys, time, math
from pathlib import Path
import numpy as np
W = Path("/home/wahln/aiprojects/ion-mc-worktrees/v3-005b")
S = Path(sys.argv[1]); NP = int(sys.argv[2]); NS = int(sys.argv[3])
sys.path[:0] = [str(W / "validation/scripts/transport"), str(W / "src")]
import v7r
spec = importlib.util.spec_from_file_location("tv7r", W / "tests/ionmc/test_v7r_coverage.py")
T = importlib.util.module_from_spec(spec); spec.loader.exec_module(T)
st = json.loads((W / ".ionmc-cache/transport-samples/dd976c5/78b/samples/v7r-calibration-status.json").read_text())
fit = st["surrogates"]["archive_fitted"]["fit"]
pool = np.concatenate([np.load(S / f"v7r-data/dd976c5/v7r-s{i}-escaped_neutral.npz")["escaped_neutral_block_sums"] for i in (0, 1)])
sams = {"archive_fitted": T.ZeroInflatedGamma(fit["k"], fit["pi"]), "rehearsal_block_sums": T.PooledBlocks(pool)}
SEEDS = {"archive_fitted": (20989901, 20989902, 20989903), "rehearsal_block_sums": (20989911, 20989912, 20989913)}
out = {}
for name, sam in sams.items():
    rec = st["surrogates"][name]["kappas"]
    t0 = time.time(); rep = T._margin_kappas(v7r, sam); t1 = time.time()
    print(name, "reproduced kappas", {m: (rep[m], rec[m], rep[m] == rec[m]) for m in rec}, f"{t1-t0:.0f}s", flush=True)
    sd, sp, ss = SEEDS[name]
    x = sam(np.random.default_rng(sd), 2 * NP)
    ip = v7r.pair_ingredients(x, sp, reps_per_shard=NP)
    t2 = time.time(); print(name, f"pairs {NP} {t2-t1:.0f}s", flush=True)
    xs = sam(np.random.default_rng(sd + 100), NS * 4 // 3)
    isg = v7r.single_ingredients(xs, ss, reps_per_shard=NS * 4 // 3)
    t3 = time.time(); print(name, f"singles {isg['m_hat'].size} {t3-t2:.0f}s", flush=True)
    np.savez(S / f"v7r-data/direct-{name}.npz", pa=ip["d_hat"] / ip["se_d"], plo=ip["t_lo"], phi=ip["t_hi"], pdeg=ip["n_degenerate"],
             sa=(isg["m_hat"] - sam.mu) / isg["se"], slo=isg["t_lo"], shi=isg["t_hi"], sdeg=isg["n_degenerate"])
    out[name] = {"kappas_recorded": rec, "kappas_reproduced": rep}
(S / "v7r-data/direct-kappas.json").write_text(json.dumps(out, indent=1))
print("done", flush=True)
