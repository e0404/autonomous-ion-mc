"""Reviewer spot check: reproduce the recorded rehearsal paired_low kappa (n=20000, SEED_KAPPA+1) with the repo code."""
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np
W = Path("/home/wahln/aiprojects/ion-mc-worktrees/v3-005b"); S = Path(sys.argv[1])
sys.path[:0] = [str(W / "validation/scripts/transport"), str(W / "src")]
import v7r
spec = importlib.util.spec_from_file_location("tv7r", W / "tests/ionmc/test_v7r_coverage.py")
T = importlib.util.module_from_spec(spec); spec.loader.exec_module(T)
st = json.loads((W / ".ionmc-cache/transport-samples/dd976c5/78b/samples/v7r-calibration-status.json").read_text())
fx = W / "tests/ionmc/fixtures/v7r/rehearsal"
print(sorted(p.name for p in fx.iterdir()))
pool = np.concatenate([np.load(fx / f"v7r-s{i}-escaped_neutral.npz")["escaped_neutral_block_sums"] for i in (0, 1)])
sam = T.PooledBlocks(pool)
t = time.time(); k = T.margin_kappas(v7r, sam, "pair", v7r.V7R_REGION[0], T.SEED_KAPPA + 1)
rec = st["surrogates"]["rehearsal_block_sums"]["kappas"]["paired_low"]
print("reproduced", repr(k), "recorded", repr(rec), "equal", k == rec, f"{time.time()-t:.0f}s")
