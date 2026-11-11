"""EXPLORATORY planner measurement (not evidence), Revision 4 (after Codex REVIEW-bb2a8151): the COMPLETE two-gate
joint pass probability of the V7-R2 `escaped_neutral` bootstrap-t gates at the nominal (kappa = 1), measured directly
with the repository's own calibration machinery: `tests/ionmc/test_v7r_coverage.py::run_calibration` (both gates on
the SAME synthetic replicates and the same bootstrap seed base, exactly as the calibration steps do) and
`validation/scripts/transport/v7r.py` (pair_gate, single_gate with the 5000-resample box). The only change is the
gate level: `steps_v5b.V7_ALPHA_TOST` is set to the V7-R2 value 0.005 in this process (both gates read it at call time).
Layout: 48 shards x 225 replicates (10 800 replicates, 5400 pairs; pairs rows i / i + 5400, H = last quarter) for
design 5400, 64 x 225 (14 400) for design 7200.
Surrogates: archive-fitted zero-inflated Gamma (k, pi of record) and the pooled rehearsal block sums (committed fixture,
sha256 5d28a9b3... / d4f7d749..., identical at dd976c5a).
Seeds: roots 20989920/20989921 (archive data/boot), 20989922/20989923 (rehearsal data/boot) for design 5400;
+10 for design 7200; replicate r uses data seed SeedSequence([root_data, r]).generate_state(1)[0] and bootstrap seed
base SeedSequence([root_boot, r]).generate_state(1)[0] (planner block 20989900-20989999, never evidence).
Usage: python joint_power.py <design 5400|7200> <reps per surrogate> <out.jsonl>  (appends one line per replicate;
re-running resumes after the last completed replicate)."""
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np
W = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(W / "validation/scripts/transport"), str(W / "src")]
import steps_v5b, v7r  # noqa: E402
design, reps, out = int(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
steps_v5b.V7_ALPHA_TOST = 0.005
spec = importlib.util.spec_from_file_location("tv7r", W / "tests/ionmc/test_v7r_coverage.py")
T = importlib.util.module_from_spec(spec); spec.loader.exec_module(T)
FIT_K, FIT_PI = 5.070278501518333, 0.21666666666666667  # Amendment 14a (b)1 fit of record (status 8451f463...)
fx = W / "tests/ionmc/fixtures/v7r/rehearsal"
pool = np.concatenate([np.load(fx / f"v7r-s{i}-escaped_neutral.npz")["escaped_neutral_block_sums"] for i in (0, 1)])
sams = {"archive_fitted": T.ZeroInflatedGamma(FIT_K, FIT_PI), "rehearsal_block_sums": T.PooledBlocks(pool)}
nrep = {5400: 10_800, 7200: 14_400}[design]
off = {5400: 0, 7200: 10}[design]
roots = {"archive_fitted": (20989920 + off, 20989921 + off), "rehearsal_block_sums": (20989922 + off, 20989923 + off)}
done = {}
if out.exists():
    for line in out.read_text().splitlines():
        d = json.loads(line); done.setdefault(d["surrogate"], set()).add(d["r"])
with out.open("a") as fh:
    for name, sam in sams.items():
        rd, rb = roots[name]
        for r in range(reps):
            if r in done.get(name, ()):  # resume
                continue
            ds = int(np.random.SeedSequence([rd, r]).generate_state(1)[0]); bs = int(np.random.SeedSequence([rb, r]).generate_state(1)[0])
            t0 = time.perf_counter()
            res = T.run_calibration(v7r, steps_v5b, sam, {"nominal": 1.0}, 1, replicates=nrep, reps_per_shard=225,
                                    data_seed0=ds, boot_seed0=bs)["nominal"][0]
            fh.write(json.dumps({"design": design, "alpha_tost": steps_v5b.V7_ALPHA_TOST, "surrogate": name, "r": r,
                                 "paired": bool(res[0]), "single": bool(res[1]), "joint": bool(res[0] and res[1]),
                                 "data_seed": ds, "boot_seed_base": bs, "wall_s": round(time.perf_counter() - t0, 2)}) + "\n")
            fh.flush()
