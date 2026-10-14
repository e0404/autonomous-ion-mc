"""Summary of scripts/joint_power.py output against the pre-declared Revision 4 rule (section 3.8): per surrogate the
joint / paired / single pass counts, the one-sided 95 % CP lower bound p_low of the joint rate, P(Bin(3200, p_low) >= 2909),
the Frechet cross-check p_paired + p_single - 1 (on the CP lower bounds), and the runtime. Usage: joint_power_summary.py <jsonl>"""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from fastbin import cp_lower, cdf
rows = [json.loads(x) for x in Path(sys.argv[1]).read_text().splitlines()]
ok_all = True
for s in ("archive_fitted", "rehearsal_block_sums"):
    r = [x for x in rows if x["surrogate"] == s]; n = len(r)
    if not n: print(s, "no rows"); ok_all = False; continue
    j = sum(x["joint"] for x in r); pa = sum(x["paired"] for x in r); si = sum(x["single"] for x in r)
    pl = cp_lower(j, n, 0.05); pp = 1 - cdf(3200, 2908, pl)
    fre = cp_lower(pa, n, 0.05) + cp_lower(si, n, 0.05) - 1
    sel = n >= 400 and j >= 378
    ok_all &= sel
    print(f"{s}: n {n}, joint {j} ({j/n:.4f}), paired {pa} ({pa/n:.4f}), single {si} ({si/n:.4f}); joint CP95 lower {pl:.4f}; "
          f"P(power clause | p_low) {pp:.4f}; Frechet lower (CP bounds) {fre:.4f}; rule (>= 378/400) {'MET' if sel else 'NOT met'}; "
          f"wall {sum(x['wall_s'] for x in r)/3600:.2f} h ({sum(x['wall_s'] for x in r)/n:.1f} s/rep)")
print("selection:", "5400/0.005" if ok_all else "rule not met for 5400 (measure 7200 next)")
