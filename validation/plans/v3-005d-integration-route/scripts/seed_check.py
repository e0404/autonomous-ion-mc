"""Seed-by-seed check of the declared fresh bases against the rehearsal seeds actually used (exploratory).
used = transport seeds recorded in validation/generated/RUN-*/ (used_seeds.py output) + E-shape 20509500+k (declared,
Amendment 15 (b)6 / 14a (c)4; not recorded under a 'seed' key) + the rehearsal fixture-halving seeds 20521000..20521003."""
import json, sys
rec = {int(k) for k in json.load(open(sys.argv[1]))}
used = (rec - {20505000}) | {20509500 + k for k in range(10)}   # 20505000 is a base value, recorded as seed_base
planned_fixture = {20521000 + k for k in range(4)}
print("recorded used seeds:", sorted(rec - {20505000}))
def seeds(base, rs, kmax): return {base + 1000*r + k: (r, k) for r in rs for k in range(kmax)}
for name, base, rs in (("lv 20511004", 20511004, range(0, 21)), ("hr 20521004", 20521004, range(0, 21)),
                       ("lv 20531004", 20531004, range(0, 21)), ("hr 20541004", 20541004, range(0, 21))):
    for kmax in (100, 1000):
        s = seeds(base, rs, kmax)
        hit = sorted((x, s[x]) for x in s if x in used | planned_fixture | {20505000})
        print(f"{name} r 0..20 k<{kmax}: shared with used/fixture seeds: {hit}")
lv = seeds(20511004, range(0, 21), 100); hr = seeds(20521004, range(0, 21), 100)
print("lv20511004 ∩ hr20521004 (cross-base identities, r_lv = r_hr + 10):",
      sorted({(lv[x], hr[x]) for x in lv if x in hr})[:3], "...", len([x for x in lv if x in hr]))
v7r2 = {20511004 + 16000 + k for k in range(49)}
print("V7-R2 transport (r16, k 0..48):", min(v7r2), max(v7r2), "∩ used:", sorted(v7r2 & used), "= hr r_index 6 seeds:", 20521004 + 6000 in v7r2)
v6f = {20531004 + 1000*r + k for r in (7, 8, 9) for k in range(16)}
print("V6-F B=20531004 r7-9:", min(v6f), max(v6f), "∩ lv/hr 2051/2052 (k<100):", len(v6f & (set(lv) | set(hr))), "∩ used:", len(v6f & used))
