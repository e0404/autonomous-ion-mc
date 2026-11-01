"""Reviewer round 2: V7-R2 transport seeds 20531004 + 16000 + k (k=0..48) vs known used/reserved integer seeds
(conservative envelopes where only ranges are documented)."""
v = set(20531004 + 16000 + k for k in range(49))
sets = {
 "consumed (20441004/20451004/20471004, r0-11, k0-32)": {b + 1000*r + k for b in (20441004, 20451004, 20471004) for r in range(12) for k in range(33)},
 "lv5c 20481004 r0-21 k0-99": {20481004 + 1000*r + k for r in range(22) for k in range(100)},
 "hr5c 20491004 r0-13 k0-99": {20491004 + 1000*r + k for r in range(14) for k in range(100)},
 "rehearsal envelope 20505000-20525999": set(range(20505000, 20526000)),
 "V6 never-run": {b + k for b in (20493004, 20494004, 20501004) for k in range(100)},
 "skipped lv 20511004 r0-21 k0-99": {20511004 + 1000*r + k for r in range(22) for k in range(100)},
 "skipped hr 20521004 r0-13 k0-99": {20521004 + 1000*r + k for r in range(14) for k in range(100)},
 "V6-F candidate B=20551004 r7-9 k0-15": {20551004 + 1000*r + k for r in (7, 8, 9) for k in range(16)},
 "hr 20541004 r1-3 k0-99": {20541004 + 1000*r + k for r in (1, 2, 3) for k in range(100)},
 "hr 20541004 r6 (reserved)": {20541004 + 6000 + k for k in range(100)},
 "lv 20531004 other r (0-15,17-21) k0-99": {20531004 + 1000*r + k for r in list(range(16)) + list(range(17, 22)) for k in range(100)},
 "2098xxxx/2099xxxx": set(range(20980000, 21000000)),
}
print("V7-R2 range", min(v), max(v))
for n, s in sets.items():
    print(f"{n:45s} overlap {len(v & s)}")
