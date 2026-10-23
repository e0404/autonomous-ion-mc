"""Collect every integer seed recorded in host-run outputs (validation/generated/RUN-*) in 20500000..20600000.
Reads JSON documents (keys containing 'seed') and environment.txt lines with 'seed'. Exploratory, read-only."""
import json, os, re, sys
root = sys.argv[1]
LO, HI = 20500000, 20600000
found = {}  # seed -> set of (run, file, key)
def walk(o, path, run, f):
    if isinstance(o, dict):
        for k, v in o.items():
            if 'seed' in k.lower():
                for s in ints(v):
                    if LO <= s < HI: found.setdefault(s, set()).add((run, os.path.basename(f), k))
            if isinstance(v, (dict, list)) and k not in ('dose_batches',):
                walk(v, path + [k], run, f)
    elif isinstance(o, list):
        for x in o[:5000]:
            if isinstance(x, (dict, list)): walk(x, path, run, f)
def ints(v):
    if isinstance(v, bool): return []
    if isinstance(v, int): return [v]
    if isinstance(v, list): return [x for y in v for x in ints(y)]
    if isinstance(v, dict): return [x for y in v.values() for x in ints(y)]
    return []
for run in sorted(os.listdir(root)):
    rd = os.path.join(root, run)
    for dp, dn, fn in os.walk(rd):
        for n in fn:
            f = os.path.join(dp, n)
            if n.endswith('.json') and os.path.getsize(f) < 50_000_000:
                try: walk(json.load(open(f)), [], run, f)
                except Exception: pass
            elif n in ('environment.txt',) or n.endswith('.txt') and os.path.getsize(f) < 2_000_000:
                for line in open(f, errors='replace'):
                    if 'seed' in line.lower():
                        for m in re.findall(r'\b(20[5]\d{6})\b', line):
                            s = int(m)
                            if LO <= s < HI: found.setdefault(s, set()).add((run, n, 'txt'))
json.dump({str(k): sorted(map(list, v))[:6] for k, v in sorted(found.items())}, open(sys.argv[2], 'w'), indent=0)
print(len(found), "distinct seeds in", LO, HI)
