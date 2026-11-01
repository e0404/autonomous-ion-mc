"""Measured n+natC / n+natO total cross sections read at run time from the ionmc acquisition cache
(EXFOR 13569 = Finlay 1993, 13753 = Abfalterer 2001). No measured value is stored in this file."""
import hashlib, json, os, re, sys
from pathlib import Path
import numpy as np
CACHE = Path(os.environ.get("IONMC_CACHE", "/home/wahln/aiprojects/ion-mc-worktrees/v3-005b/.ionmc-cache"))
ENTRIES = ("13569", "13753")
TARGETS = {"C-12": "6-C-0", "O-16": "8-O-0"}   # natural-element n,TOT datasets used for the A=12/16 comparison
WINDOW = 0.03  # +-3 % in energy

def _load(entry):
    man = CACHE / "ionmc-data" / "manifests" / f"exfor-{entry}.json"
    if not man.is_file():
        sys.exit(f"exfor_tot: acquisition manifest {man} not found - acquire EXFOR {entry} first (set IONMC_CACHE)")
    m = json.loads(man.read_text())
    src = Path(m["source_path"])
    if not src.is_file():
        sys.exit(f"exfor_tot: cached EXFOR {entry} file {src} missing")
    raw = src.read_bytes()
    if hashlib.sha256(raw).hexdigest() != m["sha256"]:
        sys.exit(f"exfor_tot: sha256 mismatch for EXFOR {entry} ({src})")
    return raw.decode("ascii", "replace"), m["sha256"]

def _series(txt, tgt):
    out = []
    for s in re.split(r"\nSUBENT\s+", txt):
        if not re.search(r"REACTION\s+\(" + re.escape(tgt) + r"\(N,TOT\),,SIG\)", s):
            continue
        d = s.split("\nDATA", 1)
        if len(d) < 2:
            continue
        lines = d[1].splitlines()
        heads, units = lines[1].split(), lines[2].split()
        rows = []
        for ln in lines[3:]:
            if ln.startswith("ENDDATA"):
                break
            try:
                rows.append([float(x) for x in ln[:66].split()])
            except ValueError:
                pass
        a = np.array([r for r in rows if len(r) == len(rows[0])])
        ie, idd = heads.index("EN"), heads.index("DATA")
        e = a[:, ie] * (1e-6 if units[ie] == "EV" else 1.0)
        v = a[:, idd] * (1e3 if units[idd].startswith("B") else 1.0)
        out.append((e, v))
    return out

def measured_sigma_tot(energies=(100, 150, 200)):
    """{target: {E: mean over entries of the per-entry mean of points within +-3 % of E}} [mb]."""
    data = {k: _load(k) for k in ENTRIES}
    res = {}
    for name, tgt in TARGETS.items():
        res[name] = {}
        for e0 in energies:
            per = []
            for k, (txt, _) in data.items():
                vals = [v[np.abs(e - e0) <= WINDOW * e0] for e, v in _series(txt, tgt)]
                vals = np.concatenate(vals) if vals else np.array([])
                if vals.size:
                    per.append(float(vals.mean()))
            if not per:
                sys.exit(f"exfor_tot: no {tgt}(N,TOT) points within 3 % of {e0} MeV")
            res[name][e0] = float(np.mean(per))
    return res
