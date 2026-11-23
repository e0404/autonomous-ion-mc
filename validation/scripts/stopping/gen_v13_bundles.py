"""Generate the V3-013 TOPAS case bundles of row V-C1 (He 75 / 125 / 200 MeV/u in water): 27 bundles.

Single source of ``validation/reference_cases/topas/he-water-*mevu-v13-idd-em-*`` and of their frozen
manifest ``validation/reference_cases/topas/v13-bundles-manifest.json`` (plan
``validation/plans/v3-013-acceptance.md``, section A.6 row V-C1, section C and section D step 3). It must
not be changed without a dated plan amendment (a change after the seal is a new seal). Nothing is run;
only ``case.json`` and ``input.txt`` are written.

Per case (plan sizing table, rules of ``validation/plans/v3-008/r2-topas-target.md`` section 3): box
length L, IDD bin h_z, total beam energy 4 x E/u. Physics, cuts, source, step cap and scorer
conventions are those of the TD-EM EM-only He bundles of plan v3-008 (``he-water-100mevu-r2-idd-em-part0``):
``g4em-standard_opt4`` only, ``CutForAllParticles`` 0.05 mm, ``CutForElectron`` 1 km, ``EMRangeMax`` 10 GeV,
``Ts/MaxStepNumber`` 1e7, ``b:Ph/ListProcesses`` True, monoenergetic zero-size pencil beam from z = -0.001 mm,
a 120 x 120 mm lateral extent, and the IDD scorer only (no lateral scorers).

Parts (plan section C, family ``topas13``: ``seed = 20610000 + 100*part + case``, case 1 = He 75, 2 = He
125, 3 = He 200 MeV/u): parts 0-3 evidence (2.5e4 primaries each); parts 4-7 frozen repeats (part + 4,
same configuration; drawn only on a run-validity repeat trigger of r2-topas-target section 6); part 19
execution smoke (2000 primaries, never evidence).

Usage::

    python -I validation/scripts/stopping/gen_v13_bundles.py [--out DIR]   # write bundles
    python -I validation/scripts/stopping/gen_v13_bundles.py --check       # regenerate in a temp dir, compare
                                                                           # byte for byte, print the digest
    python -I validation/scripts/stopping/gen_v13_bundles.py --manifest    # write the frozen manifest
"""

import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "validation" / "reference_cases" / "topas"
REL_PREFIX = "validation/reference_cases/topas"
DEFAULT_MANIFEST = DEFAULT_OUT / "v13-bundles-manifest.json"
MANIFEST_SCHEMA = "ionmc-v13-bundles-1"
BUNDLE_NAME_TAG = "-v13-"
BUNDLE_SET_DIGEST = "aa2f24e436f5cfb22b6bc0920367bfe7d2677bb00ea05cd9a05d407cd1466b83"

TOPAS13_BASE = 20610000  # seed = TOPAS13_BASE + 100*part + case (plan section C)
EVIDENCE_PARTS = (0, 1, 2, 3)
REPEAT_OFFSET = 4  # parts 4-7 = frozen repeats of parts 0-3
REPEAT_PARTS = tuple(p + REPEAT_OFFSET for p in EVIDENCE_PARTS)
SMOKE_PART = 19
N_EVIDENCE = 25000  # primaries per evidence (and repeat) part
N_SMOKE = 2000
# case number, MeV/u, planning CSDA range R_plan of the default model (mm, rounded to 0.1 mm),
# box length L (mm), IDD bin h_z (mm); A = 4, z = 2, particle alpha. Plan section A.6 sizing table.
A_NUCLEONS = 4
FIELD_MM = 120
CASES = [
    (1, 75, 46.7, 100, 0.05),
    (2, 125, 115.7, 200, 0.05),
    (3, 200, 261.9, 350, 0.1),
]
EM_MODULES = 'sv:Ph/Default/Modules = 1 "g4em-standard_opt4"'
UNITS = (
    "IDD: TOPAS CSV with '#' header, dose in Gy (run total; divide by histories); converted to "
    "laterally integrated energy per depth (MeV/mm per primary) with the slab mass 120 x 120 mm x bin "
    "depth at 1.0 g/cm3. Fluence 1/mm2 run total. Bins 0-based, z from the entrance face"
)


def fmt(x):
    return f"{x:g}"


def seed_of(part, case):
    return TOPAS13_BASE + 100 * part + case


def name_of(c, tail):
    return f"he-water-{c[1]}mevu-v13-idd-em-{tail}"


def build(c, part, histories, kind):
    """One bundle; kind is "evidence", "repeat" or "smoke"."""
    num, e, r_plan, length, hz = c
    seed = seed_of(part, num)
    total_e = A_NUCLEONS * e
    nz_idd = round(length / hz)
    hdr = {
        "evidence": f"V3-013 plan A.6 row V-C1, case {num}, part {part} (seed = 20610000 + 100*part + case).",
        "repeat": (
            f"V3-013 plan A.6 row V-C1, case {num}, part {part} = FROZEN REPEAT of part {part - REPEAT_OFFSET} "
            "(seed = 20610000 + 100*part + case; drawn only on a run-validity repeat trigger of "
            "r2-topas-target section 6)."
        ),
        "smoke": (
            f"V3-013 plan A.6 row V-C1, case {num}, part {part} = EXECUTION SMOKE (seed = 20610000 + "
            "100*part + case; execution only, 2000 primaries, never evidence)."
        ),
    }[kind]
    inp = [
        f"# Native OpenTOPAS 4.3 input. {hdr} No validation claim.",
        "i:Ts/NumberOfThreads = 1",
        f"i:Ts/Seed = {seed}",
        'b:Ts/PauseBeforeQuit = "False"',
        'b:Ts/ShowCPUTime = "True"',
        "i:Ts/TrackingVerbosity = 0",
        "# Step cap raised from the TOPAS default (1e6); any killed track (stderr 'stuck'/'killed') makes the run invalid (plan v3-008 revision 2, run validity).",
        "i:Ts/MaxStepNumber = 10000000",
        "",
        f"# Geometry: {FIELD_MM} x {FIELD_MM} x {length} mm G4_WATER box, entrance face at z = 0, beam axis = +z through the face centre.",
        's:Ge/World/Material = "Vacuum"',
        "d:Ge/World/HLX = 20 cm",
        "d:Ge/World/HLY = 20 cm",
        "d:Ge/World/HLZ = 60 cm",
        's:Ge/Water/Parent = "World"',
        's:Ge/Water/Type = "TsBox"',
        's:Ge/Water/Material = "G4_WATER"',
        f"d:Ge/Water/HLX = {FIELD_MM // 2} mm",
        f"d:Ge/Water/HLY = {FIELD_MM // 2} mm",
        f"d:Ge/Water/HLZ = {fmt(length / 2)} mm",
        f"d:Ge/Water/TransZ = {fmt(length / 2)} mm",
        's:Ge/Source/Parent = "World"',
        's:Ge/Source/Type = "Group"',
        "d:Ge/Source/TransZ = -0.001 mm",
        "",
        f"# Source: monoenergetic, zero-size, zero-divergence pencil beam along +z. BeamEnergy = TOTAL kinetic energy ({A_NUCLEONS} x {e} MeV/u = {total_e} MeV; P4 convention).",
        's:So/Beam/Type = "Beam"',
        's:So/Beam/Component = "Source"',
        's:So/Beam/BeamParticle = "alpha"',
        f"d:So/Beam/BeamEnergy = {total_e} MeV",
        "u:So/Beam/BeamEnergySpread = 0",
        's:So/Beam/BeamPositionDistribution = "None"',
        's:So/Beam/BeamAngularDistribution = "None"',
        f"i:So/Beam/NumberOfHistoriesInRun = {histories}",
        "",
        "# Physics: EM only, matched electron cut: CutForElectron = 1 km so that delta electrons deposit locally (as ionmc; P4 verbose smoke: e- threshold 9.55 GeV in G4_WATER).",
        's:Ph/Default/Type = "Geant4_Modular"',
        EM_MODULES,
        "d:Ph/Default/CutForAllParticles = 0.05 mm",
        "d:Ph/Default/CutForElectron = 1 km",
        "# Prints the per-particle process list so that the executed module set, g4decay included, is verifiable from the log (Codex REVIEW-6b4731a7).",
        'b:Ph/ListProcesses = "True"',
        "# EM table upper edge above the highest total kinetic energy (default 600 MeV is below He 150-250 MeV/u totals).",
        "d:Ph/Default/EMRangeMax = 10 GeV",
        "",
        f"# Laterally integrated depth dose: {nz_idd} depth bins of {fmt(hz)} mm over the full {FIELD_MM} x {FIELD_MM} mm field.",
        's:Sc/IDD/Quantity = "DoseToMedium"',
        's:Sc/IDD/Component = "Water"',
        "i:Sc/IDD/XBins = 1",
        "i:Sc/IDD/YBins = 1",
        f"i:Sc/IDD/ZBins = {nz_idd}",
        's:Sc/IDD/OutputType = "csv"',
        's:Sc/IDD/OutputFile = "idd_dose"',
        's:Sc/IDD/IfOutputFileAlreadyExists = "Overwrite"',
        'sv:Sc/IDD/Report = 2 "Sum" "Standard_Deviation"',
        "",
    ]
    where = f"V3-013 plan A.6 row V-C1, case {num}, part {part}; seal section D step 3; "
    rationale = {
        "evidence": (
            where + f"case {num} = He {e} MeV/u, evidence part {part} of parts 0-3, "
            "reference side of row V-C1 (TOPAS 4.3 / Geant4 11.4.2, g4em-standard_opt4 only; shared-lineage "
            "transport agreement, never a validation of the stopping physics); "
            "not run before the reviewed seal and before the ionmc samples are committed"
        ),
        "repeat": (
            where + f"case {num} = He {e} MeV/u, frozen repeat of evidence part {part - REPEAT_OFFSET} "
            "(parts 4-7), drawn only on a run-validity repeat trigger of r2-topas-target section 6; "
            "not run before the reviewed seal and before the ionmc samples are committed"
        ),
        "smoke": (
            where + f"case {num} = He {e} MeV/u, execution smoke (part 19, 2000 primaries; "
            "execution only, outputs read for run-validity and timing diagnostics only, never evidence); "
            "not run before the reviewed seal and before the ionmc samples are committed"
        ),
    }[kind]
    case = {
        "input": "input.txt",
        "physics": [
            'Geant4_Modular (OpenTOPAS 4.3 / Geant4 11.4.2), EM ONLY: sv:Ph/Default/Modules = 1 "g4em-standard_opt4"; no hadronic, ion, elastic or decay modules'
        ],
        "source": f"alpha (he {e} MeV/u), BeamEnergy {total_e} MeV total kinetic energy (P4 convention), monoenergetic (BeamEnergySpread 0), zero-size, zero-divergence pencil beam from z=-0.001 mm along +z at the face centre",
        "geometry": f"G4_WATER box {FIELD_MM} x {FIELD_MM} x {length} mm (x,y centred on the beam axis; z from 0 to {length} mm) in a vacuum world 400 x 400 x 1200 mm; L = 1.3 x planning CSDA range ({r_plan} mm, default Bethe model) rounded up to 50 mm",
        "materials": "G4_WATER (Geant4 NIST material, I=78 eV, density 1.0 g/cm3)",
        "cuts": "CutForAllParticles = 0.05 mm; CutForElectron = 1 km (accepted, P4 verbose smoke: e- production threshold 9.55 GeV in G4_WATER, no delta electrons); EMRangeMax = 10 GeV; Ts/MaxStepNumber = 1e7",
        "histories": histories,
        "seeds": [seed],
        "scorers": [
            f"IDD: DoseToMedium, 1x1x{nz_idd} bins (full {FIELD_MM}x{FIELD_MM} mm lateral extent, {fmt(hz)} mm depth bins); Report Sum, Standard_Deviation"
        ],
        "units": UNITS,
        "rationale": rationale,
        "expected_outputs": ["idd_dose.csv"],
        "arguments": [],
    }
    return case, "\n".join(inp)


def write(out, name, case, inp):
    d = out / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "case.json").write_text(json.dumps(case, indent=2, ensure_ascii=False) + "\n")
    (d / "input.txt").write_text(inp)


def generate(out):
    """Write the 27 bundles under `out`; return [(name, seed)]."""
    out = Path(out)
    made = []
    for c in CASES:
        for part in EVIDENCE_PARTS:
            case, inp = build(c, part, N_EVIDENCE, "evidence")
            nm = name_of(c, f"part{part}")
            write(out, nm, case, inp)
            made.append((nm, seed_of(part, c[0])))
        for part in REPEAT_PARTS:
            case, inp = build(c, part, N_EVIDENCE, "repeat")
            nm = name_of(c, f"part{part}")
            write(out, nm, case, inp)
            made.append((nm, seed_of(part, c[0])))
        case, inp = build(c, SMOKE_PART, N_SMOKE, "smoke")
        nm = name_of(c, f"part{SMOKE_PART}-smoke")
        write(out, nm, case, inp)
        made.append((nm, seed_of(SMOKE_PART, c[0])))
    return made


def _bundle_dirs(out):
    out = Path(out)
    return sorted(d for d in out.iterdir() if d.is_dir() and BUNDLE_NAME_TAG in d.name)


def bundle_set_digest(out):
    """sha256 of the sha256sum-style listing (sorted by path, final newline) of all v13 bundle files."""
    lines = []
    for d in _bundle_dirs(out):
        for f in ("case.json", "input.txt"):
            h = hashlib.sha256((d / f).read_bytes()).hexdigest()
            lines.append((f"{REL_PREFIX}/{d.name}/{f}", f"{h}  {REL_PREFIX}/{d.name}/{f}"))
    lines.sort()
    text = "".join(line + "\n" for _, line in lines)
    return hashlib.sha256(text.encode()).hexdigest(), len(lines)


def manifest_text(out):
    """Frozen bundle manifest (JSON text) of the v13 bundles under `out`.

    Per bundle the sha256 of case.json and input.txt, the counts and the bundle-set digest.
    Deterministic: sorted keys, indent 1, trailing newline. Fails if the bundle-set digest is not
    BUNDLE_SET_DIGEST.
    """
    bundles = {
        d.name: {f: hashlib.sha256((d / f).read_bytes()).hexdigest() for f in ("case.json", "input.txt")}
        for d in _bundle_dirs(out)
    }
    digest, nfiles = bundle_set_digest(out)
    if digest != BUNDLE_SET_DIGEST:
        raise ValueError(f"bundle-set digest {digest} != recorded {BUNDLE_SET_DIGEST}")
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "generated_from": "validation/scripts/stopping/gen_v13_bundles.py",
        "bundle_count": len(bundles),
        "file_count": nfiles,
        "bundle_set_digest": digest,
        "bundles": bundles,
    }
    return json.dumps(manifest, sort_keys=True, indent=1) + "\n"


def write_manifest(out, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(manifest_text(out))
    return path


def check(committed, manifest_path=DEFAULT_MANIFEST):
    committed = Path(committed)
    manifest_path = Path(manifest_path)
    with tempfile.TemporaryDirectory() as tmp:
        made = generate(tmp)
        try:
            fresh_manifest = manifest_text(tmp)
        except ValueError as exc:
            fresh_manifest = None
            print("MANIFEST NOT GENERATED:", exc)
        bad = []
        for nm, _ in made:
            for f in ("case.json", "input.txt"):
                p = committed / nm / f
                if not p.is_file() or p.read_bytes() != (Path(tmp) / nm / f).read_bytes():
                    bad.append(f"{nm}/{f}")
        names = {nm for nm, _ in made}
        extra = [d.name for d in _bundle_dirs(committed) if d.name not in names] if committed.is_dir() else []
        digest, nfiles = bundle_set_digest(tmp)
    print(f"generated {len(made)} bundles, {nfiles} files; bundle-set digest {digest}")
    ok = digest == BUNDLE_SET_DIGEST
    print("digest matches recorded value" if ok else f"DIGEST MISMATCH (recorded {BUNDLE_SET_DIGEST})")
    for b in bad:
        print("DIFFERS or MISSING:", b)
    for e in extra:
        print("EXTRA committed v13 bundle:", e)
    manifest_ok = (
        fresh_manifest is not None and manifest_path.is_file() and manifest_path.read_text() == fresh_manifest
    )
    print(
        f"manifest {manifest_path} matches a fresh regeneration"
        if manifest_ok
        else f"MANIFEST MISSING or DIFFERS: {manifest_path}"
    )
    return ok and not bad and not extra and manifest_ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output directory (default: validation/reference_cases/topas)")
    ap.add_argument("--check", action="store_true", help="regenerate into a temp dir and compare byte for byte with --out")
    ap.add_argument("--manifest", action="store_true", help="write the frozen bundle manifest computed from --out")
    ap.add_argument("--manifest-path", type=Path, default=DEFAULT_MANIFEST, help="manifest path (default: validation/reference_cases/topas/v13-bundles-manifest.json)")
    args = ap.parse_args(argv)
    if args.check:
        return 0 if check(args.out, args.manifest_path) else 1
    if args.manifest:
        print("wrote", write_manifest(args.out, args.manifest_path))
        return 0
    for nm, s in generate(args.out):
        print(s, nm)
    return 0


if __name__ == "__main__":
    sys.exit(main())
