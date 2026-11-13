"""Generate the V3-008 V5-SS (branch B) TOPAS bundles: 16 evidence parts + 1 scorer smoke.

Single source of ``validation/reference_cases/topas/v5ss-*`` and of their frozen manifest
``validation/reference_cases/topas/v5ss-frozen-manifest.json`` (plan row V5-SS; dated "V5-SS
execution note (2026-11-13)" of ``validation/plans/v3-008-acceptance.md``). The revision-2
generator ``gen_r2_bundles.py`` and its manifest stay byte-identical (the r2 evaluator rejects a
run whose r2 manifest at ``code_sha`` differs from HEAD's; ``r2_io.frozen_manifest``). Nothing is
run; only ``case.json`` and ``input.txt`` are written.

Every line of ``input.txt`` is that of the committed 1 mm SS smoke
(``carbon-water-290mevu-ss-slab1mm-smoke``) except seed, histories, particle, energy and the
scorer: the only scorer is a PhaseSpace scorer on ``Water/ZPlusSurface`` (primaries, outgoing,
ASCII, no time/run/event/track ids). He uses ``BeamParticle = "alpha"`` as the committed r2 helium
bundles do; ``BeamEnergy`` is the TOTAL kinetic energy (A x MeV/u; P4 convention).

Usage::

    python -I validation/scripts/ion_em/gen_ss_bundles.py [--out DIR] [--manifest-path P]
    python -I validation/scripts/ion_em/gen_ss_bundles.py --check
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "validation" / "reference_cases" / "topas"
DEFAULT_MANIFEST = DEFAULT_OUT / "v5ss-frozen-manifest.json"
REL_PREFIX = "validation/reference_cases/topas"
MANIFEST_SCHEMA = "ionmc-v5ss-frozen-bundles-1"
PREFIX = "v5ss-"

TOPAS_BASE = 20619000  # seed = 20619000 + 100 * part + case (plan section C)
PARTS = (0, 1)
HISTORIES_PER_PART = 30_000  # 6e4 primaries per case in two parts (row V5-SS, 2026-11-01)
# case -> (ion tag, MeV/u, A, BeamParticle); He as the committed r2 helium bundles ("alpha")
CASES: dict[int, tuple[str, int, int, str]] = {
    22: ("he", 150, 4, "alpha"),
    25: ("carbon", 150, 12, "GenericIon(6,12)"),
    26: ("oxygen", 150, 16, "GenericIon(8,16)"),
    27: ("carbon", 250, 12, "GenericIon(6,12)"),
    28: ("he", 250, 4, "alpha"),
    29: ("carbon", 430, 12, "GenericIon(6,12)"),
    30: ("oxygen", 430, 16, "GenericIon(8,16)"),
    31: ("oxygen", 250, 16, "GenericIon(8,16)"),
}
SMOKE_NAME = "v5ss-helium-water-250mevu-scorer-smoke"
SMOKE_SEED = 20600438  # rehearsal family; 20600439 reserved for one re-smoke
SMOKE_HISTORIES = 200
SMOKE_CASE = 28  # He 250 configuration

SCORER_LINES = [
    's:Sc/SSExit/Quantity = "PhaseSpace"',
    's:Sc/SSExit/Surface = "Water/ZPlusSurface"',
    's:Sc/SSExit/OutputType = "ASCII"',
    's:Sc/SSExit/OutputFile = "ss_exit"',
    's:Sc/SSExit/IfOutputFileAlreadyExists = "Overwrite"',
    'b:Sc/SSExit/OutputToConsole = "False"',
    's:Sc/SSExit/OnlyIncludeParticlesGoing = "Out"',
    's:Sc/SSExit/OnlyIncludeParticlesOfGeneration = "Primary"',
    'b:Sc/SSExit/IncludeTOPASTime = "False"',
    'b:Sc/SSExit/IncludeTimeOfFlight = "False"',
    'b:Sc/SSExit/IncludeRunID = "False"',
    'b:Sc/SSExit/IncludeEventID = "False"',
    'b:Sc/SSExit/IncludeTrackID = "False"',
]
EXPECTED_OUTPUTS = ["ss_exit.phsp", "ss_exit.header"]


def topas_seed(case: int, part: int) -> int:
    return TOPAS_BASE + 100 * part + case


def bundle_name(case: int, part: int) -> str:
    ion, e, _a, _p = CASES[case]
    return f"{PREFIX}{ion}-water-{e}mevu-case{case:02d}-part{part}"


def build(case: int, seed: int, histories: int, kind: str, part: int | None = None):
    """(case.json dict, input.txt text) of one bundle; ``kind`` "evidence" or "smoke"."""
    ion, e, a, particle = CASES[case]
    total = a * e
    if kind == "evidence":
        head = (
            f"V3-008 row V5-SS (branch B), case {case} ({ion} {e} MeV/u), part {part} of parts 0-1 "
            f"(seed = 20619000 + 100*part + case)."
        )
        rationale = (
            f"V3-008 row V5-SS (branch B; plan validation/plans/v3-008-acceptance.md, V5-SS "
            f"execution note 2026-11-13): reference side, case {case} ({ion} {e} MeV/u), part "
            f"{part} of parts 0-1, seed 20619000 + 100*{part} + {case} = {seed}; exit angles of "
            "the primaries at the far face of the 1 mm slab under g4em-standard_SS; run only after "
            "the reviewed SS seal commit (SS-SEAL-1)"
        )
    else:
        head = (
            "V5-SS SCORER SMOKE (execution only): parameter acceptance and phase-space "
            "header/column layout; not evidence."
        )
        rationale = (
            "execution-only scorer smoke: parameter acceptance, header/column layout and "
            "precision; angles never used; not evidence"
        )
    inp = [
        f"# Native OpenTOPAS 4.3 input. {head} No validation claim.",
        "i:Ts/NumberOfThreads = 1",
        f"i:Ts/Seed = {seed}",
        'b:Ts/PauseBeforeQuit = "False"',
        "# Step cap raised far above the per-track single-scattering step count of a 1 mm slab "
        "(the 300 mm SS smoke lost every carbon track at the 1e6 default); any killed track "
        "invalidates the run.",
        "i:Ts/MaxStepNumber = 100000000",
        'b:Ts/ShowCPUTime = "True"',
        "# Tracking verbosity 0: no per-step output.",
        "i:Ts/TrackingVerbosity = 0",
        "",
        "# Geometry (V5-SS design, P10): 120 x 120 x 1 mm G4_WATER slab (z 0..1 mm), entrance face "
        "at z = 0, beam axis = +z through the face centre.",
        's:Ge/World/Material = "Vacuum"',
        "d:Ge/World/HLX = 20 cm",
        "d:Ge/World/HLY = 20 cm",
        "d:Ge/World/HLZ = 40 cm",
        's:Ge/Water/Parent = "World"',
        's:Ge/Water/Type = "TsBox"',
        's:Ge/Water/Material = "G4_WATER"',
        "d:Ge/Water/HLX = 60 mm",
        "d:Ge/Water/HLY = 60 mm",
        "d:Ge/Water/HLZ = 0.5 mm",
        "d:Ge/Water/TransZ = 0.5 mm",
        's:Ge/Source/Parent = "World"',
        's:Ge/Source/Type = "Group"',
        "d:Ge/Source/TransZ = -0.001 mm",
        "",
        "# Source: monoenergetic, zero-size, zero-divergence pencil beam along +z. BeamEnergy = "
        f"TOTAL kinetic energy ({a} x {e} MeV/u = {total} MeV; P4 convention).",
        's:So/Beam/Type = "Beam"',
        's:So/Beam/Component = "Source"',
        f's:So/Beam/BeamParticle = "{particle}"',
        f"d:So/Beam/BeamEnergy = {total} MeV",
        "u:So/Beam/BeamEnergySpread = 0",
        's:So/Beam/BeamPositionDistribution = "None"',
        's:So/Beam/BeamAngularDistribution = "None"',
        f"i:So/Beam/NumberOfHistoriesInRun = {histories}",
        "",
        "# Physics: EM only, single scattering (g4em-standard_SS; alias, IonCoulombScattering for "
        "GenericIon and alpha and no msc confirmed by the P10 SS smokes).",
        's:Ph/Default/Type = "Geant4_Modular"',
        'sv:Ph/Default/Modules = 1 "g4em-standard_SS"',
        "# Production cut: single cut for all particles, electron cut 1 km so that delta electrons "
        "deposit locally (e- energy cut >= T_max of the fastest case, about 1.16 MeV; checked in "
        "the couple table of the log).",
        "d:Ph/Default/CutForAllParticles = 0.05 mm",
        "d:Ph/Default/CutForElectron = 1 km",
        "# EM table upper edge above the highest total kinetic energy (O 430 MeV/u: 6880 MeV).",
        "d:Ph/Default/EMRangeMax = 10 GeV",
        "# Printouts: per-particle process list (no msc may be registered).",
        'b:Ph/ListProcesses = "True"',
        "",
        "# The only scorer: phase space of the primaries leaving the far face (z = 1 mm), ASCII; "
        "no time, run, event or track ids. Header and column layout: see ss_evaluate.py.",
        *SCORER_LINES,
        "",
    ]
    case_doc = {
        "input": "input.txt",
        "physics": [
            "Geant4_Modular (OpenTOPAS 4.3 / Geant4 11.4.2), EM ONLY, single scattering: "
            'sv:Ph/Default/Modules = 1 "g4em-standard_SS"; no hadronic, ion, elastic or decay '
            "modules"
        ],
        "source": f"{particle} ({ion} {e} MeV/u), BeamEnergy {total} MeV total kinetic energy (P4 "
        "convention), monoenergetic (BeamEnergySpread 0), zero-size, zero-divergence pencil beam "
        "from z=-0.001 mm along +z at the face centre",
        "geometry": "G4_WATER slab 120 x 120 x 1 mm (V5-SS design; x,y centred on the beam axis; "
        "z from 0 to 1 mm) in a vacuum world 400 x 400 x 800 mm",
        "materials": "G4_WATER (Geant4 NIST material, I=78 eV, density 1.0 g/cm3)",
        "cuts": "CutForAllParticles = 0.05 mm; CutForElectron = 1 km (e- energy cut >= 1.16 MeV "
        "required, checked in the couple table); EMRangeMax = 10 GeV; Ts/MaxStepNumber = 1e8",
        "histories": histories,
        "seeds": [seed],
        "scorers": ["PhaseSpace SSExit at Water/ZPlusSurface, primaries, outgoing, ASCII"],
        "units": "cm; direction cosines; MeV",
        "rationale": rationale,
        "expected_outputs": list(EXPECTED_OUTPUTS),
        "arguments": [],
    }
    return case_doc, "\n".join(inp)


def write(out: Path, name: str, case_doc: dict, inp: str) -> None:
    d = out / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "case.json").write_text(json.dumps(case_doc, indent=2, ensure_ascii=False) + "\n")
    (d / "input.txt").write_text(inp)


def generate(out: Path | str) -> list[tuple[str, int]]:
    """Write the 16 evidence bundles and the scorer smoke under ``out``; return [(name, seed)]."""
    out = Path(out)
    made = []
    for case in CASES:
        for part in PARTS:
            seed = topas_seed(case, part)
            case_doc, inp = build(case, seed, HISTORIES_PER_PART, "evidence", part)
            write(out, bundle_name(case, part), case_doc, inp)
            made.append((bundle_name(case, part), seed))
    case_doc, inp = build(SMOKE_CASE, SMOKE_SEED, SMOKE_HISTORIES, "smoke")
    write(out, SMOKE_NAME, case_doc, inp)
    made.append((SMOKE_NAME, SMOKE_SEED))
    return made


def bundle_names() -> list[str]:
    return sorted([bundle_name(c, p) for c in CASES for p in PARTS] + [SMOKE_NAME])


def manifest_text(out: Path | str) -> str:
    """Frozen SS manifest: sha256 of every file per bundle plus one overall digest (sha256 of the
    sorted sha256sum-style listing). Deterministic (sorted keys, indent 1, trailing newline)."""
    out = Path(out)
    bundles, lines = {}, []
    for name in bundle_names():
        entry = {}
        for f in ("case.json", "input.txt"):
            h = hashlib.sha256((out / name / f).read_bytes()).hexdigest()
            entry[f] = h
            lines.append(f"{h}  {REL_PREFIX}/{name}/{f}")
        bundles[name] = entry
    listing = "".join(line + "\n" for line in sorted(lines, key=lambda s: s.split("  ", 1)[1]))
    doc = {
        "schema": MANIFEST_SCHEMA,
        "generated_from": "gen_ss_bundles.py",
        "bundle_count": len(bundles),
        "file_count": len(lines),
        "bundle_set_digest": hashlib.sha256(listing.encode()).hexdigest(),
        "bundles": bundles,
    }
    return json.dumps(doc, sort_keys=True, indent=1) + "\n"


def check(committed: Path, manifest_path: Path) -> bool:
    with tempfile.TemporaryDirectory() as tmp:
        made = generate(tmp)
        bad = [
            f"{nm}/{f}"
            for nm, _ in made
            for f in ("case.json", "input.txt")
            if not (committed / nm / f).is_file()
            or (committed / nm / f).read_bytes() != (Path(tmp) / nm / f).read_bytes()
        ]
        fresh = manifest_text(tmp)
    ok_manifest = manifest_path.is_file() and manifest_path.read_text() == fresh
    for b in bad:
        print("DIFFERS or MISSING:", b)
    print("manifest matches" if ok_manifest else f"MANIFEST MISSING or DIFFERS: {manifest_path}")
    return not bad and ok_manifest


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--manifest-path", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    if a.check:
        return 0 if check(a.out, a.manifest_path) else 1
    for nm, s in generate(a.out):
        print(s, nm)
    a.manifest_path.write_text(manifest_text(a.out))
    print("wrote", a.manifest_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
