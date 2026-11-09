"""Generate the V3-008 plan revision 2 TOPAS case bundles (219 bundles, 438 files).

This is the single source of the revision-2 bundles under validation/reference_cases/topas/*-r2-*
(plan: validation/plans/v3-008/r2-topas-target.md, appendix A). It must not be changed without a plan
amendment. Nothing is run; only case.json and input.txt are written.

Amendment R2-SIZE-1 (adopted 2026-11-09; plan section 8, v3-008-acceptance.md "Amendments"):
histories per (ion, configuration) (table N_PART), CutForElectron = 1 km in the LAT-FULL evidence
bundles (seven modules unchanged; IDD-FULL keeps 0.05 mm), three LAT-FULL ecut1km timing smokes
and a LAT-EM step-length smoke.

Resolution R2-SMOKE-LP (2026-11-09): ten execution-only ListProcesses variants (seeds
20600428-20600437) of the ten smokes whose signed re-runs lacked the process list; the originals
stay byte-identical.

Usage:
    python -I validation/scripts/ion_em/gen_r2_bundles.py [--out DIR]   # write bundles (default: committed dir)
    python -I validation/scripts/ion_em/gen_r2_bundles.py --check       # regenerate in a temp dir, compare byte
                                                                        # for byte, print the bundle-set digest
    python -I validation/scripts/ion_em/gen_r2_bundles.py --manifest

--check also compares the committed frozen bundle manifest
(validation/results/transport/v8r2-topas/frozen-bundles.json) with a fresh regeneration;
--manifest writes that manifest from the bundles under --out.
"""

import json
import argparse
import hashlib
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "validation" / "reference_cases" / "topas"
REL_PREFIX = "validation/reference_cases/topas"
DEFAULT_MANIFEST = (
    REPO_ROOT / "validation" / "results" / "transport" / "v8r2-topas" / "frozen-bundles.json"
)
MANIFEST_SCHEMA = "ionmc-v8r2-frozen-bundles-1"
BUNDLE_SET_DIGEST = "76d785667211bb39e96d517f46d636ff4b60be212b0d7fda9d968993ef47ed66"

TOPAS_BASE = 20619000
# case number, ion tag, MeV/u, A, z, particle, planning CSDA range (project Bethe table, mm), L, h_z, h_x, set
CASES = [
    (1, "he", 50, 4, 2, "alpha", 22.5, 50, 0.02, 0.01, "main"),
    (2, "he", 150, 4, 2, "alpha", 159.3, 250, 0.1, 0.02, "main"),
    (3, "he", 250, 4, 2, "alpha", 382.6, 500, 0.1, 0.02, "main"),
    (4, "carbon", 80, 12, 6, "GenericIon(6,12)", 17.5, 50, 0.02, 0.01, "main"),
    (5, "carbon", 290, 12, 6, "GenericIon(6,12)", 163.2, 250, 0.1, 0.02, "main"),
    (6, "carbon", 430, 12, 6, "GenericIon(6,12)", 309.3, 450, 0.1, 0.02, "main"),
    (7, "oxygen", 100, 16, 8, "GenericIon(8,16)", 19.6, 50, 0.02, 0.01, "main"),
    (8, "oxygen", 350, 16, 8, "GenericIon(8,16)", 166.6, 250, 0.1, 0.02, "main"),
    (9, "oxygen", 430, 16, 8, "GenericIon(8,16)", 231.9, 350, 0.1, 0.02, "main"),
    (10, "he", 100, 4, 2, "alpha", 78.0, 150, 0.05, 0.02, "held-out"),
    (11, "carbon", 150, 12, 6, "GenericIon(6,12)", 53.1, 100, 0.05, 0.02, "held-out"),
    (12, "oxygen", 200, 16, 8, "GenericIon(8,16)", 65.5, 100, 0.05, 0.02, "held-out"),
]
# configuration -> first part index of its block
CONFIGS = {
    "idd-em": 0,
    "lat-em": 4,
    "idd-full": 8,
    "lat-full": 12,
}
# histories per part, N_part, per (configuration, ion); every case of an ion uses the upper-end
# value of its ion (dated sizing note 2026-11-09, plan section 8). LAT-FULL is re-sized after the
# ecut1km timing smokes.
N_PART = {
    "idd-em": {"he": 25000, "carbon": 25000, "oxygen": 25000},
    "lat-em": {"he": 25000, "carbon": 25000, "oxygen": 25000},
    "idd-full": {"he": 5000, "carbon": 3333, "oxygen": 2529},
    "lat-full": {"he": 5000, "carbon": 5000, "oxygen": 5000},
}
K_PARTS = 4
SMOKE_SEEDS = {  # rehearsal family 2060xxxx (execution only)
    (3, "idd-em"): 20600411, (3, "lat-em"): 20600412, (3, "idd-full"): 20600413, (3, "lat-full"): 20600414,
    (6, "idd-em"): 20600415, (6, "lat-em"): 20600416, (6, "idd-full"): 20600417, (6, "lat-full"): 20600418,
    (9, "idd-em"): 20600419, (9, "lat-em"): 20600420, (9, "idd-full"): 20600421, (9, "lat-full"): 20600422,
}
VERBOSE_SMOKE = (5, "idd-em", 20600423)
# amendment R2-SIZE-1 (rehearsal family, execution only)
ECUT_SMOKE_SEEDS = {(3, "lat-full"): 20600424, (6, "lat-full"): 20600425, (9, "lat-full"): 20600426}
LAT_VERBOSE_SMOKE = (5, "lat-em", 20600427)
# resolution R2-SMOKE-LP (2026-11-09, rehearsal family, execution only): ListProcesses variants of
# the ten smokes whose signed re-runs failed run validity for lack of the process list; each is
# identical to its original except the process-list request and the seed
LP_SMOKE_SEEDS = {  # (case, configuration): seed; originals: SMOKE_SEEDS 20600411-13, 15-17, 19-21
    (3, "idd-em"): 20600428, (3, "lat-em"): 20600429, (3, "idd-full"): 20600430,
    (6, "idd-em"): 20600431, (6, "lat-em"): 20600432, (6, "idd-full"): 20600433,
    (9, "idd-em"): 20600434, (9, "lat-em"): 20600435, (9, "idd-full"): 20600436,
}
LP_VERBOSE_SMOKE = (5, "idd-em", 20600437)  # original: VERBOSE_SMOKE (20600423)

EM_MODULES = 'sv:Ph/Default/Modules = 1 "g4em-standard_opt4"'
FULL_MODULES = (
    'sv:Ph/Default/Modules = 7 "g4em-standard_opt4" "g4h-phy_QGSP_BIC_HP" "g4decay" '
    '"g4h-elastic_HP" "g4stopping" "g4ion-QMD" "g4radioactivedecay"'
)
# LAT-FULL under amendment R2-SIZE-1: the seven modules with the matched electron cut of the EM runs
ECUT_CONFIG = "seven-module TOPAS configuration, delta electrons deposited locally"
ECUT_COMMENT = (
    "# Physics: the operator's final reference configuration (directive 2026-11-07), seven "
    "modules, delta electrons deposited locally (plan amendment R2-SIZE-1): matched electron cut "
    "CutForElectron = 1 km so that delta electrons deposit locally (as ionmc and the EM "
    "configuration; P4 verbose smoke: e- threshold 9.55 GeV in G4_WATER); CutForAllParticles "
    "0.05 mm for all other particles."
)
ECUT_PHYSICS = (
    f"Geant4_Modular (OpenTOPAS 4.3 / Geant4 11.4.2), FULL, {ECUT_CONFIG} (operator's final "
    "reference configuration, directive 2026-11-07, with the matched electron cut "
    "CutForElectron = 1 km of plan amendment R2-SIZE-1): sv:Ph/Default/Modules = 7 "
)
ECUT_CUTS = (
    "CutForAllParticles = 0.05 mm (all particles other than e-); CutForElectron = 1 km (plan "
    "amendment R2-SIZE-1, matched electron cut as in the EM configuration: e- production threshold "
    "9.55 GeV in G4_WATER, delta electrons deposit locally); EMRangeMax = 10 GeV; "
    "Ts/MaxStepNumber = 1e7"
)
ECUT_TIMING_RATIONALE = (
    "V3-008 plan revision 2 timing smoke for amendment R2-SIZE-1 (execution only): LAT-FULL with "
    f"CutForElectron = 1 km ({ECUT_CONFIG}); wall time and printouts (process list, module "
    "acceptance, killed tracks) re-size the LAT-FULL parts; run only after the sample seal commit; "
    "outputs archived by the controlled runner; read for sizing/step diagnostics only; no evidence"
)
FULL_LIST = (
    "g4em-standard_opt4, g4h-phy_QGSP_BIC_HP, g4decay, g4h-elastic_HP, g4stopping, g4ion-QMD, "
    "g4radioactivedecay"
)


def fmt(x):
    s = f"{x:g}"
    return s


def name_of(c, cfg, tail):
    _, ion, e = c[0], c[1], c[2]
    return f"{ion}-water-{e}mevu-r2-{cfg}-{tail}"


def scorer_block(name, quantity, xb, yb, zb, out, outtype, report, gen=False, atomic=None):
    lines = [
        f's:Sc/{name}/Quantity = "{quantity}"',
        f's:Sc/{name}/Component = "Water"',
        f"i:Sc/{name}/XBins = {xb}",
        f"i:Sc/{name}/YBins = {yb}",
        f"i:Sc/{name}/ZBins = {zb}",
    ]
    if gen:
        lines.append(f's:Sc/{name}/OnlyIncludeParticlesOfGeneration = "Primary"')
    if atomic is not None:
        lines.append(f"i:Sc/{name}/OnlyIncludeParticlesOfAtomicNumber = {atomic}")
    lines += [
        f's:Sc/{name}/OutputType = "{outtype}"',
        f's:Sc/{name}/OutputFile = "{out}"',
        f's:Sc/{name}/IfOutputFileAlreadyExists = "Overwrite"',
        f"sv:Sc/{name}/Report = {len(report)} " + " ".join(f'"{r}"' for r in report),
    ]
    return lines


def build(
    c, cfg, seed, histories, kind, part=None, ecut_1km=None, list_processes=True, lp_variant_of=None
):
    """One bundle. ecut_1km: CutForElectron = 1 km in FULL (default: LAT-FULL, R2-SIZE-1).

    list_processes: request the per-particle process list (False only for the 13 smokes
    already run).
    lp_variant_of: seed of the original smoke of which this is the ListProcesses variant
    (R2-SMOKE-LP); only the rationale text records it.
    """
    num, ion, e, a, z, particle, r_plan, L, hz, hx, cset = c
    full = cfg.endswith("full")
    lat = cfg.startswith("lat")
    if ecut_1km is None:
        ecut_1km = full and lat
    full_ecut = full and ecut_1km
    total_e = a * e
    nz_idd = round(L / hz)
    nx_lat = round(120 / hx)
    nz_lat = L
    hdr_kind = {
        "evidence": f"V3-008 plan revision 2, row {'TL' if lat else 'TD'}-{'FULL' if full else 'EM'}"
        f"{' (held-out set)' if cset == 'held-out' else ''}, case {num}, part {part} (seed = 20619000 + 100*part + case).",
        "timing": "V3-008 plan revision 2 TIMING SMOKE (execution only): timing/step smoke; run only after the sample seal commit; outputs archived by the controlled runner; read for sizing/step diagnostics only; no evidence.",
        "verbose": "V3-008 plan revision 2 STEP-LENGTH SMOKE (execution only): TrackingVerbosity 1 on 3 histories, to record whether the scorer bins limit the primary step length; run only after the sample seal commit; outputs archived by the controlled runner; read for sizing/step diagnostics only; no evidence.",
    }[kind]
    inp = [
        f"# Native OpenTOPAS 4.3 input. {hdr_kind} No validation claim.",
        "i:Ts/NumberOfThreads = 1",
        f"i:Ts/Seed = {seed}",
        'b:Ts/PauseBeforeQuit = "False"',
        'b:Ts/ShowCPUTime = "True"',
        f"i:Ts/TrackingVerbosity = {1 if kind == 'verbose' else 0}",
        "# Step cap raised from the TOPAS default (1e6); any killed track (stderr 'stuck'/'killed') makes the run invalid (plan revision 2, run validity).",
        "i:Ts/MaxStepNumber = 10000000",
        "",
        f"# Geometry: 120 x 120 x {L} mm G4_WATER box, entrance face at z = 0, beam axis = +z through the face centre.",
        's:Ge/World/Material = "Vacuum"',
        "d:Ge/World/HLX = 20 cm",
        "d:Ge/World/HLY = 20 cm",
        "d:Ge/World/HLZ = 60 cm",
        's:Ge/Water/Parent = "World"',
        's:Ge/Water/Type = "TsBox"',
        's:Ge/Water/Material = "G4_WATER"',
        "d:Ge/Water/HLX = 60 mm",
        "d:Ge/Water/HLY = 60 mm",
        f"d:Ge/Water/HLZ = {fmt(L / 2)} mm",
        f"d:Ge/Water/TransZ = {fmt(L / 2)} mm",
        's:Ge/Source/Parent = "World"',
        's:Ge/Source/Type = "Group"',
        "d:Ge/Source/TransZ = -0.001 mm",
        "",
        f"# Source: monoenergetic, zero-size, zero-divergence pencil beam along +z. BeamEnergy = TOTAL kinetic energy ({a} x {e} MeV/u = {total_e} MeV; P4 convention).",
        's:So/Beam/Type = "Beam"',
        's:So/Beam/Component = "Source"',
        f's:So/Beam/BeamParticle = "{particle}"',
        f"d:So/Beam/BeamEnergy = {total_e} MeV",
        "u:So/Beam/BeamEnergySpread = 0",
        's:So/Beam/BeamPositionDistribution = "None"',
        's:So/Beam/BeamAngularDistribution = "None"',
        f"i:So/Beam/NumberOfHistoriesInRun = {histories}",
        "",
    ]
    if full_ecut:
        inp += [
            ECUT_COMMENT,
            's:Ph/Default/Type = "Geant4_Modular"',
            FULL_MODULES,
            "d:Ph/Default/CutForAllParticles = 0.05 mm",
            "d:Ph/Default/CutForElectron = 1 km",
        ]
    elif full:
        inp += [
            "# Physics: the operator's final reference configuration (directive 2026-11-07), seven modules, default production cut of the proton full-physics cases.",
            's:Ph/Default/Type = "Geant4_Modular"',
            FULL_MODULES,
            "d:Ph/Default/CutForAllParticles = 0.05 mm",
        ]
    else:
        inp += [
            "# Physics: EM only (intermediate configuration), matched electron cut: CutForElectron = 1 km so that delta electrons deposit locally (as ionmc; P4 verbose smoke: e- threshold 9.55 GeV in G4_WATER).",
            's:Ph/Default/Type = "Geant4_Modular"',
            EM_MODULES,
            "d:Ph/Default/CutForAllParticles = 0.05 mm",
            "d:Ph/Default/CutForElectron = 1 km",
        ]
    if list_processes:
        inp += [
            "# Prints the per-particle process list so that the executed module set, g4decay "
            "included, is verifiable from the log (Codex REVIEW-6b4731a7).",
            'b:Ph/ListProcesses = "True"',
        ]
    inp += [
        "# EM table upper edge above the highest total kinetic energy (default 600 MeV is below He 150-250 MeV/u totals).",
        "d:Ph/Default/EMRangeMax = 10 GeV",
        "",
    ]
    outputs = []
    scorers = []
    if not lat:
        inp.append(f"# Laterally integrated depth dose: {nz_idd} depth bins of {fmt(hz)} mm over the full 120 x 120 mm field.")
        inp += scorer_block("IDD", "DoseToMedium", 1, 1, nz_idd, "idd_dose", "csv", ["Sum", "Standard_Deviation"])
        outputs.append("idd_dose.csv")
        scorers.append(
            f"IDD: DoseToMedium, 1x1x{nz_idd} bins (full 120x120 mm lateral extent, {fmt(hz)} mm depth bins); Report Sum, Standard_Deviation"
        )
        if full:
            inp += [""] + scorer_block(
                "IDDPrimary", "DoseToMedium", 1, 1, nz_idd, "idd_primary_dose", "csv", ["Sum"], gen=True
            )
            inp += [""] + scorer_block(
                "PrimaryFluence", "Fluence", 1, 1, L, "primary_fluence", "csv", ["Sum"], gen=True, atomic=z
            )
            outputs += ["idd_primary_dose.csv", "primary_fluence.csv"]
            scorers += [
                f"IDDPrimary: as IDD with OnlyIncludeParticlesOfGeneration = Primary (report-only attribution diagnostic; excludes delta-electron and all secondary dose); Report Sum",
                f"PrimaryFluence: Fluence (1/mm2), 1x1x{L} bins (1 mm), generation Primary and atomic number {z} (primary-ion survival, report-only); Report Sum",
            ]
    else:
        inp.append(
            f"# Projected lateral dose profiles: x profile ({nx_lat} x 1 x {nz_lat}) and y profile (1 x {nx_lat} x {nz_lat}), {fmt(hx)} mm lateral bins over the full 120 mm field, 1 mm depth bins. Binary, Sum only."
        )
        inp += scorer_block("LatX", "DoseToMedium", nx_lat, 1, nz_lat, "lat_x", "Binary", ["Sum"])
        inp += [""] + scorer_block("LatY", "DoseToMedium", 1, nx_lat, nz_lat, "lat_y", "Binary", ["Sum"])
        outputs += ["lat_x.bin", "lat_x.binheader", "lat_y.bin", "lat_y.binheader"]
        scorers += [
            f"LatX: DoseToMedium, {nx_lat}x1x{nz_lat} bins ({fmt(hx)} mm in x over 120 mm, full y, 1 mm depth); OutputType Binary; Report Sum",
            f"LatY: DoseToMedium, 1x{nx_lat}x{nz_lat} bins (full x, {fmt(hx)} mm in y over 120 mm, 1 mm depth); OutputType Binary; Report Sum",
        ]
        if full:
            inp += [""] + scorer_block("LatXPrimary", "DoseToMedium", nx_lat, 1, nz_lat, "lat_x_primary", "Binary", ["Sum"], gen=True)
            inp += [""] + scorer_block("LatYPrimary", "DoseToMedium", 1, nx_lat, nz_lat, "lat_y_primary", "Binary", ["Sum"], gen=True)
            outputs += ["lat_x_primary.bin", "lat_x_primary.binheader", "lat_y_primary.bin", "lat_y_primary.binheader"]
            scorers += [
                "LatXPrimary: as LatX with OnlyIncludeParticlesOfGeneration = Primary (the TL-FULL gating core-width observable); Binary; Sum",
                "LatYPrimary: as LatY with OnlyIncludeParticlesOfGeneration = Primary (the TL-FULL gating core-width observable); Binary; Sum",
            ]
    inp.append("")
    physics = (
        [ECUT_PHYSICS + FULL_LIST]
        if full_ecut
        else [f"Geant4_Modular (OpenTOPAS 4.3 / Geant4 11.4.2), FULL (operator's final reference configuration, directive 2026-11-07): sv:Ph/Default/Modules = 7 {FULL_LIST}"]
        if full
        else ['Geant4_Modular (OpenTOPAS 4.3 / Geant4 11.4.2), EM ONLY (intermediate configuration): sv:Ph/Default/Modules = 1 "g4em-standard_opt4"; no hadronic, ion, elastic or decay modules']
    )
    cuts = (
        ECUT_CUTS
        if full_ecut
        else "CutForAllParticles = 0.05 mm (default production cut of the proton full-physics cases); EMRangeMax = 10 GeV; Ts/MaxStepNumber = 1e7"
        if full
        else "CutForAllParticles = 0.05 mm; CutForElectron = 1 km (accepted, P4 verbose smoke: e- production threshold 9.55 GeV in G4_WATER, no delta electrons); EMRangeMax = 10 GeV; Ts/MaxStepNumber = 1e7"
    )
    if kind == "evidence":
        rationale = (
            f"V3-008 plan revision 2 (validation/plans/v3-008/r2-topas-target.md): reference side of row "
            f"{'TL' if lat else 'TD'}-{'FULL' if full else 'EM'}{' held-out set' if cset == 'held-out' else ''}, "
            f"case {num} ({ion} {e} MeV/u), part {part} of 4 "
            f"(parts {CONFIGS[cfg]}-{CONFIGS[cfg] + 3}); "
            "agreement with the specified TOPAS configuration is the acceptance reference "
            "(operator directive 2026-11-07"
            f"{f'; configuration: {ECUT_CONFIG}, plan amendment R2-SIZE-1' if full_ecut else ''}); "
            "not run before the reviewed freeze of revision 2; outputs not materialized before the ionmc samples of the row are recorded"
        )
    elif kind == "timing" and full_ecut:
        rationale = ECUT_TIMING_RATIONALE
    elif kind == "timing":
        rationale = "V3-008 plan revision 2 timing smoke (execution only): wall time and printouts (process list, module acceptance, killed tracks) size the parts; run only after the sample seal commit; outputs archived by the controlled runner; read for sizing/step diagnostics only; no evidence"
    else:
        rationale = "V3-008 plan revision 2 step-length smoke (execution only): TrackingVerbosity 1 records whether scorer-bin boundaries limit the primary step length (relevant to the Urban step dependence); run only after the sample seal commit; outputs archived by the controlled runner; read for sizing/step diagnostics only; no evidence"
    if lp_variant_of is not None:
        rationale += (
            f"; ListProcesses variant (R2-SMOKE-LP, 2026-11-09): identical to seed {lp_variant_of} "
            "except the process-list printout and the seed; execution only, not evidence"
        )
    units = (
        "IDD: TOPAS CSV with '#' header, dose in Gy (run total; divide by histories); converted to laterally integrated energy per depth (MeV/mm per primary) with the slab mass 120 x 120 mm x bin depth at 1.0 g/cm3. Fluence 1/mm2 run total. Bins 0-based, z from the entrance face"
        if not lat
        else "Lat*: TOPAS binary (dose per bin, run total Gy; float64 and x-fastest order as established for the proton lateral cases, to be confirmed against the .binheader); bins 0-based, x,y centred on the beam axis, z from the entrance face in 1 mm bins"
    )
    case = {
        "input": "input.txt",
        "physics": physics,
        "source": f"{particle} ({ion} {e} MeV/u), BeamEnergy {total_e} MeV total kinetic energy (P4 convention), monoenergetic (BeamEnergySpread 0), zero-size, zero-divergence pencil beam from z=-0.001 mm along +z at the face centre",
        "geometry": f"G4_WATER box 120 x 120 x {L} mm (x,y centred on the beam axis; z from 0 to {L} mm) in a vacuum world 400 x 400 x 1200 mm; L = 1.3 x planning CSDA range ({r_plan} mm, project Bethe table) rounded up to 50 mm",
        "materials": "G4_WATER (Geant4 NIST material, I=78 eV, density 1.0 g/cm3)",
        "cuts": cuts,
        "histories": histories,
        "seeds": [seed],
        "scorers": scorers,
        "units": units,
        "rationale": rationale,
        "expected_outputs": outputs,
        "arguments": [],
    }
    return case, "\n".join(inp)


def write(out, name, case, inp):
    d = out / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "case.json").write_text(json.dumps(case, indent=2, ensure_ascii=False) + "\n")
    (d / "input.txt").write_text(inp)


def generate(out):
    """Write all 219 r2 bundles under `out`; return [(name, seed)]."""
    out = Path(out)
    made = []
    for c in CASES:
        for cfg, p0 in CONFIGS.items():
            n = N_PART[cfg][c[1]]
            for k in range(K_PARTS):
                part = p0 + k
                seed = TOPAS_BASE + 100 * part + c[0]
                case, inp = build(c, cfg, seed, n, "evidence", part)
                nm = name_of(c, cfg, f"part{part}")
                write(out, nm, case, inp)
                made.append((nm, seed))
    bynum = {c[0]: c for c in CASES}
    for (num, cfg), seed in SMOKE_SEEDS.items():
        n = 1000 if cfg.endswith("em") else 200
        # the 12 timing smokes ran 2026-11-09 and stay as recorded
        # (0.05 mm LAT-FULL, no process list)
        case, inp = build(bynum[num], cfg, seed, n, "timing", ecut_1km=False, list_processes=False)
        nm = name_of(bynum[num], cfg, "timing-smoke")
        write(out, nm, case, inp)
        made.append((nm, seed))
    num, cfg, seed = VERBOSE_SMOKE
    case, inp = build(bynum[num], cfg, seed, 3, "verbose", list_processes=False)  # already run
    nm = name_of(bynum[num], cfg, "steps-verbose-smoke")
    write(out, nm, case, inp)
    made.append((nm, seed))
    for (num, cfg), seed in ECUT_SMOKE_SEEDS.items():
        case, inp = build(bynum[num], cfg, seed, 200, "timing", ecut_1km=True)
        nm = name_of(bynum[num], cfg, "ecut1km-timing-smoke")
        write(out, nm, case, inp)
        made.append((nm, seed))
    num, cfg, seed = LAT_VERBOSE_SMOKE
    case, inp = build(bynum[num], cfg, seed, 3, "verbose")
    nm = name_of(bynum[num], cfg, "steps-verbose-smoke")
    write(out, nm, case, inp)
    made.append((nm, seed))
    # R2-SMOKE-LP: ListProcesses variants of the ten smokes above (re-runs lacked the process list)
    for (num, cfg), seed in LP_SMOKE_SEEDS.items():
        n = 1000 if cfg.endswith("em") else 200
        orig = SMOKE_SEEDS[(num, cfg)]
        case, inp = build(bynum[num], cfg, seed, n, "timing", ecut_1km=False, lp_variant_of=orig)
        nm = name_of(bynum[num], cfg, "lp-timing-smoke")
        write(out, nm, case, inp)
        made.append((nm, seed))
    num, cfg, seed = LP_VERBOSE_SMOKE
    assert (num, cfg) == VERBOSE_SMOKE[:2]
    case, inp = build(bynum[num], cfg, seed, 3, "verbose", lp_variant_of=VERBOSE_SMOKE[2])
    nm = name_of(bynum[num], cfg, "lp-steps-verbose-smoke")
    write(out, nm, case, inp)
    made.append((nm, seed))
    return made


def bundle_set_digest(out):
    """sha256 of the sha256sum-style listing (sorted by path, final newline) of all r2 bundle files."""
    out = Path(out)
    lines = []
    for d in out.iterdir():
        if d.is_dir() and "-r2-" in d.name:
            for f in ("case.json", "input.txt"):
                h = hashlib.sha256((d / f).read_bytes()).hexdigest()
                lines.append((f"{REL_PREFIX}/{d.name}/{f}", f"{h}  {REL_PREFIX}/{d.name}/{f}"))
    lines.sort()
    text = "".join(line + "\n" for _, line in lines)
    return hashlib.sha256(text.encode()).hexdigest(), len(lines)


def manifest_text(out):
    """Frozen bundle manifest (JSON text) of the r2 bundles under `out`.

    Every TOPAS run of the set is bound to this manifest (Codex review requirement): per bundle the
    sha256 of case.json and input.txt, the counts and the bundle-set digest. Deterministic: sorted
    keys, indent 1, trailing newline. Fails if the bundle-set digest is not BUNDLE_SET_DIGEST.
    """
    out = Path(out)
    bundles = {}
    for d in out.iterdir():
        if d.is_dir() and "-r2-" in d.name:
            bundles[d.name] = {
                f: hashlib.sha256((d / f).read_bytes()).hexdigest()
                for f in ("case.json", "input.txt")
            }
    digest, nfiles = bundle_set_digest(out)
    if digest != BUNDLE_SET_DIGEST:
        raise ValueError(f"bundle-set digest {digest} != recorded {BUNDLE_SET_DIGEST}")
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "generated_from": "gen_r2_bundles.py",
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
        extra = sorted(d.name for d in committed.iterdir() if d.is_dir() and "-r2-" in d.name) if committed.is_dir() else []
        extra = [e for e in extra if e not in {nm for nm, _ in made}]
        digest, nfiles = bundle_set_digest(tmp)
    print(f"generated {len(made)} bundles, {nfiles} files; bundle-set digest {digest}")
    ok = digest == BUNDLE_SET_DIGEST
    print("digest matches recorded value" if ok else f"DIGEST MISMATCH (recorded {BUNDLE_SET_DIGEST})")
    for b in bad:
        print("DIFFERS or MISSING:", b)
    for e in extra:
        print("EXTRA committed r2 bundle:", e)
    manifest_ok = (
        fresh_manifest is not None
        and manifest_path.is_file()
        and manifest_path.read_text() == fresh_manifest
        and json.loads(manifest_path.read_text())["bundle_set_digest"] == BUNDLE_SET_DIGEST
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
    ap.add_argument(
        "--manifest",
        action="store_true",
        help="write the frozen bundle manifest computed from --out (default: committed dir)",
    )
    ap.add_argument(
        "--manifest-path",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="manifest path (default: validation/results/transport/v8r2-topas/frozen-bundles.json)",
    )
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
