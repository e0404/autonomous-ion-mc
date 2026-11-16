"""Generate the V3-008 amendment I3-FU-1 step 2 TOPAS bundles: 40 + 40 repeat parts + 4 smokes.

Single source of ``validation/reference_cases/topas/i3fu2-*`` and of their frozen manifest
``validation/reference_cases/topas/i3fu2-frozen-manifest.json`` (plan
``validation/plans/v3-008-acceptance.md``, amendment I3-FU-1 step 2 and the section-C note of
2026-11-16). Nothing is run; only ``case.json`` and ``input.txt`` are written. Diagnostic bundles:
no validation claim.

Thin-slab mean energy loss of 12C, 16O (and a 4He control) in G4_WATER at E0 = 2.4-50 MeV/u
(TOPAS cases 13-21, :data:`CASES`; amendment I3-FU-1 revision 1). Evidence parts (:data:`PARTS`):
0 = C slab t, 1 = C slab 1.6t, 2 = O slab t, 3 = O slab 1.6t (all nine cases), 4 = He slab t (cases
17, 21), 5 = C slab t with energy-loss fluctuations disabled (cases 13, 21): 40 bundles. The
pre-declared repeat is frozen now: part + 8 (:data:`REPEAT_OFFSET`) for every evidence bundle,
same configuration (8 = C t, 9 = C 1.6t, 10 = O t, 11 = O 1.6t, 12 = He t, 13 = C fluct. off;
kind "repeat", run only on the repeat trigger of the evaluator's repeat rule): 40 bundles. Parts
6, 7 and 14-19 reserved (Codex REVIEW-22f384ed finding 5: repeats renumbered to leave room).
Seed = 20619000 + 100*part + case (unchanged formula). Four execution smokes
(:data:`SMOKES`): C 10 MeV/u, slab t, 2000 primaries at rehearsal seed 20600440 (standard lines)
and 20600441 (the fluctuation-off line), so that every parameter line of an evidence bundle is
exercised by a smoke; and two step-count smokes, C 2.4 and 50 MeV/u, slab 1.6t, 50 primaries,
``i:Ts/TrackingVerbosity = 1`` (seeds 20600443, 20600444), from whose step printout the
orchestrator counts the ionIoni steps per track inside the slab (expected: one) before the
evidence runs. 20600442 is reserved by the plan for one re-smoke (not generated).

Slab thickness from ionmc's own default Bethe model (``bethe_mass_stopping`` default options,
``water(78.0)``): t_mm = 10 * 0.03 * E0 * A / S(E0) (Delta E / E = 3 % at the entrance energy),
rounded to 3 significant digits; the second slab is exactly 1.6 x t (``thickness_factor``
0.048, Delta E / E ~ 4.8 %; the factors are exact so that the zero-thickness extrapolation
(f2 r1 - f1 r2) / (f2 - f1) of the evaluator removes a bias linear in the thickness). Both slab
losses exceed Geant4's GenericIon ``linLossLimit`` (2 %) and stay below the ~5.8 % energy loss
of one step allowed by the GenericIon ionisation step function (0.1, 1 um) (~0.1 R + 1.8 um),
so that each slab is meant to be crossed in one range-table-branch step (plan revision 1,
SR-3 finding 1; 4.8 % rather than 6 % keeps this even if TOPAS's stopping were 13 % above
ionmc's; checked by the step-count smokes). The model
exit energy of every slab must stay >= 2.0 MeV/u (:func:`exit_energy_per_u`, asserted for every
bundle). These are model-only numbers. The
physics lines are those of the TD-EM bundles (g4em-standard_opt4, CutForAllParticles 0.05 mm,
CutForElectron 1 km, EMRangeMax 10 GeV, ListProcesses), the fluctuation-off line for parts 5
and 13 and the fluctuation-off smoke, plus the printout ``i:Ph/Verbosity = 1``. Scorers: the
V5-SS ASCII phase space ``SSExit`` of the primaries leaving the far face and an
``EnergyDeposit`` Sum/Standard_Deviation scorer on the slab (cross check). ``BeamEnergy`` is
the TOTAL kinetic energy (A x MeV/u; P4 convention).

Usage::

    python validation/scripts/stopping/gen_i3fu2_bundles.py [--out DIR] [--manifest-path P]
    python validation/scripts/stopping/gen_i3fu2_bundles.py --check
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import sys
import tempfile
from decimal import Decimal
from pathlib import Path
from typing import Any

from ionmc.materials import water
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, Projectile
from ionmc.physics.stopping import bethe_mass_stopping

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "validation" / "reference_cases" / "topas"
DEFAULT_MANIFEST = DEFAULT_OUT / "i3fu2-frozen-manifest.json"
REL_PREFIX = "validation/reference_cases/topas"
MANIFEST_SCHEMA = "ionmc-i3fu2-frozen-bundles-1"
PREFIX = "i3fu2-"

TOPAS_BASE = 20619000  # seed = 20619000 + 100 * part + case (plan section C)
#: TOPAS case -> E0 [MeV/u] (section-C note of 2026-11-16); kept as decimal strings.
CASES: dict[int, str] = {
    13: "2.4",
    14: "3",
    15: "5",
    16: "7",
    17: "10",
    18: "15",
    19: "20",
    20: "30",
    21: "50",
}
#: ion tag -> (projectile, Z, A, TOPAS BeamParticle); He as the committed r2/V5-SS He bundles.
IONS: dict[str, tuple[Projectile, int, int, str]] = {
    "carbon": (CARBON12, 6, 12, "GenericIon(6,12)"),
    "oxygen": (OXYGEN16, 8, 16, "GenericIon(8,16)"),
    "he": (ALPHA, 2, 4, "alpha"),
}
ALL_CASES = tuple(CASES)
#: evidence part -> (ion tag, slab index (1: t, 2: 1.6 t), energy-loss fluctuations on, cases)
PARTS: dict[int, tuple[str, int, bool, tuple[int, ...]]] = {
    0: ("carbon", 1, True, ALL_CASES),
    1: ("carbon", 2, True, ALL_CASES),
    2: ("oxygen", 1, True, ALL_CASES),
    3: ("oxygen", 2, True, ALL_CASES),
    4: ("he", 1, True, (17, 21)),
    5: ("carbon", 1, False, (13, 21)),
}
REPEAT_OFFSET = 8  # repeat part = evidence part + 8 (frozen, same configuration)
REPEAT_PARTS = tuple(p + REPEAT_OFFSET for p in PARTS)  # 8-13
RESERVED_PARTS = (6, 7, *range(14, 20))  # not generated
HISTORIES = 50_000
DELTA_E_FRACTION = 0.03  # Delta E / E at E0 by the ionmc model (slab t)
#: slab index -> (multiple of t, thickness factor f = Delta E / E by design); exact decimals
SLABS: dict[int, tuple[Decimal, Decimal]] = {
    1: (Decimal("1"), Decimal("0.03")),
    2: (Decimal("1.6"), Decimal("0.048")),
}
I_EV_MODEL = 78.0
E_EXIT_MIN_PER_U = 2.0  # every slab's model exit energy must stay >= this [MeV/u]
N_RK4_GEN = 50  # RK4 steps of the generator's exit-energy guard (the evaluator uses >= 200)
#: smoke name -> (seed, evidence part whose configuration it runs, case, histories,
#: Ts/TrackingVerbosity); the step-count smokes print every step (TrackingVerbosity 1)
SMOKES: dict[str, tuple[int, int, int, int, int]] = {
    "i3fu2-smoke-carbon-water-10mevu": (20600440, 0, 17, 2000, 0),
    "i3fu2-smoke-carbon-water-10mevu-fluct-off": (20600441, 5, 17, 2000, 0),
    "i3fu2-smoke-carbon-water-2p4mevu-steps": (20600443, 1, 13, 50, 1),
    "i3fu2-smoke-carbon-water-50mevu-steps": (20600444, 1, 21, 50, 1),
}
SMOKE_NAME = "i3fu2-smoke-carbon-water-10mevu"
RESERVED_RESMOKE_SEED = 20600442  # plan: one re-smoke if a parameter name is wrong; not generated

#: Physics lines, byte-identical to the TD-EM bundles (e.g. carbon-water-150mevu-r2-idd-em-part0).
PHYSICS_LINES = [
    's:Ph/Default/Type = "Geant4_Modular"',
    'sv:Ph/Default/Modules = 1 "g4em-standard_opt4"',
    "d:Ph/Default/CutForAllParticles = 0.05 mm",
    "d:Ph/Default/CutForElectron = 1 km",
    'b:Ph/ListProcesses = "True"',
    "d:Ph/Default/EMRangeMax = 10 GeV",
]
PH_VERBOSITY_LINE = "i:Ph/Verbosity = 1"
FLUCTUATIONS_OFF_LINE = (
    'b:Ph/Default/ForRegion/DefaultRegionForTheWorld/EnableEnergyLossFluctuations = "False"'
)
#: Phase-space scorer, identical to the V5-SS bundles (gen_ss_bundles.SCORER_LINES).
SSEXIT_LINES = [
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
EDEP_LINES = [
    's:Sc/EDep/Quantity = "EnergyDeposit"',
    's:Sc/EDep/Component = "Water"',
    "i:Sc/EDep/XBins = 1",
    "i:Sc/EDep/YBins = 1",
    "i:Sc/EDep/ZBins = 1",
    's:Sc/EDep/OutputType = "csv"',
    's:Sc/EDep/OutputFile = "edep"',
    's:Sc/EDep/IfOutputFileAlreadyExists = "Overwrite"',
    'sv:Sc/EDep/Report = 2 "Sum" "Standard_Deviation"',
]
EXPECTED_OUTPUTS = ["ss_exit.header", "ss_exit.phsp", "edep.csv"]


def topas_seed(case: int, part: int) -> int:
    return TOPAS_BASE + 100 * part + case


def e_tag(case: int) -> str:
    """E0 for bundle names: "2.4" -> "2p4"."""
    return CASES[case].replace(".", "p")


def base_part(part: int) -> int:
    """Evidence part whose configuration ``part`` (evidence 0-5 or repeat 8-13) runs."""
    if part in PARTS:
        return part
    if part in REPEAT_PARTS:
        return part - REPEAT_OFFSET
    raise ValueError(f"part {part} is neither an evidence nor a repeat part")


def bundle_name(case: int, part: int) -> str:
    ion = PARTS[base_part(part)][0]
    return f"{PREFIX}{ion}-water-{e_tag(case)}mevu-case{case:02d}-part{part}"


def evidence_keys() -> list[tuple[int, int]]:
    """(case, part) of the 40 evidence bundles."""
    return [(c, p) for p, spec in PARTS.items() for c in spec[3]]


def repeat_keys() -> list[tuple[int, int]]:
    """(case, part + 8) of the 40 frozen repeat bundles."""
    return [(c, p + REPEAT_OFFSET) for c, p in evidence_keys()]


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


def slab_mm(ion: str, case: int, multiple: int = 1) -> float:
    """Slab thickness [mm] (water, 1.0 g/cm3): t with Delta E / E = 3 % at E0 by the ionmc
    default Bethe model (I = 78 eV), 3 significant digits; slab index ``multiple`` 2: exactly
    1.6 t (decimal arithmetic)."""
    if multiple not in SLABS:
        raise ValueError(f"slab index {multiple} is not 1 or 2")
    proj, _z, a, _p = IONS[ion]
    e0 = float(CASES[case])
    s = float(bethe_mass_stopping(e0, proj, water(I_EV_MODEL))[0])  # MeV cm2/g
    t = Decimal(f"{10.0 * DELTA_E_FRACTION * e0 * a / s:.3g}")
    return float(SLABS[multiple][0] * t)


@functools.cache
def exit_energy_per_u(ion: str, case: int, multiple: int) -> float:
    """Model exit energy [MeV/u] of the slab (RK4, :data:`N_RK4_GEN` steps, ionmc default
    Bethe); guard only (the evaluator's predictions use their own >= 200-step integration)."""
    proj, _z, a, _p = IONS[ion]
    mat = water(I_EV_MODEL)
    h = slab_mm(ion, case, multiple) / 10.0 / N_RK4_GEN  # cm, rho = 1 g/cm3

    def f(e: float) -> float:
        if e / a <= 1.0:  # far below the guard; keeps the model inside its domain
            raise ValueError(f"{ion} case {case} x{multiple}: slab leaves the model domain")
        return -float(bethe_mass_stopping(e / a, proj, mat)[0])

    e = a * float(CASES[case])
    for _ in range(N_RK4_GEN):
        k1 = f(e)
        k2 = f(e + 0.5 * h * k1)
        k3 = f(e + 0.5 * h * k2)
        k4 = f(e + h * k3)
        e += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return e / a


def fmt_mm(x: float) -> str:
    """Length in mm with at most 6 decimals; must be exact at that precision."""
    s = f"{x:.6f}".rstrip("0").rstrip(".")
    if abs(float(s) - x) > 1e-12:
        raise ValueError(f"{x} mm is not representable with 6 decimals")
    return s


def beam_energy(ion: str, case: int) -> str:
    """Total kinetic energy A x E0 [MeV] as a clean decimal string (12 x 2.4 -> "28.8")."""
    a = IONS[ion][2]
    return format((Decimal(a) * Decimal(CASES[case])).normalize(), "f")


def build(
    case: int, part: int, seed: int, histories: int, kind: str, tracking: int = 0
) -> tuple[dict[str, Any], str]:
    """(case.json dict, input.txt text) of one bundle; ``kind`` "evidence", "repeat" (``part``
    8-13, the configuration of ``part`` - 8) or "smoke" (``part`` the evidence part whose
    configuration it runs)."""
    ion, mult, fluct, _cases = PARTS[base_part(part)]
    if (kind == "repeat") != (part in REPEAT_PARTS):
        raise ValueError(f"kind {kind!r} does not fit part {part}")
    _proj, z, a, particle = IONS[ion]
    e0 = CASES[case]
    t = slab_mm(ion, case, mult)
    e_exit = exit_energy_per_u(ion, case, mult)
    if not e_exit >= E_EXIT_MIN_PER_U:
        raise ValueError(f"{ion} case {case} x{mult}: model exit {e_exit} MeV/u < 2.0 MeV/u")
    hlz = fmt_mm(round(t / 2.0, 6))
    t_txt = fmt_mm(t)
    total = beam_energy(ion, case)
    slab_word = "1.6t" if mult == 2 else "t"
    fl_word = "on" if fluct else "OFF"
    if kind == "repeat":
        orig = part - REPEAT_OFFSET
        head = (
            f"V3-008 amendment I3-FU-1 step 2 FROZEN REPEAT of part {orig} (run only on the "
            "evaluator's repeat trigger for this point), "
            f"case {case} ({ion} {e0} MeV/u), part {part} "
            f"(slab {slab_word} = {t_txt} mm, energy-loss fluctuations {fl_word}; seed = "
            f"20619000 + 100*part + case = {seed})."
        )
        rationale = (
            "V3-008 amendment I3-FU-1 step 2 (plan validation/plans/v3-008-acceptance.md): "
            f"pre-declared repeat (part {part} = part {orig} + 8, same configuration) of the "
            f"thin-slab mean energy loss of {ion} at E0 = {e0} MeV/u in G4_WATER, case {case} "
            f"(slab {slab_word} = {t_txt} mm, fluctuations {fl_word}), seed 20619000 + "
            f"100*{part} + {case} = {seed}; used only if the original run or its point is "
            "invalid; diagnostic, no validation claim"
        )
    elif kind == "evidence":
        head = (
            f"V3-008 amendment I3-FU-1 step 2 (thin-slab mean energy loss), case {case} "
            f"({ion} {e0} MeV/u), part {part} (slab {slab_word} = {t_txt} mm, energy-loss "
            f"fluctuations {fl_word}; seed = 20619000 + 100*part + case = {seed})."
        )
        rationale = (
            "V3-008 amendment I3-FU-1 step 2 (plan validation/plans/v3-008-acceptance.md): "
            "what stopping does the acceptance reference carry; thin-slab mean energy loss "
            f"of {ion} at E0 = {e0} MeV/u in G4_WATER, case {case}, part {part} (slab "
            f"{slab_word} = {t_txt} mm, fluctuations {fl_word}), seed 20619000 + 100*{part} + "
            f"{case} = {seed}; diagnostic, no validation claim"
        )
    elif kind == "smoke" and tracking:
        head = (
            "I3-FU-1 step 2 STEP-COUNT SMOKE (TrackingVerbosity 1: count the ionIoni steps per "
            f"track inside the slab, expected one; not evidence): case {case} configuration ({ion} "
            f"{e0} MeV/u, slab {slab_word} = {t_txt} mm), rehearsal seed {seed}."
        )
        rationale = (
            "V3-008 amendment I3-FU-1 step 2 step-count smoke: Ts/TrackingVerbosity 1 prints "
            "every step; the orchestrator counts the ionIoni steps per primary inside the slab "
            f"(expected one, range-table branch); {ion} {e0} MeV/u slab {slab_word} = {t_txt} mm, "
            f"fluctuations {fl_word}, rehearsal seed {seed} (20600442 reserved for one "
            "re-smoke); energies never used; not evidence; diagnostic, no validation claim"
        )
    elif kind == "smoke":
        head = (
            "I3-FU-1 step 2 EXECUTION SMOKE (parameter acceptance, scorer outputs; not "
            f"evidence): case {case} configuration ({ion} {e0} MeV/u, slab t = {t_txt} mm), "
            f"fluctuations {fl_word}, rehearsal seed {seed}."
        )
        rationale = (
            "V3-008 amendment I3-FU-1 step 2 execution smoke: parameter acceptance "
            "(Ph/Verbosity"
            + ("" if fluct else ", the energy-loss-fluctuation switch")
            + ", the phase-space and energy-deposit scorers), output layout; "
            f"{ion} {e0} MeV/u slab t = {t_txt} mm, fluctuations {fl_word}, rehearsal seed "
            f"{seed} (20600442 reserved for one re-smoke); energies never used; not evidence; "
            "diagnostic, no validation claim"
        )
    else:
        raise ValueError(f"unknown kind {kind!r}")
    physics = list(PHYSICS_LINES)
    if not fluct:
        physics.append(FLUCTUATIONS_OFF_LINE)
    inp = [
        f"# Native OpenTOPAS 4.3 input. {head} Diagnostic; no validation claim.",
        "i:Ts/NumberOfThreads = 1",
        f"i:Ts/Seed = {seed}",
        'b:Ts/PauseBeforeQuit = "False"',
        'b:Ts/ShowCPUTime = "True"',
        f"i:Ts/TrackingVerbosity = {tracking}",
        "# Step cap raised from the TOPAS default (1e6); any killed track invalidates the run.",
        "i:Ts/MaxStepNumber = 10000000",
        "# Run verbosity 2: production-cut couple table in the log; printout only.",
        "i:Ts/RunVerbosity = 2",
        "",
        f"# Geometry: 120 x 120 x {t_txt} mm G4_WATER slab (z 0..{t_txt} mm), entrance face at "
        "z = 0, beam axis = +z through the face centre; vacuum world.",
        's:Ge/World/Material = "Vacuum"',
        "d:Ge/World/HLX = 20 cm",
        "d:Ge/World/HLY = 20 cm",
        "d:Ge/World/HLZ = 40 cm",
        's:Ge/Water/Parent = "World"',
        's:Ge/Water/Type = "TsBox"',
        's:Ge/Water/Material = "G4_WATER"',
        "d:Ge/Water/HLX = 60 mm",
        "d:Ge/Water/HLY = 60 mm",
        f"d:Ge/Water/HLZ = {hlz} mm",
        f"d:Ge/Water/TransZ = {hlz} mm",
        's:Ge/Source/Parent = "World"',
        's:Ge/Source/Type = "Group"',
        "d:Ge/Source/TransZ = -0.001 mm",
        "",
        "# Source: monoenergetic, zero-size, zero-divergence pencil beam along +z. BeamEnergy = "
        f"TOTAL kinetic energy ({a} x {e0} MeV/u = {total} MeV; P4 convention).",
        's:So/Beam/Type = "Beam"',
        's:So/Beam/Component = "Source"',
        f's:So/Beam/BeamParticle = "{particle}"',
        f"d:So/Beam/BeamEnergy = {total} MeV",
        "u:So/Beam/BeamEnergySpread = 0",
        's:So/Beam/BeamPositionDistribution = "None"',
        's:So/Beam/BeamAngularDistribution = "None"',
        f"i:So/Beam/NumberOfHistoriesInRun = {histories}",
        "",
        "# Physics: EM only, the physics lines of the TD-EM bundles (g4em-standard_opt4; "
        "CutForElectron 1 km so that delta electrons deposit locally)"
        + ("" if fluct else "; energy-loss fluctuations DISABLED (control (iii))")
        + ".",
        *physics,
        "# Physics printout (verbosity 1); printout only.",
        PH_VERBOSITY_LINE,
        "",
        "# Scorer 1: phase space of the primaries leaving the far face, ASCII (as V5-SS); "
        "no time, run, event or track ids.",
        *SSEXIT_LINES,
        "",
        "# Scorer 2: energy deposited in the slab (cross check: Sum / histories = A E0 - "
        "<E_exit>).",
        *EDEP_LINES,
        "",
    ]
    case_doc: dict[str, Any] = {
        "input": "input.txt",
        "physics": [
            "Geant4_Modular (OpenTOPAS 4.3 / Geant4 11.4.2), EM ONLY: sv:Ph/Default/Modules = 1 "
            '"g4em-standard_opt4" (physics lines of the TD-EM bundles); no hadronic, ion, '
            "elastic or decay modules"
            + (
                ""
                if fluct
                else "; energy-loss fluctuations disabled "
                "(Ph/Default/ForRegion/DefaultRegionForTheWorld/EnableEnergyLossFluctuations)"
            )
        ],
        "source": f"{particle} ({ion} {e0} MeV/u), BeamEnergy {total} MeV total kinetic energy "
        "(P4 convention), monoenergetic (BeamEnergySpread 0), zero-size, zero-divergence pencil "
        "beam from z=-0.001 mm along +z at the face centre",
        "geometry": f"G4_WATER slab 120 x 120 x {t_txt} mm (x,y centred on the beam axis; z from "
        f"0 to {t_txt} mm) in a vacuum world 400 x 400 x 800 mm; slab {slab_word}, t from the "
        "ionmc default Bethe model (78 eV) for Delta E / E = 3 % at E0, 3 significant digits"
        + ("; 1.6t = exactly 1.6 x t" if mult == 2 else ""),
        "materials": "G4_WATER (Geant4 NIST material, I=78 eV, density 1.0 g/cm3)",
        "cuts": "CutForAllParticles = 0.05 mm; CutForElectron = 1 km (no delta electrons); "
        "EMRangeMax = 10 GeV; Ts/MaxStepNumber = 1e7",
        "histories": histories,
        "seeds": [seed],
        "scorers": [
            "PhaseSpace SSExit at Water/ZPlusSurface, primaries, outgoing, ASCII",
            "EDep: EnergyDeposit in Water, 1x1x1 bin; Report Sum, Standard_Deviation; csv",
        ],
        "units": "phase space: cm; direction cosines; kinetic energy MeV (total per ion); "
        "EDep: MeV, run total (divide by histories)",
        "rationale": rationale,
        "expected_outputs": list(EXPECTED_OUTPUTS),
        "arguments": [],
        "i3fu2": {
            "case": case,
            "part": part,
            "ion": ion,
            "z": z,
            "a": a,
            "e0_per_u_mev": float(e0),
            "slab_mm": t,
            "thickness_factor": float(SLABS[mult][1]),
            "kind": kind,
            "fluctuations": fluct,
        },
    }
    return case_doc, "\n".join(inp)


def write(out: Path, name: str, case_doc: dict[str, Any], inp: str) -> None:
    d = out / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "case.json").write_text(json.dumps(case_doc, indent=2, ensure_ascii=False) + "\n")
    (d / "input.txt").write_text(inp)


def generate(out: Path | str, histories: int = HISTORIES) -> list[tuple[str, int]]:
    """Write the 40 evidence bundles, the 40 repeat bundles and the four smokes under ``out``;
    return [(name, seed)]. ``histories`` (evidence and repeat bundles) differs from
    :data:`HISTORIES` only in tests."""
    out = Path(out)
    made = []
    for keys, kind in ((evidence_keys(), "evidence"), (repeat_keys(), "repeat")):
        for case, part in keys:
            seed = topas_seed(case, part)
            case_doc, inp = build(case, part, seed, histories, kind)
            write(out, bundle_name(case, part), case_doc, inp)
            made.append((bundle_name(case, part), seed))
    for name, (seed, part, case, hist, tracking) in SMOKES.items():
        case_doc, inp = build(case, part, seed, hist, "smoke", tracking)
        write(out, name, case_doc, inp)
        made.append((name, seed))
    return made


def bundle_names() -> list[str]:
    keys = evidence_keys() + repeat_keys()
    return sorted([bundle_name(c, p) for c, p in keys] + list(SMOKES))


def manifest_text(out: Path | str) -> str:
    """Frozen manifest: sha256 of every file per bundle plus one overall digest (sha256 of the
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
        "generated_from": "gen_i3fu2_bundles.py",
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
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
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
