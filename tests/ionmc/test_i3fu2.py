"""Amendment I3-FU-1 step 2 (revision 1): bundle generator, frozen manifest and evaluator.

Bundles ``validation/reference_cases/topas/i3fu2-*`` and their frozen manifest come from
``validation/scripts/stopping/gen_i3fu2_bundles.py``; the evaluator
``validation/scripts/stopping/i3fu2_evaluate.py`` is exercised in temporary git repositories
holding a synthetic bundle set (1000 histories), a committed review pin and Ed25519-signed
synthetic reference-runner records (test key from a fixed seed, committed as the trusted key of
the temporary repository, as ``test_v5ss_evaluate.py`` does), the sealed decision-producing
sources copied from this checkout, two synthetic passed review records (evaluator and markers
pins), a synthetic log-marker file and a plan holding the review-pin note. Synthetic TOPAS
outputs (phase space with the ``ss_evaluate`` header layout, a TOPAS CSV ``edep.csv``) carry
energy losses equal to the ionmc prediction times planted factors (by default including the
model's one-step algorithm bias 1 + g(f)); synthetic step-count smoke logs use a small
SteppingVerbose-like format parsed with the test marker regexes; the synthetic libdEdx
archive's ICRU73_NEW carbon/oxygen water tables are the ionmc model times a known factor. No
real libdEdx, ICRU or TOPAS data are used.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import re
import shutil
import subprocess
import sys
import tarfile
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from ionmc.materials import water
from ionmc.physics.projectiles import CARBON12, OXYGEN16
from ionmc.physics.stopping import bethe_mass_stopping

ROOT = Path(__file__).resolve().parents[2]
TOPAS = ROOT / "validation" / "reference_cases" / "topas"
MANIFEST = TOPAS / "i3fu2-frozen-manifest.json"
TD_EM = TOPAS / "carbon-water-150mevu-r2-idd-em-part0" / "input.txt"
SCRIPTS = ROOT / "validation" / "scripts"
GRID = np.geomspace(0.025, 1000.0, 53)
SYN_HISTORIES = 1000
ROW_NOISE = 0.004  # relative spread of the synthetic exit energies
TABLE_FACTOR = 1.05  # synthetic ICRU73_NEW = ionmc model x this factor (s ~ +5 %)
KEY_SEED = bytes(range(32))  # synthetic reference-runner key
VM_SEED = bytes(range(1, 33))  # synthetic validation-manager key (trusted in the temp repos)
ROGUE_SEED = bytes(range(2, 34))  # a validation-manager key nobody trusts
ATT_A_SIGNED = "2026-11-16T12:02:00+00:00"  # pin A attestation (before the smokes)
ATT_B_SIGNED = "2026-11-16T12:35:00+00:00"  # pin B attestation (before the evidence runs)
PIN_A_STARTED = "2026-11-16T11:50:00+00:00"  # evaluator review (pin A)
PIN_A_FINISHED = "2026-11-16T12:00:00+00:00"
SMOKE_STARTED = "2026-11-16T12:05:00+00:00"  # smokes: after pin A, before pin B starts
SMOKE_FINISHED = "2026-11-16T12:10:00+00:00"
PIN_B_STARTED = "2026-11-16T12:20:00+00:00"  # markers review (pin B)
PIN_B_FINISHED = "2026-11-16T12:30:00+00:00"
STARTED = "2026-11-16T13:00:00+00:00"  # evidence runs
FINISHED = "2026-11-16T14:00:00+00:00"
#: synthetic marker file of the temporary repositories (the real one is pinned from smoke logs)
TEST_MARKERS = {
    "fluct_on_marker": r"^Fluctuations of dE/dx are enabled\s+1\s*$",
    "fluct_off_marker": r"^Fluctuations of dE/dx are enabled\s+0\s*$",
    "step_line_re": r"^\s*(?P<step>\d+)\s+\S+ \S+\s+\S+ \S+\s+(?P<z>-?[0-9.eE+-]+) "
    r"(?P<zunit>fm|Ang|nm|um|mm|cm|m|km)\s",
    "track_header_re": r"^\* G4Track Information:\s+Particle = \S+,\s+Track ID = "
    r"(?P<track>\d+),\s+Parent ID = (?P<parent>\d+)",
}


def _load(name: str, path: Path) -> ModuleType:
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


gen = _load("gen_i3fu2_bundles", SCRIPTS / "stopping" / "gen_i3fu2_bundles.py")
ev = _load("i3fu2_evaluate", SCRIPTS / "stopping" / "i3fu2_evaluate.py")
gen_ss = _load("gen_ss_bundles", SCRIPTS / "ion_em" / "gen_ss_bundles.py")
seeds = _load("_v3_008_seed_guard", ROOT / "tests" / "ionmc" / "test_v3_008_seed_families.py")
IO = ev.rio

BUNDLES = sorted(d for d in TOPAS.iterdir() if d.is_dir() and d.name.startswith("i3fu2-"))
SMOKES = [d for d in BUNDLES if d.name in gen.SMOKES]
STEP_SMOKES = {  # step-count smokes: name -> case (configuration of part 1, slab 1.6 t)
    "i3fu2-smoke-carbon-water-2p4mevu-steps": 13,
    "i3fu2-smoke-carbon-water-50mevu-steps": 21,
}
RUNS = [d for d in BUNDLES if d.name not in gen.SMOKES]


def _native(d: Path, key: str) -> str:
    m = re.findall(rf"^{re.escape(key)} = (.+)$", (d / "input.txt").read_text(), re.M)
    assert len(m) == 1, (d.name, key)
    return str(m[0]).strip()


def _case(d: Path) -> dict[str, Any]:
    return dict(json.loads((d / "case.json").read_text()))


# =============================================================================================
# bundles
# =============================================================================================


def test_bundle_set() -> None:
    assert len(BUNDLES) == 84 and len(RUNS) == 80 and len(SMOKES) == 4
    assert [d.name for d in BUNDLES] == gen.bundle_names()
    evidence = sorted(
        [(c, p) for p in (0, 1, 2, 3) for c in range(13, 22)]
        + [(c, 4) for c in (17, 21)]
        + [(c, 5) for c in (13, 21)]
    )
    assert sorted(gen.evidence_keys()) == evidence and len(evidence) == 40
    assert gen.REPEAT_OFFSET == 8
    assert sorted(gen.repeat_keys()) == sorted((c, p + 8) for c, p in evidence)
    assert gen.REPEAT_PARTS == (8, 9, 10, 11, 12, 13)  # Codex REVIEW-22f384ed: part + 8
    assert gen.RESERVED_PARTS == (6, 7, 14, 15, 16, 17, 18, 19)
    assert not set(gen.PARTS) & set(gen.RESERVED_PARTS)
    assert not set(gen.REPEAT_PARTS) & set(gen.RESERVED_PARTS)
    energies = {c: float(e) for c, e in gen.CASES.items()}
    assert energies == {13: 2.4, 14: 3, 15: 5, 16: 7, 17: 10, 18: 15, 19: 20, 20: 30, 21: 50}
    kinds = {"evidence": 0, "repeat": 0}
    for d in RUNS:
        m = re.fullmatch(
            r"i3fu2-(carbon|oxygen|he)-water-([0-9p]+)mevu-case(\d\d)-part(\d+)", d.name
        )
        assert m, d.name
        info = _case(d)["i3fu2"]
        part = int(m.group(4))
        assert (info["ion"], info["case"], info["part"]) == (m.group(1), int(m.group(3)), part)
        kind = "repeat" if part >= 8 else "evidence"
        assert info["kind"] == kind
        kinds[kind] += 1
        base = part - 8 if part >= 8 else part
        assert info["thickness_factor"] == (0.048 if base in (1, 3) else 0.03)
        if part >= 8:  # a repeat is byte-identical to its original except seed and wording
            orig = _case(TOPAS / gen.bundle_name(info["case"], base))
            for key in ("physics", "source", "geometry", "histories", "scorers", "units"):
                assert _case(d)[key] == orig[key], (d.name, key)
            assert {k: v for k, v in info.items() if k not in ("part", "kind")} == {
                k: v for k, v in orig["i3fu2"].items() if k not in ("part", "kind")
            }
    assert kinds == {"evidence": 40, "repeat": 40}
    for d in SMOKES:
        info = _case(d)["i3fu2"]
        assert info["kind"] == "smoke" and info["ion"] == "carbon"
        if d.name in STEP_SMOKES:
            assert (info["case"], info["part"], info["thickness_factor"]) == (
                STEP_SMOKES[d.name],
                1,
                0.048,
            )
        else:
            assert info["case"] == 17 and info["thickness_factor"] == 0.03


def test_regeneration_is_byte_identical_and_manifest(tmp_path: Path) -> None:
    made = gen.generate(tmp_path)
    assert sorted(nm for nm, _ in made) == [d.name for d in BUNDLES]
    for nm, _ in made:
        for f in ("case.json", "input.txt"):
            assert (TOPAS / nm / f).read_bytes() == (tmp_path / nm / f).read_bytes(), nm
    assert gen.manifest_text(tmp_path) == MANIFEST.read_text()
    assert gen.check(TOPAS, MANIFEST)
    m = json.loads(MANIFEST.read_text())
    assert m["schema"] == "ionmc-i3fu2-frozen-bundles-1"
    assert m["bundle_count"] == 84 and m["file_count"] == 168
    assert ev.manifest_digest(m["bundles"]) == m["bundle_set_digest"]
    for d in BUNDLES:
        for f in ("case.json", "input.txt"):
            assert m["bundles"][d.name][f] == hashlib.sha256((d / f).read_bytes()).hexdigest()
    for p in [MANIFEST] + [d / f for d in BUNDLES for f in ("case.json", "input.txt")]:
        b = p.read_bytes()
        assert b.endswith(b"\n") and not b.endswith(b"\n\n"), p


def test_seeds_formula_and_reservation() -> None:
    topas = seeds.topas_seeds()
    used = []
    for d in RUNS:
        c = _case(d)
        case_no, part = c["i3fu2"]["case"], c["i3fu2"]["part"]
        s = seeds.TOPAS_8_BASE + 100 * part + case_no
        assert c["seeds"] == [s] == [gen.topas_seed(case_no, part)]
        assert int(_native(d, "i:Ts/Seed")) == s
        assert s in topas and s in seeds.V3_008_RESERVED
        assert case_no in seeds.TOPAS_CASES and part in seeds.TOPAS_PARTS
        used.append(s)
    assert len(set(used)) == 80 and min(used) == 20619013 and max(used) == 20620321
    smoke_seeds = {}
    for d in SMOKES:
        c = _case(d)
        assert c["histories"] == (50 if d.name in STEP_SMOKES else 2000)
        assert "not evidence" in c["rationale"]
        assert "20600442 reserved" in c["rationale"]
        smoke_seeds[d.name] = c["seeds"]
        assert int(_native(d, "i:Ts/Seed")) == c["seeds"][0]
    assert smoke_seeds == {
        "i3fu2-smoke-carbon-water-10mevu": [20600440],
        "i3fu2-smoke-carbon-water-10mevu-fluct-off": [20600441],
        "i3fu2-smoke-carbon-water-2p4mevu-steps": [20600443],
        "i3fu2-smoke-carbon-water-50mevu-steps": [20600444],
    }
    for s in (20600440, 20600441, 20600442, 20600443, 20600444):
        assert s in seeds.REHEARSAL_8 and s in seeds.V3_008_RESERVED
    assert gen.RESERVED_RESMOKE_SEED == 20600442
    all_text = "".join((d / "input.txt").read_text() for d in BUNDLES)
    assert "Ts/Seed = 20600442" not in all_text


def test_thickness_rule_and_exit_energy() -> None:
    proj = {"carbon": CARBON12, "oxygen": OXYGEN16}
    slabs: dict[tuple[str, int], dict[float, float]] = {}
    for d in BUNDLES:
        info = _case(d)["i3fu2"]
        ion, a, e0, t = info["ion"], info["a"], info["e0_per_u_mev"], info["slab_mm"]
        p = gen.IONS[ion][0]
        if ion in proj:
            assert p is proj[ion]
        s = float(bethe_mass_stopping(e0, p, water(78.0))[0])
        t_full = float(f"{10 * 0.03 * e0 * a / s:.3g}")
        factor = info["thickness_factor"]
        assert factor in (0.03, 0.048)
        assert t == pytest.approx(t_full * factor / 0.03, abs=1e-12), d.name  # 1.6t = 1.6 x t
        assert abs(10 * 0.03 * e0 * a / s / t_full - 1.0) <= 0.005, d.name  # 3-digit rounding
        de = ev.predicted_loss(f"model:{ion}", ev.model_stopping(p), a, e0, t)  # raises < 2 MeV/u
        assert (a * e0 - de) / a >= 2.0, d.name
        assert gen.exit_energy_per_u(ion, info["case"], 2 if factor == 0.048 else 1) > 2.0
        assert abs(de / (a * e0) - factor) < 0.003, d.name
        slabs.setdefault((ion, info["case"]), {})[factor] = t
        hlz = float(_native(d, "d:Ge/Water/HLZ").removesuffix(" mm"))
        assert _native(d, "d:Ge/Water/TransZ") == _native(d, "d:Ge/Water/HLZ")
        assert abs(2 * hlz - t) < 1e-12, d.name
        assert re.fullmatch(r"\d+\.\d{1,6} mm", _native(d, "d:Ge/Water/HLZ"))
    for ion in ("carbon", "oxygen"):  # every point has the t / 1.6 t pair
        for c in gen.CASES:
            pair = slabs[(ion, c)]
            assert pair[0.048] == pytest.approx(1.6 * pair[0.03], abs=1e-12)
    # the lowest point: 2.4 MeV/u, 4.8 % -> ~2.28 MeV/u, above the 2.0 MeV/u floor
    assert 2.25 < gen.exit_energy_per_u("carbon", 13, 2) < 2.31
    assert 2.25 < gen.exit_energy_per_u("oxygen", 13, 2) < 2.31


def test_physics_block_identical_to_td_em() -> None:
    def ph(text: str) -> list[str]:
        return [ln for ln in text.splitlines() if re.match(r"^\w+:Ph/", ln)]

    td = ph(TD_EM.read_text())
    assert td == gen.PHYSICS_LINES
    fluct_off = {gen.bundle_name(c, p) for c in (13, 21) for p in (5, 13)}
    fluct_off.add("i3fu2-smoke-carbon-water-10mevu-fluct-off")
    for d in BUNDLES:
        lines = ph((d / "input.txt").read_text())
        extra = [gen.PH_VERBOSITY_LINE]
        if d.name in fluct_off:
            extra.append(gen.FLUCTUATIONS_OFF_LINE)
            assert _case(d)["i3fu2"]["fluctuations"] is False
        else:
            assert _case(d)["i3fu2"]["fluctuations"] is True
        assert sorted(ln for ln in lines if ln not in td) == sorted(extra), d.name
        assert [ln for ln in lines if ln in td] == td, d.name
    # every parameter line of an evidence bundle is exercised by a smoke
    smoke_lines = {ln for d in SMOKES for ln in (d / "input.txt").read_text().splitlines()}
    for d in RUNS:
        for ln in (d / "input.txt").read_text().splitlines():
            if re.match(r"^\w+:(Ph|Sc|Ts)/", ln) and not ln.startswith("i:Ts/Seed"):
                assert ln in smoke_lines, (d.name, ln)
    # step-count smokes: the part-1 configuration of their case; only the seed, the histories
    # and i:Ts/TrackingVerbosity (smoke-only, 1) differ
    skip = ("i:Ts/Seed", "i:So/Beam/NumberOfHistoriesInRun", "i:Ts/TrackingVerbosity", "#")
    for name, case in STEP_SMOKES.items():
        smoke = (TOPAS / name / "input.txt").read_text().splitlines()
        ref = (TOPAS / gen.bundle_name(case, 1) / "input.txt").read_text().splitlines()
        assert [ln for ln in smoke if not ln.startswith(skip)] == [
            ln for ln in ref if not ln.startswith(skip)
        ]
        assert _native(TOPAS / name, "i:Ts/TrackingVerbosity") == "1"
        assert _native(TOPAS / name, "i:So/Beam/NumberOfHistoriesInRun") == "50"


def test_native_lines_scorers_and_case_keys() -> None:
    r2_keys = list(_case(TOPAS / "he-water-150mevu-r2-lat-em-part4"))
    assert gen.SSEXIT_LINES == gen_ss.SCORER_LINES
    for d in BUNDLES:
        c = _case(d)
        assert list(c) == r2_keys + ["i3fu2"]
        assert c["expected_outputs"] == ["ss_exit.header", "ss_exit.phsp", "edep.csv"]
        n_hist = (50 if d.name in STEP_SMOKES else 2000) if d in SMOKES else 50000
        assert c["histories"] == n_hist, d.name
        assert int(_native(d, "i:So/Beam/NumberOfHistoriesInRun")) == c["histories"]
        info = c["i3fu2"]
        _p, z, a, particle = gen.IONS[info["ion"]]
        assert (info["z"], info["a"]) == (z, a)
        assert _native(d, "s:So/Beam/BeamParticle") == f'"{particle}"'
        energy = _native(d, "d:So/Beam/BeamEnergy")
        assert abs(float(energy.removesuffix(" MeV")) - a * info["e0_per_u_mev"]) < 1e-9
        assert re.fullmatch(r"\d+(\.\d+)? MeV", energy) and not energy.startswith("0")
        for key, want in (
            ("i:Ts/MaxStepNumber", "10000000"),
            ("i:Ts/RunVerbosity", "2"),
            ("i:Ph/Verbosity", "1"),
            ("i:Ts/TrackingVerbosity", "1" if d.name in STEP_SMOKES else "0"),
            ("i:Ts/NumberOfThreads", "1"),
            ("b:Ts/PauseBeforeQuit", '"False"'),
            ("b:Ts/ShowCPUTime", '"True"'),
            ("d:Ge/Source/TransZ", "-0.001 mm"),
            ("u:So/Beam/BeamEnergySpread", "0"),
            ("d:Ge/World/HLZ", "40 cm"),
            ("d:Ge/Water/HLX", "60 mm"),
        ):
            assert _native(d, key) == want, (d.name, key)
        text = (d / "input.txt").read_text()
        assert [
            ln for ln in text.splitlines() if ln.startswith(("s:Sc", "b:Sc", "i:Sc", "sv:Sc"))
        ] == (gen.SSEXIT_LINES + gen.EDEP_LINES)
        head = text.splitlines()[0]
        assert "I3-FU-1 step 2" in head and "no validation claim" in head
        if d not in SMOKES:
            assert "20619000 + 100*part + case" in head
            assert "diagnostic, no validation claim" in c["rationale"]
            if c["i3fu2"]["kind"] == "repeat":
                assert "FROZEN REPEAT" in head and "only if the original" in c["rationale"]


# =============================================================================================
# evaluator: pure functions
# =============================================================================================

HEADER_TEMPLATE = """TOPAS ASCII Phase Space

Number of Original Histories: {n}
Number of Original Histories that Reached Phase Space: {reached}
Number of Scored Particles: {scored}

Columns of data are as follows:
 1: Position X [cm]
 2: Position Y [cm]
 3: Position Z [cm]
 4: Direction Cosine X
 5: Direction Cosine Y
 6: Energy [MeV]
 7: Weight
 8: Particle Type (in PDG Format)
 9: Flag to tell if Third Direction Cosine is Negative (1 means true)
10: Flag to tell if this is the First Scored Particle from this History (1 means true)

"""
EDEP_TEMPLATE = """# TOPAS Version: 4.3
# Parameter File: input.txt
# Results for scorer: EDep
# Scored in component: Water
# X in 1 bin  of 12 cm
# Y in 1 bin  of 12 cm
# Z in 1 bin  of {z_cm} cm
# EnergyDeposit ( MeV ) : Sum   Standard_Deviation
0, 0, 0, {total!r}, {sd:.6g}
"""


def test_synthetic_header_matches_fixture_layout() -> None:
    fx = (ROOT / "tests" / "ionmc" / "data" / "v5ss" / "ss_exit.header").read_text()
    syn = HEADER_TEMPLATE.format(n=3, reached=3, scored=3)
    assert ev.ss.parse_header(fx)["columns"] == ev.ss.parse_header(syn)["columns"]


def test_classify_rule() -> None:
    sig = 0.0025
    k = 1.645 * sig
    assert ev.classify(0.0, 0.0, 0.003, sig) == "U"  # X = 0.001 <= k
    assert ev.classify(0.05, 0.0, 0.05, sig) == "I"
    assert ev.classify(0.08, 0.03, 0.05, sig) == "H"  # TOPAS above the table
    assert ev.classify(0.0, -0.05, 0.05, sig) == "M"
    assert ev.classify(-0.03, -0.076, 0.05, sig) == "L"
    assert ev.classify(0.02, -0.03, 0.05, sig) == "N"
    assert ev.classify(0.06, 0.012, 0.05, sig) == "N"  # above the table, but not significantly
    assert ev.classify(0.0, 0.0, 3 * k, sig) == "U"  # X == k: U checked first
    assert ev.classify(0.0, 0.0, 0.0, sig) == "U"  # s = 0 is always U
    assert ev.band_of(2.4) == "[2,5)" and ev.band_of(5.0) == "[5,10)"
    assert ev.band_of(20.0) == "[20,50]" and ev.band_of(50.0) == "[20,50]"
    with pytest.raises(ValueError):
        ev.band_of(50.1)


def test_classify_relative_to_the_sign_of_s() -> None:
    """Reviewer pass 2 item (d): with the table BELOW ionmc (s_b < 0), H is TOPAS further below
    the table and L is TOPAS above ionmc; the mirrored inputs give the mirrored classes."""
    sig = 0.0025
    assert ev.classify(-0.05, 0.0, -0.05, sig) == "I"
    assert ev.classify(-0.08, -0.03, -0.05, sig) == "H"  # beyond the table, away from ionmc
    assert ev.classify(0.0, 0.05, -0.05, sig) == "M"
    assert ev.classify(0.03, 0.076, -0.05, sig) == "L"  # beyond ionmc, away from the table
    assert ev.classify(-0.02, 0.03, -0.05, sig) == "N"
    # the unsigned rule would have called these H / L: with s < 0 they are not
    assert ev.classify(0.08, 0.13, -0.05, sig) == "L"
    assert ev.classify(-0.03, 0.02, -0.05, sig) == "N"
    for r_m, r_i, s in ((0.08, 0.03, 0.05), (-0.03, -0.076, 0.05), (0.02, -0.03, 0.05)):
        assert ev.classify(-r_m, -r_i, -s, sig) == ev.classify(r_m, r_i, s, sig)


def _bands(classes: dict[str, tuple[str, str, str]]) -> dict[str, Any]:
    names = ("[2,5)", "[5,10)", "[10,20)")
    return {
        ion: {b: {"class": c} for b, c in zip(names, cl, strict=True)}
        for ion, cl in classes.items()
    }


@pytest.mark.parametrize(
    "c, o, want",
    [
        (("I", "I", "I"), ("I", "H", "U"), "ICRU"),
        (("H", "H", "H"), ("I", "I", "I"), "ICRU"),
        (("M", "L", "U"), ("M", "M", "M"), "MODEL"),
        (("L", "L", "L"), ("M", "M", "L"), "MODEL"),
        (("I", "I", "N"), ("I", "I", "I"), "MIXED"),
        (("I", "I", "I"), ("M", "M", "M"), "MIXED"),
        (("M", "M", "I"), ("M", "M", "M"), "MIXED"),
        (("U", "I", "I"), ("I", "I", "I"), "MIXED"),
        (("I", "N", "I"), ("I", "I", "I"), "MIXED"),
        (("none", "I", "I"), ("I", "I", "I"), "MIXED"),
    ],
)
def test_pattern_rule(c: tuple[str, str, str], o: tuple[str, str, str], want: str) -> None:
    assert ev.pattern(_bands({"carbon": c, "oxygen": o})) == want


def test_richardson_cancels_a_linear_bias() -> None:
    for r0, b in ((0.05, 0.1), (-0.03, -0.2), (0.0, 0.13)):
        assert abs(ev.richardson(r0 + b * 0.03, r0 + b * 0.048) - r0) < 1e-12
    assert ev.FACTORS == (0.03, 0.048)
    se = ev.richardson_se(1e-3, 2e-3)
    assert se == pytest.approx(np.hypot(0.048e-3, 0.03 * 2e-3) / 0.018, rel=1e-12)


def test_rk4_against_constant_stopping_and_domain_guard() -> None:
    assert ev.slab_energy_loss(120.0, 12, 0.1, lambda e: 500.0) == pytest.approx(5.0, rel=1e-12)
    with pytest.raises(ValueError, match="200"):
        ev.slab_energy_loss(120.0, 12, 0.1, lambda e: 500.0, n_steps=100)
    with pytest.raises(ValueError, match="domain"):
        ev.slab_energy_loss(28.8, 12, 1.0, lambda e: 500.0)  # 28.8 - 50 MeV < 2 MeV/u
    s = ev.table_stopping(GRID, 100.0 * GRID**-0.8)
    assert s(10.0) == pytest.approx(100.0 * 10.0**-0.8, rel=1e-12)


# --- stepping model (Codex REVIEW-22f384ed finding 4) ----------------------------------------


def test_stepping_model_residual_within_limit_at_every_point() -> None:
    """R_model <= 0.001 for all 18 C/O points and both He points (real model); the frozen-charge
    algorithm loses more than the exact path (q^2 falls with the energy along the slab, so
    q^2(E_pre) >= q^2(E)) and g grows with the thickness."""
    table = ev.stepping_model_check()
    assert sorted(table) == ["carbon", "he", "oxygen"]
    assert len(table["carbon"]) == len(table["oxygen"]) == 9 and sorted(table["he"]) == [17, 21]
    worst = {ion: max(abs(v["r_model"]) for v in pts.values()) for ion, pts in table.items()}
    assert all(w <= ev.MODEL_RESIDUAL_LIMIT for w in worst.values())
    # P13 record: R_model ~ -0.033 % (O 2.4 MeV/u), -0.024 % (C 2.4 MeV/u)
    assert 2e-4 < worst["carbon"] < 3e-4 and 3e-4 < worst["oxygen"] < 4e-4 and worst["he"] < 1e-6
    for ion in ("carbon", "oxygen"):
        assert all(v["r_model"] < 0.0 for v in table[ion].values())
    d_max = max(abs(v["d_pred"]) for v in table["oxygen"].values())
    assert d_max == pytest.approx(0.0033, abs=1e-4)  # P13: Delta_pred up to ~0.33 % (O 2.4)
    for ion in ("carbon", "oxygen"):
        for case, v in table[ion].items():
            assert v["d_pred"] == pytest.approx(v["g2"] - v["g1"])
            if case <= 20:  # <= 30 MeV/u: the algorithm loses more, more so in the thicker slab
                assert 0.0 < v["g1"] < v["g2"]
        assert table[ion][13]["g1"] > 3e-3  # ~ +0.36 % (C) / +0.52 % (O) at 2.4 MeV/u, f = 3 %
        assert abs(table[ion][21]["g2"]) < 1e-4  # 50 MeV/u: negligible either way


def test_algorithm_loss_decomposition() -> None:
    """P13 record: the C/O algorithm loss = table part with q^2(E_pre) S_base + the correction
    term at the mid-step energy; it exceeds the earlier whole-stopping frozen-charge emulation
    at 2.4 MeV/u, and He uses the frozen whole-stopping path."""
    for ion in ("carbon", "oxygen"):
        proj, z, a, _p = gen.IONS[ion]
        full = ev.model_stopping(proj)
        base = ev.base_stopping_per_q2(proj)
        unit = bethe_mass_stopping(10.0, proj, water(78.0), barkas=False, bloch=False,
                                   mott=False, effective_charge=False)[0]  # fmt: skip
        assert base(10.0) * z**2 == pytest.approx(float(unit), rel=1e-14)
        t = gen.slab_mm(ion, 13, 1)
        whole = ev.slab_energy_loss(a * 2.4, a, t, ev.frozen_charge_stopping(proj, full, 2.4))
        assert ev.algorithm_loss(proj, a, 2.4, t, ion) > whole
    proj, _z, a, _p = gen.IONS["he"]
    full = ev.model_stopping(proj)
    t = gen.slab_mm("he", 17, 1)
    assert ev.algorithm_loss(proj, a, 10.0, t, "he") == ev.slab_energy_loss(
        a * 10.0, a, t, ev.frozen_charge_stopping(proj, full, 10.0)
    )


def test_frozen_charge_stopping_is_exact_at_the_pre_step_energy() -> None:
    proj = gen.IONS["carbon"][0]
    s = ev.model_stopping(proj)
    alg = ev.frozen_charge_stopping(proj, s, 2.4)
    assert alg(2.4) == pytest.approx(s(2.4), rel=1e-14)
    assert alg(2.3) > s(2.3)  # q^2(2.4) > q^2(2.3): the frozen charge overestimates S
    assert ev.effective_charge_sq(proj, 2.4) < 36.0


def test_planted_algorithm_bias_is_consistent_and_extrapolates() -> None:
    """A TOPAS bias equal to the model's g(f) passes the consistency check and the
    extrapolation recovers the planted r0 up to exactly the modelled residual (1 + r0) R_model
    (|R_model| <= 3.5e-4, P13 record); a bias of the opposite sign fails the check at the low
    points."""
    fails = []
    for ion in ("carbon", "oxygen"):
        for case in gen.ALL_CASES:
            m = ev.stepping_bias(ion, case)
            for r0 in (0.05, -0.04):
                r1 = (1 + r0) * (1 + m["g1"]) - 1
                r2 = (1 + r0) * (1 + m["g2"]) - 1
                assert abs((r2 - r1) - m["d_pred"]) <= ev.CONSISTENCY_TOL
                got = ev.richardson(r1, r2)
                assert abs(got - (r0 + (1 + r0) * m["r_model"])) < 1e-12
                assert abs(got - r0) < 3.5e-4
                pt = ev.resolve_point(
                    {gen.bundle_name(case, 0): _fake(r=r1), gen.bundle_name(case, 1): _fake(r=r2)},
                    case,
                    (0, 1),
                    _ratios,
                    m["d_pred"],
                )
                assert pt.valid
            anti = (0.05 - m["g1"], 0.05 - m["g2"])
            if abs((anti[1] - anti[0]) - m["d_pred"]) > ev.CONSISTENCY_TOL:
                fails.append((ion, case))
    assert ("carbon", 13) in fails and ("oxygen", 13) in fails and ("oxygen", 14) in fails


def test_model_residual_above_limit_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(ev, "MODEL_RESIDUAL_LIMIT", 1e-5)
    with pytest.raises(ev.RecordError, match="design premise"):
        ev.stepping_model_check()


# --- repeat rule (Codex REVIEW-22f384ed finding 5) -------------------------------------------


def _fake(valid: bool = True, reason: str | None = None, r: float = 0.0) -> dict[str, Any]:
    return {"valid": valid, "reason": reason, "r": r}


def _ratios(run: dict[str, Any]) -> dict[str, Any]:
    return {"r_M": run["r"]}


def test_repeat_rule_unit() -> None:
    n = gen.bundle_name
    d = 0.001  # predicted difference

    def res(runs: dict[str, Any], case: int = 15, parts: tuple[int, ...] = (0, 1)) -> Any:
        return ev.resolve_point(runs, case, parts, _ratios, d if len(parts) == 2 else None)

    ok = {n(15, 0): _fake(r=0.01), n(15, 1): _fake(r=0.012)}
    pt = res(ok)
    assert pt.valid and pt.used == {} and pt.thickness_difference() == pytest.approx(0.002)
    assert pt.consistency_residual() == pytest.approx(0.001)
    for extra in (n(15, 8), n(15, 9)):  # a repeat without a recorded failure
        with pytest.raises(ev.RecordError, match="without a recorded failure"):
            res({**ok, extra: _fake()})
    # validity failure of one run: exactly that run is repeated
    bad_t = {n(15, 0): _fake(False, "deficit"), n(15, 1): _fake(r=0.012)}
    pt = res(bad_t)
    assert not pt.valid and pt.reason == "invalid run; no repeat"
    pt = res({**bad_t, n(15, 8): _fake(r=0.011)})
    assert pt.valid and pt.used == {n(15, 8): f"original {n(15, 0)} invalid: deficit"}
    with pytest.raises(ev.RecordError, match="is valid"):  # repeat of the valid 1.6t run
        res({**bad_t, n(15, 8): _fake(r=0.011), n(15, 9): _fake()})
    with pytest.raises(ev.RecordError, match="is valid"):
        res({**bad_t, n(15, 9): _fake()})
    pt = res({**bad_t, n(15, 8): _fake(False, "killed")})
    assert not pt.valid and "second failure" in str(pt.reason)
    pt = res({**bad_t, n(15, 8): _fake(r=0.0)})  # repeat valid but the new pair inconsistent
    assert not pt.valid and str(pt.reason).startswith(ev.CONSISTENCY_REASON)
    assert "second failure" in str(pt.reason)
    both_bad = {n(15, 0): _fake(False, "deficit"), n(15, 1): _fake(False, "killed")}
    with pytest.raises(ev.RecordError, match="partial repeat set"):
        res({**both_bad, n(15, 8): _fake(r=0.01)})
    assert res({**both_bad, n(15, 8): _fake(r=0.01), n(15, 9): _fake(r=0.011)}).valid
    # consistency failure (|0.01 - 0.001| > 0.003): BOTH slabs are repeated, the pair replaced
    incons = {n(15, 0): _fake(r=0.0), n(15, 1): _fake(r=0.01)}
    pt = res(incons)
    assert not pt.valid and pt.thickness_difference() == pytest.approx(0.01)
    assert pt.reason == f"{ev.CONSISTENCY_REASON}; no repeat"
    for rep in (n(15, 8), n(15, 9)):
        with pytest.raises(ev.RecordError, match="partial repeat set"):
            res({**incons, rep: _fake(r=0.0)})
    pt = res({**incons, n(15, 8): _fake(r=0.02), n(15, 9): _fake(r=0.021)})
    assert pt.valid and set(pt.used) == {n(15, 8), n(15, 9)}
    assert pt.thickness_difference() == pytest.approx(0.001)
    pt = res({**incons, n(15, 8): _fake(r=0.0), n(15, 9): _fake(r=0.02)})
    assert not pt.valid and "second failure" in str(pt.reason)
    pt = res({**incons, n(15, 8): _fake(False, "deficit"), n(15, 9): _fake(r=0.0)})
    assert not pt.valid and "second failure" in str(pt.reason)
    # single-run points (He, fluctuations off)
    with pytest.raises(ev.RecordError, match="missing"):
        res({n(17, 12): _fake()}, 17, (4,))
    assert res({}, 17, (4,)).reason == "missing run"
    assert res({n(17, 4): _fake(False, "x"), n(17, 12): _fake(r=0.0)}, 17, (4,)).valid
    with pytest.raises(ev.RecordError, match="without a recorded failure"):
        res({n(17, 4): _fake(), n(17, 12): _fake()}, 17, (4,))
    with pytest.raises(ValueError, match="predicted"):
        ev.resolve_point(ok, 15, (0, 1), _ratios)


# --- step-count parser (Codex REVIEW-22f384ed finding 3) -------------------------------------

STEP_RE = re.compile(TEST_MARKERS["step_line_re"], re.M)
TRACK_RE = re.compile(TEST_MARKERS["track_header_re"], re.M)


def _track(
    slab_mm: float,
    *,
    split: bool = False,
    primary: bool = True,
    tid: int = 1,
    stop: bool = False,
) -> str:
    """One synthetic SteppingVerbose block (G4BestUnit style values); ``stop``: the track ends
    inside the slab (half way) instead of crossing it."""
    parent = 0 if primary else 1
    lines = [
        f"* G4Track Information:   Particle = C12,   Track ID = {tid},   Parent ID = {parent}",
        "Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume Proc",
        "    0      0 fm      0 fm     -1 um   28.8 MeV   0 eV   0 fm   0 fm   World initStep",
        "    1      0 fm      0 fm      0 fm   28.8 MeV   0 eV   1 um   1 um   Water Transp",
    ]
    zs = [slab_mm / 2] if stop else ([slab_mm / 2, slab_mm] if split else [slab_mm])
    for i, z in enumerate(zs, start=2):
        lines.append(
            f"    {i}      0 fm      0 fm   {z * 1000:.6g} um   27 MeV   1 MeV   1 um   1 um   "
            "World Transportation"
        )
    if not stop:
        k = len(zs) + 2
        lines.append(
            f"    {k}      0 fm      0 fm     40 cm   27 MeV   0 eV   1 cm   1 cm   OutOfWorld T"
        )
    lines.append("    :----- List of 2ndaries - #SpawnInStep=  0(Rest= 0,Along= 0,Post= 0)")
    return "\n".join(lines) + "\n"


def _step_stdout(
    n: int, slab_mm: float, *, split: int = -1, stop: int = -1, primaries: int | None = None
) -> str:
    blocks = [_stdout(n)]
    for i in range(n if primaries is None else primaries):
        blocks.append(_track(slab_mm, split=i == split, stop=i == stop))
        blocks.append(_track(slab_mm / 3, primary=False, tid=2, split=True))  # ignored
    return "".join(blocks)


def _steps(counts: dict[int, dict[str, Any]]) -> dict[int, int]:
    return {k: int(v["steps"]) for k, v in counts.items()}


def test_count_slab_steps_contract() -> None:
    slab = 0.0657
    counts = ev.count_slab_steps(_step_stdout(5, slab), slab, STEP_RE, TRACK_RE)
    assert _steps(counts) == dict.fromkeys(range(5), 1)  # secondaries ignored, IDs repeat
    assert all(c["traverses"] for c in counts.values())
    counts = ev.count_slab_steps(_step_stdout(5, slab, split=3), slab, STEP_RE, TRACK_RE)
    assert _steps(counts) == {0: 1, 1: 1, 2: 1, 3: 2, 4: 1}  # a step count of 2
    counts = ev.count_slab_steps(_step_stdout(5, slab, stop=2), slab, STEP_RE, TRACK_RE)
    assert _steps(counts) == dict.fromkeys(range(5), 1)  # truncated track: one step ...
    assert [c["traverses"] for c in counts.values()] == [True, True, False, True, True]
    assert ev.count_slab_steps("no steps\n", slab, STEP_RE, TRACK_RE) == {}
    bad = _track(slab).replace("    2      0 fm", "    5      0 fm")
    with pytest.raises(ev.RecordError, match="consecutive"):
        ev.count_slab_steps(bad, slab, STEP_RE, TRACK_RE)
    with pytest.raises(ev.RecordError, match="step 0"):
        ev.count_slab_steps(_track(slab).replace("    0      0 fm", "    7      0 fm"), slab,
                            STEP_RE, TRACK_RE)  # fmt: skip
    plain = re.compile(r"^\s*(?P<step>\d+)\s+\S+ \S+\s+\S+ \S+\s+(?P<z>-?[0-9.eE+-]+) mm", re.M)
    mm_text = (
        "* G4Track Information:   Particle = C12,   Track ID = 1,   Parent ID = 0\n"
        "    0      0 fm      0 fm   -0.001 mm\n    1      0 fm      0 fm   0 mm\n"
        "    2      0 fm      0 fm   0.0657 mm\n    3      0 fm      0 fm   400 mm\n"
    )
    counts = ev.count_slab_steps(mm_text, slab, plain, TRACK_RE)  # no zunit group: mm
    assert counts == {0: {"steps": 1, "traverses": True}}


def test_count_slab_steps_printed_precision_tolerance() -> None:
    """Reviewer pass 3 item C3: the far face 0.06574 mm printed as "65.7 um" (3 digits, below the
    thickness) still ends the counted step at the far face, and the next pre-step z (within
    half a printed digit of the far face) is outside the slab: one step, traversing."""
    slab = 0.06574
    text = (
        "* G4Track Information:   Particle = C12,   Track ID = 1,   Parent ID = 0\n"
        "    0      0 fm      0 fm     -1 um   28.8 MeV\n"
        "    1      0 fm      0 fm      0 fm   28.8 MeV\n"
        "    2      0 fm      0 fm   65.7 um   27 MeV\n"
        "    3      0 fm      0 fm     40 cm   27 MeV\n"
    )
    assert ev.count_slab_steps(text, slab, STEP_RE, TRACK_RE) == {
        0: {"steps": 1, "traverses": True}
    }
    exact = text.replace("65.7 um", "65.74 um")  # printed exactly: the same outcome
    assert ev.count_slab_steps(exact, slab, STEP_RE, TRACK_RE)[0] == {
        "steps": 1,
        "traverses": True,
    }
    short = text.replace("65.7 um", "65.6 um")  # more than half a digit short: inside
    out = ev.count_slab_steps(short, slab, STEP_RE, TRACK_RE)[0]
    assert out == {"steps": 2, "traverses": True}  # the second step crosses the face
    stopped = text.replace("65.7 um", "65.6 um").replace("    3      0 fm      0 fm     40 cm"
                                                          "   27 MeV\n", "")  # fmt: skip
    assert ev.count_slab_steps(stopped, slab, STEP_RE, TRACK_RE)[0] == {
        "steps": 1,
        "traverses": False,
    }
    assert ev.printed_tolerance("65.7", 1e-3) == pytest.approx(0.05e-3)
    assert ev.printed_tolerance("40", 10.0) == pytest.approx(5.0)
    assert ev.printed_tolerance("1.25e-02", 1.0) == pytest.approx(0.005e-2)
    with pytest.raises(ev.RecordError, match="not a number"):
        ev.printed_tolerance("abc", 1.0)


def test_marker_file_parsing() -> None:
    ok = ev.parse_markers(json.dumps(TEST_MARKERS).encode(), "x")
    assert set(ok) == set(ev.MARKER_KEYS)
    for key in ev.MARKER_KEYS:
        with pytest.raises(ev.RecordError, match="non-empty"):
            ev.parse_markers(json.dumps({**TEST_MARKERS, key: ""}).encode(), "x")
    with pytest.raises(ev.RecordError, match="exactly the keys"):
        ev.parse_markers(json.dumps({**TEST_MARKERS, "extra": "a"}).encode(), "x")
    with pytest.raises(ev.RecordError, match="does not compile"):
        ev.parse_markers(json.dumps({**TEST_MARKERS, "fluct_on_marker": "("}).encode(), "x")
    with pytest.raises(ev.RecordError, match="named groups"):
        ev.parse_markers(json.dumps({**TEST_MARKERS, "step_line_re": r"^\s*\d+"}).encode(), "x")


def test_markers_discriminate_unit() -> None:
    assert ev.markers_discriminate([(True, True, False), (False, False, True)])
    assert not ev.markers_discriminate([(True, True, True), (False, False, True)])  # off hits on
    assert not ev.markers_discriminate([(True, True, False), (False, True, True)])  # on hits off
    assert not ev.markers_discriminate([(True, False, False), (False, False, True)])
    assert not ev.markers_discriminate([(True, True, False)])  # no fluctuation-off run


def _phase_space(
    case_doc: dict[str, Any], de: float, rng: np.random.Generator, *, noise: float, deficit: int = 0
) -> tuple[str, str, float]:
    """(header, rows, sum of the energy losses) of a synthetic run."""
    info = case_doc["i3fu2"]
    n = int(case_doc["histories"])
    a, e0, t = info["a"], info["e0_per_u_mev"], info["slab_mm"]
    de_rows = (
        de * (1.0 + noise * rng.standard_normal(n - deficit)) if noise else np.full(n - deficit, de)
    )
    rows = np.zeros((n - deficit, 10))
    rows[:, 2] = t / 10.0
    rows[:, 5] = a * e0 - de_rows
    rows[:, 6] = 1.0
    rows[:, 7] = ev.PDG[info["ion"]]
    rows[:, 9] = 1.0
    buf = io.StringIO()
    np.savetxt(buf, rows, fmt=["%.6g"] * 5 + ["%.17g", "%g", "%d", "%d", "%d"])
    header = HEADER_TEMPLATE.format(n=n, reached=n - deficit, scored=n - deficit)
    return header, buf.getvalue(), float(np.sum(de_rows))


def _edep(case_doc: dict[str, Any], total: float) -> str:
    return EDEP_TEMPLATE.format(z_cm=case_doc["i3fu2"]["slab_mm"] / 10.0, total=total, sd=total)


def test_measure_outputs_validity() -> None:
    doc = _case(TOPAS / gen.bundle_name(17, 0))
    doc["histories"] = 400
    de = 7.0
    h, rows, tot = _phase_space(doc, de, np.random.default_rng(1), noise=0.01)
    m = ev.measure_outputs(h, rows, _edep(doc, tot), doc)
    assert m["de"] == pytest.approx(de, rel=3e-3) and m["scored"] == 400
    h2, rows2, tot2 = _phase_space(doc, de, np.random.default_rng(1), noise=0.01, deficit=1)
    with pytest.raises(ev.RunInvalid, match="deficit"):
        ev.measure_outputs(h2, rows2, _edep(doc, tot2), doc)
    with pytest.raises(ev.RunInvalid, match="cross check"):
        ev.measure_outputs(h, rows, _edep(doc, tot * 1.002), doc)
    with pytest.raises(ev.RunInvalid, match="PDG"):
        ev.measure_outputs(
            h, rows, _edep(doc, tot), {**doc, "i3fu2": {**doc["i3fu2"], "ion": "oxygen"}}
        )
    thick = {**doc, "i3fu2": {**doc["i3fu2"], "slab_mm": 2 * doc["i3fu2"]["slab_mm"]}}
    with pytest.raises(ev.RunInvalid, match="far face"):
        ev.measure_outputs(h, rows, _edep(doc, tot), thick)
    assert (
        ev.validity_problems("Particle source Beam: Total number of histories: 5\n", b"", 5) == []
    )
    assert ev.validity_problems("Particle source Beam: Total number of histories: 5\n", b"x", 5)
    assert ev.validity_problems(
        "track killed\nParticle source Beam: Total number of histories: 5\n", b"", 5
    )
    assert ev.validity_problems("Particle source Beam: Total number of histories: 4\n", b"", 5)


# =============================================================================================
# evaluator: synthetic signed runs in temporary git repositories
# =============================================================================================


def _stdout(n: int, particle: str = "GenericIon", fluct: bool = True) -> str:
    return (
        "Registered Physics Processes:\n"
        "CoupledTransportation,     msc,     eIoni,     ionIoni,     CoulombScat\n\n"
        "Use ICRU90 data                                    1\n"
        f"Fluctuations of dE/dx are enabled                  {int(fluct)}\n\n"
        "msc:  for e-  SubType= 10\n"
        "      ===== EM models for the G4Region  DefaultRegionForTheWorld ======\n"
        "  GoudsmitSaunderson : Emin=    0 eV  Emax=  100 MeV Nbins=120 100 eV  - 100 MeV\n\n"
        f"ionIoni:  for {particle}  XStype:3  SubType=2\n"
        f"      StepFunction=(0.1, 0.001 mm), integ: 3, fluct: {int(fluct)}, linLossLim= 0.02\n\n"
        "Begin processing of event 0\n"
        f"Particle source Beam: Total number of histories: {n}\n"
    )


def _git(repo: Path, *args: str) -> str:
    env = {
        "GIT_AUTHOR_NAME": "t",
        "GIT_AUTHOR_EMAIL": "t@x",
        "GIT_COMMITTER_NAME": "t",
        "GIT_COMMITTER_EMAIL": "t@x",
        "PATH": "/usr/bin:/bin:/usr/local/bin",
        "HOME": str(repo),
    }
    return subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=True, env=env
    ).stdout.strip()


def _commit(repo: Path, msg: str) -> str:
    _git(repo, "add", "-A")
    _git(repo, "-c", "commit.gpgsign=false", "commit", "-q", "-m", msg)
    return _git(repo, "rev-parse", "HEAD")


def _archive(tmp: Path, factor: float) -> tuple[Path, str, dict[str, str]]:
    """Synthetic libdEdx archive: ICRU73_NEW water C/O = ionmc model x ``factor`` (E >= 1)."""
    energy_text = f"{GRID.size}\n" + "".join(f"{e:.10g}\n" for e in GRID)
    blocks = []
    for target in (104, 276):
        for z in range(1, 10):
            if target == 276 and z in (6, 8):
                proj = {6: CARBON12, 8: OXYGEN16}[z]
                s = bethe_mass_stopping(np.maximum(GRID, 1.0), proj, water(78.0)) * factor
            else:
                s = np.full(GRID.size, 1.0e6 + z)  # decoys
            blocks.append(f"#{target}:{z:02d}:{GRID.size}\n" + "".join(f"{v:.12E}\n" for v in s))
    files = {
        "libdedx-main/data/raw/icru73_newEng.dat": energy_text.encode(),
        "libdedx-main/data/raw/ICRU73_NEW.dat": "\n\n".join(blocks).encode(),
    }
    path = tmp / "archive.tar.gz"
    with tarfile.open(path, "w:gz") as tar:
        for name, data in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    members = {
        m: hashlib.sha256(files[f"libdedx-main/{m}"]).hexdigest()
        for m in ("data/raw/ICRU73_NEW.dat", "data/raw/icru73_newEng.dat")
    }
    return path, hashlib.sha256(path.read_bytes()).hexdigest(), members


def _review(review_id: str, head: str, started: str, finished: str, **over: Any) -> dict[str, Any]:
    """A synthetic passed codex-worker review record (schema of ``ss-seal-review.json``)."""
    doc = {
        "base_sha": "0" * 40,
        "experiment_id": "experiment-v3",
        "finished_at": finished,
        "head_sha": head,
        "report": {"findings": [], "head_sha": head, "summary": "synthetic", "verdict": "pass"},
        "report_sha256": hashlib.sha256(review_id.encode()).hexdigest(),
        "review_id": review_id,
        "started_at": started,
        "status": "passed",
        "task_id": "V3-008",
        "tool": "codex exec",
    }
    doc.update(over)
    return doc


def _attestation(
    review: dict[str, Any],
    signed_at: str,
    *,
    seed: bytes = VM_SEED,
    details: dict[str, Any] | None = None,
    **over: Any,
) -> dict[str, Any]:
    """A validation-manager local-validation record carrying ``review`` under
    ``details.review_record``, signed at ``signed_at`` (the layout of ``VAL-*.json``)."""
    at = IO.attest_module()
    rec: dict[str, Any] = {
        "schema_version": 1,
        "validation_id": f"VAL-syn-{review['review_id']}",
        "task_id": "V3-008",
        "recorded_at_utc": signed_at,
        "suite": "I3-FU-1 step-2 review pin",
        "status": "passed",
        "summary": "synthetic pin record",
        "details": {"review_record": review} if details is None else details,
        "git_branch": "task/v3-008",
        "git_sha": review["head_sha"],
        "worktree": "/synthetic",
    }
    rec.update(over)
    key = at.keypair_from_seed(seed)
    message = at.signed_message(rec, "validation-manager", signed_at, key.key_id)
    rec["attestation"] = {
        "schema": at.SCHEMA,
        "alg": at.ALG,
        "signer": "validation-manager",
        "key_id": key.key_id,
        "public_key": key.public_key.hex(),
        "signed_at_utc": signed_at,
        "payload_sha256": hashlib.sha256(at.canonical_payload(rec)).hexdigest(),
        "signature": at.sign(message, key.seed).hex(),
    }
    return rec


def _note(eval_id: str, eval_sha: str, markers_id: str, markers_sha: str) -> str:
    return (
        "**Dated note (2026-11-16, orchestrator; I3-FU-1 step-2 review pins).** evaluator review "
        f"{eval_id} at {eval_sha}; markers review {markers_id} at {markers_sha}.\n"
    )


def _write_json(path: Path, doc: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=1) + "\n")


def _sealed_copy_list() -> list[str]:
    """Files copied from this checkout into the temporary repositories: the executable closure,
    every ionmc source file (so that a later lazy import is present) and the sealed files."""
    rels = set(ev.executable_closure()) | set(ev.SEALED_FILES)
    rels |= {
        p.relative_to(ROOT).as_posix()
        for p in (ROOT / "src" / "ionmc").rglob("*.py")
        if "__pycache__" not in p.parts
    }
    return sorted(rels)


@pytest.fixture(scope="module")
def base(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """Temporary repository: the sealed closure + synthetic bundles (1000 histories) + manifest +
    trust (commit "reviewed" = the evaluator review's head_sha); pin A (commit "pin_a": the
    smokes run here); the marker file (commit "markers" = the markers review's head_sha); pin B
    + the plan note (commit "code_sha")."""
    root = tmp_path_factory.mktemp("i3fu2")
    repo = root / "repo"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "main")
    for rel in _sealed_copy_list():
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, repo / rel)
    bundles = repo / gen.REL_PREFIX
    gen.generate(bundles, histories=SYN_HISTORIES)
    (repo / ev.MANIFEST).write_text(gen.manifest_text(bundles))
    kp = IO.attest_module().keypair_from_seed(KEY_SEED)
    vm = IO.attest_module().keypair_from_seed(VM_SEED)
    _write_json(
        repo / IO.DEFAULT_TRUST,
        {
            "reference-runner": [{"key_id": kp.key_id, "public_key": kp.public_key.hex()}],
            "validation-manager": [{"key_id": vm.key_id, "public_key": vm.public_key.hex()}],
        },
    )
    seal = repo / ev.ENV_SEAL
    seal.parent.mkdir(parents=True, exist_ok=True)
    seal.write_text(ev.env_seal.seal_text(ev.env_seal.current_environment()))
    (repo / ev.PLAN_PATH).parent.mkdir(parents=True, exist_ok=True)
    (repo / ev.PLAN_PATH).write_text("# synthetic acceptance plan\n")
    reviewed = _commit(repo, "sealed evaluator and bundles")
    pin_a = _review("REVIEW-syn-a", reviewed, PIN_A_STARTED, PIN_A_FINISHED)
    _write_json(repo / ev.REVIEW_PIN, pin_a)
    _write_json(repo / ev.REVIEW_ATTESTATION, _attestation(pin_a, ATT_A_SIGNED))
    pin_a_sha = _commit(repo, "pin A and its signed validation record")
    _write_json(repo / ev.MARKERS, TEST_MARKERS)
    markers_sha = _commit(repo, "markers from the smoke logs")
    pin_b = _review("REVIEW-syn-b", markers_sha, PIN_B_STARTED, PIN_B_FINISHED)
    _write_json(repo / ev.MARKERS_PIN, pin_b)
    _write_json(repo / ev.MARKERS_ATTESTATION, _attestation(pin_b, ATT_B_SIGNED))
    (repo / ev.PLAN_PATH).write_text(
        "# synthetic acceptance plan\n\n"
        + _note("REVIEW-syn-a", reviewed, "REVIEW-syn-b", markers_sha)
    )
    code_sha = _commit(repo, "pin B and plan note")
    archive, sha, members = _archive(root, TABLE_FACTOR)
    cases, de_m = {}, {}
    for name in [gen.bundle_name(c, p) for c, p in gen.evidence_keys() + gen.repeat_keys()] + list(
        ev.ALL_SMOKES
    ):
        doc = json.loads((bundles / name / "case.json").read_text())
        info = doc["i3fu2"]
        stop = ev.model_stopping(gen.IONS[info["ion"]][0])
        cases[name] = doc
        de_m[name] = ev.predicted_loss(
            f"model:{info['ion']}", stop, info["a"], info["e0_per_u_mev"], info["slab_mm"]
        )
    return {
        "root": root,
        "repo": repo,
        "pin_a": pin_a,
        "pin_b": pin_b,
        "reviewed": reviewed,
        "pin_a_sha": pin_a_sha,
        "markers_sha": markers_sha,
        "code_sha": code_sha,
        "archive": archive,
        "sha": sha,
        "members": members,
        "cases": cases,
        "de_m": de_m,
    }


def _g(info: dict[str, Any]) -> float:
    """The model's one-step algorithm bias g(f) of this run's slab (planted as TOPAS's)."""
    m = ev.stepping_bias(info["ion"], info["case"])
    return float(m["g2"] if info["thickness_factor"] == 0.048 else m["g1"])


def _write_run(
    base: dict[str, Any],
    repo: Path,
    label: str,
    bundle: str,
    r: float,
    rng: np.random.Generator,
    *,
    noise: float = ROW_NOISE,
    deficit: int = 0,
    stdout: str | None = None,
    stderr: str = "",
    started: str = STARTED,
    finished: str = FINISHED,
    code_sha: str | None = None,
    succeeded: bool = True,
    exit_code: int = 0,
    tamper_record: bool = False,
    plant_g: bool = True,
) -> str:
    """Materialized outputs under ``<root>/out-<label>/<run_id>/`` and the signed record +
    listing under ``<repo>/records/<label>/<run_id>/`` (not committed here). The planted loss is
    the model loss x (1 + r) x (1 + g(f)) (``plant_g``: TOPAS's one-step algorithm)."""
    doc = base["cases"][bundle]
    info = doc["i3fu2"]
    code_sha = code_sha or base["code_sha"]
    bdir = repo / gen.REL_PREFIX / bundle
    files = {}
    for name in ("case.json", "input.txt"):
        blob = (bdir / name).read_bytes()
        files[name] = {"sha256": hashlib.sha256(blob).hexdigest(), "bytes": len(blob)}
    request = {
        "schema_version": 2,
        "code_sha": code_sha,
        "input_files": files,
        "case": doc,
        "dirty": False,
        "engine": {"engine": "topas"},
    }
    cfg = IO.canonical_digest(request)
    run_id = f"REF-{cfg[:20]}-{hashlib.sha256(label.encode()).hexdigest()[:8]}"
    odir = base["root"] / f"out-{label.removesuffix('-smokes')}" / run_id
    (odir / "work").mkdir(parents=True)
    (odir / "request.json").write_text(json.dumps(request))
    particle = ev.BEAM_PARTICLE[info["ion"]]
    n = int(doc["histories"])
    (odir / "stdout.txt").write_text(
        stdout if stdout is not None else _stdout(n, particle, bool(info["fluctuations"]))
    )
    (odir / "stderr.txt").write_text(stderr)
    for name in ("case.json", "input.txt"):
        shutil.copy(bdir / name, odir / "work" / name)
    g = _g(info) if plant_g and info["kind"] != "smoke" else 0.0
    h, rows, tot = _phase_space(
        doc, base["de_m"][bundle] * (1.0 + r) * (1.0 + g), rng, noise=noise, deficit=deficit
    )
    (odir / "work" / "ss_exit.header").write_text(h)
    (odir / "work" / "ss_exit.phsp").write_text(rows)
    (odir / "work" / "edep.csv").write_text(_edep(doc, tot))
    listing = [
        {
            "path": p.relative_to(odir).as_posix(),
            "bytes": p.stat().st_size,
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        }
        for p in sorted(odir.rglob("*"))
        if p.is_file()
    ]
    sha = {e["path"]: e["sha256"] for e in listing}
    work = {
        e["path"][5:]: {"sha256": e["sha256"], "bytes": e["bytes"]}
        for e in listing
        if e["path"].startswith("work/")
    }
    at = IO.attest_module()
    result = at.sign_record(
        {
            "schema_version": 3,
            "run_id": run_id,
            "configuration_id": cfg,
            "code_sha": code_sha,
            "task_id": "V3-008",
            "engine": "topas",
            "request_sha256": sha["request.json"],
            "started_at": started,
            "finished_at": finished,
            "seconds": 60.0,
            "exit_code": exit_code,
            "timed_out": False,
            "missing_outputs": [],
            "succeeded": succeeded,
            "output_count": len(work),
            "output_files_sha256": IO.canonical_digest(work),
            "stdout_sha256": sha["stdout.txt"],
            "stderr_sha256": sha["stderr.txt"],
        },
        at.keypair_from_seed(KEY_SEED),
        "reference-runner",
    )
    if tamper_record:  # changed after signing
        result["seconds"] = 61.0
    rdir = repo / "records" / label / run_id
    rdir.mkdir(parents=True)
    (rdir / "record.json").write_text(json.dumps(result))
    rb = (rdir / "record.json").read_bytes()
    listing.append(
        {"path": "record.json", "sha256": hashlib.sha256(rb).hexdigest(), "bytes": len(rb)}
    )
    (rdir / "artifacts.json").write_text(
        json.dumps({"run_id": run_id, "sha": code_sha, "files": listing})
    )
    return run_id


def _write_index(repo: Path, label: str, entries: dict[str, str]) -> None:
    keys = ev.run_bundles()
    runs = [
        {"case": keys[b][0], "part": keys[b][1], "run_id": rid, "bundle": b}
        for b, rid in sorted(entries.items())
    ]
    (repo / "records" / label).mkdir(parents=True, exist_ok=True)
    (repo / "records" / label / "index.json").write_text(json.dumps({"runs": runs}, indent=1))


def _write_smokes(
    base: dict[str, Any],
    repo: Path,
    label: str,
    *,
    split: str = "",
    stop: str = "",
    deficit: str = "",
    only: tuple[str, ...] = ev.ALL_SMOKES,
    over: dict[str, dict[str, Any]] | None = None,
) -> None:
    """The four smokes (``<label>-smokes``), run at the pin-A commit between the two reviews.
    Step-count smokes: ``split`` / ``stop``: that smoke's primary 7 takes two steps / stops in
    the slab; ``deficit``: that smoke's phase space misses a row. ``over``: per-smoke keywords."""
    entries: dict[str, str] = {}
    for name in only:
        info = base["cases"][name]["i3fu2"]
        n = int(base["cases"][name]["histories"])
        kw: dict[str, Any] = {
            "code_sha": base["pin_a_sha"],
            "started": SMOKE_STARTED,
            "finished": SMOKE_FINISHED,
        }
        if name in ev.STEP_SMOKES:
            kw["stdout"] = _step_stdout(
                n, float(info["slab_mm"]), split=7 if name == split else -1,
                stop=7 if name == stop else -1,
            )  # fmt: skip
            kw["deficit"] = 1 if name == deficit else 0
        kw.update((over or {}).get(name, {}))
        rng = np.random.default_rng(len(entries))
        entries[name] = _write_run(base, repo, f"{label}-smokes", name, 0.0, rng, **kw)
    (repo / "records" / f"{label}-smokes").mkdir(parents=True, exist_ok=True)
    (repo / "records" / f"{label}-smokes" / "index.json").write_text(
        json.dumps({"runs": [{"bundle": b, "run_id": r} for b, r in sorted(entries.items())]})
    )


def _scenario(
    base: dict[str, Any],
    label: str,
    r_of: Callable[[int, int], float],
    *,
    noise: float = ROW_NOISE,
    run_kw: dict[str, dict[str, Any]] | None = None,
    repeats: dict[str, float] | None = None,
    plant_g: bool = True,
    smoke_kw: dict[str, Any] | None = None,
) -> str:
    """All 40 evidence runs (r_of(case, part) planted relative to the model slab loss), plus
    the ``repeats`` {repeat bundle: r} and the four smokes; ``run_kw`` per bundle."""
    repo, entries = base["repo"], {}
    run_kw = run_kw or {}
    keys = ev.run_bundles()
    todo = [(gen.bundle_name(c, p), r_of(c, p)) for c, p in gen.evidence_keys()]
    todo += sorted((repeats or {}).items())
    for b, r in todo:
        c, p = keys[b]
        rng = np.random.default_rng([20261116, c, p])  # per run: identical across scenarios
        kw = {"noise": noise, "plant_g": plant_g, **run_kw.get(b, {})}
        entries[b] = _write_run(base, repo, label, b, r, rng, **kw)
    _write_index(repo, label, entries)
    _write_smokes(base, repo, label, **(smoke_kw or {}))
    return label


def _evaluate(base: dict[str, Any], label: str, repo: Path | None = None) -> dict[str, Any]:
    return dict(
        ev.evaluate(
            repo or base["repo"],
            f"records/{label}",
            base["root"] / f"out-{label}",
            ev.REVIEW_PIN,
            ev.MANIFEST,
            base["archive"],
            base["code_sha"],
            smokes=f"records/{label}-smokes",
            archive_sha256=base["sha"],
            member_sha256=base["members"],
        )
    )


S = TABLE_FACTOR - 1.0


def _r_icru(c: int, p: int) -> float:
    """TOPAS = the table, except C 10-15 MeV/u 3 % above it (class H) and He = model."""
    if p == 4:
        return 0.0
    if p in (0, 1, 5) and c in (17, 18):
        return TABLE_FACTOR * 1.03 - 1.0
    return S


def _r_model(c: int, p: int) -> float:
    return 0.0


def _r_mixed(c: int, p: int) -> float:
    return S if p in (0, 1, 5) else 0.0  # C sides with the table, O with the model


def _r_he_fail(c: int, p: int) -> float:
    return 0.02 if p == 4 else _r_icru(c, p)


def _r_fragile(c: int, p: int) -> float:
    """TOPAS 0.4 % above the table everywhere: I at the declared systematic, N when doubled."""
    return 0.0 if p == 4 else TABLE_FACTOR * 1.004 - 1.0


R0, B_LIN = 0.05, 0.1  # Richardson scenario: r = R0 + B_LIN f exactly (noise free, no g)


def _r_linear(c: int, p: int) -> float:
    if p == 4:
        return 0.0
    return R0 + B_LIN * (0.048 if p in (1, 3) else 0.03)


def _r_inconsistent(c: int, p: int) -> float:
    if p == 4:
        return 0.0
    return S + (0.01 if (c, p) == (19, 3) else 0.0)  # O 20 MeV/u: Delta_meas off by ~0.01


BAD = {"deficit": 1}


def _both_lines(n: int, fluct: bool) -> str:
    """A stdout showing BOTH fluctuation dump lines (a marker set that cannot discriminate)."""
    return _stdout(n, fluct=fluct) + f"Fluctuations of dE/dx are enabled   {int(not fluct)}\n"


RESULT_LABELS = (
    "icru", "model", "mixed", "he_fail", "linear", "invalid", "second_invalid", "repeat_ok",
    "split", "stopped", "smoke_deficit", "off_hits_on", "on_hits_off", "fragile",
)  # fmt: skip


@pytest.fixture(scope="module")
def results(base: dict[str, Any]) -> dict[str, dict[str, Any]]:
    c15, o19_2t = gen.bundle_name(15, 0), gen.bundle_name(19, 3)
    _scenario(base, "icru", _r_icru)
    _scenario(base, "model", _r_model)
    _scenario(base, "mixed", _r_mixed)
    _scenario(base, "he_fail", _r_he_fail)
    _scenario(base, "linear", _r_linear, noise=0.0, plant_g=False)
    _scenario(base, "invalid", _r_icru, run_kw={c15: BAD})
    _scenario(  # C 5 MeV/u invalid + an invalid repeat; O 20 MeV/u inconsistent, no repeat
        base,
        "second_invalid",
        _r_inconsistent,
        run_kw={c15: BAD, gen.bundle_name(15, 8): BAD},
        repeats={gen.bundle_name(15, 8): S},
    )
    _scenario(  # both restored: C part 8 for the invalid run; O 20 MeV/u BOTH slabs (10, 11)
        base,
        "repeat_ok",
        _r_inconsistent,
        run_kw={c15: BAD},
        repeats={gen.bundle_name(15, 8): S, gen.bundle_name(19, 10): S, gen.bundle_name(19, 11): S},
    )
    _scenario(base, "split", _r_icru, smoke_kw={"split": ev.STEP_SMOKES[1]})
    _scenario(base, "stopped", _r_icru, smoke_kw={"stop": ev.STEP_SMOKES[0]})
    _scenario(base, "smoke_deficit", _r_icru, smoke_kw={"deficit": ev.STEP_SMOKES[1]})
    _scenario(  # the fluctuation-off marker also matches a fluctuation-ON run
        base,
        "off_hits_on",
        _r_icru,
        run_kw={gen.bundle_name(17, 0): {"stdout": _both_lines(SYN_HISTORIES, True)}},  # fmt: skip
    )
    _scenario(  # the fluctuation-on marker also matches a fluctuation-OFF run
        base,
        "on_hits_off",
        _r_icru,
        run_kw={gen.bundle_name(21, 5): {"stdout": _both_lines(SYN_HISTORIES, False)}},  # fmt: skip
    )
    _scenario(base, "fragile", _r_fragile)
    _commit(base["repo"], "records")
    out = {lab: _evaluate(base, lab) for lab in RESULT_LABELS}
    assert o19_2t in out["second_invalid"]["provenance"]["runs"]
    return out


def _cls(doc: dict[str, Any], ion: str) -> dict[str, str]:
    return {b: v["class"] for b, v in doc["bands"][ion].items()}


def test_icru_scenario_with_h_band(results: dict[str, dict[str, Any]], base: dict) -> None:
    doc = results["icru"]
    assert doc["pattern"] == doc["pattern_before_controls"] == "ICRU" and doc["conclusive"]
    assert _cls(doc, "carbon") == {"[2,5)": "I", "[5,10)": "I", "[10,20)": "H", "[20,50]": "I"}
    assert set(_cls(doc, "oxygen").values()) == {"I"}
    for ion in ("carbon", "oxygen"):
        for b, v in doc["bands"][ion].items():
            assert v["mean_s"] == pytest.approx(S, abs=4e-3)
            assert v["x_b"] == 0.01
            assert v["n"] == {"[2,5)": 2, "[5,10)": 2, "[10,20)": 2, "[20,50]": 3}[b]
            assert v["sigma_b"] == pytest.approx(np.hypot(v["se_b"], 0.0025))
            assert 1e-4 < v["se_b"] < 6e-4  # 0.4 % row noise, 1000 rows, extrapolated points
            assert v["max_abs_consistency_residual"] < 2e-3  # planted g: Delta_meas ~ Delta_pred
            assert 0.0 <= v["max_abs_model_residual"] < 3.5e-4
            assert v["max_abs_predicted_difference"] < 3.5e-3
            if (ion, b) != ("carbon", "[10,20)"):
                assert v["mean_r_I"] == pytest.approx(0.0, abs=4e-3)
    assert doc["bands"]["carbon"]["[10,20)"]["mean_r_I"] == pytest.approx(0.03, abs=4e-3)
    assert doc["robust"] is True and doc["pattern_sigma_doubled"] == "ICRU"
    assert doc["classes_sigma_doubled"] == {ion: _cls(doc, ion) for ion in ("carbon", "oxygen")}
    c = doc["controls"]
    assert c == {
        "helium": {"pass": True, "max_abs_r_m": c["helium"]["max_abs_r_m"]},
        "fluctuations_off": {
            "pass": True,
            "max_abs_difference": c["fluctuations_off"]["max_abs_difference"],
            "markers_discriminate": True,
        },
        "controls_pass": True,
    }
    assert c["helium"]["max_abs_r_m"] <= 2e-3 and c["fluctuations_off"]["max_abs_difference"] < 2e-3
    assert doc["invalid_runs"] == {} and doc["invalid_points"] == {} and doc["repeats_used"] == {}
    assert doc["inconclusive_reasons"] == []
    one = {
        "primaries": 50,
        "max_steps_in_slab": 1,
        "traversing_primaries": 50,
        "phase_space_complete": True,
    }
    assert doc["one_step_premise"] == {
        "checked": True,
        "smokes": dict.fromkeys(ev.STEP_SMOKES, one),
        "pass": True,
    }
    prov = doc["provenance"]
    runs = prov["runs"]
    assert len(runs) == 40 and set(runs) == {gen.bundle_name(c, p) for c, p in gen.evidence_keys()}
    assert all(
        v["scored"] == SYN_HISTORIES and v["run_id"].startswith("REF-") for v in runs.values()
    )
    assert set(prov["smokes"]) == set(ev.ALL_SMOKES)
    assert {v["code_sha"] for v in prov["smokes"].values()} == {base["pin_a_sha"]}
    pins = prov["review_pins"]
    assert pins["evaluator"]["review_id"] == "REVIEW-syn-a"
    assert pins["evaluator"]["head_sha"] == base["reviewed"]
    assert pins["markers"]["head_sha"] == base["markers_sha"]
    assert pins["markers"]["finished_at"] == PIN_B_FINISHED
    assert prov["manifest"]["path"] == ev.MANIFEST
    assert prov["log_markers"]["path"] == ev.MARKERS
    assert "I3-FU-1 step-2 review pins" in prov["plan_note"]["line"]
    sealed = prov["sealed_sources"]
    assert set(sealed) == set(ev.executable_closure()) | set(ev.SEALED_FILES) | {ev.ENV_SEAL}
    assert "validation/scripts/stopping/env_seal.py" in sealed
    for rel, digest in sealed.items():
        source = (base["repo"] if rel == ev.ENV_SEAL else ROOT) / rel
        assert digest == hashlib.sha256(source.read_bytes()).hexdigest()
    env = prov["environment"]
    assert env["lock"]["numpy"] == env["lock"]["numpy_locked"] == np.__version__
    assert env["lock"]["python"] == env["lock"]["python_locked"]
    seal = env["seal"]
    assert seal["path"] == ev.ENV_SEAL and seal["match"] is True
    assert seal["numpy_version"] == np.__version__ and seal["numpy_record_files_verified"] > 100
    assert seal["seal_sha256"] == sealed[ev.ENV_SEAL]
    for name, signed in (("evaluator", ATT_A_SIGNED), ("markers", ATT_B_SIGNED)):
        att = pins[name]["attestation"]
        assert att["signed_at_utc"] == signed and att["validation_id"].startswith("VAL-syn-")
    assert doc["amendment"] == "I3-FU-1 step 2 (revision 2)"


def test_model_scenario(results: dict[str, dict[str, Any]]) -> None:
    doc = results["model"]
    assert doc["pattern"] == "MODEL" and doc["conclusive"] and doc["robust"]
    for ion in ("carbon", "oxygen"):
        for v in doc["bands"][ion].values():
            assert v["class"] == "M"
            assert v["mean_r_M"] == pytest.approx(0.0, abs=2e-3)
            assert v["mean_r_I"] == pytest.approx(-S / TABLE_FACTOR, abs=4e-3)


def test_mixed_scenario(results: dict[str, dict[str, Any]]) -> None:
    doc = results["mixed"]
    assert doc["pattern"] == "MIXED" and doc["conclusive"] and doc["robust"]
    assert set(_cls(doc, "carbon").values()) == {"I"}
    assert set(_cls(doc, "oxygen").values()) == {"M"}


def test_pattern_not_robust_to_the_doubled_systematic_is_inconclusive(
    results: dict[str, dict[str, Any]],
) -> None:
    """REVIEW-aa4e7b1b item 4: TOPAS 0.4 % above the table is I at sigma_sys = 0.25 % but N at
    0.5 %: the pattern ICRU is not robust and the result is INCONCLUSIVE."""
    doc = results["fragile"]
    assert doc["pattern_before_controls"] == "ICRU" and doc["pattern_sigma_doubled"] == "MIXED"
    assert doc["robust"] is False and doc["pattern"] == "INCONCLUSIVE"
    assert doc["inconclusive_reasons"] == ["not_robust_to_systematic"] == [ev.ROBUST_REASON]
    for ion in ("carbon", "oxygen"):
        for b in ("[2,5)", "[5,10)", "[10,20)"):
            assert doc["bands"][ion][b]["class"] == "I"
            assert doc["classes_sigma_doubled"][ion][b] == "N"


def test_failed_control_is_inconclusive(results: dict[str, dict[str, Any]]) -> None:
    doc = results["he_fail"]
    assert doc["pattern_before_controls"] == "ICRU" and doc["pattern"] == "INCONCLUSIVE"
    assert not doc["conclusive"] and doc["inconclusive_reasons"] == [
        "a control failed or is incomplete"
    ]
    he = doc["controls"]["helium"]
    assert he["pass"] is False and he["max_abs_r_m"] == pytest.approx(0.02, abs=2e-3)
    assert doc["controls"]["fluctuations_off"]["pass"] is True


def test_richardson_removes_a_planted_linear_bias(results: dict[str, dict[str, Any]]) -> None:
    doc = results["linear"]
    assert doc["conclusive"]
    for ion in ("carbon", "oxygen"):
        for v in doc["bands"][ion].values():
            assert v["mean_r_M"] == R0  # written rounded; the planted value is exact at 0.1 pp
            assert v["se_b"] < 1e-9
            assert v["max_abs_thickness_difference"] == pytest.approx(B_LIN * 0.018, abs=1e-9)
    # fluctuation control: same planted r as C slab t -> zero difference
    assert doc["controls"]["fluctuations_off"]["max_abs_difference"] < 1e-9


def test_invalid_run_without_repeat_is_inconclusive(results: dict[str, dict[str, Any]]) -> None:
    doc = results["invalid"]
    c15 = gen.bundle_name(15, 0)
    assert doc["pattern"] == "INCONCLUSIVE" and doc["pattern_before_controls"] == "ICRU"
    assert list(doc["invalid_runs"]) == [c15] and "deficit" in doc["invalid_runs"][c15]
    assert doc["invalid_points"] == {"carbon case 15": "invalid run; no repeat"}
    assert doc["bands"]["carbon"]["[5,10)"]["n"] == 1
    assert doc["inconclusive_reasons"] == ["an evidence point is invalid without a valid repeat"]
    assert doc["provenance"]["runs"][c15]["scored"] is None


def test_second_invalid_and_consistency_failure(results: dict[str, dict[str, Any]]) -> None:
    doc = results["second_invalid"]
    assert doc["pattern"] == "INCONCLUSIVE" and not doc["conclusive"]
    assert set(doc["invalid_runs"]) == {gen.bundle_name(15, 0), gen.bundle_name(15, 8)}
    assert doc["invalid_points"] == {
        "carbon case 15": "repeat run invalid (second failure of the point)",
        "oxygen case 19": f"{ev.CONSISTENCY_REASON}; no repeat",
    }
    assert list(doc["repeats_used"]) == [gen.bundle_name(15, 8)]
    o = doc["bands"]["oxygen"]["[20,50]"]  # 20 MeV/u
    assert o["n"] == 2 and o["max_abs_consistency_residual"] == pytest.approx(0.01, abs=2e-3)


def test_valid_repeats_restore_conclusiveness(results: dict[str, dict[str, Any]]) -> None:
    doc = results["repeat_ok"]
    assert doc["conclusive"] and doc["pattern"] == "ICRU"
    assert doc["invalid_points"] == {} and list(doc["invalid_runs"]) == [gen.bundle_name(15, 0)]
    cons = f"{ev.CONSISTENCY_REASON}"
    assert doc["repeats_used"] == {
        gen.bundle_name(15, 8): f"original {gen.bundle_name(15, 0)} invalid: "
        + doc["invalid_runs"][gen.bundle_name(15, 0)],
        gen.bundle_name(19, 10): f"original {gen.bundle_name(19, 2)} invalid: {cons}",
        gen.bundle_name(19, 11): f"original {gen.bundle_name(19, 3)} invalid: {cons}",
    }
    assert (
        doc["bands"]["carbon"]["[5,10)"]["n"] == 2 and doc["bands"]["oxygen"]["[20,50]"]["n"] == 3
    )
    assert doc["bands"]["oxygen"]["[20,50]"]["max_abs_consistency_residual"] < 2e-3
    assert len(doc["provenance"]["runs"]) == 43


@pytest.mark.parametrize(
    "label, smoke, want",
    [
        ("split", 1, {"max_steps_in_slab": 2, "traversing_primaries": 49}),  # a step count of 2
        ("stopped", 0, {"max_steps_in_slab": 1, "traversing_primaries": 49}),  # truncated track
        ("smoke_deficit", 1, {"phase_space_complete": False}),  # a missing phase-space row
    ],
)
def test_one_step_premise_failures_are_inconclusive(
    results: dict[str, dict[str, Any]], label: str, smoke: int, want: dict[str, Any]
) -> None:
    doc = results[label]
    assert doc["pattern"] == "INCONCLUSIVE" and doc["pattern_before_controls"] == "ICRU"
    assert doc["inconclusive_reasons"] == [ev.PREMISE_REASON] == ["one_step_premise"]
    prem = doc["one_step_premise"]
    assert prem["checked"] and not prem["pass"]
    bad = prem["smokes"][ev.STEP_SMOKES[smoke]]
    assert {k: bad[k] for k in want} == want and bad["primaries"] == 50
    good = prem["smokes"][ev.STEP_SMOKES[1 - smoke]]
    assert good["traversing_primaries"] == 50 and good["phase_space_complete"]


@pytest.mark.parametrize("label", ["off_hits_on", "on_hits_off"])
def test_non_discriminating_markers_fail_the_fluctuation_control(
    results: dict[str, dict[str, Any]], label: str
) -> None:
    doc = results[label]
    fl = doc["controls"]["fluctuations_off"]
    assert fl["markers_discriminate"] is False and fl["pass"] is False
    assert fl["max_abs_difference"] < 2e-3  # the values alone would have passed
    assert doc["pattern"] == "INCONCLUSIVE" and doc["pattern_before_controls"] == "ICRU"
    assert doc["inconclusive_reasons"] == ["a control failed or is incomplete"]


def _walk(x: Any, path: str = "") -> list[tuple[str, Any]]:
    out = [(path, x)]
    if isinstance(x, dict):
        for k, v in x.items():
            out += _walk(v, f"{path}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            out += _walk(v, f"{path}[{i}]")
    return out


def test_aggregates_only_and_table_relative_values_rounded(
    results: dict[str, dict[str, Any]],
) -> None:
    """The result holds per-ion/band aggregates, control maxima and counts only (no per-point
    list or value) and every table-relative aggregate is written rounded to 0.1 pp."""
    band_keys = {
        "n", "mean_r_M", "mean_r_I", "mean_s", "se_b", "sigma_b", "x_b", "class",
        "max_abs_thickness_difference", "max_abs_consistency_residual",
        "max_abs_predicted_difference", "max_abs_model_residual",
    }  # fmt: skip
    for doc in results.values():
        assert doc["rounding_note"] == ev.ROUNDING_NOTE
        for path, v in _walk(doc):
            if isinstance(v, list):
                assert len(v) <= 5, path
        assert set(doc["bands"]) == {"carbon", "oxygen"}
        for ion in doc["bands"].values():
            assert set(ion) == {"[2,5)", "[5,10)", "[10,20)", "[20,50]"}
            for band in ion.values():
                assert set(band) == band_keys
                for key in ("mean_r_M", "mean_r_I", "mean_s", "x_b"):
                    v = band[key]
                    assert v is None or round(v, 3) == v, (key, v)
        assert set(doc["classes_sigma_doubled"]) == {"carbon", "oxygen"}
        he = doc["controls"]["helium"]["max_abs_r_m"]
        assert he is None or round(he, 3) == he
        for run in doc["provenance"]["runs"].values():
            assert set(run) == {"run_id", "seed", "histories", "scored"}
        assert set(doc["controls"]["helium"]) == {"pass", "max_abs_r_m"}
        assert set(doc["controls"]["fluctuations_off"]) == {
            "pass", "max_abs_difference", "markers_discriminate",
        }  # fmt: skip
        prem = doc["one_step_premise"]
        assert set(prem) == {"checked", "smokes", "pass"}
        for s in prem["smokes"].values():
            assert set(s) == {
                "primaries", "max_steps_in_slab", "traversing_primaries", "phase_space_complete",
            }  # fmt: skip
        for reason in [*doc["invalid_runs"].values(), *doc["invalid_points"].values()]:
            rest = reason.replace(ev.CONSISTENCY_REASON, "")
            assert not re.search(r"\d\.\d{3,}", rest), reason  # reasons carry no values
        text = json.dumps(doc)
        for bad in ("r_i_per", 'delta_e"', "de_m", "de_i", "exit_energy", "differences", '"g1"'):
            assert bad not in text.lower()


def test_classes_use_full_precision_values() -> None:
    """A band whose rounded mean would flip the class keeps the full-precision class."""
    pts = [
        {"e0": 2.4, "valid": True, "r_M": 0.0, "r_I": 0.0, "s": 0.0304, "se_rel": 0.0,
         "d_thick": None, "d_resid": None, "d_pred": 0.0, "r_model": 0.0},
    ]  # fmt: skip
    band = ev.band_aggregates(pts)["[2,5)"]
    # X = 0.0304 / 3 = 0.010133 -> capped 0.01; |r_I| + k = 0.0041 <= X: I
    assert band["class"] == "I" and band["mean_s"] == 0.03 and band["x_b"] == 0.01
    pts[0]["r_I"] = 0.0056  # full: 0.0097 <= 0.01 (I); rounded 0.006: 0.0101 > 0.01 (N)
    band = ev.band_aggregates(pts)["[2,5)"]
    assert band["mean_r_I"] == 0.006 and band["class"] == "I"
    assert ev.band_classes(pts, 0.0025)["[2,5)"] == "I"
    assert ev.band_classes(pts, 0.005)["[2,5)"] == "M"  # the robustness re-classification
    assert ev.band_classes([], 0.0025) == dict.fromkeys(("[2,5)", "[5,10)", "[10,20)", "[20,50]"),
                                                        "none")  # fmt: skip


# --- authentication, sealing and binding ------------------------------------------------------


def _world(
    base: dict[str, Any],
    tmp: Path,
    bundle: str = "",
    *,
    r: float = 0.0,
    commit: bool = True,
    smoke_kw: dict[str, Any] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """A copy of the base repository with ONE run of ``bundle`` (default C 10 MeV/u slab t) and
    the four smokes."""
    repo = tmp / "repo"
    shutil.copytree(base["repo"], repo)
    bundle = bundle or gen.bundle_name(17, 0)
    label = f"t-{tmp.name}"
    w = {**base, "root": tmp, "repo": repo}
    rid = _write_run(w, repo, label, bundle, r, np.random.default_rng(3), **kw)
    _write_index(repo, label, {bundle: rid})
    _write_smokes(w, repo, label, **(smoke_kw or {}))
    if commit:
        _commit(repo, "records")
    return {"repo": repo, "label": label, "bundle": bundle, "run_id": rid, "root": tmp}


def _ctx(base: dict[str, Any], w: dict[str, Any], code_sha: str | None = None) -> Any:
    return ev.context(
        w["repo"],
        f"records/{w['label']}",
        w["root"] / f"out-{w['label']}",
        ev.REVIEW_PIN,
        ev.MANIFEST,
        code_sha or base["code_sha"],
        smokes=f"records/{w['label']}-smokes",
    )


def _check(base: dict[str, Any], w: dict[str, Any]) -> dict[str, Any]:
    return dict(ev.check_run(_ctx(base, w), w["bundle"], w["run_id"]))


def test_authenticated_run_is_valid(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path)
    out = _check(base, w)
    assert out["valid"] and out["scored"] == SYN_HISTORIES and out["run_id"] == w["run_id"]
    g = _g(base["cases"][w["bundle"]]["i3fu2"])
    assert out["de"] == pytest.approx(base["de_m"][w["bundle"]] * (1 + g), rel=3e-3)
    assert out["markers"] == (True, False)


@pytest.mark.parametrize(
    "kw, match",
    [
        ({"tamper_record": True}, "attestation"),
        ({"started": "2026-11-16T11:59:59+00:00"}, "not after both review pins"),
        ({"started": "2026-11-16T12:15:00+00:00"}, "not after both review pins"),  # A < t < B
        ({"succeeded": False}, "not a successful topas run"),
        ({"exit_code": 1}, "not a successful topas run"),
        (
            {"stdout": _stdout(SYN_HISTORIES).replace("Registered Physics Processes:", "")},
            "ListProcesses",
        ),
        (
            {"stdout": _stdout(SYN_HISTORIES).replace("Use ICRU90 data", "Use ICRU90")},
            "g4em-standard_opt4",
        ),
        ({"stdout": _stdout(SYN_HISTORIES, particle="alpha")}, "ionIoni"),
    ],
)
def test_tampering_raises(
    base: dict[str, Any], tmp_path: Path, kw: dict[str, Any], match: str
) -> None:
    w = _world(base, tmp_path, **kw)
    with pytest.raises(ev.RecordError, match=match) as exc:
        _check(base, w)
    assert w["bundle"] in str(exc.value) and w["run_id"] in str(exc.value)  # names the run


def test_altered_artifact_bytes_raise(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path)
    edep = w["root"] / f"out-{w['label']}" / w["run_id"] / "work" / "edep.csv"
    edep.write_text(edep.read_text().replace("0, 0, 0, ", "0, 0, 0, 1"))
    with pytest.raises(ev.RecordError, match="sha256 mismatch"):
        _check(base, w)


def test_markers_matching_the_native_input_raise(base: dict[str, Any], tmp_path: Path) -> None:
    """Reviewer pass 3 item C2: a fluctuation marker that matches the TOPAS parameter echo (the
    bundle's input text) does not identify the G4EmParameters printout."""
    w = _world(base, tmp_path)
    for key, rx in (
        ("fluct_on_marker", r"Ph/Default/Modules"),
        ("fluct_off_marker", r"EnergyLossFluctuations|Ph/Default/Type"),
    ):
        ctx = _ctx(base, w)
        ctx.markers = {**ctx.markers, key: re.compile(rx, re.M)}
        with pytest.raises(ev.RecordError, match=f"{key} matches the bundle's native input"):
            ev.check_run(ctx, w["bundle"], w["run_id"])
        smoke = ev.ALL_SMOKES[0]
        with pytest.raises(ev.RecordError, match=f"{key} matches the bundle's native input"):
            ev.check_smoke(ctx, smoke, ctx.smoke_index[smoke])


def test_wrong_code_sha_raises(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path, code_sha=base["markers_sha"])  # run executed before pin B
    with pytest.raises(ev.RecordError, match="code_sha"):
        _check(base, w)
    with pytest.raises(ev.RecordError, match="found 0|is not committed at"):
        _ctx(base, w, code_sha=base["reviewed"])  # no plan note / pins at the reviewed commit
    with pytest.raises(ev.RecordError, match="40-hex"):
        _ctx(base, w, code_sha="abc")


def test_manifest_and_pin_mismatch_raise(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path / "m")
    man = w["repo"] / ev.MANIFEST
    man.write_text(man.read_text() + " ")
    _commit(w["repo"], "manifest changed after the runs")
    with pytest.raises(ev.RecordError, match="frozen manifest .* differs"):
        _ctx(base, w)
    for rel in (ev.REVIEW_PIN, ev.MARKERS_PIN):
        w = _world(base, tmp_path / Path(rel).stem)
        pin = w["repo"] / rel
        pin.write_text(pin.read_text().replace("2026-11-16T12", "2026-11-16T10"))
        _commit(w["repo"], "pin forged later")
        with pytest.raises(ev.RecordError, match="review pin .* differs"):
            _ctx(base, w)
    w = _world(base, tmp_path / "mk")
    mk = w["repo"] / ev.MARKERS
    mk.write_text(json.dumps({**TEST_MARKERS, "fluct_off_marker": "x"}))
    _commit(w["repo"], "markers changed after the runs")
    with pytest.raises(ev.RecordError, match="log-marker file .* differs"):
        _ctx(base, w)
    w = _world(base, tmp_path / "u", commit=False)  # uncommitted records: fail closed
    with pytest.raises(ev.RecordError, match="not committed|differs from HEAD"):
        _check(base, w)


def test_executable_closure_is_computed_from_the_imports() -> None:
    closure = ev.executable_closure()
    for rel in (
        "validation/scripts/stopping/i3fu2_evaluate.py",
        "validation/scripts/stopping/gen_i3fu2_bundles.py",
        "validation/scripts/stopping/compare_icru73.py",
        "validation/scripts/ion_em/r2_io.py",
        "validation/scripts/ion_em/ss_evaluate.py",
        "src/ionmc/physics/stopping.py",
        "src/ionmc/data/nist_star.py",  # imported by stopping.py, never hand-listed
        "src/ionmc/__init__.py",
        "infrastructure/experiment_v3/attest.py",
    ):
        assert rel in closure, rel
    assert all(p.is_file() and ROOT.resolve() in p.parents for p in closure.values())
    assert not any(rel.startswith("tests/") for rel in closure)


def test_sealed_sources_committed_or_uncommitted_change_raises(
    base: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A later committed change of a closure file (also one no list names: nist_star.py), of
    uv.lock or pyproject.toml, an uncommitted change, and a loaded module that differs from the
    worktree each refuse."""
    for rel in (
        "src/ionmc/physics/stopping.py",
        "validation/scripts/stopping/i3fu2_evaluate.py",
        "src/ionmc/data/nist_star.py",
        "uv.lock",
        "pyproject.toml",
    ):
        w = _world(base, tmp_path / ("c-" + Path(rel).name))
        (w["repo"] / rel).write_text((w["repo"] / rel).read_text() + "# later change\n")
        _commit(w["repo"], "sealed source changed after the runs")
        # a closure module is caught by the loaded-vs-worktree check, uv.lock / pyproject.toml
        # by the commit comparison
        want = f"sealed source {re.escape(rel)}(: the loaded module)? differs"
        with pytest.raises(ev.RecordError, match=want):
            _ctx(base, w)
    w = _world(base, tmp_path / "dirty")
    rel = "src/ionmc/physics/projectiles.py"
    (w["repo"] / rel).write_text((w["repo"] / rel).read_text() + "# uncommitted\n")
    with pytest.raises(ev.RecordError, match="not clean"):
        _ctx(base, w)
    w = _world(base, tmp_path / "loaded")
    other = tmp_path / "other.py"
    other.write_text("# not the sealed bytes\n")
    real = ev.executable_closure()
    monkeypatch.setattr(
        ev, "executable_closure", lambda: {**real, "src/ionmc/data/nist_star.py": other}
    )
    with pytest.raises(ev.RecordError, match="nist_star.py: the loaded module differs"):
        _ctx(base, w)


def test_environment_must_match_the_lock(
    base: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lock = (ROOT / "uv.lock").read_text()
    env = ev.environment_check(lock)
    assert env["numpy"] == np.__version__ and env["python_locked"] == "3.12"
    with pytest.raises(ev.RecordError, match="not found exactly once"):
        ev.environment_check(lock.replace('name = "numpy"', 'name = "numpyx"'))
    with pytest.raises(ev.RecordError, match="requires-python"):
        ev.environment_check(
            lock.replace('requires-python = ">=3.12"', 'requires-python = ">=3.13"')
        )
    w = _world(base, tmp_path)
    monkeypatch.setattr(np, "__version__", "0.0.1")
    # caught first by the environment seal (the installed distribution is not the imported numpy)
    with pytest.raises(ev.RecordError, match="imported numpy 0.0.1|numpy 0.0.1 differs"):
        _ctx(base, w)
    with pytest.raises(ev.RecordError, match="numpy 0.0.1 differs from the uv.lock pin"):
        ev.environment_check(lock)


def _repin(base: dict[str, Any], tmp: Path, *, pin_a: Any = None, pin_b: Any = None,
           note: str | None = None, att_a: Any = None,
           att_b: Any = None) -> dict[str, Any]:  # fmt: skip
    """A copy of the base repository whose pins / attestations / plan note are rewritten in ONE
    new commit, which then serves as code_sha."""
    w = _world(base, tmp)
    for rel, doc in (
        (ev.REVIEW_PIN, pin_a),
        (ev.MARKERS_PIN, pin_b),
        (ev.REVIEW_ATTESTATION, att_a),
        (ev.MARKERS_ATTESTATION, att_b),
    ):
        if doc is not None:
            _write_json(w["repo"] / rel, doc)
    if note is not None:
        (w["repo"] / ev.PLAN_PATH).write_text("# synthetic acceptance plan\n\n" + note)
    w["code_sha"] = _commit(w["repo"], "rewritten pins")
    return w


@pytest.mark.parametrize(
    "which, over, match",
    [
        ("a", {"status": "failed"}, "status 'failed' is not 'passed'"),
        ("b", {"status": "running"}, "status 'running' is not 'passed'"),
        ("a", {"report": {"verdict": "changes_required"}}, "verdict 'changes_required'"),
        ("b", {"report": {"verdict": "changes_required"}}, "verdict 'changes_required'"),
        ("a", {"report": {"summary": "x"}}, "missing record field"),
        ("a", {"report": None}, "missing record field"),
        ("b", {"finished_at": None}, "missing record field"),
        ("a", {"task_id": "V3-005"}, "task_id"),
        ("b", {"finished_at": "yesterday"}, "ISO-8601"),
        ("a", {"head_sha": "abc"}, "40-hex"),
        ("a", {"report_sha256": "abc"}, "report_sha256 is not a 64-hex"),
        ("b", {"report_sha256": None}, "missing record field"),
        ("a", {"started_at": PIN_A_FINISHED}, "started_at is not before finished_at"),
        ("b", {"started_at": "2026-11-16T13:00:00+00:00"}, "started_at is not before"),
    ],
)
def test_review_record_must_be_a_passed_review(
    base: dict[str, Any], tmp_path: Path, which: str, over: dict[str, Any], match: str
) -> None:
    a = _review("REVIEW-syn-a", base["reviewed"], PIN_A_STARTED, PIN_A_FINISHED)
    b = _review("REVIEW-syn-b", base["markers_sha"], PIN_B_STARTED, PIN_B_FINISHED)
    if which == "a":
        a.update(over)
    else:
        b.update(over)
    w = _repin(base, tmp_path, pin_a=a if which == "a" else None, pin_b=b if which == "b" else None)
    with pytest.raises(ev.RecordError, match=match):
        _ctx(base, w, code_sha=w["code_sha"])


def test_pins_must_be_named_by_the_plan_note(base: dict[str, Any], tmp_path: Path) -> None:
    ok = ("REVIEW-syn-a", base["reviewed"], "REVIEW-syn-b", base["markers_sha"])
    for i, (bad, match) in enumerate(
        (
            (("REVIEW-other", *ok[1:]), "not the evaluator review named"),
            ((*ok[:3], base["reviewed"]), "not the markers review named"),
            ((ok[0], base["markers_sha"], *ok[2:]), "not the evaluator review named"),
        )
    ):
        w = _repin(base, tmp_path / f"n{i}", note=_note(*bad))
        with pytest.raises(ev.RecordError, match=match):
            _ctx(base, w, code_sha=w["code_sha"])
    w = _repin(base, tmp_path / "none", note="no note\n")
    with pytest.raises(ev.RecordError, match="found 0"):
        _ctx(base, w, code_sha=w["code_sha"])
    w = _repin(base, tmp_path / "two", note=_note(*ok) + _note(*ok))
    with pytest.raises(ev.RecordError, match="found 2"):
        _ctx(base, w, code_sha=w["code_sha"])
    # the note at code_sha differs from HEAD's (rewritten after the runs)
    w = _world(base, tmp_path / "later")
    (w["repo"] / ev.PLAN_PATH).write_text(
        "# synthetic acceptance plan\n\n" + _note("REVIEW-syn-a", ok[1], "REVIEW-syn-b", ok[1])
    )
    _commit(w["repo"], "plan note rewritten after the runs")
    with pytest.raises(ev.RecordError, match="plan note differs between code_sha"):
        _ctx(base, w)


def test_bundle_entry_not_in_manifest_raises(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path)
    idx = w["repo"] / "records" / w["label"] / "index.json"
    doc = json.loads(idx.read_text())
    doc["runs"][0]["bundle"] = gen.SMOKE_NAME
    idx.write_text(json.dumps(doc))
    _commit(w["repo"], "index points at a smoke")
    with pytest.raises(ev.RecordError, match="not i3fu2"):
        _ctx(base, w)


def test_log_validity_failure_is_invalid_not_error(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path, stdout=_stdout(SYN_HISTORIES) + "G4Exception: track stuck\n")
    out = _check(base, w)
    assert out["valid"] is False and "killed/stuck" in out["reason"]


def test_repeat_without_recorded_failure_raises(base: dict[str, Any], tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    shutil.copytree(base["repo"], repo)
    w = {**base, "root": tmp_path, "repo": repo}
    _scenario(w, "extra", _r_model, repeats={gen.bundle_name(13, 8): 0.0})
    _commit(repo, "records")
    with pytest.raises(ev.RecordError, match="without a recorded failure"):
        _evaluate(w, "extra", repo)


# --- smokes: completeness and chronology (Codex REVIEW-aa4e7b1b finding 3) --------------------


def _smoke_check(base: dict[str, Any], tmp: Path, name: str, **kw: Any) -> dict[str, Any]:
    """A copy of the base repository with the four smokes (``kw`` for smoke ``name``); checks
    that smoke."""
    w = _world(base, tmp, smoke_kw={"over": {name: kw}})
    ctx = _ctx(base, w)
    return dict(ev.check_smoke(ctx, name, ctx.smoke_index[name]))


def test_smokes_authenticated_with_markers(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path)
    ctx = _ctx(base, w)
    out = {n: ev.check_smoke(ctx, n, rid) for n, rid in ctx.smoke_index.items()}
    assert set(out) == set(ev.ALL_SMOKES)
    assert out["i3fu2-smoke-carbon-water-10mevu"]["markers"] == (True, False)
    assert out["i3fu2-smoke-carbon-water-10mevu-fluct-off"]["markers"] == (False, True)
    assert out["i3fu2-smoke-carbon-water-10mevu-fluct-off"]["fluctuations"] is False
    for n in ev.STEP_SMOKES:
        assert out[n]["ok"] and out[n]["traversing_primaries"] == 50
        assert out[n]["code_sha"] == base["pin_a_sha"]


@pytest.mark.parametrize(
    "kw, match",
    [
        ({"finished": "2026-11-16T12:25:00+00:00"}, "not before the markers review"),
        ({"finished": PIN_B_FINISHED}, "not before the markers review"),
        ({"started": "2026-11-16T11:55:00+00:00"}, "not after the evaluator review pin"),
        ({"code_sha": "reviewed"}, "smoke .* not committed"),  # does not carry pin A
        ({"code_sha": "code_sha"}, "not an ancestor of the markers review"),  # after pin B
        ({"tamper_record": True}, "attestation"),
        ({"stdout": _stdout(2000).replace("Use ICRU90 data", "Use ICRU90")}, "opt4"),
    ],
)
def test_smoke_chronology_and_binding_raise(
    base: dict[str, Any], tmp_path: Path, kw: dict[str, Any], match: str
) -> None:
    if kw.get("code_sha") in ("reviewed", "code_sha"):
        kw = {**kw, "code_sha": base[kw["code_sha"]]}
    with pytest.raises(ev.RecordError, match=match):
        _smoke_check(base, tmp_path, "i3fu2-smoke-carbon-water-10mevu", **kw)


def test_smoke_at_the_markers_review_head_is_accepted(base: dict[str, Any], tmp_path: Path) -> None:
    out = _smoke_check(
        base, tmp_path, ev.STEP_SMOKES[0], code_sha=base["markers_sha"],
        started="2026-11-16T12:03:00+00:00", finished="2026-11-16T12:19:00+00:00",
    )  # fmt: skip
    assert out["ok"] and out["code_sha"] == base["markers_sha"]


def test_missing_or_late_smokes_raise_in_the_evaluation(
    base: dict[str, Any], tmp_path: Path
) -> None:
    """Regression (REVIEW-aa4e7b1b finding 3): smokes finishing after the markers review began,
    or a missing smoke, stop the evaluation (nothing graded)."""
    repo = tmp_path / "repo"
    shutil.copytree(base["repo"], repo)
    w = {**base, "root": tmp_path, "repo": repo}
    late = {n: {"finished": "2026-11-16T13:30:00+00:00"} for n in ev.ALL_SMOKES}
    _scenario(w, "late", _r_icru, smoke_kw={"over": late})
    _scenario(w, "three", _r_icru, smoke_kw={"only": ev.ALL_SMOKES[:3]})
    _commit(repo, "records")
    with pytest.raises(ev.RecordError, match="not before the markers review"):
        _evaluate(w, "late", repo)
    with pytest.raises(ev.RecordError, match="lacks .*all four smokes are required"):
        _evaluate(w, "three", repo)


# --- review pins bound to signed validation records (Codex REVIEW-c4f3fdd1 finding 1) --------


def _att_variant(base: dict[str, Any], which: str, kind: str) -> dict[str, Any]:
    pin = base["pin_a"] if which == "a" else base["pin_b"]
    signed = ATT_A_SIGNED if which == "a" else ATT_B_SIGNED
    if kind == "altered":  # the review inside the details changed after signing
        doc = _attestation(pin, signed)
        doc["details"]["review_record"] = {**pin, "summary": "forged"}
        return doc
    if kind == "untrusted":
        return _attestation(pin, signed, seed=ROGUE_SEED)
    if kind == "git_sha":
        return _attestation(pin, signed, git_sha="f" * 40)
    if kind == "status":
        return _attestation(pin, signed, status="failed")
    if kind == "task":
        return _attestation(pin, signed, task_id="V3-005")
    if kind == "other_review":  # validly signed, but for a different review record
        return _attestation(pin, signed, details={"review_record": {**pin, "summary": "other"}})
    if kind == "early":  # signed before the review finished
        return _attestation(pin, "2026-11-16T11:00:00+00:00")
    raise ValueError(kind)


@pytest.mark.parametrize("which", ["a", "b"])
@pytest.mark.parametrize(
    "kind, match",
    [
        ("altered", "attestation invalid: payload hash mismatch"),
        ("untrusted", "attestation unverified"),
        ("git_sha", "git_sha is not the pin's head_sha"),
        ("status", "status 'failed' is not passed"),
        ("task", "task_id 'V3-005'"),
        ("other_review", "does not carry the committed review pin"),
        ("early", "not after the review's finished_at"),
    ],
)
def test_pin_attestation_tampering_raises(
    base: dict[str, Any], tmp_path: Path, which: str, kind: str, match: str
) -> None:
    doc = _att_variant(base, which, kind)
    w = _repin(base, tmp_path, att_a=doc) if which == "a" else _repin(base, tmp_path, att_b=doc)
    with pytest.raises(ev.RecordError, match=match):
        _ctx(base, w, code_sha=w["code_sha"])


def test_pin_attestation_layouts_and_missing(base: dict[str, Any], tmp_path: Path) -> None:
    """The listed layout (``{"validation": record}`` with the service's derived fields) is
    accepted; a missing attestation file is a RecordError."""
    rec = _attestation(base["pin_b"], ATT_B_SIGNED)
    wrapped = {"validation": {**rec, "attestation_status": "verified", "attestation_reason": "ok"}}
    w = _repin(base, tmp_path / "wrapped", att_b=wrapped)
    ctx = _ctx(base, w, code_sha=w["code_sha"])
    assert ctx.pins["markers"]["attestation"]["signed_at_utc"] == ATT_B_SIGNED
    w = _world(base, tmp_path / "missing")
    (w["repo"] / ev.MARKERS_ATTESTATION).unlink()
    w["code_sha"] = _commit(w["repo"], "attestation removed")
    with pytest.raises(ev.RecordError, match="not committed"):
        _ctx(base, w, code_sha=w["code_sha"])
    w = _world(base, tmp_path / "moved")  # pin A's record rewritten after the markers review
    _write_json(
        w["repo"] / ev.REVIEW_ATTESTATION, _attestation(base["pin_a"], ATT_A_SIGNED, summary="x")
    )
    w["code_sha"] = _commit(w["repo"], "pin A record re-signed later")
    with pytest.raises(ev.RecordError, match="differs between the markers review"):
        _ctx(base, w, code_sha=w["code_sha"])


def test_runs_must_start_after_the_pin_attestations(base: dict[str, Any], tmp_path: Path) -> None:
    """Signed after the runs: smokes must start after pin A's signed_at_utc, evidence runs after
    pin B's (the pins' signed times are set later than the synthetic runs)."""
    w = _world(base, tmp_path)
    ctx = _ctx(base, w)
    ctx.pins["evaluator"]["signed"] = IO._parse_time("2026-11-16T12:07:00+00:00", "t")
    smoke = ev.STEP_SMOKES[0]
    with pytest.raises(ev.RecordError, match="not after the evaluator review attestation"):
        ev.check_smoke(ctx, smoke, ctx.smoke_index[smoke])
    ctx = _ctx(base, w)
    ctx.pins["markers"]["signed"] = IO._parse_time("2026-11-16T13:30:00+00:00", "t")
    with pytest.raises(ev.RecordError, match="not after the markers review attestation"):
        ev.check_run(ctx, w["bundle"], w["run_id"])
    assert ev.check_smoke(ctx, smoke, ctx.smoke_index[smoke])["ok"]  # smokes depend on pin A


def test_environment_seal_mismatch_raises(
    base: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    w = _world(base, tmp_path)
    real = ev.env_seal.current_environment()
    for key in ev.env_seal.SEAL_FIELDS:
        monkeypatch.setattr(ev.env_seal, "current_environment", lambda k=key: {**real, k: "x"})
        with pytest.raises(
            ev.RecordError, match=f"environment seal .*differs from the seal.*{key}"
        ):
            _ctx(base, w)
    monkeypatch.setattr(ev.env_seal, "current_environment", lambda: real)
    assert _ctx(base, w).environment["seal"]["match"] is True
