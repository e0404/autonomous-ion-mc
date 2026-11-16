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

import functools
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
CW_SEED = bytes(range(3, 35))  # synthetic codex-worker key (trusted in the temp repos)
ROGUE_SEED = bytes(range(2, 34))  # a codex-worker key nobody trusts
RID_A = "REVIEW-" + "a" * 32  # pin A (evaluator review)
RID_B = "REVIEW-" + "b" * 32  # pin B (markers review)
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
    "step_limiter_marker": r"^Step limiter process: StepLimiter$",
}
LIMITER_LINE = "Step limiter process: StepLimiter\n"
_SYN_PATH = ["/synthetic/src", "/synthetic/scripts", "/synthetic/site-packages"]
VALID_LAUNCH = {  # a launch record as sealed_launch.py + env_seal.capture_launch produce it
    "dont_write_bytecode_flag": 1,
    "pycache_prefix": "/synthetic/empty-prefix",
    "pycache_prefix_exists": True,
    "pycache_prefix_files_at_start": 0,
    "isolated_flag": 1,
    "no_site_flag": 1,
    "meta_path": ["BuiltinImporter", "FrozenImporter", "PathFinder"],
    "machinery_problems": [],
    "sys_path_at_capture": list(_SYN_PATH),
    "launcher_record": {
        "origin": "validation/scripts/stopping/sealed_launch.py",
        "launcher_sha256": hashlib.sha256(
            (ROOT / "validation/scripts/stopping/sealed_launch.py").read_bytes()
        ).hexdigest(),
        "sys_path": list(_SYN_PATH),
        "stdlib_modules_at_launch": {},
        "target": "validation/scripts/stopping/i3fu2_evaluate.py",
        "target_path": str((ROOT / "validation/scripts/stopping/i3fu2_evaluate.py").resolve()),
        "target_sha256": hashlib.sha256(
            (ROOT / "validation/scripts/stopping/i3fu2_evaluate.py").read_bytes()
        ).hexdigest(),
        "sys_modules_at_launch": ["builtins", "sys"],
    },
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


@pytest.fixture(autouse=True, scope="module")
def _bytecode_stub() -> Any:
    """The test interpreter's venv and worktree hold ``__pycache__`` files (the orchestrator
    removes them only before the real run), so the bytecode-cache control is replaced by its
    ``sys.dont_write_bytecode`` part here; the full control is unit-tested in
    ``test_i3_band_diagnostic.py`` and :func:`test_runtime_checks_fail_closed`."""
    real = ev.env_seal.bytecode_check

    def stub(modules: Any) -> dict[str, Any]:
        return {
            "dont_write_bytecode": sys.dont_write_bytecode,
            "modules_without_bytecode_cache": len(modules),
            "test_stub": True,
        }

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(ev.env_seal, "bytecode_check", stub)
        mp.setattr(ev, "_REAL_BYTECODE_CHECK", real, raising=False)
        # the test interpreter is not launched with -B -X pycache_prefix (unit-tested launch
        # control in test_i3_band_diagnostic.py and test_launch_record_is_enforced)
        mp.setattr(ev, "LAUNCH", dict(VALID_LAUNCH))
        yield


BUNDLES = sorted(d for d in TOPAS.iterdir() if d.is_dir() and d.name.startswith("i3fu2-"))
SMOKES = [d for d in BUNDLES if d.name in gen.SMOKES]
STEP_SMOKES = {  # TrackingVerbosity-1 smokes: name -> (case, part whose configuration it runs)
    "i3fu2-smoke-carbon-water-2p4mevu-steps": (13, 1),  # slab 1.6 t: one step expected
    "i3fu2-smoke-carbon-water-50mevu-steps": (21, 1),
    "i3fu2-smoke-carbon-water-2p4mevu-finesteps": (13, 6),  # slab t, MaxStepSize t/20
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
    assert len(BUNDLES) == 121 and len(RUNS) == 116 and len(SMOKES) == 5
    assert [d.name for d in BUNDLES] == gen.bundle_names()
    evidence = sorted(
        [(c, p) for p in (0, 1, 2, 3, 6, 7) for c in range(13, 22)]
        + [(c, 4) for c in (17, 21)]
        + [(c, 5) for c in (13, 21)]
    )
    assert sorted(gen.evidence_keys()) == evidence and len(evidence) == 58
    assert gen.REPEAT_OFFSET == 8
    assert sorted(gen.repeat_keys()) == sorted((c, p + 8) for c, p in evidence)
    assert gen.REPEAT_PARTS == tuple(range(8, 16))  # part + 8 for parts 0-7
    assert gen.RESERVED_PARTS == (16, 17, 18, 19)
    assert gen.FINE_PARTS == (6, 7) and gen.FINE_DIVISOR == 20
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
        if base in gen.FINE_PARTS:  # MaxStepSize = t / 20 in um, exactly
            t_um = info["slab_mm"] * 1000.0
            assert info["max_step_um"] == pytest.approx(t_um / 20, rel=1e-12), d.name
            assert (
                _native(d, "d:Ge/Water/MaxStepSize")
                == f"{gen.max_step_size_um(info['ion'], info['case'])} um"
            )
            assert info["slab_mm"] == gen.slab_mm(info["ion"], info["case"], 1)
        else:
            assert info["max_step_um"] is None
            assert "MaxStepSize" not in (d / "input.txt").read_text()
        if part >= 8:  # a repeat is byte-identical to its original except seed and wording
            orig = _case(TOPAS / gen.bundle_name(info["case"], base))
            for key in ("physics", "source", "geometry", "histories", "scorers", "units"):
                assert _case(d)[key] == orig[key], (d.name, key)
            assert {k: v for k, v in info.items() if k not in ("part", "kind")} == {
                k: v for k, v in orig["i3fu2"].items() if k not in ("part", "kind")
            }
    assert kinds == {"evidence": 58, "repeat": 58}
    for d in SMOKES:
        info = _case(d)["i3fu2"]
        assert info["kind"] == "smoke" and info["ion"] == "carbon"
        if d.name in STEP_SMOKES:
            case, part = STEP_SMOKES[d.name]
            factor = 0.048 if part == 1 else 0.03
            assert (info["case"], info["part"], info["thickness_factor"]) == (case, part, factor)
            assert (info["max_step_um"] is not None) == (part == 6)
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
    assert m["bundle_count"] == 121 and m["file_count"] == 242
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
    assert len(set(used)) == 116 and min(used) == 20619013 and max(used) == 20620521
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
        "i3fu2-smoke-carbon-water-2p4mevu-finesteps": [20600445],
    }
    for s in (20600440, 20600441, 20600442, 20600443, 20600444, 20600445):
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
    for name, (case, part) in STEP_SMOKES.items():
        smoke = (TOPAS / name / "input.txt").read_text().splitlines()
        ref = (TOPAS / gen.bundle_name(case, part) / "input.txt").read_text().splitlines()
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


def _emulated_g4_losses(
    ion: str, case: int, mult: int, rho0: float = 0.0, dlam: float = 0.0
) -> tuple[float, float]:
    """Independent emulation of Geant4 (reviewer pass 4): (one-step loss over slab(mult), 20-step
    linear-branch loss over slab t) with the Geant4 stopping = ionmc's with the base part times
    (1 + rho0 E0 / E) and the effective charge q^2 (E / E0)^dlam (a different charge slope)."""
    proj, _z, a, _p = gen.IONS[ion]
    e0 = float(gen.CASES[case])
    full, base = ev.model_stopping(proj), ev.base_stopping_per_q2(proj)

    def q2(e: float) -> float:
        return float(ev.effective_charge_sq(proj, e) * (e / e0) ** dlam)

    def corr(e: float) -> float:
        return float(full(e) - ev.effective_charge_sq(proj, e) * base(e))

    def base_g4(e: float) -> float:
        return float(base(e) * (1.0 + rho0 * e0 / e))

    t = gen.slab_mm(ion, case, mult)
    q2p = q2(e0)
    table = ev.slab_energy_loss(a * e0, a, t, lambda e: q2p * base_g4(e))
    mid = max(e0 - table / a / 2.0, e0 / 2.0)
    one = table + (t / 10.0) * corr(mid)
    tf, e = gen.slab_mm(ion, case, 1), a * e0
    step = tf / 10.0 / gen.FINE_DIVISOR
    for _ in range(gen.FINE_DIVISOR):
        eu = e / a
        lt = step * q2(eu) * base_g4(eu)
        e -= lt + step * corr(max(eu - lt / a / 2.0, eu / 2.0))
    return one, a * e0 - e


def _model(ion: str, case: int) -> dict[str, Any]:
    return {
        **ev.stepping_bias(ion, case),
        "expect": functools.partial(ev.expected_g_meas, ion, case),
    }


@pytest.mark.parametrize("ion", ["carbon", "oxygen"])
@pytest.mark.parametrize("rho0", [0.0, 0.05, 0.10, 0.15])
def test_same_slab_gate_is_not_confounded_by_the_stopping_ratio(ion: str, rho0: float) -> None:
    """Reviewer pass 4 (c): TOPAS built from the emulated one-step and 20-step algorithms with a
    Geant4 stopping ionmc x (1 + rho0 E0 / E) on the base part (ICRU-like up to +15 %) passes the
    same-slab gate e1 at C/O 2.4 MeV/u; the report-only e2 is not part of the gate."""
    case, proj = 13, gen.IONS[ion][0]
    a, e0 = gen.IONS[ion][2], 2.4
    de_m = ev.predicted_loss(
        f"model:{ion}", ev.model_stopping(proj), a, e0, gen.slab_mm(ion, case, 1)
    )
    one, fine = _emulated_g4_losses(ion, case, 1, rho0)
    rs = [{"r_M": one / de_m - 1.0}, {"r_M": 0.0}, {"r_M": fine / de_m - 1.0}]
    e1 = ev.emulation_residuals(rs, _model(ion, case))[0]
    assert abs(e1) <= 2e-4 < ev.EMULATION_TOL


@pytest.mark.parametrize("ion", ["carbon", "oxygen"])
@pytest.mark.parametrize("dlam", [0.2, -0.2])
def test_same_slab_gate_detects_a_wrong_charge_slope(ion: str, dlam: float) -> None:
    case, proj = 13, gen.IONS[ion][0]
    a, e0 = gen.IONS[ion][2], 2.4
    de_m = ev.predicted_loss(
        f"model:{ion}", ev.model_stopping(proj), a, e0, gen.slab_mm(ion, case, 1)
    )
    one, fine = _emulated_g4_losses(ion, case, 1, 0.1, dlam)
    rs = [{"r_M": one / de_m - 1.0}, {"r_M": 0.0}, {"r_M": fine / de_m - 1.0}]
    assert abs(ev.emulation_residuals(rs, _model(ion, case))[0]) > ev.EMULATION_TOL


def test_fine_step_linear_branch_bias() -> None:
    """b_fine (20 forward-Euler steps): -0.04 to -0.07 %, negative at every C/O point."""
    for ion in ("carbon", "oxygen"):
        for case in gen.ALL_CASES:
            b = ev.stepping_bias(ion, case)["b_fine"]
            assert -7e-4 < b < -3.5e-4, (ion, case, b)


def test_model_residual_above_limit_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(ev, "MODEL_RESIDUAL_LIMIT", 1e-5)
    with pytest.raises(ev.RecordError, match="design premise"):
        ev.stepping_model_check()


# --- repeat rule (deterministic; revision 5: the triple, gate e1) ----------------------------


def _fake(valid: bool = True, reason: str | None = None, r: float = 0.0) -> dict[str, Any]:
    return {"valid": valid, "reason": reason, "r": r}


def _ratios(run: dict[str, Any]) -> dict[str, Any]:
    return {"r_M": run["r"]}


FIXED_MODEL = {"g1": 0.003, "g2": 0.005, "d_pred": 0.002, "r_model": 0.0,
               "expect": lambda r_fine: (0.003, 0.005)}  # fmt: skip


def test_check_formulas() -> None:
    m = FIXED_MODEL
    rs = [{"r_M": 0.0531}, {"r_M": 0.05525}, {"r_M": 0.05}]
    g1, g2 = 1.0531 / 1.05 - 1, 1.05525 / 1.05 - 1
    assert ev.emulation_residuals(rs, m) == pytest.approx((g1 - 0.003, g2 - 0.005))
    assert ev.emulation_residual(rs, m) == pytest.approx(abs(g1 - 0.003))
    r0 = (0.048 * 0.0531 - 0.03 * 0.05525) / 0.018
    assert ev.extrapolation_difference(rs) == pytest.approx(r0 - 0.05)
    assert ev.EMULATION_TOL == 0.0015
    g_star = ev.expected_g_meas("oxygen", 13, 0.1)  # g* at the measured loss fraction
    proj, a = gen.IONS["oxygen"][0], 16
    t = gen.slab_mm("oxygen", 13, 1) * 1.1
    exact = ev.slab_energy_loss(a * 2.4, a, t, ev.model_stopping(proj))
    b_s = ev.fine_algorithm_loss(proj, a, 2.4, t, "oxygen", 20) / exact - 1
    g_s = ev.algorithm_loss(proj, a, 2.4, t, "oxygen") / exact - 1
    assert g_star[0] == pytest.approx((1 + g_s) / (1 + b_s) - 1, rel=1e-12)


def test_planted_algorithm_losses_pass_the_gate() -> None:
    """TOPAS = the emulated algorithms with stopping (1 + r0) x ionmc's: |e1| < 1e-5 at every
    point; the opposite one-step bias fails the gate at the low-energy points."""
    fails = []
    for ion in ("carbon", "oxygen"):
        proj, _z, a, _p = gen.IONS[ion]
        for case in gen.ALL_CASES:
            e0, m = float(gen.CASES[case]), _model(ion, case)
            t = gen.slab_mm(ion, case, 1)
            de_m = ev.predicted_loss(f"model:{ion}", ev.model_stopping(proj), a, e0, t)
            for r0 in (0.05, -0.04):
                one = ev.algorithm_loss(proj, a, e0, t * (1 + r0), ion)
                fine = ev.fine_algorithm_loss(proj, a, e0, t * (1 + r0), ion, 20)
                rs = [{"r_M": one / de_m - 1}, {"r_M": 0.0}, {"r_M": fine / de_m - 1}]
                assert abs(ev.emulation_residuals(rs, m)[0]) < 1e-5, (ion, case, r0)
            g1 = m["g1"]
            anti = [{"r_M": 0.05 - g1}, {"r_M": 0.0}, {"r_M": 0.05}]
            if ev.emulation_residual(anti, m) > ev.EMULATION_TOL:
                fails.append((ion, case))
    assert ("carbon", 13) in fails and ("oxygen", 13) in fails and ("oxygen", 15) in fails
    assert ("carbon", 21) not in fails


def test_repeat_rule_unit() -> None:
    n = gen.bundle_name
    m = FIXED_MODEL
    parts = (0, 1, 6)
    reps = (8, 9, 14)

    def res(runs: dict[str, Any], case: int = 15, prt: tuple[int, ...] = parts) -> Any:
        return ev.resolve_point(runs, case, prt, _ratios, m if len(prt) == 3 else None)

    def triple(r1: float, r2: float, rf: float, ps: tuple[int, ...] = parts) -> dict[str, Any]:
        return {n(15, p): _fake(r=r) for p, r in zip(ps, (r1, r2, rf), strict=True)}

    ok = triple(0.003, 0.005, 0.0)  # e1 = e2 = 0
    pt = res(ok)
    assert pt.valid and pt.used == {} and pt.emulation_residual() == pytest.approx(0.0, abs=1e-12)
    assert pt.emulation_residual_f2() == pytest.approx(0.0, abs=1e-12)
    assert pt.extrapolation_difference() == pytest.approx(-1.0 / 3000.0, rel=1e-9)
    assert pt.thickness_difference() == pytest.approx(0.002)
    for extra in reps:  # a repeat without a recorded failure
        with pytest.raises(ev.RecordError, match="without a recorded failure"):
            res({**ok, n(15, extra): _fake()})
    # validity failure of one run: exactly that run is repeated
    bad = {**ok, n(15, 6): _fake(False, "deficit")}
    pt = res(bad)
    assert not pt.valid and pt.reason == "invalid run; no repeat"
    pt = res({**bad, n(15, 14): _fake(r=0.0)})
    assert pt.valid and pt.used == {n(15, 14): f"original {n(15, 6)} invalid: deficit"}
    with pytest.raises(ev.RecordError, match="is valid"):  # repeat of a valid run
        res({**bad, n(15, 14): _fake(r=0.0), n(15, 8): _fake(r=0.003)})
    pt = res({**bad, n(15, 14): _fake(False, "killed")})
    assert not pt.valid and "second failure" in str(pt.reason)
    pt = res({**bad, n(15, 14): _fake(r=0.004)})  # valid repeat, but the triple fails e1
    assert not pt.valid and str(pt.reason).startswith(ev.EMULATION_REASON)
    # gate failure (e1): the WHOLE triple is repeated and replaced
    emul = triple(0.0055, 0.005, 0.0)  # e1 = 0.0025
    pt = res(emul)
    assert not pt.valid and pt.reason == f"{ev.EMULATION_REASON}; no repeat"
    assert pt.emulation_residual() == pytest.approx(0.0025, rel=1e-6)
    for drop in reps:
        partial = {**emul, **{n(15, p): _fake(r=0.0) for p in reps if p != drop}}
        with pytest.raises(ev.RecordError, match="partial repeat set"):
            res(partial)
    good = {**emul, **triple(0.013, 0.015, 0.01, reps)}
    pt = res(good)
    assert pt.valid and set(pt.used) == {n(15, p) for p in reps}
    assert all(ev.EMULATION_REASON in v for v in pt.used.values())
    again = {**emul, **triple(0.0055, 0.005, 0.0, reps)}
    pt = res(again)
    assert not pt.valid and "second failure" in str(pt.reason)
    # e2 and the extrapolation difference are reported only: no failure
    pt = res(triple(0.003, 0.0125, 0.0))  # e2 = 0.0075, r_0 - r_fine = -0.0167
    assert pt.valid and pt.emulation_residual_f2() == pytest.approx(0.0075, rel=1e-9)
    assert abs(pt.extrapolation_difference() or 0.0) > 0.01
    # single-run points (He, fluctuations off)
    with pytest.raises(ev.RecordError, match="missing"):
        res({n(17, 12): _fake()}, 17, (4,))
    assert res({}, 17, (4,)).reason == "missing run"
    assert res({n(17, 4): _fake(False, "x"), n(17, 12): _fake(r=0.0)}, 17, (4,)).valid
    with pytest.raises(ev.RecordError, match="without a recorded failure"):
        res({n(17, 4): _fake(), n(17, 12): _fake()}, 17, (4,))
    with pytest.raises(ValueError, match="triple"):
        ev.resolve_point(ok, 15, parts, _ratios)
    with pytest.raises(ValueError, match="triple"):
        ev.resolve_point(ok, 15, (0, 1), _ratios, m)


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


def _fine_track(slab_mm: float, n_steps: int = 20) -> str:
    """A primary crossing the slab in ``n_steps`` equal limiter steps (fine-step smoke)."""
    lines = [
        "* G4Track Information:   Particle = C12,   Track ID = 1,   Parent ID = 0",
        "    0      0 fm      0 fm     -1 um   28.8 MeV   0 eV   0 fm   0 fm   World initStep",
        "    1      0 fm      0 fm      0 fm   28.8 MeV   0 eV   1 um   1 um   Water Transp",
    ]
    for k in range(1, n_steps + 1):
        z_um = slab_mm * 1000.0 * k / n_steps
        lines.append(
            f"    {k + 1}      0 fm      0 fm   {z_um:.6g} um   27 MeV   1 keV   1 nm   1 um   "
            + ("World Transportation" if k == n_steps else "Water StepLimiter")
        )
    lines.append(
        f"    {n_steps + 2}      0 fm      0 fm     40 cm   27 MeV   0 eV   1 cm   1 cm   Out T"
    )
    return "\n".join(lines) + "\n"


def _fine_stdout(n: int, slab_mm: float, *, coarse: int = -1) -> str:
    """Fine-step smoke stdout: every primary 20 steps in the slab, primary ``coarse`` 1 step."""
    return _stdout(n) + "".join(_fine_track(slab_mm, 1 if i == coarse else 20) for i in range(n))


def test_count_slab_steps_fine_steps() -> None:
    slab = 0.00208
    counts = ev.count_slab_steps(_fine_stdout(4, slab, coarse=2), slab, STEP_RE, TRACK_RE)
    assert [c["steps"] for c in counts.values()] == [20, 20, 1, 20]
    assert all(c["traverses"] for c in counts.values())


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


def _review(
    review_id: str,
    head: str,
    started: str,
    finished: str,
    *,
    seed: bytes | None = CW_SEED,
    verdict: str = "pass",
    **over: Any,
) -> dict[str, Any]:
    """A codex-worker-signed review record in the layout the review service produces
    (``review.build_record`` + ``attest.sign_record``, as ``inspect_codex_review`` returns it);
    ``over`` changes the payload before signing; ``seed`` None: unsigned."""
    report = {"head_sha": head, "verdict": verdict, "summary": "synthetic", "findings": []}
    payload = {
        "schema_version": 1,
        "review_id": review_id,
        "task_id": "V3-008",
        "head_sha": head,
        "base_sha": "0" * 40,
        "status": "passed",
        "started_at": started,
        "finished_at": finished,
        "model": "gpt-5.6-sol",
        "provider": "openai",
        "tool": "codex exec",
        "experiment_id": ev.review_records_module().current_experiment_id(),
        "report_sha256": hashlib.sha256(json.dumps(report).encode()).hexdigest(),
        "report": report,
    }
    payload.update(over)
    if seed is None:
        return payload
    at = IO.attest_module()
    return dict(at.sign_record(payload, at.keypair_from_seed(seed), "codex-worker"))


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
        for top in ("src/ionmc", "infrastructure")
        for p in (ROOT / top).rglob("*.py")
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
    cw = IO.attest_module().keypair_from_seed(CW_SEED)
    _write_json(
        repo / IO.DEFAULT_TRUST,
        {
            "reference-runner": [{"key_id": kp.key_id, "public_key": kp.public_key.hex()}],
            "codex-worker": [{"key_id": cw.key_id, "public_key": cw.public_key.hex()}],
        },
    )
    seal = repo / ev.ENV_SEAL
    seal.parent.mkdir(parents=True, exist_ok=True)
    seal.write_text(ev.env_seal.seal_text(ev.env_seal.current_environment()))
    (repo / ev.PLAN_PATH).parent.mkdir(parents=True, exist_ok=True)
    (repo / ev.PLAN_PATH).write_text("# synthetic acceptance plan\n")
    reviewed = _commit(repo, "sealed evaluator and bundles")
    pin_a = _review(RID_A, reviewed, PIN_A_STARTED, PIN_A_FINISHED)
    _write_json(repo / ev.REVIEW_PIN, pin_a)
    pin_a_sha = _commit(repo, "pin A and its signed validation record")
    _write_json(repo / ev.MARKERS, TEST_MARKERS)
    markers_sha = _commit(repo, "markers from the smoke logs")
    pin_b = _review(RID_B, markers_sha, PIN_B_STARTED, PIN_B_FINISHED)
    _write_json(repo / ev.MARKERS_PIN, pin_b)
    (repo / ev.PLAN_PATH).write_text(
        "# synthetic acceptance plan\n\n" + _note(RID_A, reviewed, RID_B, markers_sha)
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


_PLANTED: dict[tuple[str, int, int, str, float], float] = {}


def _planted_loss(info: dict[str, Any], de_m: float, r: float) -> float:
    """Synthetic TOPAS mean loss: Geant4's emulated algorithms (one-step for slab t / 1.6 t and
    He, 20-step linear branch for the fine-step runs) with stopping (1 + r) x ionmc's, i.e. on
    the slab scaled by (1 + r); smokes: the model loss x (1 + r)."""
    if info["kind"] == "smoke":
        return de_m * (1.0 + r)
    ion, case, a = info["ion"], int(info["case"]), int(info["a"])
    fine = info.get("max_step_um") is not None
    key = (ion, case, 1 if info["thickness_factor"] == 0.03 else 2, "fine" if fine else "one", r)
    if key not in _PLANTED:
        proj = gen.IONS[ion][0]
        t = float(info["slab_mm"]) * (1.0 + r)
        e0 = float(info["e0_per_u_mev"])
        if fine:
            _PLANTED[key] = ev.fine_algorithm_loss(proj, a, e0, t, ion, gen.FINE_DIVISOR)
        else:
            _PLANTED[key] = ev.algorithm_loss(proj, a, e0, t, ion)
    return _PLANTED[key]


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
        stdout
        if stdout is not None
        else _stdout(n, particle, bool(info["fluctuations"]))
        + (LIMITER_LINE if info.get("max_step_um") is not None else "")
    )
    (odir / "stderr.txt").write_text(stderr)
    for name in ("case.json", "input.txt"):
        shutil.copy(bdir / name, odir / "work" / name)
    de = _planted_loss(info, base["de_m"][bundle], r) if plant_g else base["de_m"][bundle] * (1 + r)
    h, rows, tot = _phase_space(doc, de, rng, noise=noise, deficit=deficit)
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
    coarse: bool = False,
    only: tuple[str, ...] = ev.ALL_SMOKES,
    over: dict[str, dict[str, Any]] | None = None,
) -> None:
    """The five smokes (``<label>-smokes``), run at the pin-A commit between the two reviews.
    Step-count smokes: ``split`` / ``stop``: that smoke's primary 7 takes two steps / stops in
    the slab; ``deficit``: that smoke's phase space misses a row; ``coarse``: primary 7 of the
    fine-step smoke crosses the slab in one step. ``over``: per-smoke keywords."""
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
        if name == ev.FINE_SMOKE:
            kw["stdout"] = _fine_stdout(n, float(info["slab_mm"]), coarse=7 if coarse else -1)
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
    if p in (0, 1, 5, 6) and c in (17, 18):
        return TABLE_FACTOR * 1.03 - 1.0
    return S


def _r_model(c: int, p: int) -> float:
    return 0.0


def _r_mixed(c: int, p: int) -> float:
    return S if p in (0, 1, 5, 6) else 0.0  # C sides with the table, O with the model


def _r_he_fail(c: int, p: int) -> float:
    return 0.02 if p == 4 else _r_icru(c, p)


def _r_fragile(c: int, p: int) -> float:
    """TOPAS 0.5 % above the table everywhere: I at the declared systematic, N when doubled."""
    return 0.0 if p == 4 else TABLE_FACTOR * 1.005 - 1.0


def _r_extrapolation(c: int, p: int) -> float:
    """O 20 MeV/u: the one-step pair off by +0.09 % (t) and -0.09 % (1.6 t) relative to the
    model's one-step bias: check (i) passes (<= 0.15 %), check (ii) fails (|r_0 - r_fine| ~
    0.37 %)."""
    if p == 4:
        return 0.0
    return S + {(19, 2): 0.0009, (19, 3): -0.0009}.get((c, p), 0.0)


def _r_inconsistent(c: int, p: int) -> float:
    if p == 4:
        return 0.0
    return S + (0.01 if (c, p) == (19, 2) else 0.0)  # O 20 MeV/u slab t off by ~1 %: e1 fails


BAD = {"deficit": 1}


def _both_lines(n: int, fluct: bool) -> str:
    """A stdout showing BOTH fluctuation dump lines (a marker set that cannot discriminate)."""
    return _stdout(n, fluct=fluct) + f"Fluctuations of dE/dx are enabled   {int(not fluct)}\n"


RESULT_LABELS = (
    "icru", "model", "mixed", "he_fail", "extrapolation", "invalid", "second_invalid",
    "repeat_ok", "split", "stopped", "smoke_deficit", "off_hits_on", "on_hits_off", "fragile",
    "coarse",
)  # fmt: skip


@pytest.fixture(scope="module")
def results(base: dict[str, Any]) -> dict[str, dict[str, Any]]:
    c15, o19_2t = gen.bundle_name(15, 0), gen.bundle_name(19, 3)
    _scenario(base, "icru", _r_icru)
    _scenario(base, "model", _r_model)
    _scenario(base, "mixed", _r_mixed)
    _scenario(base, "he_fail", _r_he_fail)
    _scenario(base, "extrapolation", _r_extrapolation)
    _scenario(base, "invalid", _r_icru, run_kw={c15: BAD})
    _scenario(  # C 5 MeV/u invalid + an invalid repeat; O 20 MeV/u inconsistent, no repeat
        base,
        "second_invalid",
        _r_inconsistent,
        run_kw={c15: BAD, gen.bundle_name(15, 8): BAD},
        repeats={gen.bundle_name(15, 8): S},
    )
    _scenario(  # both restored: C part 8 for the invalid run; O 20 MeV/u the whole triple
        base,
        "repeat_ok",
        _r_inconsistent,
        run_kw={c15: BAD},
        repeats={gen.bundle_name(15, 8): S, **{gen.bundle_name(19, p): S for p in (10, 11, 15)}},
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
    _scenario(base, "coarse", _r_icru, smoke_kw={"coarse": True})
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
            assert v["sigma_b"] == pytest.approx(np.hypot(v["se_b"], 0.002))
            assert 3e-5 < v["se_b"] < 3e-4  # 0.4 % row noise, 1000 rows, fine-step points
            assert v["max_abs_consistency_residual"] < 2e-3  # planted g: Delta_meas ~ Delta_pred
            assert v["max_abs_emulation_residual"] < 1e-3  # planted g: g_meas ~ g_model
            assert v["max_abs_extrapolation_difference_report_only"] < 1.5e-3
            assert v["max_abs_emulation_residual_f2_report_only"] < 1e-3
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
    assert doc["fine_step_premise"] == {
        "primaries": 50,
        "min_steps_in_slab": 20,
        "traversing_primaries": 50,
        "phase_space_complete": True,
        "pass": True,
    }
    assert doc["one_step_premise"] == {
        "checked": True,
        "smokes": dict.fromkeys(ev.STEP_SMOKES, one),
        "pass": True,
    }
    prov = doc["provenance"]
    runs = prov["runs"]
    assert len(runs) == 58 and set(runs) == {gen.bundle_name(c, p) for c, p in gen.evidence_keys()}
    assert all(
        v["scored"] == SYN_HISTORIES and v["run_id"].startswith("REF-") for v in runs.values()
    )
    assert set(prov["smokes"]) == set(ev.ALL_SMOKES)
    assert {v["code_sha"] for v in prov["smokes"].values()} == {base["pin_a_sha"]}
    pins = prov["review_pins"]
    assert pins["evaluator"]["review_id"] == RID_A
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
    assert seal["stdlib_modules_in_seal"] > 500 and seal["stdlib_modules_verified"] > 20
    assert seal["sys_prefix"] == sys.prefix and seal["sys_base_prefix"] == sys.base_prefix
    assert env["runtime"]["stdlib_modules_verified"] > 20
    assert (
        env["runtime_after"]["stdlib_modules_verified"] >= env["runtime"]["stdlib_modules_verified"]
    )
    assert "validation/scripts/stopping/sealed_launch.py" in sealed
    cw_key = IO.attest_module().keypair_from_seed(CW_SEED).key_id
    for name in ("evaluator", "markers"):
        att = pins[name]["attestation"]
        assert att["signer"] == "codex-worker" and att["key_id"] == cw_key
        assert len(att["signed_at_utc"]) > 10
    assert doc["amendment"] == f"I3-FU-1 step 2 (revision {ev.AMENDMENT_REVISION})"
    assert f"I3-FU-1 step 2 (revision {ev.AMENDMENT_REVISION})" in (ev.__doc__ or "")
    assert prov["plan_revision_records"] == 0  # the synthetic plan holds no revision record
    conv = doc["convention_sensitivity"]
    assert conv["invariant"] is True and conv["pattern_alt"] == "ICRU"
    assert conv["classes_alt"] == {ion: _cls(doc, ion) for ion in ("carbon", "oxygen")}
    for ion in ("carbon", "oxygen"):
        for v in conv["max_abs_delta_r_i_per_band"][ion].values():
            assert v is not None and round(v, 3) == v and v <= 0.001
    assert "E/A, A = mass number" in doc["method"]["energy_convention"]
    runtime = prov["environment"]["runtime"]
    assert runtime["launch"]["launch_ok"] is True
    assert (
        prov["environment"]["runtime_after"]["numpy_modules_verified"]
        >= runtime["numpy_modules_verified"]
    )


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


def test_extrapolation_difference_is_reported_only(results: dict[str, dict[str, Any]]) -> None:
    """Reviewer pass 4: r_0 - r_fine (~0.37 % here) and e2 are reported, not gates (they carry
    the slope of the stopping ratio); the point stays valid and the result conclusive."""
    doc = results["extrapolation"]
    assert doc["pattern"] == "ICRU" and doc["conclusive"] and doc["invalid_points"] == {}
    o = doc["bands"]["oxygen"]["[20,50]"]
    assert o["n"] == 3 and o["max_abs_emulation_residual"] <= ev.EMULATION_TOL
    assert o["max_abs_extrapolation_difference_report_only"] > 0.002


def test_points_are_the_fine_step_runs(results: dict[str, dict[str, Any]]) -> None:
    """The band means are those of the fine-step runs (no one-step bias, no extrapolation): in
    the ICRU scenario r_I of the fine runs is 0 up to noise, also where |R_model| ~ 3e-4."""
    doc = results["icru"]
    assert abs(doc["bands"]["oxygen"]["[2,5)"]["mean_r_I"]) <= 0.001
    assert doc["method"]["method_systematic"] == 0.002
    assert doc["method"]["method_systematic_robustness"] == 0.004
    assert doc["method"]["emulation_tolerance"] == 0.0015
    assert "not gates" in doc["method"]["report_only"]


def test_fine_step_premise_failure_is_inconclusive(results: dict[str, dict[str, Any]]) -> None:
    doc = results["coarse"]
    assert doc["pattern"] == "INCONCLUSIVE" and doc["inconclusive_reasons"] == ["fine_step_premise"]
    fine = doc["fine_step_premise"]
    assert fine["pass"] is False and fine["min_steps_in_slab"] == 1
    assert fine["traversing_primaries"] == 49 and fine["primaries"] == 50


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
        "oxygen case 19": f"{ev.EMULATION_REASON}; no repeat",
    }
    assert list(doc["repeats_used"]) == [gen.bundle_name(15, 8)]
    o = doc["bands"]["oxygen"]["[20,50]"]  # 20 MeV/u
    assert o["n"] == 2 and o["max_abs_consistency_residual"] == pytest.approx(0.01, abs=2e-3)
    assert o["max_abs_emulation_residual"] == pytest.approx(0.01, abs=2e-3)


def test_valid_repeats_restore_conclusiveness(results: dict[str, dict[str, Any]]) -> None:
    doc = results["repeat_ok"]
    assert doc["conclusive"] and doc["pattern"] == "ICRU"
    assert doc["invalid_points"] == {} and list(doc["invalid_runs"]) == [gen.bundle_name(15, 0)]
    cons = ev.EMULATION_REASON
    assert doc["repeats_used"] == {
        gen.bundle_name(15, 8): f"original {gen.bundle_name(15, 0)} invalid: "
        + doc["invalid_runs"][gen.bundle_name(15, 0)],
        gen.bundle_name(19, 10): f"original {gen.bundle_name(19, 2)} invalid: {cons}",
        gen.bundle_name(19, 11): f"original {gen.bundle_name(19, 3)} invalid: {cons}",
        gen.bundle_name(19, 15): f"original {gen.bundle_name(19, 7)} invalid: {cons}",
    }
    assert (
        doc["bands"]["carbon"]["[5,10)"]["n"] == 2 and doc["bands"]["oxygen"]["[20,50]"]["n"] == 3
    )
    assert doc["bands"]["oxygen"]["[20,50]"]["max_abs_consistency_residual"] < 2e-3
    assert len(doc["provenance"]["runs"]) == 62


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
        "max_abs_emulation_residual", "max_abs_emulation_residual_f2_report_only",
        "max_abs_extrapolation_difference_report_only",
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
        assert set(doc["fine_step_premise"]) == {
            "primaries", "min_steps_in_slab", "traversing_primaries", "phase_space_complete",
            "pass",
        }  # fmt: skip
        for reason in [*doc["invalid_runs"].values(), *doc["invalid_points"].values()]:
            rest = reason.replace(ev.EMULATION_REASON, "")
            assert not re.search(r"\d\.\d{3,}", rest), reason  # reasons carry no values
        text = json.dumps(doc)
        for bad in ('delta_e"', "de_m", "de_i", "exit_energy", "differences", '"g1"', "r_i_alt"):
            assert bad not in text.lower()
        conv = doc["convention_sensitivity"]
        assert set(conv) == {
            "max_abs_delta_r_i_per_band",
            "pattern_alt",
            "classes_alt",
            "invariant",
        }
        for ion in conv["max_abs_delta_r_i_per_band"].values():
            assert all(v is None or round(v, 3) == v for v in ion.values())


def test_classes_use_full_precision_values() -> None:
    """A band whose rounded mean would flip the class keeps the full-precision class."""
    pts = [
        {"e0": 2.4, "valid": True, "r_M": 0.0, "r_I": 0.0, "s": 0.0304, "se_rel": 0.0,
         "d_thick": None, "d_resid": None, "d_pred": 0.0, "r_model": 0.0, "d_emul": None,
         "d_extra": None, "d_emul2": None},
    ]  # fmt: skip
    band = ev.band_aggregates(pts)["[2,5)"]
    # X = 0.0304 / 3 = 0.010133 -> capped 0.01; k = 1.645 x 0.002 = 0.00329; |r_I| + k <= X: I
    assert band["class"] == "I" and band["mean_s"] == 0.03 and band["x_b"] == 0.01
    pts[0]["r_I"] = 0.0066  # full: 0.00989 <= 0.01 (I); rounded 0.007: 0.01029 > 0.01 (N)
    band = ev.band_aggregates(pts)["[2,5)"]
    assert band["mean_r_I"] == 0.007 and band["class"] == "I"
    assert ev.band_classes(pts, 0.002)["[2,5)"] == "I"
    assert ev.band_classes(pts, 0.004)["[2,5)"] == "M"  # the robustness re-classification
    assert ev.band_classes([], 0.002) == dict.fromkeys(("[2,5)", "[5,10)", "[10,20)", "[20,50]"),
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
           note: str | None = None) -> dict[str, Any]:  # fmt: skip
    """A copy of the base repository whose pins / plan note are rewritten in ONE new commit,
    which then serves as code_sha."""
    w = _world(base, tmp)
    for rel, doc in ((ev.REVIEW_PIN, pin_a), (ev.MARKERS_PIN, pin_b)):
        if doc is not None:
            _write_json(w["repo"] / rel, doc)
    if note is not None:
        (w["repo"] / ev.PLAN_PATH).write_text("# synthetic acceptance plan\n\n" + note)
    w["code_sha"] = _commit(w["repo"], "rewritten pins")
    return w


def _tampered(doc: dict[str, Any]) -> dict[str, Any]:
    return {**doc, "finished_at": "2026-11-16T10:00:00+00:00"}  # changed after signing


@pytest.mark.parametrize("which", ["a", "b"])
@pytest.mark.parametrize(
    "kind, match",
    [
        ("untrusted", "codex-worker review record unverified"),
        ("unsigned", "codex-worker review record unverified"),
        ("tampered", "codex-worker review record invalid"),
        ("changes_required", "status 'changes_required' is not 'passed'"),
        ("verdict_mismatch", "codex-worker review record invalid: status 'passed' inconsistent"),
        ("wrong_task", "signed identity mismatch: task_id"),
        ("wrong_experiment", "signed identity mismatch: experiment_id"),
    ],
)
def test_review_pin_must_be_a_verified_passed_codex_worker_record(
    base: dict[str, Any], tmp_path: Path, which: str, kind: str, match: str
) -> None:
    """Rev13: the pins are codex-worker-signed review records verified with
    ``review_records.verify_review_record`` against the committed trust map; anything but
    ``verified`` (untrusted key, unsigned, tampered, inconsistent, foreign identity) and any
    non-passed review fails closed."""
    rid, head, st, fin = (
        (RID_A, base["reviewed"], PIN_A_STARTED, PIN_A_FINISHED)
        if which == "a"
        else (RID_B, base["markers_sha"], PIN_B_STARTED, PIN_B_FINISHED)
    )
    doc = {
        "untrusted": lambda: _review(rid, head, st, fin, seed=ROGUE_SEED),
        "unsigned": lambda: _review(rid, head, st, fin, seed=None),
        "tampered": lambda: _tampered(_review(rid, head, st, fin)),
        "changes_required": lambda: _review(
            rid, head, st, fin, verdict="changes_required", status="changes_required"
        ),
        "verdict_mismatch": lambda: _review(rid, head, st, fin, verdict="changes_required"),
        "wrong_task": lambda: _review(rid, head, st, fin, task_id="V3-005"),
        "wrong_experiment": lambda: _review(rid, head, st, fin, experiment_id="experiment-v2"),
    }[kind]()
    w = _repin(base, tmp_path, pin_a=doc) if which == "a" else _repin(base, tmp_path, pin_b=doc)
    with pytest.raises(ev.RecordError, match=match):
        _ctx(base, w, code_sha=w["code_sha"])


def test_review_pin_chronology_uses_the_signed_times(base: dict[str, Any], tmp_path: Path) -> None:
    """Runs starting before a pin's signed finished_at are refused (the times are inside the
    signed payload, so they cannot be moved without breaking the signature)."""
    w = _world(base, tmp_path, started="2026-11-16T12:29:59+00:00")  # before pin B finished
    with pytest.raises(ev.RecordError, match="not after both review pins"):
        _check(base, w)
    pin = json.loads((w["repo"] / ev.MARKERS_PIN).read_text())
    assert pin["finished_at"] == PIN_B_FINISHED and "attestation" in pin


def test_pins_must_be_named_by_the_plan_note(base: dict[str, Any], tmp_path: Path) -> None:
    ok = (RID_A, base["reviewed"], RID_B, base["markers_sha"])
    for i, (bad, match) in enumerate(
        (
            (("REVIEW-" + "c" * 32, *ok[1:]), "not the evaluator review named"),
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
        "# synthetic acceptance plan\n\n" + _note(RID_A, ok[1], RID_B, ok[1])
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
    fine = out[ev.FINE_SMOKE]
    assert fine["ok"] and fine["min_steps_in_slab"] == 20 and fine["traversing_primaries"] == 50


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
    with pytest.raises(ev.RecordError, match="lacks .*all five smokes are required"):
        _evaluate(w, "three", repo)


# --- review pins bound to signed validation records (Codex REVIEW-c4f3fdd1 finding 1) --------


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


def test_runtime_checks_fail_closed(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Revision 5 (Codex REVIEW-183057ce finding 3): a numpy module loaded from outside the
    verified distribution, and an existing bytecode cache file of a sealed module, refuse."""
    out = ev.runtime_checks()
    assert out["numpy_modules_verified"] > 10 and out["dont_write_bytecode"] is True
    shadow = ModuleType("numpy._shadow_for_test")
    shadow.__file__ = str(tmp_path / "shadow.py")
    monkeypatch.setitem(sys.modules, "numpy._shadow_for_test", shadow)
    with pytest.raises(ev.RecordError, match="numpy._shadow_for_test is loaded from outside"):
        ev.runtime_checks()
    monkeypatch.delitem(sys.modules, "numpy._shadow_for_test")
    real = ev._REAL_BYTECODE_CHECK
    cached = tmp_path / "gen.cpython-312.pyc"
    cached.write_bytes(b"stale")
    monkeypatch.setattr(ev.gen, "__cached__", str(cached))
    monkeypatch.setattr(ev.env_seal, "bytecode_check", real)
    monkeypatch.setattr(ev.env_seal, "numpy_modules", lambda: {})
    monkeypatch.setattr(ev, "closure_modules", lambda: {"gen_i3fu2_bundles": ev.gen})
    with pytest.raises(ev.RecordError, match="bytecode cache files exist for .*gen_i3fu2_bundles"):
        ev.runtime_checks()


# --- rev8: energy convention, launch, post-computation recheck, limiter marker -------------


def _conv_point(e0: float, r_i: float, r_i_alt: float) -> dict[str, Any]:
    return {"e0": e0, "valid": True, "r_M": 0.06, "r_I": r_i, "s": 0.05, "se_rel": 0.0,
            "r_I_alt": r_i_alt, "s_alt": 0.05}  # fmt: skip


def test_convention_sensitivity_flips_at_a_class_boundary() -> None:
    """X = 0.01, k = 1.645 x 0.002: I iff |r_I| <= 0.00671; r_I 0.0067 (I) vs the alternative
    reading 0.0070 (N) flips the pattern ICRU -> MIXED: not invariant."""
    pts = {
        ion: [_conv_point(e, 0.0, 0.0) for e in (2.4, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0)]
        for ion in ("carbon", "oxygen")
    }
    out = ev.convention_sensitivity(pts, "ICRU")
    assert out["invariant"] and out["pattern_alt"] == "ICRU"
    pts["carbon"][0] = _conv_point(2.4, 0.0067, 0.0070)
    pts["carbon"][1] = _conv_point(3.0, 0.0067, 0.0070)
    out = ev.convention_sensitivity(pts, "ICRU")
    assert not out["invariant"] and out["pattern_alt"] == "MIXED"
    assert out["classes_alt"]["carbon"]["[2,5)"] == "N"
    assert out["max_abs_delta_r_i_per_band"]["carbon"]["[2,5)"] == 0.0  # 0.03 pp -> 0.0 at 0.1 pp
    assert ev.mass_convention_factor("carbon") == pytest.approx(1.000274, abs=1e-6)
    assert ev.mass_convention_factor("oxygen") == pytest.approx(1.000592, abs=1e-6)


def test_launch_record_is_enforced(
    base: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    w = _world(base, tmp_path)
    record = VALID_LAUNCH["launcher_record"]
    assert isinstance(record, dict)
    for bad, match in (
        ({"dont_write_bytecode_flag": 0}, "python -B"),
        ({"pycache_prefix": None}, "pycache_prefix=<fresh empty directory>"),
        ({"pycache_prefix_exists": False}, "does not exist"),
        ({"pycache_prefix_files_at_start": 3}, "was not empty"),
        ({"isolated_flag": 0}, "python -I -S"),
        ({"no_site_flag": 0}, "python -I -S"),
        ({"launcher_record": None}, "not launched through sealed_launch.py"),
        ({"machinery_problems": ["foreign finder"]}, "import machinery"),
        ({"sys_path_at_capture": ["/elsewhere"]}, "explicit list"),
        ({"launcher_record": {**record, "launcher_sha256": "0" * 64}}, "not the sealed"),
        ({"launcher_record": {**record, "target_sha256": "0" * 64}}, "differs from this file"),
        ({"launcher_record": {**record, "target_path": "/elsewhere/x.py"}}, "the launcher ran"),
        (
            {"launcher_record": {**record, "sys_modules_at_launch": ["sys", "numpy.linalg"]}},
            "modules imported before the launch",
        ),
    ):
        monkeypatch.setattr(ev, "LAUNCH", {**VALID_LAUNCH, **bad})
        with pytest.raises(ev.RecordError, match=match):
            _ctx(base, w)


def test_post_computation_recheck_covers_lazy_numpy_imports(
    base: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A numpy module imported lazily during the computation (simulated: a shadow module from
    outside the verified distribution appears in sys.modules) is caught by the recheck."""
    repo = tmp_path / "repo"
    shutil.copytree(base["repo"], repo)
    w = {**base, "root": tmp_path, "repo": repo}
    _scenario(w, "lazy", _r_icru)
    _commit(repo, "records")
    real = ev.stepping_model_check
    shadow = ModuleType("numpy._lazy_shadow")
    shadow.__file__ = str(tmp_path / "lazy_shadow.py")

    def lazy() -> Any:
        sys.modules["numpy._lazy_shadow"] = shadow
        return real()

    monkeypatch.setattr(ev, "stepping_model_check", lazy)
    try:
        with pytest.raises(ev.RecordError, match="numpy._lazy_shadow is loaded from outside"):
            _evaluate(w, "lazy", repo)
    finally:
        sys.modules.pop("numpy._lazy_shadow", None)


def test_fine_step_run_requires_the_step_limiter_marker(
    base: dict[str, Any], tmp_path: Path
) -> None:
    w = _world(base, tmp_path, gen.bundle_name(13, 6), stdout=_stdout(SYN_HISTORIES))
    with pytest.raises(ev.RecordError, match="step limiter marker is missing"):
        _check(base, w)
    w2 = _world(base, tmp_path / "ok", gen.bundle_name(13, 6))
    assert _check(base, w2)["valid"]


def test_extreme_r_fine_is_an_emulation_domain_failure() -> None:
    """rev9: g* is not computable for r_fine ~ 3 (the slab scaled by 4 leaves the emulation
    domain): the point is INVALID with the reason ``emulation_domain`` (no uncaught ValueError);
    the reported residuals are None."""
    n = gen.bundle_name
    runs = {n(13, p): _fake(r=r) for p, r in zip((0, 1, 6), (3.01, 3.02, 3.0), strict=True)}
    pt = ev.resolve_point(runs, 13, (0, 1, 6), _ratios, _model("carbon", 13))
    assert not pt.valid and pt.reason == f"{ev.EMULATION_DOMAIN_REASON}; no repeat"
    assert ev.EMULATION_DOMAIN_REASON == "emulation_domain"
    assert pt.emulation_residual() is None and pt.emulation_residual_f2() is None


def test_amendment_revision_matches_the_plan() -> None:
    """The evaluator's AMENDMENT_REVISION equals the highest "**Revision N of amendment
    I3-FU-1" record of the acceptance plan: a plan revision without a matching constant fails."""
    plan = (ROOT / "validation" / "plans" / "v3-008-acceptance.md").read_text(encoding="utf-8")
    found = [int(n) for n in re.findall(r"^\*\*Revision (\d+) of amendment I3-FU-1\b", plan, re.M)]
    assert found and max(found) == ev.AMENDMENT_REVISION == ev.plan_revision_records(plan)
    assert ev.plan_revision_records("no records") == 0


def test_lazily_imported_stdlib_module_must_be_in_the_seal(
    base: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A stdlib module loaded during the computation but missing from (or differing from) the
    environment seal stops the evaluation at the post-computation recheck."""
    w = _world(base, tmp_path)
    ctx = _ctx(base, w)
    import wave  # noqa: F401 - any stdlib module; present in the seal

    seal = dict(ctx.stdlib_seal)
    rel = "wave.py"
    assert rel in seal
    seal.pop(rel)
    with pytest.raises(
        ev.RecordError, match="stdlib module wave.py is not in the environment seal"
    ):
        ev.runtime_checks(seal)
    with pytest.raises(ev.RecordError, match="stdlib module wave.py differs"):
        ev.runtime_checks({**ctx.stdlib_seal, rel: "0" * 64})
    assert ev.runtime_checks(ctx.stdlib_seal)["stdlib_modules_verified"] > 20
