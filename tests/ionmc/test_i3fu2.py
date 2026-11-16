"""Amendment I3-FU-1 step 2 (revision 1): bundle generator, frozen manifest and evaluator.

Bundles ``validation/reference_cases/topas/i3fu2-*`` and their frozen manifest come from
``validation/scripts/stopping/gen_i3fu2_bundles.py``; the evaluator
``validation/scripts/stopping/i3fu2_evaluate.py`` is exercised in temporary git repositories
holding a synthetic bundle set (1000 histories), a committed review pin and Ed25519-signed
synthetic reference-runner records (test key from a fixed seed, committed as the trusted key of
the temporary repository, as ``test_v5ss_evaluate.py`` does). Synthetic TOPAS outputs (phase space
with the ``ss_evaluate`` header layout, a TOPAS CSV ``edep.csv``) carry energy losses equal to the
ionmc prediction times planted factors; the synthetic libdEdx archive's ICRU73_NEW carbon/oxygen
water tables are the ionmc model times a known factor. No real libdEdx, ICRU or TOPAS data are
used.
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
ROW_NOISE = 0.01
TABLE_FACTOR = 1.05  # synthetic ICRU73_NEW = ionmc model x this factor (s ~ +5 %)
KEY_SEED = bytes(range(32))
PIN_FINISHED = "2026-11-16T12:00:00+00:00"
STARTED = "2026-11-16T13:00:00+00:00"
FLUCT_OFF_TEST_MARKER = r"^Fluctuations of dE/dx are enabled\s+0\s*$"


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
    assert sorted(gen.repeat_keys()) == sorted((c, p + 6) for c, p in evidence)
    assert gen.REPEAT_PARTS == (6, 7, 8, 9, 10, 11)
    assert gen.RESERVED_PARTS == tuple(range(12, 20))
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
        kind = "repeat" if part >= 6 else "evidence"
        assert info["kind"] == kind
        kinds[kind] += 1
        base = part - 6 if part >= 6 else part
        assert info["thickness_factor"] == (0.048 if base in (1, 3) else 0.03)
        if part >= 6:  # a repeat is byte-identical to its original except seed and wording
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
    assert len(set(used)) == 80 and min(used) == 20619013 and max(used) == 20620121
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
    fluct_off = {gen.bundle_name(c, p) for c in (13, 21) for p in (5, 11)}
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
    assert ev.band_of(2.4) == "[2,5)" and ev.band_of(5.0) == "[5,10)"
    assert ev.band_of(20.0) == "[20,50]" and ev.band_of(50.0) == "[20,50]"
    with pytest.raises(ValueError):
        ev.band_of(50.1)


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


def _fake(valid: bool = True, reason: str | None = None, r: float = 0.0) -> dict[str, Any]:
    return {"valid": valid, "reason": reason, "r": r}


def _ratios(run: dict[str, Any]) -> dict[str, Any]:
    return {"r_M": run["r"]}


def test_repeat_rule_unit() -> None:
    n = gen.bundle_name
    ok = {n(15, 0): _fake(r=0.01), n(15, 1): _fake(r=0.012)}
    pt = ev.resolve_point(ok, 15, (0, 1), _ratios)
    assert pt.valid and pt.used == {} and pt.thickness_difference() == pytest.approx(0.002)
    with pytest.raises(ev.RecordError, match="without an invalid original"):
        ev.resolve_point({**ok, n(15, 6): _fake()}, 15, (0, 1), _ratios)
    bad_t = {n(15, 0): _fake(False, "deficit"), n(15, 1): _fake(r=0.012)}
    pt = ev.resolve_point(bad_t, 15, (0, 1), _ratios)
    assert not pt.valid and "no valid repeat" in str(pt.reason)
    pt = ev.resolve_point({**bad_t, n(15, 6): _fake(r=0.011)}, 15, (0, 1), _ratios)
    assert pt.valid and pt.used == {n(15, 6): f"original {n(15, 0)} invalid: deficit"}
    with pytest.raises(ev.RecordError, match="is valid"):  # repeat of the valid 1.6t run
        ev.resolve_point({**bad_t, n(15, 7): _fake()}, 15, (0, 1), _ratios)
    pt = ev.resolve_point({**bad_t, n(15, 6): _fake(False, "killed")}, 15, (0, 1), _ratios)
    assert not pt.valid and "second invalid" in str(pt.reason)
    incons = {n(15, 0): _fake(r=0.0), n(15, 1): _fake(r=0.01)}
    pt = ev.resolve_point(incons, 15, (0, 1), _ratios)
    assert not pt.valid and pt.thickness_difference() == pytest.approx(0.01)
    assert str(pt.reason).startswith(ev.CONSISTENCY_REASON)
    for rep in (n(15, 6), n(15, 7)):  # either slab may be repeated after a consistency failure
        r = 0.001 if rep.endswith("part6") else 0.0
        pt = ev.resolve_point({**incons, rep: _fake(r=r)}, 15, (0, 1), _ratios)
        assert not pt.valid if rep.endswith("part6") else pt.valid
    pt = ev.resolve_point({**incons, n(15, 7): _fake(r=0.02)}, 15, (0, 1), _ratios)
    assert not pt.valid and "second invalid" in str(pt.reason)
    with pytest.raises(ev.RecordError, match="missing"):
        ev.resolve_point({n(17, 10): _fake()}, 17, (4,), _ratios)
    assert ev.resolve_point({}, 17, (4,), _ratios).reason == "missing run"


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


@pytest.fixture(scope="module")
def base(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """Temporary repository: synthetic bundles (1000 histories) + manifest + trust (commit
    "reviewed"), then the review pin (commit "code_sha"); archive and model predictions."""
    root = tmp_path_factory.mktemp("i3fu2")
    repo = root / "repo"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "main")
    bundles = repo / gen.REL_PREFIX
    gen.generate(bundles, histories=SYN_HISTORIES)
    (repo / ev.MANIFEST).write_text(gen.manifest_text(bundles))
    kp = IO.attest_module().keypair_from_seed(KEY_SEED)
    trust = repo / IO.DEFAULT_TRUST
    trust.parent.mkdir(parents=True)
    trust.write_text(
        json.dumps({"reference-runner": [{"key_id": kp.key_id, "public_key": kp.public_key.hex()}]})
    )
    reviewed = _commit(repo, "bundles")
    pin = repo / ev.REVIEW_PIN
    pin.parent.mkdir(parents=True)
    pin.write_text(
        json.dumps({"review_id": "REVIEW-syn", "head_sha": reviewed, "finished_at": PIN_FINISHED})
        + "\n"
    )
    code_sha = _commit(repo, "review pin")
    archive, sha, members = _archive(root, TABLE_FACTOR)
    cases, de_m = {}, {}
    for c, p in gen.evidence_keys() + gen.repeat_keys():
        name = gen.bundle_name(c, p)
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
        "reviewed": reviewed,
        "code_sha": code_sha,
        "archive": archive,
        "sha": sha,
        "members": members,
        "cases": cases,
        "de_m": de_m,
    }


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
    code_sha: str | None = None,
    succeeded: bool = True,
    exit_code: int = 0,
    tamper_record: bool = False,
) -> str:
    """Materialized outputs under ``<root>/out-<label>/<run_id>/`` and the signed record +
    listing under ``<repo>/records/<label>/<run_id>/`` (not committed here)."""
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
    odir = base["root"] / f"out-{label}" / run_id
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
    h, rows, tot = _phase_space(
        doc, base["de_m"][bundle] * (1.0 + r), rng, noise=noise, deficit=deficit
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
            "finished_at": "2026-11-16T14:00:00+00:00",
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
    (repo / "records" / label / "index.json").write_text(json.dumps({"runs": runs}, indent=1))


def _scenario(
    base: dict[str, Any],
    label: str,
    r_of: Callable[[int, int], float],
    *,
    noise: float = ROW_NOISE,
    run_kw: dict[str, dict[str, Any]] | None = None,
    repeats: dict[str, float] | None = None,
) -> str:
    """All 40 evidence runs (r_of(case, part) planted relative to the model slab loss), plus
    the ``repeats`` {repeat bundle: r}; ``run_kw`` per bundle (e.g. a deficit)."""
    repo, entries = base["repo"], {}
    run_kw = run_kw or {}
    keys = ev.run_bundles()
    todo = [(gen.bundle_name(c, p), r_of(c, p)) for c, p in gen.evidence_keys()]
    todo += sorted((repeats or {}).items())
    for b, r in todo:
        c, p = keys[b]
        rng = np.random.default_rng([20261116, c, p])  # per run: identical across scenarios
        entries[b] = _write_run(base, repo, label, b, r, rng, noise=noise, **run_kw.get(b, {}))
    _write_index(repo, label, entries)
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


R0, B_LIN = 0.05, 0.1  # Richardson scenario: r = R0 + B_LIN f exactly (noise free)


def _r_linear(c: int, p: int) -> float:
    if p == 4:
        return 0.0
    return R0 + B_LIN * (0.048 if p in (1, 3) else 0.03)


def _r_inconsistent(c: int, p: int) -> float:
    if p == 4:
        return 0.0
    return S + (0.01 if (c, p) == (19, 3) else 0.0)  # O 20 MeV/u: |r(f2) - r(f1)| = 0.01


BAD = {"deficit": 1}


@pytest.fixture(scope="module")
def results(base: dict[str, Any]) -> dict[str, dict[str, Any]]:
    c15, o19_2t = gen.bundle_name(15, 0), gen.bundle_name(19, 3)
    _scenario(base, "icru", _r_icru)
    _scenario(base, "model", _r_model)
    _scenario(base, "mixed", _r_mixed)
    _scenario(base, "he_fail", _r_he_fail)
    _scenario(base, "linear", _r_linear, noise=0.0)
    _scenario(base, "invalid", _r_icru, run_kw={c15: BAD})
    _scenario(  # C 5 MeV/u invalid + an invalid repeat; O 20 MeV/u inconsistent, no repeat
        base,
        "second_invalid",
        _r_inconsistent,
        run_kw={c15: BAD, gen.bundle_name(15, 6): BAD},
        repeats={gen.bundle_name(15, 6): S},
    )
    _scenario(  # both restored by valid repeats (C part 6; O 1.6t part 9)
        base,
        "repeat_ok",
        _r_inconsistent,
        run_kw={c15: BAD},
        repeats={gen.bundle_name(15, 6): S, gen.bundle_name(19, 9): S},
    )
    _commit(base["repo"], "records")
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(ev, "FLUCT_OFF_MARKER", FLUCT_OFF_TEST_MARKER)
        labels = ("icru", "model", "mixed", "he_fail", "linear", "invalid")
        out = {lab: _evaluate(base, lab) for lab in (*labels, "second_invalid", "repeat_ok")}
    assert o19_2t in out["second_invalid"]["provenance"]["runs"]
    return out


def _cls(doc: dict[str, Any], ion: str) -> dict[str, str]:
    return {b: v["class"] for b, v in doc["bands"][ion].items()}


def test_icru_scenario_with_h_band(results: dict[str, dict[str, Any]]) -> None:
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
            assert 1e-4 < v["se_b"] < 1.5e-3  # 1 % row noise, 1000 rows, extrapolated points
            assert v["max_abs_thickness_difference"] < 3e-3
            if (ion, b) != ("carbon", "[10,20)"):
                assert v["mean_r_I"] == pytest.approx(0.0, abs=4e-3)
    assert doc["bands"]["carbon"]["[10,20)"]["mean_r_I"] == pytest.approx(0.03, abs=4e-3)
    c = doc["controls"]
    assert c == {
        "helium": {"pass": True, "max_abs_r_m": c["helium"]["max_abs_r_m"]},
        "fluctuations_off": {
            "pass": True,
            "max_abs_difference": c["fluctuations_off"]["max_abs_difference"],
        },
        "controls_pass": True,
    }
    assert c["helium"]["max_abs_r_m"] < 2e-3 and c["fluctuations_off"]["max_abs_difference"] < 2e-3
    assert doc["invalid_runs"] == {} and doc["invalid_points"] == {} and doc["repeats_used"] == {}
    assert doc["inconclusive_reasons"] == []
    runs = doc["provenance"]["runs"]
    assert len(runs) == 40 and set(runs) == {gen.bundle_name(c, p) for c, p in gen.evidence_keys()}
    assert all(
        v["scored"] == SYN_HISTORIES and v["run_id"].startswith("REF-") for v in runs.values()
    )
    pin = doc["provenance"]["review_pin"]
    assert pin["review_id"] == "REVIEW-syn" and pin["finished_at"] == PIN_FINISHED
    assert doc["provenance"]["manifest"]["path"] == ev.MANIFEST


def test_model_scenario(results: dict[str, dict[str, Any]]) -> None:
    doc = results["model"]
    assert doc["pattern"] == "MODEL" and doc["conclusive"]
    for ion in ("carbon", "oxygen"):
        for v in doc["bands"][ion].values():
            assert v["class"] == "M"
            assert v["mean_r_M"] == pytest.approx(0.0, abs=2e-3)
            assert v["mean_r_I"] == pytest.approx(-S / TABLE_FACTOR, abs=4e-3)


def test_mixed_scenario(results: dict[str, dict[str, Any]]) -> None:
    doc = results["mixed"]
    assert doc["pattern"] == "MIXED" and doc["conclusive"]
    assert set(_cls(doc, "carbon").values()) == {"I"}
    assert set(_cls(doc, "oxygen").values()) == {"M"}


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
    for ion in ("carbon", "oxygen"):
        for v in doc["bands"][ion].values():
            assert abs(v["mean_r_M"] - R0) < 1e-9
            assert v["se_b"] < 1e-9
            assert v["max_abs_thickness_difference"] == pytest.approx(B_LIN * 0.018, abs=1e-9)
    # fluctuation control: same planted r as C slab t -> zero difference
    assert doc["controls"]["fluctuations_off"]["max_abs_difference"] < 1e-9


def test_invalid_run_without_repeat_is_inconclusive(results: dict[str, dict[str, Any]]) -> None:
    doc = results["invalid"]
    c15 = gen.bundle_name(15, 0)
    assert doc["pattern"] == "INCONCLUSIVE" and doc["pattern_before_controls"] == "ICRU"
    assert list(doc["invalid_runs"]) == [c15] and "deficit" in doc["invalid_runs"][c15]
    assert doc["invalid_points"] == {"carbon case 15": "invalid run; no valid repeat"}
    assert doc["bands"]["carbon"]["[5,10)"]["n"] == 1
    assert doc["inconclusive_reasons"] == ["an evidence point is invalid without a valid repeat"]
    assert doc["provenance"]["runs"][c15]["scored"] is None


def test_second_invalid_and_consistency_failure(results: dict[str, dict[str, Any]]) -> None:
    doc = results["second_invalid"]
    assert doc["pattern"] == "INCONCLUSIVE" and not doc["conclusive"]
    assert set(doc["invalid_runs"]) == {gen.bundle_name(15, 0), gen.bundle_name(15, 6)}
    assert doc["invalid_points"] == {
        "carbon case 15": "repeat run invalid (second invalid of the point)",
        "oxygen case 19": f"{ev.CONSISTENCY_REASON}; no valid repeat",
    }
    assert list(doc["repeats_used"]) == [gen.bundle_name(15, 6)]
    o = doc["bands"]["oxygen"]["[20,50]"]  # 20 MeV/u
    assert o["n"] == 2 and o["max_abs_thickness_difference"] == pytest.approx(0.01, abs=2e-3)


def test_valid_repeats_restore_conclusiveness(results: dict[str, dict[str, Any]]) -> None:
    doc = results["repeat_ok"]
    assert doc["conclusive"] and doc["pattern"] == "ICRU"
    assert doc["invalid_points"] == {} and list(doc["invalid_runs"]) == [gen.bundle_name(15, 0)]
    assert doc["repeats_used"] == {
        gen.bundle_name(15, 6): f"original {gen.bundle_name(15, 0)} invalid: "
        + doc["invalid_runs"][gen.bundle_name(15, 0)],
        gen.bundle_name(
            19, 9
        ): f"original {gen.bundle_name(19, 3)} invalid: {ev.CONSISTENCY_REASON}",
    }
    assert (
        doc["bands"]["carbon"]["[5,10)"]["n"] == 2 and doc["bands"]["oxygen"]["[20,50]"]["n"] == 3
    )
    assert doc["bands"]["oxygen"]["[20,50]"]["max_abs_thickness_difference"] < 3e-3
    assert len(doc["provenance"]["runs"]) == 42


def _walk(x: Any, path: str = "") -> list[tuple[str, Any]]:
    out = [(path, x)]
    if isinstance(x, dict):
        for k, v in x.items():
            out += _walk(v, f"{path}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            out += _walk(v, f"{path}[{i}]")
    return out


def test_no_per_point_leak(results: dict[str, dict[str, Any]]) -> None:
    band_keys = {
        "n", "mean_r_M", "mean_r_I", "mean_s", "se_b", "sigma_b", "x_b", "class",
        "max_abs_thickness_difference",
    }  # fmt: skip
    for doc in results.values():
        for path, v in _walk(doc):
            if isinstance(v, list):
                assert len(v) <= 3, path
        assert set(doc["bands"]) == {"carbon", "oxygen"}
        for ion in doc["bands"].values():
            assert set(ion) == {"[2,5)", "[5,10)", "[10,20)", "[20,50]"}
            for band in ion.values():
                assert set(band) == band_keys
        for run in doc["provenance"]["runs"].values():
            assert set(run) == {"run_id", "seed", "histories", "scored"}
        assert set(doc["controls"]["helium"]) == {"pass", "max_abs_r_m"}
        assert set(doc["controls"]["fluctuations_off"]) == {"pass", "max_abs_difference"}
        for reason in [*doc["invalid_runs"].values(), *doc["invalid_points"].values()]:
            rest = reason.replace(ev.CONSISTENCY_REASON, "")
            assert not re.search(r"\d\.\d{3,}", rest), reason  # reasons carry no values
        text = json.dumps(doc)
        for bad in ("r_i_per", 'delta_e"', "de_m", "de_i", "exit_energy", "differences"):
            assert bad not in text.lower()


# --- authentication and binding (Codex REVIEW-cfdff398 finding 1) ----------------------------


def _world(
    base: dict[str, Any],
    tmp: Path,
    bundle: str = "",
    *,
    r: float = 0.0,
    commit: bool = True,
    **kw: Any,
) -> dict[str, Any]:
    """A copy of the base repository with ONE run of ``bundle`` (default C 10 MeV/u slab t)."""
    repo = tmp / "repo"
    shutil.copytree(base["repo"], repo)
    bundle = bundle or gen.bundle_name(17, 0)
    label = f"t-{tmp.name}"
    w = {**base, "root": tmp, "repo": repo}
    rid = _write_run(w, repo, label, bundle, r, np.random.default_rng(3), **kw)
    _write_index(repo, label, {bundle: rid})
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
    )


def _check(base: dict[str, Any], w: dict[str, Any]) -> dict[str, Any]:
    return dict(ev.check_run(_ctx(base, w), w["bundle"], w["run_id"]))


def test_authenticated_run_is_valid(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path)
    out = _check(base, w)
    assert out["valid"] and out["scored"] == SYN_HISTORIES and out["run_id"] == w["run_id"]
    assert out["de"] == pytest.approx(base["de_m"][w["bundle"]], rel=3e-3)


@pytest.mark.parametrize(
    "kw, match",
    [
        ({"tamper_record": True}, "attestation"),
        ({"started": "2026-11-16T11:59:59+00:00"}, "not after the review pin"),
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
        ({"stdout": _stdout(SYN_HISTORIES, fluct=False)}, "fluctuations enabled"),
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


def test_wrong_code_sha_raises(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path, code_sha=base["reviewed"])  # run executed at the bundles commit
    with pytest.raises(ev.RecordError, match="code_sha"):
        _check(base, w)
    with pytest.raises(ev.RecordError, match="i3fu2-review.json is not committed"):
        _ctx(base, w, code_sha=base["reviewed"])
    with pytest.raises(ev.RecordError, match="40-hex"):
        _ctx(base, w, code_sha="abc")


def test_manifest_and_pin_mismatch_raise(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path / "m")
    man = w["repo"] / ev.MANIFEST
    man.write_text(man.read_text() + " ")
    _commit(w["repo"], "manifest changed after the runs")
    with pytest.raises(ev.RecordError, match="frozen manifest .* differs"):
        _ctx(base, w)
    w = _world(base, tmp_path / "p")
    pin = w["repo"] / ev.REVIEW_PIN
    pin.write_text(pin.read_text().replace(PIN_FINISHED, "2026-11-16T10:00:00+00:00"))
    _commit(w["repo"], "pin forged later")
    with pytest.raises(ev.RecordError, match="review pin .* differs"):
        _ctx(base, w)
    w = _world(base, tmp_path / "u", commit=False)  # uncommitted records: fail closed
    with pytest.raises(ev.RecordError, match="not committed|differs from HEAD"):
        _check(base, w)


def test_bundle_entry_not_in_manifest_raises(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path)
    idx = w["repo"] / "records" / w["label"] / "index.json"
    doc = json.loads(idx.read_text())
    doc["runs"][0]["bundle"] = gen.SMOKE_NAME
    idx.write_text(json.dumps(doc))
    _commit(w["repo"], "index points at a smoke")
    with pytest.raises(ev.RecordError, match="not i3fu2"):
        _ctx(base, w)


def test_fluct_off_marker_must_be_pinned(
    base: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    w = _world(base, tmp_path, gen.bundle_name(13, 5))
    assert ev.FLUCT_OFF_MARKER == ""  # left for the orchestrator to pin from the smoke log
    with pytest.raises(ev.RecordError, match="not pinned"):
        _check(base, w)
    monkeypatch.setattr(ev, "FLUCT_OFF_MARKER", r"^this marker is absent$")
    with pytest.raises(ev.RecordError, match="marker is missing"):
        _check(base, w)
    monkeypatch.setattr(ev, "FLUCT_OFF_MARKER", FLUCT_OFF_TEST_MARKER)
    assert _check(base, w)["valid"]


def test_log_validity_failure_is_invalid_not_error(base: dict[str, Any], tmp_path: Path) -> None:
    w = _world(base, tmp_path, stdout=_stdout(SYN_HISTORIES) + "G4Exception: track stuck\n")
    out = _check(base, w)
    assert out["valid"] is False and "killed/stuck" in out["reason"]


def test_repeat_without_invalid_original_raises(base: dict[str, Any], tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    shutil.copytree(base["repo"], repo)
    w = {**base, "root": tmp_path, "repo": repo}
    _scenario(w, "extra", _r_model, repeats={gen.bundle_name(13, 6): 0.0})
    _commit(repo, "records")
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(ev, "FLUCT_OFF_MARKER", FLUCT_OFF_TEST_MARKER)
        with pytest.raises(ev.RecordError, match="without an invalid original"):
            _evaluate(w, "extra", repo)
