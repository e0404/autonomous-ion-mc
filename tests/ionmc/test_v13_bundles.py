"""Pure-python guard of the V3-013 TOPAS bundles (plan A.6 row V-C1, sections C and D step 3).

Regeneration of the committed bundles and of the frozen manifest is byte-identical, the seeds follow
the section-C formula and lie inside the ``topas13`` family of ``test_v3_013_seed_families.py``, and
the sizing literals equal the plan's sizing table (exact strings of the plan file).
"""

import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PLAN = REPO / "validation" / "plans" / "v3-013-acceptance.md"
TOPAS = REPO / "validation" / "reference_cases" / "topas"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


GEN = _load("_gen_v13", REPO / "validation" / "scripts" / "stopping" / "gen_v13_bundles.py")
FAM = _load("_v3_013_seeds", HERE / "test_v3_013_seed_families.py")

# the plan's sizing table as exact row strings: case | MeV/u | total | R_plan | L | h_z | window | g
PLAN_ROWS = {
    1: "| 1 | 75 | 300 | 46.7 | 100 | 0.05 | [1, 3] | 2 |",
    2: "| 2 | 125 | 500 | 115.7 | 200 | 0.05 | [1, 3] | 2 |",
    3: "| 3 | 200 | 800 | 261.9 | 350 | 0.1 | [1, 3] | 5 |",
}
PARTS = (*range(0, 8), 19)


def _case_dirs():
    return sorted(d for d in TOPAS.iterdir() if d.is_dir() and "-v13-" in d.name)


def test_plan_sizing_table_literals():
    a6 = PLAN.read_text().split("### A.6", 1)[1].split("## B.", 1)[0]
    for row in PLAN_ROWS.values():
        assert row in a6
    for c in GEN.CASES:
        num, e, r_plan, length, hz = c
        total = GEN.A_NUCLEONS * e
        row = f"| {num} | {e} | {total} | {r_plan} | {length} | {GEN.fmt(hz)} | [1, 3] | "
        assert PLAN_ROWS[num].startswith(row)
    assert (GEN.N_EVIDENCE, GEN.N_SMOKE) == (25000, 2000)
    assert "4 parts × 2.5e4 primaries per case" in a6
    assert "(2 000 primaries; never evidence)" in PLAN.read_text().split("## C. Seeds", 1)[1]


def test_bundle_set_is_complete():
    names = {d.name for d in _case_dirs()}
    expect = {
        f"he-water-{e}mevu-v13-idd-em-part{p}" + ("-smoke" if p == 19 else "")
        for e in (75, 125, 200)
        for p in PARTS
    }
    assert names == expect and len(names) == 27


def test_regeneration_is_byte_identical(tmp_path):
    made = GEN.generate(tmp_path)
    assert len(made) == 27
    for nm, _ in made:
        for f in ("case.json", "input.txt"):
            assert (TOPAS / nm / f).read_bytes() == (tmp_path / nm / f).read_bytes(), (nm, f)
    assert GEN.manifest_text(tmp_path) == (TOPAS / "v13-bundles-manifest.json").read_text()
    man = json.loads((TOPAS / "v13-bundles-manifest.json").read_text())
    assert man["bundle_count"] == 27 and man["file_count"] == 54
    assert man["bundle_set_digest"] == GEN.BUNDLE_SET_DIGEST
    assert GEN.bundle_set_digest(TOPAS) == (GEN.BUNDLE_SET_DIGEST, 54)


def test_files_end_with_one_newline():
    for d in _case_dirs():
        for f in ("case.json", "input.txt"):
            b = (d / f).read_bytes()
            assert b.endswith(b"\n") and not b.endswith(b"\n\n"), (d.name, f)


def test_seeds_follow_the_formula_and_the_topas13_family():
    fam = FAM.topas13_seeds()
    assert GEN.TOPAS13_BASE == FAM.TOPAS13_BASE
    for e, case in ((75, 1), (125, 2), (200, 3)):
        for p in PARTS:
            tail = f"part{p}" + ("-smoke" if p == 19 else "")
            d = TOPAS / f"he-water-{e}mevu-v13-idd-em-{tail}"
            seed = 20610000 + 100 * p + case
            cj = json.loads((d / "case.json").read_text())
            assert cj["seeds"] == [seed]
            assert f"i:Ts/Seed = {seed}\n" in (d / "input.txt").read_text()
            assert seed in fam
            assert GEN.seed_of(p, case) == seed


def test_bundle_content_matches_the_sizing_and_conventions():
    for c in GEN.CASES:
        num, e, r_plan, length, hz = c
        for p in PARTS:
            tail = f"part{p}" + ("-smoke" if p == 19 else "")
            d = TOPAS / f"he-water-{e}mevu-v13-idd-em-{tail}"
            inp = (d / "input.txt").read_text()
            cj = json.loads((d / "case.json").read_text())
            n = 2000 if p == 19 else 25000
            assert cj["histories"] == n and f"i:So/Beam/NumberOfHistoriesInRun = {n}\n" in inp
            assert f"d:So/Beam/BeamEnergy = {4 * e} MeV\n" in inp
            assert f"i:Sc/IDD/ZBins = {round(length / hz)}\n" in inp
            assert f"d:Ge/Water/HLZ = {length / 2:g} mm\n" in inp
            assert 'b:Ph/ListProcesses = "True"\n' in inp
            assert "i:Ts/MaxStepNumber = 10000000\n" in inp
            assert "d:Ph/Default/CutForElectron = 1 km\n" in inp
            assert 'sv:Ph/Default/Modules = 1 "g4em-standard_opt4"\n' in inp
            assert (
                "Sc/LatX" not in inp and "Sc/LatY" not in inp and inp.count("Sc/IDD/Quantity") == 1
            )
            assert (
                f"V3-013 plan A.6 row V-C1, case {num}, part {p}; seal section D step 3; "
                in cj["rationale"]
            )
            assert (
                "not run before the reviewed seal and before the ionmc samples are committed"
                in cj["rationale"]
            )
            assert cj["expected_outputs"] == ["idd_dose.csv"]
