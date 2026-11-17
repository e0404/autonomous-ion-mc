"""Seed-family guard for V3-013 (pure python, no transport imports).

The V3-013 families are declared in validation/plans/v3-013-acceptance.md section C and lie in the
gap 20610000-20618999 of the V3-008 reservation (transferred by a dated note in plan v3-008 C).
The V3-008 and V3-005 families are taken from the existing guard of
tests/ionmc/test_v3_008_seed_families.py.
"""

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLAN = HERE.parents[1] / "validation" / "plans" / "v3-013-acceptance.md"

_spec = importlib.util.spec_from_file_location(
    "_v3_008_seeds", HERE / "test_v3_008_seed_families.py"
)
assert _spec is not None and _spec.loader is not None
V8 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V8)

GAP = range(20610000, 20619000)  # 20610000-20618999 inclusive
TOPAS13_BASE = 20610000  # seed = base + 100*part + case
TOPAS13_PARTS = range(0, 20)  # parts 0-3 evidence, 4-7 repeats, 19 smokes, others reserved
TOPAS13_CASES = range(1, 41)
ION13_BASE = 20612004  # seed = base + 100*r_index + k
ION13_R_MAX = 60
ION13_K_MAX = 99


def topas13_seeds():
    return {TOPAS13_BASE + 100 * p + c for p in TOPAS13_PARTS for c in TOPAS13_CASES}


def ion13_seeds():
    return set(range(ION13_BASE, ION13_BASE + 100 * ION13_R_MAX + ION13_K_MAX + 1))


def families():
    return {"topas13": topas13_seeds(), "ion13": ion13_seeds()}


def test_families_lie_inside_the_transferred_gap_and_the_v3_008_reservation():
    for name, fam in families().items():
        assert min(fam) in GAP and max(fam) in GAP, name
        assert min(fam) in V8.V3_008_RESERVED and max(fam) in V8.V3_008_RESERVED, name


def test_declared_extents():
    assert (min(topas13_seeds()), max(topas13_seeds())) == (20610001, 20611940)
    assert (min(ion13_seeds()), max(ion13_seeds())) == (20612004, 20618103)


def test_families_mutually_disjoint():
    assert not topas13_seeds() & ion13_seeds()


def test_disjoint_from_every_v3_008_and_v3_005_family():
    mine = topas13_seeds() | ion13_seeds()
    for name, fam in V8.v8_families().items():
        assert not mine & set(fam), name
    for name, ranges in V8.V3_005_FAMILIES.items():
        for r in ranges:
            assert not mine & set(r), (name, r)


def test_plan_section_c_matches_constants():
    c = PLAN.read_text().split("## C. Seeds", 1)[1].split("\n## D.", 1)[0]
    assert f"{GAP[0]}–{GAP[-1]}" in c
    assert f"seed = {TOPAS13_BASE} + 100·part + case" in c
    assert f"Extent {min(topas13_seeds())}–{max(topas13_seeds())}" in c
    assert f"seed = {ION13_BASE} + 100·r_index + k" in c
    assert f"k ≤ {ION13_K_MAX}" in c and f"r_index ≤ {ION13_R_MAX}" in c
    assert f"Extent {min(ion13_seeds())}–{max(ion13_seeds())}" in c
