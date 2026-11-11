"""Seed-family guard for V3-008 (pure python, no transport imports).

The V3-008 families are declared in validation/plans/v3-008-acceptance.md section C.
The V3-005 families are restated as literals because that branch is unmerged.
Source of the V3-005 literals: tests/ionmc/test_v3_005d_seed_families.py on the
V3-005B branch, commit 06974d07 (V3_008_RESERVED, CONSUMED_BASES, FRESH_BASES,
REHEARSAL, TOPAS, E_SHAPE, ELASTIC, VR_DIAG). Re-check at integration.
"""

import re
from pathlib import Path

PLAN = Path(__file__).resolve().parents[2] / "validation" / "plans" / "v3-008-acceptance.md"

# ---- V3-008 (plan section C) ----
V3_008_RESERVED = range(20572004, 20632004)  # 20572004-20632003 inclusive
LV8_BASE = 20581004
HR8_BASE = 20591004
R_MAX = 60  # r_index in 0..60
K_MAX = 99  # k in 0..99
REHEARSAL_8 = range(20600000, 20610000)
TOPAS_8_BASE = 20619000  # seed = base + 100*part + case
TOPAS_PARTS = range(0, 20)
TOPAS_CASES = range(1, 41)
RERUN_BASE_1 = 20621004  # full r range
RERUN_BASE_2 = 20631004  # r_index <= 9 only
RERUN_BASE_3 = 20572004  # full r range; the unassigned leading gap of the block, below LV8_BASE
RERUN_R_MAX = {RERUN_BASE_1: R_MAX, RERUN_BASE_2: 9, RERUN_BASE_3: R_MAX}


def row_seeds(base, r_max=R_MAX):
    return range(base, base + 100 * r_max + K_MAX + 1)


def topas_seeds():
    return {TOPAS_8_BASE + 100 * p + c for p in TOPAS_PARTS for c in TOPAS_CASES}


def v8_families():
    fam = {
        "lv8": row_seeds(LV8_BASE),
        "hr8": row_seeds(HR8_BASE),
        "rehearsal8": REHEARSAL_8,
        "rerun1": row_seeds(RERUN_BASE_1, RERUN_R_MAX[RERUN_BASE_1]),
        "rerun2": row_seeds(RERUN_BASE_2, RERUN_R_MAX[RERUN_BASE_2]),
        "rerun3": row_seeds(RERUN_BASE_3, RERUN_R_MAX[RERUN_BASE_3]),
    }
    t = sorted(topas_seeds())
    fam["topas8"] = t
    return fam


# ---- V3-005 (restated; source above) ----
ENVELOPE = 21000  # base + 0..20999
V3_005_NEXT_BASE = 20641004
CONSUMED_BASES = range(20261004, 20491004 + 1, 10000)  # superset: gaps not excluded
FRESH_BASES = (20511004, 20521004, 20531004, 20541004, 20551004)
V3_005_FAMILIES = {
    "consumed": [range(b, b + ENVELOPE) for b in CONSUMED_BASES],
    "fresh": [range(b, b + ENVELOPE) for b in FRESH_BASES],
    "rehearsal": [range(20505000, 20526000)],
    "topas": [range(20270000 + 1000 * row + 1, 20270000 + 1000 * row + 4) for row in range(1, 24)],
    "e_shape": [range(20509500, 20509510)],
    "elastic": [range(20990000, 21000000)],
    "vr_diag": [range(20980000, 20990000)],
    "continuation": [range(V3_005_NEXT_BASE, V3_005_NEXT_BASE + ENVELOPE)],
}


def _bounds(fam):
    return min(fam), max(fam)


def test_every_v8_seed_inside_reservation():
    for name, fam in v8_families().items():
        lo, hi = _bounds(fam)
        assert lo in V3_008_RESERVED and hi in V3_008_RESERVED, name


def test_declared_extents_match_plan_arithmetic():
    fam = v8_families()
    assert _bounds(fam["lv8"]) == (20581004, 20587103)
    assert _bounds(fam["hr8"]) == (20591004, 20597103)
    assert _bounds(fam["topas8"]) == (20619001, 20620940)
    assert _bounds(fam["rerun2"])[1] == V3_008_RESERVED[-1]
    assert _bounds(fam["rerun3"]) == (20572004, 20578103)


def test_rerun3_in_leading_gap_of_reservation():
    """rerun3 starts at the first value of the block and ends below every other V3-008 family."""
    fam = v8_families()
    assert RERUN_BASE_3 == V3_008_RESERVED[0]
    hi = _bounds(fam["rerun3"])[1]
    assert hi < LV8_BASE
    assert all(hi < _bounds(v)[0] for k, v in fam.items() if k != "rerun3")


def test_v8_families_mutually_disjoint():
    fam = {k: set(v) for k, v in v8_families().items()}
    names = sorted(fam)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            assert not (fam[a] & fam[b]), (a, b)


def test_topas_formula_below_first_rerun_base():
    assert max(topas_seeds()) < RERUN_BASE_1
    assert len(topas_seeds()) == len(TOPAS_PARTS) * len(TOPAS_CASES)  # injective


def test_reservation_disjoint_from_v3_005():
    lo, hi = V3_008_RESERVED[0], V3_008_RESERVED[-1]
    for name, ranges in V3_005_FAMILIES.items():
        for r in ranges:
            assert r[-1] < lo or r[0] > hi, (name, r)


def test_v8_seeds_disjoint_from_v3_005_families():
    v8 = set()
    for fam in v8_families().values():
        v8.update(fam)
    for name, ranges in V3_005_FAMILIES.items():
        for r in ranges:
            assert not (v8 & set(r)), (name, r)


def test_no_2098_2099_seeds():
    for name, fam in v8_families().items():
        for s in (min(fam), max(fam)):
            assert not str(s).startswith(("2098", "2099")), name
    assert V3_008_RESERVED[-1] < 20980000


def test_plan_section_c_matches_constants():
    text = PLAN.read_text()
    c = text.split("## C. Seeds", 1)[1].split("\n## D.", 1)[0]
    assert f"{V3_008_RESERVED[0]}–{V3_008_RESERVED[-1]}" in c
    assert f"range({V3_008_RESERVED[0]}, {V3_008_RESERVED[-1] + 1})" in c
    assert f"lv8 qualification **{LV8_BASE}**" in c
    assert f"hr8 qualification **{HR8_BASE}**" in c
    assert f"continues at **{V3_005_NEXT_BASE}**" in c
    assert f"{REHEARSAL_8[0]}–{REHEARSAL_8[-1]}" in c
    assert f"`{TOPAS_8_BASE} + 100·part + case`" in c
    assert f"part {TOPAS_PARTS[0]}–{TOPAS_PARTS[-1]}" in c
    assert f"case {TOPAS_CASES[0]}–{TOPAS_CASES[-1]}" in c
    assert f"r_index ≤ {R_MAX}" in c and "k < 100" in c
    assert f"from **{RERUN_BASE_1}**, then {RERUN_BASE_2}" in c
    assert f"third re-run base **{RERUN_BASE_3}**" in c
    assert f"{RERUN_BASE_3}–{_bounds(row_seeds(RERUN_BASE_3))[1]}" in c
    assert f"{LV8_BASE}–{_bounds(row_seeds(LV8_BASE))[1]}" in c
    assert f"{HR8_BASE}–{_bounds(row_seeds(HR8_BASE))[1]}" in c
    assert re.search(rf"highest is {max(topas_seeds())}", c)
    assert "Never 2098xxxx or 2099xxxx" in c
