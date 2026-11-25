"""Seed-family guard for V3-006 (pure python, no transport imports).

The ``out6`` family is declared in validation/plans/v3-006-acceptance.md section C. The V3-008 and
V3-005 families come from tests/ionmc/test_v3_008_seed_families.py, the V3-013 ones from
tests/ionmc/test_v3_013_seed_families.py.
"""

import importlib.util
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLAN = HERE.parents[1] / "validation" / "plans" / "v3-006-acceptance.md"


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


V8 = _load("_v3_008_seeds", "test_v3_008_seed_families.py")
V13 = _load("_v3_013_seeds", "test_v3_013_seed_families.py")

OUT6_BASE = 20700004  # seed = base + 100*r_index + k
OUT6_R_MAX = 49
OUT6_K_MAX = 99
EXCLUDED_V3_005_BASES = (20681004, 20691004, 20701004)


def out6_seeds() -> set[int]:
    return {OUT6_BASE + 100 * r + k for r in range(OUT6_R_MAX + 1) for k in range(OUT6_K_MAX + 1)}


def test_declared_extent():
    s = out6_seeds()
    assert (min(s), max(s)) == (20700004, 20705003)
    assert len(s) == 5000


def test_disjoint_from_v3_008_and_v3_005_families():
    mine = out6_seeds()
    for name, fam in V8.v8_families().items():
        assert not mine & set(fam), name
    for name, ranges in V8.V3_005_FAMILIES.items():
        for r in ranges:
            assert not mine & set(r), (name, r)


def test_disjoint_from_v3_013_families():
    mine = out6_seeds()
    for name, fam in V13.families().items():
        assert not mine & set(fam), name


def test_excluded_v3_005_bases_would_collide():
    """The exclusion is necessary: the 21 000-seed windows of these bases reach out6, so the
    guard documents that they must not be used (they are not members of V3_005_FAMILIES)."""
    mine = out6_seeds()
    listed = {r.start for rs in V8.V3_005_FAMILIES.values() for r in rs}
    for base in EXCLUDED_V3_005_BASES:
        assert base not in listed
        assert mine & set(range(base, base + V8.ENVELOPE)), base


def test_plan_section_c_matches_constants():
    c = PLAN.read_text(encoding="utf-8").split("## C. Seeds", 1)[1].split("\n## D.", 1)[0]
    assert f"seed = {OUT6_BASE} + 100·r_index + k" in c
    assert f"k ≤ {OUT6_K_MAX}" in c and f"r_index ≤ {OUT6_R_MAX}" in c
    assert f"Extent {min(out6_seeds())}–{max(out6_seeds())}" in c
    for base in EXCLUDED_V3_005_BASES:
        assert str(base) in c
    nums = {int(n) for n in re.findall(r"\b(2070\d{4})\b", c)}
    assert nums <= {OUT6_BASE, min(out6_seeds()), max(out6_seeds())} | set(EXCLUDED_V3_005_BASES)
