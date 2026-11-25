"""Seed-family guard for V3-007 (pure python, no transport imports).

The ``inf7`` family is declared in validation/plans/v3-007-acceptance.md section C. The V3-006
family comes from test_v3_006_seed_families.py, V3-008/V3-005 from test_v3_008_seed_families.py and
V3-013 from test_v3_013_seed_families.py.
"""

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLAN = HERE.parents[1] / "validation" / "plans" / "v3-007-acceptance.md"


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


V6 = _load("_v3_006_seeds", "test_v3_006_seed_families.py")
V8 = _load("_v3_008_seeds", "test_v3_008_seed_families.py")
V13 = _load("_v3_013_seeds", "test_v3_013_seed_families.py")

INF7_BASE = 20800000  # seed = base + 100*r + k
INF7_R_MAX = 49
INF7_K_MAX = 99
EXCLUDED_V3_005_BASES = (20781004, 20791004, 20801004)


def inf7_seeds() -> set[int]:
    return {INF7_BASE + 100 * r + k for r in range(INF7_R_MAX + 1) for k in range(INF7_K_MAX + 1)}


def test_declared_extent():
    s = inf7_seeds()
    assert (min(s), max(s)) == (20800000, 20804999)
    assert len(s) == 5000


def test_disjoint_from_v3_006_family():
    assert not inf7_seeds() & V6.out6_seeds()


def test_disjoint_from_v3_008_and_v3_005_families():
    mine = inf7_seeds()
    for name, fam in V8.v8_families().items():
        assert not mine & set(fam), name
    for name, ranges in V8.V3_005_FAMILIES.items():
        for r in ranges:
            assert not mine & set(r), (name, r)


def test_disjoint_from_v3_013_families():
    mine = inf7_seeds()
    for name, fam in V13.families().items():
        assert not mine & set(fam), name


def test_excluded_v3_005_bases_would_collide():
    mine = inf7_seeds()
    listed = {r.start for rs in V8.V3_005_FAMILIES.values() for r in rs}
    for base in EXCLUDED_V3_005_BASES:
        assert base not in listed
        assert mine & set(range(base, base + V8.ENVELOPE)), base


def test_plan_section_c_matches_constants():
    c = PLAN.read_text(encoding="utf-8").split("## C. Seeds", 1)[1].split("\n## D.", 1)[0]
    assert f"seed = {INF7_BASE} + 100·r + k" in c
    assert f"k ≤ {INF7_K_MAX}" in c and f"r ≤ {INF7_R_MAX}" in c
    assert f"{min(inf7_seeds())}–{max(inf7_seeds())}" in c
    for base in EXCLUDED_V3_005_BASES:
        assert str(base) in c
