"""V3-008 revision 2: the registered (not run) ionmc sample steps ``steps_v8.py``: seeds, cases,
grids and configuration (validated only; nothing is transported). Seed guard vs the plan."""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "validation" / "scripts" / "transport"
sys.path.insert(0, str(SCRIPTS))
_spec = importlib.util.spec_from_file_location("steps_v8_under_test", SCRIPTS / "steps_v8.py")
assert _spec is not None and _spec.loader is not None
V8 = importlib.util.module_from_spec(_spec)
sys.modules["steps_v8_under_test"] = V8
_spec.loader.exec_module(V8)
sys.path.remove(str(SCRIPTS))

RESERVED = range(20572004, 20632004)
LV8_BASE, TOPAS_BASE = 20581004, 20619000


def test_seed_formula_and_r_index_30_to_41() -> None:
    assert sorted(V8.CASES) == list(range(1, 13))
    assert [V8.r_index_of(c) for c in range(1, 13)] == list(range(30, 42))
    assert V8.seed_of(1) == 20581004 + 3000
    assert V8.seed_of(12, 99) == 20581004 + 4100 + 99
    seeds = {V8.seed_of(c, k) for c in V8.CASES for k in range(100)}
    assert len(seeds) == 1200
    assert all(s in RESERVED for s in seeds)
    assert all(LV8_BASE <= s <= LV8_BASE + 100 * 60 + 99 for s in seeds)  # lv8 family r <= 60
    topas = {TOPAS_BASE + 100 * p + c for p in range(20) for c in range(1, 41)}
    assert not seeds & topas
    rerun = {V8.seed_of(c, k, V8.RERUN_BASE_1) for c in V8.CASES for k in range(100)}
    assert all(s in RESERVED for s in rerun) and not rerun & seeds
    for bad in (-1, 100):
        with pytest.raises(ValueError):
            V8.seed_of(1, bad)
    with pytest.raises(ValueError):
        V8.seed_of(13)


def test_cases_match_the_plan_table() -> None:
    text = (ROOT / "validation/plans/v3-008/r2-topas-target.md").read_text()
    rows = re.findall(
        r"^\| (\d+) \| (He|C|O) \| (\d+) \| [^|]+\| \d+ \| [\d.]+ \| (\d+) \| ([\d.]+) \| "
        r"([\d.]+) \| (main|held-out)",
        text,
        re.M,
    )
    assert len(rows) == 12
    names = {"He": "he", "C": "carbon", "O": "oxygen"}
    for num, ion, e_u, length, hz, hx, kind in rows:
        got = V8.CASES[int(num)]
        assert got == (names[ion], float(e_u), float(length), float(hz), float(hx), kind.split()[0])
    assert V8.N_HISTORIES == 1_000_000 and V8.N_BATCHES == 40
    assert "20581004 + 100" in (ROOT / "validation/plans/v3-008/r2-topas-target.md").read_text()


@pytest.mark.parametrize("case", range(1, 13))
def test_sample_config_validates_with_three_topas_grids(case: int) -> None:
    from ionmc.config import validate

    cfg = V8.sample_config(case)
    assert cfg.run.backend == "warp-cpu" and cfg.run.precision == "float64"
    assert cfg.run.n_histories == 1_000_000 and cfg.run.n_batches == 40
    assert cfg.run.seed == V8.seed_of(case)
    ion, e_u, length, hz, hx, _ = V8.CASES[case]
    idd, latx, laty = cfg.scoring
    assert (idd.name, latx.name, laty.name) == ("idd", "latx", "laty")
    assert idd.shape == (1, 1, round(length / hz)) and idd.spacing_mm == (120.0, 120.0, hz)
    assert latx.shape == (round(120 / hx), 1, round(length)) and latx.spacing_mm == (hx, 120.0, 1.0)
    assert laty.shape == (1, round(120 / hx), round(length))
    for g in cfg.scoring:  # the same edges as the TOPAS scorers (full 120 mm field from -60 mm)
        assert g.origin_mm[:2] == (-60.0, -60.0)
    eff = validate(cfg)
    assert eff.requested.physics.nuclear is False
    assert eff.summary()["physics"]["ls_requested"] == "off"
    assert eff.requested.source.kinetic_energy_mev == pytest.approx(
        e_u * eff.requested.source.projectile.a
    )
