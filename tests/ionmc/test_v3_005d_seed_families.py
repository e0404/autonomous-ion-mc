"""V3-005D seed-family disjointness (pure python, no transport).

Evidence seed formulas (validation/scripts/transport/steps_v5*.py, steps_v6.py:170, v7r.py:85):
per-row seed = base + 1000 * r_index + k (r_index 0..20, k < 1000 is the conservative
envelope), V7-R shards base + 16000 + k. The formulas are restated here, not imported, to keep
the test side-effect free.

The planner's one-off numpy seed 20990001 was used in a scratch calculation, below the family start
(20990100 is the first seed used there); it consumes nothing.
"""

# ruff: noqa: E501

from __future__ import annotations

CONSUMED_BASES = (
    20261004,
    20271004,
    20281004,
    20291004,
    20301004,
    20311004,
    20321004,
    20331004,
    20341004,  # fmt: skip
    20351004,
    20361004,
    20371004,
    20381004,
    20391004,
    20401004,
    20421004,
    20441004,
    20451004,  # fmt: skip
    20471004,
    20481004,
    20491004,  # fmt: skip
)
FRESH_BASES = (20511004, 20521004, 20531004, 20541004, 20551004)
ENVELOPE = range(0, 21000)  # base + 1000 * r_index (0..20) + k (k < 1000), and base + 16000 + k

REHEARSAL_FAMILY = range(20505000, 20510000)
REHEARSAL_USED = range(
    20505000, 20526000
)  # fixed seeds 20505000-20508000, shards base + 16000 + k, largest 20525999
E_SHAPE = {20509500 + k for k in range(10)}
TOPAS = {20270000 + 1000 * row + k for row in range(1, 24) for k in (1, 2, 3)}

ELASTIC = range(20990000, 21000000)
VR_DIAG = range(20980000, 20990000)
ELASTIC_USED = {
    20990100,
    20990101,
    20990200,
    20990350,
    20990400,
    *(20991000 + c for c in range(10)),
    *(20991100 + c for c in range(10)),  # Amendment 20 (c): M2 diagnostic D4 seeds
    *(20991200 + c for c in range(10)),  # Amendment 21 (e): V-M2 D4 arm (report-only)
    *(20991300 + c for c in range(10)),  # Amendment 21 (d): M2' diagnostic D4 seeds
}  # 20990200: D-E6 CDF bootstrap base; 20990350 / 20990400: per-energy seeds (base + 150 / 200)
VR_CAL = range(20980000, 20985000)
VR_SANDBOX = range(20989900, 20990000)


def _evidence_seeds() -> set[int]:
    seeds: set[int] = set(TOPAS) | E_SHAPE | set(REHEARSAL_USED) | set(REHEARSAL_FAMILY)
    for base in (*CONSUMED_BASES, *FRESH_BASES):
        seeds.update(base + d for d in ENVELOPE)
    return seeds


def test_evidence_set_is_below_investigation_families() -> None:
    evidence = _evidence_seeds()
    assert max(evidence) == 20551004 + 20999 < min(VR_DIAG)
    assert not evidence & set(ELASTIC)
    assert not evidence & set(VR_DIAG)


def test_investigation_families_are_mutually_disjoint() -> None:
    assert not set(ELASTIC) & set(VR_DIAG)
    assert not set(VR_CAL) & set(VR_SANDBOX)
    assert set(VR_CAL) | set(VR_SANDBOX) <= set(VR_DIAG)


def test_elastic_seeds_inside_family() -> None:
    assert ELASTIC_USED <= set(ELASTIC)
    assert {20990100, 20990101, 20990200, 20990350, 20990400} | set(
        range(20991000, 20991010)
    ) | set(range(20991100, 20991110)) | set(range(20991200, 20991210)) | set(
        range(20991300, 20991310)
    ) == ELASTIC_USED
    assert 20990001 not in ELASTIC_USED and 20990001 < 20990100  # scratch-only seed, unused


def test_v7r2_transport_seeds_vs_rehearsal_and_hr_row6_block() -> None:
    r2 = {20531004 + 16000 + k for k in range(49)}
    assert not r2 & set(REHEARSAL_USED)
    # Computed finding: base 20531004 + 16000 == 20541004 + 6000, so the V7-R2 seeds coincide
    # exactly with the first 49 seeds of the hr base 20541004 r_index 6 block (20547004..20547999).
    hr_row6 = set(range(20541004 + 6000, 20541004 + 7000))
    assert hr_row6 == set(range(20547004, 20548004))
    assert r2 & hr_row6 == r2 == set(range(20547004, 20547053))


# -- Amendment 20 (f)4-(f)7: fixture, timing, V6-F seeds and the reserved hr r_index blocks ---------
LV_BASE, HR_BASE, V6F_BASE = 20531004, 20541004, 20551004
LV5D_USED_R = {
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    13,
    14,
    15,
    16,
    18,
    19,
}  # r_index rows that draw seeds (f)3
HR5D_USED_R = {1, 2, 3}  # (f)4: hr5d uses r_index 1-3 only (r_index 3 draws no seed)
HR_RESERVED = {0, 4, 5, 6, 8, 9, 17, 18, 19}  # (f)4 table; 1-3 follow the Amendment 14 (j) pattern
V6F = {
    "cpu": [V6F_BASE + 7000 + k for k in range(16)],
    "sens": [V6F_BASE + 8000 + k for k in range(8)],
    "cuda": [V6F_BASE + 9000 + k for k in range(16)],
}
FIXTURE_SEEDS = [20505000 + 16000 + k for k in range(4)]
TIMING_SEED = 20505000 + 9000 + 0
ORIGINAL_V6 = {b + k for b in (20493004, 20494004, 20501004) for k in range(1000)}


def _row_seeds(base: int, rows: set[int], k_max: int = 1000) -> set[int]:
    return {base + 1000 * r + k for r in rows for k in range(k_max)}


def test_v6f_seed_ranges_as_registered() -> None:
    assert V6F["cpu"][0] == 20558004 and V6F["cpu"][-1] == 20558019
    assert (V6F["sens"][0], V6F["sens"][-1]) == (20559004, 20559011)
    assert (V6F["cuda"][0], V6F["cuda"][-1]) == (20560004, 20560019)
    assert V6F_BASE + 7000 > max(
        REHEARSAL_USED
    )  # B + 7000 above the largest rehearsal seed used (20525999)


def test_v6f_seeds_satisfy_checks_a_to_e() -> None:
    every = [s for g in V6F.values() for s in g]
    assert len(set(every)) == len(every)  # (e) pairwise disjoint, no repeats
    assert not set(every) & set(REHEARSAL_USED)  # (b)
    assert not set(every) & ORIGINAL_V6  # (c) the never-run original V6 seeds
    consumed = {b + 1000 * r + k for b in CONSUMED_BASES for r in range(21) for k in range(1000)}
    assert not set(every) & consumed  # (a), over the conservative envelope
    other = _row_seeds(LV_BASE, set(range(21))) | _row_seeds(HR_BASE, HR5D_USED_R)
    assert not set(every) & other  # (d) the other rows of the executing task
    assert not set(every) & (set(ELASTIC) | set(VR_DIAG))
    assert not set(every) & (TOPAS | E_SHAPE)


def test_hr_reserved_r_index_blocks_coincide_only_with_declared_rows() -> None:
    """Bases 10 000 apart: lv r == hr r - 10 at equal k; V6-F r_index 7-9 == hr r_index 17-19 (Amendment 20 (f)4)."""
    for r in HR_RESERVED:
        hr = _row_seeds(HR_BASE, {r}, 100)
        lv_rows = {q for q in range(21) if hr & _row_seeds(LV_BASE, {q}, 100)}
        v6f_rows = {
            r_v for r_v, key in ((7, "cpu"), (8, "sens"), (9, "cuda")) if hr & set(V6F[key])
        }
        if r in (0, 4, 5, 6, 8, 9):
            assert lv_rows == {r + 10} and not v6f_rows
        else:  # 17-19
            assert not lv_rows and v6f_rows == {r - 10}
    assert not HR_RESERVED & HR5D_USED_R and HR5D_USED_R == {1, 2, 3}
    # the lv rows that draw seeds and their hr aliases: only the unused lv 11-13 alias hr 1-3 (Amendment 14 (j))
    assert {r for r in LV5D_USED_R if r - 10 in HR5D_USED_R} == {
        13
    }  # lv 13 (E1-B) aliases hr 3, which draws none


def test_v6f_hr_17_to_19_blocks_equal_the_v6f_seeds() -> None:
    assert set(V6F["cpu"]) <= _row_seeds(HR_BASE, {17}, 100)
    assert set(V6F["sens"]) <= _row_seeds(HR_BASE, {18}, 100)
    assert set(V6F["cuda"]) <= _row_seeds(HR_BASE, {19}, 100)


def test_rehearsal_fixture_and_timing_seeds() -> None:
    assert FIXTURE_SEEDS == [
        20521000,
        20521001,
        20521002,
        20521003,
    ]  # 4 x 2.25e6 at base 20505000, r_index 16
    assert TIMING_SEED == 20514000  # V6-F CUDA timing step: rehearsal base, r_index 9, k = 0
    assert set(FIXTURE_SEEDS) | {TIMING_SEED} <= set(REHEARSAL_USED)
    assert TIMING_SEED not in FIXTURE_SEEDS and not (set(FIXTURE_SEEDS) | {TIMING_SEED}) & (
        set(V6F["cpu"])
        | set(V6F["sens"])
        | set(V6F["cuda"])
        | {20531004 + 16000 + k for k in range(65)}
    )
    assert not (set(FIXTURE_SEEDS) | {TIMING_SEED}) & (
        set(ELASTIC) | set(VR_DIAG) | E_SHAPE | TOPAS
    )
    # never in a fresh evidence base: the rehearsal family lies below every fresh base envelope
    assert max(REHEARSAL_USED) < LV_BASE


def test_v7r2_calibration_roots_inside_the_calibration_family() -> None:
    roots = {20980000, 20980001, 20982000, 20982001, 20984000}
    assert roots <= set(VR_CAL) and not roots & (set(VR_SANDBOX) | ELASTIC_USED)
    assert not roots & _evidence_seeds()
