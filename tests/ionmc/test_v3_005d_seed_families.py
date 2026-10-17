"""V3-005D seed-family disjointness (pure python, no transport).

Evidence seed formulas (validation/scripts/transport/steps_v5*.py, steps_v6.py:170, v7r.py:85):
per-row seed = base + 1000 * r_index + k (r_index 0..20, k < 1000 is the conservative
envelope), V7-R shards base + 16000 + k. The formulas are restated here, not imported, to keep
the test side-effect free.

The planner's one-off numpy seed 20990001 was used in a scratch calculation, below the family start
(20990100 is the first seed used there); it consumes nothing.
"""

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
