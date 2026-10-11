"""V3-005C C8b (plan Amendment 14a): calibration sharding and the combine consumer of the V7-R ``escaped_neutral``
calibration, the 3200-replicate constant and its propagation, the E1-B shard registration and the V6 fail-closed
state. CI scale: synthetic records, no Monte Carlo."""

# ruff: noqa: E501

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from tests.ionmc.test_run_suite import REPO, _load

RS = _load("run_suite")
V7R = _load("v7r")


def _seal(doc: dict[str, Any]) -> dict[str, Any]:
    norm = json.loads(json.dumps(doc))
    norm.pop("content_sha256", None)
    text = json.dumps(norm, sort_keys=True, separators=(",", ":"))
    return {**doc, "content_sha256": hashlib.sha256(text.encode()).hexdigest()}


FIT = {"family": "zero-inflated-gamma", "k": 1.5, "pi": 0.1}


def _share(total: int, sizes: tuple[int, ...], part: int) -> int:
    """Part ``part``'s share of ``total`` accepted replicates, filled greedily (never above the part size)."""
    before = sum(sizes[:part])
    return max(0, min(total - before, sizes[part]))


def _records(
    rates: dict[str, tuple[int, int, int]] | None = None, reps: int = 3200
) -> list[dict[str, Any]]:
    """One record per (surrogate, part); ``rates[margin] = (joint, paired, single)`` accepted in total
    (spread greedily over the parts)."""
    rates = rates or {}
    sizes = V7R.cal_part_sizes(reps)
    out = []
    for si, s in enumerate(V7R.CAL_SURROGATES):
        for p in range(V7R.CAL_PARTS):
            margins = {}
            for mi, m in enumerate(V7R.CAL_MARGINS):
                tot = rates.get(m, (2950, 2950, 2950) if m == "nominal" else (0, 0, 0))
                j, pa, sg = (_share(t, sizes, p) for t in tot)
                margins[m] = {"kappa": 1.0 + 0.01 * mi, "accepted_joint": j, "accepted_paired": pa,
                              "accepted_single": sg}  # fmt: skip
            out.append({"schema": V7R.CAL_SCHEMA, "surrogate": s, "part": p, "n": sizes[p], "margins": margins,
                        "fit": FIT, "reps_env": reps, "seeds": V7R.cal_part_seeds(si, p), "run_sha": "abc",
                        "suite": "lv5c"})  # fmt: skip
    return out


def test_constants_agree_between_run_suite_and_v7r() -> None:
    assert RS.V7R_CAL_REPS == V7R.V7R_CAL_REPS == 3200
    assert RS.V7R_CAL_SURROGATES == V7R.CAL_SURROGATES and RS.V7R_CAL_MARGINS == V7R.CAL_MARGINS
    assert RS.V7R_CAL_PARTS == V7R.CAL_PARTS == 15
    assert (
        V7R.cal_part_sizes(3200) == (214,) * 5 + (213,) * 10 and sum(V7R.cal_part_sizes(200)) == 200
    )


def test_part_seeds_are_derived_deterministically_and_distinct() -> None:
    seen = set()
    for si in range(2):
        for p in range(15):
            sd = V7R.cal_part_seeds(si, p)
            assert sd == V7R.cal_part_seeds(si, p)
            seen.add((sd["data"], sd["boot"]))
    assert len(seen) == 30


def test_partition_shape_order_and_registration() -> None:
    part = RS.V7R_CALIBRATION_PARTITION
    assert RS.V7R_CALIBRATION_SHARDS == len(part) == 31
    cases = [(s, p) for s in RS.V7R_CAL_SURROGATES for p in range(15)]
    for i, (s, p) in enumerate(cases):
        targets, k, want = part[i]
        assert targets == (RS.V7R_FILE,) and want == 1
        assert k == f"test_v7r_calibration_shard and {s}-part{p:02d}" == (
            f"test_v7r_calibration_shard and {RS.cal_case_id(s, p)}")  # fmt: skip
    assert part[14][1].endswith("archive_fitted-part14") and part[15][1].endswith(
        "rehearsal_block_sums-part00"
    )
    assert part[30] == RS.GAUSSIAN_SHARD and RS.GAUSSIAN_SHARD[2] == 11
    assert sum(w for _, _, w in part) == sum(RS.EXPECTED_CALIBRATION_CASES.values()) == 41
    names = [n.split("-", 1)[1] for n, _, _ in RS.suite_steps("lv5c", 1, 1.0)]
    i0 = names.index("pytest-v7r-calibration-s0")
    assert names[i0 : i0 + 33] == [*(f"pytest-v7r-calibration-s{k}" for k in range(31)),
                                   "v7r-calibration-combine", "v7r-combine"]  # fmt: skip


def test_reps_and_out_dir_are_exported_to_every_calibration_step_only() -> None:
    steps = RS.suite_steps("lv5c", 1, 1.0, out=Path("OUT"))
    for name, _, env in steps:
        if "pytest-v7r-calibration" in name:
            assert env["IONMC_V7R_CAL_REPS"] == "3200" and env["IONMC_V7R_CAL_OUT"] == "OUT/samples"
        else:
            assert "IONMC_V7R_CAL_REPS" not in env
    env = RS.step_environment({"IONMC_V7R_CAL_REPS": "7"}, "20-pytest-v7r-calibration-s3", {})
    assert env["IONMC_V7R_CAL_REPS"] == "3200"  # the caller's value cannot lower the count


def test_combine_passes_and_pools_the_counts() -> None:
    out = V7R.calibration_combine(_records({"paired_low": (3, 20, 5), "single_high": (1, 4, 30)}))
    assert (
        out["pass"] is True and out["reps"] == 3200 and out["part_sizes"] == [214] * 5 + [213] * 10
    )
    for s in V7R.CAL_SURROGATES:
        v = out["surrogates"][s]
        assert v["margins"]["paired_low"]["joint"] == 3 and v["margins"]["paired_low"]["n"] == 3200
        assert v["margins"]["nominal"]["joint"] == 2950 and v["kappas"]["paired_low"] == 1.01
        assert v["margins"]["single_high"]["gate_alone"]["single"]["rate"] == 30 / 3200


@pytest.mark.parametrize(
    ("mutate", "match"),
    [
        (lambda r: r.pop(7), "missing"),
        (lambda r: r.append(dict(r[0])), "duplicate"),
        (lambda r: r[4]["margins"]["paired_low"].update(kappa=2.0), "kappa differs"),
        (lambda r: r[4].update(fit={"k": 9}), "fit differs"),
        (lambda r: r[4].update(reps_env=200), "reps_env"),
        (lambda r: r[4].update(n=5), "is not"),
        (lambda r: r[4].update(seeds={"data": 1, "boot": 2}), "seeds"),
        (lambda r: r[4]["margins"]["nominal"].update(accepted_joint=-1), "invalid"),
        (lambda r: r[4]["margins"].pop("nominal"), "margins are not exactly"),
    ],
)
def test_combine_fails_closed_on_inconsistent_records(mutate, match: str) -> None:  # type: ignore[no-untyped-def]
    recs = _records()
    mutate(recs)
    with pytest.raises(SystemExit, match=match):
        V7R.calibration_combine(recs)


def test_combine_n_must_sum_to_3200() -> None:
    with pytest.raises(SystemExit, match="reps_env"):
        V7R.calibration_combine(_records(reps=200))
    assert V7R.calibration_combine(_records(reps=200), reps=200)["reps"] == 200


def _failures(rates: dict[str, tuple[int, int, int]]) -> list[str]:
    out = V7R.calibration_combine(_records(rates))
    return [f for v in out["surrogates"].values() for f in v["failures"]]


def test_combine_requirement_verdicts() -> None:
    assert _failures({}) == []
    assert any("joint" in f for f in _failures({"single_low": (51, 51, 51)}))  # CP95 upper > 0.02
    assert _failures({"single_low": (50, 50, 50)}) == []
    assert any(
        "power" in f for f in _failures({"nominal": (2908, 2908, 2908)})
    )  # one short of 2909
    assert _failures({"nominal": (2909, 2909, 2909)}) == []
    # gate alone by construction: refuted only if the CP95 lower bound exceeds 0.04 at the gate's own margin
    assert any(
        "refuted" in f and "paired gate" in f for f in _failures({"paired_high": (0, 160, 0)})
    )
    assert _failures({"paired_high": (0, 140, 0)}) == []  # rate 0.044, lower bound < 0.04: stands
    assert (
        _failures({"single_low": (0, 400, 0)}) == []
    )  # the paired gate at a single margin: not gated
    out = V7R.calibration_combine(_records({"paired_high": (0, 140, 0)}))
    alone = out["surrogates"]["archive_fitted"]["margins"]["paired_high"]["gate_alone"]["paired"]
    assert alone["construction_refuted"] is False
    assert alone["contrary_reading_cp95_upper_le_0.04"] is False and alone["cp95_upper"] > 0.04


def _write(d: Path, recs: list[dict[str, Any]]) -> None:
    d.mkdir(exist_ok=True)
    for r in recs:
        (d / f"v7r-calibration-{r['surrogate']}-p{r['part']:02d}.json").write_text(
            json.dumps(_seal(r), sort_keys=True)
        )


@pytest.fixture(scope="module")
def v5c() -> ModuleType:
    return _load("steps_v5c")


def test_loader_verifies_seal_run_identity_and_uniqueness(
    v5c: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("IONMC_RUN_SHA", "abc")
    monkeypatch.setenv("IONMC_RUN_SUITE", "lv5c")
    d = tmp_path / "samples"
    _write(d, _records())
    a = argparse.Namespace(dirs=[str(tmp_path)], partials_manifest=None)
    recs = v5c.load_calibration_records(a)
    assert len(recs) == 30 and len(a.attested_partials) == 30
    assert all(p["origin"] == "current-run" for p in a.attested_partials)
    assert V7R.calibration_combine(recs)["pass"] is True
    victim = d / "v7r-calibration-archive_fitted-p00.json"
    doc = json.loads(victim.read_text())
    victim.write_text(json.dumps({**doc, "accepted_joint": 0}))  # altered after sealing
    with pytest.raises(SystemExit, match="content_sha256"):
        v5c.load_calibration_records(
            argparse.Namespace(dirs=[str(tmp_path)], partials_manifest=None)
        )
    victim.write_text(json.dumps(_seal({**doc, "run_sha": "other"})))
    with pytest.raises(SystemExit, match="run_sha"):
        v5c.load_calibration_records(
            argparse.Namespace(dirs=[str(tmp_path)], partials_manifest=None)
        )
    victim.unlink()
    with pytest.raises(SystemExit, match="found 0 copies"):
        v5c.load_calibration_records(
            argparse.Namespace(dirs=[str(tmp_path)], partials_manifest=None)
        )
    victim.write_text(json.dumps(_seal(doc)))
    (tmp_path / "v7r-calibration-archive_fitted-p00.json").write_text(victim.read_text())
    with pytest.raises(SystemExit, match="found 2 copies"):
        v5c.load_calibration_records(
            argparse.Namespace(dirs=[str(tmp_path)], partials_manifest=None)
        )


def test_e1b_registration_three_shards_plus_compare_and_v6_stays_unsized() -> None:
    d = _load("steps_v5d")
    assert RS.E1B_SHARDS == d.E1B_SHARDS == 3
    assert d.e1b_shard_n(3) == 33340 and 3 * 33340 >= 100_000 and 33340 % 20 == 0
    names = [n.split("-", 1)[1] for n, _, _ in RS.suite_steps("lv5c", 1, 1.0)]
    assert [n for n in names if n.startswith("e1-b")] == [
        "e1-b-s0",
        "e1-b-s1",
        "e1-b-s2",
        "e1-b-compare",
    ]
    assert "v6-unsized" in names and RS.V6_SHARDS is None


def test_summarize_tag_and_timeout_floor_of_the_consumer() -> None:
    su = _load("summarize")
    assert su.expected_tag("50-v7r-calibration-combine") == "v7r-calibration-combine"
    assert RS.step_timeout_s("lv5c", "50-v7r-calibration-combine", 1500) >= 1500
    assert REPO.exists()
