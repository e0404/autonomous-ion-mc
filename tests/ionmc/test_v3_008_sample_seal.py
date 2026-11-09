"""V3-008 revision 2: structural consistency of the committed ionmc sample seal record.

``validation/results/transport/v8r2-samples/sample-seal.json`` is the sealed record of the
twelve 1e6-history TD/TL samples (plan ``validation/plans/v3-008/r2-topas-target.md``, sections
"Order" and 7, dated chunk-plan note of 2026-11-08). It lists every chunk's identity (case, k, seed,
histories), the signed host-runner run it came from (run id, output-manifest digest) and the sha256
of its npz; the arrays themselves are not in Git. This test checks the record against the frozen
chunk plan only (counts, seeds, totals, uniqueness, SHA binding); it reads no physics value and the
signed run records are verified by the validation-manager gate, not here.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
RECORD = REPO / "validation/results/transport/v8r2-samples/sample-seal.json"
PLAN_N = {1: 4, 2: 5, 3: 10, 4: 4, 5: 8, 6: 20, 7: 4, 8: 8, 9: 10, 10: 4, 11: 4, 12: 5}
SEED_BASE, R_OFF = 20581004, 29
HEX64 = re.compile(r"^[0-9a-f]{64}$")
SHA40 = re.compile(r"^[0-9a-f]{40}$")
RUN_ID = re.compile(r"^RUN-\d{8}T\d{6}Z-[0-9a-f]{8}$")

pytestmark = pytest.mark.skipif(not RECORD.is_file(), reason="sample seal record not committed yet")


@pytest.fixture(scope="module")
def record() -> dict:
    return json.loads(RECORD.read_text(encoding="utf-8"))


def test_header(record: dict) -> None:
    assert record["schema"] == "ionmc-v8r2-sample-seal-1"
    assert record["task_id"] == "V3-008" and record["revision"] == 2
    assert SHA40.match(record["implementation_sha"])
    assert record["suite"] == "lv8" and record["single_process"] is True
    assert {int(k): v for k, v in record["chunk_plan"]["n"].items()} == PLAN_N
    assert record["local_verification"]["ok"] is True
    assert record["local_verification"]["problems"] == []


def test_chunks_follow_the_frozen_plan(record: dict) -> None:
    chunks = record["chunks"]
    assert len(chunks) == sum(PLAN_N.values()) == record["chunk_count"] == 86
    seen = Counter()
    for c in chunks:
        case, k = c["case"], c["k"]
        assert 1 <= k <= PLAN_N[case]
        seen[(case, k)] += 1
        assert c["r_index"] == R_OFF + case
        assert c["seed"] == SEED_BASE + 100 * (R_OFF + case) + k
        assert c["histories"] == round(1e6 / PLAN_N[case])
        assert c["npz"] == f"ionmc-case{case}-k{k}.npz"
        assert HEX64.match(c["npz_sha256"]) and HEX64.match(c["step_stdout_sha256"])
        assert HEX64.match(c["output_manifest_sha256"]) and HEX64.match(c["wrapper_stdout_sha256"])
        assert RUN_ID.match(c["run_id"])
        assert c["npz_bytes"] > 0 and c["wall_s"] > 0
        assert c["counters_all_zero"] is True
        att = c["host_runner_attestation"]
        assert att["signer"] == "host-runner" and re.match(r"^[0-9a-f]{16}$", att["key_id"])
        assert HEX64.match(att["payload_sha256"]) and re.match(r"^[0-9a-f]{128}$", att["signature"])
    assert all(n == 1 for n in seen.values())
    assert set(seen) == {(case, k) for case, n in PLAN_N.items() for k in range(1, n + 1)}


def test_totals_hashes_and_run_binding(record: dict) -> None:
    chunks = record["chunks"]
    per_case = Counter()
    for c in chunks:
        per_case[c["case"]] += c["histories"]
    assert all(per_case[case] == 1_000_000 for case in PLAN_N)
    assert {int(k): v["histories"] for k, v in record["per_case"].items()} == dict(per_case)
    npz_hashes = [c["npz_sha256"] for c in chunks]
    assert len(set(npz_hashes)) == len(npz_hashes), "two chunks share an npz digest"
    assert record["total_npz_bytes"] == sum(c["npz_bytes"] for c in chunks)
    runs = {(r["run_id"], r["output_manifest_sha256"]) for r in record["host_runs"]}
    assert len(runs) == len(record["host_runs"]) == len({r["run_id"] for r in record["host_runs"]})
    assert {(c["run_id"], c["output_manifest_sha256"]) for c in chunks} == runs
    shared = Counter(c["run_id"] for c in chunks)
    for c in chunks:
        if shared[c["run_id"]] > 1:
            assert "run_shared_with" in c, f"case{c['case']} k{c['k']} shares a run without saying so"
    assert len(record["superseded_runs"]) >= 1
    for failed in record["failed_runs"]:
        assert failed["succeeded"] is False and RUN_ID.match(failed["run_id"])
        assert failed["run_id"] not in {r["run_id"] for r in record["host_runs"]}
