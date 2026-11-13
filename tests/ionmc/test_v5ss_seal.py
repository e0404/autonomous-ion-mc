"""V3-008 V5-SS: inventory of the committed thin-slab sample seal and its signed host-run records.

``validation/results/transport/v8-ss/ss-sample-seal.json`` seals the eight 4e5-history V5-SS
thin-slab samples (cases 22, 25-31, one run each, k = 0). ``v8-ss/host-records/<run_id>.json`` holds
the complete signed host-runner record behind every entry. Cache-free: reads the committed JSON,
the trust file and the two step modules only (no npz, no run directory, no network). The test
re-derives every seal field from the frozen step constants (seed formula, window, npz name), binds
each entry to its signed record (signature, commit, exit state, attestation, output manifest,
stdout digest, the JSON document printed by the wrapper, the exact rendered wrapper ``argv``) and
requires the record set to equal the sealed run set. Everything fails closed; a negative check
shows that an altered seal is rejected.
"""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
SS_DIR = REPO / "validation/results/transport/v8-ss"
SEAL = SS_DIR / "ss-sample-seal.json"
HOST_RECORDS = SS_DIR / "host-records"
TRUST = REPO / "infrastructure/keys/trusted_attestation_keys.json"
STEPS_V8 = REPO / "validation/scripts/transport/steps_v8.py"
WRAPPER = REPO / "validation/scripts/ion_em/ss_wrapper_template.py"
KEY_ID = "4d02f521f3fb9344"
SHA40 = re.compile(r"^[0-9a-f]{40}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def _load(name: str, path: Path, extra_path: Path | None = None) -> Any:
    if extra_path is not None:
        sys.path.insert(0, str(extra_path))
    try:
        spec = importlib.util.spec_from_file_location(name, path)
        assert spec is not None and spec.loader is not None, path
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
    finally:
        if extra_path is not None:
            sys.path.pop(0)
    return mod


@pytest.fixture(scope="module")
def steps() -> Any:
    return _load("steps_v8_v5ss_seal", STEPS_V8, STEPS_V8.parent)


@pytest.fixture(scope="module")
def wrapper() -> Any:
    return _load("ss_wrapper_v5ss_seal", WRAPPER)


@pytest.fixture(scope="module")
def seal() -> dict:
    return json.loads(SEAL.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def trusted() -> Any:
    from infrastructure.experiment_v3 import attest

    return attest.load_trusted(TRUST)


def _json_doc(stdout: str) -> dict:
    assert stdout.count("#JSON-BEGIN") == 1 and stdout.count("#JSON-END") == 1
    body = stdout.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0]
    return json.loads(body)


def check_seal(seal: dict, steps: Any, wrapper: Any) -> None:
    """Header, per-case derivation and run list of the seal (raises AssertionError)."""
    assert seal["schema"] == "ionmc-v8-ss-sample-seal-1"
    assert seal["task_id"] == "V3-008" and seal["row"] == "V5-SS"
    assert seal["backend"] == "warp-cpu" and seal["precision"] == "float64"
    assert seal["single_process"] is True
    assert seal["histories"] == steps.SS_HISTORIES == 400_000
    assert seal["seed_base"] == steps.SS_SEED_BASE == 20631004
    assert seal["window_k"] == 2.5
    assert seal["wrapper_template_sha256"] == wrapper.SS_WRAPPER_SHA256
    wrapper.check_digest()
    assert SHA40.match(seal["implementation_sha"])
    assert seal["npz_subdir"] == "v5ss-samples"
    cases = seal["cases"]
    assert sorted(int(c["case"]) for c in cases) == sorted(steps.SS_CASES)
    assert len({int(c["case"]) for c in cases}) == len(cases)
    runs = [c["run_id"] for c in cases]
    assert len(set(runs)) == len(runs) == len(steps.SS_CASES)
    assert sorted(r["run_id"] for r in seal["host_runs"]) == sorted(runs)
    manifests = {c["run_id"]: c["output_manifest_sha256"] for c in cases}
    assert {r["run_id"]: r["output_manifest_sha256"] for r in seal["host_runs"]} == manifests
    for c in cases:
        case = int(c["case"])
        assert c["r_index"] == steps.SS_R_INDEX[case]
        assert c["k"] == 0
        assert c["seed"] == seal["seed_base"] + 100 * c["r_index"] + c["k"]
        assert c["histories"] == 400_000
        assert c["npz"] == f"ionmc-ss-case{case}-k0.npz"
        assert c["window_k"] == 2.5
        assert abs(c["window_rad"] - 2.5 * c["theta_anchor_rad"]) <= 1e-12 * c["window_rad"]
        assert 0 < c["theta_w_rad"] < c["theta_anchor_rad"]
        assert 0.98 < c["p_window"] < 0.995
        assert abs(c["n_in_window"] - round(c["p_window"] * 2 * c["histories"])) <= 1
        assert c["run_started_at"] < c["run_finished_at"]
        assert 0 < c["wall_s"] < 1500
        for key in ("npz_sha256", "step_stdout_sha256", "output_manifest_sha256"):
            assert HEX64.match(c[key]), key
        assert HEX64.match(c["wrapper_stdout_sha256"])
        assert HEX64.match(c["host_record_sha256"])
        assert c["npz_bytes"] > 0


def check_records(seal: dict, trusted: Any, wrapper: Any, records_dir: Path) -> None:
    """Bind every sealed case to its signed host-run record; no extra record files."""
    from infrastructure.experiment_v3 import attest

    impl = seal["implementation_sha"]
    for c in seal["cases"]:
        path = records_dir / f"{c['run_id']}.json"
        assert path.is_file(), path
        blob = path.read_bytes()
        rec = json.loads(blob.decode("utf-8"))
        assert attest.classify_record(rec, trusted, "host-runner")[0] == "verified"
        assert rec["attestation"]["key_id"] == KEY_ID
        assert rec["attestation"] == c["host_runner_attestation"]
        assert rec["run_id"] == c["run_id"]
        assert rec["sha"] == impl
        assert rec["task_id"] == "V3-008"
        assert rec["succeeded"] is True and rec["exit_code"] == 0
        assert rec["committed_state_unchanged"] is True
        assert rec["dirty_before"] is False and rec["dirty_after"] is False
        assert rec["timed_out"] is False
        assert rec["output_manifest_sha256"] == c["output_manifest_sha256"]
        assert rec["started_at"] == c["run_started_at"]
        assert rec["finished_at"] == c["run_finished_at"]
        assert hashlib.sha256(blob).hexdigest() == c["host_record_sha256"]
        out_sha = hashlib.sha256(rec["stdout"].encode("utf-8")).hexdigest()
        assert out_sha == rec["stdout_sha256"] == c["wrapper_stdout_sha256"]
        doc = _json_doc(rec["stdout"])
        for key in ("case", "k", "seed", "r_index", "histories", "npz", "npz_sha256"):
            assert doc[key] == c[key], key
        for key in ("theta_anchor_rad", "theta_w_rad"):
            assert doc[key] == c[key], key
        assert doc["seed_base"] == seal["seed_base"]
        assert doc["git_sha"] == impl
        assert doc["pass"] is True and doc["valid"] is True
        assert doc["suite"] == "v5ss" and doc["n_escaped"] == 400_000
        assert doc["counters"] and all(v == 0 for v in doc["counters"].values())
        want = ["python", "-c", wrapper.render_wrapper(int(c["case"]), impl, "v5ss-samples", "1")]
        assert rec["argv"] == want
    on_disk = {p.stem for p in records_dir.glob("*.json")}
    assert on_disk == {c["run_id"] for c in seal["cases"]}
    assert len(list(records_dir.iterdir())) == len(on_disk)


def test_seal_header_and_cases(seal: dict, steps: Any, wrapper: Any) -> None:
    check_seal(seal, steps, wrapper)


def test_records_bind_the_seal(seal: dict, trusted: Any, wrapper: Any) -> None:
    check_records(seal, trusted, wrapper, HOST_RECORDS)


def test_altered_seal_is_rejected(seal: dict, steps: Any, wrapper: Any, trusted: Any) -> None:
    anchor = copy.deepcopy(seal)
    anchor["cases"][0]["theta_anchor_rad"] *= 1.001
    with pytest.raises(AssertionError):
        check_seal(anchor, steps, wrapper)
    # the signed record binding alone also catches an altered anchor (window kept consistent)
    anchor["cases"][0]["window_rad"] *= 1.001
    check_seal(anchor, steps, wrapper)
    with pytest.raises(AssertionError):
        check_records(anchor, trusted, wrapper, HOST_RECORDS)
    seed = copy.deepcopy(seal)
    seed["cases"][0]["seed"] += 1
    with pytest.raises(AssertionError):
        check_seal(seed, steps, wrapper)
