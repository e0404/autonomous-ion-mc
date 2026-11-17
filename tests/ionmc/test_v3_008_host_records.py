"""V3-008 revision 2: inventory of the committed signed host-run records of the sample chunks.

``validation/results/transport/v8r2-samples/host-records/<run_id>.json`` holds the complete signed
host-runner record of every sample chunk run of the two seals (``sample-seal.json``,
``sample-seal-heldout-r2stop3.json``); the evaluator (``r2_io.verify_host_records``) binds each
graded chunk to its record. Cache-free: reads the committed JSON and trust file only (no npz, no
published run directory). The records were delivered by the operator (IR-20261111-101840-EDCAE6).
Every signed ``argv`` must be one of the four reviewed wrapper templates
(``validation/scripts/ion_em/r2_wrapper_templates.py``) rendered for the run's sealed chunks.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
SAMPLES = REPO / "validation/results/transport/v8r2-samples"
HOST_RECORDS = SAMPLES / "host-records"
SEALS = (SAMPLES / "sample-seal.json", SAMPLES / "sample-seal-heldout-r2stop3.json")
TRUST = REPO / "infrastructure/keys/trusted_attestation_keys.json"
TEMPLATES = REPO / "validation/scripts/ion_em/r2_wrapper_templates.py"


def _attest():
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    from infrastructure.experiment_v3 import attest

    return attest


def _wrappers():
    spec = importlib.util.spec_from_file_location("r2_wrapper_templates_t", TEMPLATES)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses resolve their module
    spec.loader.exec_module(mod)
    return mod


def _runs() -> dict[str, dict[str, Any]]:
    """run_id -> {implementation_sha, output_manifest_sha256} of both seals' ``host_runs``."""
    runs: dict[str, dict[str, str]] = {}
    for path in SEALS:
        seal = json.loads(path.read_text(encoding="utf-8"))
        for h in seal["host_runs"]:
            assert h["run_id"] not in runs, f"{h['run_id']} listed twice"
            chunks = [c for c in seal["chunks"] if c["run_id"] == h["run_id"]]
            assert len({(c["case"], c["scale"]) for c in chunks}) == 1, h["run_id"]
            runs[h["run_id"]] = {
                "implementation_sha": seal["implementation_sha"],
                "output_manifest_sha256": h["output_manifest_sha256"],
                "kind": "heldout" if "heldout" in seal["schema"] else "main",
                "case": chunks[0]["case"],
                "ks": sorted(c["k"] for c in chunks),
                "scale": chunks[0]["scale"],
            }
    return runs


def test_seal_host_runs_count() -> None:
    assert len(_runs()) == 76 + 13


def test_committed_host_records_verify_and_bind() -> None:
    """Every file present must be ``<run_id>.json`` of a sealed run, signed by the trusted
    host-runner key, and carry that run id, implementation SHA and output-manifest digest."""
    runs = _runs()
    wrappers = _wrappers()
    matched: Counter = Counter()
    files = sorted(HOST_RECORDS.iterdir()) if HOST_RECORDS.is_dir() else []
    trusted = _attest().load_trusted(TRUST)
    for path in files:
        assert path.is_file() and not path.is_symlink(), path
        assert path.suffix == ".json" and path.stem in runs, f"{path.name}: not a sealed run"
        rec = json.loads(path.read_text(encoding="utf-8"))
        status, why = _attest().classify_record(rec, trusted, "host-runner")
        assert status == "verified", f"{path.name}: {status} ({why})"
        want = runs[path.stem]
        assert rec["run_id"] == path.stem
        assert rec["sha"] == want["implementation_sha"]
        assert rec["output_manifest_sha256"] == want["output_manifest_sha256"]
        argv = rec["argv"]
        assert isinstance(argv, list) and len(argv) == 3 and argv[:2] == ["python", "-c"]
        ids = [
            t.id
            for t in wrappers.TEMPLATES
            if t.kind == want["kind"] and _render(wrappers, t, want) == argv[2]
        ]
        assert len(ids) == 1, f"{path.name}: wrapper source matches templates {ids}"
        matched[ids[0]] += 1
    if files:
        assert matched == {"T1": 1, "T2": 30, "T3": 45, "T4": 13}


def _render(wrappers, t, want) -> str | None:
    try:
        return wrappers.render(t, want["case"], want["ks"], want["scale"])
    except ValueError:
        return None


def test_wrapper_template_digests() -> None:
    wrappers = _wrappers()
    wrappers.check_digests()
    assert [t.id for t in wrappers.TEMPLATES] == ["T1", "T2", "T3", "T4"]
    for t in wrappers.TEMPLATES:
        assert hashlib.sha256(t.source.encode("utf-8")).hexdigest() == t.sha256
        placeholder = wrappers.PLACEHOLDER_KS if t.multi_k else wrappers.PLACEHOLDER_K
        assert t.source.split("\n").count(placeholder) == 1


def test_every_sealed_run_has_a_host_record() -> None:
    present = {p.stem for p in HOST_RECORDS.glob("*.json")} if HOST_RECORDS.is_dir() else set()
    missing = sorted(set(_runs()) - present)
    assert not missing, f"{len(missing)} sealed runs without a committed host record"
