"""V3-008: inventory of the committed signed host-run records of the LV/HR suites.

``validation/results/transport/v8-lvhr/host-records/<run_id>.json`` holds the complete signed
host-runner record of every lv8 / hr8 step run behind the combined summaries
(``lv8-712f484a-steps22-135-single-process.json``, ``hr8-712f484a-single-process.json``) and of
the preserved first I9 run, the abandoned I12 shard, the superseded duplicate of step 70 and the
rejected hr8 launch at 712f484a; ``v8-lvhr/index.json`` is a derived index (only ``role`` and
``note`` are curated). Cache-free: reads the committed JSON and the trust file only (no npz, no run
directory, no network). The test re-renders the signed wrapper ``argv`` of every record, binds the
seed base to the declared one, re-derives the index and binds both combined summaries to the
evidence records (partial digests and step verdicts). Everything fails closed.
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
LVHR = REPO / "validation/results/transport/v8-lvhr"
HOST_RECORDS = LVHR / "host-records"
INDEX = LVHR / "index.json"
TRUST = REPO / "infrastructure/keys/trusted_attestation_keys.json"
TRANSPORT = REPO / "validation/results/transport"
SUMMARIES = {
    "hr8": TRANSPORT / "hr8-712f484a-single-process.json",
    "lv8": TRANSPORT / "lv8-712f484a-steps22-135-single-process.json",
}
BRANCH = "task/v3-008-helium-carbon-and-oxygen-electromagnetic-transport-effective"
KEY_ID = "4d02f521f3fb9344"
SHA_FINAL = "712f484a13dbee4a181d15addf127b1c58dc1cbc"
SHA_FIRST = "e602a58383e30036ab2ed17fa64fffc99d29d3af"
SEED_BASES = {
    ("lv8", SHA_FIRST): 20581004,
    ("lv8", SHA_FINAL): 20572004,
    ("hr8", SHA_FINAL): 20591004,
}
ROLES = {
    "evidence",
    "superseded_duplicate",
    "rejected_launch",
    "preserved_fail",
    "preserved_fail_shard",
    "abandoned_partial",
}
ROLE_COUNTS = {
    ("lv8", "712f484a", "evidence"): 114,
    ("hr8", "712f484a", "evidence"): 12,
    ("lv8", "e602a583", "preserved_fail_shard"): 6,
    ("lv8", "e602a583", "preserved_fail"): 1,
    ("lv8", "e602a583", "abandoned_partial"): 1,
    ("lv8", "712f484a", "superseded_duplicate"): 1,
    ("hr8", "712f484a", "rejected_launch"): 1,
}
CURATED = {"role", "note"}
CACHE = "/workspace/.ionmc-cache/transport-samples"


def _attest():
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    from infrastructure.experiment_v3 import attest

    return attest


def render(suite: str, sha: str, step: str, seed_base: int, import_steps: list[str]) -> str:
    """The reviewed one-line wrapper source of a single-process ``run_suite.py --only`` run."""
    sha8 = sha[:8]
    imports = ""
    if import_steps:
        dirs = ", ".join(f"'{CACHE}/{sha8}/{s}'" for s in import_steps)
        imports = (
            f", '--import-dirs', {dirs}, '--partials-manifest', '{CACHE}/{sha8}/manifest.json'"
        )
    return (
        "import os, sys, runpy; os.environ['PYTHONPATH'] = '/workspace/src'; "
        "os.environ['PYTHONUNBUFFERED'] = '1'; os.environ['PYTHONFAULTHANDLER'] = '1'; "
        "sys.path.insert(0, '/workspace/validation/scripts/transport'); "
        "sys.path.insert(0, '/workspace/src'); "
        f"sys.argv = ['run_suite.py', '--suite', '{suite}', "
        f"'--out', 'validation/generated/transport/{suite}-{sha8}-{step}', "
        f"'--expected-sha', '{sha}', '--single-process', '--workers', '1', "
        f"'--seed-base', '{seed_base}', '--only', '{step}'{imports}]; "
        "runpy.run_path('validation/scripts/transport/run_suite.py', run_name='__main__')"
    )


def _records() -> dict[str, dict[str, Any]]:
    assert HOST_RECORDS.is_dir(), f"{HOST_RECORDS} is missing"
    out: dict[str, dict[str, Any]] = {}
    for path in sorted(HOST_RECORDS.iterdir()):
        assert path.is_file() and not path.is_symlink(), f"{path.name}: not a regular file"
        assert path.suffix == ".json", f"{path.name}: not a .json file"
        out[path.stem] = json.loads(path.read_text(encoding="utf-8"))
    assert out, "no host records"
    return out


def _index() -> dict[str, Any]:
    return json.loads(INDEX.read_text(encoding="utf-8"))


def _flag_after(argv_list: list[str], flag: str) -> str:
    assert argv_list.count(flag) == 1, f"{flag} must occur exactly once"
    return argv_list[argv_list.index(flag) + 1]


def _sys_argv(source: str) -> list[str]:
    """The ``sys.argv`` list literal of a wrapper source (all items are plain quoted strings)."""
    m = re.search(r"sys\.argv = \[(.*?)\]; runpy", source)
    assert m, "no sys.argv in the wrapper source"
    return re.findall(r"'([^']*)'", m.group(1))


def _derive(rec: dict[str, Any]) -> dict[str, Any]:
    """Every index field except ``role`` and ``note``, re-derived from the record."""
    source = rec["argv"][2]
    args = _sys_argv(source)
    stdout = rec["stdout"]
    name = re.search(r"^== (.+): exit=(\d+)$", stdout, re.MULTILINE)
    imports = (
        args[args.index("--import-dirs") + 1 : args.index("--partials-manifest")]
        if ("--import-dirs" in args)
        else []
    )
    return {
        "run_id": rec["run_id"],
        "suite": _flag_after(args, "--suite"),
        "step": _flag_after(args, "--only"),
        "step_name": name.group(1) if name else None,
        "sha": rec["sha"],
        "seed_base": int(_flag_after(args, "--seed-base")),
        "import_steps": [i.rsplit("/", 1)[1] for i in imports],
        "started_at": rec["started_at"],
        "finished_at": rec["finished_at"],
        "duration_s": round(rec["duration_seconds"], 1),
        "exit_code": rec["exit_code"],
        "succeeded": rec["succeeded"],
        "committed_state_unchanged": rec["committed_state_unchanged"],
        "stdout_sha256": rec["stdout_sha256"],
        "output_manifest_sha256": rec["output_manifest_sha256"],
        "partials": [
            {"name": m.group(1), "sha256": m.group(2)}
            for m in re.finditer(r"^PARTIAL (\S+) ([0-9a-f]{64})$", stdout, re.MULTILINE)
        ],
        "attestation_key_id": rec["attestation"]["key_id"],
    }


def _by_id(index: dict[str, Any]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for e in index["records"]:
        assert e["run_id"] not in out, f"{e['run_id']} listed twice in the index"
        out[e["run_id"]] = e
    return out


def test_records_verify_and_are_clean() -> None:
    trusted = _attest().load_trusted(TRUST)
    recs = _records()
    assert len(recs) == 136, len(recs)
    for run_id, rec in recs.items():
        assert rec["run_id"] == run_id, run_id
        assert rec["task_id"] == "V3-008" and rec["schema_version"] == 1, run_id
        assert rec["branch"] == BRANCH, run_id
        assert rec["dirty_before"] is False and rec["dirty_after"] is False, run_id
        assert rec["committed_state_unchanged"] is True, run_id
        assert rec["timed_out"] is False, run_id
        assert rec["sha_after"] == rec["sha"], run_id
        assert rec["stdout_truncated"] is False, f"{run_id}: stdout truncated"
        status, why = _attest().classify_record(rec, trusted, "host-runner")
        assert status == "verified", f"{run_id}: {status} ({why})"
        assert rec["attestation"]["key_id"] == KEY_ID, run_id


def test_signed_argv_is_the_rendered_wrapper_at_the_declared_seed_base() -> None:
    index = _by_id(_index())
    for run_id, rec in _records().items():
        assert run_id in index, f"{run_id}: not in the index"
        d = _derive(rec)
        suite, sha = d["suite"], rec["sha"]
        assert (suite, sha) in SEED_BASES, f"{run_id}: undeclared (suite, sha) {suite} {sha[:8]}"
        assert d["seed_base"] == SEED_BASES[(suite, sha)], f"{run_id}: seed base {d['seed_base']}"
        assert isinstance(rec["argv"], list) and len(rec["argv"]) == 3, run_id
        assert rec["argv"][:2] == ["python", "-c"], run_id
        want = render(suite, sha, d["step"], d["seed_base"], d["import_steps"])
        assert rec["argv"][2] == want, f"{run_id}: argv differs from the rendered template"
        if suite == "lv8":
            assert re.fullmatch(r"\d{2,3}", d["step"]) and 22 <= int(d["step"]) <= 135, run_id
        elif d["step"] != "1":  # the one rejected launch ran with the unpadded step
            assert re.fullmatch(r"\d{2}", d["step"]) and 1 <= int(d["step"]) <= 12, run_id


def test_index_is_derived_from_the_records() -> None:
    index = _index()
    assert index["schema"] == "ionmc-v8-lvhr-host-record-index-1"
    recs = _records()
    entries = _by_id(index)
    assert set(entries) == set(recs), "index and record files differ"
    assert len(index["records"]) == len(recs)
    counts = Counter((e["suite"], e["sha"][:8], e["role"]) for e in index["records"])
    assert all(e["role"] in ROLES for e in index["records"])
    assert counts == Counter(ROLE_COUNTS), counts
    for run_id, rec in recs.items():
        want = _derive(rec)
        got = {k: v for k, v in entries[run_id].items() if k not in CURATED}
        assert got == want, f"{run_id}: index differs from the record"
    missing = index["runs_without_record"]
    assert len(missing) == 1 and missing[0]["run_id"] == "RUN-20261111T183304Z-bd523546", missing
    assert missing[0]["run_id"] not in recs


def test_every_evidence_record_succeeded() -> None:
    recs = _records()
    for e in _index()["records"]:
        if e["role"] != "evidence":
            continue
        rec = recs[e["run_id"]]
        assert rec["succeeded"] is True and rec["exit_code"] == 0, e["run_id"]
        assert "(failed_steps=0)" in rec["stdout"], f"{e['run_id']}: failed_steps not 0"


def test_combined_summaries_bind_to_the_evidence_records() -> None:
    recs = _records()
    entries = list(_index()["records"])
    for suite, path in SUMMARIES.items():
        s = json.loads(path.read_text(encoding="utf-8"))
        assert s["suite"] == suite and s["git_sha"] == SHA_FINAL, path.name
        assert s["pass"] is True, path.name
        evidence = {
            e["step_name"]: e
            for e in entries
            if e["role"] == "evidence" and e["suite"] == suite and e["sha"] == SHA_FINAL
        }
        n_ev = sum(
            1
            for e in entries
            if e["role"] == "evidence" and e["suite"] == suite and e["sha"] == SHA_FINAL
        )
        assert n_ev == len(evidence) == len(s["steps"]), f"{suite}: evidence/step count"
        for key in s["steps"]:
            assert key in evidence, f"{suite}: no evidence record for step {key}"
            stdout = recs[evidence[key]["run_id"]]["stdout"]
            assert f"== {key}: exit=0" in stdout.splitlines(), f"{suite}: {key} did not exit 0"
        partials = s["attestation"]["partials"]
        for p in partials:
            rec = recs.get(p["host_run_id"])
            assert rec is not None, f"{suite}: unknown host_run_id {p['host_run_id']}"
            entry = next(e for e in entries if e["run_id"] == p["host_run_id"])
            assert rec["sha"] == s["git_sha"] and entry["role"] == "evidence", p
            line = f"PARTIAL {p['name']} {p['content_sha256']}"
            assert line in rec["stdout"].splitlines(), f"{suite}: {line} not in the record"
            assert re.sub(r"^\d+-", "", entry["step_name"]) == p["name"].removesuffix(".json"), p
        if suite == "hr8":
            assert len(partials) == 10 and s["complete"] is True and n_ev == 12
        else:
            assert len(partials) == 110 and s["complete"] is False and n_ev == 114
            assert s["subset"] is True and len(s["missing_steps"]) == 21


def test_specific_records() -> None:
    recs = _records()
    entries = _index()["records"]

    def role(r: str) -> list[dict[str, Any]]:
        return [e for e in entries if e["role"] == r]

    (dup,) = role("superseded_duplicate")
    assert dup["run_id"] == "RUN-20261111T222703Z-13954013" and dup["step"] == "70"
    assert recs[dup["run_id"]]["exit_code"] == 0
    (ev70,) = [
        e for e in role("evidence")
        if e["suite"] == "lv8" and e["step"] == "70" and e["sha"] == SHA_FINAL
    ]  # fmt: skip
    assert ev70["run_id"] == "RUN-20261112T140322Z-964c767c"
    assert len(dup["partials"]) == 1 and dup["partials"] == ev70["partials"], "step 70 differs"

    (rej,) = role("rejected_launch")
    r = recs[rej["run_id"]]
    assert rej["run_id"] == "RUN-20261113T015423Z-d1d4986b"
    assert r["exit_code"] == 1 and r["stdout"] == "" and "does not select exactly" in r["stderr"]
    assert r["output_manifest_sha256"] is None and r["succeeded"] is False

    (fail,) = role("preserved_fail")
    r = recs[fail["run_id"]]
    assert fail["run_id"] == "RUN-20261111T132917Z-098c9e9a" and fail["step"] == "28"
    assert fail["suite"] == "lv8" and fail["sha"] == SHA_FIRST
    assert r["exit_code"] == 1 and '"pass": false' in r["stdout"]

    shards = role("preserved_fail_shard") + role("abandoned_partial")
    assert len(shards) == 7
    for e in shards:
        assert recs[e["run_id"]]["exit_code"] == 0, e["run_id"]
