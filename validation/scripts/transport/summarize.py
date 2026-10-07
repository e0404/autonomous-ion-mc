"""Verify ``run_suite.py`` archives, attest the source, combine subsets (fail closed).

Usage::

    python summarize.py DIR --expected-sha SHA [--attest-sha SHA]
    python summarize.py --combine DIR [DIR ...] --expected-sha SHA --out FILE

Single archive: ``environment.txt`` carries the expected SHA; ``manifest.txt`` lists exactly the
step files present; every step file has the four header lines with the same SHA and an ``# exit=``
trailer as its last line. A step passes iff its exit code is 0 (steps that print a JSON document
between ``#JSON-BEGIN`` and ``#JSON-END`` must also have ``"pass": true``). The summary is a
``subset`` if the manifest is not the complete step list of the suite (``--only``); a subset, a
run with reduced history counts or a run whose source is neither from a clean git tree
(``tree_dirty=no``, ``sha_source=git``) nor attested is never ``conformant``. ``--attest-sha``
(run where git exists) compares the ``source_hashes`` of the archive with the blobs of
``git ls-tree -r SHA`` and records the attestation in ``summary.json``.

Single-process diagnostic archives (``execution_mode=single-process-diagnostic`` in
``environment.txt``) may contain steps with ``# status: deferred`` (only those listed in
``run_suite.DEFERRED_STEPS``): they are reported in ``deferred_steps``, ``pass`` is judged over the
executed steps, and the archive is never ``conformant`` (reasons: the execution mode and
``deferred multiprocessing checks``); ``--combine`` carries the deferred list through.

``--combine``: every directory must verify, with the same suite, SHA, scale and
``source_hashes`` (the *identity*); the union of the manifests must be exactly the complete step
list of the suite, with no step duplicated and none missing. Only then is the combined summary
``conformant`` (if every part is clean or attested and none is reduced).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
QUALIFICATION_SEED_BASE = 20391004
"""Only an archive made with this seed base can be conformant. The bases 20261004 (rehearsal),
20271004 (T9 investigation), 20281004, 20291004, 20301004, 20311004, 20321004, 20331004 and
20341004 (observed at the head before V3-003D)
(first to sixth qualification attempts, consumed) are recorded as non-qualification evidence
and never qualify."""
V4_QUALIFICATION_SEED_BASE = 20401004
"""Qualification base of the suites ``lv4`` and ``hr4`` (V3-004); 20351004 is their rehearsal base
and 20361004, 20371004 and 20381004 are consumed (observed before amendments 4 and 5 of the plan). A base whose full-scale
results were observed is consumed (plan, section Seeds)."""
V5_QUALIFICATION_SEED_BASE = 20421004
"""Qualification base of the suite ``lv5`` (V3-005A); the 2043xxxx family are rehearsals."""
IDENTITY_KEYS = ("git_sha", "suite", "scale", "python_parts", "seed_base")


def run_suite_suites() -> tuple[str, ...]:
    sys.path.insert(0, str(HERE))
    import run_suite

    return tuple(run_suite.SUITES)


def parse_env(text: str) -> dict[str, Any]:
    env: dict[str, Any] = {}
    hashes: dict[str, str] = {}
    in_hashes = False
    for line in text.splitlines():
        if line == "source_hashes:":
            in_hashes = True
        elif in_hashes and line.startswith("  "):
            digest, path = line.strip().split("  ", 1)
            hashes[path] = digest
        elif "=" in line:
            k, v = line.split("=", 1)
            env[k] = v
    env["source_hashes"] = hashes
    return env


def identity(env: dict[str, Any]) -> str:
    """sha256 of the parts of an environment that identify the code and suite under test."""
    blob = json.dumps(
        {k: env.get(k) for k in IDENTITY_KEYS} | {"source_hashes": env["source_hashes"]},
        sort_keys=True,
    )
    return hashlib.sha256(blob.encode()).hexdigest()


STEP_TAGS = (
    ("pytest", None),
    ("lv5-throughput", "lv5-throughput"),
    ("n1-", "n1"),
    ("v2-combine", "v2-combine"),
    ("v2-probe-combine", "v2-probe-combine"),
    ("v2-probe-", "v2-probe-shard"),
    ("v2-", "v2-shard"),
    ("v3-lv", "v3-lv"),
    ("v3-workers-", "v3-workers"),
    ("v4-v4b", "v4-v4b"),
    ("x1", "x1"),
    ("e1", "e1"),
    ("r1-", "r1"),
    ("a9-part", "a9-part"),
    ("a9-compare", "a9-compare"),
    ("a7-", "a7"),
    ("a8-", "a8"),
    ("a11-", "a11"),
    ("a13-", "a13"),
    ("a15-", "a15"),
    ("a16-", "a16"),
    ("t12-python-sample", "t12-sample"),
    ("t12-accelerated-samples", "t12-sample"),
    ("t12-compare", "t12-compare"),
    ("t-r1-", "t-r1"),
    ("t13-", "t13"),
    ("t14-", "t14"),
    ("t10-", "t10"),
    ("t1-", "t1"),
    ("t2-", "t2"),
    ("t8-", "t8"),
    ("t9-", "t9"),
)


def expected_tag(name: str) -> str | None | bool:
    """Document tag a step must print (``None``: a pytest step has no document; ``False``: the
    step name is unknown, which is a failure)."""
    base = name.split("-", 1)[1] if "-" in name else name
    for prefix, tag in STEP_TAGS:
        if base.startswith(prefix):
            return tag
    return False


def parse_step(
    path: Path, sha: str, suite: str | None = None, seed_base: str | None = None
) -> dict[str, Any]:
    lines = path.read_text().splitlines()
    problems: list[str] = []
    if "# status: deferred" in lines[:8]:
        reason = next((x[10:] for x in lines[:8] if x.startswith("# reason: ")), "")
        if len(lines) < 5 or not lines[0].startswith("# command: ") or lines[-1] != "# exit=0":
            problems.append("malformed deferred step")
        if f"# git_sha: {sha}" not in "\n".join(lines[:4]):
            problems.append("SHA header missing or different")
        return {
            "file": path.name,
            "exit": 0,
            "pass": False,
            "status": "deferred",
            "reason": reason,
            "problems": problems,
            "reduced": False,
            "histories": None,
            "frozen_histories": None,
        }
    if len(lines) < 5 or not lines[0].startswith("# command: "):
        problems.append("missing header")
    header = "\n".join(lines[:4])
    if f"# git_sha: {sha}" not in header:
        problems.append("SHA header missing or different")
    m = re.fullmatch(r"# exit=(\d+)", lines[-1]) if lines else None
    if not m:
        problems.append("missing exit trailer")
    code = int(m.group(1)) if m else -1
    doc = None
    text = "\n".join(lines)
    if "#JSON-BEGIN" in text:
        blob = text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0]
        try:
            doc = json.loads(blob)
        except json.JSONDecodeError:
            problems.append("unparseable JSON document")
    tag = expected_tag(path.stem)
    if tag is False:
        problems.append("unknown step name")
    elif tag is not None:
        # every non-pytest step must print its result document, naming itself, the suite and the SHA
        if doc is None:
            problems.append("missing JSON result document")
        else:
            checks = [("step", tag), ("suite", suite), ("git_sha", sha)]
            if seed_base is not None:  # every statistical step used the archive's seed base
                checks.append(("seed_base", int(seed_base)))
            for key, want in checks:
                if doc.get(key) != want:
                    problems.append(f"document {key} is {doc.get(key)!r}, expected {want!r}")
    verdict = (
        code == 0
        and not problems
        and (tag is None or (doc is not None and doc.get("pass") is True))
    )
    sel = re.search(r"(\d+) deselected", text) if tag is None else None
    return {
        "file": path.name,
        "exit": code,
        "pass": bool(verdict),
        "status": "executed",
        "deselected_tests": int(sel.group(1)) if sel else None,
        "problems": problems,
        "reduced": bool(doc and doc.get("reduced")),
        "histories": doc.get("histories") if doc else None,
        "frozen_histories": doc.get("frozen_histories") if doc else None,
        "attestation": doc.get("attestation") if doc else None,
    }


def attestation_problems(att: dict[str, Any] | None, sha: str, env: dict[str, Any]) -> list[str]:
    """Defects of the partials attestation block of a combine step (None: not a combine step).
    An imported partial without a ``host_run_id``, a manifest digest that is not the one recorded
    in ``environment.txt`` or a run SHA other than the archive's makes the archive non-conformant;
    the protected host-runner records themselves cannot be verified here."""
    if att is None:
        return []
    out = []
    imported = [p for p in att.get("partials", []) if p.get("origin") == "imported"]
    if imported and not att.get("manifest_sha256"):
        out.append("imported partials without a recorded manifest sha256")
    out += [f"imported partial {p.get('name')} lacks a host_run_id"
            for p in imported if not str(p.get("host_run_id") or "").strip()]  # fmt: skip
    recorded = env.get("partials_manifest_sha256")
    if att.get("manifest_sha256") and recorded != att["manifest_sha256"]:
        out.append("combine manifest sha256 differs from the one recorded in environment.txt")
    if att.get("run_sha") != sha:
        out.append("combine attestation run_sha differs from the archive SHA")
    return out


def _git(*args: str) -> bytes | None:
    try:
        r = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout if r.returncode == 0 else None


def attest(env: dict[str, Any], sha: str) -> dict[str, Any]:
    """Compare the archive's ``source_hashes`` with the blobs of commit ``sha``."""
    sys.path.insert(0, str(HERE))
    import run_suite

    prefixes = run_suite.SOURCE_PREFIXES
    suite_files = run_suite.source_file_list(env.get("suite"))
    listing = _git("ls-tree", "-r", "--name-only", sha)
    result: dict[str, Any] = {"attested_sha": sha, "valid": False, "mismatches": []}
    if listing is None:
        result["mismatches"].append("git cannot list the commit")
        return result
    tracked = {
        p
        for p in listing.decode().splitlines()
        if (p.startswith(tuple(prefixes)) or p in suite_files)
        and not p.endswith(".pyc")
        and "__pycache__" not in p
    }
    recorded = env["source_hashes"]
    if env.get("git_sha") != sha:
        result["mismatches"].append(f"archive SHA {env.get('git_sha')} differs from {sha}")
    for path in sorted(tracked - set(recorded)):
        result["mismatches"].append(f"tracked but not hashed: {path}")
    for path in sorted(set(recorded) - tracked):
        result["mismatches"].append(f"hashed but not tracked: {path}")
    for path in sorted(tracked & set(recorded)):
        blob = _git("cat-file", "blob", f"{sha}:{path}")
        if blob is None or hashlib.sha256(blob).hexdigest() != recorded[path]:
            result["mismatches"].append(f"content differs: {path}")
    result["valid"] = not result["mismatches"]
    return result


def seed_blockers(seed_base: Any, suite: str | None = None) -> list[str]:
    """Reasons why an archive's seed base cannot qualify (empty for the qualification base)."""
    if seed_base is None:
        return ["seed_base not recorded in environment.txt"]
    if suite == "lv5":
        if int(seed_base) != V5_QUALIFICATION_SEED_BASE:
            return [
                f"seed_base {int(seed_base)} is not the qualification base "
                f"{V5_QUALIFICATION_SEED_BASE} (the 2043xxxx family are rehearsals; any other base "
                "is non-qualification evidence)"
            ]
        return []
    if suite in ("lv4", "hr4"):
        if int(seed_base) != V4_QUALIFICATION_SEED_BASE:
            return [
                f"seed_base {int(seed_base)} is not the qualification base "
                f"{V4_QUALIFICATION_SEED_BASE} (20351004 is the rehearsal base, 20361004, 20371004 and 20381004 are consumed, 2041xxxx are V3-003D rehearsals; any other "
                "base is non-qualification evidence)"
            ]
        return []
    if int(seed_base) != QUALIFICATION_SEED_BASE:
        return [
            f"seed_base {int(seed_base)} is not the qualification base {QUALIFICATION_SEED_BASE} "
            "(20261004 is the rehearsal base, 20271004 was used by the T9 investigation, "
            "20281004, 20291004, 20301004, 20311004, 20321004, 20331004 and 20341004 by the first to seventh "
            "qualification attempts: all are "
            "non-qualification evidence)"
        ]
    return []


def source_ok(env: dict[str, Any], attestation: dict[str, Any] | None) -> bool:
    clean = env.get("tree_dirty") == "no" and str(env.get("sha_source", "")).startswith("git")
    return bool(clean or (attestation and attestation.get("valid")))


def verify(d: Path, sha: str, attest_sha: str | None = None) -> dict[str, Any]:
    """Verify one archive and return its summary (does not write it)."""
    problems: list[str] = []
    env_path = d / "environment.txt"
    env = parse_env(env_path.read_text()) if env_path.exists() else {"source_hashes": {}}
    if f"git_sha={sha}" not in (env_path.read_text() if env_path.exists() else ""):
        problems.append("environment.txt lacks the expected git_sha")
    names = (d / "manifest.txt").read_text().split() if (d / "manifest.txt").exists() else []
    if not names:
        problems.append("manifest.txt missing or empty")
    if len(set(names)) != len(names):
        problems.append("manifest lists a step twice")
    present = sorted(p.stem for p in d.glob("[0-9][0-9]-*.txt"))
    if sorted(names) != present:
        problems.append(f"manifest {sorted(names)} differs from step files {present}")
    steps: dict[str, Any] = {}
    for n in names:
        p = d / f"{n}.txt"
        steps[n] = (
            parse_step(p, sha, env.get("suite"), env.get("seed_base"))
            if p.exists()
            else {"pass": False, "problems": ["missing"]}
        )
    suite = env.get("suite")
    subset = True
    if suite in run_suite_suites():
        sys.path.insert(0, str(HERE))
        import run_suite

        full = run_suite.full_step_names(suite, int(env.get("python_parts", 2)))
        subset = sorted(names) != sorted(full)
    else:
        problems.append("environment.txt lacks a valid suite")
    reduced = any(s.get("reduced") for s in steps.values())
    att_problems = [
        f"{n}: {m}"
        for n, s in steps.items()
        for m in attestation_problems(s.get("attestation"), sha, env)
    ]
    attestation = attest(env, attest_sha) if attest_sha else None
    mode = env.get("execution_mode", "standard")
    deferred = sorted(n for n, s in steps.items() if s.get("status") == "deferred")
    if deferred and suite in run_suite_suites():
        import run_suite as _rs

        allowed = set(_rs.deferred_step_names(suite, int(env.get("python_parts", 2))))
        if mode != _rs.EXECUTION_SINGLE_PROCESS:
            problems.append("deferred steps in an archive that is not single-process diagnostic")
        for n in deferred:
            if n not in allowed:
                problems.append(f"step {n} may not be deferred")
    ok = not problems and all(s["pass"] for s in steps.values() if s.get("status") != "deferred")
    src_ok = source_ok(env, attestation)
    blockers = seed_blockers(env.get("seed_base"), suite)
    if deferred and suite not in run_suite_suites():
        problems.append("deferred steps in an archive of an unknown suite")
    if mode != "standard":
        blockers = [*blockers, f"execution mode {mode} (diagnostic, not the qualification mode)"]
    if deferred:
        blockers = [*blockers, "deferred multiprocessing checks"]
    return {
        "execution_mode": mode,
        "single_process_env": env.get("single_process_env") or None,
        "deferred_steps": deferred,
        "git_sha": sha,
        "suite": suite,
        "pass": ok,
        "subset": subset,
        "conformant": bool(
            ok and not subset and not reduced and src_ok and not blockers and not att_problems
        ),
        "non_conformant_reasons": [*blockers, *att_problems],
        "partials_attestation_problems": att_problems,
        "seed_base": env.get("seed_base"),
        "reduced_history_counts": reduced,
        "tree_dirty": env.get("tree_dirty"),
        "sha_source": env.get("sha_source"),
        "source_ok": src_ok,
        "attestation": attestation,
        "identity": identity(env),
        "problems": problems,
        "steps": steps,
    }


def combine(dirs: list[Path], sha: str, attest_sha: str | None = None) -> dict[str, Any]:
    parts = [verify(d, sha, attest_sha) for d in dirs]
    problems: list[str] = []
    for d, p in zip(dirs, parts, strict=True):
        if not p["pass"]:
            problems.append(f"{d}: archive does not verify ({p['problems']})")
    if len({p["identity"] for p in parts}) != 1:
        problems.append("archives differ in suite, SHA, scale, python_parts or source hashes")
    seen: dict[str, str] = {}
    for d, p in zip(dirs, parts, strict=True):
        for name in p["steps"]:
            if name in seen:
                problems.append(f"step {name} appears in {seen[name]} and {d}")
            seen[name] = str(d)
    suite = parts[0]["suite"] if parts else None
    missing: list[str] = []
    if suite in run_suite_suites():
        sys.path.insert(0, str(HERE))
        import run_suite

        env = parse_env((dirs[0] / "environment.txt").read_text())
        full = run_suite.full_step_names(suite, int(env.get("python_parts", 2)))
        missing = sorted(set(full) - set(seen))
        extra = sorted(set(seen) - set(full))
        if missing:
            problems.append(f"missing steps: {missing}")
        if extra:
            problems.append(f"unexpected steps: {extra}")
    ok = not problems and all(p["pass"] for p in parts)
    att_problems = [m for p in parts for m in p["partials_attestation_problems"]]
    reduced = any(p["reduced_history_counts"] for p in parts)
    deferred = sorted({n for p in parts for n in p["deferred_steps"]})
    return {
        "execution_modes": sorted({p["execution_mode"] for p in parts}),
        "deferred_steps": deferred,
        "git_sha": sha,
        "suite": suite,
        "pass": ok,
        "complete": not missing and not problems,
        "conformant": bool(
            ok
            and not reduced
            and all(p["source_ok"] for p in parts)
            and all(not p["non_conformant_reasons"] for p in parts)
            and not att_problems
        ),
        "attestation": {
            "run_sha": sha,
            "manifests": sorted(
                {
                    a["manifest_sha256"]
                    for p in parts
                    for s in p["steps"].values()
                    if (a := s.get("attestation")) and a.get("manifest_sha256")
                }
            ),
            "partials": [
                {"step": n, **q}
                for p in parts
                for n, s in p["steps"].items()
                if s.get("attestation")
                for q in s["attestation"].get("partials", [])
            ],
            "protected_host_records_verified_by_code": False,
            "statement": "the manifest was built by the orchestrator from the protected "
            "host-runner PARTIAL stdout lines; this code cannot verify those records, compare "
            "the digests and host_run_ids above with them and with record_local_validation",
        },
        "non_conformant_reasons": sorted({r for p in parts for r in p["non_conformant_reasons"]}),
        "reduced_history_counts": reduced,
        "problems": problems,
        "missing_steps": missing,
        "directories": [str(d) for d in dirs],
        "parts": {
            str(d): {k: p[k] for k in ("pass", "subset", "source_ok")}
            for d, p in zip(dirs, parts, strict=True)
        },
        "steps": {n: s for p in parts for n, s in p["steps"].items()},
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("directory", nargs="?")
    ap.add_argument("--expected-sha", required=True)
    ap.add_argument("--attest-sha")
    ap.add_argument("--combine", nargs="+", metavar="DIR")
    ap.add_argument("--out", help="combined summary file (must not exist)")
    args = ap.parse_args(argv)
    if args.combine:
        if not args.out:
            raise SystemExit("--combine needs --out")
        out = Path(args.out)
        if out.exists():
            raise SystemExit(f"refusing to overwrite {out}")
        summary = combine([Path(d) for d in args.combine], args.expected_sha, args.attest_sha)
        out.write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    else:
        if not args.directory:
            raise SystemExit("a directory is required")
        d = Path(args.directory)
        summary = verify(d, args.expected_sha, args.attest_sha)
        (d / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("pass", "conformant")}))
    return 0 if summary["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
