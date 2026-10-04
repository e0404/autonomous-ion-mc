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
IDENTITY_KEYS = ("git_sha", "suite", "scale", "python_parts")


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


def parse_step(path: Path, sha: str, suite: str | None = None) -> dict[str, Any]:
    lines = path.read_text().splitlines()
    problems: list[str] = []
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
            for key, want in (("step", tag), ("suite", suite), ("git_sha", sha)):
                if doc.get(key) != want:
                    problems.append(f"document {key} is {doc.get(key)!r}, expected {want!r}")
    verdict = (
        code == 0
        and not problems
        and (tag is None or (doc is not None and doc.get("pass") is True))
    )
    return {
        "file": path.name,
        "exit": code,
        "pass": bool(verdict),
        "problems": problems,
        "reduced": bool(doc and doc.get("reduced")),
        "histories": doc.get("histories") if doc else None,
        "frozen_histories": doc.get("frozen_histories") if doc else None,
    }


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
    listing = _git("ls-tree", "-r", "--name-only", sha)
    result: dict[str, Any] = {"attested_sha": sha, "valid": False, "mismatches": []}
    if listing is None:
        result["mismatches"].append("git cannot list the commit")
        return result
    tracked = {
        p
        for p in listing.decode().splitlines()
        if (p.startswith(tuple(prefixes)) or p in run_suite.SOURCE_FILES)
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
            parse_step(p, sha, env.get("suite"))
            if p.exists()
            else {"pass": False, "problems": ["missing"]}
        )
    suite = env.get("suite")
    subset = True
    if suite in ("lv", "hr"):
        sys.path.insert(0, str(HERE))
        import run_suite

        full = run_suite.full_step_names(suite, int(env.get("python_parts", 2)))
        subset = sorted(names) != sorted(full)
    else:
        problems.append("environment.txt lacks a valid suite")
    reduced = any(s.get("reduced") for s in steps.values())
    attestation = attest(env, attest_sha) if attest_sha else None
    ok = not problems and all(s["pass"] for s in steps.values())
    src_ok = source_ok(env, attestation)
    return {
        "git_sha": sha,
        "suite": suite,
        "pass": ok,
        "subset": subset,
        "conformant": bool(ok and not subset and not reduced and src_ok),
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
    if suite in ("lv", "hr"):
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
    reduced = any(p["reduced_history_counts"] for p in parts)
    return {
        "git_sha": sha,
        "suite": suite,
        "pass": ok,
        "complete": not missing and not problems,
        "conformant": bool(ok and not reduced and all(p["source_ok"] for p in parts)),
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
