"""Verify a ``run_suite.py`` archive and write ``summary.json`` (fail closed).

Checks: ``environment.txt`` carries the expected SHA; ``manifest.txt`` lists exactly the step
files present; every step file has the four header lines with the same SHA and an ``# exit=``
trailer as its last line. A step passes iff its exit code is 0 (steps that print a JSON document
between ``#JSON-BEGIN`` and ``#JSON-END`` must also have ``"pass": true`` and are marked
``reduced`` when run with fewer histories than frozen). The overall verdict is the conjunction;
the exit status is non-zero for any failure or inconsistency.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


def parse_step(path: Path, sha: str) -> dict[str, Any]:
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
    verdict = code == 0 and not problems and (doc is None or doc.get("pass") is True)
    return {
        "file": path.name,
        "exit": code,
        "pass": bool(verdict),
        "problems": problems,
        "reduced": bool(doc and doc.get("reduced")),
        "histories": doc.get("histories") if doc else None,
        "frozen_histories": doc.get("frozen_histories") if doc else None,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("directory")
    ap.add_argument("--expected-sha", required=True)
    args = ap.parse_args(argv)
    d = Path(args.directory)
    problems: list[str] = []
    env = (d / "environment.txt").read_text() if (d / "environment.txt").exists() else ""
    if f"git_sha={args.expected_sha}" not in env:
        problems.append("environment.txt lacks the expected git_sha")
    names = (d / "manifest.txt").read_text().split() if (d / "manifest.txt").exists() else []
    if not names:
        problems.append("manifest.txt missing or empty")
    present = sorted(p.stem for p in d.glob("[0-9][0-9]-*.txt"))
    if sorted(names) != present:
        problems.append(f"manifest {sorted(names)} differs from step files {present}")
    steps = {}
    for n in names:
        p = d / f"{n}.txt"
        steps[n] = (
            parse_step(p, args.expected_sha)
            if p.exists()
            else {"pass": False, "problems": ["missing"]}
        )
    reduced = any(s.get("reduced") for s in steps.values())
    ok = not problems and all(s["pass"] for s in steps.values())
    summary = {
        "git_sha": args.expected_sha,
        "pass": ok,
        "conformant": bool(ok and not reduced),
        "reduced_history_counts": reduced,
        "problems": problems,
        "steps": steps,
    }
    (d / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(json.dumps({"pass": ok, "conformant": summary["conformant"]}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
