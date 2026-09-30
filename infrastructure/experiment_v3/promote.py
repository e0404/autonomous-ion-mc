"""Trusted, retriable release promotion after re-evaluating exact-SHA evidence."""

import argparse
import json
import re
import subprocess
import tempfile
from pathlib import Path

from infrastructure.experiment_v3.bootstrap import IDENTITY
from infrastructure.experiment_v3.common import ROOT, event, git
from infrastructure.experiment_v3.release import evaluate


def gh(root, *args):
    return subprocess.check_output(["gh", *args], cwd=root, text=True).strip()


def check_ci(pr, sha, required=("pre-commit", "tests", "docs")):
    if pr.get("headRefOid") != sha:
        raise ValueError("PR head differs from qualified SHA")
    checks = pr.get("statusCheckRollup") or []
    names = set()
    for check in checks:
        name = check.get("name", check.get("context"))
        names.add(name)
        if check.get("conclusion", check.get("state")) not in ("SUCCESS",):
            raise ValueError(f"CI is not passing: {name}")
    if not set(required) <= names:
        raise ValueError("Required lightweight CI is missing")


def promote(root, plan, freeze_sha, report, artifacts, tag):
    if not re.fullmatch(r"ionmc-v3-[0-9]+\.[0-9]+\.[0-9]+", tag):
        raise ValueError("Tag must be ionmc-v3-X.Y.Z")
    root = Path(root)
    condition = json.loads((root / ".ionmc-condition.json").read_text())
    base, release = condition["integration_branch"], condition["release_branch"]
    if git(root, "branch", "--show-current") != base:
        raise ValueError("Promote from the clean v3 integration checkout")
    qualification = evaluate(root, plan, freeze_sha, report, artifacts)
    if not qualification["release_ready"]:
        return qualification
    sha = qualification["code_sha"]
    subprocess.run(
        ["git", "-C", str(root), "fetch", "origin", base, release], check=True
    )
    if git(root, "rev-parse", "origin/" + base) != sha:
        raise ValueError("Push/integrate exact qualified SHA before promotion")
    prs = json.loads(
        gh(
            root,
            "pr",
            "list",
            "--head",
            base,
            "--base",
            release,
            "--state",
            "all",
            "--json",
            "number,state,headRefOid",
            "--limit",
            "100",
        )
    )
    pr = next(
        (p for p in prs if p["headRefOid"] == sha and p["state"] in ("OPEN", "MERGED")),
        None,
    )
    with tempfile.TemporaryDirectory(prefix="ionmc-release-") as tmp:
        notes = Path(tmp) / "notes.md"
        notes.write_text(
            "IonMC v3 exact-SHA qualification\n\n"
            + json.dumps(qualification, indent=2)
            + "\n"
        )
        if pr is None:
            gh(
                root,
                "pr",
                "create",
                "--base",
                release,
                "--head",
                base,
                "--title",
                tag,
                "--body-file",
                str(notes),
            )
            prs = json.loads(
                gh(
                    root,
                    "pr",
                    "list",
                    "--head",
                    base,
                    "--base",
                    release,
                    "--json",
                    "number,state,headRefOid",
                )
            )
            pr = next(p for p in prs if p["headRefOid"] == sha)
        info = json.loads(
            gh(
                root,
                "pr",
                "view",
                str(pr["number"]),
                "--json",
                "headRefOid,statusCheckRollup,state",
            )
        )
        check_ci(info, sha)
        if info["state"] != "MERGED":
            gh(
                root,
                "pr",
                "merge",
                str(pr["number"]),
                "--merge",
                "--match-head-commit",
                sha,
            )
        subprocess.run(["git", "-C", str(root), "fetch", "origin", release], check=True)
        release_sha = git(root, "rev-parse", "origin/" + release)
        if git(root, "rev-parse", release_sha + "^{tree}") != git(
            root, "rev-parse", sha + "^{tree}"
        ):
            raise ValueError(
                "Promoted tree differs from the qualified tree; do not tag"
            )
        tags = git(root, "tag", "--list", tag)
        if tags:
            if git(root, "rev-parse", tag + "^{}") != release_sha:
                raise ValueError("Existing immutable tag points elsewhere")
        else:
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(root),
                    *IDENTITY,
                    "tag",
                    "-a",
                    tag,
                    release_sha,
                    "-m",
                    f"V3 qualified source {sha}; "
                    f"report {qualification['report_sha256']}",
                ],
                check=True,
            )
        subprocess.run(
            ["git", "-C", str(root), "push", "origin", "refs/tags/" + tag], check=True
        )
        existing = subprocess.run(
            ["gh", "release", "view", tag], cwd=root, capture_output=True
        )
        if existing.returncode:
            gh(
                root,
                "release",
                "create",
                tag,
                "--verify-tag",
                "--title",
                tag,
                "--notes-file",
                str(notes),
                str(report),
            )
    result = qualification | {
        "promoted": True,
        "release_sha": release_sha,
        "tag": tag,
        "pull_request": pr["number"],
    }
    event("release_promoted", sha=sha, details=result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=ROOT)
    p.add_argument("--plan", default="validation/release-plan.json")
    p.add_argument("--freeze-sha", required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--artifacts", type=Path, required=True)
    p.add_argument("--tag", required=True)
    a = p.parse_args()
    result = promote(a.root, a.plan, a.freeze_sha, a.report, a.artifacts, a.tag)
    print(json.dumps(result, indent=2))
    return 0 if result.get("promoted") else 1


if __name__ == "__main__":
    raise SystemExit(main())
