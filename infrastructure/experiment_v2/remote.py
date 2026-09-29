"""Read-only launch checks for v2 branch protection and exact-head lightweight CI."""

import json
import subprocess
import urllib.parse

from infrastructure.experiment_v2.common import exact_state, git

REPOSITORY = "e0404/autonomous-ion-mc"
REQUIRED = {"pre-commit", "tests", "docs"}


def protection_ok(data, *, linear):
    return bool(
        data.get("required_pull_request_reviews") is not None
        and data.get("enforce_admins", {}).get("enabled") is True
        and data.get("allow_force_pushes", {}).get("enabled") is False
        and data.get("allow_deletions", {}).get("enabled") is False
        and data.get("required_status_checks", {}).get("strict") is True
        and REQUIRED <= set(data.get("required_status_checks", {}).get("contexts", []))
        and data.get("required_linear_history", {}).get("enabled") is linear
    )


def inspect(root):
    condition = json.loads((root / ".ionmc-condition.json").read_text())
    sha = exact_state(root)
    remote = git(root, "remote", "get-url", "origin")
    if remote not in (
        f"https://github.com/{REPOSITORY}",
        f"https://github.com/{REPOSITORY}.git",
    ):
        raise ValueError("Unexpected experiment origin")

    def api(path):
        return json.loads(
            subprocess.check_output(
                ["gh", "api", f"repos/{REPOSITORY}/{path}"],
                cwd=root,
                text=True,
                stderr=subprocess.DEVNULL,
            )
        )

    protections = {}
    for key, linear in [("integration_branch", True), ("release_branch", False)]:
        branch = condition[key]
        data = api("branches/" + urllib.parse.quote(branch, safe="") + "/protection")
        protections[branch] = protection_ok(data, linear=linear)
    pushed = git(
        root, "ls-remote", "origin", "refs/heads/" + condition["integration_branch"]
    ).split()[0]
    checks = api(f"commits/{sha}/check-runs?per_page=100&filter=latest")["check_runs"]
    passed = {
        c["name"]
        for c in checks
        if c["head_sha"] == sha
        and c["status"] == "completed"
        and c["conclusion"] == "success"
    }
    return {
        "sha": sha,
        "exact_sha_pushed": pushed == sha,
        "protected_branches": protections,
        "passed_checks": sorted(passed),
        "ready": pushed == sha and all(protections.values()) and REQUIRED <= passed,
    }
