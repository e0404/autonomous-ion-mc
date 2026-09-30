"""Check the practical v3 history boundary without retrieving excluded content."""

import argparse
import json
import os
import subprocess
from pathlib import Path

from infrastructure.experiment_v3.common import ROOT, git


def inspect(root=ROOT):
    """Inspect the integration clone, not an individual task worktree.

    This catches configuration mistakes, not a determined bypass through public
    network access. Excluded commit probes inspect object existence only.
    """
    root = Path(root).resolve()
    checks = {}
    try:
        condition = json.loads((root / ".ionmc-condition.json").read_text())
        policy = condition.get("history_access", {})
        checks["policy"] = policy.get("mode") == "v3-only"
        common = Path(
            git(root, "rev-parse", "--path-format=absolute", "--git-common-dir")
        )
        checks["independent_git_directory"] = (
            (root / ".git").is_dir()
            and not (root / ".git").is_symlink()
            and common.resolve() == root / ".git"
        )
        checks["no_shared_objects"] = not (
            (common / "objects/info/alternates").exists()
            or (common / "objects/info/http-alternates").exists()
            or os.environ.get("GIT_ALTERNATE_OBJECT_DIRECTORIES")
            or os.environ.get("GIT_OBJECT_DIRECTORY")
            or (common / "objects").is_symlink()
        )
        checks["complete_v3_ancestry"] = (
            git(root, "rev-parse", "--is-shallow-repository") == "false"
        )
        checks["only_origin"] = git(root, "remote").splitlines() == ["origin"]
        branches = {condition["integration_branch"], condition["release_branch"]}
        checks["v3_branch_names"] = branches == {"v3/develop", "v3/main"}
        expected = {f"+refs/heads/{b}:refs/remotes/origin/{b}" for b in branches}
        checks["restricted_fetch"] = (
            set(git(root, "config", "--get-all", "remote.origin.fetch").splitlines())
            == expected
        )
        checks["no_automatic_tags"] = (
            git(root, "config", "--get", "remote.origin.tagOpt") == "--no-tags"
        )
        checks["v3_branches_available"] = all(
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(root),
                    "show-ref",
                    "--verify",
                    "--quiet",
                    "refs/remotes/origin/" + b,
                ],
                check=False,
            ).returncode
            == 0
            for b in branches
        )
        refs = git(root, "for-each-ref", "--format=%(refname)").splitlines()
        allowed = {f"refs/heads/{b}" for b in branches} | {
            f"refs/remotes/origin/{b}" for b in branches
        }
        allowed.add("refs/remotes/origin/HEAD")
        checks["only_v3_refs"] = all(
            ref in allowed
            or ref.startswith(
                (
                    "refs/heads/task/",
                    "refs/remotes/origin/task/",
                    "refs/tags/ionmc-v3-",
                    "refs/tags/experiment-v3-",
                )
            )
            for ref in refs
        )
        excluded = policy.get("excluded_commit_ids", [])
        checks["excluded_objects_absent"] = bool(excluded) and all(
            subprocess.run(
                ["git", "-C", str(root), "cat-file", "-e", sha + "^{commit}"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            ).returncode
            != 0
            for sha in excluded
        )
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        return {"ready": False, "checks": checks, "error_type": type(exc).__name__}
    return {"ready": all(checks.values()), "checks": checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    result = inspect(parser.parse_args().root)
    print(json.dumps(result, indent=2))
    return 0 if result["ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
