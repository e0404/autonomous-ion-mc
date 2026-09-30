"""Bind external workers to registered tasks in the active independent clone."""

import re
from pathlib import Path

from infrastructure.experiment_v3.common import ROOT, git


def validate_task(task_id, worktree, *, root=ROOT):
    task_id = task_id.strip().upper()
    if not re.fullmatch(r"V3-[A-Z0-9][A-Z0-9._-]{0,60}", task_id):
        raise ValueError("Expected a safe V3- task ID")
    root = Path(root).resolve()
    expected = root.parent / (root.name + "-worktrees") / task_id.lower()
    if expected.is_symlink() or Path(worktree).resolve() != expected:
        raise ValueError("Worker path must equal the registered V3 task worktree")
    registered = git(root, "worktree", "list", "--porcelain").splitlines()
    if "worktree " + str(expected) not in registered:
        raise ValueError("Worktree is not registered in the active clone")
    common = Path(
        git(expected, "rev-parse", "--path-format=absolute", "--git-common-dir")
    )
    if common.resolve() != root / ".git":
        raise ValueError("Worktree shares an unexpected Git database")
    branch = git(expected, "branch", "--show-current")
    prefix = "task/" + task_id.lower()
    if branch != prefix and not branch.startswith(prefix + "-"):
        raise ValueError("Task ID and worktree branch disagree")
    return task_id, expected
