"""Create fresh v2 from pre-scientific v1 plus tracked infrastructure."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from infrastructure.experiment_v2.common import (
    ROOT,
    exact_state,
    file_hash,
    git,
    write_json,
)

# Explicit infrastructure boundaries. Never overlay src, scientific tests,
# validation, benchmarks, root scientific package config or v1 science decisions.
OVERLAY = (
    "infrastructure/",
    "tests/infrastructure/",
    "experiment/v1/",
    "experiment/v2/",
    "experiment/prompts/prepare-v2.md",
    "experiment/prompts/kickoff-v2.md",
    ".mcp.json",
    ".pre-commit-config.yaml",
    "AGENTS.md",
    "decisions/0035-v2-experiment-condition.md",
)
EXCLUDE = {"infrastructure/diagnostics/warp_execution_model_probe.py"}
IDENTITY = [
    "-c",
    "user.name=Autonomous IonMC Agent",
    "-c",
    "user.email=autonomous-ionmc-agent@users.noreply.github.com",
]


def create(destination, *, root=ROOT):
    root, destination = Path(root).resolve(), Path(destination).resolve()
    setup_sha = exact_state(root)
    provenance = json.loads((root / "experiment/v1/provenance.json").read_text())
    baseline = provenance["start_sha"]
    if destination.exists():
        raise ValueError(
            "Destination already exists; bootstrap never overwrites a checkout"
        )
    tracked = git(root, "ls-files").splitlines()
    selected = [
        p
        for p in tracked
        if p not in EXCLUDE
        and any(p == x or (x.endswith("/") and p.startswith(x)) for x in OVERLAY)
    ]
    if not selected or "experiment/v2/PROTOCOL.md" not in selected:
        raise ValueError("No committed v2 setup found")
    subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "worktree",
            "add",
            "-b",
            "v2/develop",
            str(destination),
            baseline,
        ],
        check=True,
    )
    for name in selected:
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(
            subprocess.check_output(
                ["git", "-C", str(root), "show", f"{setup_sha}:{name}"]
            )
        )
    # Executable modes (e.g. existing telemetry scripts) remain from baseline.
    settings = json.loads(
        (destination / "experiment/v2/claude-settings.json").read_text()
    )
    sibling = str(destination.parent / (destination.name + "-worktrees"))
    settings["permissions"]["additionalDirectories"] = [sibling]
    for key in ("allowRead", "allowWrite"):
        values = settings["sandbox"]["filesystem"][key]
        settings["sandbox"]["filesystem"][key] = [
            sibling if v.endswith("ion-mc-worktrees") else v for v in values
        ]
    write_json(destination / ".claude/settings.json", settings)
    # Baseline docs/CI remain infrastructure-only. Extend them with v2 pages/checks.
    ci = destination / ".github/workflows/ci.yml"
    text = ci.read_text().replace(
        "      - main\n", '      - main\n      - "v2/develop"\n      - "v2/main"\n'
    )
    ci.write_text(text)
    nav = destination / "zensical.toml"
    nav.write_text(
        nav.read_text().replace(
            '    "experiment/index.md",',
            '    "experiment/index.md",\n    { "V2 setup" = '
            '"experiment/v2/SETUP.md" },\n    { "V2 protocol" = '
            '"experiment/v2/PROTOCOL.md" },\n    { "V2 requirements" '
            '= "experiment/v2/REQUIREMENTS.md" },\n    { "V2 '
            'operations" = "experiment/v2/AGENTS.md" },\n    { "V2 '
            'telemetry" = "experiment/v2/TELEMETRY.md" },\n    { "V2 '
            'kickoff" = "experiment/prompts/kickoff-v2.md" },',
        )
    )
    with (destination / ".gitignore").open("a") as out:
        out.write(
            "\n# V2 generated docs and task data caches\n"
            "docs/experiment/v2/\n.ionmc-cache/\n.uv-cache/\n"
        )
    condition = {
        "schema_version": 2,
        "experiment_id": "experiment-v2",
        "baseline_sha": baseline,
        "setup_sha": setup_sha,
        "integration_branch": "v2/develop",
        "release_branch": "v2/main",
        "scientific_implementation_inherited": False,
        "overlay": {name: file_hash(destination / name) for name in selected},
    }
    write_json(destination / ".ionmc-condition.json", condition)
    if (destination / "src").exists() or (destination / "tests/ionmc").exists():
        raise RuntimeError("Contaminated scientific starting tree")
    subprocess.run(["git", "-C", str(destination), "add", "-A"], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(destination),
            *IDENTITY,
            "commit",
            "-m",
            "Freeze fresh v2 experiment infrastructure condition",
        ],
        check=True,
    )
    sha = exact_state(destination)
    subprocess.run(
        ["git", "-C", str(destination), "branch", "v2/main", baseline], check=True
    )
    subprocess.run(
        [
            "git",
            "-C",
            str(destination),
            *IDENTITY,
            "tag",
            "-a",
            "experiment-v2-start",
            sha,
            "-m",
            "Fresh scientific starting tree; inherited v1 "
            "infrastructure plus v2 condition, no v1 science",
        ],
        check=True,
    )
    return {
        "path": str(destination),
        "sha": sha,
        "setup_sha": setup_sha,
        "next": "Run preflight, publish v2 branches through maintainer "
        "setup, then use kickoff-v2.md",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--destination", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(create(a.destination), indent=2))


if __name__ == "__main__":
    main()
