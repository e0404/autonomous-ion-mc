#!/usr/bin/env python3

from __future__ import annotations

import shutil
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = REPO_ROOT / "docs"
GENERATED_ROOT = DOCS_ROOT / "generated"

DECISIONS_SOURCE = REPO_ROOT / "decisions"
DECISIONS_DEST = GENERATED_ROOT / "decisions"

CANONICAL_EXPERIMENT_FILES = {
    "EXPERIMENT.md": DOCS_ROOT / "experiment" / "protocol.md",
    "REQUIREMENTS.md": DOCS_ROOT / "experiment" / "requirements.md",
}


def reset_generated_decisions() -> None:
    if DECISIONS_DEST.exists():
        shutil.rmtree(DECISIONS_DEST)

    DECISIONS_DEST.mkdir(parents=True, exist_ok=True)


def copy_decisions() -> list[Path]:
    reset_generated_decisions()

    copied: list[Path] = []

    if not DECISIONS_SOURCE.exists():
        return copied

    for source in sorted(DECISIONS_SOURCE.glob("*.md")):
        destination = DECISIONS_DEST / source.name
        shutil.copyfile(source, destination)
        copied.append(destination)

    return copied


def copy_experiment_files() -> list[Path]:
    copied: list[Path] = []

    for source_name, destination in CANONICAL_EXPERIMENT_FILES.items():
        source = REPO_ROOT / source_name

        if not source.is_file():
            raise RuntimeError(
                f"required canonical documentation source missing: {source}"
            )

        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        copied.append(destination)

    return copied


def write_decision_index(decisions: list[Path]) -> Path:
    index = DOCS_ROOT / "experiment" / "decisions.md"

    lines = [
        "# Design decisions",
        "",
        "Consequential design decisions are recorded canonically under "
        "`decisions/` in the repository and mirrored here automatically "
        "during documentation builds.",
        "",
    ]

    if not decisions:
        lines.append("No design decisions have been recorded yet.")
        lines.append("")
    else:
        for decision in decisions:
            title = decision.stem

            for line in decision.read_text(encoding="utf-8").splitlines():
                if line.startswith("# "):
                    title = line[2:].strip()
                    break

            relative = Path("../generated/decisions") / decision.name
            lines.append(f"- [{title}]({relative.as_posix()})")

        lines.append("")

    index.write_text("\n".join(lines), encoding="utf-8")
    return index


def main() -> int:
    decisions = copy_decisions()
    copied = copy_experiment_files()
    index = write_decision_index(decisions)

    print(f"generated {len(decisions)} decision page(s)")
    print(f"copied {len(copied)} canonical experiment document(s)")
    print(f"generated {index.relative_to(REPO_ROOT)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
