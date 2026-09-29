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

CANONICAL_PROMPT_FILES = {
    "experiment/prompts/kickoff-v1.md":
        DOCS_ROOT / "experiment" / "prompts" / "kickoff-v1.md",
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


def copy_prompt_files() -> list[Path]:
    copied: list[Path] = []

    for source_name, destination in CANONICAL_PROMPT_FILES.items():
        source = REPO_ROOT / source_name

        if not source.is_file():
            raise RuntimeError(
                f"required canonical prompt missing: {source}"
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


def write_requirement_ledger() -> Path | None:
    """Render validation/requirement-ledger.json as a documentation page."""
    import json

    source = REPO_ROOT / "validation" / "requirement-ledger.json"
    if not source.is_file():
        return None
    ledger = json.loads(source.read_text(encoding="utf-8"))
    lines = [
        "# Requirement ledger",
        "",
        "Generated from `validation/requirement-ledger.json`, the canonical "
        "machine-readable record mapping every MUST requirement to its "
        "implementation status, planned evidence suites and tasks. Status "
        "vocabulary: " + ", ".join(f"`{s}`" for s in ledger["status_vocabulary"]) + ".",
        "",
        f"Ledger revision: {ledger['revision']} ({ledger['updated']}).",
        "",
        "| ID | Title | Status | Suites | Tasks | Notes |",
        "|---|---|---|---|---|---|",
    ]
    for item in ledger["requirements"]:
        lines.append(
            "| {id} | {title} | `{status}` | {suites} | {tasks} | {notes} |".format(
                id=item["id"],
                title=item["title"],
                status=item["status"],
                suites=", ".join(item["suites"]) or "-",
                tasks=", ".join(item["tasks"]) or "-",
                notes=item.get("notes", "").replace("|", "/"),
            )
        )
    lines.append("")
    destination = DOCS_ROOT / "validation" / "requirement-ledger.md"
    destination.write_text("\n".join(lines), encoding="utf-8")
    return destination


def main() -> int:
    decisions = copy_decisions()
    copied = copy_experiment_files()
    prompts = copy_prompt_files()
    v2 = REPO_ROOT / "experiment" / "v2"
    if v2.exists():
        destination = DOCS_ROOT / "experiment" / "v2"
        destination.mkdir(parents=True, exist_ok=True)
        for source in v2.glob("*.md"):
            shutil.copyfile(source, destination / source.name)
        for name in ("kickoff-v2.md", "kickoff-v2-isolated.md", "prepare-v2.md"):
            shutil.copyfile(REPO_ROOT / "experiment/prompts" / name, DOCS_ROOT / "experiment/prompts" / name)
    index = write_decision_index(decisions)
    write_requirement_ledger()

    print(f"generated {len(decisions)} decision page(s)")
    print(f"copied {len(copied)} canonical experiment document(s)")
    print(f"copied {len(prompts)} canonical prompt document(s)")
    print(f"generated {index.relative_to(REPO_ROOT)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
