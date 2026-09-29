"""Command-line entry point.

Only ``ionmc version`` exists at this stage; scientific commands are added by
later tasks together with their documentation and tests.
"""

from __future__ import annotations

import argparse
import json
import sys

from ionmc.provenance import code_identity


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ionmc", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("version", help="print the package version and code identity")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "version":
        json.dump(code_identity(), sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return 0
    return 2  # pragma: no cover - argparse enforces the subcommand set


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
