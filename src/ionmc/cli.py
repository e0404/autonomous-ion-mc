"""Command-line interface for IonMC."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence

from ionmc._version import __version__


def build_parser() -> argparse.ArgumentParser:
    """Create the ``ionmc`` argument parser."""
    parser = argparse.ArgumentParser(prog="ionmc", description="IonMC command-line interface.")
    parser.add_argument("--version", action="version", version=f"ionmc {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="command")
    sub.add_parser("version", help="print the ionmc version")
    sub.add_parser("info", help="print the runtime environment as JSON")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI and return a process exit code."""
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return exc.code if isinstance(exc.code, int) else 2
    if args.command == "version":
        print(__version__)
        return 0
    if args.command == "info":
        from ionmc.environment import describe_environment

        print(json.dumps(describe_environment(), indent=2, sort_keys=True))
        return 0
    parser.print_usage(sys.stderr)
    return 2
