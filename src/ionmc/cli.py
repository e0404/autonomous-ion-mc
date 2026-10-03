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
    data = sub.add_parser("data", help="manage cached external datasets")
    data_sub = data.add_subparsers(dest="data_command", metavar="action")
    for action, help_text in (
        ("list", "list registered datasets and their cache status"),
        ("fetch", "download a dataset into the cache and verify its SHA-256"),
        ("verify", "re-hash a cached dataset"),
        ("path", "print the path of a cached dataset"),
    ):
        p = data_sub.add_parser(action, help=help_text)
        if action != "list":
            p.add_argument("dataset_id")
        p.add_argument("--cache-dir", default=None, help="cache directory")
        p.add_argument("--offline", action="store_true", help="never use the network")
    return parser


def _run_data(args: argparse.Namespace) -> int:
    """Execute ``ionmc data <action>``."""
    from ionmc.data import cache
    from ionmc.data.acquire import OfflineError, fetch
    from ionmc.data.registry import DATASETS

    action = args.data_command
    if action is None:
        return 2
    cdir = cache.resolve_cache_dir(args.cache_dir)
    if action == "list":
        for ds in DATASETS.values():
            cached = cache.read_manifest(ds.id, cdir) is not None
            print(f"{ds.id}\t{'cached' if cached else 'missing'}\t{ds.description}")
        return 0
    try:
        if action == "fetch":
            print(fetch(args.dataset_id, cdir, offline=args.offline))
        elif action == "verify":
            print(f"OK {cache.verify(args.dataset_id, cdir)}")
        else:  # path
            print(fetch(args.dataset_id, cdir, offline=True))
    except (OfflineError, FileNotFoundError, cache.IntegrityError, KeyError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


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
    if args.command == "data":
        code = _run_data(args)
        if code == 2:
            parser.print_usage(sys.stderr)
        return code
    parser.print_usage(sys.stderr)
    return 2
