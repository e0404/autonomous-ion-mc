"""Command-line entry point: ``python -m ionmc.data <command>``.

Commands
--------
``acquire NAME [NAME ...]`` or ``acquire --all``
    Download, verify and cache the named registered datasets (needs network).
``status``
    Show which registered datasets are present in the cache.
``cache-dir``
    Print the resolved cache directory.
"""

from __future__ import annotations

import argparse
import sys

from ionmc.data import cache, registry


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m ionmc.data")
    parser.add_argument(
        "--cache-dir", default=None, help="override the cache directory"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    acq = sub.add_parser("acquire", help="download and cache datasets")
    acq.add_argument("names", nargs="*")
    acq.add_argument("--all", action="store_true")
    acq.add_argument("--force", action="store_true")
    sub.add_parser("status", help="show cache status of registered datasets")
    sub.add_parser("cache-dir", help="print the resolved cache directory")
    args = parser.parse_args(argv)

    if args.command == "cache-dir":
        print(cache.cache_dir(args.cache_dir))
        return 0
    if args.command == "status":
        for name, spec in registry.DATASETS.items():
            state = "cached" if cache.is_cached(spec, args.cache_dir) else "missing"
            print(f"{name}\t{spec.version[:12]}\t{state}")
        return 0
    names = list(registry.DATASETS) if args.all else args.names
    if not names:
        parser.error("acquire needs dataset names or --all")
    unknown = [n for n in names if n not in registry.DATASETS]
    if unknown:
        parser.error(
            f"unknown dataset(s): {unknown}; known: {sorted(registry.DATASETS)}"
        )
    for name in names:
        path = cache.acquire(registry.DATASETS[name], args.cache_dir, force=args.force)
        print(f"{name}\t{path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
