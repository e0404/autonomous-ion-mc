"""Command-line entry point ``ionmc``.

Implemented commands:

* ``ionmc version`` — package version and exact code identity (JSON).
* ``ionmc data fetch [--material NAME ...] [--offline]`` — acquire (or verify
  from the cache) the external stopping-power datasets needed for the given
  materials; prints provenance records.
* ``ionmc data list`` — list cached dataset records.
* ``ionmc data cache-dir`` — print the cache location.
* ``ionmc stopping --species S --material M --energy T [T ...]`` — electronic
  mass stopping power (MeV cm²/g) and CSDA range (g/cm²) from the default
  tables at kinetic energies per nucleon T (MeV/u), with table provenance.

Scientific transport commands are added by later tasks.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys

from ionmc.provenance import code_identity


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ionmc", description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("version", help="print the package version and code identity")

    data = sub.add_parser("data", help="external data cache")
    data_sub = data.add_subparsers(dest="data_command", required=True)
    fetch = data_sub.add_parser("fetch", help="acquire or verify datasets")
    fetch.add_argument(
        "--material",
        action="append",
        default=None,
        help="ionmc material name (repeatable); default: all with tabulated data",
    )
    fetch.add_argument("--offline", action="store_true", help="never use the network")
    data_sub.add_parser("list", help="list cached dataset records")
    data_sub.add_parser("cache-dir", help="print the cache directory")

    stop = sub.add_parser("stopping", help="stopping power and range from the tables")
    stop.add_argument("--species", required=True)
    stop.add_argument("--material", required=True)
    stop.add_argument("--energy", type=float, nargs="+", required=True, help="MeV/u")
    stop.add_argument("--offline", action="store_true")
    return parser


def _cmd_data(args: argparse.Namespace) -> int:
    from ionmc.data import DataCache, nist_star
    from ionmc.materials import get_material, list_materials

    cache = DataCache()
    if args.data_command == "cache-dir":
        print(cache.root)
        return 0
    if args.data_command == "list":
        records = [dataclasses.asdict(r) for r in cache.list_records()]
        json.dump(records, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return 0
    names = args.material or [
        m for m in list_materials() if nist_star.material_key(get_material(m))
    ]
    out = []
    for name in names:
        key = nist_star.material_key(get_material(name))
        if key is None:
            out.append(
                {"material": name, "status": "no tabulated dataset (analytic layer)"}
            )
            continue
        for program in ("PSTAR", "ASTAR"):
            _, record = nist_star.load_star_table(
                cache, program, key, offline=args.offline
            )
            out.append(
                {"material": name, "program": program, **dataclasses.asdict(record)}
            )
    json.dump(out, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


def _cmd_stopping(args: argparse.Namespace) -> int:
    from ionmc.data import DataCache
    from ionmc.physics.tables import build_stopping_table, default_low_energy_source

    source = default_low_energy_source(DataCache(), offline=args.offline)
    table = build_stopping_table(args.species, args.material, source)
    rows = [
        {
            "t_mev_per_u": t,
            "electronic_mev_cm2_g": float(table.stopping_at(t)),
            "csda_range_g_cm2": float(table.range_at(t)),
        }
        for t in args.energy
    ]
    json.dump({"provenance": table.provenance, "values": rows}, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "version":
        json.dump(code_identity(), sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return 0
    if args.command == "data":
        return _cmd_data(args)
    if args.command == "stopping":
        return _cmd_stopping(args)
    return 2  # pragma: no cover


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
