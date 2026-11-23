"""Command-line interface for IonMC."""

from __future__ import annotations

import argparse
import contextlib
import json
import sys
import traceback
from collections.abc import Sequence
from typing import Any, TextIO

from ionmc._version import __version__


def build_parser() -> argparse.ArgumentParser:
    """Create the ``ionmc`` argument parser."""
    parser = argparse.ArgumentParser(prog="ionmc", description="IonMC command-line interface.")
    parser.add_argument("--version", action="version", version=f"ionmc {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="command")
    sub.add_parser("version", help="print the ionmc version")
    sub.add_parser("info", help="print the runtime environment as JSON")
    sub.add_parser("notices", help="print the packaged third-party notices")
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
    imp = data_sub.add_parser(
        "import", help="import a local file as a registered dataset (hash and size verified)"
    )
    imp.add_argument("path", help="local file to import")
    imp.add_argument("--dataset", required=True, dest="dataset_id", help="registered dataset id")
    imp.add_argument("--cache-dir", default=None, help="cache directory")
    bld = data_sub.add_parser("build", help="build a derived table from cached datasets")
    bld.add_argument("table", choices=["nuclear-proton"])
    bld.add_argument("--cache-dir", default=None, help="cache directory")
    bld.add_argument("--points-per-decade", type=int, default=50)
    bld.add_argument("--diagnostic-events", type=int, default=20000)
    bld.add_argument(
        "--strict", action="store_true", help="fail when a lambda node does not converge to 1e-3"
    )
    run = sub.add_parser("run", help="run a TOML configuration and write a result bundle")
    run.add_argument("config", metavar="CONFIG.toml")
    run.add_argument("--out", default=None, help="bundle directory (default: [output].path)")
    run.add_argument("--overwrite", action="store_true", help="replace an existing bundle")
    chk = sub.add_parser(
        "check", help="load, validate and reconcile a TOML configuration (no transport)"
    )
    chk.add_argument("config", metavar="CONFIG.toml")
    cap = sub.add_parser("capabilities", help="print the capability report as JSON")
    cap.add_argument("--nuclear", action="store_true", help="the report for physics.nuclear = true")
    ver = sub.add_parser("verify", help="verify a result bundle (tier 1; --rederive adds tier 2)")
    ver.add_argument("bundle", metavar="DIR")
    ver.add_argument(
        "--rederive", action="store_true", help="tier 2: re-derive the effective block"
    )
    ver.add_argument("--cache-dir", default=None, help="data cache (required with --rederive)")
    ver.add_argument(
        "--inputs", action="append", default=[], metavar="ROLE=PATH", help="trusted input file"
    )
    rep = sub.add_parser("report", help="print the regenerated report of a bundle")
    rep.add_argument("bundle", metavar="DIR")
    rep.add_argument(
        "--check", action="store_true", help="compare the bytes on disk, write nothing"
    )
    return parser


def _build_nuclear(args: argparse.Namespace) -> int:
    from ionmc.data import cache
    from ionmc.nuclear.build import BuildError, BuildOptions, build_nuclear_proton
    from ionmc.nuclear.tables import qualification_failures

    opts = BuildOptions(
        points_per_decade=args.points_per_decade,
        diagnostic_events=args.diagnostic_events,
        strict=args.strict,
    )
    try:
        res = build_nuclear_proton(cache.resolve_cache_dir(args.cache_dir), opts, log=print)
    except (BuildError, FileNotFoundError, cache.IntegrityError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    gate = res.info["gate_d6"]
    mult = res.info["multiplicity"]
    pmin = mult["p_accept_min"]
    print(f"table id {res.table_id}")
    print(f"build seconds {res.info['timing_s']:.1f}")
    print(f"grid size {res.info['grid']['n_points']}")
    print(f"P_accept min {pmin['p_accept']:.4f} ({pmin['target']}, {pmin['e_mev']:.4g} MeV)")
    print(f"non-converged lambda nodes {len(mult['non_converged_nodes'])}")
    print(f"transport_energy_bound_mev {res.info['transport_energy_bound_mev']:.2f}")
    print(f"history_energy_bound_mev {res.info['history_energy_bound_mev']:.2f}")
    print(f"recoil_t_max_mev {res.info['recoil_t_max_mev']:.4f}")
    print(f"transport_path_bound_terms {res.info['transport_path_bound_terms']}")
    for key, num in gate["numbers"].items():
        print(
            f"D6 {key} MeV: G={num['G']:.4f} D={num['D']:.3e} "
            f"R99.9={num['percentile_range_g_cm2']:.4f} g/cm2"
        )
    print(
        f"tier1_pass={gate['tier1_pass']} tier2_pass={gate['tier2_pass']} "
        f"ceiling_pass={gate['ceiling_pass']}"
    )
    reasons = qualification_failures(res.info)
    print("qualified for transport: " + ("yes" if not reasons else "NO: " + "; ".join(reasons)))
    return 0


def _run_data(args: argparse.Namespace) -> int:
    """Execute ``ionmc data <action>``."""
    from ionmc.data import cache
    from ionmc.data.acquire import OfflineError, fetch, import_file
    from ionmc.data.registry import DATASETS

    action = args.data_command
    if action is None:
        return 2
    if action == "build":
        return _build_nuclear(args)
    cdir = cache.resolve_cache_dir(args.cache_dir)
    if action == "list":
        for ds in DATASETS.values():
            cached = cache.read_manifest(ds.id, cdir) is not None
            print(f"{ds.id}\t{'cached' if cached else 'missing'}\t{ds.description}")
        return 0
    try:
        if action == "import":
            print(import_file(args.path, args.dataset_id, cdir))
        elif action == "fetch":
            print(fetch(args.dataset_id, cdir, offline=args.offline))
        elif action == "verify":
            print(f"OK {cache.verify(args.dataset_id, cdir)}")
        else:  # path
            print(fetch(args.dataset_id, cdir, offline=True))
    except (OfflineError, FileNotFoundError, cache.IntegrityError, KeyError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


# Exception -> exit code (decision 0045 section 7, frozen; the first matching row wins, so
# subclasses precede their bases). The "run completed with valid = false" row (exit 3, bundle
# written) and success (0) are not exceptions and are handled in ``_cmd_run``. ``_UsageError`` is
# the pre-check form of row 1 (argparse usage error); ``VerificationError`` is a tier-2 failure of
# ``verify --rederive`` (amendment A-2) and shares row 8's code.


class _UsageError(Exception):
    """A command-line misuse found after parsing (for example ``run`` without an output path)."""


def _cache_cause(exc: BaseException) -> bool:
    from ionmc.data.acquire import OfflineError
    from ionmc.data.cache import IntegrityError

    seen: set[int] = set()
    cur: BaseException | None = exc
    while cur is not None and id(cur) not in seen:
        seen.add(id(cur))
        if isinstance(cur, IntegrityError | OfflineError | FileNotFoundError):
            return True
        cur = cur.__cause__
    return False


def _exit_code(exc: BaseException) -> int:
    """The frozen exception -> exit-code table of decision 0045 section 7 (first match wins)."""
    from ionmc.configio import ConfigError, DatasetError
    from ionmc.errors import (
        BackendUnavailableError,
        TransportLimitError,
        UnsupportedCombinationError,
    )
    from ionmc.results import ResultFormatError, VerificationError

    table: tuple[tuple[Any, int], ...] = (
        (lambda e: isinstance(e, _UsageError), 2),  # 1: usage error
        (lambda e: isinstance(e, ConfigError), 2),  # 2: before its base class
        (lambda e: isinstance(e, DatasetError), 1),  # 3
        (lambda e: isinstance(e, UnsupportedCombinationError) and _cache_cause(e), 1),  # 3
        (lambda e: isinstance(e, UnsupportedCombinationError), 4),  # 4
        (lambda e: isinstance(e, BackendUnavailableError), 5),  # 5
        (lambda e: isinstance(e, TransportLimitError), 3),  # 6
        (lambda e: isinstance(e, ResultFormatError | VerificationError), 6),  # 8
    )
    for matches, code in table:
        if matches(exc):
            return int(code)
    return 70  # 9: any other exception (a defect)


def _canonical(obj: Any) -> str:
    from ionmc.configio import canonical_json

    return canonical_json(obj)


def _cmd_run(args: argparse.Namespace, out_stream: TextIO) -> int:
    from pathlib import Path

    from ionmc.configio import load_toml_full, verify_datasets
    from ionmc.results import save
    from ionmc.simulation import Simulation

    loaded = load_toml_full(args.config)
    out = args.out if args.out is not None else loaded.output.path
    if out is None:
        raise _UsageError("no output directory: pass --out DIR or set [output].path")
    dest = Path(out)
    if dest.parent.exists() and not dest.parent.is_dir():
        raise _UsageError(f"parent of {dest.name!r} is not a directory")
    if dest.exists() and not args.overwrite:
        raise _UsageError(f"{dest.name!r} exists (pass --overwrite)")
    verify_datasets(loaded.config)
    sim = Simulation(loaded.config)  # validate and reconcile; raises before any transport
    result = sim.run()
    dest.parent.mkdir(parents=True, exist_ok=True)  # nothing is created before this point
    bundle = save(
        result,
        dest,
        manifest=loaded.manifest,
        persist_batches=loaded.output.persist_batches,
        overwrite=args.overwrite,
    )
    print(f"bundle {bundle}", file=out_stream)
    if not result.valid:
        print(
            "run completed with valid = false (allow_invalid_result); see the banner in report.md",
            file=sys.stderr,
        )
        return 3
    return 0


def _cmd_check(args: argparse.Namespace, out_stream: TextIO) -> int:
    from ionmc.configio import load_toml_full, verify_datasets
    from ionmc.simulation import Simulation

    loaded = load_toml_full(args.config)
    verify_datasets(loaded.config)
    sim = Simulation(loaded.config)
    print(_canonical(sim.effective.summary()), file=out_stream)
    return 0


def _cmd_capabilities(args: argparse.Namespace, out_stream: TextIO) -> int:
    from ionmc import capabilities

    print(_canonical(capabilities(args.nuclear)), file=out_stream)
    return 0


def _cmd_verify(args: argparse.Namespace, out_stream: TextIO) -> int:
    import dataclasses

    from ionmc.results import load, verify_effective

    if args.rederive and args.cache_dir is None:
        raise _UsageError("--rederive needs --cache-dir")
    if (args.cache_dir is not None or args.inputs) and not args.rederive:
        raise _UsageError("--cache-dir and --inputs apply to --rederive only")
    files: dict[str, str] = {}
    for item in args.inputs:
        role, sep, path = item.partition("=")
        if not sep or not role or not path:
            raise _UsageError(f"--inputs expects ROLE=PATH, got {item!r}")
        files[role] = path
    loaded = load(args.bundle)
    out: dict[str, Any] = {"tier1": "ok", "binding_limits": list(loaded.binding_limits)}
    code = 0
    if args.rederive:
        rep = verify_effective(loaded, cache_dir=args.cache_dir, input_files=files or None)
        out["tier2"] = dataclasses.asdict(rep)
        code = 0 if rep.ok else 6
    print(_canonical(out), file=out_stream)
    return code


def _cmd_report(args: argparse.Namespace, out_stream: TextIO) -> int:
    from pathlib import Path

    from ionmc.report import regenerate_report
    from ionmc.results import ResultFormatError

    files = regenerate_report(args.bundle)
    if not args.check:
        out_stream.write(files["report.md"].decode("utf-8"))
        return 0
    bad = []
    for name, data in sorted(files.items()):
        try:
            on_disk = (Path(args.bundle) / name).read_bytes()
        except OSError:
            on_disk = None
        if on_disk != data:
            bad.append(name)
    if bad:
        raise ResultFormatError(f"derived files differ from their regeneration: {bad}")
    print(f"report ok ({len(files)} files)", file=out_stream)
    return 0


_COMMANDS = {
    "run": _cmd_run,
    "check": _cmd_check,
    "capabilities": _cmd_capabilities,
    "verify": _cmd_verify,
    "report": _cmd_report,
}


def _dispatch(args: argparse.Namespace) -> int:
    """Run a V3-006 subcommand; every failure maps through :func:`_exit_code`."""
    out_stream = sys.stdout
    try:
        # Warp prints its banner and kernel-load messages on stdout: keep stdout for the result
        # (JSON, report text) by sending everything else to stderr.
        with contextlib.redirect_stdout(sys.stderr):
            return _COMMANDS[args.command](args, out_stream)
    except Exception as exc:
        code = _exit_code(exc)
        if code == 70:
            traceback.print_exc()
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return code


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
    if args.command == "notices":
        from importlib.resources import files

        print((files("ionmc") / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8"))
        return 0
    if args.command == "data":
        code = _run_data(args)
        if code == 2:
            parser.print_usage(sys.stderr)
        return code
    if args.command in _COMMANDS:
        return _dispatch(args)
    parser.print_usage(sys.stderr)
    return 2
