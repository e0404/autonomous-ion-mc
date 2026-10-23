"""Command-line interface for IonMC."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

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
    bld.add_argument("table", choices=["nuclear-proton", "elastic-proton"])
    bld.add_argument("--cache-dir", default=None, help="cache directory")
    bld.add_argument("--points-per-decade", type=int, default=50)
    bld.add_argument("--diagnostic-events", type=int, default=20000)
    bld.add_argument(
        "--model",
        choices=["pinned", "m2", "m2prime", "m3"],
        default="pinned",
        help="elastic-proton only: m2 = exploratory optical-limit Glauber shapes, "
        "m2prime = m2 with the Bertulani-De Conti Pauli-blocked in-medium sigma_NN, "
        "m3 = ENDF/B-VIII.0 neutron MF4/MT2 mirror (C-12, O-16)",
    )
    bld.add_argument(
        "--pauli-variant",
        choices=["clementel-villi", "bdc-pn"],
        default=None,
        help="m2prime: reading of rho1, rho2 in the in-medium factor (required for m2prime)",
    )
    bld.add_argument(
        "--pauli-emin-mev",
        type=float,
        default=None,
        help="m2prime: nodes below this lab energy keep the black-disk node (required)",
    )
    bld.add_argument(
        "--measured-sigma-tot",
        default=None,
        help="m2/m2prime: JSON table of measured sigma_tot for the report-only comparison",
    )
    bld.add_argument("--endf-dir", default=None, help="neutron ENDF files (--model m3)")
    bld.add_argument(
        "--above-emax",
        choices=["error", "hold-last", "black-disk"],
        default="error",
        help="m3: nodes above the end of the ENDF evaluation (20 MeV)",
    )
    bld.add_argument("--beta", choices=["unitarity", "horiuchi-table"], default="unitarity")
    bld.add_argument(
        "--pp-level", choices=["unchanged", "said-anchored-ge150"], default="unchanged"
    )
    bld.add_argument("--unfold", choices=["formfactor", "rms-matched"], default="rms-matched")
    bld.add_argument("--said-dir", default=None, help="SAID SGT/ALFA text outputs (--model m2)")
    bld.add_argument(
        "--dsg-dir",
        default=None,
        help="SAID DSG segments + dsg_manifest.json (--pp-level said-anchored-ge150)",
    )
    bld.add_argument("--out-dir", default=None, help="output directory (default cache/derived)")
    bld.add_argument(
        "--allow-unpinned-inputs",
        action="store_true",
        help="m2: allow empty file_sha256 pins (marks exploratory.inputs_unpinned)",
    )
    bld.add_argument(
        "--strict", action="store_true", help="fail when a lambda node does not converge to 1e-3"
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


def _build_elastic(args: argparse.Namespace) -> int:
    from ionmc.data import cache
    from ionmc.nuclear import elastic_glauber as glauber
    from ionmc.nuclear import elastic_m3
    from ionmc.nuclear.build import BuildError
    from ionmc.nuclear.elastic_build import ElasticBuildOptions, build_elastic_proton

    opts = ElasticBuildOptions(points_per_decade=args.points_per_decade)
    m2 = None
    try:
        if args.model in ("m2", "m2prime"):
            if not args.said_dir:
                raise BuildError(f"--model {args.model} needs --said-dir")
            m2 = glauber.M2Config(
                Path(args.said_dir),
                dsg_files=glauber.load_dsg_dir(Path(args.dsg_dir)) if args.dsg_dir else None,
                beta=args.beta,
                pp_level=args.pp_level,
                unfold=args.unfold,
                allow_unpinned_inputs=args.allow_unpinned_inputs,
                pauli=args.model == "m2prime",
                pauli_variant=args.pauli_variant,
                pauli_emin_mev=args.pauli_emin_mev,
                measured_sigma_tot=(
                    Path(args.measured_sigma_tot) if args.measured_sigma_tot else None
                ),
            )
        m3 = None
        if args.model == "m3":
            if not args.endf_dir:
                raise BuildError("--model m3 needs --endf-dir")
            m3 = elastic_m3.M3Config(
                Path(args.endf_dir),
                above_emax=args.above_emax,
                allow_unpinned_inputs=args.allow_unpinned_inputs,
            )
        res = build_elastic_proton(
            cache.resolve_cache_dir(args.cache_dir),
            opts,
            log=print,
            m2=m2,
            m3=m3,
            out_dir=Path(args.out_dir) if args.out_dir else None,
        )
    except (
        BuildError,
        glauber.M2Error,
        elastic_m3.M3Error,
        FileNotFoundError,
        cache.IntegrityError,
    ) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    info = res.info
    print(f"table id {res.table_id}")
    print(f"npz sha256 {info['npz_sha256']}")
    print(f"build seconds {info['timing_s']:.1f}")
    print(f"grid size {info['grid']['n_points']}")
    for t in info["targets"]:
        key = "e_min_pp_mev" if t["kind"] == "law5-ltp1" else "e_min_shape_mev"
        print(f"{t['name']}: {key} {t[key]:.4g}")
    print(f"O-16 MT2 is a C-12 copy: {info['o16_finding']['o16_mt2_is_c12_copy']}")
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
        return _build_elastic(args) if args.table == "elastic-proton" else _build_nuclear(args)
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
    parser.print_usage(sys.stderr)
    return 2
