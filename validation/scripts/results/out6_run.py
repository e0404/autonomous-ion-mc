"""The ``out6`` workflow driver of plan v3-006 (step S8): LV/HR evidence runs and their checks.

Subcommands: ``run`` (config, ``ionmc run``, clean-process ``verify`` and ``report --check``,
summary), ``compare`` (CLI-5), ``identities`` (CLI-6), ``reproduce`` (OUT-3), ``prov7`` (PROV-7)
and ``host`` (the controlled host-runner entry). Each check prints one line
``RESULT <row> <PASS|FAIL> <json path>``; the host stdout is the only channel readable afterwards.

Every ``ionmc`` call is a child process, one at a time, with the thread limits at 1; the driver
reads no environment variable itself (children inherit the environment, with
``IONMC_CACHE_DIR`` set from ``--cache-dir``). The summaries carry identities and digests only:
no table value, no array and no path of the path set P (decision 0045 section 2). Not a clinical
claim; synthetic tiny runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import out6_cases  # noqa: E402

# the in-process equivalent of the ``ionmc`` console script (the package has no __main__)
MAIN = "import sys; from ionmc.cli import main; sys.exit(main(sys.argv[1:]))"
CLEAN_MAIN = "import sys; sys.path.insert(0, sys.argv[1]); from ionmc.cli import main; sys.exit(main(sys.argv[2:]))"  # noqa: E501
THREADS = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "IONMC_SINGLE_PROCESS": "1",
    "PYTHONFAULTHANDLER": "1",
}
#: a path of the set P as text (the same pattern as the DOC-1 page hygiene test)
P_TEXT = re.compile(r"(?<![\w.])(?:/home/|/Users/|/root/|/tmp/|/mnt/|/workspace)\S*|[A-Za-z]:\\\S*")
DATASET_ICRU90 = "geant4-icru90-stopping-11.4.2"
#: the asserted-identity tests of row REG-1 (the list of tests/ionmc/test_v3_006_regression.py)
REG1_TESTS = (
    "tests/ionmc/test_run_suite.py::test_a16_has_no_recorded_exception_and_the_baseline_is_the_v3_003d_merge",  # noqa: E501
    "tests/ionmc/test_run_suite.py::test_a16_gate_accepts_the_intended_change_and_rejects_everything_else",  # noqa: E501
    "tests/ionmc/test_scoring_transport.py::test_a16_tallies_do_not_change_the_qualified_outputs",
    "tests/ionmc/test_stopping_models.py::test_id_def_default_table_identity_equals_the_base_commit",
)
TRACE_TEST = "tests/ionmc/test_transport_diagnostic.py::test_default_trace_is_unchanged_from_the_stored_baseline"  # noqa: E501


# ---------------------------------------------------------------------------------------------
# helpers


#: set by ``cmd_host``: the cache-backed tests must run, never skip (IONMC_REQUIRE_DATA, LV tier)
_EXTRA_ENV: dict[str, str] = {}


def _env(cache_dir: str | Path) -> dict[str, str]:
    env = dict(os.environ)
    env.update(THREADS)
    env["IONMC_CACHE_DIR"] = str(cache_dir)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env.update(_EXTRA_ENV)
    return env


def _package_parent() -> str:
    import ionmc

    return str(Path(ionmc.__file__).resolve().parent.parent)


def _child(
    argv: Sequence[str],
    cwd: Path,
    log_dir: Path,
    cache_dir: str | Path,
    timeout: float,
    log: str | None = None,
) -> tuple[int, float]:
    """Run one child to completion (one at a time). Its output goes to ``log_dir/logs``, the
    writable ``--out`` directory, never to ``cwd`` (which may be a read-only snapshot)."""
    logs = log_dir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    label = re.sub(r"[^\w.-]+", "_", "_".join(str(a) for a in argv[3:6]))[:60]
    t0 = time.monotonic()
    with open(logs / (log or f"{label}.log"), "ab") as fh:
        proc = subprocess.run(  # noqa: S603 - argv, no shell
            list(argv), cwd=cwd, env=_env(cache_dir), stdout=fh, stderr=subprocess.STDOUT,
            timeout=timeout, check=False,
        )  # fmt: skip
    return proc.returncode, time.monotonic() - t0


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _dump(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=1, sort_keys=True) + "\n", encoding="utf-8")


def _result(row: str, ok: bool, path: Path) -> None:
    print(f"RESULT {row} {'PASS' if ok else 'FAIL'} {path}", flush=True)


def _flatten(obj: Any, ptr: str = "") -> dict[str, Any]:
    if isinstance(obj, Mapping):
        out: dict[str, Any] = {}
        for k, v in obj.items():
            out.update(_flatten(v, f"{ptr}/{k}"))
        return out or {ptr: {}}
    if isinstance(obj, list):
        out = {}
        for i, v in enumerate(obj):
            out.update(_flatten(v, f"{ptr}/{i}"))
        return out or {ptr: []}
    return {ptr: obj}


def differing_pointers(a: Any, b: Any) -> list[str]:
    fa, fb = _flatten(a), _flatten(b)
    return sorted(p for p in set(fa) | set(fb) if fa.get(p, _MISSING) != fb.get(p, _MISSING))


_MISSING = object()


def _strip_cache_dir(obj: Any) -> Any:
    if isinstance(obj, Mapping):
        return {k: _strip_cache_dir(v) for k, v in obj.items() if k != "cache_dir"}
    if isinstance(obj, list):
        return [_strip_cache_dir(v) for v in obj]
    return obj


NUCLEAR_IDENTITY_KEYS = (
    "table_id", "npz_sha256", "schema", "builder_version", "source_sha256",
    "deuteron_tables_sha256", "grid_size", "elements",
)  # fmt: skip


def _nuclear_identity(eff: Mapping[str, Any]) -> dict[str, Any] | None:
    """The identity fields of ``effective.nuclear`` (no table value, no cache path)."""
    nuc = eff.get("nuclear")
    if nuc is None:
        return None
    return {key: nuc[key] for key in NUCLEAR_IDENTITY_KEYS if key in nuc}


def _longest_numeric_list(obj: Any) -> int:
    if isinstance(obj, list):
        n = len(obj) if obj and all(isinstance(x, int | float) for x in obj) else 0
        return max([n, *(_longest_numeric_list(x) for x in obj)])
    if isinstance(obj, Mapping):
        return max([0, *(_longest_numeric_list(v) for v in obj.values())])
    return 0


def check_summary_text(text: str, p_values: Sequence[str]) -> None:
    """Raise when a summary holds a path of P (pattern or one of the actual values)."""
    if P_TEXT.search(text):
        raise AssertionError(f"a path pattern of P in the summary: {P_TEXT.search(text)}")
    for v in p_values:
        if v and len(v) > 3 and v in text:
            raise AssertionError("a value of the path set P in the summary")
    if len(text.encode("utf-8")) >= 64 * 1024:
        raise AssertionError("summary larger than 64 KiB")


def _bundle_of(summary: Path) -> Path:
    return summary.with_name(summary.name.removesuffix(".summary.json"))


def _load_doc(bundle: Path) -> dict[str, Any]:
    return json.loads((bundle / "result.json").read_text(encoding="utf-8"))  # type: ignore[no-any-return]


# ---------------------------------------------------------------------------------------------
# run


def bundle_facts(loaded: Any) -> dict[str, Any]:
    """Every summary field that is derived from the verified bundle (``results.LoadedResult``).

    ``build_summary`` writes these fields and ``summary_binding_failures`` regenerates them from
    the sibling bundle, so a summary cannot claim more than its bundle holds."""
    doc = loaded.doc
    eff = doc["effective"]
    run = doc["run"]
    out_path = run["output_path"]
    if os.path.isabs(out_path) or ".." in Path(out_path).parts:
        raise AssertionError("the persisted output_path is not relative")
    stopping = _strip_cache_dir(eff["stopping"])
    nuclear = _nuclear_identity(eff)
    for name, block in (("stopping", stopping), ("nuclear", nuclear)):
        if _longest_numeric_list(block) > 8:
            raise AssertionError(f"effective.{name} holds a numeric table")
    return {
        "seed": run["seed"],
        "backend": run["backend"],
        "device": run["device"],
        "precision": run["precision"],
        "production": run["production"],
        "n_histories": run["n_histories"],
        "n_batches": run["n_batches"],
        "valid": doc["valid"],
        "output_path_persisted": out_path,
        "result_json_sha256": _sha256_file(Path(loaded.path) / "result.json"),
        "files": doc["files"],
        "stopping": stopping,
        "nuclear": nuclear,
        "requested_nuclear_table_id": doc["requested"]["physics"].get("nuclear_table_id"),
        "code": doc["provenance"]["code"],
        "array_content_sha256": {e["name"]: e["content_sha256"] for e in doc["arrays"]},
    }


def build_summary(
    bundle: Path, r: int, k: int, codes: Mapping[str, int], durations: Mapping[str, float]
) -> dict[str, Any]:
    from ionmc import results

    loaded = results.load(bundle)  # tier 1: every digest recomputed
    return {
        "schema": "ionmc.out6-summary/1",
        "r_index": r,
        "k": k,
        "exit_codes": dict(codes),
        "durations_s": {n: round(v, 3) for n, v in durations.items()},
        **bundle_facts(loaded),
    }


def run_case(
    r: int,
    k: int,
    out: Path,
    cache_dir: Path,
    *,
    backend: str | None = None,
    precision: str | None = None,
    timeout: float = 1800.0,
) -> Path:
    """One case: config, ``ionmc run``, clean-process ``verify`` and ``report --check``, summary."""
    out.mkdir(parents=True, exist_ok=True)
    name = f"r{r}k{k}"
    if (out / name).exists() or (out / f"{name}.summary.json").exists():
        raise SystemExit(f"{name} exists under --out: every run writes to a fresh path")
    cfg = out_cases_write(r, k, out / "cfg", cache_dir, backend, precision, name)
    py = sys.executable
    codes: dict[str, int] = {}
    dur: dict[str, float] = {}
    codes["run"], dur["run"] = _child(
        [py, "-c", MAIN, "run", cfg, "--out", name], out, out, cache_dir, timeout
    )
    pkg = _package_parent()
    clean = [py, "-I", "-c", CLEAN_MAIN, pkg]
    codes["verify"], dur["verify"] = _child([*clean, "verify", name], out, out, cache_dir, timeout)
    codes["report_check"], dur["report_check"] = _child(
        [*clean, "report", name, "--check"], out, out, cache_dir, timeout
    )
    summary = build_summary(out / name, r, k, codes, dur)
    text = json.dumps(summary, indent=1, sort_keys=True) + "\n"
    check_summary_text(text, [str(out), str(cache_dir), str(Path.home()), pkg, sys.prefix])
    path = out / f"{name}.summary.json"
    path.write_text(text, encoding="utf-8")
    ok = all(c == 0 for c in codes.values())
    _result(f"RUN-r{r}", ok, path)
    return path


def out_cases_write(
    r: int, k: int, cfg_dir: Path, cache_dir: Path, backend: str | None, precision: str | None,
    name: str,
) -> str:  # fmt: skip
    """Write the config; the file is named for the case (``name``) and returned relative to out."""
    p = out6_cases.write_config(r, k, cfg_dir, cache_dir, backend=backend, precision=precision)
    target = cfg_dir / f"{name}.toml"
    if target != p:
        p.replace(target)
    return f"cfg/{name}.toml"


# ---------------------------------------------------------------------------------------------
# compare (CLI-5)

#: pointers that may differ between the backends/precisions of one case (CLI-5): the backend,
#: device, precision and production fields, what a backend changes without being a physics input
#: (timings, environment, the configuration file's digest), the seed (the out6 seed depends on
#: r_index by construction, checked against the formula), the path set P, the array contents and
#: everything derived from them (digests, file hashes, accounting counters, statistics, plots)
_ALLOWED = tuple(
    re.compile(p)
    for p in (
        r"^/run/(backend|device|precision|production|seed|cpu_workers|chunk_histories)$",
        r"^/run/timings(/.*)?$",
        r"^/run/(output_path|cache_dir|cache_dir_source)$",
        r"^/requested/run/(backend|precision|seed)$",
        r"^/effective/(backend|precision|production|seed|cpu_workers|chunk_histories)$",
        r"^/provenance/environment(/.*)?$",
        r"^/input_manifest/(base_dir|config_file|config_sha256)(/.*)?$",
        r"^/(run|effective)/rng/u01$",  # the uniform mapping follows the precision (24/23 bit)
        r"^/effective_extra/resolutions/(run|chunk)/(class|rule)$",  # the rows of those fields
        r"^/arrays/\d+/(content_sha256|dtype|sha256|nbytes|bytes)$",
        r"^/files(/.*)?$",
        r"^/accounting(/.*)?$",
        r"^/statistics(/.*)?$",
        r"^/effective/(stopping|nuclear)(/.*)?/cache_dir$",
    )
)


def _by_row(doc: Mapping[str, Any]) -> dict[str, Any]:
    """``doc`` with the reconciliation list keyed by its row id (a stable pointer per row)."""
    out = dict(doc)
    extra = dict(doc["effective_extra"])
    extra["resolutions"] = {r["row"]: r for r in extra["resolutions"]}
    out["effective_extra"] = extra
    return out


def compare_docs(
    docs: Sequence[dict[str, Any]], seeds_by_r: Sequence[tuple[int, int, int]]
) -> dict[str, Any]:
    """Pointers differing between ``docs[0]`` and each later document; which are not allowed."""
    for r, k, seed in seeds_by_r:
        if seed != out6_cases.seed_of(r, k):
            return {"verdict": "FAIL", "reason": f"seed of r{r}k{k} is not the out6 value"}
    pairs = []
    bad_total = 0
    for i, doc in enumerate(docs[1:], start=1):
        diff = differing_pointers(_by_row(docs[0]), _by_row(doc))
        bad = [p for p in diff if not any(a.match(p) for a in _ALLOWED)]
        bad_total += len(bad)
        pairs.append({"pair": [0, i], "differing_pointers": diff, "not_allowed": bad})
    return {"verdict": "FAIL" if bad_total else "PASS", "pairs": pairs}


#: the case identity of CLI-5: r_index -> (backend, precision)
CLI5_CASES = {0: ("python", "float64"), 1: ("warp-cpu", "float32"), 2: ("warp-cuda", "float32")}


def compare_identity_failures(sms: Sequence[Mapping[str, Any]], tier: str) -> list[str]:
    """The violated conditions of the case identity (empty when the summaries are the cases)."""
    want = [0, 1, 2] if tier == "hr" else [0, 1]
    bad: list[str] = []
    got = sorted(int(sm["r_index"]) for sm in sms)
    if tier not in ("hr", "lv"):
        return [f"unknown tier {tier!r}"]
    if got != want:
        bad.append(f"r_index set {got} is not exactly {want} for tier {tier}")
    if len({sm["k"] for sm in sms}) != 1:
        bad.append("the summaries do not share one k")
    if len({sm["result_json_sha256"] for sm in sms}) != len(sms):
        bad.append("result.json digests are not pairwise different")
    for sm in sms:
        exp = CLI5_CASES.get(int(sm["r_index"]))
        if exp is not None and (sm["backend"], sm["precision"]) != exp:
            bad.append(
                f"r{sm['r_index']}: backend/precision {(sm['backend'], sm['precision'])} "
                f"is not {exp}"
            )
    return bad


def _fact_mismatches(sm: Mapping[str, Any], loaded: Any) -> list[str]:
    """Every bundle-derived field of the summary that differs from the regenerated value."""
    return [
        f"summary {key} differs from the bundle's {key}"
        for key, value in bundle_facts(loaded).items()
        if sm.get(key, _MISSING) != value
    ]


def summary_binding_failures(sm: Mapping[str, Any], bundle: Path) -> list[str]:
    """Conditions under which a summary does not describe the verified bundle next to it."""
    from ionmc import results

    bad: list[str] = []
    codes = sm.get("exit_codes")
    if not isinstance(codes, Mapping) or set(codes) != {"run", "verify", "report_check"}:
        bad.append("exit_codes are missing or incomplete")
    elif any(c != 0 for c in codes.values()):
        bad.append(f"an exit code is not 0: {dict(codes)}")
    try:
        loaded = results.load(bundle)  # tier 1: every digest of the bundle is recomputed
    except Exception as exc:  # noqa: BLE001 - any failure of the reader is a violated condition
        return [*bad, f"the bundle does not load: {type(exc).__name__}"]
    doc = loaded.doc
    for field in ("backend", "precision", "seed"):
        if doc["requested"]["run"].get(field) != sm.get(field):
            bad.append(f"summary {field} differs from requested.run.{field} of the bundle")
    bad += _fact_mismatches(sm, loaded)
    return bad


def cmd_compare(summaries: Sequence[Path], out: Path | None = None, tier: str = "hr") -> Path:
    sms = [json.loads(s.read_text(encoding="utf-8")) for s in summaries]
    ident = compare_identity_failures(sms, tier)
    docs: list[dict[str, Any]] = []
    for s, sm in zip(summaries, sms, strict=True):
        bound = summary_binding_failures(sm, _bundle_of(s))
        ident += [f"{s.name}: {f}" for f in bound]
        docs.append(
            {}
            if bound and not (_bundle_of(s) / "result.json").is_file()
            else _load_doc(_bundle_of(s))
        )
    seeds = [(sm["r_index"], sm["k"], sm["seed"]) for sm in sms]
    if ident:
        verdict: dict[str, Any] = {"verdict": "FAIL", "identity_failures": ident}
    else:
        verdict = compare_docs(docs, seeds)
    for pair in verdict.get("pairs", []):
        print("differing pointers", pair["pair"], len(pair["differing_pointers"]))
        for p in pair["differing_pointers"]:
            print("  ", p)
    path = out or summaries[0].with_name("compare.json")
    _dump(path, verdict)
    _result("CLI-5-LV" if tier == "lv" else "CLI-5", verdict["verdict"] == "PASS", path)
    return path


# ---------------------------------------------------------------------------------------------
# identities (CLI-6)


def identity_failures(
    summary: Mapping[str, Any], doc: Mapping[str, Any], expect: str, cache_dir: Path | None
) -> list[str]:
    bad: list[str] = []
    eff = doc["effective"]
    stopping = _strip_cache_dir(eff["stopping"])
    if summary["stopping"] != stopping:
        bad.append("summary stopping block differs from the bundle's effective.stopping")
    if summary["nuclear"] != _nuclear_identity(eff):
        bad.append("summary nuclear block differs from the bundle's effective.nuclear")
    if expect == "carbon":
        if eff["source"].get("projectile") != "carbon12":
            bad.append(f"projectile {eff['source'].get('projectile')!r} is not carbon12")
        if doc["requested"]["source"].get("projectile") != "carbon12":
            bad.append("requested projectile is not carbon12")
    elif expect == "icru90":
        if stopping.get("model_id_effective") != "icru90-water-v1":
            bad.append(f"stopping model {stopping.get('model_id_effective')!r}")
        pins = [u["pin"] for u in stopping.get("data_used", [])]
        if not any(p.get("dataset_id") == DATASET_ICRU90 for p in pins):
            bad.append("the ICRU 90 dataset is not among data_used")
        from ionmc.data.registry import DATASETS

        reg = DATASETS.get(DATASET_ICRU90)
        for p in pins:
            if (
                p.get("dataset_id") == DATASET_ICRU90
                and reg is not None
                and p.get("archive_sha256") != reg.sha256
            ):
                bad.append("archive_sha256 differs from the registry pin")
        if not re.fullmatch(r"[0-9a-f]{64}", str(stopping.get("content_sha256"))):
            bad.append("content_sha256 is not a digest")
        if not re.fullmatch(r"[0-9a-f]{64}", str(stopping.get("spec_sha256"))):
            bad.append("spec_sha256 is not a digest")
    elif expect == "nuclear":
        nuc = _nuclear_identity(eff) or {}
        tid = nuc.get("table_id")
        if not tid or tid != summary["requested_nuclear_table_id"]:
            bad.append("nuclear table_id is missing or differs from the requested id")
        sha = nuc.get("npz_sha256")
        if not re.fullmatch(r"[0-9a-f]{64}", str(sha)):
            bad.append("npz_sha256 is missing")
        elif cache_dir is not None and tid:
            f = cache_dir / "derived" / f"nuclear-proton-{tid}.npz"
            if not f.is_file() or _sha256_file(f) != sha:
                bad.append("npz_sha256 differs from the cached table file")
    else:
        bad.append(f"unknown expectation {expect!r}")
    return bad


def cmd_identities(
    summary: Path, expect: str, cache_dir: Path | None, out: Path | None = None
) -> Path:
    sm = json.loads(summary.read_text(encoding="utf-8"))
    bad = identity_failures(sm, _load_doc(_bundle_of(summary)), expect, cache_dir)
    path = out or summary.with_name(f"identities-{expect}.json")
    _dump(path, {"expect": expect, "verdict": "FAIL" if bad else "PASS", "failures": bad})
    _result(f"CLI-6-{expect}", not bad, path)
    return path


# ---------------------------------------------------------------------------------------------
# reproduce (OUT-3)


def _reproduce_child(bundle: Path, out: Path, cache_dir: Path) -> int:
    """In the child: ``from_dict(requested)``, the same seed, a fresh bundle (API, not a TOML)."""
    from ionmc import configio, results
    from ionmc.simulation import Simulation

    doc = _load_doc(bundle)
    cfg = configio.from_dict(doc["requested"])
    persisted = doc["statistics"]["batch_arrays"]["status"] == "persisted"
    results.save(Simulation(cfg).run(), out, persist_batches=persisted)
    return 0


def reproduction_verdict(orig: Mapping[str, Any], again: Mapping[str, Any]) -> dict[str, Any]:
    """OUT-3: every array ``content_sha256`` and the ``effective`` block (pointer-wise) agree."""
    da = {e["name"]: e["content_sha256"] for e in orig["arrays"]}
    db = {e["name"]: e["content_sha256"] for e in again["arrays"]}
    eff = differing_pointers(orig["effective"], again["effective"])
    eff = [p for p in eff if not p.endswith("/cache_dir")]
    equal = da == db and bool(da)
    same_requested = orig["requested"] == again["requested"]
    return {
        "arrays_equal": equal,
        "differing_array_names": sorted(n for n in set(da) | set(db) if da.get(n) != db.get(n)),
        "effective_differing_pointers": eff,
        "requested_equal": same_requested,
        "ok": equal and not eff and same_requested,
    }


def cmd_reproduce(bundle: Path, out: Path, cache_dir: Path, timeout: float = 1800.0) -> Path:
    from ionmc import results

    out.mkdir(parents=True, exist_ok=True)
    name = f"{bundle.name}-repro"
    if (out / name).exists():
        raise SystemExit(f"{name} exists under --out: every run writes to a fresh path")
    # the child runs with cwd = out (its bundle path is then the relative name, so the persisted
    # output_path carries no absolute path): the ORIGINAL bundle must be passed resolved
    code, dur = _child(
        [sys.executable, str(Path(__file__).resolve()), "_reproduce-child", "--bundle",
         str(bundle.resolve()), "--out", name, "--cache-dir", str(cache_dir)],
        out, out, cache_dir, timeout,
    )  # fmt: skip
    verdict: dict[str, Any] = {"child_exit_code": code, "method": "configio.from_dict + API"}
    if code == 0:
        a, b = results.load(bundle), results.load(out / name)
        verdict.update(reproduction_verdict(a.doc, b.doc))
        ok = verdict["ok"]
    else:
        ok = False
    verdict.pop("ok", None)
    verdict["verdict"] = "PASS" if ok else "FAIL"
    verdict["duration_s"] = round(dur, 3)
    path = out / f"reproduce-{bundle.name}.json"
    _dump(path, verdict)
    _result("OUT-3", ok, path)
    return path


# ---------------------------------------------------------------------------------------------
# prov7


def prov7_failures(code: Mapping[str, Any], sha: str, tree_oid: str | None) -> list[str]:
    from ionmc import provenance

    if not tree_oid:
        return ["no tree oid given: the package_tree_oid cannot be checked"]
    exp = provenance.host_runner_identity_expectation(sha, tree_oid)
    bad = provenance.host_runner_identity_mismatches(code, exp)
    if code.get("source") != "git":
        bad.append(f"source: {code.get('source')!r} is not 'git'")
    if code.get("git_dirty") is not False:
        bad.append(f"git_dirty: {code.get('git_dirty')!r} is not False")
    return sorted(set(bad))


def cmd_prov7(summary: Path, sha: str, tree_oid: str | None, out: Path | None = None) -> Path:
    """PROV-7 on the code block of the verified sibling bundle; the summary copy must equal it."""
    from ionmc import results

    sm = json.loads(summary.read_text(encoding="utf-8"))
    bad: list[str] = []
    code: Mapping[str, Any] = {}
    try:
        loaded = results.load(_bundle_of(summary))
        code = loaded.doc["provenance"]["code"]
        if sm.get("code") != code:
            bad.append("summary code differs from the bundle's provenance.code")
        bad += prov7_failures(code, sha, tree_oid)
    except Exception as exc:  # noqa: BLE001 - an unreadable bundle is a FAIL, never a skip
        bad.append(f"the bundle does not load: {type(exc).__name__}")
    path = out or summary.with_name("prov7.json")
    shown = {
        k: code.get(k) for k in ("source", "status", "git_commit", "git_dirty", "package_tree_oid")
    }
    _dump(
        path,
        {"verdict": "FAIL" if bad else "PASS", "failures": sorted(set(bad)), "bundle_code": shown},
    )
    _result("PROV-7", not bad, path)  # every check echoes one RESULT line (the host stdout)
    return path


# ---------------------------------------------------------------------------------------------
# host


_SUMMARY = re.compile(r"^=*\s*(\d+ \w+(?:, \d+ \w+)*) in [\d.]+s", re.MULTILINE)
_COUNT = re.compile(r"(\d+) (\w+)")


def parse_pytest_summary(text: str) -> dict[str, int] | None:
    """Counts of the final pytest summary line (``None`` when it is missing)."""
    found = _SUMMARY.findall(text)
    if not found:
        return None
    counts: dict[str, int] = {}
    for n, word in _COUNT.findall(found[-1]):
        counts[word] = counts.get(word, 0) + int(n)
    return counts


def pytest_ok(code: int, counts: dict[str, int] | None) -> bool:
    """Exit 0, a summary line, at least one pass and no skip, failure or error."""
    if code != 0 or counts is None:
        return False
    bad = ("skipped", "failed", "error", "errors", "xfailed", "xpassed")
    return counts.get("passed", 0) > 0 and not any(counts.get(w, 0) for w in bad)


def _pytest(
    targets: Sequence[str],
    label: str,
    out: Path,
    cache_dir: Path,
    timeout: float,
    repo: Path = HERE.parents[2],
) -> tuple[int, bool, dict[str, int] | None]:
    log = f"pytest-{re.sub(r'[^A-Za-z0-9]+', '_', label)}.log"
    code, _dur = _child(
        [sys.executable, "-m", "pytest", *targets, "-p", "no:cacheprovider", "-rs", "-q"],
        repo, out, cache_dir, timeout, log,
    )  # fmt: skip
    try:
        text = (out / "logs" / log).read_text(encoding="utf-8", errors="replace")
    except OSError:
        text = ""
    counts = parse_pytest_summary(text)
    ok = pytest_ok(code, counts)
    print(f"PYTEST {label} {code} {'ok' if ok else 'NOT-OK'} {counts}", flush=True)
    return code, ok, counts


def cmd_host(out: Path, cache_dir: Path, sha: str, tree_oid: str | None) -> Path:
    if not tree_oid:
        raise SystemExit("host needs a non-empty --tree-oid")
    _EXTRA_ENV.update({"IONMC_REQUIRE_DATA": "1", "IONMC_LV_TIER": "1"})
    out.mkdir(parents=True, exist_ok=True)
    s2 = run_case(2, 0, out, cache_dir)
    s0 = run_case(0, 0, out, cache_dir)
    s1 = run_case(1, 0, out, cache_dir)
    summaries = [json.loads(p.read_text(encoding="utf-8")) for p in (s2, s0, s1)]
    cmp_path = cmd_compare([s0, s1, s2], out / "compare.json")
    prov_path = cmd_prov7(s2, sha, tree_oid, out / "prov7.json")
    T_OUT2 = "tests/ionmc/test_v3_006_results.py::test_out2_clean_process_round_trip"
    T_CAP = "tests/ionmc/test_v3_006_contract.py"
    pt = {
        "out2": _pytest([T_OUT2], "out2", out, cache_dir, 1500),
        "cap1": _pytest([T_CAP], "contract", out, cache_dir, 1500),
        "trace": _pytest([TRACE_TEST], "trace", out, cache_dir, 1500),
        "a16": _pytest(list(REG1_TESTS), "a16", out, cache_dir, 1500),
    }
    pytest_codes = {n: v[0] for n, v in pt.items()}
    host = {
        "schema": "ionmc.out6-host-summary/1",
        "sha": sha,
        "cases": {f"r{s['r_index']}k{s['k']}": {"device": s["device"], "backend": s["backend"],
                  "precision": s["precision"], "exit_codes": s["exit_codes"]} for s in summaries},
        "compare": json.loads(cmp_path.read_text(encoding="utf-8"))["verdict"],
        "prov7": json.loads(prov_path.read_text(encoding="utf-8"))["verdict"],
        "pytest_exit_codes": pytest_codes,
        "pytest_counts": {n: v[2] for n, v in pt.items()},
        "pytest_zero_skip_ok": {n: v[1] for n, v in pt.items()},
        "no_network": "every run is a child of the controlled host runner (network denied)",
    }  # fmt: skip
    path = out / "host-summary.json"
    text = json.dumps(host, indent=1, sort_keys=True) + "\n"
    check_summary_text(text, [str(out), str(cache_dir), str(Path.home()), _package_parent()])
    path.write_text(text, encoding="utf-8")
    runs_ok = all(c == 0 for s in summaries for c in s["exit_codes"].values())
    _result("OUT-2-HR", runs_ok and pt["out2"][1], path)
    _result("CLI-4-HR", runs_ok, path)
    _result("REG-1-HR", pt["trace"][1] and pt["a16"][1], path)
    _result("CAP-1b-HR", pt["cap1"][1], path)
    return path


# ---------------------------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="out6_run", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run")
    r.add_argument("--r", type=int, required=True)
    r.add_argument("--k", type=int, default=0)
    r.add_argument("--out", required=True)
    r.add_argument("--cache-dir", required=True)
    r.add_argument("--backend-override")
    r.add_argument("--precision-override")
    c = sub.add_parser("compare")
    c.add_argument("--summaries", nargs="+", required=True)
    c.add_argument("--tier", choices=["hr", "lv"], default="hr")
    i = sub.add_parser("identities")
    i.add_argument("--summary", required=True)
    i.add_argument("--expect", choices=["carbon", "nuclear", "icru90"], required=True)
    i.add_argument("--cache-dir")
    q = sub.add_parser("reproduce")
    q.add_argument("--bundle", required=True)
    q.add_argument("--out", required=True)
    q.add_argument("--cache-dir", required=True)
    v = sub.add_parser("prov7")
    v.add_argument("--summary", required=True)
    v.add_argument("--sha", required=True)
    v.add_argument("--tree-oid", required=True)
    h = sub.add_parser("host")
    h.add_argument("--out", required=True)
    h.add_argument("--cache-dir", required=True)
    h.add_argument("--sha", required=True)
    h.add_argument("--tree-oid", required=True)
    rc = sub.add_parser("_reproduce-child")
    rc.add_argument("--bundle", required=True)
    rc.add_argument("--out", required=True)
    rc.add_argument("--cache-dir", required=True)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    a = build_parser().parse_args(argv)
    if a.command == "run":
        run_case(a.r, a.k, Path(a.out), Path(a.cache_dir),
                 backend=a.backend_override, precision=a.precision_override)  # fmt: skip
    elif a.command == "compare":
        cmd_compare([Path(s) for s in a.summaries], tier=a.tier)
    elif a.command == "identities":
        cmd_identities(Path(a.summary), a.expect, Path(a.cache_dir) if a.cache_dir else None)
    elif a.command == "reproduce":
        cmd_reproduce(Path(a.bundle), Path(a.out), Path(a.cache_dir))
    elif a.command == "prov7":
        cmd_prov7(Path(a.summary), a.sha, a.tree_oid)
    elif a.command == "host":
        cmd_host(Path(a.out), Path(a.cache_dir), a.sha, a.tree_oid)
    else:
        return _reproduce_child(Path(a.bundle), Path(a.out), Path(a.cache_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
