"""Rows CLI-1, CLI-2, CLI-3, CLI-6 (CI part), CAP-7 and REG-5 of plan v3-006 (step S6): the
``ionmc run/check/capabilities/verify/report`` subcommands and the frozen exception -> exit-code
table of decision 0045 section 7.

Tiny runs only (at most 8 histories, python backend). ``main([...])`` is driven in process; the
exit codes are also observed through one child process at a time (``python -I``, the ``ionmc``
entry function, ``OMP_NUM_THREADS=1``). Cache-backed cases use the ``_require_datasets``
convention of the other V3-006 tests. No test spawns worker processes (REG-5).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from types import MappingProxyType
from typing import Any

import pytest

import ionmc
from ionmc import cli
from ionmc.cli import main
from ionmc.configio import ConfigError, DatasetError, load_toml_full
from ionmc.data import cache
from ionmc.data.acquire import OfflineError
from ionmc.data.registry import DATASETS
from ionmc.errors import (
    BackendUnavailableError,
    TransportLimitError,
    UnsupportedCombinationError,
)
from ionmc.physics import stopping_models
from ionmc.physics.stopping_models import STOPPING_MODELS, DataPin, StoppingModelSpec
from ionmc.results import ResultFormatError, VerificationError

DATA = Path(__file__).resolve().parent / "data" / "v3_006"
CARBON = DATA / "carbon_box.toml"
PROTON = DATA / "proton_box.toml"
ICRU90 = DATA / "icru90_box.toml"
NUCLEAR = DATA / "nuclear_box.toml"
EXAMPLES = (PROTON, CARBON, DATA / "proton_box_warp_material.toml", ICRU90)
ICRU_DS = "geant4-icru90-stopping-11.4.2"
NIST_DS = "nist-pstar-water-2005"
CLEAN_ENV = {"OMP_NUM_THREADS": "1", "IONMC_SINGLE_PROCESS": "1"}


# ---------------------------------------------------------------------------------------------
# helpers


def _require_datasets(*dataset_ids: str) -> None:
    """Skip a cache-backed case without the dataset (CI tier); fail when the dataset is absent
    under IONMC_REQUIRE_DATA=1 or IONMC_LV_TIER=1. A corrupt object is never a skip."""
    cdir = cache.resolve_cache_dir(None)
    required = "1" in (os.environ.get("IONMC_REQUIRE_DATA"), os.environ.get("IONMC_LV_TIER"))
    for dataset_id in dataset_ids:
        try:
            cache.verify(dataset_id, cdir)
        except FileNotFoundError as exc:
            if required:
                pytest.fail(f"{dataset_id} is required here but not cached: {exc}")
            pytest.skip(f"{dataset_id} is not in the cache (IONMC_CACHE_DIR)")


def _nuclear_table_id() -> str:
    """The id of a built, qualified nuclear table of the cache (skip/fail as for datasets)."""
    from ionmc.nuclear.tables import NuclearTable, NuclearTableError

    required = "1" in (os.environ.get("IONMC_REQUIRE_DATA"), os.environ.get("IONMC_LV_TIER"))
    try:
        for p in sorted((cache.resolve_cache_dir(None) / "derived").glob("nuclear-proton-*.json")):
            tid = p.stem.removeprefix("nuclear-proton-")
            try:
                if "transport_path_bound_terms" in NuclearTable.load(None, tid).info:
                    return tid
            except NuclearTableError:
                continue
    except OSError:
        pass
    if required:
        pytest.fail("no built nuclear table in the cache")
    pytest.skip("no built nuclear table in the cache (IONMC_CACHE_DIR)")


def _v(x: Any) -> str:
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, str):
        return json.dumps(x)
    if isinstance(x, list | tuple):
        return "[" + ", ".join(_v(i) for i in x) + "]"
    if isinstance(x, dict):
        return "{ " + ", ".join(f"{k} = {_v(i)}" for k, i in x.items()) + " }"
    return repr(x)


def _base() -> dict[str, Any]:
    return {
        "schema": "ionmc.config",
        "schema_version": "1.0",
        "source": {
            "projectile": "proton",
            "position_mm": [0.0, 0.0, -5.0],
            "direction": [0.0, 0.0, 1.0],
            "kinetic_energy_mev": 20.0,
        },
        "geometry": {
            "type": "box",
            "lower_mm": [-10.0, -10.0, 0.0],
            "size_mm": [20.0, 20.0, 20.0],
            "material": "water",
        },
        "scoring": [
            {
                "name": "g",
                "origin_mm": [-10.0, -10.0, 0.0],
                "spacing_mm": [10.0, 10.0, 10.0],
                "shape": [2, 2, 2],
            }
        ],
        "physics": {"nuclear": False, "stopping_model": "bethe-ionmc-v1"},
        "run": {
            "backend": "python",
            "precision": "float64",
            "seed": 3,
            "n_histories": 4,
            "n_batches": 2,
        },
    }


def _write(tmp: Path, edits: dict[str, Any] | None = None, name: str = "cfg.toml") -> Path:
    """Write the base configuration with ``edits`` ({table: {key: value | None (delete)}} or
    {table: list | None}) as a TOML file (the one input format of the CLI)."""
    doc = _base()
    for tbl, change in (edits or {}).items():
        if change is None:
            del doc[tbl]
        elif isinstance(change, dict) and isinstance(doc.get(tbl), dict):
            for k, val in change.items():
                if val is None:
                    doc[tbl].pop(k, None)
                else:
                    doc[tbl][k] = val
        else:
            doc[tbl] = change
    lines: list[str] = []
    for key, val in doc.items():
        if not isinstance(val, dict | list):
            lines.append(f"{key} = {_v(val)}")
    for key, val in doc.items():
        for item in val if isinstance(val, list) else [val] if isinstance(val, dict) else []:
            lines.append(f"\n[[{key}]]" if isinstance(val, list) else f"\n[{key}]")
            lines.extend(f"{k} = {_v(i)}" for k, i in item.items())
    path = tmp / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _nuc_edits(tid: str, **source: Any) -> dict[str, Any]:
    src = {"kinetic_energy_mev": 100.0, **source}
    return {
        "source": src,
        "physics": {"nuclear": True, "nuclear_table_id": tid},
        "geometry": {"size_mm": [20.0, 20.0, 100.0]},
    }


def _main(argv: list[str], capsys: pytest.CaptureFixture[str]) -> tuple[int, str, str]:
    code = main(argv)
    out = capsys.readouterr()
    return code, out.out, out.err


def _child(
    args: list[str], env_extra: dict[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    """One child process (a tiny CLI call), isolated mode, single-threaded."""
    env = {**os.environ, **CLEAN_ENV, **(env_extra or {})}
    code = "import sys; from ionmc.cli import main; sys.exit(main(sys.argv[1:]))"
    return subprocess.run(
        [sys.executable, "-I", "-c", code, *args],
        capture_output=True, text=True, env=env, timeout=120, check=False,
    )  # fmt: skip


def _leftovers(parent: Path, keep: set[str]) -> list[str]:
    """Anything in ``parent`` besides ``keep``: a bundle, a temporary or an old directory."""
    return sorted(p.name for p in parent.iterdir() if p.name not in keep)


@pytest.fixture(scope="module")
def carbon_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    d = tmp_path_factory.mktemp("carbon")
    assert main(["run", str(CARBON), "--out", str(d / "bundle")]) == 0
    return d / "bundle"


def _doc(bundle: Path) -> dict[str, Any]:
    return json.loads((bundle / "result.json").read_text(encoding="utf-8"))  # type: ignore[no-any-return]


# ---------------------------------------------------------------------------------------------
# CLI-1


def test_cli1_run_verify_report_check_in_process(
    carbon_bundle: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    capsys.readouterr()
    code, out, _ = _main(["verify", str(carbon_bundle)], capsys)
    assert code == 0 and json.loads(out)["tier1"] == "ok"
    code, out, _ = _main(["report", str(carbon_bundle), "--check"], capsys)
    assert code == 0 and out.startswith("report ok")
    code, out, _ = _main(["report", str(carbon_bundle)], capsys)
    assert code == 0 and out == (carbon_bundle / "report.md").read_text(encoding="utf-8")
    assert _doc(carbon_bundle)["input_manifest"]["source"] == "toml"


def test_cli1_run_verify_report_check_child_processes(tmp_path: Path) -> None:
    """The same chain through child processes (one at a time): exit codes 0, 0, 0."""
    bundle = tmp_path / "b"
    r = _child(["run", str(CARBON), "--out", str(bundle)])
    assert r.returncode == 0, r.stderr[-800:]
    assert r.stdout.strip().endswith(str(bundle))
    assert _child(["verify", str(bundle)]).returncode == 0
    assert _child(["report", str(bundle), "--check"]).returncode == 0


def test_cli1_run_uses_output_table_and_refuses_overwrite(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """``[output].path`` and ``persist_batches`` are honoured; an existing bundle is refused (2)
    before any transport unless --overwrite."""
    cfg = _write(tmp_path, {"output": {"path": "out/b", "persist_batches": False}})
    assert main(["run", str(cfg)]) == 0
    bundle = tmp_path / "out" / "b"
    assert _doc(bundle)["requested"] is not None
    assert "batch_energy_mev" not in {e["origin"]["array"] for e in _doc(bundle)["arrays"]}
    assert main(["run", str(_write(tmp_path, name="d.toml")), "--out", str(tmp_path / "d")]) == 0
    assert "batch_energy_mev" in {e["origin"]["array"] for e in _doc(tmp_path / "d")["arrays"]}
    capsys.readouterr()
    code, _, err = _main(["run", str(cfg)], capsys)
    assert code == 2 and "--overwrite" in err
    assert main(["run", str(cfg), "--overwrite"]) == 0
    assert _leftovers(tmp_path / "out", {"b"}) == []
    shutil.rmtree(tmp_path / "d")
    # neither --out nor [output].path
    code, _, err = _main(["run", str(_write(tmp_path, name="nopath.toml"))], capsys)
    assert code == 2 and "--out" in err


# ---------------------------------------------------------------------------------------------
# CLI-2: every line of the exit-code table, provoked separately; no partial bundle


def _clean(parent: Path, keep: set[str]) -> None:
    assert _leftovers(parent, keep) == [], "a bundle or temporary was left behind"


def test_cli2_row1_argparse_usage_error() -> None:
    assert main(["run"]) == 2  # missing positional
    assert main(["verify"]) == 2
    assert main(["capabilities", "--bogus"]) == 2
    assert main(["report", "x", "--bogus"]) == 2
    assert _child(["run"]).returncode == 2


def test_cli2_row1_pre_checks_are_usage_errors(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    cfg = _write(tmp_path)
    assert main(["verify", str(tmp_path), "--rederive"]) == 2  # no --cache-dir
    assert main(["verify", str(tmp_path), "--cache-dir", str(tmp_path)]) == 2
    assert main(["verify", str(tmp_path), "--rederive", "--cache-dir", "c", "--inputs", "x"]) == 2
    capsys.readouterr()
    assert main(["run", str(cfg)]) == 2
    _clean(tmp_path, {"cfg.toml"})


def _row2_cases(tmp: Path) -> dict[str, Path]:
    bad = tmp / "syntax.toml"
    bad.write_text("this is = = not toml\n", encoding="utf-8")
    return {
        "toml-syntax": bad,
        "missing-file": tmp / "absent.toml",
        "unknown-key": _write(tmp, {"physics": {"bogus_key": 1}}, "k.toml"),
        "unknown-table": _write(tmp, {"extra": {"a": 1}}, "t.toml"),
        "wrong-type": _write(tmp, {"run": {"seed": "three"}}, "w.toml"),
        "missing-key": _write(tmp, {"run": {"seed": None}}, "m.toml"),
        "missing-stopping-model": _write(tmp, {"physics": {"stopping_model": None}}, "s.toml"),
        "unknown-stopping-model": _write(
            tmp, {"physics": {"stopping_model": "no-such-model"}}, "u.toml"
        ),
        "stopping-object": _write(tmp, {"physics": {"stopping": {"model_id": "x"}}}, "o.toml"),
    }


def test_cli2_row2_config_error_is_2_and_names_the_key(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """ConfigError (a subclass of UnsupportedCombinationError, which is row 4) exits 2."""
    assert issubclass(ConfigError, UnsupportedCombinationError)
    names = {
        "unknown-key": "bogus_key",
        "wrong-type": "run.seed",
        "missing-key": "run.seed",
        "missing-stopping-model": "physics.stopping_model",
        "unknown-stopping-model": "physics.stopping_model",
        "stopping-object": "physics.stopping",
    }
    cases = _row2_cases(tmp_path)
    keep = {p.name for p in cases.values()}
    for label, cfg in cases.items():
        for command in ("check", "run"):
            argv = [command, str(cfg)] + (["--out", str(tmp_path / "bundle")] * (command == "run"))
            code, _, err = _main(argv, capsys)
            assert code == 2, (label, command, err)
            assert names.get(label, "") in err, (label, err)
    _clean(tmp_path, keep)
    assert _child(["check", str(cases["unknown-key"])]).returncode == 2


def test_cli2_row3_missing_cache_object_is_1(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A missing cache object (DatasetError) exits 1, also through a child process."""
    empty = tmp_path / "empty-cache"
    empty.mkdir()
    cfg = _write(tmp_path, {"physics": {"stopping_model": "icru90-water-v1"}})
    monkeypatch.setenv(cache.CACHE_ENV, str(empty))
    for command in ("check", "run"):
        argv = [command, str(cfg)] + (["--out", str(tmp_path / "bundle")] * (command == "run"))
        code, _, err = _main(argv, capsys)
        assert code == 1 and ICRU_DS in err
    _clean(tmp_path, {"cfg.toml", "empty-cache"})
    r = _child(["check", str(cfg)], {cache.CACHE_ENV: str(empty)})
    assert r.returncode == 1 and ICRU_DS in r.stderr


def test_cli2_row3_pin_mismatch_is_1(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    _require_datasets(NIST_DS)
    cfg = _write(tmp_path, {"physics": {"stopping_model": "nist-pstar-astar-water-2005"}})
    assert main(["check", str(cfg)]) == 0  # positive control
    pin = DATASETS[NIST_DS]
    monkeypatch.setitem(DATASETS, NIST_DS, type(pin)(**{**pin.__dict__, "sha256": "f" * 64}))
    code, _, err = _main(["check", str(cfg)], capsys)
    assert code == 1 and NIST_DS in err


@pytest.mark.parametrize("cause", [cache.IntegrityError, OfflineError, FileNotFoundError])
def test_cli2_row3_cause_chain_beats_row4(
    cause: type[Exception],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An UnsupportedCombinationError whose __cause__ chain holds a cache error exits 1 (a
    dataset parse failure reached through select()); without such a cause it exits 4."""
    chained = UnsupportedCombinationError("stopping data unreadable")
    chained.__cause__ = RuntimeError("wrapper")
    chained.__cause__.__cause__ = cause("cause")
    assert cli._exit_code(chained) == 1
    assert cli._exit_code(UnsupportedCombinationError("plain")) == 4
    assert cli._exit_code(DatasetError("x")) == 1 and cli._exit_code(ConfigError("k", "m")) == 2

    def refuse(_config: Any) -> Any:
        raise chained

    monkeypatch.setattr("ionmc.simulation.validate", refuse)
    cfg = _write(tmp_path)
    for command in ("check", "run"):
        argv = [command, str(cfg)] + (["--out", str(tmp_path / "bundle")] * (command == "run"))
        assert _main(argv, capsys)[0] == 1
    _clean(tmp_path, {"cfg.toml"})


def test_cli2_row4_unsupported_combination_is_4(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    cfg = _write(tmp_path, {"source": {"kinetic_energy_mev": 900.0}})  # above the model domain
    for command in ("check", "run"):
        argv = [command, str(cfg)] + (["--out", str(tmp_path / "bundle")] * (command == "run"))
        code, _, err = _main(argv, capsys)
        assert code == 4 and "UnsupportedCombinationError" in err
    _clean(tmp_path, {"cfg.toml"})
    assert _child(["check", str(cfg)]).returncode == 4


def test_cli2_row5_backend_unavailable_is_5(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("ionmc.config.cuda_available", lambda: False)
    cfg = _write(tmp_path, {"run": {"backend": "warp-cuda", "precision": "float32"}})
    for command in ("check", "run"):
        argv = [command, str(cfg)] + (["--out", str(tmp_path / "bundle")] * (command == "run"))
        code, _, err = _main(argv, capsys)
        assert code == 5 and "BackendUnavailableError" in err
    _clean(tmp_path, {"cfg.toml"})
    assert cli._exit_code(BackendUnavailableError("x")) == 5


def test_cli2_rows6_7_invalid_result_not_allowed_and_allowed(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """max_steps = 1 exhausts the step budget: TransportLimitError (3, no bundle) unless the
    configuration allows an invalid result (3, bundle written with valid = false, verify 0)."""
    refused = _write(tmp_path, {"run": {"max_steps": 1}}, "refused.toml")
    code, _, err = _main(["run", str(refused), "--out", str(tmp_path / "no")], capsys)
    assert code == 3 and "TransportLimitError" in err
    _clean(tmp_path, {"refused.toml"})
    allowed = _write(tmp_path, {"run": {"max_steps": 1, "allow_invalid_result": True}}, "ok.toml")
    code, out, err = _main(["run", str(allowed), "--out", str(tmp_path / "bundle")], capsys)
    assert code == 3 and "valid = false" in err
    assert _doc(tmp_path / "bundle")["valid"] is False
    assert "INVALID" in (tmp_path / "bundle" / "report.md").read_text(encoding="utf-8").upper()
    assert main(["verify", str(tmp_path / "bundle")]) == 0
    assert main(["report", str(tmp_path / "bundle"), "--check"]) == 0
    assert _leftovers(tmp_path, {"refused.toml", "ok.toml", "bundle"}) == []
    assert cli._exit_code(TransportLimitError("x")) == 3


def test_cli2_row3_child_process_invalid_result(tmp_path: Path) -> None:
    cfg = _write(tmp_path, {"run": {"max_steps": 1}})
    r = _child(["run", str(cfg), "--out", str(tmp_path / "bundle")])
    assert r.returncode == 3 and not (tmp_path / "bundle").exists()


def _tamper(bundle: Path, dest: Path, what: str) -> Path:
    shutil.copytree(bundle, dest)
    if what == "arrays":
        raw = bytearray((dest / "arrays.npz").read_bytes())
        raw[len(raw) // 2] ^= 0xFF
        (dest / "arrays.npz").write_bytes(bytes(raw))
    elif what == "report":
        (dest / "report.md").write_text("changed\n", encoding="utf-8")
    elif what == "result":
        (dest / "result.json").write_bytes((dest / "result.json").read_bytes()[:-5])
    return dest


@pytest.mark.parametrize("what", ["arrays", "report", "result", "missing"])
def test_cli2_row8_tampered_bundle_is_6(
    what: str, carbon_bundle: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    b = tmp_path / "b" if what == "missing" else _tamper(carbon_bundle, tmp_path / "b", what)
    for argv in (["verify", str(b)], ["report", str(b)], ["report", str(b), "--check"]):
        code, _, err = _main(argv, capsys)
        assert code == 6 and "ResultFormatError" in err, (argv, err)
    _clean(tmp_path, {"b"} if what != "missing" else set())


def test_cli2_row8_child_process_tampered(carbon_bundle: Path, tmp_path: Path) -> None:
    b = _tamper(carbon_bundle, tmp_path / "b", "arrays")
    assert _child(["verify", str(b)]).returncode == 6


def test_cli2_row8_report_check_detects_a_regenerated_difference(
    carbon_bundle: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """When the writer's output would differ from the bytes on disk (a code change after the
    bundle was written), ``report --check`` exits 6; the bundle is never modified."""
    from ionmc import report

    real = report.regenerate_report
    monkeypatch.setattr(
        "ionmc.report.regenerate_report",
        lambda b: {n: d + b"x" for n, d in real(b).items()},
    )
    before = (carbon_bundle / "report.md").read_bytes()
    code, _, err = _main(["report", str(carbon_bundle), "--check"], capsys)
    assert code == 6 and "report.md" in err
    assert (carbon_bundle / "report.md").read_bytes() == before


def test_cli2_row8_verify_rederive_mismatch_is_6(
    carbon_bundle: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    cdir = str(cache.resolve_cache_dir(None))
    code, out, _ = _main(["verify", str(carbon_bundle), "--rederive", "--cache-dir", cdir], capsys)
    assert code == 0
    t2 = json.loads(out)["tier2"]
    assert t2["ok"] is True and t2["status"] == "reproduced"
    # a tier-2 mismatch: re-derivation reports a difference -> exit 6 (tier 1 stays ok)
    import ionmc.results as results

    real = results.verify_effective

    def differs(*a: Any, **k: Any) -> Any:
        rep = real(*a, **k)
        return type(rep)(False, "differs", ("/effective/seed",), rep.exclusions, (), {})

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(results, "verify_effective", differs)
        code, out, _ = _main(
            ["verify", str(carbon_bundle), "--rederive", "--cache-dir", cdir], capsys
        )
    assert code == 6 and json.loads(out)["tier2"]["differences"] == ["/effective/seed"]
    assert (
        cli._exit_code(VerificationError("x")) == 6 and cli._exit_code(ResultFormatError("x")) == 6
    )


def test_cli2_row9_unexpected_exception_is_70_with_traceback(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    def boom(self: Any) -> Any:
        raise ZeroDivisionError("injected defect")

    monkeypatch.setattr("ionmc.simulation.Simulation.run", boom)
    code, _, err = _main(["run", str(_write(tmp_path)), "--out", str(tmp_path / "bundle")], capsys)
    assert code == 70 and "Traceback" in err and "injected defect" in err
    _clean(tmp_path, {"cfg.toml"})
    assert cli._exit_code(KeyError("k")) == 70 and cli._exit_code(FileExistsError("e")) == 70


def test_cli2_failure_inside_save_leaves_no_partial_bundle(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """save fails once its temporary directory exists: it is removed; exit 6."""
    import ionmc.results as results

    def fail(*_a: Any, **_k: Any) -> Any:
        raise ResultFormatError("injected save failure")

    monkeypatch.setattr(results.report_mod, "render_all", fail)
    code, _, _ = _main(["run", str(_write(tmp_path)), "--out", str(tmp_path / "bundle")], capsys)
    assert code == 6
    _clean(tmp_path, {"cfg.toml"})


def test_cli2_exit_code_table_precedence() -> None:
    """The table as data of behaviour: subclasses before bases (decision 0045 section 7)."""
    cfg_err = ConfigError("physics.nuclear", "bad")
    cfg_err.__cause__ = FileNotFoundError("x")  # row 2 beats the cache-cause line of row 3
    assert cli._exit_code(cfg_err) == 2
    assert cli._exit_code(ValueError("v")) == 70
    assert cli._exit_code(cli._UsageError("u")) == 2


# ---------------------------------------------------------------------------------------------
# CLI-3


def test_cli3_capabilities_json_equals_capabilities(capsys: pytest.CaptureFixture[str]) -> None:
    for nuclear in (False, True):
        argv = ["capabilities"] + (["--nuclear"] if nuclear else [])
        code, out, _ = _main(argv, capsys)
        assert code == 0
        assert json.loads(out) == json.loads(json.dumps(ionmc.capabilities(nuclear)))
    assert json.loads(_main(["capabilities"], capsys)[1]) != json.loads(
        _main(["capabilities", "--nuclear"], capsys)[1]
    )


def test_cli3_capabilities_child_process() -> None:
    r = _child(["capabilities"])
    assert r.returncode == 0
    assert json.loads(r.stdout) == json.loads(json.dumps(ionmc.capabilities()))


def test_check_prints_the_effective_summary_and_runs_no_transport(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    def forbidden(self: Any) -> Any:
        raise AssertionError("check must not transport")

    monkeypatch.setattr("ionmc.simulation.Simulation.run", forbidden)
    code, out, _ = _main(["check", str(CARBON)], capsys)
    assert code == 0
    from ionmc.simulation import Simulation

    eff = Simulation(load_toml_full(CARBON).config).effective.summary()
    assert json.loads(out) == json.loads(json.dumps(eff))
    assert main(["check", str(_write(tmp_path))]) == 0
    assert not list(tmp_path.glob("*bundle*"))


# ---------------------------------------------------------------------------------------------
# CLI-6 (CI part): stopping / nuclear identities in the bundles of the TOML examples


def _identity_of(config: Path) -> dict[str, Any]:
    from ionmc.physics.stopping_models import spec_sha256

    mid = load_toml_full(config).config.physics.stopping_model
    assert mid is not None
    return {"id": mid, "sha": spec_sha256(STOPPING_MODELS[mid])}


def test_cli6_carbon_bundle_carries_its_stopping_identity(carbon_bundle: Path) -> None:
    ident = _identity_of(CARBON)
    st = _doc(carbon_bundle)["effective"]["stopping"]
    assert st["model_id_effective"] == st["model_id_requested"] == ident["id"] == "bethe-ionmc-v1"
    assert st["spec_sha256"] == ident["sha"]
    assert _doc(carbon_bundle)["requested"]["source"]["projectile"] == "carbon12"


def test_cli6_icru90_bundle_carries_its_stopping_identity(tmp_path: Path) -> None:
    _require_datasets(ICRU_DS)
    assert main(["run", str(ICRU90), "--out", str(tmp_path / "b")]) == 0
    ident = _identity_of(ICRU90)
    doc = _doc(tmp_path / "b")
    st = doc["effective"]["stopping"]
    assert st["model_id_effective"] == ident["id"] == "icru90-water-v1"
    assert st["spec_sha256"] == ident["sha"]
    assert ICRU_DS in json.dumps(doc["provenance"])
    assert main(["verify", str(tmp_path / "b")]) == 0
    assert main(["report", str(tmp_path / "b"), "--check"]) == 0


def test_cli6_nuclear_bundle_carries_its_table_id(tmp_path: Path) -> None:
    tid = _nuclear_table_id()
    cfg = tmp_path / "nuclear.toml"
    cfg.write_text(NUCLEAR.read_text(encoding="utf-8").replace("@TABLE_ID@", tid), encoding="utf-8")
    assert main(["run", str(cfg), "--out", str(tmp_path / "b")]) == 0
    doc = _doc(tmp_path / "b")
    assert doc["effective"]["nuclear"]["table_id"] == tid
    assert doc["provenance"]["data"]["nuclear"]["table_id"] == tid
    assert doc["requested"]["physics"]["nuclear_table_id"] == tid
    assert main(["verify", str(tmp_path / "b")]) == 0
    assert main(["report", str(tmp_path / "b"), "--check"]) == 0


# ---------------------------------------------------------------------------------------------
# CAP-7: the V3-005A C1 list and the V3-013 DOM/ROLE rows through check and run

CaseEdit = Callable[[str], dict[str, Any]]
ICRU = {"physics": {"stopping_model": "icru90-water-v1"}}
NIST = {"physics": {"stopping_model": "nist-pstar-astar-water-2005"}}


def _merge(*parts: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for part in parts:
        for k, v in part.items():
            out[k] = {**out.get(k, {}), **v} if isinstance(v, dict) else v
    return out


def _alpha(e: float) -> dict[str, Any]:
    return {"source": {"projectile": "alpha", "kinetic_energy_mev": e}}


def _ion(name: str) -> dict[str, Any]:
    return {"source": {"projectile": name, "kinetic_energy_mev": 1200.0}}


def _tally(species: str) -> dict[str, Any]:
    return {"tallies": [{"name": "t", "grid": "g", "quantity": "dose", "species": [species]}]}


# (case id, datasets needed, nuclear table needed, edits builder(tid), expected exit code)
CAP7: tuple[tuple[str, tuple[str, ...], bool, CaseEdit, int], ...] = (
    # --- V3-005A row C1 (domain / role / combination 4; missing table 1) ---------------------
    ("c1-warp-cpu-nuclear", (), True,
     lambda t: _merge(_nuc_edits(t), {"run": {"backend": "warp-cpu", "precision": "float32"}}), 4),
    ("c1-warp-cuda-nuclear", (), True,
     lambda t: _merge(_nuc_edits(t), {"run": {"backend": "warp-cuda", "precision": "float32"}}), 4),
    ("c1-e0-plus-6sigma-above-250", (), True,
     lambda t: _nuc_edits(t, kinetic_energy_mev=245.0, energy_sigma_mev=1.0), 4),
    ("c1-e0-above-250", (), True, lambda t: _nuc_edits(t, kinetic_energy_mev=260.0), 4),
    ("c1-non-proton-source", (), True, lambda t: _nuc_edits(t, projectile="deuteron"), 4),
    ("c1-nist-star-with-nuclear", (NIST_DS,), True,
     lambda t: _merge(_nuc_edits(t), NIST), 4),  # also DOM-5 (existing rule)
    ("c1-unsupported-element", (), True,
     lambda t: _merge(_nuc_edits(t), {"geometry": {"material": "lead"}}), 4),
    ("c1-no-table-id", (), False,
     lambda t: {"physics": {"nuclear": True}, "source": {"kinetic_energy_mev": 100.0}}, 4),
    ("c1-unproducible-alpha", (), True, lambda t: _merge(_nuc_edits(t), _tally("alpha")), 4),
    ("c1-unproducible-triton", (), True, lambda t: _merge(_nuc_edits(t), _tally("triton")), 4),
    ("c1-missing-table", (), False, lambda t: _nuc_edits("0" * 64), 1),
    # --- V3-013 DOM rows -----------------------------------------------------------------------
    ("dom1-icru90-carbon12", (ICRU_DS,), False,
     lambda t: _merge(ICRU, _ion("carbon12")), 4),
    ("dom1-icru90-oxygen16", (ICRU_DS,), False,
     lambda t: _merge(ICRU, _ion("oxygen16")), 4),
    ("dom1-icru90-deuteron", (ICRU_DS,), False,
     lambda t: _merge(ICRU, {"source": {"projectile": "deuteron"}}), 4),
    ("dom1-nist-carbon12", (NIST_DS,), False,
     lambda t: _merge(NIST, _ion("carbon12")), 4),
    ("dom2-icru90-pmma", (ICRU_DS,), False,
     lambda t: _merge(ICRU, {"geometry": {"material": "pmma"}}), 4),
    ("dom2-nist-pmma", (NIST_DS,), False,
     lambda t: _merge(NIST, {"geometry": {"material": "pmma"}}), 4),
    ("dom3-icru90-alpha-above-250-per-u", (ICRU_DS,), False,
     lambda t: _merge(ICRU, _alpha(1100.0)), 4),
    ("dom3-icru90-proton-above-500", (ICRU_DS,), False,
     lambda t: _merge(ICRU, {"source": {"kinetic_energy_mev": 520.0}}), 4),
    ("dom4-icru90-water-I-not-78", (ICRU_DS,), False,
     lambda t: _merge(ICRU, {
         "materials": [{"name": "water80", "density_g_cm3": 1.0,
                        "mass_fractions": {"H": 0.111894, "O": 0.888106}, "I_eV": 80.0}],
         "geometry": {"material": "water80"}}), 4),
    ("dom5-icru90-with-nuclear", (ICRU_DS,), True, lambda t: _merge(_nuc_edits(t), ICRU), 4),
    ("dom6-missing-cache-object", (), False,
     lambda t: ICRU, 1),  # run with an empty cache directory (see the test)
    ("dom7-unknown-model-id", (), False,
     lambda t: {"physics": {"stopping_model": "no-such-model"}}, 2),
    ("dom7-neither-id-nor-object", (), False, lambda t: {"physics": {"stopping_model": None}}, 2),
    ("dom7-both-id-and-object", (), False,
     lambda t: {"physics": {"stopping": {"model_id": "bethe-ionmc-v1"}}}, 2),
    ("dom7-object-without-spec", (), False,
     lambda t: {"physics": {"stopping_model": None, "stopping": {"kind": "analytic"}}}, 2),
)  # fmt: skip


def test_cap7_case_count_covers_the_c1_list_and_the_dom_rows() -> None:
    ids = [c[0] for c in CAP7]
    assert len(ids) == len(set(ids)) == 26
    for row in ("c1-", "dom1-", "dom2-", "dom3-", "dom4-", "dom5-", "dom6-", "dom7-"):
        assert any(i.startswith(row) for i in ids), row


@pytest.mark.parametrize("case", CAP7, ids=[c[0] for c in CAP7])
def test_cap7_check_and_run_give_the_table_code_and_no_bundle(
    case: tuple[str, tuple[str, ...], bool, CaseEdit, int],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    label, datasets, needs_table, build, expected = case
    if label != "dom6-missing-cache-object":
        _require_datasets(*datasets)
    tid = _nuclear_table_id() if needs_table else "unused"
    if label == "dom6-missing-cache-object":
        empty = tmp_path / "empty-cache"
        empty.mkdir()
        monkeypatch.setenv(cache.CACHE_ENV, str(empty))
    cfg = _write(tmp_path, build(tid))
    keep = {cfg.name, "empty-cache"}
    for command in ("check", "run"):
        argv = [command, str(cfg)] + (["--out", str(tmp_path / "bundle")] * (command == "run"))
        code, _, err = _main(argv, capsys)
        assert code == expected, (label, command, code, err[-400:])
        assert "Traceback" not in err
    _clean(tmp_path, keep)


def test_cap7_positive_controls_and_role3(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The base configuration and ROLE-3 (z = 1 with ls-eff: requested != effective id) pass."""
    assert main(["check", str(_write(tmp_path))]) == 0
    capsys.readouterr()
    cfg = _write(tmp_path, {"physics": {"stopping_model": "bethe-ionmc-v1+ls-eff"}}, "ls.toml")
    code, out, _ = _main(["check", str(cfg)], capsys)
    st = json.loads(out)["stopping"]
    assert code == 0
    assert st["model_id_requested"] == "bethe-ionmc-v1+ls-eff"
    assert st["model_id_effective"] == "bethe-ionmc-v1"


def test_cap7_role1_non_construction_pin_is_4(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    default = STOPPING_MODELS["bethe-ionmc-v1"]
    for role in ("evaluation-only", "lineage-only"):
        pin = DataPin("synthetic", "0" * 64, None, None, "none", "synthetic", role)
        synth = StoppingModelSpec(
            "synthetic-role-test", "1", "tabulated", {}, default.interpolation,
            default.range_start, (pin,), default.domain, "test",
        )  # fmt: skip
        with monkeypatch.context() as mp:
            mp.setattr(
                stopping_models, "STOPPING_MODELS",
                MappingProxyType({**STOPPING_MODELS, synth.model_id: synth}),
            )  # fmt: skip
            cfg = _write(tmp_path, {"physics": {"stopping_model": "synthetic-role-test"}})
            for command in ("check", "run"):
                argv = [command, str(cfg)] + (
                    ["--out", str(tmp_path / "bundle")] * (command == "run")
                )
                code, _, err = _main(argv, capsys)
                assert code == 4 and role in err, (role, command, code, err[-300:])
    _clean(tmp_path, {"cfg.toml"})


def test_cap7_warp_plus_nuclear_exits_4_until_v3_005b(tmp_path: Path) -> None:
    """One child process: warp-cpu with nuclear = true is refused (4) until V3-005B."""
    tid = _nuclear_table_id()
    cfg = _write(
        tmp_path, _merge(_nuc_edits(tid), {"run": {"backend": "warp-cpu", "precision": "float32"}})
    )
    r = _child(["run", str(cfg), "--out", str(tmp_path / "bundle")])
    assert r.returncode == 4 and "V3-005B" in r.stderr
    _clean(tmp_path, {"cfg.toml"})


# ---------------------------------------------------------------------------------------------
# REG-5: single process only


def test_reg5_cli_runs_single_process() -> None:
    """Every example configuration has cpu_workers = 1 (the default), the CLI never raises it and
    imports no process machinery, and the tests set IONMC_SINGLE_PROCESS for their children."""
    for path in EXAMPLES:
        assert load_toml_full(path).config.run.cpu_workers == 1, path.name
    src = Path(cli.__file__).read_text(encoding="utf-8")
    assert "multiprocessing" not in src and "cpu_workers" not in src and "ProcessPool" not in src
    assert CLEAN_ENV["IONMC_SINGLE_PROCESS"] == "1" and CLEAN_ENV["OMP_NUM_THREADS"] == "1"
