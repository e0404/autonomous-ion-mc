"""Rows OUT-2, OUT-3, CLI-5, CLI-6, PROV-7 (the checks of the S8 driver) of plan v3-006.

The ``out6`` cases and the driver ``validation/scripts/results/out6_run.py`` are exercised on tiny
runs with the CI seed k = 99 of the out6 family (reserved for tests; the LV/HR evidence uses k = 0
and 1). Cache-backed cases (r_index 5 and 7) skip without the data, and fail under
IONMC_REQUIRE_DATA=1 or IONMC_LV_TIER=1. Each ``run`` is three child processes (one at a time).
Not a clinical claim.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from ionmc import provenance as prov
from ionmc import results
from ionmc.configio import load_toml_full
from ionmc.data import cache
from ionmc.simulation import Simulation

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parents[1] / "validation" / "scripts" / "results"
ICRU_DS = "geant4-icru90-stopping-11.4.2"
K_TEST = 99


def _load(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CASES = _load("out6_cases", SCRIPTS / "out6_cases.py")
RUN = _load("out6_run", SCRIPTS / "out6_run.py")
SEEDS = _load("_v3_006_seed_families_for_out6", HERE / "test_v3_006_seed_families.py")


def _cache_dir() -> Path:
    return cache.resolve_cache_dir(None)


def _require(dataset_id: str) -> None:
    required = "1" in (os.environ.get("IONMC_REQUIRE_DATA"), os.environ.get("IONMC_LV_TIER"))
    try:
        cache.verify(dataset_id, _cache_dir())
    except FileNotFoundError as exc:
        if required:
            pytest.fail(f"{dataset_id} is required here but not cached: {exc}")
        pytest.skip(f"{dataset_id} is not in the cache (IONMC_CACHE_DIR)")


def _require_nuclear() -> None:
    try:
        CASES.nuclear_table_id(_cache_dir())
    except RuntimeError as exc:
        if "1" in (os.environ.get("IONMC_REQUIRE_DATA"), os.environ.get("IONMC_LV_TIER")):
            pytest.fail(str(exc))
        pytest.skip(f"{exc} (IONMC_CACHE_DIR)")


# ---- the eight templates and their seeds


def test_cases_are_the_eight_out6_cases_with_formula_seeds() -> None:
    cases = CASES.cases()
    assert len(cases) == 8
    family = SEEDS.out6_seeds()
    assert CASES.OUT6_BASE == SEEDS.OUT6_BASE
    for r, (path, seed, backend, precision, notes) in enumerate(cases):
        assert path.is_file() and notes
        assert seed == CASES.OUT6_BASE + 100 * r and seed in family
        assert (backend, precision) in {
            ("python", "float64"), ("warp-cpu", "float32"), ("warp-cuda", "float32"),
            ("warp-cpu", "float64"),
        }  # fmt: skip
        assert CASES.seed_of(r, K_TEST) in family
        assert f"seed = {seed}\n" in path.read_text(encoding="utf-8")
    with pytest.raises(ValueError, match="out of range"):
        CASES.seed_of(0, 100)
    with pytest.raises(ValueError, match="out of range"):
        CASES.seed_of(50, 0)


@pytest.mark.parametrize("r", range(8))
def test_templates_load_and_validate(r: int, tmp_path: Path) -> None:
    if r == 5:
        _require_nuclear()
    if r == 7:
        _require(ICRU_DS)
    path = CASES.write_config(r, K_TEST, tmp_path, _cache_dir())
    assert "@TABLE_ID@" not in path.read_text(encoding="utf-8")
    loaded = load_toml_full(path)
    cfg = loaded.config
    assert cfg.run.seed == CASES.seed_of(r, K_TEST) and cfg.run.seed in SEEDS.out6_seeds()
    _name, _seed, backend, precision, _notes = CASES.cases()[r]
    assert (cfg.run.backend, cfg.run.precision) == (backend, precision)
    assert cfg.run.n_histories <= 400 and cfg.run.n_batches >= 4
    for grid in cfg.scoring:
        assert len(grid.shape) == 3 and max(grid.shape[:2]) <= 8 and grid.shape[2] <= 16
    if backend != "warp-cuda":  # no CUDA needed to validate the others
        Simulation(cfg)
    if r == 3:
        assert {t.quantity for t in cfg.tallies} == {"let_d", "lookup_sum"}
        assert loaded.manifest.files  # the slab files and the lookup are referenced by digest


def test_write_config_overrides_and_missing_table(tmp_path: Path) -> None:
    p = CASES.write_config(0, 1, tmp_path, tmp_path, backend="warp-cpu", precision="float32")
    run = load_toml_full(p).config.run
    assert (run.backend, run.precision, run.seed) == ("warp-cpu", "float32", CASES.seed_of(0, 1))
    with pytest.raises(RuntimeError, match="no built nuclear table"):
        CASES.write_config(5, K_TEST, tmp_path, tmp_path)  # an empty cache fails, never guesses


# ---- the driver on tiny runs (k = 99)


@pytest.fixture(scope="module")
def py_run(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    out = tmp_path_factory.mktemp("py")
    return out, RUN.run_case(0, K_TEST, out, out / "cache")


@pytest.fixture(scope="module")
def wp_run(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    out = tmp_path_factory.mktemp("wp")
    return out, RUN.run_case(
        0, K_TEST, out, out / "cache", backend="warp-cpu", precision="float64"
    )  # fmt: skip


@pytest.fixture(scope="module")
def r1_run(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    out = tmp_path_factory.mktemp("r1")
    return out, RUN.run_case(1, K_TEST, out, out / "cache")


def _summary(p: Path) -> dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def test_out2_end_to_end_run_is_verified_and_clean_process_read(
    py_run: tuple[Path, Path],
) -> None:
    out, sp = py_run
    s = _summary(sp)
    assert s["exit_codes"] == {"run": 0, "verify": 0, "report_check": 0}
    assert s["seed"] == CASES.seed_of(0, K_TEST) and s["r_index"] == 0 and s["k"] == K_TEST
    assert (s["backend"], s["precision"]) == ("python", "float64")
    assert s["output_path_persisted"] == f"r0k{K_TEST}"  # relative: no path of P
    bundle = out / f"r0k{K_TEST}"
    loaded = results.load(bundle)
    assert s["array_content_sha256"] == {
        e["name"]: e["content_sha256"] for e in loaded.doc["arrays"]
    }
    assert set(s["files"]) >= {"arrays.npz", "report.md"} and len(s["result_json_sha256"]) == 64
    assert s["stopping"]["model_id_effective"] == "bethe-ionmc-v1" and s["nuclear"] is None
    assert s["code"]["package_tree_oid"] == loaded.doc["provenance"]["code"]["package_tree_oid"]


def test_summaries_hold_no_path_of_p_and_no_array(py_run: tuple[Path, Path]) -> None:
    out, sp = py_run
    text = sp.read_text(encoding="utf-8")
    RUN.check_summary_text(text, [str(out), str(Path.home()), str(sys.prefix)])
    assert RUN.P_TEXT.search(text) is None and text.endswith("}\n")
    assert prov.privacy_violations(_summary(sp)) == []  # the privacy rule of the provenance tests
    assert RUN._longest_numeric_list(_summary(sp)) <= 8 and len(text) < 64 * 1024
    for bad in (f"{out}/x", "/home/ann/x", "see /tmp/y"):
        with pytest.raises(AssertionError):
            RUN.check_summary_text(json.dumps({"p": bad}), [str(out)])


def test_run_refuses_an_existing_path(py_run: tuple[Path, Path]) -> None:
    out, _ = py_run
    with pytest.raises(SystemExit, match="exists"):
        RUN.run_case(0, K_TEST, out, out / "cache")


def test_cli5_compare_passes_on_backends_only_and_fails_on_a_mutated_pointer(
    py_run: tuple[Path, Path],
    wp_run: tuple[Path, Path],
    r1_run: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    (_, spy), (_, swp) = py_run, wp_run
    docs = [json.loads((p.parent / f"r0k{K_TEST}" / "result.json").read_text()) for p in (spy, swp)]
    seeds = [(0, K_TEST, CASES.seed_of(0, K_TEST))] * 2
    verdict = RUN.compare_docs(docs, seeds)
    assert verdict["verdict"] == "PASS", verdict
    diff = set(verdict["pairs"][0]["differing_pointers"])
    assert {"/run/backend", "/requested/run/backend", "/run/device"} <= diff
    assert not any(p.startswith(("/requested/source", "/effective/geometry")) for p in diff)
    s1 = r1_run[1]
    path = RUN.cmd_compare([spy, s1], tmp_path / "compare.json", tier="lv")  # r0 and r1, k = 99
    assert json.loads(path.read_text())["verdict"] == "PASS", path.read_text()
    # the case identity binds: three r0 summaries, the wrong tier, a wrong precision all FAIL
    sms = [_summary(spy), _summary(s1)]
    assert RUN.compare_identity_failures(sms, "lv") == []
    three = RUN.compare_identity_failures([sms[0]] * 3, "hr")
    assert any("r_index set" in f for f in three) and any("pairwise" in f for f in three)
    assert RUN.compare_identity_failures(sms, "hr")  # {0, 1} is not the HR set
    wrong = copy.deepcopy(sms)
    wrong[1]["precision"] = "float64"
    assert any("r1" in f and "precision" in f for f in RUN.compare_identity_failures(wrong, "lv"))
    other_k = copy.deepcopy(sms)
    other_k[1]["k"] = 0
    assert any("one k" in f for f in RUN.compare_identity_failures(other_k, "lv"))
    f3 = RUN.cmd_compare([spy, spy, spy], tmp_path / "three.json")
    assert json.loads(f3.read_text())["verdict"] == "FAIL"
    mutated = copy.deepcopy(docs[1])
    mutated["effective"]["source"]["kinetic_energy_mev"] += 1.0
    bad = RUN.compare_docs([docs[0], mutated], seeds)
    assert bad["verdict"] == "FAIL"
    assert bad["pairs"][0]["not_allowed"] == ["/effective/source/kinetic_energy_mev"]
    wrong_seed = RUN.compare_docs(docs, [(0, K_TEST, 12345)] * 2)
    assert wrong_seed["verdict"] == "FAIL"  # a seed that is not the out6 formula value


def test_cli6_identities_detect_a_wrong_model_id(py_run: tuple[Path, Path], tmp_path: Path) -> None:
    out, sp = py_run
    s = _summary(sp)
    doc = json.loads((out / f"r0k{K_TEST}" / "result.json").read_text())
    assert RUN.identity_failures(s, doc, "icru90", None)  # the default model is not icru90
    assert RUN.identity_failures(s, doc, "carbon", None)  # a proton run is not carbon-12
    assert RUN.identity_failures(s, doc, "nuclear", None)  # no nuclear table was used
    wrong = copy.deepcopy(s)
    wrong["stopping"]["model_id_effective"] = "icru90-water-v1"
    assert any("differs" in f for f in RUN.identity_failures(wrong, doc, "icru90", None))
    path = RUN.cmd_identities(sp, "icru90", None, tmp_path / "id.json")
    assert json.loads(path.read_text())["verdict"] == "FAIL"


def test_cli6_carbon_identity_on_a_run(tmp_path: Path) -> None:
    sp = RUN.run_case(4, K_TEST, tmp_path, tmp_path / "cache")
    s = _summary(sp)
    assert s["exit_codes"] == {"run": 0, "verify": 0, "report_check": 0}
    doc = json.loads((tmp_path / f"r4k{K_TEST}" / "result.json").read_text())
    assert RUN.identity_failures(s, doc, "carbon", None) == []
    assert doc["effective"]["source"]["projectile"] == "carbon12"


def test_cli6_icru90_alpha_identity_on_a_run(tmp_path: Path) -> None:
    _require(ICRU_DS)
    sp = RUN.run_case(7, K_TEST, tmp_path, _cache_dir())
    s = _summary(sp)
    assert s["exit_codes"] == {"run": 0, "verify": 0, "report_check": 0}
    doc = json.loads((tmp_path / f"r7k{K_TEST}" / "result.json").read_text())
    assert RUN.identity_failures(s, doc, "icru90", _cache_dir()) == []
    assert s["stopping"]["model_id_effective"] == "icru90-water-v1"


def test_cli6_nuclear_identity_on_a_run(tmp_path: Path) -> None:
    _require_nuclear()
    sp = RUN.run_case(5, K_TEST, tmp_path, _cache_dir())
    s = _summary(sp)
    assert s["exit_codes"] == {"run": 0, "verify": 0, "report_check": 0}
    doc = json.loads((tmp_path / f"r5k{K_TEST}" / "result.json").read_text())
    assert RUN.identity_failures(s, doc, "nuclear", _cache_dir()) == []
    assert s["nuclear"]["table_id"] == s["requested_nuclear_table_id"]
    bad = copy.deepcopy(s)
    bad["requested_nuclear_table_id"] = "0" * 64
    assert RUN.identity_failures(bad, doc, "nuclear", _cache_dir())


def test_out3_reproduction_passes_and_fails_on_a_mutated_digest(
    py_run: tuple[Path, Path], tmp_path: Path
) -> None:
    out, _ = py_run
    bundle = out / f"r0k{K_TEST}"
    path = RUN.cmd_reproduce(bundle, tmp_path, tmp_path / "cache")
    v = json.loads(path.read_text())
    assert v["verdict"] == "PASS", v
    assert v["arrays_equal"] and v["effective_differing_pointers"] == []
    a, b = results.load(bundle).doc, results.load(tmp_path / f"r0k{K_TEST}-repro").doc
    assert RUN.reproduction_verdict(a, b)["ok"]
    mutated = copy.deepcopy(b)
    mutated["arrays"][0]["content_sha256"] = "0" * 64
    bad = RUN.reproduction_verdict(a, mutated)
    assert not bad["ok"] and bad["differing_array_names"] == [mutated["arrays"][0]["name"]]
    mutated = copy.deepcopy(b)
    mutated["effective"]["seed"] += 1
    assert RUN.reproduction_verdict(a, mutated)["effective_differing_pointers"] == ["/seed"]


def test_prov7_fails_on_a_non_git_code_block() -> None:
    sha, oid = "a" * 40, "b" * 40
    good = {"source": "git", "status": "ok", "git_dirty": False, "git_commit": sha,
            "package_tree_oid": oid}  # fmt: skip
    assert RUN.prov7_failures(good, sha, oid) == []
    assert RUN.prov7_failures({**good, "source": "stamp"}, sha, oid)
    assert RUN.prov7_failures({**good, "source": "unknown", "git_commit": None}, sha, oid)
    assert RUN.prov7_failures({**good, "git_dirty": True}, sha, oid)
    assert RUN.prov7_failures({**good, "git_commit": "c" * 40}, sha, oid)
    assert RUN.prov7_failures({**good, "package_tree_oid": "d" * 40}, sha, oid)
    assert RUN.prov7_failures({**good, "git_dirty": None}, sha, oid)
    assert RUN.prov7_failures(good, sha, None) and RUN.prov7_failures(
        good, sha, ""
    )  # never skipped


def test_host_pytest_rows_require_a_summary_and_zero_skips() -> None:
    ok = RUN.parse_pytest_summary("..\n==== 3 passed in 1.20s ====\n")
    assert ok == {"passed": 3} and RUN.pytest_ok(0, ok)
    skipped = RUN.parse_pytest_summary("3 passed, 1 skipped in 2.00s\n")
    assert skipped == {"passed": 3, "skipped": 1} and not RUN.pytest_ok(0, skipped)
    assert RUN.parse_pytest_summary("no summary here") is None
    assert not RUN.pytest_ok(0, None) and not RUN.pytest_ok(1, ok)
    assert not RUN.pytest_ok(0, {"passed": 2, "failed": 1}) and not RUN.pytest_ok(0, {"skipped": 1})
