"""Regression rows REG-1..REG-3 of plan v3-006 (S1 part).

REG-1: the A16 digest and T1 trace-fixture tests are not edited (they run unchanged in the same
logged run; this module guards that they still exist). REG-2: ``EffectiveConfig.summary()`` of a
fixed configuration set is byte-identical to the fixture generated at the base commit 95937553
(``data/v3_006/gen_summaries.py``). REG-3: the ``Result`` fields and the top-level keys of
``capabilities()`` of the base commit are still present (additions are allowed). No table values
and no seed from the ``out6`` family appear here; nothing runs transport or worker processes.
"""

from __future__ import annotations

import importlib.util
import json
import os
from dataclasses import fields
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from ionmc.errors import UnsupportedCombinationError
from ionmc.simulation import Result, capabilities

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "v3_006"


def _gen() -> ModuleType:
    spec = importlib.util.spec_from_file_location("v3_006_gen_summaries", DATA / "gen_summaries.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GEN = _gen()
FIXTURE: dict[str, Any] = json.loads((DATA / "summaries-95937553.json").read_text("utf-8"))
SURFACE: dict[str, Any] = json.loads((DATA / "api-surface-95937553.json").read_text("utf-8"))


# REG-1 -----------------------------------------------------------------------------------------

_IDENTITY_TESTS = (
    (
        "test_run_suite.py",
        "test_a16_has_no_recorded_exception_and_the_baseline_is_the_v3_003d_merge",
    ),
    ("test_run_suite.py", "test_a16_gate_accepts_the_intended_change_and_rejects_everything_else"),
    ("test_scoring_transport.py", "test_a16_tallies_do_not_change_the_qualified_outputs"),
    ("test_transport_diagnostic.py", "test_default_trace_is_unchanged_from_the_stored_baseline"),
    ("test_stopping_models.py", "test_id_def_default_table_identity_equals_the_base_commit"),
)


@pytest.mark.parametrize(("module", "name"), _IDENTITY_TESTS)
def test_reg1_a16_and_t1_identity_tests_still_exist(module: str, name: str) -> None:
    """The asserted-identity tests are run unchanged in the same logged run; a deleted or renamed
    one fails here, so REG-1 cannot pass by omission."""
    assert f"def {name}(" in (HERE / module).read_text("utf-8")


# REG-2 -----------------------------------------------------------------------------------------


def test_reg2_fixture_is_canonical_and_holds_no_table_values() -> None:
    text = (DATA / "summaries-95937553.json").read_text("utf-8")
    assert text == GEN.dumps(FIXTURE)  # byte-exact canonical form, one final newline
    assert FIXTURE["base_commit"] == "95937553"
    assert set(FIXTURE["cases"]) >= {
        "proton_box_python",
        "proton_box_warp_cpu_float32",
        "voxel_slab_let_lookup_warp_cpu",
        "carbon_box_warp_cpu",
    }
    assert "/home/" not in text  # machine paths are normalised away

    def longest_numeric_list(obj: Any) -> int:
        if isinstance(obj, list):
            n = len(obj) if obj and all(isinstance(x, int | float) for x in obj) else 0
            return max([n, *(longest_numeric_list(x) for x in obj)])
        if isinstance(obj, dict):
            return max([0, *(longest_numeric_list(v) for v in obj.values())])
        return 0

    assert longest_numeric_list(FIXTURE) <= 8  # identities and digests, never table arrays


@pytest.mark.parametrize("name", sorted(FIXTURE["cases"]))
def test_reg2_effective_summary_is_byte_identical_to_the_base_commit(name: str) -> None:
    case = FIXTURE["cases"][name]
    try:
        summary = GEN.summary_of(name, case["params"])
    except UnsupportedCombinationError as exc:
        if case["requires"]:
            if os.environ.get("IONMC_REQUIRE_DATA") == "1":
                pytest.fail(f"{name}: required data are missing: {exc}")
            pytest.skip(f"{name}: needs {case['requires']} in the cache (IONMC_CACHE_DIR)")
        raise
    if os.environ.get("IONMC_LV_TIER") == "1":
        assert GEN.dumps(summary) == GEN.dumps(case["summary"])  # byte-identical on the LV machine
    else:
        # GitHub-runner floating point differs in the last bits of table-derived values (and so do
        # the table digests computed from them): outside the LV tier compare structurally with a
        # relative tolerance and compare digests of table-derived blocks by presence only.
        _assert_platform_equal(summary, case["summary"], "")


_TABLE_DIGEST_KEYS = {"sha256", "content_sha256", "source_sha256"}


def _assert_platform_equal(a: Any, b: Any, path: str) -> None:
    if isinstance(a, dict) and isinstance(b, dict):
        assert set(a) == set(b), (path, set(a) ^ set(b))
        for k in a:
            if k in _TABLE_DIGEST_KEYS and isinstance(a[k], str) and isinstance(b[k], str):
                assert len(a[k]) == len(b[k]) == 64, f"{path}/{k}"
                continue
            _assert_platform_equal(a[k], b[k], f"{path}/{k}")
    elif isinstance(a, list) and isinstance(b, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _assert_platform_equal(x, y, f"{path}[{i}]")
    elif isinstance(a, float) and isinstance(b, float):
        assert a == b or abs(a - b) <= 1e-9 * max(abs(a), abs(b)), (path, a, b)
    else:
        assert a == b, (path, a, b)


# REG-3 -----------------------------------------------------------------------------------------


def test_reg3_result_fields_of_the_base_commit_are_still_present() -> None:
    names = {f.name for f in fields(Result)}
    missing = [n for n in SURFACE["result_fields"] if n not in names]
    assert not missing, missing


def test_reg3_capabilities_top_level_keys_of_the_base_commit_are_still_present() -> None:
    for key, report in (
        ("capabilities_keys", capabilities()),
        ("capabilities_nuclear_keys", capabilities(nuclear=True)),
    ):
        missing = [k for k in SURFACE[key] if k not in report]
        assert not missing, (key, missing)
