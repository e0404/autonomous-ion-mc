import copy

import pytest

from infrastructure.experiment_v2.common import file_hash
from infrastructure.experiment_v2.release import (
    BACKENDS,
    CATEGORIES,
    METRICS,
    WORKLOADS,
    qualify,
)


@pytest.fixture
def bundle(tmp_path):
    artifact = tmp_path / "evidence.json"
    artifact.write_text('{"actual_test_fixture":true}')
    ref = {"path": artifact.name, "sha256": file_hash(artifact)}
    requirements = {"MUST-1", "MUST-2"}
    plan = {
        "schema_version": 2,
        "status": "frozen",
        "hardware": "fixture CPU/GPU",
        "calibration_artifacts": ["fixture-calibration"],
        "suites": [],
        "performance_targets": {
            w: {m: {"direction": "max", "value": 10} for m in METRICS}
            for w in WORKLOADS
        },
    }
    report = {
        "code_sha": "a" * 40,
        "dirty": False,
        "scope": "full",
        "filters": [],
        "plan_sha256": "b" * 64,
        "requirements": dict.fromkeys(requirements, "satisfied"),
        "suites": [],
        "performance": {
            w: {
                "hardware": plan["hardware"],
                "effective_physics_verified": True,
                "measurement_artifacts": [ref],
                "metrics": dict.fromkeys(METRICS, 1),
            }
            for w in WORKLOADS
        },
    }
    for category in sorted(CATEGORIES):
        plan["suites"].append(
            {
                "id": category,
                "category": category,
                "requirements": sorted(requirements),
                "backends": sorted(BACKENDS),
                "acceptance": "predeclared fixture",
                "rationale": "software test",
                "independent_required": category == "physics",
                "species": ["proton", "helium", "carbon", "oxygen"],
            }
        )
        report["suites"].append(
            {
                "id": category,
                "status": "passed",
                "code_sha": "a" * 40,
                "tests_executed": 1,
                "skipped": 0,
                "backends": sorted(BACKENDS),
                "artifacts": [ref],
                "independence": {
                    "category": "measured",
                    "rationale": "fixture",
                    "shared_lineage": [],
                    "role": "evaluation",
                    "source_ids": ["fixture"],
                    "sufficiency_reviewed": True,
                },
            }
        )
    return (
        plan,
        report,
        {
            "sha": "a" * 40,
            "requirement_ids": requirements,
            "artifacts": tmp_path,
            "plan_hash": "b" * 64,
        },
    )


def test_complete_fixture_passes(bundle):
    plan, report, kw = bundle
    assert qualify(plan, report, **kw)["release_ready"]


@pytest.mark.parametrize(
    "defect",
    [
        "zero_plan",
        "zero_results",
        "subset",
        "unknown_filter",
        "missing_scope",
        "stale_sha",
        "dirty",
        "missing_must",
        "unresolved_must",
        "missing_backend",
        "skipped",
        "zero_tests",
        "duplicate",
        "missing_artifact",
        "hash_changed",
        "internal_only",
        "calibration_reused",
        "missing_lineage",
        "no_targets",
        "performance_failure",
        "missing_workload",
        "nonfinite",
        "bool_tests",
        "unfrozen",
        "different_plan",
        "no_ions",
        "missing_physics",
        "missing_product",
        "no_measurement",
    ],
)
def test_fail_closed(bundle, defect):
    plan, report, kw = copy.deepcopy(bundle)
    physics = next(r for r in report["suites"] if r["id"] == "physics")
    if defect == "zero_plan":
        plan["suites"] = []
    if defect == "zero_results":
        report["suites"] = []
    if defect == "subset":
        report["scope"] = "diagnostic"
    if defect == "unknown_filter":
        report["filters"] = ["NO_SUCH_SUITE"]
    if defect == "missing_scope":
        del report["scope"]
    if defect == "stale_sha":
        report["code_sha"] = "c" * 40
    if defect == "dirty":
        report["dirty"] = True
    if defect == "missing_must":
        del report["requirements"]["MUST-1"]
    if defect == "unresolved_must":
        report["requirements"]["MUST-1"] = "deferred"
    if defect == "missing_backend":
        physics["backends"] = ["python"]
    if defect == "skipped":
        physics["skipped"] = 1
    if defect == "zero_tests":
        physics["tests_executed"] = 0
    if defect == "bool_tests":
        physics["tests_executed"] = True
    if defect == "duplicate":
        report["suites"].append(physics)
    if defect == "missing_artifact":
        physics["artifacts"] = []
    if defect == "hash_changed":
        physics["artifacts"][0]["sha256"] = "0" * 64
    if defect == "internal_only":
        physics["independence"]["category"] = "backend-parity"
    if defect == "calibration_reused":
        physics["independence"]["role"] = "calibration"
    if defect == "missing_lineage":
        del physics["independence"]["shared_lineage"]
    if defect == "no_targets":
        plan["performance_targets"] = {}
    if defect == "performance_failure":
        report["performance"]["influence"]["metrics"]["wall_seconds"] = 100
    if defect == "missing_workload":
        del report["performance"]["carbon-fragments"]
    if defect == "nonfinite":
        report["performance"]["influence"]["metrics"]["gpu_peak_mib"] = float("nan")
    if defect == "unfrozen":
        plan["status"] = "unratified"
    if defect == "different_plan":
        report["plan_sha256"] = "wrong"
    if defect == "no_ions":
        next(s for s in plan["suites"] if s["category"] == "physics")["species"] = [
            "proton"
        ]
    if defect == "missing_physics":
        plan["suites"] = [s for s in plan["suites"] if s["category"] != "physics"]
    if defect == "missing_product":
        plan["suites"] = [s for s in plan["suites"] if s["category"] != "workflow"]
    if defect == "no_measurement":
        report["performance"]["influence"]["measurement_artifacts"] = []
    assert not qualify(plan, report, **kw)["release_ready"]


@pytest.mark.parametrize(
    "value",
    [
        {},
        None,
        [],
        {"suites": None},
        {"schema_version": 2, "status": "frozen", "suites": [None]},
    ],
)
def test_malformed_inputs_never_pass(bundle, value):
    _, report, kw = bundle
    assert not qualify(value, report, **kw)["release_ready"]


def test_promotion_requires_exact_check_level_evidence():
    from infrastructure.experiment_v2.promote import check_ci

    pr = {"headRefOid": "sha", "mergeStateStatus": "CLEAN", "statusCheckRollup": []}
    with pytest.raises(ValueError):
        check_ci(pr, "sha")
    pr["statusCheckRollup"] = [
        {"name": n, "conclusion": "SUCCESS"} for n in ("pre-commit", "tests", "docs")
    ]
    check_ci(pr, "sha")
    with pytest.raises(ValueError):
        check_ci(pr, "other-sha")
    pr["statusCheckRollup"][0]["conclusion"] = "SKIPPED"
    with pytest.raises(ValueError):
        check_ci(pr, "sha")
