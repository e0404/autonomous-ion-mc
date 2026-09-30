"""Fail-closed release evidence qualification, independent of v1's validator."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path

from infrastructure.experiment_v3.common import (
    ROOT,
    event,
    exact_state,
    file_hash,
    git,
    safe_path,
    write_json,
)

BACKENDS = {"python", "warp-cpu", "warp-cuda"}
CATEGORIES = {
    "physics",
    "numerical",
    "capability",
    "provenance",
    "performance",
    "workflow",
}
INDEPENDENT = {
    "measured",
    "ion-specific-tabulated",
    "independent-monte-carlo",
    "independent-theory",
}
TAXONOMY = INDEPENDENT | {"related-model", "backend-parity", "self-consistency"}
WORKLOADS = {"proton-3d", "influence", "carbon-fragments"}
METRICS = {
    "wall_seconds",
    "histories_per_second",
    "cold_start_seconds",
    "compilation_seconds",
    "host_seconds",
    "gpu_kernel_seconds",
    "gpu_peak_mib",
    "host_peak_mib",
}


def positive_number(value):
    return type(value) in (float, int) and math.isfinite(value) and value > 0


def validate_plan(plan, requirement_ids):
    if not isinstance(plan, dict):
        raise ValueError("Release plan must be an object")
    if plan.get("schema_version") != 2 or plan.get("status") != "frozen":
        raise ValueError("Release plan is not frozen")
    suites = plan.get("suites")
    if not isinstance(suites, list) or not suites:
        raise ValueError("Nonempty required release membership is mandatory")
    names, coverage, backends, categories = set(), set(), set(), set()
    for suite in suites:
        name = suite["id"]
        if not isinstance(name, str) or not name or name in names:
            raise ValueError("Empty/duplicate suite ID")
        names.add(name)
        if (
            not suite.get("requirements")
            or not set(suite["requirements"]) <= requirement_ids
        ):
            raise ValueError("Suite must map to known MUST requirements")
        coverage.update(suite["requirements"])
        backends.update(suite["backends"])
        categories.add(suite["category"])
        if not suite.get("acceptance") or not suite.get("rationale"):
            raise ValueError("Acceptance criteria and rationale must be frozen")
    if (
        coverage != requirement_ids
        or not BACKENDS <= backends
        or not CATEGORIES <= categories
    ):
        raise ValueError("Incomplete requirement/backend/category membership")
    physics = [s for s in suites if s["category"] == "physics"]
    if not physics or any(s.get("independent_required") is not True for s in physics):
        raise ValueError(
            "Each major physical qualification requires independent evidence"
        )
    if not {"proton", "helium", "carbon", "oxygen"} <= {
        x for s in physics for x in s.get("species", [])
    }:
        raise ValueError("Ion-specific release evidence is incomplete")
    targets = plan.get("performance_targets", {})
    if (
        set(targets) != WORKLOADS
        or not plan.get("calibration_artifacts")
        or not plan.get("hardware")
    ):
        raise ValueError(
            "All workloads need ratified targets and calibration/hardware provenance"
        )
    for target in targets.values():
        if not METRICS <= set(target):
            raise ValueError("Incomplete performance/resource targets")
        for value in target.values():
            if value["direction"] not in ("min", "max") or not positive_number(
                value["value"]
            ):
                raise ValueError("Invalid numeric performance target")
    return suites


def qualify(plan, report, *, sha, requirement_ids, artifacts, plan_hash):
    """Return reasons on failure; malformed input must never return ready."""
    errors = []
    try:
        suites = validate_plan(plan, requirement_ids)
        if report.get("code_sha") != sha or report.get("dirty") is not False:
            errors.append("Evidence is not for the exact clean current SHA")
        if report.get("scope") != "full" or report.get("filters") != []:
            errors.append("Filtered/subset/missing-scope runs are diagnostic only")
        if report.get("plan_sha256") != plan_hash:
            errors.append("Evidence does not match the frozen plan")
        records = report["suites"]
        by_id = {s["id"]: s for s in records}
        if len(by_id) != len(records) or set(by_id) != {s["id"] for s in suites}:
            errors.append("Missing, extra or duplicate suite results")
        statuses = report.get("requirements", {})
        if set(statuses) != requirement_ids or any(
            v != "satisfied" for v in statuses.values()
        ):
            errors.append("Unresolved or missing MUST requirements")
        for suite in suites:
            result = by_id.get(suite["id"], {})
            if (
                result.get("status") != "passed"
                or result.get("code_sha") != sha
                or type(result.get("tests_executed")) is not int
                or result["tests_executed"] <= 0
                or result.get("skipped") != 0
            ):
                errors.append(
                    f"{suite['id']}: failed/missing/stale/empty/skipped evidence"
                )
            if not set(suite["backends"]) <= set(result.get("backends", [])):
                errors.append(f"{suite['id']}: missing required backend evidence")
            refs = result.get("artifacts", [])
            if not refs:
                errors.append(f"{suite['id']}: no evidence artifacts")
            for artifact in refs:
                path = safe_path(artifacts, artifact["path"])
                if not path.is_file() or file_hash(path) != artifact["sha256"]:
                    errors.append(f"{suite['id']}: missing or changed artifact")
            if suite.get("independent_required"):
                evidence = result.get("independence", {})
                if (
                    evidence.get("category") not in INDEPENDENT
                    or not evidence.get("rationale")
                    or "shared_lineage" not in evidence
                    or evidence.get("role") != "evaluation"
                    or not evidence.get("source_ids")
                    or evidence.get("sufficiency_reviewed") is not True
                ):
                    errors.append(
                        f"{suite['id']}: missing independent evaluation/lineage review"
                    )
        measures = report.get("performance", {})
        if set(measures) != WORKLOADS:
            errors.append("Missing performance workload evidence")
        for name, targets in plan["performance_targets"].items():
            values = measures.get(name, {})
            if (
                values.get("effective_physics_verified") is not True
                or values.get("hardware") != plan["hardware"]
            ):
                errors.append(f"{name}: physics/hardware not verified")
            if not values.get("measurement_artifacts"):
                errors.append(f"{name}: missing measurement artifacts")
            for artifact in values.get("measurement_artifacts", []):
                path = safe_path(artifacts, artifact["path"])
                if not path.is_file() or file_hash(path) != artifact["sha256"]:
                    errors.append(f"{name}: invalid measurement artifact")
            for metric, target in targets.items():
                value = values.get("metrics", {}).get(metric)
                if not positive_number(value):
                    errors.append(f"{name}: invalid/missing {metric}")
                elif (target["direction"] == "max" and value > target["value"]) or (
                    target["direction"] == "min" and value < target["value"]
                ):
                    errors.append(f"{name}: {metric} target failed")
    except (KeyError, TypeError, ValueError, OSError, AttributeError) as exc:
        errors.append(f"Invalid qualification input: {exc}")
    return {
        "schema_version": 2,
        "code_sha": sha,
        "plan_sha256": plan_hash,
        "release_ready": not errors,
        "errors": errors,
    }


def evaluate(root, plan_path, freeze_sha, report_path, artifact_root):
    sha = exact_state(root)
    frozen_sha = git(root, "rev-parse", freeze_sha + "^{commit}")
    subprocess.run(
        ["git", "-C", str(root), "merge-base", "--is-ancestor", frozen_sha, sha],
        check=True,
    )
    plan_file = safe_path(root, plan_path)
    frozen = subprocess.check_output(
        ["git", "-C", str(root), "show", f"{frozen_sha}:{plan_path}"]
    )
    if plan_file.read_bytes() != frozen:
        raise ValueError(
            "Release plan changed after freeze; explicit amendment "
            "and requalification required"
        )
    ids = {
        r["id"]
        for r in json.loads(
            (Path(root) / "experiment/v3/requirements-index.json").read_text()
        )["requirements"]
    }
    result = qualify(
        json.loads(frozen),
        json.loads(Path(report_path).read_text()),
        sha=sha,
        requirement_ids=ids,
        artifacts=artifact_root,
        plan_hash=file_hash(plan_file),
    )
    result["freeze_sha"] = frozen_sha
    result["report_sha256"] = file_hash(report_path)
    event("release_evaluated", sha=sha, details=result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=ROOT)
    p.add_argument("--plan", default="validation/release-plan.json")
    p.add_argument("--freeze-sha", required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--artifacts", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    try:
        result = evaluate(a.root, a.plan, a.freeze_sha, a.report, a.artifacts)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        result = {"release_ready": False, "errors": [str(exc)]}
    write_json(a.output, result)
    print(json.dumps(result, indent=2))
    return 0 if result["release_ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
