"""Choose economical explicit subagent models; log exceptional top-tier use."""

import json
import re
import sys

from infrastructure.experiment_v3.common import event

SCIENTIFIC = {
    "physics-researcher",
    "scientific-planner",
    "performance-specialist",
    "reviewer",
    "Plan",
}


def route(payload):
    if payload.get("tool_name") not in ("Agent", "Task"):
        return {}
    args = dict(payload.get("tool_input", {}))
    role = args.get("subagent_type", "general-purpose")
    model = args.get("model")
    if not model or model == "inherit":
        model = "opus" if role in SCIENTIFIC else "sonnet"
    family = next(
        (
            m
            for m in ("sonnet", "opus", "haiku", "fable")
            if model == m or model.startswith("claude-" + m + "-")
        ),
        None,
    )
    reason = None
    if family not in ("sonnet", "opus", "haiku"):
        match = re.search(
            r"(?m)^TOP_TIER_JUSTIFICATION:\s*(.{40,})$", args.get("prompt", "")
        )
        if family != "fable" or not match:
            return {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "deny",
                    "permissionDecisionReason": (
                        "Use Opus for scientific research/planning and Sonnet for "
                        "planned implementation/retrieval. Fable requires "
                        "TOP_TIER_JUSTIFICATION: with a concrete explanation of "
                        "why lower tiers are insufficient in the delegated prompt."
                    ),
                }
            }
        reason = match.group(1)[:600]
    args["model"] = model
    return {
        "hookSpecificOutput": {"hookEventName": "PreToolUse", "updatedInput": args},
        "_routing": {"role": role, "model": model, "top_tier_justification": reason},
    }


def main():
    result = route(json.load(sys.stdin))
    routing = result.pop("_routing", None)
    if routing:
        event("agent_routed", details=routing)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
