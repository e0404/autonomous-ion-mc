"""Operational evidence-blocker context required for v2 scientific interventions."""


def validate_context(context):
    if not isinstance(context, dict):
        raise ValueError("V2 intervention requires structured investigation context")
    for key in ("blocked_claim", "budget", "alternatives", "attempts"):
        if not context.get(key):
            raise ValueError(f"Missing intervention context: {key}")
    if not isinstance(context["attempts"], list):
        raise ValueError("attempts must be a list")
    for attempt in context["attempts"]:
        if not isinstance(attempt, dict) or not all(
            attempt.get(k) for k in ("approach", "outcome", "artifact")
        ):
            raise ValueError("Each attempt requires approach, outcome and artifact")
    return context
