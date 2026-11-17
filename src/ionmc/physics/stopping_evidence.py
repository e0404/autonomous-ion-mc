"""Append-only evidence registry of the stopping models (decision 0044, section 1).

Evidence is kept apart from the immutable model spec: a grading appends one record keyed by
``spec_sha256`` and never changes a spec, a ``spec_sha256`` or a table hash. Statuses:
``unvalidated``, ``consistency-only``, ``topas-agreement-shared-lineage:<row>:PASS``,
``failed:<row>`` and ``as-graded:<decision>``. No status is named ``validated``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from ionmc.physics.stopping_models import STOPPING_MODELS, spec_sha256


@dataclass(frozen=True)
class EvidenceRecord:
    """One evidence record of the model with hash ``spec_sha256``."""

    spec_sha256: str
    status: str
    row: str | None
    plan: str
    graded_at_sha: str | None


def _initial() -> dict[str, tuple[EvidenceRecord, ...]]:
    ids = STOPPING_MODELS
    default = spec_sha256(ids["bethe-ionmc-v1"])
    out: dict[str, tuple[EvidenceRecord, ...]] = {
        default: (
            EvidenceRecord(
                default, "as-graded:0043", None, "validation/plans/v3-008-acceptance.md", None
            ),
        )
    }
    for mid in ("bethe-ionmc-v1+ls-eff", "bethe-ionmc-v1+ls-bare_mixed"):
        h = spec_sha256(ids[mid])
        out[h] = (
            EvidenceRecord(h, "unvalidated", None, "validation/plans/v3-013-acceptance.md", None),
        )
    return out


# Append-only by convention and by type: records are added in this module's source, never mutated.
STOPPING_EVIDENCE: Mapping[str, tuple[EvidenceRecord, ...]] = MappingProxyType(_initial())
