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
    for mid in ("nist-pstar-astar-water-2005", "icru90-water-v1"):
        h = spec_sha256(ids[mid])
        out[h] = (
            EvidenceRecord(
                h, "consistency-only", None, "validation/plans/v3-013-acceptance.md", None
            ),
        )
    return out


# Append-only by convention and by type: records are added in this module's source, never mutated.
_RECORDS: dict[str, tuple[EvidenceRecord, ...]] = _initial()
STOPPING_EVIDENCE: Mapping[str, tuple[EvidenceRecord, ...]] = MappingProxyType(_RECORDS)


# Dataset ids that are also used as evaluation references of the model's rows (decision 0044
# section 1: the NIST model is circular for the PSTAR/ASTAR rows), keyed by spec_sha256. For
# ``icru90-water-v1`` the dataset id is followed by the names of the rows that use the ICRU 90
# arrays as an evaluation reference (row PER-2): the ICRU 90 part of I2 of plan v3-008 and the
# ICRU 90 comparisons of decision 0038; all of them are circular for a run of that model.
STOPPING_CIRCULAR_ROWS: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        spec_sha256(STOPPING_MODELS["nist-pstar-astar-water-2005"]): (
            "nist-pstar-water-2005",
            "nist-astar-water-2005",
        ),
        spec_sha256(STOPPING_MODELS["icru90-water-v1"]): (
            "geant4-icru90-stopping-11.4.2",
            "I2 (ICRU 90 part)",
            "decision-0038 ICRU 90 comparisons",
        ),
    }
)


def add_test_record(record: EvidenceRecord) -> None:
    """TEST ONLY (``register_test_model``): add the record of a temporarily registered model."""
    if record.spec_sha256 in _RECORDS:
        raise ValueError(f"evidence for {record.spec_sha256} already exists")
    _RECORDS[record.spec_sha256] = (record,)


def remove_test_record(digest: str) -> None:
    """TEST ONLY (``register_test_model``): remove the record added by :func:`add_test_record`."""
    _RECORDS.pop(digest, None)
