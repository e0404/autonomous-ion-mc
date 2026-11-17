"""Helper for test stopping sources that declare themselves (decision 0044 section 2, DOM-7).

A ``stopping`` object of an unsupported type is refused; the one documented hook is a
``stopping_spec()`` method returning a :class:`StoppingModelSpec` (decision ``"user"``). Test
sources use :func:`user_spec` or the :class:`UserSpecMixin`; the spec declares no window narrower
than the source's own ``e_min_per_u``/``e_max_per_u`` attributes (when it has them) and otherwise
0 to the largest float per species, so that the behaviour of the source is unchanged.
"""

from __future__ import annotations

import sys
from typing import Any

from ionmc.physics.stopping_models import MATERIAL_RULE_ANY, StoppingDomain, StoppingModelSpec

SPECIES = ("proton", "deuteron", "alpha", "carbon12", "oxygen16")
UNBOUNDED = sys.float_info.max


def user_spec(source: Any, kind: str = "analytic", **options: Any) -> StoppingModelSpec:
    """The ``StoppingModelSpec`` of a test source (``decision = "user"``, ``data = ()``)."""
    lo = float(getattr(source, "e_min_per_u", 0.0))
    hi = float(getattr(source, "e_max_per_u", UNBOUNDED))
    return StoppingModelSpec(
        "unregistered",
        "1",
        kind,
        {"source_name": str(getattr(source, "name", type(source).__name__)), **options},
        "test source",
        "test source",
        (),
        StoppingDomain(
            frozenset(SPECIES),
            MATERIAL_RULE_ANY,
            dict.fromkeys(SPECIES, lo),
            dict.fromkeys(SPECIES, hi),
        ),
        "user",
    )


class UserSpecMixin:
    """Gives a test source ``stopping_spec()`` (kind ``analytic``)."""

    def stopping_spec(self) -> StoppingModelSpec:
        return user_spec(self)
