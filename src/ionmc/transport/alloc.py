"""The one instrumentable allocation point of the accumulator payload (plan section B.3, row MEM-1).

Every array of the accumulator-payload class that a transport backend allocates (the grid
accumulators ``edep``, the per-history tally and counter rows, the influence accumulators) is
created through :func:`alloc_zeros` and, when it is dropped before the end of the run, announced
through :func:`release`. Without an active recorder both are plain ``numpy`` calls; inside
``with record_allocations() as rec`` every call appends an event ``(kind, tag, nbytes)`` with
``kind`` in ``("alloc", "release")`` so that the audit can take the peak of the running sum
(amendment A-8: lifetimes across chunks are observed, a lingering array fails the audit).
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import DTypeLike, NDArray

__all__ = ["AllocationRecorder", "alloc_zeros", "record_allocations", "release"]


@dataclass
class AllocationRecorder:
    """Allocation and release events of the payload class, in call order."""

    events: list[tuple[str, str, int]] = field(default_factory=list)

    def running(self) -> list[int]:
        """The running sum of live bytes after every event."""
        total, out = 0, []
        for kind, _tag, nbytes in self.events:
            total += nbytes if kind == "alloc" else -nbytes
            out.append(total)
        return out

    @property
    def peak_bytes(self) -> int:
        """The peak of the running sum (0 without events)."""
        return max(self.running(), default=0)

    @property
    def live_bytes(self) -> int:
        """Bytes allocated and not released at the end of the recording."""
        r = self.running()
        return r[-1] if r else 0


_ACTIVE: AllocationRecorder | None = None


@contextmanager
def record_allocations() -> Iterator[AllocationRecorder]:
    """Record the payload allocations of the code run inside the ``with`` block."""
    global _ACTIVE
    previous, rec = _ACTIVE, AllocationRecorder()
    _ACTIVE = rec
    try:
        yield rec
    finally:
        _ACTIVE = previous


def alloc_zeros(tag: str, shape: tuple[int, ...], dtype: DTypeLike) -> NDArray[Any]:
    """``numpy.zeros`` of a payload array, recorded as an allocation of its ``nbytes``."""
    arr = np.zeros(shape, dtype=dtype)
    if _ACTIVE is not None:
        _ACTIVE.events.append(("alloc", tag, int(arr.nbytes)))
    return arr


def release(tag: str, arr: NDArray[Any]) -> None:
    """Announce that the caller drops the payload array ``arr`` now (recorded with its nbytes)."""
    if _ACTIVE is not None:
        _ACTIVE.events.append(("release", tag, int(arr.nbytes)))
