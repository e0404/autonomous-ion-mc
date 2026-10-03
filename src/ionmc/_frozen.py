"""Truly read-only numpy arrays: backed by an immutable ``bytes`` buffer.

An owning array that is only marked non-writeable can be made writable again with
``setflags(write=True)``; an array that views a ``bytes`` object cannot.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import DTypeLike, NDArray


def freeze_array(a: Any, dtype: DTypeLike, name: str = "array") -> NDArray[Any]:
    """Return a read-only copy of ``a`` (C-contiguous, ``dtype``) backed by immutable bytes.

    The caller's array is left untouched. Raises ``ValueError`` if the result is not frozen.
    """
    arr = np.asarray(a)
    raw = np.ascontiguousarray(arr, dtype=dtype).tobytes()
    frozen: NDArray[Any] = np.frombuffer(raw, dtype=dtype).reshape(arr.shape)
    root: Any = frozen
    while isinstance(root, np.ndarray):
        if root.flags.writeable:
            raise ValueError(f"{name} could not be made read-only")
        root = root.base
    if not isinstance(root, bytes):
        raise ValueError(f"{name} could not be backed by an immutable buffer")
    return frozen
