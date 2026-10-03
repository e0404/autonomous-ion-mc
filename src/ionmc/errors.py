"""Exception types of the transport engine (decision 0039).

``UnsupportedCombinationError`` is raised for every configuration the engine refuses
(fail-closed capability contract); it derives from ``ValueError`` so that callers
validating input may catch either. ``BackendUnavailableError`` signals a backend that is
not implemented or not usable on this machine; the engine never falls back to another
backend. ``TransportLimitError`` signals that a transport limit counter (step truncation,
stall, rejection limit, genealogy overflow, source energy out of range) is nonzero after a
run; the invalid result is attached as ``.result``.
"""

from __future__ import annotations

from typing import Any


class UnsupportedCombinationError(ValueError):
    """The requested configuration is invalid or not supported by the engine."""


class BackendUnavailableError(RuntimeError):
    """The requested backend is not implemented or has no usable device."""


class TransportLimitError(RuntimeError):
    """A transport limit counter is nonzero; ``result`` holds the invalid result (or None)."""

    def __init__(self, message: str, result: Any = None) -> None:
        super().__init__(message)
        self.result = result


class CounterOverflowError(OverflowError):
    """A random-stream counter field exceeded its documented bound (decision 0037)."""
