"""IonMC: Monte Carlo framework for ion therapy built on Python and NVIDIA Warp."""

from typing import Any

from ionmc._version import __version__

__all__ = ["__version__", "capabilities"]


def capabilities() -> dict[str, Any]:
    """What the transport engine supports (see :func:`ionmc.simulation.capabilities`).

    The simulation module (and with it Warp) is imported on first use.
    """
    from ionmc.simulation import capabilities as _capabilities

    return _capabilities()
