"""Backend dispatch: run the transport of a validated configuration and reduce it.

``python``, ``warp-cpu`` and ``warp-cuda`` all produce per-history exact tallies that are
reduced by :func:`ionmc.transport.tally.merge_partials`. ``cpu_workers > 1`` splits the
histories over spawned worker processes (:mod:`ionmc.transport.pool`).
"""

from __future__ import annotations

from ionmc.config import EffectiveConfig
from ionmc.transport.tally import PartialTransport, RawTransport, build_diagnostics, merge_partials


def run_range(eff: EffectiveConfig, h0: int, h1: int) -> PartialTransport:
    """Transport histories ``[h0, h1)`` in this process on the backend of ``eff``."""
    if eff.backend == "python":
        from ionmc.transport.reference import run_reference_range

        return run_reference_range(eff, h0, h1)
    from ionmc.transport.warp_driver import device_for_backend, run_warp_range

    return run_warp_range(eff, h0, h1, device_for_backend(eff.backend))


def run_transport(eff: EffectiveConfig) -> RawTransport:
    """Run all histories of ``eff`` and return the reduced raw result."""
    run = eff.requested.run
    n = run.n_histories
    if run.cpu_workers > 1:
        from ionmc.transport.pool import run_pool

        parts = run_pool(eff)
    else:
        parts = [run_range(eff, 0, n)]
    diag = eff.requested.diagnostics
    raw = merge_partials(parts, n, len(eff.requested.scoring))
    raw.diagnostics = build_diagnostics(
        parts, diag.track_end_positions, diag.escape_records, diag.trace_histories
    )
    raw.meta = {
        "workers": run.cpu_workers,
        "partials": [dict(p.meta) for p in sorted(parts, key=lambda q: q.h0)],
    }
    return raw
