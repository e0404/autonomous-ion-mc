"""Backend-parity measurements of the scoring channels (rows A11 and A15 of
``validation/plans/v3-004-acceptance.md``; the frozen T12 machinery of
:mod:`ionmc.transport.parity` is reused unchanged).

* :func:`channel_t12_observables` builds, per batch and per primary, the **linear** depth profiles
  of
  the channels ``L``, ``LS``, ``LS2``, ``ES``, ``E_step`` and ``FE`` (a synthetic lookup, when
  requested), the depth-integrated fluence spectrum, and the channel totals of ``L``, ``LS``,
  ``LS2`` and ``ES`` as scalars. Linear channels are compared because the frozen T12 grouping and
  batch permutation are calibrated for additive batch sums (a merged group of ratio bins is not the
  ratio of the merged sums, and per-batch ratios carry an O(1/n_b) Jensen bias); equality in
  distribution of numerator and denominator implies parity of the ratio.
* :func:`compare_channel_runs` applies :func:`~ionmc.transport.parity.t12_compare` to them (profiles
  at p > 0.001 and max-T p > 0.001, inconclusive fails; scalars ``|z| < 3.5``; the deterministic
  total-deposit rule is *not* used) and reports the delta-method z of the ``LET_t`` and ``LET_d``
  ratio profiles (:func:`ratio_z_report`), which is informative and never gates.
* :func:`compare_channel_partition` (A15) requires bit-identical accumulators (N and N_local
  included), residuals, out-of-domain counts, counters and deposits of two runs that differ only
  in the partition.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.scoring import reduce_ratio
from ionmc.simulation import Result
from ionmc.transport.channels import ChannelPlan
from ionmc.transport.parity import T12Observables, t12_compare

PROFILE_QUANTITIES = {
    "let_t": ("LS", "L"),
    "let_d": ("LS2", "LS"),
    "let_d_eps": ("ES", "E_step"),
}
"""Ratio quantity -> (numerator, denominator) channel labels."""


def _quantity(plan: ChannelPlan, quantity: str, grid: str) -> Any:
    """First quantity definition of kind ``quantity`` on grid ``grid`` (None if absent)."""
    for q in plan.quantities:
        if q.request.quantity == quantity and q.request.grid == grid:
            return q
    return None


def _batch_values(result: Result, ci: int) -> NDArray[np.float64]:
    """Per-batch, per-primary values ``[B, size]`` of channel ``ci`` in channel units."""
    plan = result.effective_config.channels
    assert plan is not None
    hpb = result.n_histories // result.n_batches
    q = plan.channels[ci].quantum
    return result.channel_batches(ci).astype(np.float64) * q / hpb


def piece_counts_for(result: Result, ci: int) -> NDArray[np.float64]:
    """Per-batch piece counts ``[B, n_voxels]`` that bound the rounding of channel ``ci``: the
    step count N, the local count N_local or their sum according to the channel's class mask
    (``ChannelPlan.piece_count_indices``); a spectrum channel's bins share its voxel counts."""
    plan = result.effective_config.channels
    assert plan is not None
    out: NDArray[np.float64] | None = None
    for ni in plan.piece_count_indices(ci):
        v = result.channel_batches(ni).astype(np.float64) * plan.channels[ni].quantum
        out = v if out is None else out + v
    assert out is not None
    return out


def channel_labels(result: Result, grid: str) -> dict[str, int]:
    """Channel index of each label (``L``, ``LS``, ``LS2``, ``ES``, ``E_step``, ``FE``,
    ``spectrum``) that the run has on ``grid`` (labels without a request are absent)."""
    plan = result.effective_config.channels
    if plan is None:
        raise ValueError("this result has no scoring channels")
    out: dict[str, int] = {}
    for qname, (num, den) in PROFILE_QUANTITIES.items():
        q = _quantity(plan, qname, grid)
        if q is not None:
            out[num], out[den] = q.numerator, q.denominator
    for quantity, label in (("lookup_sum", "FE"), ("fluence_spectrum", "spectrum")):
        q = _quantity(plan, quantity, grid)
        if q is not None:
            out[label] = q.numerator
    if "L" not in out:
        q = _quantity(plan, "fluence", grid)
        if q is not None:
            out["L"] = q.numerator
    return out


def channel_t12_observables(result: Result, grid: str = "idd") -> T12Observables:
    """Batch-wise linear channel profiles (depth bins of ``grid``) and totals of one run (see the
    module docstring)."""
    labels = channel_labels(result, grid)
    arrays: dict[str, NDArray[np.float64]] = {}
    scalars: dict[str, NDArray[np.float64]] = {}
    for label, ci in labels.items():
        v = _batch_values(result, ci)
        if label == "spectrum":
            nvox = next(
                g for g in result.effective_config.requested.scoring if g.name == grid
            ).n_voxels
            v = v.reshape(v.shape[0], nvox, -1).sum(axis=1)  # depth-integrated bins
        arrays[label] = v
        if label in ("L", "LS", "LS2", "ES"):
            scalars[f"{label}_total"] = v.sum(axis=1)
    return T12Observables(
        arrays,
        scalars,
        result.n_batches,
        result.n_histories // result.n_batches,
        result.precision,
    )


def ratio_z_report(a: Result, b: Result, grid: str = "idd") -> dict[str, Any]:
    """Delta-method z of the ``LET_t`` and ``LET_d`` (and ``LET_d^eps``) profiles of two
    independent samples, over the bins defined in both: ``z = (R_a - R_b) / sqrt(V_a + V_b)``.
    Informative only (non-gating)."""
    la, lb = channel_labels(a, grid), channel_labels(b, grid)
    out: dict[str, Any] = {}
    for qname, (num, den) in PROFILE_QUANTITIES.items():
        if num not in la or den not in la or num not in lb or den not in lb:
            continue
        sa = reduce_ratio(_batch_values(a, la[num]), _batch_values(a, la[den]))
        sb = reduce_ratio(_batch_values(b, lb[num]), _batch_values(b, lb[den]))
        ok = sa.defined_mask & sb.defined_mask
        var = sa.variance_of_mean + sb.variance_of_mean
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (sa.mean - sb.mean) / np.sqrt(var)
        z = np.where(ok & (var > 0.0), z, np.nan)
        finite = z[np.isfinite(z)]
        out[qname] = {
            "n_defined_both": int(ok.sum()),
            "max_abs_z": float(np.abs(finite).max()) if finite.size else None,
            "rms_z": float(np.sqrt(np.mean(finite**2))) if finite.size else None,
            "gating": False,
        }
    return out


def compare_channel_runs(
    a: Result, b: Result, grid: str = "idd", **t12_kwargs: Any
) -> dict[str, Any]:
    """A11 (HR) verdict for two independent samples: ``t12_compare`` of the linear channel
    observables (gating) plus the reported ratio z (non-gating). ``pass`` is the T12 verdict."""
    t12 = t12_compare(
        channel_t12_observables(a, grid), channel_t12_observables(b, grid), **t12_kwargs
    )
    return {"t12": t12, "ratio_z": ratio_z_report(a, b, grid), "pass": bool(t12["pass"])}


def compare_channel_partition(a: Result, b: Result) -> dict[str, Any]:
    """A15: accumulators, residual tallies, out-of-domain count, counters and the qualified
    deposit grids of two runs of the same configuration but different partitions (chunk size or
    worker count) must be bit-identical."""
    ca, cb = a.channel_raw, b.channel_raw
    if ca is None or cb is None:
        raise ValueError("both results need scoring channels")
    pa, pb = a.effective_config.channels, b.effective_config.channels
    assert pa is not None and pb is not None
    counts = [i for i, c in enumerate(pa.channels) if c.kind == "N"]  # N (step) and N_local
    checks = {
        "acc": bool(np.array_equal(ca.acc, cb.acc)),
        "piece_counts": all(
            np.array_equal(a.channel_batches(i), b.channel_batches(i)) for i in counts
        )
        and len(counts) >= 2,
        "residual": bool(np.array_equal(ca.residual, cb.residual)),
        "lookup_out_of_domain": ca.lookup_out_of_domain == cb.lookup_out_of_domain,
        "path_bound_exceeded": ca.path_bound_exceeded == cb.path_bound_exceeded,
        "counters": a.counters == b.counters,
        "edep": all(
            np.array_equal(ga.batch_energy_mev, gb.batch_energy_mev)
            for ga, gb in zip(a.grids, b.grids, strict=True)
        ),
    }
    return {"checks": checks, "pass": all(checks.values())}
