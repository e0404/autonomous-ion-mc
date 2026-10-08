"""V3-005B C11: the nuclear Warp kernel variant against the Python reference (row V8 CI).

K = 8 histories at 150 MeV and the nuclear-dense configuration (every Sigma x 40, as in the Python
loop tests): the float64 warp-cpu kernel follows the Python trajectory (identical events,
species, genealogy ids; continuous columns <= 1e-10), closes the energy balance, fires the
fail-closed counters on forced inputs, is bitwise chunk-invariant, and the float32 variant stays
within the 0039 class tolerances (maxima reported). Data-backed (``IONMC_CACHE_DIR``)."""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Any

import numpy as np
import pytest

import ionmc.transport.kernels_nuclear as kn
import ionmc.transport.nuclear_device as nuclear_device
from ionmc.config import DiagnosticsOptions
from ionmc.simulation import Simulation
from ionmc.transport.reference import run_reference_range
from ionmc.transport.tally import COUNTER_NAMES, NUCLEAR_COUNTER_NAMES
from ionmc.transport.warp_driver import run_warp_range
from tests.ionmc.test_nuclear_transport_loop import _config, _scaled_rows, _table_id

K = 8
TOL = 1e-10


@pytest.fixture(scope="module")
def tid() -> str:
    return _table_id()


def _eff(tid: str, backend: str, precision: str, n: int = K, energy: float = 150.0,
         trace: bool = True, **kw: Any) -> Any:  # fmt: skip
    cfg = _config(tid, energy=energy, n=n, track_end=True, **kw)
    cfg = replace(
        cfg,
        run=replace(cfg.run, backend=backend, precision=precision),
        diagnostics=DiagnosticsOptions(track_end_positions=True, trace_histories=n if trace else 0),
    )
    return Simulation(cfg).effective


def _sums(part: Any) -> tuple[list[float], list[int]]:
    return [math.fsum(c) for c in part.tally_components], list(part.counter_sums)


def _compare(py: Any, wp_: Any) -> dict[str, float]:
    """Parity of two partials; returns the maxima; asserts identical discrete columns."""
    out: dict[str, float] = {}
    dp, dw = py.diagnostics, wp_.diagnostics
    assert dp is not None and dw is not None
    np.testing.assert_array_equal(dp.trace_int, dw.trace_int)  # steps, voxels, reasons, blocks
    np.testing.assert_array_equal(dp.end_code, dw.end_code)
    out["trace"] = float(np.max(np.abs(dp.trace_float - dw.trace_float), initial=0.0))
    out["end"] = max(
        float(np.max(np.abs(dp.end_position_mm - dw.end_position_mm))),
        float(np.max(np.abs(dp.end_energy_mev - dw.end_energy_mev))),
    )
    tp, tw = py.meta["nuclear_trace"], wp_.meta["nuclear_trace"]
    for key in ("events", "secondaries"):
        a, b = tp[key], tw[key]
        assert a.shape == b.shape, (key, a.shape, b.shape)
        disc = [0, 1, 2] if key == "events" else [0, 1, 2, 3, 4]
        if key == "events":
            disc += list(range(3, 11))  # counts, Z_r, A_r, attempts
        np.testing.assert_array_equal(a[:, disc], b[:, disc], err_msg=key)  # ids, species, gen
        cont = [11] if key == "events" else list(range(5, 12))
        out[key] = float(np.max(np.abs(a[:, cont] - b[:, cont]), initial=0.0))
    sp_, cp = _sums(py)
    sw, cw = _sums(wp_)
    assert cp == cw, (cp, cw)
    out["tallies"] = max(abs(x - y) for x, y in zip(sp_, sw, strict=True))
    out["edep"] = max(
        float(np.max(np.abs(a - b))) for a, b in zip(py.edep, wp_.edep, strict=True)
    ) * 2.0**-30  # fmt: skip
    return out


@pytest.mark.parametrize("dense", [False, True])
def test_v8_ci_trajectory_parity_float64(
    tid: str, monkeypatch: pytest.MonkeyPatch, dense: bool
) -> None:
    n = 64 if dense else K
    if dense:
        _scaled_rows(monkeypatch, 40.0)
    ep = _eff(tid, "python", "float64", n=n, energy=100.0 if dense else 150.0)
    ew = _eff(tid, "warp-cpu", "float64", n=n, energy=100.0 if dense else 150.0)
    py = run_reference_range(ep, 0, n)
    wr = run_warp_range(ew, 0, n, "cpu")
    m = _compare(py, wr)
    print(
        "V8 parity",
        "dense" if dense else "ci",
        m,
        "events",
        len(wr.meta["nuclear_trace"]["events"]),
    )
    if dense:
        assert len(wr.meta["nuclear_trace"]["events"]) > 0
        assert len(wr.meta["nuclear_trace"]["secondaries"]) > 0
    assert max(m.values()) <= TOL, m


def test_energy_balance_and_chunk_invariance_warp_cpu(
    tid: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    _scaled_rows(monkeypatch, 40.0)
    n = 40
    cfg = _config(tid, energy=100.0, n=n, track_end=False)
    cfg = replace(cfg, run=replace(cfg.run, backend="warp-cpu", precision="float64"))
    res = Simulation(cfg).run()
    b = res.energy_balance
    print("warp-cpu nuclear closure", b.relative_residual, res.counters.as_dict())
    assert res.valid and b.relative_residual <= 1e-12
    assert res.energy_balance.nuclear_mev["nuclear_local"] > 0.0
    ref = None
    for chunk in (1, 4, 8):
        eff = Simulation(cfg).effective
        object.__setattr__(eff.requested.run, "chunk_histories", chunk)
        part = run_warp_range(eff, 0, 16, "cpu")
        assert part.meta["n_chunks"] == 16 // chunk
        cur = (_sums(part), [e.tobytes() for e in part.edep])  # exact column sums
        if ref is None:
            ref = cur
        assert cur == ref, f"chunk {chunk} differs"


def _forced(tid: str, monkeypatch: pytest.MonkeyPatch) -> Any:
    _scaled_rows(monkeypatch, 40.0)
    cfg = _config(tid, energy=100.0, n=32, track_end=False, allow_invalid=True)
    cfg = replace(cfg, run=replace(cfg.run, backend="warp-cpu", precision="float64"))
    return Simulation(cfg).run()


def _closed(res: Any, counter: str) -> None:
    assert res.counters.as_dict()[counter] > 0 and not res.valid
    assert res.energy_balance.relative_residual <= 1e-12, res.energy_balance.relative_residual


@pytest.mark.parametrize(
    "counter", ["queue_overflow", "genealogy_overflow", "nuclear_conservation"]
)
def test_forced_runtime_limits_fail_closed_warp_cpu(
    tid: str, monkeypatch: pytest.MonkeyPatch, counter: str
) -> None:
    if counter == "queue_overflow":
        monkeypatch.setattr(kn, "STACK_CAPACITY", 0)
    elif counter == "genealogy_overflow":
        monkeypatch.setattr(kn, "CHILD_LIMIT", 0)
    else:
        monkeypatch.setattr(kn, "LEDGER_TOL", -1.0)
    _closed(_forced(tid, monkeypatch), counter)


@pytest.mark.parametrize("counter", ["nuclear_rejection_limit", "majorant_violation"])
def test_forced_table_inputs_fail_closed_warp_cpu(
    tid: str, monkeypatch: pytest.MonkeyPatch, counter: str
) -> None:
    orig = nuclear_device.pack_nuclear

    def packed(*a: Any, **k: Any) -> Any:
        h = orig(*a, **k)
        if counter == "nuclear_rejection_limit":
            h.arrays["m_res"] = np.full_like(h.arrays["m_res"], np.inf)  # no residual exists
        else:
            h.arrays["sigma"] = h.arrays["sigma"] * 1.0e9  # Sigma(E1) > S^(E0)
        return h

    monkeypatch.setattr(nuclear_device, "pack_nuclear", packed)
    _closed(_forced(tid, monkeypatch), counter)


def test_float32_parity_within_class_tolerance(tid: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """float32 against the float64 reference: histories whose events agree are compared; the
    maxima are reported (decision 0039 classes: momenta 5e-5, mu 5e-3, ledger 0.1 of the event
    energy)."""
    _scaled_rows(monkeypatch, 40.0)
    n = 64
    e64 = _eff(tid, "warp-cpu", "float64", n=n, energy=100.0, trace=False)
    e32 = _eff(tid, "warp-cpu", "float32", n=n, energy=100.0, trace=False)
    p64 = run_warp_range(e64, 0, n, "cpu")
    p32 = run_warp_range(e32, 0, n, "cpu")
    s64, c64 = _sums(p64)
    s32, c32 = _sums(p32)
    names = COUNTER_NAMES + NUCLEAR_COUNTER_NAMES
    print("float32 counters", dict(zip(names, c32, strict=True)))
    initial = s64[0]
    rel = max(abs(a - b) for a, b in zip(s64[:6], s32[:6], strict=True)) / initial
    print("float32 vs float64 max tally diff / initial", rel)
    assert c32[names.index("nuclear_conservation")] == 0
    assert rel <= 1e-3
