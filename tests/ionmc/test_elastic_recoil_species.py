"""V3-005C C6: the pseudo-species ``elastic_recoil`` (id 65), plan Amendment 17 (b)5.

The elastic p + A recoil is deposited with its own species id so that a species-scoring channel
resolves it in depth (rung 3 of the V5 attribution ladder). ``nuclear_local`` is the local-deposit
class and selects ids {64, 65}, so every ``nuclear_local`` channel is unchanged; the new channel is
a subset of it. Python and warp-cpu float64 agree on it; with ``elastic=False`` the request fails
closed (not producible)."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
import pytest

from ionmc.errors import UnsupportedCombinationError
from ionmc.scoring import TallyRequest
from ionmc.simulation import Simulation
from ionmc.species import ELASTIC_RECOIL_ID, N_SPECIES_SLOTS, PROTON, producible, species_by_name
from ionmc.transport.run import run_range
from ionmc.transport.tally import QUANTUM_MEV, tally_column_index
from tests.ionmc.test_elastic_transport import _cfg, _scale_elastic, tid  # noqa: F401

NUC = TallyRequest("nuc_local", "dose", "edep", species=("nuclear_local",))
REC = TallyRequest("recoil", "dose", "edep", species=("elastic_recoil",))
ALL = TallyRequest("all", "dose", "edep")


def _cfgt(table_id: str, n: int, tallies: tuple[TallyRequest, ...], **phys: Any) -> Any:
    return replace(_cfg(table_id, 100.0, n, **phys), tallies=tallies)


def _channel(eff: Any, part: Any, name: str) -> np.ndarray:
    plan = eff.channels
    c = plan.channels[{q.name: q.numerator for q in plan.quantities}[name]]
    return np.asarray(part.channel_acc[:, c.offset : c.offset + c.size], dtype=np.int64)


def _run(
    table_id: str, backend: str, tallies: tuple[TallyRequest, ...], n: int = 60
) -> tuple[Any, Any]:
    cfg = _cfgt(table_id, n, tallies)
    cfg = replace(cfg, run=replace(cfg.run, backend=backend, precision="float64"))
    eff = Simulation(cfg).effective
    return eff, run_range(eff, 0, n)


def test_registry_entry_and_width() -> None:
    assert species_by_name("elastic_recoil").id == ELASTIC_RECOIL_ID == 65 < N_SPECIES_SLOTS == 66
    assert species_by_name("elastic_recoil").transported is False
    on, off = producible(PROTON, True, True), producible(PROTON, True, False)
    assert ("elastic_recoil", "primary") in on and ("elastic_recoil", "secondary") in on
    assert not {p for p in off if p[0] == "elastic_recoil"}
    assert producible(PROTON, False, True) == producible(PROTON, False, False)


def test_request_fails_closed_without_elastic(tid: str) -> None:  # noqa: F811
    with pytest.raises(UnsupportedCombinationError):
        _ = Simulation(_cfgt(tid, 10, (REC,), elastic=False)).effective


@pytest.mark.parametrize("backend", ["python", "warp-cpu"])
def test_recoil_channel_is_the_elastic_recoil_subset(tid: str, backend: str) -> None:  # noqa: F811
    with pytest.MonkeyPatch.context() as mp:
        _scale_elastic(mp, 30.0)
        eff, part = _run(tid, backend, (NUC, REC, ALL))
    el = eff.nuclear is not None and eff.nuclear.elastic is not None
    assert el
    rec, nuc = _channel(eff, part, "recoil"), _channel(eff, part, "nuc_local")
    comps = part.tally_components
    ev = math_fsum(comps[tally_column_index(len(comps), True, "elastic_recoil_local")])
    assert ev > 0.0, "no elastic p + A event in the sample"
    assert abs(rec.sum() * QUANTUM_MEV - ev) <= 1e-9 * ev + 2 * QUANTUM_MEV * rec.size
    # the local-deposit class contains the elastic recoil: nuc_local >= recoil bin by bin
    assert np.all(nuc - rec >= 0)


def math_fsum(x: Any) -> float:
    import math

    return math.fsum(x)


def test_python_and_warp_cpu_agree_on_the_recoil_channel(tid: str) -> None:  # noqa: F811
    with pytest.MonkeyPatch.context() as mp:
        _scale_elastic(mp, 30.0)
        ep, pp = _run(tid, "python", (NUC, REC), 40)
        ew, pw = _run(tid, "warp-cpu", (NUC, REC), 40)
    for name in ("recoil", "nuc_local"):
        a, b = _channel(ep, pp, name), _channel(ew, pw, name)
        assert np.max(np.abs(a - b)) * QUANTUM_MEV <= 1e-10
    assert _channel(ep, pp, "recoil").sum() > 0


def test_nuclear_local_channel_is_unchanged_by_requesting_the_recoil(tid: str) -> None:  # noqa: F811
    with pytest.MonkeyPatch.context() as mp:
        _scale_elastic(mp, 30.0)
        e1, p1 = _run(tid, "python", (NUC, ALL), 40)
        e2, p2 = _run(tid, "python", (NUC, ALL, REC), 40)
    for name in ("nuc_local", "all"):
        assert np.array_equal(_channel(e1, p1, name), _channel(e2, p2, name))
    assert all(np.array_equal(a, b) for a, b in zip(p1.edep, p2.edep, strict=True))
