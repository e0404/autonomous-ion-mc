"""Offline LET from fluence spectra (helper of row A8): closed forms on synthetic spectra, edge
handling, and a small transport run whose online LET it must reproduce (CI size)."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace

import numpy as np
import pytest

from ionmc.config import SimulationConfig
from ionmc.geometry import BoxPhantom
from ionmc.let_offline import bin_centres, let_from_spectrum, log_edges
from ionmc.materials import WATER
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import Simulation


def test_log_edges_are_uniform_with_the_requested_density() -> None:
    e = log_edges(1.0, 200, 2.0)
    assert e.size == 401 and e[0] == 1.0 and e[-1] == pytest.approx(100.0, rel=1e-12)
    r = np.log(e[1:] / e[:-1])
    assert np.ptp(r) < 1e-12 and r[0] == pytest.approx(np.log(10.0) / 200, rel=1e-12)
    with pytest.raises(ValueError):
        log_edges(0.0, 200, 1.0)


def test_bin_centres_geometric_and_arithmetic() -> None:
    assert bin_centres([1.0, 4.0, 16.0]) == pytest.approx([2.0, 8.0])
    assert bin_centres([1.0, 3.0, 5.0], log_axis=False) == pytest.approx([2.0, 4.0])
    with pytest.raises(ValueError):
        bin_centres([2.0, 1.0])


def test_closed_form_two_bins_and_outside_weights() -> None:
    edges = [1.0, 2.0, 4.0]  # log centres sqrt(2), sqrt(8)
    s = lambda e: 10.0 / e  # noqa: E731
    sp = np.array([[0.5, 3.0, 1.0, 0.5], [0.0, 0.0, 0.0, 0.0]])
    r = let_from_spectrum(sp, edges, s)
    c = np.array([2.0**0.5, 8.0**0.5])
    sb = 10.0 / c
    w = np.array([3.0, 1.0])
    assert r.let_t[0] == pytest.approx((w * sb).sum() / w.sum(), rel=1e-14)
    assert r.let_d[0] == pytest.approx((w * sb * sb).sum() / (w * sb).sum(), rel=1e-14)
    assert r.outside_fraction[0] == pytest.approx(1.0 / 5.0)
    assert np.isnan(r.let_t[1]) and np.isnan(r.let_d[1]) and np.isnan(r.outside_fraction[1])
    # a nucleon number scales the energy at which S is evaluated
    r4 = let_from_spectrum(sp, edges, s, a_nucleon=4)
    assert r4.let_t[0] == pytest.approx(r.let_t[0] / 4.0, rel=1e-14)


def test_input_validation() -> None:
    s = lambda e: 1.0  # noqa: E731
    with pytest.raises(ValueError, match="entries"):
        let_from_spectrum(np.ones((1, 3)), [1.0, 2.0, 4.0], s)
    with pytest.raises(ValueError, match="non-negative"):
        let_from_spectrum(np.array([[0.0, -1.0, 1.0, 0.0]]), [1.0, 2.0, 4.0], s)
    with pytest.raises(ValueError, match="positive"):
        let_from_spectrum(np.ones((1, 4)), [1.0, 2.0, 4.0], lambda e: 0.0)


def test_reproduces_the_online_let_of_a_transport_run(
    make_config: Callable[..., SimulationConfig],
) -> None:
    edges = tuple(log_edges(1.0, 200, 1.4))  # 1 .. 25 MeV
    grid = ScoringGrid((-10.0, -10.0, 0.0), (20.0, 20.0, 1.0), (1, 1, 14), name="dose")
    cfg = replace(
        make_config(
            energy=20.0, n=40, n_batches=4, position=(0.0, 0.0, 0.0), scoring=(grid,),
            geometry=BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, 14.0), WATER), max_step=1.0,
        ),
        tallies=(
            TallyRequest("lt", "dose", "let_t"),
            TallyRequest("ld", "dose", "let_d"),
            TallyRequest("sp", "dose", "fluence_spectrum", energy_edges_mev_per_u=edges),
        ),
    )  # fmt: skip
    res = Simulation(cfg).run()
    q = res.grid("dose").quantities
    off = let_from_spectrum(
        q["sp"].mean.reshape(14, -1), edges, res.effective_config.tables.s_water
    )
    ok = q["lt"].defined_mask.reshape(-1) & q["ld"].defined_mask.reshape(-1)
    assert ok.sum() >= 5 and np.all(off.outside_fraction[ok] == 0.0)
    np.testing.assert_allclose(off.let_t[ok], q["lt"].mean.reshape(-1)[ok], rtol=5e-3)
    np.testing.assert_allclose(off.let_d[ok], q["ld"].mean.reshape(-1)[ok], rtol=5e-3)
