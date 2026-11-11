"""V3-008 R2-STOP3: window-cycle rule of ``core_sigma`` in
``validation/scripts/ion_em/r2_statistics.py`` (section 4.2, amendment R2-STOP3)."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

_PATH = Path(__file__).resolve().parents[2] / "validation" / "scripts" / "ion_em"
_spec = importlib.util.spec_from_file_location(
    "r2_statistics_cycle_under_test", _PATH / "r2_statistics.py"
)
assert _spec is not None and _spec.loader is not None
ST = importlib.util.module_from_spec(_spec)
sys.modules["r2_statistics_cycle_under_test"] = ST
_spec.loader.exec_module(ST)

_cdf = np.vectorize(lambda x: 0.5 * (1.0 + math.erf(x / math.sqrt(2.0))))
EDGES = 0.01 * np.arange(-3000, 3001)  # h = 0.01 mm, +-30 mm
CENTRES = 0.5 * (EDGES[:-1] + EDGES[1:])


def gauss(sigma: float, n: float = 1e6, edges: np.ndarray = EDGES) -> np.ndarray:
    return n * np.diff(_cdf(edges / sigma))


def old_core_sigma(w: np.ndarray, edges: np.ndarray) -> tuple[float, int]:
    """Frozen copy of the pre-amendment ``core_sigma`` (sigma, iterations)."""
    w = np.asarray(w, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.float64)
    c = 0.5 * (edges[:-1] + edges[1:])
    s0 = ST.second_moment_sigma(w, edges)
    s = s0
    for it in range(1, ST.FIT_MAX_ITER + 1):
        sel = np.nonzero(np.abs(c) <= ST.FIT_WINDOW_SIGMAS * s + ST._EDGE_EPS)[0]
        if sel.size < 3 or not np.sum(w[sel]) > 0.0 or np.any(np.diff(sel) != 1):
            raise ST.StatisticsError("fewer than 3 bins (or no signal) inside +-2 sigma")
        s_new = ST._window_ml(w[sel], edges[sel], edges[sel + 1], s)
        if abs(s_new - s) < ST.FIT_REL_TOL * s:
            return s_new, it
        s = s_new
    raise ST.StatisticsError(f"core fit did not converge in {ST.FIT_MAX_ITER} iterations")


def cycling_profile() -> np.ndarray:
    """Deterministic perturbation (found by search) of bins near |x| ~ 0.1-0.2 mm of a
    sigma = 0.0825 mm Gaussian for which the old rule 2-cycles."""
    rng = np.random.default_rng(23)
    p = gauss(0.0825)
    m = (np.abs(CENTRES) > 0.1) & (np.abs(CENTRES) < 0.2)
    p[m] *= 1.0 + 0.02 * rng.standard_normal(int(m.sum()))
    return p


def test_t1_two_cycle_old_fails_new_returns_mean() -> None:
    p = cycling_profile()
    with pytest.raises(ST.StatisticsError, match="did not converge"):
        old_core_sigma(p, EDGES)
    fit = ST.core_sigma(p, EDGES)
    assert fit.cycle_period == 2 and 0.0 <= fit.cycle_amplitude < 0.02
    assert fit.sigma == pytest.approx(0.0825, rel=0.02)
    mirrored = ST.core_sigma(p[::-1], EDGES)
    # x -> -x only reorders floating-point sums: equal to rounding, same cycle
    assert mirrored.sigma == pytest.approx(fit.sigma, rel=1e-7) and mirrored.cycle_period == 2
    assert mirrored.iterations == fit.iterations
    again = ST.core_sigma(p, EDGES)
    assert (again.sigma, again.iterations, again.cycle_amplitude) == (
        fit.sigma, fit.iterations, fit.cycle_amplitude,
    )  # fmt: skip
    # the cycle mean is the mean of the two ML scales of the two windows
    assert ST.core_xy(p, p, EDGES, EDGES) == fit.sigma
    assert ST.core_xy_detail(p, p[::-1], EDGES, EDGES)[1] == 2


def test_t2_bit_identity_on_converging_profiles() -> None:
    rng = np.random.default_rng(2024)
    n = 0
    for i in range(600):
        sigma = float(rng.uniform(0.06, 0.4))
        shift = float(rng.uniform(0.0, 0.01))  # bin-phase of the profile relative to the grid
        edges = EDGES + shift
        p = gauss(sigma, 1e6, edges)
        if i % 2:
            p = p + gauss(float(rng.uniform(3.0, 9.0)), 2e4, edges)
        if i % 3 == 0:
            p = p + 5.0  # flat pedestal
        try:
            sig_old, it_old = old_core_sigma(p, edges)
        except ST.StatisticsError:
            continue
        fit = ST.core_sigma(p, edges)
        assert fit.sigma == sig_old and fit.iterations == it_old and fit.cycle_period == 0
        n += 1
    assert n >= 200


def _fake_ml(monkeypatch: pytest.MonkeyPatch, scales: list[float], s_start: float) -> None:
    """``_window_ml`` returns ``scales[k % len]`` at call k (a pure function of the call index, so
    the windows visited follow the scales); the starting scale is ``s_start``."""
    calls = {"n": 0}

    def fake(w, lo, hi, s0):  # noqa: ANN001, ARG001
        k = calls["n"]
        calls["n"] += 1
        return scales[k % len(scales)]

    monkeypatch.setattr(ST, "_window_ml", fake)
    monkeypatch.setattr(ST, "second_moment_sigma", lambda w, e: s_start)


def test_t3_large_amplitude_cycle_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    # s alternates 0.08 <-> 0.09: windows |x| <= 0.16 and |x| <= 0.18 (h = 0.01 mm)
    _fake_ml(monkeypatch, [0.09, 0.08], 0.08)
    with pytest.raises(ST.StatisticsError, match=r"window cycle amplitude .* > 0\.02"):
        ST.core_sigma(gauss(0.0825), EDGES)


def test_t3b_small_amplitude_two_cycle_returns_mean(monkeypatch: pytest.MonkeyPatch) -> None:
    # 2 s = 0.1746 / 0.1754 straddles the bin centre 0.175: two different windows
    _fake_ml(monkeypatch, [0.0877, 0.0873], 0.0873)
    fit = ST.core_sigma(gauss(0.0825), EDGES)
    assert fit.cycle_period == 2 and fit.sigma == pytest.approx(0.0875)
    assert fit.cycle_amplitude == pytest.approx(0.0004 / 0.0875)


def test_t4_period_three_cycle(monkeypatch: pytest.MonkeyPatch) -> None:
    # 2 s = 0.5952 / 0.6000 / 0.6052 select three different windows (bin centres at 0.595, 0.605)
    cyc = [0.3000, 0.3026, 0.2976]
    _fake_ml(monkeypatch, cyc, 0.2976)
    fit = ST.core_sigma(gauss(0.3), EDGES)
    assert fit.cycle_period == 3
    assert fit.sigma == pytest.approx(sum(cyc) / 3)
    assert fit.cycle_amplitude == pytest.approx((0.3026 - 0.2976) / (sum(cyc) / 3))


def test_t5_non_repeating_non_convergence_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    state = {"n": 0}

    def growing(w, lo, hi, s0):  # noqa: ANN001, ARG001
        state["n"] += 1
        return s0 * 1.05  # every window is new until the cap

    monkeypatch.setattr(ST, "_window_ml", growing)
    monkeypatch.setattr(ST, "second_moment_sigma", lambda w, e: 0.05)
    with pytest.raises(ST.StatisticsError, match="did not converge"):
        ST.core_sigma(gauss(0.0825), EDGES)


def test_t6_continuity_across_bin_edge_thresholds() -> None:
    # 2 sigma crosses a bin centre at sigma = (k + 1/2) h / 2 ... sweep finely around several
    for centre in (0.0825, 0.0875, 0.0925, 0.1025):
        for sigma in centre + np.linspace(-0.0008, 0.0008, 17):
            sigma = float(sigma)
            fit = ST.core_sigma(gauss(sigma), EDGES)
            assert abs(fit.sigma / sigma - 1.0) < 1e-3, (sigma, fit)


def test_t0_cycle_fields_default() -> None:
    f = ST.CoreFit(1.0, 1.0, 3)
    assert f.cycle_period == 0 and f.cycle_amplitude == 0.0
