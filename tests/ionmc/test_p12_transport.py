"""Regression coverage of the P12 transport reproducer (validation/plans/v3-008/p12-scripts/
transport.py; Codex REVIEW-c43a2e52): the cutoff point deposit is normalised exactly once, a
double-normalised term fails the check, and the reported peak-to-entrance / range changes of
the band-average profile at the standard points are pinned."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "validation/plans/v3-008/p12-scripts/transport.py"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("p12_transport", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    sys.modules["p12_transport"] = m
    spec.loader.exec_module(m)
    return m


def test_point_deposit_normalised_once(mod):
    x = np.arange(0, 200, 0.02)
    term = mod.cutoff_point_deposit(x, 100.0, 0.4, 12)
    integ = mod.check_point_deposit(x, term, 12)
    assert integ == pytest.approx(2.0 * 12, rel=1e-4)
    with pytest.raises(RuntimeError):  # the former double normalisation fails closed
        mod.check_point_deposit(x, term / math.sqrt(2 * math.pi), 12)


def test_band_profile_pinned_at_the_standard_points(mod):
    names = {"band-avg profile"}
    _, c = mod.run_case(mod.P.CARBON12, 290, 0.10, names)
    _, o = mod.run_case(mod.P.OXYGEN16, 350, 0.135, names)
    c, o = c["band-avg profile"], o["band-avg profile"]
    assert c["pe_change_pct"] == pytest.approx(1.80, abs=0.15)
    assert o["pe_change_pct"] == pytest.approx(2.77, abs=0.15)
    assert c["dR_mm"] == pytest.approx(-0.058, abs=0.01)
    assert o["dR_mm"] == pytest.approx(-0.092, abs=0.01)
    assert 0.4 < c["sigma_r_mm"] < 0.5 and 0.35 < o["sigma_r_mm"] < 0.45
