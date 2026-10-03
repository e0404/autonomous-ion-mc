"""Lateral-profile estimator and TOPAS binary parser on synthetic data."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ionmc.reference import ParseError, parse_topas_binary
from ionmc.reference.metrics import (
    bin_centres,
    lateral_sigma2_at_depth_fractions,
    lateral_sigma2_slab,
    lateral_variance_1d,
)

_ERF = np.vectorize(math.erf)


def gauss_binned(n: int, h: float, sigma: float, centre: float = 0.0) -> np.ndarray:
    """Bin-integrated unit Gaussian on ``n`` bins of width ``h`` centred on the grid centre."""
    edges = (np.arange(n + 1) - 0.5 * n) * h
    cdf = 0.5 * (1.0 + _ERF((edges - centre) / (sigma * math.sqrt(2.0))))
    return np.diff(cdf)


@pytest.mark.parametrize("sigma", [2.0, 3.0, 4.0, 5.0, 6.0])
def test_sigma_recovered_within_half_percent(sigma: float) -> None:
    h, n = 0.5, 240
    x = bin_centres(n, h)
    p = gauss_binned(n, h, sigma) * 1e6 + 3.0  # constant background
    v = lateral_variance_1d(x, p, h)
    assert math.sqrt(v) == pytest.approx(sigma, rel=5e-3)
    if sigma <= 4.0:  # a narrow window needs its background zone (0.9 W) >= 4.5 sigma
        v20 = lateral_variance_1d(x, p, h, window_mm=20.0)
        assert math.sqrt(v20) == pytest.approx(sigma, rel=5e-3)


def test_sheppard_correction_needed_and_offset_beam() -> None:
    h, n, sigma = 2.0, 60, 3.0  # coarse bins: the correction is large
    x = bin_centres(n, h)
    p = gauss_binned(n, h, sigma, centre=1.3)
    v = lateral_variance_1d(x, p, h)
    assert v == pytest.approx(sigma**2, rel=1e-6)
    assert v + h * h / 12 > 1.03 * sigma**2  # without Sheppard the bias would be material


def test_point_sampled_gaussian() -> None:
    h, n, sigma = 0.5, 240, 2.0
    x = bin_centres(n, h)
    p = np.exp(-0.5 * (x / sigma) ** 2)
    assert math.sqrt(lateral_variance_1d(x, p, h)) == pytest.approx(sigma, rel=5e-3)


def test_slab_projections_and_depth_fractions() -> None:
    h, n = 0.5, 120
    gx = gauss_binned(n, h, 3.0)
    gy = gauss_binned(n, h, 4.0)
    slab = np.outer(gx, gy)
    res = lateral_sigma2_slab(slab, (h, h))
    assert res["sigma2_x_mm2"] == pytest.approx(9.0, rel=1e-5)
    assert res["sigma2_y_mm2"] == pytest.approx(16.0, rel=1e-5)
    assert res["sigma2_mm2"] == pytest.approx(12.5, rel=1e-5)
    d3 = np.stack([slab * (1 + z) for z in range(10)], axis=2)
    out = lateral_sigma2_at_depth_fractions(d3, (h, h, 1.0), 8.0, (0.5, 0.9))
    assert out[0.5]["slab_centre_mm"] == 4.5 and out[0.9]["slab_centre_mm"] == 7.5
    assert out[0.9]["slab_centre_over_r80"] == pytest.approx(7.5 / 8.0)
    assert out[0.9]["sigma2_mm2"] == pytest.approx(12.5, rel=1e-5)
    with pytest.raises(ValueError):
        lateral_sigma2_at_depth_fractions(d3, (h, h, 1.0), 20.0)


def test_estimator_errors() -> None:
    x = bin_centres(40, 0.5)
    with pytest.raises(ValueError):
        lateral_variance_1d(x, np.zeros(40), 0.5)  # no signal
    with pytest.raises(ValueError):
        lateral_variance_1d(x, np.ones(39), 0.5)
    with pytest.raises(ValueError):
        lateral_variance_1d(x, np.full(40, np.nan), 0.5)
    with pytest.raises(ValueError):
        lateral_variance_1d(x, np.ones(40), 0.5, window_mm=1.0)


HEADER = (
    "# X in 2 bins of 0.5 mm\n# Y in 3 bins of 0.5 mm\n# Z in 4 bins of 1 mm\n"
    "# DoseToMedium ( Gy ) : Sum\n"
)


def test_topas_binary_layouts_and_precision() -> None:
    arr = np.arange(24, dtype="<f8")
    s = parse_topas_binary(HEADER, arr.tobytes())
    assert s.bins == (2, 3, 4) and s.statistics == ["Sum"] and s.meta["dtype"] == "<f8"
    assert s.values["Sum"][1, 2, 3] == 23 and s.values["Sum"][0, 0, 1] == 1  # z fastest
    f = parse_topas_binary(HEADER, arr.astype("<f4").tobytes(), order="F")
    assert f.meta["dtype"] == "<f4" and f.values["Sum"][1, 0, 0] == 1  # x fastest
    assert s.bin_width == (0.5, 0.5, 1.0) and s.bin_unit == ("mm", "mm", "mm")
    nohash = "\n".join(ln[2:] for ln in HEADER.splitlines())
    assert parse_topas_binary(nohash, arr.tobytes()).bins == (2, 3, 4)


def test_topas_binary_fail_closed() -> None:
    arr = np.arange(24, dtype="<f8")
    with pytest.raises(ParseError):
        parse_topas_binary(HEADER, arr.tobytes()[:-1])
    with pytest.raises(ParseError):
        parse_topas_binary(HEADER, arr.tobytes() + b"\0" * 8)
    with pytest.raises(ParseError):
        parse_topas_binary("# nothing\n", arr.tobytes())
    with pytest.raises(ParseError):
        parse_topas_binary(HEADER.replace("Sum", "Sum Standard_Deviation"), arr.tobytes())
    bad = arr.copy()
    bad[3] = np.nan
    with pytest.raises(ParseError):
        parse_topas_binary(HEADER, bad.tobytes())
    with pytest.raises(ValueError):
        parse_topas_binary(HEADER, arr.tobytes(), order="X")
