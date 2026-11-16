"""Tests of validation/scripts/stopping/i3_band_diagnostic.py (row I3 diagnostic) on synthetic data.

The synthetic libdEdx archive has the ICRU73_NEW block format, the 53-point grid and liquid-water
blocks Z = 1-18 (plus air decoys); the Z = 1-10 water values are the ionmc Bethe model (78 eV,
default options, the script's own projectile definitions) times known factors, so every band
mean and z-scan outcome is known in advance. No real libdEdx or ICRU data are used.

The script rounds every table-relative band mean to 0.1 percentage point (``round(d, 3)``)
before writing it and before any derived quantity; 0.1 pp is the precision already present in Git
(plan results note of 2026-11-15). The plants below therefore use band deviations that are
multiples of 0.001 where an exact z-scan fit is asserted, and the derived quantities are checked
to be reproducible from the rounded means.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import re
import sys
import tarfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc.data.cache import IntegrityError
from ionmc.materials import water
from ionmc.physics.stopping import bethe_mass_stopping

SCRIPT = Path(__file__).resolve().parents[2] / "validation/scripts/stopping/i3_band_diagnostic.py"
# rounded to the 10 significant digits written to the energy file, so the parsed grid is exact
GRID = np.array([float(f"{e:.10g}") for e in np.geomspace(0.025, 1000.0, 53)])
SEL = GRID >= 2.0
I3_MAXIMA = ("max_abs_relative_deviation_2_to_10", "max_abs_relative_deviation_ge_10")

Factor = Callable[[int, np.ndarray], np.ndarray]


@pytest.fixture(scope="module")
def mod():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location("i3_band_diagnostic", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["i3_band_diagnostic"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def model(mod) -> dict[int, np.ndarray]:  # type: ignore[no-untyped-def]
    """Default-option Bethe 78 eV on the grid (constant below 1 MeV/u, unused below 2 MeV/u)."""
    e = np.maximum(GRID, 1.0)
    return {z: bethe_mass_stopping(e, p, water(78.0)) for z, p in mod.PROJECTILES.items()}


def _build(
    tmp_path: Path,
    model: dict[int, np.ndarray],
    factor: Factor | None = None,
    skip: tuple[int, ...] = (),
) -> tuple[Path, str, dict[str, str]]:
    """Write a synthetic archive; ``factor(z, E)`` scales the Z = 1-10 water tables."""
    energy_text = f"{GRID.size}\n" + "".join(f"{e:.10g}\n" for e in GRID)
    blocks = []
    for target in (104, 276):
        for z in range(1, 19):
            if target == 276 and z in skip:
                continue
            if target == 276 and z in model:
                s = model[z] * (factor(z, GRID) if factor else 1.0)
            else:
                s = np.full(GRID.size, 1.0e6 + z)  # decoys: must never be selected
            blocks.append(f"#{target}:{z:02d}:{GRID.size}\n" + "".join(f"{v:.12E}\n" for v in s))
    files = {
        "libdedx-main/data/raw/icru73_newEng.dat": energy_text.encode(),
        "libdedx-main/data/raw/ICRU73_NEW.dat": "\n\n".join(blocks).encode(),
        "libdedx-main/README.md": b"synthetic\n",
    }
    path = tmp_path / "archive.tar.gz"
    with tarfile.open(path, "w:gz") as tar:
        for name, data in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    members = {
        m: hashlib.sha256(files[f"libdedx-main/{m}"]).hexdigest()
        for m in ("data/raw/ICRU73_NEW.dat", "data/raw/icru73_newEng.dat")
    }
    return path, hashlib.sha256(path.read_bytes()).hexdigest(), members


def _run(mod, tmp_path: Path, model, factor: Factor | None = None) -> dict[str, Any]:  # type: ignore[no-untyped-def]
    path, sha, members = _build(tmp_path, model, factor)
    out: dict[str, Any] = mod.run(path, "deadbeef", archive_sha256=sha, member_sha256=members)
    return out


def _band_plant(mod, values: Callable[[int, int], float]) -> Factor:  # type: ignore[no-untyped-def]
    """Factor that plants the band-constant deviation ``values(z, band index)``."""

    def factor(z: int, e: np.ndarray) -> np.ndarray:
        b = _band_index(mod, e)
        dev = np.array([0.0 if i < 0 else values(z, int(i)) for i in b])
        return 1.0 / (1.0 + dev)

    return factor


def _band_index(mod, e: np.ndarray) -> np.ndarray:  # type: ignore[no-untyped-def]
    """Index of the band of each energy (-1 below 2 MeV/u)."""
    idx = np.full(e.shape, -1)
    for i, (lo, hi) in enumerate(mod.BANDS):
        idx[(e >= lo) & ((e < hi) | (e <= hi) & (i == len(mod.BANDS) - 1))] = i
    return idx


def test_band_means_recovered(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    def g(z: int, b: int) -> float:  # 0.7 of the 0.1 pp quantum off a multiple: no rounding tie
        return 0.001 * z * (b + 1) - 0.0023

    def factor(z: int, e: np.ndarray) -> np.ndarray:
        b = _band_index(mod, e)
        return np.where(b < 0, 1.0, 1.0 / (1.0 + np.array([g(z, i) for i in b])))

    out = _run(mod, tmp_path, model, factor)
    assert out["row"] == "I3-diagnostic" and "not re-graded" in out["status_note"]
    assert len(out["ions"]) == 10
    for ion in out["ions"].values():
        z = ion["z"]
        assert set(ion) - set(I3_MAXIMA) == {"z", "a", "mass_mev", "bands"}
        assert all((k in ion) == (z in (6, 8)) for k in I3_MAXIMA)
        for i, (key, band) in enumerate(ion["bands"].items()):
            assert set(band) == {"n_points", "band_per_u_mev", "mean_relative_deviation"}
            assert band["band_per_u_mev"] == list(mod.BANDS[i])
            assert band["n_points"] == int(np.count_nonzero(_band_index(mod, GRID) == i)), key
            # written rounded to 0.1 pp, not the planted full-precision value
            assert band["mean_relative_deviation"] == round(g(z, i), 3)
            assert band["mean_relative_deviation"] != pytest.approx(g(z, i), rel=1e-6)


# band coefficients that are multiples of 0.001: c z^2 and c are exact after rounding to 0.1 pp
BLOCH_COEFF = (0.002, -0.001, 0.001, 0.003, 0.001)
CONSTANT_COEFF = (0.002, -0.001, 0.003, 0.001, 0.001)
# charge plant: per-energy K (1 - q^2/z^2); band means K f_charge are not multiples of 0.001
CHARGE_K = (0.5, 1.0, 6.0, 10.0, 1.0)


def _zscan_plant(mod, model, tmp_path: Path, shape: str) -> dict[str, Any]:  # type: ignore[no-untyped-def]
    def factor(z: int, e: np.ndarray) -> np.ndarray:
        b = _band_index(mod, e)
        if shape == "bloch":
            dev = np.where(b < 0, 0.0, np.array([BLOCH_COEFF[i] for i in b])) * z**2
        elif shape == "constant":
            dev = np.where(b < 0, 0.0, np.array([CONSTANT_COEFF[i] for i in b]))
        else:
            k = np.where(b < 0, 0.0, np.array([CHARGE_K[i] for i in b]))
            dev = k * mod.charge_shape(e, mod.PROJECTILES[z])
        return 1.0 / (1.0 + dev)

    return _run(mod, tmp_path, model, factor)


def test_z_scan_identifies_bloch_like(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    scan = _zscan_plant(mod, model, tmp_path, "bloch")["z_scan"]["bands"]
    assert len(scan) == 4
    for i, band in enumerate(scan.values()):
        assert band["ranking"][0] == "bloch_like"
        fit = band["hypotheses"]["bloch_like"]
        assert fit["residual_rms"] < 1e-12
        assert fit["coefficient"] == pytest.approx(BLOCH_COEFF[i], rel=1e-10)
        assert fit["r_squared"] == pytest.approx(1.0, abs=1e-12)
        assert band["shape_ratio_to_z6"]["10"] == pytest.approx(100.0 / 36.0, rel=1e-10)
        assert set(band["f_charge_band_mean"]) == {str(z) for z in range(3, 11)}
        # rule 1b: clean plant, attributed
        verdict = band["verdict"]
        assert verdict["attributed_to"] == "bloch_like" and verdict["reason"] is None
        assert verdict["common_sign"] is True and verdict["rms_ratio_1_to_2"] < 1e-9


def test_z_scan_identifies_constant(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    out = _zscan_plant(mod, model, tmp_path, "constant")
    for i, (key, band) in enumerate(out["z_scan"]["bands"].items()):
        assert band["ranking"][0] == "constant" and len(band["ranking"]) == 4
        fit = band["hypotheses"]["constant"]
        assert fit["residual_rms"] < 1e-12
        assert fit["coefficient"] == pytest.approx(CONSTANT_COEFF[i], rel=1e-10)
        assert band["verdict"]["attributed_to"] == "constant"
        z1 = out["informative_z1_z2"]["bands"][key]["z1"]
        assert z1["prediction_from_z3_10_fit"]["constant"] == pytest.approx(fit["coefficient"])
        assert z1["mean_relative_deviation"] == CONSTANT_COEFF[i]


def test_z_scan_identifies_charge(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    out = _zscan_plant(mod, model, tmp_path, "charge")
    for i, band in enumerate(out["z_scan"]["bands"].values()):
        assert band["ranking"][0] == "charge"
        fit = band["hypotheses"]["charge"]
        # the rounding to 0.1 pp leaves a residual (<= 0.05 pp per ion) and biases c slightly
        assert fit["residual_rms"] < 5e-4
        assert fit["coefficient"] == pytest.approx(CHARGE_K[i], rel=0.02)
        verdict = band["verdict"]
        if i < 3:
            assert verdict["attributed_to"] == "charge" and verdict["reason"] is None
        else:
            # [20, 50): the planted d_z of Z = 3, 4 round to 0 (the quantisation limit of this
            # band), so the common-sign condition fails whatever r^2 says
            assert verdict["attributed_to"] is None and "sign_change" in verdict["reason"]
            assert verdict["common_sign"] is False
    info = out["informative_z1_z2"]
    assert "table lineage for Z = 1, 2 OPEN" in info["lineage_note"]
    for band in info["bands"].values():
        assert band["z1"]["mean_relative_deviation"] == 0.0
        assert band["z1"]["f_charge_band_mean"] == 0.0


def _z_plant(mod, model, tmp_path: Path, d: dict[int, float]) -> dict[str, Any]:  # type: ignore[no-untyped-def]
    """The same band-constant d_z (multiples of 0.001) in every band below 50 MeV/u."""
    return _run(
        mod, tmp_path, model, _band_plant(mod, lambda z, b: d.get(z, 0.0) if b < 4 else 0.0)
    )


def test_z_scan_not_attributed_on_sign_change(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    # bloch-like for Z = 4-10, but Z = 3 negative: r^2 >= 0.99, yet the sign rule refuses
    d = {z: 0.001 * z**2 for z in range(4, 11)} | {3: -0.001}
    scan = _z_plant(mod, model, tmp_path, d)["z_scan"]["bands"]
    for key, band in scan.items():
        verdict = band["verdict"]
        best = band["hypotheses"][band["ranking"][0]]
        assert best["r_squared"] >= 0.99
        assert verdict["attributed_to"] is None and verdict["common_sign"] is False
        assert "sign_change" in verdict["reason"] and "r_squared" not in verdict["reason"]
        if key in ("10_to_20_MeV_u", "20_to_50_MeV_u"):
            # there the RMS ratio would allow attribution: the sign rule alone decides
            assert verdict["reason"] == ["sign_change"] and verdict["rms_ratio_1_to_2"] <= 0.5


def test_z_scan_not_attributed_for_collinear_shapes(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    # d_z = 0.0005 z (z + 5): half z, half z^2 (multiples of 0.001). Every shape is monotone in
    # z, so the two best shapes fit about equally well (RMS ratio > 0.5) while r^2 > 0.99 and
    # the sign is common: only the RMS-ratio condition fails
    d = {z: 0.0005 * z * (z + 5) for z in range(3, 11)}
    scan = _z_plant(mod, model, tmp_path, d)["z_scan"]["bands"]
    for band in scan.values():
        verdict = band["verdict"]
        fits = band["hypotheses"]
        rms_1, rms_2 = (fits[n]["residual_rms"] for n in band["ranking"][:2])
        assert 0.5 < rms_1 / rms_2 < 1.0
        assert fits[band["ranking"][0]]["r_squared"] >= 0.99
        assert verdict["rms_ratio_1_to_2"] == pytest.approx(rms_1 / rms_2, rel=1e-12)
        assert verdict["attributed_to"] is None and verdict["reason"] == ["rms_ratio"]
        assert verdict["common_sign"] is True


def _hypothesis_plant(mod, overrides: dict[tuple[str, int], float]) -> Factor:  # type: ignore[no-untyped-def]
    """C/O band deviations at the quoted hypothesis values (0 at >= 50 MeV/u), plus overrides."""
    names = {6: "carbon12", 8: "oxygen16"}

    def values(z: int, b: int) -> float:
        if z not in names:
            return 0.0
        planted: list[float] = [*mod.HYPOTHESIS[names[z]].values(), 0.0]
        return overrides.get((names[z], b), planted[b])

    return _band_plant(mod, values)


@pytest.mark.parametrize(
    ("overrides", "verdict", "failed"),
    [
        ({}, "CONFIRMED", None),
        # within the 0.1 pp quantum: |delta| = 0.04 pp; |mean >= 50| = 0.04 % rounds to 0
        ({("carbon12", 0): -0.0864, ("oxygen16", 4): 0.0004}, "CONFIRMED", None),
        ({("carbon12", 2): -0.017}, "WITHDRAWN", ("carbon12", "10_to_20_MeV_u")),
        ({("oxygen16", 4): 0.001}, "WITHDRAWN", ("oxygen16", "50_to_1000_MeV_u")),
    ],
)
def test_hypothesis_check(mod, model, tmp_path: Path, overrides, verdict, failed) -> None:  # type: ignore[no-untyped-def]
    out = _run(mod, tmp_path, model, _hypothesis_plant(mod, overrides))
    check = out["hypothesis_check"]
    assert check["verdict"] == verdict and "rule 1a" in check["rule"]
    assert set(check["per_ion"]) == {"carbon12", "oxygen16"}
    for name, entry in check["per_ion"].items():
        assert list(entry) == [*(f"{k}_MeV_u" for k in mod.HYPOTHESIS[name]), "50_to_1000_MeV_u"]
        for band, item in entry.items():
            assert item["committed"] == out["ions"][name]["bands"][band]["mean_relative_deviation"]
            assert item["match"] is (failed != (name, band))
            if band == "50_to_1000_MeV_u":
                assert item["quoted_abs_bound"] == mod.HYPOTHESIS_HIGH_ABS[name]
                assert set(item) == {"quoted_abs_bound", "committed", "match"}
            else:
                assert item["quoted"] == mod.HYPOTHESIS[name][band.removesuffix("_MeV_u")]
                assert set(item) == {"quoted", "committed", "match"}


def test_missing_block_and_wrong_hashes_fail_closed(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    path, sha, members = _build(tmp_path, model, skip=(5,))
    with pytest.raises(ValueError, match="ion Z = 5"):
        mod.run(path, None, archive_sha256=sha, member_sha256=members)
    path, sha, members = _build(tmp_path, model)
    with pytest.raises(IntegrityError):
        mod.run(path, None, archive_sha256="0" * 64, member_sha256=members)
    bad = dict(members)
    bad["data/raw/icru73_newEng.dat"] = "0" * 64
    with pytest.raises(IntegrityError):
        mod.run(path, None, archive_sha256=sha, member_sha256=bad)
    with pytest.raises(IntegrityError):  # default pins are the real archive's
        mod.run(path, None)
    assert mod.ARCHIVE_SHA256.startswith("b941a0c3")


def test_band_with_fewer_than_two_points_fails_closed(mod) -> None:  # type: ignore[no-untyped-def]
    e = np.array([3.0, 4.0, 6.0, 12.0, 15.0, 30.0, 40.0, 60.0, 70.0])
    with pytest.raises(ValueError, match="per-energy"):
        mod.band_means(e, np.ones_like(e))


def _floats(obj: Any) -> list[float]:
    if isinstance(obj, dict):
        return [x for v in obj.values() for x in _floats(v)]
    if isinstance(obj, list):
        return [x for v in obj for x in _floats(v)]
    return [float(obj)] if isinstance(obj, float) else []


def _lists(obj: Any) -> list[list[Any]]:
    if isinstance(obj, dict):
        return [x for v in obj.values() for x in _lists(v)]
    if isinstance(obj, list):
        return [obj, *[x for v in obj for x in _lists(v)]]
    return []


def test_no_per_energy_values_in_output(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    def factor(z: int, e: np.ndarray) -> np.ndarray:
        return 1.0 + 0.01 * np.sin(np.log(e))

    out = _run(mod, tmp_path, model, factor)
    text = json.dumps(out)
    assert out["redistribution_note"] and "TO-VERIFY (P5b)" in out["i_value_note"]
    # the carbon/oxygen band maxima (already published by row I3) are the only floats allowed to
    # coincide with a per-energy |deviation|; with a Z-independent factor they also coincide with
    # the |deviation| of other ions at the same energy, so they are removed before the check
    maxima = {
        z: [out["ions"][n].pop(k) for k in I3_MAXIMA] for z, n in ((6, "carbon12"), (8, "oxygen16"))
    }
    floats = np.array(_floats(out))
    for z in range(1, 11):
        table = model[z] * factor(z, GRID)
        dev = model[z][SEL] / table[SEL] - 1.0
        if z in (6, 8):
            low = GRID[SEL] < 10.0
            assert np.max(np.abs(dev[low])) == pytest.approx(maxima[z][0], rel=1e-9)
            assert np.max(np.abs(dev[~low])) == pytest.approx(maxima[z][1], rel=1e-9)
        for value in np.concatenate((table, dev, np.abs(dev))):
            assert not np.any(np.isclose(floats, value, rtol=1e-9, atol=0.0)), (z, value)
    # no numeric list longer than 3; the only longer lists are the 4-name hypothesis rankings
    for lst in _lists(out):
        assert len(lst) <= 3 or (len(lst) == 4 and all(isinstance(v, str) for v in lst)), lst
    assert "deadbeef" in text
    prov = out["provenance"]
    assert prov["model"]["options"]["ls"] == "off" and prov["model"]["I_eV"] == 78.0
    assert prov["projectiles"]["nitrogen14"]["mass"] == "A u - Z m_e"
    assert 0.0 < prov["mass_sensitivity_max_rel"] < 1e-3
    assert "A u - Z m_e" in prov["mass_note"]


def test_i3b_variants(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    out = _run(mod, tmp_path, model)
    i3b = out["i3b"]["ions"]
    assert set(i3b) == {"carbon12", "oxygen16"}
    assert i3b["carbon12"]["standard_energies_per_u_mev"] == [80.0, 290.0, 430.0]
    for ion in i3b.values():
        variants = ion["variants"]
        assert list(variants) == list(mod.VARIANTS)
        for band in variants["default"]["bands"].values():
            assert abs(band["mean_relative_deviation"]) < 1e-11
            assert band["model_only_stopping_ratio_to_default_band_mean"] == 1.0
        assert set(variants["default"]["csda_range_shift_mm_vs_default"].values()) == {0.0}
        for name, variant in variants.items():
            for band in variant["bands"].values():
                # only the default variant carries a table-relative value
                assert ("mean_relative_deviation" in band) == (name == "default"), name
                assert set(band) - {"mean_relative_deviation"} == {
                    "n_points",
                    "band_per_u_mev",
                    "model_only_stopping_ratio_to_default_band_mean",
                }
        low = variants["effective_charge_off"]["bands"]["2_to_5_MeV_u"]
        assert low["model_only_stopping_ratio_to_default_band_mean"] > 1.0  # bare charge > q
        shifts = variants["effective_charge_off"]["csda_range_shift_mm_vs_default"]
        assert len(shifts) == 3 and all(v < 0.0 for v in shifts.values())  # more stopping
        assert variants["ls_eff"]["options"]["ls"] == "eff"


TABLE_DEPENDENT_WHITELIST = (
    # first-hand table aggregates: one band mean per (ion, band); the C/O band maxima that the
    # row I3 file already publishes
    r"ions/[a-z0-9]+/bands/[^/]+/mean_relative_deviation",
    r"ions/(carbon12|oxygen16)/max_abs_relative_deviation_(2_to_10|ge_10)",
    # exact duplicates of the band means (asserted below)
    r"informative_z1_z2/bands/[^/]+/z[12]/mean_relative_deviation",
    r"i3b/ions/(carbon12|oxygen16)/variants/default/bands/[^/]+/mean_relative_deviation",
    r"hypothesis_check/per_ion/(carbon12|oxygen16)/[^/]+/committed",
    # functions of the rounded Z = 3-10 band means only (z-scan fits, ratios, fit predictions,
    # verdicts; reproduced from the rounded means below)
    r"z_scan/bands/[^/]+/hypotheses/[a-z_]+/(coefficient|residual_rms|r_squared)",
    r"z_scan/bands/[^/]+/shape_ratio_to_z6/[0-9]+",
    r"z_scan/bands/[^/]+/verdict/rms_ratio_1_to_2",
    r"informative_z1_z2/bands/[^/]+/z[12]/prediction_from_z3_10_fit/[a-z_]+",
)
MEAN_PATTERN = r".*/mean_relative_deviation|hypothesis_check/.*/committed"


def _paths(obj: Any, prefix: str = "") -> dict[str, float]:
    if isinstance(obj, dict):
        return {k: v for key, val in obj.items() for k, v in _paths(val, f"{prefix}{key}/").items()}
    if isinstance(obj, list):
        return {k: v for i, val in enumerate(obj) for k, v in _paths(val, f"{prefix}{i}/").items()}
    return {prefix.rstrip("/"): float(obj)} if isinstance(obj, float) else {}


def test_only_whitelisted_fields_depend_on_the_table(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    """Two different tables: every float that changes must be on the whitelist.

    With the duplicates asserted equal, each (ion, band) carries one table quantity, the band
    mean ROUNDED to 0.1 pp, plus the two full-precision band maxima for carbon and oxygen that
    row I3 already publishes. A full-precision mean would not be safe: on the real 53-point grid
    the [10, 20) band holds two points and its 10 MeV/u point carries the published >= 10 MeV/u
    maximum (monotone deficit), so d(15) = 2 mean - (+-max) and the table value at 15 MeV/u would
    follow from the public model. The rounded mean fixes d(15) only to +-0.1 pp, the precision
    already present in Git (plan results note of 2026-11-15). Every derived quantity (z-scan fits,
    ratios, predictions, verdicts, hypothesis check) is reproduced here from the rounded means, so
    the output carries no information beyond them.
    """
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    out_a = _run(mod, tmp_path / "a", model, lambda z, e: 1.0 + 0.01 * (1 + z) * np.sin(np.log(e)))
    out_b = _run(mod, tmp_path / "b", model, lambda z, e: 1.0 + 0.02 * np.cos(np.log(e)))
    a, b = _paths(out_a), _paths(out_b)
    assert set(a) == set(b)
    changed = {p for p in a if a[p] != b[p]}
    assert changed  # the check is not vacuous
    for path in changed:
        assert any(re.fullmatch(w, path) for w in TABLE_DEPENDENT_WHITELIST), path
    n_duplicates = 0
    for run in (a, b):
        for path, value in run.items():
            m = re.fullmatch(
                r"informative_z1_z2/bands/([^/]+)/z([12])/mean_relative_deviation", path
            )
            if m:
                name = {"1": "proton", "2": "alpha"}[m.group(2)]
                assert value == run[f"ions/{name}/bands/{m.group(1)}/mean_relative_deviation"]
                n_duplicates += 1
            m = re.fullmatch(
                r"i3b/ions/([a-z0-9]+)/variants/default/bands/([^/]+)/mean_relative_deviation",
                path,
            )
            if m:
                assert value == run[f"ions/{m.group(1)}/bands/{m.group(2)}/mean_relative_deviation"]
                n_duplicates += 1
    assert n_duplicates == 2 * (2 * 4 + 2 * 5)  # both runs: z1/z2 x 4 bands, C/O x 5 bands
    for run, out in ((a, out_a), (b, out_b)):
        means = [v for p, v in run.items() if re.fullmatch(MEAN_PATTERN, p)]
        assert len(means) == 10 * 5 + 2 * 4 + 2 * 5 + 2 * 5
        assert all(v == round(v, 3) for v in means)  # at most 3 decimals (0.1 pp)
        _assert_reproducible_from_rounded_means(mod, out)


def _assert_reproducible_from_rounded_means(mod, out: dict[str, Any]) -> None:  # type: ignore[no-untyped-def]
    """Recompute every z-scan quantity from the output's rounded means (1e-12)."""
    by_z = {ion["z"]: ion for ion in out["ions"].values()}
    zs = np.arange(3.0, 11.0)
    for key, band in out["z_scan"]["bands"].items():
        d = np.array([by_z[z]["bands"][key]["mean_relative_deviation"] for z in range(3, 11)])
        f_charge = np.array([band["f_charge_band_mean"][str(z)] for z in range(3, 11)])
        shapes = {
            "charge": f_charge,
            "barkas_like": zs,
            "bloch_like": zs**2,
            "constant": np.ones_like(zs),
        }
        coef = {}
        for name, f in shapes.items():
            c = float(f @ d / (f @ f))
            res = d - c * f
            fit = band["hypotheses"][name]
            assert fit["coefficient"] == pytest.approx(c, rel=0.0, abs=1e-12)
            assert fit["residual_rms"] == pytest.approx(np.sqrt(np.mean(res**2)), abs=1e-12)
            assert fit["r_squared"] == pytest.approx(1.0 - res @ res / (d @ d), abs=1e-12)
            coef[name] = c
        rms = sorted(band["hypotheses"][n]["residual_rms"] for n in shapes)
        assert [band["hypotheses"][n]["residual_rms"] for n in band["ranking"]] == rms
        for z in range(3, 11):
            assert band["shape_ratio_to_z6"][str(z)] == pytest.approx(d[z - 3] / d[3], abs=1e-12)
        r2_1 = band["hypotheses"][band["ranking"][0]]["r_squared"]
        common = bool(np.all(d > 0) or np.all(d < 0))
        ok = rms[0] <= 0.5 * rms[1] and r2_1 >= 0.9 and common
        assert band["verdict"]["rms_ratio_1_to_2"] == pytest.approx(rms[0] / rms[1], abs=1e-12)
        assert band["verdict"]["attributed_to"] == (band["ranking"][0] if ok else None)
        assert band["verdict"]["common_sign"] is common
        for z in (1, 2):
            info = out["informative_z1_z2"]["bands"][key][f"z{z}"]
            pred = info["prediction_from_z3_10_fit"]
            assert pred["charge"] == pytest.approx(
                coef["charge"] * info["f_charge_band_mean"], abs=1e-12
            )
            assert pred["barkas_like"] == pytest.approx(coef["barkas_like"] * z, abs=1e-12)
            assert pred["bloch_like"] == pytest.approx(coef["bloch_like"] * z**2, abs=1e-12)
            assert pred["constant"] == pytest.approx(coef["constant"], abs=1e-12)
            assert (
                info["mean_relative_deviation"] == by_z[z]["bands"][key]["mean_relative_deviation"]
            )


def _i3_file(mod, tmp_path: Path, path: Path, sha: str, members: dict[str, str]) -> Path:  # type: ignore[no-untyped-def]
    """Row I3 JSON of the same synthetic archive (compare_icru73.run, as the CLI writes it)."""
    i3 = mod.I3.run(path, "i3sha", archive_sha256=sha, member_sha256=members)
    out = tmp_path / "icru73_comparison.json"
    out.write_text(json.dumps(i3, indent=2) + "\n", encoding="utf-8")
    return out


def test_cli_expect_i3(mod, model, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:  # type: ignore[no-untyped-def]
    path, sha, members = _build(tmp_path, model, lambda z, e: 1.0 + 0.01 * np.sin(np.log(e)))
    i3_path = _i3_file(mod, tmp_path, path, sha, members)
    monkeypatch.setattr(
        mod,
        "run",
        lambda a, c, _run=mod.run: _run(a, c, archive_sha256=sha, member_sha256=members),
    )

    def cli(out_dir: Path, expect: Path | None) -> int:
        argv = ["i3_band_diagnostic", "--archive", str(path), "--output-dir", str(out_dir)]
        argv += ["--code-sha", "abc"]
        monkeypatch.setattr(
            sys, "argv", argv if expect is None else [*argv, "--expect-i3", str(expect)]
        )
        result: int = mod.main()
        return result

    # --expect-i3 is required: no run, no output without it
    with pytest.raises(SystemExit) as excinfo:
        cli(tmp_path / "none", None)
    assert excinfo.value.code == 2
    assert not (tmp_path / "none").exists()

    assert cli(tmp_path / "ok", i3_path) == 0
    raw = (tmp_path / "ok" / "icru73_band_diagnostic.json").read_text()
    assert raw.endswith("}\n") and not raw.endswith("\n\n")
    data = json.loads(raw)
    assert data["provenance"]["code_sha"] == "abc" and data["row"] == "I3-diagnostic"
    assert "0.1 percentage point" in data["rounding_note"]
    record = data["provenance"]["expect_i3"]
    assert record["matched"] is True and record["path"] == str(i3_path)
    assert record["sha256"] == hashlib.sha256(i3_path.read_bytes()).hexdigest()
    i3 = json.loads(i3_path.read_text())
    for name in ("carbon12", "oxygen16"):
        stopping = i3["ions"][name]["stopping"]
        ion = data["ions"][name]
        assert ion["max_abs_relative_deviation_2_to_10"] == pytest.approx(
            stopping["2_to_10_MeV_u"]["max_abs_relative_deviation"], rel=1e-9
        )
        assert ion["max_abs_relative_deviation_ge_10"] == pytest.approx(
            stopping["ge_10_MeV_u"]["max_abs_relative_deviation"], rel=1e-9
        )

    stopping = i3["ions"]["oxygen16"]["stopping"]["ge_10_MeV_u"]
    stopping["max_abs_relative_deviation"] *= 1.0 + 1e-6
    tampered = tmp_path / "tampered.json"
    tampered.write_text(json.dumps(i3), encoding="utf-8")
    with pytest.raises(RuntimeError, match="oxygen16"):
        cli(tmp_path / "bad", tampered)
    assert not (tmp_path / "bad" / "icru73_band_diagnostic.json").exists()
