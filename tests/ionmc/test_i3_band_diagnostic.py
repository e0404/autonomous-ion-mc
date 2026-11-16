"""Tests of validation/scripts/stopping/i3_band_diagnostic.py (row I3 diagnostic) on synthetic data.

The synthetic libdEdx archive has the ICRU73_NEW block format, the 53-point grid and liquid-water
blocks Z = 1-18 (plus air decoys); the Z = 1-10 water values are the ionmc Bethe model (78 eV,
default options, the script's own projectile definitions) times known factors, so every band
mean and z-scan outcome is known in advance. No real libdEdx or ICRU data are used.
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


def _band_index(mod, e: np.ndarray) -> np.ndarray:  # type: ignore[no-untyped-def]
    """Index of the band of each energy (-1 below 2 MeV/u)."""
    idx = np.full(e.shape, -1)
    for i, (lo, hi) in enumerate(mod.BANDS):
        idx[(e >= lo) & ((e < hi) | (e <= hi) & (i == len(mod.BANDS) - 1))] = i
    return idx


def test_band_means_recovered(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    def g(z: int, b: int) -> float:
        return 0.001 * z * (b + 1) - 0.0025

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
            assert band["mean_relative_deviation"] == pytest.approx(g(z, i), rel=1e-9)


def _zscan_plant(mod, model, tmp_path: Path, shape: str) -> dict[str, Any]:  # type: ignore[no-untyped-def]
    coeff = (2e-4, -1e-4, 5e-5, 3e-5, 1e-5)

    def factor(z: int, e: np.ndarray) -> np.ndarray:
        b = _band_index(mod, e)
        c = np.where(b < 0, 0.0, np.array([coeff[i] for i in b]))
        if shape == "bloch":
            dev = c * z**2
        elif shape == "constant":
            dev = c
        else:
            dev = 100.0 * c * mod.charge_shape(e, mod.PROJECTILES[z])
        return 1.0 / (1.0 + dev)

    return _run(mod, tmp_path, model, factor)


def test_z_scan_identifies_bloch_like(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    scan = _zscan_plant(mod, model, tmp_path, "bloch")["z_scan"]["bands"]
    assert len(scan) == 4
    for i, band in enumerate(scan.values()):
        assert band["ranking"][0] == "bloch_like"
        fit = band["hypotheses"]["bloch_like"]
        assert fit["residual_rms"] < 1e-11
        assert fit["coefficient"] == pytest.approx((2e-4, -1e-4, 5e-5, 3e-5)[i], rel=1e-8)
        assert fit["r_squared"] == pytest.approx(1.0, abs=1e-12)
        assert band["shape_ratio_to_z6"]["10"] == pytest.approx(100.0 / 36.0, rel=1e-8)
        assert set(band["f_charge_band_mean"]) == {str(z) for z in range(3, 11)}


def test_z_scan_identifies_constant(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    out = _zscan_plant(mod, model, tmp_path, "constant")
    for i, (key, band) in enumerate(out["z_scan"]["bands"].items()):
        assert band["ranking"][0] == "constant" and len(band["ranking"]) == 4
        fit = band["hypotheses"]["constant"]
        assert fit["residual_rms"] < 1e-11
        assert fit["coefficient"] == pytest.approx((2e-4, -1e-4, 5e-5, 3e-5)[i], rel=1e-8)
        z1 = out["informative_z1_z2"]["bands"][key]["z1"]
        assert z1["prediction_from_z3_10_fit"]["constant"] == pytest.approx(fit["coefficient"])
        assert z1["mean_relative_deviation"] == pytest.approx(fit["coefficient"], rel=1e-8)


def test_z_scan_identifies_charge(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    out = _zscan_plant(mod, model, tmp_path, "charge")
    for i, band in enumerate(out["z_scan"]["bands"].values()):
        assert band["ranking"][0] == "charge"
        fit = band["hypotheses"]["charge"]
        assert fit["residual_rms"] < 1e-11
        assert fit["coefficient"] == pytest.approx(100.0 * (2e-4, -1e-4, 5e-5, 3e-5)[i], rel=1e-8)
    info = out["informative_z1_z2"]
    assert "table lineage for Z = 1, 2 OPEN" in info["lineage_note"]
    for band in info["bands"].values():
        assert band["z1"]["mean_relative_deviation"] == pytest.approx(0.0, abs=1e-12)
        assert band["z1"]["f_charge_band_mean"] == 0.0


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
    # functions of the Z = 3-10 band means only (z-scan fits, ratios, fit predictions)
    r"z_scan/bands/[^/]+/hypotheses/[a-z_]+/(coefficient|residual_rms|r_squared)",
    r"z_scan/bands/[^/]+/shape_ratio_to_z6/[0-9]+",
    r"informative_z1_z2/bands/[^/]+/z[12]/prediction_from_z3_10_fit/[a-z_]+",
)


def _paths(obj: Any, prefix: str = "") -> dict[str, float]:
    if isinstance(obj, dict):
        return {k: v for key, val in obj.items() for k, v in _paths(val, f"{prefix}{key}/").items()}
    if isinstance(obj, list):
        return {k: v for i, val in enumerate(obj) for k, v in _paths(val, f"{prefix}{i}/").items()}
    return {prefix.rstrip("/"): float(obj)} if isinstance(obj, float) else {}


def test_only_whitelisted_fields_depend_on_the_table(mod, model, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    """Two different tables: every float that changes must be on the whitelist.

    With the duplicates asserted equal, each (ion, band) carries one table equation (the band
    mean), plus the overall maximum for carbon and oxygen; nothing allows a linear solve for the
    per-energy table values of a band.
    """
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    a = _paths(
        _run(mod, tmp_path / "a", model, lambda z, e: 1.0 + 0.01 * (1 + z) * np.sin(np.log(e)))
    )
    b = _paths(_run(mod, tmp_path / "b", model, lambda z, e: 1.0 + 0.02 * np.cos(np.log(e))))
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


def test_cli_writes_json(mod, model, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:  # type: ignore[no-untyped-def]
    path, sha, members = _build(tmp_path, model)
    monkeypatch.setattr(
        mod,
        "run",
        lambda a, c, _run=mod.run: _run(a, c, archive_sha256=sha, member_sha256=members),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "i3_band_diagnostic",
            "--archive",
            str(path),
            "--output-dir",
            str(tmp_path / "o"),
            "--code-sha",
            "abc",
        ],
    )
    assert mod.main() == 0
    raw = (tmp_path / "o" / "icru73_band_diagnostic.json").read_text()
    assert raw.endswith("}\n") and not raw.endswith("\n\n")
    data = json.loads(raw)
    assert data["provenance"]["code_sha"] == "abc" and data["row"] == "I3-diagnostic"
    assert data["provenance"]["expect_i3"] is None


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

    def cli(out_dir: Path, expect: Path) -> int:
        argv = ["i3_band_diagnostic", "--archive", str(path), "--output-dir", str(out_dir)]
        monkeypatch.setattr(sys, "argv", [*argv, "--expect-i3", str(expect)])
        result: int = mod.main()
        return result

    assert cli(tmp_path / "ok", i3_path) == 0
    data = json.loads((tmp_path / "ok" / "icru73_band_diagnostic.json").read_text())
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
