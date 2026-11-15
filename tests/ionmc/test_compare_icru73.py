"""Tests of validation/scripts/stopping/compare_icru73.py (row I3) on a synthetic archive.

The synthetic libdEdx archive has the ICRU73_NEW block format (``#target:ion:n`` headers, one
value per line, blank lines between blocks), a 53-point grid and values generated from the ionmc
Bethe model (78 eV) times known factors, so every gate outcome is known in advance. No real
libdEdx or ICRU data are used.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import sys
import tarfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc.data.cache import IntegrityError
from ionmc.materials import water
from ionmc.physics.projectiles import CARBON12, OXYGEN16
from ionmc.physics.stopping import bethe_mass_stopping

SCRIPT = Path(__file__).resolve().parents[2] / "validation/scripts/stopping/compare_icru73.py"
GRID = np.geomspace(0.025, 1000.0, 53)


@pytest.fixture(scope="module")
def mod():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location("compare_icru73", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["compare_icru73"] = module
    spec.loader.exec_module(module)
    return module


def _bethe78(z: int) -> np.ndarray:
    """Synthetic table: Bethe 78 eV at E >= 1 MeV/u, constant below (unused below 2 MeV/u)."""
    proj = {6: CARBON12, 8: OXYGEN16}[z]
    e = np.maximum(GRID, 1.0)
    return bethe_mass_stopping(e, proj, water(78.0))


Factor = Callable[[np.ndarray], np.ndarray]


def _build(tmp_path: Path, factor: Factor | None = None) -> tuple[Path, str, dict[str, str]]:
    """Write a synthetic archive; ``factor(E)`` scales the Z = 6 and 8 water tables."""
    energy_text = f"{GRID.size}\n" + "".join(f"{e:.10g}\n" for e in GRID)
    blocks = []
    for target in (104, 276):
        for z in range(1, 10):
            if target == 276 and z in (6, 8):
                s = _bethe78(z) * (factor(GRID) if factor else 1.0)
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


def _run(mod, tmp_path: Path, factor: Factor | None = None) -> dict[str, Any]:  # type: ignore[no-untyped-def]
    path, sha, members = _build(tmp_path, factor)
    return mod.run(path, "deadbeef", archive_sha256=sha, member_sha256=members)


def test_parser_selects_water_carbon_and_oxygen(mod, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    path, sha, members = _build(tmp_path)
    tables = mod.load_icru73_water(path, archive_sha256=sha, member_sha256=members)
    assert sorted(tables) == [6, 8]
    for z, t in tables.items():
        assert t.z == z and t.energy_per_u_mev.size == 53
        np.testing.assert_allclose(t.energy_per_u_mev, GRID, rtol=1e-9)
        np.testing.assert_allclose(t.stopping_mev_cm2_g, _bethe78(z), rtol=1e-11)


def test_real_member_hashes_pinned(mod) -> None:  # type: ignore[no-untyped-def]
    assert mod.ARCHIVE_SHA256.startswith("b941a0c3")
    assert mod.MEMBER_SHA256["data/raw/ICRU73_NEW.dat"].startswith("749768c1")
    assert mod.MEMBER_SHA256["data/raw/icru73_newEng.dat"].startswith("8b5b809f")


def test_wrong_hashes_fail_closed(mod, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    path, sha, members = _build(tmp_path)
    with pytest.raises(IntegrityError):
        mod.load_icru73_water(path, archive_sha256="0" * 64, member_sha256=members)
    bad = dict(members)
    bad["data/raw/ICRU73_NEW.dat"] = "0" * 64
    with pytest.raises(IntegrityError):
        mod.load_icru73_water(path, archive_sha256=sha, member_sha256=bad)
    with pytest.raises(IntegrityError):  # default pins are the real archive's
        mod.load_icru73_water(path)


def test_parser_rejects_malformed_blocks(mod) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(ValueError, match="declared 3"):
        mod.parse_stopping_file("#276:06:3\n1.0\n2.0\n")
    with pytest.raises(ValueError, match="duplicate"):
        mod.parse_stopping_file("#276:06:1\n1.0\n#276:6:1\n1.0\n")
    with pytest.raises(ValueError, match="before the first"):
        mod.parse_stopping_file("1.0\n#276:06:1\n1.0\n")
    with pytest.raises(ValueError, match="declares"):
        mod.parse_energy_file("3\n1\n2\n")


def test_unscaled_table_passes_and_garbage_below_2_mev_is_ignored(mod, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    out = _run(mod, tmp_path, lambda e: np.where(e < 2.0, 50.0, 1.0))
    assert out["pass"] is True
    for name, ion in out["ions"].items():
        assert all(ion["gates"].values()), name
        assert ion["stopping"]["ge_10_MeV_u"]["max_abs_relative_deviation"] < 1e-9
        assert ion["csda_range"]["max_abs_relative_deviation"] < 2e-3


def test_windows_respected(mod, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    out = _run(mod, tmp_path)
    for ion in out["ions"].values():
        low, high = ion["stopping"]["2_to_10_MeV_u"], ion["stopping"]["ge_10_MeV_u"]
        assert low["n_points"] == int(np.count_nonzero((GRID >= 2.0) & (GRID < 10.0)))
        assert high["n_points"] == int(np.count_nonzero(GRID >= 10.0))
        assert low["energy_per_u_range_mev"][0] >= 2.0 and low["energy_per_u_range_mev"][1] < 10
        assert high["energy_per_u_range_mev"][0] >= 10.0
        assert high["energy_per_u_range_mev"][1] == pytest.approx(1000.0)
        info = ion["stopping"]["informative_10_to_500_MeV_u"]
        assert info["energy_per_u_range_mev"][1] <= 500.0
    assert out["ions"]["carbon12"]["csda_range"]["energies_per_u_mev"] == [80.0, 290.0, 430.0]
    assert out["ions"]["oxygen16"]["csda_range"]["energies_per_u_mev"] == [100.0, 350.0, 430.0]


def test_scaled_table_fails_the_right_gates(mod, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    # S_table = 1.01 S_model: |dev| = 0.99 % everywhere -> max gate passes, RMS gate (0.8 %) and
    # range gate (0.5 %) fail, the 2-10 MeV/u gate (6 %) passes.
    out = _run(mod, tmp_path, lambda e: np.full_like(e, 1.01))
    assert out["pass"] is False
    for ion in out["ions"].values():
        g = ion["gates"]
        assert g["stopping_ge_10_MeV_u_max_abs_le_0.015"] is True
        assert g["stopping_ge_10_MeV_u_rms_le_0.008"] is False
        assert g["stopping_2_to_10_MeV_u_max_abs_le_0.06"] is True
        assert g["csda_range_max_abs_le_0.005"] is False
        assert 0.009 < ion["csda_range"]["max_abs_relative_deviation"] < 0.0115


def test_low_band_perturbation_fails_only_low_gate(mod, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    out = _run(mod, tmp_path, lambda e: np.where(e < 10.0, 1.08, 1.0))
    for ion in out["ions"].values():
        g = ion["gates"]
        assert g["stopping_2_to_10_MeV_u_max_abs_le_0.06"] is False
        assert g["stopping_ge_10_MeV_u_max_abs_le_0.015"] is True
        assert g["stopping_ge_10_MeV_u_rms_le_0.008"] is True
        assert g["csda_range_max_abs_le_0.005"] is True
        assert ion["pass"] is False


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


def test_no_per_energy_values_in_output(mod, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    out = _run(mod, tmp_path, lambda e: 1.0 + 0.01 * np.sin(np.log(e)))
    text = json.dumps(out)
    assert "redistribution_note" in out and "TO-VERIFY (P5b)" in out["i_value_note"]
    assert "75.0 eV" in out["i_value_note"] and "78 eV" in out["i_value_note"]
    floats = np.array(_floats(out))
    maxima = np.array(
        [
            band["max_abs_relative_deviation"]
            for ion in out["ions"].values()
            for band in ion["stopping"].values()
        ]
    )
    sel = GRID >= 2.0
    for z in (6, 8):
        table = _bethe78(z) * (1.0 + 0.01 * np.sin(np.log(GRID)))
        dev = _bethe78(z)[sel] / table[sel] - 1.0
        for value in np.concatenate((table, dev, np.abs(dev))):
            if np.any(np.isclose(floats, value, rtol=1e-9, atol=0.0)):
                # only the reported band maxima may coincide with a per-energy |deviation|
                assert np.any(np.isclose(maxima, value, rtol=1e-9, atol=0.0))
    assert max(len(v) for v in _lists(out)) <= 3
    assert "deadbeef" in text
    assert out["provenance"]["model"]["options"]["ls"] == "off"
    assert out["provenance"]["model"]["I_eV"] == 78.0


def test_cli_writes_json(mod, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:  # type: ignore[no-untyped-def]
    path, sha, members = _build(tmp_path)
    # the pinned defaults are the real archive's; route the CLI through the synthetic hashes
    monkeypatch.setattr(
        mod,
        "run",
        lambda a, c, _run=mod.run: _run(a, c, archive_sha256=sha, member_sha256=members),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "compare_icru73",
            "--archive",
            str(path),
            "--output-dir",
            str(tmp_path / "o"),
            "--code-sha",
            "abc",
        ],
    )
    assert mod.main() == 0
    data = json.loads((tmp_path / "o" / "icru73_comparison.json").read_text())
    assert data["provenance"]["code_sha"] == "abc" and data["row"] == "I3"
