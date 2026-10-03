"""Tests of the reference-output parsers on small synthetic files."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ionmc.reference import ParseError, read_metaimage, read_topas_csv

TOPAS_CSV = """# TOPAS Version: 4.3
# Parameter File: input.txt
# Results for scorer IDD
# Scored in component: Water
# X in 1 bin  of 120 mm
# Y in 1 bin  of 120 mm
# Z in 3 bins of 1 mm
# DoseToMedium ( Gy ) : Sum  Standard_Deviation
0, 0, 0, 1.5, 0.1
0, 0, 2, 3.5, 0.3
0, 0, 1, 2.5, 0.2
"""


def test_topas_csv(tmp_path: Path) -> None:
    p = tmp_path / "idd.csv"
    p.write_text(TOPAS_CSV)
    s = read_topas_csv(p)
    assert s.bins == (1, 1, 3)
    assert s.bin_width == (120.0, 120.0, 1.0)
    assert s.bin_unit == ("mm", "mm", "mm")
    assert (s.quantity, s.unit) == ("DoseToMedium", "Gy")
    assert s.statistics == ["Sum", "Standard_Deviation"]
    assert s.meta["scorer"] == "IDD" and s.meta["component"] == "Water"
    np.testing.assert_allclose(s.values["Sum"][0, 0], [1.5, 2.5, 3.5])
    np.testing.assert_allclose(s.values["Standard_Deviation"][0, 0], [0.1, 0.2, 0.3])


def test_topas_csv_missing_bins_are_nan(tmp_path: Path) -> None:
    p = tmp_path / "a.csv"
    p.write_text(TOPAS_CSV.replace("0, 0, 1, 2.5, 0.2\n", ""))
    s = read_topas_csv(p)
    assert np.isnan(s.values["Sum"][0, 0, 1])


def test_topas_csv_errors(tmp_path: Path) -> None:
    p = tmp_path / "bad.csv"
    p.write_text(TOPAS_CSV.replace("0, 0, 2,", "0, 0, 7,"))
    with pytest.raises(ParseError):
        read_topas_csv(p)
    p.write_text("# Z in 3 bins of 1 mm\n0,0,0,1\n")
    with pytest.raises(ParseError):
        read_topas_csv(p)
    p.write_text(TOPAS_CSV.replace("0, 0, 0, 1.5, 0.1", "0, 0, 0, 1.5"))
    with pytest.raises(ParseError):
        read_topas_csv(p)


def _write_mhd(tmp_path: Path, arr: np.ndarray, msb: bool = False) -> Path:
    nz, ny, nx = arr.shape
    dt = ">f4" if msb else "<f4"
    arr.astype(dt).tofile(tmp_path / "v.raw")
    mhd = tmp_path / "v.mhd"
    mhd.write_text(
        "ObjectType = Image\nNDims = 3\n"
        f"DimSize = {nx} {ny} {nz}\nElementSpacing = 2 2 1\nOffset = 1 2 3\n"
        f"ElementType = MET_FLOAT\nElementByteOrderMSB = {msb}\nElementDataFile = v.raw\n"
    )
    return mhd


@pytest.mark.parametrize("msb", [False, True])
def test_metaimage_roundtrip(tmp_path: Path, msb: bool) -> None:
    arr = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)  # (nz, ny, nx)
    img = read_metaimage(_write_mhd(tmp_path, arr, msb))
    assert img.dims == (4, 3, 2)
    assert img.spacing == (2.0, 2.0, 1.0)
    assert img.offset == (1.0, 2.0, 3.0)
    np.testing.assert_array_equal(img.data, arr)
    assert img.data[1, 2, 3] == 23  # x fastest: index = x + nx*(y + ny*z)


def test_metaimage_errors(tmp_path: Path) -> None:
    mhd = _write_mhd(tmp_path, np.zeros((2, 3, 4), dtype=np.float32))
    (tmp_path / "v.raw").write_bytes(b"\x00" * 8)
    with pytest.raises(ParseError):
        read_metaimage(mhd)
    mhd.write_text("NDims = 3\nDimSize = 1 1 1\n")
    with pytest.raises(ParseError):
        read_metaimage(mhd)


def test_metaimage_local(tmp_path: Path) -> None:
    arr = np.arange(2 * 3 * 4, dtype="<f4").reshape(2, 3, 4)
    head = (
        "ObjectType = Image\nNDims = 3\nDimSize = 4 3 2\nBinaryData = True\n"
        "BinaryDataByteOrderMSB = False\nCompressedData = False\nOffset = -1 -2 0.5\n"
        "ElementSpacing = 1 1 1\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n"
    )
    p = tmp_path / "local.mhd"
    p.write_bytes(head.encode() + arr.tobytes() + b"")
    img = read_metaimage(p)
    np.testing.assert_array_equal(img.data, arr)
    assert img.offset == (-1.0, -2.0, 0.5)
    p.write_bytes(head.encode() + arr.tobytes()[:-4])
    with pytest.raises(ParseError):
        read_metaimage(p)


def test_topas_csv_count_scorer_without_unit(tmp_path: Path) -> None:
    p = tmp_path / "count.csv"
    p.write_text(
        "# TOPAS Version: 4.3\n# Results for scorer: PrimaryCount\n"
        '# Filtered by: OnlyIncludeParticlesGoing = "In"\n'
        "# Scored on surface: Water/ZMinusSurface\n"
        "# X in 1 bin  of 12 cm\n# Y in 1 bin  of 12 cm\n# Z in 2 bins of 0.1 cm\n"
        "# SurfaceTrackCount : Sum   \n0, 0, 0, 200\n0, 0, 1, 197\n"
    )
    s = read_topas_csv(p)
    assert s.quantity == "SurfaceTrackCount" and s.unit == ""
    assert s.statistics == ["Sum"]
    assert s.meta["scorer"] == "PrimaryCount"
    np.testing.assert_array_equal(s.values["Sum"][0, 0], [200, 197])


def test_topas_csv_nested_parentheses_in_unit(tmp_path: Path) -> None:
    p = tmp_path / "let.csv"
    p.write_text(
        "# Results for scorer: LETd\n# Warning: bins set to 0: 119\n"
        "# X in 1 bin  of 12 cm\n# Y in 1 bin  of 12 cm\n# Z in 1 bin of 0.1 cm\n"
        "# ProtonLET ( MeV/mm/(g/cm3) ) : Sum   \n0, 0, 0, 0.52\n"
    )
    s = read_topas_csv(p)
    assert (s.quantity, s.unit) == ("ProtonLET", "MeV/mm/(g/cm3)")
