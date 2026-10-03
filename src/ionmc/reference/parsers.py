"""Parsers for reference-engine output files (numpy only).

TOPAS CSV scorer output
    '#' header comment lines followed by rows ``ix, iy, iz, v1, v2, ...``. The header gives the
    bin structure per axis ("X in 60 bins of 2 mm"), the quantity with its unit and the ordered
    statistics ("DoseToMedium ( Gy ) : Sum  Standard_Deviation"). Only those header facts are
    relied upon; row indices (0-based) place values, so row order does not matter.

MetaImage (.mhd with raw data file)
    ``ElementType = MET_FLOAT`` (and other fixed-width numeric types), little- or big-endian
    as stated by ``ElementByteOrderMSB``. Arrays are returned in numpy order ``(nz, ny, nx)``
    (x fastest in the file); ``MetaImage.dims`` keeps the header ``(nx, ny, nz)`` order.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


class ParseError(ValueError):
    """Raised when a reference output file does not have the expected structure."""


_AXIS = re.compile(
    r"^#\s*([XYZ]|R|Phi|Theta)\s+in\s+(\d+)\s+bins?\s+of\s+([-+0-9.eE]+)\s*(\S*)\s*$"
)
_QUANTITY = re.compile(
    r"^#\s*(?P<name>[^\s(][^(]*?)\s*\(\s*(?P<unit>[^)]*?)\s*\)\s*:\s*(?P<stats>.*)$"
)
_SCORER = re.compile(r"^#\s*Results for scorer\s+(\S+)", re.IGNORECASE)
_COMPONENT = re.compile(r"^#\s*Scored in component:\s*(\S+)", re.IGNORECASE)


@dataclass
class TopasScorer:
    """Dense TOPAS scorer result. ``values[stat]`` has shape ``(nx, ny, nz)``."""

    values: dict[str, np.ndarray]
    bins: tuple[int, int, int]
    bin_width: tuple[float, float, float]
    bin_unit: tuple[str, str, str]
    quantity: str
    unit: str
    statistics: list[str]
    meta: dict[str, Any] = field(default_factory=dict)


def read_topas_csv(path: str | Path) -> TopasScorer:
    """Read a TOPAS CSV scorer file; missing bins (not written by TOPAS) are NaN."""
    header: list[str] = []
    rows: list[list[float]] = []
    for line_no, raw in enumerate(Path(path).read_text().splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            header.append(line)
            continue
        try:
            rows.append([float(tok) for tok in line.split(",")])
        except ValueError as exc:
            raise ParseError(f"{path}:{line_no}: non-numeric data row {line!r}") from exc

    axes: dict[str, tuple[int, float, str]] = {}
    quantity = unit = ""
    statistics: list[str] = []
    meta: dict[str, Any] = {"header": header}
    for line in header:
        m = _AXIS.match(line)
        if m:
            axes[m.group(1).upper()] = (int(m.group(2)), float(m.group(3)), m.group(4))
            continue
        m = _SCORER.match(line)
        if m:
            meta["scorer"] = m.group(1)
            continue
        m = _COMPONENT.match(line)
        if m:
            meta["component"] = m.group(1)
            continue
        m = _QUANTITY.match(line)
        if m and not statistics:
            quantity, unit = m.group("name").strip(), m.group("unit").strip()
            statistics = m.group("stats").replace(",", " ").split()
    if not all(a in axes for a in "XYZ"):
        raise ParseError(f"{path}: header lacks X/Y/Z bin structure")
    if not statistics:
        raise ParseError(f"{path}: header lacks 'Quantity ( unit ) : statistics' line")
    if not rows:
        raise ParseError(f"{path}: no data rows")

    nx, ny, nz = (axes[a][0] for a in "XYZ")
    ncol = 3 + len(statistics)
    if any(len(r) != ncol for r in rows):
        raise ParseError(f"{path}: every row must have {ncol} columns (ix,iy,iz + {statistics})")
    data = np.asarray(rows, dtype=np.float64)
    idx = data[:, :3].astype(np.int64)
    if (
        np.any(idx < 0)
        or np.any(idx[:, 0] >= nx)
        or np.any(idx[:, 1] >= ny)
        or np.any(idx[:, 2] >= nz)
    ):
        raise ParseError(f"{path}: bin index outside the declared bin structure")
    values: dict[str, np.ndarray] = {}
    for k, name in enumerate(statistics):
        arr = np.full((nx, ny, nz), np.nan)
        arr[idx[:, 0], idx[:, 1], idx[:, 2]] = data[:, 3 + k]
        values[name] = arr
    return TopasScorer(
        values=values,
        bins=(nx, ny, nz),
        bin_width=tuple(axes[a][1] for a in "XYZ"),  # type: ignore[arg-type]
        bin_unit=tuple(axes[a][2] for a in "XYZ"),  # type: ignore[arg-type]
        quantity=quantity,
        unit=unit,
        statistics=statistics,
        meta=meta,
    )


@dataclass
class MetaImage:
    """MetaImage volume: ``data`` has numpy shape ``(nz, ny, nx)``."""

    data: np.ndarray
    dims: tuple[int, int, int]
    spacing: tuple[float, float, float]
    offset: tuple[float, float, float]
    header: dict[str, str]


_MET_TYPES = {
    "MET_FLOAT": "f4",
    "MET_DOUBLE": "f8",
    "MET_SHORT": "i2",
    "MET_USHORT": "u2",
    "MET_INT": "i4",
    "MET_UINT": "u4",
    "MET_UCHAR": "u1",
    "MET_CHAR": "i1",
}


def _floats(text: str, n: int, default: float) -> tuple[float, ...]:
    if not text:
        return (default,) * n
    vals = tuple(float(t) for t in text.split())
    if len(vals) != n:
        raise ParseError(f"expected {n} numbers, got {text!r}")
    return vals


def read_metaimage(path: str | Path) -> MetaImage:
    """Read a 3D MetaImage ``.mhd`` and its raw data file (same directory unless absolute)."""
    path = Path(path)
    header: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, _, value = line.partition("=")
            header[key.strip()] = value.strip()
    for key in ("NDims", "DimSize", "ElementType", "ElementDataFile"):
        if key not in header:
            raise ParseError(f"{path}: missing {key}")
    if header["NDims"] != "3":
        raise ParseError(f"{path}: only NDims = 3 is supported")
    if header["ElementType"] not in _MET_TYPES:
        raise ParseError(f"{path}: unsupported ElementType {header['ElementType']}")
    if header.get("CompressedData", "False").lower() == "true":
        raise ParseError(f"{path}: compressed MetaImage data are not supported")
    if int(header.get("ElementNumberOfChannels", "1")) != 1:
        raise ParseError(f"{path}: multi-channel MetaImage is not supported")
    dims = tuple(int(t) for t in header["DimSize"].split())
    if len(dims) != 3:
        raise ParseError(f"{path}: DimSize must have 3 entries")
    msb = header.get("ElementByteOrderMSB", header.get("BinaryDataByteOrderMSB", "False"))
    dtype = np.dtype(("<" if msb.lower() != "true" else ">") + _MET_TYPES[header["ElementType"]])
    raw_name = header["ElementDataFile"]
    if raw_name == "LOCAL":
        raise ParseError(f"{path}: ElementDataFile = LOCAL is not supported")
    raw_path = Path(raw_name)
    if not raw_path.is_absolute():
        raw_path = path.parent / raw_path
    count = dims[0] * dims[1] * dims[2]
    data = np.fromfile(raw_path, dtype=dtype)
    if data.size != count:
        raise ParseError(f"{raw_path}: expected {count} elements, found {data.size}")
    spacing = _floats(header.get("ElementSpacing", header.get("ElementSize", "")), 3, 1.0)
    offset = _floats(header.get("Offset", header.get("Position", "")), 3, 0.0)
    return MetaImage(
        data=data.reshape(dims[2], dims[1], dims[0]),
        dims=(dims[0], dims[1], dims[2]),
        spacing=(spacing[0], spacing[1], spacing[2]),
        offset=(offset[0], offset[1], offset[2]),
        header=header,
    )
