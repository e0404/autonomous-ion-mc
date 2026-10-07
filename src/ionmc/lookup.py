"""External lookup tables for lookup tallies (decision 0040, section 7).

A :class:`LookupTable` is a per-species function ``f(x)`` on a **uniform** grid (linear or
logarithmic spacing) of an energy-per-nucleon axis [MeV/u] or a LET-in-water axis [keV/um]. The
transport evaluates it per scoring piece (linear interpolation) and accumulates ``eps f`` into the
channel FE. Values must be finite and non-negative so every channel stays non-negative.

Non-uniform source data are rejected on load. :func:`resample_uniform` is the explicit helper that
turns such data into a uniform table; what it did is recorded in the table provenance
(``resampling``) and ends up in the ``Result``.

**No clinical tables and no biological (RBE) formulas ship with the code.** The only fixture is the
mathematical test function ``tests/data/synthetic_lookup.json`` (``synthetic: true``,
``f = 1 + 0.1 L``) with no biological meaning.

File format (JSON)::

    {"name": ..., "quantity": ..., "units": ..., "axis": "energy_per_nucleon_mev" |
     "let_water_kev_um", "axis_spacing": "linear" | "log", "axis_values": [...],
     "species": {"proton": [...], ...}, "citation": ..., "license": ..., "source": ...,
     "synthetic": true | false}
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ionmc._frozen import freeze_array
from ionmc._validate import fail
from ionmc.species import species_by_name

AXES = ("energy_per_nucleon_mev", "let_water_kev_um")
AXIS_SPACINGS = ("linear", "log")
UNIFORM_RTOL = 1.0e-9
MIN_POINTS = 2
REQUIRED_TEXT_FIELDS = ("name", "quantity", "units", "citation", "license", "source")


def is_uniform(axis: NDArray[np.float64], spacing: str) -> bool:
    """True if ``axis`` is uniform in ``x`` (``linear``) or in ``ln x`` (``log``)."""
    a = np.log(axis) if spacing == "log" else axis
    d = np.diff(a)
    return bool(np.all(np.abs(d - d.mean()) <= UNIFORM_RTOL * abs(float(d.mean()))))


@dataclass(frozen=True, eq=False)
class LookupTable:
    """A validated lookup table (see the module docstring).

    ``values[species_name]`` is a read-only array with one value per axis point. ``file_sha256`` is
    the hash of the file bytes (None for a table built in memory), ``content_sha256`` the hash of
    the canonical content (axis, spacing, values, metadata); ``resampling`` records
    :func:`resample_uniform` (None if the data were used as given).
    """

    name: str
    quantity: str
    units: str
    axis: str
    axis_spacing: str
    axis_values: NDArray[np.float64]
    values: Mapping[str, NDArray[np.float64]]
    citation: str
    license: str
    source: str
    synthetic: bool
    file_sha256: str | None = None
    resampling: Mapping[str, Any] | None = None
    content_sha256: str = field(init=False, default="")

    def __post_init__(self) -> None:
        for f in REQUIRED_TEXT_FIELDS:
            v = getattr(self, f)
            if not isinstance(v, str) or not v.strip():
                raise fail(f"lookup table field {f!r} must be a non-empty string")
        if self.axis not in AXES:
            raise fail(f"lookup axis must be one of {AXES}, got {self.axis!r}")
        if self.axis_spacing not in AXIS_SPACINGS:
            raise fail(f"axis_spacing must be one of {AXIS_SPACINGS}, got {self.axis_spacing!r}")
        if not isinstance(self.synthetic, bool):
            raise fail("synthetic must be a bool")
        axis = np.asarray(self.axis_values, dtype=np.float64)
        if axis.ndim != 1 or axis.size < MIN_POINTS or not np.all(np.isfinite(axis)):
            raise fail(f"axis_values must be a finite 1-D sequence of >= {MIN_POINTS} numbers")
        if np.any(np.diff(axis) <= 0.0):
            raise fail("axis_values must be strictly increasing")
        if axis[0] < 0.0 or (self.axis_spacing == "log" and axis[0] <= 0.0):
            raise fail("axis_values must be positive (log spacing) or non-negative (linear)")
        if not is_uniform(axis, self.axis_spacing):
            raise fail(
                f"lookup table {self.name!r}: the axis is not uniform in "
                f"{'ln x' if self.axis_spacing == 'log' else 'x'} (rtol {UNIFORM_RTOL}); "
                "non-uniform tables are rejected, use resample_uniform() and record it"
            )
        if not self.values:
            raise fail(f"lookup table {self.name!r} has no species")
        frozen_values: dict[str, NDArray[np.float64]] = {}
        for sp, v in self.values.items():
            try:
                transported = species_by_name(sp).transported
            except ValueError as exc:
                raise fail(f"lookup table {self.name!r}: {exc}") from exc
            if not transported:
                raise fail(f"lookup table {self.name!r}: species {sp!r} is not transported")
            arr = np.asarray(v, dtype=np.float64)
            if arr.shape != axis.shape:
                raise fail(f"species {sp!r}: {arr.size} values for {axis.size} axis points")
            if not np.all(np.isfinite(arr)):
                raise fail(f"species {sp!r}: lookup values must be finite")
            if np.any(arr < 0.0):
                raise fail(f"species {sp!r}: lookup values must be non-negative")
            frozen_values[sp] = freeze_array(arr, np.float64, f"values[{sp}]")
        object.__setattr__(self, "axis_values", freeze_array(axis, np.float64, "axis_values"))
        object.__setattr__(self, "values", frozen_values)
        if self.resampling is not None:
            object.__setattr__(self, "resampling", json.loads(json.dumps(dict(self.resampling))))
        object.__setattr__(self, "content_sha256", self._content_hash())

    def _content_hash(self) -> str:
        doc = {
            "name": self.name,
            "quantity": self.quantity,
            "units": self.units,
            "axis": self.axis,
            "axis_spacing": self.axis_spacing,
            "axis_values": self.axis_values.tolist(),
            "species": {k: self.values[k].tolist() for k in sorted(self.values)},
            "synthetic": self.synthetic,
        }
        blob = json.dumps(doc, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(blob).hexdigest()

    @property
    def axis_range(self) -> tuple[float, float]:
        """First and last axis value."""
        return float(self.axis_values[0]), float(self.axis_values[-1])

    def covers(self, lo: float, hi: float) -> bool:
        """True if ``[lo, hi]`` lies inside the axis (tolerance 1e-12 relative)."""
        a, b = self.axis_range
        tol = 1e-12 * max(abs(a), abs(b))
        return lo >= a - tol and hi <= b + tol

    def evaluate(self, species: str, x: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
        """Reference (numpy) evaluation: linear interpolation of ``values[species]`` at ``x`` and a
        mask that is True where ``x`` is inside the axis (outside, the end value is returned and the
        transport would count ``lookup_out_of_domain``)."""
        if species not in self.values:
            raise fail(f"lookup table {self.name!r} has no species {species!r}")
        xs = np.asarray(x, dtype=np.float64)
        a, b = self.axis_range
        inside = (xs >= a) & (xs <= b)
        return np.interp(xs, self.axis_values, self.values[species]), inside

    def provenance(self) -> dict[str, Any]:
        """JSON-serialisable provenance recorded in the ``Result`` (decision 0040, section 7)."""
        return {
            "name": self.name,
            "quantity": self.quantity,
            "units": self.units,
            "axis": self.axis,
            "axis_spacing": self.axis_spacing,
            "axis_range": list(self.axis_range),
            "n_points": int(self.axis_values.size),
            "species": sorted(self.values),
            "file_sha256": self.file_sha256,
            "content_sha256": self.content_sha256,
            "citation": self.citation,
            "license": self.license,
            "source": self.source,
            "synthetic": self.synthetic,
            "resampling": None if self.resampling is None else dict(self.resampling),
        }

    @classmethod
    def from_file(cls, path: str | Path, expected_sha256: str | None = None) -> LookupTable:
        """Read a JSON table; fail closed on a hash mismatch, a missing field, a non-uniform axis,
        an unknown species or a negative or non-finite value."""
        raw = Path(path).read_bytes()
        sha = hashlib.sha256(raw).hexdigest()
        if expected_sha256 is not None and sha != expected_sha256.lower():
            raise fail(
                f"lookup file {str(path)!r} has sha256 {sha}, expected {expected_sha256.lower()}"
            )
        try:
            doc = json.loads(raw)
        except ValueError as exc:
            raise fail(f"lookup file {str(path)!r} is not valid JSON: {exc}") from exc
        if not isinstance(doc, dict):
            raise fail("lookup file must contain a JSON object")
        required = (*REQUIRED_TEXT_FIELDS, "axis", "axis_spacing", "axis_values", "species")
        missing = [k for k in (*required, "synthetic") if k not in doc]
        if missing:
            raise fail(f"lookup file lacks required fields {missing}")
        if not isinstance(doc["species"], dict):
            raise fail("lookup file 'species' must be an object of value arrays")
        return cls(
            name=doc["name"],
            quantity=doc["quantity"],
            units=doc["units"],
            axis=doc["axis"],
            axis_spacing=doc["axis_spacing"],
            axis_values=np.asarray(doc["axis_values"], dtype=np.float64),
            values={k: np.asarray(v, dtype=np.float64) for k, v in doc["species"].items()},
            citation=doc["citation"],
            license=doc["license"],
            source=doc["source"],
            synthetic=doc["synthetic"],
            file_sha256=sha,
        )


def resample_uniform(
    axis: ArrayLike,
    values: Mapping[str, ArrayLike],
    *,
    n_points: int,
    spacing: str,
) -> tuple[NDArray[np.float64], dict[str, NDArray[np.float64]], dict[str, Any]]:
    """Resample non-uniform source data onto a uniform grid (linear in x for ``linear`` spacing,
    linear in ``ln x`` for ``log``) by linear interpolation of the values in x.

    Returns ``(new_axis, new_values, record)``; pass the record as ``resampling=`` of the
    :class:`LookupTable` so the provenance states that, how and from what the data were resampled.
    The new axis spans exactly the source range, so no value is extrapolated.
    """
    if spacing not in AXIS_SPACINGS:
        raise fail(f"spacing must be one of {AXIS_SPACINGS}, got {spacing!r}")
    if isinstance(n_points, bool) or not isinstance(n_points, int) or n_points < MIN_POINTS:
        raise fail(f"n_points must be an integer >= {MIN_POINTS}")
    x = np.asarray(axis, dtype=np.float64)
    if x.ndim != 1 or x.size < MIN_POINTS or not np.all(np.isfinite(x)) or np.any(np.diff(x) <= 0):
        raise fail("source axis must be finite, 1-D and strictly increasing")
    if spacing == "log":
        if x[0] <= 0.0:
            raise fail("log resampling needs a positive axis")
        new = np.exp(np.linspace(math.log(x[0]), math.log(x[-1]), n_points))
        new[0], new[-1] = x[0], x[-1]
    else:
        new = np.linspace(x[0], x[-1], n_points)
    out: dict[str, NDArray[np.float64]] = {}
    for sp, v in values.items():
        arr = np.asarray(v, dtype=np.float64)
        if arr.shape != x.shape:
            raise fail(f"species {sp!r}: {arr.size} values for {x.size} axis points")
        out[sp] = np.interp(new, x, arr)
    record = {
        "method": "linear interpolation of the values in x onto a uniform grid",
        "spacing": spacing,
        "n_points_in": int(x.size),
        "n_points_out": int(n_points),
        "source_axis_sha256": hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest(),
    }
    return new, out, record
