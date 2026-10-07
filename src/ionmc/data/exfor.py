"""Parser for EXFOR master-entry text (one ``ENTRY`` with its ``SUBENT`` subentries).

EXFOR lines are 80 columns: 66 data columns (keyword in columns 1-10, text from column 12),
then the entry/subentry id and line number. Supported blocks:

* ``BIB``: keyword blocks with continuation lines (keyword field blank). ``REACTION`` entries may
  carry a pointer character in column 11 (``REACTION  1(6-C-12(P,NON),,SIG)``);
* ``COMMON``: ``N1`` field head line, unit line and value line(s), 11-character fields, at most
  six fields (fail closed otherwise);
* ``DATA``: ``N1`` fields per row, ``N2`` rows, a head and a unit line, then rows of fixed
  11-character fields (blank is ``None``). A pointer character in column 11 of a head field is
  kept in ``DataTable.pointers``.

Numbers may be written with an exponent without ``E`` (``1.5+3``). Units handled by the helpers:
energy ``MEV``, ``KEV``, ``MEV/A`` (per nucleon, flagged), cross section ``B``, ``MB``,
``MUB``, uncertainty ``PER-CENT``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

_FIELD = 11
_DATA_COLS = 66
_NUMBER = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+))(?:[eE]?([+-]?\d+))?$")
_ENERGY_TO_MEV = {"MEV": 1.0, "KEV": 1.0e-3, "MEV/A": 1.0}
_XS_TO_MB = {"MB": 1.0, "B": 1.0e3, "MUB": 1.0e-3}


class ExforError(ValueError):
    """Malformed or unsupported EXFOR text."""


def parse_exfor_number(text: str) -> float | None:
    """Parse one EXFOR numeric field (``12.5``, ``1.5+3``, ``2.E-2``); blank is ``None``."""
    s = text.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        pass
    m = _NUMBER.match(s)
    if m is None or m.group(2) is None:
        raise ExforError(f"cannot parse EXFOR number {text!r}")
    return float(f"{m.group(1)}e{m.group(2)}")


@dataclass(frozen=True)
class DataTable:
    """DATA block: ``heads``, ``units``, ``pointers`` (column-11 characters) and ``rows``."""

    heads: tuple[str, ...]
    units: tuple[str, ...]
    pointers: tuple[str, ...]
    rows: tuple[tuple[float | None, ...], ...]

    def column(self, index: int) -> NDArray[np.float64]:
        """Column ``index`` as float64 (``None`` becomes NaN)."""
        return np.asarray(
            [np.nan if r[index] is None else r[index] for r in self.rows], dtype=np.float64
        )

    def columns_named(self, head: str) -> list[int]:
        """Indices of the columns whose head is ``head`` (several with pointers)."""
        return [i for i, h in enumerate(self.heads) if h == head]


@dataclass(frozen=True)
class CommonBlock:
    """COMMON block: ``heads``, ``units`` and ``values`` of the quantities shared by all rows."""

    heads: tuple[str, ...]
    units: tuple[str, ...]
    values: tuple[float | None, ...]


@dataclass(frozen=True)
class Reaction:
    """One REACTION entry: ``pointer`` (column-11 character, blank if none) and ``text``."""

    pointer: str
    text: str


@dataclass(frozen=True)
class Subentry:
    """One subentry: id, BIB keys (continuation lines kept), COMMON and DATA."""

    subentry_id: str
    bib: dict[str, tuple[str, ...]]
    reactions: tuple[Reaction, ...]
    common: CommonBlock | None
    data: DataTable | None

    @property
    def title(self) -> str:
        return " ".join(s.strip() for s in self.bib.get("TITLE", ()))

    @property
    def first_author(self) -> str | None:
        """First name of the AUTHOR list (parentheses removed), or None."""
        lines = self.bib.get("AUTHOR")
        if not lines:
            return None
        text = " ".join(s.strip() for s in lines).strip().lstrip("(")
        return text.split(",")[0].strip().rstrip(")") or None

    @property
    def reference_year(self) -> int | None:
        """Four-digit year ending the first REFERENCE group (``(J,PR/C,71,064606,2005)``)."""
        lines = self.bib.get("REFERENCE")
        if not lines:
            return None
        text = " ".join(s.strip() for s in lines)
        group = text[text.find("(") + 1 :].split(")")[0]
        tail = group.split(",")[-1].strip()
        return int(tail) if re.fullmatch(r"\d{4}", tail) else None


@dataclass(frozen=True)
class ExforEntry:
    """A whole entry: ``entry_id`` and its subentries in file order."""

    entry_id: str
    subentries: tuple[Subentry, ...] = field(default_factory=tuple)

    @property
    def first_author(self) -> str | None:
        """First author from the first subentry that lists AUTHOR (normally subentry 001)."""
        return next((s.first_author for s in self.subentries if s.first_author), None)

    @property
    def reference_year(self) -> int | None:
        """Publication year from the first subentry with a parseable REFERENCE."""
        return next((s.reference_year for s in self.subentries if s.reference_year), None)

    @property
    def title(self) -> str:
        return next((s.title for s in self.subentries if s.title), "")


def _fields(line: str, n: int) -> list[str]:
    line = line.ljust(_FIELD * 6)
    return [line[i * _FIELD : (i + 1) * _FIELD] for i in range(n)]


def _block_counts(header: str) -> tuple[int, int]:
    parts = header.split()
    if len(parts) < 3:
        raise ExforError(f"malformed block header {header!r}")
    return int(parts[1]), int(parts[2])


def _parse_common(lines: list[str], i: int) -> tuple[CommonBlock, int]:
    n_fields, n_lines = _block_counts(lines[i])
    if n_fields > 6 or n_lines != 3:
        raise ExforError("COMMON blocks of more than six fields are not supported")
    heads, units, values = (_fields(lines[i + k], n_fields) for k in (1, 2, 3))
    return (
        CommonBlock(
            tuple(h.strip() for h in heads),
            tuple(u.strip() for u in units),
            tuple(parse_exfor_number(v) for v in values),
        ),
        i + 1 + n_lines,
    )


def _parse_data(lines: list[str], i: int) -> tuple[DataTable, int]:
    n_fields, n_rows = _block_counts(lines[i])
    if n_fields > 6:
        raise ExforError("DATA blocks of more than six fields are not supported")
    head_f = _fields(lines[i + 1], n_fields)
    unit_f = _fields(lines[i + 2], n_fields)
    rows = tuple(
        tuple(parse_exfor_number(v) for v in _fields(lines[i + 3 + r], n_fields))
        for r in range(n_rows)
    )
    return (
        DataTable(
            heads=tuple(h[: _FIELD - 1].strip() for h in head_f),
            units=tuple(u.strip() for u in unit_f),
            pointers=tuple(h[_FIELD - 1].strip() for h in head_f),
            rows=rows,
        ),
        i + 3 + n_rows,
    )


def _parse_bib(
    lines: list[str], i: int
) -> tuple[dict[str, tuple[str, ...]], tuple[Reaction, ...], int]:
    n_keys, n_lines = _block_counts(lines[i])
    bib: dict[str, list[str]] = {}
    reactions: list[Reaction] = []
    key = ""
    for line in lines[i + 1 : i + 1 + n_lines]:
        line = line.ljust(_DATA_COLS)
        head = line[:10].strip()
        if head:
            key = head
            bib.setdefault(key, [])
        elif not key:
            raise ExforError("BIB continuation line before any keyword")
        text = line[11:_DATA_COLS].rstrip()
        bib[key].append(text)
        if key == "REACTION":
            pointer = line[10]
            if head or pointer.strip() or not reactions:
                reactions.append(Reaction(pointer.strip(), text.strip()))
            else:
                last = reactions[-1]
                reactions[-1] = Reaction(last.pointer, last.text + " " + text.strip())
    if len(bib) != n_keys:
        raise ExforError(f"BIB announces {n_keys} keywords, found {len(bib)}")
    return {k: tuple(v) for k, v in bib.items()}, tuple(reactions), i + 1 + n_lines


def parse_entry(text: str) -> ExforEntry:
    """Parse one EXFOR entry; the first line must be ``ENTRY``."""
    lines = [ln[:_DATA_COLS].rstrip() for ln in text.splitlines()]
    lines = [ln for ln in lines if ln.strip()]
    if not lines or not lines[0].startswith("ENTRY"):
        raise ExforError("text does not start with an ENTRY line")
    entry_id = lines[0].split()[1]
    subs: list[Subentry] = []
    i = 1
    while i < len(lines):
        word = lines[i].split()[0]
        if word == "ENDENTRY":
            break
        if word != "SUBENT":
            raise ExforError(f"expected SUBENT, found {lines[i]!r}")
        sub_id = lines[i].split()[1]
        i += 1
        bib: dict[str, tuple[str, ...]] = {}
        reactions: tuple[Reaction, ...] = ()
        common: CommonBlock | None = None
        data: DataTable | None = None
        while True:
            if i >= len(lines):
                raise ExforError(f"subentry {sub_id} lacks ENDSUBENT")
            word = lines[i].split()[0]
            if word == "BIB":
                bib, reactions, i = _parse_bib(lines, i)
            elif word == "COMMON":
                common, i = _parse_common(lines, i)
            elif word == "DATA":
                data, i = _parse_data(lines, i)
            elif word == "ENDSUBENT":
                i += 1
                break
            elif word in {"ENDBIB", "ENDCOMMON", "ENDDATA", "NOCOMMON", "NODATA"}:
                i += 1
            else:
                raise ExforError(f"unexpected line in subentry {sub_id}: {lines[i]!r}")
        subs.append(Subentry(sub_id, bib, reactions, common, data))
    return ExforEntry(entry_id, tuple(subs))


def energy_to_mev(values: NDArray[np.float64], unit: str) -> tuple[NDArray[np.float64], bool]:
    """Convert energies to MeV; returns ``(values, per_nucleon)``. ``MEV/A`` is NOT multiplied
    by the projectile mass number here: the flag tells the caller to do so."""
    try:
        factor = _ENERGY_TO_MEV[unit]
    except KeyError:
        raise ExforError(f"unsupported energy unit {unit!r}") from None
    return np.asarray(values, dtype=np.float64) * factor, unit == "MEV/A"


def energy_total_mev(
    values: NDArray[np.float64], unit: str, mass_number: int | None
) -> NDArray[np.float64]:
    """Total projectile kinetic energy [MeV]; ``MEV/A`` needs ``mass_number`` (else ExforError)."""
    out, per_nucleon = energy_to_mev(values, unit)
    if per_nucleon:
        if mass_number is None:
            raise ExforError("MEV/A energies need the projectile mass number")
        out = out * mass_number
    return out


def xs_to_mb(values: NDArray[np.float64], unit: str) -> NDArray[np.float64]:
    """Convert cross sections (``B``, ``MB``, ``MUB``) to millibarn."""
    try:
        factor = _XS_TO_MB[unit]
    except KeyError:
        raise ExforError(f"unsupported cross-section unit {unit!r}") from None
    return np.asarray(values, dtype=np.float64) * factor


def percent_to_fraction(values: NDArray[np.float64], unit: str) -> NDArray[np.float64]:
    """Relative uncertainty as a fraction; the unit must be ``PER-CENT``."""
    if unit != "PER-CENT":
        raise ExforError(f"expected PER-CENT, got {unit!r}")
    return np.asarray(values, dtype=np.float64) * 0.01
