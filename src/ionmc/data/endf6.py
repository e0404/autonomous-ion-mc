"""Parser for ENDF-6 evaluated nuclear data files (proton sublibrary, MF1/MF3/MF6 subset).

Basis: ENDF-6 Formats Manual (BNL-203218-2018-INRE), chapters 0-3 (record types, MF1, MF3) and
chapter 6 (MF6 product energy-angle distributions). Units are those of the format: energies in
eV, cross sections in barns.

Records are 80 columns: six 11-character fields (columns 1-66), MAT (67-70), MF (71-72), MT
(73-75) and the line number NS (76-80). Floats are written without the ``E`` (``9.986200-1``);
integers are right-aligned in 11-character fields. The record types read here are

* HEAD/CONT  ``C1 C2 L1 L2 N1 N2`` (HEAD: ``ZA AWR L1 L2 N1 N2``),
* LIST       CONT followed by ``N1`` floats (six per line),
* TAB1       CONT (``N1 = NR``, ``N2 = NP``), ``NR`` (NBT, INT) integer pairs, ``NP`` (x, y) pairs,
* TAB2       CONT (``N1 = NR``, ``N2 = NZ``), ``NR`` (NBT, INT) integer pairs.

Interpolation laws (INT): 1 histogram, 2 lin-lin, 3 lin-log (y linear in ln x), 4 log-lin
(ln y linear in x), 5 log-log. MF6 is supported for LAW=1 (continuum energy-angle distribution)
with LANG=1 (Legendre coefficients) and LANG=2 (Kalbach-Mann, NA=1: ``b_0 = f_0``, ``b_1 = r``);
LAW=5 (charged-particle elastic scattering) is recognised and stored opaque; any other LAW or LANG
raises :class:`UnsupportedEndfError` (fail closed). Files are untrusted data: they are parsed only,
never executed.
"""

from __future__ import annotations

import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray

_FIELD = 11
_NFIELDS = 6
_MAX_MEMBER_BYTES = 256 * 1024 * 1024
MEMBER_PATTERN = re.compile(r"^ENDF-B-VIII\.0_protons/p-\d{3}_[A-Za-z]{1,2}_\d{3}(?:m\d)?\.endf$")
_MANTISSA_EXPONENT = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$")


class EndfError(ValueError):
    """Malformed ENDF-6 text."""


class UnsupportedEndfError(EndfError):
    """A valid ENDF-6 construct outside the supported subset (fail closed)."""


def parse_endf_float(field_text: str) -> float:
    """Parse one 11-character ENDF float field: ``1.234567+3``, ``1.2E+03``, ``0.5``; blank is 0."""
    s = field_text.strip()
    if not s:
        return 0.0
    try:
        return float(s.replace("D", "E").replace("d", "e"))
    except ValueError:
        pass
    match = _MANTISSA_EXPONENT.match(s)
    if match is None:
        raise EndfError(f"cannot parse ENDF float {field_text!r}")
    return float(f"{match.group(1)}e{match.group(2)}")


def parse_endf_int(field_text: str) -> int:
    """Parse one 11-character ENDF integer field (blank is 0)."""
    s = field_text.strip()
    if not s:
        return 0
    try:
        return int(s)
    except ValueError:
        raise EndfError(f"cannot parse ENDF integer {field_text!r}") from None


class Cont(NamedTuple):
    """CONT/HEAD record: ``c1, c2`` floats and ``l1, l2, n1, n2`` integers."""

    c1: float
    c2: float
    l1: int
    l2: int
    n1: int
    n2: int


def _cont_from_line(line: str) -> Cont:
    f = [line[i * _FIELD : (i + 1) * _FIELD] for i in range(_NFIELDS)]
    return Cont(
        parse_endf_float(f[0]),
        parse_endf_float(f[1]),
        parse_endf_int(f[2]),
        parse_endf_int(f[3]),
        parse_endf_int(f[4]),
        parse_endf_int(f[5]),
    )


def _lines_for(n_values: int) -> int:
    return -(-n_values // _NFIELDS)


class _Reader:
    """Sequential record reader over the 66-character data parts of one section."""

    def __init__(self, lines: tuple[str, ...]) -> None:
        self._lines = lines
        self._pos = 0

    @property
    def exhausted(self) -> bool:
        return self._pos >= len(self._lines)

    def _next_line(self) -> str:
        if self._pos >= len(self._lines):
            raise EndfError("unexpected end of section")
        line = self._lines[self._pos].ljust(_FIELD * _NFIELDS)
        self._pos += 1
        return line

    def cont(self) -> Cont:
        return _cont_from_line(self._next_line())

    def text(self) -> str:
        return self._next_line().rstrip()

    def _floats(self, n: int) -> list[float]:
        out: list[float] = []
        for _ in range(_lines_for(n)):
            line = self._next_line()
            out.extend(
                parse_endf_float(line[i * _FIELD : (i + 1) * _FIELD]) for i in range(_NFIELDS)
            )
        return out[:n]

    def _ints(self, n: int) -> list[int]:
        out: list[int] = []
        for _ in range(_lines_for(n)):
            line = self._next_line()
            out.extend(parse_endf_int(line[i * _FIELD : (i + 1) * _FIELD]) for i in range(_NFIELDS))
        return out[:n]

    def list_body(self, head: Cont) -> NDArray[np.float64]:
        """The ``N1`` floats of a LIST record whose CONT was already read."""
        return np.asarray(self._floats(head.n1), dtype=np.float64)

    def interpolation(self, nr: int) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
        if nr < 1:
            raise EndfError(f"NR must be >= 1, got {nr}")
        pairs = self._ints(2 * nr)
        return (
            np.asarray(pairs[0::2], dtype=np.int64),
            np.asarray(pairs[1::2], dtype=np.int64),
        )

    def tab1(self) -> tuple[Cont, Tab1]:
        head = self.cont()
        return head, self.tab1_body(head)

    def tab1_body(self, head: Cont) -> Tab1:
        nbt, law = self.interpolation(head.n1)
        flat = np.asarray(self._floats(2 * head.n2), dtype=np.float64)
        return Tab1(nbt=nbt, interp=law, x=flat[0::2].copy(), y=flat[1::2].copy())

    def tab2(self) -> tuple[Cont, NDArray[np.int64], NDArray[np.int64]]:
        head = self.cont()
        nbt, law = self.interpolation(head.n1)
        return head, nbt, law


@dataclass(frozen=True)
class Tab1:
    """A TAB1 function y(x): ``nbt`` (1-based last point of each region), ``interp`` (INT 1-5)."""

    nbt: NDArray[np.int64]
    interp: NDArray[np.int64]
    x: NDArray[np.float64]
    y: NDArray[np.float64]

    def __post_init__(self) -> None:
        if self.x.shape != self.y.shape or self.x.ndim != 1 or self.x.size < 2:
            raise EndfError("TAB1 needs matching 1-D x and y with at least two points")
        if np.any(np.diff(self.x) < 0.0):
            raise EndfError("TAB1 abscissae must be non-decreasing")
        if int(self.nbt[-1]) != self.x.size:
            raise EndfError("last NBT must equal the number of points NP")
        if np.any(np.diff(self.nbt) <= 0) or np.any((self.interp < 1) | (self.interp > 5)):
            raise EndfError("invalid TAB1 interpolation table")

    def interpolate(self, x: float | NDArray[np.float64]) -> NDArray[np.float64]:
        """Evaluate y at ``x`` (vectorised) with the law of the region containing each interval.

        The interval between points ``i`` and ``i+1`` (0-based) belongs to the first region whose
        NBT is at least ``i + 2``. At a repeated abscissa (a discontinuity) the value of the
        upper interval is returned. ``x`` outside ``[x_0, x_{NP-1}]`` raises ``ValueError``.
        """
        q = np.atleast_1d(np.asarray(x, dtype=np.float64))
        if np.any(q < self.x[0]) or np.any(q > self.x[-1]) or np.any(np.isnan(q)):
            raise ValueError("argument outside the tabulated range of the TAB1")
        n = self.x.size
        i = np.clip(np.searchsorted(self.x, q, side="right") - 1, 0, n - 2)
        region = np.searchsorted(self.nbt, i + 2, side="left")
        law = self.interp[region]
        x0, x1 = self.x[i], self.x[i + 1]
        y0, y1 = self.y[i], self.y[i + 1]
        out = np.empty_like(q)
        with np.errstate(divide="ignore", invalid="ignore"):
            dx = x1 - x0
            lin = np.where(dx > 0.0, (q - x0) / dx, 1.0)
            out = np.where(law == 1, y0, out)
            out = np.where(law == 2, y0 + (y1 - y0) * lin, out)
            if np.any(law == 3):
                if np.any(q[law == 3] <= 0.0) or np.any(x0[law == 3] <= 0.0):
                    raise ValueError("lin-log interpolation needs positive abscissae")
                t3 = np.where(dx > 0.0, np.log(q / x0) / np.log(x1 / x0), 1.0)
                out = np.where(law == 3, y0 + (y1 - y0) * t3, out)
            if np.any(law == 4):
                if np.any(y0[law == 4] <= 0.0) or np.any(y1[law == 4] <= 0.0):
                    raise ValueError("log-lin interpolation needs positive ordinates")
                out = np.where(law == 4, np.exp(np.log(y0) + lin * (np.log(y1) - np.log(y0))), out)
            if np.any(law == 5):
                if np.any((y0[law == 5] <= 0.0) | (y1[law == 5] <= 0.0) | (q[law == 5] <= 0.0)):
                    raise ValueError("log-log interpolation needs positive abscissae and ordinates")
                t5 = np.where(dx > 0.0, np.log(q / x0) / np.log(x1 / x0), 1.0)
                out = np.where(law == 5, np.exp(np.log(y0) + t5 * (np.log(y1) - np.log(y0))), out)
        # tabulated nodes are returned exactly (exp/log round trips lose an ulp)
        out = np.where((law > 1) & (q == x0), y0, out)
        out = np.where((law > 1) & (q == x1), y1, out)
        return out


@dataclass(frozen=True)
class EnergyDistribution:
    """LAW=1 data for one incident energy: ``rows[k] = (E', b_0, ..., b_NA)`` for ``k < NEP``."""

    energy: float
    nd: int
    na: int
    nw: int
    nep: int
    rows: NDArray[np.float64]


@dataclass(frozen=True)
class Product:
    """One MF6 product. ``lang``, ``lep`` and ``distributions`` are set for LAW=1 only;
    ``opaque_records`` holds the raw 66-character data lines of a LAW=5 distribution."""

    zap: int
    awp: float
    lip: int
    law: int
    yield_: Tab1
    lang: int | None = None
    lep: int | None = None
    nr: int | None = None
    ne: int | None = None
    energy_nbt: NDArray[np.int64] | None = None
    energy_interp: NDArray[np.int64] | None = None
    distributions: tuple[EnergyDistribution, ...] = ()
    opaque_records: tuple[str, ...] = ()


@dataclass(frozen=True)
class Mf6Section:
    """MF6 section: ``lct`` (reference frame: 1 lab, 2 CM, 3 light products CM, recoils lab)
    and products."""

    za: float
    awr: float
    lct: int
    products: tuple[Product, ...]


@dataclass(frozen=True)
class EndfSection:
    """Raw data (columns 1-66) of the records of one (MF, MT) section, SEND excluded."""

    mf: int
    mt: int
    lines: tuple[str, ...]


@dataclass
class EndfMaterial:
    """A parsed material: ``za``, ``awr``, ``emax`` [eV, from MF1/MT451 if present], sections."""

    za: float
    awr: float
    emax: float | None
    mat: int
    sections: dict[tuple[int, int], EndfSection] = field(default_factory=dict)

    def cross_section(self, mt: int) -> Tab1:
        """MF3 section ``mt`` as a :class:`Tab1` (E in eV, sigma in b)."""
        return parse_mf3(self.sections[(3, mt)])

    def products(self, mt: int) -> Mf6Section:
        """MF6 section ``mt``.

        Raises :class:`UnsupportedEndfError` outside the supported subset.
        """
        return parse_mf6(self.sections[(6, mt)])


def parse_mf3(section: EndfSection) -> Tab1:
    """MF3 section: HEAD (ZA, AWR) then TAB1 (QM, QI, 0, LR, NR, NP)."""
    if section.mf != 3:
        raise EndfError("not an MF3 section")
    reader = _Reader(section.lines)
    reader.cont()
    _, tab = reader.tab1()
    return tab


def _parse_law1(reader: _Reader) -> dict[str, object]:
    head, nbt, interp = reader.tab2()
    lang, lep, ne = head.l1, head.l2, head.n2
    if lang not in (1, 2):
        raise UnsupportedEndfError(f"MF6 LAW=1 LANG={lang} is not supported (only 1 and 2)")
    if lep not in (1, 2):
        raise UnsupportedEndfError(f"MF6 LAW=1 LEP={lep} is not supported (only 1 and 2)")
    dists = []
    for _ in range(ne):
        lhead = reader.cont()
        nd, na, nw, nep = lhead.l1, lhead.l2, lhead.n1, lhead.n2
        if nw != nep * (na + 2):
            raise EndfError(f"LIST NW={nw} != NEP*(NA+2) = {nep * (na + 2)}")
        if lang == 2 and na != 1:
            raise UnsupportedEndfError(f"Kalbach-Mann (LANG=2) needs NA=1, got NA={na}")
        body = reader.list_body(lhead)
        dists.append(
            EnergyDistribution(
                energy=lhead.c2, nd=nd, na=na, nw=nw, nep=nep, rows=body.reshape(nep, na + 2)
            )
        )
    return {
        "lang": lang,
        "lep": lep,
        "nr": head.n1,
        "ne": ne,
        "energy_nbt": nbt,
        "energy_interp": interp,
        "distributions": tuple(dists),
    }


def _parse_law5(reader: _Reader) -> tuple[str, ...]:
    """LAW=5 (charged-particle elastic): TAB2 over energies, a LIST per energy; stored opaque."""
    start = reader._pos
    head, _, _ = reader.tab2()
    for _ in range(head.n2):
        lhead = reader.cont()
        reader.list_body(lhead)
    return reader._lines[start : reader._pos]


def parse_mf6(section: EndfSection) -> Mf6Section:
    """MF6 section: HEAD (ZA, AWR, JP, LCT, NK) and NK products (TAB1 yield + law data)."""
    if section.mf != 6:
        raise EndfError("not an MF6 section")
    reader = _Reader(section.lines)
    head = reader.cont()
    products = []
    for _ in range(head.n1):
        phead, yield_tab = reader.tab1()
        law = phead.l2
        common = {
            "zap": int(round(phead.c1)),
            "awp": phead.c2,
            "lip": phead.l1,
            "law": law,
            "yield_": yield_tab,
        }
        if law == 1:
            products.append(Product(**common, **_parse_law1(reader)))  # type: ignore[arg-type]
        elif law == 5:
            products.append(Product(**common, opaque_records=_parse_law5(reader)))  # type: ignore[arg-type]
        else:
            raise UnsupportedEndfError(f"MF6 LAW={law} is not supported (only 1; 5 opaque)")
    if not reader.exhausted:
        raise EndfError("unexpected records after the last MF6 product")
    return Mf6Section(za=head.c1, awr=head.c2, lct=head.l2, products=tuple(products))


def parse_endf(text: str) -> EndfMaterial:
    """Parse ENDF-6 text into an :class:`EndfMaterial` (sections kept as raw records).

    The TPID line (MF 0), FEND/MEND/TEND lines and SEND lines (MT 0) are consumed. Exactly one
    MAT is allowed. ``emax`` is the upper energy limit of MF1/MT451 (C2 of its third record).
    """
    sections: dict[tuple[int, int], list[str]] = {}
    mats: set[int] = set()
    for raw in text.splitlines():
        if len(raw.rstrip()) < 72:
            if not raw.strip():
                continue
            raise EndfError(f"short record: {raw!r}")
        try:
            mat, mf, mt = int(raw[66:70]), int(raw[70:72]), int(raw[72:75])
        except ValueError:
            raise EndfError(f"cannot read MAT/MF/MT from {raw!r}") from None
        if mf == 0 or mt == 0 or mat <= 0:
            continue
        mats.add(mat)
        sections.setdefault((mf, mt), []).append(raw[:66])
    if not sections:
        raise EndfError("no data sections")
    if len(mats) != 1:
        raise EndfError(f"expected exactly one MAT, found {sorted(mats)}")
    first = _cont_from_line(next(iter(sections.values()))[0])
    emax = None
    if (1, 451) in sections:
        r = _Reader(tuple(sections[(1, 451)]))
        r.cont()
        r.cont()
        emax = r.cont().c2
    return EndfMaterial(
        za=first.c1,
        awr=first.c2,
        emax=emax,
        mat=next(iter(mats)),
        sections={k: EndfSection(k[0], k[1], tuple(v)) for k, v in sections.items()},
    )


def _checked_member(name: str) -> str:
    if not MEMBER_PATTERN.match(name):
        raise ValueError(f"not a valid ENDF proton-sublibrary member name: {name!r}")
    return name


def list_zip_members(zip_path: str | Path) -> list[str]:
    """Names of the ``p-ZZZ_El_AAA.endf`` members of the sublibrary zip (others are ignored)."""
    with zipfile.ZipFile(zip_path) as zf:
        return sorted(n for n in zf.namelist() if MEMBER_PATTERN.match(n))


def read_member(zip_path: str | Path, name: str) -> str:
    """Read one validated member as ASCII text (no path traversal; size-capped)."""
    _checked_member(name)
    with zipfile.ZipFile(zip_path) as zf:
        info = zf.getinfo(name)
        if info.file_size > _MAX_MEMBER_BYTES:
            raise EndfError(f"member {name} is larger than {_MAX_MEMBER_BYTES} bytes")
        return zf.read(info).decode("ascii")


__all__ = [
    "Cont",
    "EndfError",
    "EndfMaterial",
    "EndfSection",
    "EnergyDistribution",
    "Mf6Section",
    "Product",
    "Tab1",
    "UnsupportedEndfError",
    "list_zip_members",
    "parse_endf",
    "parse_endf_float",
    "parse_endf_int",
    "parse_mf3",
    "parse_mf6",
    "read_member",
]
