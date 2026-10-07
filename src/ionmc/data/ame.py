"""Parser for the AME2020 atomic mass table (``mass_1.mas20.txt``, Wang et al. 2021).

Fixed-width Fortran layout (the file's own header gives it)::

    a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,f10.5,1x,a2,f13.5,f11.5,1x,i3,1x,f13.6,f12.6
    cc NZ  N  Z  A    el  o     mass  unc binding unc      B  beta  unc    atomic_mass   unc

Used here: ``Z``, ``A``, the element symbol, the mass excess [keV] and its uncertainty, and the
atomic mass [micro-u, written as an integer millionths part ``i3`` followed by ``f13.6``] with
its uncertainty. A ``#`` replaces the decimal point of values estimated from systematic trends
(not measured); such entries are flagged ``estimated``. The header (about 36 lines) is skipped;
after the first data line every non-blank line must parse (fail closed).
"""

from __future__ import annotations

from dataclasses import dataclass

U_MEV = 931.49410242  # atomic mass unit [MeV/c^2], CODATA 2018 (NIST SP 961 / PDG 2020)
ELECTRON_MASS_MEV = 0.51099895  # electron mass [MeV/c^2], CODATA 2018
_MIN_LINE = 123


class AmeError(ValueError):
    """Malformed AME2020 text."""


@dataclass(frozen=True)
class AmeEntry:
    """One nuclide: mass excess [keV], atomic mass [u] (neutral atom), uncertainties, flag."""

    z: int
    a: int
    element: str
    mass_excess_kev: float
    err_kev: float
    atomic_mass_u: float
    err_u: float
    estimated: bool


def _number(text: str) -> tuple[float, bool]:
    s = text.strip()
    estimated = "#" in s
    return float(s.replace("#", ".")), estimated


def _parse_line(line: str) -> AmeEntry:
    line = line.ljust(135)
    n, z, a = int(line[4:9]), int(line[9:14]), int(line[14:19])
    if z + n != a:
        raise AmeError(f"Z + N != A in {line!r}")
    me, est1 = _number(line[28:42])
    me_err, est2 = _number(line[42:54])
    micro_hi = int(line[106:109])
    micro_lo, est3 = _number(line[110:123])
    amu_err, est4 = _number(line[123:135])
    return AmeEntry(
        z=z,
        a=a,
        element=line[20:23].strip(),
        mass_excess_kev=me,
        err_kev=me_err,
        atomic_mass_u=(micro_hi * 1.0e6 + micro_lo) * 1.0e-6,
        err_u=amu_err * 1.0e-6,
        estimated=est1 or est2 or est3 or est4,
    )


def load_ame2020(text: str) -> dict[tuple[int, int], AmeEntry]:
    """Parse the table text into ``{(Z, A): AmeEntry}``; duplicates or bad lines raise ``AmeError``."""
    entries: dict[tuple[int, int], AmeEntry] = {}
    started = False
    for line in text.splitlines():
        if not line.strip():
            continue
        if len(line.rstrip()) < _MIN_LINE and not started:
            continue
        try:
            entry = _parse_line(line)
        except ValueError as exc:
            if not started:
                continue
            raise AmeError(f"cannot parse AME2020 line {line!r}: {exc}") from None
        started = True
        key = (entry.z, entry.a)
        if key in entries:
            raise AmeError(f"duplicate nuclide {key}")
        entries[key] = entry
    if not entries:
        raise AmeError("no AME2020 data lines found")
    return entries


def nuclear_mass_mev(table: dict[tuple[int, int], AmeEntry], z: int, a: int) -> float:
    """Bare nuclear mass [MeV/c^2] = atomic mass x u - Z m_e.

    The total binding energy of the Z electrons (about 14.4 Z^2.39 eV, 2 keV for oxygen, 0.1 MeV
    for lead) is ignored: it is far below the uncertainties of the physics it feeds.
    """
    entry = table[(z, a)]
    return entry.atomic_mass_u * U_MEV - z * ELECTRON_MASS_MEV
