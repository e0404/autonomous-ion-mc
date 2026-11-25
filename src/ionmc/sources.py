"""Particle sources."""

from __future__ import annotations

import hashlib
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc._validate import fail, integer, real, triple
from ionmc.physics.projectiles import Projectile

UNIT_DIRECTION_STATUS_BEAMLETS = "not-applicable: beamlets"
"""Status recorded in place of the single unit direction of a :class:`BeamletSet` (V3-007,
plan section B.8): a set has one direction per beamlet (:meth:`BeamletSet.unit_directions`)."""
CANONICAL_TABLE_COLUMNS = (
    "position_x_mm",
    "position_y_mm",
    "position_z_mm",
    "unit_direction_x",
    "unit_direction_y",
    "unit_direction_z",
    "kinetic_energy_mev",
    "energy_sigma_mev",
    "lateral_sigma_mm",
    "n_histories",
    "has_weight",
    "weight",
)
"""Columns of :meth:`BeamletSet.canonical_table`."""
SIGMA_ENVELOPE = 6.0
"""Half-width, in energy standard deviations, of the energy range of a beamlet (E +- 6 sigma_E),
the same bound the pencil source rules use."""


def unit_vector(d: tuple[float, float, float]) -> tuple[float, float, float]:
    """The normalised vector ``d / |d|`` (float64); the one definition every source uses."""
    n = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
    return (d[0] / n, d[1] / n, d[2] / n)


@dataclass(frozen=True)
class PencilBeamSource:
    """A pencil beam: one position, one direction, Gaussian energy and lateral spread.

    ``position_mm`` is the point of emission [mm]; the beam travels in vacuum until it enters
    the geometry. ``direction`` is any nonzero vector; the normalised value is recorded in
    the effective configuration. ``kinetic_energy_mev`` is the mean total kinetic energy of
    the ion [MeV]; ``energy_sigma_mev`` its Gaussian standard deviation [MeV];
    ``lateral_sigma_mm`` the standard deviation of the position offset in each of the two
    directions perpendicular to the beam [mm]. No axis is privileged.
    """

    projectile: Projectile
    position_mm: tuple[float, float, float]
    direction: tuple[float, float, float]
    kinetic_energy_mev: float
    energy_sigma_mev: float = 0.0
    lateral_sigma_mm: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.projectile, Projectile):
            raise fail(f"projectile must be a Projectile, got {self.projectile!r}")
        object.__setattr__(self, "position_mm", triple("position_mm", self.position_mm))
        d = triple("direction", self.direction)
        if math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) <= 1e-12:
            raise fail("direction must be a nonzero vector")
        object.__setattr__(self, "direction", d)
        real("kinetic_energy_mev", self.kinetic_energy_mev, positive=True)
        real("energy_sigma_mev", self.energy_sigma_mev, nonnegative=True)
        real("lateral_sigma_mm", self.lateral_sigma_mm, nonnegative=True)

    def unit_direction(self) -> tuple[float, float, float]:
        """The normalised direction (float64)."""
        return unit_vector(self.direction)


@dataclass(frozen=True)
class Beamlet:
    """One beamlet of a :class:`BeamletSet`: the parameters of a pencil beam plus its number of
    primaries.

    ``position_mm``, ``direction``, ``kinetic_energy_mev``, ``energy_sigma_mev`` and
    ``lateral_sigma_mm`` mean exactly what they mean for :class:`PencilBeamSource`.
    ``n_histories`` is the number of primaries transported for this beamlet. ``weight`` is the
    planning weight (unit: primaries delivered by the beamlet in the plan), metadata for
    the planning helpers: it never changes sampling, columns or statistics; ``None`` or, checked
    by ``validate``, a finite float >= 0. ``label`` is an optional name (grammar and uniqueness are
    checked by ``validate``). Construction checks types and finiteness only; every range rule of
    the set is a fail-closed rule of ``validate`` (decision 0046).
    """

    position_mm: tuple[float, float, float]
    direction: tuple[float, float, float]
    kinetic_energy_mev: float
    energy_sigma_mev: float = 0.0
    lateral_sigma_mm: float = 0.0
    n_histories: int = 1
    weight: float | None = None
    label: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "position_mm", triple("position_mm", self.position_mm))
        d = triple("direction", self.direction)
        if math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) <= 1e-12:
            raise fail("direction must be a nonzero vector")
        object.__setattr__(self, "direction", d)
        real("kinetic_energy_mev", self.kinetic_energy_mev, positive=True)
        real("energy_sigma_mev", self.energy_sigma_mev, nonnegative=True)
        real("lateral_sigma_mm", self.lateral_sigma_mm, nonnegative=True)
        integer("n_histories", self.n_histories, minimum=1)
        if self.label is not None and not isinstance(self.label, str):
            raise fail(f"label must be a string or None, got {self.label!r}")

    def unit_direction(self) -> tuple[float, float, float]:
        """The normalised direction (float64)."""
        return unit_vector(self.direction)

    @property
    def e_lo_mev(self) -> float:
        """Lower end ``E - 6 sigma_E`` of the energy range of the beamlet [MeV]."""
        return self.kinetic_energy_mev - SIGMA_ENVELOPE * self.energy_sigma_mev

    @property
    def e_hi_mev(self) -> float:
        """Upper end ``E + 6 sigma_E`` of the energy range of the beamlet [MeV]."""
        return self.kinetic_energy_mev + SIGMA_ENVELOPE * self.energy_sigma_mev


@dataclass(frozen=True)
class BeamletEnvelope:
    """What code written for one source sees of a :class:`BeamletSet` (plan section B.1):
    ``e_lo_mev`` / ``e_hi_mev`` the smallest ``E - 6 sigma_E`` and largest ``E + 6 sigma_E`` over
    the beamlets, ``e_cap_mev`` the energy bound of the fixed-point capacity check
    (``max_j (E_j + 6 sigma_j)``, equal to ``e_hi_mev``), ``e0_min_mev`` / ``e0_max_mev`` the
    extreme beamlet energies, ``sigma_e_max_mev`` and ``sigma_lat_max_mm`` the largest spreads."""

    e_lo_mev: float
    e_hi_mev: float
    e_cap_mev: float
    e0_min_mev: float
    e0_max_mev: float
    sigma_e_max_mev: float
    sigma_lat_max_mm: float

    def to_dict(self) -> dict[str, float]:
        """JSON-serialisable form (the ``envelope`` of the effective summary)."""
        return {
            "e_lo_mev": self.e_lo_mev,
            "e_hi_mev": self.e_hi_mev,
            "e_cap_mev": self.e_cap_mev,
            "e0_min_mev": self.e0_min_mev,
            "e0_max_mev": self.e0_max_mev,
            "sigma_e_max_mev": self.sigma_e_max_mev,
            "sigma_lat_max_mm": self.sigma_lat_max_mm,
        }


@dataclass(frozen=True)
class BeamletSet:
    """A spot list: ``beamlets`` (in spot-list order) of one ``projectile`` (decision 0046).

    Beamlet ``j`` owns the global history indices ``[s_j, s_j + n_j)`` with ``s_j`` the running
    sum of the earlier ``n_i`` (:meth:`offsets`); the batch of history ``h`` is ``h mod B`` as for
    every source, so a set whose every ``n_j`` is a multiple of ``B`` puts exactly ``n_j / B``
    histories of every beamlet into each batch. ``RunOptions.n_histories`` must equal the sum of
    the ``n_j``. The set has no single energy, position or direction: those attributes do not
    exist, and code that needs one uses :meth:`envelope` or the per-beamlet values.
    Not a clinical claim.
    """

    projectile: Projectile
    beamlets: tuple[Beamlet, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.projectile, Projectile):
            raise fail(f"projectile must be a Projectile, got {self.projectile!r}")
        if isinstance(self.beamlets, list):
            object.__setattr__(self, "beamlets", tuple(self.beamlets))
        if not isinstance(self.beamlets, tuple) or not all(
            isinstance(b, Beamlet) for b in self.beamlets
        ):
            raise fail("beamlets must be a tuple of Beamlet")

    @property
    def n_beamlets(self) -> int:
        """The number of beamlets J."""
        return len(self.beamlets)

    @property
    def n_histories_total(self) -> int:
        """The sum of the ``n_histories`` of all beamlets (the history count of the run)."""
        return sum(b.n_histories for b in self.beamlets)

    def offsets(self) -> NDArray[np.int64]:
        """``s_j``: int64 array of length J + 1, ``s_0 = 0``, ``s_J`` the total history count."""
        out = np.zeros(len(self.beamlets) + 1, dtype=np.int64)
        np.cumsum([b.n_histories for b in self.beamlets], out=out[1:])
        return out

    def beamlet_index(self, history: Any) -> Any:
        """Beamlet id of a global history index (scalar or array): ``searchsorted(s, h, "right")
        - 1`` on :meth:`offsets`. Raises ``ValueError`` outside ``[0, s_J)`` (never wraps)."""
        s = self.offsets()
        h = np.asarray(history)
        if h.size and (h.min() < 0 or h.max() >= s[-1]):
            raise ValueError(f"history index outside [0, {int(s[-1])})")
        j = np.searchsorted(s, h, side="right") - 1
        return int(j) if j.ndim == 0 else j

    def histories_per_batch(self, n_batches: int) -> NDArray[np.int64]:
        """``n_j / B`` per beamlet (exact when every ``n_j`` is a multiple of ``B``, which
        ``validate`` requires)."""
        n = np.diff(self.offsets())
        if np.any(n % n_batches):
            raise ValueError(f"a beamlet n_histories is not a multiple of n_batches {n_batches}")
        return n // n_batches

    def envelope(self) -> BeamletEnvelope:
        """The envelope of the set (see :class:`BeamletEnvelope`); needs at least one beamlet."""
        return envelope_of(self.beamlets)

    def unit_direction(self) -> None:
        """A set has no single direction: ``None`` (the status is
        :data:`UNIT_DIRECTION_STATUS_BEAMLETS`); see :meth:`unit_directions`."""
        return None

    def unit_directions(self) -> NDArray[np.float64]:
        """The normalised direction of every beamlet, float64 ``[J, 3]``."""
        return np.array([b.unit_direction() for b in self.beamlets], dtype=np.float64).reshape(
            -1, 3
        )

    def canonical_table(self) -> NDArray[np.float64]:
        """The numeric canonical table (see :func:`canonical_table_of`)."""
        return canonical_table_of(self.beamlets)

    def canonical_sha256(self) -> str:
        """sha256 of the canonical encoding (see :func:`canonical_table_sha256`)."""
        return canonical_table_sha256(self.canonical_table(), [b.label for b in self.beamlets])


def envelope_of(beamlets: Sequence[Beamlet]) -> BeamletEnvelope:
    """The :class:`BeamletEnvelope` of ``beamlets`` (at least one)."""
    if not beamlets:
        raise ValueError("the envelope of an empty beamlet set is undefined")
    e_hi = max(b.e_hi_mev for b in beamlets)
    return BeamletEnvelope(
        e_lo_mev=min(b.e_lo_mev for b in beamlets),
        e_hi_mev=e_hi,
        e_cap_mev=e_hi,
        e0_min_mev=min(b.kinetic_energy_mev for b in beamlets),
        e0_max_mev=max(b.kinetic_energy_mev for b in beamlets),
        sigma_e_max_mev=max(b.energy_sigma_mev for b in beamlets),
        sigma_lat_max_mm=max(b.lateral_sigma_mm for b in beamlets),
    )


def canonical_table_of(beamlets: Sequence[Beamlet]) -> NDArray[np.float64]:
    """The numeric canonical table, float64 ``[J, 12]``, columns :data:`CANONICAL_TABLE_COLUMNS`:
    position, UNIT direction, energy, sigmas, ``n_histories`` (exact as float64, <= 2**32),
    ``has_weight`` (0 or 1) and ``weight`` (0.0 when absent). Labels are not numeric and enter
    only :func:`canonical_table_sha256`."""
    rows = [
        [
            *b.position_mm,
            *b.unit_direction(),
            b.kinetic_energy_mev,
            b.energy_sigma_mev,
            b.lateral_sigma_mm,
            float(b.n_histories),
            0.0 if b.weight is None else 1.0,
            0.0 if b.weight is None else float(b.weight),
        ]
        for b in beamlets
    ]
    return np.array(rows, dtype=np.float64).reshape(len(beamlets), len(CANONICAL_TABLE_COLUMNS))


def canonical_table_sha256(table: NDArray[np.float64], labels: Any) -> str:
    """The canonical digest of a spot list. Bytes hashed, in order: the ASCII text
    ``ionmc.beamlets/1`` and a NUL; J and the column count as little-endian uint64; the table as
    little-endian float64 in C order (``table.astype('<f8').tobytes()``); a NUL; the labels
    UTF-8 encoded and joined by LF, an absent label being the empty string."""
    t = np.ascontiguousarray(table, dtype="<f8")
    h = hashlib.sha256()
    h.update(b"ionmc.beamlets/1\0")
    h.update(np.array(t.shape, dtype="<u8").tobytes())
    h.update(t.tobytes())
    h.update(b"\0")
    h.update("\n".join("" if x is None else x for x in labels).encode("utf-8"))
    return h.hexdigest()
