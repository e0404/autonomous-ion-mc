"""Append-only species registry for species-resolved scoring (decision 0040, section 6).

Species ids are stable integers: an id is never reused, reordered or renumbered. Ions occupy ids
0 to ``MAX_ION_SPECIES - 1`` (V3-009 appends further ions); pseudo-species (energy that no
transported particle carries, for example the local deposit of nuclear fragments, V3-005) start at
``PSEUDO_BASE``. Channels select species with an ``int8`` match matrix of ``N_SPECIES_SLOTS``
columns indexed by id (not with bitmasks), so appending a species never changes an existing
column.

Which species the engine can *produce* is an engine capability, not a property of the registry:
:func:`producible` returns the set of ``(species name, generation)`` pairs a given source can
create in the current engine. A request for anything else fails closed before transport.
"""

from __future__ import annotations

from dataclasses import dataclass

from ionmc.physics.projectiles import (
    ALPHA,
    DEUTERON,
    HELIUM3,
    PROTON,
    TRITON,
    Projectile,
)

MAX_ION_SPECIES = 64
PSEUDO_BASE = 64
N_SPECIES_SLOTS = 65  # ids 0..63 ions, 64 nuclear_local (room for more pseudo-species by growth)
GENERATIONS = ("all", "primary", "secondary")


@dataclass(frozen=True)
class Species:
    """A registry entry: ``id`` (stable), ``name`` (the projectile name or a pseudo-species
    name), ``transported`` (carries a track through the engine; False for pseudo-species)."""

    id: int
    name: str
    transported: bool
    projectile: Projectile | None = None


def _ion(species_id: int, projectile: Projectile) -> Species:
    return Species(species_id, projectile.name, True, projectile)


# Append-only. Never reorder or reuse an id; append new ions at the end of the ion block.
_ION_ENTRIES: tuple[Species, ...] = (
    _ion(0, PROTON),
    _ion(1, DEUTERON),
    _ion(2, TRITON),
    _ion(3, HELIUM3),
    _ion(4, ALPHA),
)
_PSEUDO_ENTRIES: tuple[Species, ...] = (Species(PSEUDO_BASE, "nuclear_local", False),)

REGISTRY: tuple[Species, ...] = _ION_ENTRIES + _PSEUDO_ENTRIES
NUCLEAR_LOCAL = "nuclear_local"
_BY_NAME = {s.name: s for s in REGISTRY}
_BY_ID = {s.id: s for s in REGISTRY}

if len(_BY_NAME) != len(REGISTRY) or len(_BY_ID) != len(REGISTRY):  # pragma: no cover
    raise RuntimeError("species registry has duplicate names or ids")
if any(s.id >= N_SPECIES_SLOTS for s in REGISTRY):  # pragma: no cover
    raise RuntimeError("species id exceeds N_SPECIES_SLOTS")


def species_by_name(name: str) -> Species:
    """The registry entry called ``name``; ``ValueError`` for an unknown species."""
    try:
        return _BY_NAME[name]
    except (KeyError, TypeError):
        raise ValueError(
            f"unknown species {name!r}; the registry knows {[s.name for s in REGISTRY]}"
        ) from None


def species_by_id(species_id: int) -> Species:
    """The registry entry with id ``species_id``."""
    try:
        return _BY_ID[species_id]
    except KeyError:
        raise ValueError(f"unknown species id {species_id!r}") from None


def species_of_projectile(projectile: Projectile) -> Species:
    """The registry entry of a projectile (identified by its name and exact fields)."""
    s = _BY_NAME.get(projectile.name)
    if s is None or s.projectile != projectile:
        raise ValueError(f"projectile {projectile.name!r} is not in the species registry")
    return s


def transported_ids() -> tuple[int, ...]:
    """Ids of every transported (charged-ion) species."""
    return tuple(s.id for s in REGISTRY if s.transported)


def producible(projectile: Projectile) -> frozenset[tuple[str, str]]:
    """Engine capability: the ``(species name, generation)`` pairs the engine can create from a
    source of ``projectile``. V3-004 (no secondary transport): the source projectile as a
    primary. V3-005A extends this set; nothing else may be requested until then."""
    return frozenset({(species_of_projectile(projectile).name, "primary")})
