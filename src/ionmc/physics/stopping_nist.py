"""The registered model ``nist-pstar-astar-water-2005`` (decision 0044, section 1; V3-013 step S3).

:class:`NistStarCombinedSource` joins the two pre-existing single-program sources
(:class:`~ionmc.physics.stopping.NistStarStoppingSource`, one PSTAR or ASTAR table each) into the
one model of the registry: ``table(material, projectile)`` delegates to the program of the
projectile, so every table, its identity and its hash are those of the old single-program source
(row ID-NIST); the model id is *not* put into that metadata (the model existed before the freeze,
decision 0044 section 4), the single-model assertion recognises the tables by their dataset pin
(:func:`~ionmc.physics.stopping_selection.table_model_id`). Any other projectile fails closed
(no fallback, no automatic composite).

:func:`load_nist_model` loads both cached objects through the cache and registry API: a missing
object, a hash that differs from the registry pin or a parse failure raises (row DOM-6); nothing
is downloaded. The upper window ends are those of the grids the tables are built on
(``min(e_max_per_u, table end / a)``) and are read from the cached tables at resolve time, never
written into Git.

Canonical serialisation of the parsed arrays (:meth:`NistStarCombinedSource.content_sha256_for`):
SHA-256 over, for the program name (ASCII, then a zero byte), the total kinetic energy
``energy_mev``, the electronic mass stopping power ``s_electronic`` and the CSDA range
``csda_range``, each as C-contiguous little-endian float64 bytes in this order; for the
whole source the programs in the order PSTAR, ASTAR.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

import numpy as np

from ionmc.data import cache as _cache
from ionmc.data.nist_star import StarTable, load_star_table
from ionmc.data.registry import DATASETS
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import Material
from ionmc.physics.projectiles import Projectile
from ionmc.physics.stopping import NistStarStoppingSource, StoppingTable
from ionmc.physics.stopping_models import (
    DATASET_OF,
    MATERIAL_RULE_WATER,
    NIST_CONTENT_SHA256,
    NIST_INTERPOLATION,
    NIST_MODEL_ID,
    NIST_RANGE_START,
    PROGRAM_OF,
    StoppingDomain,
    StoppingModelSpec,
    nist_options,
    nist_pins,
)


def _arrays_digest(tables: list[StarTable]) -> str:
    h = hashlib.sha256()
    for t in tables:
        h.update(t.program.encode("ascii") + b"\0")
        for a in (t.energy_mev, t.s_electronic, t.csda_range):
            h.update(np.ascontiguousarray(a, dtype="<f8").tobytes())
    return h.hexdigest()


@dataclass(frozen=True, eq=False)
class NistStarCombinedSource:
    """PSTAR (proton) and ASTAR (alpha) liquid-water tables as one stopping model."""

    pstar: NistStarStoppingSource
    astar: NistStarStoppingSource
    name: str = "nist-star"

    def __post_init__(self) -> None:
        for src, prog in ((self.pstar, "PSTAR"), (self.astar, "ASTAR")):
            if not isinstance(src, NistStarStoppingSource) or src.star_table.program != prog:
                raise ValueError(f"NistStarCombinedSource needs a {prog} NistStarStoppingSource")

    def _program_source(self, projectile: Projectile) -> NistStarStoppingSource:
        if projectile.name == "proton":
            return self.pstar
        if projectile.name == "alpha":
            return self.astar
        raise UnsupportedCombinationError(
            f"stopping model {NIST_MODEL_ID!r} has no table for projectile {projectile.name!r} "
            "(species proton and alpha only; no fallback model)"
        )

    def table(self, material: Material, projectile: Projectile) -> StoppingTable:
        """The table of the single-program source of ``projectile`` (metadata unchanged)."""
        return self._program_source(projectile).table(material, projectile)

    # -- identity --

    def content_sha256_for(self, projectile: Projectile) -> str:
        """Digest of the parsed arrays of the program of ``projectile`` (see the module doc)."""
        return _arrays_digest([self._program_source(projectile).star_table])

    @property
    def content_sha256(self) -> str:
        """Digest of the parsed arrays of both programs (PSTAR, then ASTAR)."""
        return _arrays_digest([self.pstar.star_table, self.astar.star_table])

    def datasets_used(self, projectile: Projectile) -> tuple[str, ...]:
        """The dataset id read for the tables of ``projectile`` (PSTAR for protons only)."""
        sid = self._program_source(projectile).star_table.dataset_id
        return () if sid is None else (sid,)

    def window_per_u(self, projectile: Projectile) -> tuple[float, float]:
        """The MeV/u window of the grid the table of ``projectile`` is built on (resolved from the
        cached table: the grid cap, or the table end if that is lower)."""
        src = self._program_source(projectile)
        end = float(src.star_table.energy_mev[-1]) / projectile.a
        return float(src.e_min_per_u), min(float(src.e_max_per_u), end)

    def is_exact_implementation(self) -> bool:
        """True only for the exact classes (a subclass may override ``table``). The start range of
        this model is read from the table's ``csda_range``, which is part of the content digest,
        so no separate start-range digest is needed."""
        return (
            type(self) is NistStarCombinedSource
            and type(self.pstar) is NistStarStoppingSource
            and type(self.astar) is NistStarStoppingSource
        )

    def instance_spec(self) -> StoppingModelSpec:
        """The spec built from the fields of this instance, for the structural match of
        ``spec_of``: registered pins only for tables loaded from the verified cache."""
        pins = []
        live = {  # recomputed from the live arrays: the identity never trusts the metadata alone
            "proton": _arrays_digest([self.pstar.star_table]),
            "alpha": _arrays_digest([self.astar.star_table]),
        }
        for sp, src in (("proton", self.pstar), ("alpha", self.astar)):
            ds = src.star_table.dataset_id
            if (
                ds == DATASET_OF[sp]
                and src.star_table.sha256 == DATASETS[ds].sha256
                and live[sp] == NIST_CONTENT_SHA256[sp]
            ):
                pins.append(nist_pins()[0 if sp == "proton" else 1])
        species = frozenset(PROGRAM_OF)
        lo = {"proton": float(self.pstar.e_min_per_u), "alpha": float(self.astar.e_min_per_u)}
        hi = {"proton": float(self.pstar.e_max_per_u), "alpha": float(self.astar.e_max_per_u)}
        ppd = self.pstar.points_per_decade
        options = nist_options(str(self.name), ppd, live)
        if self.astar.points_per_decade != ppd:
            options["points_per_decade_astar"] = int(self.astar.points_per_decade)
        return StoppingModelSpec(
            NIST_MODEL_ID,
            "1",
            "tabulated",
            options,
            NIST_INTERPOLATION,
            NIST_RANGE_START,
            tuple(pins),
            StoppingDomain(species, MATERIAL_RULE_WATER, lo, hi),
            "0044",
        )


def load_nist_model(cache_dir: str | None = None) -> NistStarCombinedSource:
    """Load PSTAR and ASTAR from the verified cache and return the combined source (row DOM-6).

    Raises ``FileNotFoundError`` (object or manifest absent), ``IntegrityError`` (bytes differ
    from the manifest or the registry pin) or ``ValueError`` (parse failure, wrong program); never
    downloads and never falls back to anything else.
    """
    cdir = _cache.resolve_cache_dir(cache_dir)
    tables: dict[str, StarTable] = {}
    for sp, dataset_id in DATASET_OF.items():
        path = _cache.verify(dataset_id, cdir)
        table = load_star_table(path)  # re-hashes the bytes against the registry pin
        if table.dataset_id != dataset_id or table.sha256 != DATASETS[dataset_id].sha256:
            raise _cache.IntegrityError(f"{dataset_id}: loaded table is not the pinned dataset")
        if table.program != PROGRAM_OF[sp]:
            raise ValueError(f"{dataset_id}: program {table.program!r}, expected {PROGRAM_OF[sp]}")
        tables[sp] = table
    return NistStarCombinedSource(
        NistStarStoppingSource(tables["proton"]), NistStarStoppingSource(tables["alpha"])
    )
