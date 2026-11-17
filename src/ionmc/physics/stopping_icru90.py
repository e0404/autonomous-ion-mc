"""The registered model ``icru90-water-v1`` (decision 0044, section 1; V3-013 step S4).

:class:`Icru90WaterSource` builds the liquid-water stopping tables of protons and alpha particles
from the ICRU 90 electronic stopping power arrays that Geant4 embeds in ``G4ICRU90StoppingData.cc``
(parsed by :mod:`ionmc.data.icru90`). One generic
:class:`~ionmc.physics.stopping_tabulated.TabulatedStoppingSource` per projectile does the table
construction (log-log interpolation of S onto 200 points per decade from 1 MeV/u, exact log-log
range quadrature, ``model_id`` and ``version`` in every table's metadata); this module adds what
is specific to the model: the per-projectile windows (proton 1 to 500 MeV, alpha 1 to 250 MeV/u,
the alpha array ending at 1000 MeV total energy), the rule that the requested water I-value must
be 78 eV (row DOM-4), the start range and the identity.

Start range ``icru90-integrated-from-lowest-node-v1``: the range at the first grid node
(1 MeV/u) is ``R(1 MeV/u) = a E_0 / S_0 + int_{E_0}^{1 MeV/u} dE_total / S``, with ``E_0`` the
lowest node of the array of the projectile (total energy ``a E_0`` per u times the mass number
``a``), ``S_0`` the stopping power there, and the integral the exact integral of the log-log
interpolant of the array's own nodes (the same closed form as every other range of the model,
:func:`~ionmc.physics.stopping.exact_loglog_range_increments`). Only below the lowest node the
range is the constant-S approximation ``a E_0 / S_0`` (the range of the Bethe tables below their
grid, ``_BETHE_RANGE_START``); no NIST range enters. Every later range is the start range plus the
exact increments of the grid, so the node arrays ``R(E_i) - R(E_0)`` do not depend on the start
range at all (row CC-I2).

Canonical serialisation of the parsed arrays (:attr:`Icru90WaterSource.content_sha256`): SHA-256
over, for the projectiles in the order ``alpha``, ``proton`` (sorted names), the name (ASCII, then
a zero byte), the energy per nucleon [MeV/u] and the electronic mass stopping power
[MeV cm2/g], each as C-contiguous little-endian float64 bytes. The alpha energies are the array's
total energies divided by 4.

:func:`load_icru90_model` loads the cached object through the cache and registry API: a missing
object, a hash that differs from the registry pin or a parse failure raises (row DOM-6); nothing
is downloaded. No array value is stored in Git.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.data import cache as _cache
from ionmc.data.icru90 import Icru90Water, load_icru90_water
from ionmc.data.registry import DATASETS
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import Material
from ionmc.physics.projectiles import Projectile
from ionmc.physics.stopping import (
    StoppingTable,
    _require_nist_water,
    exact_loglog_range_increments,
)
from ionmc.physics.stopping_models import (
    ICRU_DATASET_ID,
    ICRU_I_VALUE_EV,
    ICRU_INTERPOLATION,
    ICRU_MODEL_ID,
    ICRU_RANGE_START,
    ICRU_SPECIES,
    MATERIAL_RULE_WATER,
    MATERIAL_RULE_WATER_I78,
    STOPPING_MODELS,
    StoppingDomain,
    StoppingModelSpec,
    icru_options,
    icru_pins,
)
from ionmc.physics.stopping_tabulated import TabulatedStoppingSource, _digest

_PROJECTILE_A = {"proton": 1.0, "alpha": 4.0}
_START_PER_U = 1.0  # MeV/u, the first node of every table


def icru90_start_range(
    energy_per_u: NDArray[np.float64],
    s_el_mass: NDArray[np.float64],
    a: float,
    e_start_per_u: float = _START_PER_U,
) -> float:
    """Range [g/cm2] at ``e_start_per_u`` MeV/u from the array's own nodes (see the module doc).

    ``a E_0 / S_0`` below the lowest node plus the exact log-log integral from the lowest node to
    ``e_start_per_u``; raises ``ValueError`` if the array starts above ``e_start_per_u``.
    """
    e = np.asarray(energy_per_u, dtype=np.float64)
    s = np.asarray(s_el_mass, dtype=np.float64)
    if e[0] > e_start_per_u * (1.0 + 1e-12):
        raise ValueError(
            f"the array starts at {e[0]} MeV/u, above the start energy {e_start_per_u} MeV/u"
        )
    r_below = a * float(e[0]) / float(s[0])
    k = int(np.searchsorted(e, e_start_per_u * (1.0 + 1e-12), side="left"))  # nodes below start
    ln_e = np.log(e[:k])
    ln_s = np.log(s[:k])
    if k == 0 or e[k - 1] < e_start_per_u * (1.0 - 1e-12):
        ln_e = np.append(ln_e, math.log(e_start_per_u))
        ln_s = np.append(ln_s, float(np.interp(math.log(e_start_per_u), np.log(e), np.log(s))))
    if ln_e.size < 2:
        return r_below
    return r_below + float(np.sum(exact_loglog_range_increments(ln_e, ln_s, a)))


@dataclass(frozen=True, eq=False)
class Icru90WaterSource:
    """ICRU 90 liquid-water tables of protons and alpha particles as one stopping model."""

    proton: TabulatedStoppingSource
    alpha: TabulatedStoppingSource
    name: str = "icru90"

    def __post_init__(self) -> None:
        for sp, src in (("proton", self.proton), ("alpha", self.alpha)):
            if not isinstance(src, TabulatedStoppingSource) or set(src._arrays) != {sp}:
                raise ValueError(f"Icru90WaterSource needs a TabulatedStoppingSource of {sp} only")

    def _source(self, projectile: Projectile) -> TabulatedStoppingSource:
        if projectile.name == "proton":
            return self.proton
        if projectile.name == "alpha":
            return self.alpha
        raise UnsupportedCombinationError(
            f"stopping model {ICRU_MODEL_ID!r} has no table for projectile {projectile.name!r} "
            "(species proton and alpha only; no fallback model)"
        )

    def _refuse(
        self, projectile: Projectile, material: Material, what: str
    ) -> UnsupportedCombinationError:
        """The fail-closed error naming model, species, material, energy and domain (DOM-1/4)."""
        from ionmc.physics.stopping_selection import _domain_text  # noqa: PLC0415

        return UnsupportedCombinationError(
            f"stopping model {ICRU_MODEL_ID!r} does not support species {projectile.name!r}, "
            f"material {material.name!r}{what}, energy any; its domain: "
            f"{_domain_text(STOPPING_MODELS[ICRU_MODEL_ID])} (no extrapolation, no fallback "
            "model, no automatic composite)"
        )

    def table(self, material: Material, projectile: Projectile) -> StoppingTable:
        """The table of ``projectile`` in ``material`` (row DOM-1, DOM-4; water only)."""
        if projectile.name not in ICRU_SPECIES:
            raise self._refuse(projectile, material, "")
        src = self._source(projectile)
        try:
            _require_nist_water(material)  # structural liquid water (row DOM-2)
        except ValueError as exc:
            raise self._refuse(projectile, material, "") from exc
        if not math.isclose(
            material.mean_excitation_eV, ICRU_I_VALUE_EV, rel_tol=0.0, abs_tol=1e-9
        ):
            raise self._refuse(
                projectile,
                material,
                f" with requested I = {material.mean_excitation_eV} eV (the tabulated data are "
                f"for I = {ICRU_I_VALUE_EV} eV only)",
            )
        return src.table(material, projectile)

    # -- identity --

    def content_sha256_for(self, projectile: Projectile) -> str:
        """Digest of the parsed arrays of ``projectile`` (see the module doc)."""
        return self._source(projectile).content_sha256

    @property
    def content_sha256(self) -> str:
        """Digest of the parsed arrays of both projectiles."""
        return _digest({**self.alpha._arrays, **self.proton._arrays})

    def datasets_used(self, projectile: Projectile) -> tuple[str, ...]:
        """The dataset id read for the table of ``projectile`` (the one ICRU 90 file)."""
        sid = self._source(projectile).metadata.get("dataset_id")
        return () if sid is None else (str(sid),)

    def window_per_u(self, projectile: Projectile) -> tuple[float, float]:
        """The MeV/u window of the grid the table of ``projectile`` is built on."""
        src = self._source(projectile)
        return float(src.e_min_per_u), float(src.e_max_per_u)

    def instance_spec(self) -> StoppingModelSpec:
        """The spec built from the fields of this instance, for the structural match of
        ``spec_of``: the registered pin only for arrays loaded from the verified cache."""
        meta = self.proton.metadata
        pins = (
            icru_pins()
            if meta.get("dataset_id") == ICRU_DATASET_ID
            and meta.get("dataset_sha256") == DATASETS[ICRU_DATASET_ID].sha256
            and self.alpha.metadata.get("dataset_sha256") == meta.get("dataset_sha256")
            else ()
        )
        species = frozenset(ICRU_SPECIES)
        lo = {"proton": float(self.proton.e_min_per_u), "alpha": float(self.alpha.e_min_per_u)}
        hi = {"proton": float(self.proton.e_max_per_u), "alpha": float(self.alpha.e_max_per_u)}
        ppd = self.proton.points_per_decade
        options = icru_options(str(self.name), ppd)
        if self.alpha.points_per_decade != ppd:
            options["points_per_decade_alpha"] = int(self.alpha.points_per_decade)
        return StoppingModelSpec(
            ICRU_MODEL_ID,
            "1",
            "tabulated",
            options,
            ICRU_INTERPOLATION,
            ICRU_RANGE_START,
            pins,
            StoppingDomain(species, MATERIAL_RULE_WATER_I78, lo, hi),
            "0044",
        )


def build_icru90_source(
    data: Icru90Water, *, name: str = "icru90", points_per_decade: int = 200
) -> Icru90WaterSource:
    """The model source from parsed ICRU 90 arrays (energies in total MeV, converted to MeV/u).

    The windows are the registered ones, capped at the array end of each projectile. The table
    metadata carry ``model_id``, ``version``, ``decision`` and the provenance fields of ``data``.
    """
    spec = STOPPING_MODELS[ICRU_MODEL_ID]
    meta: dict[str, Any] = {
        "model_id": ICRU_MODEL_ID,
        "version": spec.version,
        "decision": spec.decision,
        "dataset_id": data.dataset_id,
        "dataset_version": data.version,
        "dataset_sha256": data.sha256,
        "retrieved_at": data.retrieved_at,
        "file_content_sha256": data.content_sha256,
        "I_eV_model": ICRU_I_VALUE_EV,
    }
    sources: dict[str, TabulatedStoppingSource] = {}
    for sp, e_tot, s in (
        ("proton", data.proton_energy_mev, data.proton_stopping),
        ("alpha", data.alpha_energy_mev, data.alpha_stopping),
    ):
        a = _PROJECTILE_A[sp]
        e_u = np.asarray(e_tot, dtype=np.float64) / a
        s_arr = np.asarray(s, dtype=np.float64)
        lo = spec.domain.e_min_per_u[sp]
        hi = min(spec.domain.e_max_per_u[sp], float(e_u[-1]))

        def start(
            _p: Projectile,
            _m: Material,
            _e: NDArray[np.float64],
            _s: NDArray[np.float64],
            _e_u: NDArray[np.float64] = e_u,
            _s_arr: NDArray[np.float64] = s_arr,
            _a: float = a,
        ) -> float:
            return icru90_start_range(_e_u, _s_arr, _a)

        sources[sp] = TabulatedStoppingSource(
            name,
            meta,
            {sp: (e_u, s_arr)},
            lo,
            hi,
            ICRU_RANGE_START,
            start,
            MATERIAL_RULE_WATER,
            points_per_decade,
        )
    return Icru90WaterSource(sources["proton"], sources["alpha"], name)


def load_icru90_model(cache_dir: str | None = None) -> Icru90WaterSource:
    """Load the ICRU 90 object from the verified cache and return the model (row DOM-6).

    Raises ``FileNotFoundError`` (object or manifest absent), ``IntegrityError`` (bytes differ
    from the manifest or the registry pin) or ``ValueError`` (parse failure); never downloads and
    never falls back to anything else.
    """
    cdir = _cache.resolve_cache_dir(cache_dir)
    path = _cache.verify(ICRU_DATASET_ID, cdir)
    data = load_icru90_water(path)  # re-hashes the bytes against the registry pin
    if data.dataset_id != ICRU_DATASET_ID or data.sha256 != DATASETS[ICRU_DATASET_ID].sha256:
        raise _cache.IntegrityError(f"{ICRU_DATASET_ID}: loaded object is not the pinned dataset")
    return build_icru90_source(data)


__all__ = [
    "Icru90WaterSource",
    "build_icru90_source",
    "icru90_start_range",
    "load_icru90_model",
]
