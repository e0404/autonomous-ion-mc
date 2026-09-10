"""Python-facing API for tabulated stopping powers and CSDA ranges.

:class:`TabulatedStoppingPower` evaluates a prepared
:class:`~ionmc.data.stopping_tables.StoppingTable` through the shared-source
functions of :mod:`ionmc.physics.tabulated` on the ``python`` (float64
reference), ``numpy`` or ``warp`` (float32 kernels) execution path, with the
same method names and units as :class:`ionmc.stopping_power.AnalyticStoppingPower`
so that transport code does not care which layer supplies the quantity.

Tables come from the data layer (:func:`from_dataset`), which requires the
dataset to have been acquired into the local cache beforehand (offline use
phase, decision ``0007``).
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np

from ionmc.backend import mathlib, reference
from ionmc.data import cache as data_cache
from ionmc.data.cache import DatasetSpec
from ionmc.data.stopping_tables import StoppingTable, load_stopping_table
from ionmc.materials import Material
from ionmc.particles import Particle


class TabulatedStoppingPower:
    """Stopping power and CSDA range of one projectile in one material from a table.

    Parameters
    ----------
    table:
        Prepared table (energies in MeV, mass stopping power in MeV cm^2/g).
    material, particle:
        What the table describes; recorded in the provenance and used for
        consistency checks only (the numbers come from the table).
    path:
        ``"python"``, ``"numpy"`` or ``"warp"``.
    device:
        Warp device for ``path="warp"``.
    allow_extrapolation:
        If false (default), energies outside ``[E_min, E_max]`` of the table
        raise ``ValueError``; if true, the boundary segments' power laws are
        extrapolated.
    """

    def __init__(
        self,
        table: StoppingTable,
        material: Material,
        particle: Particle,
        path: str = "numpy",
        device: str = "cpu",
        allow_extrapolation: bool = False,
    ) -> None:
        if path not in ("python", "numpy", "warp"):
            raise ValueError(f"unknown execution path {path!r}")
        self.table = table
        self.material = material
        self.particle = particle
        self.path = path
        self.device = device
        self.allow_extrapolation = allow_extrapolation
        if path == "warp":
            from ionmc.backend import warp_kernels

            self._warp = warp_kernels.TabulatedKernels(table, device)
        else:
            self._module = reference.load_bound_module("ionmc.physics.tabulated", path)

    # -- construction helpers -------------------------------------------------

    @classmethod
    def from_dataset(
        cls,
        spec: DatasetSpec,
        material: Material,
        particle: Particle,
        cache_root: str | os.PathLike[str] | None = None,
        **kwargs: Any,
    ) -> TabulatedStoppingPower:
        """Load a registered dataset from the local cache (no network)."""
        path = data_cache.load_path(spec, cache_root)
        manifest = data_cache.read_manifest(spec, cache_root)
        provenance = {
            "dataset": spec.name,
            "version": spec.version,
            "sha256": spec.sha256,
            "license": spec.license,
            "source": dict(spec.provenance),
            "retrieved_at_utc": manifest.get("retrieved_at_utc"),
            "file": str(path),
        }
        table = load_stopping_table(path, provenance)
        return cls(table, material, particle, **kwargs)

    # -- evaluation -----------------------------------------------------------

    @property
    def energy_range_mev(self) -> tuple[float, float]:
        return float(self.table.energy_mev[0]), float(self.table.energy_mev[-1])

    def _check(self, energies: np.ndarray) -> None:
        if np.any(energies <= 0.0):
            raise ValueError("kinetic energies must be positive")
        if not self.allow_extrapolation:
            lo, hi = self.energy_range_mev
            if np.any(energies < lo) or np.any(energies > hi):
                raise ValueError(
                    f"energies outside the tabulated range [{lo}, {hi}] MeV; "
                    "pass allow_extrapolation=True to extrapolate"
                )

    def _table_args(self) -> tuple[Any, ...]:
        t = self.table
        return (
            t.energy_mev,
            t.stopping_mev_cm2_per_g,
            t.slope,
            t.size,
            t.bisection_steps,
        )

    def mass_stopping_power(self, kinetic_energy_mev: Any) -> np.ndarray:
        """Mass stopping power [MeV cm^2/g] at the given energies [MeV]."""
        energies = np.atleast_1d(np.asarray(kinetic_energy_mev, dtype=np.float64))
        self._check(energies)
        if self.path == "warp":
            return self._warp.mass_stopping_power(energies)
        if self.path == "numpy":
            values = self._module.tabulated_mass_stopping_power(
                energies, *self._table_args()
            )
            return np.asarray(values, dtype=np.float64)
        args = self._table_args()
        return np.array(
            [
                self._module.tabulated_mass_stopping_power(float(e), *args)
                for e in energies
            ],
            dtype=np.float64,
        )

    def csda_range(self, kinetic_energy_mev: Any) -> np.ndarray:
        """CSDA range [g/cm^2] from the table's first energy up to the given energies.

        The range below ``table.range_floor_energy_mev`` is not included
        (decision ``0008`` quantifies the neglected residual).
        """
        energies = np.atleast_1d(np.asarray(kinetic_energy_mev, dtype=np.float64))
        self._check(energies)
        t = self.table
        if self.path == "warp":
            return self._warp.csda_range(energies)
        args = (
            t.energy_mev,
            t.stopping_mev_cm2_per_g,
            t.slope,
            t.csda_range_g_per_cm2,
            t.size,
            t.bisection_steps,
        )
        if self.path == "numpy":
            values = self._module.tabulated_csda_range(energies, *args)
            return np.asarray(values, dtype=np.float64)
        return np.array(
            [self._module.tabulated_csda_range(float(e), *args) for e in energies],
            dtype=np.float64,
        )

    # -- provenance -----------------------------------------------------------

    def provenance(self) -> dict[str, Any]:
        """Metadata identifying exactly which table and evaluation path were used."""
        return {
            "model": "tabulated-pchip",
            "module": "ionmc.physics.tabulated",
            "path": self.path,
            "device": self.device if self.path == "warp" else None,
            "binding_default": mathlib.default_binding(),
            "material": self.material.name,
            "particle": self.particle.name,
            "table": dict(self.table.provenance),
            "table_size": self.table.size,
            "energy_range_mev": list(self.energy_range_mev),
            "range_floor_energy_mev": self.table.range_floor_energy_mev,
            "allow_extrapolation": self.allow_extrapolation,
        }
