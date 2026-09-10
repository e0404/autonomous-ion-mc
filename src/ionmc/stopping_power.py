"""Python-facing API for the analytical electronic stopping power.

:class:`AnalyticStoppingPower` binds a :class:`~ionmc.materials.Material` and a
:class:`~ionmc.particles.Particle` to the shared-source physics of
:mod:`ionmc.physics.stopping` and evaluates it on the requested execution
path:

* ``"python"``  - the float64 reference path (pure Python, scalar loop);
* ``"numpy"``   - the same source vectorised with numpy (fast reference);
* ``"warp"``    - Warp kernels on a CPU or CUDA device (float32).

The material-dependent constants that the shared functions take as scalars are
collected once in :class:`StoppingPowerParameters`, which is also what the
Warp kernels receive. Every result can be accompanied by :meth:`provenance`,
which records the mean excitation energy, its source, the density-effect
parameter set and the enabled corrections, as required for reproducibility.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from ionmc.backend import mathlib, reference
from ionmc.materials import ELEMENTS, Material
from ionmc.particles import Particle
from ionmc.physics.barkas_table import BARKAS_F, BARKAS_W, barkas_b_parameter
from ionmc.units import ev_to_mev

#: Default number of Simpson sub-intervals per CSDA-range integration.
DEFAULT_RANGE_STEPS: int = 200

#: Lower limit of the analytic CSDA-range integration, MeV. Below this the
#: analytic model is not trusted; the residual range of a 1 MeV proton in
#: water is about 0.0025 g/cm^2 and is neglected (documented in decision 0006).
RANGE_ENERGY_FLOOR_MEV: float = 1.0


@dataclass(frozen=True)
class Corrections:
    """Which Bethe corrections are enabled."""

    density_effect: bool = True
    shell: bool = True
    barkas: bool = True
    bloch: bool = True


#: Default correction set (all corrections enabled).
DEFAULT_CORRECTIONS: Corrections = Corrections()


@dataclass(frozen=True)
class StoppingPowerParameters:
    """Scalar/table inputs of the shared stopping-power functions."""

    rest_energy_mev: float
    charge: float
    za_ratio: float
    i_mev: float
    cbar: float
    x0: float
    x1: float
    a: float
    mexp: float
    delta0: float
    shell_m2: float
    shell_m3: float
    elem_z: tuple[float, ...]
    elem_b: tuple[float, ...]
    elem_f: tuple[float, ...]
    table_w: tuple[float, ...]
    table_f: tuple[float, ...]
    use_density: float
    use_shell: float
    use_barkas: float
    use_bloch: float

    @property
    def n_elem(self) -> int:
        return len(self.elem_z)

    @property
    def n_table(self) -> int:
        return len(self.table_w)

    def as_call_args(self, arrays: Any = None) -> tuple[Any, ...]:
        """Arguments after ``kinetic_energy`` for ``mass_stopping_power``.

        ``arrays`` optionally maps the five table names to array objects to
        substitute for the tuples (used for the numpy and Warp paths).
        """
        tabs = {
            "elem_z": self.elem_z,
            "elem_b": self.elem_b,
            "elem_f": self.elem_f,
            "table_w": self.table_w,
            "table_f": self.table_f,
        }
        if arrays:
            tabs.update(arrays)
        return (
            self.rest_energy_mev,
            self.charge,
            self.za_ratio,
            self.i_mev,
            self.cbar,
            self.x0,
            self.x1,
            self.a,
            self.mexp,
            self.delta0,
            self.shell_m2,
            self.shell_m3,
            tabs["elem_z"],
            tabs["elem_b"],
            tabs["elem_f"],
            self.n_elem,
            tabs["table_w"],
            tabs["table_f"],
            self.n_table,
            self.use_density,
            self.use_shell,
            self.use_barkas,
            self.use_bloch,
        )


def build_parameters(
    material: Material,
    particle: Particle,
    corrections: Corrections = DEFAULT_CORRECTIONS,
) -> StoppingPowerParameters:
    """Derive the shared-function parameters for a material/projectile pair."""
    za = material.electrons_per_gram_ratio
    shell_m2 = 0.0
    shell_m3 = 0.0
    elem_z = []
    elem_b = []
    elem_f = []
    for symbol in material.mass_fractions:
        elem = ELEMENTS[symbol]
        atoms_per_electron = material.atoms_per_gram(symbol) / za
        i_ev = elem.mean_excitation_energy_ev
        shell_m2 += atoms_per_electron * i_ev * i_ev
        shell_m3 += atoms_per_electron * i_ev * i_ev * i_ev
        elem_z.append(float(elem.atomic_number))
        elem_b.append(barkas_b_parameter(elem.atomic_number))
        elem_f.append(material.electron_fraction(symbol))
    de = material.density_effect
    if de is None:
        cbar, x0, x1, a, mexp, delta0 = 0.0, 0.0, 0.0, 0.0, 1.0, 0.0
        use_density = 0.0
    else:
        cbar, x0, x1, a, mexp, delta0 = de.cbar, de.x0, de.x1, de.a, de.m, de.delta0
        use_density = 1.0 if corrections.density_effect else 0.0
    return StoppingPowerParameters(
        rest_energy_mev=particle.rest_energy_mev,
        charge=particle.charge,
        za_ratio=za,
        i_mev=ev_to_mev(material.mean_excitation_energy.value_ev),
        cbar=cbar,
        x0=x0,
        x1=x1,
        a=a,
        mexp=mexp,
        delta0=delta0,
        shell_m2=shell_m2,
        shell_m3=shell_m3,
        elem_z=tuple(elem_z),
        elem_b=tuple(elem_b),
        elem_f=tuple(elem_f),
        table_w=BARKAS_W,
        table_f=BARKAS_F,
        use_density=use_density,
        use_shell=1.0 if corrections.shell else 0.0,
        use_barkas=1.0 if corrections.barkas else 0.0,
        use_bloch=1.0 if corrections.bloch else 0.0,
    )


class AnalyticStoppingPower:
    """Analytical electronic stopping power for one material and projectile.

    Parameters
    ----------
    material, particle:
        Target and projectile.
    corrections:
        Enabled Bethe corrections.
    path:
        ``"python"`` (float64 reference), ``"numpy"`` (vectorised reference) or
        ``"warp"`` (float32 Warp kernels; see ``device``).
    device:
        Warp device alias for ``path="warp"`` (``"cpu"`` or ``"cuda:0"``).
    """

    def __init__(
        self,
        material: Material,
        particle: Particle,
        corrections: Corrections = DEFAULT_CORRECTIONS,
        path: str = "numpy",
        device: str = "cpu",
    ) -> None:
        if path not in ("python", "numpy", "warp"):
            raise ValueError(f"unknown execution path {path!r}")
        self.material = material
        self.particle = particle
        self.corrections = corrections
        self.path = path
        self.device = device
        self.parameters = build_parameters(material, particle, corrections)
        if path == "warp":
            from ionmc.backend import warp_kernels

            self._warp = warp_kernels.StoppingPowerKernels(self.parameters, device)
        else:
            self._module = reference.load_bound_module("ionmc.physics.stopping", path)

    # -- evaluation ---------------------------------------------------------

    def mass_stopping_power(self, kinetic_energy_mev: Any) -> np.ndarray:
        """Electronic mass stopping power [MeV cm^2/g] at the given energies [MeV]."""
        energies = np.atleast_1d(np.asarray(kinetic_energy_mev, dtype=np.float64))
        if self.path == "warp":
            return self._warp.mass_stopping_power(energies)
        if self.path == "numpy":
            args = self.parameters.as_call_args(_numpy_tables(self.parameters))
            values = self._module.mass_stopping_power(energies, *args)
            return np.asarray(values, dtype=np.float64)
        args = self.parameters.as_call_args()
        return np.array(
            [self._module.mass_stopping_power(float(e), *args) for e in energies],
            dtype=np.float64,
        )

    def csda_range(
        self,
        kinetic_energy_mev: Any,
        energy_floor_mev: float = RANGE_ENERGY_FLOOR_MEV,
        n_steps: int = DEFAULT_RANGE_STEPS,
    ) -> np.ndarray:
        """CSDA range [g/cm^2] from ``energy_floor_mev`` up to the given energies.

        The contribution below the floor is neglected (about 0.0025 g/cm^2 for
        protons in water at the default floor of 1 MeV).
        """
        if n_steps % 2:
            raise ValueError("n_steps must be even for Simpson integration")
        energies = np.atleast_1d(np.asarray(kinetic_energy_mev, dtype=np.float64))
        if np.any(energies <= energy_floor_mev):
            raise ValueError("energies must exceed the integration floor")
        if self.path == "warp":
            return self._warp.csda_range(energies, energy_floor_mev, n_steps)
        if self.path == "numpy":
            args = self.parameters.as_call_args(_numpy_tables(self.parameters))
            values = self._module.csda_range_increment(
                energy_floor_mev, energies, n_steps, *args
            )
            return np.asarray(values, dtype=np.float64)
        args = self.parameters.as_call_args()
        return np.array(
            [
                self._module.csda_range_increment(
                    energy_floor_mev, float(e), n_steps, *args
                )
                for e in energies
            ],
            dtype=np.float64,
        )

    # -- provenance ---------------------------------------------------------

    def provenance(self) -> dict[str, Any]:
        """Metadata identifying exactly which model and constants were used."""
        de = self.material.density_effect
        return {
            "model": "bethe-analytic",
            "module": "ionmc.physics.stopping",
            "path": self.path,
            "device": self.device if self.path == "warp" else None,
            "binding_default": mathlib.default_binding(),
            "material": self.material.name,
            "material_source": self.material.source,
            "mean_excitation_energy_ev": self.material.mean_excitation_energy.value_ev,
            "mean_excitation_energy_source": (
                self.material.mean_excitation_energy.source
            ),
            "density_effect": asdict(de) if de is not None else None,
            "corrections": asdict(self.corrections),
            "particle": self.particle.name,
            "range_energy_floor_mev": RANGE_ENERGY_FLOOR_MEV,
        }


def _numpy_tables(params: StoppingPowerParameters) -> dict[str, np.ndarray]:
    return {
        name: np.asarray(getattr(params, name), dtype=np.float64)
        for name in ("elem_z", "elem_b", "elem_f", "table_w", "table_f")
    }


def available_paths() -> Sequence[str]:
    """Execution paths usable in this environment."""
    paths = ["python", "numpy"]
    if mathlib.HAVE_WARP:
        paths.append("warp")
    return tuple(paths)
