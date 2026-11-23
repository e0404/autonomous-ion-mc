"""Generic tabulated stopping source (decision 0044, sections 1 and 4; V3-013 step S2b).

:class:`TabulatedStoppingSource` turns per-projectile arrays of the mass electronic stopping power
into :class:`~ionmc.physics.stopping.StoppingTable` objects. It knows nothing about where the
arrays come from (no data file, no parser): a model that reads data (S3, S4) or a test model builds
the arrays and hands them over together with its table metadata, its start-range rule and its
material rule.

Construction of a table (the one rule of every tabulated model):

* the energy grid is ``log_grid(e_min_per_u, e_max_per_u, points_per_decade)``, both ends inside the
  span of the array of the projectile (no extrapolation: a request outside the span, or for a
  projectile without an array, fails closed with ``ValueError``);
* ``S`` on the grid is the log-log interpolant of the array (piecewise linear in ``ln E`` and
  ``ln S``), the same interpolation the table itself uses between its nodes;
* the range at the first node is the ``start_range`` rule; the range at the other nodes is the exact
  integral of the log-log interpolated ``S`` (:func:`~ionmc.physics.stopping.build_table`,
  ``range_construction = "exact-loglog-quadrature-v1"``).

Every table carries ``model_id`` and ``version`` in its metadata (decision 0044 section 4), so the
single-model assertion of ``validate`` can name the model a table identity belongs to, and the
metadata (hence ``TransportTables.sha256``) of two models never coincide.
"""

from __future__ import annotations

import hashlib
import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ionmc.materials import Material
from ionmc.physics.projectiles import Projectile
from ionmc.physics.stopping import (
    RANGE_CONSTRUCTION,
    StoppingTable,
    _require_nist_water,
    build_table,
    log_grid,
)
from ionmc.physics.stopping_models import (
    MATERIAL_RULE_ANY,
    MATERIAL_RULE_WATER,
    StoppingDomain,
    StoppingModelSpec,
    snapshot,
)

# start-range rule: (projectile, material, energy grid [MeV/u], S on the grid) -> range [g/cm2]
StartRange = Callable[[Projectile, Material, NDArray[np.float64], NDArray[np.float64]], float]

_TABULATED_INTERPOLATION = f"log-log interpolation of S onto the grid; range {RANGE_CONSTRUCTION}"
_REQUIRED_METADATA = ("model_id", "version")
_SPEC_ONLY_KEYS = ("model_id", "version", "decision")
_SPAN_REL = 1.0e-12


def _digest(arrays: Mapping[str, tuple[NDArray[np.float64], NDArray[np.float64]]]) -> str:
    h = hashlib.sha256()
    for name in sorted(arrays):
        e, s = arrays[name]
        h.update(name.encode("ascii") + b"\0")
        h.update(np.ascontiguousarray(e, dtype="<f8").tobytes())
        h.update(np.ascontiguousarray(s, dtype="<f8").tobytes())
    return h.hexdigest()


@dataclass(frozen=True, slots=True, eq=False)
class TabulatedStoppingSource:
    """Stopping source from per-projectile arrays (see the module docstring).

    ``arrays`` maps a projectile name to ``(energy_per_u [MeV/u], s_el_mass [MeV cm2/g])``
    (strictly increasing positive energies, positive finite stopping powers). ``metadata`` must
    carry ``model_id`` and ``version`` (strings) and may carry ``decision`` (default ``"user"``)
    and any further JSON-ready entries; all of it is copied into every table's metadata.
    ``start_range`` is the range at the first grid node: a callable
    ``(projectile, material, energy_per_u, s_el_mass) -> g/cm2`` or a mapping projectile name to
    a value; ``start_range_id`` names the rule (it is part of the model's identity; a callable
    itself cannot be hashed). ``material_rule`` is the model's
    :attr:`~ionmc.physics.stopping_models.StoppingDomain.material_rule`: structural liquid water
    (checked by composition and density, not the name) or any material.
    """

    name: str
    metadata: Mapping[str, Any]
    arrays: Mapping[str, tuple[ArrayLike, ArrayLike]]
    e_min_per_u: float
    e_max_per_u: float
    start_range_id: str
    start_range: StartRange | Mapping[str, float]
    material_rule: str = MATERIAL_RULE_WATER
    points_per_decade: int = 200
    _arrays: dict[str, tuple[NDArray[np.float64], NDArray[np.float64]]] = field(
        init=False, repr=False
    )

    def __post_init__(self) -> None:
        for key in _REQUIRED_METADATA:
            if not isinstance(self.metadata.get(key), str) or not self.metadata[key]:
                raise ValueError(f"TabulatedStoppingSource metadata must carry a string {key!r}")
        if self.material_rule not in (MATERIAL_RULE_ANY, MATERIAL_RULE_WATER):
            raise ValueError(f"unknown material rule {self.material_rule!r}")
        if not callable(self.start_range) and not isinstance(self.start_range, Mapping):
            raise ValueError("start_range must be a callable or a mapping of projectile to value")
        if not 0.0 < self.e_min_per_u < self.e_max_per_u:
            raise ValueError("need 0 < e_min_per_u < e_max_per_u")
        if self.points_per_decade < 1:
            raise ValueError("points_per_decade must be positive")
        clean: dict[str, tuple[NDArray[np.float64], NDArray[np.float64]]] = {}
        for pname, (e_raw, s_raw) in self.arrays.items():
            e = np.array(e_raw, dtype=np.float64)
            s = np.array(s_raw, dtype=np.float64)
            if e.ndim != 1 or e.shape != s.shape or e.size < 2:
                raise ValueError(f"{pname}: energy and stopping arrays must be 1-D, equal, >= 2")
            if not (np.all(np.isfinite(e)) and np.all(np.isfinite(s))):
                raise ValueError(f"{pname}: arrays contain NaN or infinity")
            if e[0] <= 0.0 or np.any(np.diff(e) <= 0.0):
                raise ValueError(f"{pname}: energy_per_u must be positive and strictly increasing")
            if np.any(s <= 0.0):
                raise ValueError(f"{pname}: s_el_mass must be positive")
            if self.e_min_per_u < e[0] * (1 - _SPAN_REL) or self.e_max_per_u > e[-1] * (
                1 + _SPAN_REL
            ):
                raise ValueError(
                    f"{pname}: grid [{self.e_min_per_u}, {self.e_max_per_u}] MeV/u is outside the "
                    f"array span [{e[0]}, {e[-1]}] MeV/u (no extrapolation)"
                )
            e.setflags(write=False)
            s.setflags(write=False)
            clean[pname] = (e, s)
        if not clean:
            raise ValueError("at least one projectile array is required")
        object.__setattr__(self, "_arrays", clean)

    # -- identity --

    @property
    def content_sha256(self) -> str:
        """SHA-256 of the arrays (names, energies, stopping powers), sorted by projectile."""
        return _digest(self._arrays)

    def instance_spec(self) -> StoppingModelSpec:
        """The spec of this source built from its fields (``data = ()``: no pins at this level).

        The model id, version and decision are the source's own metadata entries; whether they are
        a registered identity is decided by ``spec_of`` (structural comparison, content included).
        """
        names = frozenset(self._arrays)
        extra = {k: v for k, v in self.metadata.items() if k not in _SPEC_ONLY_KEYS}
        return StoppingModelSpec(
            str(self.metadata["model_id"]),
            str(self.metadata["version"]),
            "tabulated",
            {
                "points_per_decade": int(self.points_per_decade),
                "source_name": str(self.name),
                "content_sha256": self.content_sha256,
                "metadata": extra,
            },
            _TABULATED_INTERPOLATION,
            self.start_range_id,
            (),
            StoppingDomain(
                names,
                self.material_rule,
                dict.fromkeys(names, float(self.e_min_per_u)),
                dict.fromkeys(names, float(self.e_max_per_u)),
            ),
            str(self.metadata.get("decision", "user")),
        )

    # -- tables --

    def _start(
        self,
        projectile: Projectile,
        material: Material,
        e: NDArray[np.float64],
        s: NDArray[np.float64],
    ) -> float:
        if isinstance(self.start_range, Mapping):
            if projectile.name not in self.start_range:
                raise ValueError(f"no start range for projectile {projectile.name!r}")
            return float(self.start_range[projectile.name])
        return float(self.start_range(projectile, material, e, s))

    def table(self, material: Material, projectile: Projectile) -> StoppingTable:
        """Build the table of ``projectile`` in ``material`` (fail closed outside the domain)."""
        if projectile.name not in self._arrays:
            raise ValueError(
                f"source {self.name!r} has no data for projectile {projectile.name!r} "
                f"(has {sorted(self._arrays)})"
            )
        if self.material_rule == MATERIAL_RULE_WATER:
            try:
                _require_nist_water(material)
            except ValueError as exc:
                raise ValueError(
                    f"source {self.name!r} describes liquid water only (composition and density); "
                    f"material {material.name!r} does not match"
                ) from exc
        # the digest is taken from the very copies consumed below, inside this call (no cached
        # digest): the recorded content identity is that of the data the table is built from
        consumed = {
            n: (np.array(ea, dtype=np.float64), np.array(sa, dtype=np.float64))
            for n, (ea, sa) in self._arrays.items()
        }
        content = _digest(consumed)
        e_a, s_a = consumed[projectile.name]
        e = log_grid(self.e_min_per_u, self.e_max_per_u, self.points_per_decade)
        s = np.exp(np.interp(np.log(e), np.log(e_a), np.log(s_a)))
        r0 = self._start(projectile, material, e, s)
        if not math.isfinite(r0) or r0 <= 0.0:
            raise ValueError(f"start range {r0!r} is not positive and finite")
        meta: dict[str, Any] = {
            "source": "tabulated",
            "source_name": self.name,
            "material": material.name,
            "material_source": material.source,
            "material_density_g_cm3": material.density_g_cm3,
            "material_mass_fractions": dict(material.mass_fractions),
            "projectile": projectile.name,
            "I_eV": material.mean_excitation_eV,
            "content_sha256": content,
            "start_range": self.start_range_id,
            **self.metadata,
        }
        return build_table(projectile, material, e, s, r0, meta)


_SNAP = snapshot(TabulatedStoppingSource)
