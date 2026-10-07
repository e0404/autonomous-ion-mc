"""Transport tables: log-log stopping, range and inverse-range tables per material.

Built (float64) from :class:`ionmc.physics.stopping.StoppingTable` objects:

* ``ln S(E)`` on a uniform grid in ``ln E`` (at least 200 points per decade), log-log interpolated
  from the source table;
* the CSDA range as the exact integral of that interpolated ``S``
  (``range_construction = exact-loglog-quadrature-v1``, decision 0039, V3-003D). On interval ``i``
  ``f = E / S`` is ``f_i exp(d_i (u - u_i) / h)`` (``u = ln E``, ``d_i = ln f_{i+1} - ln f_i``,
  ``h`` the grid step), so the node ranges are ``R_{i+1} = R_i + h f_i (e^{d_i} - 1) / d_i``
  starting from the source table's range at the first grid energy, and the range at the fraction
  ``phi`` of a bin is ``R_i + h f_i (e^{d_i phi} - 1) / d_i`` (shared function ``range_in_bin``).
  The arrays ``r_mass``, ``f_mass`` and ``d_f`` hold ``R_i``, ``f_i`` and ``d_i``; ``ln_r_mass``
  is ``ln R_i`` (checks only);
* ``ln E(R)`` on a uniform grid in ``ln R`` (800 points per decade by default), the exact
  inverse of that closed-form ``R(E)``, ``phi = log1p(d_i (r - R_i) / (h f_i)) / d_i`` on the
  bin that contains ``r``, sampled at those ranges (the end points are pinned).

Reading the stopping power is the shared bin location (``log_bin_index``), two array reads and
the shared interpolation (``interp_exp``); reading the range is the same bin and ``range_in_bin``;
the Python and Warp backends differ only in the memory access. The projectile must be a proton so
that MeV per nucleon equals MeV.

Units: energies MeV, mass stopping power MeV cm2/g, ranges g/cm2, ``inv_rho_xs`` cm2/g.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from types import MappingProxyType, SimpleNamespace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc._frozen import freeze_array
from ionmc.materials import Material
from ionmc.physics import stopping as _stopping
from ionmc.physics.projectiles import Projectile
from ionmc.physics.scattering import inverse_scattering_length_cm2_per_g
from ionmc.physics.stopping import RANGE_CONSTRUCTION, StoppingTable

MIN_POINTS_PER_DECADE = 200
DEFAULT_RANGE_POINTS_PER_DECADE = 800


@dataclass(frozen=True, eq=False)
class TransportTables:
    """Per-material transport tables (float64 numpy arrays; see module docstring).

    Arrays have a leading material axis. ``ln_s_mass`` / ``ln_r_mass`` have ``n_e`` columns,
    ``ln_e_of_r`` has ``n_r`` columns. ``ln_e0``/``inv_dln_e`` and ``ln_r0``/``inv_dln_r``
    locate the uniform grids. ``identity`` lists the source and I-value per material and
    ``sha256`` hashes all numerical content. The hashed arrays are float64 results of ``exp``,
    ``expm1`` and ``log`` whose last bits depend on the platform libm, so ``sha256`` is machine
    specific: it is the provenance of a run on one machine (equal for equal inputs there), not
    an identity across machines. The portable identity is the ``identity`` dict (source,
    content and material hashes, grid limits, range construction).
    """

    projectile: Projectile
    materials: tuple[Material, ...]
    e_min_mev: NDArray[np.float64]
    e_max_mev: NDArray[np.float64]
    ln_e0: NDArray[np.float64]
    inv_dln_e: NDArray[np.float64]
    ln_s_mass: NDArray[np.float64]
    ln_r_mass: NDArray[np.float64]
    r_mass: NDArray[np.float64]
    f_mass: NDArray[np.float64]
    d_f: NDArray[np.float64]
    r_min_g_cm2: NDArray[np.float64]
    r_max_g_cm2: NDArray[np.float64]
    ln_r0: NDArray[np.float64]
    inv_dln_r: NDArray[np.float64]
    ln_e_of_r: NDArray[np.float64]
    z_over_a: NDArray[np.float64]
    inv_rho_xs_cm2_g: NDArray[np.float64]
    nominal_density_g_cm3: NDArray[np.float64]
    identity: tuple[dict[str, Any], ...]
    sha256: str
    water_ln_s_mass: NDArray[np.float64] | None = None
    water_ln_e0: float = 0.0
    water_inv_dln_e: float = 0.0
    water_density_g_cm3: float = 0.0
    water_identity: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        for name in _ARRAY_FIELDS:
            a = getattr(self, name)
            if not isinstance(a, np.ndarray):
                raise ValueError(f"{name} must be a numpy array")
            object.__setattr__(self, name, freeze_array(a, np.float64, name))
        object.__setattr__(self, "identity", _freeze(self.identity))
        if self.water_ln_s_mass is not None:
            w = np.asarray(self.water_ln_s_mass, dtype=np.float64)
            if (
                w.ndim != 1
                or w.size < 2
                or not np.all(np.isfinite(w))
                or not (math.isfinite(self.water_ln_e0) and self.water_inv_dln_e > 0.0)
                or not self.water_density_g_cm3 > 0.0
                or self.water_identity is None
            ):
                raise ValueError("the water row must be a finite 1-D table with identity")
            object.__setattr__(self, "water_ln_s_mass", freeze_array(w, np.float64, "water row"))
            object.__setattr__(self, "water_identity", _freeze(self.water_identity))
        nm = len(self.materials)
        if nm < 1:
            raise ValueError("tables need at least one material")
        for name in (
            "e_min_mev",
            "e_max_mev",
            "ln_e0",
            "inv_dln_e",
            "r_min_g_cm2",
            "r_max_g_cm2",
            "ln_r0",
            "inv_dln_r",
            "z_over_a",
            "inv_rho_xs_cm2_g",
            "nominal_density_g_cm3",
        ):
            a = getattr(self, name)
            if a.shape != (nm,) or not np.all(np.isfinite(a)):
                raise ValueError(f"{name} must be finite with shape ({nm},)")
        for name in ("ln_s_mass", "ln_r_mass", "r_mass", "f_mass", "d_f", "ln_e_of_r"):
            a = getattr(self, name)
            if a.ndim != 2 or a.shape[0] != nm or a.shape[1] < 2 or not np.all(np.isfinite(a)):
                raise ValueError(f"{name} must be finite with shape ({nm}, n >= 2)")
        if np.any(np.diff(self.ln_r_mass, axis=1) <= 0.0) or np.any(
            np.diff(self.ln_e_of_r, axis=1) <= 0.0
        ):
            raise ValueError("range and inverse-range tables must be strictly increasing")
        if np.any(self.inv_dln_e <= 0.0) or np.any(self.inv_dln_r <= 0.0):
            raise ValueError("grid steps must be positive")

    @property
    def has_water(self) -> bool:
        """True if the water stopping row of the LET definition (decision 0040) is present."""
        return self.water_ln_s_mass is not None

    @property
    def n_water(self) -> int:
        """Grid points of the water row (0 without a row)."""
        return 0 if self.water_ln_s_mass is None else int(self.water_ln_s_mass.size)

    @property
    def n_materials(self) -> int:
        """Number of materials."""
        return len(self.materials)

    @property
    def n_e(self) -> int:
        """Grid points of the energy grid."""
        return int(self.ln_s_mass.shape[1])

    @property
    def n_r(self) -> int:
        """Grid points of the range grid."""
        return int(self.ln_e_of_r.shape[1])

    @classmethod
    def from_stopping_tables(
        cls,
        tables: Sequence[StoppingTable],
        points_per_decade: int = MIN_POINTS_PER_DECADE,
        range_points_per_decade: int = DEFAULT_RANGE_POINTS_PER_DECADE,
        water: StoppingTable | None = None,
    ) -> TransportTables:
        """Build transport tables from one proton :class:`StoppingTable` per material.

        ``water`` (the proton table of liquid water from the *same* ``StoppingSource``) adds the
        water stopping row of the LET definition: ``ln S`` on a uniform ``ln E`` grid of at least
        ``points_per_decade``, resampled like a material row; its identity is recorded and enters
        the content hash (without it the hash is that of the tables without a water row)."""
        if not tables:
            raise ValueError("at least one stopping table is required")
        if points_per_decade < MIN_POINTS_PER_DECADE or range_points_per_decade < 200:
            raise ValueError(f"tables need at least {MIN_POINTS_PER_DECADE} points per decade")
        projectile = tables[0].projectile
        if projectile.a != 1 or projectile.z != 1:
            raise ValueError("transport tables are implemented for protons only (MeV/u = MeV)")
        for t in tables:
            if t.projectile != projectile:
                raise ValueError("all stopping tables must describe the same projectile")
        n_e = 0
        n_r = 0
        for t in tables:
            n_e = max(n_e, _count(t.energy_per_u[-1] / t.energy_per_u[0], points_per_decade))
            n_r = max(
                n_r, _count(t.csda_range_g_cm2[-1] / t.csda_range_g_cm2[0], range_points_per_decade)
            )
        nm = len(tables)
        e_min = np.empty(nm)
        e_max = np.empty(nm)
        ln_e0 = np.empty(nm)
        inv_dln_e = np.empty(nm)
        ln_s = np.empty((nm, n_e))
        ln_r = np.empty((nm, n_e))
        r_mass = np.empty((nm, n_e))
        f_mass = np.empty((nm, n_e))
        d_f = np.zeros((nm, n_e))
        r_min = np.empty(nm)
        r_max = np.empty(nm)
        ln_r0 = np.empty(nm)
        inv_dln_r = np.empty(nm)
        ln_e_of_r = np.empty((nm, n_r))
        for m, t in enumerate(tables):
            lo, hi = math.log(t.energy_per_u[0]), math.log(t.energy_per_u[-1])
            grid = np.linspace(lo, hi, n_e)
            e = np.exp(grid)
            e[0], e[-1] = t.energy_per_u[0], t.energy_per_u[-1]
            lns = np.interp(grid, np.log(t.energy_per_u), np.log(t.s_el_mass))
            ln_f = grid - lns  # ln(E / S): dR/d ln E is f = E / S for a proton
            inc = _stopping.exact_loglog_range_increments(grid, lns, 1.0)
            r = float(t.csda_range_g_cm2[0]) + np.concatenate(([0.0], np.cumsum(inc)))
            r_mass[m] = r
            f_mass[m] = np.exp(ln_f)
            d_f[m, :-1] = np.diff(ln_f)
            e_min[m], e_max[m] = e[0], e[-1]
            ln_e0[m] = grid[0]
            inv_dln_e[m] = (n_e - 1) / (hi - lo)
            ln_s[m] = lns
            ln_r[m] = np.log(r)
            r_min[m], r_max[m] = r[0], r[-1]
            rlo, rhi = float(np.log(r[0])), float(np.log(r[-1]))
            rgrid = np.linspace(rlo, rhi, n_r)
            ln_r0[m] = rgrid[0]
            inv_dln_r[m] = (n_r - 1) / (rhi - rlo)
            ln_e_of_r[m] = _exact_inverse(
                np.exp(rgrid), r, f_mass[m], d_f[m], 1.0 / inv_dln_e[m], grid
            )
            ln_e_of_r[m, 0], ln_e_of_r[m, -1] = grid[0], grid[-1]
        materials = tuple(t.material for t in tables)
        z_over_a = np.array([mat.z_over_a for mat in materials])
        inv_xs = np.array([inverse_scattering_length_cm2_per_g(mat) for mat in materials])
        density = np.array([mat.density_g_cm3 for mat in materials])
        identity = tuple(
            {
                "material": mat.name,
                "projectile": projectile.name,
                "source": t.metadata.get("source"),
                "I_eV_requested": mat.mean_excitation_eV,
                # the I value of the data actually used (a NIST table is at 75 eV whatever
                # the requested material says); falls back to the material's value
                "I_eV_effective": float(t.metadata.get("I_eV", mat.mean_excitation_eV)),
                "dataset_id": t.metadata.get("dataset_id"),
                "source_sha256": t.metadata.get("sha256"),
                "content_sha256": t.metadata.get("content_sha256"),
                "material_sha256": material_fingerprint(mat),
                "e_min_mev": float(e_min[i]),
                "e_max_mev": float(e_max[i]),
                "range_construction": RANGE_CONSTRUCTION,
                "metadata": _jsonable(t.metadata),
            }
            for i, (mat, t) in enumerate(zip(materials, tables, strict=True))
        )
        h = hashlib.sha256()
        for arr in (
            e_min,
            e_max,
            ln_s,
            ln_r,
            r_mass,
            f_mass,
            d_f,
            ln_e_of_r,
            z_over_a,
            inv_xs,
            density,
        ):
            h.update(np.ascontiguousarray(arr).tobytes())
        h.update(json.dumps(identity, sort_keys=True).encode())
        water_kwargs: dict[str, Any] = {}
        if water is not None:
            if water.projectile != projectile:
                raise ValueError("the water table must describe the same projectile")
            n_w = _count(water.energy_per_u[-1] / water.energy_per_u[0], points_per_decade)
            wlo, whi = math.log(water.energy_per_u[0]), math.log(water.energy_per_u[-1])
            wgrid = np.linspace(wlo, whi, n_w)
            w_ln_s = np.interp(wgrid, np.log(water.energy_per_u), np.log(water.s_el_mass))
            w_identity = {
                "material": water.material.name,
                "projectile": projectile.name,
                "role": "LET medium (decision 0040)",
                "source": water.metadata.get("source"),
                "I_eV_effective": float(
                    water.metadata.get("I_eV", water.material.mean_excitation_eV)
                ),
                "dataset_id": water.metadata.get("dataset_id"),
                "source_sha256": water.metadata.get("sha256"),
                "content_sha256": water.metadata.get("content_sha256"),
                "material_sha256": material_fingerprint(water.material),
                "e_min_mev": float(water.energy_per_u[0]),
                "e_max_mev": float(water.energy_per_u[-1]),
                "n_points": int(n_w),
                "metadata": _jsonable(water.metadata),
            }
            h.update(np.ascontiguousarray(w_ln_s).tobytes())
            h.update(json.dumps(w_identity, sort_keys=True).encode())
            water_kwargs = {
                "water_ln_s_mass": w_ln_s,
                "water_ln_e0": float(wgrid[0]),
                "water_inv_dln_e": float((n_w - 1) / (whi - wlo)),
                "water_density_g_cm3": float(water.material.density_g_cm3),
                "water_identity": w_identity,
            }
        return cls(
            projectile=projectile,
            materials=materials,
            e_min_mev=e_min,
            e_max_mev=e_max,
            ln_e0=ln_e0,
            inv_dln_e=inv_dln_e,
            ln_s_mass=ln_s,
            ln_r_mass=ln_r,
            r_mass=r_mass,
            f_mass=f_mass,
            d_f=d_f,
            r_min_g_cm2=r_min,
            r_max_g_cm2=r_max,
            ln_r0=ln_r0,
            inv_dln_r=inv_dln_r,
            ln_e_of_r=ln_e_of_r,
            z_over_a=z_over_a,
            inv_rho_xs_cm2_g=inv_xs,
            nominal_density_g_cm3=density,
            identity=identity,
            sha256=h.hexdigest(),
            **water_kwargs,
        )

    def _locate(
        self, material: int, x: float, ln0: float, inv_dl: float, n: int
    ) -> tuple[int, float]:
        t = (math.log(max(x, 1.0e-30)) - ln0) * inv_dl
        i = max(min(int(math.floor(t)), n - 2), 0)
        return i, min(max(t - i, 0.0), 1.0)

    def stopping_mass(self, material: int, energy_mev: float) -> float:
        """Mass stopping power [MeV cm2/g] at ``energy_mev`` (log-log interpolation)."""
        i, f = self._locate(
            material, energy_mev, self.ln_e0[material], self.inv_dln_e[material], self.n_e
        )
        row = self.ln_s_mass[material]
        return math.exp(row[i] * (1.0 - f) + row[i + 1] * f)

    def s_water(self, energy_mev: float, species: str = "proton") -> float:
        """Linear unrestricted electronic stopping power of ``species`` in water [MeV/mm]
        (numerically keV/um) at the kinetic energy ``energy_mev`` (log-log interpolation of the
        water row; clamped to the table range). Only the proton row exists (V3-004); the species
        argument is the interface for V3-005A (``S_w(E, species)``)."""
        if self.water_ln_s_mass is None:
            raise ValueError("these transport tables have no water stopping row")
        if species != self.projectile.name:
            raise ValueError(
                f"the water row describes {self.projectile.name!r}, not species {species!r}"
            )
        i, f = self._locate(
            0, energy_mev, self.water_ln_e0, self.water_inv_dln_e, int(self.water_ln_s_mass.size)
        )
        row = self.water_ln_s_mass
        mass = math.exp(row[i] * (1.0 - f) + row[i + 1] * f)
        return mass * self.water_density_g_cm3 / 10.0

    def range_g_cm2(self, material: int, energy_mev: float) -> float:
        """CSDA range [g/cm2] at ``energy_mev``: the exact integral of the log-log interpolated
        ``S`` (closed form in the bin, the same arithmetic as the shared ``range_in_bin``)."""
        i, f = self._locate(
            material, energy_mev, self.ln_e0[material], self.inv_dln_e[material], self.n_e
        )
        d = float(self.d_f[material, i])
        x = d * f
        g = f * (1.0 + x * 0.5 + x * x / 6.0)
        if abs(x) >= 1.0e-5:
            g = (math.exp(x) - 1.0) / d
        h = 1.0 / float(self.inv_dln_e[material])
        return float(self.r_mass[material, i]) + h * float(self.f_mass[material, i]) * g

    def energy_from_range(self, material: int, range_g_cm2: float) -> float:
        """Energy [MeV] whose CSDA range is ``range_g_cm2``; clamped to the table range."""
        i, f = self._locate(
            material, range_g_cm2, self.ln_r0[material], self.inv_dln_r[material], self.n_r
        )
        row = self.ln_e_of_r[material]
        return math.exp(row[i] * (1.0 - f) + row[i + 1] * f)

    def to_warp(self, device: str, dtype: Any) -> SimpleNamespace:
        """Copy the tables to Warp arrays of precision ``dtype`` (``wp.float32``/``float64``)
        on ``device``; table values are stored as logarithms."""
        import warp as wp

        np_dtype = np.float32 if dtype is wp.float32 else np.float64
        if dtype is not wp.float32 and dtype is not wp.float64:
            raise ValueError("dtype must be wp.float32 or wp.float64")

        def arr2(a: NDArray[np.float64]) -> Any:
            return wp.array(a.astype(np_dtype), dtype=dtype, device=device)

        water: dict[str, Any] = {}
        if self.water_ln_s_mass is not None:
            water = {
                "ln_s_water": wp.array(
                    self.water_ln_s_mass.astype(np_dtype), dtype=dtype, device=device
                ),
                "ln_e0_water": float(self.water_ln_e0),
                "inv_dln_e_water": float(self.water_inv_dln_e),
                "rho_water": float(self.water_density_g_cm3),
                "n_water": self.n_water,
            }
        return SimpleNamespace(
            **water,
            ln_s_mass=arr2(self.ln_s_mass),
            ln_r_mass=arr2(self.ln_r_mass),
            r_mass=arr2(self.r_mass),
            f_mass=arr2(self.f_mass),
            d_f=arr2(self.d_f),
            ln_e_of_r=arr2(self.ln_e_of_r),
            ln_e0=arr2(self.ln_e0),
            inv_dln_e=arr2(self.inv_dln_e),
            ln_r0=arr2(self.ln_r0),
            inv_dln_r=arr2(self.inv_dln_r),
            z_over_a=arr2(self.z_over_a),
            inv_rho_xs=arr2(self.inv_rho_xs_cm2_g),
            n_e=self.n_e,
            n_r=self.n_r,
        )


_ARRAY_FIELDS = (
    "e_min_mev",
    "e_max_mev",
    "ln_e0",
    "inv_dln_e",
    "ln_s_mass",
    "ln_r_mass",
    "r_mass",
    "f_mass",
    "d_f",
    "r_min_g_cm2",
    "r_max_g_cm2",
    "ln_r0",
    "inv_dln_r",
    "ln_e_of_r",
    "z_over_a",
    "inv_rho_xs_cm2_g",
    "nominal_density_g_cm3",
)


def _exact_inverse(
    r: NDArray[np.float64],
    r_nodes: NDArray[np.float64],
    f_nodes: NDArray[np.float64],
    d_nodes: NDArray[np.float64],
    h: float,
    grid: NDArray[np.float64],
) -> NDArray[np.float64]:
    """``ln E`` at the ranges ``r``: the exact inverse of the closed-form range on its bin,
    ``phi = log1p(d (r - R_i) / (h f_i)) / d`` (``(r - R_i) / (h f_i)`` for ``|d| < 1e-12``)."""
    n = r_nodes.size
    i = np.clip(np.searchsorted(r_nodes, r, side="right") - 1, 0, n - 2)
    d = d_nodes[i]
    x = (r - r_nodes[i]) / (h * f_nodes[i])
    flat = np.abs(d) < 1.0e-12
    phi = np.where(flat, x, np.log1p(d * x) / np.where(flat, 1.0, d))
    return grid[i] + np.clip(phi, 0.0, 1.0) * h


def _freeze(value: Any) -> Any:
    """Recursively immutable copy: dict -> MappingProxyType, list/tuple -> tuple."""
    if isinstance(value, dict | MappingProxyType):
        return MappingProxyType({k: _freeze(v) for k, v in value.items()})
    if isinstance(value, list | tuple):
        return tuple(_freeze(v) for v in value)
    return value


def thaw(value: Any) -> Any:
    """Deep mutable copy (dict/list) of a frozen structure, for summaries."""
    if isinstance(value, MappingProxyType | dict):
        return {k: thaw(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [thaw(v) for v in value]
    return value


def material_fingerprint(material: Material) -> str:
    """SHA-256 of the complete physics definition of a material (composition, density, I,
    density-effect parameters)."""
    sp = material.sternheimer
    payload = {
        "name": material.name,
        "density_g_cm3": material.density_g_cm3,
        "mass_fractions": {k: material.mass_fractions[k] for k in sorted(material.mass_fractions)},
        "I_eV": material.I_eV,
        "sternheimer": None
        if sp is None
        else {"x0": sp.x0, "x1": sp.x1, "cbar": sp.cbar, "a": sp.a, "m": sp.m},
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _jsonable(value: Any) -> Any:
    """Canonical JSON-compatible copy of table metadata (recorded in identity, hash, summary)."""
    return json.loads(json.dumps(value, sort_keys=True, default=str))


def _count(ratio: float, points_per_decade: int) -> int:
    return max(2, int(math.ceil(points_per_decade * math.log10(ratio))) + 1)
