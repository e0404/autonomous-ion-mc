"""Stopping-power and range tables consumed by the transport layer.

A :class:`StoppingTable` holds, on a fixed logarithmic grid of kinetic energy
per nucleon, the electronic mass stopping power (MeV cm²/g) and the CSDA
range per unit density (g/cm²) of one species in one material, together
with the provenance of every source that contributed. Transport kernels only
see these arrays (decision 0038); whether a value came from the analytic
Bethe layer (``ionmc.physics.stopping``) or from external tabulated data is a
property of the table's provenance, not of the kernels.

Construction (:func:`build_stopping_table`):

* above the blend window: analytic Bethe layer with corrections and the
  material's mean excitation energy;
* below the blend window: an external low-energy source — ICRU 90 tables for
  water (protons, alphas, carbon ions; ``ionmc.data.icru90``) or NIST
  PSTAR/ASTAR for other materials (``ionmc.data.nist_star``) — scaled for
  other ions by the ratio of effective charges squared at equal velocity;
* inside the window (8–16 MeV/u for z ≤ 2, 16–32 MeV/u for z ≥ 3): a smooth
  cosine blend in ln t, so tables are continuous; the two sources agree
  within ≈ 0.25 % (p, He) and ≈ 0.6 % (C) there at equal mean excitation
  energy (validated against ICRU 90).

When no external source is available the caller may request
``low_energy="bethe-extrapolated"``; that table is marked as such in its
provenance and documented as inaccurate below 8 MeV/u.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ionmc.materials import Material, get_material
from ionmc.physics import stopping
from ionmc.species import Species, get_species

GRID_T_MIN = 1.0e-3  # MeV/u
GRID_T_MAX = 1.0e3  # MeV/u
GRID_POINTS = 721  # 120 points per decade
BLEND_WINDOWS = {  # MeV/u: (low, high) by projectile charge class
    "light": (8.0, 16.0),  # z <= 2: Bethe layer within 0.25 % of ICRU 90 above 8 MeV/u
    "heavy": (16.0, 32.0),  # z >= 3: carbon within 0.6 % above 16 MeV/u, 0.3 % above 30
}


def blend_window(species: Species) -> tuple[float, float]:
    return BLEND_WINDOWS["light" if species.z <= 2 else "heavy"]


LowEnergySource = Callable[
    [Species, Material, np.ndarray], tuple[np.ndarray, dict[str, Any]]
]
"""Returns proton-equivalent-scaled electronic mass stopping power on the given
energies-per-nucleon grid and a provenance dictionary."""


def log_grid(
    t_min: float = GRID_T_MIN, t_max: float = GRID_T_MAX, n: int = GRID_POINTS
) -> np.ndarray:
    return np.geomspace(t_min, t_max, n)


@dataclass(frozen=True)
class StoppingTable:
    species: Species
    material: Material
    t_mev_per_u: np.ndarray  # log-spaced grid, MeV/u
    electronic_mev_cm2_g: np.ndarray  # mass electronic stopping power
    csda_range_g_cm2: np.ndarray  # CSDA range for the ion (total), g/cm²
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def log_t(self) -> np.ndarray:
        return np.log(self.t_mev_per_u)

    def stopping_at(self, t: np.ndarray | float) -> np.ndarray:
        """Mass stopping power at t (MeV/u), log-log interpolated."""
        lt = np.log(
            np.clip(
                np.asarray(t, dtype=np.float64),
                self.t_mev_per_u[0],
                self.t_mev_per_u[-1],
            )
        )
        return np.exp(np.interp(lt, self.log_t, np.log(self.electronic_mev_cm2_g)))

    def range_at(self, t: np.ndarray | float) -> np.ndarray:
        """CSDA range (g/cm²) at t (MeV/u), log-log interpolated."""
        lt = np.log(
            np.clip(
                np.asarray(t, dtype=np.float64),
                self.t_mev_per_u[0],
                self.t_mev_per_u[-1],
            )
        )
        return np.exp(np.interp(lt, self.log_t, np.log(self.csda_range_g_cm2)))

    def energy_at_range(self, r: np.ndarray | float) -> np.ndarray:
        """Inverse range: energy per nucleon (MeV/u) at CSDA range r (g/cm²)."""
        lr = np.log(
            np.clip(
                np.asarray(r, dtype=np.float64),
                self.csda_range_g_cm2[0],
                self.csda_range_g_cm2[-1],
            )
        )
        return np.exp(np.interp(lr, np.log(self.csda_range_g_cm2), self.log_t))


def csda_range(species: Species, t: np.ndarray, s_mass: np.ndarray) -> np.ndarray:
    """CSDA range (g/cm²) by trapezoidal integration in ln t.

    R(t) = A ∫_0^t dt'/S(t'). The contribution below the first grid point is
    approximated by assuming S ∝ t^{-1/2}... more precisely by extrapolating
    S constant, which for t_min = 1 keV/u is < 1e-5 g/cm² and irrelevant.
    """
    t = np.asarray(t, dtype=np.float64)
    integrand = t / s_mass  # dt / S = t dln t / S
    lnt = np.log(t)
    r = np.zeros_like(t)
    r[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(lnt))
    r += t[0] / s_mass[0]  # residual below t_min with constant S
    return species.a * r


def bethe_extrapolated_source(
    species: Species, material: Material, t: np.ndarray
) -> tuple[np.ndarray, dict[str, Any]]:
    """Fallback low-energy source: the Bethe layer itself (inaccurate < 8 MeV/u)."""
    s = stopping.mass_stopping_power(species, material, t)
    # Guard against the Bethe logarithm turning negative at very low velocity.
    s = np.where(np.isfinite(s) & (s > 0), s, np.nan)
    if np.isnan(s).any():
        first = int(np.argmax(~np.isnan(s)))
        s[:first] = (
            s[first] * (t[:first] / t[first]) ** 0.5
        )  # ~ velocity-proportional stopping
    return s, {
        "low_energy_source": "bethe-extrapolated",
        "accuracy_note": "not validated below 8 MeV/u",
    }


def build_stopping_table(
    species: Species | str,
    material: Material | str,
    low_energy: LowEnergySource | None = None,
    *,
    grid: np.ndarray | None = None,
) -> StoppingTable:
    species = get_species(species)
    material = get_material(material)
    t = log_grid() if grid is None else np.asarray(grid, dtype=np.float64)
    s_high = stopping.mass_stopping_power(species, material, t)
    source = low_energy or bethe_extrapolated_source
    s_low, low_prov = source(species, material, t)
    blend_low, blend_high = blend_window(species)
    w = np.clip(
        (np.log(t) - np.log(blend_low)) / (np.log(blend_high) - np.log(blend_low)),
        0.0,
        1.0,
    )
    w = 0.5 - 0.5 * np.cos(np.pi * w)  # smooth cosine blend in ln t
    s = (1.0 - w) * s_low + w * s_high
    if not np.all(np.isfinite(s)) or np.any(s <= 0):
        raise ValueError(
            f"invalid stopping table for {species.name} in {material.name}"
        )
    r = csda_range(species, t, s)
    prov = {
        "species": species.name,
        "material": material.name,
        "mean_excitation_ev": material.mean_excitation_ev,
        "high_energy_source": "bethe-corrected (ionmc.physics.stopping)",
        "high_energy_domain_mev_per_u": [blend_high, float(t[-1])],
        "blend_window_mev_per_u": [blend_low, blend_high],
        "grid": {
            "t_min": float(t[0]),
            "t_max": float(t[-1]),
            "points": int(t.size),
            "spacing": "log",
        },
        **low_prov,
    }
    return StoppingTable(species, material, t, s, r, prov)


def scaled_low_energy_source(
    proton_like: Callable[[Material, np.ndarray], tuple[np.ndarray, dict[str, Any]]],
    alpha_like: Callable[[Material, np.ndarray], tuple[np.ndarray, dict[str, Any]]]
    | None = None,
) -> LowEnergySource:
    """Build a low-energy source from proton (and optionally alpha) tabulations.

    For species other than the tabulated ones the stopping power at equal
    velocity is scaled by (z_eff,ion / z_eff,proton)²; helium uses the alpha
    table directly when provided.
    """

    def source(
        species: Species, material: Material, t: np.ndarray
    ) -> tuple[np.ndarray, dict[str, Any]]:
        if species.name == "he4" and alpha_like is not None:
            s, prov = alpha_like(material, t)
            return s, {**prov, "low_energy_scaling": "none (alpha table)"}
        s_p, prov = proton_like(material, t)
        if species.name == "proton":
            return s_p, {**prov, "low_energy_scaling": "none (proton table)"}
        beta, _ = stopping.kinematics(species, t)
        beta_p, _ = stopping.kinematics(get_species("proton"), t)
        zeff = stopping.effective_charge(species.z, beta)
        zeff_p = stopping.effective_charge(1, beta_p)
        return s_p * (zeff / zeff_p) ** 2, {
            **prov,
            "low_energy_scaling": (
                "effective-charge squared (Pierce-Blann) at equal velocity"
            ),
        }

    return source


def icru90_water_source(
    species: Species, material: Material, t: np.ndarray
) -> tuple[np.ndarray, dict[str, Any]]:
    """Low-energy source for liquid water from the ICRU 90 tables (p, He-4, C-12 direct;
    other ions by effective-charge scaling of the proton table)."""
    from ionmc.data import icru90

    if material.name != "water":
        raise ValueError("icru90_water_source applies to water only")
    prov_base = {
        "low_energy_source": "ICRU 90 liquid water (NIST addendum)",
        "low_energy_dataset": icru90.provenance()["dataset_id"],
        "low_energy_version": icru90.provenance()["version"],
    }
    if species.name in icru90.available_species():
        return icru90.water_table(species.name).electronic_at(t), {
            **prov_base,
            "low_energy_scaling": "none (species table)",
        }
    table = icru90.water_table("proton")
    beta, _ = stopping.kinematics(species, t)
    beta_p, _ = stopping.kinematics(get_species("proton"), t)
    scale = (
        stopping.effective_charge(species.z, beta)
        / stopping.effective_charge(1, beta_p)
    ) ** 2
    return table.electronic_at(t) * scale, {
        **prov_base,
        "low_energy_scaling": (
            "effective-charge squared (Pierce-Blann) of the proton table"
        ),
    }


def default_low_energy_source(
    cache: Any | None = None, *, offline: bool = False
) -> LowEnergySource:
    """Compose the default low-energy source.

    Water uses the shipped ICRU 90 tables. Materials with a NIST STAR code use
    PSTAR (and ASTAR for helium-4) through the data cache; anything else falls
    back to the Bethe-extrapolated layer, which the table provenance records.
    """
    from ionmc.data import DataCache, nist_star

    cache = cache if cache is not None else DataCache()

    def star_source(
        program: str,
    ) -> Callable[[Material, np.ndarray], tuple[np.ndarray, dict[str, Any]]]:
        def source(
            material: Material, t: np.ndarray
        ) -> tuple[np.ndarray, dict[str, Any]]:
            code_name = nist_star.material_key(material)
            if code_name is None:
                raise ValueError(f"no NIST STAR table for {material.name}")
            table, record = nist_star.load_star_table(
                cache, program, code_name, offline=offline
            )
            t_tab = table.energy_per_nucleon_mev()
            tt = np.clip(t, t_tab.min(), t_tab.max())
            values = np.asarray(table.electronic_at(tt), dtype=np.float64)
            return values, {
                "low_energy_source": f"NIST {program} {code_name} (ICRU 49)",
                "low_energy_dataset": record.dataset_id,
                "low_energy_version": record.version,
                "low_energy_sha256": record.sha256,
            }

        return source

    scaled = scaled_low_energy_source(star_source("PSTAR"), star_source("ASTAR"))

    def source(
        species: Species, material: Material, t: np.ndarray
    ) -> tuple[np.ndarray, dict[str, Any]]:
        if material.name == "water":
            return icru90_water_source(species, material, t)
        if nist_star.material_key(material) is not None:
            return scaled(species, material, t)
        return water_shape_source(species, material, t)

    return source


def water_shape_source(
    species: Species, material: Material, t: np.ndarray
) -> tuple[np.ndarray, dict[str, Any]]:
    """Low-energy source for materials without tabulated data (e.g. ICRP muscle, lung).

    Uses the ICRU 90 water shape scaled by the Bethe-layer ratio
    S_material/S_water evaluated at the low edge of the blend window, i.e. the
    material's electron density and mean excitation energy set the level while
    the energy dependence below the window follows water. Documented as an
    approximation valid for water-like tissues; the affected residual range is
    below ~1 mm.
    """
    water = get_material("water")
    s_water, prov = icru90_water_source(species, water, t)
    t_ref = blend_window(species)[0]
    ratio = float(
        stopping.mass_stopping_power(species, material, np.array([t_ref]))[0]
        / stopping.mass_stopping_power(species, water, np.array([t_ref]))[0]
    )
    return s_water * ratio, {
        **prov,
        "low_energy_source": (
            f"ICRU 90 water shape scaled by Bethe ratio {ratio:.4f} at "
            f"{t_ref} MeV/u (no tabulated data for {material.name})"
        ),
    }
