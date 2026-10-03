"""Electronic stopping power of ions: analytic Bethe-Bloch model and tables.

Model (mass electronic stopping power, units MeV cm2/g), after ICRU Report 49 (1993),
ICRU Report 90 (J. ICRU 14(1), 2014/2016) and the PDG review "Passage of particles through
matter"::

    S/rho = K q^2 (Z/A) beta^-2 [ L0 + q L1 + q^2 L2' + M ]

with K = 4 pi N_A r_e^2 m_e c^2 = 0.307075 MeV cm2/mol, q the (effective) projectile
charge and

* ``L0 = 1/2 ln(2 m_e c^2 beta^2 gamma^2 Tmax / I^2) - beta^2 - delta/2 - C/Z``
  with the exact ``Tmax = 2 m_e c^2 beta^2 gamma^2 / (1 + 2 gamma m_e/M + (m_e/M)^2)``;
* ``delta``: Sternheimer density effect (Sternheimer, Berger, Seltzer, ADNDT 30 (1984) 261),
  zero when the material has no parameters;
* ``C/Z``: shell correction in the Bichsel parameterisation of ICRU 49 as given by the PDG
  (older editions): with eta = beta*gamma and I in eV,
  ``C = (0.422377 eta^-2 + 0.0304043 eta^-4 - 0.00038106 eta^-6) 1e-6 I^2
  + (3.858019 eta^-2 - 0.1667989 eta^-4 + 0.00157955 eta^-6) 1e-9 I^3``,
  valid for eta >= 0.13 (about 8 MeV/u for protons); below, C is held at its value at
  eta = 0.13. For compounds C/Z is the electron-weighted sum over elements of C(I_i)/Z_i
  with the elemental I_i (additivity per electron);
* ``q L1``: Barkas term, the Ashley-Ritchie-Brandt (1972) form as implemented in Geant4
  ``G4EmCorrections::BarkasCorrection`` (v11.4.2): per element
  ``1.29 q F(b/sqrt(X)) / (sqrt(Z) X^(3/2))`` with ``X = beta^2/(alpha^2 Z)``, tabulated
  function F (47 points from the Geant4 source, linear interpolation, scaled by Wmax/W
  beyond W = 10) and shell-dependent ``b`` (1.8 for Z <= 10 and Ar, 1.4 for 11 <= Z <= 17
  and 19 <= Z <= 25), atom-number weighted over the elements;
* ``q^2 L2'``: Bloch term ``-y^2 sum_{n>=1} 1/(n (n^2 + y^2))``, ``y = q alpha / beta``;
* ``M``: Mott term ``0.5 pi alpha beta q`` for projectiles with z >= 2 (as in Geant4);
* ``q``: Pierce-Blann effective charge ``z (1 - exp(-125 beta z^(-2/3)))`` for z >= 2
  (Barkas/Pierce-Blann form, ICRU 49 and G4/Ziegler-type effective charge); q = z for z = 1.

Domain: energy per nucleon >= 1 MeV/u; no value is returned below (ValueError).
All arithmetic is float64.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ionmc.data.nist_star import StarTable
from ionmc.materials import Material
from ionmc.physics.projectiles import ELECTRON_MASS_MEV, Projectile

K_MEV_CM2_MOL = 0.307075
ALPHA_FS = 1.0 / 137.035999084
E_MIN_PER_U_MEV = 1.0
ETA_SHELL_MIN = 0.13

# Ashley-Ritchie-Brandt function F(W) as tabulated in G4EmCorrections.cc (Geant4 v11.4.2).
# This table derives from Geant4 source; see src/ionmc/THIRD_PARTY_NOTICES.md
# (Geant4 Software License).
_BARKAS_F = np.array(
    [
        (0.02, 21.5), (0.03, 20.0), (0.04, 18.0), (0.05, 15.6), (0.06, 15.0), (0.07, 14.0),
        (0.08, 13.5), (0.09, 13.0), (0.1, 12.2), (0.2, 9.25), (0.3, 7.0), (0.4, 6.0),
        (0.5, 4.5), (0.6, 3.5), (0.7, 3.0), (0.8, 2.5), (0.9, 2.0), (1.0, 1.7), (1.2, 1.2),
        (1.3, 1.0), (1.4, 0.86), (1.5, 0.7), (1.6, 0.61), (1.7, 0.52), (1.8, 0.5), (1.9, 0.43),
        (2.0, 0.42), (2.1, 0.3), (2.4, 0.2), (3.0, 0.13), (3.08, 0.1), (3.1, 0.09),
        (3.3, 0.08), (3.5, 0.07), (3.8, 0.06), (4.0, 0.051), (4.1, 0.04), (4.8, 0.03),
        (5.0, 0.024), (5.1, 0.02), (6.0, 0.013), (6.5, 0.01), (7.0, 0.009), (7.1, 0.008),
        (8.0, 0.006), (9.0, 0.0032), (10.0, 0.0025),
    ],
    dtype=np.float64,
)  # fmt: skip
_W_MAX = 10.0


def _barkas_b(z_el: int) -> float:
    """Shell parameter b of the Barkas function for target element number ``z_el``."""
    if z_el == 1:
        return 1.8
    if z_el == 2:
        return 0.6
    if z_el <= 10:
        return 1.8
    if z_el <= 17:
        return 1.4
    if z_el == 18:
        return 1.8
    if z_el <= 25:
        return 1.4
    if z_el <= 50:
        return 1.35
    return 1.3


@dataclass(frozen=True)
class BetheOptions:
    """Switches of the analytic model (all on by default); recorded in table metadata."""

    shell: bool = True
    barkas: bool = True
    bloch: bool = True
    mott: bool = True
    effective_charge: bool = True


def kinematics(
    energy_per_u_mev: ArrayLike, projectile: Projectile
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Return ``(beta^2, (beta*gamma)^2, Tmax [MeV])`` for kinetic energy per nucleon [MeV/u].

    The total kinetic energy is ``energy_per_u * projectile.a``.
    """
    e = np.asarray(energy_per_u_mev, dtype=np.float64) * projectile.a
    gamma = 1.0 + e / projectile.mass_mev
    bg2 = gamma * gamma - 1.0
    beta2 = bg2 / (gamma * gamma)
    ratio = ELECTRON_MASS_MEV / projectile.mass_mev
    tmax = 2.0 * ELECTRON_MASS_MEV * bg2 / (1.0 + 2.0 * gamma * ratio + ratio * ratio)
    return beta2, bg2, tmax


def density_effect(bg2: NDArray[np.float64], material: Material) -> NDArray[np.float64]:
    """Sternheimer density-effect delta (dimensionless) at ``(beta*gamma)^2``; 0 if no data."""
    sp = material.sternheimer
    if sp is None:
        return np.zeros_like(bg2)
    x = 0.5 * np.log10(bg2)
    ln10 = math.log(10.0)
    delta = np.where(x >= sp.x1, 2.0 * ln10 * x - sp.cbar, 0.0)
    mid = (x >= sp.x0) & (x < sp.x1)
    delta_mid = 2.0 * ln10 * x - sp.cbar + sp.a * np.maximum(sp.x1 - x, 0.0) ** sp.m
    return np.where(mid, delta_mid, delta)


def _shell_c(eta: NDArray[np.float64], i_ev: float) -> NDArray[np.float64]:
    """Bichsel/ICRU 49 shell correction C (dimensionless) for one atom, I in eV; eta >= 0.13."""
    e = np.maximum(eta, ETA_SHELL_MIN)
    t1 = 0.422377 * e**-2 + 0.0304043 * e**-4 - 0.00038106 * e**-6
    t2 = 3.858019 * e**-2 - 0.1667989 * e**-4 + 0.00157955 * e**-6
    return t1 * 1e-6 * i_ev**2 + t2 * 1e-9 * i_ev**3


def shell_correction_over_z(bg2: NDArray[np.float64], material: Material) -> NDArray[np.float64]:
    """Shell correction C/Z (dimensionless) of ``material`` at ``(beta*gamma)^2``.

    Electron-weighted sum over elements of C(I_i)/Z_i (Bichsel form, held at eta = 0.13
    below its validity limit); see module docstring.
    """
    eta = np.sqrt(bg2)
    za = material.z_over_a
    out = np.zeros_like(bg2)
    for el, w in material._fractions:
        electron_fraction = w * el.Z / el.A_g_mol / za
        out += electron_fraction * _shell_c(eta, el.I_eV) / el.Z
    return out


def barkas_term(
    beta2: NDArray[np.float64], q: NDArray[np.float64], material: Material
) -> NDArray[np.float64]:
    """Barkas term ``q L1`` (dimensionless, added to L0) from the Ashley-Ritchie-Brandt form."""
    total = np.zeros_like(beta2)
    n_atoms = 0.0
    for el, w in material._fractions:
        n_i = w / el.A_g_mol  # proportional to atoms per volume
        n_atoms += n_i
        x = beta2 / (ALPHA_FS**2 * el.Z)
        wv = _barkas_b(el.Z) / np.sqrt(x)
        f = np.interp(wv, _BARKAS_F[:, 0], _BARKAS_F[:, 1])
        f = np.where(wv > _W_MAX, f * (_W_MAX / wv), f)
        total += n_i * f / (np.sqrt(el.Z * x) * x)
    return 1.29 * q * total / n_atoms


def bloch_term(beta2: NDArray[np.float64], q: NDArray[np.float64]) -> NDArray[np.float64]:
    """Bloch term ``-y^2 sum_n 1/(n(n^2+y^2))`` (dimensionless), ``y = q alpha/beta``.

    The series is summed to n = 2000 with an integral estimate of the tail.
    """
    y2 = (q * ALPHA_FS) ** 2 / beta2
    n_terms = 2000
    n = np.arange(1, n_terms + 1, dtype=np.float64)
    shape = y2.shape
    flat = y2.reshape(-1)
    s = np.zeros_like(flat)
    for start in range(0, flat.size, 256):
        chunk = flat[start : start + 256, None]
        s[start : start + 256] = np.sum(1.0 / (n * (n * n + chunk)), axis=1)
    s += 0.5 / (n_terms + 0.5) ** 2
    return (-flat * s).reshape(shape)


def pierce_blann_charge(beta: NDArray[np.float64], projectile: Projectile) -> NDArray[np.float64]:
    """Pierce-Blann effective charge ``z (1 - exp(-125 beta z^(-2/3)))`` in units of e (z >= 2)."""
    z = float(projectile.z)
    q: NDArray[np.float64] = z * (1.0 - np.exp(-125.0 * beta * z ** (-2.0 / 3.0)))
    return q


def bethe_mass_stopping(
    energy_per_u_mev: ArrayLike,
    projectile: Projectile,
    material: Material,
    *,
    shell: bool = True,
    barkas: bool = True,
    bloch: bool = True,
    mott: bool = True,
    effective_charge: bool = True,
) -> NDArray[np.float64]:
    """Mass electronic stopping power [MeV cm2/g] at kinetic energy per nucleon [MeV/u].

    Implements the model in the module docstring. Raises ``ValueError`` for energies
    below 1 MeV/u. The flags switch the shell, Barkas, Bloch and Mott terms and the
    Pierce-Blann effective charge (the latter only acts for z >= 2).
    """
    e = np.atleast_1d(np.asarray(energy_per_u_mev, dtype=np.float64))
    if not np.all(np.isfinite(e)):
        raise ValueError("energies must be finite")
    if np.any(e < E_MIN_PER_U_MEV):
        raise ValueError(f"energy below the {E_MIN_PER_U_MEV} MeV/u domain limit of the model")
    beta2, bg2, tmax = kinematics(e, projectile)
    beta = np.sqrt(beta2)
    if effective_charge and projectile.z >= 2:
        q = pierce_blann_charge(beta, projectile)
    else:
        q = np.full_like(e, float(projectile.z))
    i_ev = material.mean_excitation_eV
    i_mev = i_ev * 1e-6
    bracket = (
        0.5 * np.log(2.0 * ELECTRON_MASS_MEV * bg2 * tmax / i_mev**2)
        - beta2
        - 0.5 * density_effect(bg2, material)
    )
    if shell:
        bracket = bracket - shell_correction_over_z(bg2, material)
    if barkas:
        bracket = bracket + barkas_term(beta2, q, material)
    if bloch:
        bracket = bracket + bloch_term(beta2, q)
    if mott and projectile.z >= 2:
        bracket = bracket + 0.5 * math.pi * ALPHA_FS * beta * q
    out: NDArray[np.float64] = K_MEV_CM2_MOL * q**2 * material.z_over_a / beta2 * bracket
    return out


@dataclass(eq=False)
class StoppingTable:
    """Electronic stopping power and CSDA range of one projectile in one material.

    Arrays (float64, one entry per grid point): ``energy_per_u`` [MeV/u], log-spaced;
    ``s_el_mass`` [MeV cm2/g]; ``s_el_linear`` [MeV/mm] = S rho / 10;
    ``csda_range_g_cm2`` [g/cm2]; ``range_mm`` [mm]. ``metadata`` records source, I-value,
    options and the start-range approximation. Interpolation is log-log (piecewise linear
    in ln E, ln S and ln R); range is monotone increasing so it can be inverted exactly.
    """

    projectile: Projectile
    material: Material
    energy_per_u: NDArray[np.float64]
    s_el_mass: NDArray[np.float64]
    s_el_linear: NDArray[np.float64]
    csda_range_g_cm2: NDArray[np.float64]
    range_mm: NDArray[np.float64]
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        rho = self.material.density_g_cm3
        names = ("energy_per_u", "s_el_mass", "s_el_linear", "csda_range_g_cm2", "range_mm")
        arrays = [getattr(self, n) for n in names]
        for name, arr in zip(names, arrays, strict=True):
            if arr.ndim != 1:
                raise ValueError(f"{name} must be one-dimensional")
            if arr.size != arrays[0].size:
                raise ValueError(f"{name} has length {arr.size}, expected {arrays[0].size}")
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"{name} contains NaN or infinity")
        if self.energy_per_u.size < 2:
            raise ValueError("table needs at least two grid points")
        if np.any(self.energy_per_u <= 0.0) or np.any(np.diff(self.energy_per_u) <= 0.0):
            raise ValueError("energy_per_u must be positive and strictly increasing")
        if np.any(self.s_el_mass <= 0.0):
            raise ValueError("s_el_mass must be positive")
        if np.any(self.csda_range_g_cm2 <= 0.0) or np.any(np.diff(self.csda_range_g_cm2) <= 0.0):
            raise ValueError("csda_range_g_cm2 must be strictly positive and strictly increasing")
        if np.any(self.range_mm <= 0.0) or np.any(np.diff(self.range_mm) <= 0.0):
            raise ValueError("range_mm must be strictly positive and strictly increasing")
        if not np.allclose(self.s_el_linear, self.s_el_mass * rho / 10.0, rtol=1e-12, atol=0.0):
            raise ValueError("s_el_linear is inconsistent with s_el_mass * density / 10")
        if not np.allclose(self.range_mm, self.csda_range_g_cm2 / rho * 10.0, rtol=1e-12, atol=0.0):
            raise ValueError("range_mm is inconsistent with csda_range_g_cm2 / density * 10")

    def _check(self, e: NDArray[np.float64]) -> None:
        if not np.all(np.isfinite(e)):
            raise ValueError("energies must be finite")
        if np.any(e < self.energy_per_u[0] * (1 - 1e-12)) or np.any(
            e > self.energy_per_u[-1] * (1 + 1e-12)
        ):
            raise ValueError(
                f"energy outside table [{self.energy_per_u[0]}, {self.energy_per_u[-1]}] MeV/u"
            )

    def stopping_at(self, energy_per_u_mev: ArrayLike) -> NDArray[np.float64]:
        """Mass electronic stopping power [MeV cm2/g] at energy per nucleon [MeV/u]."""
        e = np.asarray(energy_per_u_mev, dtype=np.float64)
        self._check(e)
        return np.exp(np.interp(np.log(e), np.log(self.energy_per_u), np.log(self.s_el_mass)))

    def range_at(self, energy_per_u_mev: ArrayLike) -> NDArray[np.float64]:
        """CSDA range [g/cm2] at energy per nucleon [MeV/u] (log-log interpolation)."""
        e = np.asarray(energy_per_u_mev, dtype=np.float64)
        self._check(e)
        return np.exp(
            np.interp(np.log(e), np.log(self.energy_per_u), np.log(self.csda_range_g_cm2))
        )

    def energy_from_range(self, range_g_cm2: ArrayLike) -> NDArray[np.float64]:
        """Energy per nucleon [MeV/u] for a CSDA range [g/cm2]; exact inverse of ``range_at``."""
        r = np.asarray(range_g_cm2, dtype=np.float64)
        if not np.all(np.isfinite(r)):
            raise ValueError("ranges must be finite")
        lo, hi = self.csda_range_g_cm2[0], self.csda_range_g_cm2[-1]
        if np.any(r < lo * (1 - 1e-12)) or np.any(r > hi * (1 + 1e-12)):
            raise ValueError(f"range outside table [{lo}, {hi}] g/cm2")
        return np.exp(
            np.interp(np.log(r), np.log(self.csda_range_g_cm2), np.log(self.energy_per_u))
        )


def log_grid(e_min: float, e_max: float, points_per_decade: int) -> NDArray[np.float64]:
    """Log-spaced energy grid [MeV/u] from ``e_min`` to ``e_max`` (inclusive)."""
    n = max(2, int(math.ceil(points_per_decade * math.log10(e_max / e_min))) + 1)
    return np.geomspace(e_min, e_max, n)


def build_table(
    projectile: Projectile,
    material: Material,
    energy_per_u: NDArray[np.float64],
    s_el_mass: NDArray[np.float64],
    start_range_g_cm2: float,
    metadata: dict[str, Any],
) -> StoppingTable:
    """Assemble a :class:`StoppingTable` from stopping powers on a log grid.

    ``R(E) = R(E_min) + int dE_total / S`` with ``dE_total = a dE_u``, evaluated by the
    trapezoid rule of ``E_total/S`` over ln E on the given grid; ``start_range_g_cm2`` is
    the range at the first grid energy.
    """
    e = np.asarray(energy_per_u, dtype=np.float64)
    s = np.asarray(s_el_mass, dtype=np.float64)
    if not (np.all(np.isfinite(e)) and np.all(np.isfinite(s)) and math.isfinite(start_range_g_cm2)):
        raise ValueError("energies, stopping powers and start range must be finite")
    if np.any(np.diff(e) <= 0.0) or np.any(s <= 0.0) or e[0] <= 0.0 or start_range_g_cm2 <= 0.0:
        raise ValueError(
            "energies must be positive and increase; stopping powers and start range positive"
        )
    f = projectile.a * e / s  # dR/dlnE [g/cm2]
    dln = np.diff(np.log(e))
    r = start_range_g_cm2 + np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * dln)))
    rho = material.density_g_cm3
    if not np.all(np.isfinite(r)):
        raise ValueError("computed ranges are not finite")
    return StoppingTable(
        projectile=projectile,
        material=material,
        energy_per_u=e,
        s_el_mass=s,
        s_el_linear=s * rho / 10.0,
        csda_range_g_cm2=r,
        range_mm=r / rho * 10.0,
        metadata=metadata,
    )


class StoppingSource(Protocol):
    """Provider of :class:`StoppingTable` objects for a material and projectile."""

    @property
    def name(self) -> str:
        """Short source identifier recorded in results."""
        ...

    def table(self, material: Material, projectile: Projectile) -> StoppingTable:
        """Return the stopping table of ``projectile`` in ``material``."""
        ...


@dataclass(frozen=True)
class BetheStoppingSource:
    """Analytic Bethe-Bloch source (construction data: none; see module docstring).

    ``e_min_per_u`` and ``e_max_per_u`` in MeV/u (e_min >= 1), ``points_per_decade`` grid
    density. The range at ``e_min_per_u`` is approximated by ``a E_min / S(E_min)``, the
    integral of dE/S below E_min for a constant stopping power equal to S(E_min). This is
    an approximation, not a bound: the electronic stopping power rises below E_min only
    down to its maximum (around 0.1 MeV/u for protons in water) and falls again at lower
    velocity, so the true residual range below 1 MeV/u can be smaller or larger. The offset
    is of the order of 0.02 mm of water for protons and is recorded in the table metadata;
    transport deposits the residual energy locally below E_cut, so it never enters a
    transported range.
    """

    options: BetheOptions = field(default_factory=BetheOptions)
    e_min_per_u: float = 1.0
    e_max_per_u: float = 500.0
    points_per_decade: int = 200
    name: str = "bethe"

    def table(self, material: Material, projectile: Projectile) -> StoppingTable:
        """Build the table (float64) for ``projectile`` in ``material``."""
        o = self.options
        e = log_grid(self.e_min_per_u, self.e_max_per_u, self.points_per_decade)
        s = bethe_mass_stopping(
            e,
            projectile,
            material,
            shell=o.shell,
            barkas=o.barkas,
            bloch=o.bloch,
            mott=o.mott,
            effective_charge=o.effective_charge,
        )
        r0 = projectile.a * e[0] / s[0]
        meta = {
            "source": "bethe",
            "I_eV": material.mean_excitation_eV,
            "material": material.name,
            "material_source": material.source,
            "material_density_g_cm3": material.density_g_cm3,
            "material_mass_fractions": dict(material.mass_fractions),
            "projectile": projectile.name,
            "options": {
                "shell": o.shell,
                "barkas": o.barkas,
                "bloch": o.bloch,
                "mott": o.mott,
                "effective_charge": o.effective_charge,
            },
            "start_range": "constant-S approximation a*E_min/S(E_min)",
        }
        return build_table(projectile, material, e, s, r0, meta)


_WATER_FRACTIONS = {"H": 0.111894, "O": 0.888106}


def _require_nist_water(material: Material) -> None:
    """Raise ``ValueError`` unless ``material`` is structurally liquid water (H2O, 1 g/cm3).

    The check uses composition (only H and O, mass fractions within 1e-4 of water's) and
    density (within 1e-6 g/cm3 of 1), never the material name. The NIST table is for
    I = 75 eV; the requested I is recorded in the table metadata but not used.
    """
    fr = material.mass_fractions
    ok = (
        set(fr) == set(_WATER_FRACTIONS)
        and all(abs(fr[k] - v) <= 1e-4 for k, v in _WATER_FRACTIONS.items())
        and abs(material.density_g_cm3 - 1.0) <= 1e-6
    )
    if not ok:
        raise ValueError(
            f"the NIST STAR source describes liquid water only (H and O with mass fractions "
            f"{_WATER_FRACTIONS} within 1e-4, density 1.0 g/cm3 within 1e-6); material "
            f"{material.name!r} does not match. Use BetheStoppingSource for other materials."
        )


_WATER_DATASETS = {"proton": "nist-pstar-water-2005", "alpha": "nist-astar-water-2005"}


def _require_water_table(star: StarTable, projectile: Projectile) -> None:
    """Raise ``ValueError`` unless ``star`` is the NIST liquid-water table for ``projectile``.

    A table loaded from the verified cache must carry the pinned dataset id of the water
    table of that projectile; an unverified table (``dataset_id`` None) must name liquid
    water in its parsed header material label (case-insensitive).
    """
    if star.dataset_id is not None:
        if star.dataset_id != _WATER_DATASETS.get(projectile.name):
            raise ValueError(
                f"dataset {star.dataset_id!r} is not the NIST liquid-water table for "
                f"{projectile.name}"
            )
        return
    label = star.material.upper()
    if "WATER" not in label or "LIQUID" not in label:
        raise ValueError(
            f"STAR table material {star.material!r} is not liquid water; "
            "refusing to use it as water"
        )


@dataclass(frozen=True)
class NistStarStoppingSource:
    """Tabulated source from a NIST PSTAR (protons) or ASTAR (alpha) liquid-water table.

    The table's electronic stopping power is interpolated log-log onto the grid
    (energy per nucleon = table energy / a); the range starts from the NIST CSDA range
    at ``e_min_per_u`` and then follows the same integral as the analytic source. When
    used, the NIST table is construction data of the run (decision 0038).
    """

    star_table: StarTable
    e_min_per_u: float = 1.0
    e_max_per_u: float = 500.0
    points_per_decade: int = 200
    name: str = "nist-star"

    def table(self, material: Material, projectile: Projectile) -> StoppingTable:
        """Build the table for ``projectile`` in liquid water from the STAR data."""
        st = self.star_table
        expected = {"PSTAR": "proton", "ASTAR": "alpha"}.get(st.program)
        if expected != projectile.name:
            raise ValueError(f"{st.program} table cannot describe projectile {projectile.name}")
        _require_water_table(st, projectile)
        _require_nist_water(material)
        e_u = st.energy_mev / projectile.a
        e_max = min(self.e_max_per_u, float(e_u[-1]))
        if self.e_min_per_u < e_u[0]:
            raise ValueError("requested e_min below the table")
        e = log_grid(self.e_min_per_u, e_max, self.points_per_decade)
        s = np.exp(np.interp(np.log(e), np.log(e_u), np.log(st.s_electronic)))
        r0 = float(np.exp(np.interp(math.log(e[0]), np.log(e_u), np.log(st.csda_range))))
        meta = {
            "source": "nist-star",
            "I_eV": 75.0,
            "projectile": projectile.name,
            "program": st.program,
            **st.provenance(),
            "requested_material": {
                "name": material.name,
                "density_g_cm3": material.density_g_cm3,
                "mass_fractions": dict(material.mass_fractions),
                "I_eV": material.mean_excitation_eV,
            },
            "effective_material": {
                "name": "NIST liquid water (matno 276)",
                "density_g_cm3": 1.0,
                "mass_fractions": dict(_WATER_FRACTIONS),
                "I_eV": 75.0,
            },
            "start_range": "NIST CSDA range at e_min (log-log interpolation)",
        }
        return build_table(projectile, material, e, s, r0, meta)


def default_source() -> StoppingSource:
    """Return the default analytic source (decision 0038)."""
    return BetheStoppingSource()
