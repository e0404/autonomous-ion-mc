# Electronic stopping power

Implemented in `ionmc.physics.stopping` (model and tables), `ionmc.materials` and
`ionmc.physics.projectiles`. The design decision is
[0038](../generated/decisions/0038-electronic-stopping-power-and-data-roles.md). The
module docstring of `stopping.py` is the authoritative statement of each term.

## Model

The mass electronic stopping power of an ion with charge number z in a material is

    S/rho = K q^2 (Z/A) beta^-2 [ L0 + q L1 + L2 + M ]      [MeV cm^2/g]

with K = 0.307075 MeV cm^2/mol, (Z/A) in mol/g and q the effective charge. All
arithmetic is float64. Energies are kinetic energy per nucleon (MeV/u); the total
kinetic energy of the ion is E_u times the mass number.

| Term | Implementation |
|---|---|
| L0 | `1/2 ln(2 m_e c^2 beta^2 gamma^2 Tmax / I^2) - beta^2 - delta/2 - C/Z`, with the exact `Tmax = 2 m_e c^2 beta^2 gamma^2 / (1 + 2 gamma m_e/M + (m_e/M)^2)` |
| delta | Sternheimer density effect (Sternheimer, Berger, Seltzer, At. Data Nucl. Data Tables 30 (1984) 261); zero if the material has no parameters. Parameters are provided for water (x0 0.2400, x1 2.8004, Cbar 3.5017, a 0.09116, m 3.4773) |
| C/Z | Shell correction in the Bichsel parameterisation of ICRU Report 49 as written in the PDG review (polynomial in eta = beta gamma, I in eV, valid for eta >= 0.13, about 8 MeV/u for protons). Below that limit the value at eta = 0.13 is held. For compounds, C(I_i)/Z_i of each element (elemental I) is weighted by its electron fraction |
| q L1 | Barkas term in the Ashley-Ritchie-Brandt form (Phys. Rev. B 5 (1972) 2393; ICRU 49) as implemented in Geant4 `G4EmCorrections::BarkasCorrection` (v11.4.2); the tabulated function F(W) is the 47-point table of that source, interpolated linearly |
| L2 | Bloch term `-y^2 sum_{n>=1} 1/(n (n^2 + y^2))`, `y = q alpha / beta` |
| M | Mott term `1/2 pi alpha beta q` for projectiles with z >= 2 (as in Geant4) |
| q | Pierce-Blann effective charge `z (1 - exp(-125 beta z^(-2/3)))` for z >= 2; q = z for z = 1 |

Because Barkas, Bloch and Mott depend on the actual z, helium, carbon and oxygen
stopping powers are not z^2-scaled proton stopping powers. Each correction can be
switched off through `BetheOptions` (or the keyword arguments of
`bethe_mass_stopping`); the options are stored in the table metadata.

## Domain

Energies below 1 MeV/u raise `ValueError`: the model does not return stopping
powers there, and transport is expected to deposit the residual energy locally.
Between 1 MeV/u and about 8 MeV/u (protons) the shell correction is held at its
boundary value, so the model is less accurate in that interval; the range effect is
confined to the last part of the track.

## Tables

`BetheStoppingSource(options, e_min_per_u, e_max_per_u, points_per_decade).table(material,
projectile)` returns a `StoppingTable` on a log-spaced grid (defaults 1 to 500 MeV/u,
200 points per decade) with `s_el_mass` [MeV cm^2/g], `s_el_linear` [MeV/mm] =
S rho / 10, `csda_range_g_cm2` and `range_mm`. The CSDA range is the trapezoid integral
of `a E_u / S` over ln E_u starting from `a E_min / S(E_min)`, which is the integral of
dE/S below E_min for a constant stopping power. Since the true stopping power is larger
at lower energy, this start value is an upper bound of the range below E_min. Stopping
power, range and the inverse `energy_from_range` use one log-log piecewise-linear
interpolant, so `energy_from_range(range_at(E))` returns E to rounding error.
`NistStarStoppingSource` produces the same structure from a cached PSTAR (protons) or
ASTAR (alpha) liquid-water table, starting the range from the NIST CSDA range at E_min.
Both implement the `StoppingSource` protocol; the table metadata records the source and
the I value.

## Mean excitation energy of water

`ionmc.materials.WATER` uses I = 78 eV (ICRU Report 90, 2016). `water(I_eV)` builds a
variant, for example `water(75.0)` for the ICRU 49 / NIST PSTAR/ASTAR value. Material I
is configurable per material; if it is `None` the Bragg-additivity value is used.
Comparisons with other engines must state each engine's I.

## Data roles

| Dataset | Role |
|---|---|
| NIST PSTAR / ASTAR liquid water (I = 75 eV) | evaluation of the analytic model at I = 75 eV; construction data when `NistStarStoppingSource` is selected |
| ICRU 90 water arrays (from the Geant4 source) | evaluation of the analytic model at I = 78 eV |

PSTAR, ASTAR and ICRU 90 share the Bethe-theory lineage above about 1 MeV/u; agreement
with them is tabulated, not measured, evidence. A run that uses the NIST tables as its
stopping source cannot be evaluated against the same tables.

The energy grid of the ICRU 90 alpha array is interpreted as total alpha kinetic
energy (same span as ASTAR); `compare_nist.py` reports the ratio of the two datasets
on their shared energies so that this reading is checked on every run.

## Comparison script

After fetching the three datasets (see [data layer](../architecture/data-layer.md)):

```bash
uv run python validation/scripts/stopping/compare_nist.py --output-dir OUT_DIR
```

It reads only the cache (no network unless `--online`) and writes
`stopping_comparison.json` containing relative deviations on the reference energy grids
for E >= 2 MeV/u (Bethe at I = 75 eV against PSTAR/ASTAR; Bethe at I = 78 eV against
the ICRU 90 arrays), CSDA range differences at 100/150/200/250 MeV protons and
100/150/200 MeV/u alpha particles, and the range shift between I = 78 and 75 eV.
Results are produced by running the script; this page does not quote them.

## Sources

- ICRU Report 49 (1993), Stopping Powers and Ranges for Protons and Alpha Particles.
- ICRU Report 90 (2016), Key Data for Ionizing-Radiation Dosimetry (J. ICRU 14(1)).
- Particle Data Group, Passage of particles through matter.
- M.J. Berger, J.S. Coursey, M.A. Zucker, J. Chang (2005), ESTAR, PSTAR, and ASTAR
  (version 2.0.1), NIST, https://doi.org/10.18434/T4NC7P.
- Geant4 11.4.2: `G4EmCorrections.cc`, `G4ICRU90StoppingData.cc`, `G4NistMaterialBuilder.cc`.
