# Physics

This section documents implemented transport and interaction physics, model
assumptions and validity ranges. Transport itself is not yet implemented;
the following data and physics layers are.

## Species and materials

`ionmc.species` defines the transported nuclear species (p, d, t, ³He, ⁴He,
⁶⁻⁷Li, ⁷⁻¹⁰Be, ¹⁰⁻¹¹B, ¹⁰⁻¹²C, ¹³⁻¹⁵N, ¹⁴⁻¹⁶O) with nuclear masses from
AME2020; ion kinetic energies are handled per nucleon (MeV/u) in the physics
layer. `ionmc.materials` defines elemental compositions (mass fractions),
densities and mean excitation energies: elements from ICRU 37/49, water
(78 eV), air (85.7 eV) and graphite (81 eV) from ICRU Report 90, tissues
from the ICRU 44/46 substitutes as tabulated by Geant4's NIST database, and
Bragg additivity otherwise.

## Electronic stopping power

Two layers feed one table type (decision 0038):

1. **Analytic layer** (`ionmc.physics.stopping`): Bethe stopping number with
   shell correction (ICRU 37 parameterization, valid for βγ ≥ 0.13), Barkas
   (approximate Ashley–Ritchie–Brandt form), Bloch (exact series), a
   leading-order Mott term, the Sternheimer density effect for water, and the
   Pierce–Blann effective charge. Validated against the ICRU 90 water tables
   at equal mean excitation energy: protons within 0.12 % (16–500 MeV),
   helium within 0.25 % (16–250 MeV/u), carbon within 0.6 % above 16 MeV/u
   and 0.3 % above 30 MeV/u; against PSTAR/ASTAR (ICRU 49, I = 75 eV) within
   0.1 % (protons, 8–500 MeV) and 0.25 % (helium, 8–250 MeV/u) when
   evaluated at I = 75 eV. Below ≈ 10 MeV/u the layer is *not* accurate for
   heavy ions (carbon −1.5 % at 10 MeV/u, −9 % at 2 MeV/u) and is not used
   there.
2. **Tabulated layer** (`ionmc.data`): the ICRU 90 liquid-water tables for
   protons, alpha particles and carbon ions (I = 78 eV; shipped with the
   package, ≈ 8 kB, provenance in the file header) and NIST PSTAR/ASTAR
   (ICRU 49) for other materials with verified NIST material codes,
   downloaded on demand into the provenance cache.

`ionmc.physics.tables.build_stopping_table` combines them on a logarithmic
grid (1 keV/u – 1 GeV/u, 120 points per decade): the tabulated source below
the blend window, the analytic layer above it, and a smooth cosine blend
inside (8–16 MeV/u for z ≤ 2, 16–32 MeV/u for z ≥ 3). Ions without their
own table use the proton table scaled by the squared effective-charge ratio
at equal velocity (within ±2 % of the ICRU 90 helium and carbon tables above
1 MeV/u). Materials without any tabulated data use the ICRU 90 water shape
scaled by the analytic ratio at the window edge, which the table provenance
records as an approximation. The CSDA range is integrated on the same grid;
constructed water ranges reproduce the ICRU 90 range column within 0.4 %
(e.g. 290 MeV/u carbon: 16.3 g/cm²).

Mean excitation energy lineage matters at the 0.5 % level in range: the
default tables use ICRU 90 values (78 eV for water), so ranges are ≈ 0.5 %
longer than PSTAR/ASTAR (75 eV). Every table carries its I value, sources,
dataset identifiers/hashes and domains in `StoppingTable.provenance`.

## Condensed-history transport (electromagnetic)

Implemented on the reference Python backend and the Warp CPU/CUDA backends
(decisions 0041, 0042) with the shared step physics of
`ionmc.transport.step_physics`:

| Element | Model | Domain / limitation |
|---|---|---|
| Step limits | voxel face, `max_step_mm` (1 mm), 10 % CSDA energy loss | configurable; step-size independence tested |
| Mean energy loss | exact inversion of the log-log CSDA range table | any step length |
| Straggling | Gamma with CSDA mean and Bohr variance (relativistic factor, z_eff) | integrated quantities; not Landau-shaped per step |
| Multiple scattering | Gaussian core, Gottschalk differential Molière scattering power T_dM, z² scaling, random hinge | no single-scattering tail |
| Cutoff | 0.5 MeV/u, residual energy deposited locally and counted | — |
| Nuclear interactions | not implemented; requesting them raises `UnsupportedConfigurationError` | — |

Verified on the reference backend: energy conservation to 1e-9 with escaped
and cutoff energies accounted separately; 100 MeV protons in water R80 =
77.57 mm versus the ICRU 90 CSDA range 77.59 mm, unchanged between 1 mm/10 %
and 0.2 mm/2 % steps; scoring grids of different resolution or alignment
conserve the deposited energy.
