# Physics

This section documents implemented transport and interaction physics, model
assumptions, and applicable validity ranges. Implemented so far (task
`DEV-002`): the analytical electronic stopping power and CSDA range. No
transport exists yet.

## Electronic stopping power (analytical layer)

Module: ``ionmc.physics.stopping``; API: ``ionmc.stopping_power.AnalyticStoppingPower``;
decision: `0006`.

The mass electronic stopping power is the Bethe formula with corrections:

```
-(1/rho) dE/dx = K z^2 (Z/A) / beta^2 * L(beta)

L = 1/2 ln(2 m_e c^2 beta^2 gamma^2 T_max / I^2) - beta^2 - delta/2 - C/Z + z L1 + L2

T_max = 2 m_e c^2 beta^2 gamma^2 / (1 + 2 gamma m_e/M + (m_e/M)^2)
```

with `K = 4 pi N_A r_e^2 m_e c^2` (0.307075 MeV cm^2/mol from CODATA 2018
constants), `z` the projectile charge, `Z/A` the electrons per unit mass of
the target, `I` the mean excitation energy and `M` the projectile rest energy.

| term | model | source |
|---|---|---|
| density effect `delta` | Sternheimer–Berger–Seltzer parameterisation in `x = log10(beta gamma)` with per-material `Cbar, x0, x1, a, m, delta0` | SBS 1984 via Geant4 `G4DensityEffectData` |
| shell correction `C/Z` | Barkas–Berger empirical `C(I, eta) = A(eta) I^2 + B(eta) I^3` per target element with ICRU 37 elemental `I`, combined as `sum_k n_k C_k / N_e`; frozen below the velocity of an 8 MeV proton and tapered to zero at that of a 2 MeV proton | Leo eq. 2.33; Geant4 `G4IonisParamElm`, `ShellCorrectionSTD` |
| Barkas `z L1` | Ashley–Ritchie–Brandt, `L1 = 1.29 F(b/sqrt(X)) / (sqrt(Z X) X)`, `X = beta^2/(alpha^2 Z)`, tabulated `F` (47 points), element parameter `b`; electron-weighted over elements | Ashley & Ritchie 1972; ICRU 49; Geant4 `G4EmCorrections` |
| Bloch `L2` | exact series `-y^2 sum_j 1/(j (j^2 + y^2))`, `y = z alpha/beta`, 16 terms + tail | Geant4 `BlochCorrection` |

Every correction can be switched off individually (``Corrections``), which
serves the configurable-fidelity requirement and makes each term's magnitude
testable. For protons in water at 10 MeV the terms are: shell −0.78 %, Barkas
+0.30 %, Bloch −0.05 %, density 0 (the density effect is identically zero in
water below about 900 MeV).

### Mean excitation energy policy

`I` is a per-material, provenance-tagged parameter (``MeanExcitationEnergy``).
``WATER`` uses 75.0 eV (ICRU 37/49, the value underlying NIST PSTAR) so that
validation against PSTAR-derived tables is I-consistent; ``WATER_ICRU90``
uses 78 eV (ICRU 90). The ICRU 90 value lowers the stopping power by
0.4–1.1 % over 1–400 MeV. Density-effect parameter sets are bound to the `I`
they were fitted with. Elemental `I` values enter only the shell correction;
the compound value for liquid water is never derived by Bragg additivity
(which would give about 67 eV).

### Validity

Claimed accurate (within the 1 % criterion of decision `0006`) for proton
kinetic energies of 10–1000 MeV; measured deviation from PSTAR over
10–400 MeV is at most 0.075 %. Below 10 MeV the analytic layer degrades
(+0.9 % at 5 MeV, +3.5 % at 2 MeV, +5.2 % at 1 MeV against PSTAR) and is
documented as such; the tabulated layer (task `DEV-003`) is intended to
cover that region.

### Units

Energies in MeV, `I` converted to MeV inside the formula, stopping power in
MeV cm²/g (decision `0004`).

## CSDA range

``AnalyticStoppingPower.csda_range`` integrates `dE/S(E)` over `ln E`
by composite Simpson quadrature (200 intervals; quadrature error below
10⁻⁶) from a 1 MeV floor upward. The neglected residual range below 1 MeV is
about 0.0025 g/cm² for protons in water. Protons in water: 7.717 g/cm² at
100 MeV, 15.776 at 150 MeV, 25.962 at 200 MeV, 37.944 at 250 MeV.

## Not yet implemented

Energy-loss straggling, multiple Coulomb scattering, nuclear interactions,
tabulated stopping powers and material libraries, and all transport are
future stages of the [roadmap](../development/roadmap.md).
