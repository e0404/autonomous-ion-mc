# Physics

This section documents implemented transport and interaction physics, model
assumptions, and applicable validity ranges. Implemented so far: the analytical
electronic stopping power and CSDA range (task `DEV-002`), a tabulated
stopping-power layer from external data (task `DEV-003`), continuous-slowing-down
proton transport in homogeneous water with energy-loss straggling and multiple
Coulomb scattering (milestone V1, tasks `DEV-004`/`DEV-005`/`DEV-006`),
nonelastic nuclear attenuation of primaries (task `DEV-007`), and secondary
charged-particle transport (task `DEV-008`, closing milestone V2).

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

## Tabulated stopping power (data layer)

Module: ``ionmc.physics.tabulated``; API:
``ionmc.tabulated_stopping_power.TabulatedStoppingPower``; data layer:
``ionmc.data``; decisions `0007`, `0008`.

The same interface as the analytical layer (`mass_stopping_power`,
`csda_range`, same units) is served from an external table. The reference
table is the proton stopping power of liquid water from the MCsquare
repository (``Materials/Water/PSTAR_Stop_Pow.dat``, Apache-2.0, UCLouvain;
commit ``211eefe6``): NIST PSTAR (SRD 124) *total* mass stopping power on a
0.5 MeV grid, 0-400 MeV. Nuclear stopping is at most 0.1 % of the total above
0.5 MeV, so the table is used as the electronic stopping power at therapeutic
energies.

| aspect | choice | reason |
|---|---|---|
| interpolation | monotone cubic Hermite (PCHIP), slopes precomputed at load | reproduces the data generator; shape-preserving; no transcendentals; error against the model ≲ 1e-5 above 1 MeV |
| segment lookup | fixed-count bisection (``ceil(log2 n)`` steps) | grid-agnostic; runs identically on the Python, numpy and Warp paths |
| CSDA range | cumulative range precomputed per segment by 4-point Gauss-Legendre of ``1/S``; a query integrates the partial top segment with the same rule | continuity across grid nodes; range from the 0.5 MeV floor, residual below it added by the caller |
| out of range | raises unless ``allow_extrapolation=True`` | avoids MCsquare's out-of-bounds read |

The layer separates an **acquire** phase (needs network; verifies a SHA-256
and writes a manifest) from an offline **use** phase (cache-only, fails with
an actionable message when data are missing), because the host runner has no
network (decision `0007`). The cache location is configurable
(``IONMC_CACHE_DIR`` or ``$XDG_CACHE_HOME/ionmc``). A second table
(``mcsquare-g4-water``, Geant4-derived) is registered as a selectable
alternative source; it lies 0.4-0.5 % below the PSTAR table, matching the
ICRU 90 versus ICRU 49 analytic offset.

### Validity and accuracy

The tabulated layer agrees with the analytic model (decision `0006`) to within
0.25 % over 10-400 MeV and reproduces NIST PSTAR CSDA ranges at
100/150/200/250 MeV to within 1.8e-4 (with the 8.9 µm residual range below
0.5 MeV added). It is defined over the tabulated energy range 0.5-400 MeV.

## Continuous-slowing-down transport (Stage 1, longitudinal)

Modules: ``ionmc.transport`` (``TransportEngine``, ``PencilBeamSource``,
``WaterSlab``, ``DepthDoseGrid``); shared-source step physics
``ionmc.physics.transport``; decision `0009`.

A monoenergetic proton pencil beam is transported through homogeneous water by
the continuous-slowing-down approximation with **no straggling and no
scattering yet**: one history per thread, an in-kernel step loop with a
midpoint (RK2) energy-loss integration, and energy deposited along the axis
into a 1-D integral-depth-dose grid by atomics. The step is limited by a
fractional energy loss and by the depth-bin boundary; below a 0.5 MeV cutoff
the residual energy is deposited locally. The stopping power comes from the
tabulated layer (task DEV-003).

The reference Python path is a scalar per-history loop over the *same*
shared-source step functions the Warp CPU/CUDA kernel runs. Energy is conserved
exactly (float64) and the distal 80 % depth (R80) reproduces the tabulated CSDA
range to within 0.15 % at 100-200 MeV. Because the distal edge is a
near-discontinuity without straggling, cross-backend agreement is measured on
the *cumulative* depth dose (reference vs Warp ≤ 1e-4 of the total), not
bin-by-bin.

### Energy-loss straggling and the Bragg peak (task DEV-005, decision 0010)

Each transport step adds a Gaussian fluctuation to the mean energy loss,
``dE = clamp(dE_mean + xi sigma, 0, E)`` with the Bohr variance
``sigma^2 = K m_e c^2 (Z/A) rho z^2 f(beta) dx`` (relativistic factor
``f(beta) = (1 - beta^2/2)/(1 - beta^2)``) and ``xi`` a standard normal from
the per-history counter-based stream (drawn in the execution layer, not the
shared source). This gives a realistic Bragg peak. Transport now steps by a
physics-sized step (fractional energy loss) and distributes the deposited
energy across the depth bins the step spans, so the mean range is
grid-independent and unbiased.

Validated on the **range straggling** (std of stopping depths), which isolates
straggling from the peak shape: the mean range reproduces the CSDA range to
<0.1 %, and sigma_R matches the model's own analytic Bohr integral to <1 % and
Bortfeld's ``0.012 R^0.935`` to +5-7 % (the pure-Bohr value is slightly high
because it omits the electron-binding reduction), staying at 1.0-1.1 % of the
range over 100-200 MeV. Statistical uncertainty is estimated by batches
(``TransportEngine.run_batched``); two independent seeds are consistent within
their combined standard error (decision 0010).

### Multiple Coulomb scattering and lateral spread (task DEV-006, decision 0011)

Each step deflects the proton by a small angle sampled from the Highland /
Lynch-Dahl projected RMS ``theta0 = (13.6 MeV / (beta c p)) z sqrt(x/X0)`` (the
scattering-power form, no log term; water radiation length ``X0 = 36.08
g/cm^2``), applied with a **random hinge** (the full-step deflection at a
uniform random point along the step, which gives the correct lateral
displacement) as a two-plane tilt of the direction. Transport is now 3-D and
energy is scored into a 2-D (depth, lateral-``x``) grid; the depth dose is the
lateral marginal and ``sigma_x(z)`` is the lateral second moment per depth
slice.

Validated against the **Fermi-Eyges** analytic lateral spread built from the
same scattering power (``ionmc.physics.fermi_eyges``): the Monte Carlo
``sigma_x(z)`` agrees to within 1 % at 0.5R and 0.8R for 150 and 200 MeV, and
the lateral ``sigma_x`` at 0.8R (2.45 mm at 150 MeV, 3.9 mm at 200 MeV) matches
published values within ~2 %. Multiple scattering shortens the projected range
by the detour factor (< 0.1 %). A single Gaussian central-scattering model is
used (no single-scattering tail), which is standard for fast therapy Monte
Carlo and adequate for the lateral-``sigma`` core.

This closes the physics of milestone **V1** (proton transport in homogeneous
water).

## Nonelastic nuclear attenuation (Stage 2, task DEV-007, decision 0012)

Module: ``ionmc.physics.nuclear``; wired into ``TransportEngine`` and the Warp
depth-dose kernel behind a ``nuclear`` flag.

Real protons are removed along the track by **nonelastic nuclear reactions**,
so about 20 % (150 MeV) to 27 % (200 MeV) never reach the Bragg peak. Each
alive step draws one uniform and, with probability ``Sigma(E) x s`` (the
thin-step macroscopic nonelastic rate over the step length ``s``), removes the
primary: a local fraction ``f_local = 0.30`` of its energy is deposited at the
vertex (short-range recoils and fragments) and the remainder is booked to an
**escaping/deferred channel** for the secondary transport of Stage 2 (DEV-008).
The macroscopic cross section is oxygen-only, ``Sigma = n_O sigma_nonel(E)``;
the hydrogen channel is proton-proton *elastic* scattering, which deflects
rather than removes the primary and is excluded here. ``sigma_nonel(E)`` on
oxygen is an analytic parameterization of the ICRU-63 shape (7 MeV threshold,
~550 mb near 20 MeV, ~340-400 mb plateau over 100-250 MeV), written in the
shared-source math namespace so the Warp kernel and the float64 reference
evaluate it identically.

The nuclear uniform is drawn after the straggling normal and unconditionally
per alive step, so the reference and Warp counter streams stay bit-aligned and
the two backends remove the identical primary set (decision 0001). With
``nuclear=False`` no extra draw is made and the DEV-004/005/006 results are
reproduced exactly. Energy is fully accounted: ``deposited + escaped =
energy_in``. Validated on the primary survival to the Bragg peak (0.813 /
0.716 at 150 / 200 MeV, within the published tolerances) and on the analytic
reaction fraction (see the validation index).

The **absolute** peak-to-entrance ratio still reads high until the escaping
secondary protons are transported (below).

## Secondary charged-particle transport (Stage 2, task DEV-008, decision 0013)

Module: ``ionmc.physics.secondaries`` (host-side sampling), wired into
``TransportEngine`` as a second transport pass.

The escaping energy booked by DEV-007 is not all lost: about half is carried by
**secondary protons** that travel on and deposit a broad low-level dose (the
nuclear plateau/halo). Each nonelastic reaction of residual energy ``E`` is now
partitioned into local heavy fragments (``f_heavy = 0.12``, deposited at the
vertex), transported secondary protons (``f_p = 0.50``), and truly escaping
neutrons, gammas and binding energy (``f_esc = 0.38``, booked as escaping). The
secondary-proton multiplicity is ``nu_p(E) = 0.5 + 0.004 E`` (Poisson); their
energies are sampled from a two-component spectrum (an evaporation Maxwellian
``E exp(-E/T)`` with ``T = 2`` MeV mixed with a forward cascade component uniform
on ``[10 MeV, E]``) and **renormalised so their per-reaction sum equals
``f_p E`` exactly**, which keeps the energy budget closed. Because the
multiplicity is low and the energies are renormalised, this budget-closing
renormalisation dominates the effective spectrum (a single-secondary reaction
emits one proton of ``f_p E`` regardless of the sampled shape); the evaporation/
cascade parameters are a second-order influence. This is an intentional depth-
dose surrogate; a faithful differential spectrum is deferred with the tabulated-
data follow-up (decision 0013). Secondaries are emitted forward from the vertex
and transported by the same CSDA + straggling proton engine (their own nuclear
removal off, a documented sub-percent simplification).

Because each primary reacts at most once, the transport drivers emit one
reaction record per history (vertex depth, residual energy); the secondaries are
generated **host-side** from those records and transported in a second pass.
Since DEV-007 already makes the reference and Warp paths remove the same primary
set, and host-side generation is deterministic (counter RNG keyed by the history
index), the two backends produce the same secondary set and dose in practice,
within the decision-0001 tolerances (not strictly bit-identical: the reaction
energy fed to the sampler is float32 on the Warp path and float64 on the
reference path, so a rare boundary case could shift a count; the cumulative
depth-dose agreement is the real gate). Validated: secondary protons contribute
~1-2 % of the local dose at entrance and ~4-7 % of the total dose at 150/200 MeV;
their dose *fraction* rises to a ~5-10 % plateau proximal to the peak, then
collapses at the sharp Bragg peak (~0.1 %) where the primary dose dominates. The
energy budget ``deposited + escaped = energy_in`` stays exact for any geometry.
This **closes milestone V2**.

Deferred with quantitative justification (Paganetti 2002): explicit deuteron/
triton/alpha and recoil transport (deposited locally, < 0.1 % of dose), neutron
and prompt-gamma transport (dropped as escaping, < 0.05 %), tertiary reactions,
and the ICRU-63/TENDL tabulated double-differential path. Not yet: nuclear
removal and the lateral halo on the 3-D scattering path, heterogeneous voxel
geometry (Stage 3), and treatment-planning scoring and influence matrices
(Stage 4).
