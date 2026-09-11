# Physics

This section documents implemented transport and interaction physics, model
assumptions, and applicable validity ranges. Implemented so far: the analytical
electronic stopping power and CSDA range (task `DEV-002`), a tabulated
stopping-power layer from external data (task `DEV-003`), continuous-slowing-down
proton transport in homogeneous water with energy-loss straggling and multiple
Coulomb scattering (milestone V1, tasks `DEV-004`/`DEV-005`/`DEV-006`),
nonelastic nuclear attenuation of primaries (task `DEV-007`), secondary
charged-particle transport (task `DEV-008`, closing milestone V2), and 1-D
voxelized density heterogeneity (task `DEV-009`, opening milestone V3), and
per-voxel tissue materials via stopping-power ratios (task `DEV-010`), and
density-heterogeneous 3-D multiple scattering (task `DEV-011`), and non-water
materials on the 3-D scattering path (task `DEV-012`).

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
and the ICRU-63/TENDL tabulated double-differential path.

## Voxelized density heterogeneity (Stage 3, task DEV-009, decision 0014)

Module: ``ionmc.transport.geometry`` (``VoxelSlab``), wired into the depth-dose
engine and Warp kernel as a per-step density lookup.

Real treatment planning runs on a voxelized patient geometry where the dose
depends on the **water-equivalent thickness (WET)** the beam traverses, not the
geometric depth. A ``VoxelSlab`` gives the beam axis a 1-D stack of voxels, each
with its own mass density (built from ``(thickness, density)`` layers with exact
interfaces, or a uniform grid). The tabulated stopping power is a *mass*
stopping power, so the linear energy loss is ``S/rho * rho * dl``: the transport
looks up the local voxel density each step and limits the step to the voxel
boundary so the density is unambiguous, tracking a **non-decreasing** voxel index
(forward transport never backscatters). The nonelastic nuclear rate uses a
per-voxel oxygen number density that scales with the voxel mass density.

Consecutive equal-density voxels are merged, so a uniform slab (any voxel count)
and a homogeneous ``WaterSlab`` both collapse to a single voxel and reproduce the
homogeneous transport **bit-for-bit** (a validation gate). Validated: the CSDA
range scales as ``1/rho`` (R80 at ``R_water/rho`` to ~1e-4 on a 0.1 mm grid); a
dense layer shifts the Bragg peak proximally by exactly its extra water-
equivalent thickness (a 20 mm, 1.85 g/cm^3 layer shifts the peak 17 mm); the
energy budget stays exact with nuclear and secondary transport across a density
interface; and the reference and Warp paths agree across the interface. This
**opens milestone V3**.

## Per-voxel tissue materials (Stage 3, task DEV-010, decision 0015)

Module: ``ionmc.materials`` (tissue library) and
``ionmc.stopping_power.mass_stopping_power_ratio``, wired into the engine's
per-voxel arrays.

Different tissues stop protons differently *per gram* than water (through their
electron density ``<Z/A>`` and mean excitation energy ``I``), so density alone
(DEV-009) is not enough. An ICRU-44 tissue library is added (cortical bone,
adipose, soft tissue, skeletal muscle, lung, air), and each voxel is transported
as **water at its water-equivalent density** ``rho_we = SPR(material) x
rho_phys``. The stopping-power ratio ``SPR`` is the ratio of the analytic Bethe
mass stopping powers (material / water) at a 150 MeV reference energy, so the
water table and the DEV-009 per-voxel-density transport are reused **unchanged**;
the material physics enters entirely through the per-voxel density and nuclear
arrays the engine assembles at construction. Because ``SPR ~ <Z/A>_mat/
<Z/A>_wat`` up to a small ``I``-value term, transporting in this water-equivalent
frame (water table, water ``<Z/A>``) reproduces the material's linear straggling
to a sub-percent term. The nonelastic nuclear rate uses a **composition-scaled**
oxygen-equivalent number density ``sum_{Z>1} (w_i/A_i)(A_i/A_O)^(2/3)`` (geometric
``A^(2/3)`` scaling of the oxygen cross-section shape), which reduces to ``n_O``
for water and is ~2x larger in bone / ~3.4x in adipose, removing the large error
of oxygen-only removal in bony and fatty tissue.

Validated: the water-equivalent ratios ``SPR x rho`` fall in the published
Schneider/ICRU bands (cortical bone 1.70, adipose 0.97, muscle 1.04, lung 1.04);
a tissue slab puts R80 at ``R80_water / WER``; a water ``VoxelSlab`` reproduces
the homogeneous baseline bit-for-bit; energy is conserved across material
interfaces with nuclear and secondary transport; and the reference and Warp paths
agree. Radiation lengths (Tsai's formula) are populated per material for the
future heterogeneous scattering path.

Deferred: energy-dependent SPR / per-material stopping tables (cortical bone's
2.5 % energy dependence; soft tissues are already < 0.25 % with the 150 MeV
scalar), element-specific measured nonelastic cross sections (C, Ca), per-tissue
density-effect parameters, and non-water/heterogeneous-material multiple
scattering on the 3-D path (the water-only scattering path rejects both a
heterogeneous slab and a homogeneous non-water material with a
``NotImplementedError`` rather than silently giving a wrong Bragg depth).

## Density-heterogeneous 3-D scattering (Stage 3, task DEV-011, decision 0016)

The 3-D multiple-scattering path (DEV-006) is brought to **1-D voxelized density
heterogeneity** for water, reusing the DEV-009 per-voxel-density-by-depth
pattern: the scattering reference driver and Warp kernel look up the local water
density by depth ``pz``, limit the step so its depth advance ``s d_z`` stays
within the voxel, and advance a non-decreasing voxel index. Since the Highland
scattering power uses the areal thickness ``rho s / X0``, the same per-voxel
density that drives the energy loss also drives the lateral spread; for water
``X0`` and ``<Z/A>`` are constant, so no new physics is needed. ``run_scattering``
now accepts a density-heterogeneous water ``VoxelSlab`` and still rejects
non-water materials (whose MCS needs the physical density and material radiation
length, distinct from the water-equivalent density used for stopping — a later
task).

Validated: a uniform water ``VoxelSlab`` reproduces the homogeneous ``WaterSlab``
scattering bit-for-bit; the lateral ``sigma_x`` matches the Fermi-Eyges oracle at
a uniform density and a **piecewise-density** Fermi-Eyges variant across a water/
dense/water interface to < 1 %; energy is conserved; and the reference and Warp
paths agree (``sigma_x`` tightly, the depth dose within the 3-D float32 budget).

## Non-water materials on the 3-D scattering path (Stage 3, task DEV-012, decision 0017)

DEV-012 extends the scattering path to per-voxel non-water **materials**,
mirroring DEV-010 for the scattering kernel. The subtlety is that stopping and
scattering use *different* densities: the energy loss uses the per-voxel
water-equivalent density ``rho_we = SPR(material) rho_phys`` (with water's
``<Z/A>``), while the Highland scattering power uses the per-voxel **physical**
density and **material** radiation length, ``theta0^2 = (13.6/pv)^2 (rho_phys s /
X0_mat)``. The engine builds and merges a per-voxel physical density and radiation
length alongside the water-equivalent and oxygen-equivalent densities; the
scattering reference driver and Warp kernel look both up by depth. The
``is_water_only`` guard is lifted; only a known radiation length per voxel is
required. A uniform single-material slab still collapses to one voxel, so water
reproduces the prior result bit-for-bit.

Validated against a **material-aware** Fermi-Eyges oracle
(``lateral_sigma_x_material_mm``: energy vs depth from the integrated
water-equivalent thickness, scattering power from the physical density and
material radiation length): a homogeneous bone slab (R80 at ``R_water/WER``,
``sigma_x`` with bone's smaller X0) and a water/bone/water interface reproduce the
oracle to < 1 %, energy is conserved, and the backends agree. The depth-dose and
scattering paths now share the same per-voxel material capability.

## Arbitrary beam incidence via a beam frame (Stage 3, task DEV-013, decision 0018)

DEV-013 lets the pencil beam enter at an arbitrary position and unit
``direction`` and the whole scene (beam + phantom) be rigidly rotated, on the
scattering path. Transport runs in a **canonical beam frame** (origin at the
entry point, ``+z'`` along the beam), built per history from the beam direction
with the same transverse-frame construction the scattering sampler uses. Because
the beam is always ``+z'`` in this frame, the (non-rotation-equivariant)
scattering sampler is unchanged and cross-backend parity is preserved; multiple
Coulomb scattering is isotropic, so it is automatically rotation-covariant.

The plane-layered phantom carries a lab unit ``normal`` (default ``+z``); the
voxel is looked up by the **material coordinate** ``u = normal . position = u0 +
m_hat . (beam position)`` with ``m_hat = R^T normal``, and the step is limited so
its advance in ``u`` (rate ``m_hat . direction``) stays within the layer. When
``normal = direction = +z`` this reduces bit-for-bit to the DEV-012 path. Scoring
is in the beam frame, so an axis-aligned and a rotated run share the identical
grid and discretisation is not confused with physics. A ``scattering`` toggle
adds a deterministic straight-ray mode (used by the deterministic
rotation-equivalence check and for scattering-only vs no-scattering studies).

The frame-invariant observables (integral depth dose in beam depth, lateral
``sigma_x'``, R80, total deposited energy, energy balance) are guaranteed
invariant under a rigid rotation because energy loss depends only on the
water-equivalent path length and MCS only on ``rho s / X0`` and ``pv``. Validated
(decision 0018): a rigidly rotated scene reproduces the axis-aligned beam-frame
depth dose to round-off with scattering off and within statistics with it on; an
oblique beam through a slab of thickness ``D`` traverses the same
water-equivalent path as a normal beam through ``D/cos(theta)``; and the
reference, Warp CPU and CUDA paths agree.

## Scoring grids decoupled from the transport grid (Stage 3, task DEV-014, decision 0019)

DEV-014 makes the scoring grid's **alignment** free, completing the decoupling of
scoring from transport and **closing milestone V3**. The scoring grids were
already resolution-independent of the transport voxel grid (deposition splits
each step by depth overlap); DEV-014 adds an `origin_mm` to `DepthDoseGrid` and a
`depth_origin_mm`/`lateral_center_mm` to `DepthLateralGrid`, so the scored window
can be placed anywhere along and across the beam. The deposition (reference and
Warp) indexes bins as ``floor((z - z_origin)/dz)`` with edges ``z_origin +
b*dz``; because ``floor((z - z_origin)/dz) = b`` guarantees ``z_origin + (b+1)*dz
> z`` for any sign of ``b``, the depth-overlap loop stays well-formed for a grid
that starts before, at, or after the entrance. The defaults (origin 0, centre 0)
reproduce the prior grids bit-for-bit.

Validated (decision 0019): the total deposited energy and the energy per depth
interval are invariant under resolution changes (a 4x-finer grid summed 4:1
equals the coarse, to round-off) and under depth-origin/lateral-centre shifts; a
grid starting past the entrance captures strictly less energy; the lateral origin
is checked discriminatingly (a lateral window shifted off the on-axis beam
captures no dose, and the energy-weighted mean lateral position is the beam axis,
not the grid centre -- which the translation-invariant ``sigma_x`` alone cannot
show, though ``sigma_x`` invariance under a lateral shift is also confirmed); and
the reference, Warp CPU and CUDA paths agree on a shifted grid. This closes **V3**
(voxelized heterogeneous geometry and materials).

## 3-D voxel grid with ray/voxel DDA traversal (Stage 3, task DEV-015, decision 0020)

DEV-015 adds a true 3-D voxel geometry (`VoxelGrid3D`: lab-axis-aligned,
`Nx*Ny*Nz`, per-voxel mass density, a single material — the density-only cut) so
the scattering path transports through a patient-like grid at arbitrary
incidence. Traversal uses an **Amanatides-Woo voxel DDA** (not Siddon: the
random-hinge step changes direction every step and is already physics-bounded, so
a whole-ray integral is the wrong amortisation). It is phrased as a strict
generalisation of the beam-frame 1-D limiter: the single material coordinate
``u = n̂·p`` and its step limit ``(voxel_z[v+1]−u)/mproj`` become **three lab-axis
coordinates** ``u_k = p0[k] + mᵏ·q`` with ``mᵏ`` the rows of the beam-frame
rotation ``R = (e1|e2|d̂)`` (the 1-D path already computes the ``z`` row),
``rateₖ = mᵏ·D``, and a per-axis distance to the next voxel face; the step is
clipped to the nearest face over the three axes. Geometry state is **recomputed
from the position and direction each step** (rather than carried), so no float
state drifts between the float64 reference and the float32 kernel. Axis-parallel
rays take an infinite face distance; on-face ties are resolved by a direction-aware
``s_tie = 1e-4 mm`` nudge with a fixed axis priority. The grid is not rotated — an
oblique beam through an axis-aligned grid already exercises full 3-D traversal.

Parallel drivers keep the 1-D `VoxelSlab` path untouched, so a single-column grid
(`Nx=Ny=1`) reproduces it bit-for-bit when the uniform spacing reconstructs the
slab's z-faces exactly (the x/y face distances are infinite and the z-arithmetic
mirrors the 1-D path; the source must enter within the grid bounding box). Validated (decision 0020): the reduction is
bit-exact; a homogeneous box reproduces `WaterSlab` at the step-partition
discretization level; an oblique beam through an off-axis dense insert stops where
an independent Siddon `∫ρ dl` oracle reaches the water CSDA range (a genuine
interior-crossing physics check); energy is conserved for a contained beam; and
the reference, Warp CPU and CUDA paths agree (statistically, since rare near-corner
face flips decorrelate float32/float64 histories per decision 0001). Scoring stays
beam-frame marginal.

Deferred: a per-voxel material map (a pure data extension of the flat arrays),
grid rotation, non-uniform spacing, and CT/Hounsfield ingestion.

## Lab-frame 3-D dose scoring (Stage 4, task DEV-016, decision 0021)

DEV-016 begins Stage 4 with a **lab-frame 3-D dose grid** (`DoseGrid3D`), the
volumetric-dose prerequisite for beamlet-resolved influence matrices. Scoring so
far was a beam-frame 2-D marginal (depth × one transverse axis); treatment
planning needs dose per voxel in a fixed patient frame, so beamlets entering at
different positions/directions accumulate into a common grid. The 3-D voxel-grid
transport path already computes the lab position `x = p0 + R·q` every step, so the
scorer needs only a deposition mapping: when a `DoseGrid3D` is supplied, each step
additionally deposits its energy at the step's **lab midpoint** into the
containing dose voxel (one accumulation per step, so energy is conserved exactly
for a beam contained in the dose grid). Point-at-midpoint deposition (the physics
step ~1 mm ≤ a dose voxel) is accurate to O(step/voxel) while conserving energy; a
path-length-splitting scoring DDA is a later refinement.

Validated (decision 0021, reference + Warp CPU/CUDA): the summed 3-D dose equals
the deposited energy; the 3-D dose projected onto the beam axis reproduces the
beam-frame depth-dose R80 within a dose voxel; the total dose is invariant to the
dose grid's resolution and alignment; and the backends agree. Deferred: a 3-D
scorer on the 1-D-geometry path, path-length-splitting deposition, and
density-driven dose (energy/mass) units.

## Beamlet-resolved scoring and sparse influence matrices (Stage 4, task DEV-017, decision 0022)

DEV-017 adds the treatment-planning core: a **dose-influence matrix** with one row
per beamlet (pencil-beam spot) and one column per dose voxel.
`TransportEngine.run_scattering_multi` transports a list of beamlets **together**
in one batched launch (each beamlet `i` seeded `seed + i`, keeping its `beamlet`
id and per-history RNG streams, without reinitialising the immutable
physics/material data); single- and multi-beamlet runs share the same
`_scatter_state` dispatch. `assemble_influence_matrix` transports each beamlet with
a `DoseGrid3D` scorer (decision 0021), thresholds its dose (dropping voxels below a
fraction of that beamlet's peak), and stores the rows as a CSR
`SparseInfluenceMatrix` (flat voxel index `(i·ny+j)·nz+k`, dose in MeV), exported
to a documented `.npz`.

The key property is **exact additivity**: because each beamlet's separate run and
its contribution to the batched broad-field run use the identical
seed-`(seed+i)` histories and dose is a linear sum of per-history deposits, the sum
of the beamlet rows equals the batched broad-field dose to round-off — no
statistical tolerance needed. Validated (decision 0022): the summed influence dose
matches `run_scattering_multi` to ~1e-13 (reference; the first V4 gate), a 1 % peak
threshold keeps > 99 % of the energy while making the matrix sparse, energy is
conserved, and the backends agree. Deferred: a GPU (voxel, beamlet) hash-table
assembly, per-beamlet scoring inside the kernel (currently one launch per beamlet),
LET/fluence/species-resolved scorers, and beamlet-resolved uncertainty (remaining
Stage-4 items).
