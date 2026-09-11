# Validation

This section documents software, numerical, physical, statistical, and
cross-backend validation. The pages under this section describe the
diagnostic instruments; this page lists the scientific validation performed
so far and the rules it follows.

## Rules

- Acceptance quantities, metrics, tolerances and reference data are fixed in
  the decision that opens a stage, before the comparison is run (roadmap,
  *Rules applied to every validation milestone*).
- Decision `0001` covers deterministic Warp CPU-vs-CUDA parity by expectation
  class; decision `0005` adds the criterion for the float64 reference path
  versus the float32 Warp kernels; stochastic transport output will need its
  own criterion (Stage 1).
- Every result is recorded for the exact committed SHA with
  ``record_local_validation``. Scripts under ``validation/`` run on the
  controlled host runner and depend only on numpy, Warp and the standard
  library.

## V0 — analytical stopping power (task DEV-002)

Script: ``validation/v0_stopping_power.py``. Reference: the 20 native PSTAR
grid values for protons in liquid water kept in
``ionmc.reference_data.pstar_water`` (NIST PSTAR via the Apache-2.0 MCsquare
repository, commit ``211eefe6``; two independent fetches agreed). Warp
results below are from host run ``RUN-20260910T080455Z-e666b469`` (RTX A6000,
CUDA 12.9, Warp 1.17.0) at task SHA ``d23c78a``, recorded in the DEV-002
local validation record; the reference-path results were reproduced in the
sandbox as well.

| comparison | criterion (decision 0006) | result |
|---|---|---|
| analytic (Python, float64) vs PSTAR, 17 points from 10 to 400 MeV | ≤ 1.0 % each | max 0.075 %, all points low by 0.01–0.07 % |
| analytic vs PSTAR at 5 / 2 / 1 MeV | regression guards 2 / 5 / 10 % | +0.90 / +3.46 / +5.20 % |
| CSDA range 100–250 MeV vs recalled PSTAR values | 0.5 % guard (informational) | within 0.04 % |
| R(200) − R(100) vs recalled 18.242 g/cm² | 0.2 % | +0.02 % |
| numpy binding vs Python binding | bitwise equal | equal |
| Warp CPU / CUDA (float32) vs Python reference, 400 energies 2–400 MeV | S: rtol 1e-5; range: rtol 2e-5 | S 4.2e-7 on both devices (normalized 0.042); range 2.4e-7 (normalized 0.012) |
| Warp CPU vs CUDA | S: rtol 4e-6, atol 1e-6 (transcendental class); range: rtol 1e-5, atol 1e-6 (iterative) | S max abs 1.5e-5 MeV cm²/g = 1 float32 ULP at 2 MeV (normalized 0.043); range 9.5e-7 g/cm² (normalized 0.006) |

## V0 — tabulated stopping power and range (task DEV-003)

Script: ``validation/v0_tabulated_stopping_power.py``. Data: the MCsquare
PSTAR and Geant4 water tables through the data layer (decision `0007`),
acquired into the host-runner cache. Reference: the analytic model
(decision `0006`) and NIST PSTAR CSDA ranges for liquid water (retrieved
2026-09-10). Criteria fixed in decision `0008` before the comparison.

| comparison | criterion | result (reference paths) |
|---|---|---|
| tabulated (PCHIP) reproduces the committed 20-point PSTAR subset | rtol 1e-12 | exact |
| tabulated vs analytic, 10-400 MeV | ≤ 1.0 % | max +0.25 % (at 11 MeV) |
| tabulated CSDA range + 0.5 MeV residual vs NIST PSTAR (100/150/200/250 MeV) | ≤ 0.1 % | ≤ 0.018 % |
| analytic CSDA range + 1 MeV residual vs NIST PSTAR | ≤ 0.5 % | ≤ 0.053 % |
| numpy binding vs Python binding | bitwise equal | equal |
| Warp CPU / CUDA (float32) vs Python reference | S rtol 1e-5; range rtol 2e-5 | S 5.8e-8, range 4.5e-8 (both devices) |
| Warp CPU vs CUDA | transcendental class (rtol 4e-6, atol 1e-6) | bitwise identical (no transcendentals) |

Warp results are from host run `RUN-20260910T100641Z-80a76b4a` (RTX A6000, CUDA 12.9, Warp 1.17.0) at
task SHA `dfd9995`, recorded in the DEV-003 local validation record.

This closes milestone V0: proton stopping power in water is now available from
both an analytical model and an I-value-consistent external table, on all three
execution paths, validated against reference data with pre-fixed tolerances.

## V1 (partial) - CSDA proton depth dose (task DEV-004)

Script: ``validation/v1_depth_dose_csda.py``. Reference: the tabulated CSDA
range (task DEV-003). Criteria fixed in decision `0009` before the comparison.
This is the deterministic, no-straggling foundation of milestone V1; the
straggling/Bragg-shape and scattering/lateral parts follow in later Stage-1
tasks.

| check | criterion | result (reference paths) |
|---|---|---|
| energy conservation (reference, float64) | \|balance\| ≤ 1e-9 | ~1e-16 |
| R80 vs tabulated CSDA range, 100/150/200 MeV | ≤ 0.3 % | ≤ 0.15 % |
| step-size convergence (fixed grid, 0.02 vs 0.002) | drift ≤ 0.05 % | ≤ 3e-4 % |
| range cross-check (stopping depth vs CSDA range) | ≤ 0.2 % | ≤ 0.15 % |
| reference vs Warp CPU/CUDA, cumulative depth dose | ≤ 1e-4 of total | ≤ 2.7e-5 (both devices) |
| Warp CPU vs CUDA, cumulative depth dose | ≤ 1e-5 of total | ≤ 2.3e-8 (near bitwise) |

Warp results are from host run `RUN-20260910T135739Z-8166c2f4` (RTX A6000, CUDA 12.9, Warp 1.17.0) at
SHA `b726507`, recorded in the DEV-004 local validation record.

## V1 (straggling) - Bragg peak and range straggling (task DEV-005)

Script: ``validation/v1_bragg_straggling.py``. References: the analytic Bohr
range-straggling integral and Bortfeld's ``sigma = 0.012 R^0.935`` (decision
0010). The metric is the range straggling ``sigma_R`` (std of stopping depths),
which isolates straggling from the peak shape.

| check | criterion | result (reference paths) |
|---|---|---|
| mean range vs CSDA range, 100/150/200 MeV | within 0.1 % | <0.01 % |
| sigma_R vs analytic Bohr integral | within 3 % | <1 % |
| sigma_R vs Bortfeld 0.012 R^0.935 | within 10 % | +5-7 % (Bohr electron-binding gap) |
| sigma_R / R band | 0.9-1.2 % | 1.04-1.11 % |
| energy conservation (reference float64 / Warp float32) | 1e-9 / 5e-5 | ~1e-16 / ~1e-5 |
| Bragg peak present (peak > 3x entrance, deep) | yes | yes |
| two independent seeds statistically consistent | RMS(t) in [0.5, 1.6] | ~1.1 |
| reference vs Warp / Warp CPU vs CUDA cumulative | 1e-4 / 1e-5 | 5.8e-7 / 8.1e-8 |

Warp results are from host run `RUN-20260910T144031Z-8dd34007` (RTX A6000, CUDA 12.9, Warp 1.17.0) at
SHA `100653d`, recorded in the DEV-005 local validation record.

## V1 (scattering) - lateral spread from multiple Coulomb scattering (task DEV-006)

Script: ``validation/v1_lateral_scattering.py``. Reference: the Fermi-Eyges
lateral-spread oracle (``ionmc.physics.fermi_eyges``) built from the same
scattering power, and published values (decision 0011). Metric: lateral
``sigma_x(z)`` from the 2-D depth-lateral dose.

| check | criterion | result (reference paths) |
|---|---|---|
| theta0 / pv vs Highland hand values | within 0.5 % | matches |
| MC sigma_x(z) vs Fermi-Eyges at 0.5R, 0.8R (150/200 MeV) | within 3 % | within 1 % |
| MC sigma_x(0.8R) vs published (2.4 / 3.9 mm) | within 8 % | within ~2 % |
| energy conservation (reference / Warp) | 1e-9 / 5e-5 | ~1e-16 / ~1e-6 |
| depth-dose Bragg peak (marginal) | present | peak/entrance ~6.8 |
| detour factor (mean projected range vs CSDA) | shortened ~0.1 % | -0.10 % |
| Warp CPU vs CUDA sigma_x(0.8R) | within 0.02 mm | agree to 5 digits |

Warp results are from host run `RUN-20260910T151223Z-286dcc42` (RTX A6000, CUDA 12.9, Warp 1.17.0) at
SHA `69d5401`, recorded in the DEV-006 local validation record.

This closes the physics of milestone **V1** (proton transport in homogeneous water).

The ICRU 90 mean excitation energy (78 eV) is reported as a separate model
difference (−1.14 % at 1 MeV to −0.43 % at 400 MeV), not absorbed into any
tolerance.

## V2 (attenuation) - proton nonelastic nuclear removal (task DEV-007)

Script: ``validation/v2_nuclear_attenuation.py``. References: the published
primary survival to the Bragg peak (Paganetti 2002; Gottschalk), the analytic
``1 - exp(-integral Sigma/S dE)`` reaction fraction, and cross-backend parity
(decision 0012). Primaries are removed catastrophically at the macroscopic
nonelastic rate on oxygen; a local fraction (``f_local = 0.30``) is deposited
at the vertex and the remainder is booked to an audited escaping channel.

| check | criterion | result (reference / Warp CPU) |
|---|---|---|
| primary survival to peak, 150 MeV | 0.80 ± 0.03 | 0.813 |
| primary survival to peak, 200 MeV | 0.73 ± 0.04 | 0.716 |
| MC reaction fraction vs analytic | within 5σ + 0.01 | 0.187 vs 0.189; 0.284 vs 0.285 |
| energy budget deposited + escaped = in (reference / Warp) | 1e-9 / 1e-5 | ~1e-16 / ~4e-7 |
| ``nuclear=False`` regression | zero reactions, exact EM budget | identical |
| reference vs Warp CPU (matched N, seed) | same reactions, cumulative ≤ 1e-4 | identical set, ~8e-7 |

The absolute Bragg peak-to-entrance ratio still reads high until DEV-008
transports the escaping secondary protons; the V2 attenuation gate therefore
checks primary removal and the energy budget, not the peak-to-entrance ratio.
Milestone **V2** opens here and closes with DEV-008. Warp CPU/CUDA results are
recorded in the DEV-007 local validation record (host run below).

## V2 (secondaries) - secondary charged-particle transport (task DEV-008)

Script: ``validation/v2_secondary_transport.py``. References: the published
secondary-proton dose contribution (Paganetti 2002) and cross-backend parity
(decision 0013). Each nonelastic reaction is partitioned into local heavy
fragments, transported secondary protons, and truly escaping neutrals; the
secondary protons are generated host-side from the per-history reaction records
and transported in a second pass through the same proton engine.

| check | criterion | result (reference / Warp CPU) |
|---|---|---|
| secondary dose fraction at entrance (150/200 MeV) | 0.5-4 % | 2.1 % / 2.5 % |
| secondary dose fraction of total dose | 2-12 % | 4.0 % / 6.9 % |
| secondary fraction plateau mean (rises entrance→plateau, drops at peak) | 3-12 %, > entrance, > peak | 5.4 % / 8.6 % |
| energy budget deposited + escaped = in (reference / Warp) | 1e-9 / 1e-5 | ~0 / ~4e-7 |
| secondaries-off regression | zero secondaries, no second pass | identical |
| reference vs Warp CPU (matched N, seed) | same secondaries, cumulative ≤ 1e-4 | identical set, ~1e-6 |

This **closes milestone V2** (nuclear interactions for protons in water). The
absolute peak-to-entrance ratio is not gated (it is geometry- and convention-
dependent); the acceptance is the secondary-dose fraction and plateau shape plus
the unchanged DEV-007 survival gate. Warp CPU/CUDA results are recorded in the
DEV-008 local validation record (host run below).

The V0 milestone of the roadmap is closed by task ``DEV-003`` (tabulated
layer), which will replace the recalled range values by table-integrated
ones from the acquired dataset.

## V3 (opening) - voxelized density heterogeneity (task DEV-009)

Script: ``validation/v3_density_heterogeneity.py``. References: the analytic
water-equivalent-thickness relation (CSDA range scales as ``1/rho``) and cross-
backend parity (decision 0014). A ``VoxelSlab`` gives each voxel its own mass
density; the transport looks up the local density per step.

| check | criterion | result |
|---|---|---|
| homogeneous equivalence (uniform voxels vs WaterSlab) | bit-exact (reference) | max diff 0 |
| WET scaling: R80 at R_water/rho (rho 0.5, 1.2; 150/200 MeV) | within 0.3 % | ~1e-4 |
| layered interface: peak shift = extra WET (20 mm x 0.85) | within 0.5 mm | 17.0 mm |
| energy conservation with nuclear+secondaries across interface | 1e-9 (reference) | ~0 |
| reference vs Warp CPU across the interface | same reactions, cumulative <= 1e-4 | ~6e-7 |

This **opens milestone V3** (voxelized heterogeneous geometry). Deferred to later
Stage-3 tasks: 3-D voxel geometry and arbitrary incidence, per-voxel material
composition, and decoupled scoring grids. Warp CPU/CUDA results are recorded in
the DEV-009 local validation record (host run below).

## V3 (materials) - per-voxel tissue materials (task DEV-010)

Script: ``validation/v3_material_composition.py``. References: the published
Schneider/ICRU water-equivalent ratios and the analytic stopping-power ratio
(decision 0015). Each voxel transports as water at its water-equivalent density
``SPR(material) x rho`` with a composition-scaled nuclear rate.

| check | criterion | result |
|---|---|---|
| water VoxelSlab vs homogeneous WaterSlab | bit-exact (reference) | max diff 0 |
| water-equivalent ratio in published band | bone 1.60-1.72, adipose 0.95-0.98, ... | 1.70 / 0.97 / 1.04 |
| R80 at R_water/WER (bone, adipose, muscle; 150/200 MeV) | within 0.3 % | ~1e-4 |
| nuclear composition scaling (bone / adipose vs oxygen-only) | ~2x / ~3.4x, n_O for water | 2.08 / 3.40 |
| energy conservation across material interfaces | 1e-9 (reference) | ~0 |
| reference vs Warp CPU across material interfaces | same reactions, cumulative <= 1e-4 | ~4e-7 |

The **independent physics check** is the water-equivalent-ratio band membership;
the R80-at-``R_water/WER`` check is a *self-consistency* check (it uses the same
analytic SPR to build both the transported density and the expectation), so it
confirms the engine applies the SPR correctly, not the SPR value itself.

Deferred: energy-dependent SPR / per-material stopping tables (bone's 2.5 % energy
dependence), element-specific nuclear cross sections, and non-water/heterogeneous-
material scattering (the 3-D path rejects both). Warp CPU/CUDA results are
recorded in the DEV-010 local validation record (host run below).

## V3 (scattering) - density-heterogeneous 3-D scattering (task DEV-011)

Script: ``validation/v3_scattering_heterogeneity.py``. References: the Fermi-Eyges
lateral-spread oracle (uniform density) and a new piecewise-density Fermi-Eyges
variant (layered), plus cross-backend parity (decision 0016). The 3-D scattering
path now looks up the local water density by depth.

| check | criterion | result |
|---|---|---|
| uniform water VoxelSlab vs WaterSlab scattering | bit-exact (reference) | max diff 0 |
| uniform-density lateral sigma_x vs Fermi-Eyges (0.5R, 0.8R) | within 3 % | < 1 % |
| layered sigma_x vs piecewise-density Fermi-Eyges (60/90/120 mm) | within 3 % | < 1 % |
| energy conservation (layered scattering) | 1e-9 (reference) | ~0 |
| reference vs Warp CPU (layered) | sigma_x <= 0.05 mm, depth dose cumulative <= 5e-4 | ~1e-4 |
| CUDA sigma_x vs Fermi-Eyges (its own physics gate) | within 3 % | < 1 % |
| CUDA vs CPU (layered depth dose) | cumulative <= 5e-4 | ~4e-6 |

The depth-dose cross-backend tolerance (5e-4) is looser than the 1-D path's 1e-4
because the 3-D scattering kernel accumulates more float32 rounding at voxel
boundaries; the lateral ``sigma_x`` is the tight cross-backend metric. Warp
CPU/CUDA results are recorded in the DEV-011 local validation record (host run
below).

## V3 (material scattering) - non-water materials on the scattering path (task DEV-012)

Script: ``validation/v3_material_scattering.py``. References: a **material-aware**
Fermi-Eyges oracle (``lateral_sigma_x_material_mm``: energy vs depth from the
water-equivalent thickness, scattering power from the physical density and
material radiation length) and cross-backend parity (decision 0017).

| check | criterion | result |
|---|---|---|
| water VoxelSlab vs WaterSlab scattering | bit-exact (reference) | max diff 0 |
| homogeneous bone sigma_x vs material Fermi-Eyges | within 3 % | < 1 % |
| water/bone/water interface sigma_x vs material Fermi-Eyges | within 3 % | < 1 % |
| energy conservation (material scattering) | 1e-9 (reference) | ~0 |
| reference vs Warp CPU (material interface) | sigma_x <= 0.05 mm, depth dose <= 5e-4 | within budget |
| CUDA sigma_x vs material oracle / CUDA vs CPU depth dose | < 3 % / <= 5e-4 | < 1 % / tight |

The depth-dose and scattering paths now share the same per-voxel material
capability. Deferred: 3-D voxel geometry with arbitrary incidence and decoupled
scoring grids (V3 closure). Warp CPU/CUDA results are recorded in the DEV-012
local validation record (host run below).

## V3 (arbitrary incidence) - beam-frame transport, rotated-vs-axis-aligned (task DEV-013)

Script: ``validation/v3_arbitrary_incidence.py``. References: the axis-aligned
run itself (rigid-rotation invariance), the homogeneous Fermi-Eyges oracle, an
analytic oblique water-equivalent-path relation, and cross-backend parity
(decision 0018). Tolerances are pre-registered before running.

| check | criterion | result |
|---|---|---|
| axis-aligned energy conservation | 1e-9 (reference) | ~0 |
| deterministic rotation equivalence (scattering off): energy / depth-dose cumulative / R80 | 1e-12 / 1e-9 / 1e-4 mm | round-off |
| statistical rotation covariance (scattering on): rotated sigma_x' vs Fermi-Eyges, axis-vs-rot R80 | within 3 % / <= 1 mm | < 1 % / within stats |
| oblique WET traversal (heterogeneous, interior boundaries): tilted vs cos-stretched-normal deposited energy + profile | within 1e-6 / 1e-9 | round-off |
| reference vs Warp CPU (rotated): sigma_x', depth dose cumulative | <= 0.05 mm / <= 5e-4 | within budget |
| CUDA sigma_x' vs Fermi-Eyges / CUDA vs CPU depth dose | < 3 % / <= 5e-4 | within budget |

The deterministic gate isolates the frame/geometry transform (a straight ray has
no azimuthal gauge, so a rotated scene is the same history in different
coordinates, agreeing to round-off); the statistical gate confirms the physics
is genuinely rotation-covariant with scattering on; the oblique gate exercises
the interior oblique voxel-boundary traversal through a heterogeneous slab that a
stopping beam in a thick homogeneous slab hides. A unit test additionally asserts
the beam-frame material-coordinate identity ``u = normal . (p0 + R . p_beam)``
directly, so a wrong transverse cross-term is falsifiable.
Deferred: arbitrary 3-D per-voxel maps with Siddon traversal and decoupled
scoring grids (V3 closure, DEV-014). Warp CPU/CUDA results are recorded in the
DEV-013 local validation record (host run below).

## V3 (scoring grid) - grid-independence of integral dose (task DEV-014, closes V3)

Script: ``validation/v3_scoring_grid.py``. References: the transport run itself
(grid-independence under scoring-grid changes) and cross-backend parity (decision
0019). Tolerances are pre-registered before running.

| check | criterion | result |
|---|---|---|
| integral dose vs resolution (total; fine summed 4:1 vs coarse) | 1e-12 / 1e-12 | round-off |
| integral dose vs alignment (depth-origin + lateral-centre shift) | 1e-12 | round-off |
| lateral origin discriminated (off-beam window captures none; mean_x = beam axis) | 0 MeV / < 0.2 mm | holds |
| lateral sigma_x invariant under a lateral shift | <= 1e-6 mm | round-off |
| partial coverage (grid past the entrance captures strictly less) | 0 < E_down < E_full | holds |
| reference vs Warp CPU (shifted grid): depth dose cumulative | <= 5e-4 | within budget |
| CUDA sigma_x vs Fermi-Eyges / CUDA vs CPU (shifted grid) | < 3 % / <= 5e-4 | within budget |

The scoring grid is now decoupled from the transport grid in resolution *and*
alignment, satisfying the last **V3** gate: **milestone V3 is closed**. Warp
CPU/CUDA results are recorded in the DEV-014 local validation record (host run
below).

## V3 (voxel grid) - 3-D voxel geometry with ray/voxel DDA (task DEV-015)

Script: ``validation/v3_voxel_grid_3d.py``. References: the 1-D VoxelSlab / homogeneous
WaterSlab (reduction), an independent Siddon water-equivalent-path oracle, and
cross-backend parity (decision 0020). Tolerances are pre-registered before running.

| check | criterion | result |
|---|---|---|
| single-column grid reduces to VoxelSlab (deterministic) | depth-dose cumulative <= 1e-9 | bit-exact |
| homogeneous box vs WaterSlab (step-partition discretization) | energy <= 1e-4 / R80 <= 1e-2 mm | within budget |
| oblique beam through a dense insert vs Siddon WET range | within 1e-2 | within budget |
| energy conservation (contained beam, scattering on) | 1e-9 (reference) | ~0 |
| reference vs Warp CPU (grid, scattering): depth dose / sigma_x | <= 3e-3 (DDA face-flip) / <= 0.05 mm | within budget |
| CUDA sigma_x vs Fermi-Eyges / CUDA vs CPU depth dose (same precision) | < 3 % / <= 5e-4 | within budget |

This adds the general 3-D voxel geometry (beyond the slab/layered phantoms the V3
gates used); the single-column reduction is the strict-generalization safety anchor
and the Siddon oracle is an implementation-independent interior-crossing physics
check. Warp CPU/CUDA results are recorded in the DEV-015 local validation record
(host run below).

## V4 (dose scoring) - lab-frame 3-D dose on the voxel-grid path (task DEV-016)

Script: ``validation/v4_dose3d.py``. References: energy conservation, the
beam-frame depth dose (R80 consistency), and cross-backend parity (decision 0021).

| check | criterion | result |
|---|---|---|
| energy conservation (dose sum == deposited; deterministic + scattering) | 1e-9 | ~0 |
| depth-dose consistency (3-D z-marginal R80 vs beam-frame R80) | <= 1 mm (a dose voxel) | within budget |
| grid-independence (total dose vs dose-grid resolution/alignment) | 1e-12 | round-off |
| reference vs Warp CPU dose (deterministic): total / per-voxel | <= 1e-5 / <= 5e-3 | float32 point deposition |
| CUDA vs CPU total dose (same precision) | <= 1e-4 | within budget |

This is the volumetric-dose prerequisite for beamlet-resolved influence matrices.
Warp CPU/CUDA results are recorded in the DEV-016 local validation record (host
run below).

## V4 (influence) - beamlet-resolved scoring and sparse influence matrices (task DEV-017)

Script: ``validation/v4_influence.py``. References: the batched broad-field dose
(beamlet-sum additivity), the dense matrix (sparse-vs-dense), energy conservation,
and cross-backend parity (decision 0022).

| check | criterion | result |
|---|---|---|
| beamlet sum == broad-field dose (the first V4 gate) | max rel diff <= 1e-9 | ~1e-18 |
| sparse vs dense (1 % threshold keeps energy) | >= 99 % kept, nnz shrinks | ~99.4 % |
| energy conservation (summed influence == deposited) | 1e-9 | ~0 |
| reference vs Warp CPU broad-field dose (deterministic) | total <= 1e-5 / per-voxel <= 5e-3 | 1.7e-7 / 2.7e-5 |
| CUDA vs CPU broad-field dose (deterministic, per voxel) | total <= 1e-5 / per-voxel <= 5e-3 | 2.0e-9 / 1.4e-7 |
| CUDA vs CPU broad-field total (scattering on) | total <= 1e-4 | ~1.9e-10 |
| CUDA vs CPU per voxel (scattering on) | diagnostic only (decorrelates) | ~3 % of peak |

The beamlet-sum gate is exact (round-off) because the per-beamlet and batched runs
share the identical seed-partitioned histories. Warp CPU/CUDA results are recorded
in the DEV-017 local validation record (host run below).

## V4 (LET) - dose-averaged LET scoring (task DEV-018)

Script: ``validation/v4_let.py``. References: NIST PSTAR water electronic stopping
power (thin-voxel analytic limit), the scored dose energy (dose-weight
consistency), and cross-backend parity (decisions 0023, 0022).

| check | criterion | result |
|---|---|---|
| thin-voxel LET_d == PSTAR S(E) (250/150/100 MeV) | rel <= 2 % and LET_d >= S(E) | ~0.3 % |
| dose-weight consistency (den == dose, num>0 where dose>0) | exact | pass |
| distal LET peak (entrance ~0.5445 keV/um, peak distal to dose) | entrance rel <= 5 %, peak > 3 keV/um | 0.545 / ~14 |
| voxel-size convergence (plateau LET_d vs dose-voxel z-size) | rel spread < 3 % over 0.5/1/2 mm | ~0.3 % |
| reference vs Warp CPU LET_d (deterministic, per voxel) | per-voxel <= 5e-3 | ~4.8e-5 |
| CUDA vs CPU LET_d (deterministic, per voxel) | per-voxel <= 5e-3 | tight (float32) |
| CUDA vs CPU LET numerator total (scattering on) | total <= 1e-4 | recorded (host) |
| CUDA vs CPU LET_d per voxel (scattering on) | diagnostic only (decorrelates) | recorded (host) |

LET_d "Method C" uses the tabulated `S(E_mid)`, so it is immune to step
face-clipping and straggling; the denominator equals the scored dose energy by
construction. The LET_d ratio inherits the dose scorer's cross-backend story
(decision 0022): the deterministic per-voxel comparison is the tight, discriminating
spatial gate, while under scattering only the totals stay tight and the per-voxel
ratio is a recorded diagnostic.

## V4 (uncertainty) - batch-based statistical uncertainty for 3-D dose (task DEV-019)

Script: ``validation/v4_uncertainty.py``. `run_scattering_batched` runs independent
history batches into a `DoseGrid3D` and reports the per-voxel mean and standard
error of the mean (decision 0024).

| check | criterion | result |
|---|---|---|
| se_scaling (4x histories -> ~0.5x rel SEM, shared mask) | ratio in [0.35, 0.71] | ~0.53 |
| mean energy conservation (batch-mean total == input) | rel <= 1e-9 | ~0 |
| uncertainty sanity (SEM>0 in dose, ==0 outside) | exact | pass |
| reference vs Warp CPU batched-mean total dose | rel <= 1e-4 | ~1.4e-9 |
| CUDA vs CPU batched-mean total dose | rel <= 1e-4 | recorded (host) |

The standard error of the mean is `std(batches, ddof=1)/√n_batches`; the `1/√N`
scaling is checked over a **shared** high-dose voxel mask so both runs average over
the identical voxel set (per-voxel SEM estimates are individually noisy at a
practical batch count). Under scattering the per-voxel batched *mean* decorrelates
across backends like dose (decision 0022), so the cross-backend gate is on the tight
total.

## V4 (beamlet uncertainty) - beamlet-resolved uncertainty (task DEV-020)

Script: ``validation/v4_beamlet_uncertainty.py``. `assemble_influence_matrix_batched`
stores each beamlet's per-voxel mean dose and its standard error (`data_sigma`) in
the sparse influence matrix (decision 0025).

| check | criterion | result |
|---|---|---|
| sigma alignment (data_sigma vs data, finite, >= 0) | exact | pass |
| per-beamlet SE scaling (4x histories -> ~0.5x rel SEM) | ratio in [0.35, 0.71] | ~0.53 |
| beamlet sum within statistics (vs high-stat broad field) | >= 95 % of hi voxels within 4 sigma | ~99.9 % |
| reference vs Warp CPU per-beamlet dose (deterministic, per voxel) | per-voxel <= 5e-3 | ~1.8e-5 |
| CUDA vs CPU per-beamlet dose (deterministic, per voxel) | per-voxel <= 5e-3 | recorded (host) |

The **beamlet_sum_within_statistics** gate is the V4 milestone gate ("sum of beamlet
doses equals the broad-field dose within statistics") in its full statistical form:
the summed independent-batch beamlet means agree per voxel with a high-statistics
broad field within the combined SEM. It complements DEV-017's exact same-partition
round-off check. Per-beamlet variances add for the broad-field SEM (`√Σσ²`) because
the beamlets are transported as independent batched estimates. The per-voxel
`z = |Δ|/σ` uses a ~16-batch SEM, so it is t-distributed (heavier tails than
normal); the loose 4σ / 95 %-coverage band accommodates that. The cross-backend
gate validates the per-beamlet *mean*; the per-beamlet *SEM* inherits its
cross-backend soundness from the single-source batched estimator (decision 0024).

## Unit and regression tests

``tests/ionmc/`` covers units and constants, materials, the RNG mirror
(pinned to a kernel dump from Warp 1.17.0 on the workstation), the binding
mechanism, and the stopping-power model on the reference paths; Warp tests
are skipped where Warp or CUDA is unavailable. GitHub CI runs the suite
without a GPU.
