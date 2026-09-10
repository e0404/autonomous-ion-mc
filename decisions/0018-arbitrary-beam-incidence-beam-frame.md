# 0018 — Arbitrary beam incidence via a canonical beam frame

- Status: accepted
- Date: 2026-09-11
- Task: DEV-013
- Affects: architecture, physical accuracy, validation strategy

## Problem

Through DEV-012 the 3-D multiple-scattering path transports a proton pencil
beam that always enters at the origin along **+z**, with the voxel geometry
looked up by the depth coordinate `pz` and scored in a fixed depth/lateral grid
anchored to that axis. Stage 3 and milestone **V3** require *arbitrary beam
incidence* and *beam-frame transforms* (`REQUIREMENTS.md`, *Arbitrary beam
incidence*; roadmap Stage 3), validated by the **rotated-vs-axis-aligned beam
equivalence** gate. This task lets the source enter at an arbitrary position and
unit direction, and the whole scene (beam + phantom) be rigidly rotated, without
changing any interaction physics.

Scope is deliberately bounded to a **plane-layered** phantom: interfaces are
parallel planes with a (possibly rotated) lab normal `n̂`, so geometry traversal
is *distance along the ray to the next interface plane*. A full arbitrary 3-D
per-voxel map with Siddon/DDA ray traversal and a decoupled 3-D scoring grid is
the remaining V3 work (DEV-014).

## Context and evidence

- **The state is already 3-D.** `ParticleState` carries `(n,3)` position and
  direction; the scattering reference loop and Warp kernel already move a 3-D
  position and direction. Only the *source* (hardcoded +z) and the *geometry
  lookup* (`searchsorted(voxel_z, pz)`, step limit `(voxel_z[v+1]-pz)/dz`) and
  the *scoring anchor* were axis-locked.
- **Frame choice determines whether MCS covariance is automatic.** The
  scattering sampler (`_scatter_direction`) builds its transverse frame from the
  *least-aligned world axis*; that construction is **not** rotation-equivariant.
  Transporting in the lab frame with a tilted direction would therefore require
  proving `F(R·d) = R·F(d)`, which fails at the sample level. Transporting in a
  **canonical beam frame** (origin at the entry point, `+z'` along the beam)
  side-steps this entirely: the sampler only ever sees a beam-frame direction,
  exactly as before, so its gauge choice is harmless and the existing bit-exact
  cross-backend parity (decision `0001`, `0011`) is preserved.
- **Frame-invariant observables are guaranteed invariant.** Energy loss depends
  only on the water-equivalent path length `∫ρ dl` (a rigid-motion invariant);
  Highland MCS is azimuthally isotropic and depends only on `ρ s / X0` and `pv`.
  So integral depth dose (in beam depth / WET), lateral spread, R80, total
  deposited energy and energy balance are rotation invariants; absolute lab
  positions are covariant, not invariant, so scoring is done in the beam frame.

## Decision

1. **Canonical beam-frame transport.** For each history the engine builds an
   orthonormal frame `R = (e1, e2, d̂)` from the unit beam direction `d̂` with the
   **same** transverse-frame construction the scattering sampler uses (factored
   as `_transverse_frame`, shared by the reference path and mirrored in the Warp
   kernel). Transport runs in the beam frame: position starts at `0` (relative to
   the entry point `p0`), direction at `+z'`. Scoring uses the beam-frame
   position components (depth = `z'`, lateral = `x'`), so an axis-aligned and a
   rotated run score into the *identical* beam-frame grid and discretisation is
   not confused with physics.

2. **Material-coordinate plane traversal.** The layer stack has a lab unit normal
   `n̂` (default `+z`) and boundaries `voxel_z` measured along it. The engine
   looks up the voxel by the **material coordinate** `u = n̂·position = u0 +
   m̂·(beam position)`, with `m̂ = Rᵀ n̂` and `u0 = n̂·p0`. The step is limited so
   its advance in `u` (rate `m̂·direction`) stays within the current layer and the
   geometry, replacing the former `pz`/`dz` depth limit. When `n̂ = d̂ = +z` this
   reduces **bit-for-bit** to the DEV-012 path.

3. **Arbitrary source.** `PencilBeamSource` gains a `direction` (default `+z`,
   normalised on use); `WaterSlab`/`VoxelSlab` gain a `normal` (default `+z`).
   The pure depth-dose path (`TransportEngine.run`) is longitudinal and rejects a
   non-+z beam with an actionable error; arbitrary incidence is a scattering-path
   (`run_scattering`) capability.

4. **Scattering toggle.** A `scattering` engine flag (default on, symmetric with
   `straggling`) allows a deterministic straight-ray transport. This enables the
   tight deterministic rotation-equivalence gate below (and scattering-only vs
   no-scattering studies) without perturbing the default behaviour.

## Validation (milestone V3, rotated-vs-axis-aligned equivalence)

Pre-registered before running (decision `0001` methodology), in
`validation/v3_arbitrary_incidence.py` on the host runner (CPU + CUDA):

- **Axis-aligned reduction** — the beam-frame path with `n̂ = d̂ = +z` reproduces
  the DEV-012 `run_scattering` result **bit-for-bit** (`max|Δ| == 0`).
- **Deterministic rotation equivalence (scattering off)** — a rigidly rotated
  scene reproduces the axis-aligned beam-frame depth dose to round-off:
  total deposited energy `|ΔE|/E ≤ 1e-12`, depth-dose cumulative
  `max|Δcumsum|/ΣE ≤ 1e-9`, R80 `≤ 1e-4 mm` (reference float64).
- **Statistical rotation covariance (scattering on)** — with the frame
  reconstructed from `d̂` the azimuthal gauge differs, so histories agree only in
  distribution: `σ_x'` vs the Fermi-Eyges oracle `≤ 3 %`, R80 within the batch
  SEM, energy balance `≤ 1e-9` (reference) / `≤ 1e-5` (Warp float32).
- **Oblique WET traversal (analytic)** — a beam tilted by `θ` through a lab-fixed
  slab of thickness `D` deposits the same energy as a normal beam through
  `D/cosθ` (same water-equivalent path), to `≤ 1e-3`; this exercises the
  `m̂ ≠ ẑ'` plane-traversal code that a stopping beam in a thick slab would hide.
- **Cross-backend** — reference vs Warp CPU vs CUDA on a rotated config:
  depth-dose cumulative `≤ 5e-4`, `σ_x'` `≤ 0.05 mm`, consistent with
  decisions `0001`/`0011`.

## Consequences

- Arbitrary beam incidence is available on the scattering path for plane-layered
  phantoms, satisfying the V3 rotated-vs-axis-aligned equivalence gate. Scoring
  is in the beam frame (a pencil beam is centred on its own axis).
- The depth-dose path stays longitudinal (+z only) and is unchanged; its guard
  prevents silent misuse with an oblique beam.
- **Deferred to DEV-014 (remaining V3 work):** arbitrary 3-D per-voxel material
  maps with robust ray/voxel (Siddon/DDA) traversal, and scoring grids decoupled
  from the transport grid in resolution *and* alignment (the grid-independence
  gate). The beam-frame machinery introduced here is the prerequisite for both.
