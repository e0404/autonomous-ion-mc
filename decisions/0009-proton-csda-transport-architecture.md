# 0009 — Proton transport architecture and longitudinal CSDA depth-dose

- Status: accepted (Stage 1, first task; the execution model of decision 0005 is re-examined here under real transport as promised)
- Date: 2026-09-10
- Task: DEV-004
- Affects: architecture, performance, numerical accuracy, validation strategy, maintainability, scientific interpretation

## Problem

Stage 1 of the roadmap builds proton transport in homogeneous water. Its first
bounded task must establish the transport *architecture* — the particle-state
layout, the GPU-oriented transport-loop organisation, the step-control scheme,
and the scoring mechanism — on a physics model simple enough to validate
deterministically, before energy-loss straggling and multiple Coulomb
scattering (which make the output stochastic) are added. This decision fixes
that architecture and the acceptance criteria for a **continuous-slowing-down**
(CSDA) proton pencil beam depositing energy along its path in water, with **no
straggling and no scattering** yet.

## Context

- `REQUIREMENTS.md`: GPU-oriented architecture must avoid a Python-level loop
  per particle or step on accelerated backends; particle state must carry
  species and beamlet identity from the start (no proton-only assumption); the
  reference Python path is the correctness oracle; scoring grids may differ
  from the transport geometry.
- Decision 0005 (shared-source execution model) is *provisional until Stage 1*.
  This task is where a real transport kernel — a per-thread ``while`` loop with
  register-resident state, data-dependent trip count, and atomic scoring —
  first exercises it. The DEV-002 range integral and DEV-003 bisection already
  showed dynamic loops and per-thread state compile on Warp CPU and CUDA;
  atomics on multi-dimensional float arrays were shown by the DEV-002 probe.
- The stopping power is available from DEV-003 as a tabulated shared-source
  function (``ionmc.physics.tabulated``) and from DEV-002 analytically. This
  task transports on the **tabulated** layer: it is the production data source,
  its interface takes a few flat arrays (kernel-friendly), and it is defined
  from 0.5 MeV.
- Without scattering the beam stays on its axis, so all energy is deposited on
  the central ray. The physically meaningful 1-D quantity is then the
  **integral depth dose** (energy deposited per depth bin, MeV), not a
  volumetric dose (which needs a lateral extent introduced with scattering).

## Candidate approaches

**Transport loop.** (1) One history per thread with an in-kernel step loop
(selected): the design every published GPU proton engine uses (decision 0002
evidence). (2) One step per thread with global synchronisation between steps:
rejected — needs a Python-level loop over steps on the accelerated backend,
which the requirements forbid, and repeatedly reloads state from global memory.

**Step integration of the energy loss.** (1) Forward Euler with a small
fractional-energy-loss step: simplest, first-order, biased in range. (2)
**Midpoint (RK2) with a fractional-energy-loss step, additionally limited to
the next depth-bin boundary (selected):** second-order in the step, so the
range is accurate at a modest step count, and the boundary limit makes each
step deposit into exactly one bin (exact binning). (3) Transport-by-range using
the precomputed CSDA range table: accurate but bypasses the stepping loop that
is the architectural point of this task; kept only as an independent
cross-check in the validation.

**Scorer.** (1) 1-D depth grid of deposited energy (selected for this task):
the natural quantity for an on-axis pencil beam; the scoring grid is decoupled
from the (homogeneous) transport geometry, as the requirements demand. (2) Full
3-D voxel dose: deferred to Stage 3 (heterogeneous geometry) — premature here,
and volumetric dose is ill-defined for a zero-width beam.

## Selected approach

### Particle state (`ionmc.transport.state.ParticleState`)

Structure-of-arrays over the history batch, one entry per history, float32 on
the device and float64 in the reference path:

| field | unit | note |
|---|---|---|
| ``position`` (x, y, z) | mm | start of the current step |
| ``direction`` (x, y, z) | — | unit vector; constant in this task |
| ``energy`` | MeV | kinetic |
| ``weight`` | — | statistical weight (1.0 here) |
| ``species`` | int id | proton = 0; present from the start |
| ``beamlet`` | int id | for influence matrices (Stage 4); present from the start |
| ``rng_state`` | uint32 | per-history stream (decision 0005); unused until straggling |
| ``status`` | int | 0 alive, 1 stopped (energy cutoff), 2 escaped geometry |

The per-thread kernel copies its history's state into registers, runs the step
loop, and writes nothing back except through the scorer (state write-back is
added when secondaries/queues arrive).

### Transport loop (one history per thread)

```
while status == alive and step_count < max_steps:
    dl = min(fractional-energy-loss step, distance to next depth bin, distance to geometry exit)
    dE = midpoint energy loss over dl                     # shared-source function
    atomic_add(depth_edep[bin(z)], weight * dE)           # kernel; mirrored in the reference loop
    z += dl;  E -= dE;  step_count += 1
    if E <= E_cut:  atomic_add(depth_edep[bin(z)], weight * E);  E = 0;  status = stopped
    if z outside geometry:  status = escaped
```

``max_steps`` is a hard cap with a per-run truncation counter (a history that
hits it is recorded, not silently dropped). ``E_cut`` defaults to the 0.5 MeV
table floor; the residual range below it (~9 µm in water) is deposited locally.

### Shared-source step physics (`ionmc.physics.transport`)

Written against the math namespace ``m`` and decorated with ``@func`` so the
identical text runs as a Warp function per thread, as float64 Python (the
reference oracle) and vectorised with numpy (unit-tested): the linear stopping
power ``S_lin = S_mass · ρ / 10`` [MeV/mm], the fractional-energy-loss step
length, and the midpoint energy loss over a step. The kernel and the reference
Python driver are the thin layer around these functions.

### Execution paths

The **reference Python** path runs the loop as a scalar per-history Python loop
calling the float64-bound step functions — the oracle. The **Warp** path runs
the kernel on CPU and CUDA. The numpy binding remains the elementwise validator
of the step physics (as in DEV-003); the transport *loop* itself is reference
Python and Warp only, because a per-history data-dependent loop does not
vectorise — this is an execution choice, not a second physics model, and is the
Stage-1 outcome of decision 0005's "provisional" status.

### Geometry and source

Homogeneous water slab (`ionmc.transport.geometry`), beam along +z from z = 0.
Monoenergetic pencil source (`ionmc.transport.source`); N identical histories
are supported (the batch/beamlet machinery) but the physics is deterministic
here, so the depth dose scales linearly with N and one history determines it.

## Acceptance criteria (fixed before the comparison)

Deterministic output, so decision 0001 (Warp CPU vs CUDA) and decision 0005
(float64 reference vs float32 Warp) apply directly; no statistical criterion is
needed until straggling (its own later decision).

| check | criterion |
|---|---|
| energy conservation (reference) | \|Σ deposited − E₀·N\| / (E₀·N) ≤ 1e-9 |
| energy conservation (Warp float32) | ≤ 1e-5 |
| distal edge (R80 of the integral depth dose) vs tabulated CSDA range R(E₀) | ≤ 0.3 % for E₀ = 100, 150, 200 MeV |
| step-size convergence (R80 on a fixed fine grid, max_fraction 0.02 vs 0.002) | drift ≤ 0.05 % |
| independent range cross-check: Σ step lengths vs R(E₀) − R(E_cut) from the CSDA table | ≤ 0.2 % |
| reference Python vs Warp CPU/CUDA: cumulative depth dose ∫₀ᶻ edep | max_z \|ΔC(z)\| / C_total ≤ 1e-4 |
| reference Python vs Warp CPU/CUDA: R80 | within one bin |
| reference Python vs Warp CPU/CUDA: per-bin (edge-discontinuity-aware) | normalized \|Δ\|/(1e-3·peak + 2e-3·\|ref\|) ≤ 1 |
| Warp CPU vs Warp CUDA: cumulative depth dose | max_z \|ΔC(z)\| / C_total ≤ 1e-5 |
| Warp CPU vs Warp CUDA: per-bin | normalized \|Δ\|/(1e-5·peak + 1e-4·\|ref\|) ≤ 1 |

**Why the cross-backend metric is cumulative, not per-bin relative.** Without
straggling the distal edge of the depth dose is a near-discontinuity (it falls
from the maximum to zero within about one bin). A per-bin *relative* difference
there measures the sub-bin position of that edge under float32 versus float64
arithmetic, not the agreement of the physics — the same zero-crossing inflation
decision 0001 identified. The **cumulative** depth dose ∫₀ᶻ edep is insensitive
to which side of a bin boundary the edge energy lands on, so it isolates the
physical agreement; measured reference-vs-Warp-CPU cumulative disagreement is
≤ 2.7e-5 of the total (a ~4× margin under the 1e-4 gate), while a genuine
transport-model divergence would be orders of magnitude larger. The per-bin
gate is kept but written in decision 0001's mixed abs+rel form with the
absolute term scaled to the peak, so it is not inflated at the edge. DEV-005
(straggling) smears the edge over millimetres and removes this sensitivity.

## Rationale

One-history-per-thread with an in-kernel loop is the validated GPU transport
pattern and is the only one consistent with the no-Python-per-step requirement.
Doing energy loss first, deterministically, lets the whole architecture
(register-resident state, dynamic step loop, atomic depth scoring, geometry and
scorer separation, beamlet/species fields) be built and validated against
exact, tolerance-free physics — energy conservation and the CSDA range — before
stochastic straggling and scattering are layered on. The midpoint step with a
bin-boundary limit gives an accurate range and exact binning at a modest step
count. Carrying species, beamlet, RNG and status from the first version avoids a
state-layout rewrite when secondaries and beamlets arrive.

## Expected tradeoffs

- No Bragg-peak *shape* physics yet: without straggling the distal edge is a
  near-step, so R80 ≈ R100 ≈ R_CSDA. The peak-to-entrance ratio is the CSDA
  stopping-power ratio, not the measured Bragg ratio; DEV-005 (straggling) adds
  the realistic shape and its own peak-width validation.
- Integral depth dose, not volumetric dose: dose needs a lateral extent, which
  arrives with multiple scattering (DEV-005/006) or a broad beam.
- The step loop is sequential per thread; occupancy and divergence tuning are a
  Stage 6 concern, not now.
- Transport uses the tabulated stopping power only; the analytic layer is an
  alternative source and an independent range cross-check.

## Validation strategy

Unit tests of the shared-source step physics (all three bindings) and of the
state/source/geometry/scorer; a reference-vs-Warp and CPU-vs-CUDA depth-dose
comparison and the energy-conservation and range checks in
`tests/ionmc/test_transport.py`; and `validation/v1_depth_dose_csda.py` on the
host runner for the exact SHA (Warp CPU and CUDA), recorded below.

## Later validation outcome

To be filled in from the DEV-004 host run.
