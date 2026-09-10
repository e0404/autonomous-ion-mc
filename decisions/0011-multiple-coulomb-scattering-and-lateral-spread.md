# 0011 — Multiple Coulomb scattering, lateral spread, and 2-D scoring

- Status: accepted
- Date: 2026-09-10
- Task: DEV-006
- Affects: physical accuracy, architecture, validation strategy, scientific interpretation, performance

## Problem

DEV-004/005 transported protons along the beam axis with energy loss and
straggling (a Bragg peak, but no lateral spread). This task adds **multiple
Coulomb scattering (MCS)** so protons deflect and the pencil beam acquires a
growing lateral width `sigma_x(z)`, turning the on-axis integral depth dose
into a genuine 2-D (depth x lateral) distribution. It fixes, before the first
comparison: the scattering-angle model, how the deflection is placed within a
step (which controls lateral displacement), the direction-update geometry, the
2-D scoring, and the acceptance targets. This closes the physics of milestone
**V1**.

## Context and evidence

Research (Claude physics-researcher subagent, 2026-09-10; sources below) and
in-repo checks established:

- **Highland / Lynch-Dahl** projected RMS angle over thickness `x`:
  `theta0 = (13.6 MeV / (beta c p)) z sqrt(x/X0) [1 + 0.038 ln(x/(X0 beta^2))]`,
  with `beta c p = pv = T(T+2M)/(T+M)` [MeV], `X0(water) = 36.08 g/cm^2` (PDG).
  `theta0` is the **projected** (one-plane) RMS; the space-angle variance is
  `2 theta0^2`. Verified: pv = 190.37 / 364.86 MeV and theta0(1 mm, bracket 1)
  = 3.76 / 1.96 mrad at 100 / 200 MeV.
- **The log bracket is pathological per 1 mm step** (it is not additive across
  sub-steps; for a 1 mm water step it is ~0.84, and its argument goes negative
  for thin steps). The scattering-power form (bracket = 1) is the clean choice:
  it makes the per-step model and the Fermi-Eyges oracle use the *identical*
  scattering power `T = (13.6/pv)^2 / X0` and reproduces published lateral
  spreads within a few percent.
- **Random hinge** (Kawrakow/Fippel): the full-step deflection applied at a
  uniform random point along the step recovers the correct lateral variance
  `<y^2> = theta0^2 s^2 / 3` per plane; applying it at the step end gives zero
  displacement, at the start 3x too much.
- **Fermi-Eyges**: `sigma_x^2(z) = integral_0^z (z-u)^2 T(u) du` (projected) is
  the analytic lateral-spread oracle.
- A Gaussian central-MCS model (no single-scattering tail) is standard for
  fast therapy MC (MCsquare, FRED, MOQUI) and adequate for the lateral-`sigma`
  core targeted here; the far penumbra tail is out of scope.

## Selected approach

### Scattering model (shared source `ionmc.physics.transport`)

`highland_theta0(E, M, z, step_mm, rho, X0)` returns the projected RMS angle
with **bracket = 1** (scattering-power form). Per step, two independent
projected angles `theta_x, theta_y ~ N(0, theta0)` are drawn in the execution
layer (`wp.randn` / the `ionmc.rng` mirror) and applied as a two-plane tilt of
the direction, giving space-angle variance `2 theta0^2` and uniform azimuth
automatically. Scattering is switched off below the same 2 MeV floor as
straggling.

### Random hinge and 3-D step

Each step of geometric length `s`: draw `a = xi s`, `xi ~ U(0,1)`; move straight
`a` along the current direction; tilt the direction; move straight `s - a`; then
deposit the step's energy. The energy-loss physics step is the fractional
-energy-loss step (decision 0010), *not* bin-limited. The step is capped so its
z-projection does not exceed the geometry.

### Direction-update geometry (mirrored, not shared source)

`d' = normalize(d + theta_x e1 + theta_y e2)` with `{e1, e2}` an orthonormal
transverse frame built from the least-aligned world axis (avoids degeneracy).
This vector geometry is written with identical plain-float arithmetic in the
reference driver (`engine._scatter_direction`) and the Warp kernel
(`_scatter_dir`, using `wp.vec3`), the same "thin mirrored layer" discipline as
the deposition loop (decision 0005): the physics stays shared, the geometry is
byte-identical between paths and checked by the cross-backend comparison.

### 2-D scoring (`ionmc.transport.DepthLateralGrid`)

Energy is scored into a (depth `z` x lateral `x`) grid, summed over the second
transverse axis `y` (a marginal). Each step's energy is split across the depth
bins it spans (as in decision 0010) at the lateral bin of the step's mid-`x`
(the within-step lateral motion is ~micrometres and negligible). The depth dose
is the `x`-marginal; `sigma_x(z)` is the second `x`-moment per depth slice. The
scoring accumulator is float64 on the device (decision 0005).

## Acceptance criteria (fixed before the comparison), milestone V1 (lateral part)

The primary check is against the **Fermi-Eyges oracle** built from the *same*
scattering power (a tight sampling/hinge-correctness check); a secondary check
is against published measured/Geant4 spreads.

| check | criterion |
|---|---|
| theta0 and pv vs the Highland/Lynch-Dahl hand values | within 0.5 % |
| MC sigma_x(z) vs the Fermi-Eyges oracle at 0.5R and 0.8R, 150/200 MeV | within 3 % (primary) |
| MC sigma_x(0.8R) vs published (~2.4 mm at 150 MeV, ~3.9 mm at 200 MeV) | within 8 % (secondary) |
| energy conservation (reference float64 / Warp float32) | 1e-9 / 5e-5 |
| projected range vs the CSDA range (detour factor) | shortened by < 0.1 % (not a bug) |
| reference vs Warp CPU/CUDA: sigma_x(0.8R) | within the batch statistical error |
| depth dose (x-marginal) with scattering on vs the DEV-005 Bragg peak | consistent (R80, peak) |

## Rationale

The bracket-1 scattering power resolves the thin-step log pathology and, being
identical in the step model and the Fermi-Eyges oracle, makes the primary gate
a clean test of the *sampling and hinge* implementation (the classic
factor-2 projected-vs-space or missing-1/3-hinge bugs would show as tens of
percent). The random hinge is required for the correct lateral displacement;
the two-plane sampler is the cheapest correct GPU form (one Box-Muller, no
azimuth trig). Validating on `sigma_x(0.8R)` follows the research
recommendation: it is sensitive to the scattering power but robust against the
end-of-range divergence and single-scattering tail.

## Expected tradeoffs

- Bracket-1 Highland (no log term): the log correction is ~few percent over the
  full path; dropping it is covered by the 8 % physics tolerance and is the
  clean per-step choice. A total-pathlength log term is a possible later
  refinement.
- Gaussian central MCS only: correct for the lateral-`sigma` core; the far
  penumbra tail (large-angle single scattering) is not modelled and is out of
  scope for V1.
- The lateral position of a step's deposition uses the step mid-`x`
  (micrometre error); `sigma_x` accumulates over many steps from the
  hinge-correct positions, so this is negligible.
- MCS shortens the projected range by the detour factor (< 0.1 %); documented,
  not a bug.
- 2-D scoring adds an `n_depth x n_lateral` accumulator; fine at these sizes.

## Validation strategy

`tests/ionmc/test_scattering.py` (theta0/pv values, MC sigma_x vs Fermi-Eyges
and vs targets, energy conservation, detour factor, cross-backend, depth-dose
marginal consistency) and `validation/v1_lateral_scattering.py` on the host
runner for the exact SHA (Warp CPU and CUDA), recorded below.

## Later validation outcome

To be filled in from the DEV-006 host run.
