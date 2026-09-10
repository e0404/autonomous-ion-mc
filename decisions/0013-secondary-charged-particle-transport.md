# 0013 — Secondary charged-particle transport from nonelastic reactions

- Status: accepted
- Date: 2026-09-10
- Task: DEV-008
- Affects: physical accuracy, architecture, validation strategy, scientific interpretation

## Problem

DEV-007 (decision `0012`) removes primaries at the macroscopic nonelastic
nuclear rate, deposits a local fraction of the reacting primary's energy at the
vertex, and books the rest to an abstract **escaping** channel. That escaping
energy is not all lost: about half of it is carried by **secondary protons**
that travel on and deposit a broad low-level dose (the "nuclear halo/plateau"),
lifting the entrance/plateau dose by ~1-2 % and the pre-peak dose to ~5-10 %
of the local dose (Paganetti 2002). This task replaces the single escaping
bucket with a physical three-way partition and transports the secondary
protons, **closing milestone V2** (nuclear interactions for protons in water).

This decision fixes, before the first comparison: the per-reaction energy
partition, the secondary-proton multiplicity/spectrum/angle model, the
transport architecture (how secondaries are generated and tracked while keeping
cross-backend bit-parity), what is deferred, and the acceptance targets.

## Context and evidence

Research (Claude physics-researcher subagent, 2026-09-10; sources below) and the
existing DEV-007 machinery established:

- **Energy partition of a p+16O nonelastic reaction** (roughly constant over
  100-250 MeV). Of the reacting proton kinetic energy `E`: secondary protons
  ~0.50 (0.45-0.60 across sources), heavy charged fragments (alphas, d, t, 3He,
  recoils) ~0.10-0.12, neutrons ~0.22, prompt gammas ~0.03, binding/Q ~0.15.
  Charged (local + transported) ~0.60; neutral + binding ~0.40. Whole-track
  cross-check: at 177 MeV ~20 % of protons react and the nonelastic channel
  takes ~12 % of the beam energy (Gottschalk 2015). Reactions are weighted to
  the high-energy entrance/plateau (representative reaction energy ~96 MeV on a
  177 MeV track), not the peak.
- **Dominant channels**: proton knockout `16O(p,2p)15N` (Q = -12.1 MeV),
  neutron knockout `16O(p,pn)15O` (Q = -15.7 MeV) (Gottschalk 2015).
- **Secondary-proton multiplicity** `nu_p(E)` (transportable protons above a
  ~1-2 MeV cut): ~0.6-0.8 at 100 MeV rising to ~1.2-1.6 at 250 MeV. A
  defensible linear form `nu_p(E) = 0.5 + 0.004 (E/MeV)`.
- **Secondary-proton spectrum**: two components. Evaporation (compound nucleus),
  `dN/dE ~ E exp(-E/T)`, nuclear temperature `T ~ 2 MeV`, isotropic, short
  range; cascade/knockout, ~10 MeV to ~E, forward-peaked, broad (a quasi-elastic
  bump near E and near E/2 for the two (p,2p) protons).
- **Fippel-Soukup (2004) fast-MC model family** (VMCpro; and as used by
  MCsquare 2016 and FRED 2017): on a nonelastic event the primary is removed;
  **secondary protons are transported as primaries**; **deuterons, tritons,
  alphas and heavy recoils are deposited locally** (range << voxel); **neutrons
  and prompt gammas are not transported, their energy is treated as locally
  lost** (exactly the escaping channel). Outgoing protons are "almost isotropic
  in the CM system", sampled in CM and boosted to lab, or sampled from the ICRU
  63 double-differential tables.
- **Dose magnitudes** (acceptance anchors, Paganetti 2002): secondary protons
  deliver up to ~10 % of the total dose proximal to the Bragg peak and ~1-2 % at
  entrance; d/t/3He/alpha together < 0.1 % of total dose; neutron dose distal to
  an SOBP < 0.05 %. These quantify both the acceptance target and the safe
  deferrals.

## Decision

1. **Three-way per-reaction partition** (replacing the DEV-007 lumped local
   fraction, only when secondary transport is enabled). Of the reacting
   primary's kinetic energy `E` at the vertex:
   - **heavy fragments, local**: `f_heavy = 0.12` deposited at the vertex bin;
   - **secondary protons, transported**: `f_p = 0.50` carried by sampled
     secondary protons;
   - **truly escaping**: `f_esc = 1 - f_heavy - f_p = 0.38` (neutrons, gammas,
     binding) booked to `escaped_mev`.
   These are named engine constants, overridable. With secondary transport
   *off*, DEV-007 is reproduced exactly (`f_local = 0.30`, escaping `0.70`, no
   second pass, no extra RNG draws).

2. **Secondary sampling** (`ionmc.physics.secondaries`, host-side numpy, used
   identically by both backends). Per reaction of weighted energy `w E`:
   - draw the count `N` from the multiplicity `nu_p(E)` (rounded/│Poisson-like
     via the counter RNG); if `N = 0` the whole `f_p w E` folds into the
     escaping channel (documented);
   - sample `N` proton energies from the evaporation+cascade mixture, then
     **renormalize their sum to `f_p w E`** so the event energy balance stays
     closed to round-off (the decision-0012 `deposited + escaped = energy_in`
     gate is preserved exactly);
   - energies below a transport cut (`1 MeV`, or the table floor) are deposited
     locally at the vertex rather than transported (folds the short-range
     evaporation tail into the local channel).

3. **Forward emission for the 1-D depth dose.** Secondaries are emitted along
   `+z` from the vertex and transported by the existing CSDA depth-dose engine.
   The dose-dominant cascade protons are genuinely forward; the low-energy
   isotropic evaporation protons are short-range and deposit near the vertex
   regardless of direction, so the depth distribution is insensitive to their
   angle. Full 3-D angular emission and the lateral nuclear halo belong to the
   3-D scattering path and are deferred (they do not affect the integral depth
   dose materially).

4. **Two-pass architecture preserving cross-backend parity.** Each primary
   reacts at most once (catastrophic removal), so the transport drivers emit one
   reaction record per history (vertex depth, weighted residual energy) with no
   extra RNG draws. Pass 1 transports primaries (reference or Warp) and returns
   the records; the host generates the secondary `ParticleState` deterministically
   from the records using counter-based streams keyed by a stable reaction
   ordinal (primary history index); pass 2 transports the secondaries through the
   same driver and its dose is added to the grid. Because DEV-007 already
   guarantees the reference and Warp paths remove the **identical** primary set
   (decision `0001`), and secondary generation is host-side and deterministic,
   the secondary dose is bit-reproducible across backends.

5. **Secondaries are transported without their own nuclear removal** (pure EM:
   CSDA + energy-loss straggling; scattering off on the 1-D path). Secondary-
   secondary (tertiary) reactions are a sub-percent effect and are deferred; this
   also keeps the secondary energy fully deposited, so the budget closes cleanly.

## Deferred (with quantitative justification)

- Explicit deuteron/triton/alpha/recoil transport: deposited locally
  (< 0.1 % of dose; Paganetti 2002).
- Neutron transport: dropped as escaping (distal neutron dose < 0.05 %).
- Prompt-gamma transport: dropped as escaping (needed for range verification,
  not dose).
- Tertiary (secondary-secondary) nuclear reactions.
- Full intranuclear-cascade / pre-compound modelling and the ICRU-63 double-
  differential *tabulated* sampling: the analytic evaporation+cascade surrogate
  reproduces the pristine Bragg curve within a few percent; the tabulated path
  is deferred to the same openly-licensed TENDL-2021 follow-up already planned
  in decision `0012` (ICRU 63 is copyrighted and is not vendored).

## Consequences

- The simulated depth dose now carries the secondary-proton plateau, so the
  absolute peak-to-entrance ratio drops toward measured pristine Bragg curves.
  The peak-to-entrance ratio is geometry- and convention-dependent, so it is not
  hard-coded as a tight gate; the acceptance is the **secondary-dose fraction vs
  depth** and the lifted entrance plateau, plus the unchanged survival gate.
- `DepthDoseResult` gains the secondary dose (its own array) and a secondary
  count, so the secondary contribution is auditable per depth.
- The DEV-007 depth-dose engine gains a second transport pass; the 3-D
  scattering path is unchanged (nuclear transport there remains future work).

## Acceptance targets (validation `v2_secondary_transport.py`)

- **Secondary-dose fraction**: ~1-2 % of the local dose at entrance and a few
  percent (~2-10 %) of the total dose at 150 and 200 MeV (Paganetti 2002).
- **Secondary plateau shape**: the secondary dose is a broad plateau, not a
  peak; its fraction at the Bragg peak is well below the entrance value (the
  primary peak dominates there). This is the clean, baseline-independent check
  (a comparison against the DEV-007 curve is confounded because that curve uses
  the lumped local fraction `f_local = 0.30` rather than `f_heavy = 0.12`).
- **Energy bookkeeping**: `deposited + escaped = energy_in` to < 1e-6 (Warp
  float32) and < 1e-9 (reference), now with the secondary pass included.
- **Secondaries-off regression**: identical to the DEV-007 baseline (no second
  pass, no extra draws, unchanged depth dose and escaped energy).
- **Cross-backend**: reference vs Warp CPU and (host) CUDA reproduce the same
  secondary set and a cumulative total-depth-dose difference within the
  decision-`0001` tolerances.

## Sources

- Paganetti, Phys. Med. Biol. 47 (2002) 747 (secondary-particle dose fractions).
- Fippel & Soukup, Med. Phys. 31 (2004) 2263 (VMCpro fast-MC nuclear model).
- Souris, Lee & Sterpin, Med. Phys. 43 (2016) 1700 (MCsquare).
- Schiavi et al., Phys. Med. Biol. 62 (2017) 7482 (FRED).
- Gottschalk et al., Phys. Med. Biol. 60 (2015) 5627 / arXiv:1412.0045 (halo,
  energy accounting).
- Holmes, Fast Monte Carlo Dose Calculation in Proton Therapy review,
  arXiv:2404.14543 (2024) (practical parameterizations, deferrals).
