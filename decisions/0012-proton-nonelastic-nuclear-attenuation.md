# 0012 — Proton nonelastic nuclear attenuation and local deposition

- Status: accepted
- Date: 2026-09-10
- Task: DEV-007
- Affects: physical accuracy, architecture, validation strategy, scientific interpretation, data-layer policy

## Problem

DEV-004/005/006 transport a proton by continuous slowing down with energy-loss
straggling and multiple Coulomb scattering. Every primary survives to its CSDA
range, so the simulated integral depth dose carries **too much dose in the
Bragg peak relative to the entrance**: real protons are removed along the track
by **nonelastic nuclear reactions** on target nuclei (mostly oxygen-16 in
water), which convert a primary into short-range recoils and heavy fragments
(deposited locally) plus secondary protons, neutrons and gammas (transported or
carried away). About **20 % of 150 MeV protons and 27 % of 200 MeV protons
react before reaching the Bragg peak** (Paganetti 2002; Gottschalk). This task
adds catastrophic nonelastic removal of primaries with a local-deposition
fraction, and books the remainder to an escaping/deferred channel that DEV-008
(secondary transport) will consume. It opens milestone **V2**.

This decision fixes, before the first comparison: the cross-section model and
its data provenance, which reaction channels remove a primary, the per-step
removal algorithm and its RNG-stream placement, the local/escaping energy
split, and the acceptance targets.

## Context and evidence

Research (Claude physics-researcher subagent, 2026-09-10) and in-repo numerical
checks established:

- **Nonelastic cross section on oxygen-16.** The proton-oxygen nonelastic cross
  section rises from a ~7 MeV threshold, peaks near ~550 mb around 20 MeV, and
  settles to a mild plateau of ~340–400 mb over 100–250 MeV (ICRU Report 63
  shape; Geant4/ICRU compilations). Water contains `n_O = 3.343e22` oxygen
  atoms/cm^3 (Avogadro x mass-fraction/atomic-weight x density), so the
  macroscopic removal cross section is `Sigma = n_O sigma_nonel(E) ~ 0.012–
  0.013 /cm` over the therapeutic plateau, consistent with Gottschalk's
  `delta ~ 0.012 /cm` rule of thumb.
- **Only the oxygen channel removes a primary.** The hydrogen channel at these
  energies is proton-proton *elastic* scattering, which **deflects** the primary
  (a large-angle event that belongs with lateral spread and secondary-proton
  production) rather than **removing** it from the primary population. Counting
  it as catastrophic removal would roughly double the attenuation and overshoot
  the measured survival. It is therefore excluded here and deferred to secondary
  production.
- **Published primary survival to the Bragg peak.** `~0.80` at 150 MeV and
  `~0.73` at 200 MeV (Paganetti, Phys. Med. Biol. 47 (2002) 747; Gottschalk,
  "Passive Beam Spreading"). These are the primary-fluence attenuation factors
  and are the acceptance target.
- **Local vs. escaping energy split.** Of the reacting primary's kinetic energy,
  a fraction is deposited within ~a millimetre of the vertex (nuclear recoils,
  alphas and heavier fragments, short-range secondary protons), and the rest is
  carried by longer-range secondary protons, neutrons and photons. Published
  fragment-energy accounting puts the promptly-local fraction near ~0.3 of the
  primary energy at therapeutic energies; the remainder is exactly the pool a
  secondary-transport stage must conserve. We take **`f_local = 0.30`** as the
  DEV-007 constant, with the remaining `0.70` booked to an escaping/deferred
  channel. Both `f_local` and the escaping pool are explicit so DEV-008 can
  replace the escaping term with transported secondaries without changing the
  bookkeeping contract.

In-repo verification (float64 reference, 4000 histories, MCsquare/PSTAR water
table): the analytic `sigma_nonel(E)` gives `Sigma = 0.0127 /cm` at 150 MeV and
`0.0120 /cm` at 200 MeV, integral reaction fractions **0.183 / 0.279**, i.e.
primary survival **0.817 / 0.721**, matching the published 0.80 / 0.73 within
tolerance. Energy conservation `deposited + escaped = energy_in` holds to
float64 round-off (reference) and to ~4e-7 (Warp float32). The Warp CPU path
reproduces the reference reaction count and escaped energy **exactly** (aligned
counter-based RNG) with a cumulative depth-dose difference of 1.2e-6.

## Decision

1. **Analytic cross-section model, not a vendored table.** The nonelastic
   cross section on oxygen is an analytic parameterization of the ICRU-63 shape,
   written in the shared-source math namespace so the Warp kernel and the
   float64 reference evaluate it identically
   (`ionmc.physics.nuclear.nonelastic_cross_section_oxygen`):
   a linear threshold turn-on from 7 to 20 MeV, a plateau declining from 0.40 b
   at 100 MeV to ~0.34 b at 250 MeV, and a Gaussian resonance bump of 0.15 b
   centred at 20 MeV. This satisfies the REQUIREMENTS "analytical physics layer"
   MUST and avoids depending on a copyrighted tabulated report; a TENDL-2021 /
   ICRU-63 **tabulated** path is deferred (see "Data-layer policy" below).

2. **Oxygen-only macroscopic removal.** `Sigma(E) = n_O sigma_nonel(E)`, with
   `n_O` derived from the medium's composition and density in the engine
   (`AVOGADRO x atoms_per_gram('O') x density`; zero, disabling nuclear removal,
   for oxygen-free media). Hydrogen is excluded (elastic, not removal).

3. **Per-step catastrophic removal.** Behind a `nuclear: bool = False` flag on
   `TransportEngine`, each *alive* step draws **one** uniform and tests
   `xi < P`, with `P = Sigma(E_step-entry) x (s/10)` the thin-step reaction
   probability over the actual (post-geometry-clamp) step length `s` [mm]. The
   step-entry energy is used (the cross section varies slowly over one step). On
   a reaction the primary is removed: it deposits `w x f_local x E` locally at
   the vertex bin, books `w x (1 - f_local) x E` to the escaping accumulator,
   increments the reaction counter, is marked `Status.REACTED`, and stops. The
   thin-step form `Sigma x s` (vs. `1 - exp(-Sigma s)`) is accurate to better
   than 0.1 % for the ~1 mm steps used here.

4. **RNG-stream placement fixes cross-backend parity.** The nuclear uniform is
   drawn **after** the straggling normal and **unconditionally per alive step**
   when `nuclear=True` (the probability is zero below the 7 MeV threshold, so
   the draw never triggers there). This keeps the reference and Warp counter
   streams bit-aligned, so the two backends produce the identical reaction set
   (decision 0001). With `nuclear=False` no extra draw is made, so DEV-004/005/
   006 streams and results are reproduced exactly.

5. **Energy bookkeeping is closed.** `DepthDoseResult` gains `escaped_mev` and
   `n_reactions`; `energy_balance` becomes `(deposited + escaped - in) / in`.
   The escaping pool is a first-class output, not discarded, so it is auditable
   now and consumable by DEV-008.

6. **`f_local = 0.30`, escaping = 0.70**, both as named engine defaults
   (`DEFAULT_NUCLEAR_LOCAL_FRACTION`), overridable per engine.

## Data-layer policy

The DEV-007 research recommended against vendoring the ICRU-63-derived MCsquare
oxygen cross-section table (it is derived from a copyrighted ICRU report and its
redistribution terms are not Apache-compatible). The oxygen dataset registration
that had been drafted was reverted; DEV-007 ships the analytic model only. A
future tabulated path should use an openly-licensed evaluation (e.g. TENDL-2021,
which is public-domain) registered through the decision-0007 `DatasetSpec`
mechanism, and validate the analytic model against it. This is deferred, not
rejected.

## Consequences

- The simulated Bragg peak-to-entrance ratio drops toward measured values, but
  will still read **too high in absolute terms until DEV-008 transports the
  escaping secondary protons** (which deposit a broad low-level dose plateau).
  DEV-007's depth-dose gate therefore checks *primary attenuation and energy
  bookkeeping*, not the absolute peak-to-entrance ratio.
- `Status.REACTED = 4` is added (3 already marks truncation). Reacted primaries
  are excluded from the stopped-range statistics (they did not reach their CSDA
  range).
- Nuclear removal is wired into the **depth-dose** path (1-D) only. The 3-D
  scattering path (`run_scattering`) keeps EM-only transport for now; nuclear
  removal there is a small follow-up once secondary lateral transport exists.

## Acceptance targets (validation `v2_nuclear_attenuation.py`)

- **Primary survival to the Bragg peak**: `0.80 +- 0.03` at 150 MeV,
  `0.73 +- 0.04` at 200 MeV.
- **Energy bookkeeping**: `deposited + escaped = energy_in` to <1e-6 (Warp
  float32) and <1e-12 (reference).
- **Statistical self-consistency**: the Monte Carlo reaction fraction agrees
  with the analytic `1 - exp(-integral Sigma dl)` along the mean track within
  the batch statistical uncertainty.
- **`nuclear=False` regression**: identical results to the DEV-005 baseline
  (zero reactions, unchanged depth dose).
- **Cross-backend**: reference vs. Warp CPU and (host) CUDA reproduce the same
  reaction set and a cumulative depth-dose difference within the decision-0001
  tolerances.
