# V3-003 acceptance criteria (frozen before any LV/HR result)

Status: frozen with task V3-003A, before any local-validation (LV) or
host-runner (HR) result of the transport engine exists. Later changes
require a new dated entry below with the scientific reason (decision 0039).
Criteria apply to the proton electromagnetic transport engine of decision
0039; nuclear interactions, LET and persistence are covered by later plans.
Tiers: CI = GitHub (small, Warp CPU); LV = local exact-SHA validation on the
CPU; HR = controlled host runner with CUDA. Evidence classes follow
`experiment/v3/PROTOCOL.md`.

Lead amendments at freeze time: `E_cut` is 2 MeV (table floor 1 MeV/u, see
decision 0039), so T2/T3 use `R_csda(E_cut)` with that value; the default
batch count is 20.

## Criteria
Tiers:
- **CI:** GitHub, small, warp CPU.
- **LV:** local exact-SHA validation on the CPU.
- **HR:** host runner with CUDA. It has only numpy, warp and pytest and no network, so every test uses the offline analytic `StoppingSource`. Chi-square and normal quantiles come from `statistics.NormalDist` plus Wilson–Hilferty, without scipy.

| # | Check | Criterion (pass iff) | Evidence class | Tier |
|---|---|---|---|---|
| U1 | Every shared func, Python scope vs Warp CPU kernel (f32, f64), 10⁴ args incl. edges | Decision 0001 classes: f64 rtol 1e-12 / atol 1e-14, f32 1e-6 / 1e-6; `dda_next` reasons and indices identical; NaN fails | backend parity | CI |
| U2 | Table round trip \|Rinv(R(E))−E\|/E and monotonicity of R and Rinv | ≤ 1e-5 (f64 build), strictly monotone | self-consistency | CI |
| U3 | X_S from the formula for water, Be, Al, Cu, Pb | within 0.3% of 46.88 / 92.60 / 28.75 / 14.62 / 6.62 g/cm². This also settles ln vs log10 in eq. XS. | source-model reproduction⁴ (Gottschalk's own X_S table) | CI |
| U4 | θ_dM(x) by numpy quadrature of T_dM along CSDA E(x), 158.6 MeV, x/R1 ∈ {0.01, 0.1, 0.5, 0.9}, for the materials V3-002 provides | within 1.5% of Gottschalk's θ_Hanson·(1+dM%) (theta0Single table). Residual is from the range-table I-value. | source-model reproduction³ (shared T_dM model and coefficients; not independent) | LV¹ |
| U4b³ | Same θ_dM(x) as U4 compared directly with Gottschalk's tabulated θ_Hanson (Molière/Fano/Hanson theory, the quantity T_dM was fitted to; independent of the dM fit) | \|θ_dM/θ_Hanson − 1\| ≤ 4.5 % at every frozen point, derived before measurement from the table's own max \|dM%\| = 2.74 % (x/R1 ≥ 0.01) compounded with the U4 bound 1.5 % | source-model reproduction⁴ (θ_Hanson is the calibration target of T_dM and the bound is implied by U4; not independent) | LV¹ |
| T15³ | **Independent-MC lateral check:** 150 MeV protons in water, all physics on, lateral σ of the dose profile at z/R ∈ {0.5, 0.9} (0.2 mm lateral bins, 1 mm slabs) versus a TOPAS/Geant4 11 run of the same geometry (3-D dose grid, ≥ 10⁵ primaries, ≥ 2 seeds) from the V3-010 reference service⁶ | \|σ_ionmc/σ_TOPAS − 1\| ≤ 3 % and \|z\| < 3 using both standard errors; result labelled by TOPAS physics list and I-value (78 eV) | independent Monte Carlo | LV (after V3-003B and the V3-010 3-D dose case) |
| U5 | **Step independence of the MCS integrator (deterministic):** Σ_steps T(E_mid)·s along CSDA in water, 150 MeV, s ∈ {0.01, 0.1, 0.5, 1, 2, 5} mm, compared with adaptive quadrature at x/R1 ∈ {0.05, 0.25, 0.5, 0.9} | \|Δθ²/θ²\| ≤ 2e-3 for every s at x/R1 ≥ 0.05. Negative control computed analytically: per-step Highland varies ≥ 5% over the same s range, so the test has power. | self-consistency | CI |
| U6 | RNG tests (a)–(e) of §3 | exact | self-consistency | CI |
| T1 | **Trajectory parity**, python vs warp-cpu f64 trace: K = 16 (CI) and K = 256 (LV) histories, 100 MeV (CI) and 150 MeV (LV) in a water box, all physics on | Discrete sequences (idx, step reason, c2 count, attempt counts) identical for 100% of steps. Continuous state per step within rtol = atol = 1e-10. Any branch flip fails. | backend parity | CI/LV |
| T2 | **Deterministic CSDA** (straggling off, MCS off): end depth for 100/150/200 MeV | R_csda(E0) − R_csda(e_cut) − 1e-4·R ≤ z_end ≤ R_csda(E0) + 1e-4·R (f32), using the project table | self-consistency | CI (100 MeV) / LV |
| T3 | Finite-slab escape, 100 MeV through 30 mm water, deterministic | E_escaped per primary = Rinv(R(100) − 3.0 g/cm²) within 1e-4 relative | self-consistency | CI |
| T4 | Energy balance: initial = Σgrid + Σoutside + cutoff + escaped + truncated + unaccounted, all accumulated independently | rel ≤ 1e-5 (f32), 1e-12 (f64 and python); counters 0 | conservation | CI |
| T5 | R80 of the IDD (0.1 mm depth bins, straggling on, MCS off), 10⁵ histories | \|R80/R_csda − 1\| ≤ 0.2% at 100/150/200 MeV | related model (Bortfeld R80 ≈ R0) | LV |
| T6 | σ_R from end-depth diagnostics (straggling on, MCS off) | within 5% of σ_R² = ∫(dσ²/dx)/S³ dE (numpy, same Bohr variance); within 10% of 0.012·R^0.935 cm | self-consistency + related model | LV |
| T7 | Lateral variance of deposits in 1 mm slabs (0.2 mm lateral bins, Sheppard-corrected) at z/R ∈ {0.25, 0.5, 0.75, 0.9}, MCS on, straggling off, 150 MeV | within 2% of Fermi-Eyges A2(z) = ∫(z−z')²T_dM dz'; within 3% of the generalised-Highland y_rms (related model) | self-consistency⁴ (Fermi-Eyges with the implemented T_dM) + related model | LV |
| T8 | **MCS step independence in transport:** exit θ_rms of a 0.5·R1 water slab and T7 σ at 0.9R, s_max ∈ {0.1, 0.5, 1, 5} mm, 10⁶ histories warp-cpu | pairwise θ_rms ≤ 0.5%, σ ≤ 1%, and θ_rms within 0.5% of the U5 quadrature | self-consistency (falsification) | LV |
| T9 | IDD step independence, s_max ∈ {0.1, 0.5, 1} mm, f_E ∈ {0.005, 0.02} | \|ΔR80\| ≤ 0.1 mm; IDD χ² (dose > 1% of max) p > 0.001 | falsification | LV |
| T10 | Rotation invariance: beam along +x, +y, +z, −z, and (1,1,0)/√2 and (1,1,1)/√3 in a 200 mm cube; IDD versus projected depth | permutations: χ² p > 0.001; obliques: \|ΔR80\| ≤ 0.3 mm and total energy equal within 3σ; zero counters | falsification | LV |
| T11 | DDA adversarial cases: source exactly on planes, edges and corners, directions with a zero component, grazing entry | no stall or truncation; T4 holds | falsification | CI |
| T12 | **Statistical parity:** python (4·10³, spawn pool), warp-cpu (10⁶), warp-cuda (10⁶) at 150 MeV with all physics on; also the warp f32 vs f64 bias probe | per-voxel z = Δ/√(σa²+σb²) on the IDD and on lateral profiles at 3 depths with dose > 1% of max: χ² p > 0.001 and max\|z\| < Bonferroni Φ⁻¹(1−0.001/2n); scalars (R80, total deposit, σ_lat at 0.5R) \|z\| < 3.5. CI runs a reduced version: python 200 vs warp-cpu 2·10⁴ at 70 MeV. | backend parity | CI (reduced) / LV / HR (CUDA) |
| T13 | Partition invariance: 1 vs 3 workers, and CUDA chunk sizes 2¹⁰ vs 2¹⁸ | counters and tallies identical; edep within rel 1e-5 (f32) or 1e-12 (f64) | self-consistency | CI (CPU) / HR |
| C1 | Fail-closed dispatch: each rule in §1 raises before transport; requested and effective configs are present | | capability | CI |
| T14² | **Voxel-boundary scattering bias (grid size and alignment):** 150 MeV in water, MCS on, straggling off, 10⁶ histories warp-cpu; transport voxel size {0.5, 1, 2, 5} mm, each also with the phantom grid shifted by half a voxel along the beam and laterally; observables: exit θ_rms of a 0.5·R1 slab and the T7 lateral σ at z/R ∈ {0.5, 0.9} (scoring grid fixed at 0.2 mm lateral bins, independent of the transport voxels) | pairwise θ_rms ≤ 0.5 %, σ_lat ≤ 1 % across all voxel sizes and shifts; θ_rms within 0.5 % of the U5 quadrature; σ_lat within 2 % of Fermi-Eyges A2(z) (as T7). Negative control: forcing the hinge angle to be sampled for the truncated length only (diagnostic switch) must change θ_rms by < 0.5 % at 1 mm voxels — if it changes more, the approximation is material and the default must switch | falsification (boundary approximation) | LV (V3-003B) |

CI ≤ ~3 min. This file is that frozen table.


## Amendments

² Added 2026-10-03 after the first Codex review of V3-003A (head de212db) and before any T14
measurement exists: decision 0039 acknowledges that the hinge angle is sampled for the planned
step while the second leg may be truncated at a voxel plane, which overestimates scattering on
boundary steps; U5, T8 and T11 do not measure this. T14 bounds it with a grid-size and
half-voxel-alignment refinement using angular and lateral observables (requirement V2-NUM
boundary alignment). The numbers are frozen now; the test runs on the Warp CPU backend
(V3-003B), where 10⁶ histories are feasible.

¹ Venue amendment (2026-10-03, after the table was materialized, before the test was first run against it):
the theta0Single table is parsed at run time from the materialized arXiv:0908.1413 source
(provenance cache, sha256 67fb1478…), which is not in Git; the U4 test therefore skips in CI and is
executed locally (LV) and on the host runner. The acceptance number (1.5 %) is unchanged. The table
contains Be, Al, Cu and Pb blocks only (no water), so U4 covers those four V3-002 materials; water
is covered by the cache-independent generalised-Highland cross-check (bounded at 6 %, not an
acceptance criterion) and by the X_S comparison with Gottschalk's Table tbl:LS (within 0.5 %).

³ Classification amendment (2026-10-03, after the second Codex review of V3-003A, head f4aad1f):
U4 compares against θ_Hanson·(1+dM%), which reconstructs Gottschalk's own differential-Molière
result — the model, coefficients and calibration paper that the engine implements — so U4 is a
source-model reproduction check, not independent theory. U4b (direct comparison with the
Molière/Fano/Hanson column, which T_dM was fitted to) and T15 (independent Monte Carlo lateral
spreading) are added as the independent scattering evidence. The U4b tolerance is derived from
the paper's table (max |dM%| 2.74 % for x/R1 ≥ 0.01, Pb at 0.9 R1 = −2.51 %) and the frozen U4
bound, not from any ionmc result; the U4b numbers had not been computed when this was frozen.

⁴ Classification amendment (2026-10-03, after the third Codex review of V3-003A, head dc88004):
U3, U4, U4b and the X_S comparison all reproduce Gottschalk (2010) — the paper whose formulae and
coefficients the engine implements — from that paper's own tables, and T7 integrates the
implemented T_dM in Fermi-Eyges theory; none of them is independent evidence for the scattering
physics. They are retained as source-model reproduction and self-consistency checks with their
frozen numbers unchanged. The only independent scattering qualification in this plan is T15
(independent Monte Carlo), still pending; the requirement ledger must not mark scattering
requirements as independently evidenced before T15 passes.

⁶ Clarification (2026-10-04, task V3-010B, before any T15 reference run was launched): the 0.2 mm
lateral bins of T15 refer to ionmc's own scoring grid; the TOPAS reference scores the 3-D dose on
0.5 mm lateral × 1 mm depth bins over the full ±60 mm field, and both σ estimators apply the
Sheppard bin-width correction on the same lateral window (full field and a ±20 mm window are both
reported). The TOPAS case carries two 3-D scorers: all particles, and primary-generation only
(TOPAS generation filter). EM-only backends (V3-003B, no nuclear interactions) are compared with the
primary-generation scorer — the residual difference is Geant4's hadron-elastic deflection of
primaries, which ionmc lacks until V3-005 and which is therefore part of the measured deviation,
not an allowance; the all-particle scorer is the T15 reference once nuclear interactions exist
(V3-005). Tolerances are unchanged.

⁸ Correction of clarification ⁶ (2026-10-04, task V3-010B, after Codex review REVIEW-87d6aaa0e15a4866916dfefac62485c4
and before any T15 grading): the primary-generation-filtered TOPAS scorer excludes dose deposited by
electromagnetic secondaries (delta electrons), whereas ionmc deposits the entire electronic energy
loss locally — a different observable. The T15 reference for EM-only backends is therefore a
matched EM-only, all-particle TOPAS run (electromagnetic physics modules only, same geometry,
scorers and seeds: case `topas/proton-water-150mev-lateral-emonly`); the generation-filtered scorer
of the all-physics lateral case is informative only. The remaining known difference between that
reference and ionmc is Geant4's transport of delta electrons above its production cut versus local
deposition; it is part of the measured deviation, not an allowance. The all-physics, all-particle
scorer is the T15 reference once nuclear interactions exist (V3-005). Tolerances are unchanged.

¹⁰ Window designation for T15 (2026-10-04, task V3-010B, after Codex review REVIEW-0aa3fa48f2f54b4f81e5c3a2cab13da3
and before any T15 grading): the T15 acceptance observable is the Sheppard-corrected lateral σ of
the dose profile within the ±20 mm window around the beam axis (1 mm slabs at z/R80 = 0.5 and
0.9), computed identically for ionmc and the reference; the full-field (±60 mm) second moment is
reported alongside as informative only. Reason: the full-field moment is dominated by the far
tails, and ionmc's Gaussian multiple-scattering model has no single-scattering tail, so its
full-field σ is expected to fall below the reference — that expected deficit is a documented
model limitation (to be probed in V3-011), not a T15 pass/fail quantity. The frozen tolerance
(3 % and |z| < 3) applies to the ±20 mm window at both depths; both depths must pass.
