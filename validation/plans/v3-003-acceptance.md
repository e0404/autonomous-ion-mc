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
| T9-CI⁷ | **Scoring aliasing probe (deterministic):** single-voxel water box, 150 MeV, straggling and MCS off, 2·10⁴ histories, IDD with 1 mm bins for s_max ∈ {0.33, 0.5, 0.9, 1.0} mm and 2 mm bins for s_max = 2.0 mm, each compared with the 0.1 mm-step run | max \|IDD deviation\| ≤ 2e-3 in 20–120 mm and ≤ 1e-2 in 125–140 mm for every s_max; pre-fix midpoint scoring gave 32–80 % (plateau) and 62 % (125–140 mm) — those numbers stay in the test as the negative control | falsification (numerical) | CI |
| T13⁷ (edep) | Partition and chunk invariance of the deposit grid with int64 fixed-point accumulators | edep bit-identical across 1 vs 3 workers (CPU) and chunk sizes 2¹⁰ vs 2¹⁸ (CUDA) within each precision; the former relative bounds (1e-5 f32 / 1e-12 f64) are reported but no longer the criterion | self-consistency | CI (CPU) / HR |
| T10 | Rotation invariance: beam along +x, +y, +z, −z, and (1,1,0)/√2 and (1,1,1)/√3 in a 200 mm cube; IDD versus projected depth | permutations: χ² p > 0.001; obliques: \|ΔR80\| ≤ 0.3 mm and total energy equal within 3σ; zero counters | falsification | LV |
| T11 | DDA adversarial cases: source exactly on planes, edges and corners, directions with a zero component, grazing entry | no stall or truncation; T4 holds | falsification | CI |
| T12 | **Statistical parity:** python (4·10³, spawn pool), warp-cpu (10⁶), warp-cuda (10⁶) at 150 MeV with all physics on; also the warp f32 vs f64 bias probe | per-voxel z = Δ/√(σa²+σb²) on the IDD and on lateral profiles at 3 depths⁵ with dose > 1% of max: χ² p > 0.001 and max\|z\| < Bonferroni Φ⁻¹(1−0.001/2n); scalars (R80, total deposit, σ_lat at 0.5R) \|z\| < 3.5. CI runs a reduced version: python 200 vs warp-cpu 2·10⁴ at 70 MeV. | backend parity | CI (reduced) / LV / HR (CUDA) |
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

⁵ Clarification (2026-10-04, before any T12 run, task V3-003B): the three lateral-profile depths of
T12 are z/R ∈ {0.25, 0.5, 0.9} (1 mm slabs, 0.2 mm lateral bins, a subset of the T7 depths), and
each statistical sample uses a distinct seed (shared seeds would correlate the samples and
invalidate the z statistic). Tolerances are unchanged.

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

⁷ Amendment (2026-10-04, task V3-003B, before any T9/T13 result is graded): two frozen probes
falsified implementation choices — (a) the T9 step-independence probe exposed point-sampling
aliasing of the midpoint deposit with the IDD bin edges (up to 80 % deviations for steps that are
not an integer fraction of the bin; details and numbers in decision 0039, "Later validation
outcome"); scoring is changed to path-length-proportional apportioning along both hinge legs and
the deterministic T9-CI probe is added with its tolerances frozen here; T9 itself is unchanged.
(b) The T13 chunk-invariance check failed on CUDA for float32 deposit grids (tallies and counters
identical, grid outside 1e-5; float32 atomic_add order depends on the chunking); decision 0037 is
amended to int64 fixed-point accumulators, and the T13 deposit criterion becomes bit-identity
within each precision (stronger than the former bounds). The failed float32 result is preserved
in the validation ledger (VAL-20261003-234609-14DE3F).

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

⁹ Statistical clarification (2026-10-04, task V3-003B, after Codex review REVIEW-23d3b5d7946b473d967409284fa816d4
and before any T9/T10/T12 LV/HR result is graded): the per-bin z statistics use the batch standard
errors of each sample and are marginally valid; the profile-wide χ² = Σ z² treats bins as
independent although bins are correlated (a history deposits in several bins and all bins share
batch fluctuations), so its p-value is calibrated by a batch-level studentized permutation test:
per-batch per-primary means m_b with batch sizes n_b from both samples are pooled, residuals
r_b = (m_b − μ̂)·√n_b are permuted across all batches (≥ 2000 permutations, recorded seed), the
samples are reconstructed as m*_b = μ̂ + r*_b/√n_b and χ² is recomputed; p = (1 + #{χ²* ≥ χ²})/(1 + n_perm).
The frozen thresholds (p > 0.001; Bonferroni max|z|) are unchanged; the Wilson–Hilferty p-value is
reported as informative only. The permutation assumes equal per-primary variance under the null.
Also: the T12 "reduced" label is evaluated per sample against its own frozen count (python 4·10³,
accelerated 10⁶); the T10 total-energy criterion uses the batch standard error (|z| < 3); runner
outputs produced with step selection are labelled subset and are conformant only through a
combined summary covering every step of the suite at the same SHA and environment; a run on a
dirty tree is refused, and a plain-snapshot run must be attested against the committed blobs.

¹⁰ Window designation for T15 (2026-10-04, task V3-010B, after Codex review REVIEW-0aa3fa48f2f54b4f81e5c3a2cab13da3
and before any T15 grading): the T15 acceptance observable is the Sheppard-corrected lateral σ of
the dose profile within the ±20 mm window around the beam axis (1 mm slabs at z/R80 = 0.5 and
0.9), computed identically for ionmc and the reference; the full-field (±60 mm) second moment is
reported alongside as informative only. Reason: the full-field moment is dominated by the far
tails, and ionmc's Gaussian multiple-scattering model has no single-scattering tail, so its
full-field σ is expected to fall below the reference — that expected deficit is a documented
model limitation (to be probed in V3-011), not a T15 pass/fail quantity. The frozen tolerance
(3 % and |z| < 3) applies to the ±20 mm window at both depths; both depths must pass.

¹¹ Definition of the T14 negative control (2026-10-04, task V3-003B, after Codex reviews
REVIEW-e9f5410b7d9c4cd0ad4540756901ae4d and REVIEW-94ca934293af49ceb0633ff55ec0d959, before any
T14 result is graded): the diagnostic switch applies "truncate first": before the angle is sampled,
the straight-line path with the pre-hinge direction over the planned step is tested against the
voxel planes; if it crosses one, the step is shortened to end exactly on that plane
(s_trunc = straight-line distance to the plane), otherwise s_trunc is the planned step. The hinge
angle is then sampled with the variance for s_trunc, the hinge is placed within s_trunc, both legs
are travelled without re-cutting, and the end point is snapped onto the plane along the cut axis
(lateral displacement kept); energy loss and scoring use s_trunc. By construction the travelled
hinge-path length equals the length used for the variance and the direction is untouched (tested);
the only approximation is the snap displacement δ of the end point along the cut axis, which is
second order in θ for planes normal to the beam but first order for grazing lateral planes. It is
recorded per step, summed vectorially per history, and the control verdict requires the RMS over
histories of |Σδ| at the slab exit to be below 1 % of the control sample's lateral σ at that exit
plane (positional bias of the control < 1 % of the lateral observable; quantiles of |Σδ| and of
per-step δ/s are reported as information). The frozen control quantity is unchanged:
|Δθ_rms(default − control)| < 0.5 % at 1 mm voxels. T14's angular observable is the raw pairwise
exit θ_rms of the same 0.5·R₁ slab (exact depth, no rounding) for every voxel size and shift,
obtained with a world exit plane at the slab thickness that is independent of the voxel grid
(grid shifts overhang the world instead of lengthening it); the frozen bound "θ_rms within 0.5 %
of the U5 quadrature" remains part of the T14 verdict; the quadrature ratio across shifts is
informative only.

¹² Statistical clarifications for T12 (2026-10-04, task V3-003B, written after the first full-scale
HR T12 rehearsal at f243b15 revealed two degenerate cases, before any T12 result is graded):
(a) a scalar whose batch standard errors are both below 1e-8 relative (the total deposit is fixed
by energy conservation) cannot be judged by a z statistic; it is compared with the deterministic
precision bound of the less precise sample (T4: 1e-5 relative for float32, 1e-12 for float64 and
python) and the rule applied is recorded; (b) the Bonferroni max|z| statistic includes only bins
where both samples have at least max(2, ⌈B/2⌉) batches with non-zero content (the same
defined-value rule as the ratio estimators of V3-004), because a batch standard error estimated
from a handful of non-zero batches is unreliable in profile tails; bins excluded by this rule are
listed in the result document. The χ² permutation statistic is unchanged. The float32-versus-
float64 probe of T12 keeps its full criteria: the rehearsal showed a significant float32 bias
(|ΔR80| ≈ 0.02 mm, IDD χ² at the permutation floor), which is treated as a finding to be fixed in
the engine (float64 energy bookkeeping inside the float32 kernel), not as an allowance.

¹³ Held-out qualification seeds (2026-10-04, task V3-003B, after Codex review REVIEW-2d89ea43c6784baea45cc8946c62850d):
the full-scale rehearsal of the LV/HR suites at f243b15 (seed base 20261004; validation record
VAL-20261004-064044-9EEB7A) exposed runner defects and the two degenerate T12 cases of ¹², and the
statistical rules were clarified after those outcomes were observed. The rehearsal is preserved as
contrary evidence and is not acceptance evidence. The qualification runs of every statistical step
(T8, T9, T10, T12 samples, T14, R1) use the independent seed base 20271004, recorded in the
environment, the step documents and each sample's metadata; samples with different seed bases are
never compared. All tolerances and the frozen history counts are unchanged.

¹⁴ Qualification seed base re-frozen (2026-10-04, task V3-003B, after Codex review REVIEW-bfc1d5ece8d74615ba1dea70f8c3aa72):
the base 20271004 of ¹³ was consumed by the T9 investigation (its full-scale results were
inspected and used to choose the Gamma straggling default and the stopping-power-ramp
apportioning), so it is no longer held out. The qualification seed base is now 20281004 and is
used only for the final qualification runs of this task after the engine is frozen; 20261004
(rehearsal) and 20271004 (T9 investigation) are preserved as non-qualification evidence and any
further diagnostic run uses a base from the 2027xxxx family. The runner enforces conformance only
for base 20281004. Tolerances and history counts are unchanged.

¹⁵ Profile χ² over supported bins, and the qualification seed base re-frozen again (2026-10-04, task V3-003B):
the T12 qualification run at 0d7b36e (base 20281004) passed every comparison except the lateral
profile χ² of the 4·10³-history python sample against the accelerated samples at 0.9·R (and
marginally 0.5·R). A control comparison of the python sample against the float64 Warp-CPU sample —
algorithmically identical code paths, bit-identical in T1 — failed in exactly the same way
(permutation p 0.0005 at 0.9·R, 0.001 at 0.5·R; IDD, 0.25·R, R80, σ and total deposit all pass),
while the Warp-CPU float64 vs CUDA float32 pair passed everything. The χ² is therefore
miscalibrated for sparse profile bins of the small python sample (discrete per-batch means break
the exchangeability of the studentized residuals), not evidence of a physics difference. Test
correction: the profile χ² includes only the supported bins of ¹²(b) (both samples with at least
max(2, ⌈B/2⌉) non-zero batches; unsupported bins are listed in the result), and for pairs with
unequal batch structures (the 4·10³-history python sample in 40 batches against accelerated samples
in 100 batches) its p-value is calibrated by a within-sample studentized bootstrap-t (batches
resampled with replacement within each sample, ≥ 2000 replicates, recorded seed, null formed by
recentring the bootstrapped difference at the observed one) — the studentized permutation of ⁹
remained anti-conservative there (≈ 15 % of sparse-null trials below p = 0.05, because the small
sample's batch means are skewed) while it is calibrated for equal batch structures (T9, T10 and
the accelerated-vs-accelerated pairs: χ² 150–160 for same-configuration controls), where it is
kept. A synthetic sparse-null calibration test (≤ 8 % below 0.05, ≤ 2 % below 0.01 over 200
trials) guards the method. The python-vs-float64-CPU pair joins the T12 comparisons as the
same-algorithm control. Because the 20281004 samples were observed, the qualification seed base is re-frozen to
20291004; 20261004, 20271004 and 20281004 are preserved as non-qualification evidence. Tolerances
and history counts are unchanged.
