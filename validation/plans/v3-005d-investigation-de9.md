# V3-005D investigation D-E9 — the 200 MeV gain residual: archived-baseline sensitivity and the p-p contribution (prospective, report-only; revision 2, independently reviewed in rounds D-E9-1 (REVISE) and D-E9-2 (ACCEPT) of `validation/plans/v3-005d-elastic-attribution/review.md`; operator authorisation: response to IR-20261018-091819-903058, 2026-10-19)

## Review history

| Revision | Review round | Changes |
|---|---|---|
| 1 | D-E9-1: REVISE, nothing blocking | (1) A1 calibration ranges set to [0.556, 1.418] (K = 8) and [0.342, 1.614] (K = 4). (2) Lead A sensitivity rows labelled "sensitivity; frozen d = 0.00462 stands", never compared with the band; K fixed. (3) B0: independence, the "dip persists" rule and the relative 1e-3 tolerance defined before computing. (4) B2 restricted to the existing scoring channels, no `src/` change, dropped if they cannot separate. (5) B3 no-inversion rule; the monotone variant is never a candidate. (6) Physics-reason criteria: the "sufficient for the residual" criterion removed; criterion 2 recast as a level discrepancy beyond combined uncertainties. (7) Absolute cross sections replaced by ratios; the reviewer's A20-3 σ_anchor ratios and the PP-LEVEL ratios disclosed; the list of what needs an amendment added |
| 2 | D-E9-2: ACCEPT, with three inserted statements | (a) The PP-LEVEL verdict is the researcher's, not adopted as a plan conclusion; B0 decides. (c) Table/SGT quoted only at its minimum (0.821) and at the step (1.021). (d) Both table/anchor value sets recorded (1.2–2.2 % apart; cause unresolved; not B0-relevant), and the deviation's span (125–150 MeV) stated as covered by B0's dip test |

**Status: exploratory and report-only (Amendment 19 class), never evidence.**
- Written 2026-10-19 under the operator's response to IR-20261018-091819-903058: "Continue bounded autonomous investigation of the remaining discrepancy. Assess the documented proton–proton contribution and the sensitivity of the screening result to uncertainty in its archived baseline. These are investigative leads, not instructions to tune to TOPAS, substitute a favorable baseline, or repeat seeds until a test passes."
- **Nothing here changes a recorded verdict.**
  - The D-E8 STOP (failing (i), (i-b) and (ii) at 200 MeV) stands.
  - V5, V11 and V10-A stay failed.
  - D-E8 is never re-run or re-evaluated.
  - The ≈ 43.5 h closure stays deferred.

## 0. Disclosure

**Both leads were identified after the D-E8 STOP** (review round D-E8, item 3).

**Numbers already seen before this plan:**
- the archived V11 ionmc gain 0.02672 at 200 MeV against 0.0249 at the D-E3 CRN seed;
- D-E6 σ_pp arm Δ_g −0.00287;
- the table/Geant4 p-p rate 1.165 / 1.159;
- the non-monotonic Hale-NI table level of the Amendment 21 (b′)4 rev 4 lead, here as a ratio: level(150 MeV)/level(130 MeV) = 1.18. The absolute values are in the committed Amendment 20 / 21 text and are not repeated in de9 records;
- the reviewer's σ_anchor computation in round A20-3, as ratios table/σ_anchor at 100 / 125 / 150 / 200 MeV: **1.001 / 0.936 / 1.131 / 1.142** (from the values committed in Amendment 21 (b′)4 rev 4).
- **the physics researcher's report-only quantification of the pinned table's p-p level** (PP-LEVEL record, ratios only; committed with this plan as a projected record), **obtained after the STOP.** No table of record changes under this plan:

  | T_lab [MeV] | 100 | 125 | 130 | 150 | 175 | 200 | 250 |
  |---|---|---|---|---|---|---|---|
  | table/SAID, like-for-like (same \|μ_CM\| ≤ 0.96 cut, same Mott subtraction) | 0.979 | 0.922 | – | 1.119 | 1.128 | 1.131 | 1.129 |
  | table/BGG | 0.989 | 0.920 | – | 1.101 | 1.101 | 1.101 | 1.101 |

  - Above 150 MeV, table/BGG is 1.1008 at every node **by construction**: the 150 MeV node value is carried upward by the S(E)/S(150) scaling (`elastic_build.py` docstring, lines 33–44).
  - Table/SGT is quoted only at its minimum, **0.821 at 130 MeV**, and at the step, **1.021 at 150 MeV**: about 20 % between 130 and 150 MeV. Neither SAID reference has a physical counterpart: both fall smoothly, and there is no threshold below about 280 MeV. The cut and interference conventions are equally present at 100–125 MeV, where the table sits 2–8 % **below** SAID.
  - The D-E4 table/Geant4 p-p ratios 1.165 / 1.159 are this same step measured against Geant4, whose p-p rate behaves as BGG σ_pp. The expected table/(BGG × accepted fraction) is 1.153 / 1.148, against the measured 1.165(15) / 1.159(15).
  - **The researcher's verdict:** a defect of the table construction entering at the 150 MeV node. It is undetermined whether it sits in the LA150 H-1 LAW=5 data at the top nodes or in the reconstruction.
  - This verdict is the researcher's. It is not adopted as a conclusion of this plan; B0 decides.
  - **Two sets of table/anchor values.** The researcher's values sit 1.2–2.2 % below the reviewer's (A20-3): 0.979 / 0.922 / 1.119 / 1.131 against 1.001 / 0.936 / 1.131 / 1.142 at 100 / 125 / 150 / 200 MeV. The likely cause is the Mott η / k definitions or the angle interpolation. Both sets are recorded, with the cause unresolved; the difference is not B0-relevant.
  - **The deviation spans 125–150 MeV**: table/SAID ≈ 0.92 at 125–130 MeV, 0.95 at 140 MeV and 1.12 at 150 MeV. B0's dip test covers this range.
  - **Open observation, report-only and not a lead:** the D-E4 30 MeV p-p cell differs by about 6σ. It has not been investigated.

**Computed now by the planner from the generated D-E3 documents of record** (sha256 bound in `de6/de6-report.md`; plateau 20–60 mm, ratio of batch means; reported here and not yet in a committed record):

| Quantity | 150 MeV | 200 MeV |
|---|---|---|
| Base gain g_all | 0.02207 (jackknife SE 0.00079) | 0.02491 (SE 0.00099) |
| pp-only | 0.01368 | **0.01982** |
| pA-only | 0.00809 | **0.00475** |

**At 200 MeV the p-p channel carries about 80 % of ionmc's elastic plateau gain.**

## 1. Lead A — sensitivity to the archived baseline

**Facts.**
- The frozen D2 is d = g_arch + (g_M2′ − g_base) − g_TOPAS,arch.
- g_arch = 0.02672 is the archived ionmc V11 gain at 200 MeV. Its source is lv5c **step 16** (V11 compare, ionmc seed 20495005), bound through `de0h5/h5.json`, not step 21. Its 20-batch SE is 0.00077 (`plateau_gain/ionmc/se`).
- The D-E3 base realisation at the same configuration (same pinned table; 1e5 histories, 20 batches, V5 geometry) gives 0.02491 ± 0.00099.
- The difference is 0.00181, z = 1.44.
- Δ = g_M2′ − g_base = −0.00369 (common random numbers; SE 0.00039).

**Report-only sensitivity, computed now** (not a re-evaluation; not a substitution). **No row except the first is compared with the 0.00365 band. Each row is labelled "sensitivity; frozen d = 0.00462 stands".**

| Reading | d (report-only) |
|---|---|
| The prediction as frozen, with prediction SE √(0.00077² + 0.00039²) = 0.00087 | 0.00462 (the frozen value; the STOP stands) |
| Archived baseline ± 2 SE | 0.00307 – 0.00616. Sensitivity; frozen d = 0.00462 stands |
| D-E3 baseline | 0.00281. Sensitivity; frozen d = 0.00462 stands |
| Inverse-variance pooled baseline, 0.02604 ± 0.00061 | 0.00393. Sensitivity; frozen d = 0.00462 stands |

**Conclusion of the sensitivity:** the frozen prediction of D2 carries the archived baseline's own SE (0.00077), which is comparable to its margin. The screen's band accounts for the future run's SE (Amendment 20 (c)3) but not for the SE of the archived term inside D2. **This is a planning finding about the screen's design. It does not change the STOP.**

**Diagnostic A1, pre-registered: is 0.02672 a fluctuation?**
- **Runs:** K = 8 independent common-random-number pairs (pinned table: EM+elastic and EM-only) at 200 MeV, plus K = 4 at 150 MeV as control.
  - Harness `de3_probes.py`: elastic-only configuration, 1e5 histories, 20 batches, lateral off.
  - New seeds: 20992000 + k at 200 MeV (k = 0..7) and 20992010 + k at 150 MeV (k = 0..3).
- **Reading rule, fixed now:**
  1. μ̂ = mean gain over the K pairs; s = sample SD; SE_μ = s/√K.
  2. **SE calibration:** ρ = s / mean(per-run jackknife SE). At 200 MeV (K = 8), ρ outside **[0.556, 1.418]** (90 % range of √(χ²₇/7)) is recorded as "per-run SE mis-calibrated". At 150 MeV (K = 4 control), the range is **[0.342, 1.614]** (√(χ²₃/3)). This matters because the frozen bands use per-run SEs.
  3. z_arch = (0.02672 − μ̂)/√(0.00077² + SE_μ²), and z_DE3 the same way for 0.02491. |z| ≤ 2 → "consistent with the realisation distribution"; |z| > 2 → "outlier at 2σ". Recorded for both values.
  4. Report-only d(μ̂) = μ̂ + Δ − g_TOPAS,arch with its SE, labelled "sensitivity; frozen d = 0.00462 stands". **It is never compared with the 0.00365 band.** K is fixed: 8 + 4 pairs, never extended.
- **Supports:** whether the archived baseline is representative; the prediction SE of D2; the design lesson that a single archived realisation needs its SE in the band. **Does not support:** any change of the D-E8 verdict, substituting μ̂ for g_arch, a prediction of a fresh lv5d V11 run, or re-running D-E8.

## 2. Lead B — the p-p contribution (open finding H3)

**Facts.**
- The table's Hale-NI partial σ_pp is **non-monotonic** (PP-LEVEL: table/SGT minimum 0.821 at 130 MeV and 1.021 at the 150 MeV step; table/SAID like-for-like 1.119–1.131 above 150 MeV; table/BGG 1.1008 by construction above 150 MeV). Level(150)/level(130) = 1.18. The S(E) rule above 150 MeV inherits the 150 MeV level. Table/σ_anchor (SAID, half sphere) is 1.131 / 1.142 at 150 / 200 MeV.
- The p-p rate ratio table/Geant4 is 1.165 / 1.159 at 150 / 200 MeV.
- The D-E6 σ_pp arm moved the 200 MeV gain by −0.00287.
- The D7 recoil share at 200 MeV is 0.13, below the 20 % bound (`de8`), so the recoil handling is not the leading term.

**Diagnostics, in order (cheapest first; no transport table of record is changed):**
- **B0, deterministic, minutes: P6-D-type check of the Hale NI integration, 100–200 MeV.**
  - **Independent** means that the check does **not** import `elastic_build`, `endf6` or the LAW=5 reconstruction code. It reads the raw LA150 H-1 file directly (the pinned acquisition) and implements ENDF-102 eqs 6.9 / 6.10 / 6.14 anew.
  - Grid: 1 MeV in energy over 100–200 MeV; the half sphere 0 ≤ μ_CM ≤ 0.96, with a 4× finer μ grid than the builder.
  - **Definitions, fixed before computing:**
    - "match" = relative difference ≤ 1e-3 at every table node and node midpoint;
    - "dip persists" = an interior minimum between 100 and 150 MeV that lies more than 1 % below both endpoint values.
  - **Reading, data against reconstruction** (fixed before computing). The independent reconstruction of the cut integral at the 130, 140 and 150 MeV nodes is the discriminating step:
    - **the independent reconstruction reproduces the step from the raw data** (it matches the table within 1e-3 and the dip persists) → an **evaluation feature at the LA150 top nodes**. The question then becomes whether the S(E) carry-up of a top-node value is a construction choice. That choice is recorded, and any change to it needs its own amendment (§3);
    - **it does not** (a mismatch beyond 1e-3, or no dip) → a **builder defect**, which takes the defect route and an amendment (§3).
- **B1, zero cost: decomposition already recorded.** The pp-only / pA-only / all gains above enter the de9 record, with their additivity (0.01982 + 0.00475 against 0.02491).
- **B2, existing scoring channels only (orchestrator decision; no `src/` change).** B2 uses only the existing species and genealogy scoring channels (decision 0040; the harness tallies) to separate the plateau dose of p-p secondary protons (generation 1) from primaries, in elastic-only runs with the pinned table.
  - Seeds 20992100 (150 MeV) and 20992101 (200 MeV).
  - **If the existing channels cannot separate them, B2 is dropped, and that is recorded.**
  - **Reading:** the share of the elastic plateau gain carried by p-p secondary dose.
- **B3, report-only sensitivity family (exploratory tables; `make_variant_tables --scratch-cache`; loader-refused outside this step).**
  - Variants:
    - (a) σ_pp × s at E ≥ 150 MeV, with s ∈ {0.90, 0.95};
    - (b) a **monotone-level** variant: the 130–150 MeV upturn replaced by the table's own 125 MeV level carried by the table's own S(E) ratio. No external data and no Geant4 input.
  - The Hale shape is unchanged.
  - Common random numbers: seed 20990101 at 200 MeV and 20990100 at 150 MeV (the declared D-E3 CRN seeds).
  - **Reading:** the slope dg/d ln s at 150 and 200 MeV, and the Δg of each variant. These are reported as gains only.
  - **No-inversion rule.** The record **never computes or reports the s, or the level, that would close the residual**, and never relates a variant's Δg to the residual.
  - **The monotone variant was designed after the dip was seen.** It stays a sensitivity and is **never a candidate model**. No variant is adopted.
- **What would constitute a physics reason for a prospective p-p amendment** (any one, independently reviewed):
  - B0 finds a builder or interpolation defect: the defect route;
  - B0 confirms the dip is in the evaluation, **and** the table's level differs from an independent-lineage level (SAID-type; not Geant4 / TOPAS) **beyond their combined uncertainties**, established by deterministic data comparison. Monotonicity alone is a shape statement and does not qualify.

  Even then, a σ_pp change is a new model change. It needs its own decision, pre-registration and re-qualification. The previously applied SAID re-anchoring rule failed its own continuity check (Amendment 21 (b′)4 rev 4) and stays not applied.

## 3. What needs an amendment rather than an exploratory record

- Any `src/` change. B2 is restricted so that it needs none.
- Any B0 defect: the defect route, with an intended-change record and a new reviewed amendment.
- Any screen-design change learned from Lead A, such as baseline SE in a band, for a future screen.
- Any σ_pp change.

## 3b. What neither lead can justify

- Tuning σ_pp, or any parameter, to TOPAS or Geant4.
- Substituting a baseline in D2.
- Re-running or re-evaluating D-E8.
- Repeating seeds, or adding realisations until a criterion passes.
- Changing a band.
- Any passing record.

D-E9 outcomes are recorded **whichever way they fall**, and are distinguished in every record from the historical qualification archives (lv5c / hr5c dd976c5a) and from the D-E8 screen.

## 4. Seeds, budget, records

**Seeds.** New registered block **20992000–20992199** in the 2099xxxx elastic family:

| Seeds | Use |
|---|---|
| 20992000–20992007 | A1, 200 MeV |
| 20992010–20992013 | A1, 150 MeV |
| 20992100–20992101 | B2 |
| all other seeds of 20992000–20992199 | reserved |

B3 reuses the D-E3 CRN seeds 20990100 / 20990101, declared. The seed-family test is extended to register and assert the block before any run.

**Budget: sandbox only, single process, one heavy job, no host job.**

| Item | Size | Wall time |
|---|---|---|
| A1 | 24 runs | ≤ 43 s each (D-E3 `wall_s` 15–43 s), ≈ 17 min |
| B2 | 2 runs | – |
| B3 | 6 runs (3 variants × 2 energies) | – |
| B2 + B3 together | – | ≈ 6 min |
| B0 | – | minutes |
| **Total** | – | **under 40 min** |

**Records.** `validation/results/elastic_attribution/de9/` under the `common.py` projection:
- gains and ratios, z values, flags and SE ratios;
- **no absolute cross sections**: the B0 result is projected as ratios to the 125 MeV node and as monotonicity flags;
- provenance: HEAD, clean flag, script sha256, seeds, an attempt log under the D-E8 crash and invalid-run rules.

The header states "exploratory; leads identified after the D-E8 STOP; not a re-evaluation".

**Order:**
1. this plan, independently reviewed;
2. the seed block registered;
3. harness options (generation tally; the B3 variants) committed and Codex-reviewed;
4. B0;
5. B1;
6. A1;
7. B2 and B3;
8. one synthesis record;
9. any next step only through a new prospective amendment.

## 5. Independent-review checklist

1. No step changes a recorded verdict, D-E8, a band, or a seed or table of record.
2. Lead A is report-only, with its reading rule and K fixed before A1. No baseline substitution; both z values are reported.
3. g_arch provenance (step 16 via h5) and its SE (0.00077) are correct.
4. B3 uses only the table's own data and rules; s is never selected.
5. B0's thresholds and its feature/defect definition are fixed before computing.
6. The "physics reason" criteria are evidence, not V5/V11 outcomes.
7. The seed block is disjoint from all used, reserved and rehearsal seeds.
8. The Git projection contains no absolute cross sections.
9. The disclosure of seen numbers is complete.
10. Single process, one heavy job, no host use.

*Dated note (2026-10-19, after Codex REVIEW-a701b2c8 and before any B0 number exists): domain of B0.*

The raw LA150 H-1 LAW=5 data end at 150 MeV. Above 150 MeV, the table's level is the builder's S(E)/S(150) carry-up of the 150 MeV node. B0 is therefore evaluated in two blocks.

1. **Raw-data domain, 100–150 MeV.**
   - The ≤ 1e-3 relative match is checked at table nodes and node midpoints.
   - The dip test ("an interior minimum between 100 and 150 MeV more than 1 % below both endpoints") uses a 1 MeV grid.
2. **Carry-up check, 150–200 MeV (separate block).**
   - The level above 150 MeV is reconstructed independently as level(150, reconstructed) × S(E)/S(150), with S = `bgg.pp_bgg_mb`. The block imports only `ionmc.nuclear.bgg`, nothing else from the package.
   - It is compared, as ratios, with the pinned nodes at 175, 200 and 250 MeV, and with the table's interpolated level on a 1 MeV grid over 150–200 MeV.
   - Ratios within 1e-3 mean "the carry-up rule reproduces". Any other result is reported as found.

The carry-up is a construction choice that is already disclosed (PP-LEVEL: table/BGG = 1.1008 by construction above 150 MeV). It is **not part of the "dip persists" rule**.

All other reading rules are unchanged.

*Dated note (2026-10-19, after Codex REVIEW-6b2120c4 and before any B0 number exists): B0 angular quadrature.*

The pre-registered phrase "a μ grid four times finer than the builder" is corrected to: **Gauss–Legendre quadrature (1024 nodes) on |μ| ≤ 0.96.** This is a different rule from the builder's 4001-point uniform Simpson grid.

The B0 record carries a **convergence block** at the 100, 125 and 150 MeV nodes, comparing four rules:
- GL 1024;
- GL 2048;
- Simpson 4001;
- Simpson 20001.

**Pre-declared acceptance:** the quadrature error is below 1e-6 relative, negligible against the 1e-3 match tolerance.

*Dated note (2026-10-19, after Codex REVIEW-1b6174d3): chronology of the quadrature note.*

Correction: before the previous note, the implementer had observed the 100 MeV convergence check and the 100 MeV level match, both at about 2e-12 (disclosed in `PRE_RUN_LOOKS`).

The 1e-6 acceptance is therefore data-informed for that node.

**Independent rationale.**
- 1024-node Gauss–Legendre quadrature integrates polynomials up to degree 2047 exactly.
- The LAW=5 LTP=1 integrand is a Legendre sum of degree NL ≤ 64, times Mott and interference factors analytic on |μ| ≤ 0.96, for which convergence is geometric.
- The expected error is near machine precision, far below 1e-3.

**Reading rules unaffected.** The 100 MeV node is an endpoint, not interior to the dip test. Its table value was already known.

*Dated note (2026-10-19; after Codex review and B0's result; before any B3 number).*

**(a) B3 monotone variant.** The table has no 125 MeV node (the cache rebuild failed before any B3 run). "The table's own 125 MeV level" therefore means the table's p-p level at 125 MeV, interpolated between bracketing nodes by the table's convention. It is carried upward by S(E)/S(125), with S = `bgg.pp_bgg_mb`. This makes a report-only sensitivity precise; B0's result does not bear on it.

**(b) B0 binding.** B0 ran once, at clean, reviewed HEAD 0dfca566, before the B3 defect was found. It stands and is not re-run. The harness fix changes the frozen `make_variant_tables.py`. The evaluator therefore binds B0 by ancestry of its HEAD to the reviewed SHA, the `de9_pp_level.py` hash, a clean tree and the attempt-log digest. A1, B2 and B3 stay bound by HEAD equality.

*Dated note (2026-10-19; before any synthesis number: no A1, B2 or B3 statistic computed or seen).*

The single production sequence ran at HEAD db7d84e2 (reviewed REVIEW-50f285f6). All 32 run documents were published once and are published-valid in the attempt log: A1 24 (8 + 4 pairs), B2 2, B3 6.

The synthesis refused before writing anything, on an evaluator defect: the monotone variant's upper bracket is exactly 130 MeV, against a "< 5 MeV" check.

The fix is confined to `de9_record.py`. The producers (`de3_probes.py`, `make_variant_tables.py`, `de9_pp_level.py`) are unchanged, and their hashes are bound in every document.

Once-only forbids re-producing the documents, and the evaluator runs at a later reviewed HEAD. The documents are therefore bound like B0: HEAD an ancestor of the reviewed SHA, a clean tree, producer hashes equal to the reviewed commit's, and the attempt-log digest. No document, seed, setting or reading rule changes.
