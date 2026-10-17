# IR-20261017-223039-9DA564 — batched V3-005D operator intervention (category experiment-boundary)

Status: open (raised 2026-10-17T22:30:39Z at commit 804729ee; notification sent). Request record: `request.json` in this directory. The full decision memo follows.


- **Category:** experiment-boundary (ambiguous experimental boundary, EXPERIMENT.md §5).
- **Task id:** V3-005B (scope name V3-005D). PR #82 stays unmerged (disposition A of IR-20261013-164616-DEFAF4).
- **Route of record:** `validation/plans/v3-005d-integration-route.md`, revision 6.
- **Commit of record:** `804729ee7609c54acbcdcefeab82dce321466140` on the V3-005B task branch (PR #82, unmerged); every path below resolves at that commit. The decision 0041 and Amendment 19 texts use the record's ratio direction (ionmc table / Geant4).

## Summary (read first)

- The elastic over-correction of V5/V11 at dd976c5a is now **partially attributed** under the frozen rule (§5 of `validation/plans/v3-005d-elastic-attribution.md`; Amendment 19 addendum (d); decision 0041, section 2026-10-17).
- The p+A angular shape carries almost all of it.
- Two minor hypotheses could not be excluded: the p-p rate (H3) and the recoil share (H4).
- The plan authorises a model change only for "established", so **no model change has been made**.
- **What is asked:**
  - item 0: whether to change the model now;
  - items 1–5: how the failed merge gates are carried;
  - items 6–8: statements only, no decision needed.
- **Nothing here relaxes a tolerance, re-labels a verdict or reuses a seed.**
- **Recorded failures that stay as they are:**
  - V5, V11, V10-A, V7-R, V7-scan-el (dd976c5a);
  - V6 / V6-sens, failed — not sized;
  - lv5b V5 and V7;
  - P6 for H-1 below 12.53 MeV;
  - the MF6 copy rule.
- **Default without an answer:** no model change, no re-qualification, PR #82 stays open (Amendment 18 (c)).

## Facts carried (all exploratory, report-only; records under `validation/results/elastic_attribution/`)

| Fact | Value | Record |
|---|---|---|
| §5 clauses 1–4 | hold at both energies | `de6/de6-report.md`, `de6/de6-summary.json` |
| Clause 1: base reproduces archived V11 F | +0.0013 / +0.0015 (2 SE 0.0042 / 0.0075) | same |
| Clause 2: frozen 2 × 2 closure residual r | **+0.0014 (SE_r 0.0027) at 150 MeV, +0.0051 (0.0037) at 200 MeV**; threshold 0.010 | same |
| Clause 3: H2 (p+A shape) share of Δ(both) | 100 % / 99 % | same |
| swap-pA share of the observed archived V11 peak gap | 96 % / 88 % | same; review round D-E6 in `validation/plans/v3-005d-elastic-attribution/review.md` |
| Clause 4: D-E5 signature | holds (exceedance ratio ≥ 1.15 at 5 of 5 energies; tail moment ratio 2.4–25) | `de4/de4-ratios-report.md` |
| **Clause 5 fails: H3 (p-p) not excluded** | p-p rate ratio **ionmc table / Geant4 = 1.165 / 1.159** at 150 / 200 MeV (outside 1 ± 0.05; direction as defined in `de4/de4-ratios-report.md`, "ratio = ionmc table / Geant4"). (decision 0041 and Amendment 19 addendum (d) use this direction; the review record carries an orchestrator note); shape ratios ≈ 1.00; table σ_pp = 1.10 × BGG at 150 / 200 MeV | `de4/de4-ratios.json`; `de2/de2-report.md` |
| **Clause 5 fails: H4 (recoil) not excluded at 150 MeV** | recoil part 25 % of Δ_g against the 20 % bound (14 % at 200 MeV). It arises only inside the shape-swap arms; that is a contrary reading, and the verdict is unchanged | `de6/de6-report.md`; review round D-E6 |
| H1 (σ level), H5 (MCS overlap), H6 (case diff) | excluded (p+A rates 1.005–1.031; EM-only within ±2 %) | `de4/`, `de0h5/`, `de6/` |
| Post-hoc p-p σ arm (addendum (c)3, declared after D-E4) | 6 % / 13 % of Δ(all), r_all −0.0007 / −0.0012; changes no label | `de6/de6-report.md` |
| Mechanism | The peak lever is the **core width**, not the tail. `pA-gauss` leaves F unchanged; `pA-s0.9` closes most of the gap; CHIPS has the narrower core. The original sharp-tail mechanism (H2a) is **refuted for the peak** | `de3/de3-summary.md`; review rounds D-E3, D-E6 |
| Contrary evidence | r = +0.0051 at 200 MeV (≈ 12 % of the gap left); the gain closure overshoots at 150 MeV and undershoots at 200 MeV; the peak-versus-gain ordering flips between the energies | review round D-E6, item 5 |
| Plan conflict | §5 authorises a model change only for "established"; §4 says "M2 adopted if §5 attributes H2". Not read as an authorisation; recorded | decision 0041, 2026-10-17 |

**What is NOT known.**
- The swap test shows that Geant4 CHIPS's shape closes the gap. It does **not** show that the pre-committed replacement M2 (optical-limit Glauber) reproduces the CHIPS core width. No M2 table exists yet.
- Whether the p-p rate excess (H3) matters at V5's tolerance after a shape change. Estimate: the 6 % / 13 % of the σ arm.
- Whether H4 is an independent bookkeeping defect or the recoil part of H2.

## Item 0 (new): model-change authorisation

| Option | Content | Cost | Consequences |
|---|---|---|---|
| **0A (recommended)** | Adopt **M2** on the partial attribution, with H3 and H4 carried as **open findings** into the pre-registered re-qualification. The inputs are pre-committed in plan §4.1: SAID NN σ_NN/α_NN (fallback SP07), Horiuchi 2007 β_NN, de Vries 1987 densities, r_p 0.8414 fm; σ_el stays Barashenkov. A new decision 00xx records the lineage and the §5 contradiction verbatim. The new-model V10-A is report-only (data seen). The M3 triggers of plan §4.2 apply unchanged. **Added safeguard, pre-registered before the M2 table is built:** a report-only M2 diagnostic (D-E3 harness F at 150 / 200 MeV, and the D-E5 comparison against the recorded D-E4 data). It may stop host spending and escalate, and it never selects another model | **Acquisitions:** three `acquire_reference_data` calls (SAID live service may be unreachable; then SP07; if neither, M2 is "construction impossible" → M3). **Calendar:** decision, builder, tests, Codex. **Sandbox:** the diagnostic at 15–43 s per 1e5-history run (`wall_s` values in `de3/de3-summary.json`), plus ≈ 2.2 h for the V7-R2 joint-power re-measurement. **Host:** the re-qualification only, ≈ 37.2 h (item 8) | The fastest route to a V5/V11 re-test. Risk: M2 may not reproduce the CHIPS core. V5/V11 would then fail again at the one granted base (route §6 S2). The diagnostic shows this before the 37 h are spent |
| 0B | First pre-register and run a follow-up that resolves H3 and H4, then decide. H3 needs a p-p σ question answered: the table is 1.10 × BGG above 100 MeV and 1.16 × Geant4 in rate. H4 needs a recoil-bookkeeping test | **Sandbox:** hours of D-E3-type runs. **Host:** 0, unless new TOPAS rows are needed (each ≤ 1500 s). Then 0A's costs follow anyway | Delays the re-test by one review cycle. **A declared σ_pp correction is itself a model change.** Moving σ_pp toward BGG/Geant4 would give the p-p channel TOPAS lineage, and that channel is currently independent of TOPAS (Amendment 14 (c)3). The H3 effect measured so far is 6 % / 13 % of Δ(all) |
| 0C | No model change. V5, V11 and V10-A stay failed. Integrate (if at all) via Option I of item 1, with the elastic depth-dose claim withdrawn | **Host:** 0 with item 3 = 3B (dd976c5a archives plus a truthful head-check record at the PR head; `src/` is unchanged since dd976c5a). ≈ 37 h with 3A (V7-R2 must run at an evidence SHA) | The elastic channel ships without a validated depth-dose claim, or stays off by default |

**Recommendation: 0A with the diagnostic safeguard.**
1. The frozen closure holds, and the named component carries 100 % / 99 %.
2. The clause-5 failures are small: the H3 σ arm is 6 % / 13 % of Δ(all). H4's excess arises only in the shape-swap arms.
3. 0B's natural fix for H3 would move an independent channel onto TOPAS lineage.
4. The diagnostic exposes the main unknown, M2's core width, before any host time is spent.

**Fallback:** 0B if the operator requires clause 5 to be resolved first. 0C if neither is acceptable.

## Items 1–8 (route §4.1 (iv); unchanged in substance)

1. **Merge gate (a)6.** It cannot be met as written: V6 is failed — not sized; V7-R's frozen design cannot pass; V10-A (c) fails on the shared σ_el.
   - **Option I (recommended):** amend (a)6. V6 and V10-A are carried as recorded failures, V7-R2 replaces V7-R, and the claims are narrowed: no lateral-halo or measured-dose claim; the p+A shape is gated by V11.
   - **Option II:** as I, plus V6-F executed once as a binding gate (+≈ 6.2 h CPU plus CUDA).
   - **Option III:** literal (a)6, i.e. never merge.
   - Options I and II **reduce the gate set** and override the Amendment 18 (c) default. The archive will read `pass: false`, unlike the V3-005A precedent.
2. **V10-A after a shape change.**
   - **2B (recommended):** the dd976c5a failure is carried; the new-model V10-A is report-only (data seen; (e)9 precedent); V11 is added to (a)6.
   - 2A: count the data-seen result toward (a)6.
   - 2C: literal, cannot pass.
   - Under 0C this item is moot.
3. **V7-R slot.**
   - **3A (recommended):** the new row V7-R2 takes V7-R's (a)6 slot; V7-R stays failed.
   - 3B: V7-R carried; `escaped_neutral` coverage claim withheld.
   - 3C: literal.
4. **V7-R2 calibration variant.** **As written, recommended** (`validation/plans/v3-005d-v7r-calibration-repair.md`, revision 4).
   - α_tost 0.005, 48 × 2.25e6, R = 3200.
   - The size rule uses the measured two-gate joint pass rate at κ = 1, which needs ≥ 378/400. Today: 383/400 and 384/400 → 5400/0.005 selected. It is re-measured on the rebuilt pool.
   - Cost ≈ 33.9 h, or ≈ 44 h at 7200.
   - Alternatives: (B) exact paired certification with R = 800, ≈ 16.5 h; (C) deterministic certification, ≈ 12.5 h. Both change how a clause is evaluated.
5. **V5 fresh base (decision).**
   - Amendment 17 (b)4: once the stop rule fires, "no fresh base is taken for V5, no tolerance changes, and the orchestrator raises an intervention". It fired at dd976c5a.
   - **Recommended:** grant **one** fresh V5 base (lv 20531004) after the item-0 model decision, with tolerances unchanged.
   - Default if refused: V5 stays failed (route §6 S1). Moot under 0C.
6. **Bases (statement).** lv 20531004 / hr 20541004; V6-F 20551004 (Option II). The pre-declared 20511004 / 20521004 are retired unused, declared before observation, so no seed is repeated. Seed-by-seed record: route §4.4.
7. **Conformity (statement).** Single-process directive. Archives `conformant: false` (imported partials; deferred partition step 86). The operator may lift the directive; it is not requested.
8. **Host budget (measured facts; route §4.2).** D-E4 is done (18 OpenTOPAS runs, `de4/de4-report.md`); its host time is spent and not re-budgeted.

   | Plan | Re-qualification host hours |
   |---|---|
   | 0A with V7-R2 at 5400 (as selected today) | **≈ 37.2** |
   | V7-R2 at 7200 | ≈ 47.3 |
   | with calibration (B) | ≈ 19.8 |
   | with calibration (C) | ≈ 15.8 |
   | Option II | + ≈ 6.2 CPU + CUDA |

   - Sandbox, not host: ≈ 2.2 h joint-power re-measurement, plus D-E3-type diagnostics at under a minute per 1e5-history run.
   - These supersede the ≈ 9.5 h quoted for V7-R in IR-DEFAF4.
   - Under 0C with 3B: 0 h.

## Requested input

Choose for item 0 (0A / 0B / 0C) and for items 1–5, or state another disposition. Items 6–8 need no answer. The response is recorded verbatim with `resolve_human_intervention` and changes no row verdict.

## Basis of the request and the orchestrator's own disposition (added 2026-10-18 after Codex REVIEW-4996f6dd)

Codex REVIEW-4996f6dd (on commit f3c3d9b9) objected that this request "misclassifies a contradiction between two internally authored validation-plan clauses as an ambiguous experimental boundary" and "asks the operator to select the physics model, statistical calibration design, acceptance-gate changes, and fresh evaluation base", which EXPERIMENT.md assigns to the orchestrator. The objection is preserved here verbatim in substance. The request is **not withdrawn**, for these recorded reasons:

1. **The operator required it.** The resolution of IR-20261013-164616-DEFAF4 (2026-10-13, verbatim): "Before expensive qualification, document and independently review the route to integration, including remaining requirement coverage, V6/V6-F treatment, and execution conformity under the existing stability restrictions. V6-F passing does not authorize replacing V6 or merging. Bring back any concrete proposal that changes required validation scope; routine scientific investigation and implementation should continue autonomously." Items 1–5 change the required validation scope (the (a)6 merge-gate set, the V10-A and V7-R treatment, the V7-R2 design that replaces a frozen row, a fresh V5 base). Item 0 changes the re-qualification scope (≈ 37 host hours) and rests on a frozen pre-registration whose two clauses conflict.
2. **The frozen pre-registration requires it.** Amendment 17 (b)4: once the V5 stop rule fires, "no fresh base is taken for V5, no tolerance changes, and the orchestrator raises an intervention". It fired at dd976c5a. Taking a fresh base without the intervention would break the pre-registration.
3. **Category.** Of the five permitted categories, the request falls under *experiment-boundary* (the frozen plan's §4/§5 conflict on what authorises a model change, and the operator's standing "bring back" boundary) and equally under *authorization* (host budget ≈ 37 h; a fresh base after a stop rule; reducing a frozen merge-gate set). The tool records one category; the request carries both grounds.
4. **No scientific judgment is delegated.** The attribution verdict ("partially attributed"), the mechanism, the V7-R2 design and the recommendation for every item are the orchestrator's own and are independently reviewed in the cited records. An operator response is recorded as a scope and authorization decision, never as scientific validation; it changes no row verdict, tolerance or seed. The orchestrator's disposition, stated in the memo, is: 0A with the pre-registered diagnostic safeguard; Option I; 2B; 3A; V7-R2 as written; one fresh V5 base at lv 20531004. **Default without an answer** stays as stated: no model change, no re-qualification, PR #82 open.

## Correction of the request (2026-10-18, after Codex REVIEW-787e49bb; request.json is the tool's record and is not edited)

Codex REVIEW-787e49bb maintains that the model selection, the calibration design and the fresh evaluation base are orchestrator decisions and that the ≈ 37 h host run is not an authorization-category action. Applied correction and recorded disagreement:

- **Item 4 is withdrawn as a question.** The V7-R2 calibration variant is a design choice: the orchestrator decides **as written** (revision 4 of the repair plan; 5400/0.005 by the measured rule) and only informs the operator. Variants (B)/(C) are not proposed unless the operator raises a budget constraint.
- **Items 0, 1, 2, 3 and 5 remain operator items**, on these grounds, each a change of required validation scope or of a frozen pre-registration — exactly what the operator's resolution of IR-20261013-164616-DEFAF4 asked to be brought back ("Bring back any concrete proposal that changes required validation scope"): item 0 adopts a model although the frozen §5 clause authorises a model change only for "established" (a change of the pre-registered rule's consequence; the orchestrator may not relax a frozen rule); items 1–3 reduce the frozen (a)6 merge-gate set (Options I/II override the Amendment 18 (c) default) and re-assign gate slots; item 5 takes a fresh V5 base after the Amendment 17 (b)4 stop rule fired, which that amendment makes conditional on an intervention.
- **Unresolved disagreement, recorded:** Codex holds that these are autonomous decisions; the orchestrator holds that the operator directive and the frozen amendments reserve them. The operator's response resolves the disagreement either way. Until it arrives the conservative default applies: no model change, no re-qualification, no host run, PR #82 open. Routine investigation and implementation that change no scope continue autonomously.
