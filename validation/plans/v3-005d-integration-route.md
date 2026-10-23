# V3-005 route to integration (V3-005D planning document; planner output, not evidence)

Status: revision 7 (2026-10-18; one-line correction: the V6-F r_index-7 seed range ends at 20558019, k = 0..15, not 20558015 — found by the Amendment 20 review); revision 6 (2026-10-14), after two independent review rounds and Codex REVIEW-bb2a8151 (`validation/plans/v3-005d-integration-route/review.md`); to be committed as `validation/plans/v3-005d-integration-route.md`, the "route to integration" that the operator requires. Branch state: V3-005B worktree `/home/wahln/aiprojects/ion-mc-worktrees/v3-005b`, HEAD 749f768c (PR #82, unmerged by operator decision). Evidence SHA of the last qualification: dd976c5a (bases lv5c 20481004, hr5c 20491004).

This document changes no frozen row, tolerance, seed base or verdict. Every row result quoted here is preserved as archived. Proposals that change the required validation scope are collected in §4 (iv) for ONE operator intervention; nothing in them is applied before the operator answers.

## Revision 6 (2026-10-14)

Triggered by Codex **REVIEW-bb2a8151** on 9a7b142; synchronises the route with revision 4 of the V7-R calibration-repair plan. Nothing else changes: Options I/II/III, the intervention item list, the seed bases, and every other section are unchanged.

- **Size-selection rule.** The V7-R2 rule is now the **directly measured complete two-gate joint pass rate at κ = 1**: the repository `run_calibration`, with both gates on the same replicates, α_tost = 0.005, 48 × 225 replicates, R_meas = 400 per surrogate, seed roots 20989920–20989923 (planner block). Freeze 5400/0.005 (48 × 2.25e6) iff both surrogates reach **≥ 378/400** joint passes, i.e. the CP95 lower bound p_low satisfies P(Bin(3200, p_low) ≥ 2909) ≥ 0.99. Otherwise measure 7200/0.005 (64 × 2.25e6; roots + 10) and freeze it iff both reach ≥ 378/400. Otherwise no V7-R2 design is frozen and the finding goes to the operator (repair plan revision 4, §3.8 item 2).
  - It replaces the revision-3 independence product OC_paired × OC_single with a fitted box shrink, everywhere the route restated it: §2.2 Option I description, §2.4, §4.1 (i), §4.1 (vi) step 2.
  - Measured on the current pool (repair plan revision 4, `validation/plans/v3-005d-v7r-calibration-repair/data/joint_power_5400.*`): archive-fitted 383/400 (CP95 lower 0.9369), rehearsal block sums 384/400 (0.9399), so **5400/0.005 is selected today**.
  - The revision-3 powers 0.954 / 0.962 survive only as labelled model values.
- **Rebuilt pool.** The rule is re-applied on the rebuilt pool before the sizes are frozen. That re-measurement is sandbox work: single process, one heavy job, never overlapping the D-E3 probes or any host run. It takes ≈ 2.2 h for 5400/0.005 and ≈ 3 h more if 7200/0.005 has to be measured. It is added to the critical path (§4.3).
- **Budgets.** No host hours are added, so the host totals are unchanged: ≈ 41.3 / 51.5 / 23.9 / 19.9 h (§4.2). The V7-R2 block stays 33.9 h / 44 h, (B) 16.5 h, (C) 12.5 h.

## Revision 5 (2026-10-14)

Answers review round 2 (REVISE, minor, no blocking).

- **R2-1 (selection rule).** The revision-2 window rule is replaced in the revision-4 header, §2.4 and §4.1 (vi) step 2 by the computable rule of the repair plan revision 3: 5400/0.005 (48 × 2.25e6) iff P(Bin(3200, pw) ≥ 2909) ≥ 0.99 for both surrogates, where pw = OC_paired(c_p) × OC_single(c_s) is the exact binomial OC at the measured κ = 1 coverages of the rebuilt pool (the measured sandbox single-gate power replaces the modelled factor when lower); otherwise 7200/0.005 (64 × 2.25e6) iff the same computation gives ≥ 0.99 there; otherwise no V7-R2 design is frozen and the finding goes to the operator (repair plan revision 3, §3.8 item 2). *(Superseded by revision 6: the measured joint-rate rule of repair plan revision 4.)*
- **R2-2 (budgets).** Updated to the repair plan revision 3 (§6.2 there):
  - V7-R2 33.9 h (5400/0.005) or 44 h (7200/0.005); (B) 16.5 h; (C) 12.5 h.
  - Re-qualification 37.2 / 47.3 / 19.8 / 15.8 h.
  - Totals with D-E4 41.3 / 51.5 / 23.9 / 19.9 h (§4.2; exact sums in Appendix A).
- **R2-3 (layout).** The stale "24 × 4.5e6" in §2.2 and §2.4 is replaced by 48 × 2.25e6.
- **R2-4 (fallback seeds).** The 7200/0.005 fallback seeds are stated: 64 shards k = 0..63 and the reference k = 64, i.e. 20547004–20547068. The reservation of hr 20541004 r_index 6 covers k ≤ 64.
- **R2-5 (citations).** The repair plan is cited as revision 3.
- **R2-6 (fixture).** The halved rehearsal fixture is declared in the (v) amendment text.
- **Paths.** Cross-references use the committed paths:
  - `validation/plans/v3-005d-integration-route.md`;
  - `validation/plans/v3-005d-elastic-attribution.md`;
  - `validation/plans/v3-005d-v7r-calibration-repair.md`;
  - their `review.md` and `scripts/` under the same-named directories.

## Revision 4 (2026-10-14)

- **Bases.** lv **20531004** / hr **20541004** are adopted as the re-qualification bases, and V6-F (Option II) goes to **20551004**. This is consistent with the V7-R repair plan revision 2. They satisfy the frozen range rule of Amendment 18 (b)5 without a proof.
- **Skipped bases.** 20511004 / 20521004 are retired unused: declared skipped before any observation, no seed drawn, no seed repetition. The revision-3 seed-by-seed result (no shared seed for k ≤ 995) is kept as the recorded reason that they would also have been admissible with a proof (§4.4).
- **Seed numbers aligned** with the V7-R plan revision 2 and the elastic plan revision 4 (§4.4 table):
  - V7-R2 transport 20547004–20547052, layout 48 × 2.25e6;
  - hr r_index 6 reserved;
  - V6-F r_index 7–9;
  - calibration and κ 20980000–20984999;
  - 2099xxxx from 20990100.
- **V7-R2 power re-check rule** taken from the repair plan revision 2. *(Superseded by revision 5, R2-1: the computable rule of repair plan revision 3.)*

## Revision 3 (2026-10-14)

Answers the independent review `validation/plans/v3-005d-integration-route/review.md` (verdict REVISE: B1, M1–M5, m1–m7) and incorporates the conclusions of the V7-R plan review (`validation/plans/v3-005d-v7r-calibration-repair/review.md`, F1–F4). Text, planning arithmetic and one read-only seed scan only. No row verdict, tolerance or recorded failure changes.

- **B1 (V5 fresh base is a decision, not a statement).** Amendment 17 (b)4 says that once the V5 stop rule fires, "no fresh base is taken for V5". It fired at dd976c5a. Re-running V5 at a fresh base is therefore operator decision item 5 of (iv), quoted verbatim, with a default if refused (§1.2, §4.1 (iv), §6 S1).
- **M1 (Options I/II reduce (a)6).** §2.2 now states plainly that Options I and II amend (a)6: V6 and V10-A are removed as gates, and V7-R2 is substituted for V7-R. These are scope reductions only the operator can make. They override the merge default of Amendment 18 (c), which leaves that to the operator. The V3-005A precedent (#81) had archive `pass: true`; this run's archive will be `pass: false` (§2.2, §3.1).
- **M2 (suite and seed-base plumbing).** §4.1 (iii) lists the code changes the re-qualification needs:
  - suite identity;
  - qualification-base constants and `seed_blockers`;
  - step list, preflight and single-process force list;
  - the V7-R2 reference index (`V7R_REF_INDEX = 16` collides with shard 16 of a 48-shard layout);
  - `declared_plans` entries or a timing rule for the calibration and κ parts.

  The CLI suite names are corrected (§3.2).
- **M3 (surrogate pool regenerated).** The V7-R2 surrogate *definitions* are unchanged. The rehearsal-pool *content* is regenerated at the evidence digest after the model change. A pre-declared sandbox re-check of κ = 1 coverage and power on the rebuilt pool precedes the freezing of the V7-R2 sizes (§2.4, §4.1 (vi)). The fixture shards are halved, with rehearsal seeds 20521000–20521003 (m5).
- **M4 (Option II binding).** Once the (v) amendment is frozen, the Option II choice is binding. A V6-F failure then blocks integration. The rev-2 S7 fallback is removed (§2.2, §6 S7).
- **M5 (seed-by-seed evidence).** The check was performed (`validation/plans/v3-005d-integration-route/scripts/used_seeds.py`, `seed_check.py`; §4.4).
  - Recorded rehearsal transport seeds: 20506000, 20517000, 20521000, 20521001, 20525000; plus the declared E-shape 20509500 + k.
  - **The pre-declared bases lv 20511004 / hr 20521004 share no used seed for any k ≤ 995.** The first shared values need k = 996–999.
  - Revision 2's claim that they "collide" overstated the evidence; it is withdrawn. Only the V6-F range criterion (B + 7000 > largest used) fails for B = 20511004, and V6-F does not use that base.
  - Recommendation: keep the Amendment 14 (j) bases 20511004 / 20521004 with the seed-by-seed proof and a test. The rev-2 instruction that Amendment 19 must not name them is withdrawn.
  - Allocation: lv 20511004, hr 20521004, then V6-F (Option II) 20531004; the sequence skips 2098xxxx / 2099xxxx. *(Superseded by revision 4: lv 20531004 / hr 20541004, V6-F 20551004.)*
- **m1.** V7-R2 design, pre-registration, implementation and sandbox work are autonomous. Its host execution (≈ 33.7 h at revision 3 of this route; ≈ 33.9 h per repair plan revision 3) happens only in (vi), after (iv) and (v), with its budget in (iv). IR-DEFAF4 had quoted 9.5 h for V7-R.
- **m2.** The throughput step runs first in (vi).
- **m3.** §3.1: `conformant` also requires `pass`, and missing documents make `pass` false.
- **m4.** One ordered base allocation (§4.4, §6 S8).
- **m5.** See M3.
- **m6.** The calibration is 39 parts per surrogate × 82/83 repetitions = 3200, 78 parts in total; the cost is unchanged.
- **m7.** Under Option II the V6-F code enters the hashed source set. It must land before the fixture rebuild (§4.1 (iii), (vi)).
- **V7-R plan review.** V7-R2 is a new row, and the repair plan's citation of 14 (h)4 for the shard increase is void, because that permission expired at the evidence run. Transport shards are re-sized: 48 × 2.25e6 (≈ 631 s) here, or the reviewer's 36 × 3e6 (≈ 841 s); either is fixed in the V7-R2 amendment. The power clause is re-measured on the regenerated surrogate.

## Revision 2 (2026-10-14)

Changes from revision 1. They are text and planning arithmetic only. No row verdict, tolerance or recorded failure changes.

1. **V7-R calibration-repair plan integrated** (`validation/plans/v3-005d-v7r-calibration-repair.md`).
   - Diagnosis in §1.2. Design V7-R2, the cheaper variants (B) and (C), and their costs in §2.3, §2.4 and §4.
   - New §2.4 assesses that plan's claim that V7-R2 "needs no operator approval".
   - Option I's conditional ("only if a sandbox study shows …") is replaced by the plan's quantitative answer.
   - Finding: the plan's 4.5e6-history V7-R2 transport shards (≈ 1262 s) exceed the 1200 s `shard_table_ok` limit of Amendment 14a (c)5. They must be halved (§3.2, §4.2).
2. **Elastic attribution plan revision 4** (accepted after three review rounds, `validation/plans/v3-005d-elastic-attribution/review.md`) is referenced:
   - rows D-E0–D-E6;
   - 10 TOPAS thin-target cases (rows 14–23; 2e7 primaries with 2 mm targets at 30/50 MeV, 4e6 with 1 cm targets at 100–200 MeV), each ≤ 1500 s, after a 1e5-primary smoke case;
   - investigation seeds from 20990100;
   - M2 (optical-limit Glauber) pre-committed.

   §4 phase costs are updated accordingly.
3. *(Superseded by revision 3, M5: the seed-by-seed check shows no shared seed; the bases 20511004 / 20521004 are kept.)* **Fresh-base collision stated explicitly** with the arithmetic (§4.4). The declared bases lv 20511004 / hr 20521004 overlap the rehearsal seeds actually used. The re-qualification amendment must declare **20531004 / 20541004 (or later)** with a seed-by-seed check. The closing line and seed list of the elastic plan's Amendment 19 draft, and the V7-R plan's transport-seed paragraph, must not name 20511004 / 20521004 as re-qualification bases. The V6-F base under Option II moves to 20551004 (candidate).
4. Budget (§4.2), critical path (§4.3), operator items (§2.2, §2.3, §4.1 (iv)) and stop rules (§6) are revised to match.

## 0. Binding frame

- **Operator disposition A**, IR-20261013-164616-DEFAF4 (`reference_requests/interventions/IR-20261013-164616-DEFAF4/response.json`, verbatim): "Proceed with option A: keep PR #82 unmerged and continue dependent follow-up work from its unmerged branch state. Preserve all recorded failures and do not create an implementation-only passing record. Prioritize elastic-model attribution and V7-R calibration repair before another full qualification run. Before expensive qualification, document and independently review the route to integration, including remaining requirement coverage, V6/V6-F treatment, and execution conformity under the existing stability restrictions. V6-F passing does not authorize replacing V6 or merging. Bring back any concrete proposal that changes required validation scope; routine scientific investigation and implementation should continue autonomously."
- **Earlier operator conditions still binding.**
  - IR-20261009-032741-C590CF (Amendment 14 (a)): interpretation B; pre-register before results; no threshold relaxation; no seed repetition to obtain a pass; preserve failed and inconclusive records; no implementation-only passing record; merge only after the applicable scientific and implementation gates pass.
  - IR-20261011-031924-C205C8 (Amendment 18 (a)): V6/V6-sens stay failed at the sizing stage; applicable merge gates stay enforced.
  - Single-process directive (acceptance plan, "Execution under the single-process directive"; Amendment 14 (k)): `cpu_workers=1`, one compute job globally, OMP/OPENBLAS/MKL/NUMEXPR = 1, `IONMC_SINGLE_PROCESS=1`, CUDA through one controlling process, multiprocessing checks DEFERRED.
  - Operator diagnostic rules of 2026-10-12 (memory note `v3-environment-constraints.md`, "Operator operational rules"): one heavy local job at a time; complete external logs under `telemetry/test-logs/<UTC>-<label>/`; instrumented pytest; after a second crash, narrow first.
- **Plan rules used throughout.**
  - Later changes need a new dated amendment (acceptance plan header).
  - A failed base is consumed (Seeds; Amendment 14 (a)4).
  - Fresh bases after a reviewed fix were pre-declared: lv **20511004**, hr **20521004** (Amendment 14 (j)). *Revision 4:* the re-qualification uses lv **20531004** / hr **20541004**; 20511004 / 20521004 are retired unused (§4.4).
  - Exactly one fresh base per reviewed fix (Amendment 14 (f)).

## 1. Requirement coverage at dd976c5a

### 1.1 Merge gate (a)6, exact wording

Amendment 14 (a)6: "The task is merged only if the scientific gates (V5, V6 with the status fixed by V6-sens, V7-R including its calibration step, V8-nuc, V10-A) and the implementation gates (P6, P7, P5-ext, C1-ext, R1-D2, V8-LV, CI) pass. A summary that passes only the implementation rows is not a passing record of this task, and none is created."

Modifications already in force:
- P6 is replaced by **P6-D** in the gate list. P6 stays FAILED for H-1 below 12.53 MeV (Amendment 15 (a)5).
- Amendment 17 (c)7: if V6-sens fails its condition, (a)6 treats V6 as "status report-only, recorded".
- Amendment 18 (c): "A row recorded as 'failed — not sized under the pre-registered rule' is a failed scientific gate: it is not a pass and not report-only, and V6-F cannot stand in for it. V3-005C therefore cannot be merged on V6 under (a)6." Default without an operator decision: NOT merged.

### 1.2 Gate entries that cannot be satisfied as written

| (a)6 entry | Why it cannot pass as written | Source |
|---|---|---|
| **V6 (with V6-sens)** | Recorded "failed — not sized under the pre-registered rule". `V6_SHARDS` stays `None`. The r = 0 clause needs ≥ 8910 shards at z = 22.4 cm, and ≥ 1.3e10 histories at z = 22.6 cm (≥ 19 and ≥ 28 days of host time at 5437 hist/s). V6-F cannot stand in (Amendment 18 (b), (c)). The operator also stated that V6-F passing does not authorise merging. **No future run can turn V6 into a pass.** | Amendment 14a (a)2–4; 18 (b)–(c); IR-DEFAF4 |
| **V10-A clause (c)** | It fails on the σ_el level, independent of the angular shape. σ_el(p+C, 100 MeV) = 253 mb (Barashenkov, shared with TOPAS; P6 (5) value) against the post-1997 proxy 216 ± 21 mb, i.e. +17 % > 12 %. A shape change cannot repair (c). Changing σ_el at a seen point would be tuning, and is excluded (Amendment 14 (a)3). | elastic-attribution-plan §0, §7.2; Amendment 14 (d) V10-A, (e)2 |
| **V10-A (a)/(b) after a shape change** | Any new p+A shape is chosen after V10-A was observed at dd976c5a (step 03, and already at 37e87ecf). The pre-registered (e)9 logic for the analogous (ii′) case says: "All V10-A results for (ii′) are exploratory, because the data have been seen. Gating for the shape then rests on V6 and V11." V6 is permanently failed. So a new-model V10-A cannot be a blind gate. | Amendment 14 (e)9; elastic-attribution-plan §7.1 |
| **V7-R (including its calibration step)** | The calibration of record failed at dd976c5a. The frozen outcome rule reads "Fail: a recorded scientific failure **without a further base**" (Amendment 14 (h), Outcome). The calibration is deterministic (Amendment 14a (b)6), so re-running the same design reproduces the same counts. Amendment 14a (b)4 also records that the gate-alone clause (CP95 upper ≤ 0.04 for a gate built at α_tost = 0.04) "may be unprovable by Monte Carlo at any n". The rehearsal block sums give 11.5 % / 11.6 % gate-alone acceptance at the margins. **The frozen V7-R design cannot pass, and the frozen text forbids a further base for it.** A repaired design is a new row. *Revision 2 (repair plan §1–§2):* the defect is quantified. The paired gate at α = 0.04 with 3600 pairs has an exact false acceptance of 0.0383 at the margin. The ≤ 109/3200 gate-alone clause therefore passes with P ≈ 0.11 per margin even with a perfect κ (≈ 3·10⁻⁴ for all four margins). The κ from a 2·10⁴ sample has a coverage SE of 0.0034–0.0039. The '11.5 %' is the gate's correct acceptance at a true coverage of ≈ 0.658, where the κ of record landed. All κ of record are reproduced bit for bit. The real data, 2467/3600 = 0.6853 (CP [0.6715, 0.6988]), are consistent with a correct uncertainty, but that is not a pass. | Amendment 14 (h); 14a (b)3–6; README "V3-005C results"; v7r calibration-repair plan |; `validation/plans/v3-005d-v7r-calibration-repair.md`
| V5 | **Cannot be re-run without an operator decision.** The (f) stop rule fired at dd976c5a (peak/plateau low at both energies). Amendment 17 (b)4 then says "no fresh base is taken for V5, no tolerance changes, and the orchestrator raises an intervention". Re-running V5 at a fresh base after the attribution is therefore operator decision item 5 of §4.1 (iv), not a statement. Default if refused: V5 stays failed and §6 S1 applies. | Amendment 14 (f); 17 (b)4 |

**Consequence.** Amendment 14 (a)6 cannot be met by any future run as written. Its V6 entry is permanently failed. Its V7-R entry cannot pass under the frozen outcome rule. Its V10-A entry fails (c) on a shared quantity that the plan forbids tuning. Integration of V3-005 therefore requires an operator decision on how (a)6 is carried (§2, §4 (iv)). Until that decision is made, the default stays: PR #82 unmerged (Amendment 18 (c)).

### 1.3 Row table (status at dd976c5a unless stated)

Status vocabulary:
- **pass**
- **FAIL**
- **failed-not-sized**
- **deferred**
- **report-only** (executed, no verdict, or a verdict that does not gate)
- **pending** (pre-registered, not executed)

"Integrable" means what each row needs at the re-qualification evidence SHA.

| Row | (a)6 role | Status | Evidence pointer | Needed for an integrable state |
|---|---|---|---|---|
| P1–P5, N1, V1, V3, V4, V4b, V9, R1, X1, C1, D6 (slice A) | merged with V3-005A (#81, 087a70a) | pass at b84fdf3 (D6: both tiers failed; local α deposition under a declared ceiling, decision 0041 §7) | `lv5-b84fdf3-single-process.json`; acceptance archive table | Nothing new. R1 is re-checked by R1-D2 (A) at the new SHA. |
| V1b, E1 | – | report-only | `lv5-b84fdf3…` | – |
| V2b (slice B, pre-elastic) | – | pass at 7aae5bb | `lv5b-7aae5bb…` steps 09–11 | Superseded in scope by V2b-el; stays archived. |
| V7 replicate coverage (slice B) | – | **FAIL** at 7aae5bb (`escaped_neutral` paired gate, lower bound 0.6445 < 0.6464) | `lv5b-7aae5bb…` steps 14–24 | Preserved. Re-tested by V7-R. |
| V8 `nuc_local` (slice B) | – | inconclusive at 7aae5bb | `hr5-7aae5bb…` | Preserved. Answered by V8-nuc. |
| **P6** (frozen) | replaced by P6-D | **FAIL** (H-1 item (4) below 12.53 MeV); items (1)–(3), (5)–(7) pass | Amendment 15 (a)1; ledger S-PHYS-PROTON-NUCLEAR | Preserved; disclosed in the PR and the merge-gate report. |
| **P6-D** | implementation gate | pass (CI) | IR-DEFAF4 summary; CI run 38431680322 at bb8c198f | Re-pass in CI at the new head. Re-check if the builder changes (the table id changes with a new shape). |
| X-MT5-SPEC | qualification flag | pass (KS D 0.042–0.230, 9 of 9) | ledger; Amendment 16 | Re-pass if the table is rebuilt. |
| **P5-ext**, **C1-ext** | implementation gates | pass (CI) | IR-DEFAF4 | Re-pass in CI. New twins or fail-closed cases are needed for any new shape code (the U1 harness covers every `@wp.func`). |
| **P7** | implementation gate | pass (step 10) | lv5c archive `rows.P7` | Re-run at the fresh base. The tables change with the shape. |
| **R1-D2** | implementation gate | pass (step 08) | lv5c `rows.R1-D2` | Re-run. (A) and (B) must stay bit-identical. A new (C) intended-change record is needed for the shape change (directions re-declared before the run). |
| **V8-LV** | implementation gate | pass (step 09) | lv5c `rows.V8` | Re-run. |
| **CI** | implementation gate | green at bb8c198f (after the data-availability skips, Amendment 18 (e)) | CI run 38431680322 | Green at the final head. |
| **V5** | scientific gate | **FAIL** (peak/plateau −3.66 % / −4.89 %; ΔIDD +10.35 % at 150 MeV); stop rule triggered | lv5c steps 21–22; README | Preserved. Re-run at lv 20531004 only if the operator grants a fresh V5 base (§4.1 (iv) item 5), after the reviewed model change, with tolerances unchanged. |
| **V6, V6-sens** | scientific gate | **failed-not-sized** | lv5c steps 28–29, hr5c 10–11; Amendment 14a (a)4 | Cannot pass (§1.2). Operator item 1 (§2). |
| **V7-R** | scientific gate | **FAIL (calibration)**; gates pass (0.685 / 0.674 / 0.673) | lv5c steps 30–80, 78b; README | Cannot pass as frozen (§1.2). New row V7-R2 (repair plan §3.4; autonomous as a row, §2.4); its (a)6 slot is operator item 3. |
| **V8-nuc** | scientific gate | pass (T12 conclusive, p 0.50, 160 groups) | hr5c step 08 | Re-run at hr 20541004 (§4.4). Under its unchanged rule an inconclusive verdict is not a pass, and no larger sample follows (Amendment 14 (d)). |
| **V10-A** | scientific gate | **FAIL** ((a) at 15°/20°, (b), (c)); trigger did not fire | lv5c steps 03, 11 | (c) cannot pass; a new-model (a)/(b) is data-seen (§1.2). Operator item 2. |
| V11 | not in (a)6 ("report with tolerance") | **FAIL** (F 0.857 / 0.844 vs 0.892 / 0.887; gain +0.006 / +0.008 above TOPAS) | lv5c step 16 | Re-run with the new model. Under the recommended option 2B, it becomes a gate (an addition to (a)6, operator item 2). |
| V2b-el, V7-shift-el | not in (a)6 | pass | lv5c 25, 27 | Re-run (the model changes the step machinery inputs). |
| V7-scan-el | not in (a)6 | **FAIL** (`sec_p_dose` slope −0.5526 vs [−0.55, −0.45]) | lv5c 26 | Re-run at the fresh base. The failure stays recorded. A second fail narrows the V2-NUM claim for `sec_p_dose` (§6). |
| V10, X-ENDF, D9-BGG, E1-B | report-only | executed | lv5c 02, 05, 06, 84 | Re-run (deterministic or cheap). |
| E-shape | exploratory | executed (85b direct; archive row `pass: false` from wrapper run 85) | README "Incidents" | Re-run through the fixed wrapper (C11) so that the archive row is complete. |
| V6-F, V6-sens-F | not in (a)6 (Amendment 18 (b)) | pending (pre-registered, independently reviewed; not executed) | `validation/plans/v3-005c-v6f-followup.md` | Executed once, gating, in the task that makes the lateral-halo claim. See §2. |
| `v3-workers-partition` (1-vs-N) | conformance | **deferred** | lv5c step 86 | Needs the directive lifted (operator); see §3. |

### 1.4 Requirement IDs touched by V3-005

Source: `validation/requirement-ledger.md`, rows quoted at dd976c5a.

| Requirement | V3-005 contribution | Status at dd976c5a | Gap to an integrable claim |
|---|---|---|---|
| S-PHYS-PROTON-NUCLEAR (suite) | the whole row set above | non-elastic (A) pass; elastic channel (C): V5, V11, V10-A FAIL; V6 failed-not-sized | The elastic claim needs §4. The halo claim is withheld until V6-F. |
| V1-MUST-005 Ion transport | proton nuclear (non-elastic + elastic) on python, warp-cpu, warp-cuda | "in progress"; parity rows pass; V5, V11, V10-A fail | V5 pass at the fresh base. The p+A elastic claim scope follows operator item 2. |
| V1-MUST-013 Absorbed dose | nuclear dose contribution; V5 and V6 dose comparisons | V5 FAIL on peak/plateau, ΔIDD 150 MeV; plateau, R80 and total deposit pass; V6 failed-not-sized | V5 pass. **Measured-dose (V6) evidence gap stays** under every option of §2. |
| V1-MUST-014 Energy deposition | nuclear bookkeeping tallies (V3), recoil `nuclear_local` | pass (slice A V3; R1-D2) | Re-run R1-D2 / V3-relevant rows. |
| V1-MUST-029 Cross-backend statistical consistency | V8, V8-nuc | pass | Re-run hr5c'. |
| V1-MUST-030 Statistical uncertainty estimation | V7, V7-R, V7-scan-el | V7-R calibration FAIL; V7-scan-el FAIL | V7-R2 (operator item 3), or a narrowed claim. |
| V1-MUST-040 Backend consistency | V8-LV, V7 f32/f64 CUDA | pass | Re-run. |
| V2-CAP | C1, C1-ext, the elastic-domain counters | pass | New fail-closed cases for the new shape table. |
| V2-NUM | V2b-el, V7-scan/shift-el, V7-R | V2b-el and V7-shift-el pass; V7-scan-el and V7-R fail | As above. |
| V2-EVID | lineage tables (Amendment 14 (c)); "no fully independent evidence for the elastic channel exists in this scope" (14 (c)4) | elastic evidence: V5/V11 share σ_el with TOPAS; V10-A related lineage; V6 failed | The claim must be narrowed to the lineage that remains (§2, option 2B). |
| V2-RELEASE, V1-MUST-012 (conformance) | – | archives `conformant: false` by code | Release-level, not task-level (§3). |

## 2. V6 / V6-F treatment and the operator decision

### 2.1 Fixed facts (no decision needed)

1. **V6 and V6-sens stay "failed — not sized under the pre-registered rule"** in the ledger, the results README, the decision 0041 outcome and the merge-gate report. This holds at dd976c5a and at every later evidence SHA: the steps `v6-unsized` are deterministic and re-record the same verdict (Amendment 14a (a)4; Amendment 18 (b), (c); v6f-followup §1).
2. **V6-F is a separate, new row.**
   - It is pre-registered (revision 2, independently reviewed: `validation/plans/v3-005c-v6f-followup.md`).
   - It is executed **once**, as a **gating** row, "in a later task that makes the lateral-halo claim", on that task's fresh base B with seeds B + 7000/8000/9000 + k (Amendment 18 (b)5; v6f-followup §7, §8.2).
   - Its pass is reported only as "V6-F (reduced scope, item 4) passed", never as "V6 passed" or "halo validated" (v6f-followup §1.4).
3. **V6-F passing never replaces V6 and never authorises merging** (operator, IR-DEFAF4; Amendment 18 (b)). P6-D is explicitly not a precedent (v6f-followup §1.2).
4. **Base constraint for V6-F.** V6-F needs its own base B with B + 7000 > 20525999 (the block bound of the rehearsal seeds; largest recorded 20525000), or a seed-by-seed proof, and checks (a)–(e) (v6f-followup §7). It cannot share the lv base, because its r_index 7–9 coincide with lv rows. Under Option II, B = **20551004** (B + 7000 = 20558004; §4.4).
5. **What V6-F can establish.**
   - It can establish r = 0 absolute dose for z ≤ 21.6 cm, and the charged-secondary halo at r = 1–5 cm, at core-insensitive points with aura share ≤ 5 %.
   - It cannot establish the frozen V6 claim, the aura or any neutral dose, r ≥ 6 cm, the tail beyond R0 + σ_sem, independence from the beam fit, or any other energy (v6f-followup §10).

### 2.2 What the operator must decide

Merge gate (a)6 cannot be met as written (§1.2). The operator decides **how (a)6 is carried for V3-005**. The question covers three entries together, because each one blocks integration:
- V6 (permanently failed);
- V10-A (clause (c) on the shared σ_el level; a data-seen evaluation after a shape change);
- V7-R (the frozen design cannot pass; "without a further base").

None of the options relaxes a tolerance, re-labels a recorded verdict or reuses a seed. **But Options I and II reduce the merge-gate set of (a)6 itself** (review-route M1): V6 and V10-A stop being gates and become carried failures, and V7-R2 replaces V7-R. They override the Amendment 18 (c) default ("NOT merged, PR #82 stays open"), which leaves exactly that decision to the operator. They are scope reductions that **only the operator can make**. Unlike the V3-005A precedent (#81, archive `lv5-b84fdf3` with `pass: true`), the re-qualification archive will carry `pass: false`, because `v6-unsized` re-records V6 as failed.

**Option I — integrate with the recorded failures carried and the claims narrowed (recommended).**
- *V6.* Carried in (a)6 as "failed — not sized, recorded". V3-005 makes **no lateral-halo claim and no measured-dose claim** against Gottschalk 2015. Decision 0041 limitation 7 ("No lateral-halo claim …") is re-stated with "until V6-F is executed in the halo-claim task". V6-F stays scheduled there, as Amendment 18 (b)5 already says.
- *V10-A.* The dd976c5a failure (including (c)) is carried as recorded.
  - The new-model V10-A is computed with the frozen tolerances and **reported** with a "data seen" disclosure. This follows the pre-registered (e)9 precedent ("All V10-A results for (ii′) are exploratory").
  - The p+A shape is then gated by **V11 added to (a)6**. V11 is currently "report with tolerance"; adding it is a tightening, unchanged in tolerance.
  - The claim is narrowed to: "elastic angular effect agrees with Geant4 `g4h-elastic_HP` within the V11 tolerances; σ_el(p+A) level is the Barashenkov compilation shared with TOPAS; p+C σ_el at 100 MeV is +17 % above the post-1997 proxy (recorded deficiency, analogous to the pre-declared LA150 p+C deficiency of V1b)".
- *V7-R.* Carried as recorded failed. The repaired, separately pre-registered **new row V7-R2** stands in its (a)6 slot. The repair plan has now answered the feasibility question offline (§2.4): 48 shards × 2.25e6 (5400 pairs), α_tost = 0.005 on both `escaped_neutral` gates, κ from 1e6 samples, R = 3200. Expected false acceptance 0.0050. The measured joint pass rate at κ = 1 on the current pool is 383/400 and 384/400 (CP95 lower 0.9369 / 0.9399), so the power clause holds with P ≈ 1.000. It is re-measured on the regenerated pool (§2.4); ≈ 33.9 host-hours (repair plan revision 3). Design, pre-registration, implementation and sandbox work are autonomous. Host execution happens only in (vi), with its budget in (iv) item 8. The slot substitution is the proposal (§2.4).
- *Consequences.*
  - Integration becomes reachable with one re-qualification (§4.1 (vi): ≈ 37.2 h host with V7-R2 as written, ≈ 15.8 h with calibration variant (C)).
  - The ledger shows V6, V10-A and V7-R as failed next to the passes.
  - V1-MUST-013 keeps its measured-dose evidence gap.
  - V2-EVID records that the elastic channel has model-to-model evidence (V11, shared σ_el lineage) and related-lineage data consistency only (Amendment 14 (c)4 already states that no fully independent evidence exists in scope).
  - Risk: V11 is the only gating shape test, and it is a Geant4 comparison. The D-E4 thin-target study also uses Geant4, so it adds no new lineage.

**Option II — as Option I, plus V6-F executed inside V3-005D as an additional gate.**
- V3-005D declares itself the task that makes the lateral-halo claim. It executes V6-F once, gating, on its own base B = 20551004, which passes the v6f-followup §7 checks (§4.4).
- V6 is still carried as failed: V6-F does not replace it.
- *Consequences.*
  - +72 CPU shards ≈ 3.68 h at the Amendment 18 rate, i.e. ≈ 6.2 h at the slower host rate (§4.2), plus the CUDA replicate (rate not assumed, v6f-followup §8.3).
  - It adds measured-data evidence for the halo, which is sensitive to elastic and non-elastic secondaries and MCS tails (v6f-followup §10). This is the only measured test of the new elastic shape in scope.
  - New risks: V6-F is a new failure surface. Its blindness is weaker than that of a row frozen before the Amendment 14a (a)5 display (Amendment 18 (b)5). The new shape changes the halo, so a V6-F fail would block integration again.
  - Needs the operator's explicit agreement, because the response tied V6-F to "a later task", and V6-F would gate here.
  - **Binding:** once the (v) amendment is frozen, the Option II gate set is fixed. A V6-F failure then blocks integration. No fall-back to Option I after seeing the result; a later change would be a new, post-hoc operator decision (review-route M4).
  - The V6-F code enters the hashed source set, and so the rehearsal-fixture digest. It must land before the fixture rebuild (m7).

**Option III — keep (a)6 literal.**
- PR #82 stays unmerged until every listed entry passes. V6 can never pass (§1.2), so under this option V3-005 never integrates.
- *Consequences.*
  - Proton nuclear physics (non-elastic V3-005A is already on develop; elastic and Warp nuclear are not) stays off `v3/develop`.
  - Every dependent task must branch from the V3-005B branch state. `create_task_worktree` cannot do that (§5).
  - The release cannot carry the elastic channel.
- Option (B) of IR-DEFAF4 (non-qualified integration) is not re-proposed, because the operator chose A.

**Recommendation: Option I.**
- It honours Amendment 18 (b) and the operator's sentence on V6-F: no V6-F result is used to merge. It requires the operator to amend (a)6, which 18 (c) leaves to the operator.
- It spends no host time on rows whose outcome cannot change a gate.
- It narrows each claim to the evidence that exists (V2-EVID) instead of inventing a substitute gate.
- It keeps V6-F where Amendment 18 put it.

Option II is the better scientific choice **only if** the operator wants measured-data evidence for the elastic channel before integration. Its marginal cost (≈ 6.2 h host CPU) is about a sixth of the re-qualification with V7-R2 (≈ 37.2 h). If the operator selects II, the V6-F rows enter the pre-registration amendment of §4 (v) with B fixed there.

### 2.3 Items 2 and 3 in isolation (if the operator prefers to decide per row)

| Item | Option A | Option B | Option C (literal) | Recommended |
|---|---|---|---|---|
| 2: V10-A | Count the new-model V10-A (a)/(b) toward (a)6 at frozen tolerances with a "data seen" disclosure; (c) carried failed, with a σ_el claim restriction | dd976c5a V10-A carried failed; new-model V10-A reported (data seen, (e)9 precedent); V11 added to (a)6 as the shape gate; σ_el claim restricted | cannot pass: (c) fails on shared σ_el; tuning is excluded | **B**. Revision 4 of the attribution plan pre-commits M2. But it declares every later use of D-E1 'data-informed with respect to V10-A' (plan §4), so a new-model (a)/(b) would not be blind. |
| 3: V7-R slot | V7-R2 (pre-registered repaired design, §2.4) stands in the (a)6 slot; V7-R stays failed | V7-R carried failed; `escaped_neutral` coverage claim withheld (V1-MUST-030 narrowed: SE reported, coverage not qualified); V7-R2 may still run and is reported | cannot pass ("without a further base"; deterministic calibration) | **A**. The repair plan shows P(all clauses) ≈ 1.000 for V7-R2 with clauses as written, computed on the old pool and re-checked on the regenerated pool before the sizes are frozen (§2.4). The calibration verdict comes before the transport. |
| 3′: calibration variant (only if the operator wants a cheaper run) | V7-R2 as written (design autonomous; host budget ≈ 33.9 h is an operator item) | (B) exact paired certification + R = 800 (≈ 16.5 h; decreases R, analytic clause) | (C) deterministic certification (≈ 12.5 h; analytic false-acceptance clauses, power by Monte Carlo at R = 400) | **V7-R2 as written**: no clause change. (B)/(C) only by operator decision. |

### 2.4 What is autonomous about V7-R2? (assessment of the repair plan's claim; the V7-R plan review agrees)

**The repair plan's original claim** (revision 1, §3.4, §4 there; revision 3 has since adopted the classification below): V7-R2 needs no operator approval, because it changes only:
- the α_tost of the two `escaped_neutral` gates (0.04 → 0.005, stricter);
- the transport size (16 shards of 4.5e6 → 48 × 2.25e6, the same total as 24 × 4.5e6);
- n_κ (2·10⁴ → 10⁶, full sample);
- the step layout and the seeds.

The margins [0.6527, 0.7127] and the clauses (joint ≤ 0.02, gate alone ≤ 0.04, power ≥ 0.9, R = 3200, B = 1999, the surrogate definitions) are unchanged.

**Frozen wording it has to be read against:**
- Amendment 14 (h), Outcome: "Fail: a recorded scientific failure without a further base."
- 14 (h)4: increases are allowed "in Amendment 14a before the evidence run, never after and never decreased. If the calibration step at the evidence SHA misses a requirement, V7-R fails …, whatever the gates give."
- 14a (b)6: "a miss there fails V7-R … and no further count increase follows."
- The plan header: later changes need a dated amendment with the scientific reason.
- Operator C590CF: pre-register follow-up validation; do not repeat seeds to obtain a pass.
- Operator DEFAF4: "Prioritize … V7-R calibration repair"; "Bring back any concrete proposal that changes required validation scope".

**Assessment:**
1. **Autonomous:** the diagnosis record, and the design, pre-registration (dated amendment), implementation, review and sandbox work of V7-R2 **as a new row**.
   - The frozen sentences bind V7-R, which stays failed with its counts.
   - V7-R2 has a stated scientific reason: a defect in the design of the calibration. The finite-sample size 0.0383 of the α = 0.04 test at 3600 pairs makes the ≤ 109/3200 gate-alone clause pass with P ≈ 0.11 per margin even with an exact κ (P(≤ 109 | 0.0374) = 0.17, per the V7-R review F8). The κ-sample SE of 0.0034–0.0039 makes every margin clause a lottery.
   - No tolerance changes, and the gate is stricter.
   - Precedent: V7 → V7-R (Amendments 11, 14), now explicitly prioritised by DEFAF4.
   - The repair plan's citation of 14 (h)4 as authority for the shard increase is void (V7-R review F1), and is removed. V7-R2's sizes are new-row parameters, not increases of V7-R.
2. **Not autonomous (operator, in (iv)):**
   - **V7-R2's role in merge gate (a)6 in place of V7-R** (item 3). Read literally, "without a further base" forbids re-testing the same coverage claim on a new base. That is the substance of V7-R2, even though the row is new. Amendment 18 and the operator refused substitution for V6.
   - **Its ≈ 33.9 h host budget** (item 8). IR-DEFAF4 quoted ≈ 9.5 h for the frozen design at the earlier rate. Host execution happens only in (vi), after (iv) and (v) (review-route m1; V7-R review F1).
3. **Operator proposals:** the cheaper variants (B) and (C). (B) decreases R below 3200, which "never decreased" forbids. Both evaluate a false-acceptance clause analytically, the "construction reading" that 14a (b)4 withdrew. Thresholds are kept, but the clause text changes.
4. **Not available at any price:**
   - re-evaluating the observed dd976c5a gate data under a new α;
   - repeating the dd976c5a calibration with more repetitions.

   The real-data result stays contrary evidence next to the failure, not a pass: 2467/3600 = 0.6853, CP bounds at α = 0.04 [0.6715, 0.6988].

**V7-R2 at a glance** (repair plan revision 3, §3.4, §5, §6.2; its review F2–F4, F8):
- **Transport** (repair plan rev 3 §3.4): **48 × 2.25e6** histories (225 replicates per shard, 10 800 replicates; pairs s ↔ s + 24, s = 0..23 → 5400 pairs; H = shards 36–47, E = 0–35; ≈ 631 s per shard, so `shard_table_ok` holds at the 1200 s limit); reference 1e6 at k = 48 (never k = 16, M2). Seeds 20547004–20547052 (§4.4).
- **Gates:** α_tost = 0.005 on both `escaped_neutral` gates; paired acceptance band 3615–3762 of 5400. `sec_p` / `nuclear_local` stay under the Amendment 13 rule.
- **κ and calibration:** κ from 1e6 samples in a sealed `v7r-kappa` step; R = 3200; 39 parts per surrogate × 82/83 repetitions, 78 parts (m6).
- **Seeds:** κ root 20980000, calibration root 20982000, A13 shard 20984000. The reserved retry block 20985000–20988999 is usable only under the V7-R review F6 conditions.
- **Surrogates** (M3, V7-R review F3):
  - The **definitions** are unchanged: the archive-fitted zero-inflated Gamma, re-fitted with seed 20505000 to k = 5.0703, π = 0.2167; and the rehearsal block-sum pool.
  - The pool **content** is regenerated at the evidence digest after the model change, so the pinned pool hashes 5d28a9b3… / d4f7d749… are dropped for evidence use.
  - The expected figures were computed on the old pool: false acceptance 0.0050 (paired, exact) and ≤ 0.0050 + 0.01 rigorous for the single gate (V7-R review F7); power 0.954 / 0.962 (repair plan revision 3 model values from the independence product; labelled, superseded by the measured joint rate below); P(all clauses) ≈ 1.000. They are **re-checked on the regenerated pool before the sizes are frozen** (§4.1 (vi) step 2).
  - The binding risk is the power clause. The pre-declared size-selection rule uses the **directly measured complete two-gate joint pass rate at κ = 1**: the repository `run_calibration`, with both gates on the same replicates, α_tost = 0.005, 48 × 225 replicates, R_meas = 400 per surrogate, seed roots 20989920–20989923 (planner block). Freeze 5400/0.005 (48 × 2.25e6) iff both surrogates reach **≥ 378/400** joint passes, i.e. the CP95 lower bound p_low satisfies P(Bin(3200, p_low) ≥ 2909) ≥ 0.99. Otherwise measure 7200/0.005 (64 × 2.25e6; roots + 10) and freeze it iff both reach ≥ 378/400. Otherwise no V7-R2 design is frozen and the finding goes to the operator (repair plan revision 4, §3.8 item 2). Measured on the current pool (repair plan revision 4, `validation/plans/v3-005d-v7r-calibration-repair/data/joint_power_5400.*`): archive-fitted 383/400 (CP95 lower 0.9369), rehearsal block sums 384/400 (0.9399), so **5400/0.005 is selected today**. The rule and R_meas were fixed before that measurement. The Fréchet bound (CP lower 0.9272 / 0.9310) is a cross-check only. The rule is re-applied on the rebuilt pool before the sizes are frozen (§4.1 (vi) step 2).
- **Host time:** ≈ 33.9 h (fixture 0.7, κ 1.1, calibration 23.2, A13 0.2, transport 8.7). 7200/0.005 (64 × 2.25e6): ≈ 44 h, with transport seeds 20547004–20547068 (k = 0..63, reference k = 64). (B): ≈ 16.5 h (22 h at R = 1600). (C): ≈ 12.5 h (+ a small end-to-end Monte Carlo).
- **Estimator risk:** a correct estimator at true paired coverage 0.69–0.695 passes the stricter gate with only 0.86 / 0.61.

## 3. Execution conformity under the stability restrictions

### 3.1 What the flags of the dd976c5a archives mean (verified with `json` on both files)

| Flag | lv5c-dd976c5 | hr5c-dd976c5 | Meaning |
|---|---|---|---|
| `pass` | false | false | Failed steps exist: lv5c `03, 16, 21, 22, 26, 28, 29, 78, 79, 85`; hr5c `10, 11`. `pass` is the conjunction of step verdicts. |
| `conformant` | false | false | `summarize` requires `pass` for `conformant`, so the failed steps alone make it false. Missing documents (`problems`) also make `pass` false. In addition: `non_conformant_reasons` (the next three rows). |
| ↳ imported partials | steps `11, 16, 21, 22, 25, 79, 80, 84`: `imported_partials_unverified_by_code` | step `08` | Amendment 5. A combine step that reads partials produced by other host calls cannot verify the protected host-runner records by code. The attestation is the orchestrator's validation record, which lists the host run ids whose `PARTIAL <name> <sha256>` stdout lines match the manifest. The archive's `attestation.protected_host_records_verified_by_code` is false. |
| ↳ deferred multiprocessing | `86-v3-workers-partition` (in `deferred_steps`) | – | 1-vs-N-worker partition invariance and any `cpu_workers > 1` row; "not passed and not part of a conformant qualification" (plan, "Deferred multiprocess checks"). |
| ↳ execution mode | `single-process-diagnostic` | same | "diagnostic, not the qualification mode". |
| `complete` | **false** | true | lv5c `problems`: steps 78 and 85 have no JSON result document. Wrapper defects 1 and 2: 78 lacked `--out-dir`, 85 built into the read-only cache. Both were executed directly as 78b / 85b at the same SHA; their documents are in the protected host records, not in the archive. Both defects were fixed in C11 (README "Incidents"). |
| `reduced_history_counts` | false | false | No count was reduced. |

Further deviations recorded in the README:
- The calibration records are sealed in-file with no `PARTIAL` line (a C8b omission). The manifest was built from the internal `content_sha256` and the host run ids (ledger `sealed_records`).
- One interpreter anomaly in calibration shard s27 (74 → 74b), with a clean deterministic re-run.
- Precedent: the V3-005A archive `lv5-b84fdf3` was equally non-conformant by code and was integrated (#81, 087a70a). It had `pass: true`; a future V3-005 archive will not (V6 is re-recorded as failed). So `conformant: false` under the directive has not blocked task integration before. It does block a release claim (V2-RELEASE; release report requires every suite "passed").

### 3.2 How the re-qualification run is executed so that its archive is reproducible and honest

1. **One host call per step**, through the CLI only (`run_suite.py --suite lv5c|hr5c --only <step>` today, or the successor suites `lv5d|hr5d` decided in §4.1 (iii); Amendment 14 (k)). One heavy job globally. No sandbox pytest or prototypes while a host step runs. CUDA through one controlling process. Thread limits 1 and `IONMC_SINGLE_PROCESS=1` are set by `run_suite --single-process` for every step subprocess.
2. **Step sizing at the slower host rate.**
   - Hard limit 1500 s per call: the MCP idle abort is at 1800 s (memory note; Amendment 14 (k)).
   - Effective planning limit 1200 s = 1500/1.25 (Amendment 14a (c)5). The throughput step at the evidence head must report `shard_table_ok` for every shard plan, or the run stops before the first evidence step.
   - Measured slowdown: identical V7-R rehearsal shards took 752–771 s up to dd976c5a and **≈ 1262 s** at a6ca56a4 (RUN-20261013T151818Z-9aa05382 and its sibling; orchestrator progress log 2026-10-13, after the operator disabled CPU idle states). The factor is 1262/755 = **1.67** in wall time, i.e. throughput −40 %.
   - Consequences, computed (§4.2):
     - V7-R(2) transport shards of 4.5e6 histories: 755 s → ≈ 1262 s, **above the 1200 s limit**. The repair plan's 24 × 4.5e6 layout (its §6.2 costs 24 × 1262 s) would therefore fail `shard_table_ok`. V7-R2 must use 48 × 2.25e6 (≈ 631 s each; 225 replicates per shard; pairs s ↔ s + 24; same 10 800 replicates, 5400 pairs and total histories), declared in the V7-R2 amendment before observation. Alternatively, a change of the 1.25 planning factor would have to be pre-registered; this plan does not propose one.
     - Frozen-design calibration shards ≈ 1080 s → ≈ 1805 s, **above 1500 s**. V7-R2 re-sizes its calibration into parts of ≤ 1100 s (repair plan §6.1–6.2: 78 parts × 84/83 repetitions ≈ 1095 s at 12.2 s per repetition for 10 800 replicates; 6 κ parts ≈ 645 s).
     - E1-B shards 543 s → ≈ 908 s and hr5c python shards ≈ 546 s → ≈ 913 s: they fit.
     - V2b-el shards 376 s → ≈ 628 s: fits.
   - Re-sizing only adds shards. History counts and criteria never decrease (plan, "Execution under the single-process directive").
3. **Logging.** Every local pytest and every sandbox compute job is logged under `/home/wahln/.local/share/ionmc-experiment/experiment-v3/telemetry/test-logs/<UTC>-<label>/{command.txt, meta.json, log.txt}`, unfiltered, with `python -X faulthandler -u -m pytest -v --capture=tee-sys --tb=long`. The full suite is run in two halves until a full run completes cleanly (crash #3 rule). Host steps carry `PYTHONFAULTHANDLER=1` (STEP_BASE_ENV).
4. **PARTIAL ledger.** Every producer step prints `PARTIAL <name>.json <sha256>`, **including the calibration shards**: the C8b omission must be fixed before the run, so the sealed-record deviation does not recur. The orchestrator copies each line from the protected stdout into the run ledger (`partials-<sha>.json`) immediately after each call, because a lost call (MCP abort) loses that stdout.
5. **Sealed records and attestation.** Consumers verify manifest digests before loading. `summarize.py --combine … --expected-sha <sha> --attest-sha <sha>` builds the archive. `record_local_validation` at the exact evidence SHA lists every host run id. A truthful head-check record at the final PR head references the evidence record and lists what was not re-run (memory note, V3-005A lesson).
6. **No direct executions.** With the C11 wrapper fixes, every step must run through its registration, so `complete: true` is achievable. Any direct execution is an incident and makes `complete: false`; it is recorded, not hidden.
7. **Expected flags of a clean future run:** `complete: true`, `reduced_history_counts: false`, `conformant: false` (single-process mode, imported partials, deferred step 86). `pass` follows the step verdicts. V6 is re-recorded as failed-not-sized by the deterministic `v6-unsized` steps, so under every option of §2 `pass` stays false at archive level, and the merge decision rests on the operator-confirmed (a)6 reading, never on the archive's `pass` field.

### 3.3 Lifting the deferral (operator only)

- **Deferred multiprocessing checks.** Running step 86 and any `cpu_workers > 1` row needs the operator to lift the single-process directive.
- **Imported-partial non-conformance.** This is structural (Amendment 5). It disappears only when all shards and their combine run in one invocation. That needs both the directive lifted and either a call budget above the 1800 s MCP idle limit or code-verifiable host records (`protected_host_records_verified_by_code`).
- **Proposal content (informational).** No lifting is requested now: the host stability problem that caused the directive (crashes #1–#6 and WSL crashes on 2026-10-11/12/13) has no recorded resolution. Task integration under `conformant: false` follows the V3-005A precedent. A release (V2-RELEASE) will need the deferral lifted, or an operator decision, at release time.

## 4. Sequence to an integrable state

All work stays on the V3-005B task branch under task id V3-005B (Amendment 14 (a)1; §5).

Compute jobs are serial:
- one sandbox or host job at a time;
- no sandbox compute while a host job runs;
- never two compute-heavy agents concurrently.

**No expensive qualification step runs before (iv) is answered and (v) is frozen.**

### 4.1 Phases and gates

**(i) Sandbox investigation and pre-registrations (host 0 h).**
- **Elastic: commit Amendment 19 first.** Commit the draft text of `validation/plans/v3-005d-elastic-attribution.md` revision 4, §6, before D-E0 runs.
  - Its seed list (§6 (f)) and closing line currently name lv 20511004 / hr 20521004. They must be changed before the commit to lv 20531004 / hr 20541004 (maximum 20562003), with 20511004 / 20521004 retired unused (§4.4).
  - Its 2099xxxx family (from 20990100) is disjoint from both bases.
- **Elastic: run D-E0 → D-E1 → D-E2 → D-E3** (attribution plan §3).
  - D-E3 is warp-cpu f64, 9 runs × 2 energies × 1e5 histories, with common-random-number seeds 20990100 / 20990101.
  - All of it uses the exploratory-table route without a `src/` change (plan §3.1).
- **V7-R2: commit the dated design amendment** (content of repair plan §3.7, corrected per V7-R review F1–F8):
  - V7-R2 is a new row;
  - its design;
  - the pre-declared size-selection rule of §2.4 (measured joint pass rate, repair plan revision 4);
  - seeds and transport layout.

  The sizes are frozen later, in a 14a-style addendum after the pool re-check of (vi), and never decrease.
- **V7-R2: joint-power measurement on the current pool (done).** 383/400 and 384/400 select 5400/0.005 (repair plan revision 4, §3.8). It informs the design only; it is not the check of record, which is the re-measurement on the rebuilt pool in (vi) step 2.
- Gate to (ii): D-E0–D-E3 reported (exploratory) and both amendments committed and reviewed.

**(ii) Host reference calculations (host ≤ 4.2 h plus a smoke case).**
- **D-E4.** Ten TOPAS/Geant4 thin-target cases, rows 14–23, through `run_reference_calculation` (seeds 20270000 + 1000·row + k).
  - Primaries: 2e7 with 2 mm targets at 30 / 50 MeV; 4e6 with 1 cm targets at 100–200 MeV.
  - A 1e5-primary smoke case comes first.
  - Each case runs ≤ 1500 s. If the smoke timing extrapolates above that, the case is split into k = 1, 2, as declared there.
  - Bound: 10 × 1500 s = 4.17 h.
  - D-E4 and the M2 input acquisitions run before the (iv) record is written into the main clone (§5 F5).
- **D-E5 and D-E6 (sandbox).**
- **Synthesis** under attribution plan §5, recorded as exploratory. V5, V11 and V10-A stay failed.
- Gate to (iii): "established" for H2 authorises the pre-committed M2 decision. M3 is allowed only on a pre-declared M2 failure. Other outcomes go to §6 S1.

**(iii) Model decision, implementation, CI, Codex (host: CUDA test calls only, each ≤ 1500 s; not budgeted).**
- **Elastic model.**
  - The M2 decision, with its inputs pinned.
  - Builder and table changes, with P6 (6), P6-D, X-MT5-SPEC and the flags re-checked.
  - P5-ext twins and C1-ext cases.
  - R1-D2 (C) directions re-declared; (A)/(B) bit-identical.
  - The `run_suite.py` 2099xxxx seed guard.
- **V7-R2** (repair plan §6.5):
  - a per-gate α parameter (default 0.04 keeps lv5b / dd976c5a bit-identical);
  - the `v7r-kappa` step and its sealed record;
  - calibration parts that consume it, each ≤ 1100 s;
  - the transport layout of §2.4;
  - a **V7-R2 reference index distinct from every shard index** (the constant `V7R_REF_INDEX = 16` collides with shard 16 of a 48-shard layout);
  - calibration parts that print `PARTIAL` lines (C8b omission).
- **Suite and seed-base plumbing** (review-route M2). Today these are hard-wired to lv5c/hr5c and 20481004/20491004. A decision on the suite identity is needed. Recommended: new suites `lv5d` / `hr5d`, following lv5 → lv5b → lv5c, so that the lv5c / hr5c archives keep re-verifying under their own definitions. Items:
  - qualification-base constants and `summarize.seed_blockers` (today it accepts only 20481004 for lv5c and 20491004 for hr5c; any other base becomes "non-qualification evidence");
  - `DEFAULT_SEED_BASES`;
  - the step list (`run_suite.full_step_names`: κ step, 78 calibration parts, 48 transport shards (64 under the 7200/0.005 fallback), reference k = 48);
  - `rehearsal_preflight` (keyed to suite `lv5c`);
  - the single-process workers = 1 force list (`run_suite` l.1318);
  - row mapping and the pending-rows accounting;
  - the cross-base reservation of hr r_index 6 (it equals the V7-R2 transport seeds, §4.4);
  - a **timing rule for rep-based steps.** `declared_plans()` covers only history-based plans, so the calibration and κ parts get either `declared_plans` entries (repetitions × measured seconds per repetition) or a stated separate rule. One short calibration part is measured on the slow host before the part sizes are frozen (V7-R review F8).
- **Option II only:** V6-F code (rows V6-F / V6-sens-F, seeds B + 7000 / 8000 / 9000 + k) enters the hashed source set. It lands in this phase, before the fixture rebuild (m7).
- Reviews and tests: Codex rounds until pass; CI green; the sandbox full suite in two logged halves.
- Gate to (v): Codex pass, CI green, and a complete suite record.

**(iv) ONE operator intervention (batched; category "ambiguous experimental boundary").** Raised once the attribution outcome is known, in parallel with (iii). It must be answered before (v) is frozen. Contents:
1. **(a)6 carrying.** Options I / II / III (§2.2). State plainly that I and II **reduce** the (a)6 gate set: V6 and V10-A are removed as gates, and V7-R2 replaces V7-R. They override the Amendment 18 (c) default "NOT merged". Unlike #81, the archive `pass` will be false. Recommendation I.
2. **V10-A.** 2A / 2B / 2C (§2.3). Recommendation 2B (V11 added to (a)6).
3. **V7-R slot.** 3A / 3B / 3C (§2.3). Recommendation 3A.
4. **V7-R2 calibration variant.** As written (≈ 33.9 h, or ≈ 44 h if the selection rule picks 7200/0.005; no clause change), or (B) (≈ 16.5 h), or (C) (≈ 12.5 h). Recommendation: as written.
5. **V5 fresh base (decision; review-route B1).**
   - Amendment 17 (b)4 reads verbatim: "Stop rule (Amendment 14 (f)): when `v5-attribution` detects the trigger it writes its document, exits non-zero and the remaining steps of the C9 order still run (they are independent evidence and are archived); no fresh base is taken for V5, no tolerance changes, and the orchestrator raises an intervention. 'The task stops' refers to fix-and-rerun, not to the archiving of the remaining rows."
   - The stop rule fired at dd976c5a.
   - Question: does disposition A authorise **one** fresh V5 base after the reviewed attribution and the M2 decision, with tolerances unchanged and the dd976c5a V5 / V11 failures kept?
   - **Default if refused:** V5 stays failed and is not re-run; V11 still runs as the model-comparison row. (a)6 cannot then pass on V5, and §6 S1 applies: carry V5 as failed and withdraw the elastic depth-dose claim, or keep PR #82 open.
6. **Bases.** Statement: lv 20531004 / hr 20541004 (V6-F 20551004 under Option II). The pre-declared 20511004 / 20521004 are retired unused, declared before observation, with no seed repetition. The seed-by-seed record (no shared seed for k ≤ 995) is attached (§4.4).
7. **Conformity.** Statement: single-process directive, `conformant: false`, imported partials, deferred step 86 (§3). The operator may lift the directive.
8. **Host budget** (§4.2): ≈ 41.3 h with V7-R2 at 5400/0.005, ≈ 51.5 h at 7200/0.005, ≈ 23.9 h with (B), ≈ 19.9 h with (C); Option II + ≈ 6.2 h CPU plus CUDA. This supersedes the ≈ 9.5 h quoted for V7-R in IR-DEFAF4.

The response is recorded verbatim and changes no row verdict.

**(v) Pre-registration amendment for re-qualification (dated, before any observation; independently reviewed).**
- **Suites and bases:** lv5d / hr5d (or the decided identity) at lv **20531004** / hr **20541004**. 20511004 / 20521004 are retired unused. Under Option II, V6-F runs at **20551004**. A committed test asserts the disjointness checks (a)–(e) and the reservations of §4.4.
- **Unchanged:** every frozen tolerance.
- **Rehearsal fixture (R2-6):** the amendment text declares the halved fixture, 4 × 2.25e6 histories at rehearsal seeds 20505000 + 16000 + k = 20521000–20521003 (k = 0..3). It is never evidence and is regenerated once at the evidence digest.
- **Rows:**
  - the §1.3 row set, with V6 / V6-sens re-recorded by `v6-unsized`;
  - V5 only if item 5 is granted;
  - V7-R2 (V7-R itself stays failed);
  - V11 into (a)6 if 2B;
  - V6-F if Option II.
- **Option II is binding once this amendment is frozen** (review-route M4). A V6-F failure then blocks integration. A later change of the gate set would be a new, post-hoc operator decision, never a pass of this route.
- **Seed checks:** every seed checked against the consumed list, the recorded rehearsal seeds, 2098xxxx / 2099xxxx, the TOPAS seeds, and one another. Known cross-base identities stay unused:
  - lv r_index 11–13 = hr r_index 1–3 (lv r_index 11 stays unused);
  - lv r_index 16 = hr r_index 6;
  - lv r_index 10 + k = hr r_index 0 + k (the hr base value itself; hr r_index 0 stays unused).
- **Crash rule:** a step that died without a published document is re-executed once at the same seeds, with both attempts recorded. A step with a published verdict is never re-run.

**(vi) Rehearsal and re-qualification run (≈ 37.2 h with V7-R2 at 5400/0.005, ≈ 47.3 h at 7200/0.005; ≈ 19.8 h with (B); ≈ 15.8 h with (C); Option II + ≈ 6.2 h).**
1. **Fixture rebuild at the evidence digest.**
   - The rehearsal digest is the hashed set minus `validation/plans/`, `decisions/`, `docs/`, README files and `tests/`. Only `src/` and `validation/scripts/` changes invalidate it.
   - Freeze both, including the Option II V6-F code and the M2 plumbing, and rebuild once.
   - Halved shards: 4 × 2.25e6 at rehearsal seeds 20505000 + 16000 + k = **20521000–20521003** (≈ 0.70 h total; m5).
2. **Pool re-check (sandbox, pre-declared; M3).** On the rebuilt pool, re-measure the complete two-gate joint pass rate at κ = 1 (R_meas = 400 per surrogate).
   - Cost: ≈ 2.2 h, and ≈ 3 h more if 7200/0.005 has to be measured.
   - Execution: single process, one heavy job at a time, never overlapping the D-E3 probes or any host run.

   Apply the pre-declared size-selection rule of §2.4: the **directly measured complete two-gate joint pass rate at κ = 1**: the repository `run_calibration`, with both gates on the same replicates, α_tost = 0.005, 48 × 225 replicates, R_meas = 400 per surrogate, seed roots 20989920–20989923 (planner block). Freeze 5400/0.005 (48 × 2.25e6) iff both surrogates reach **≥ 378/400** joint passes, i.e. the CP95 lower bound p_low satisfies P(Bin(3200, p_low) ≥ 2909) ≥ 0.99. Otherwise measure 7200/0.005 (64 × 2.25e6; roots + 10) and freeze it iff both reach ≥ 378/400. Otherwise no V7-R2 design is frozen and the finding goes to the operator (repair plan revision 4, §3.8 item 2). Then freeze the V7-R2 sizes addendum. Nothing decreases, and the rule is not changed after the check.
3. **Throughput step first** (review-route m2; 14a (c)5). `shard_table_ok` covers every history-based shard plan, and the rep-based rule of (iii) covers the κ and calibration parts. A non-fitting plan stops the run before the first evidence step.
4. **Evidence order** (one step per call; PARTIAL lines into the ledger immediately):
   1. deterministic steps (v6-e0, V10, V10-A, X-ENDF, D9-BGG);
   2. R1-D2, V8-LV, P7, v10-a-trigger;
   3. V11, then V5 (if granted) and its attribution step;
   4. V2b-el, V7-scan-el, V7-shift-el, `v6-unsized` / `v6-compare`;
   5. V6-F block (Option II);
   6. **V7-R2 block:** `v7r-kappa` → calibration parts → calibration combine (verdict before any V7-R2 transport; on a miss, V7-R2 fails and no transport is spent) → transport shards and reference → `v7r-combine`, diagnostic;
   7. E1-B, E-shape, the deferred partition step;
   8. hr: V8-stat (python × 4, cpu64, cuda32, cuda64, compare), V7 f32/f64 CUDA, `v6-unsized`, `v6-cuda-compare`, V6-F CUDA (Option II).
5. If the V5 stop rule fires again, the remaining steps still run and are archived (Amendment 17 (b)4).

**(vii) Closure and merge conditions.**
- Records: archives, README, ledger and the decision 0041 outcome; `record_local_validation` at the evidence SHA; a head-check at the PR head; Codex; CI.
- **Merge iff** every entry of (a)6, as carried by the operator's answer to (iv), passes at the evidence SHA.
- Every recorded failure is listed next to the passes. No implementation-only passing record.

### 4.2 Host-hour budget (computed; s = 1262/755 = 1.672 from the identical V7-R rehearsal shard, §3.2)

| Item | Basis (source) | Host hours at s |
|---|---|---|
| D-E4: 10 TOPAS cases + smoke | 10 × ≤ 1500 s (attribution plan §3) | ≤ 4.17 + smoke |
| lv/hr suite without V7-R | 14.45 h × s = 24.16 h, minus V7-R transport 5.69 h and calibration 15.19 h | 3.28 (includes orchestration gaps) |
| **V7-R2 block as written (5400/0.005)** | repair plan rev 3 §6.2 (fixture 0.7, κ 1.1, calibration 23.2, A13 0.2, transport 8.7) | **33.9** |
| V7-R2 7200/0.005 (selection-rule fallback) | repair plan rev 3 §6.2 | 44 |
| (B) / (C) (operator) | repair plan rev 3 §6.2 | 16.5 (22 at R = 1600) / 12.5 |
| V6-F (Option II) | 72 × 184 s × s (v6f-followup §8.3), CPU | 6.15 + CUDA |

**Re-qualification:** 37.2 h (V7-R2 at 5400/0.005), 47.3 h (7200/0.005), 19.8 h ((B)), 15.8 h ((C)).

**Total host time to an integrable state:**

| Path | Host hours |
|---|---|
| Option I + V7-R2 (5400/0.005) | ≈ **41.3** |
| V7-R2 at 7200/0.005 | ≈ 51.5 |
| with (B) | ≈ 23.9 |
| with (C) | ≈ 19.9 |
| Option II | + ≈ 6.2 (+ CUDA) |

Revision 6: the V7-R2 joint-power re-measurement on the rebuilt pool is sandbox time (≈ 2.2 h; ≈ 3 h more for 7200/0.005). It adds no host hours, so these totals are unchanged.

If item 5 is refused, the V5 steps drop out: four short V5 shards of 10–25 s each at the dd976c5a rate, which is negligible.

### 4.3 Critical path

1. Amendment 19 and the V7-R2 design amendment.
2. D-E0 → D-E3 (sandbox).
3. D-E4 (host ≤ 4.2 h).
4. D-E5, D-E6 and the synthesis.
5. M2 decision → implementation (elastic, V7-R2, suite plumbing, Option II code) → Codex / CI. This is the longest calendar item, with no host hours. Intervention (iv) runs in parallel and must be answered here.
6. Re-qualification amendment (v) and its review.
7. Fixture (0.7 h host), then the sandbox joint-power re-measurement on the rebuilt pool (≈ 2.2 h sandbox; ≈ 3 h more if 7200/0.005 is needed). It runs single-process, as the only heavy job: no host run and no D-E3 probe at the same time.
8. Throughput, then the re-qualification ≈ 37.2 h serial host time (≈ 47.3 h at 7200/0.005) (about 1.5 days of back-to-back calls). The V7-R2 calibration block (≈ 24.5 h) is its largest item and cannot start before the model and fixture are final.
9. Closure.

(C) would shorten the run to ≈ 15.8 h; that is the operator's lever. The revision-6 re-measurement adds only sandbox time, so the host totals of §4.2 are unchanged.

### 4.4 Fresh bases (revision 4: 20531004 / 20541004 adopted; seed-by-seed record kept)

**Adopted allocation** (declared in the (v) amendment before any observation; consistent with the V7-R repair plan revision 2, §5 there):

| Base | Use |
|---|---|
| lv **20531004** | elastic re-qualification rows and V7-R2 |
| hr **20541004** | – |
| **20551004** | V6-F (Option II only) |
| 20561004, … | next; the sequence skips 2098xxxx / 2099xxxx |

**Why these bases.** Amendment 18 (b)5 / v6f-followup §7 (b): "B + 7000 must exceed the largest rehearsal seed actually used (20525999) or a seed-by-seed disjointness proof is recorded". Every seed of 20531004 and 20541004 (base + 1000·r + k ≥ 20531004) exceeds 20525999, so the simple range rule holds without a proof.

**The skipped bases.** The pre-declared Amendment 14 (j) bases **20511004 / 20521004 are retired unused**. They are declared skipped before any observation; no seed of them is ever drawn, so there is no seed repetition.

**Seed-by-seed record (kept as the reason the skipped bases would also have been admissible with a proof):**
- Method: `validation/plans/v3-005d-integration-route/scripts/used_seeds.py` scans every JSON document and seed-bearing text line under `validation/generated/RUN-*` (337 runs, read-only). `validation/plans/v3-005d-integration-route/scripts/seed_check.py` compares the result with `base + 1000·r_index + k` for r_index 0–20 and k < 1000. Output: `validation/plans/v3-005d-integration-route/scripts/seed_check.out`.
- Recorded rehearsal transport seeds:

  | Seed | Use |
  |---|---|
  | 20506000 | throughput and `v8-stat-python-s0` rehearsals |
  | 20517000 | `v6-s0` |
  | 20521000, 20521001 | V7-R fixture shards |
  | 20525000 | V6-sens nominal and st-hi |

  Added from the plan: E-shape 20509500 + k and the planned halved fixture seeds 20521000–20521003.
- **20511004 / 20521004 share none of them for any k ≤ 995.** The first shared values need k = 996–999.
- Revision 2's "collision" claim is withdrawn. Those bases fail only the range rule above. The V7-R repair plan revision 2 still says they "collide" (§5 there); that wording should read "fail the range rule".
- The same scan finds nothing for 20531004 / 20541004 at any k < 1000.

**Seeds of record:**

| Item | Seeds | Source |
|---|---|---|
| V7-R2 transport (lv r_index 16) | **20547004 … 20547051** (k = 0..47); reference k = 48 → **20547052**. 7200/0.005 fallback: k = 0..63, reference k = 64 → **20547004–20547068** | repair plan rev 3 §5 |
| Reserved hr r_index 6 at 20541004 | numerically equal to the V7-R2 seeds (k ≤ 64); must stay unused | repair plan rev 3 §5 |
| Other cross-base identities (lv 20531004 r ≡ hr 20541004 (r − 10) at equal k) | lv r_index 10 ≡ hr r_index 0, unused; lv 11–13 ≡ hr 1–3, as in Amendment 14 (j), lv 11 unused | – |
| V6-F (Option II, B = 20551004) | r_index 7 / 8 / 9: 20558004 … 20558019, 20559004 … 20559011, 20560004 … 20560019. B + 7000 = 20558004 > 20525999. Above every lv 20531004 seed (≤ 20552003 for r_index ≤ 20). Numerically equal only to hr 20541004 r_index 17–19, which no hr row uses; reserved | v6f-followup §7 |
| V7-R2 κ and calibration | roots 20980000 (κ) / 20982000 (calibration) / 20984000 (A13 shard), within **20980000–20984999**. Retry block 20985000–20988999 under the F6 conditions; planner diagnostics 20989900–20989999 | repair plan rev 3 §5 |
| Elastic investigation (Amendment 19) | family **2099xxxx from 20990100**: D-E3 / D-E6 20990100 / 20990101; D-E5 20991000 + case | elastic plan rev 4 §6 (f) |
| TOPAS | 20270000 + 1000·row + k, rows 14–23 | elastic plan rev 4 §6 (f) |

**Required edit to the elastic plan (before Amendment 19 is committed).** Its disjointness list (§6 (f): "the fresh bases 20511004 / 20521004 and their r_index blocks (maximum 20542003)") and its closing line ("… at the next declared fresh bases (lv 20511004, hr 20521004) …") must name **20531004 / 20541004** (maximum 20562003) and state that 20511004 / 20521004 are retired unused. The 2099xxxx family is disjoint from both pairs either way.

## 5. Open infrastructure items

| # | Item | Consequence for the route | Source |
|---|---|---|---|
| F1 | **The PR body of #82 cannot be updated.** `create_task_pull_request` refuses an existing PR, there is no update tool, and direct use of the GitHub CLI is forbidden. | The results README section and decision 0041 are the description of record. The closure (vii) points there. | IR-DEFAF4 `reason_autonomous_resolution_not_possible` |
| F2 | **`create_task_worktree` branches only from `v3/develop`.** | All V3-005D work stays on the V3-005B branch under task id V3-005B. "V3-005D" names a scope only, as "V3-005C" did (Amendment 14 (a)1). If `v3/develop` moves, `start_codex_review` fails until develop is merged into the branch (`git -c user.name=… -c user.email=… merge --no-commit --no-ff v3/develop`, tests, then `commit_task_changes`). A merge that touches hashed files invalidates the rehearsal fixture. | memory `v3-environment-constraints.md` |
| F3 | **CI has no data cache.** Data-dependent tests skip in CI and fail under `IONMC_REQUIRE_DATA=1`. The 36 slice-C tests lacked the guard and failed CI at 4ce33832 until fixed. | Green CI is not evidence for data-dependent rows. Every new test that opens a nuclear or elastic table needs the guard and must run locally in required mode. The local logged suite is the record. | Amendment 18 (e) |
| F4 | **WSL and host instability.** Crashes during full pytest runs (2026-10-11, 2026-10-12) and crash #3 on 2026-10-13 (hypervisor CPU power states; the operator disabled CPU idle states, after which the host is ≈ 40 % slower in throughput). Earlier: six crashes up to 2026-10-07, which led to the single-process directive. | The two-half pytest protocol continues until a full run completes cleanly. After any crash: one instrumented re-run, and after a second crash, narrow before any heavy run. A host step that died without publishing a document may be re-executed once at the same seeds (precedent 74 → 74b), with both attempts recorded; this is declared in the (v) amendment. A step with a published verdict is never re-run. | README "Incidents"; memory |
| F5 | **The v3-research MCP requires a clean main clone** (`/home/wahln/aiprojects/ion-mc`). `run_reference_calculation`, `acquire_reference_data` and the related tools fail with a bare error on any untracked file there. The intervention tools write IR records into that clone. `materialize_reference_data` also fails on a dirty task worktree. | Order: D-E4 and any new acquisitions (for example the M2 inputs: the SAID NN solution, the Horiuchi et al. 2007 β_NN and the de Vries 1987 densities, attribution plan §4.1) run **before** the (iv) intervention is raised, or the IR record is committed on the task branch and its untracked copy is removed from the main clone (backup in scratchpad) before the next research call. Materialise data before workers dirty the worktree. | memory |
| F6 | **The MCP idle abort is at 1800 s.** A silent host call loses its protected stdout (the `PARTIAL` lines). | Steps ≤ 1500 s hard and ≤ 1200 s planned; ledger every `PARTIAL` line immediately (§3.2). | Amendment 14 (k); memory |
| F7 | **Codex review and host runs need a frozen worktree.** "Worktree changed during review" aborts a review, and host runs require a clean tree. | Sequence: worker idle → commit → review and host runs → fixes. No edits during the re-qualification. | memory |
| F8 | **Merge tooling.** `merge_task_pull_request` needs a passing `record_local_validation` at the exact PR head. | A truthful head-check record at the final head references the evidence-SHA record (§3.2.5). | memory (V3-005A lesson) |

## 6. Risks and stop rules

Each stop below means: record the failure, archive the remaining steps, and take no further base without a new reviewed physics reason (Amendment 14 (f)). What goes to the operator is a **proposal**, never a relaxation.

| # | Result | Why integration becomes impossible without a scope change | What goes to the operator |
|---|---|---|---|
| S1 | D-E5 gives **"unexplained"** or only **"partially attributed"**, **or the operator refuses a fresh V5 base** ((iv) item 5; Amendment 17 (b)4) | No reviewed physics reason or no V5 base, so V5 / V11 stay failed under every (a)6 reading | A proposal to carry V5 / V11 as recorded failures and withdraw the elastic depth-dose claim (elastic channel shipped as an option without a validated depth-dose claim, or off by default), or keep PR #82 open |
| S2 | Attribution established, M2 implemented, V5 granted, but **V5 or V11 fails again** at lv 20531004 | The granted base is used. A second fix needs a new reviewed physics reason and a new operator decision under 17 (b)4; M3 only on its pre-declared triggers | As S1, with both failures listed |
| S3 | The **V5 stop rule fires again** (peak/plateau low, or > 5 %) | "The task stops" (fix-and-rerun) per Amendment 14 (f) and 17 (b)4 | Escalation with the attribution record |
| S4 | **V7-R2 fails.** The calibration block misses a clause at the evidence SHA (repair plan risk: the modelled single-gate box power at α = 0.005), or a gate fails on the new data (a correct estimator at true paired coverage 0.69–0.695 passes with only 0.86 / 0.61) | V7-R2 is a recorded failure with no further count increase and no further base (repair plan §3.7 item 5). A new design would need a new reviewed reason | Option 3B (coverage claim withheld) as the fallback. Any clause change (for example (B)/(C) after the fact) is excluded post hoc |
| S5 | D-E4 shows a **σ_el rate difference** (H1: Geant4/ionmc rate ratio outside 1 ± 0.03) | The Barashenkov transcription would be an implementation defect (rung 0). The fix is autonomous, but V10-A (c) and the "shared lineage" statements change | Reported in (iv). The (c) verdict at dd976c5a stays recorded |
| S6 | **V8-nuc inconclusive** or failing at the hr re-qualification base | Under its unchanged rule an inconclusive verdict is not a pass and no larger sample follows (Amendment 14 (d)) | A proposal on how to carry it |
| S7 | Option II only: **V6-F fails**, or its coverage is "not established (beam-sensitive)" | Option II is binding once (v) is frozen; one execution only (v6f-followup §8.2). Integration is blocked | A recorded failure. There is no fall-back to Option I; any change would be a new, post-hoc operator decision, never a pass of this route |
| S8 | The committed seed test finds a shared seed with a newly recorded run (none at present: §4.4) | – (not a scope change) | Recorded. The next base in the declared allocation (20561004, …, skipping 2098xxxx / 2099xxxx) is declared before observation |
| S9 | **Host slowdown worsens**, or the throughput step reports a shard plan > 1200 s (V7-R2 parts are planned at ≤ 1100 s, transport at ≈ 631 s) | – | Re-shard (counts unchanged) before the first evidence step. The budget in (iv) is updated |
| S10 | **V7-scan-el fails again** (not in (a)6) | – | The V2-NUM / V1-MUST-030 claim for `sec_p_dose` slope is narrowed in the ledger. No gate change |

The planner sees no route in which (a)6 **as written** is met. Integration of V3-005 therefore depends on the operator's answer to §4 (iv), item 1, in every scenario.

## Appendix A. Computations used (planner arithmetic; reproducible)

```
s = 1262/755 = 1.672           # identical V7-R rehearsal shard: 752–771 s (≤ dd976c5a) vs ≈ 1262 s at a6ca56a4
V7-R shards: (16*755 + 166) s = 3.40 h -> 5.69 h; one shard 755 s -> 1262 s (> 1200 s effective limit)
Calibration: (30*1080 + 308) s = 9.09 h -> 15.19 h; one shard 1080 s -> 1805 s (> 1500 s)
  parts per surrogate p: t = (3200/p * 4.86 + 42) * s -> p=15: 1803 s, 20: 1370 s, 25: 1110 s, 30: 937 s
V2b-el 2*376 s -> 0.35 h; E1-B 3*543 s -> 0.76 h (908 s per shard); hr5c python shard 546 s -> 913 s
hr5c compute (546+541+540+539+173+22+26) s = 0.66 h -> 1.11 h
lv5c span 2026-10-12T14:23:31 -> 2026-10-13T04:09:17 = 13.76 h; hr5c 04:20:34 -> 05:02:02 = 0.69 h
  sum 14.45 h -> 24.16 h
Fixture 2*1262 s = 0.70 h; V6-F 72*184 s = 3.68 h -> 6.15 h; (rev 1 D-E4 4*1200 s -> 2.23 h, superseded)
Revision 2:
  non-V7-R remainder of lv/hr = 24.16 - 5.69 - 15.19 = 3.28 h
  (rev 2 figures, superseded by revision 5) re-qualification = 3.28 + 33.7 (V7-R2) = 36.98 h; + 16.3 (B) = 19.58 h; + 22.0 (B, R 1600) = 25.28 h; + 12.3 (C) = 15.58 h
  D-E4 upper bound 10 * 1500 s = 4.17 h; totals 4.17 + 36.98 = 41.15 h; 4.17 + 19.58 = 23.75 h; 4.17 + 15.58 = 19.75 h
  V7-R2 transport shard 4.5e6: 755 s * 1.672 = 1262 s > 1200 s; 2.25e6: 631 s
  lv 20511004 + 1000 r for r = 10..14: 20521004, 20522004, 20523004, 20524004, 20525004 (all <= 20525999)
  hr 20521004 (r = 1..3: 20522004 .. 20524004) inside 20505000..20525999
  (rev 2, superseded) lv 20531004 + 16000 = 20547004; V6-F 20551004
Revision 3: seed-by-seed scan validation/plans/v3-005d-integration-route/scripts/used_seeds.py -> used {20506000, 20517000, 20521000, 20521001, 20525000};
  validation/plans/v3-005d-integration-route/scripts/seed_check.py -> lv 20511004 / hr 20521004 share none for k <= 995 (first at k = 996..999);
  (rev 3, superseded) V7-R2 transport at 20511004 would have been 20527004 .. 20527052;
  (rev 3) V6-F B = 20531004: r7-9 seeds 20538004 .. 20540019
Revision 5: re-qualification 3.28 + 33.9 = 37.18 h; 3.28 + 44 = 47.28 h; 3.28 + 16.5 = 19.78 h; 3.28 + 12.5 = 15.78 h;
  with D-E4 (4.17 h): 41.35 / 51.45 / 23.95 / 19.95 h; 7200 transport seeds 20547004 + k, k = 0..64 -> 20547068
Revision 4: adopted lv 20531004 / hr 20541004: every seed >= 20531004 > 20525999; scan finds no shared seed at any k < 1000;
  V7-R2 transport 20531004 + 16000 + k = 20547004 .. 20547052 (k <= 48) = hr 20541004 r_index 6 numbers (reserved);
  V6-F B = 20551004: r7 20558004..20558019, r8 20559004..20559011, r9 20560004..20560019; B + 7000 = 20558004 > 20525999;
  lv 20531004 max seed for r_index <= 20, k <= 999: 20552003 < 20558004
  calibration parts 39 per surrogate * (82 or 83) = 3200; 78 parts
Seed algebra: 20511004 + 1000 r + k = 20505000 + 1000 (r+6) + (k+4);
              20521004 + 1000 r + k = 20505000 + 1000 (r+16) + (k+4)
```

Assumptions, stated: the factor s is assumed constant across workloads. It was measured on one workload, warp-cpu V7-R shards, and the throughput step at the evidence head re-measures every rate. The lv5c span includes orchestration gaps, so scaling it by s is conservative for compute and ignores gap changes. The CUDA rates after the slowdown are not measured.

## Appendix B. Files read for this document

- Worktree root: `REQUIREMENTS.md`, `EXPERIMENT.md` (Permitted Human Intervention), `AGENTS.md`; `experiment/v3/REQUIREMENTS.md` (V2-CAP, V2-NUM, V2-EVID, V2-RELEASE); `experiment/v3/requirements-index.json`.
- `validation/plans/v3-005-acceptance.md`: criteria, Seeds, single-process and deferred sections, evidence archive table, Amendments 5, 14 (a)–(m), 14 addendum, 15, 15 addenda, 16, 17, 14a, 18.
- `validation/plans/v3-005c-v6f-followup.md` §§1, 7, 8, 10.
- `validation/results/transport/README.md`: "V3-005C results", Incidents, Summarizer corrections.
- `validation/results/transport/{lv5c,hr5c}-dd976c5-single-process.json`: summarised with `json`.
- `validation/requirement-ledger.md`.
- `decisions/0041-proton-nuclear-interactions.md`: limitations, V3-005C outcome.
- `reference_requests/interventions/IR-20261013-164616-DEFAF4/{request,response}.json`.
- `tests/ionmc/fixtures/v7r/rehearsal/PROVENANCE.json`.
- Reviews: `validation/plans/v3-005d-integration-route/review.md` (revision-2 review), `validation/plans/v3-005d-v7r-calibration-repair/review.md`. Scripts: `validation/plans/v3-005d-integration-route/scripts/used_seeds.py`, `validation/plans/v3-005d-integration-route/scripts/seed_check.py` (+ `.out`), which read `validation/generated/RUN-*` read-only.
- Sibling plans: `validation/plans/v3-005d-elastic-attribution.md` revision 4 (with `validation/plans/v3-005d-elastic-attribution/review.md`, round-3 verdict ACCEPT conditional on R3-1/R3-2, applied in revision 4) and `validation/plans/v3-005d-v7r-calibration-repair.md` (§§1–6). In revision 1 the V7-R plan did not yet exist.
- Memory notes `v3-release-plan-schema.md`, `v3-environment-constraints.md`, `v3-operator-directives.md`, `v3-progress-log.md` (timing of the a6ca56a4 fixture shards).
