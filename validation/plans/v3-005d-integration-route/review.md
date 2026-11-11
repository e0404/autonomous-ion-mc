# Independent review: V3-005 route to integration (integration-route.md, revision 2)

Reviewer: independent sub-agent, 2026-10-14. Repository: V3-005B worktree at HEAD 749f768c (read only).
Scope: requirement coverage, V6/V6-F treatment, execution conformity, sequence and budget, seeds,
infrastructure and stop rules, and the classification of V7-R2.

## Checks performed (all verified, nothing re-derived from the route itself)

- The archives lv5c/hr5c-dd976c5 were summarised with python json:
  - lv5c: pass/conformant/complete = false/false/false; failed steps 03,16,21,22,26,28,29,78,79,85; imported-partial steps 11,16,21,22,25,79,80,84; deferred 86; problems = missing JSON in 78 and 85.
  - hr5c: false/false/true; imported 08; no deferred step.
  - Row verdicts and step numbers in route §1.3 and §3.1 match. The V7-R row has gates_pass true and calibration_pass false.
- Spans from the run ids: lv5c 14:23:31 to 04:07:58 (last start) = 13.74 h, hr5c 0.69 h. These are consistent with the route's 13.76 h (end time) and 0.69 h.
  - Per-step gaps: V7-R transport about 763 s, calibration about 1088 s, s30 about 318 s, 73 to 74b 1753 s (the anomaly). They agree with Appendix A.
- Arithmetic reproduced:
  - s = 1262/755 = 1.6715; 14.45 h × s = 24.16 h; 3.40 h × s = 5.69 h; 9.09 h × s = 15.19 h; remainder 3.28 h.
  - 3.28 + 33.7 = 36.98 h; with (B) 19.58 h; with (C) 15.58 h. With D-E4: 41.15 / 23.75 / 19.75 h.
  - V6-F: 72 × 184 × s = 6.15 h. V7-R2 calibration: 6400 reps × 12.2 s = 21.7 h, plus 78 × 70 s = 1.5 h.
- 1200 s limit: `steps_v5d.plan_table` sets `shard_plan_c(..., HOST_STEP_WALL_LIMIT_S=1500, v5.MARGIN=1.25)`, and `step_throughput` sets `shard_table_ok`. 14a (c)5 makes a plan that does not fit a stop.
  - The plan time is n / block rate. 14a (c)3 recorded 813 s planned against 758 s measured, so the block-rate plan is more conservative than wall time. The 4.5e6-shard failure claim therefore holds a fortiori (about 813 × 1.67 ≈ 1359 s planned).
- Seeds:
  - Amendment 14 (j): base + 1000·r + k; r_index table; declared lv 20511004 / hr 20521004, then 20531004 …
  - v6f-followup §7 (b): B + 7000 > 20525999, or a seed-by-seed proof.
  - Code: `v7r.V7R_R_INDEX = 16`, `V7R_REF_INDEX = 16` (a constant). The `v7r_rehearsal.py` docstring says the rehearsal seeds coincide with 20521004 "only as a number, never as a used seed".
  - Route algebra checked: 20511004 + 1000r + k = 20505000 + 1000(r+6) + (k+4), and 20521004 + 1000r + k = 20505000 + 1000(r+16) + (k+4). 20531004 + 16000 = 20547004. 20551004 + 7000 = 20558004.
- Code paths:
  - `summarize.verify/combine`: conformant requires pass; complete = no missing steps and no problems.
  - `seed_blockers` hard-codes 20481004 / 20491004 for lv5c/hr5c.
  - `rehearsal_preflight` applies only to suite == "lv5c".
  - `run_suite` l.1318 forces workers = 1 only for the listed suites.
  - `SINGLE_PROCESS_ENV` and `STEP_BASE_ENV` are as the route describes.
- Precedent: lv5-b84fdf3 has pass true and conformant false (integrated in #81).
- The numbers in §1.2 and §1.3 trace to the README V3-005C section, Amendments 14a and 18, IR-DEFAF4 request.json (CI 38431680322), and the elastic plan §0 / table l.94 (−3.66 / −4.89 / +10.35, from the step-21 document; the README rounds them to −3.7 / −4.9 / +10.3).

## Findings

### B1 (blocking): V5 fresh-base re-run is presented as "statement, no decision requested", but a frozen rule forbids it
- **Evidence.**
  - Amendment 17 (b)4: when the V5 stop rule fires, "**no fresh base is taken for V5**, no tolerance changes, and the orchestrator raises an intervention".
  - Amendment 14 (f): "the task stops" (fix-and-rerun).
  - The stop rule fired at dd976c5a.
  - Route §1.2 (V5 row), §4.1 (iv) item 5 and (v) nonetheless schedule V5 at a fresh base pair as a statement only.
  - IR-DEFAF4 ("prioritize elastic-model attribution … before another full qualification run") implies a further qualification. It does not explicitly lift the 17 (b)4 prohibition for V5.
- **Why this blocks.** This is the most consequential row of the re-qualification. Taking a fresh V5 base contrary to a frozen sentence is a change to the required validation rules, which the operator asked to have brought back.
- **Fix.** Make item 5 a decision item in (iv). Quote 17 (b)4 verbatim, and ask whether disposition A authorises one fresh V5 base after a reviewed attribution. Default if refused: V5 stays failed, and S1/S3 handling applies. Remove "no decision is requested" in §1.2.

### M1 (major): Option I is a reduction of the (a)6 gate set, and the route does not say so; it also misstates 18 (c)
- **Evidence.**
  - Option I merges while V6, V10-A and V7-R (the slot) are recorded as failed scientific gates.
  - Amendment 18 (c): "V3-005C therefore cannot be merged on V6 under (a)6".
  - The recommendation says Option I "honours Amendment 18 (b)/(c) … literally". It honours 18 (b) (no V6-F substitution). It overrides 18 (c)'s merge default by amending (a)6.
  - "None of the options … relaxes a tolerance" is true, but it hides that Option I relaxes the merge gate itself: two gating rows become carried failures.
  - The V3-005A precedent (lv5 archive pass = true) does not cover merging with archive pass = false.
- **Fix.** State plainly that Options I/II amend (a)6 by removing V6 and V10-A as gates (and substituting V7-R2) and are therefore scope reductions only the operator can make. Replace "honours 18 (c) literally" with "requires the operator to amend (a)6, which 18 (c) leaves to the operator". Note that, unlike #81, the archive pass flag will be false.

### M2 (major): The implementation list omits the suite and seed-base plumbing; the expected flags and the throughput claim are wrong without it
- **Evidence.**
  - `summarize.seed_blockers` accepts only 20481004 (lv5c) and 20491004 (hr5c). Any other base adds "not the qualification base … non-qualification evidence".
  - `run_suite.full_step_names` fixes the lv5c step list (16 V7-R shards, 31 calibration shards). V7-R2 changes it (κ step, 78 parts, 48 shards, reference index), and old-archive re-verification depends on it.
  - `rehearsal_preflight` and the workers = 1 force (l.1318) are keyed to suite names.
  - `V7R_REF_INDEX = 16` would collide with shard 16 of a 48-shard layout.
  - `declared_plans()` has no entry for calibration or κ parts: they are rep-based, not history-based. So "shard_table_ok … including the V7-R2 transport and calibration parts" (§4.1 (vi)) is not true of the code.
  - §3.2 item 1 uses `--suite <lv|hr>`. Those names are V3-003D suites; this task uses lv5c/hr5c.
- **Fix.** Add to §4.1 (iii) a decision on the suite identity (e.g. new suites lv5d/hr5d, as lv5 → lv5b → lv5c did), with:
  - qualification base constants;
  - seed_blockers;
  - DEFAULT_SEED_BASES;
  - the single-process force list;
  - the preflight;
  - row mapping;
  - V7R2 reference index;
  - a timing entry for calibration and κ parts in the throughput table, or a stated separate rule.

  Correct the CLI names.

### M3 (major): The V7-R2 surrogate pool is stated both as fixed and as rebuilt
- **Evidence.**
  - The repair plan (§3.4, §5) holds the rehearsal pool fixed (sha 5d28a9b3…/d4f7d749…, "identical at dd976c5a and at HEAD").
  - Route §4.1 (vi) and §4.3 say the calibration consumes the fixture rebuilt after the M2 model change. This is correct under 14a (b)2 and 17 (a)7, and `rehearsal_preflight` enforces it.
  - After M2 the block sums, and so the pool hashes, change. The repair plan's OC numbers (power 0.954, P(all) ≈ 1.000) were computed on the old pool.
  - The route quotes them without this caveat. §2.4 lists "the surrogates" as unchanged.
- **Fix.**
  - State that the surrogate *definition* is unchanged and the pool *content* is regenerated at the evidence digest.
  - Require the repair plan to drop the pinned pool hashes for evidence use.
  - Add a pre-evidence sandbox re-check of the OC/power on the rebuilt pool. This is a pre-declared smoke, not a gate change.
  - If the fixture shards are halved, declare the k range (seeds 20521000–20521003).

### M4 (major): The S7 fallback lets a gate be dropped after it fails
- **Evidence.** Under Option II, V6-F is a gate. S7 says that if V6-F fails, "fall back to Option I treatment by operator decision". That is choosing the gate set after seeing the result. It conflicts with C590CF (pre-register; no reinterpretation to obtain a pass) and with "V6-F passing does not authorize … merging".
- **Fix.** State that the Option II choice is binding once (v) is frozen. A V6-F failure blocks integration under Option II. Any later change is a new operator decision, recorded as post hoc and never presented as a pass of the Option II route. Add this consequence to Option II's risk list so that the operator chooses knowingly.

### M5 (major): The collision claim is overstated; seed-by-seed evidence is missing
- **Evidence.**
  - Route header item 3 and §4.4 say 20511004 / 20521004 "overlap the rehearsal seeds actually used" and that an exact collision exists "whenever" r′ = r + 6, k′ = k + 4. The route never shows a used rehearsal seed with k′ ≥ 4 in r′ = 16–20.
  - The known used seeds (20517000, 20521000/1, 20525000) have k′ ≤ 1, and `v7r_rehearsal.py` says the coincidence is "only as a number".
  - What actually fails is the range criterion (v6f-followup §7 (b): B + 7000 > largest used, or a seed-by-seed proof).
  - "Largest used seed 20525999" is a block bound; the largest listed used seed is about 20525000.
- **Fix.** Either:
  - perform the seed-by-seed check from the C8 rehearsal run records (cheap) and report actual intersections; or
  - reword to "fails the range-disjointness criterion; skipped conservatively". In that case keep the skip as an amendment item that changes a pre-declared Amendment 14 (j) base: bookkeeping, but declared before observation, with the reason stated accurately.

### m1 (minor): V7-R2 is called "autonomous to run"
- **Evidence.** §2.2 and §2.4 say "designing and running V7-R2 is autonomous". The repair plan's own table and the V7-R plan review say host execution (≈ 33.7 h) is scheduled under disposition A after the reviewed route, and that the budget is an operator item.
- **Fix.** Say: "design, pre-registration, implementation and sandbox work are autonomous; host execution occurs only in (vi), after (iv) and (v), with its budget in (iv) item 8."
- Otherwise item 6 is consistent with 14 (h)4 and 14a (b)6. V7-R stays failed; V7-R2 is a new row like V7 → V7-R; the slot substitution is an operator item.

### m2 (minor): The throughput gate is ordered after evidence steps
- **Evidence.** 14a (c)5 says a non-fitting plan "stops C9 before the first evidence step". The route's own evidence order (vi) runs V10-A, X-ENDF and the other deterministic steps before throughput, as at dd976c5a (step 07).
- **Fix.** Run throughput first, or treat the rehearsal throughput at the evidence digest as the binding gate and say so.

### m3 (minor): §3.1 understates why `conformant` is false
- **Evidence.** Code: conformant requires pass, so failed steps alone make it false. Missing documents (problems) also make `pass` false in `combine`.
- **Fix.** Add both to the table.

### m4 (minor): S8 reuses the V6-F candidate base
- **Evidence.** S8 names "20551004 / 20561004" as the next lv/hr pair, but 20551004 is the Option II V6-F candidate.
- **Fix.** Declare one ordered allocation (e.g. V6-F 20551004, then lv/hr 20561004 / 20571004).

### m5 (minor): Rehearsal fixture shards are over the planning target
- **Evidence.** The fixture shards (2 × 4.5e6, ≈ 1262 s) exceed the 1200 s / 1100 s targets, though not the 1500 s hard limit.
- **Fix.** Halve them by default (as the route allows) and declare the seeds.

### m6 (minor): Calibration repetition count does not add up
- **Evidence.** "78 parts × 84/83 repetitions" gives about 6513 repetitions, not 6400. 39 parts per surrogate × 82/83 = 3200.
- **Fix.** Correct this here and in the repair plan; the cost is unchanged (21.7 h).

### m7 (minor): V6-F under Option II enters the hashed source set
- **Evidence.** 18 (b)5 kept V6-F code out of the lv5c hashed set. Under Option II it enters the set and so affects the fixture.
- **Fix.** Note the fixture/digest consequence and the sequencing.

## Items confirmed correct
- §1.3 statuses and pointers.
- §3.1 flags.
- 1200 s limit and the need to halve the 4.5e6 transport shards (48 × 2.25e6, ≈ 631 s).
- Budget arithmetic.
- 18 (b) handling: no V6-F substitution.
- Phase gating: no expensive qualification step before (iv) and (v). D-E4 is investigation under disposition A.
- F1–F8 are accurate, and the requirement IDs are complete against the ledger rows that reference V3-005.

## Verdict
**REVISE.** One blocking finding (B1) and five major findings (M1–M5). After B1 and M1–M5 are fixed, I expect no blocking issue to remain. The route's structure, arithmetic and phase gating are sound.

---

# Round 2 (2026-10-14): integration-route.md revision 4

Cross-checked against `v7r-calibration-repair-plan.md` **revision 3** and `elastic-attribution-plan.md` revision 4. Also checked: `route-scripts/used_seeds.py`, `seed_check.py`, and the round-1 code facts.

## Round-1 findings: status

| Finding | Status | Evidence |
|---|---|---|
| **B1** | **Fixed** | (iv) item 5 quotes 17 (b)4 verbatim, asks the question and gives the refusal default. V5 runs in (v)/(vi) only "if granted". §1.2 row, S1 and S2 are updated. |
| **M1** | **Fixed** | §2.2 states that Options I/II reduce the (a)6 gate set and override the 18 (c) default (operator only). It notes the archive will have `pass: false`, unlike #81. The recommendation no longer claims to honour 18 (c) literally. |
| **M2** | **Fixed** | §4.1 (iii) lists: lv5d/hr5d, base constants, `seed_blockers`, step list, preflight, single-process force list, reference index ≠ shard index, rep-based timing rule. §3.2 uses the right CLI names. |
| **M3** | **Fixed in substance** | Definitions are unchanged, the pool is regenerated, the pinned hashes are dropped, and a pre-freeze re-check is declared. The rule quoted is stale (R2-1). |
| **M4** | **Fixed** | Option II is binding once (v) is frozen, and the S7 fallback is removed. |
| **M5** | **Fixed** | The scan of 337 runs gives used seeds 20506000, 20517000, 20521000/1 and 20525000. I re-derived "no shared seed for k ≤ 995" for 20511004/20521004: the first shared seeds are at r = 5/9/13, k = 996/997. The collision claim is withdrawn. Retiring the bases unused, declared before observation, is a legitimate conservative choice and is stated honestly. |
| m1–m7 | **Fixed** | – |

## Remaining findings

### R2-1 (major): the V7-R2 size-selection rule is the superseded revision-2 window rule and contradicts the route's own baseline
- **Evidence.**
  - The route (header; §2.4 "Surrogates"; §4.1 (vi) step 2) says: freeze 5400/0.005 if both κ = 1 coverages lie in [0.680, 0.686], else 7200/0.005 if [0.677, 0.689] holds.
  - Repair plan rev 3 (R2-1, §3.8 item 2) replaced this with: **5400/0.005 iff P(power clause) = P(Bin(3200, pw) ≥ 2909) ≥ 0.99 for both surrogates, using the exact binomial OC at the measured κ = 1 coverages; else 7200/0.005 if the same computation gives ≥ 0.99; else no design, and the finding goes to the operator.** The plan states that the windows are illustrations, not the rule.
  - The current archive single coverage, 0.68664, lies outside [0.680, 0.686]. The route's rule would therefore select 7200/0.005 (≈ 44 h; re-qualification 47.3 h) on today's numbers. That contradicts its 5400 baseline and every budget derived from it.
- **Fix.** Quote the rev-3 rule verbatim in all three places and cite "repair plan rev 3 §3.8". Keep the windows only as illustrations.

### R2-2 (minor): budget figures lag repair plan rev 3
- **Evidence.** Rev 3 gives 5400/0.005 ≈ 33.9 h (transport 8.7 h including the extra shard start-ups), (B) ≈ 16.5 h and (C) ≈ 12.5 h. The route still uses 33.7 / 16.3 / 12.3 h in §2.2, §2.3, §2.4, (iv) items 4 and 8, (vi), §4.2, §4.3 and Appendix A.
- **Fix.** Recompute:
  - re-qualification 3.28 + 33.9 = 37.2 h; (B) 19.8 h; (C) 15.8 h; 7200 fallback 47.3 h;
  - totals with D-E4: 41.3 / 23.9 / 19.9 h.

### R2-3 (minor): stale transport layout text
- **Evidence.** §2.2 Option I still says "24 shards × 4.5e6 (5400 pairs; halved shards …)". The §2.4 claim list still says "the transport size (16 → 24 shards of 4.5e6)".
- **Fix.** Say 48 × 2.25e6 (64 × 2.25e6 under the fallback) in both places.

### R2-4 (minor): seeds for the 7200 fallback are not stated
- **Evidence.** §4.4 gives only the 48-shard seeds. The fallback layout uses k = 0..63, i.e. 20547004–20547067, with the reference at k = 64, i.e. 20547068.
- **Fix.** State this range. Note that the hr r_index 6 reservation at 20541004 covers k ≤ 64.

### R2-5 (minor): the route cites the wrong revision of the repair plan
- **Evidence.** The route cites "repair plan revision 2" throughout; the current revision is 3. Repair plan rev 3 §5 still says lv 20511004 "lie[s] inside the used range (exact collisions wherever …)". That is conditional and correct, but the route asks for "fails the range rule" wording, and it should be aligned at commit.
- **Note.** The elastic plan's Amendment 19 seed list (§6 (f): "fresh bases 20511004 / 20521004 … maximum 20542003") and its closing line still need the edit the route specifies: 20531004 / 20541004, maximum 20562003, which I verified as hr 20541004 + 20999. The coordinator has already noted this.

### R2-6 (minor): fixture composition change is not on the amendment list
- **Evidence.** The halved fixture (4 × 2.25e6, seeds 20521000–20521003) changes the composition named in Amendment 14a (b)2 ("shards v7r-s0 and v7r-s1").
- **Fix.** Declare it in the V7-R2 amendment or the (v) amendment, not only in (vi).

## Consistency checks passed
- lv 20531004 r ≤ 20 gives at most 20552003; hr 20541004 r ≤ 3 gives at most 20544010.
- V6-F seeds 20558004–20558015, 20559004–20559011 and 20560004–20560019 are disjoint from both bases' used r_index blocks, from 20547004–20547068, from 2098xxxx/2099xxxx and from TOPAS 2028xxxx–2029xxxx. They equal only the unused hr r_index 17–19 and lv r_index 27–29.
- Cross-base identities stated correctly: lv 10 ≡ hr 0, lv 11–13 ≡ hr 1–3, lv 16 ≡ hr 6.
- Sequence: no V7-R2 host step (including the fixture rebuild) runs before (iv) is answered, which is consistent with repair plan rev 3. Throughput comes first, and the sizes addendum is frozen by a pre-declared rule before any V7-R2 result.

## Verdict (round 2)
**REVISE, minor.** No blocking issue remains. R2-1 should be corrected before commit because it changes the selected design and budget on current numbers. R2-2 to R2-6 are text edits. If R2-1 is applied as stated, no further review round is needed.
