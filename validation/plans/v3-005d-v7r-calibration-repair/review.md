# Independent review — V3-005D V7-R calibration repair plan (v7r-calibration-repair-plan.md)

Reviewer: independent statistical/scientific review, 2026-10-14. Worktree v3-005b HEAD 749f768c (read-only).
Reviewer scripts (exploratory, not evidence): `validation/plans/v3-005d-v7r-calibration-repair/review-scripts/oc.py` (independent CP-TOST OC via the binomial
tail identity: accept iff P(X>=k|0.6527) <= a and P(X<=k|0.7127) <= a, which avoids CP bisection), `consistency.py`,
`sens.py`, `pall.py` (seed 20989995, inside the planner-diagnostic block), `kappa_spot.py`.

## What was verified (and agrees with the plan)

* Gate in code (`v7r.pair_gate`, `single_gate`): two one-sided Clopper–Pearson bounds at
  `steps_v5b.V7_ALPHA_TOST = 0.04` per side; pass iff lower >= 0.6527 and upper <= 0.7127. Single gate: the same
  on K_min/K_max over the H box. The plan's model (exact binomial TOST) is the code's gate.
* n = 3600, a = 0.04: acceptance set 2401–2517; FA 0.0374 / 0.0383; power 0.9636 at 0.6827. P(<=109/3200 | 0.0383)
  = 0.113 (confirmed); at 0.0374 it is **0.172** (plan: "≈ 0.15"). P(all four paired margin clauses) = 3.7e-4 with exact κ.
* κ-SE: sqrt(c(1−c)/n) = 0.0034 (pairs, 2e4), 0.0039/0.0037 (single, 15 000 — `single_ingredients` does apply the
  3:1 split inside `margin_kappas`; confirmed in code). The κ of record reproduce bit for bit (reviewer spot check:
  rehearsal paired_low 0.9493916971221448 from the HEAD fixture, 20 s; planner log shows all ten equal).
* Direct-test table in §2.2 matches `validation/plans/v3-005d-v7r-calibration-repair/data/analyse_direct.out`. Consistency of the observed paired counts with the
  exact OC at the directly measured coverage: chi2 = 14.5 on 10 rows (p ≈ 0.15); one row in tension
  (rehearsal paired_high: predicted 0.0095, observed 69/3200 = 0.0216, z ≈ +3.1; implied coverage 0.7146 vs direct
  0.7170, 2.4 SE of the direct sample). The pattern across the rehearsal rows is consistent with one ≈2σ common-mode
  error of the single direct sample. The diagnosis — κ-solve sampling error of the frozen 2e4 sample (rehearsal
  paired +1.6σ, single −1.8σ; archive ≈ 0.4σ) plus defect (a) for the archive — is sound; surrogate misspecification
  cannot act inside the calibration (same sampler for κ and MC). 11.5 % at true coverage ≈ 0.658 is what the gate
  should do there.
* Recommended design 5400 pairs / a = 0.005: acceptance 3615–3762 (0.6694–0.6967); FA 0.0050/0.0050; paired pass
  probability 0.969 at 0.6827, 0.9595 at the observed 0.6853 (frozen 3600/0.04: 0.964/0.957). Independent P(all
  clauses, both surrogates, r = 1, single gate unshrunk) = 1.000 at n_κ = 1e6 (0.077 at n_κ = 2e4 — the κ step is the
  crux). Option (B) R = 800: 0.86 (plan 0.83); R = 1600: 0.993. κ error of 3 SE (1e6) raises FA only to 0.0092.
* Lowering a_tost is not a relaxation of the false-acceptance guarantee (0.04 → 0.005) and the real-data pass
  probability of a correct estimator is preserved by the larger n. It is, however, a change of a frozen gate
  parameter (14 (h): "Clopper–Pearson at α_tost 0.04 per side"), i.e. a criterion change — see F1.
* Host arithmetic: 24 × 1262 s + 280 s = 8.5 h; MC 6400 × 12.2 s = 21.7 h + 78 × 70 s = 1.5 h; κ 1.1 h; total 33.7 h
  (consistent).

## Findings

**F1 (blocking) — "needs operator: no" for V7-R2 is not defensible as written.**
Evidence: 14 (h) Outcome "Fail: a recorded scientific failure without a further base"; 14 (h)4 allows count increases
"in Amendment 14a before the evidence run, never after"; 14a (b)6 "the C9 run … is the calibration of record; a miss
there fails V7-R and no further count increase follows". The plan cites 14 (h)4 as authority for 16 → 24 shards (§3.6,
§4) — that permission expired at the evidence run. V7-R2 is a new gating row (new α, new n, new κ procedure, new base)
for the same estimator claim. The operator's disposition A reserves scope changes ("bring back any concrete proposal
that changes required validation scope") and, for the analogous V6-F, ruled "V6-F passing does not authorize
replacing V6 or merging". Whether a V7-R2 pass discharges the V7-R / V2-NUM requirement in the merge gate 14 (a)6 is
exactly such a scope decision. Also the IR told the operator "V7-R calibration ≈ 9.5 h"; the plan needs 33.7 h.
Fix: classify as: diagnosis record, design, pre-registration (Amendment 19), implementation and review = autonomous
(the operator-prioritised "V7-R calibration repair"); **V7-R2's role relative to the failed V7-R in the merge gate /
requirement ledger and the 33.7 h budget = operator proposal**, sent before host execution (host execution is anyway
gated by the independently reviewed integration route). State V7-R2 as a new row; remove the 14 (h)4 citation.

**F2 (important) — transport shards violate the plan's own step-time target.** 4.5e6-history shards take 1262 s on the
slow host (84 % of the 1500 s hard limit, 16 % headroom), against the plan's 1100 s planning target; 14 (h) halved the
lv5b shards precisely because 13 % headroom was judged insufficient, and an elastic-model change from the sibling plan
may raise the per-history cost. Fix: 36 shards × 3e6 (300 replicates; pairs s ↔ s + 18; H = shards 27–35, 2700; E =
0–26, 8100; ≈ 841 s each), same totals, or justify the exception explicitly in the amendment.

**F3 (important) — "both surrogates unchanged" and P(all) ≈ 1 are conditional on the rehearsal pool, which will be
regenerated.** The rehearsal fixture must be at the evidence digest (14a (b)2, 17 (a)7); if the elastic repair (same
priority under disposition A) lands first, the escaped-neutral block sums change and pool sha 5d28a9b3…/d4f7d749…
no longer holds. The power clause is the binding one: joint power ≥ ≈ 0.92 (needed for P(≥ 2909/3200) ≈ 0.99) requires
the surrogate's κ = 1 paired coverage within ≈ [0.679, 0.687] at 5400/0.005 (paired OC 0.911 at 0.678, 0.890 at
0.689) — no wider than the frozen design's window, and the single-gate box shrink is a fit at α = 0.04. Fix: fix the
order (elastic decision → rehearsal regeneration → measure κ = 1 coverage and single-gate power on the regenerated
pool in the sandbox → freeze Amendment 19 numbers), or pre-register the 7200/0.005 variant (window ≈ ±0.007), and
state the power risk explicitly.

**F4 (important) — fresh-base/seed statement is under-specified.** V7-R2 transport seeds 20511004 + 16000 + k =
20527004–20527028 are disjoint from the rehearsal seeds actually used (base 20505000 + 1000 r + k, r ≤ 20, ≤ 20525999),
from 2098xxxx and from 2099xxxx — correct. But "20511004 if it is still unobserved" is the wrong criterion: the V6-F
follow-up already found 20511004 colliding with the rehearsal envelope (r_index 0–14 lie inside 20505000–20525999;
lv r ≡ rehearsal r + 6, k + 4), and lv 20511004 r16 ≡ hr 20521004 r6 (cross-base identity; hr 20521004 r_index 0–4 also
lie inside the rehearsal envelope). Also V7-R2 shares the lv base with the elastic re-qualification (different r_index),
so "contention" is a sequencing, not an ownership question. Fix: adopt the V6-F checks (a)–(e) verbatim (seed-by-seed
disjointness from consumed and rehearsal seeds, cross-base identities reserved: hr r_index 6 at 20521004 must stay
unused), and a test asserting them.

**F5 (important) — reproducibility of the design tables.** `designs3.out`, `designs4.out`, `designs5.out` have no
scripts in the scratchpad (only `designs.py`, `designs2.py`); `designs5.out` labels the (B) R-reduction block "with the
recommended 7200/alpha .005" while §3.5 attributes the numbers to 5400/0.005 (reviewer: 0.86 at R = 800 for 5400/0.005,
close but not identical). Fix: commit the generating scripts with the amendment, regenerate every table cited in
Amendment 19 from them, pin the OC numbers in tests (as §6.7-3 proposes).

**F6 (minor) — reserved retry block 20985000–20988999.** A pre-declared second attempt "only if an implementation
defect voids the first" invites re-labelling a statistical miss. Fix: require an independently reviewed defect with a
reproducing test that does not depend on the counts, the first attempt's counts reported, and an operator notice.

**F7 (minor) — single-gate bound wording.** §1.2 says the single gate's FA is "at most" the binomial OC; it is at most
OC + P(box misses δ*) ≤ OC + 0.01 (the 0.99 box level rests on the A13 Gaussian assumption). §3.5 (C) states it
correctly; align §1.2/§3.2 ("FA single (bound)" 0.0050 should read ≤ 0.0050 + 0.01 rigorous, 0.0050 given box coverage).

**F8 (minor) — numeric/consistency nits.** P(≤109 | 0.0374) is 0.172, not ≈ 0.15. 3200/39 parts gives 83/82
repetitions, not 84/83. "1100 s leaves 36 % margin" vs "27 % headroom" (1 − 1100/1500 = 27 %). Acceptance test 2
(bit-identical κ at 2e4 with old seeds) conflicts with the new full-sample single κ unless a legacy mode is kept —
say so. The 1.67× host factor is from transport; the per-repetition calibration rate (12.2 s) is extrapolated, not
measured on the slow host — measure one short calibration part before freezing part sizes. "the observed rate is
exactly what the gate should do" (§1.4) — predicted 0.137 ± 0.02 vs observed 0.115; say "consistent within 1σ".
Report the rehearsal paired_high 2.4σ tension rather than "the OC reproduces the observed counts".

**F9 (minor, options B/C).** Their OCs are honestly stated and correctly routed to the operator. Deterministic
certification of the paired gate is mathematically legitimate (hits are i.i.d. Bernoulli by construction, OC is
exact and κ-/surrogate-free) and strictly more informative than the κ-perturbed MC for the paired FA; but the MC also
served as an end-to-end implementation check — (C) should keep one (e.g. a small MC at fixed known coverage, or
pinned acceptance-set tests). The single-gate bound in (C) is assumption-bearing (box level) and must be labelled so.

## Verdict

**Revise.** The diagnosis and the statistical design are sound and independently reproduced; no frozen tolerance is
relaxed. One blocking issue remains (F1: governance classification of V7-R2 and its budget); F2–F5 should be fixed
before Amendment 19 is frozen. Once F1 is reclassified as an operator proposal for V7-R2's merge-gate role and budget,
no blocking issue remains.

---

# Round 2 (Revision 2 of the plan; rev-1 backup `v7r-calibration-repair-plan.rev1.bak`)

Reviewer checks (scripts in `validation/plans/v3-005d-v7r-calibration-repair/review-scripts/`; `seeds_r2.py` new):

* **F1 — resolved.** V7-R2 is declared a new row; V7-R stays failed; the 14 (h)4 citation is removed; merge-gate
  role (14 (a)6 / V2-NUM / V1-MUST-030) and the ≈ 33.9 h budget (vs 9.5 h in IR-DEFAF4) are operator-proposal items
  sent before any V7-R2 host step (§3.4, §4 rows added, §6.2, §6.4 step 3). The V6-F precedent is cited.
* **F2 — resolved.** 48 × 2.25e6 (225 replicates/shard, pairs s ↔ s + 24, H 36–47, E 0–35, reference k = 48).
  Checked against the code: `steps_v5d.shard_plan_c` uses limit 1500 s / `steps_v5.MARGIN` 1.25 = 1200 s;
  1.08e8 / 3566 hist/s (4.5e6 / 1262 s) / 1200 = 25.2 → 26 needed ≤ 48 declared (ok); 2.25e6 / 3566 = 631 s.
  Revision 1's 24 × 4.5e6 would indeed have failed (26 > 24). 2.25e6 = 225 × 20 × 500 is block-consistent.
* **F4 — resolved.** 20531004 + 16000 + k, k = 0..48 → 20547004–20547052. Numerically disjoint from: consumed
  bases 20441004/20451004/20471004 (r 0–11, k 0–32), lv5c 20481004 and hr5c 20491004 envelopes, the rehearsal envelope
  20505000–20525999, the never-run V6 seeds, the skipped 20511004/20521004 blocks, hr 20541004 r 1–3, other lv
  20531004 rows, the V6-F candidate B = 20551004 (r 7–9), 2098xxxx/2099xxxx. The only overlap is the declared
  cross-base identity hr 20541004 r_index 6 (49 seeds), which the plan reserves and tests. V6-F checks (a)–(e) adopted.
  Wording nit: §5 says lv 20511004 "r_index 10–14 lie inside the used range"; numerically r_index 0–14 do
  (20511004 itself is inside) — irrelevant now that the base is skipped, but correct it.
* **F5 — resolved.** `designs3.py`–`designs6.py` re-run: all four outputs byte-identical to `validation/plans/v3-005d-v7r-calibration-repair/data/designs*.out`
  (42 s, 11 s, 21 s, 5 s). (B) 5400/0.005 R 800 / 1600: 0.859 / 0.995 (r = 1), matching my 0.86 / 0.993.
* **F6–F9 — resolved** (retry-block conditions, "+ 0.01 box miss" wording, 0.172, 83/82 parts, 27 % headroom, legacy
  κ mode, slow-host calibration timing to be measured, §1.4/§2.2 wording with the 2.4σ residual, (C) end-to-end MC).

## Remaining findings

**R2-1 (important) — the §3.8 pre-freeze rule is computable and pre-registered, but as written it already selects
7200/0.005, contradicting the recommendation and the budget sent to the operator.** Rule item 2 requires "for both
surrogates, both κ = 1 coverages [paired and single] in [0.680, 0.686]". The archive-fitted surrogate is fixed (not
regenerated) and its κ = 1 single coverage is already measured at **0.68664 ± 0.0012** (`analyse_direct.out`), outside
0.686. So the 5400 branch fails today and the rule switches to 7200/0.005 (≈ 44 h); the "33.9 h" recommendation and the
operator budget item are inconsistent with the plan's own rule. The window is also the wrong criterion for the single
coverage: `designs6.py` derives it with paired = single coverage, but the single gate runs on n_E = 8100 and is far
more powerful; at the measured (paired 0.6840, single 0.6866) the 5400/0.005 joint power is 0.954 and P(power clause)
≈ 1.000. Fix: replace the coverage-window condition by the direct criterion already in item 2 — P(power clause) ≥ 0.99
computed from the measured paired and single κ = 1 coverages, each moved one measurement SE toward the worse side, and
the measured single-gate power taken at its CP95 lower bound (R ≈ 200 has SE ≈ 0.01). Or keep the windows but apply
them to the paired coverage only. Then state which branch the current archive measurement selects, and give the
operator both budgets (33.9 h / 44 h) with the rule that picks between them. The rule is not post hoc: it uses
surrogate inputs only and is fixed before any V7-R2 result. But it must be stated before the regenerated-pool
measurement, and the archive measurement is already known — say so.

**R2-2 (minor).** The §5 wording nit above (r_index 0–14, not 10–14). §3.4 says executing V7-R2 is "autonomous (when
scheduled)", while §3.4/§6.4 make host execution wait for the operator's budget decision. One sentence should say that
host execution waits for that decision.

## Round 2 verdict

**Revise (minor revision).** No blocking issue remains: F1–F9 are fixed and verified (the scripts reproduce, the shard
arithmetic matches `shard_plan_c`, the seeds are disjoint except the reserved identity). Before Amendment 19 is frozen,
R2-1 must be corrected so that the selection rule, the recommended size and the budget sent to the operator agree.

---

# Round 3 (Revision 3; backup `v7r-calibration-repair-plan.rev2.bak`)

* **R2-1 — resolved.** §3.8 item 2 now reads: 5400/0.005 iff P(Bin(3200, pw) ≥ 2909) ≥ 0.99 for both surrogates,
  where pw = OC_paired(c_p) × OC_single(c_s). The OCs are exact binomial with the conservative box shrink, and the
  sandbox-measured single-gate power is used where it is lower. Otherwise 7200/0.005 if that passes, else operator. The
  windows are labelled as illustrations only. `designs7.py` re-run: output byte-identical to `designs7.out`. pw =
  0.9537 / 0.9618 (5400) and 0.9937 / 0.9954 (7200); P(power clause) = 1.0000 everywhere. So 5400 is selected at the
  current coverages, matching my round-1 independent value 0.954 / 0.962.
  The rule is computable and pre-specified. Its inputs are surrogate measurements only; it is fixed before the
  regenerated-pool measurement and before any V7-R2 result; it is re-applied to the rebuilt pool. The fact that the
  archive inputs are already known is stated openly.
* **R2-2 — resolved.** §5 now says r_index 0–14. Host execution of any V7-R2 step waits for the operator's decision on
  role and budget (header, §3.4, §6.4 step 3). Both budgets go to the operator: 33.9 h / 44 h, plus (B) 16.5 h and
  (C) 12.5 h. Seeds 20547004–20547052 at lv 20531004 are as verified in Round 2; V6-F at 20551004 r 7–9 is disjoint.

**R3-1 (minor, non-blocking).** The rule uses point estimates. The measured single-gate power from R ≈ 200 has SE
≈ 0.01, and the coverages have SE ≈ 0.001. For a borderline rebuilt pool, P(power clause) should be computed with the
coverages moved one SE toward the worse side and the MC power at its CP95 lower bound. This is immaterial at the
current margin (pw 0.954 vs the ≈ 0.915 needed). Optional: the §3.8 header still says "(Revision 2, F3)".

## Final verdict

**Accept.** No blocking issue remains. The diagnosis, the design (5400/0.005, n_κ = 1e6, R = 3200, 48 × 2.25e6), the
governance split (design autonomous; V7-R2's merge-gate role and budget are operator proposals sent before any host
step), the seeds and the reproducibility have been checked independently. R3-1 can be folded into Amendment 19 when
it is written.
