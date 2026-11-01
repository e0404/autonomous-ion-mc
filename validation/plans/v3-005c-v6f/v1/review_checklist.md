# V6-F — independent review checklist

The reviewer answers each item yes, no or finding, with the evidence used. Any "no" blocks committing
Amendment 18.

## A. The original failure is preserved

1. Does any text, code path or ledger entry change the verdict "failed — not sized under the
   pre-registered rule" of V6 and V6-sens (Amendment 14a (a)4), `V6_SHARDS = None`, or the `v6-unsized`
   steps? The required answer is no.
2. Is V6-F kept out of the Amendment 14 (a)6 merge-gate list of V3-005C, and does every report show the V6
   failure next to any V6-F verdict (design 1.2)?
3. Does the text avoid the P6-D replacement pattern (Amendment 15 (a)5)? If it does not, is the deviation
   justified against the operator response to IR-20261011-031924-C205C8?
4. Are seeds lv5c r_index 12 and 20 and hr5c r_index 3 left unrun and reserved?

## B. Pre-registration completeness

5. Is every gating quantity defined without free choices left for after the data? Check each of the
   following.
   - The point set: z_F = 21.6 cm from R0 + σ_sem; the 1e-4 floor; the measured r = 0 maximum,
     5.95 MeV/g/p at z = 21.0 cm.
   - The observable, which is Amendment 17 (c)1 unchanged.
   - E₀ = 177.84075737 MeV, unchanged.
   - Tolerances, aggregation and sign-test groups, with ties as in Amendment 17 (c)3.
   - The screen metric.
   - The coverage minimum: ≥ 40 points, ≥ 3 radii, ≥ 10 sign-test depths.
   - N_F = 16 and k = 8.
   - Seeds, backends and the disagreement rule.
6. Is the outcome "not established (beam-sensitive)" defined and stated never to be a pass?
7. Are the contrary readings (design 10 a–c) fixed now and declared verdict-neutral?
8. Is the zero-mean rule of the screen (a zero nominal or variant mean makes the point sensitive) stated?
9. Is the rule "realized SE is reported and never acted on" (no top-up, no point dropped for precision)
   stated?

## C. No post-hoc tuning

10. Contamination check. Amendment 14a (a)5 records that the orchestrator saw measured and model mean
    doses side by side for the 1e6 rehearsal shard.
    - Can each selection rule (z_F, the floor, the screen, the coverage minimum) be derived from frozen
      quantities and the physics argument alone?
    - Would a reviewer who never saw that display arrive at the same rule?
    - Check in particular whether 21.6 cm rather than 21.8 cm or 22.0 cm could have been chosen because
      of known agreement at 21.8 cm.
    - Confirm that contrary reading (b) reports the 21.8 cm verdict.
11. Does `se_table.py` compute only precision quantities (relative SE, non-empty batches, paired σ_Δ)?
    Read it. It must have no model mean, no residual and no ratio |D_var/D_nom − 1| in its outputs.
12. Was the rehearsal sensitivity ratio (st-hi vs nominal) ever computed or displayed? If yes, the screen
    design must be re-justified.
13. Was the core-sensitivity estimate (`core_sens_estimate.py`) computed from beam parameters and an
    analytic MCS model only, with no ionmc output and no measured dose?
14. Is the per-point screen (design D3) a change of consequence only (point leaves the gating set), with
    the frozen metric and the 0.15 tolerance unchanged? Could it select points by agreement? It must not:
    it is model against model.
15. Does any retained gating point get a tolerance other than the frozen 2 % (r = 0) or 15 % / 90 % / 10 %
    and p > 0.01 (halo)?

## D. Feasibility arithmetic

Recompute these with `python3 -I se_table.py <worktree>`.

16. Reproduction of Amendment 14a: halo clause 88 shards at (22.0, 3); r = 0 clause 116 / 651 / 8910 / ∞;
    V6-sens 32 shards on the 153-point subset at σ_Δ ≤ 0.075.
17. V6-F numbers:
    - P_F,0: 29 points, 8 shards (12 with the × 1.25 safety factor).
    - H_F: 145 points (30 / 29 / 28 / 20 / 21 / 17 at r = 1–6 cm), 10 shards (15 with safety).
    - 24 sign-test depths.
    - Worst σ_Δ = 0.0988 per 1e6 → k = 7 for 2σ_Δ ≤ 0.075; k = 8 gives 0.070.
18. Budget: 72 CPU shards × 184 s = 13 248 s = 3.68 h ≤ 4 h. Each step is 184 s against the 1200 s
    effective limit. CUDA is 16 shards (planning 10× → 294 s; worst case 1× → 2944 s), and the 10× is
    stated as an unmeasured assumption.
19. Is SE ∝ N^−1/2 scaling from 20 batches of 5e4 histories adequate at the worst H_F point (1.533 cm,
    3 cm; 20/20 batches non-empty)? Is the × 1.25 safety factor sufficient? A 20-batch SE has about 16 %
    relative uncertainty for Gaussian batches, and more for skewed ones.
20. Is the single rehearsed variant pair (st-hi) an adequate proxy for the σ_x and ρ variants' σ_Δ? What
    happens if it is not? Points drop from H_g, and the coverage minimum then decides.
21. Expected H_g size: does the core estimate suggest the coverage minimum is reachable? Points at r ≥ 3 cm,
    z ≲ 15 cm, and r = 1 cm at z ≳ 18 cm are expected beam-insensitive. If the estimate suggests
    coverage fails, say so now.

## E. Seeds

22. Are the seeds 20502004–20502019 (lv5c r21), 20503004–20503011 (lv5c r22, CRN across seven runs) and
    20504004–20504019 (hr5c r13) disjoint from all of the following?
    - `consumed_seed_set()`, whose maximum is 20482036.
    - The original V6 reservations: 20493004+, 20494004+, 20501004+.
    - Every other lv5c/hr5c r_index of Amendment 14 (j).
    - The rehearsal family 20505000–20509999, including the E-shape seeds 20509500+.
    - The fresh bases 20511004, 20521004, ….
23. Are lv5c r_index 23 and hr5c r_index 11 and 12 declared unused, and is lv5c r_index ≥ 24 forbidden?
24. Do the compare steps fail closed on any overlap? Do the CPU and CUDA samples have disjoint seeds, so
    that the paired z is valid?

## F. Scientific scope and claims

25. Is the aura argument (no neutral transport; the source's l.45, l.74, l.250) correct for the current
    transport? Check whether any neutral energy is deposited locally. If it is, the scope statement must
    change.
26. Is the r = 0 distal-edge risk disclosed? Retained points at 21.2–21.6 cm have measured log-gradients
    of 12–38 % per mm, so 2 % corresponds to 0.05–0.17 mm of range, against the source's ~1 mm
    positioning uncertainty (l.321). This is the price of respecting Amendment 14 revision 1.
27. Does design 12 ("can and cannot establish") avoid any claim about the aura, r ≥ 7 cm, depths beyond
    21.6 cm, beam-fit independence, normalisation better than ±3 %, other energies, or the frozen V6
    claim?
28. Is the shared beam lineage (p1, p5–p7 from the same `tbl:dmlg` fit) disclosed, together with the fact
    that E₀ comes from R0 = d80?

## G. Process

29. Is "report-only for V3-005C, gating for the later halo-claim task" argued? Is the single pre-declared
    gating execution in the later task (fresh base, criteria unchanged, earlier results reported next to
    it) stated so that it cannot become a retry?
30. Hashed-source digest: adding V6-F steps changes the lv5c digest. Is the V7-R rehearsal fixture of
    Amendment 14a (b)2 then stale? Is its re-run (about 0.45 h) planned before the C9 evidence head?
31. Is Amendment 18 committed and reviewed before any V6-F shard? Does the throughput step at the
    evidence head list the V6-F shard plans as fitting?
