# V6-F — independent review checklist (revision 2)

The reviewer answers each item yes, no or finding, with the evidence used. Any "no" blocks committing the
V6-F items of Amendment 18. Items marked [r2] are new or changed in revision 2. The revision-1 review is
`review.md`.

## A. The original failure is preserved

1. Does any text, code path or ledger entry change the verdict "failed — not sized under the
   pre-registered rule" of V6 and V6-sens (Amendment 14a (a)4), `V6_SHARDS = None`, or the `v6-unsized`
   steps? The required answer is no.
2. Is V6-F kept out of the Amendment 14 (a)6 list? Does every report show the V6 failure next to any V6-F
   verdict (design 1.2)? Is the P6-D replacement pattern rejected?
3. [r2] Does design 1.3 leave the (a)6 treatment of "failed — not sized" to the separate orchestrator item
   of Amendment 18, with no ruling inside V6-F?
4. [r2] Are the never-run original V6 seeds left unused? These are lv5c r_index 12 and 20 and hr5c
   r_index 3 of 20481004 and 20491004, i.e. 20493004+, 20501004+ and 20494004+.

## B. Pre-registration completeness

5. [r2] Is every gating quantity fixed? Check each of the following.
   - z_F = 21.6 cm.
   - The 1e-4 floor on the measured r = 0 maximum, 5.95 MeV/g/p.
   - The aura share ≤ 0.05, with D_au from p23–p25 by the stated formula.
   - The observable, Amendment 17 (c)1.
   - E₀ = 177.84075737 MeV.
   - The tolerances and sign-test groups.
   - The screen metric with measured denominator and σ′ as the paired SE of the difference.
   - The coverage minimum: ≥ 40 points, ≥ 3 radii, ≥ 12 depths after ties.
   - N_F = 16 and k = 8.
   - The r_index mapping 7 / 8 / 9 and the backend rule.
6. Is "not established (beam-sensitive)" defined and stated never to be a pass?
7. [r2] Are the contrary readings (a)–(e) fixed now and verdict-neutral? Is the bias of (a) declared
   beforehand?
8. [r2] Is the zero-mean rule dropped (R2-m2), with the reason stated? Is the rule that realized SE is
   reported and never acted on stated?

## C. No post-hoc tuning

9. [r2] Were the D7 aura rule (5 %), the D3 measured denominator and the D8 coverage numbers fixed by
   orchestrator decision before the revision-2 script was run? Check the timestamps of `review.md`, the
   orchestrator message and `se_table.json`.
10. [r2] Contrary-evidence notes. Two are recorded:
    - Amendment 14a (a)5 (a side-by-side display);
    - the revision-2 note in design.md: revision-1 model-normalised σ_Δ (0.0988 at (21.6, 3)) and
      revision-2 measured-normalised σ′ (≤ 0.037) were both seen, which implies the 1e6 model is low there.

    Is any selection rule plausibly influenced by either? Is the D1 disclosure (convention, weak budget
    discriminator, high-risk points kept, reading (b)) adequate?
11. [r2] Read `se_table.py`.
    - Does it compute only precision quantities, the measured fraction and the aura share?
    - Are `ratio_paired` and any model mean absent?
    - Are the model-normalised sensitivity SE and the measured-normalised σ′ not both present?
    - Is the remaining weak leak (rel_se of v6-s0 next to σ′/D_meas of the sensitivity sample) acceptable?
12. [r2] Is the residual dependence of the D3 screen on the model (σ′ ∝ √D_model) correctly stated as
    second order and opposite in direction to the revision-1 bias?
13. Does any retained gating point get a tolerance other than 2 % (r = 0) or 15 % / 90 % / 10 % and
    p > 0.01?

## D. Arithmetic

Recompute with `python3 -I se_table.py <worktree>` and `python3 -I core_sens_estimate.py`.

14. [r2] Point counts.
    - P_F,0 = 29.
    - H_F before D7: 145 points, 24 depths.
    - Aura share by r: medians 0.03 / 0.30 / 0.91 / 3.2 / 6.0 / 10.1 % at r = 1–6 cm.
    - D7 excludes 41 points (r = 3 / 4 / 5 / 6 cm: 2 / 4 / 18 / 17).
    - H_F = 104 points (r = 1–5 cm: 30 / 29 / 26 / 16 / 3), 16 sign-test depths.
15. [r2] Sizing.
    - r = 0: rule 8 (12 with safety).
    - Halo: rule 6 (9 with safety).
    - N rule 12, declared 16.
    - σ′/D_meas worst 0.0369 at (20.4, 5) → k rule 1, declared 8 (2σ′ worst 0.026).
    - CPU planning 72 × 184 s = 3.68 h.
    - A2 at 21.8 cm: N 27, 4.24 h.
16. [r2] Sign-test boundaries: at 12 depths ≥ 11/12 fails (p = 0.0063); at 16 depths ≥ 14/16 fails
    (p = 0.0042). Power at q = 0.7 / 0.8 / 0.9 is 0.085 / 0.275 / 0.659 (12 depths) and
    0.099 / 0.352 / 0.789 (16 depths).
17. [r2] Coverage expectation (planning): ≈ 58 points, r = 1 / 2 / 3 / 4 / 5 cm with 17 / 0 / 22 / 16 / 3,
    and 16 depths.
    - Is the method acceptable for a planning statement? It uses a core-only Gottschalk T_dM with
      X_S = 46.88 cm, sphere averaging, the core scaled to the measured r = 0 dose and the st-hi proxy for
      σ′.
    - Is the margin to the minimum (58 vs 40; 16 vs 12) adequate?
18. Is the st-hi proxy for σ′ of the σ_x and ρ variants acceptable, given that errors only shrink H_g?

## E. Seeds (single later execution)

19. [r2] Is the scheme fully declared?
    - Base B is fixed in the executing task's amendment before observation. It is the next unused base
      that passes checks (a)–(e); skipped bases are recorded and never reused.
    - Seeds are B + 7000 + k (k = 0..15, CPU), B + 8000 + k (k = 0..7, all seven CRN runs) and
      B + 9000 + k (k = 0..15, CUDA).
    - r_index ≤ 9, which avoids the cross-base collision at r_index ≥ 10.
20. [r2] Does the disjointness check use the rehearsal **seeds actually used**, up to 20525999 for the C8
    rehearsal, rather than the base range? Is B = 20511004 correctly identified as colliding
    (B + 7000 = 20518004)?
21. Are the CPU and CUDA samples disjoint, so that the paired z is valid?

## F. Scope and claims

22. [r2] Aura: the transport deposits no neutral energy. Does D7 remove the points where the source's
    own aura term exceeds 5 %? Is the role extension of `tbl:parms` (p23–p25 used for selection) stated,
    and its recording assigned to the executing task?
23. Is the r = 0 distal-edge risk (0.17 / 0.08 / 0.05 mm per 2 % at 21.2 / 21.4 / 21.6 cm) stated in D1?
24. [r2] Does design 10 avoid any claim about:
    - the aura;
    - r ≥ 6 cm or points with aura share > 5 %;
    - depths beyond 21.6 cm;
    - beam-fit independence (the beam, E₀ and the aura selection all have shared lineage);
    - normalisation better than ±3 %;
    - other energies;
    - the frozen V6 claim?

    Is the r = 0 statement qualified "given the fitted beam"?

## G. Process

25. [r2] Is it stated that V3-005C runs no V6-F step, adds no V6-F code to the lv5c/hr5c hashed set and
    spends no V6-F host time?
26. [r2] Is it stated that V6-F is executed exactly once, gating, in the later task, with no second
    execution without a new reviewed physics reason and a new pre-registration?
27. Will the V6-F items of Amendment 18 be committed and reviewed before any V6-F shard? Is
    `shard_table_ok` required at the executing task's evidence head?
28. [r2] Is `v1/se_table.py` neutralised (no per-point precision column, no σ_Δ)? Is every revision-1
    per-point table quarantined (never committed or archived next to revision 2)?
