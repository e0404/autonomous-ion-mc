# V6-F — draft pre-registration text (proposed Amendment 18)

Status of this text: a design draft for independent review. It is not committed and has not been executed.
Numbers marked (T) come from `se_table.json` / `se_table.py` in this directory. Those numbers were computed
from the three C8 rehearsal partials at base 20505000 (head c847cc75; never evidence) and from the pinned
`tbl:dmlg`. Only precision quantities were computed from the rehearsal. No model mean, residual or
sensitivity ratio was computed.

---

**Amendment 18 (draft, 2026-10-11; written after the C8 rehearsal and after Amendment 14a (a)4; before any
V6, V6-F or V6-sens-F evidence shard) — V6-F, a separately pre-registered and computationally feasible
follow-up to the failed row V6/V6-sens.** This amendment does not change the frozen row V6, the row
V6-sens, Amendments 14 (g), 17 (c) or 14a (a), or any tolerance, seed base or verdict recorded there. It
adds two new rows, V6-F and V6-sens-F, with their own scope and a documented list of differences
(item 4). It answers the operator response to IR-20261011-031924-C205C8, recorded verbatim in that
intervention's `response.json`.

## 1. The original failure stays recorded

1. Rows **V6** and **V6-sens** stay recorded as **failed — not sized under the pre-registered rule**
   (Amendment 14a (a)4). Steps `v6-unsized` (lv5c, hr5c) keep that verdict, and `V6_SHARDS` stays `None`.
   Seeds lv5c r_index 12 and 20 and hr5c r_index 3 stay reserved for the original rows and are never
   run.
2. V6-F never replaces V6 in the requirement ledger, the results README, the decision 0041 outcome or
   the merge-gate list of Amendment 14 (a)6. Every report that shows a V6-F verdict also shows, next to
   it: "V6/V6-sens: failed — not sized under the pre-registered rule (Amendment 14a (a)4)". Unlike P6-D
   (Amendment 15 (a)5), V6-F does **not** take the original row's place in the merge-gate list.
3. Lateral-halo claim. Decision 0041 (l.201) makes no lateral-halo claim until V6. That stays true for
   V3-005C whatever V6-F gives. A V6-F pass is reported as "V6-F (reduced scope, item 4) passed", never
   as "V6 passed" or as "halo validated".

## 2. Why the original rule could not be sized (recorded facts, not new findings)

- The frozen history-count rule (Amendment 14 (g)) requires 0.5 % relative SE at every r = 0 point and
  3 % at every r ≥ 1 cm point with measured dose ≥ 1e-4 of the r = 0 maximum. The r = 0 maximum is
  5.95 MeV/g/p at z = 21.0 cm (T).
- At r = 0 it applies at z = 22.0 / 22.2 / 22.4 / 22.6 cm. Measured dose there is 3.6e-2 / 6.0e-3 /
  6.5e-4 / 7.0e-5 of the maximum (T). The rule needs 116 / 651 / 8910 shards, and an unbounded number at
  22.6 cm, where the model dose was zero in 1e6 histories (Amendment 14a (a)2, reproduced exactly in (T)).
- The halo clause needs 88 shards because of one point, (z 22.0 cm, r 3 cm), at 28.1 % SE per 1e6 (T).
- Of the 261 halo points (Amendment 17 (c)6), 33 had zero model dose in v6-s0 and 32 had zero model dose
  in both V6-sens rehearsal shards (T).
- Every unsizable point lies beyond the Bragg-peak region (z ≥ 21.8 cm), or below 1e-4 of the maximum
  at large r. There the measured dose includes the **aura**: neutral secondaries and the charged
  particles they set in motion. The source says the aura dominates at large r and beyond range
  (halo.tex l.45, l.74, l.250). The transport has no neutral-particle transport: neutral energy is
  scored as `escaped_neutral` and is never deposited. So the model has no physics for these points.
  They are outside model scope, not merely short of statistics.

## 3. Rows

| # | Check | Pass iff | Class / tier | Justification |
|---|---|---|---|---|
| V6-F | Gottschalk 2015 `tbl:dmlg` (177 MeV, the Amendment 14 (g) beam, the E₀ of Amendment 14a (a)1, the observable of Amendment 17 (c)1) on the reduced point set P_F of item 4. Two components: r = 0 and halo | (i) r = 0: \|D/D_meas − 1\| ≤ 0.02 at **every** point of P_F,0. (ii) Halo, on the beam-insensitive subset H_g ⊆ H_F (V6-sens-F) and only if the coverage minimum of item 6 holds: per point 0.15; ≥ 90 % of points pass; median \|residual\| ≤ 0.10; no r-trend (two-sided binomial sign test over depths with points in both the 1–3 cm and the ≥ 4 cm groups, p > 0.01, ties excluded and counted). (iii) Both backends give the same verdict on every component (Amendment 17 (c)8 rule). V6-F passes iff (i), (ii) and (iii) all hold. If the coverage minimum fails, V6-F is **"not established (beam-sensitive)"**, which is never a pass | measured; evaluation with shared beam lineage (as V6); LV + HR; **report-only for V3-005C** (item 9) | The tolerances are the frozen V6 tolerances, unchanged. Only the point set and the sensitivity semantics differ (item 4) |
| V6-sens-F | Beam sensitivity on H_F with common random numbers: σ_x, σ_θ, ρ_xθ each × 0.8 and × 1.2 (six variants, as Amendment 17 (c)5), k shards per variant plus a CRN nominal | Per point, not global: a point of H_F belongs to H_g iff, for each of the six variants, \|D_var/D_nom − 1\| + 2σ_Δ < 0.15 (the frozen per-point metric). A point where the nominal or a variant mean is zero is beam-sensitive by definition. No row-level pass. The output is H_g and the per-point metrics | –, LV | Item 5 |

## 4. Point set and the documented differences from the frozen row

**Definitions.** These use only the frozen beam inputs and the measured table, never model output:

- z_F = **21.6 cm**: the deepest data depth ≤ R0 + σ_sem = 21.321 + 0.293 = 21.614 cm, with p1 and p2 as
  frozen in Amendment 14 (g).
- P_F,0: every non-missing r = 0 point with z ≤ z_F. That is **29 points** (T).
- H_F: every non-missing r ≥ 1 cm point with z ≤ z_F **and** measured dose ≥ 1e-4 of the measured r = 0
  maximum. The floor is the frozen sizing-subset floor of Amendment 14 (g) / 17 (c)6, not a new number.
  That is **145 points** (T): r = 1, 2, 3, 4, 5, 6 cm with 30 / 29 / 28 / 20 / 21 / 17 points, 34
  depths, and **24** sign-test depths with both groups.
- H_g: the subset of H_F that V6-sens-F finds beam-insensitive.

**Differences from the frozen V6 / V6-sens, each with its reason.** No tolerance is changed. The r = 0
2 %, halo 15 % / 90 % / 10 % and sign-test p > 0.01 criteria are retained unchanged wherever a point gates.

| # | Frozen row | V6-F | Reason |
|---|---|---|---|
| D1 | Every r = 0 point gates (34) | r = 0 points with z ≥ 21.8 cm are excluded from gating and reported: 21.8, 22.0, 22.2, 22.4, 22.6 cm (5 points) | (a) Unsizable or very expensive under the 0.5 % SE rule: 18 / 116 / 651 / 8910 / ∞ shards (T). (b) Beyond R0 + σ_sem the measured r = 0 log-gradient is 56–111 % per mm at 21.8–22.4 cm (central differences of `tbl:dmlg`), so 2 % corresponds to ≤ 0.04 mm of range. The source states positioning errors "of order 1 mm in both z and r" (halo.tex l.321), so these points test positioning, not transport. (c) z ≥ 22.2 cm contains aura dose, which the model does not transport. The frozen annotation (Amendment 14 (g), "Informative annotation") already records the gradient issue. The withdrawal of the 0.95·R0 exclusion (Amendment 14, revision 1) is **respected**: every r = 0 point up to 21.6 cm, including the proximal distal edge 21.0–21.6 cm, gates at 2 % |
| D2 | Every non-missing r ≥ 1 cm point gates (261) | H_F (145). Excluded (116): 42 points beyond z_F (r = 1 / 2 / 3 / 4 / 5 / 6 / 8 / 10 cm: 5 / 9 / 8 / 5 / 5 / 4 / 2 / 4) and 74 points below 1e-4 within z_F (r = 4 / 5 / 6 / 7 / 8 / 10 cm: 3 / 6 / 12 / 23 / 16 / 14). All r = 7, 8 and 10 cm points are excluded | Aura scope (item 2). The 33 points with zero model dose in 1e6 are all in the excluded set. Selection uses the measured value and depth only, so it cannot preferentially drop points where the model is low. A model-count floor could (alternatives.md A3) |
| D3 | V6-sens all-or-nothing: one failing point at one variant makes the whole V6 report-only (Amendment 17 (c)7) | Per-point screen: a beam-sensitive point leaves the gating set and is reported. The halo gates only if the coverage minimum (item 6) holds | A Fermi–Eyges core-only planning estimate (`core_sens_estimate.py`; no ionmc output, no measured dose) gives a core change of −36 % to −45 % at r = 1 cm, z ≤ 10 cm (σ_x × 0.8), and +31 % to +39 % at r = 2 cm, z ≈ 20 cm (σ_θ × 1.2). These points are core-dominated, so under the frozen all-or-nothing rule the halo would probably be report-only from physics alone, even with infinite statistics. The screen keeps the frozen per-point metric and tolerance (0.15, 2σ_Δ). It changes only the consequence: the point leaves the gating set instead of the whole halo becoming report-only. The screen depends on model-vs-model output only, never on residuals |
| D4 | R0 ± 0.5 mm variants run, report-only | Not run | Budget. They were report-only in the frozen row. The r = 0 range-shift equivalent of each r = 0 residual is still reported (Amendment 14 (g)) |
| D5 | History-count rule: rehearsal-estimated SE, no safety factor, N never decreased | Same targets (3 % halo, 0.5 % r = 0) on P_F, with the rehearsal relative SE inflated × 1.25 to allow for the uncertainty of a 20-batch SE estimate | Gives N_F = 15 → declared **16** (item 5) |
| D6 | V6-sens on the frozen point set, k from the frozen point set (32 or 345) | k = 8 per variant on H_F | Item 5 |

## 5. History counts and the SE table that justifies them

All counts are 1e6-history shards, 20 batches each, with SE ∝ N^−1/2 scaling from the rehearsal (T).

| Set | Points | Worst rehearsal relative SE per 1e6 | Shards for the target | With × 1.25 safety |
|---|---|---|---|---|
| P_F,0 (target 0.5 %) | 29 | 1.35 % at z 21.6 cm | 8 | 12 |
| H_F (target 3 %) | 145 | 9.09 % at (1.533 cm, 3 cm); median 1.67 % | 10 | 15 |
| Frozen V6 halo clause (for comparison) | 153 | 28.1 % at (22.0, 3) | 88 | 138 |
| Frozen V6 r = 0 clause (for comparison) | 34 | ∞ at 22.6 cm | ∞ | ∞ |

- **N_F = 16 shards = 1.6e7 histories per backend.** Expected realized worst relative SE: 2.3 % (halo)
  and 0.34 % (r = 0). Every H_F point had all 20 batches non-empty in the rehearsal (T).
- **V6-sens-F, k = 8 shards per run, 7 runs** (nominal + six variants, CRN): 5.6e7 histories.
  - The sizing proxy is the paired σ_Δ of the only rehearsed pair (st-hi vs nominal). Its worst value on
    H_F is 0.0988 per 1e6 at (21.6 cm, 3 cm) (T).
  - Target 2σ_Δ ≤ 0.075, half the 0.15 tolerance, which leaves ≥ 0.075 for \|R − 1\|. That needs k = 7;
    k = 8 is declared, giving a worst 2σ_Δ of 0.070.
  - No safety factor is applied here. An under-estimated σ_Δ only removes points from H_g, and the
    coverage minimum protects the row. This is disclosed.
- **Realized precision is reported, never acted on.** After evidence, no shard is added and no point is
  dropped for precision. A gating point whose realized relative SE exceeds twice its target is
  flagged in the report and still evaluated. There is no sequential top-up (no optional stopping).

## 6. Coverage minimum for the halo component (fixed now)

The halo component gates only if H_g has:

- ≥ **40 points**;
- points at ≥ **3 distinct radii**;
- ≥ **10 sign-test depths** (depths with H_g points in both the 1–3 cm and the ≥ 4 cm groups). With 10
  depths the smallest attainable two-sided p is 0.002 < 0.01, so the sign test can fail. With fewer than
  8 depths it is vacuous.

Otherwise V6-F is "not established (beam-sensitive)". The r = 0 verdict and every halo residual are still
reported.

## 7. Backends, steps and seeds

**Seed bases and scheme.** The bases are lv5c 20481004 and hr5c 20491004 (Amendment 14 (j)). The seed is
base + 1000·r_index + k.

| Step | Suite / backend | r_index | Seeds | Shards |
|---|---|---|---|---|
| `v6f-shard` | lv5c, warp-cpu float64 | **21** | 20502004 + k, k = 0..15 | 16 |
| `v6f-sens-shard` (nominal, sx-lo, sx-hi, st-lo, st-hi, rho-lo, rho-hi) | lv5c, warp-cpu float64 | **22** | 20503004 + k, k = 0..7, the same seeds for all seven runs (CRN, Amendment 17 (c)4) | 56 |
| `v6f-cuda-shard` | hr5c, warp-cuda float32, chunk ≤ 2^14 | **13** | 20504004 + k, k = 0..15 | 16 |
| `v6f-compare`, `v6f-cuda-compare` | consumers, no seed | – | – | – |

**Seed checks.**

- **Declared unused** because of the equalities of the scheme: lv5c r_index 23 (equals hr5c 13), and hr5c
  r_index 11 and 12 (equal lv5c 21 and 22).
- **Disjointness**, checked in the compare steps (fail closed, Amendment 14 (j)):
  - all V6-F seeds (20502004–20504019) lie above the largest consumed seed of `consumed_seed_set()`
    (20482036);
  - they are below the rehearsal family 20505000–20509999;
  - they are disjoint from the reserved original V6 seeds (20493004+, 20494004+, 20501004+);
  - they are disjoint from the fresh-base sequence 20511004, 20521004, …
- lv5c r_index 24 (20505004) would enter the rehearsal family and must never be used.

**V6-sens-F evaluation.** V6-sens-F is evaluated once, on lv5c. The CUDA verdict inherits H_g, as with the
existing `fixed_status` mechanism.

**Disagreement rule.** This is Amendment 17 (c)8 unchanged. V6-F fails if lv5c and hr5c differ on any
component (r = 0, halo, overall, coverage). The per-point paired z is reported, and any |z| > 4 is
listed as a finding. The CUDA rate is **not measured** for this design. **≥ 10× the CPU rate is assumed
as a planning figure only**. Per-step safety does not depend on it, because a 1e6 shard takes ≤ 184 s
even at 1×.

## 8. Host-time budget (warp-cpu float64 at 184 s per 1e6, Amendment 14a (a)2)

| Item | Shards | CPU s | Hours |
|---|---|---|---|
| V6-F gating | 16 | 2 944 | 0.82 |
| V6-sens-F (7 × 8) | 56 | 10 304 | 2.86 |
| **CPU total** | **72** | **13 248** | **3.68** (≤ 4 h) |
| V6-F CUDA (planning, 10×) | 16 | ≈ 294 | 0.08 (worst case 1×: 0.82) |
| Possible V7-R rehearsal re-run (item 11, risk R1) | 2 | ≈ 1 600 | 0.45 |

- Each step is one 1e6 shard: 184 s, against the 1500 s step limit and the 1200 s effective limit of
  Amendment 14a (c)5. The throughput step at the evidence head must list the V6-F shard plans as fitting
  (`shard_table_ok`).
- Total host time with the CUDA planning figure is ≈ 3.8 h, and ≈ 4.5 h at CUDA 1×. Both are under 5 %
  of the 88 + 288-shard plan that IR-20261011-031924-C205C8 put at 20–41 h.

## 9. Merge-gate status (recommendation)

V6-F is **report-only for V3-005C**: its verdict is computed and archived but is not part of (a)6. It is
recommended as **gating for the later task that makes the lateral-halo claim** (decision 0041 l.201).

Reasons:

1. The operator directive keeps V6/V6-sens failed and the applicable gates enforced. Making a redesigned
   row gate the same task would in effect replace the failed row, which item 1.2 excludes.
2. V6-F was designed after the rehearsal. The design was made after the contrary-evidence event of
   Amendment 14a (a)5: measured and model mean doses of the 1e6 rehearsal shard were displayed side by
   side to the orchestrator. Mitigations:
   - every selection rule is stated in frozen quantities (p1, p2, the frozen 1e-4 floor);
   - every rule is independently reviewed;
   - the contrary readings of item 10 are reported.

   Even so, its blindness is weaker than a row frozen before any model output. A gate in the task where
   that event happened invites the reading that the scope was chosen to pass.
3. V3-005C has its own gating evidence for the elastic channel (V10-A, V11, V5). The halo is not needed
   to decide the V3-005C physics question.
4. In the later task, V6-F gates exactly as pre-registered here. It is run on that task's declared fresh
   base at its own evidence SHA, without changes to these criteria.
   - Any V3-005C V6-F result, pass or fail, is reported next to the later one, never in its place
     (Amendment 14 (a)2).
   - A V3-005C failure stays recorded. The later execution is the single pre-declared gating execution,
     not a retry.

If the orchestrator prefers not to spend the 3.7 h in V3-005C, the alternative is that V6-F runs only once,
in the later task, as gating. The design is unchanged in that case.

## 10. Reporting (no verdict effect)

The following are reported next to the verdict and never change it:

- Residuals with 90 % intervals at all 295 present points, and the range-shift equivalent at r = 0.
- The ±3 % monitor normalisation, reported and not added.
- Contrary readings:
  - (a) the halo verdict on H_F without screening (frozen all-or-nothing semantics);
  - (b) the V6-F verdict with z_F = 21.8 cm, flagged as under-sized at its deepest points;
  - (c) the V6-sens-F screen evaluated with 1σ_Δ instead of 2σ_Δ.
- The r = 6 cm scan, which the source names as an outlier (l.321), is **retained** and not post-hoc
  excluded. Its share of halo failures is reported.

## 11. Implementation conditions (before any V6-F evidence)

- New steps in `steps_v6.py` or a new module, the point-set masks of item 4, the screen and coverage of
  item 6, and the seed guards of item 7. These go through Codex review.
- The existing `evaluate`, `_sens_status` and `combine_backends` stay unchanged for the frozen rows.
  V6-F uses new functions.
- Adding these files changes the lv5c hashed-source digest. Under Amendment 14a (b)2 that may require
  re-running the two V7-R rehearsal shards (≈ 0.45 h) before the C9 evidence head. This is risk R1.

## 12. What V6-F can and cannot establish

**Can establish:**

- With the frozen tolerances, absolute dose per primary (no fitted scale) agrees at r = 0 for z ≤ 21.6 cm,
  including the proximal distal edge.
- The charged-secondary halo (core-insensitive points, r = 1–6 cm, z ≤ 21.6 cm, dose ≥ 1e-4 of the
  maximum) agrees within 15 % per point, ≥ 90 % pass, median ≤ 10 %, with no radial trend.
- This is a test of non-elastic plus elastic nuclear secondaries and of MCS tails, made where calibration
  inputs cannot move the model by the tolerance.

**Cannot establish:**

- The aura or any neutral-particle dose.
- The distal tail beyond R0 + σ_sem.
- The halo at r ≥ 7 cm or below 1e-4 of the maximum.
- The core-dominated points (screened out).
- Independence from the beam fit: the r = 0 points and the beam stay shared lineage, and E₀ comes from
  R0 = d80.
- Absolute normalisation better than the source's ±3 %.
- Any energy other than 177 MeV.
- The full frozen V6 claim. That stays failed.
