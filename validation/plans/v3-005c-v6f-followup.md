# V3-005C row V6-F — follow-up validation design (pre-registration, revision 2) and its independent review record

Governing summary: validation/plans/v3-005-acceptance.md, Amendment 18 (b). This file governs in any difference. Scripts of the design (se_table.py, core_sens_estimate.py) and the SE table are archived with the C9 results; revision-1 per-point tables are quarantined and never committed.

---

# V6-F — pre-registration text (Amendment 18, V6-F items), revision 2 — FROZEN 2026-10-11 at the commit that adds this file

**Status.**

- Revision 2 answers the independent review `review.md` and the orchestrator's decisions on it. They are
  listed in "Changes after review" at the end. Revision 1 is kept in `v1/`.
- Status: frozen with Amendment 18 (this file is committed; the row has not been executed and is not executed in V3-005C). The wording "draft" inside the review record below is historical (it reviewed revisions 1 and 2 before this freeze) and is kept verbatim.
- Numbers marked (T) are from `se_table.json` / `se_table.py`, revision 2. The planning estimates (E)
  are from `core_sens_estimate.py`.
- Neither script reads any ionmc dose output other than precision quantities of the three C8 rehearsal
  partials. The rehearsal seeds are 20517000 and 20525000 at rehearsal base 20505000, head c847cc75,
  and are never evidence.

---

**Amendment 18 — V6-F items (frozen 2026-10-11; after the C8 rehearsal and Amendment 14a (a)4; before any
V6, V6-F or V6-sens-F evidence).** These items do not change:

- the frozen rows V6 and V6-sens;
- Amendments 14 (g), 17 (c) or 14a (a);
- any tolerance, seed base or verdict recorded there.

They pre-register two new rows, **V6-F** and **V6-sens-F**, with a documented list of scope and criteria
differences (item 4). They answer the operator response to IR-20261011-031924-C205C8, recorded verbatim
in its `response.json`. V6-F is **pre-registered now and executed once, as a gating row, in a later task**
(item 8). It is not executed in V3-005C.

## 1. The original failure stays recorded

1. Rows **V6** and **V6-sens** stay recorded as **failed — not sized under the pre-registered rule**
   (Amendment 14a (a)4).
   - `V6_SHARDS` stays `None` and the steps `v6-unsized` keep their verdict.
   - The original V6 seeds are lv5c r_index 12 and 20 and hr5c r_index 3, of bases 20481004 and
     20491004. They were never run, and no V6-F seed reuses them (item 7).
2. V6-F never replaces V6 in the requirement ledger, the results README, the decision 0041 outcome or the
   merge-gate list of Amendment 14 (a)6. P6-D (Amendment 15 (a)5) is not a precedent here.
   - Every report that shows a V6-F verdict shows, next to it: "V6/V6-sens: failed — not sized under the
     pre-registered rule (Amendment 14a (a)4)".
3. **Merge gate (a)6.** The treatment of the "failed — not sized" V6 under the Amendment 14 (a)6 merge
   gate of V3-005C is fixed by the orchestrator in a separate item of Amendment 18. These V6-F items do
   not rule on it.
4. **Claims.** A V6-F pass is reported as "V6-F (reduced scope, item 4) passed", never as "V6 passed" or
   "halo validated".

## 2. Recorded facts behind the redesign

- **The frozen rule could not be sized** (Amendment 14a (a)2–4, reproduced in revision 1 of (T)). The
  rule needs 0.5 % SE at every r = 0 point and 3 % at every r ≥ 1 cm point ≥ 1e-4 of the maximum.
  - At r = 0 it needs 116 / 651 / 8910 shards at z = 22.0 / 22.2 / 22.4 cm. At 22.6 cm the model dose
    was zero in 1e6 histories, so no shard count suffices.
  - The halo needs 88 shards.
  - 32 halo points had zero model dose in both V6-sens rehearsal shards.
- **The transport deposits no neutral energy.** Neutrons and γ are scored as escaped neutral energy only.
  The measured dose contains the **aura**: neutral secondaries and the charged particles they set in
  motion (halo.tex l.45, l.74, l.250).
- **Gottschalk's MD fit gives an empirical aura term** (l.797–798, `tbl:parms` p23 = 0.0008981,
  p24 = 0.8092, p25 = −0.1431):

  ```
  D_au = p23 · E_2D(r_n, r_au),   r_au = p24 · (1 + p25 · r_n),   r_n = r / (R0/3),
  E_2D(r, r_e) = exp(−r / r_e) / (2π r_e²)
  ```

  - D_au is independent of z: 1.90 / 1.65 / 1.42 / 1.21 / 1.03 / 0.86 × 10⁻⁴ MeV/g/p at r = 1–6 cm (T).
  - Inside the revision-1 halo set H_F the aura share D_au/D_meas is not negligible at r ≥ 4 cm.
    Median / max (T):

    | r | median | max |
    |---|---|---|
    | 3 cm | 0.9 % | 7.4 % |
    | 4 cm | 3.2 % | 16 % |
    | 5 cm | 6.0 % | 15 % |
    | 6 cm | 10.1 % | 13.5 % |

  - So a correct charged-halo model reads low there by physics this transport omits. That biases the
    per-point test and the r-trend sign test (review B1).
- **The core dominates the core/halo transition.** A planning estimate (E: core-only Fermi–Eyges,
  Gottschalk's T_dM with X_S = 46.88 cm, averaged over the 0.4 cm sphere) gives:
  - at r = 1 cm, z ≤ 10 cm: −33 % to −42 % (σ_x × 0.8);
  - at r = 2 cm, z = 20–21 cm: +43 % to +49 % (σ_θ × 1.2).

  Under the frozen all-or-nothing V6-sens rule the halo would therefore probably be report-only from
  physics alone. The bound |R_total − 1| ≤ |R_core − 1| holds only if the halo itself is
  beam-insensitive.

## 3. Rows

| # | Check | Pass iff | Class / tier | Justification |
|---|---|---|---|---|
| V6-F | Gottschalk 2015 `tbl:dmlg` (177 MeV; the Amendment 14 (g) beam; E₀ = 177.84075737 MeV of Amendment 14a (a)1; the observable of Amendment 17 (c)1) on the point sets P_F,0 and H_F of item 4. Two components: r = 0 and halo | (i) r = 0: \|D/D_meas − 1\| ≤ 0.02 at **every** point of P_F,0. (ii) Halo, on the beam-insensitive subset H_g ⊆ H_F (V6-sens-F), and only if the coverage minimum of item 6 holds: per point 0.15; ≥ 90 % of points pass; median \|residual\| ≤ 0.10; no r-trend (two-sided binomial sign test over depths with H_g points in both the 1–3 cm and the ≥ 4 cm groups, p > 0.01; ties excluded and counted, Amendment 17 (c)3). (iii) Backend agreement, Amendment 17 (c)8 rule (item 7). V6-F passes iff (i), (ii) and (iii) all hold. If the coverage minimum fails, V6-F is **"not established (beam-sensitive)"**, never a pass | measured; evaluation with shared beam lineage (as V6); LV + HR; **gating in the executing later task** (item 8) | The tolerances are the frozen V6 tolerances, unchanged. The differences are in item 4 |
| V6-sens-F | Beam sensitivity on H_F, with common random numbers. σ_x, σ_θ and ρ_xθ each × 0.8 and × 1.2 (six variants, Amendment 17 (c)5), k shards each, plus a CRN nominal | Per point. A point of H_F belongs to H_g iff, for each of the six variants, **\|D_var − D_nom\| / D_meas + 2σ′ / D_meas < 0.15**. σ′ is the paired batch SE of the difference D_var − D_nom. With the measured denominator the metric is defined even where a model mean is zero, so no zero-mean rule applies (revision 2, R2-m2). There is no row-level pass; the output is H_g and the per-point metrics | –, LV | Item 4 D3, item 5 |

## 4. Point sets and the documented differences from the frozen row

**Definitions.** They use only frozen beam inputs, the measured table and the MD-fit aura term, never
model output.

- **z_F = 21.6 cm**: the deepest data depth ≤ R0 + σ_sem = 21.321 + 0.2925 = 21.614 cm (p1 and p2 as
  frozen).
- **P_F,0**: every non-missing r = 0 point with z ≤ z_F. That is **29 points** (T).
- **H_F**: every non-missing r ≥ 1 cm point with z ≤ z_F, measured dose ≥ 1e-4 of the measured r = 0
  maximum (5.95 MeV/g/p at z = 21.0 cm), **and** aura share D_au/D_meas ≤ 0.05.
  - That is **104 points** (T): r = 1 / 2 / 3 / 4 / 5 cm with 30 / 29 / 26 / 16 / 3 points; none at
    r = 6 cm.
  - There are **16 sign-test depths** with points in both groups before screening.
- **H_g**: the subset of H_F that V6-sens-F finds beam-insensitive.

**Differences.** No tolerance is changed. The r = 0 2 %, halo 15 % / 90 % / 10 % and sign-test p > 0.01
criteria are retained wherever a point gates.

| # | Frozen row | V6-F | Reason |
|---|---|---|---|
| D1 | Every r = 0 point gates (34) | r = 0 points at z = 21.8, 22.0, 22.2, 22.4, 22.6 cm are excluded from gating and reported | See "D1 in detail" below |
| D2 | Every non-missing r ≥ 1 cm point gates (261) | Points beyond z_F (42) and points below 1e-4 within z_F (74, including every r = 7, 8, 10 cm point) are excluded and reported | Aura-dominated or unsizable (item 2). All 33 halo points with zero model dose in 1e6 are among them. Selection uses measured value and depth only. A model-count floor would preferentially drop points where the model is low (alternatives A3) |
| D3 | V6-sens all-or-nothing on the model ratio \|D_var/D_nom − 1\| + 2σ_Δ < 0.15 (Amendment 17 (c)7) | Per-point screen with a **measured denominator**: \|D_var − D_nom\|/D_meas + 2σ′/D_meas < 0.15. A beam-sensitive point leaves the gating set and is reported. The halo gates only under the coverage minimum (item 6) | See "D3 in detail" below |
| D4 | R0 ± 0.5 mm variants run, report-only | Not run | Planning budget. They were report-only. The r = 0 range-shift equivalents stay reported |
| D5 | History-count rule: rehearsal SE, no safety factor | Same targets on P_F,0 and H_F, rehearsal relative SE × 1.25 | Rule value 12 → declared **N_F = 16** (item 5) |
| D6 | V6-sens k from the frozen point set (32 or 345) | k = 8 per variant on H_F | Item 5 |
| **D7** | Every non-missing r ≥ 1 cm point gates; the aura is not considered | **Aura-share exclusion.** Points of H_F with D_au/D_meas > 0.05 are excluded and reported. That is 41 points: r = 3 / 4 / 5 / 6 cm with 2 / 4 / 18 / 17. D_au uses `tbl:parms` p23–p25 and the measured table only | See "D7 in detail" below |
| **D8** | No coverage minimum (all-or-nothing) | Coverage minimum of item 6: ≥ 40 points, ≥ 3 radii, ≥ 12 sign-test depths after tie exclusion | It stops the per-point screen from shrinking the halo test to triviality. It is a new criterion, declared before evidence |

**D1 in detail.**

- (a) The 0.5 % SE rule needs 18 / 116 / 651 / 8910 / ∞ shards at these depths (T).
- (b) The measured r = 0 log-gradient is 56–111 % per mm at 21.8–22.4 cm, so 2 % is ≤ 0.04 mm of range,
  against the ~1 mm positioning error that the source states (halo.tex l.321).
- (c) At z ≥ 22.2 cm the dose contains aura.
- **Choice of z_F, disclosed (review I3).**
  - "+1 σ_sem" is a convention. The gradient argument also applies at 21.4–21.6 cm, and the aura argument
    only from about 22.0–22.2 cm.
  - **Budget discriminates 21.6 from 21.8 cm.** At 21.8 cm the r = 0 point alone raises the N rule from 12 to
    27 shards. The total becomes 27 + 7 × 8 = 83 shards, 4.24 h, above the 4 h target (alternatives A2).
    This is a weak discriminator, and the text says so.
  - The influence of the Amendment 14a (a)5 display, where measured and model mean doses of the
    rehearsal shard were shown side by side, **cannot be excluded** by inspection.
  - Against a favourable choice: the cut keeps the high-risk distal points 21.2 / 21.4 / 21.6 cm. Their
    measured log-gradients are 11.7 / 23.7 / 38.3 % per mm, so 2 % is 0.17 / 0.08 / 0.05 mm of range. A
    CSDA-vs-projected-range offset of a few tenths of a millimetre fails them.
  - The 21.8 cm verdict is reported as contrary reading (b).
  - Every r = 0 point ≤ 21.6 cm gates at 2 % (Amendment 14 revision 1 respected).

**D3 in detail.**

- Review B2: the model-normalised ratio has the model's own total dose in the denominator. A point where
  the model under-predicts the halo has a larger model core fraction, so it would be preferentially
  screened out.
- The measured denominator removes that dependence from the **numerator's normalisation**.
- **The revision-1 claim that the screen "cannot preferentially drop points where the model is low" is
  withdrawn.**
- **Residual dependence, stated.** σ′ scales with the model's own fluctuations, roughly √D_model for a
  given measured dose. A model that is low at a point therefore has a smaller σ′/D_meas, which makes it
  slightly *more* likely to stay in H_g. So the dependence runs opposite to the revision-1 bias, and it
  is second order (only the 2σ′ term).
- The frozen-metric screen (\|D_var/D_nom − 1\| + 2σ_Δ < 0.15) is reported as contrary reading (d).
- The core estimate (item 2) explains why some screen is needed: the core/halo transition is beam-
  sensitive at ±20 %.

**D7 in detail.**

- The transport deposits no neutral energy. Where the aura exceeds 5 % of the measurement, the
  comparison tests the missing physics, not the charged halo, and it biases the ≥ 4 cm group of the sign
  test toward negative residuals.
- 5 % is one third of the 15 % per-point tolerance.
- **Uncertainty of the classification.** D_au is a crude, depth-independent empirical term of an MD fit
  whose global rms ratio of measurement to fit is 15 % (halo.tex l.49). A point near the 5 % boundary may
  therefore be classified either way, with that uncertainty. The per-point aura share and contrary
  reading (a) make such points visible.
- **Role extension, disclosed.** `tbl:parms` was frozen as "construction for the beam only"
  (Amendment 14 (g), 17 (c)10). p23–p25 now also serve as a **construction input of the V6-F point
  selection**. They come from the same 25-parameter fit to `tbl:dmlg`, so this is shared lineage.
- The executing task records the role extension in the role record of the arXiv:1409.1938 acquisition
  and in the registry description before execution.
- The unexcluded verdict is contrary reading (a).

## 5. History counts and the SE table

All counts are 1e6-history shards, 20 batches each, with SE ∝ N^−1/2 from the rehearsal (T).

| Set | Points | Worst rehearsal relative SE per 1e6 | Shards for target | × 1.25 safety |
|---|---|---|---|---|
| P_F,0 (target 0.5 %) | 29 | 1.35 % at z 21.6 cm | 8 | 12 |
| H_F (target 3 %) | 104 | at (21.6 cm, 3 cm), 6 shards | 6 | 9 |
| Revision-1 H_F (no D7) | 145 | 9.09 % at (1.533, 3) | 10 | 15 |
| Frozen V6 halo clause | 153 | 28.1 % at (22.0, 3) | 88 | 138 |

**N_F = 16 shards = 1.6e7 histories per backend (unchanged from revision 1).**

- The rule value is 12. 16 is kept for margin. Precision never gates, and the count is not reduced once
  declared.
- Expected realized worst relative SE: 1.74 % (halo) and 0.34 % (r = 0).

**V6-sens-F: k = 8 shards per run, seven runs (unchanged).**

- Proxy: the paired SE of D_st-hi − D_nom over D_meas, the only rehearsed pair. Its worst value on H_F is
  0.0369 per 1e6 at (20.4 cm, 5 cm) (T).
- The sizing target is 2σ′/D_meas ≤ 0.075 (half of 0.15). The rule value is k = 1. k = 8 is kept, giving
  a worst 2σ′/D_meas of 0.026 and leaving ≥ 0.12 of the tolerance for the sensitivity itself.
- The σ_x and ρ variants may have a larger σ′ than st-hi. That would only shrink H_g, and the coverage
  minimum then decides.

**Precision is reported, never acted on.** No shard is added after evidence, and no point is dropped for
precision. A gating point whose realized relative SE exceeds twice its target is flagged and still
evaluated.

## 6. Coverage minimum and sign-test power (D8)

The halo component gates only if H_g has:

- ≥ **40 points**;
- points at ≥ **3 distinct radii**;
- ≥ **12 sign-test depths, counted after excluding ties**.

Otherwise V6-F is "not established (beam-sensitive)".

**This is a realistic outcome, not a remote one.** The coverage expectation below assumes:

- a beam-insensitive halo (only the core responds to the beam variants);
- the σ′ of one variant pair (st-hi) for all six variants.

Either assumption can fail. The ≥ 4 cm group of the sign test will consist mostly of r = 4 cm points
(16 of the 19 points at r ≥ 4 cm in H_F).

**Sign-test power (exact two-sided binomial, p > 0.01 passes).**

| Depths | Fails at | p at the boundary | P(fail) if a fraction q of depths share one sign (q = 0.7 / 0.8 / 0.9) |
|---|---|---|---|
| 12 | ≥ 11/12 same sign (p = 0.0063); 10/12 passes (p = 0.039) | | 0.085 / 0.275 / 0.659 |
| 16 | ≥ 14/16 (p = 0.0042); 13/16 passes (p = 0.021) | | 0.099 / 0.352 / 0.789 |

- At the floor of 12 depths the test detects only a strong, consistent radial trend.
- Correction to the orchestrator brief: at 12 depths not only 12/12 but also 11/12 same-sign fails.

**Coverage expectation (E, planning only).** This assumes the halo is beam-insensitive, and uses the core
scaled to the measured r = 0 dose, the measured denominator and σ′/√8.

- Expected H_g ≈ **58 points**:
  - r = 1 / 2 / 3 / 4 / 5 cm with 17 / 0 / 22 / 16 / 3 points;
  - r = 1 cm kept only at z = 14–21.6 cm;
  - r = 2 cm screened out entirely.
- That is 4 radii and **16 sign-test depths**.
- The minimum is expected to be reachable. The margin is 58 vs 40 points and 16 vs 12 depths, before
  ties.
- The 1–3 cm group is then deep r = 1 cm plus r = 3 cm, and the ≥ 4 cm group is mostly r = 4 cm.
- Revision 1's expectation (≈ 95 points, ≈ 20 depths) referred to H_F before D7. It is superseded.

## 7. Backends and seeds for the single execution

**Backends.**

- warp-cpu float64 (the LV gating sample and V6-sens-F) and warp-cuda float32 (HR gating replicate, chunk
  ≤ 2^14).
- V6-sens-F is evaluated once, on warp-cpu f64. The CUDA verdict uses the same H_g.
- Disagreement rule: Amendment 17 (c)8 unchanged. V6-F fails if the two backends differ on any component
  (r = 0, halo, overall, coverage). The per-point paired z is reported, and any \|z\| > 4 is listed as a
  finding.

**Seed scheme, declared now.** The executing task declares its fresh base B in an amendment before any
observation. B is the next unused value of the 10000-spaced sequence of Amendment 14 (j) (…, 20511004,
20521004, 20531004, …) **that passes checks (a)–(e) below**. A base skipped because it fails a check is
recorded in that amendment and is never used afterwards. B itself is never reused. The seeds are B + 1000·r_index + k:

| r_index | Sample | Seeds |
|---|---|---|
| **7** | V6-F gating, warp-cpu f64 | B + 7000 + k, k = 0..15 |
| **8** | V6-sens-F: nominal, sx-lo, sx-hi, st-lo, st-hi, rho-lo, rho-hi | B + 8000 + k, k = 0..7, identical for all seven runs (CRN) |
| **9** | V6-F gating replicate, warp-cuda f32 | B + 9000 + k, k = 0..15 |

r_index ≤ 9 and k ≤ 999 keep every seed below B + 10000. That avoids the cross-base collision of the
10000-spaced sequence at r_index ≥ 10 (for example 20511004 + 21000 = 20531004 + 1000).

**Required checks.** These are written in the executing amendment and fail closed in the compare steps.
Every V6-F seed must be disjoint from:

- (a) `consumed_seed_set()` and every evidence seed used up to then;
- (b) **every rehearsal seed actually used**, not the base range. The C8 rehearsal seeds are
  20505000 + 1000·r + k, with r up to 20 (20517000 and 20525000 for V6). They reach 20525999, so
  B = 20511004 would collide (B + 7000 = 20518004). The check therefore requires B + 7000 > the
  largest rehearsal seed used, or an explicit seed-by-seed disjointness proof;
- (c) the never-run original V6 seeds 20493004 + k, 20494004 + k and 20501004 + k, which are **not
  reused**;
- (d) every other row of the executing task;
- (e) one another: the CPU and CUDA samples are disjoint, so the paired z is valid.

## 8. Execution and gating status

1. **V3-005C.** No V6-F step, seed or code is added to V3-005C. No V6-F code enters the lv5c/hr5c hashed
   source set (`run_suite.SOURCE_FILES_V5C`), and no V6-F host time is spent there. The V7-R rehearsal
   fixture of Amendment 14a (b)2 is unaffected.
2. **The later task** that makes the lateral-halo claim (decision 0041 l.201) executes V6-F **exactly
   once**, as a **gating** row.
   - It runs at its own evidence SHA, on its fresh base B, with the criteria, sets and counts of these
     items unchanged.
   - A failure is a recorded scientific failure. There is no second execution on another base without a
     new, reviewed physics reason and a new pre-registration.
3. **Planning budget for that task** (warp-cpu f64 at 184 s per 1e6, Amendment 14a (a)2): 16 + 7 × 8 =
   72 shards, 13 248 s, 3.68 h. Each step is one 1e6 shard of 184 s against the 1200 s effective limit.
   - The executing task's own throughput step re-measures the rates, including CUDA, and must report
     every V6-F shard plan as fitting (`shard_table_ok`).
   - No CUDA rate is assumed here.

## 9. Reporting (verdict-neutral)

Reported next to the verdict and never changing it:

- The residual and 90 % interval at all 295 present points, and the range-shift equivalent at r = 0.
- **The aura share D_au/D_meas at every point.**
- The ±3 % monitor normalisation, reported and not added.
- The r = 6 cm scan is reported as the source's own outlier (l.321). It is not part of H_F after D7.

Contrary readings:

- **(a)** Halo verdict on the 145-point set without D7, with the bias declared beforehand: residuals
  expected low by about the aura share (median 3 / 6 / 10 % at r = 4 / 5 / 6 cm), and the sign test
  pushed toward "≥ 4 cm lower".
- **(b)** V6-F with z_F = 21.8 cm, flagged as under-sized at its deepest points.
- **(c)** The screen with 1σ′ instead of 2σ′.
- **(d)** The screen with the frozen model-normalised metric \|D_var/D_nom − 1\| + 2σ_Δ < 0.15.
- **(e)** The halo verdict on H_F \ H_g, the screened-out points.

## 10. What V6-F can and cannot establish

**Can establish, given the fitted beam.** The beam and E₀ are construction inputs from the same fit; the
r = 0 dose is itself strongly beam-sensitive (E: +57 % on axis at z = 1.5 cm for σ_x × 0.8).

- Absolute dose per primary (no fitted scale) agrees at r = 0 within 2 % for z ≤ 21.6 cm, including the
  proximal distal edge.
- The **charged-secondary halo** agrees within 15 % per point, ≥ 90 % pass, median ≤ 10 %, with no
  radial trend. This holds at core-insensitive points with r = 1–5 cm, z ≤ 21.6 cm, dose ≥ 1e-4 of the
  maximum and aura share ≤ 5 %.
- That is a test of non-elastic plus elastic nuclear secondaries and of MCS tails, where calibration
  inputs cannot move the model by the tolerance and the omitted neutral dose is small.

**Cannot establish:**

- The aura or any neutral dose.
- The halo at r ≥ 6 cm, or wherever the aura share exceeds 5 %.
- The distal tail beyond R0 + σ_sem.
- The core-dominated points (screened out, expected all of r = 2 cm and shallow r = 1 cm).
- Independence from the beam fit, which is shared lineage for the beam, E₀ (from R0 = d80) and now the
  aura selection (p23–p25).
- Normalisation better than ±3 %.
- Any other energy.
- The frozen V6 claim, which stays failed.

---

## Changes after review (revision 1 → revision 2)

| # | Review item | Change | Reason |
|---|---|---|---|
| 1 | B1 → option (i) | New difference **D7**: aura-share exclusion at 5 %. D_au from `tbl:parms` p23–p25 and the measured table only. H_F 145 → **104** points; sign-test depths 24 → **16**. Per-point aura share reported everywhere. The 145-point verdict with the declared bias becomes contrary reading (a). Role extension of `tbl:parms` stated, to be recorded by the executing task | No neutral energy is deposited, so above 5 % aura share the comparison tests the missing physics, not the charged halo. Counts recomputed: 104 points and 16 depths, matching the reviewer |
| 2 | B2 → option (a) | Screen metric \|D_var − D_nom\|/D_meas + 2σ′/D_meas < 0.15, with σ′ the paired SE of the difference. The frozen model-normalised screen becomes contrary reading (d). The revision-1 "cannot preferentially drop" claim for D3 is withdrawn, and the residual σ′ dependence is stated | Selection must not depend on the model's level through the denominator |
| 3 | B3 | Every V3-005C execution removed: revision-1 sections 7–9, 11, risk R1, the CUDA figures and the 20481004/20491004 seeds. V6-F is executed once, gating, in a later task. Seed scheme declared: base B of that task; r_index 7 / 8 / 9 for gating / sensitivity / CUDA. Disjointness checks against rehearsal seeds actually used. No original V6 seed reused | The orchestrator's disposition. Avoids a two-look history |
| 4 | I1 | Item 1.3: the (a)6 treatment of "failed — not sized" is fixed in a separate item of Amendment 18 by the orchestrator. V6-F does not rule on it. The revision-1 item 9 recommendation is removed | Scope separation |
| 5 | I2 | Coverage depths ≥ 12 (was 10), counted after ties. Exact sign-test boundaries and power at 12 and 16 depths stated. The coverage minimum is listed as difference **D8** | Power. Correction: at 12 depths, 11/12 same-sign also fails (p = 0.0063), not only 12/12 |
| 6 | I3 | D1 now states: +1 σ_sem is a convention; budget is the (weak) discriminator, 27 vs 12 shards and 4.24 vs 3.68 h; the (a)5 display cannot be excluded as an influence; the high-risk 21.2–21.6 cm points are kept; 21.8 cm is contrary reading (b) | Disclosure |
| 7 | I4 | Coverage expectation added (item 6), recomputed after D7 and with the measured-denominator screen: ≈ 58 points, 4 radii, 16 depths, r = 2 cm screened out. Distal-edge range risk 0.05–0.17 mm added to D1 | The revision-1 expectation of ≈ 95 points, ≈ 20 depths applied to the pre-D7 set |
| 8 | m1 | `se_table.py`: constant `SIGMA_SCREEN_F = 0.0375`; column `sens_k_for_2sigma_0p075`; the 0.025 columns removed | Consistency with the design |
| 9 | m2 | `se_table.py`: `ratio_paired` is no longer called. New helpers `rel_se_only` and `paired_diff_se` return no mean and no ratio. The model-normalised σ_Δ column was removed, because next to the measured-normalised σ′ it would reveal D_nom/D_meas. The revision-1 table was deleted from this directory (`v1/` keeps only the scripts and text) | Blindness by construction |
| 10 | m3 | `core_sens_estimate.py`: Gottschalk T_dM with f_dM, E_s = 15 MeV, X_S = 46.88 cm; sphere averaging; `np.trapezoid` with fallback; the "upper bound" qualified (it holds only if the halo is beam-insensitive). New part 2 gives the coverage expectation | Physics accuracy, NumPy 2 |
| 11 | m4 | Item 10 qualifies the r = 0 statement as "given the fitted beam" | r = 0 is beam-sensitive (+57 % on axis, σ_x × 0.8) |
| 12 | m5 | The rehearsal is named by its **base** (20505000) and its **seeds actually used** (20517000, 20525000); the disjointness rule uses used seeds | Precision |
| 13 | m6 | CUDA rate figures dropped. The executing task measures them | Moot without V3-005C execution |
| 14 | — | N_F = 16 and k = 8 unchanged. Rule values are now 12 and 1. Planning CPU time 3.68 h | Counts are not reduced once proposed |
| 15 | R2-I1 | The revision-1 script `v1/se_table.py` is neutralised: its per-point precision columns and the model-normalised σ_Δ computation are removed. No revision-1 table exists in this directory. **Revision-1 per-point tables (including any reviewer copy) are quarantined: they must never be committed, or archived next to the revision-2 table** | Combined with revision 2 they would give a per-point model-to-measured map |
| 16 | R2-m1 | The implied ratio is recorded as the bound ≤ 0.37 at (21.6 cm, 3 cm) only | Precision of the disclosure |
| 17 | R2-m2 | The rule that a zero nominal or variant mean makes a point beam-sensitive is **dropped** | The measured denominator makes the metric finite there. Keeping the rule would remove exactly the points where the model dose is zero, the most extreme model-low points, which is selection by agreement. Without it such a point stays in H_g and fails its residual test openly |
| 18 | R2-m3 | B is the next unused base passing checks (a)–(e). A skipped base is recorded and never reused | Fully determined base choice |
| 19 | R2-m4 | Item 6: "not established (beam-sensitive)" is stated as a realistic outcome, with the two assumptions of the expectation. The ≥ 4 cm group is mostly r = 4 cm | Honest expectation |
| 20 | R2-m5 | D7: the aura term is a crude depth-independent term from a 15 %-rms fit, so classification near 5 % carries that uncertainty | Disclosure |

**Contrary-evidence note (revision 2, disclosed by the design author).** The model-normalised screen SE
of revision 1 (worst 0.0988 per 1e6 at (21.6 cm, 3 cm), reported in the revision-1 handback) and the
measured-normalised σ′ of revision 2 (≤ 0.037 at every H_F point) were both seen by the design author.

- Their quotient approximates D_nom/D_meas at that point, up to the CRN linearisation. It bounds the 1e6
  rehearsal ratio at **(21.6 cm, 3 cm) only**: D_nom/D_meas ≲ 0.037 / 0.0988 = **0.37**.
- No other point's ratio is implied or known. Only the maximum of σ′ over H_F was seen, not the
  per-point pairing.
- That point is in the gating set H_F. No rule was changed after this observation: D7, the D3 metric,
  z_F and the coverage numbers were fixed by the orchestrator's decisions before the revision-2 script
  was run.
- It is recorded here in the manner of Amendment 14a (a)5. The reviewer's copy of the revision-1 table
  must not be combined with the revision-2 table.


## Derivation artifacts (committed with this file; sha256 at the freeze)

The numerical derivation of the point set, N_F, k and the coverage expectation is reproducible from these files (`python3 -I se_table.py <worktree>` regenerates `se_table.json` byte for byte from the three rehearsal partials named in the design; `python3 -I core_sens_estimate.py` reproduces the coverage estimate). `v1/` holds the revision-1 text and checklist and the NEUTRALISED revision-1 script (its per-point precision columns removed); no revision-1 per-point table exists or may be committed. The operator's intervention record is `reference_requests/interventions/IR-20261011-031924-C205C8/` (request.json, response.json).

| file | bytes | sha256 |
|---|---|---|
| `validation/plans/v3-005c-v6f/core_sens_estimate.py` | 5253 | `4105071f00555a84910d3b9b62de0bb7ded3f753a3bf392cde09b2ebf938324d` |
| `validation/plans/v3-005c-v6f/review_checklist.md` | 7183 | `f4a603e0002dd1ec7030ec2312fa4d5ed32444cf2fd3c33bd6b381d96aeed3a2` |
| `validation/plans/v3-005c-v6f/se_table.json` | 261735 | `b929408aed193dbd31f459e654f09c5a4afb90f5454258d131a1771f8645cd0b` |
| `validation/plans/v3-005c-v6f/se_table.py` | 14319 | `ba4d0e34f75661bd743a481e418bdfa0ea453457970ae28fd224db96196d8869` |
| `validation/plans/v3-005c-v6f/v1/design.md` | 18496 | `78408160e26ac3835e8634fbd926bf70a604fb267de1c45249d4c285eb81b439` |
| `validation/plans/v3-005c-v6f/v1/review_checklist.md` | 7251 | `577f743c1df96e873d87ec5c86d7604b6c7cfcc403611d700ee6e52ff512f713` |
| `validation/plans/v3-005c-v6f/v1/se_table.py` | 13930 | `e33a7a06a8e61c90324e7d2d423fa46ebd4f7d45e005ebfe64dc5958e9fdeaae` |

---

# Alternatives considered

# V6-F — alternatives considered (revision 2)

Revision 2 changes A1, A2, A5 and A8 and adds A10 and A11 (see design.md, "Changes after review").

All numbers come from `se_table.json`, and all costs are in 1e6-history warp-cpu f64 shards at 184 s each.
The CPU budget is 4 h, which is 78 shards in total. The chosen design (design.md) costs 16 gating shards
plus 7 × 8 sensitivity shards, 72 shards or 3.68 h.

## A0. Run the frozen rows as frozen (rejected; this was option A of IR-20261011-031924-C205C8)

- Cost: 88 gating shards plus 9 × 32 sensitivity shards, about 20 h. The intervention estimated 20–41 h.
- The r = 0 SE rule is still not met at four points.
- The V6-sens condition fails by construction: 32 points have zero model dose, and the core-dominated
  r = 1–2 cm points probably fail too (`core_sens_estimate.py`).
- So V6 is report-only whatever the data show. The operator directive excludes this option.

## A1. Chosen design (revision 2): depth cut z ≤ R0 + σ_sem, frozen 1e-4 floor, aura share ≤ 5 %, per-point screen with measured denominator

| Item | Value |
|---|---|
| r = 0 points | 29, rule 8 shards (12 with safety) |
| Halo points | 104, rule 6 shards (9 with safety) |
| N_F (rule 12) | 16 |
| Sign-test depths before screening | 16 |
| k (rule 1) | 8 |
| Expected H_g | ≈ 58 points, 16 depths |
| Planning CPU time | 3.68 h |

- Every selection quantity is frozen (p1, p2, the 1e-4 floor), a data coordinate, the measured dose, or
  the MD-fit aura term (p23–p25, a declared role extension).
- None is model output, except the screen's 2σ′ precision term, whose residual dependence is stated in
  design D3.
- The exclusions match physical scope boundaries: aura, the beyond-range tail, and the beam-sensitive
  core.

## A2. Depth cut at 21.8 cm instead of 21.6 cm (rejected on budget; reported as contrary reading (b))

| Item | Value |
|---|---|
| r = 0 points | 30, rule 18 shards (27 with safety) |
| Halo points (with D7) | 106, rule 9 with safety |
| k (rule) | 1 |
| Cost | 27 + 7 × 8 = 83 shards, 4.24 h CPU, above the 4 h target |

- It would add one r = 0 point (21.8 cm, 13.5 % of the maximum). There 2 % is 0.036 mm, against a ~1 mm
  source positioning uncertainty.
- Discriminating reason: +0.56 h for one point that tests positioning rather than transport.
- This is a weak discriminator, disclosed in design D1. The (a)5 display cannot be excluded as an
  influence.

## A3. Model-count floor instead of a measured floor (rejected)

- Rule: gate the points where the rehearsal has all 20 batches non-empty, or rehearsal-estimated N at most
  the budget.
- It is cheap: the same N_F, by construction.
- Discriminating reason: selection bias. The relative SE of a point scales roughly as (model dose)^−1/2,
  so a point where the model under-predicts the halo is more likely to be excluded. The test would then
  remove its own failures.
- Using the rehearsal for precision sizing is fine. Using it to decide which points are compared is not.

## A4. Cumulative or integral halo observables (rejected as the gate; recommended as report-only)

Examples:

- the radial dose integral ∫ D 2πr dr over r = 1–6 cm per depth, which is close to Gottschalk's S_mixed;
- the halo-to-core ratio per depth;
- depth-integrated dose per radius.

Pros:

- Variance is far lower, because averaging over 145 points gives about 12× smaller SE.
- It is less sensitive to the ~1 mm positioning error.
- The integral is less sensitive to σ_x, because it conserves energy across radius.

Discriminating reasons against gating on it:

1. The data are sparse in r: 10 radii with missing cells, 165 of 460 cells missing, and r = 7 cm present
   only at some depths. An integral needs an interpolation rule over measured points, which is a new
   construction choice with no frozen tolerance. A tolerance would have to be invented after the
   rehearsal.
2. The frozen 15 % per-point tolerance has no derived counterpart for an integral. Deriving one (for
   example from the source's ±3 % normalisation and the MD-fit residuals) would be a new criterion chosen
   after seeing model precision.
3. The integral hides radial-shape errors, which are what the r-trend sign test is for.

Recommendation: report the per-depth trapezoid integral over the measured radii, model against data, as
an annotation only.

## A5. Fixed a-priori r ≥ 3 cm halo set instead of the per-point screen (rejected; the fallback if the screen is judged too flexible)

- With D7 this is 45 points (r = 3 / 4 / 5 cm with 26 / 16 / 3), only just above the coverage minimum.
- The core estimate gives large core changes at r = 3 cm near the peak: about +177 % for σ_θ × 1.2 at
  z = 21 cm. The core fraction there is small, but a screen would still be needed.
- It would change the frozen sign-test groups (1–3 vs ≥ 4 cm).
- The screen keeps the frozen groups and lets the transport decide sensitivity.

## A6. Two-stage sequential design (rejected)

- Rule: stage 1 with 8 shards; if the realized SE misses the targets, or the verdict is within the noise
  of a tolerance, add 8 more.
- This saves about 0.4 h in the good case.
- Discriminating reason: a data-dependent extension is optional stopping. The plan forbids decreasing
  counts and forbids repeating runs to obtain a pass (Amendment 14 (a)3–4).
- A group-sequential design with an α-spending rule could be valid, but the frozen criteria are point
  criteria with no α. Inventing one is a larger criterion change than A1.
- A fixed N_F with "report, never act on, the realized SE" is simpler and has no forking paths.

## A7. Share the nominal between the gating sample and V6-sens-F (rejected)

- It saves 8 shards (0.41 h).
- It breaks Amendment 17 (c)4, which says the gating sample is not reused for the sensitivity status.
- The noise of the gating verdict and of the screen would then be correlated: a point whose nominal
  fluctuates high could both pass the screen and fail the residual test, or the other way round.
- The saving is not needed within the 4 h budget.

## A8. Run V6-sens-F on CUDA f32 and spend the CPU on a larger N_F (rejected)

- It depends on a CUDA rate that has not been measured.
- After revision 2 the sensitivity sizing rule gives k = 1, so moving it gains nothing.
- The executing task's throughput step measures both rates. No reallocation is pre-declared.

## A9. Defer V6 entirely (option C of the intervention; rejected)

- Cost: zero host time.
- It leaves the operator's "develop … a computationally feasible follow-up" unanswered, and leaves
  decision 0041's halo claim open with no plan.
- A1 costs 3.7 h and gives the r = 1–5 cm charged-halo test that motivated V6.

## A10. Aura handling: no exclusion (option (ii)), or a 10 % aura-share cut (not chosen)

- **Option (ii)**: keep all 145 points and pre-declare the bias.
  - The gate would then test "charged halo with the aura omitted".
  - A correct model is expected to read 3–10 % low at r = 4–6 cm, which uses up to two thirds of the 15 %
    per-point tolerance at r = 6 cm.
  - It pushes the sign test toward failure.
  - Kept as contrary reading (a).
- **A 10 % cut**: 133 points and 23 depths.
  - It still lets the omitted physics consume up to 10 of the 15 percentage points.
- **The 5 % cut** (one third of the tolerance) gives 104 points and 16 depths.
- Decision (orchestrator, B1): option (i) at 5 %.

## A11. Screen with the frozen model-normalised metric (revision 1; not chosen)

- The denominator is the model's own dose.
- A point where the model halo is low has a larger model core fraction, so it is preferentially screened
  out (review B2).
- Replaced by the measured denominator, which has only a second-order residual dependence through σ′.
- Kept as contrary reading (d).


---

# Independent review record (reviewer agent, read-only; three rounds)

# Independent review of draft Amendment 18 (V6-F / V6-sens-F)

Reviewer: independent Claude agent, read-only except for this file. Date 2026-10-11.
Governing texts read: EXPERIMENT.md; v3-005-acceptance.md Amendment 14 (a), (g), (j), Amendment 15 (a)5,
Amendment 17 (c), Amendment 14a (a)-(c); operator response IR-20261011-031924-C205C8 (main checkout,
`response.json`); `gottschalk_v6.py` (via se_table.py); halo.tex (extracted from the pinned tarball,
sha256 verified by se_table.py), `tbl:parms` and the MD-fit equations (l.760-800).

Reproduction: `python3 -I se_table.py <worktree>` from a scratch copy gives an se_table.json that is
byte-identical to the one in this directory. No transport was run.

**Overall verdict: approve with listed changes. Not committable as drafted: three blocking findings
(B1-B3) must be resolved in the text first. None needs new data or model output.**

---

## 1. Checklist answers

| # | Verdict | Justification |
|---|---|---|
| 1 | pass | Design 1.1 keeps "failed — not sized", `V6_SHARDS = None` and `v6-unsized` unchanged. Nothing in the draft edits Amendment 14 (g), 17 (c) or 14a (a). |
| 2 | pass (with note) | 1.2 and 9 keep V6-F out of (a)6 and require the V6 failure line next to every V6-F verdict. Note: the draft never says what the (a)6 gate does with the failed V6 for V3-005C. See I1. |
| 3 | pass | 1.2 explicitly rejects the P6-D replacement pattern. That fits the operator wording "preserving the original failure". |
| 4 | pass | 1.1 reserves lv5c r 12/20 and hr5c r 3. The V6-F seeds do not touch them (checked arithmetically). |
| 5 | fail | Point set, observable, E₀, tolerances, screen metric, coverage, N_F and k are fixed. **The seeds and base of the gating execution (later task) are not declared.** Item 7 declares only V3-005C seeds (B3). The aura question (B1) also leaves the scope incompletely specified. |
| 6 | pass | "Not established (beam-sensitive)" is defined and stated never to be a pass (rows table, item 6). |
| 7 | pass | 10 (a)-(c) are fixed and declared verdict-neutral. Two more are recommended: B1, B2. |
| 8 | pass | Rows table, V6-sens-F: a zero nominal or variant mean makes the point beam-sensitive by definition. |
| 9 | pass | 5, last bullet: realized precision is reported only, with no top-up and no dropping. A6 rejects optional stopping. |
| 10 | unclear | Cannot be verified, because the (a)5 display content is not available to me. Every rule can be written in frozen quantities. 21.6 vs 21.8 is discriminated by budget (A2), not physics. Retaining 21.2-21.6 (gradients 12-38 %/mm) is the opposite of a favourable cut. Contrary reading (b) reports 21.8. See I3. |
| 11 | pass (note) | Outputs contain only rel SE, non-empty counts and paired σ_Δ. `g.ratio_paired` does compute the ratio in memory; it is deleted, not printed or written. Rel SE indirectly encodes the model dose magnitude (rse² ∝ 1/D). That leak is weak and is not a residual. |
| 12 | pass (script) / unclear (process) | The script never outputs or prints the ratio. I cannot verify that no other session displayed it. |
| 13 | pass | `core_sens_estimate.py` uses only beam parameters (p5-p7), R0, E₀ and analytic Fermi-Eyges MCS. It reads no ionmc output and no measured dose. Physics accuracy is limited (m3), but adequate for planning. |
| 14 | **fail** | Yes, it changes only the consequence, and the metric and 0.15 are unchanged. **But it can select by agreement.** \|R_tot − 1\| ≈ f_core,model·\|R_core − 1\|, so a point where the model under-predicts the halo has a larger model core fraction. It is therefore more likely to be screened out. This is the A3 bias in weaker form (B2). |
| 15 | pass | Every retained gating point keeps 2 % (r = 0) or 15 % / 90 % / 10 % and p > 0.01. |
| 16 | pass | Reproduced: 88 at (22.0, 3) with 28.1 %; r = 0 116 / 651 / 8910 / ∞; V6-sens k = 32 on the 153-point subset. |
| 17 | pass | Reproduced: 29 pts, 8 (12) shards; 145 pts 30/29/28/20/21/17, 10 (15); 24 sign depths; σ_Δ,max = 0.0988 → k = 7; k = 8 → 2σ_Δ = 0.0698. |
| 18 | pass | 72 × 184 s = 13 248 s = 3.68 h. Each step is 184 s, below 1200 s. CUDA 294 / 2944 s, stated as unmeasured. The figure is moot if V3-005C does not execute (B3). |
| 19 | pass (note) | The worst H_F point (1.533, 3) has 20/20 batches non-empty. N = 16 against a rule value of 15. The expected worst realized SE is 2.27 % (halo) and 0.34 % (r = 0), so there is margin even without the safety factor (rule value 10 without it). 1.25 ≈ 1.5 σ of a Gaussian SE-of-SE; that is thin for skewed batches. Precision only flags and never gates, so the risk is reduced power, not a false verdict. |
| 20 | pass (note) | An under-estimated σ_Δ only shrinks H_g, which is conservative, and coverage then decides. With safety 1.25, only 1 of 145 points would exceed 2σ_Δ = 0.075 at k = 8. The proxy concern is real for σ_x variants at shallow r = 1 cm, but those points are screened out by the core sensitivity anyway. |
| 21 | unclear in draft → reachable | The draft does not answer this. My planning estimate uses core-only Fermi-Eyges with a 0.4 cm sphere average, the core scaled to the measured r = 0 dose, measured-dose denominator and rehearsal σ_Δ/√8. It gives H_g ≈ 95-97 points at all 6 radii, with ≈ 20 sign depths (r = 1: 17 points at z ≥ 14 cm; r = 2: only 2; r = 3: 19-21). Coverage looks reachable, but the 1-3 cm group is mostly r = 1 (deep) and r = 3. Add this to the design. |
| 22 | pass | lv5c r21 = 20502004-019, r22 = 20503004-011, hr5c r13 = 20504004-019. All are above max consumed 20482036. They are disjoint from 20493004+, 20494004+ and 20501004+ (r20 k ≤ 999 stays below 20502004). They are below every rehearsal seed (rehearsal seeds are base 20505000 + 1000 r + k ≥ 20506000, i.e. **outside** the "20505000-20509999" base range; still disjoint). They are below 20511004+. |
| 23 | pass | lv5c 23 (= hr5c 13) and hr5c 11/12 (= lv5c 21/22) are declared unused. lv5c 24 is forbidden. |
| 24 | unclear | Specified in the text (fail closed; CPU r21 vs CUDA hr13 disjoint, so the paired z is valid), but not implemented yet. |
| 25 | **fail** | The transport statement is correct: neutrons and gammas go to `nuclear_escaped_neutron/gamma` and nothing neutral is deposited. The scope consequence is wrong: the source's own MD-fit aura term puts a non-negligible aura share **inside** H_F at r = 4-6 cm (B1). |
| 26 | fail (minor) | The risk (0.05-0.17 mm per 2 % at 21.2-21.6; I verified gradients of 11.7/23.7/38.3 %/mm) is in the checklist, not in design.md. D1 says only that these points gate. Add it to D1 and 12. |
| 27 | pass | 12 "cannot establish" lists aura, beyond 21.6, r ≥ 7, beam-fit independence, ±3 %, other energies and the frozen V6 claim. Wording must change after B1 (r = 5-6 cm halo is partly aura). |
| 28 | pass | 12 and the rows table disclose shared beam lineage and E₀ from R0 = d80. |
| 29 | **fail** | The argument is given (9.1-9.4), but the draft plans a V3-005C report-only execution plus a later gating one. That contradicts the intended disposition (not executed in V3-005C). It also creates a two-look history that 9.4 can only assert is "not a retry" (B3). |
| 30 | pass (moot) | R1 is disclosed. The V7-R fixture was just regenerated at 1a86a08c (37e87ec). If V6-F code does not land in V3-005C (B3), R1 disappears and should be deleted. |
| 31 | unclear | Future process. The text states the order, and `shard_table_ok` is required at the evidence head. |

---

## 2. Findings

### Blocking

**B1 — The aura is not confined to the excluded points. The r-trend test is then biased by a declared
out-of-scope effect.**

- The MD fit of the source (`tbl:parms` p23-p25; D_au = p23·E_2D(r_n, r_au), r_au = p24(1 + p25 r_n),
  r_n = r/(R0/3)) gives D_au = 1.9 / 1.65 / 1.42 / 1.21 / 1.03 / 0.86 × 10⁻⁴ MeV/g/p at r = 1-6 cm,
  independent of z.
- Unit check: at r = 6, 8, 10 cm beyond the range, D_au/D_meas = 0.94-1.12, consistent with an
  aura-dominated tail, so D_au is in MeV/g/p.
- Over H_F, D_au/D_meas has median / max:

  | r | median | max |
  |---|---|---|
  | 1 cm | 0.0 % | 0.1 % |
  | 2 cm | 0.3 % | 1.5 % |
  | 3 cm | 0.9 % | 7.4 % |
  | 4 cm | 3.2 % | 16 % |
  | 5 cm | 6.0 % | 15 % |
  | 6 cm | 10.1 % | 13.5 % |

  9 of the 17 r = 6 cm points exceed 10 %.
- The model transports no neutrals, so a correct charged-halo model is expected to read 3-10 % low in
  the ≥ 4 cm group and ≈ 0-1 % low in the 1-3 cm group.
- Consequences:
  - The sign test (the sign of median(≥ 4) − median(1-3) per depth; it fails at ≥ 19/24 same sign) is
    pushed toward failure by physics declared out of scope.
  - At r = 6 cm, 8-14 % of the 15 % tolerance is consumed.
- The statement in item 2 that the excluded points are where the aura lies, and the claim in 12 that the
  test is of "the charged-secondary halo", are therefore inaccurate.

Required, choose one and declare it now:

- (i) An aura-share exclusion computed from p23-p25 and the measured table only. At 5 % this leaves 104
  points (r = 4: 16, r = 5: 3, r = 6: 0) and 16 sign depths before the screen. At 10 % it leaves 133
  points and 23 depths. This is a role extension of `tbl:parms` (frozen as "construction for the beam
  only"). It must be listed as difference D7 with its reason and recorded in the registry role text.
  The unexcluded verdict becomes a contrary reading.
- (ii) Keep H_F and state that the halo component tests "charged halo with the aura omitted". Report
  the per-point MD aura share, and pre-declare the expected bias direction.

Option (i) is preferred, because it keeps the gate interpretable. In either case the per-point
D_au/D_meas becomes a reported, verdict-neutral quantity.

**B2 — The per-point screen can select points by model–data agreement.**

- The screen ratio's denominator is the model's own total dose. Where the model's halo is low, the core
  fraction rises, |D_var/D_nom − 1| grows, and the point is more likely to be removed.
- Example: f_core 0.30 → 0.35 when the model halo is 20 % low. A core change of 0.40 then moves the
  point from 0.12 to 0.14, the edge of 0.15 − 2σ_Δ.
- This contradicts the D2/D3 rationale and the reason A3 was rejected.

Required, choose one and declare it now:

- (a) Gate the screen on |D_var − D_nom|/D_meas + 2σ_Δ′ < 0.15, which is model-vs-model in the
  numerator only, with a measured denominator and σ_Δ′ the matching paired SE. Report the frozen-ratio
  screen as a contrary reading. This must be listed as a criteria difference with its reason.
- (b) Keep the frozen metric. Pre-declare the bias, and report the halo verdict under the
  measured-denominator screen and the residual distribution of H_F \ H_g as contrary readings.

(a) is preferred. Either way, the claim "cannot preferentially drop points where the model is low" must
be removed or qualified for D3.

**B3 — The draft does not match the intended disposition, and the gating execution is under-declared.**

- Items 7, 8, 9 and 11 plan a V3-005C execution at bases 20481004/20491004 with new hashed-source code
  (risk R1), followed by a later gating run. The orchestrator's disposition is: report-only, **not
  executed** in V3-005C, and executed and gating later on a fresh base.
- A report-only first look followed by a gating second look is a two-look history. "Not a retry" in
  9.4 cannot be enforced by text once the first result is known.
- Required:
  - Remove the V3-005C execution: the steps table, the 20481004/20491004 seeds, the V3-005C budget lines
    and R1.
  - State that no V6-F code enters the V3-005C hashed set.
  - Declare the gating execution's seed scheme now: base = the later task's declared fresh base,
    r_index mapping as in item 7, and a disjointness check against consumed sets and all 2050xxxx
    rehearsal seeds actually used (20506000-20525xxx, not just 20505000-20509999).
  - Declare that the base value itself is fixed by that task's amendment before observation and is
    never re-used.
  - Declare that V6-F is executed exactly once as gating.
  - Note that the 10000-spaced fresh-base sequence collides across bases for r_index ≥ 10 (e.g.
    20511004 + 21000 = 20531004 + 1000). The later amendment must check this explicitly.
- The budget (3.68 h CPU) can stay as a planning figure for the later task, re-measured by its own
  throughput step.

### Important

- **I1 — Merge-gate consequence of the failed V6 is unstated.** Amendment 14 (a)6 lists "V6 with the
  status fixed by V6-sens". 17 (c)7 covers only "report-only", not "failed — not sized". The operator
  wrote "keeping the applicable merge gates enforced". Amendment 18 must state that it does not alter
  how (a)6 treats the failed V6 for V3-005C. If that treatment needs a ruling, it belongs in a separate,
  reviewed amendment, not in V6-F.
- **I2 — Sign-test power at the coverage floor.** At 10 depths only 10/10 same sign fails (9/10 gives
  p = 0.021). At 24 depths, ≥ 19/24 fails. The minimum is legitimate and pre-declared, but it is the
  barely non-vacuous point. Changes:
  - count coverage depths **after** tie exclusion;
  - state the power, e.g. the probability of failing for a given systematic r-trend at the realized SE;
  - consider ≥ 12 depths.
  - Note that H_g is expected to contain almost no r = 2 cm points, so the 1-3 cm group median is mostly
    r = 1 (deep) or r = 3.
- **I3 — The 21.6 cm cut.**
  - The arithmetic is right: σ_sem = 0.2925 cm and R0 + σ_sem = 21.614 cm.
  - "+1 σ" is a convention. The gradient argument applies equally to 21.4-21.6, and the aura argument
    applies only from about 22.0-22.2 cm. The discriminator is budget (A2), and the text should say so
    plainly in D1.
  - Contamination by the (a)5 display cannot be ruled out by inspection. Three things weigh against a
    favourable choice:
    - the retained 21.2-21.6 cm points are high-risk: a few tenths of a mm of CSDA-vs-projected-range
      offset uses the 2 % tolerance;
    - 21.8 cm is reported (contrary reading b, sized at about 18/19 shards, so nearly adequate at
      N = 16);
    - every alternative cut is cheaper or equally priced in the opposite direction.
  - Acceptable with the disclosure.
- **I4 — Proxy and estimate disclosure in design.** Add the H_g coverage expectation (checklist 21) and
  the r = 0 distal-edge risk (checklist 26) to design.md.

### Minor

- m1. `se_table.py` keeps `SIGMA_DELTA_F = 0.025` (2σ_Δ = 0.05) and the columns `sens_k_for_sigma_0p025*`,
  while the design uses σ_Δ ≤ 0.0375. `k_F` is computed inline with 0.0375, and `sens_worst` is keyed to
  the 0.025 column. Align the constant, or rename it and the columns, so the archived table is not
  misread.
- m2. `ratio_paired` computes the sensitivity ratio in memory before deletion. Prefer a paired-SE-only
  helper so that the blindness claim holds by construction.
- m3. `core_sens_estimate.py` uses E_s = 15 MeV with X0 = 36.08 cm. The Gottschalk differential form with
  E_s = 15 MeV uses X_S = 46.88 cm (and f_dM). This overestimates MCS by about 30 % in T. My H_g estimate
  changes from 95 to 97 points with X_S, so the conclusions are unaffected. The script also has no
  sphere averaging and uses `np.trapz`, which is deprecated in NumPy 2. Its stated "upper bound" holds
  only if the halo itself is beam-insensitive. Say so.
- m4. The r = 0 component is itself strongly beam-sensitive. The core estimate gives +59 % on axis at
  z = 1.5 cm for σ_x × 0.8, so 2 % at r = 0 is about 1 % in σ_total. This is inherited from the frozen
  V6 and disclosed as shared lineage. 12 "can establish" should qualify the r = 0 statement as "given
  the fitted beam".
- m5. "Rehearsal family 20505000-20509999" names bases. The actual rehearsal seeds (e.g. 20517000,
  20521000+k, 20525000) lie outside that range. Phrase the disjointness check against the seeds
  actually used.
- m6. The "≥ 10× CUDA" planning rate and the 4.5 h total are moot for V3-005C after B3. Keep them as
  planning figures only for the later task.

---

## 3. Specific questions of the brief

1. **Pre-registration completeness.** Every quantity is fixed from frozen inputs (p1, p2, the 1e-4
   floor), the measured table, and precision-only rehearsal quantities. No model mean enters selection
   directly. Two gaps remain: the screen depends on the model dose through f_core (B2), and the
   gating-execution seeds are undeclared (B3). The 21.6 cut equals R0 + σ_sem by correct arithmetic, is
   defensible as a declared convention with budget as the discriminator, and is mitigated against (a)5
   influence by contrary reading (b) (I3).
2. **Scope and criteria changes.** D1-D6 are listed with reasons, and no frozen tolerance is relaxed.
   Missing items:
   - the aura-scope handling (B1, needs a D7);
   - the screen-denominator choice (B2);
   - the coverage-minimum numbers are new criteria; they are declared in item 6 but not in the
     differences table, so add them as D8.

   The per-point screen (D3) is a real change of consequence, not a hidden relaxation. It cannot turn
   a failing point into a pass; it can only remove points. That is exactly why B2 matters.
3. **Feasibility.**
   - N: ⌈(1.25·0.01350/0.005)²⌉ = 12 (r = 0) and ⌈(1.25·0.09093/0.03)²⌉ = 15 (halo), so 15, declared 16.
   - k: ⌈(0.09876/0.0375)²⌉ = 7, so k = 8 gives 2σ_Δ = 0.0698.
   - CPU: 72 × 184 = 13 248 s = 3.68 h, and every step takes 184 s ≤ 1200 s.
   - All verified.
4. **Seeds.** Verified disjoint, see 22. The later-task seeds must be declared (B3).
5. **Screen and coverage.** The rule is legitimate and pre-declared. The coverage minimum is a fair
   guard, but weak at its floor (I2). The analytic core estimate is sound as an order-of-magnitude upper
   bound on core sensitivity (m3). My estimate says coverage is reachable (≈ 95 points, 6 radii,
   ≈ 20 depths).
6. **Preservation of the failure.** The failure is preserved verbatim, and V6-F cannot replace V6 in the
   ledger or (a)6 under the text. The (a)6 consequence of the failed V6 itself needs an explicit
   statement (I1).
7. **Statistics.** The sign test on 24 depths is sound and fails at ≥ 19/24. It has low power near
   10 depths and is biased by the aura (B1). Using a single rehearsed pair as the σ_Δ proxy is
   acceptable, because errors only shrink H_g (conservative). The 1.25 factor is thin for skewed
   batches, but precision never gates.

**No blocking issue remains once B1, B2 and B3 are addressed in the text. All three can be resolved
before any V6-F evidence, without new model output.**

---

## Revision 2 re-check (2026-10-11)

Re-read: design.md (rev 2), alternatives.md (diff against v1/), se_table.py (diff against v1/),
core_sens_estimate.py.

Re-run: `python3 -I se_table.py <worktree>` in a fresh scratch directory reproduces se_table.json
byte-identically. `core_sens_estimate.py` part 2 reproduces "H_g ≈ 58 points, r = 1/2/3/4/5 cm with
17/0/22/16/3, 4 radii, 16 sign-test depths". I did not combine the revision-1 and revision-2 tables.

**Verdict: no blocking issue remains. Approve with listed changes.** The changes are R2-I1 and R2-I2
(important); they change no criterion. Once they are done, the V6-F items can be committed as drafted.
"Approve as drafted" is not given yet, because of R2-I1 and because Amendment 18 is incomplete without
the (a)6 item of design 1.3.

### Status of the revision-1 findings

| Finding | Status | Evidence |
|---|---|---|
| B1 aura | **resolved** | D7 at 5 % aura share. The aura term is p23–p25 from `parse_parms`, which also asserts the frozen p1, p2, p5–p7. H_F is 145 → 104 points (r = 1–5 cm: 30/29/26/16/3), 41 excluded (r = 3/4/5/6 cm: 2/4/18/17), 16 depths, all reproduced. The role extension of `tbl:parms` is disclosed and must be recorded before execution. The per-point aura share is reported, and the 145-point verdict with the bias declared beforehand is contrary reading (a). The can/cannot claims are updated to r = 1–5 cm with aura share ≤ 5 %. |
| B2 screen bias | **resolved** | Screen \|D_var − D_nom\|/D_meas + 2σ′/D_meas < 0.15. The revision-1 claim is withdrawn. The residual dependence through σ′ (a low model gives a smaller σ′, so it stays in H_g) is stated correctly; it is second order and runs against a favourable selection. The frozen metric is contrary reading (d), and H_F \ H_g is contrary reading (e). Sizing reproduced: worst σ′/D_meas = 0.0369 per 1e6 at (20.4, 5); rule k = 1; k = 8 gives 0.026. |
| B3 disposition and seeds | **resolved** | No V6-F step, seed, code or host time in V3-005C, and the V7-R fixture is unaffected. One gating execution in a later task. Seeds B + 7000/8000/9000 + k, with r_index ≤ 9 and k ≤ 999, so every seed stays below B + 10000 and the cross-base collision is avoided. Disjointness is required against the seeds actually used, including all C8 rehearsal seeds up to 20525999. The draft notes correctly that B = 20511004 would need a seed-by-seed proof. Original V6 seeds are not reused. |
| I1 merge gate | deferred, acceptable | Design 1.3 leaves it to a separate orchestrator item of Amendment 18. **Amendment 18 must not be committed without that item.** |
| I2 sign-test power | **resolved** | ≥ 12 depths after ties. The exact boundaries are verified: n = 12 fails at ≥ 11 (p = 0.0063), 10/12 gives p = 0.039; n = 16 fails at ≥ 14 (p = 0.0042), 13/16 gives p = 0.021. Power at q = 0.8, n = 12 is 0.275 (verified). |
| I3 21.6 cut | **resolved** | Budget is disclosed as a weak discriminator (A2: 83 shards, 4.24 h). The (a)5 influence is disclosed as not excludable. |
| I4 | **resolved** | Coverage expectation and distal-edge risk are in design.md. |
| m1–m6 | **resolved** | Constants and columns are aligned; `ratio_paired` is no longer called; core estimate uses X_S = 46.88 cm with f_dM, sphere averaging and `np.trapezoid`; the r = 0 claim is qualified "given the fitted beam"; rehearsal seeds are named as used; CUDA figures are dropped. |

### Contrary-evidence note on the two normalisations

**What was seen.** The design author saw the revision-1 model-normalised σ_Δ (0.0988 at (21.6 cm, 3 cm))
and the revision-2 measured-normalised σ′ (≤ 0.037 at every H_F point). Both come from the **same**
st-hi/nominal batch pair, so their quotient at one point is ≈ D_nom/D_meas, up to the delta-method
linearisation. That is a near-residual: the model is well below the measured dose at a gating point. The
note records this correctly.

**Does any rule depend on it? No, as far as can be checked.**

- (21.6, 3) has aura share 0.0192 (reproduced), so D7 keeps it in H_F.
- The measured-denominator screen keeps points where the model is low, so it stays a likely failing
  gating point. The observation can only make a pass harder.
- D7 (5 %), the measured-denominator metric and the ≥ 12 depths are the options this review recommended
  in revision 1 (B1 (i), B2 (a), I2), and that review was written before any revision-2 number existed.
- The timing claim ("fixed before the revision-2 script was run") cannot be verified from the files.
  The file times (review.md 20:32, se_table.py 20:36, design.md 20:39) are consistent with it.

**Is the note adequate?** Adequate in substance, with two changes (R2-I1, R2-m1).

### Remaining findings (revision 2)

**R2-I1 (important) — the full residual map can be regenerated.**

- `v1/se_table.py` is kept and still writes the model-normalised σ_Δ column for every point, from the same
  st-hi/nominal pair as revision 2's σ′/D_meas.
- Running it and dividing the two columns gives D_nom/D_meas of the 1e6 rehearsal at **all** points, not
  just one. Copies of the revision-1 table also exist outside this directory: the reviewer scratch
  directory `/home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/v6rev/se_table.json`
  (revision-1 table, created by this review and never combined with revision 2).
- Required before commit:
  - Strip the σ_Δ column from `v1/se_table.py`, or do not commit/archive it next to revision 2. State
    in the note that the revision-1 table must not be regenerated.
  - Delete or quarantine the revision-1 table copies, including the reviewer copy above; the
    orchestrator decides this, since I am read-only here.
  - Never place revision-1 and revision-2 per-point precision columns in the same committed artefact.

**R2-I2 (important) — Amendment 18 is incomplete without the (a)6 item of design 1.3.** The V6-F items
are sound, but the amendment must not be committed or merged before the orchestrator's item on how
(a)6 treats "failed — not sized". That item must be reviewed, because it is the one place where the
merge gate could be weakened.

**R2-m1 (minor) — record the implied number.** The note says "well below". For a contrary-evidence
record it is better to state the implied quotient at (21.6, 3), or a bound such as "≤ 0.37 from
0.037/0.0988", so later readers can judge its weight. Record it at that one point only.

**R2-m2 (minor) — the zero-mean rule is no longer needed.** With a measured denominator, a zero nominal
mean no longer makes the metric undefined. "Zero nominal or variant mean ⇒ beam-sensitive" now removes
exactly the points where the model is lowest, which is the A3 bias in miniature. In practice it is moot:
every H_F point had 20/20 non-empty batches at 1e6, and the screen uses 8e6. Either:

- drop the rule: |0 − 0|/D_meas = 0, so the point gates and fails on its residual; or
- keep it, report every point it removes, and state the bias direction.

**R2-m3 (minor) — choice of B.** State that B is the next unused base that satisfies checks (a)–(e), so
a base that fails (b) is skipped (declared, not reused). Also state that skipping a base is recorded.

**R2-m4 (minor) — coverage margin.**

- The coverage expectation (58 vs 40 points, 16 vs 12 depths) assumes a beam-insensitive halo and uses
  the st-hi σ′ for all variants. "Not established (beam-sensitive)" is therefore a realistic outcome,
  not a remote one. Say so in item 6.
- The ≥ 4 cm group will be mostly r = 4 cm. With D7, the sign test compares r = 4 with deep r = 1 and
  r = 3.

**R2-m5 (minor) — the aura term is crude.** D_au is a z-independent empirical term of a fit with 15 % global
rms. Classification near the 5 % threshold carries that uncertainty. It is disclosed implicitly as shared
lineage; one sentence in D7 would make it explicit.

No blocking issue remains.

---

## Amendment 18 (c) check (2026-10-11)

Read: `../amendment-18-draft.md` sections (a), (b) and (c), against design.md revision 2, Amendment 14
(a)6, Amendment 17 (c)7 and the operator response. I also checked that R2-I1 is done:

- `v1/se_table.py` now writes σ_Δ = ∞ (quarantined, `ratio_paired` not called);
- the reviewer's revision-1 table copy in `.../tmp/claude-1000/v6rev/` is gone.

**Verdict: (c) does not weaken merge gate (a)6. Approve (c) with two wording changes. (b) has no number
that contradicts design.md; five minor wording fixes are listed.**

### (c) merge gate

**Consistent with the governing texts, and no weakening.**

- "Failed — not sized" is treated as a failed scientific gate: not a pass, not report-only.
- 17 (c)7 is correctly limited to the report-only status.
- V6-F is explicitly not a stand-in for V6.
- The consequence is stated plainly: V3-005C cannot merge on V6 under (a)6.
- No passing record is created.
- This matches the operator's "keeping the applicable merge gates enforced" and the binding C590CF
  condition "Merge only after the applicable scientific and implementation gates pass … Do not create
  an implementation-only passing record".

**(c)-1 (important): the planned intervention is framed as an open question.**

- The text pre-announces an "ambiguous-boundary intervention" on "integration with recorded scientific
  failures". EXPERIMENT.md category 5 is for materially incompatible interpretations of an explicit
  constraint. The operator has already stated the merge condition explicitly (C590CF (a)).
- A request phrased as an open integration question invites a gate relaxation by the operator.
- To keep (c) from becoming a route around (a)6, add the following:
  - The default without an operator decision is **not merged**; PR #82 stays unmerged.
  - The request proposes no relaxation, re-labelling or replacement of any gate. It asks only how the
    operator wants recorded failures of a task carried. The precedent is V3-005B interpretation B: keep
    unmerged and continue on the task branch.
  - Any operator answer is recorded verbatim, and it changes no row verdict.

**(c)-2 (important, provenance): the V10-A failure has no traceable record.**

- (c) states a V10-A "deterministic failure at the evidence head". That fact is in none of the texts I
  reviewed, and no record of it exists in the v3-005b worktree. Its latest commit is C9a, df59655, so no
  V10-A evidence commit has landed.
- If it is a real observation, (c) must cite its record: run id / step partial and evidence SHA.
- (c) must also state whether the (e)9 fallback (ii′) trigger was evaluated, since a V10-A (a)/(b) failure
  for p+O with (c) passing and P7 passing fires it.
- If it is not yet recorded, remove it from (c) and restrict (c) to V6.

**(c)-3 (minor): quote the gate list exactly.**

- The list should read "V7-R including its calibration step" (14 (a)6 verbatim).
- It should also name the implementation gates of (a)6, with P6 → P6-D per Amendment 15 (a)5. Otherwise
  the paraphrase reads shorter than the frozen gate.

### (b) summary versus design.md revision 2

**Checked and correct:**

| Item | Values |
|---|---|
| z_F | 21.6 / 21.614 cm |
| r = 0 | 29 gating; 5 report-only |
| Halo | 145 → 41 removed → 104 at r = 1–5 cm |
| Screen and coverage | screen formula; 40 / 3 / 12 after ties |
| Expected H_g | ≈ 58 points, 16 depths |
| Sizing | N_F = 16 (rule 12, × 1.25); 1.7 % / 0.34 %; k = 8 × 7 runs; R0 variants not run; 72 shards ≈ 3.7 h |
| Seeds | B + 7000 / 8000 / 9000; r_index ≤ 9; 20525999; original V6 seeds not reused |
| Disposition | not executed in V3-005C; once, gating, later |
| Contrary-evidence note | matches |
| B1–B3 | ↔ items 2, 3, 5 |

**Minor wording fixes:**

- (b)-1. "Point set fixed from frozen inputs and the measured table only": D7 also uses p23–p25, which are
  not among the frozen inputs (the frozen ones are p1, p2, p5–p7). Add "and the MD-fit aura term
  p23–p25 (a declared role extension of `tbl:parms`, recorded before execution)".
- (b)-2. "Replaces the all-or-nothing V6-sens condition": it could be misread as replacing row V6-sens.
  Write "V6-F's own per-point screen (instead of an all-or-nothing condition); V6-sens unchanged and
  failed".
- (b)-3. "Every step ≤ 1500 s": the binding limit is the 1200 s effective limit of 14a (c)5, and a step
  takes about 184 s. Also mention the backend-agreement rule (17 (c)8), which is part of the pass
  condition (iii) and is missing from item 1.
- (b)-4. "The file governs where this summary is shorter" should read "the file governs in any difference".
- (b)-5. Item 7 should state the re-check verdict exactly as recorded. It was "no blocking issue; approve
  with listed changes", with R2-I1 now done. R2-I2 was this (c) item. R2-m1 to m5 were minor and are
  open unless the design file shows them applied. Item 6 could carry the implied bound at
  (21.6 cm, 3 cm) (R2-m1).

No blocking issue in Amendment 18 (a)–(c) once (c)-1 and (c)-2 are addressed.
