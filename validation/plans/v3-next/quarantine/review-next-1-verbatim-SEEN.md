# QUARANTINED verbatim copy of review round NEXT-1 (2026-11-01) — contains an outcome-informed statement (section 4); SEEN; never cited for model selection

The active record `validation/plans/v3-next/selection-2026-11-01.md` carries this review with that one statement redacted (Codex REVIEW-cddf6ada: the exact statement stays solely in hash-bound quarantine). Content below is verbatim.

---

## Round NEXT-1 review (reviewer, verbatim)

**Verdict: REVISE.** The selection is careful and mostly honest, but four things need to change before it is recorded or executed:
- A′'s status and pre-registration need tightening.
- The criterion-2 rule is underdetermined and leans toward an outcome already known.
- The ranking should put V3-008 first.
- T1's step order has the outcome computed before the rule is fixed.

No binding-rule violation as drafted, if the conditions below are kept. I modified no file.

**1. A′: scientific justification and status (question (a)).**
- **(a) The physics basis is real and TOPAS-independent.** The LA150 H-1 MF1 text says the R-matrix fit covers 0–150 MeV, so 150 MeV is the fit edge. D-E9 B0 places the 130–150 MeV step inside the evaluation. Free NN cross sections have no structure below the pion threshold. SAID is the standard NN partial-wave analysis.
- **The decision to look was outcome-prompted.** The investigation followed the gain residual through the D-E2 → D-E4 → D-E6 → D-E8 chain. F3 shows the planner knows that "corrected p-p + a better p+A shape" is the closing pattern in the seen records. So A′ is one half of a seen two-component pattern, and that must be disclosed as such.
- **(b) "Physics correction, no screen, pin not moved" is an honest status only if the later use is constrained now.** Otherwise a later V5/V11 run with A′ becomes post-hoc selection. The pre-registration must state:
  1. **Adoption is irreversible.** A′ is adopted or withdrawn once, on criterion 2 and its construction checks alone, before any transport with an A′ table. If adopted, it stays the p-p component of record whatever later transport, screen or V5/V11 results show. No reversion to LA150 carry-up and no alternative p-p source may be chosen after any outcome.
  2. **No combination with seen shapes.** A′ is never combined with M2, M2′, R-B, CHIPS or any shape whose screen or swap arm has been seen. Any later p+A component is chosen only by a pre-registered rule that references no V5, V11, D-E or Geant4 outcome (for example a held-out measured-data gate, D2), and that rule is registered before any transport with the combination.
  3. **Any later qualification:** fresh bases, unchanged tolerances, the go/no-go discipline of Amendment 20 (c) carried, and an explicit statement that the two-component closing pattern was seen (F3). That qualification is therefore not blind to the pattern, only to its realisations.
  4. **Inputs.** The A′ table is built from acquired SAID data with sha256 binding, and default `la150-carryup` stays until decision 0043. Recording it as "the p-p component of record" happens only through decision 0043, not by a note.
- **(c) Disclose the (b′)4 history explicitly.** Full-domain replacement avoids the continuity check that blocked (b′)4. That is legitimate only because B0 later showed the 125–150 MeV dip is in the evaluation, so the (b′)4 check compared SAID against an evaluation feature. Say exactly this; do not let full replacement read as a way around a failed check.

**2. The criterion-2 rule is not sound as sketched (question (b)). Important.**
- The discrepancies are already seen: table/SAID like-for-like ≈ 0.92 at 125 MeV and 1.12–1.14 at 150–250 MeV. Choosing the uncertainty components and k now therefore largely determines the outcome. The rule must leave no free magnitudes.
- **Concrete problems:**
  - **(i) "LA150 χ²-scaled" has no base quantity.** The file has no MT2 covariance, and T1 step 1 is meant to establish that. A scale factor needs something to scale.
  - **(ii) SM16 − SP07 is a model-to-model difference, not an uncertainty.** At most it is a lower-bound proxy.
  - **(iii) The method residual is a computational discrepancy.** Once resolved it should be near zero, and it is not a physics uncertainty.
  - **(iv) SAID outputs gave Err = 0 in all scans so far.** A SAID uncertainty is likely not available.
- **Required rule text:**
  - Each component must come from an acquired source with stated uncertainties. Acceptable sources: documented SAID error estimates; a spread among independent partial-wave analyses (SAID SM16, SP07 and a Nijmegen-type PWA, if acquirable); or **measured pp elastic data with their stated errors** at the comparison energies.
  - The measured-data route is the most defensible. Note that the PDG points at 152–250 MeV are the V10 report-only set and have been seen, so disclose that, or use held-out EXFOR pp elastic data.
  - Fix per-energy evaluation: which energies, discrepancy defined per energy, and the logical combination (for example "at every scan energy ≥ 150 MeV").
  - Fix how the region below 150 MeV is treated. Full replacement there is a single-source consistency choice, not something criterion 2 establishes; say so.
  - **Fail closed:** if any required component cannot be formed from an acquired source, criterion 2 is "not formable", A′ is withdrawn, and this is recorded. No planner-chosen magnitude may substitute.
  - k = 3 is acceptable once fixed. State whether "discrepancy" is |table/SAID − 1| or a level difference, and at which energies.

**3. Ranking (question (c)).**
- **Reorder: rank 1 = V3-008 (He/C/O EM transport), as the main line of work.** Reasons:
  - The operator directs work "toward the original scientific and software requirements" and says not to merely repeat the named sequence. Ion transport is the largest unmet MUST block: V2-ION, "release cannot be proton-only", and the precondition for V3-009.
  - It is independent of PR #82.
  - The elastic chain has already consumed about 3 weeks.
- **A′/T1 is better placed as a bounded parallel side task (light work) with a hard stop.** Its most probable honest outcome is "criterion 2 not formable", because the uncertainty inputs are likely missing. That closes H3 cleanly at about 1 day of cost.
- Rank 3 B-wiring and rank 4 D2 are correctly placed, and D2 is correctly read-only.
- The rejected list is correct, including "A as literally posed", default-off elastic (operator-reserved), and any host qualification now.

**4. Binding rules (question (d)). No violation as drafted, with these conditions:**
- **No large data in Git:** SAID scans at 12.53–250 MeV go to the reference-data cache through acquisition records. `ppsrc.json` holds dimensionless ratios only, which matches the elastic projection rule.
- **Quarantine:** V10 (PDG) and V10-A stay report-only. If PDG pp points are used in criterion 2, declare them as seen.
- **Planner arithmetic:** the F2/F3 estimates, including "M2 + SAID p-p predicted to pass", must never be committed or cited in a record. That sentence is a tuning roadmap. Keep it out of the dated note and the amendment, and keep only the qualitative F3 hazard statement.
- **Seeds:** T1 draws none. V3-008 bases from 20561004 are clear of V6-F (maximum 20560019) and of the reserved hr r_index blocks. Register them in the seed-family test or the V3-008 guard before any run, and skip 2098xxxx and 2099xxxx.
- **History isolation:** the V3-008 worktree branches from v3/develop, and nothing is read from the v1/v2 trees.
- **Single process and one heavy job** are kept. The A16 CI fix stays separate.

**5. T1 scope and acceptance tests (question (e)). Well scoped; fix the order and add three tests.**
- **Order:** step 2 (resolve the 1.2–2.2 % method difference) produces the like-for-like ratios that criterion 2 tests. Fix the rule (step 3), including the method definition (the P6-tested `coul`, the η/k and angle-interpolation conventions), **before** the final computation. Steps 2 and 3 can be developed on fixtures. The record lists every number seen before the rule was frozen.
- **Tests to add:**
  - **(i) Freeze:** the criterion-2 rule, its components, sources and k are committed and reviewed before the computing step runs. The computing step refuses unless that commit is an ancestor of the current HEAD, the tree is clean, and the code hashes match the reviewed commit (the D-E8/D-E9 pattern), with a once-only attempt log.
  - **(ii) Provenance:** every uncertainty input is read from an acquired, sha-bound source. A test asserts that no uncertainty magnitude is a literal in the code.
  - **(iii) Fail-closed:** the "not formable" path, where a missing component leads to withdrawal, is unit-tested.
- The stop condition is correct: if criterion 2 fails or is not formable, no builder work follows.

**Fixes before recording or executing:**
1. Reorder the ranking: V3-008 first, with A′/T1 as a bounded parallel task.
2. Pre-registration clauses for A′, items 1(b)1–4 and 1(c).
3. Criterion 2 with sourced, fail-closed components, per-energy definition and treatment below 150 MeV; measured pp data with stated errors preferred, and any seen V10 points declared.
4. T1 step order and the three added tests.
5. Keep the F2/F3 arithmetic, especially "M2 + SAID predicted pass", out of every committed text.

