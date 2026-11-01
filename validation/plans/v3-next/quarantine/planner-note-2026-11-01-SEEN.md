# QUARANTINED planner note (2026-11-01) — SEEN, exploratory/tuning arithmetic, INELIGIBLE as qualification evidence

Status: this is the scientific planner's working note that produced the next-work selection of 2026-11-01. It contains outcome-informed arithmetic on seen records (sections 0 F2–F3: residual estimates and a two-component prediction). Codex REVIEW-07694336 required it to be preserved, hash-bound and quarantined so that its influence on model selection can be audited, instead of being redacted. Rules: it is never cited as evidence, never used to choose, tune or gate a model, and never read as a prediction to be confirmed; any later candidate that matches its two-component pattern is "selection by seen outcome" (reviewer round NEXT-1, finding F3) and is excluded unless chosen by a pre-registered rule referencing no seen outcome. Content below is verbatim.

---

# V3 next-task selection (planner, 2026-11-01) — for independent review

Basis: worktree v3-005b @067252e8 (PR #82 unmerged). Read: REQUIREMENTS.md, experiment/v3/AGENTS.md + REQUIREMENTS.md,
docs/development/program-plan.md, validation/requirement-ledger.md, v3-005-acceptance.md (criteria; Amendments 18-21 incl.
(b')4, (b''')4, (c), (e)-(h); dated notes of 2026-10-18 and 2026-10-23), decisions 0041 (D-1, rejected alternatives) and 0042
(§2, open findings, dated sections 10-18/10-19/10-23), D-E8 and D-E9 reports, D-E9 plan §0-1, V7-R2 plan header, steps_v6f.py
header, operator directive 2026-11-01 (+ memory of the 10-18 / 10-19 responses). Primary check done today: MF1/MT451 of the
cached ENDF/B-VIII.0 proton file p-001_H_001.endf (zip sha256 27bcafb8...).
Not read in full: integration-route.md, v6f-followup.md, the V7-R2 plan body, steps_v5d/run_suite wiring state, review.md.

## 0. Findings that drive the ranking (evidence, all from committed records)

F1. **No p+A shape alone can pass the V11 gain clause (ii) at 200 MeV with the current p-p channel.** D-E6 swap-pA with
    Geant4's own CHIPS shape leaves a residual of +0.00437 (swap-both +0.00413) against the allowed 0.00365 (Amendment 20 (c)5);
    D-E8 failed (ii) at 200 MeV for all tested arms (M2 0.00463, M2' 0.00462, R-B 0.00463). So another p+A-only candidate
    (TENDL, JENDL-5, a repaired M3, ...) is predicted STOP at the frozen screen. Repeating the candidate sequence with a new
    p+A shape would cost review cycles without a passable outcome.
F2. **A p-p change alone cannot pass either.** Under D-1, the base gain excess at 200 MeV is 0.00830 (archived; 0.0073 using
    the A1 mean). A SAID-level p-p channel (table/SAID like-for-like 1.12-1.14 at 150-250 MeV, i.e. s ~ 0.88) shifts the
    200 MeV gain by roughly -0.0025 (B3 slopes 0.016-0.021 per ln s; D-E6 sigma_pp arm -0.00287). That leaves about 0.005-0.006,
    which is above 0.00365. F closure would be about 0 (D-1 is the base), so (i)/(i-b) fail at 200 MeV (V5 peak/plateau -4.9 %
    against a band of 0.96 %). This is a planner estimate from seen records. It is not a rule and is not to be run as a test.
F3. **Only the two-component combination closes in the seen records.** The D-E6 "all" arm (Geant4 p+A shape + Geant4-level sigma_pp)
    leaves +0.00136. Another route: M2 (withdrawn at trigger 3, 0.5 % over) with a SAID p-p level is predicted to pass (0.00463
    - ~0.0025). Any candidate built now to that pattern is therefore **selection by seen outcome**. That post-hoc hazard
    cannot be removed by fresh seeds, because the screen and lv5d would only replicate a known prediction.
F4. **The p-p channel above the LA150 node has a scientific defect, independent of TOPAS.** Verified today, the MF1 of the
    LA150 H-1 file reads: "p-p elastic cross sections calculated for Ep between 0 and 150 MeV from an R-matrix analysis of
    p-p cross-section and polarization data in this energy range. The maximum nuclear partial wave allowed in the fit is l=6
    ... chi-squared per degree of freedom is 1.8" (EVAL-FEB98, Hale). So the 150 MeV node is the edge of the fit. B0 (D-E9)
    shows that the 130-150 MeV step (-18.3 % / -15.1 %) is in the evaluation, not in the builder. NN cross sections have no
    structure below the pion threshold (about 280 MeV). The carry-up S(E)/S(150) propagates the edge value to 250 MeV and is a
    project construction (decision 0041 D-1 4). SAID SM16 (world NN PWA) and BGG/PDG agree with each other at 150-250 MeV
    (table/SAID 1.12-1.13 vs table/BGG 1.101) and both lie about 10-13 % below the table. The table also lies 2-8 % *below*
    SAID at 100-125 MeV.
    Open (D-E9 criterion 2 "not established"): the SAID/LA150 uncertainties are not in hand, and the two sigma_anchor
    computations differ by 1.2-2.2 % with the cause unresolved.
F5. **Historic constraint on A.** The same idea was already pre-registered as (b')4: a level-only re-anchoring at E >= 150 MeV
    with a ramp from 125 MeV. Its pre-declared continuity check (+6.8 % > 3 % at 125 MeV) blocked it, and that outcome was seen.
    The check failed because of the dip that B0 later (10-19) attributed to the evaluation. A new p-p rule is legitimate only as
    a *different construction* motivated by new evidence (B0 + MF1). It must disclose (b')4 and must not re-use (b')4 with its
    check removed.
F6. **The requirement gap is large.** According to the ledger and the program plan, V3-006 (outputs/CLI/capability contract),
    V3-007 (beamlets and influence matrices: V1-MUST-009/010), V3-008 (He/C/O EM transport), V3-009 (ion nuclear and
    fragmentation; V2-ION "release cannot be proton-only"), V3-011, V3-012 (V2-PERF targets) and V3-013 (qualification and
    promotion) are all **not started**. Since 2026-10-10 the experiment has spent about 3 weeks on the elastic channel of one
    species. V3-008 depends only on V3-003, which is merged, so it is off PR #82's critical path.
F7. **Option B is not an integration route.** The elastic switch already exists (`PhysicsOptions.elastic`, default True,
    config.py l.110). "Default-off" is closure item (b''')4 2, which is reserved for the operator together with the deferred closure.
    Under Option II, V5 and V11 are binding. Carrying them as failures and merging *is* the deferred reduced-gate closure.
    REQUIREMENTS.md requires CI before integration and independent validation, not a particular gate set. The gate set is the
    operator's Option II selection, so B changes nothing without the operator.

## 1. Ranked options

### Rank 1 — A' (revised form of A): p-p channel source = SAID SM16, as a stand-alone physics correction (no V5/V11 claim)
- **Construction (proposal for Amendment 22 + decision 0043).** Replace the transported p-p density over the whole p-p
  domain (E >= E_min,pp = 12.53 MeV, which is unchanged) with SAID SM16: [DSG_SAID - sigma_Mott] over the same |mu_CM| <= 0.96
  half-sphere, using the P6-tested `coul`, the same angular cut and the same table schema. Level and shape are energy-dependent.
  There is no seam, no ramp, no threshold chosen from data and no carry-up. The fixed 150 MeV Hale shape above 150 MeV is
  dropped too. Fail-closed rules: NI <= 0 in any used interval, or E > 300 MeV (the SAID range in hand), raises an error.
  The LA150 reconstruction stays as a report-only cross-check (X-pp).
  - Why full replacement and not ">= 150 MeV": a threshold rule is what (b')4 tried, and any new threshold would be chosen
    knowing the dip. Single-source replacement has no free parameter.
- **Justification, without TOPAS:** F4 (the evaluation's own domain and fit edge, B0, the absence of NN structure, and SAID as
  the standard PWA).
  - Lineage, declared: SAID shares world-NN-data lineage with Hale (pre-1998 pp data) and with BGG/PDG. It shares none with
    Geant4 models. Geant4's p-p level is BGG, so the change moves the table toward TOPAS. That must be stated as a contrary
    reading, not hidden.
- **Evidence it can produce:**
  - deterministic reproduction of SAID (construction check);
  - like-for-like ratios to LA150 at <= 125 MeV (report-only);
  - V10 p-p level vs PDG (report-only, shared lineage);
  - a D-E4-type thin-target p-p rate vs Geant4 (report-only; it closes nothing, per (b)5(iv)).
  It produces **no V5/V11 evidence**. Per F2, the pre-registration must state that no screen and no lv5d V5/V11 run is
  justified by A' alone.
- **Pre-registration needs:**
  - (a) A criterion-2 rule, fixed before computing: the discrepancy must exceed k × the combined uncertainty. Components: the
    SAID PWA spread (SM16 vs SP07, or SAID error output if available), the method difference once resolved, and LA150
    chi2-scaled. Planner choice: k = 3. If the rule fails, A' is withdrawn and recorded.
  - (b) Construction checks: normalisation <= 1e-6; units/frame ratio in [0.80, 1.05] (already confirmed 0.86-0.93);
    NI > 0 over the used range; continuity of the SAID level across the scan grid.
  - (c) A full disclosure list of everything already seen: sigma_anchor values (A20-3), both table/anchor sets, PP-LEVEL
    ratios, the D-E4 rate ratios 1.165/1.159, the D-E6 sigma_pp arm, all D-E8 arms including M2, the B3 values and slopes,
    and the A1 / Lead A rows.
  - (d) An explicit ban on combining A' with M2, M2', R-B or any shape whose screen arm was seen (F3).
  - (e) No seeds. Any optional report-only harness run uses new 2099xxxx seeds asserted by the seed-family test.
  - (f) Whether the pin moves: **not now.** The table is built, checked and recorded as "the p-p component of record for any
    later elastic candidate". Moving the pin would invalidate the hashed set and the fixture without enabling any gate.
- **Cost:**
  - sandbox only: acquisitions (SAID DSG scans at the grid energies 12.53-250 MeV; SP07 tables), minutes;
  - builder extension + tests: about 1-2 days including Codex rounds;
  - host: 0 h.
- **Rule risks:**
  - It could look like tuning toward TOPAS. Mitigations: no free parameter, full-domain replacement, the MF1 and B0 basis,
    and the explicit no-V5/V11 statement.
  - It touches a STOPped chapter. Mitigation: it is a separately justified investigation (operator 10-19: a STOP "does not
    prohibit a separately justified investigation").
  - Independent review is required before commit.
- **Dependency:** none. It is the first task.

### Rank 2 — C1: start V3-008 (He / C / O electromagnetic transport) on a new worktree from v3/develop
- **Justification:** this is the largest unmet MUST block. It covers V1-MUST-005 / V2-ION and the "not proton-only" release
  condition, and it is a precondition of V3-009 (carbon fragmentation, a "central release requirement"). It is independent of
  PR #82 and of the elastic hazard. It has independent references:
  - NIST ASTAR for helium (the V3-002 pipeline has ASTAR);
  - ICRU 90 (`data/icru90.py` exists);
  - TOPAS carbon cases (working since V3-010: 290 MeV/u, about 38 ms/history);
  - measured ion Bragg-peak and range data (open item: which datasets are acquirable).
- **Evidence:**
  - ion-specific stopping and range vs tabulated data;
  - straggling and lateral spread vs theory and TOPAS;
  - backend parity across python, warp-cpu and CUDA.
- **Pre-registration:** decision 0043/0044 (effective charge, ion straggling, ion MCS), and a v3-008 acceptance plan frozen
  before any result. Also needed:
  - new seed bases from the fresh sequence (>= 20561004, skipping 2098/2099);
  - a lineage table (Geant4 ICRU73 data are shared with TOPAS and must be declared);
  - the oxygen supported domain declared.
- **Cost:** about 2-3 days plan/decision/review and 5-7 days implementation and qualification.
  - Host: reference runs one at a time. TOPAS carbon 1e5 histories takes about 1 h. A lv8 suite is estimated at a few hours
    under single-process.
- **Rule risks:** low.
  - Merge-conflict risk with PR #82 in species.py, kernels.py, funcs.py, reference.py and warp_driver.py. Mitigation: keep ion
    EM in new modules and rebase the later of the two branches.
  - The CI red on the v3-005b branch does not affect a develop-based task. The A16 digest fix is the orchestrator's, as noted.
- **Dependency:** none on #82. It can start on the same day as A' (planning is light work; there is one heavy job at a time).

### Rank 3 — B-wiring (no-regret engineering, timed after the model state settles)
- Scope:
  - the Amendment 20 (f)8 items: lv5d/hr5d constants, run_suite step names, rehearsal_preflight, the force list,
    V7R_REF_INDEX collision, the calibration timing rule and PARTIAL lines, V6-F code in the hashed set;
  - the seed-family test extension (f)7.
- No seed is drawn and nothing runs on the host.
- It is needed by *every* future route (closure or model change), but the fixture (about 0.7 h host) must follow the code
  freeze. So the wiring can be done any time, and the fixture only once the pin is final.
- Excluded: changing the elastic default (operator-reserved, F7).
- Cost: about 2-3 days with review. Host: 0 h until the freeze.
- Not a merge route (F7).

### Rank 4 — D2: research note on p+A elastic shape sources covering 200 MeV (no construction)
- Candidates:
  - TENDL-2023 p-O016 / p-C012: rejected in 0041 only because it was unreachable through the old proxy;
  - JENDL-5 proton sublibrary / JENDL/HE-2007 (open: coverage, method, and whether elastic MF6/MT2 is present to 200 MeV);
  - measured p+O / p+C elastic sets at 100-200 MeV in EXFOR, as held-out evaluation data (V10-A data are seen).
- Purpose: decide whether an independent-lineage two-component candidate (A' + new p+A shape) could ever be registered
  without the F3 hazard. It can only be justified if the shape is chosen and gated against measured data *before* any transport
  screen. Otherwise the elastic route ends at the operator's closure decision.
- Cost: hours (researcher + reviewer). Read-only.
- Risk: the F3 hazard. That is why this is a research note, not a candidate.

### Rejected (with reasons)
- **A as literally posed** ("revise p-p at E >= 150 MeV + fresh-seed screen"):
  - it is the (b')4 rule without its failed check (F5);
  - its screen outcome is predictable STOP under D-1 (F2), so a fresh-seed screen adds no information and spends a once-only
    D4/CRN design.
- **A' combined with M2, M2', R-B or CHIPS-like shapes:** selection by seen outcome (F3). M2 is withdrawn by a frozen trigger,
  M2' is STOPped, and CHIPS is TOPAS lineage (V2-EVID).
- **Changing the screen (Lead A: archived-baseline SE in the band):** allowed only prospectively for a future candidate. It
  must never re-evaluate D-E8 (M2' also failed (i) and (i-b), so it would not change that verdict anyway). Defer until a
  candidate exists.
- **New p+A-only candidate (M3 extension, TENDL, JENDL):** predicted STOP at 200 MeV (F1).
- **Default-off elastic now:** operator-reserved (F7).
- **V3-006 or V3-007 now:** both depend on V3-005. Basing them on develop collides with PR #82's config.py (+141 lines),
  simulation.py (+108) and cli.py (+119). Second in line after V3-008, or after #82's disposition.
- **Any host qualification (lv5d/hr5d, V6-F, V7-R2):** deferred by the operator, or cannot change a gate in this state.

## 2. Recommended sequence (~2 weeks; one heavy job at a time; planning and review run in parallel with light work)

| Days | Work | Host |
|---|---|---|
| 1-2 | T1 (below): p-p deterministic package + Amendment 22 / decision 0043 drafts; reviewer round; Codex. In parallel: physics-researcher (Opus) V3-008 research brief (stopping sources, effective charge, ion MCS/straggling, acquirable measured ion range data) | 0 |
| 2-3 | If criterion 2 holds: SAID scan acquisitions, builder `--pp-source said-sm16`, construction checks, tests, commit, Codex. If it fails: record the withdrawal (dated note) and close H3 as "not established" | 0 |
| 3-5 | V3-008 worktree from v3/develop; decision + acceptance plan frozen; independent review; Codex | 0 |
| 5-6 | D2 research note (read-only) on p+A sources to 200 MeV. Then decide: register an independent two-component candidate, or report to the operator that the elastic route needs the closure decision (a genuine reserved decision, an intervention only if no other route exists) | 0 |
| 6-13 | V3-008 implementation: reference backend, then Warp CPU/CUDA, CI tests, rehearsal, TOPAS He/C reference runs one at a time, lv8 qualification at a fresh base | a few hours, serial |
| 14 | Checkpoint (progress log, ledger); B-wiring scheduled when the elastic pin is final | 0 |

## 3. First concrete task (T1) — bounded; plan-amendment draft + deterministic sandbox checks; no transport

**Task: V3-005D-PP0, "p-p channel source: criterion-2 package and Amendment 22 draft"** (on the v3-005b branch; exploratory
class until adopted).

Interfaces (proposal):
- `validation/scripts/elastic_attribution/pp_source_check.py` (new; pure functions, deterministic). It outputs the
  `ppsrc.json` record: dimensionless ratios only (Git policy of the results README).
- CLI: `python -m validation.scripts.elastic_attribution.pp_source_check --said-ids <...> --la150 <cache sha>`
- Builder hook (later step, not T1): `elastic_build.py --pp-source {la150-carryup,said-sm16}`. The default stays
  la150-carryup until the decision.

Steps:
1. Record the LA150 H-1 MF1 statement (verbatim lines above) and whether the file carries MF33/MF34 for MT2. This is a
   deterministic read of the cached zip.
2. Resolve the 1.2-2.2 % method difference between the researcher's and the reviewer's sigma_anchor values with one code path:
   the P6-tested `coul`, a documented eta/k definition and the documented angle interpolation. Report both legacy sets beside it.
3. Fix the criterion-2 rule in the draft **before** computing. Uncertainty sources:
   - the SAID spread SM16 vs SP07 (acquire SP07 pp DSG/SGT at 100-250 MeV; open: availability through the CGI);
   - SAID error output, if the CGI offers one (open);
   - the method residual from step 2;
   - LA150 sqrt(chi2/dof) scaling.
   k = 3.
4. Compute criterion 2 once and record it, whichever way it falls.
5. Draft Amendment 22 + decision 0043: the construction of Rank 1, its checks and fail-closed rules, the disclosure list, the
   "no V5/V11 / no screen on A' alone" statement with the F2 reasoning, and the combination ban.
6. Independent review (reviewer agent), then Codex on the committed draft.

Acceptance tests for T1:
- the unit tests for `pp_source_check` reproduce `coul` values, half-sphere normalisation, and a units/frame ratio in
  [0.80, 1.05] on fixtures;
- the leak check (no absolute cross section in Git) passes;
- a test asserts no transport import and no seed draw;
- the record is bit-reproducible at a clean SHA;
- the amendment text passes review with all seen numbers disclosed.
- Stop condition: if criterion 2 is not met, the A' construction is withdrawn and no builder work follows.

Budget: about 1 day sandbox. Host 0 h. No heavy job.

## 4. Open items (not verified by the planner)
- Whether the SAID CGI provides uncertainties, and whether SP07 is acquirable (needed for the step-3 rule).
- Whether SAID DSG-minus-Mott is positive over |mu_CM| <= 0.96 at all p-p grid energies down to 12.53 MeV.
- SAID scans were acquired only at 100-300 MeV. Full replacement needs scans at lower grid energies and an energy-interpolation
  rule for the shape. To be defined in Amendment 22.
- TENDL-2023 / JENDL-5 proton coverage and content for C-12 / O-16 elastic to 200 MeV.
- Measured ion range datasets (C/He in water) that are open and acquirable. ICRU 73 is likely paywalled.
- The current lv5d/hr5d wiring state of run_suite / steps_v5d (not inspected).
- The F2 estimate is the planner's arithmetic from seen records. It is not a computed record and must not be committed as one.
- The orchestrator's A16 digest fix is outside this plan.
