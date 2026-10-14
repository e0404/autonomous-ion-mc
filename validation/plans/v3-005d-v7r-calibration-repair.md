# V3-005D — V7-R calibration repair: diagnosis, offline evidence, repaired gate design (planner output, NOT evidence)

**Revision 4 (2026-10-14, after Codex REVIEW-bb2a8151)** — changes_required, two findings:
- **(1, important):** the §3.8 size-selection rule no longer uses the independence product OC_paired × OC_single with a fitted box shrink. It is replaced by a pre-declared rule on the **directly measured complete two-gate joint pass rate at κ = 1** (repository `run_calibration`, both gates on the same replicates, α_tost = 0.005, 48 × 225 replicates): R_meas = 400 per surrogate, select 5400/0.005 iff ≥ 378/400 joint passes (CP95 lower bound p_low with P(Bin(3200, p_low) ≥ 2909) ≥ 0.99) for both surrogates, else 7200/0.005, else operator. Measured now (script `scripts/joint_power.py`, summary `scripts/joint_power_summary.py`, data `data/joint_power_5400.*`, seed roots 20989920–20989923, 2.15 h sandbox): archive 383/400 (CP95 lower 0.9369), rehearsal 384/400 (0.9399) → **5400/0.005 selected**; budgets unchanged (≈ 33.9 h; 7200 ≈ 44 h; (B) ≈ 16.5 h; (C) ≈ 12.5 h), all subject to the operator decision. The Revision 3 product-rule powers (0.954 / 0.962) are kept only as labelled model values.
- **(2, minor):** "24-shard V7-R2 layout" in §6.4 corrected to 48 × 2.25·10⁶; the file was grepped for every other 24-shard / 4.5·10⁶ mention (the remaining ones refer to the dd976c5a layout or to Revision 1 and say so).
- **Round-4 independent review (REVISE minor, no re-measurement):** the Codex finding is answered; three stale labels fixed afterwards — §3.1 product-power formula and the §3.2 table entry 0.954 / 0.962 labelled as model values superseded by §3.8, §3.5 (C) cites the measured rates, §6.6 risk entry rewritten (measured, residual risk on the rebuilt pool).

**Revision 3 (2026-10-14)** — after the round-2 review (REVISE, minor, no blocking):
- **R2-1:** the §3.8 window rule (which would have selected 7200/0.005 because the archive single coverage 0.68664 lies outside [0.680, 0.686]) is replaced by a computable rule: 5400/0.005 iff P(power clause) ≥ 0.99 from the exact binomial OC at the measured κ = 1 coverages of both surrogates, else 7200/0.005. The current coverages select **5400/0.005** (P = 1.0000 / 1.0000; `designs7.py`). Both budgets go to the operator: **5400: ≈ 33.9 h; 7200: ≈ 44 h; (B) ≈ 16.5 h; (C) ≈ 12.5 h.**
- **R2-2:** §5: lv 20511004 r_index **0–14** (not 10–14) lie inside the used rehearsal range.
- **Host execution of any V7-R2 step waits for the operator's decision on V7-R2's merge-gate role and on the budget** (§3.4, §4, §6.4).
- Seeds confirmed relative to the route's bases: lv **20531004** r_index 16 → shards 20547004–20547051, reference 20547052; **hr 20541004 r_index 6 reserved** (identical seeds); V6-F at 20551004 (its r_index 7–9 → 20558004+k … 20560004+k) is disjoint from 20547004–20547052, from 2098xxxx and from 2099xxxx.

**Revision 2 (2026-10-14)** — after the independent review `validation/plans/v3-005d-v7r-calibration-repair/review.md` (verdict REVISE, 1 blocking). Changes per finding:
- **F1 (blocking):** V7-R2 is now stated as a **new row** (V7-R's outcome rule says "without a further base"; the 14 (h)4 increase permission expired at the dd976c5a evidence run and is no longer cited). Its diagnosis record, design, pre-registration, implementation and review are autonomous; **its role relative to the failed V7-R in merge gate 14 (a)6 / the requirement ledger, and its host budget (≈ 33.9 h against the ≈ 9.5 h stated to the operator in IR-DEFAF4), are operator-proposal items** (§3.4, §4, §6.2). "No operator approval needed" is withdrawn.
- **F2:** transport re-sharded to **48 × 2.25·10⁶** (the route's layout; ≈ 631 s per shard, `shard_table_ok` margin verified against `steps_v5d.shard_plan_c`) (§3.4, §6.2).
- **F3:** new §3.8 — pre-freeze re-measurement rule on the regenerated rehearsal pool, with the power window [0.679, 0.687] for 5400/0.005 and the automatic switch to 7200/0.005 (window ≈ [0.676, 0.690]).
- **F4:** transport seeds moved to the route's bases lv 20531004 / hr 20541004 with the V6-F checks (a)–(e); collision arithmetic stated; hr r_index 6 at 20541004 reserved (§5).
- **F5:** the scripts behind `designs3–6.out` are now in `validation/plans/v3-005d-v7r-calibration-repair/scripts/` (`designs3.py` … `designs6.py`; re-runs reproduce every number, two labels corrected); the (B) numbers for 5400/0.005 recomputed (§3.5).
- **F6–F9:** retry-block rules (§5); single-gate bound wording "+ 0.01 box miss" (§1.2, §3.2); P(≤ 109 | 0.0374) = 0.172; parts 83/82; 27 % headroom; legacy κ mode for acceptance test 2; slow-host calibration rate to be measured; §1.4/§2.2 wording and the rehearsal paired_high 2.4σ residual; (C) keeps an end-to-end Monte Carlo check.

Status: planner document, 2026-10-14. Written in the sandbox against the V3-005B worktree
`/home/wahln/aiprojects/ion-mc-worktrees/v3-005b` (HEAD 749f768c; PR #82 unmerged). Nothing in it is
evidence, and nothing in it changes a recorded verdict. The V7-R failure at the evidence SHA dd976c5a
("uncertainty-coverage evidence not established", status partial `v7r-calibration-status.json`
content 8451f463…) **stands as recorded**. All numbers below are either read from the archived records or
computed in the sandbox with numpy (exploratory, seeds in the planner-diagnostic block 20989900–20989999,
declared in section 5).

Sources read: `validation/plans/v3-005-acceptance.md` (Amendments 12, 13, 14 (h), 17 (a), 14a (b), 18 (c));
`validation/scripts/transport/v7r.py`; `tests/ionmc/test_v7r_coverage.py` (fit, `margin_kappas`,
`run_calibration`, the clause table); `.ionmc-cache/transport-samples/dd976c5/78b/samples/v7r-calibration-status.json`
(per-surrogate counts and κ); `validation/results/transport/lv5c-dd976c5-single-process.json` (steps 30–79);
`validation/results/transport/README.md` (V3-005C table row V7-R); `decisions/0041-proton-nuclear-interactions.md`
(V7-R outcome); `tests/ionmc/fixtures/v7r/rehearsal/` (committed sidecars).

Scripts (exploratory): `validation/plans/v3-005d-elastic-attribution/scripts/binom_gate.py`, `implied_eps.py`, `candidates.py`, `candidates2.py`
(previous planner, re-run and checked), and `validation/plans/v3-005d-v7r-calibration-repair/scripts/kappa_direct.py`, `analyse_direct.py`, `fastbin.py`, `designs2.py`–`designs6.py` (this plan; every `validation/plans/v3-005d-v7r-calibration-repair/data/*.out` is regenerated by them).

---

## 1. Diagnosis

### 1.1 What the calibration measures

The V7-R calibration (Amendment 14 (h), 14a (b)3–5) runs, per surrogate (archive-fitted zero-inflated Gamma
k = 5.0703, π = 0.2167; pooled rehearsal block sums, 18 000 blocks), 3200 Monte Carlo repetitions of the
complete two-gate procedure on 7200 synthetic replicates. At the four non-nominal margins the bootstrap-t
interval is **deliberately miscalibrated** by a scale κ on its half-widths so that the true coverage equals
an edge of the region [0.6527, 0.7127]. The counts therefore measure **the gate's acceptance probability for a
hypothetical estimator whose true coverage sits at (or near) a region edge**. They are a property of the gate
procedure, not of the `escaped_neutral` estimator: the surrogates are synthetic and the κ ≠ 1 intervals are
wrong on purpose.

κ is solved once per surrogate (`margin_kappas`, `tests/ionmc/test_v7r_coverage.py`): bisection of the
**empirical** coverage of a single Monte Carlo sample — n = 20 000 pairs for the paired margins (seed
SEED_KAPPA + 1 = 20508001, one sample shared by paired_low and paired_high) and only **15 000** single
replicates for the single margins (20 000 drawn, but `single_ingredients` bootstraps only the 3:1 evaluation
split; seed 20508000, shared by single_low and single_high). The empirical coverage hits the target exactly,
but the true coverage at that κ differs from it by the sample error of the empirical coverage:
SE = √(c(1−c)/n) = 0.0034 (paired, n = 20 000) and 0.0039 / 0.0037 (single low/high, n = 15 000).
I reproduced all five recorded κ of `archive_fitted` **bit for bit** in the sandbox (worktree `.venv`,
numpy 2.5.3, the repository's own `margin_kappas`), so the κ of record are exactly this procedure.

### 1.2 The paired gate is an exact binomial test, so its operating characteristic (OC) is known exactly

The 3600 pairs are disjoint pairs of independent replicates, each pair with its own bootstrap stream; the hit
count is exactly Binomial(3600, c), where c is the true paired coverage of the bootstrap-t interval at the given
κ. The Clopper–Pearson TOST at α_tost = 0.04 per side accepts iff **2401 ≤ hits ≤ 2517**. Exact OC:

| n pairs | α_tost | acceptance set | P(accept) at 0.6527 | at 0.7127 | at 0.6827 (power) |
|---|---|---|---|---|---|
| 3600 | 0.04 (frozen) | 2401–2517 | 0.0374 | 0.0383 | 0.9636 |
| 3600 | 0.025 | 2407–2511 | 0.0231 | 0.0233 | 0.9396 |
| 3600 | 0.02 | 2409–2509 | 0.0195 | 0.0196 | 0.9292 |
| 7200 | 0.02 | 4783–5051 | 0.0197 | 0.0189 | 0.9993 |

(the single gate on 5400 evaluation replicates has the same binomial form on K_min/K_max; because the box
makes K_min ≤ K(δ*) ≤ K_max **whenever the box covers δ***, its false acceptance is at most the binomial OC of n = 5400
**plus P(box misses δ*) ≤ 0.01** (the 0.99 box level rests on the Amendment 13 Gaussian-marginal assumption, bootstrap-
checked): OC 0.0381/0.0382 at α = 0.04, 0.0242/0.0238 at 0.025, 0.0197/0.0192 at 0.02, each + 0.01 rigorously;
Amendment 13 measured 0.009/0.020 for Gaussian blocks.)

### 1.3 Why the frozen requirement cannot be met by the frozen design (a design defect)

Two independent defects, either one sufficient:

**(a) The gate-alone clause certifies a size-0.038 test against the limit 0.04.** The clause demands the
CP95 upper bound of the Monte Carlo rate ≤ 0.04, i.e. ≤ 109 of 3200 accepted (Amendment 14a (b)4). With a
perfect κ the paired gate's true rate is 0.0374 / 0.0383, i.e. 120 / 123 expected acceptances, so
P(≤ 109 | 0.0383) = **0.113** per margin (0.172 at 0.0374). The clause has to hold at both paired margins of
both surrogates: P(all four pass) ≈ 3.7·10⁻⁴ **even with an exact κ**. The limitation was recorded in advance
(Amendment 14a (b)4: "may be unprovable by Monte Carlo at any n") — it is now confirmed quantitatively. To
make it provable with α = 0.04 one would need ≳ 35 000 repetitions per margin (so that 0.0383 + 1.645·SE ≤ 0.04)
and a κ error ≲ 1·10⁻⁴ in coverage (n_κ ≳ 2·10⁷); at 8.1 s per repetition on the slower host that is
≈ 80 host-hours per surrogate, ≈ 160 h in all, and even then P(all clauses) is only 0.10 (design D0c, section 3).

**(b) The κ-solve error is of the same size as the whole margin of the test.** The OC is steep at the
edges: at n = 3600, α = 0.04, the false acceptance rises from 0.038 at the edge to ≈ 0.11 at 0.0046 inside it.
A κ solved at n = 20 000 moves the realised coverage by N(0, 0.0034), so the "false acceptance at the
margin" that the calibration records is a random quantity spread over ≈ [0.005, 0.15], depending on the
sign and size of one sampling error that is drawn once and frozen with the seed. The joint clause
(≤ 50/3200) and the single-gate clause fail for this reason (section 2).

### 1.4 What "margin acceptance 11.5 %" measures

369/3200 (rehearsal block sums, paired_low) is the acceptance probability of the paired gate for an estimator
whose true paired coverage is **≈ 0.6573**, not 0.6527: the κ of record (0.94939) lands 0.0046 inside the
region because the κ sample under-represented small studentised statistics. At that true coverage the exact
binomial OC predicts 0.115 (direct measurement of the coverage, §2.2: 0.6581 → 0.137 ± 0.02) — the observed
0.115 is consistent with what the gate should do there within ≈ 1σ. It is neither a
property of the `escaped_neutral` estimator nor a defect of the gate's finite-sample guarantee; it is the
consequence of evaluating the gate at a point that is not the margin. Likewise 370/3200 for the single gate at
single_high (section 2).

### 1.5 "Estimator correct" versus "calibration requirement satisfiable" — keep them apart

* **What the real data say (consistent with a correct uncertainty):** `escaped_neutral` paired coverage
  2467/3600 = 0.6853; the exact CP bounds at α = 0.04 are **[0.6715, 0.6988]**, ≥ 0.019 inside both edges
  (they stay inside at α = 0.025: [0.6698, 0.7004], and at 0.02: [0.6691, 0.7011]). `sec_p` 0.674 and
  `nuclear_local` 0.673 pass their Amendment 13 gates. For the paired gate the finite-sample guarantee is
  distribution-free and does not depend on the calibration Monte Carlo (an estimator with true paired coverage
  outside the region passes with probability ≤ 0.0383). The real-data passes are therefore consistent with a
  correctly estimated statistical uncertainty, and the lv5b 0.6586 paired failure was not reproduced.
* **What the pre-registered rule says:** V7-R passes only if the calibration step passes (Amendment 14 (h)4,
  14a (b)6). It did not, so **V7-R remains failed: "uncertainty-coverage evidence not established."** The
  observation above is *contrary evidence recorded next to the failure*, not a re-labelling. It cannot be
  converted into a pass by a post-hoc reading, by a re-run of the same calibration with more repetitions
  (14a (b)6: "no further count increase follows"), or by applying a new α to the already observed dd976c5 gate
  data (section 4).
* **The defect is in the calibration design**, i.e. in the combination "α_tost = 0.04 gate" × "CP95-upper
  ≤ 0.04 certification at 3200 repetitions" × "κ from one 20 000-sample". The estimator is not implicated by the
  calibration failure — and also not certified by it.

## 2. Offline evidence (sandbox, exploratory; reproducible from the committed sidecars and archives)

### 2.1 The κ of record are reproduced exactly

`validation/plans/v3-005d-v7r-calibration-repair/scripts/kappa_direct.py` re-ran the repository's `_margin_kappas` (worktree `.venv`, numpy 2.5.3 as on the
host) for both surrogates: all ten non-trivial κ equal the recorded ones **bit for bit** (archive: 0.938606/1.066375
paired, 0.929174/1.054671 single; rehearsal: 0.949392/1.073231 paired, 0.927500/1.055011 single). The committed
rehearsal sidecars at HEAD are byte-identical to those at dd976c5a (sha256 5d28a9b3… / d4f7d749…; the six
regenerations at base 20505000 produced identical block sums). Everything the calibration consumed is therefore
available offline.

### 2.2 Direct measurement of the true coverage at the recorded κ (the discriminating test)

Fresh independent samples (2·10⁵ pairs, 1.5·10⁵ single replicates per surrogate; seeds 20989901–03 / 20989911–13;
the repository's bootstrap-t functions unchanged; ≈ 5 min per surrogate). Coverage SE ≈ 0.0010 (pairs) and 0.0012
(single). "OC" is the exact binomial prediction of the gate's acceptance at that true coverage (paired: exact;
single: an upper bound, the box only lowers it).

| surrogate | margin | κ | true coverage | shift vs target | predicted acceptance | observed (3200) |
|---|---|---|---|---|---|---|
| archive | paired_low | 0.93861 | 0.65329 | +0.0006 (inside) | 0.0439 | 0.0469 (150) |
| archive | paired_high | 1.06637 | 0.71464 | +0.0019 (outside) | 0.0212 | 0.0209 (67) |
| archive | single_low | 0.92917 | 0.65047 | −0.0022 (outside) | ≤ 0.0172 | 0.0066 (21) |
| archive | single_high | 1.05467 | 0.71220 | −0.0005 (inside) | ≤ 0.0454 | 0.0241 (77) |
| rehearsal | paired_low | 0.94939 | 0.65814 | **+0.0054** (inside) | 0.137 | **0.115 (369)** |
| rehearsal | paired_high | 1.07323 | 0.71700 | +0.0043 (outside) | 0.0095 | 0.0216 (69) |
| rehearsal | single_low | 0.92750 | 0.64607 | −0.0066 (outside) | ≤ 0.0026 | 0.0013 (4) |
| rehearsal | single_high | 1.05501 | 0.70576 | **−0.0069** (inside) | ≤ 0.259 | **0.116 (370)** |
| archive / rehearsal | nominal (κ = 1) | 1 | paired 0.6840 / 0.6830; single 0.6866 / 0.6808 | | paired 0.963 / 0.964 | 0.968 / 0.960 |

Cross-checks at the *other* gate's κ also agree (paired acceptance at κ_single_high: predicted 0.092 / 0.107,
observed 0.097 / 0.135; at κ_single_low 0.0099 / 0.0070 vs 0.0109 / 0.0047). The previous planner's count-inverted
shifts (+0.0046, +0.0008 at paired_low; +0.0019/+0.0020 at paired_high) agree with the direct values within
1–2.4 σ. **Correction:** at paired_high a positive shift is *away* from the nominal (more conservative), not toward it.

### 2.3 Discrimination: κ-solve error, not surrogate misspecification

* **Signature.** Each κ sample is shared by the two margins of one kind, so a sampling error of the empirical
  distribution moves both realised coverages in the same direction (expected correlation 0.87): one margin lands
  inside the region (inflated acceptance), the other outside (deflated). Observed in all four (surrogate, kind)
  combinations: archive paired +0.0006/+0.0019, single −0.0022/−0.0005; rehearsal paired +0.0054/+0.0043, single
  −0.0066/−0.0069.
* **Size.** The shifts are 0.2–1.9 σ of the κ-solve SE (0.0034 paired, 0.0039/0.0037 single; the single solve uses
  only 15 000 of its 20 000 replicates because `single_ingredients` applies the 3:1 evaluation split — an
  implementation finding, harmless in kind but it enlarges the error).
* **Misspecification cannot act inside the calibration.** κ is solved and the Monte Carlo is run with the *same*
  sampler (truth = Gamma mean or pool mean), so a surrogate that misrepresents the real tally changes what is
  calibrated, never the agreement between target and realised coverage. The direct measurement closes the loop: the
  binomial OC at the measured coverage is consistent with the observed counts (reviewer: χ² = 14.5 on 10 rows,
  p ≈ 0.15). One row is in tension: rehearsal paired_high, predicted 0.0095 vs observed 0.0216 (z ≈ +3.1; the
  count-implied coverage 0.7146 differs from the direct 0.7170 by 2.4 SE of the direct sample). The rehearsal rows
  together are consistent with one ≈ 2σ common-mode error of the single direct sample; the residual is reported,
  not explained away.
* **Why archive 4.69 % but rehearsal 11.5 %.** Not the shape of the surrogate: the archive κ sample happened to be
  0.2 σ off at paired_low, the rehearsal sample 1.6 σ off (direct: +0.0006 vs +0.0054). At the OC slope near the edge
  (≈ 20 per unit coverage at n = 3600, α = 0.04) that is the whole difference. The archive paired_low still fails
  because of defect (a), section 1.3: even a perfect κ yields ≈ 0.038 and the clause needs ≤ 0.034 observed.
* **Representativeness (separate question, no verdict impact).** The real block sums (rehearsal pool, 18 000
  blocks of 500 histories) are unimodal, CV 0.18, range 515–2360 MeV, no zeros; the archive-fitted surrogate is a
  zero-inflated Gamma (21.7 % zero blocks) that reproduces the four studentised quantiles but not the block shape.
  Both surrogates give near-nominal bootstrap-t coverage at κ = 1 (paired 0.683–0.684, single 0.681–0.687), which is
  consistent with — not proof of — a correctly calibrated `escaped_neutral` uncertainty.

### 2.4 What an offline recalibration from the committed sidecars can and cannot establish

*Can:* reproduce the κ of record; measure the realised coverage at each κ; compute the exact paired-gate OC at any
coverage; design and size a new gate; show that the dd976c5a failure is explained by the design (sections 1–2).
*Cannot:* (1) change the V7-R verdict — the calibration of record is the C9 host step at the evidence SHA (14a (b)6);
(2) serve as calibration evidence for a new design — it would be sandbox-produced, unattested, and run after the
dd976c5a gate data were seen; (3) re-judge the observed dd976c5a gate data under a new α (post-hoc design choice);
(4) say anything about the real estimator beyond the surrogates' representativeness.

## 3. Repaired gate designs

### 3.1 Model and what "by construction" means here

Outputs of `validation/plans/v3-005d-v7r-calibration-repair/scripts/designs2.py` … `designs6.py` (helper `fastbin.py`) (logs `validation/plans/v3-005d-v7r-calibration-repair/data/designs*.out`). Paired gate: exact binomial
OC. Single gate: binomial OC on n_E with an effective box shrink of 0.0019·n_E counts per side, fitted to the ten
archived single-gate rates (shrink 0 = the rigorous upper bound on its false acceptance). κ error per kind
~ N(0, √(c(1−c)/n_κ)), same sign at both margins. Joint false acceptance = r × own-gate rate with r = 0.32 (largest
observed conditional ratio, 117/370) or **r = 1 (worst case: joint ≤ min of the two gates — the only assumption-free
bound)**. Power = OC_paired(c_p0) × OC_single(c_s0) at the measured κ = 1 coverages (**model value, independence
product; superseded for size selection by the measured joint pass rate of §3.8, Revision 4**). Clauses exactly as frozen,
at R repetitions: gate alone ≤ k_g, joint ≤ k_j, nominal ≥ k_p (R = 3200: 109 / 50 / 2909; R = 800: 22 / 9 / 735;
R = 400: 9 / 3 / 371). "By construction" = P(all clauses, both surrogates) ≈ 1 under r = 1 and the planned κ error.

### 3.2 Results

| design | pairs | α_tost | n_κ | R | FA paired / single (binomial OC given box coverage; + 0.01 rigorous) | joint power (archive / rehearsal) | P(all) r = 0.32 | P(all) r = 1 |
|---|---|---|---|---|---|---|---|---|
| D0 frozen | 3600 | 0.04 | 2·10⁴ / 1.5·10⁴ | 3200 | 0.0383 / 0.0382 | 0.945 / 0.953 | **0.000** | 0.000 |
| D0b frozen + exact-ish κ | 3600 | 0.04 | 10⁶ | 3200 | 0.0383 / 0.0382 | 0.945 / 0.953 | 0.000 | 0.000 |
| D0c frozen α, brute force | 3600 | 0.04 | 2·10⁷ | 40 000 | 0.0383 / 0.0382 | 0.945 / 0.953 | 0.10 | 0.000 |
| R1 | 3600 | 0.025 | 10⁶ | 3200 | 0.0233 / 0.0238 | 0.911 / 0.922 | 0.59 | 0.000 |
| R2 (previous planner's D2) | 3600 | 0.02 | 10⁶ | 3200 | 0.0196 / 0.0192 | **0.894** / 0.907 | 0.001 | 0.000 |
| R4 | 3600 | 0.03 | 10⁶ | 3200 | 0.0296 / 0.0295 | 0.926 / 0.937 | 0.48 | 0.000 |
| 5400 / 0.01 | 5400 | 0.01 | 10⁶ | 3200 | 0.0097 / 0.0099 | 0.975 / 0.980 | 1.000 | 0.75 |
| 5400 / 0.0075 | 5400 | 0.0075 | 10⁶ | 3200 | 0.0075 / 0.0072 | 0.966 / 0.973 | 1.000 | 0.988 |
| **5400 / 0.005 (recommended)** | 5400 | 0.005 | 10⁶ | 3200 | 0.0050 / 0.0050 | **0.954 / 0.962** (model; measured joint rate 0.9575 / 0.9600, §3.8) | 1.000 | **1.000** |
| 7200 / 0.02 | 7200 | 0.02 | 10⁶ | 3200 | 0.0197 / 0.0191 | 0.999 / 0.999 | 0.98 | 0.000 |
| 7200 / 0.01 | 7200 | 0.01 | 10⁶ | 3200 | 0.0098 / 0.0098 | 0.998 / 0.998 | 1.000 | 0.66 |
| 7200 / 0.005 (conservative variant) | 7200 | 0.005 | 10⁶ | 3200 | 0.0049 / 0.0050 | 0.995 / 0.996 | 1.000 | 1.000 |

κ accuracy is essential even for small α: 7200/0.01 drops from 0.66 to 0.14 (r = 1) at n_κ = 2·10⁵; 5400/0.005 stays
1.000 at n_κ = 10⁶ (0.77 at 2·10⁵ for 7200/0.005). **n_κ = 10⁶ (coverage SE 0.00048) is required.**

### 3.3 Which clause forces more transport (answer to "power or false acceptance?")

**Only the power clause.** The false acceptance at the margin is set by α (≈ α, almost independent of n); the
certification clauses need α ≲ 0.0075–0.01 to leave room for Monte Carlo noise and κ error. At n = 3600 pairs α ≤ 0.01
leaves a joint power of 0.71–0.81, far below the ≈ 0.915 true power that "CP95 lower ≥ 0.9 at R = 3200" needs; at
α = 0.025–0.03, where the 3600-pair power is still adequate (0.91–0.94), the certification clauses pass only with
P ≈ 0.5–0.6 (r = 0.32) and never under r = 1. **The feasible window at 3600 pairs is empty;** 5400 pairs
(10 800 replicates; 48 × 2.25·10⁶ shards) is the smallest of the evaluated sizes at which α = 0.005 keeps power ≥ 0.95. The real-data pass
probability of a correct estimator is unchanged: at true paired coverage 0.6827 / 0.6853, P(paired gate passes) =
0.969 / 0.959 for 5400/0.005, against 0.964 / 0.957 for the frozen 3600/0.04 (7200/0.005: 0.996 / 0.993). The smaller α
makes the gate *stricter*: the paired acceptance band shrinks to 3615–3762 of 5400 (0.6694–0.6967).

### 3.4 Recommendation — "V7-R2", a NEW row (design autonomous; merge-gate role and budget = operator proposal)

V7-R2 is a new gating row for the same estimator claim, not a continuation of V7-R: V7-R's outcome rule reads
"Fail: a recorded scientific failure without a further base" (14 (h)), the count-increase permission of 14 (h)4
applied only "in Amendment 14a before the evidence run, never after", and 14a (b)6 closes the calibration of record.
V7-R stays failed in every ledger. Designing, pre-registering, implementing, reviewing and (when scheduled) executing
V7-R2 is autonomous (the operator prioritised "V7-R calibration repair"). **Whether a V7-R2 pass stands in V7-R's
slot of merge gate 14 (a)6 and discharges V2-NUM / V1-MUST-030, and whether its ≈ 33.9 host-hours (against the
≈ 9.5 h stated in IR-DEFAF4) are spent, are operator-proposal items** — the operator ruled for the analogous V6-F that
its pass "does not authorize replacing V6 or merging". The proposal goes to the operator before any host step of V7-R2: **no V7-R2 host step (fixture regeneration
for V7-R2, κ step, calibration, transport) runs before the operator's decision on the role and on the budget
(5400/0.005 ≈ 33.9 h; 7200/0.005 ≈ 44 h; (B) ≈ 16.5 h; (C) ≈ 12.5 h)**; host execution is in addition gated by the
independently reviewed integration route.


* Transport: **48 shards × 2.25·10⁶** histories (225 replicates each, 10 800 replicates of 20 × 500-history blocks;
  pairs = replicate j of shard s with j of s + 24, s = 0..23 → 5400 pairs; H = shards 36–47 (2700), E = shards 0–35
  (8100)); reference 1·10⁶ (k = 48). Step time ≈ 631 s at the slow-host block rate (4.5·10⁶ / 1262 s = 3566 hist/s):
  `steps_v5d.shard_plan_c` needs ⌈1.08·10⁸ / 3566 / 1200⌉ = 26 shards, 48 declared → `shard_table_ok`, 47 % below the
  1200 s effective limit and 58 % below 1500 s (the 24 × 4.5·10⁶ of Revision 1 would have needed 26 > 24 and failed
  the check). Every other Amendment 13/14 (h) count unchanged (B = 1999, α_interval = 0.3173, 5000 box
  bootstraps, ≥ 300 pairs/replicates, the degenerate-miss rule, the bin mask rule).
* **α_tost = 0.005 per side for the two `escaped_neutral` bootstrap-t gates** (region [0.6527, 0.7127] unchanged).
  `sec_p` / `nuclear_local` keep the Amendment 13 rule (α_tost 0.04, CP / empirical Bernstein) unchanged — their
  calibration passed and their bounds are distribution-free; the Gaussian/A13 calibration shard is re-run at the new
  replicate count.
* κ: **n_κ = 10⁶ for both kinds, the full sample for the single kind** (no 3:1 split inside the κ solve), solved in
  a separate sealed `v7r-kappa` step (parts → per-row statistics with sha256 → combine bisection), consumed by every
  calibration part; the part re-derives nothing.
* Calibration Monte Carlo: **R = 3200** per margin per surrogate (unchanged, not decreased), both surrogates
  unchanged (fit seed 20505000 reproducing k = 5.0703, π = 0.2167; pool sha 5d28a9b3…/d4f7d749…), new seeds 2098xxxx.
* Expected OC: FA 0.0050 (paired, exact) / ≤ 0.0050 (single); joint ≤ 0.005 at every margin by construction;
  power 0.954 / 0.962 (modelled with the independence product; superseded by the measured joint rate of §3.8); P(all clauses) ≈ 1.000 even under r = 1; real-data pass probability of a correct estimator 0.97.
* Cost ≈ **33.9 host-hours** (section 6; the per-repetition calibration rate is extrapolated, see §6.1). Conservative
  variant 7200/0.005 (64 × 2.25·10⁶): ≈ 44 h, power 0.995. Whether this budget is spent is an operator item (above).

### 3.5 Cheaper alternatives (each is a design change that must be pre-registered; B and C also change how a
clause is evaluated and therefore need the operator)

* **(A) Same as V7-R2 at 7200 → no; at 5400 it is already the cheapest as-written design found.** Vectorising the
  per-replicate statistics (keeping the per-replicate generator draws bit-identical) could shorten the Monte Carlo;
  gain unmeasured, optional, reviewed code change.
* **(B) Exact paired certification + reduced Monte Carlo (operator proposal).** The paired gate's acceptance is
  Binomial(n, c) for *every* estimator with paired coverage c, so its false acceptance at the margins is an exact,
  κ-free, surrogate-free number (5400/0.005: 0.0050 / 0.0050). The Monte Carlo is kept only for the single gate, the
  joint and the power, at R = 800: P(all) 1.000 (r = 0.32) / 0.86 (r = 1); at R = 1600: 1.000 / 0.995 (5400/0.005,
  `designs6.py`; Revision 1 quoted the 7200/0.005 values 0.83 / 0.989). Cost ≈ 16.5 h
  (R = 800) or ≈ 22 h (R = 1600). Needs the operator because (i) "never decreased" (14 (h)4, 14a (b)3) forbids
  R < 3200 and (ii) evaluating the paired clause analytically is the "construction reading" that 14a (b)4 withdrew.
* **(C) Deterministic certification (operator proposal, cheapest acceptable).** Paired: exact OC at the margins
  (0.0050). Single: rigorous bound = binomial OC on n_E + 0.01 for box non-coverage (5400/0.005: 0.0150 ≤ 0.04);
  this bound is **assumption-bearing**: the 0.99 box level rests on the Amendment 13 Gaussian-marginal assumption,
  bootstrap-checked. (C) keeps an end-to-end implementation check: pinned acceptance-set tests plus a small Monte
  Carlo of the full two-gate procedure at a fixed known coverage (e.g. R = 100 at each paired margin). Joint ≤ min(gates)
  ≤ 0.0050 / 0.0150 ≤ 0.02 at every margin, without any dependence model. Power: the only Monte Carlo, at κ = 1,
  R = 400 per surrogate (P(CP95 lower ≥ 0.9) ≈ 0.993 per surrogate at the true 0.954 of the Revision 3 model; at the
  measured joint rates 0.9575 / 0.9600 of §3.8 it is higher). No κ step. Thresholds
  0.02 / 0.04 / 0.9 preserved; the clause *text* ("CP95 bound of the calibration Monte Carlo") is replaced by exact
  bounds for the false-acceptance clauses. Cost ≈ **12 h**.
* **(i) "Exact OC at the directly measured coverage at κ"** (coordinator's suggestion) is dominated by the κ-free
  exact OC of (B)/(C): measuring c(κ) adds uncertainty instead of removing it. With 2·10⁵ pairs the one-sided 95 %
  bound is OC(margin + 0.0018): 0.059 / 0.062 at 3600/0.04 (cannot certify ≤ 0.04), 0.038 / 0.040 at 3600/0.025,
  0.0105 / 0.0110 at 5400/0.005. It never helps the frozen α = 0.04.
* **(ii) "Fewer replicates with exact CP bounds"** is (B): workable only at small α, and it decreases R — operator.
* **Not available at any price:** re-using the dd976c5a transport data (their gates were observed under α = 0.04;
  a new α applied to them is post hoc) and repeating the dd976c5a calibration with more repetitions (14a (b)6).

### 3.6 Frozen tolerances

No frozen tolerance changes in V7-R2 or in (A)–(C): margins 0.6527 / 0.7127, joint ≤ 0.02, gate alone ≤ 0.04, power
≥ 0.9, B = 1999, R ≥ 3200 (V7-R2), surrogates, interval level. What changes in V7-R2: α_tost of the two
`escaped_neutral` gates (0.04 → 0.005, stricter), the transport layout (a new row's own size: 48 × 2.25·10⁶; 14 (h)4 is not invoked), n_κ
(2·10⁴ → 10⁶, full sample), the step layout and the seeds.

### 3.7 Required content of the new dated amendment (e.g. "Amendment 19", written before any new V7-R result)

1. Statement that the dd976c5a V7-R failure ("uncertainty-coverage evidence not established") and all its counts
   stand unchanged; the diagnosis of sections 1–2 is recorded as a finding with this document's numbers; the paired_high
   sign correction.
2. Scientific reason: the frozen α = 0.04 makes the gate-alone clause unprovable at R = 3200 (P ≈ 0.11 per margin
   even with an exact κ), and the 2·10⁴-sample κ error (SE 0.0034, measured shifts up to 0.0069) makes every margin
   clause a lottery; at 3600 pairs no α satisfies both certification and power.
3. The exact new parameters of 3.4 (48 × 2.25·10⁶ shards, pairing s ↔ s + 24, H/E split 36–47 / 0–35, α_tost = 0.005 for the `escaped_neutral` gates only,
   n_κ = 10⁶ full-sample, sealed κ step and its record schema, R = 3200, part sizes and shard count of section 6, the
   1500 s step limit and the 1100 s per-shard planning target, seeds of section 5, the fresh transport base).
4. The pre-computed OC table of 3.2 for the chosen design (exact paired OC; acceptance band 3615–3762 of 5400) and
   the expected P(all clauses), with the statement that the clauses are evaluated exactly as in 14a (b)3–4.
5. The pass rule unchanged in form: V7-R2 passes iff the calibration (κ step + all parts + combine) passes at the
   evidence SHA and all three estimators pass both gates; a miss is a recorded failure with no further count increase
   and no further base for V7-R2.
6. Contrary-evidence notes: the planner saw the dd976c5a real-data gate results (0.685 / 0.674 / 0.673) before
   choosing α; the choice is justified by the calibration operating characteristics alone, and the new transport base
   is unobserved. The real-data pass probability for a correct estimator under the new design (0.97) is stated.
7. The 2098xxxx seed family and the skip of 2098xxxx/2099xxxx in the open-ended fresh-base sequence.
8. If (B) or (C) is adopted instead, the operator's decision reference and the replaced clause text verbatim.
9. (Revision 2) V7-R2 is declared as a **new row**; V7-R stays failed; the amendment cites the operator's decision on
   V7-R2's role in merge gate 14 (a)6 / the ledger and on the budget (or states that the decision is pending and that no
   V7-R2 host step runs before it). The §3.8 pre-freeze measurements and the selected size are recorded; the transport
   layout 48 × 2.25·10⁶ with its `shard_table_ok` arithmetic; seeds lv 20531004 r_index 16 (20547004–20547052), checks
   (a)–(e) and the hr r_index 6 reservation; the generating scripts `validation/plans/v3-005d-v7r-calibration-repair/scripts/designs*.py`, `fastbin.py` are committed
   with the amendment and every cited table is regenerated from them.

### 3.8 Pre-freeze rule for the surrogates (Revision 2, F3)

The P(all) ≈ 1 of V7-R2 holds for the surrogates as measured here (κ = 1 paired coverage 0.6840 / 0.6830). The rehearsal
fixture must be at the evidence digest, and the elastic repair (same priority under disposition A) will change the
escaped-neutral block sums, so the rehearsal surrogate will be regenerated. The binding clause is power: at 5400/0.005
the joint power ≥ ≈ 0.92 that makes P(≥ 2909/3200) ≥ 0.99 requires the surrogate's κ = 1 paired (and single)
coverage within **≈ [0.679, 0.687]** (paired OC 0.911 at 0.678, 0.890 at 0.689; `designs6.py`); at 7200/0.005 the window
is **≈ [0.676, 0.690]**. Rule, fixed now and applied before the amendment's numbers are frozen and before any V7-R2
calibration or transport result:
1. Order: elastic decision → rehearsal fixture regenerated at the evidence digest → sandbox measurement on the
   regenerated pool and on the archive-fitted surrogate: κ = 1 paired and single coverage from 2·10⁵ pairs / 1.5·10⁵
   singles (≈ 5 min each), and the single-gate power at α = 0.005 by a sandbox Monte Carlo of the full gate at κ = 1
   → replaced in Revision 4 by the direct joint-power measurement of item 2 (R_meas = 400 per surrogate,
   ≈ 10 s per repetition in the sandbox: ≈ 2.2 h for 5400, ≈ 3 h for 7200) → amendment frozen.
2. **Selection rule (Revision 4, after Codex REVIEW-bb2a8151; replaces the Revision 3 product rule).** The Revision 3
   rule multiplied OC_paired × OC_single, which assumes independence of the two gates (they are evaluated on the same
   synthetic replicates and bootstrap seed base) and used a fitted box shrink for the single gate; it is withdrawn.
   The rule now uses the **directly measured complete two-gate joint pass rate at κ = 1**:
   - Measurement: `scripts/joint_power.py` — the repository's own `run_calibration` (`tests/ionmc/test_v7r_coverage.py`)
     with `v7r.pair_gate`/`v7r.single_gate` unchanged, both gates on the same replicates, the only change being
     `steps_v5b.V7_ALPHA_TOST = 0.005`; layout 48 × 225 replicates (5400 pairs; 64 × 225 for 7200); **R_meas = 400**
     repetitions per surrogate; seeds of record in the script (roots 20989920–20989923, +10 for 7200).
   - Bound: p_low = one-sided 95 % Clopper–Pearson lower bound of the joint passes. Requirement:
     P(Bin(3200, p_low) ≥ 2909) ≥ 0.99, i.e. **≥ 378 of 400 joint passes** (p_low(378/400) = 0.9224 → 0.9974;
     377 → < 0.99). P(select | true joint power 0.950 / 0.955 / 0.960) = 0.73 / 0.86 / 0.95.
   - **Freeze 5400/0.005 iff both surrogates reach ≥ 378/400;** otherwise run the same measurement for 7200/0.005 and
     freeze it iff both reach ≥ 378/400; otherwise no V7-R2 design is frozen and the finding goes to the operator.
   - The Fréchet bound p_paired + p_single − 1 from the same run is reported as a cross-check only.
   - This rule and R_meas were fixed before the measurement was started (Revision 4); it is re-applied to the rebuilt
     rehearsal pool after the elastic fix.
   *Result on the current pool* (`data/joint_power_5400.jsonl`, `.summary.txt`; 2026-10-14 10:01–12:10 UTC, one sandbox
   process, 9.7 s per repetition, 2.15 h): archive-fitted joint **383/400 = 0.9575** (paired 385, single 398), CP95 lower
   **0.9369**; rehearsal block sums joint **384/400 = 0.9600** (paired 385, single 399), CP95 lower **0.9399**; P(power
   clause | p_low) = 1.0000 for both; Fréchet cross-check (CP lower bounds) 0.9272 / 0.9310. **Both ≥ 378/400 → the rule
   selects 5400/0.005 (48 × 2.25·10⁶) today.** Disclosure: a 2-repetition timing trial (r = 0, both surrogates, both
   joint passes) ran before the rule was written; it was deleted and r = 0 was regenerated identically (deterministic
   seeds) inside the 400.
3. The measurements are recorded in the amendment as design inputs (exploratory, never evidence) with their seeds
   (planner block 20989900–20989999).

## 4. Autonomous versus operator proposal

Binding frame: recorded failures are preserved; no threshold is relaxed; no seed is repeated to obtain a pass;
a gate-design change is allowed only if pre-registered before new results, with the original failure preserved
and the change documented; scope/criteria changes go to the operator as proposals.

| Item | Class | Why |
|---|---|---|
| Record the diagnosis (section 1–2) as a dated finding next to the V7-R failure (results README, decision 0041 outcome) | **autonomous** | adds contrary evidence; changes no verdict |
| New amendment (section 3.7) defining the **new row** V7-R2: 48 × 2.25·10⁶ transport shards (5400 pairs), α_tost = 0.005 for the two `escaped_neutral` bootstrap-t gates, n_κ = 10⁶ (full sample), sealed κ step, R = 3200, re-sharding under the 1500 s limit, 2098xxxx calibration seeds — *before* any new V7-R transport or calibration result | **autonomous** as design and pre-registration (gate-design change pre-registered before new results) | margins 0.6527/0.7127, the 0.02/0.04/0.9 clauses, B = 1999, the surrogates, the replicate design and every other row stay unchanged |
| Code: α parameter for the bootstrap-t gates (default stays 0.04 for the archived suites so dd976c5a re-evaluates bit-identically), κ-solve step + sealed κ record, calibration parts consume the record, re-sharding constants, new seeds; tests pin the new OC table | **autonomous**, reviewed (Codex) before any host step | implementation of the amendment |
| Sandbox rehearsal of the new calibration at small n (smoke) and the offline OC tables | **autonomous** | exploratory, never evidence |
| Host execution: rehearsal fixture regeneration if hashed sources change, κ step, calibration shards, new V7-R transport shards at a fresh qualification base | **autonomous in principle, but scheduled by the orchestrator under the operator's disposition A**: "V7-R calibration repair before another qualification; an independently reviewed integration route precedes expensive qualification" | the host run is part of the next qualification, not a free-standing job |
| **V7-R2's role relative to the failed V7-R:** whether a V7-R2 pass stands in V7-R's slot of merge gate 14 (a)6 and discharges V2-NUM / V1-MUST-030 in the ledger | **operator proposal** (blocking review finding F1) | V7-R's outcome rule forbids a further base; the operator reserved scope changes and ruled for V6-F that a replacement row does not authorise replacing the failed one or merging |
| **V7-R2 host budget** ≈ 33.9 h (7200 variant ≈ 44 h; (B) ≈ 16.5 h; (C) ≈ 12.5 h) against the ≈ 9.5 h stated in IR-DEFAF4 | **operator proposal** | the budget told to the operator is superseded by a factor ≈ 3.6 |
| Alternatives (B) exact paired certification + R < 3200, or (C) deterministic certification (section 3.5) | **operator proposal** | (B) decreases R ("never decreased", 14 (h)4 / 14a (b)3) and both evaluate a clause analytically (the withdrawn "construction reading", 14a (b)4); thresholds unchanged |
| Applying the new α (or any new calibration) to the **already observed** dd976c5a gate data to turn V7-R into a pass | **operator proposal only — and the planner recommends against it** | it would be a design change chosen after the data were seen; the operator directive requires pre-registration before new results |
| Any change of the clauses themselves (e.g. adopting the "construction reading" of the gate-alone clause, or dropping the gate-alone clause, or certifying the paired gate analytically instead of by Monte Carlo) | **operator proposal** | criterion change; Amendment 14a (b)4 explicitly withdrew the construction reading as a post-rehearsal weakening |
| The open-ended fresh-base sequence 20511004, 20521004 (both skipped, §5), 20531004, … reaches 20981004 / 20991004 (the 2098/2099 families) after 47 steps | **autonomous** clarification in the amendment (the sequence skips 2098xxxx and 2099xxxx) | bookkeeping, prevents a future collision |
| Merge of PR #82 / the merge-gate reading of Amendment 14 (a)6 (V5, V6, V10-A also fail) | **operator** | outside this plan; a repaired V7-R alone does not make the task mergeable |

## 5. Seeds

**Consumed or reserved (from the plan file, the archives and the sibling planner):**
20261004 … 20401004 and 2041xxxx (earlier tasks); 20421004 (lv5) and 2043xxxx (V3-005A rehearsals);
20441004, 20451004, 20471004 with r_index blocks (+1000·r + k); lv5c 20481004 and hr5c 20491004 with their
r_index blocks (lv5c r_index ≤ 20 → ≤ 20501004 + k); rehearsal family 20505000–20509999 and its shard seeds
base + 16000 + k (20521000 + k), i.e. 20505000–20525999 as used; E-shape 20509500 + k; V7-R calibration constants
SEED_FIT 20505000, SEED_BOOT0 20506000 (κ-solve bootstrap base 20505999), SEED_DATA0 / CAL_SEED 20507000
(SeedSequence([20507000, s, p])), SEED_KAPPA 20508000/20508001; declared fresh bases lv 20511004 and hr 20521004
with r_index blocks ≤ +20 999 (≤ 20532003 / 20542003); elastic-attribution plan family 20990000–20999999.
A grep of plans, decisions, scripts and tests finds no integer seed in 20980000–20989999 (only digit runs inside
floating-point fixture values).

**New family for the repaired calibration — 20980000–20989999:**

| Range | Use | Notes |
|---|---|---|
| 20980000 | root of the κ-solve streams: SeedSequence([20980000, surrogate_idx, kind_idx, part]) for data, [20980001, …] for the bootstrap | kind 0 = single, 1 = pair; parts of the 10⁶ sample |
| 20982000 | root of the calibration Monte Carlo: SeedSequence([20982000, surrogate_idx, part]) (data) and [20982001, …] (bootstrap) | replaces 20507000; parts per section 6 |
| 20984000 | Gaussian / Amendment 13 shard (sec_p, nuclear_local), if re-executed | unchanged rule, new seed only because the run is new |
| 20985000–20988999 | reserved for a second calibration attempt **only if** an implementation defect voids the first: the defect must be found by an independent (Codex) review, shown by a reproducing test that does not depend on the calibration counts, the first attempt's counts are reported in full next to the second, the operator is notified, and the use is declared in an amendment before any second-attempt result | never for a statistical miss |
| 20989900–20989999 | planner/sandbox diagnostics (this document used 20989901–20989903, 20989911–20989913, 20989990) | exploratory, never evidence |

The surrogate definitions are **held fixed**, not re-drawn: the archive fit is re-computed with its frozen seed
20505000 and must reproduce k = 5.070278501518333, π = 0.21666666666666667 (I reproduced the κ of record bit for bit,
which implies the fit); the rehearsal pool is the committed sidecar pair (sha256 5d28a9b3… / d4f7d749…, identical at
dd976c5a and at HEAD). Reusing the fit seed reproduces a definition, it is not a re-draw of an observed sample.

**Transport seeds (Revision 2).** "20511004 if it is still unobserved" was the wrong criterion. The declared fresh bases
collide with the rehearsal seeds actually used (largest used rehearsal seed 20525999; the C8 rehearsal family uses
base 20505000 + 1000·r + k): lv 20511004 + 1000·r + k = 20505000 + 1000·(r + 6) + (k + 4), so lv r_index 0–14 lie
inside the used range (exact collisions wherever a used rehearsal seed has r′ = r + 6, k′ = k + 4); hr 20521004 itself and its r_index 1–4 lie inside it (= 20505000 + 1000·(r + 16) + (k + 4)).
Both are recorded as skipped and never used. V7-R2 uses the re-qualification bases of the integration route,
**lv 20531004 / hr 20541004** (20531004 > 20525999 for every r, k): V7-R2 shards (r_index 16) **20547004 … 20547051**
(k = 0..47), reference k = 48 → **20547052**. V7-R2 shares the lv base with the elastic re-qualification (different
r_index) — a sequencing question, not an ownership one.
Cross-base identities (bases 10 000 apart): lv 20531004 r ≡ hr 20541004 (r − 10) at equal k. In particular the V7-R2
seeds equal hr 20541004 r_index 6, k = 0..48 → **hr r_index 6 at 20541004 is reserved and must stay unused**; lv
r_index 11–13 equal hr r_index 1–3 at equal k (as in lv5c/hr5c, Amendment 14 (j)) and are handled there.
The amendment adopts the V6-F checks (v3-005c-v6f-followup §7) verbatim, seed by seed: V7-R2 seeds are disjoint from
(a) `consumed_seed_set()` and every evidence seed used up to then; (b) every rehearsal seed actually used (not the base
range); (c) the never-run original V6 seeds 20493004 + k, 20494004 + k, 20501004 + k; (d) every other row of the
executing task (elastic re-qualification rows, V6-F if run); (e) one another (shards, reference, calibration and κ
SeedSequence roots are distinct entropy). A test asserts (a)–(e) and the hr r_index 6 reservation.

## 6. Cost, sequence, sandbox versus host, interfaces, risks, acceptance tests

### 6.1 Rates used (slow host)

V7-R (dd976c5a layout) transport shard of 4.5·10⁶ histories: 1262 s (was 755 s; factor 1.67). Calibration Monte Carlo: 4.86 s per
repetition at 7200 replicates on the fast host (14a (b)5) → 8.1 s slow; per-repetition cost scales with the replicate
count (all bootstraps are per replicate/pair): 12.2 s at 10 800 replicates, 16.2 s at 14 400. These calibration rates
are **extrapolated** with the transport factor 1.67; one short calibration part (e.g. 20 repetitions) is measured on
the slow host before the part sizes are frozen. Per-shard fixed cost
(fit, record load) ≤ 70 s. κ solve: ≈ 1.29 ms per pair and 0.65 ms per single replicate on the slow host (sandbox
1.0 / 0.5 ms, scaled by the measured host/sandbox ratio). Planning target ≤ 1100 s per step (27 % headroom under 1500 s).

### 6.2 Host cost (one heavy job at a time, single process)

| item | V7-R2 (5400/0.005, recommended) | 7200/0.005 | (B) exact paired + R 800 | (C) deterministic |
|---|---|---|---|---|
| rehearsal fixture regeneration (hashed sources change) | 2 × 1262 s = 0.7 h | 0.7 h | 0.7 h | 0.7 h |
| `v7r-kappa` (n_κ = 10⁶ per kind, 2 surrogates) | 6 shards ≈ 645 s → 1.1 h | 1.1 h | 1.1 h | — |
| calibration Monte Carlo | R 3200: 78 shards × 83/82 reps (≈ 1080 s) → 21.7 h + 1.5 h | 102 shards × 63 reps → 28.8 h + 2.0 h | R 800: 20 shards → 5.4 h + 0.4 h | power only, R 400 at κ = 1: 10 shards → 2.7 h + 0.2 h |
| exact-OC certification step | — | — | seconds | seconds |
| `sec_p`/`nuclear_local` Gaussian (A13) shard | ≈ 0.2 h | 0.3 h | 0.2 h | 0.2 h |
| V7-R2 transport (48 or 64 shards × 2.25·10⁶ + 1·10⁶ reference) | 48 × ≈ 631 s + 280 s ≈ 8.5 h (+ ≈ 0.2 h for 24 extra shard start-ups) | 11.5 h | 8.7 h | 8.7 h |
| combines (calibration, `v7r-combine`) | minutes | minutes | minutes | minutes |
| **total** | **≈ 33.9 h** | ≈ 44 h | ≈ 16.5 h (22 h at R 1600) | **≈ 12.5 h** (+ small end-to-end MC) |
| P(all calibration clauses) r = 1 / r = 0.32 | 1.000 / 1.000 | 1.000 / 1.000 | 0.86 / 1.000 | ≈ 0.99 (power MC only) |
| operator items | merge-gate role + budget | merge-gate role + budget | + clause evaluation, R decrease | + clause evaluation |

Revision 4: the 5400/0.005 column is the design selected by the measured §3.8 rule on the current pool (measured joint
power at κ = 1: 0.9575 / 0.9600, CP95 lower 0.9369 / 0.9399 from 400 repetitions each; the power clause then holds with
P ≈ 1.000); the 7200/0.005 column applies only if the re-measurement on the rebuilt pool misses 378/400.

Sandbox cost of the preparatory work: minutes (OC tables, CI-scale smoke of the new steps at R = 20, n_κ = 2·10⁴), plus
the §3.8 joint-power measurement (single process, ≈ 2.2 h per measurement of the 5400 design, ≈ 3 h for 7200; run
once now on the current pool and once more on the rebuilt pool).

### 6.3 What can run in the sandbox and what must run on the host

* **Sandbox (planning, tests, never evidence):** the exact OC tables and the amendment's numbers; the direct-coverage
  diagnostics of section 2 and the §3.8 pre-freeze measurements (design inputs); unit/CI tests of the new code (α parameter, κ step record, part splitting, seeds, combine);
  CI-scale smoke runs of the calibration; reproduction checks (κ of record bit for bit).
* **Host, at the exact reviewed SHA (evidence):** the rehearsal-fixture regeneration (required by the digest rule
  whenever hashed sources change), the V7-R2 transport shards and reference (warp-cpu float64, attested partials),
  the `v7r-kappa` step, every calibration part and the combine, the Gaussian/A13 shard, `v7r-combine`. Even the
  deterministic certification step of (B)/(C), although its numbers can be computed anywhere, is evidence only as an
  attested host step at the evidence SHA. A sandbox computation can never replace a host step of record.

### 6.4 Sequence

1. (sandbox, autonomous) Record the diagnosis next to the V7-R failure; write Amendment 19 with section 3.7 content;
   decide V7-R2 vs (B)/(C) — (B)/(C) only after an operator decision (send the proposal with the cost table first if
   the 33.9 h budget is not available; in any case V7-R2's merge-gate role and budget go to the operator before host execution).
2. (sandbox) Implement: per-gate α parameter (default 0.04 keeps the dd976c5a/lv5b archives bit-identical), 48 × 2.25·10⁶
   V7-R2 layout, `v7r-kappa` step + sealed record (`ionmc-v7r-kappa-record-1`: surrogate, kind, part, n, seeds, sha256 of
   the per-row arrays, κ per margin after combine), calibration parts consuming the κ record, new part sizes, seeds,
   tests; Codex review; CI.
3. (host, after the operator decision on V7-R2's role and budget) Rehearsal fixture regeneration at the reviewed head (0.7 h) → **§3.8 pre-freeze measurement and design selection (sandbox, before the amendment's numbers are frozen)** → κ step (1.1 h) → calibration parts (23.2 h,
   can be interleaved with nothing else) → calibration combine. **The calibration verdict is known before any V7-R2
   transport is observed;** if it fails, V7-R2 fails and no transport is spent.
4. (host) V7-R2 transport shards + reference (48 × 2.25·10⁶, ≈ 8.7 h) → `v7r-combine` (gates) → row verdict.
5. Archives, README table row, decision 0041 outcome; the dd976c5a failure row stays as recorded.

### 6.5 Interfaces (new or changed)

* `v7r.pair_gate(ing, kappa, *, region, alpha_tost)` / `v7r.single_gate(…, alpha_tost)` — explicit α, default 0.04.
* `v7r.V7R2_SHARDS = 48`, `V7R2_SHARD_N = 2_250_000`, `V7R2_ALPHA_TOST = 0.005`, `V7R2_KAPPA_N = 1_000_000`, `V7R2_CAL_SEED = 20982000`,
  `V7R2_KAPPA_SEED = 20980000`; `cal_part_sizes(3200, parts=39)`.
* Step `v7r-kappa-s{k}` → partial `v7r-kappa-<surrogate>-<kind>-p{k}.npz` + sealed JSON; step `v7r-kappa-combine`
  → `v7r-kappa.json` (κ per margin, coverage at κ with SE, n); calibration parts refuse to run without it and record
  its content hash.
* `calibration_combine` gains a check that every part used the same κ-record hash; `calibration_evaluate` unchanged.

### 6.6 Risks

* **Joint power at α = 0.005** was modelled in Revisions 1–3 (single-gate shrink fitted to α = 0.04 data); Revision 4
  measured it directly on the current pool (§3.8: 0.9575 / 0.9600, CP95 lower 0.9369 / 0.9399, 400 repetitions per
  surrogate), so the risk that the true joint power is < 0.93 is retired for the current pool. Residual risk: the
  rebuilt pool after the elastic fix may differ; handled by the pre-declared re-measurement and rule of §3.8 (fallback
  7200/0.005, then operator).
* **Real estimator off-nominal:** a correct estimator whose true paired coverage is 0.69–0.695 passes the stricter
  gate with 0.86 / 0.61 (5400/0.005) — the price of the smaller α; the real dd976c5a value was 0.685.
* **Host time:** 33.9 h at the slow rate; a further slowdown pushes calibration parts toward 1500 s — the 1100 s target
  leaves 27 % headroom under 1500 s (and is 8 % below the 1200 s effective limit; transport shards ≈ 631 s).
* **Base sharing:** V7-R2 and the elastic re-qualification share lv 20531004 (different r_index) — sequencing only;
  hr r_index 6 at 20541004 reserved (§5).
* **Rehearsal regeneration after the elastic fix** changes the rehearsal surrogate; the power clause is the binding
  one — handled by the pre-freeze rule of §3.8.
* **Interpretive risk:** the amendment is written after the dd976c5a gate data were seen; mitigated by an unobserved
  new transport base, a stricter α, and the calibration verdict being produced before the new data.

### 6.7 Acceptance tests (for the implementation and the amendment)

1. Bit-identity: the dd976c5a calibration records and the lv5b/lv5c gate results re-evaluate identically with the
   default α = 0.04 and the old layout (archived sidecars).
2. κ reproduction: a **legacy mode** of the κ step (n = 2·10⁴, old seeds, the old 3:1 split for the single kind) reproduces
   the five κ of record of both surrogates bit for bit (as shown here in the sandbox); the V7-R2 mode (10⁶, full
   sample) is a separate code path tested for its SE.
3. Exact OC pins: acceptance band 3615–3762 of 5400 at α = 0.005; OC 0.0050 / 0.0050 at the margins; power 0.969 at
   0.6827 (paired); the R = 3200 clause table (109 / 50 / 2909) unchanged.
4. κ accuracy: the combined κ record reports coverage SE ≤ 0.0005 for every margin; a planted-error test (κ shifted
   by its 2σ) shows the expected OC change.
5. Seeds: a guard rejects any calibration/κ seed outside 20980000–20984999 and any overlap with the consumed list of
   section 5; the transport seeds pass checks (a)–(e) of §5 and hr r_index 6 at 20541004 is asserted unused.
6. Fail-closed: missing κ record, a part with a different κ hash, unequal part counts, or R ≠ 3200 → calibration fail.
7. Step timing: every new step ≤ 1100 s in the CI-scale extrapolation at the slow-host rate.
