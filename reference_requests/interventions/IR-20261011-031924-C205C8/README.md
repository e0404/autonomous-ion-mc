# Human Intervention Request IR-20261011-031924-C205C8

**Task:** V3-005B
**Category:** experiment-boundary
**Status:** resolved (2026-10-11T18:08:29Z; operator response in `response.json`, recorded through `resolve_human_intervention`; the README status line was reconciled in the task branch after the response, the request tooling does not rewrite it)

## Summary

V3-005C row V6 / V6-sens: the pre-registered history-count rule (Amendment 14 (g)) and the V6-sens point-set reading (Amendment 17 (c)6) are unsatisfiable at any feasible history count, found at the C8 rehearsal (base 20505000, never evidence). Decision requested before any V6 evidence shard is run: run as frozen (row mechanically report-only, ≈ 20 h host time), amend the two readings before evidence with the contrary reading reported, or defer V6.

## Why autonomous resolution is not possible

Rehearsal at the reviewed head c847cc7559a0daf3fefb89304a429802ca52c31e (Codex REVIEW-e1853ef7 pass), seed base 20505000, warp-cpu float64, 1e6 histories per shard, 20 batches, E0 = 177.84075737 MeV, the frozen Gottschalk 2014 beam.

(1) N_V6 rule (Amendment 14 (g): "smallest number of 1e6-history shards giving rehearsal-estimated relative SE <= 3 % at every r >= 1 cm point whose measured dose is >= 1e-4 of the r = 0 maximum, and <= 0.5 % at every gating r = 0 point"; Amendment 17 (c)2 makes every non-missing r = 0 depth a gating point). Shard v6-s0 (RUN-20261011T030353Z-0b25db97, 184 s): halo clause (153 points) worst relative SE 28.1 % at (z 22.0 cm, r 3 cm) -> 88 shards. r = 0 clause (34 points): satisfied with <= 18 shards up to z = 21.8 cm; needs 116 shards at z = 22.0 cm, 651 at 22.2 cm, 8910 at 22.4 cm, and is impossible at 22.6 cm (model dose zero in 1e6 histories; measured 4.3e-4 MeV/g). The rule therefore defines no finite N_V6.

(2) V6-sens condition (Amendment 14 (g) / Amendment 17 (c)4,6: |ratio - 1| + 2 sigma_delta < 0.15 at every non-missing r >= 1 cm point, common random numbers). Shards v6-sens-nominal-s0 (RUN-20261011T030829Z-537a9857) and v6-sens-st-hi-s0 (RUN-20261011T031141Z-88a536e0), 185 s each, same seed: of the 261 gating halo points, 32 have zero model dose at 1e6 histories (z >= 20 cm, measured 6e-6 to 4e-5 of the r = 0 maximum), so the ratio is undefined there and the condition fails at any shard count; the worst finite paired SE needs 345 shards per variant. On the sizing subset (measured >= 1e-4 of the maximum, 153 points) the worst paired SE 0.42 gives k = 32 shards per variant (9 variants, ~15 h host time); the median point needs 1 shard.

Consequence under the frozen text: V6 cannot be sized by its rule, and V6-sens cannot hold, so V6 would be report-only (Amendment 17 (c)7) whatever the data show, after ~20 h of host time (88 shards gating + 288 V6-sens shards at ~185 s, plus the hr5c CUDA shards). The operator directive of 2026-10-10 forbids relaxing acceptance thresholds; the r = 0 2 % and halo 15 %/90 %/10 % tolerances are not in question here, only the SE sizing rule and the point set on which the sensitivity ratio is defined.

Contrary-evidence note (recorded in the progress log): while selecting the rule's points the orchestrator printed measured and model mean doses side by side for the rehearsal shard; no residual, interval or verdict was computed and nothing was changed as a consequence; the evidence base 20481004 is unobserved.

No V6 evidence shard has been run. All other lv5c/hr5c rows (V5, V7-R, V8-nuc, V10/V10-A, V11, R1-D2, E1-B, E-shape) are unaffected and their rehearsal continues; the V10-A deterministic preview failure will be raised separately after C9/C10 with the recorded outcomes.

## Exact requested input or action

Choose one (or state another):
(A) Run V6 exactly as frozen: N_V6 = 88 (the halo clause; the four deepest r = 0 points recorded as "SE criterion unmet at the frozen size", their 2 % point criterion still evaluated), V6-sens with k = 32 shards per variant; the V6-sens condition then fails at the 32 zero-dose points by construction and V6 is report-only, recorded. ~20 h host time.
(B) Amend before evidence (Amendment 14a, contrary reading reported next to every verdict): N_V6 = 88 from the halo clause with the r = 0 SE clause restricted to depths where it is finite (z <= 21.8 cm), the 2 % point criterion unchanged at every r = 0 depth including the four deepest; the V6-sens condition evaluated on the sizing subset (measured >= 1e-4 of the r = 0 maximum, 153 points) with k = 32, the 32 zero-dose points reported as undetermined. Same host time; the row can gate.
(C) Defer V6 to a later task with the record "pre-registered sizing rule unsatisfiable at the rehearsal"; no V6 shards run; the merge gate (a)6 treats V6 as "deferred, recorded".
Also: if (A) or (B), confirm k = 32 for V6-sens (sized on the sizing subset) rather than k = N_V6 = 88 (which would be ~41 h).

## Decision or task depending on this request

Not specified.
