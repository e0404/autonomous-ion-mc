# Validation

This section documents software, numerical, physical, statistical, and
cross-backend validation.

## V3-013: explicit stopping-model selection (status 2026-11-18)

Decision [0044](../generated/decisions/0044-stopping-source-selection.md); acceptance plan
`validation/plans/v3-013-acceptance.md` (frozen at 7c139ca1); model documentation in
[stopping power](../physics/stopping-power.md). Not a clinical claim.

| rows | content | status |
|---|---|---|
| A.1 | regression and identity of the default (R2 digests, R1-ion, ID-DEF, ID-NIST, ID-GUARD, ID-DIAG, ID-USER) | implemented in S1 to S3; tests pass in the logged runs of the step commits |
| A.2 | consumer consistency with a synthetic test-only model (CC-R, CC-T, CC-LET, CC-RES, CC-ID, CC-NEG, CC-HASH) | implemented in S2b; tests pass in the logged runs of the step commits |
| A.3 | fail-closed domains and roles (DOM-1 to DOM-7, ROLE-1 to ROLE-3) | implemented in S2a to S4; tests pass in the logged runs of the step commits |
| A.4 | persisted identity (PER-1 clean-process round trip, PER-2, PER-3 with negative control) | implemented in S2a to S4; tests pass in the logged runs of the step commits |
| A.5 | consistency checks of `icru90-water-v1` (CC-NODE, CC-INT, CC-I2, CC-SPAN); consistency, never evidence | implemented in S4; CC-NODE, CC-I2 and CC-SPAN pass; CC-INT and V-C3 aggregates committed in `validation/results/stopping/v3-013/` |
| A.6 | V-C1 (He 75/125/200 MeV/u against TOPAS opt4 EM-only, shared-lineage agreement), V-C1-A (report-only), V-C2, V-C2-D (selection reaches transport), V-C3 (report-only) | V-C1, V-C1-A, V-C2, V-C2-D pending: no evidence run exists; V-C3 aggregates committed |

Seal: the seal material (V-C1 TOPAS bundles, ionmc sample definitions, evaluator, environment
digest record from a controlled host run, committed predictions, frozen manifest of 80 sealed paths)
is committed as the **seal candidate**; the seal is the signed review of its exact commit, which is
pending. ionmc samples of both models are produced and committed before any TOPAS run. No outcome
changes the default or any registered spec; a failed row is recorded as `failed:<row>` and the spec
is unchanged. The V3-008 verdicts (row I3 FAIL, I3-FU-1 INCONCLUSIVE) are retained.

### What is NOT claimed

- No claim that any model is TOPAS-equivalent; V-C1 is a shared-lineage agreement claim for one
  configuration.
- No physical validation of ICRU 90 stopping by V-C1, CC-NODE or CC-I2; no claim for protons with
  `icru90-water-v1` beyond consistency (no proton row is pre-registered).
- No carbon or oxygen alternative; the default's C/O low-energy limitation and the I3 FAIL stand.
- Measured He/C/O ranges, measured lateral data and measured C/O stopping at therapy energies are
  evidence gaps (experimental validation deferred by the operator directive of 2026-11-07), not
  relabelled.
- No evidence transfers between models.
