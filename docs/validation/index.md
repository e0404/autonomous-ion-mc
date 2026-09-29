# Validation

This section documents software, numerical, physical, statistical and
cross-backend validation.

- The [requirement ledger](requirement-ledger.md) maps every MUST requirement
  to its status, evidence suites and tasks (canonical source:
  `validation/requirement-ledger.json`).
- `validation/release-plan.json` defines the evidence suites; it is a
  **draft** until the calibration task freezes tolerances and performance
  targets. A draft plan is not acceptance evidence.
- Local scientific validation results are recorded under the exact-SHA
  validation gate and, for release qualification, under `validation/generated/`
  (ignored by Git) and the protected experiment state directory.
