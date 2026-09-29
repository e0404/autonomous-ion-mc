# Validation

This section documents software, numerical, physical, statistical and
cross-backend validation.

- The [requirement ledger](requirement-ledger.md) maps every MUST requirement
  to its status, evidence suites and tasks (canonical source:
  `validation/requirement-ledger.json`).
- `validation/release-plan.json` defines the 18 evidence suites, their
  acceptance criteria and the numeric performance targets. It was frozen in
  task V2-002 (decision 0039) after the performance calibration recorded in
  `benchmarks/calibration/`; the release qualification evaluates evidence
  against the plan bytes at the freeze commit.
- `validation/references/` holds the committed native input bundles for the
  TOPAS, MCsquare and FRED reference services; raw outputs are archived by
  run ID outside Git.
- `validation/research/` holds condensed research notes that precede
  consequential scientific decisions.
- Local scientific validation results are recorded under the exact-SHA
  validation gate and, for release qualification, under `validation/generated/`
  (ignored by Git) and the protected experiment state directory.
