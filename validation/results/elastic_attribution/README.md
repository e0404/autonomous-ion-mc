# V3-005D elastic attribution — exploratory records (Amendment 19; never evidence)

These directories are committed copies of the outputs that the modules under `validation/scripts/elastic_attribution/` write to `validation/generated/elastic_attribution/` (git-ignored). Each report starts with the banner "Exploratory, report-only (Amendment 19); not evidence; changes no row." and lists its inputs with sha256 and the command that regenerates it.

- `de0/` — D-E0 kernels (K_w, K* as frozen and in the corrected-sign form, G, θ_w, T), reconciliation of the planning expectation, H6 case-file check. `quarantine-data-seen.json` holds shape parameters fitted to V11; it is never cited (Amendment 19 addendum (a)3).
- `de1/` — D-E1 per-dataset shape comparison. The committed `per_dataset.json` has the experimental differential cross-section values (`data_mb_sr`, `rel_err`) stripped: no EXFOR values enter Git; the module reads them from the reference cache. Derived ratios, partial σ_tr integrals, minima diagnostics and the data-seen one-parameter R fits (diagnostic only) are kept.
- `de2/` — D-E2 p-p audit.

Independent review of these records: `validation/plans/v3-005d-elastic-attribution/review.md` (section "D-E0–D-E2 review"). The K* sign transcription error of the plan and the refuted "core 2× too wide" prior are recorded in Amendment 19 addendum (a) of `validation/plans/v3-005-acceptance.md`.
