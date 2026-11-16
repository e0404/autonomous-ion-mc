# ADDENDUM 1 to IR-20261116-204824-29B867 (2026-11-16, V3 orchestrator)

Two corrections to the request text, found by the signed Codex review REVIEW-e0d9bcd0 of the commit that recorded this
request (ee9ee73d). The options (A), (B), (C) and the declared uncertainty values are unchanged; please read the request
with these corrections.

1. **Spatial extent of the affected region (important).** The request says "The affected region is the last ≈ 0.4 mm
   (C 290 MeV/u) / ≈ 0.3 mm (O 350 MeV/u) of range". That is the region below 10 MeV/u. The region below **20 MeV/u**, where
   the declared deficit is 1.8–12 % against the revised ICRU 73 table and 0.6–5 % against the acceptance reference, is the
   last **≈ 1.45 mm of range at C 290 MeV/u (0.9 % of 163 mm) and ≈ 1.11 mm at O 350 MeV/u (0.7 % of 167 mm)**; below
   50 MeV/u (which includes the 2–20 MeV/u region above; the 20–50 MeV/u band alone has a deficit ≤ 0.8 % against the
   table) it is ≈ 7.5 mm (4.6 %) and ≈ 5.7 mm (3.4 %) [wording corrected after Codex REVIEW-21409d5d]. These are the committed P12 §4
   figures (`validation/plans/v3-008/p12-i3-low-energy-stopping.md`). Decision 0043's capability statement is corrected in
   the same commit as this addendum.
2. **Methodology wording (minor).** The TOPAS/ionmc band means are described as "one-step thin-slab mean energy loss". The
   point values entering the band means are the **fine-step runs** (20 limited steps in the slab); the one-step slab pair
   (t and 1.6 t) enters only the per-point emulation gate and the report-only checks.

The request remains open; integration of PR #84 waits for the operator's answer. Not a clinical claim.
