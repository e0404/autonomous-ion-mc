# V3-008 pre-freeze check P11: exact-Omega mapping (record, 2026-11-01, base 3d2547a)

Nothing is frozen by this note. The I11-rel / I11-abs numbers below are SEEN construction/calibration results, report only, not qualification evidence.

**Status and chronology.** Mapping, tests and the first I11-rel/I11-abs evaluation were produced together at 070b07f0; the values are SEEN construction/calibration results, not qualification evidence; a later rerun cannot make these points independent. The pre-registered order (mapping and tests committed before any evaluation) was not kept. I11-rel stays the gate as pre-registered (2 %, STOP on failure, tolerance never changed) but is evaluated pre-freeze, not blind; the freeze carries it as a disclosed gate. It is not asserted in any unit test. Genuinely held-out evidence for the charge dependence of lateral spreading is V5-SS (if its run condition holds) or measured data only, never I11-rel. (Corrected 2026-11-01 for Codex REVIEW-b469dff8. The 00af9a6c text read "V1-held/V5-SS/measured data", but V1-held has no lateral statistic and provides no lateral evidence. Measured data: EVIDENCE GAP (`p6-lateral-data-search.md`). V5-SS class: P10 §7.5.)

**Correction (Codex REVIEW-42ab68dd).** I11-rel is at EQUAL VELOCITY. The 070b07f0 driver evaluated the proton at the ion's numeric MeV/u; the corrected driver uses E_p = E_u A m_p / M (`equal_velocity_proton_energy`), for both the proton model and the proton Moliere-Hanson reference. Regression test `test_i11_rel_proton_is_at_equal_velocity`: beta_p(E_p) = beta_ion(E_u) to 1e-12 at every grid point. The tables below are regenerated.

## Files (bind in P9; sha256; regenerated after the REVIEW-42ab68dd correction)

| File | sha256 |
|---|---|
| `validation/scripts/ion_em/exact_omega.py` | `f5bc5f2c684e060f550f26b1d5e621e2ef5ecf01d4e9c35ed51d326d928c41b3` |
| `tests/ionmc/test_ion_em_exact_omega.py` | `12cbdba28530db9561ac469acba5d45f7533690f5587dd74d44b2143b542baee` |
| `validation/scripts/ion_em/moliere_hanson.py` (reference, unchanged) | `5a99b3b8748a6865584fcbe1337206222ddd7fe197e95f25656eac2324acf164` |

## Implemented equations (decision 0043 §3, quoted)

```
T_ion(pv, p1v1; material) = f_dM(pv_eq, p1v1_eq) (z E_s / pv)^2 / X_S
p1v1_eq / pv_eq = p1v1 / pv          (the memory term 1 - (pv/p1v1)^2 is unchanged)
pv_eq is the root of   b_p(pv_eq, p1v1_eq; material) = b_ion(pv, p1v1; z, M, material)
b = ln Omega + delta_b_Fano
```

b of the slab from p1v1 down to pv on the projectile's CSDA path (t = R(E1) - R(E)):
b = ln chi_c^2 - <ln chi_a^2> - ln 1.167 + <delta_b>, chi_c^2 = K z^2 sum w Z^2/A int dt/(pv)^2,
averages with weight 1/(pv)^2 (Gauss nodes of the reference, 64 nodes inside b, 400 along the path).
All formulas are the reference's (`moliere_hanson.py`).

**p1v1_eq choice.** Exactly as §3 prescribes: p1v1_eq = pv_eq * (p1v1/pv), i.e. the equivalent proton
keeps the ion's own memory ratio p1v1/pv. It is NOT the proton with the ion's initial velocity.
The proton slab is p1v1_eq -> pv_eq on the proton path.

Evaluation: per-step root solve (Brent on ln pv_eq; bracket = proton table; fail-closed ValueError
without a sign change). Tabulation variant `PvEqTable` (bilinear, axes ln pv and ln(p1v1/pv), a shear
of the decision's (ln pv, ln p1v1) plane) is implemented as the allowed alternative; its measured
f_dM error is 5.3e-6 on a He/C test domain (bound 1e-3). P9 should bind the root solve unless the
tabulation is chosen explicitly.

The proton f_dM is the numpy twin `moliere_hanson.fdm` (coefficients asserted equal to
`ionmc.physics.em.FDM_COEFFICIENTS`). There is no `src/ionmc/physics/mcs.py`; the project f_dM is
the Warp function in `em.py` (not modified, not called).

## Grid and results (water, Bethe tables, Fano reference; the exact I11 grid)

He 50/150/250, C 80/290/430, O 100/350/430 MeV/u x x/R in {0.1, 0.25, 0.5, 0.9} (36 points).

**I11-map residual |b_p(pv_eq) - b_ion|:** max 1.8e-14 on the grid (unchanged) (tests: 4.8e-14 on grid plus 7 extra
depths per energy); limit 1e-10. z = 1, M = m_p: pv_eq = pv bit-exact, f_dM equal (<= 1e-12 test, 20 points).
Omega identity (scaled synthetic path, z = 2, 6, 8, Fano and Bethe): <= 1e-10.

**I11-rel max |ratio - 1|** (equal velocity; SEEN, pre-freeze; gate 2 %):

| Ion | max | where (MeV/u, x/R) | min |
|---|---|---|---|
| He | 1.71 % | 50, 0.1 | 0.69 % |
| C | 1.82 % | 80, 0.9 | 0.84 % |
| O | 1.75 % | 100, 0.9 | 0.86 % |

(Superseded numeric-MeV/u values: He 1.73 %, C 1.83 %, O 1.76 %.) All values are positive; margin to the 2 % gate is 0.18 pp (C).

**I11-abs max |theta_model/theta_MH - 1|** (informative; unchanged by the correction): He 4.65 % (250, 0.1); C 5.31 % (430, 0.1);
O 5.32 % (430, 0.1); range +1.85 to +5.32 %.

pv_eq/pv on the grid: He 0.57-0.69, C 0.18-0.36, O 0.13-0.27. Proton p1v1_eq reaches 0.8 GeV (He),
2.5 GeV (C, O) kinetic energy: the proton table must extend to about 3 GeV (this check used 8 GeV).

## OPEN

1. Transport implementation needs proton cumulative b-tables out to ~3 GeV (above the 500 MeV/u
   project table); to be decided before P9.
2. I11-rel C margin is small (0.18 pp); STOP rule applies if freeze evaluation exceeds 2 %.
3. Warp twin of the root (I11-map "python twin vs Warp") is not implemented here.
4. Inner slab nodes 64 vs reference 400: b differs by < 1e-7.

## P11b. f_dM domain extrapolation (addendum 2026-11-01, r6, Codex REVIEW-b469dff8 finding 4)

**Chronology.** Computed by the planner on 2026-11-01 at head `b3eef90b`, after the I11-rel/I11-abs values above had
been seen. Scratch scripts are outside the repository: `domain.py` and `brange.py` in the session scratchpad
`v3-008/p11b/`. They use the committed `exact_omega.py` and `moliere_hanson.py` unchanged. Every number below is
**SEEN pre-freeze**, an aggregate, and not qualification evidence. The mapping, its scripts and tests are unchanged, so
the hashes above stand.

**What f_dM is a function of (FACT, `src/ionmc/physics/em.py` `fdm`, `moliere_hanson.fdm`).**
- f_dM = 0.5244 + 0.1975 L1 + 0.2320 L2 − 0.0098 L1 L2, with L1 = lg(1 − (pv/p1v1)²) and L2 = lg(pv / MeV).
- It is an explicit function of (pv, p1v1), not of b or Ω.
- Its qualified domain is protons with E ≤ 300 MeV (decision 0043, Context; arXiv:0908.1413), i.e. p1v1 ≤ pv(300 MeV) = 527.3 MeV.
- At fixed L1, ∂f/∂L2 = 0.232 − 0.0098 L1 ≈ 0.23–0.25 per decade.

**(i) b-range overlap.**
- At fixed memory ratio r = p1v1/pv, b_p is monotone in pv; the root bracket relies on this.
- Hence "b_ion within the proton-qualified b range at that r" is equivalent to p1v1_eq ≤ 527 MeV. Matched b does not bring the state into the domain.
- On the 36 I11 end states:
  - b_ion = 11.72–17.06;
  - the largest qualified proton b at the same r is 12.61–14.95;
  - in domain: 11/36 (He 50 and C 80 at all x/R; O 100 at three of four x/R).
- Protons 50–300 MeV over the I11 x/R cover b = 10.60–14.95.

**(ii) Extent of the extrapolation (θ²-weighted along the I11 integration paths, 60 nodes).**

| Ion, E (MeV/u) | θ² fraction from p1v1_eq > 527 MeV | max kinetic energy of p1v1_eq | max extrapolation in lg pv |
|---|---|---|---|
| He 50 / C 80 | 0 | 136 / 251 MeV | 0 |
| O 100 | 0.09–1.0 (by x/R) | 310 MeV | 0.013 |
| He 150 / He 250 | 1.0 | 460 / 807 MeV | 0.16 / 0.37 |
| C 290 / O 350 | 1.0 | 1486 / 1876 MeV | 0.59 / 0.68 |
| C 430 / O 430 | 1.0 | 2471 / 2482 MeV | 0.78 |

- f_dM on these paths: 0.38–1.25.
- Over 0.78 decades, f_dM's own L2 term changes f by ≈ +0.19. This is the physical growth of b, not the error.
- The question is whether the linear-in-L2 fit keeps tracking the Molière b-dependence beyond 300 MeV.

**(iii) Bound against Molière–Hanson proton theory (water, Fano; same reference module).**
θ_p,model/θ_MH,p − 1 for protons:

| E_p (MeV) | x/R 0.1 | 0.25 | 0.5 | 0.9 |
|---|---|---|---|---|
| 100 | +2.52 % | +2.75 % | +2.62 % | +1.25 % |
| 300 (edge) | +3.95 % | +4.01 % | +3.80 % | +2.53 % |
| 1000 | +4.80 % | +4.76 % | +4.52 % | +3.41 % |
| 3000 | +5.44 % | +5.30 % | +5.01 % | +3.95 % |

- Intermediate rows (200, 430, 600, 2000 MeV) are monotone between these.
- The drift beyond the qualified edge, up to 3 GeV, is +1.2 to +1.5 pp in θ. It is smooth and logarithmic, with a slope similar to the in-domain drift of +1.3 to +1.4 pp from 100 to 300 MeV. There is no breakdown.
- This bounds the extrapolation error *relative to Molière theory*: ≤ 1.5 % in θ for the proton-equivalent states used, ≤ 2.5 GeV.
- It is not independent evidence: the bound uses the same theory as the mapping.
- It is consistent with the seen I11-abs maxima (+5.3 % at C/O 430, x/R 0.1), which exceed the in-domain proton offset by about the drift.

**(iv) Options considered and the choice.**
- (a) In-domain reformulation in b alone: impossible for this domain. A b-only g(b, r) built from qualified protons covers b ≤ 14.95, and b_ion reaches 17.06, so it extrapolates in b instead of pv, with no better justification.
- (b) Fail-closed clamp (I0) at the qualified domain: would remove He above ≈ 100 MeV/u and C above ≈ 100 MeV/u, and almost all of O. That is a scope change against V2-ION and decision 0043 §5; not chosen. A scope reduction is the orchestrator's or operator's decision, not the planner's.
- (c) **Chosen:** keep the model unchanged (no change after seen numbers), and:
  1. declare the extrapolation;
  2. bound it as in (iii) and carry it in the declared ±5 % model uncertainty as a separate component (≤ 1.5 %; the quadrature becomes ≈ 4.9 %, decision 0043 §3);
  3. add the report-only row I11-ext (the proton drift table, SEEN) and report the extrapolated θ² fraction per I11 point;
  4. persist the domain status in metadata (I0: the p1v1_eq kinetic energy of the source and the flag `fdm_extrapolated`); the proton b-slab tables must reach ≥ 3 GeV or the run fails closed;
  5. add independent coverage through V5-SS. The 1 mm V5-SS states already map to extrapolated proton states (p1v1_eq 458 MeV for He 150, 1486 for C 290, 1875 for O 350; b_ion 10.3–11.9). The end cases He 250, C 430 and O 430 (786, 2470 and 2482 MeV; cases 28–33, ≈ +1 h REF) reach the largest extrapolation. This coverage holds only subject to the V5-SS STOP (P10 §7.4) and its class: shared screened-Coulomb lineage, independent of f_dM.
- I11-rel cannot test the extrapolation: at C/O the equal-velocity proton (up to ≈ 430 MeV) is itself extrapolated.

**OPEN after P11b.**
- Item 1 of the OPEN list above (proton b-tables to ~3 GeV) is now a fail-closed I0 requirement.
- Whether V5-SS can provide this coverage depends on P10 §7.4.
