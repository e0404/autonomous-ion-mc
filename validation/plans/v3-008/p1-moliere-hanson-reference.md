# V3-008 pre-freeze check P1: Molière–Hanson ion reference (record, 2026-11-01)

> **Status (after P1b, 2026-11-01): P1 NOT CLOSED.** With the Fano correction (section P1b) Be, Al
> and Cu pass the plan's ≤ 1.0 % proton validation; **Pb still fails** (max +1.32 % at x/R1 0.97).
> Nothing was tuned. Against the Fano-corrected water reference **no candidate passes 3 %**.
> Sections below "P1b" are the original P1 record (Bethe Z(Z+1), no Fano), kept for the record.

## Files (bind in P9; `sha256sum`, git-independent; updated for P1b)

| File | sha256 |
|---|---|
| `validation/scripts/ion_em/__init__.py` | `7962c4af19580a10f51a4bfc703ce45b67a0f96f010aaaf878aef6919b2ed4ea` |
| `validation/scripts/ion_em/moliere_hanson.py` | `9806f0cb47d00c24f89e52bea3768bddfbb6ae332c703180db012aa4ea42b5f4` |
| `tests/ionmc/test_ion_em_moliere_hanson.py` | `61403e3f3ac2f67e7e4ad8880253e99d7bacadc4d6e796bac42ead6441241133` |

## Definition implemented (θ0 convention)

θ0 is the **projected** Gaussian width of the whole slab (the Hanson width). This is the convention of
`ionmc.physics.em`, where the scattering power is in rad²/mm for the projected angle, and of Gottschalk's
θ_Hanson column. The projectile is (z, M). The slab has mass thickness t along the CSDA path, with
E(t′) taken from the project range table.

- χ_c² = K z² Σ w_i Z_i(Z_i+1)/A_i · ∫₀ᵗ dt′/(pv)², where K = 4π N_A r_e² (m_e c²)² = 0.156915 MeV² cm²/g.
  Bethe's Z(Z+1) form is used.
- χ_a,i² = χ_0,i² (1.13 + 3.76 α_i′²), with χ_0,i = α m_e c² Z_i^{1/3}/(0.885 pc) and α_i′ = z Z_i α/β.
- Compound rule: ln χ_a² = Σ W_i ln χ_a,i², with W_i ∝ w_i Z_i(Z_i+1)/A_i.
- Thick slab: ln χ_a² is averaged along the slab with weight 1/(pv)².
- b = ln(χ_c²/(1.167 χ_a²)), with 1.167 = exp(2C_E − 1). B is the root of B − ln B = b.
- θ0 = χ_c √(B − 1.2)/√2.
- Point nucleus.

The integration is Gauss–Legendre with 400 nodes, clustered at the slab end.

**Constants and sources:**
- CODATA 2018: α, r_e, m_e c², m_p c². N_A is exact (SI 2019).
- Moliere 1948 and Bethe, Phys. Rev. 89 (1953) 1256: 1.13, 3.76, 0.885 and 1.167.
- Hanson et al., Phys. Rev. 84 (1951) 634: the 1.2 offset.
- The forms as written are taken from Gottschalk, arXiv:1804.00022 (eqs. chicsq–Hanson) and arXiv:0908.1413 (eq. thetaHanson).
- Cross-check: χ_0²·1.13 equals Lynch–Dahl's 2.007e-5 Z^{2/3}/p² to within 0.2 %.
- Water composition and X_S are identical to `ionmc.materials` and `ionmc.physics.scattering` (tested).

**[Superseded by P1b] Deviation from the pre-declared definition (decision 0043 §3): the Fano atomic-electron term is NOT
implemented.** I could not source its exact Gottschalk-1993 form in this pass: NIM B 74 (1993) 467 is
paywalled, and Gottschalk's BGware site returned HTTP 502. The electron contribution therefore uses
Bethe's Z(Z+1) form. This is OPEN item 1.

## Proton validation (plan P1: ≤ 1.0 % at every point) — FAILS

**Reference and setup:**
- Reference: θ_Hanson (Molière/**Fano**/Hanson) from `tbl:theta0Single` of arXiv:0908.1413, for 158.6 MeV
  protons. The source is the hash-pinned cached file sha256 `67fb1478…`.
- The reference module uses the ionmc Bethe proton range tables. Each point is evaluated at x = (x/R1)·R1(ours).

**Results**, as θ_MH/θ_Hanson − 1 in %:

| x/R1 | 0.001 | 0.01 | 0.1 | 0.2 | 0.5 | 0.9 | 0.97 |
|---|---|---|---|---|---|---|---|
| Be | −2.73 | +0.78 | +2.97 | +3.50 | +4.22 | +5.01 | +5.40 |
| Al | −0.01 | +1.02 | +1.61 | +1.75 | +1.93 | +2.14 | +2.26 |
| Cu | +0.17 | +0.58 | +0.81 | +0.86 | +0.92 | +0.95 | +0.95 |
| Pb | +0.68 | +0.79 | +0.87 | +0.91 | +1.04 | +1.55 | +1.79 |

**Verdict:**
- Cu passes, with a maximum of 0.95 %.
- Be, Al and Pb fail. Their tests are marked `xfail(strict=True)` with the reason; the tolerance is not loosened.
- The deviation pattern is large and positive for low Z and shrinks with Z. This is the expected signature of
  Bethe's Z(Z+1) electron treatment against Fano's: Gottschalk 2010 notes the Bethe/Fano difference for Be.
- The residual +0.8–1.8 % for Pb, growing toward 0.97 R1, points to a range-table or end-of-range
  difference. Gottschalk used ICRU 49 tables with spline interpolation. This is OPEN item 2.
- **Water** is low-Z (H, O), so the missing Fano term matters for the I11 reference itself.

## Frozen candidate (d) correction (formulas; not fitted)

Notation: M is the ion rest mass, A* = M/m_p, pv* = pv m_p/M, p1v1* = p1v1 m_p/M, and β is the local velocity.

```
f_dM,ion = f_dM(pv*, p1v1*) + c_Ω · ln( A* / F_scr(β, z; material) )
ln F_scr = Σ_i W_i [ ln(1.13 + 3.76 (z Z_i α/β)²) − ln(1.13 + 3.76 (Z_i α/β)²) ],   W_i ∝ w_i Z_i(Z_i+1)/A_i
c_Ω      = (c_2 + c_3 L_1) / (ln 10 · n_Ω),      L_1 = lg(1 − (pv/p1v1)²),  c_2 = 0.2320, c_3 = −0.0098
n_Ω      = d ln R_p/d ln pv − d ln β²/d ln pv   (proton, at E* = E m_p/M, from the same material's tables;
           d ln R/d ln E = E/(S R), kinematic derivatives exact)
```

**Derivation:**
1. Ω = χ_c²/(1.167 χ_a²) ∝ z² t/(β² F). The test `test_z_scaling_identity_at_fixed_velocity` checks the
   exact identity Ω_ion/Ω_p = A*/F_scr at equal β and equal x/R (t_ion = (A*/z²) t_p).
2. f_dM depends on the proton's Ω at fixed x/R only through lg pv, because L_1 is invariant under velocity scaling.
3. Choosing pv_eq so that the proton Ω equals the ion Ω gives ln pv_eq = ln pv* + ln(A*/F_scr)/n_Ω.
4. To first order in that shift, Δf = c_Ω ln(A*/F_scr).

**Finding:** n_Ω is 0.90–1.12 for 50–300 MeV in water, so c_Ω ≈ 0.09–0.11. The decision's estimate (n_Ω ≈ 0.77,
c_Ω ≈ 0.13) was non-relativistic. Its test is `xfail(strict=True)`, so the formula stays frozen.

## Candidate predictions vs this reference (report only; water; Bethe 78 eV range tables)

Values are max |θ_cand/θ_MH − 1| over x/R ∈ {0.1, 0.25, 0.5, 0.9} at the I11 points (standard points and
ends). The full grid, including 50/100/150/200/300/430 MeV/u × {0.1, 0.5, 0.9}, comes from `python
validation/scripts/ion_em/moliere_hanson.py`.

| Candidate | He | C | O | Verdict at 3 % |
|---|---|---|---|---|
| (d) velocity-scaled f_dM + c_Ω ln(A*/F_scr) | +2.82 % | **+3.31 %** (290, x/R 0.1) | **+3.26 %** (350, 0.1) | fails marginally, only at x/R = 0.1; ≤ 2.9 % at x/R ≥ 0.25 |
| (c) slab proxy: PDG/Lynch–Dahl Highland with z² in the log | +6.8 % | +10.3 % | +11.4 % | fails |
| (b) f_dM with the ion's own pv | +3.8 % | +8.7 % | +10.5 % | fails |
| (a) velocity scaling only | −7.0 % | −6.0 % | −5.7 % | fails (as pre-declared) |
| (e) informative: Gottschalk arXiv:1204.4470 eq. fdMi, f_I = (z√A*)^−1.018 | −7.2 % | −10.0 % | −10.3 % | fails |

**Proton baseline**, f_dM vs this reference in water at 100–200 MeV:
- −0.0 to +1.1 % at x/R = 0.1;
- −1.4 to −2.4 % at x/R = 0.9.

Part of the (d) excess at x/R = 0.1 is therefore f_dM's own proton offset. The rest (≈ +2.3 %) is the
first-order linearisation of c_Ω.

**Prediction as frozen:**
- Against this provisional (no-Fano) reference, **no candidate passes 3 % at every I11 point**.
- (d) is closest and fails only at x/R = 0.1 for C/O, by ≈ 0.3 %.
- Under the plan this is STOP 2, unless the reference is completed (Fano, OPEN 1) before I11 is evaluated.
- (a) is expected to fail, as the plan requires; it misses by −5 to −7 %.

## OPEN

1. The Fano atomic-electron term, in the exact form of Gottschalk et al. 1993 (or Fano 1954 / Scott 1963
   §5), is required by decision 0043 §3 and is missing. It is needed both to pass the P1 proton validation
   and before any I11 prediction is final.
2. Pb residual (≤ 1.8 %, growing toward 0.97 R1): ICRU 49 vs the ionmc Bethe range table, to be
   attributed by repeating the Pb evaluation on ICRU 49 / PSTAR ranges.
3. The plan text says "Gottschalk's proton θ_Hanson tables (Be/Al/Cu/Pb)". There is no water block, so
   water is not validated directly.
4. Candidate (c) is a slab proxy, not the differential Kanematsu 2008 ion form (TO VERIFY).
5. The (d) ≈ +3.3 % at x/R = 0.1 comes partly from the c_Ω linearisation. An exact-Ω variant (solve for pv_eq)
   is a possible amendment, to be pre-declared before I11.
6. The proton-validation test needs the cached arXiv source. It runs with `IONMC_CACHE_DIR` set to a cache
   containing object `67fb1478…` and skips otherwise.

## P1b: Fano correction (2026-11-01)

**Source of the formula (quoted).**
- Fano, Phys. Rev. 93 (1954) 117, abstract, case (b): "For incident heavy particles, leave Z² unaltered
  but increase b by Z⁻¹{ln[1130 Z^(−4/3) (c²/v² − 1)⁻¹] − u_in − ½ v²/c²}" (abstract read on the APS
  page; the full paper was not read).
- Gottschalk's θ_Hanson uses this variant: Makarova, Gottschalk, Sauerwein, arXiv:1610.01279, main.tex
  l. 104: "using Molière's rather than Bethe's form of the theory (Z² rather than Z(Z+1)) and using the
  Fano correction for low-Z targets (see Go93)". The same is stated in arXiv:0908.1413 l. 160 ("the
  difference between the Molière/Bethe/Hanson and Molière/Fano/Hanson forms").
- u_in: "−U_in = 3.6 for hydrogen, as calculated exactly by Fano, and −U_in = 5.0 … as suggested by
  Gottschalk for other materials" (MICE, arXiv:2209.10251, apssamp.tex l. 615–617).
- Not found in any cached or open source: an explicit Gottschalk-1993 statement of the equation.
  arXiv:1204.4470 and arXiv:1804.00022 (l. 267: "Fano's approach … computes a correction to b") do not
  give the formula. The cached object `67fb1478…` is arXiv:0908.1413 (ScatPowerV2.tex), not
  arXiv:1204.4470. Sources were downloaded to the scratchpad only; nothing was added to the cache.

**Implemented (`electrons="fano"`, now the default; `electrons="bethe"` keeps the P1 variant).**
- χ_c² = K z² Σ w_i Z_i²/A_i ∫ dt′/(pv)² (Molière Z²).
- ln χ_a² compound weights W_i ∝ w_i Z_i²/A_i. χ_a,i, the 1/(pv)² slab average, 1.167, B − ln B = b,
  the Hanson θ0 and the point nucleus are unchanged.
- b = ln(χ_c²/1.167 χ_a²) + Δb, with D_i = ln[1130 Z_i^(−4/3) β²/(1 − β²)] + (−u_in,i) − β²/2.
- Δb = Σ w_i (Z_i/A_i) D_i / Σ w_i Z_i²/A_i. For an element this is D/Z (Fano). The compound weighting
  (electrons ∝ Z, nucleus ∝ Z²) and the 1/(pv)² slab average of D are **derived here, not sourced**.
- For an ion, D is taken z-independent (Fano's Born result). Higher-order (Bloch-type) electron terms
  are absent, as in Fano (Bondarenco, arXiv:2107.11409 §VI).

**Proton validation, θ_MH/θ_Hanson − 1 in % (same Gottschalk points, ionmc Bethe range tables).**

| x/R1 | 0.001 | 0.01 | 0.1 | 0.2 | 0.5 | 0.9 | 0.97 | max abs (Fano) | max abs (no Fano) |
|---|---|---|---|---|---|---|---|---|---|
| Be | +0.10 | +0.02 | +0.02 | +0.02 | +0.02 | +0.03 | +0.09 | **0.10 pass** | 5.40 |
| Al | −0.02 | −0.00 | −0.00 | −0.00 | −0.02 | −0.02 | +0.00 | **0.02 pass** | 2.26 |
| Cu | −0.11 | −0.09 | −0.09 | −0.10 | −0.11 | −0.17 | −0.21 | **0.21 pass** | 0.95 |
| Pb | +0.51 | +0.48 | +0.48 | +0.51 | +0.61 | +1.10 | +1.32 | **1.32 FAIL** | 1.79 |

- The Be/Al agreement (≤ 0.1 %) is strong evidence that the elemental Fano form is the one Gottschalk used.
- The Pb residual is a flat +0.5 % that grows to +1.3 % at x/R1 0.97. That is the end-of-range signature
  (OPEN 2). Its test stays `xfail(strict=True)`; Be and Al are no longer xfail.

**Frozen (d) correction after Fano.**
- F_scr uses the same corrected χ_a, so the weights become W_i ∝ w_i Z_i²/A_i.
- Δb depends on β only, so it cancels in Ω_ion/Ω_p at equal velocity. The identity Ω_ion/Ω_p = A*/F_scr
  holds unchanged (tested).
- c_Ω and n_Ω are unchanged.

**Candidate predictions vs the Fano-corrected reference** (water; Bethe 78 eV tables).
Grid: 50/100/150/200/300/430 MeV/u × x/R {0.1, 0.25, 0.5, 0.9}. Values are the signed extreme of
θ_cand/θ_MH − 1.

| Candidate | He | C | O | Verdict at 3 % |
|---|---|---|---|---|
| (d) | +5.78 % (200, 0.1) | +6.44 % (200, 0.1) | +6.37 % (300, 0.1) | fails; ≤ 3 % only at C/O 50 MeV/u x/R 0.9 |
| (c) slab proxy | +10.43 % | +16.53 % | +18.69 % | fails |
| (b) | +7.02 % | +13.17 % | +15.77 % | fails |
| (a) | −3.44 % (50, 0.1) | −2.93 % (100, 0.1) | −2.52 % (200, 0.1) | fails only for He at 50 MeV/u (x/R 0.1, 0.9) |
| (e) informative | −3.64 % | −7.27 % | −7.68 % | fails |

- **Proton baseline** (f_dM vs the Fano reference in water, 100–200 MeV): +2.5 to +3.5 % at x/R 0.1 and
  0.5, and +1.3 to +2.1 % at 0.9. Without Fano it was −2.4 to +1.1 %.
- f_dM therefore sits above the Fano water reference by more than Gottschalk's stated ±2 %. Either the
  derived compound rule for water (H with −u_in 3.6) is wrong, or f_dM's water accuracy is weaker than
  stated. Water cannot be checked against Gottschalk's elemental table (OPEN 3, 7).
- **Prediction as frozen:** no candidate passes 3 % at every grid point (STOP 2 under the plan).
  (d) now overshoots by ≈ +6 %. That is the proton baseline (≈ +3 %) plus the collision-number term.
  (a) misses only for He at 50 MeV/u.
- The plan pre-declared that (a) fails. Its near-pass is a finding, not a selection.

## OPEN (after P1b)

1. Fano: implemented (elemental form sourced, Be/Al pass). Not read: the full Fano 1954 paper and
   Gottschalk 1993 eqs. The source of −u_in = 5.0 is secondary (MICE quoting Gottschalk).
2. Pb fails at x/R1 ≥ 0.9 (+1.10 %, +1.32 %). Repeat the evaluation on ICRU 49 ranges (Gottschalk's
   MIXED table, spline in ln R(ln T)) before any further conclusion.
3. Water is not validated directly: there is no compound or water block in `tbl:theta0Single`.
4. Candidate (c) is a slab proxy, not the differential Kanematsu form (TO VERIFY).
5. The exact-Ω variant of (d) is still a possible pre-declared amendment.
6. The proton validation needs the cached `67fb1478…` (arXiv:0908.1413) via the repo cache or
   `IONMC_CACHE_DIR`.
7. The compound Fano rule (Z-weighted D over Z²-weighted nucleus) is derived, not sourced. The water I11
   reference and every candidate verdict depend on it. The f_dM proton offset of +2.5–3.5 % against it
   is a warning sign. A discriminating check is Gottschalk's polystyrene/Lexan/water θ_Hanson (Go93
   Table / LOOKUP), or Molière/Fano values for a compound from an independent source.

## Codex REVIEW-98315d4e fixes (2026-11-01)

Three findings fixed in the module, its tests and this note. No tolerance, candidate, formula of
the reference or plan text was changed.

**Updated file hashes** (supersede the table at the top; `__init__.py` is unchanged):

| File | sha256 |
|---|---|
| `validation/scripts/ion_em/__init__.py` | `7962c4af19580a10f51a4bfc703ce45b67a0f96f010aaaf878aef6919b2ed4ea` |
| `validation/scripts/ion_em/moliere_hanson.py` | `5a99b3b8748a6865584fcbe1337206222ddd7fe197e95f25656eac2324acf164` |
| `tests/ionmc/test_ion_em_moliere_hanson.py` | `b9dab1411e694adea406ceb82b239a0ea6661656c3eae3378b6637b465f56b41` |

### 1. n_Ω from the complete proton Ω

`n_omega` previously returned d ln(R_p/β²)/d ln(pv) and so omitted the β-dependence of the screening
factor and of Fano's Δb, both of which are in the implemented Ω (b = ln Ω + Δb). It now returns
d ln Ω_p/d ln(pv) at fixed x/R for the local form
ln Ω_p = ln R − ln β² − Σ W_i ln(1.13 + 3.76 (Z_i α)²/β²) + Δb(β), analytically:

```
n_Ω = [ E/(S R) + g (Σ_i W_i u_i/(1.13 + u_i) + k [1 + β²/(1−β²) − β²/2] − 1) ] / (d ln pv/d ln E)
u_i = 3.76 (Z_i α)²/β²,  g = d ln β²/d ln E,  k = Σ w_i Z_i/A_i / Σ w_i Z_i²/A_i  (k term only for "fano")
```

The signature gained `elements` and `electrons` (defaults water, "fano"); `theta0_candidate` passes
them. New test `test_n_omega_matches_numerical_derivative_of_full_omega` compares the analytic value
with an independent central difference in ln E (h = 1e-3) of the full ln Ω, at 30–400 MeV and at the
velocity-scaled energies of He 150, C 290, O 350 MeV/u at x/R 0.1/0.25/0.5/0.9, for "fano" and
"bethe": agreement ≤ 1e-3 relative. It also asserts that the omitted terms change n_Ω by more than
10 %.

Before/after (water, Bethe 78 eV proton tables); c_Ω = (c_2 + c_3 L_1)/(ln 10 n_Ω):

| E_p (MeV) | n_Ω before | n_Ω after | c_Ω before (L_1=0) | c_Ω after (L_1=0) | c_Ω before (L_1=lg 0.5) | c_Ω after (L_1=lg 0.5) |
|---|---|---|---|---|---|---|
| 50 | 0.899 | 1.146 | 0.1121 | 0.0879 | 0.1135 | 0.0890 |
| 80 | 0.940 | 1.155 | 0.1072 | 0.0872 | 0.1086 | 0.0883 |
| 100 | 0.963 | 1.168 | 0.1047 | 0.0863 | 0.1060 | 0.0874 |
| 150 | 1.012 | 1.204 | 0.0996 | 0.0837 | 0.1009 | 0.0847 |
| 200 | 1.053 | 1.241 | 0.0957 | 0.0812 | 0.0969 | 0.0822 |
| 250 | 1.089 | 1.275 | 0.0926 | 0.0790 | 0.0937 | 0.0800 |
| 300 | 1.119 | 1.307 | 0.0900 | 0.0771 | 0.0912 | 0.0781 |
| 430 | 1.182 | 1.375 | 0.0853 | 0.0733 | 0.0863 | 0.0742 |

So n_Ω = 1.15–1.38 (was 0.90–1.18) and c_Ω ≈ 0.073–0.089 (was 0.085–0.113); the earlier
"n_Ω 0.90–1.12, c_Ω 0.09–0.11" statements above are superseded. The decision's estimate (0.77 / 0.13)
is further away, so its xfail stays (its reason text was updated).

### 2. Exact P1/I11 grid in the driver

`I11_GRID_MEV_U` (He 50/150/250, C 80/290/430, O 100/350/430 MeV/u: the standard points and ends of
the plan's Notation section; the plan lists no extra common energies) and `I11_FRACTIONS`
(0.1, 0.25, 0.5, 0.9) are module constants; `ion_table`/`main` evaluate exactly this grid (3 energies × 4 fractions = 12
points per ion). The old driver grid
(50/100/150/200/300/430 for every ion, no 0.25) is retired; the tables above in this note that quote
it are superseded by the table below. Test `test_driver_grid_equals_plan_i11_definition` parses the
plan's Notation line and row I11 and fails on any divergence.

Candidate tables regenerated from the committed driver (ratios only; water, Bethe 78 eV tables,
Fano reference): signed range and the max |θ_cand/θ_MH − 1| per ion over the exact grid.

| Cand. | He max abs (where) | C max abs (where) | O max abs (where) | He range | C range | O range |
|---|---|---|---|---|---|---|
| (d) after fix | 4.80 % (250, 0.1) | 5.16 % (290, 0.1) | 5.13 % (350, 0.1) | +1.46…+4.80 | +1.99…+5.16 | +2.20…+5.13 |
| (d) before fix (old n_Ω, exact grid) | 5.77 % (250, 0.1) | 6.41 % (290, 0.1) | 6.33 % (350, 0.1) | +3.09…+5.77 | +3.13…+6.41 | +3.17…+6.33 |
| (c) | 10.43 % (50, 0.9) | 14.87 % (80, 0.9) | 16.05 % (100, 0.9) | +8.38…+10.43 | +10.35…+14.87 | +11.47…+16.05 |
| (b) | 6.78 % (250, 0.1) | 12.20 % (80, 0.1) | 13.98 % (100, 0.1) | +4.44…+6.78 | +8.61…+12.20 | +9.98…+13.98 |
| (a) | −3.44 % (50, 0.1) | −2.95 % (80, 0.1) | −2.45 % (100, 0.1) | −3.44…−0.57 | −2.95…−1.74 | −2.45…−1.40 |
| (e) informative | −3.64 % (50, 0.1) | −7.04 % (80, 0.1) | −7.47 % (100, 0.1) | −3.64…−0.71 | −7.04…−4.69 | −7.47…−5.40 |

Verdict unchanged in kind: no candidate is within 3 % at every point of the exact grid ((a) misses
only for He, −3.44 % at 50 MeV/u, x/R 0.1; (d) still overshoots by up to +5.2 % at x/R 0.1, about
1 percentage point less than before the n_Ω fix; the proton baseline of f_dM against the Fano
reference, +2.5–3.5 % at x/R 0.1, is part of it). The "Candidate predictions vs the Fano-corrected
reference" table and its quoted grid above are superseded by this one.

### 3. Proton-reference gate under IONMC_REQUIRE_DATA=1

`_tex()` now follows the repository convention (cf. `test_nuclear_data_parsers.py`): without the
object (sha256 `67fb1478…`, arXiv:0908.1413 e-print, repository reference cache or
`$IONMC_CACHE_DIR/objects`) it skips as before, but under `IONMC_REQUIRE_DATA=1` it fails with a
message naming the sha256, the searched directories and the materialization step. Note: the Pb test
is `xfail(strict=True)`, so under the flag without the object it reports XFAIL (the other four
proton tests FAIL, so the suite fails regardless).
