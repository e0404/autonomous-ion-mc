# V3-008 pre-freeze check P11: exact-Omega mapping (record, 2026-11-01, base 3d2547a)

Nothing is frozen by this note. The I11-rel / I11-abs numbers below are PRE-FREEZE VALUES SEEN, report only.
The gate (I11-rel <= 2 %) is evaluated once under freeze; it is not asserted in any unit test.

## Files (bind in P9; sha256)

| File | sha256 |
|---|---|
| `validation/scripts/ion_em/exact_omega.py` | `92a2feab3fb454220ad3a36736e2ed17caa9ac57ca9648d5fba4b594ed48a58c` |
| `tests/ionmc/test_ion_em_exact_omega.py` | `d4f08dc34b11e3f7f3548449501d9c9bb1b8d674c2ad5b894f73d697b9085f2e` |
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

**I11-map residual |b_p(pv_eq) - b_ion|:** max 1.8e-14 on the grid (tests: 4.8e-14 on grid plus 7 extra
depths per energy); limit 1e-10. z = 1, M = m_p: pv_eq = pv bit-exact, f_dM equal (<= 1e-12 test, 20 points).
Omega identity (scaled synthetic path, z = 2, 6, 8, Fano and Bethe): <= 1e-10.

**I11-rel max |ratio - 1|** (pre-freeze, seen; gate 2 %):

| Ion | max | where (MeV/u, x/R) | min |
|---|---|---|---|
| He | 1.73 % | 50, 0.1 | 0.70 % |
| C | 1.83 % | 80, 0.9 | 0.84 % |
| O | 1.76 % | 100, 0.9 | 0.87 % |

All values are positive (model ion/proton ratio above the MH ratio); margin to the 2 % gate is 0.17 pp (C).

**I11-abs max |theta_model/theta_MH - 1|** (informative): He 4.65 % (250, 0.1); C 5.31 % (430, 0.1);
O 5.32 % (430, 0.1); range +1.85 to +5.32 %.

pv_eq/pv on the grid: He 0.57-0.69, C 0.18-0.36, O 0.13-0.27. Proton p1v1_eq reaches 0.8 GeV (He),
2.5 GeV (C, O) kinetic energy: the proton table must extend to about 3 GeV (this check used 8 GeV).

## OPEN

1. Transport implementation needs proton cumulative b-tables out to ~3 GeV (above the 500 MeV/u
   project table); to be decided before P9.
2. I11-rel C margin is small (0.17 pp); STOP rule applies if freeze evaluation exceeds 2 %.
3. Warp twin of the root (I11-map "python twin vs Warp") is not implemented here.
4. Inner slab nodes 64 vs reference 400: b differs by < 1e-7.
