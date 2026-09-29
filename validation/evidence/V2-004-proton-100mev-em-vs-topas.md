# Evidence record: 100 MeV protons in water, electromagnetic transport vs TOPAS

- Claim/domain: reference-backend electromagnetic proton transport (energy
  loss, straggling, multiple scattering) in liquid water, 100 MeV, without
  nuclear interactions (PHY-PROTON-WATER partial; NUM step/grid probes).
- Code SHA: b77324d22232e3da7b3576a39b03f4310eb032dd (task V2-004).
- Evidence category: **independent Monte Carlo** (TOPAS 4.3 / Geant4 11.4.2,
  `g4em-standard_opt4` only, 0.1 mm cuts, electrons transported), role
  *evaluation*. Shared lineage: both use I(water) = 78 eV (ICRU 90); Geant4
  uses its own Bethe–Bloch/Bragg models and Wentzel-VI + Coulomb scattering,
  IonMC uses ICRU 90 tables below 16 MeV, its own corrected Bethe layer above,
  Bohr/Gamma straggling and Gottschalk's differential Molière core. No
  parameter of IonMC was tuned to Geant4.
- Reference run: `REF-0088173fbe7787954e15-cac6f372` (20 000 histories,
  12.0 s execution); native outputs `pdd.csv`
  (sha256 3a5f85523dd132a2d0ea9d677b2ef165479f5273f65398e0211aef42f42cece4,
  200 × 0.5 mm bins over the 60 × 60 × 100 mm box) and `lateral.csv`
  (sha256 54babfef48611d44a3c186c419b80f18cec2ebc000ac25efb60807686145b48f,
  120 × 0.5 mm bins in a 1 mm slab at 30–31 mm depth). Case bundle:
  `validation/references/cases/topas-proton-100mev-em-only`.
- IonMC run: reference backend, 4 000 histories, 8 batches, seeds 101/202,
  same geometry and scoring; ≈ 13 s each in the sandbox. Outputs under
  `validation/generated/v2-004/` (ignored) with `comparison-p100-em.json`.

| Observable (units) | TOPAS | IonMC reference | Frozen tolerance | Verdict |
|---|---|---|---|---|
| R80 (mm) | 77.69 | 77.47 | 1.0 mm | pass (Δ −0.22) |
| distal 80–20 % width (mm) | 1.17 | 1.23 | 0.5 mm | pass |
| peak / plateau(10–20 mm) | 6.36 | 6.42 | (gamma 3 %/1 mm proxy) | 1 % |
| depth-energy ratio IonMC/TOPAS, 5–70 mm | — | mean 1.004, range 0.993–1.021 | 3 % local | pass |
| σ_x at 30 mm depth (mm, second moment |x| < 5 mm) | 0.411 | 0.371 | max(5 %, 0.2 mm) | pass on the absolute term (Δ −0.04 mm, −10 %) |
| energy in the slab (MeV/primary) | 0.8935 | 0.8988 | — | +0.6 % |

Interpretation: ranges, straggling width and the depth-dose shape agree
within statistics and the frozen tolerances. The lateral second moment is
10 % narrower than TOPAS at 30 mm: the Gaussian differential-Molière core
carries no single-scattering tail, whereas Geant4's Wentzel-VI model does
(restricting to |x| < 2 mm reduces the TOPAS σ to 0.394 mm, i.e. −6 %). This
is recorded as a known limitation of decision 0041 to be re-evaluated with
full-physics profiles and the frozen PHY-PROTON-WATER criteria (three
depths, three energies) once nuclear interactions and the Warp backends
exist. The ICRU 90 CSDA range (77.59 mm) lies between the two codes.

Failed/contrary evidence retained: none for this case; earlier prototype
runs with an inverse-range table (range +2.2 % at 0.2 mm steps) are
described in decision 0041 and were superseded before this comparison.
