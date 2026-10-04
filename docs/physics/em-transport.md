# Electromagnetic transport of protons

Implemented in `ionmc.physics.em`, `ionmc.physics.kinematics`, `ionmc.physics.scattering` and
`ionmc.transport` for protons from the table floor to 500 MeV with local deposition of delta
electrons and no nuclear interactions. Decision:
[0039](../generated/decisions/0039-proton-transport-engine.md). The step algorithm and the API are in
[transport engine](../architecture/transport.md). Units: energy MeV, length mm, mass thickness and
range g/cm2, density g/cm3, scattering power rad2/mm.

## Mean energy loss

Condensed-history steps use the continuous-slowing-down (CSDA) relation. With `R(E)` the CSDA range
integrated from the same stopping power as the transport table (`dR = E d ln E / S`) and `Rinv` its
exact inverse, the mean energy after a path of mass thickness `t = rho s / 10` [g/cm2] from `E` is

    E1 = Rinv( R(E) - t )        mean loss = E - E1

which is independent of how the path is divided into steps *if the tables were exact*. They are not: every step
reads `R(E)` of the new state energy from the table, and the table round trip `R(Rinv(R)) - R` (up to the U2 bound of
1e-5, a smooth and therefore systematic error) enters once per step, so the accumulated loss depends weakly on the
number of steps (about 3e-4 of the loss between 0.1 mm and 0.01 mm steps with the telescoping form alone). For
`t < f_short R(E)` the loss is therefore the linear form at the midpoint energy, `S(E_mid) t` with
`E_mid = E - S(E) t / 2` (midpoint rule: the step-size error is second order, the loss of the same path agrees
between 1, 0.1 and 0.01 mm steps to 5e-6 and 1e-7 when every step takes this branch; the former form `S(E) t` had a
first-order bias of about `s / (2 R)` times the curvature of `S`, which showed as a step dependence of the
depth-dose at the 1e-3 level in the T9 comparison of 0.1 mm against 1 mm steps). The default is `f_short = 1e-2`, the largest
value the configuration accepts: for such steps the second-order midpoint error is smaller than the table round-trip
bias of the telescoping form at small steps, whereas where `s / R` is large (the last steps before the end of the
range) the midpoint rule is no longer accurate (its error is second order in `s / R`) and the telescoping form takes
over; the float32 cancellation that the branch once avoided is gone because the energy bookkeeping is in float64.
For `t >= R(E)` the loss is `E`. `ionmc.transport.tables.TransportTables` stores `ln S`, `ln R` on a uniform `ln E` grid
(at least 200 points per decade) and `ln E` on a uniform `ln R` grid (800 points per decade);
values are log-log linear interpolations. The test `test_u2_round_trip_and_monotonicity` requires
`|Rinv(R(E)) / E - 1| <= 1e-5` over the table range and strictly monotone `R` and `Rinv` for water and
copper, and that `R` follows the source table within 1e-4. The transport state energy itself is
carried along the track (not `Rinv(R(E))`), so the energy variable is never re-derived from the
range table; the per-step interpolation round-trip error of `E - Rinv(R(E) - t)` (≤ 1e-5 relative
per evaluation) does, however, accumulate with the number of telescoping steps — the measured
3e-4 difference between 0.01 mm and 0.1 mm steps above — which is why steps below 1e-2·R use the
midpoint branch instead.

The step length is limited by the voxel plane, by `f E / S_lin` with `f = 0.02`, by the Geant4 range
function `alpha R + rho_f (1 - alpha)(2 - rho_f / R)` for `R > rho_f` (otherwise `R`; defaults
`alpha = 0.2`, `rho_f = 0.1 mm`) and by the maximum step (default 1 mm). Energy below
`E_cut = 2 MeV` is deposited locally.

## Straggling

The Bohr variance of the energy loss over a step of length `s` at mean energy `E_m = E - mean/2`
is

    sigma^2 = (K/2) (Z/A) rho (s/10) z^2 T_max (1/beta^2 - 1/2)        [MeV^2]

with `K = 0.307075 MeV cm2/mol`, `T_max = 2 m_e c^2 beta^2 gamma^2 / (1 + 2 gamma m_e/M + (m_e/M)^2)`.
This equals the Geant4 dispersion `2 pi r_e^2 m_e c^2 n_el z^2 x T_max (1/beta^2 - 1/2)`. The loss is
sampled with the model `bohr_gamma_v1` (the default) or the older `bohr_gauss_clamped_gamma_v1` (selectable with
`PhysicsOptions.straggling_model`). The older model (ratio `r = mean / sigma`):

* `r < 3`: Gamma distribution with shape `k = r^2` and scale `sigma^2 / mean` (Marsaglia-Tsang; for
  `k < 1` the shape-plus-one sample times `u^(1/k)`), which reproduces the mean and the Bohr variance
  exactly;
* `r >= 3`: Gaussian, clamped (no resampling) to `[0, 2 mean]`. This preserves the mean and reduces
  the variance by at most about 0.5 % (0.49 % at `r = 3` from the Gaussian tails; the moment test
  allows 1 %); 0.13 % of the samples are clamped at 0 at `r = 3`.

A Gamma sample is accepted by the Marsaglia-Tsang test; a rejection draws another Philox block, at most 64
attempts (then the `straggling_rejection` counter invalidates the run). Step energies of 1 mm or more
at therapeutic energies have `r` well above 3, so the Gamma branch is used near the end of the range.
`test_straggling_gamma_moments_exact_and_gaussian_variance_loss_bounded` checks the Gamma branch for both moments and the Gaussian branch for the mean and the bounded variance deficit on a grid of `(mean, sigma)` that covers both
branches and the boundary (4e5 draws per point, three standard errors; the Gaussian branch near the boundary
is given the 1 % clamp allowance). The transport-level validation of straggling is the range
straggling `sigma_R` of the end depths (criterion T6, local validation), which is not part of the CI tests.

Default model `bohr_gamma_v1`: a Gamma distribution whose raw sampler has exactly the Bohr mean and variance for every
ratio (shape `k = r^2`, scale `sigma^2 / mean`, Marsaglia-Tsang), positive everywhere, no clamp and no Gaussian branch.
With a common scale `theta = sigma^2 / mean = kappa(E) / S(E)` along a path (`kappa` the Bohr dispersion per unit
path, `S` the stopping power) the sum of Gamma steps is exactly Gamma with the summed shape, so the whole energy-loss
distribution - not only its first two moments - is independent of the step length wherever `theta` varies slowly; the
older model switches branch with the step length (Gamma for `r < 3`, Gaussian above), so its higher moments depend on
the step. The full-scale T9 comparison of the two models is recorded in decision 0039.

Both models: the exact-moment statements above are those of the raw sampler. In transport every sampled loss is capped
at the remaining kinetic energy (`loss = min(loss, E)`). The exceedance probability `P(loss > E)` grows as the ratio
`r` decreases and as `mean / E` increases: for example `r = 0.2` at `mean = 0.02 E` gives about 2e-3, and the
demonstrated range of the moment tests is `r >= 0.7`. Which pairs `(r, mean / E)` occur in transport is fixed by the
physics: `r` scales as the square root of the step length at fixed energy, so small `r` occurs only for very short
steps, whose mean loss is a tiny fraction of `E`, and the energy-loss step limit keeps `mean / E <= 0.02` except in the
last step before the cutoff. In water (the offline Bethe tables; test
`test_transport_domain_pairs_keep_the_energy_cap_negligible`):

| E [MeV] | step 1 mm: r, mean/E | 0.1 mm | 0.01 mm | 0.001 mm |
|---|---|---|---|---|
| 10 | 48, 0.45 (beyond the limit: a last step) | 15, 0.045 | 4.8, 4.5e-3 | 1.5, 4.5e-4 |
| 50 | 13, 0.025 | 4.1, 2.5e-3 | 1.3, 2.5e-4 | 0.41, 2.5e-5 |
| 150 | 5.4, 3.6e-3 | 1.7, 3.6e-4 | 0.54, 3.6e-5 | 0.17, 3.6e-6 |

With `mean / E <= 0.02` the fraction of draws above `E` is below 3e-7 for every ratio of the table (so the cap acts only
in the last step before the cutoff, where the mean loss is a large part of `E` and the cap reduces the mean and the
variance).

## Multiple Coulomb scattering

The differential Moliere scattering power of B. Gottschalk, Med. Phys. 37 (2010) 352
(arXiv:0908.1413) is used, for the projected angle:

    T_dM = f_dM(pv, p1v1) (z E_s / pv)^2 / X_S           [rad2 per unit length], E_s = 15.0 MeV
    f_dM = 0.5244 + 0.1975 L1 + 0.2320 L2 - 0.0098 L1 L2,  L1 = lg(1 - (pv / p1v1)^2), L2 = lg(pv)

with `pv` in MeV, `p1v1` the `pv` of the particle at birth (the source energy) and `lg` the base-10
logarithm; `f_dM` is clamped at 0 (the fit diverges as `pv` approaches `p1v1`). The scattering length
`X_S` [g/cm2] follows from

    1 / X_S = alpha N_A r_e^2 (Z^2 / A) { 2 ln( 33219 (A Z)^(-1/3) ) - 1 }

per element and Bragg additivity by weight for compounds (natural logarithm). `test_u3_...` checks
water, Be, Al, Cu and Pb against the tabulated 46.88, 92.60, 28.75, 14.62 and 6.62 g/cm2 within
0.3 %. `T_dM` per mm of a medium of density `rho` is `rho / (10 X_S)` times the factor above.

Per step the variance is the midpoint rule `s T_dM(E_mid)` with `E_mid = Rinv(R - rho s / 20)`;
`p1v1` is the particle's own birth value. On the first step of a particle's life (birth step)
`1 - (pv/p1v1)^2` grows linearly from 0 and `f_dM` is logarithmically singular, so the logarithmic
term is replaced by its analytic (linearized) step average `L1 = lg(1 - (pv(E1)/p1v1)^2) - 1/ln(10)`, exact for a linear growth of `1 - (pv/p1v1)^2` along the step and accurate to about 1e-3 under CSDA slowing for steps up to 1 mm,, with `E1 = Rinv(R - rho s / 10)`
for the planned step (the shared function `scattering_variance_birth`); the smooth terms and the
prefactor are taken at `E_mid` and `f_dM` is clamped at 0 after averaging. The variance defines a two-dimensional
Gaussian polar angle `theta = sqrt(-2 var ln u)` (limited to pi) and an azimuth `2 pi u`. The
direction is rotated with the `rotateUz` construction of Geant4 and renormalised.

### Random hinge

The step of length `s` is split at a uniformly distributed fraction: the particle moves `a s` along
the old direction, is deflected, and moves `(1 - a) s` along the new direction (cut at a voxel plane
if necessary, in which case the particle snaps onto the plane). This reproduces the Fermi-Eyges
second moments of a thin layer without an extra draw; energy is deposited along both legs in
proportion to the path length in each scoring voxel (track-length apportioning). The angle is sampled for the planned `s`, so a second leg cut by a plane carries a small
overestimate of scattering on boundary steps.

### Validation status of this model (tests in `tests/ionmc/test_tables_scattering.py`)

Evidence class of the checks below (U3, U4, U4b and the `X_S` comparison): **source-model reproduction
(Gottschalk 2010 formulae and tables; not independent)**. U3 reproduces the paper's `X_S` table from its
own formula, U4 reproduces its `T_dM` column, and U4b is implied by U4 together with the paper's dM %
column (both read from the same table). The generalised-Highland comparison is a related model and the
step-independence check U5 is deterministic self-consistency. The only independent qualification of the
scattering model is T15 (independent Monte Carlo comparison), which is pending.

* `X_S` against Gottschalk's table: within 0.3 % (above; source-model reproduction).
* `theta_dM(x)` by quadrature of `T_dM` along the CSDA path at 158.6 MeV, `x/R1` = 0.01, 0.1, 0.5, 0.9,
  against `theta_Hanson (1 + dM %/100)` of the table of Gottschalk (2010), parsed at run time from the
  cached LaTeX source of arXiv:0908.1413 (not committed; the test is skipped when the cache is absent,
  so it runs locally (LV) and not in CI).
  Frozen criterion 1.5 %; maximum deviations (test `test_u4_theta_dm_against_gottschalk_table`):
  Be 0.31 %, Al 0.29 %, Cu 0.21 %, Pb 1.12 %. The table has no water block. Residuals come from the
  range tables (our Bethe `rho R1` differs from the paper's by 0.85 % for Pb and by less than 0.2 % for Be, Al and Cu).
* U4b (`test_u4b_theta_dm_against_theta_hanson`): our `theta_dM` against the paper's `theta_Hanson`
  column itself, which does not depend on the paper's dM fit; frozen tolerance 4.5 % (the paper's own
  `T_dM` is within 2.74 % of it over the frozen points and U4 bounds our reproduction at 1.5 %).
  Maximum deviations: Be 2.31 %, Al 2.26 %, Cu 1.52 %, Pb 1.41 % (the paper's own maximum `|dM %|`
  over the frozen points: Be 2.21, Al 2.06, Cu 1.41, Pb 2.51 %).
* `X_S` against the table of the same paper (`test_scattering_length_against_gottschalk_table`):
  within 0.5 % for water, Be, Al, Cu and Pb (largest deviation 0.08 %, Pb).
* Cache-independent cross-check: the quadrature against the generalised Highland formula of Gottschalk
  et al. (1993) with the radiation length from Tsai's formula agrees within 6 % for water (a related
  model with its own accuracy of a few per cent; not an acceptance criterion).
* Step independence of the integrator (test U5, frozen criterion 2e-3 for every step size at
  `x/R1 >= 0.05`, water, 150 MeV): relative deviation of the stepped sum from the quadrature:

  | x/R1 \\ s [mm] | 0.01 | 0.1 | 0.5 | 1 | 2 | 5 |
  |---|---|---|---|---|---|---|
  | 0.05 | +3.7e-6 | +4.4e-5 | +2.1e-4 | +3.7e-4 | +5.5e-4 | -3.3e-4 |
  | 0.25 | +5.8e-7 | +6.8e-6 | +3.2e-5 | +5.9e-5 | +9.5e-5 | +7.4e-5 |
  | 0.5 | +1.8e-7 | +2.6e-6 | +1.1e-5 | +1.7e-5 | +1.4e-5 | -1.1e-4 |
  | 0.9 | +7.2e-8 | -2.6e-8 | -1.6e-5 | -6.9e-5 | -2.8e-4 | -1.5e-3 |

  The midpoint rule alone is biased by the logarithmic singularity of the first step (several per
  cent for large steps); the analytic birth-step average removes it (residual of order 1e-3 from the linearization). `test_birth_variance_matches_fine_quadrature`
  checks the birth variance against quadrature (1e-3 for 0.1 to 1 mm; not tested at 0.01 mm, where
  the step changes pv by less than the 1e-5 round-trip accuracy of the tables and the quadrature's
  `1 - (pv/p1v1)^2` no longer starts from 0). `E1` is evaluated for the planned step because the
  angle is sampled before the travelled length is known. A per-step Highland formula varies
  by more than 5 % over the same step sizes (negative control).

All numbers of this section are produced by `validation/scripts/transport/mcs_checks.py` (same
functions as the tests, `ionmc.transport.mcs_checks`) and archived in
`validation/results/transport/mcs_checks.json` (aggregates only, with the hash of the consumed
source); see `validation/results/transport/README.md`. The negative control gives a relative
spread of 1.6 (160 %) of the per-step Highland sum over the step lengths, against 2e-5 for the
engine's stepped `T_dM` up to 1 mm; the worst U5 deviation is 1.5e-3.

## Sources and units

Gottschalk, "On the scattering power of radiation therapy protons", Med. Phys. 37 (2010) 352;
Gottschalk, Hanson et al., "Multiple Coulomb scattering of 160 MeV protons", Nucl. Instrum. Methods
B74 (1993) 467 (generalised Highland); Geant4 `G4IonFluctuations`, the `G4ThreeVector::rotateUz` rotation
convention and the range-limited step function; Marsaglia and Tsang,
ACM Trans. Math. Softw. 26 (2000) 363 (Gamma sampling); Particle Data Group, "Passage of particles
through matter" (Bohr variance, radiation length).
