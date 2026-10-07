# 0040 — LET definitions and extensible scoring: species, fluence, lookup tallies and channel accumulators

- Status: accepted (plan ratified 2026-10-03, plan delta ratified 2026-10-07)
- Date: 2026-10-07
- Task: V3-004 (builds on V3-003A/B; V3-005A and V3-009 extend it)
- Affects: physical accuracy, numerical accuracy, reproducibility, architecture, validation strategy, performance

## Problem

The transport engine of decision 0039 scores one quantity, the energy deposit. The requirements
ask for dose, energy deposition, track- and dose-averaged LET, fluence, species- and
energy-resolved scoring and external (biological) lookup data on the same scoring grids
(V1-MUST-012 to 016, 025, 029, 030, 037, 040). Decision 0038 defines the electronic stopping
power data roles; it does **not** define LET (a correction of an earlier assumption). The only
prior statement is the research note `docs/research/em-physics.md` sections 5 and 6, which is a
recommendation and not a decision. This decision fixes the definitions, the contributors, the
assignment rule, the accumulator and quantization policy, the ratio reduction, the species
registry and the lookup provenance policy, before any LET result exists.

## Context

- Literature on LET scoring in Monte Carlo (abstracts read and verified; no full-text claim is
  made beyond what the abstracts state):
  - Cortes-Giraldo and Carabe, Phys. Med. Biol. 60 (2015) 2645 (PubMed 25768028): the
    estimated proton LET depends strongly on the scoring method and on the voxel and step size;
    their lookup-table-S method was the most stable.
  - Guan et al., Med. Phys. 42 (2015) 6234 (PubMed 26520716): the step-limit effect is large for
    the dose-averaged LET and small for the track-averaged LET.
  - Granville and Sawakuchi, Phys. Med. Biol. 60 (2015) N283 (doi 10.1088/0031-9155/60/14/N283): the dose-averaged LET varies more between scoring techniques than the track-averaged
    LET; no single technique is recommended.
  - Grassberger and Paganetti, Phys. Med. Biol. 56 (2011) 6677 (PubMed 21965268): distal maximum of
    the dose-averaged LET of about 12 keV/um for primary protons in their beam; used only as an
    exploratory, non-gating comparison (A13).
  - Kalholm et al., Radiother. Oncol. 161 (2021) 211, doi:10.1016/j.radonc.2021.04.007: an
    unstated reference medium makes LET ill-defined, so the medium has to be part of the
    definition.
- Facts of V3-003B used here: the deposit of every step is apportioned along both hinge legs,
  piece by piece through a per-grid DDA walk, with a linear stopping-power ramp (decision 0039);
  a piece is quantized to `floor(x 2^30 + 1/2)` and added to an int64 grid; the deposit
  bookkeeping is float64 in every kernel variant; the cutoff energy is a point deposit.
- Unrestricted electronic stopping power matches the local deposition of delta electrons of
  decision 0039 (`docs/research/em-physics.md` section 5), so the LET is the unrestricted LET.

## Decision

### 1. LET definitions

Per transport step with path length l (mm, the actual step `s_act`), deposit eps (MeV), species
and generation:

- `E_mid = E0 - dE_m / 2`, with `dE_m` the CSDA mean loss that the step already computes. This is
  the path-midpoint energy to O(l^2); the sampled (straggled) eps never enters S.
- `S = S_w(E_mid)`: the unrestricted electronic stopping power **in water** of the transported
  ion, read from a water row added to `TransportTables` by the same `StoppingSource` as the
  medium tables (the source identity is recorded). MeV/mm is numerically keV/um.
- Track-averaged `LET_t = sum(l S) / sum(l)`; dose-averaged `LET_d = sum(l S^2) / sum(l S)`.
  Diagnostic `LET_d^eps = sum(eps S) / sum(eps)`, which has the same expectation as `LET_d` in
  water because `E[eps] = S l` in expectation (the Bohr/Gamma sampler reproduces the mean).
- Fluence `Phi = sum(l) / V_geom` (mm^-2 per primary); fluence spectra are `l`-weighted
  histograms of the per-nucleon energy with explicit underflow and overflow bins.
- **Why not eps/l.** The estimator `sum(eps (eps/l)) / sum(eps)` has a per-step expectation
  `S + sigma^2/(S l^2)`. With the Bohr variance `sigma^2 = kappa l` and `kappa_water` about
  0.0087 MeV^2/mm (derived in the plan from the Bohr formula) the bias is `kappa/(S l)`: about 30 %
  at l = 0.1 mm and 3 % at l = 1 mm for 150 MeV protons. This is consistent with the strong
  voxel-size dependence reported by Cortes-Giraldo and Carabe and with the larger step-limit
  effect of the dose-averaged LET in Guan et al.; the table-S estimators here are not affected.
  The negative control of acceptance row A7 reproduces the bias.

### 2. Step quantities and the midpoint-anchored ramp (assignment rule)

A step is apportioned over scoring voxels in pieces (decision 0039). Each piece p of length
`l_p` at path coordinates `[t_a, t_a + l_p]`, `tau_p = t_a + l_p/2 - s_act/2`, receives:

- From one water lookup per step at `E_mid` (bin i, fraction f): `S_mid = S_w(E_mid)` and the
  log-log slope `gamma = (ln S_{i+1} - ln S_i)/dln E` (the exact derivative of the interpolant);
  the ramp `k = dS/dt = -gamma S_mid dE_m / (E_mid s_act)` (MeV/mm^2) and the energy rate
  `Edot = dE_m / s_act`.
- `Sbar_p = S_mid + k tau_p`, `Ebar_p = E_mid - Edot tau_p` (per nucleon `Ebar_p/A`),
  `m1_p = l_p Sbar_p` and `m2_p = l_p (Sbar_p^2 + k^2 l_p^2/12)`: the exact integrals of S and S^2
  for a linear ramp.
- `eps_p`: the float64 deposit piece that V3-003B already computes (`deposit * ramp_weight`),
  passed unchanged to the edep grid and to every channel.

| Kind | Increment per piece |
|---|---|
| E | `eps_p` |
| L | `l_p` |
| LS | `m1_p` |
| LS2 | `m2_p` |
| ES | `eps_p Sbar_p` |
| FE | `eps_p f(x_p)` |
| FL | `l_p f(x_p)` |
| N | 1 (piece count, diagnostic) |

The lookup or spectrum argument `x_p` is `Ebar_p/A` for energy axes and `Sbar_p` for LET axes,
evaluated per piece. Routing:

- A step with `s_act > 0` is LET-contributing (class "step"). Its legs are walked even when
  `eps = 0`, so LET and fluence do not depend on `eps > 0`.
- Cutoff point deposits and `eps > 0` at `s_act = 0` are class "local" with `l = 0`; they reach
  only E-kind channels whose class mask includes "local" (the automatic channel
  `edep_excluded_from_let`); LET, ES, FE and FL never see them.
- `LET_d^eps = ES / E_step`, where `E_step` is the E channel of class "step" with the same species
  and generation filter, not the total edep.

**Why a midpoint-anchored ramp** rather than a constant `S_w(E_mid)` per step or a lookup per
piece: step sums are `sum_p m1_p = s S_mid` exactly (because `sum l_p tau_p = 0`), so the ratified
midpoint definition holds per step and the deterministic closed form (A4) keeps its O(l^3)
midpoint accuracy. With a constant S per step each piece is wrong to first order and the errors
cancel per voxel only for symmetric crossings. A numpy check (CSDA with S proportional to
E^-0.79, 1.5 mm voxels against 1 mm and 0.4 mm steps) gave worst per-voxel errors of `LS` against
`E(z1) - E(z2)` of 2.6e-4 (150 MeV, 1 mm) and 2.2e-4 (60 MeV, 0.4 mm) with a constant S, and 5e-8
and 7e-7 with the ramp. The cost is one water lookup per step plus three FMAs per piece.

**Water-only caveat of LET_d^eps.** `E[eps_p S_bar_p] = S_bar_p E[eps_p]` because `S_bar_p` is
deterministic given the step start. `eps_p` uses the straggled `S_m(E1)` in the ramp weight,
which shifts `E[eps_p]` by a term antisymmetric along the step (zero sum over the step) of
relative size at most `gkappa/(2 S^2)`: 0.85e-4 at 150 MeV, 1.5e-4 at 10 MeV, 1.9e-4 at 3 MeV
(`g = -dS/dE`). The identity `E[ES] = E[LS2]` therefore holds in water, or wherever `S_m/S_w` is
energy-independent over the voxel spectrum; in another medium `LET_d^eps` is a dose-to-medium
weighted mean of `S_w`. Acceptance row A7 compares the two estimators in water only.

### 3. Contributors and exclusions

- Contributors to LET and to the lookup averages: every **transported charged** particle
  (primaries now, secondaries once V3-005A transports them). Species-resolved partial sums
  (numerator and denominator per species) allow recombination.
- Excluded from LET and lookup averages: neutrals, particles below the cutoff energy and
  nuclear-local deposits (V3-005). Their energy is reported in the automatic per-voxel channel
  `edep_excluded_from_let`. A comparison with TOPAS or any other engine must apply the same
  convention; the distal LET_d is lowered by the cutoff exclusion (E_cut = 2 MeV corresponds to a
  residual range of about 0.07 mm in water).
- Dose stays dose-to-medium. A request for dose-to-water fails closed (a dedicated channel kind
  would be needed; there is no silent conversion).
- The LET medium is water. `let_medium="local"` (stopping power in the voxel's own material)
  fails closed; it is a later channel kind, not a silent default. The water table is part of the
  recorded table identity.

### 4. Accumulators and quantization (all channels)

All channels live in one int64 array `acc[B, sum_c size_c]` separate from the qualified `edep`
array, which stays untouched (bit-identical to a524f209 for `tallies = ()`; row A16).

- Increment `n = floor(x 2^k_c + 1/2)`, computed in float64 in every variant (exact power-of-two
  scale) and added with an integer atomic add. Integer sums are associative, so bit identity
  across chunk sizes, chunks and workers holds for every channel within a device kind.
- Lookup values must be non-negative, so every channel is non-negative; this keeps the existing
  wrap check (value at or above 2^62 or negative gives `accumulator_overflow`) valid.
- `k_c` is a deterministic function of the configuration and tables only; it is recorded in
  `EffectiveConfig` and in the `QuantityResult` provenance.

Per-history bound `B_c` and quantum `q_c = 2^-k_c` (`hpb = n_histories / n_batches`,
`S_max = max S_w` over `[E_cut, E_hi]`, `S_ref = S_w(E_hi)`,
`B_L = 1.25 int_0^{E_hi} dE / min_m S_lin,m(E)` in mm (25 % straggling margin; `S_lin,m = rho_min,m S_mass,m / 10` over the materials present; a runtime-checked bound, see the path bound below), `B_LS = S_bar_max B_L`, `B_LS2 = S_bar_max^2 B_L`, `B_ES = S_bar_max E_hi`,
`r_max = max S_w,lin/(rho_min S_m,lin)` over grid materials and energies, `f_max` the largest
table value):

| Kind (unit) | Quantum | B_c per history | Capacity |
|---|---|---|---|
| E (MeV), including species partials and `edep_excluded_from_let` | fixed `2^-30`, identical to edep | `E_hi` | guaranteed by the existing check; the same piece with the same rounding gives an exact integer identity |
| FE (MeV [f]) | `2^-30 2^ceil(log2 f_max)` | `f_max E_hi` | follows from the E check; with `f` a power of two `FE = f E` bitwise |
| L (mm), FL (mm [f]; spectra f = 1) | adaptive `k = min(40, floor(62 - log2(hpb B_c)))` | `B_L`; `f_max B_L` | fits by construction; spectra share L's `k`, so the bins sum to L bitwise |
| LS (MeV) | adaptive | `S_max B_L` | fits by construction (given the checked path bound) |
| LS2 (MeV^2/mm) | adaptive | `S_max^2 B_L` | fits by construction (given the checked path bound) |
| ES (MeV^2/mm) | adaptive | `S_max E_hi` | fits by construction |
| N (count) | `k = 0` | `max_steps 2 scoring_pieces` | exact |

Two automatic N channels exist per used grid (review amendment, V3-004): `N` counts the pieces of
class "step" and `N_local` the pieces of class "local" (point deposits, all species including the
pseudo-species; the same bound is used). The rounding bound of a channel uses the counts of the
classes in its class mask (edep and dose: both; the excluded-energy channel: local; all others:
step), so that the reported `rounding_bound` is a deterministic bound.

**Ramp envelope (fail closed, review amendment).** The hook scores `S_bar = S_mid + k tau` with
`S_mid = S_w(E_mid)`, `E_mid = E - dE/2`, `k = -gamma S_mid dE / (E_mid s_act)` and
`|tau| <= s_act/2`, hence `S_bar in S_mid (1 -+ |gamma| f_E / 2)`, `f_E = dE/E_mid`. The CSDA mean
loss is `E - E_1` with `E_1 >= 0`, so `dE <= E` and `E_mid >= E/2`: a step that crosses the cutoff
(`E > E_cut`) reaches `E_mid > E_cut/2` (the configuration already requires `E_cut >= 2 x` the
table floor) and `f_E <= 2`, so the amplification is at most `|gamma|`. `compile_channels`
evaluates, per bin of the runtime water row over `[E_cut/2, E_hi]`, the end values of `S_w` and the
bin slope `gamma`, and takes `S_bar_max = max S_w (1 + |gamma|)` and `S_bar_min = min S_w
(1 - |gamma|)`. `S_bar_max` replaces `S_max` in the LS, LS2 and ES bounds; `S_bar_min/max` (with the
envelope started at the table floor) give the coverage of a LET lookup axis. If any bin has
`|gamma| >= 1` (the ramp could reach zero or become negative) the configuration is rejected in
`validate()`, never clamped. A runtime guard counts any piece with `S_bar <= 0` in the
out-of-domain counter, which invalidates the result.

**Path bound B_L: heterogeneous envelope, enforced at runtime (reviews 0352ebfd, 36271bf4).**
The homogeneous range `max_m R_m(E_hi)/rho_min,m` is not a bound for a path that changes material
as the energy falls: a path that follows the lowest linear stopping power at each energy has the
length `int dE / min_m S_lin,m(E)`, which can exceed `max_m int dE / S_lin,m` by an arbitrary factor
for accepted custom tables. The runtime rows are log-log piecewise linear on the union of the
material nodes; on one interval every `ln S_m` is linear in `ln E`, so `ln min_m S_m` is concave and
lies above its chord through the node values. Replacing the minimum by that power-law chord can only
raise the integral, which is then exact per interval (`E_a/S_a (r^(1-g) - 1)/(1-g)`, `ln r` for
`g = 1`; `mixed_path_bound_mm`). Below the lowest node the rows are clamped (constant) as at
runtime; the density enters per material through the lowest density of its voxels. Water only:
`B_L = 198.3 mm`.

That is a CSDA bound. It is not a hard bound with straggling: neither sampler guarantees that a
history stops within any multiple of its CSDA path (the legacy sampler can clamp losses to zero,
Gamma losses can be very small). The rigorous alternative `max_steps max_step_mm` (a step is at most
`max_step_mm` long, a history at most `max_steps` steps; recorded as `B_L_truncation_mm`) is 3896 mm
for the water example, 20 times the CSDA path (the default `max_steps` is 1948 against about 80
steps actually made): a design on it alone would fail closed on the precision floor above about
2e6 histories per batch at 150 MeV (LS2 quantum 2^-17 at 4e6 against the floor), so it is not used
for the quanta.

**Checked assumption.** The engine accumulates each history's scored path (the sum of `s_act`, the
quantity the L channel scores, float64) in the reference and in the kernel and, when it exceeds
`B_L`, sets the per-history tally column `path_bound_exceeded` (after the residual columns and the
lookup out-of-domain count; the nine-counter layout is unchanged, the column exists only with
channels, so the qualified path and A16 are bit-identical). The reduction surfaces the count as
`ChannelRaw.path_bound_exceeded`; a nonzero count invalidates the result exactly as the lookup
out-of-domain count does. Hence either every history satisfied `L_h <= B_L`, in which case
`L, FL <= B_L`, `LS = sum l S_bar <= S_bar_max B_L` and `LS2 <= S_bar_max^2 B_L` hold for every
history and the int64 capacity proof holds exactly as compiled, or the result is invalid and no
accumulator value is ever used, so a wrap-around cannot produce a silently wrong result. `ES <=
S_bar_max E_hi` needs no assumption (`sum eps <= E_hi` by energy conservation); FE is bounded by
`f_max E_hi` likewise. The stopping-ratio candidate `1.25 r_max E_hi` of the previous LS bound is no
longer used (it is energy-derived and not covered by the checked predicate); `r_max` stays recorded.
Regressions: `test_path_bound_covers_a_path_that_changes_material_with_energy` (PMMA at 1 MeV/mm above
6 MeV, aluminium at 1 MeV/mm below, 12 MeV protons: legacy bound 8.10 mm, B_L 14.96 mm, observed
9.59 mm), `test_low_loss_sampler_violates_the_path_bound_and_invalidates` (patched near-zero-loss
sampler: flagged and invalid; ordinary sampler: 0 of 1000 histories) and the warp-cpu flag test.

**Stopping-ratio `r_max` (historical: formerly the second LS candidate, review amendments).**
SUPERSEDED for the capacity bound by the checked-path-bound design above: `B_LS = S_bar_max B_L` and
`B_LS2 = S_bar_max^2 B_L` only; the energy-derived candidate `1.25 r_max E_hi` is NOT used. `r_max`
is still computed and recorded (informative) over the complete reachable midpoint domain: `r_max = max (1 + |gamma_w(E)|) S_w(E) / (rho_min S_m(E))` over `[E_cut/2, E_hi]` for every
material row (union of material and water-row nodes plus both ends; both stopping powers are
log-log piecewise linear, so the ratio is monotone between nodes), with the ramp amplification of
the runtime water-row bin applied to `S_w` as in the envelope above. Every quantity of the proof
(`S_w`, `gamma_w`, `S_m`, `S_ref`, and the envelope and lookup-coverage extrema) is taken from the
*runtime rows* the scorer evaluates (`water_ln_s_mass` and the material rows), never from the
source-table interpolant: after resampling onto the uniform runtime grid the two agree only at
runtime nodes and differ between them (a kink inside a bin), and the candidate set is the union
of the runtime water and material nodes plus the domain ends. The previous domain
`[E_cut, E_hi]` was wrong because the scorer evaluates `S_w` down to `E_cut/2` and a custom material
table may drop sharply below the cutoff. In water `r_max = 1 + |gamma|`. (Historical: while the candidate was in use, `B_LS` grew from 188 to
339 MeV for 150 MeV protons; with the checked path bound `B_LS = S_bar_max B_L` = 8.95e3 MeV.)

**Admissible lookup magnitude (fail closed, review amendment).** The largest value of a lookup
table over the species a request selects must lie in `[2^-60, 2^60]` (or the table is all zero,
rejected separately). Then the FE exponent `k = 30 - ceil(log2 f_max)` is in `[-30, 90]`, the scale
`2^k` and the bound `f_max E_hi` are finite float64 and the accumulator capacity check is
meaningful; `validate()` rejects anything else with a message naming the range (rescale the table
and its units).

**Precision floor (fail closed).** `q_c <= 2^-16 u_c` with `u_c = 1 mm S_ref^j`, `j` the power of S
in the kind; that bounds the relative rounding of a 1 mm entrance piece by 1.5e-5. A configuration
that violates the floor raises before transport. Example (water, protons, 150 MeV, E_cut 2 MeV;
ramp envelope S_bar_max 45.1 MeV/mm (see the ramp envelope below), B_L about 198 mm, B_LS = S_bar_max B_L about
8.95e3 MeV, B_LS2 about 4.04e5): hpb 1e6 gives `q_L = 2^-34`,
`q_LS = 2^-28`, `q_LS2 = 2^-23`; hpb 1e7 gives `2^-31`, `2^-25`, `2^-20`; hpb 2.86e7 (the E limit)
gives `2^-29`, `2^-24`, `2^-18`. All are far inside the floor; 1e8 hpb already fails closed on
E by the existing rule.

**Residuals.** Every channel except N has one per-history float64 tally column `sum(x - n q_c)`,
reduced by the exact-sum expansion (`tally_rows` becomes `(n, 6 + 2G + C_res + 2)`: one residual
column per non-count channel, the lookup out-of-domain count and the path-bound-exceeded flag). They give global
closures at floating-point precision and report the quantization magnitude. Ratios have no
residual; their quantization error is approximated to first order by
`|dR| ~ (n_v q_X/2 + R n_v q_Y/2) / Y_v` (an approximation, not a bound: it substitutes the quantized R). The reported, deterministic `rounding_bound` uses the interval form `max |(X+a)/(Y-b) - R|, |(X-a)/(Y+b) - R|` with `a = n_v q_X/2`, `b = n_v q_Y/2`, `+inf` if `Y - b <= 0`, review amendment), with `n_v` the piece count of the classes of the respective channel (N, N_local; both are automatic).
Fallback if V3-012 measures a CUDA cost: residuals for E-kind channels only (this would drop the
A1/A3 global closures).

**Effect of the rounding noise.** Zero-mean, sd `q/sqrt(12)` per piece; relative per voxel
`q/(xbar sqrt(12 n_v))`, about 2e-11 (LS) and 2e-10 (LS2) on a plateau voxel at hpb 1e7 with 1e4
pieces, at most 1e-8 for a distal 10-piece voxel. The statistical error of LET_d is 1e-3 to 1e-2, so
the rounding variance is about 1e-12 of it and the induced ratio bias is about 1e-20. The batch
delta-method variance already contains the rounding noise. Pieces below `q/2` round to zero; that
bias is at most `q/2` per piece and is part of the residual. Deterministic tests use explicit
bounds, never z-scores.

**Accepted risk.** Adaptive quanta make the channel bits depend on hpb at the <= 1e-9 level.
Comparisons across different hpb are statistical anyway. `B_L` is valid for primaries only;
V3-005A must re-derive it for secondaries (the convexity argument `sum R(E_i) <= R_p(sum E_i)` for
light ions); the post-run overflow counter backstops it. Fixed quanta per kind were rejected
because `B_c` in air-gap geometries (`r_max` about 1e3) would force unnecessary fail-closed limits.

### 5. Ratio reduction

Linear channels use `reduce_batches`. A ratio `R = X/Y` of channel means uses `reduce_ratio`
with the delta method on the batch sums (`xbar`, `ybar` the means of the batch sums):

`Var R = [Var xbar - 2 R Cov(xbar, ybar) + R^2 Var ybar] / ybar^2`,
`Cov = sum(x_b - xbar)(y_b - ybar) / (B (B - 1))`.

The biased mean of per-batch ratios is not used. A voxel is **defined** iff `ybar > 0` and
`n_nonzero(y) >= max(2, ceil(B/2))`; undefined voxels are NaN with `defined_mask` False, never 0.
`QuantityResult(mean, std, defined_mask, n_nonzero, units, definition)` is stored in
`GridResult.quantities[name]`. Exact derived linear channels (for example `E_step` as the sum of
classes) are integer differences.

### 6. Species registry policy

`ionmc.species` holds an append-only registry with stable integer ids: 0 proton, 1 deuteron,
2 triton, 3 helium-3, 4 alpha; V3-009 appends ions (at most 64). An id is never reused or
reordered. Each entry has a `transported` flag; the pseudo-species `nuclear_local` is not
transported (nuclear-local deposits, V3-005). Channels select species by an `int8`
`species_match[n_ch, n_species]` matrix, not by bitmasks. The set of **producible** species is an
engine capability (in V3-004: the source projectile as a primary, generation 0); a request for a
species that the engine does not produce, for generation "secondary", or for an unknown species
fails closed before transport (no secondary transport until V3-005A). V3-005A expresses its
component arrays as `TallyRequest` objects and only supplies species and generation to the hook
(`score_piece(..., species, gen, cls)`).

### 7. Lookup tables and provenance policy

`LookupTable.from_file(path, expected_sha256=None)` reads JSON with `name`, `quantity`, `units`,
`axis` in {`energy_per_nucleon_mev`, `let_water_kev_um`}, `axis_spacing` in {linear, log},
per-species `values` on a **uniform** grid (non-uniform grids are rejected; an explicit, recorded
`resample_uniform()` helper exists), required `citation`, `license`, `source` and `synthetic`.
Values must be finite and non-negative. Interpolation is linear and reuses the shared
bin-location function. Before transport: all producible species are present, the axis covers
the reachable domain (`[E_floor, E_source,max]` or the corresponding water-S interval) and a
given `expected_sha256` matches the file; any failure raises. At run time an argument outside the
axis increments the fail-closed counter `lookup_out_of_domain`. Models of voxel LET_d (for
example a linear-quadratic RBE model) are post-processing, not lookup tallies.

`Result` carries the provenance of every table: file sha256, content sha256, citation, license,
synthetic flag and resampling record. Persistence is V3-006.

**No clinical tables and no RBE formulas ship.** The only fixture is
`tests/data/synthetic_lookup.json` (`synthetic: true`, `f = 1 + 0.1 L`), a mathematical test
function with no biological meaning.

## Alternatives considered

| Topic | Rejected | Reason |
|---|---|---|
| LET estimator | `sum(eps (eps/l))/sum(eps)` | step-size dependent through fluctuations (section 1); kept only as a negative control |
| S per step | constant `S_w(E_mid)`; lookup per piece per grid | first-order error per piece; higher cost for second-order gain |
| LET medium | local material; configurable | the medium is part of the definition (Kalholm et al.); local stays a later kind |
| Accumulators | float32 per-batch; float64 atomics | not partition-bit-identical; Warp CPU atomics are plain read-modify-write |
| Quanta | fixed per kind | air-gap geometries would fail closed unnecessarily |
| Ratio reduction | mean of per-batch ratios | O(1/n_b) bias, undefined where a batch has no denominator |
| Species selection | bitmask | ion list is open-ended; append-only ids and a match matrix scale |

## Consequences

- The no-tally path is bit-identical to a524f209 (A16); the step precompute is skipped when
  `tallies = ()`.
- Memory grows by `B (sum voxels + sum channel sizes) 8 bytes` per worker; the memory guard raises
  before transport.
- CUDA atomics, register pressure and compile time grow with the channel count: informative only
  (V3-012).
- Distal LET_d has few high-S steps and an O(1/N) ratio bias; the defined mask, `n_nonzero` and the
  coverage test A9 at the distal 80 % point guard it.
- Alignment dependence of the whole-step midpoint assignment is a V3-011 probe.

## Validation strategy (frozen before results exist)

`validation/plans/v3-004-acceptance.md` (rows A1 to A16, committed with this decision, before any
result exists). Gating Monte Carlo comparison of LET is to the TOPAS ProtonLET reference of V3-010;
A13 against Grassberger and Paganetti (2011) is exploratory and non-gating.

## Outcome

- **2026-10-07: acceptance row A4b (first result at 9f52057).** The frozen row did not define
  `E(z)`; with the table CSDA energy `Rinv(R(E0) - rho z)` A4b failed (worst 5.3e-4 / 1.04e-3 /
  1.32e-3 for a maximum step of 1 / 0.5 / 0.25 mm against 1e-4). The scoring is correct: the
  failure is the trapezoid bias of the range table (about 3.5e-5 relative per increment, about
  5 um of the range, decision 0039 outcome of the same date). `LET_t` is the path average of
  `S` of the transported particle, so the reference is that particle's energy, anchored at the
  engine step boundaries; A4b then gives 5.0e-5 / 4.6e-5 / 3.7e-5 and A4 3.8e-5 (plan footnote 1).
  The table-CSDA comparison stays as a non-gating log. The range-table inconsistency is assigned
  to the follow-up V3-003D (separate re-qualification, out of scope here).

- **2026-10-07: acceptance row A9 (observation at bases 20361004 and 20371004, plan amendment 5).**
  At 20361004 the footnote-3 coverage window [0.91, 0.99] passed (0.945 / 0.93 / 0.97); at
  20371004 (head 3c013f8) it failed at the distal 80 % depth (0.935 / 0.91 / 0.90; z sd 1.061 /
  1.119 / 1.108; z mean about -0.01; all 161 informative bins mean coverage 0.9435, min 0.900).
  Analysis of the per-batch archives: the delta-method sigma equals the batch-jackknife sigma and is
  within 0.5 % of a batch bootstrap sigma; sigma_delta / sigma_empirical has mean 0.977 and 1.026
  over informative bins (per-bin sd about 5 %, as expected for 200 seeds) and the distal 80 % ratio
  flips from 0.925 to 1.056 between the bases; pooled over 400 seeds z sd is 1.013 / 1.093 / 1.029 /
  1.023 (0.5R / peak / d80 / d50), mean coverage 0.9464, neighbouring-bin z correlation 0.95. The
  uncertainty is calibrated; the per-depth coverage of 200 seeds carries the about 5 % noise of the
  empirical sigma, which the frozen window ignored (its 1 % false-failure rate was understated).
  Contrary evidence: the Bragg-peak bin alone is suggestively low (pooled sigma ratio 0.927, about
  2 sd), an open observation for a larger-seed peak-only study. The criterion is replaced (plan
  footnote 5) by z sd in [0.87, 1.13] and |mean z| < 0.18 at the three depths (n = 200, 99 %
  intervals) and mean all-bin coverage in [0.93, 0.97]; both observed bases pass it; this
  re-calibrates the criterion's noise and does not widen a physics tolerance. Base 20371004 is
  consumed; the qualification base is 20381004.

To be appended from committed result files.

- 2026-10-07 (V3-004, known limitation recorded before merge): the water row of
  `TransportTables` is indexed by energy per nucleon (`water.energy_per_u`), while the scorer
  (`water_state`, kernel `step_state`) and the capacity-bound proof look it up with the ion's
  total kinetic energy. For protons (A = 1, the only producible species in V3-004) the two
  coincide, so nothing is affected now; when V3-009 appends ions with A > 1 the lookup must use
  E/A (bound and scorer together). Fail-closed guard: `producible()` rejects every other species.
