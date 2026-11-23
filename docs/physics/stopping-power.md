# Electronic stopping power

Implemented in `ionmc.physics.stopping` (model and tables), `ionmc.materials` and
`ionmc.physics.projectiles`. The design decision is
[0038](../generated/decisions/0038-electronic-stopping-power-and-data-roles.md). The
module docstring of `stopping.py` is the authoritative statement of each term.

## Model

The mass electronic stopping power of an ion with charge number z in a material is

    S/rho = K q^2 (Z/A) beta^-2 [ L0 + q L1 + L2 + M ]      [MeV cm^2/g]

with K = 0.307075 MeV cm^2/mol, (Z/A) in mol/g and q the effective charge. All
arithmetic is float64. Energies are kinetic energy per nucleon (MeV/u); the total
kinetic energy of the ion is E_u times the mass number.

| Term | Implementation |
|---|---|
| L0 | `1/2 ln(2 m_e c^2 beta^2 gamma^2 Tmax / I^2) - beta^2 - delta/2 - C/Z`, with the exact `Tmax = 2 m_e c^2 beta^2 gamma^2 / (1 + 2 gamma m_e/M + (m_e/M)^2)` |
| delta | Sternheimer density effect (Sternheimer, Berger, Seltzer, At. Data Nucl. Data Tables 30 (1984) 261); zero if the material has no parameters. Parameters are provided for water (x0 0.2400, x1 2.8004, Cbar 3.5017, a 0.09116, m 3.4773) |
| C/Z | Shell correction in the Bichsel parameterisation of ICRU Report 49 as written in the PDG review (polynomial in eta = beta gamma, I in eV, valid for eta >= 0.13, about 8 MeV/u for protons). Below that limit the value at eta = 0.13 is held. For compounds, C(I_i)/Z_i of each element (elemental I) is weighted by its electron fraction |
| q L1 | Barkas term in the Ashley-Ritchie-Brandt form (Phys. Rev. B 5 (1972) 2393; ICRU 49) as implemented in Geant4 `G4EmCorrections::BarkasCorrection` (v11.4.2); the tabulated function F(W) is the 47-point table of that source, interpolated linearly |
| L2 | Bloch term `-y^2 sum_{n>=1} 1/(n (n^2 + y^2))`, `y = q alpha / beta` |
| M | Mott term `1/2 pi alpha beta q` for projectiles with z >= 2 (as in Geant4) |
| q | Pierce-Blann effective charge `z (1 - exp(-125 beta z^(-2/3)))` for z >= 2; q = z for z = 1 |

Because Barkas, Bloch and Mott depend on the actual z, helium, carbon and oxygen
stopping powers are not z^2-scaled proton stopping powers. Each correction can be
switched off through `BetheOptions` (or the keyword arguments of
`bethe_mass_stopping`); the options are stored in the table metadata.

## Domain

Energies below 1 MeV/u raise `ValueError`: the model does not return stopping
powers there, and transport is expected to deposit the residual energy locally.
Between 1 MeV/u and about 8 MeV/u (protons) the shell correction is held at its
boundary value, so the model is less accurate in that interval; the range effect is
confined to the last part of the track.

## Tables

`BetheStoppingSource(options, e_min_per_u, e_max_per_u, points_per_decade).table(material,
projectile)` returns a `StoppingTable` on a log-spaced grid (defaults 1 to 500 MeV/u,
200 points per decade) with `s_el_mass` [MeV cm^2/g], `s_el_linear` [MeV/mm] =
S rho / 10, `csda_range_g_cm2` and `range_mm`. The CSDA range is the exact integral
of the log-log interpolated `a E_u / S` over ln E_u (closed form per grid interval, decision 0039
V3-003D, `range_construction = exact-loglog-quadrature-v1`; the trapezoid rule overestimated every
increment by 2e-5 to 4e-5) at the grid nodes, starting from `a E_min / S(E_min)`, which is the integral of
dE/S below E_min for a constant stopping power equal to S(E_min). This is an approximation,
not a bound: the stopping power rises below E_min only down to its low-energy maximum
(around 0.1 MeV/u for protons in water) and falls again at lower velocity, so the true
residual range can be smaller or larger. The offset is of the order of 0.02 mm of water
for protons and is recorded in the table metadata; transport deposits the energy below
E_cut locally, so the offset does not enter a transported range. Stopping
power, range and the inverse `energy_from_range` use one log-log piecewise-linear
interpolant, so `energy_from_range(range_at(E))` returns E to rounding error. `range_at` is log-log
interpolated between the nodes (within 1.2e-5 relative of the exact integral); the transport tables
evaluate the exact closed form between nodes instead.
`NistStarStoppingSource` produces the same structure from a cached PSTAR (protons) or
ASTAR (alpha) liquid-water table, starting the range from the NIST CSDA range at E_min.
`NistStarStoppingSource` accepts a material only if it is structurally liquid water
(elements H and O, mass fractions within 1e-4 of water's, density within 1e-6 g/cm^3 of 1);
names are not used. The table metadata records `requested_material` and
`effective_material` (NIST liquid water, I = 75 eV). `Material` validates its inputs
(finite positive density and I, positive finite fractions summing to 1 within 1e-6) and
exposes its mass fractions as a read-only mapping; tables reject non-finite values.
Both implement the `StoppingSource` protocol; the table metadata records the source and
the I value.

## Mean excitation energy of water

`ionmc.materials.WATER` uses I = 78 eV (ICRU Report 90, 2016). `water(I_eV)` builds a
variant, for example `water(75.0)` for the ICRU 49 / NIST PSTAR/ASTAR value. Material I
is configurable per material; if it is `None` the Bragg-additivity value is used.
Comparisons with other engines must state each engine's I.

## Data roles

| Dataset | Role |
|---|---|
| NIST PSTAR / ASTAR liquid water (I = 75 eV) | evaluation of the analytic model at I = 75 eV; construction data when `NistStarStoppingSource` is selected |
| ICRU 90 water arrays (from the Geant4 source) | evaluation of the analytic model at I = 78 eV; construction data when `icru90-water-v1` is selected (see the registry section below) |

PSTAR, ASTAR and ICRU 90 share the Bethe-theory lineage above about 1 MeV/u; agreement
with them is tabulated, not measured, evidence. A run that uses the NIST tables as its
stopping source cannot be evaluated against the same tables.

The energy grid of the ICRU 90 alpha array is interpreted as total alpha kinetic
energy (same span as ASTAR); `compare_nist.py` reports the ratio of the two datasets
on their shared energies so that this reading is checked on every run.

## Comparison script

After fetching the three datasets (see [data layer](../architecture/data-layer.md)):

```bash
uv run python validation/scripts/stopping/compare_nist.py --output-dir OUT_DIR
```

It reads only the cache (no network unless `--online`) and writes
`stopping_comparison.json`, which contains aggregates only and no per-energy values: for
each comparison (Bethe at I = 75 eV against PSTAR/ASTAR, Bethe at I = 78 eV against the
ICRU 90 arrays, ICRU 90 against NIST on shared energies) the number of points, the energy
range and the maximum and RMS relative deviation of the stopping power for E >= 2 and
E >= 10 MeV/u; for CSDA ranges at 100/150/200/250 MeV protons and 100/150/200 MeV/u alpha
particles the maximum absolute relative deviation (Bethe 75 eV against NIST, Bethe 78 eV
against the integrated ICRU 90 values); and the model-only Bethe 78 - 75 eV range shift per
energy. Results are produced by running the script; this page does not quote them.

## Stopping-model registry (V3-013)

Decision: [0044](../generated/decisions/0044-stopping-source-selection.md). Code:
`ionmc.physics.stopping_models` (registry, `resolve`, `spec_of`),
`ionmc.physics.stopping_selection` (selection point, domain checks, identity block),
`ionmc.physics.stopping_evidence` (evidence records, circular rows). Every stopping-power model
that a run can use is a registered, immutable `StoppingModelSpec` with a model id, a version, a
domain and data pins; its evidence status is kept in a separate append-only registry keyed by the
`spec_sha256` of the spec, so a grading never changes a spec or a table hash. Not a clinical claim.

### Registered models

| model id | kind | species; material; window (MeV per nucleon) | data pins (dataset id, sha256 prefix, licence, use role) | evidence record |
|---|---|---|---|---|
| `bethe-ionmc-v1` (default) | analytic | proton, deuteron, alpha, carbon12, oxygen16; any `materials.py` material (evidence: liquid water only); 1 to 500 | none (analytic; I-value from the material, water 78 eV) | `as-graded:0043` (V3-008 rows as graded, I3 FAIL retained) |
| `bethe-ionmc-v1+ls-eff` | analytic (sensitivity) | as the default; effective only for z >= 2 (for z = 1 the effective id is `bethe-ionmc-v1`) | none | `unvalidated` |
| `bethe-ionmc-v1+ls-bare_mixed` | analytic (sensitivity, attribution) | as `+ls-eff` | none | `unvalidated` |
| `nist-pstar-astar-water-2005` | tabulated | proton (PSTAR), alpha (ASTAR); structural liquid water; 1 to the table end, capped at 500 (alpha ends at 250 per nucleon); nuclear runs and C/O rejected; I = 75 eV fixed by the table, requested I recorded | `nist-pstar-water-2005` f5fc6482, `nist-astar-water-2005` dd4a2e50; NIST SRD 124, use-only, downloaded by each user; construction-permitted (the pin of the transported projectile only is used in a run) | `consistency-only`; circular for the PSTAR/ASTAR rows |
| `icru90-water-v1` | tabulated | proton (1 to 500 MeV), alpha (1 to 250 per nucleon); structural liquid water with requested I = 78 eV (otherwise fail closed); no deuteron (nuclear runs fail closed), no C/O | `geant4-icru90-stopping-11.4.2` afd14794 (single file); Geant4 Software License; construction-permitted | `consistency-only`, `topas-agreement-shared-lineage:V-C1:PASS` (V-C1 graded once at 9f1fa8ec, 2026-11-23: He 75/125/200 MeV/u vs TOPAS opt4, shared ICRU 90 lineage; transport/straggling/end-of-range agreement only, never `validated`) |

`icru90-water-v1` interpolates log-log onto 200 points per decade from 1 MeV/u and builds the range
with `exact-loglog-quadrature-v1`; the start range below 1 MeV/u is the exact log-log integral of
the array's own nodes (`icru90-integrated-from-lowest-node-v1`). The registry also lists the
excluded construction sources: G4EMLOW `ion_stopping_data` (lineage-only), libdEdx ICRU73_NEW
(evaluation-only; construction needs an operator authorization first), SRIM/MSTAR and z^2-scaled
PSTAR/ASTAR (rejected, decisions 0038 and 0043). `resolve` cannot construct a model from a pin whose
use role is `evaluation-only` or `lineage-only`. Composite models (table plus analytic with a join)
are not registered.

**Licence role of the ICRU 90 arrays.** The arrays are Geant4-licensed source-embedded data
(decision 0038; the project's assessment, not legal advice). In `icru90-water-v1` runs they are
construction data: they are read at run time from the user's verified provenance cache, never
committed and never redistributed, the ICRU report is cited and not reproduced, and the Geant4
notice is shipped in `src/ionmc/THIRD_PARTY_NOTICES.md`. The dataset's own role in the registry stays
`evaluation`; the per-use role is carried by the model's data pin. In every other run the arrays
remain an evaluation reference. No value derived node by node from them is committed.

**Circularity.** A run whose construction data are also the evaluation reference of a row makes that
row circular for that configuration. For `icru90-water-v1` these are the ICRU 90 part of row I2 and
the ICRU 90 comparisons of decision 0038; every TOPAS helium comparison above the BraggIon switch
and every TOPAS proton comparison above 2 MeV shares the same array (Geant4 11.4.2 `opt4` uses the
ICRU 90 data there) and therefore tests transport, straggling and end of range only. For
`nist-pstar-astar-water-2005` the PSTAR/ASTAR rows are circular. The rows are recorded in the
persisted `circular_rows`.

### Documented default

> Default `bethe-ionmc-v1`: analytic Bethe stopping with Barkas-125 effective charge,
> Ashley–Ritchie–Brandt Barkas term (Geant4 `G4EmCorrections` F table), Bloch, first-order Mott,
> Bichsel shell correction, Sternheimer density effect, Lindhard–Sørensen off, I_water = 78 eV
> (ICRU 90). It is chosen for consistency with the qualified proton engine and lineage alignment
> with the specified TOPAS configuration (shared I-value, ICRU 90 p/α agreement ≤ 0.2 % above
> 10 MeV/u, shared Barkas table). It is **not TOPAS-equivalent**: TOPAS (opt4) uses the ICRU 90
> alpha array for He, the ICRU 90 proton array at equal velocity × q_eff² plus Lindhard–Sørensen
> and Barkas terms for C/O, and the ICRU 90 proton array for protons above 2 MeV; measured
> TOPAS/ionmc thin-slab differences are +3.5/+5.0 % (C/O) at 2–5 MeV/u, falling to +0.1/+0.3 % at
> 20–50 MeV/u (I3-FU-1 step 2, INCONCLUSIVE by rule); row I3 (C/O vs revised ICRU 73 table) FAILED
> and is retained as a declared low-energy C/O limitation.

Not a clinical claim. The default, its options, its I-value and every V3-008 verdict are unchanged
by V3-013; the default's table identities and hashes are byte-identical to those before the registry
existed.

### Selecting a model

```python
from ionmc.config import PhysicsOptions
from ionmc.physics.stopping import BetheStoppingSource

PhysicsOptions(nuclear=False, stopping_model="icru90-water-v1")  # a registered id
PhysicsOptions(nuclear=False, stopping=BetheStoppingSource())    # a source object
```

Exactly one of `stopping_model` and `stopping` must be given (neither or both raises
`UnsupportedCombinationError`; `PhysicsOptions` has no implicit stopping choice, and `ionmc.physics.stopping.default_source()`
returns `resolve("bethe-ionmc-v1")`). A `stopping` object is identified from its instance fields, never from
its class: an exact match with a registered spec is recorded under that model id
(`BetheStoppingSource()` is `bethe-ionmc-v1`), any other source (for example
`BetheStoppingSource(e_min_per_u=1.5)` or one with `barkas=False`) is recorded as
`model_id_effective = "unregistered"` with its complete options, `decision = "user"` and the
evidence record `unvalidated`; it is never reported under a registered id. A user-written source
declares itself with a `stopping_spec()` method that returns a user spec (without data pins);
an object that neither `spec_of` supports nor declares itself this way is refused. Tabulated
models read their data from the verified cache only (no download at run time); a missing object, a
hash mismatch or a parse failure raises.

`validate` derives every stopping-dependent object from the one selected source: the material
tables, the LET water row, the deuteron tables of nuclear runs, the residual range and step bounds,
and the exact-Omega ion rows. It asserts that all table identities name one effective model.
`ExactOmegaMapping.build` and `fdm_domain_status` fail closed on a missing ion table; there is no
Bethe default in a run path. The MCS proton path (the b-matching proton path of the ion
multiple-scattering model, 1 to 4000 MeV) is fixed at `bethe-ionmc-v1` for every selected model and
is recorded as `mcs_proton_path_stopping`; for a non-default ion run this is a declared mixed
construction. `ionmc.capabilities()["physics"]["stopping_models"]` lists every registered model with domain,
data pins, licence role and evidence records.

### Fail-closed domains

Before transport, `validate` checks that the species is in the model's domain, that every geometry
material and the LET water row satisfy the material rule (structural, by composition and density,
not by name; for `icru90-water-v1` also the requested I-value), and that every energy interval used
lies inside the model's per-species window: table floor and cutoff, the LET water row, the deuteron
tables of nuclear runs and, for registered models, `E0 + 6 sigma_E`. A miss raises
`UnsupportedCombinationError` naming model, species, material, energy and domain. There is no
extrapolation, no fallback model and no automatic composite. Amendment A-1 of the plan: for an
`unregistered` user source the `E0 + 6 sigma_E` check is not made at `validate`; the pre-existing
runtime rule applies unchanged (`source_energy_out_of_range` counter, result invalid), and the
table-range, LET-row and deuteron checks still apply to it.

### Persisted identity block

`EffectiveConfig.summary()["stopping"]` is always present and lies outside the table hash. Keys:
`model_id_requested`, `model_id_effective` (equal except the z = 1 Lindhard–Sørensen no-op, where
the effective id is the default's), `version`, `kind`, `spec_sha256`, `options` (complete, with
`ls_requested`/`ls_effective` and the I-value rule), `interpolation`, `range_construction`,
`range_start`, `domain`, `data_declared` (every pin of the spec), `data_used` (only the pins actually
read for the transported projectile and the materials of this run, each with `role_in_this_run =
"construction"` and its table; a proton run of `nist-pstar-astar-water-2005` lists PSTAR only, an
alpha run ASTAR only), `content_sha256` of the parsed arrays, `evidence_at_run` (the evidence
records of the effective spec at run time), `decision`, `circular_rows`, `selected_by`
(`"stopping_model"` or `"stopping-object"`), `species_transported`, and for ion runs
`mcs_proton_path_stopping`. Tabulated sources add `domain_resolved`. The block travels through the
validation scripts' JSON step documents; there is no library-level result writer.

Models registered before decision 0044 keep byte-identical table metadata (so the T1 fixture and the
A16 digests are unchanged); every model registered later carries `model_id` and `version` in its
table metadata, so its table hash differs from every existing model.

### No silent switch

There is no fallback chain, no "best available" choice and no environment-variable override.
`DEFAULT_STOPPING_MODEL_ID` and the default's `spec_sha256` are pinned by a test; a change requires
a `STOPPING_DEFAULT_INTENDED_CHANGE` record naming a decision (the record is `None`). No validation
outcome changes the default or a registered spec: a failure is recorded by appending
`failed:<row>` to the evidence registry, and a revised alternative is a new model id or version with
fresh, held-out evidence pre-registered before it runs.

### Evidence statuses

| status | meaning |
|---|---|
| `unvalidated` | no evidence of any kind (Lindhard–Sørensen sensitivity variants; every `unregistered` user source) |
| `consistency-only` | the model reproduces its own input tables and passes implementation checks; this says nothing about physical accuracy |
| `topas-agreement-shared-lineage:<row>:PASS` | transport agreement with TOPAS for the named row; the stopping data are shared with the reference, so this tests transport, straggling and end of range only, never the stopping physics |
| `failed:<row>` | a gating row of the model failed; the spec is unchanged |
| `as-graded:<decision>` | the default: its V3-008 rows as graded in the named decision, the I3 FAIL included |

No status is named `validated`. Evidence applies to its own model id, options and domain and does not
transfer between models.

### Consistency checks of `icru90-water-v1` (not evidence)

These checks show that the implementation reproduces its input arrays. They are not validation of
the ICRU 90 stopping physics. Aggregates are committed in `validation/results/stopping/v3-013/`
(no table value is stored).

| check | result |
|---|---|
| CC-NODE: transport table against the log-log interpolant of the array at the fine-grid nodes | maximum deviation 0 (proton and alpha) |
| CC-I2: range increments R(E_i) - R(E_0) against the I2 integrated-ICRU 90 construction, node by node | maximum relative deviation <= 7e-16 |
| chord error of the fine table at the ICRU nodes themselves (aggregate, report-only) | maximum 2.2e-4 (proton), 1.3e-4 (alpha) |
| CC-SPAN: counts and coverage of both arrays | 57 and 49 nodes; windows inside the spans |
| V-C3 (report-only): `icru90-water-v1` against PSTAR (proton) and ASTAR (alpha) over the common window, stopping | maximum 1.29 % (proton), 1.23 % (alpha) |
| V-C3: CSDA range | maximum 1.49 % (proton), 2.36 % (alpha) |

The V-C3 spread between two tabulations is reported, not graded. CC-INT (interpolation error between
nodes against a monotone cubic in log-log space) is likewise report-only.

### Pending validation of `icru90-water-v1`

V-C1 (He at 75, 125 and 200 MeV/u against TOPAS 4.3 / Geant4 11.4.2 `g4em-standard_opt4`, EM only),
V-C2 and V-C2-D (the selection reaches transport as predicted from the tables alone) have **not
been run**. They are sealed before any run and are graded once. V-C1 is a shared-lineage transport
agreement: TOPAS uses the same ICRU 90 alpha array above 7.9452 MeV, so a pass supports
transport, straggling and end-of-range consistency of this configuration, never the stopping
physics. V-C1 for protons is not pre-registered; proton use of `icru90-water-v1` has consistency
checks only. Not claimed: that any model is TOPAS-equivalent; a carbon or oxygen alternative; measured
He/C/O range, lateral or stopping data (experimental validation is deferred by operator directive).
See the [validation page](../validation/index.md).

## Sources

Third-party licence texts and data attributions: `THIRD_PARTY_NOTICES.md` in the repository root (full text in `src/ionmc/THIRD_PARTY_NOTICES.md`, shipped with the package; see [third-party notices](../third-party-notices.md)).

- ICRU Report 49 (1993), Stopping Powers and Ranges for Protons and Alpha Particles.
- ICRU Report 90 (2016), Key Data for Ionizing-Radiation Dosimetry (J. ICRU 14(1)).
- Particle Data Group, Passage of particles through matter.
- M.J. Berger, J.S. Coursey, M.A. Zucker, J. Chang (2005), ESTAR, PSTAR, and ASTAR
  (version 2.0.1), NIST, https://doi.org/10.18434/T4NC7P.
- Geant4 11.4.2: `G4EmCorrections.cc`, `G4ICRU90StoppingData.cc`, `G4NistMaterialBuilder.cc`.
