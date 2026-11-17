# 0044 — Explicit selection of stopping-power models and tabulated sources

- Status: DRAFT (S0, planner, 2026-11-17; revised 2026-11-17 after the pre-freeze reviews, section
  "Pre-freeze review log"), frozen together with `validation/plans/v3-013-acceptance.md` by one
  reviewed commit before any V3-013 code or comparison output exists. Until that commit this is a freeze candidate.
- Task: V3-013 (follow-up directed by the operator in the V3-008 integration authorization IR-20261116-204824-29B867,
  option A; directive quoted in decision 0043, dated note of 2026-11-17, and in `validation/plans/v3-008-acceptance.md`,
  "Integration decision").
- Base: `v3/develop` at a9563412 (a956341291fd30ab47ad2216fe3e235db78beda0, the V3-008 squash merge, PR #84).
- Affects: architecture (stopping-source selection, identity, capability contract), data provenance and licence
  roles, validation strategy. **No physics default changes.**
- Supersedes: only the *selection mechanism* of decision 0038 ("Tabulated layer": the free-form
  `stopping_source = "bethe" | "nist-star"`) and of decision 0043 §1 ("Selected": the source object passed as
  `PhysicsOptions.stopping` with no registry, no domain, no identity), and amends decision 0038's data-role
  record for the ICRU 90 arrays ("only as a cached evaluation reference"; data-role table row `evaluation`):
  in `icru90-water-v1` runs they are construction data (section 1, licence-role extension); in every other
  run their role is unchanged. The default model, its options, its I-value,
  its domains and every V3-008 verdict are preserved (section "Preserved obligations").
- Markers: **FACT** (read in code or records at a9563412, file:line), **ESTIMATE**, **TO VERIFY**, **PROPOSAL**
  (a design choice frozen by this decision). Not a clinical claim.

## Problem

The operator directs (2026-11-17): explicit selection of stopping-power models or tabulated sources without forcing
one universal model; a documented default aligned with the TOPAS comparison objective without claiming TOPAS
equivalence; alternatives with explicit species/material/energy domains, source and licence provenance and
configuration-specific validation (reproducing an input table is a consistency check, not independent physical
validation); consistent application to range construction, transport and dependent scoring; the effective model or
table identity persisted in outputs; no silent switching and no model chosen after seeing a validation outcome to
obtain a pass. Architecture, candidate sources and sequence are chosen autonomously.

## Context (FACT unless marked)

**Selection today.**
- `PhysicsOptions.stopping: StoppingSource` has no default (`src/ionmc/config.py:159`); `__post_init__` checks only
  `.name`/`.table` and the `ls` value (`config.py:193–197`). `default_source()` returns `BetheStoppingSource()`
  (`src/ionmc/physics/stopping.py:638–640`).
- Two source classes: `BetheStoppingSource` (`stopping.py:472–533`; defaults `BetheOptions()` = shell, barkas, bloch,
  mott, effective_charge all True, `ls="off"` (`stopping.py:123–128`); `e_min_per_u=1`, `e_max_per_u=500`,
  200 points/decade, `name="bethe"`; start range constant-S `a·E_min/S(E_min)`) and `NistStarStoppingSource`
  (`stopping.py:585–`; `name="nist-star"`; structural water check `_require_nist_water`, dataset check
  `_require_water_table`; grid upper end `min(e_max_per_u, last table row)`, `stopping.py:609`).
- The I-value is a property of the material (`materials.py`, `WATER = water()`, 78 eV default), not of the source;
  `nist-star` records the requested I and uses the table's 75 eV.
- The `ls` key enters the Bethe metadata only for z ≥ 2 and `ls ≠ "off"` (`stopping.py:524–528`), so default and proton
  metadata are unchanged by the LS option (decision 0043 §1 amendment).

**Consumers of the selected source (follow the selection).** `config._stopping_tables` (`config.py:682–696`) →
`TransportTables.from_stopping_tables` (range, inverse range, step loss, `continuous_residual_v1` residual range,
`max_steps`, table energy bounds; `config.py:1059–1062`); `config._water_table` (LET water row, decision 0040;
`config.py:789–807`); deuteron tables of nuclear runs (`config.py:912–923`); the exact-Ω ion side (transport rows
packed by `_ion_mapping_data`, and `fdm_domain_status(..., ion_tables=…)` built from the same stopping tables,
`config.py:1103–1111`).

**Paths that bypass the selection (hard-wired Bethe).**
- MCS proton path: `exact_omega.proton_path()` always builds `BetheStoppingSource(e_min_per_u=1, e_max_per_u=4000,
  200/decade)` (`src/ionmc/physics/exact_omega.py:321–335`; `PROTON_PATH_EMAX_MEV = 4000`, L133); used by
  `config.py:442` and `exact_omega.py:794, 1023`.
- `ExactOmegaMapping.build(..., ion_table=None)` falls back to `BetheStoppingSource()` (`exact_omega.py:788–802`);
  `fdm_domain_status` passes `ion_tables.get(fingerprint)`, which is `None` for a missing material
  (`exact_omega.py:1022–1023`) — a latent silent fallback (not reachable from `validate`, which supplies every
  geometry material).
- Diagnostics and fixtures: `transport/parity.py:953`, `transport/ion_mcs_checks.py:42`, `nuclear/build.py:1041`
  (gate D6), `validation/scripts/transport/a16_digest.py:51, 155`, `tests/ionmc/test_transport_diagnostic.py`
  (T1 trace fixture).

**Identity and hash.** `TransportTables.identity` per material carries `source`, I requested/effective, `dataset_id`,
`source_sha256`, `content_sha256`, `range_construction` and the full `metadata` (`src/ionmc/transport/tables.py:267–288`);
`sha256` hashes the arrays and the identity (`tables.py:304`). `EffectiveConfig.summary()["tables"]` holds the sha256,
`stopping_source` (the source name), the identities and the water row (`config.py:635–644`); `ls_requested/ls_effective`
are outside the hash (`config.py:612–613`). No model id, version, domain, licence role or data-member pin exists.

**Q4 — what the R2 regression depends on.**
- A16 digests (`validation/scripts/transport/a16_digest.py`): sha256 of the run *outputs* only — per-grid batch
  energy quanta, every energy-balance field, counters, every diagnostic array (L140–148). The source is built directly
  as `BetheStoppingSource()` (L51) and enters only through the table arrays. They do **not** depend on
  `TransportTables.sha256` or on table metadata; `table_identity()` (L152–156) reads only `range_construction`, and only
  in the intended-change mode of `steps_v4.py a16`. `A16_INTENDED_CHANGE` is `None` (asserted,
  `tests/ionmc/test_run_suite.py:1014`).
- T1 trace fixture (`tests/ionmc/test_transport_diagnostic.py`, `data/trace_baseline_20mev.npz`): the portable table
  identity **is asserted** (`current == meta["table_identity"]`, L99, with `current = thaw(tables.identity)`, which
  includes the full metadata); `transport_tables_sha256` is compared softly and never asserted (L102–105; machine-
  dependent libm ulps). Consequence: any change of the default model's table metadata breaks the T1 fixture;
  the default's metadata must stay byte-identical.

**Q5 — persistence.** No on-disk result writer exists in `src/ionmc` (only the nuclear table builder writes its own
info JSON, `nuclear/build.py:1202`). Results live in `Result.effective_config` (`src/ionmc/simulation.py:319`); the
validation scripts persist `effective_config.summary()` in JSON step documents (`validation/scripts/transport/
steps_v8.py:209, 333`), and `validation/scripts/ion_em/r2_io.py:1742–1745` binds a step document's `effective_config`
to the frozen run. The identity block therefore lives in `EffectiveConfig.summary()` and travels through those JSON
writers; a library-level result writer (V2-OUTPUT) is not created by this task.

**Q7 — stopping lineage of opt4 protons in G4_WATER.** FACT at the source and printout level, one INFERRED premise:
- [BB] `G4BetheBlochModel.cc` v11.4.2 (pinned 5dc2a900…, `validation/plans/v3-008/p2b-follow-ups.md` §0–§1): L104–107
  set the ICRU 90 data for `pname == "proton"`; L273 (below `fProtonTlimit` = 10 GeV, L79) and L290–293 return
  `GetElectronicDEDXforProton` for protons in an indexed material, returning at L300 before the Bethe formula.
- Archived TOPAS 4.3 / Geant4 11.4.2 opt4 printout (run REF-e67635b8fbee212c22b2-664b209a, full-physics smoke,
  G4_WATER; `.ionmc-cache/reference-runs/…/stdout.txt`): L78 `Use ICRU90 data 1`; L231–237 `hIoni` for proton:
  `Bragg` 0–2 MeV, `BetheBloch` 2–600 MeV (table limits of that run).
- Hence opt4 protons above 2 MeV in G4_WATER use the ICRU 90 proton array. INFERRED: G4_WATER is indexed —
  `G4ICRU90StoppingData.cc` (cached, afd14794…) L82–102 registers G4_AIR, G4_WATER, G4_GRAPHITE by name; the index
  lookup `GetIndex` and the interpolation `GetDEDX` are in `G4ICRU90StoppingData.hh` (not pinned; **TO VERIFY**).
  Below 2 MeV (`G4BraggModel`): whether ICRU 90 data are used is **TO VERIFY** (`G4BraggModel.cc` v11.4.2 not pinned);
  irrelevant to transport above the 2 MeV proton cutoff.

**Q8 — ICRU 90 array spans** (cached `G4ICRU90StoppingData.cc` afd14794…, parsed by `src/ionmc/data/icru90.py`;
counts and coverage only): proton grid 57 nodes, extending beyond 1–500 MeV at both ends; alpha grid 49 nodes, ending
at 1000 MeV total (250 MeV/u, as `icru90.py` states) and starting below 1 MeV/u. Both stopping arrays have the grid
lengths. No value is reproduced here.

**Q6 — helium energies and seeds used by V3-008** (frozen bundle list `validation/results/transport/v8r2-topas/
frozen-bundles.json`; plan v3-008 §C): TOPAS He depth-dose/lateral cases at 50, 150, 250 MeV/u (TD-EM/TL-EM/TD-FULL/
TL-FULL main set, cases 1–3) and 100 MeV/u (held-out TD-H/TL-H, case 10); V5-SS He thin slabs at 150 and 250 MeV/u;
I3-FU-1 He thin slabs at 10 and 50 MeV/u; He 220 MeV/u was declared for the superseded V1-held row (never run).
TOPAS seeds 20619000 + 100·part + case; ionmc held-out samples base 20621004; i3fu2 20627200 + 100·part + case.

## Selected

### 1. Stopping-model registry (PROPOSAL; new module `src/ionmc/physics/stopping_models.py`)

```
StoppingDomain(species: frozenset[str], material_rule: str,
               e_min_per_u: Mapping[str, float], e_max_per_u: Mapping[str, float])
DataPin(dataset_id: str, archive_sha256: str, member: str | None, member_sha256: str | None,
        parser: str, licence: str, use_role: "construction-permitted" | "evaluation-only" | "lineage-only")
StoppingModelSpec(model_id: str, version: str, kind: "analytic" | "tabulated", options: Mapping,
                  interpolation: str, range_start: str, data: tuple[DataPin, ...], domain: StoppingDomain,
                  decision: str)                       # immutable; carries NO evidence status
STOPPING_MODELS: Mapping[str, StoppingModelSpec]      # immutable
DEFAULT_STOPPING_MODEL_ID = "bethe-ionmc-v1"
spec_sha256(spec) -> str                               # sha256 of the canonical JSON of the spec
EvidenceRecord(spec_sha256: str, status: str, row: str | None, plan: str, graded_at_sha: str | None)
STOPPING_EVIDENCE: Mapping[str, tuple[EvidenceRecord, ...]]   # keyed by spec_sha256; module
                                                       # stopping_evidence.py; append-only
resolve(model_id: str, *, cache_dir: str | None = None) -> StoppingSource
spec_of(source: StoppingSource) -> tuple[str, StoppingModelSpec]   # id from the instance fields (section 2)
```

- `resolve` returns the existing classes (and one new `TabulatedStoppingSource` for parsed arrays); every returned
  source carries `spec` and `identity(material, projectile) -> dict`. `resolve` of an unknown id raises; it never
  returns a source of another id. A spec's `version` changes with any numeric change of its formula, data or code path.
- Composite models (table plus analytic with a join) are allowed only as explicit, separately named models with a
  frozen join rule; none is registered by this decision.
- **Evidence status is not part of the spec** (Codex REVIEW-0f80fad8 finding 1). It lives in the append-only
  evidence registry keyed by `spec_sha256`, with the statuses `unvalidated`, `consistency-only`,
  `topas-agreement-shared-lineage:<row>:PASS`, `failed:<row>` (a FAIL of any gating row of the model) and
  `as-graded:<decision>` (the default: its V3-008 rows as graded, I3 FAIL included). A grading appends one
  record; it never changes a spec, a `spec_sha256` or a table hash. `capabilities()` joins the two registries;
  the persisted block carries the records present at run time (`evidence_at_run`). No status is ever named
  `validated`: shared-lineage agreement is recorded as such (finding 2).

**Registered models (frozen by this decision):**

| model_id | kind | options / data | domain (species; material; table window per u) | initial evidence record |
|---|---|---|---|---|
| `bethe-ionmc-v1` (**default**) | analytic | exactly `BetheOptions()` defaults: shell, barkas, bloch, mott, effective_charge on, `ls="off"`; Barkas-125 effective charge; grid 1–500 MeV/u, 200/decade; start range constant-S; I from the material (water 78 eV, ICRU 90) | proton, deuteron, alpha, carbon12, oxygen16; any `materials.py` material (evidence: water only); 1–500 MeV/u; source domains of decision 0043 §5 and the proton rules of `validate` | `as-graded:0043` (decisions 0038/0039/0043 and plan v3-008; I3 FAIL retained) |
| `bethe-ionmc-v1+ls-eff` | analytic (sensitivity) | as default with `ls="eff"` | as default; effective only for z ≥ 2 (z = 1: documented no-op, effective id `bethe-ionmc-v1`) | unvalidated (sensitivity) |
| `bethe-ionmc-v1+ls-bare_mixed` | analytic (sensitivity/attribution) | as default with `ls="bare_mixed"` | as `+ls-eff` | unvalidated (sensitivity/attribution only, decision 0043 §1) |
| `nist-pstar-astar-water-2005` (today `"nist-star"`) | tabulated | registry `nist-pstar-water-2005` f5fc6482…, `nist-astar-water-2005` dd4a2e50…; parser `nist_star`; NIST SRD 124 terms (decision 0038 legal assessment); I = 75 eV fixed by the table, requested I recorded | proton (PSTAR), alpha (ASTAR); structural liquid water; 1 MeV/u to the table end (≤ 250 MeV/u alpha; PSTAR upper end **TO VERIFY** by a loader test, counts/coverage only); nuclear and C/O rejected | consistency-only (decision 0038); circular for PSTAR/ASTAR rows |
| `icru90-water-v1` (**first new alternative**) | tabulated | registry `geant4-icru90-stopping-11.4.2` afd1479497a605d9675c1fe2a3870ddb3f305110fb2d76e13bb600fa3ab281cd (single file, `member=None`), parser `icru90` (`data/icru90.py`), Geant4 Software License; log-log interpolation of S onto 200/decade from 1 MeV/u; range by `exact-loglog-quadrature-v1`; start range `icru90-integrated-from-lowest-node-v1` (the array's own nodes below 1 MeV/u, log-log quadrature, constant-S only below the lowest node) | proton (1–500 MeV), alpha (1–250 MeV/u); structural liquid water with requested I = 78 eV (otherwise fail closed); no deuteron (nuclear runs fail closed); no C/O | `consistency-only` (a V-C1 grading appends `topas-agreement-shared-lineage:V-C1:PASS` or `failed:V-C1`; never `validated`) |

**Licence-role extension for `icru90-water-v1` (written here).** Decision 0038 records the arrays as Geant4-licensed
source-embedded data ("permits use and modification with attribution"; lead's assessment, not legal advice) and uses
them as an evaluation reference. Construction use is the same *use* under the same licence: the arrays are read at run
time from the user's verified provenance cache, never committed, never redistributed; the ICRU report is cited, not
reproduced; the Geant4 notice is already shipped (`src/ionmc/THIRD_PARTY_NOTICES.md`). The registry `Dataset.role`
stays `evaluation` (its default role); the per-use role is carried by the model's `DataPin.use_role =
"construction-permitted"`. No operator step is needed (ESTIMATE; the S6 documentation step re-states the basis). No value derived node by
node from the arrays (resampled tables, CC-NODE/CC-INT point outputs) is committed; the repository holds counts,
coverage and aggregates only.

**Circularity rule.** A run whose construction data are also the evaluation reference of a row makes that row
circular for that run configuration (decision 0038 rule, generalised by `DataPin.dataset_id`). For `icru90-water-v1`:
the ICRU 90 part of row I2 (He; plan v3-008) and the ICRU 90 comparisons of decision 0038 (`compare_nist.py`, protons
and alphas) are circular for any configuration using it; every TOPAS He comparison above the BraggIon switch, and every
TOPAS proton comparison above 2 MeV (Q7), shares the construction array and tests transport, straggling and end of
range only. Row I1 of plan v3-008 is a backend-parity row and is not affected. The default's verdicts are untouched.

**Excluded construction sources.** G4EMLOW `ion_stopping_data` files: lineage-only (README: use within Geant4 only);
`resolve` cannot construct from a `lineage-only` or `evaluation-only` pin. libdEdx ICRU73_NEW (GPL-3; project licence
undeclared): `evaluation-only`; construction use is **not in scope** and requires an operator authorization first.
SRIM/MSTAR and z²-scaled PSTAR/ASTAR stay rejected (decisions 0038, 0043 §1).

### 2. Single selection point (PROPOSAL)

- `PhysicsOptions` gains `stopping_model: str | None`; exactly one of `stopping` (object, API-compatible) and
  `stopping_model` (id, resolved by `resolve`) is given. A `stopping` object is identified from its
  instance fields by `spec_of` (never from a class constant; scientific review I4): an exact match with a registered
  spec (class, every option, grid bounds and density, start-range rule, data pins) yields that `model_id`; any other
  source (for example `BetheStoppingSource(e_min_per_u=1.5)` or `barkas=False`, as existing tests use) is recorded as
  `model_id_effective = "unregistered"` with its `kind`, class and complete options in the identity block,
  `decision = "user"` and the evidence record `unvalidated`; a source without `spec_of` support fails closed.
  `stopping=BetheStoppingSource()` therefore resolves to `bethe-ionmc-v1` (row R1-ion) and a modified source is never
  reported under a registered id (row ID-USER). `default_source()` returns `resolve(DEFAULT_STOPPING_MODEL_ID)`.
- `validate` derives every stopping-dependent object from that one source: material tables, LET water row, deuteron
  tables, exact-Ω ion rows, residual range and step bounds; it asserts that all table identities name one `model_id`.
  `ExactOmegaMapping.build` and `fdm_domain_status` fail closed on a missing ion table (no Bethe default in run paths).
- **MCS proton path (decision on Q3): fixed and recorded.** The b-matching proton path stays `bethe-ionmc-v1`
  (1–4000 MeV) for every selected model: it is part of the MCS model (I11 rows), no table reaches the proton-equivalent
  energies it needs, and following the selection would change the MCS model. Ion runs record
  `mcs_proton_path_stopping = "bethe-ionmc-v1"` in the identity block. The S difference between the ion table and the
  proton path for a non-default model enters the mapping at ≤ 0.2 % above 10 MeV/u for He (ESTIMATE from the SEEN
  0038 comparisons) and is declared.
- Diagnostics (`parity.py`, `ion_mcs_checks.py`, `nuclear/build.py` D6, `a16_digest.py`) keep their behaviour; a test
  asserts that their source equals `resolve(DEFAULT_STOPPING_MODEL_ID)` and they record the model id where they record a
  source.

### 3. Fail-closed domains (PROPOSAL)

Before transport, `validate` checks species ∈ domain, the material rule for every geometry material and the LET water
row, and that every energy interval used lies inside the model's per-species window: table floor and cutoff,
`E0 + 6σ_E`, LET water row, deuteron tables of nuclear runs. A miss raises `UnsupportedCombinationError` naming the
model, species, material, energy and domain. There is no extrapolation, no fallback model and no automatic composite.
A missing cache object, a hash mismatch or a parse failure raises. `capabilities()` lists every registered model with
domain, data pins, licence role and evidence status.

### 4. Persisted identity (PROPOSAL)

`EffectiveConfig.summary()["stopping"]` (always present, outside the table hash):
`model_id_requested`, `model_id_effective`, `version`, `kind`, `spec_sha256`, `options` (complete, incl.
`ls_requested/ls_effective` and the I-value rule), `interpolation`, `range_construction`, `range_start`, `domain`,
`data` (each pin with `role_in_this_run = "construction"`), `content_sha256` of the parsed arrays, `evidence_at_run` (the evidence records of the effective spec at run time),
`decision`, `mcs_proton_path_stopping` (ion runs), `circular_rows` (dataset ids also used as evaluation references),
`selected_by` (`"stopping_model"` or `"stopping-object"`).
- Table metadata: models registered before this decision (`bethe-ionmc-v1` and its LS variants,
  `nist-pstar-astar-water-2005`) keep **byte-identical** metadata, so their table identities and hashes, the T1 fixture
  and the A16 digests are unchanged; their existing metadata already distinguishes them (`source`, `options.ls`,
  dataset sha256). Every model registered by this decision or later carries `model_id` and `version` in its table
  metadata, so its `TransportTables.sha256` differs from every existing model. A test asserts both directions and that the map (metadata, projectile) → effective model id is injective over the
  registry (for z = 1 the LS variants share the default's metadata by design and resolve to the default's effective id).
- Persistence on disk is the validation scripts' JSON step documents (Q5); a clean-process round trip of the block is a
  gating row of the plan.

### 5. No silent switch, no post-hoc choice (PROPOSAL)

- No fallback chain, no "best available" choice, no environment-variable override. `model_id_requested ==
  model_id_effective` except the documented z = 1 LS no-op.
- `DEFAULT_STOPPING_MODEL_ID` and the default's `spec_sha256` are pinned by a test and guarded like
  `A16_INTENDED_CHANGE`: a change requires a `STOPPING_DEFAULT_INTENDED_CHANGE` record that names a decision record;
  the record is `None` in V3-013.
- Within V3-013 no validation outcome changes the default or any registered spec. A FAIL is recorded by appending `failed:<row>` to the evidence registry (the spec, its hash and every table hash are
  unchanged); a revised alternative is a new `model_id`/`version` with fresh, held-out
  evidence (new seeds, new energies), pre-registered before it runs.
- No alternative motivated by a seen FAIL of the default (I3, TD-EM low and upper ends) is introduced without stating
  that motivation, consuming the affected rows and pre-registering fresh evidence. `icru90-water-v1` is not so
  motivated: it is the first tabulated construction path for p/α and covers no C/O energy of I3.

### 6. Documented default statement (to be carried into `docs/` in step S6)

"Default `bethe-ionmc-v1`: analytic Bethe stopping with Barkas-125 effective charge, Ashley–Ritchie–Brandt Barkas
term (Geant4 `G4EmCorrections` F table), Bloch, first-order Mott, Bichsel shell correction, Sternheimer density
effect, Lindhard–Sørensen off, I_water = 78 eV (ICRU 90). It is chosen for consistency with the qualified proton
engine and lineage alignment with the specified TOPAS configuration (shared I-value, ICRU 90 p/α agreement ≤ 0.2 %
above 10 MeV/u, shared Barkas table). It is **not TOPAS-equivalent**: TOPAS (opt4) uses the ICRU 90 alpha array for
He, the ICRU 90 proton array at equal velocity × q_eff² plus Lindhard–Sørensen and Barkas terms for C/O, and the ICRU
90 proton array for protons above 2 MeV; measured TOPAS/ionmc thin-slab differences are +3.5/+5.0 % (C/O) at
2–5 MeV/u, falling to +0.1/+0.3 % at 20–50 MeV/u (I3-FU-1 step 2, INCONCLUSIVE by rule); row I3 (C/O vs revised ICRU 73
table) FAILED and is retained as a declared low-energy C/O limitation." Not a clinical claim.

## Rejected alternatives

- **Switching the default to a tabulated or TOPAS-like model now** (`icru90-water-v1` for p/He, a Geant4-like C/O
  rebuild, an ICRU 73 low-energy composite): it would change the qualified proton model and the R2 digests, make the
  TOPAS rows test only mechanics, and — after the TD-EM low-end FAILs, I3 and I3-FU-1 — be a model chosen after
  seeing outcomes.
- **A string-only selection without registry or identity** (0038 style): no domain, no licence role, no pin.
- **Implicit composites or domain fallbacks** (table where available, Bethe elsewhere): a silent switch.
- **MCS proton path following the selection:** changes the MCS model; outside this task (I11 rows).
- **Adding `model_id` to the default's metadata:** breaks the asserted T1 identity and every default hash for no
  information gain (the summary block carries the id).
- **libdEdx ICRU73_NEW as construction:** GPL-3 in a project with undeclared licence → operator authorization first;
  also I3 circularity. **G4EMLOW files:** third-party terms, excluded.
- **A G4WaterStopping.cc low-energy C/O composite in this task:** not acquired, licence/provenance review open,
  counter-indicated by I3-FU-1 (the table moves ionmc past TOPAS below 20 MeV/u); left to a later task (plan §G).

## Uncertainty

- Default: unchanged (decision 0043 §1 statement and its 2026-11-16 note).
- `icru90-water-v1`: ICRU 90 quotes ≈ 1–2 % stopping uncertainty at therapeutic energies (TO VERIFY, decision 0043);
  log-log interpolation between the coarse nodes adds an interpolation error that is reported, not bounded by theory
  (plan row CC-INT); the start-range convention below 1 MeV/u is part of the identity.
- `nist-pstar-astar-water-2005`: I = 75 eV vs the material's requested I (≈ 0.4 % range at 200 MeV, decision 0038).

## Data roles and lineage

| Item | Role in V3-013 | Shared with TOPAS? |
|---|---|---|
| Bethe model (ionmc) | construction (default) | theory, I-value and Barkas table shared (decision 0043 §6) |
| ICRU 90 p/α arrays (`geant4-icru90-stopping-11.4.2`) | evaluation (default runs); **construction** in `icru90-water-v1` runs | identical to TOPAS He stopping above 7.9452 MeV and proton stopping above 2 MeV in G4_WATER (Q7) |
| NIST PSTAR/ASTAR | evaluation (default runs); construction in `nist-pstar-astar-water-2005` runs | not used by TOPAS |
| libdEdx ICRU73_NEW | evaluation-only (I3) | largely independent above ≈ 2 MeV/u |
| G4EMLOW ion stopping files | lineage-only | Geant4 internal |
| TOPAS 4.3 / Geant4 11.4.2 opt4 EM-only | evaluation (V-C1), shared-lineage MC agreement | — |

## Supported domains

As the registry table. Evidence for every alternative is configuration-specific: it applies to its own model id,
options and domain only and never transfers to the default; the default's evidence never transfers to an alternative.
Not a clinical claim.

## Validation strategy

Frozen in `validation/plans/v3-013-acceptance.md`: gating regression rows (R2 bit-identity, R1-ion, default identity),
consumer-consistency rows with a synthetic test-only model, domain and role fail-closed rows, identity persistence,
consistency checks of `icru90-water-v1` (not evidence), its configuration-specific validation V-C1 (He, TOPAS opt4
EM-only, fresh held-out energies 75/125/200 MeV/u and a fresh seed family; shared-lineage transport agreement, never
a validation of the stopping physics), the propagation checks V-C2 (stochastic, paired, with a seal-time informativeness
precondition) and V-C2-D (deterministic straight paths) and the report-only V-C3; the ionmc samples of both models are
produced and committed before any TOPAS run (plan §D).

## Expected tradeoffs and declared limitations

- The default stays not TOPAS-equivalent; C/O low-energy limitation and I3 FAIL retained.
- `icru90-water-v1` covers p/α in liquid water only, has no deuteron table (no nuclear runs) and no C/O; its He TOPAS
  agreement is shared-construction evidence; its proton use has consistency checks only (no proton row pre-registered).
- The MCS proton path stays Bethe for every model (declared mixed construction for non-default ion runs).
- No library-level result writer; identity persistence relies on `EffectiveConfig.summary()` and the JSON step writers.

## Pre-freeze review log (2026-11-17; findings recorded, resolutions applied in place before the freeze)

- **Codex REVIEW-0f80fad83827437582b3def7bf34b9b7 at 14099895 (changes_required, 2 important).** (1) The evidence
  status was part of the immutable spec and of `spec_sha256`, so a grading would have changed the model identity →
  evidence status removed from the spec; append-only evidence registry keyed by `spec_sha256` (section 1).
  (2) `validated:V-C1` overstated shared-lineage transport agreement as validation of the stopping model →
  record `topas-agreement-shared-lineage:V-C1:PASS` or `failed:V-C1`; never `validated` (section 1; plan A.6).
- **Independent scientific review (reviewer agent, read-only, at 14099895).** B1 (blocking): the plan ran TOPAS
  before the ionmc samples with only an ancestry check, so the ionmc model could have been adjusted after seeing
  TOPAS → ionmc samples of both models produced and committed before any TOPAS run, exact seal pins (plan §D).
  I1: CC-T/CC-RES/CC-I2/CC-NODE tolerances not achievable or non-informative as written → rows rewritten (plan
  A.2, A.5). I2: V-C2 could pass trivially where |ΔR_pred| is small → seal-time informativeness precondition and
  the deterministic row V-C2-D (plan A.6). I3: V-C1 bin sizes contradicted the §3 rule of r2-topas-target.md →
  per-case sizing table frozen (plan A.6). I4: a directly constructed source with modified fields would have been
  recorded under a registered id → `spec_of` from instance fields, `unregistered` otherwise (section 2; row
  ID-USER). I5: as Codex (2). I6: the licence-role extension silently amended decision 0038 → named in
  "Supersedes"; no node-derived values committed (section 1). M1: three line citations corrected. M2: injectivity
  stated over (metadata, projectile) → effective id. M3: the seed gap is transferred from V3-008 by a dated note
  in plan v3-008 §C. M4: V-C1-A results committed whatever the outcome. M5: the V2-OUTPUT ledger row reworded
  (identity block only; no writer). M6: the ledger's "pending the remaining gates" wording marked as superseded by
  the integration note.
- **Second pass of the scientific review (verification of 92cf8929):** all findings above and both Codex findings
  RESOLVED; new findings N1–N5 (plan section H) applied in the plan: the V-C2-D prediction compares telescoped end
  points Z_m = mean(z_end + R_m(E_end)) with a per-model round-trip gate (N1), V-C2-D uses seed r_index 15 (N2), the
  evaluator checks the pinned `src/ionmc` tree hash (N3), σ_exp is fixed by literals (N4), the V-C2-D substitute
  includes CC-ID/PER-1 for `icru90-water-v1` (N5).
- Confirmed by the scientific review: directive coverage (conditional on B1), the FACT claims of "Context"
  (apart from M1), the seed families (numerically), no table values, ledger and 0043 note accuracy, internal
  consistency, the R2/T1 dependence analysis.

## Preserved obligations (from the V3-008 integration note, verbatim)

"**Unresolved qualification obligations preserved:** row I3 (C/O stopping vs the revised ICRU 73 table) FAIL;
measured-data rows are evidence gaps (experimental validation deferred, directive 2026-11-07); carbon/oxygen nuclear
interactions, fragmentation, secondary ions and distal fragment dose are V3-009; the low ends and 430 MeV/u ends of the
domains were not re-measured after the R2-STOP3 model change; lineage identity of ICRU73_NEW with the 2009 revision
OPEN. Not a clinical claim."

Also preserved: the I3-FU-1 verdict INCONCLUSIVE (code_sha 48b41981), the declared low-energy C/O stopping limitation
(decision 0043 §1/§5 note of 2026-11-16), every V3-008 verdict as graded, the seed families of plan v3-008 §C.
