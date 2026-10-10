# External data layer

Package `ionmc.data`. External datasets are never committed to Git; they are
downloaded on demand into a user-controlled cache, verified against a pinned SHA-256
and used offline afterwards. Design decision:
[0038](../generated/decisions/0038-electronic-stopping-power-and-data-roles.md).

## Registry

`ionmc.data.registry.DATASETS` maps a dataset id to a frozen `Dataset(id, version, url,
method, post_body, sha256, bytes, license, citation, parser, description, role)`. `role` is the
evidence role (`construction`, `calibration`, `evaluation` or `exploratory`, default
`exploratory`). Parser names are identifiers; the registry imports no parser. Registered:

| id | Content | Licence |
|---|---|---|
| `nist-pstar-water-2005` | NIST PSTAR liquid water (POST to `physics.nist.gov/cgi-bin/Star/apdata.pl`) | NIST SRD 124, copyright all rights reserved: use-only, downloaded by each user, never redistributed |
| `nist-astar-water-2005` | NIST ASTAR liquid water (same endpoint, `prog=ASTAR`) | NIST SRD 124, copyright all rights reserved: use-only, downloaded by each user, never redistributed |
| `geant4-icru90-stopping-11.4.2` | `G4ICRU90StoppingData.cc` of Geant4 v11.4.2 (ICRU 90 proton/alpha arrays) | Geant4 Software License |
| `endf-b8.0-protons` | ENDF/B-VIII.0 proton sublibrary zip (LA150, MF3/MF6; construction) | no licence text, free NNDC download |
| `ame2020-mass` | AME2020 `mass_1.mas20.txt` (construction) | no licence text, free AMDC download |
| `exfor-d0356` | EXFOR entry D0356, Auce et al. 2005 (evaluation) | CC BY 4.0 |
| `exfor-c1862` | EXFOR entry C1862, Slaus et al. 1975 (exploratory, report-only) | CC BY 4.0 |
| `geant-val-exfor-inelastic-7` | geant-val JSON of EXFOR-derived inelastic curves (exploratory, report-only; import only) | none stated |

The nuclear entries are the data of decision 0041. Git holds only these registry entries and
`src/ionmc/data/exfor_manifest.json` (EXFOR entry and subentry identifiers, REACTION strings,
roles, hashes; no cross-section values).

A downloaded payload whose SHA-256 or size differs from the registry is an integrity
failure (`IntegrityError`): nothing is stored.

## Cache

Directory precedence: explicit argument (`--cache-dir`), environment variable
`IONMC_CACHE_DIR`, `~/.cache/ionmc`. Layout:

    objects/<sha256>              raw bytes
    manifests/<dataset_id>.json   dataset_id, version, url, method, post_body, sha256,
                                  bytes, retrieved_at (UTC ISO-8601), license, citation
                                  (an import adds source_path and has method "import")

`ionmc.data.cache.verify` re-hashes the object and compares it with both the manifest and
the registry.

## Acquisition

`ionmc.data.acquire.fetch(dataset_id, cache_dir=None, offline=False)` returns the path of
the verified object. Cached datasets are re-hashed and returned without network access.
Otherwise the dataset is downloaded with `urllib` (https only, redirects only to https,
30 s timeout), verified, stored and a manifest written. With `offline=True` or
`IONMC_OFFLINE=1` a missing dataset raises `OfflineError`. The controlled host runner has
no network, so datasets must be fetched before such runs.

Downloads are streamed in chunks and aborted with `IntegrityError` as soon as more than
the registered `bytes` arrive; a shorter payload fails the hash check.

## Offline import

`ionmc.data.acquire.import_file(path, dataset_id, cache_dir=None)` stores a local file as a
registered dataset without any network access. The size and SHA-256 are checked against the
registry first; a mismatch raises `IntegrityError` and nothing is written. The manifest has the
shape of a download manifest with `method` `"import"` and the `source_path`.

## Parsers

`load_star_table` and `load_icru90_water` re-hash the file they parse and require the hash
to equal the pinned SHA-256 of a registered dataset (`IntegrityError` otherwise). The result
carries `dataset_id`, `version`, `sha256` (pinned), `content_sha256` (of the parsed bytes)
and `retrieved_at` (from the manifest); these are copied into the metadata of tables
built from them. `allow_unverified=True` parses any file and records only `content_sha256`
(`dataset_id` is None). `parse_*_text` functions parse text without any identity check.

- `ionmc.data.nist_star.parse_star_text` returns a `StarTable` (energy [MeV], electronic,
  nuclear, total stopping power [MeV cm^2/g], CSDA and projected range [g/cm^2],
  detour factor); energies must increase and each row has seven columns.
- `ionmc.data.icru90.parse_icru90_source` extracts the liquid-water (index 1 of
  `{G4_AIR, G4_WATER, G4_GRAPHITE}`) proton and alpha electronic stopping arrays
  [MeV cm^2/g] and their energy grids [MeV] from the C++ source.

## Nuclear-data parsers (decision 0041)

The datasets of the proton non-elastic model are read by pure-Python parsers that carry no
data values; the real files stay in the cache. Fixtures in the tests are hand-written synthetic
text. A dataset that is not cached skips the data-backed tests; `IONMC_REQUIRE_DATA=1` turns
that into a failure (the nuclear-data analogue of `IONMC_REQUIRE_NIST`).

- `ionmc.data.endf6` (basis: ENDF-6 Formats Manual, BNL-203218-2018-INRE, chapters 0-3 and 6):
  `parse_endf(text)` returns an `EndfMaterial` (`za`, `awr`, `emax` from MF1/MT451, raw
  `sections[(MF, MT)]`); `cross_section(mt)` reads MF3 into a `Tab1` (E in eV, sigma in b) whose
  `interpolate` honours NR > 1 regions and INT laws 1-5; `products(mt)` reads MF6 into
  `Product` records (ZAP, AWP, LIP, LAW, yield, LCT) with LAW=1 distributions (LANG=1 Legendre,
  LANG=2 Kalbach-Mann with NA=1: `b_0 = f_0`, `b_1 = r`; LEP 1 and 2). LAW=5 is stored opaque;
  any other LAW or LANG raises `UnsupportedEndfError` (fail closed). `list_zip_members` and
  `read_member` read the sublibrary zip with name validation (no path traversal).
- `ionmc.data.ame`: `load_ame2020(text)` returns `{(Z, A): AmeEntry}` (mass excess [keV], atomic
  mass [u], uncertainties, `estimated` for `#` values); `nuclear_mass_mev` is the atomic mass
  times u minus Z electron masses (CODATA 2018), electron binding energy ignored (about 2 keV
  for oxygen).
- `ionmc.data.exfor`: `parse_entry(text)` returns an `ExforEntry` of subentries with BIB keys,
  `REACTION` pointers, COMMON and DATA blocks (11-character fields, `None` for blanks);
  `energy_to_mev` flags `MEV/A`, `energy_total_mev` multiplies by the projectile mass number,
  `xs_to_mb` and `percent_to_fraction` convert units.
- `ionmc.data.geant_val`: `parse_geant_val(text)` returns one `GeantValCurve` per record (target,
  beam, observable, x, y, statistical and systematic y errors). Non-proton projectiles carry an
  `energy_caveat` because the file does not state whether their energy axis is total or per
  nucleon.

## Command line

    ionmc data list
    ionmc data fetch <id> [--cache-dir DIR] [--offline]
    ionmc data verify <id> [--cache-dir DIR]
    ionmc data path <id> [--cache-dir DIR]
    ionmc data import <file> --dataset <id> [--cache-dir DIR]

`path` never uses the network. Exit code 1 reports a missing, offline-unavailable or
corrupt dataset.

## Derived nuclear table (decision 0041)

`ionmc data build nuclear-proton [--cache-dir D] [--points-per-decade 50] [--diagnostic-events N]
[--strict]`
(`ionmc.nuclear.build.build_nuclear_proton`, builder `ionmc-nuclear-proton-builder-5`, schema
`ionmc-nuclear-proton-table-3`) turns the pinned sources `endf-b8.0-protons`, `ame2020-mass` and
`nist-astar-water-2005` into `<cache>/derived/nuclear-proton-<id>.npz` and `.json`;
`id = sha256(canonical JSON of the complete sidecar without its `table_id` field)` (sorted keys,
separators `,` `:`); the sidecar holds `npz_sha256`, so the id binds the npz bytes (arrays, the
`qualification` flags and the capacity `bounds`) and every sidecar field (sources, options, builder
and schema, element scales, targets, D6 numbers). The builder writes the npz, hashes it, assembles
the sidecar, computes the id and renames the files to it. The build is single-process and
deterministic (exact lambda enumeration, seeded counter-based uniforms for the diagnostics,
fixed-timestamp npz) and writes identical bytes twice. The total break-up residual
(Z_r, A_r) = (0, 0) of mass 0 is always allowed (`empty_residual_allowed: true`); the first yield
and first distribution row of every product are extended downward (constant) from the first MF6
energy to the MT5 threshold (`targets[].yield_extended_below_mev`). The build fails closed on
exact P_accept < 0.5 and on sigma > 0 without yields.

- **Cross sections.** ENDF MF3/MT5 on its native interpolation up to 150 MeV (0 below threshold);
  150-250 MeV: `sigma(150) sigma_TL(E)/sigma_TL(150)` (Tripathi light system). Union grid: 1 MeV,
  all ENDF MF3/MT5 nodes, MF6/MT5 incident energies and yield nodes of all targets in [1, 150] MeV
  and a uniform ln E grid (>= 50 points/decade, 150 MeV a node) up to the first node >= 250,
  interpolated lin-lin in E (exact for the INT=2 ENDF data; the build fails on any other law);
  lookup by the shared fixed-step bisection `grid_locate`.
  Surrogates (Na, Mg -> Al-27; S, Cl -> P-31; K, Ar -> Ca-40) are scaled by `(A/A_ref)^(2/3)`
  (recorded per element); hydrogen is 0; any other element is absent
  (`UnsupportedCombinationError` from `NuclearTable.material_rows`).
- **Product rows.** For n, p, d, alpha, gamma: ENDF yield, 64 equiprobable E'_CM bins (quantile
  interpolation between incident energies) and the Kalbach `r` per bin; above 150 MeV the 150 MeV
  rows with E' scaled by `E_avail(E)/E_avail(150)`. Kalbach separation energies use the systematics
  formula (`"kalbach_separation": "systematics-formula"`); AME2020 masses only for the residual-mass
  test and Q values. The residual gets the ENDF mean heavy-recoil energy `recoil_t_cm_mev` (formally CM under LCT=3; it
  is deposited locally without a boost, a few-MeV quantity).
- **Multiplicities.** Floor + Bernoulli per species (`multiplicity_round`), residual-existence
  acceptance (at most 64 attempts), no E*; `lam_s(E)` solved at every grid node by a fixed point on
  the exact enumeration of the 2^5 outcomes so that post-acceptance mean yields equal ENDF
  (`p_accept`, `yield_ratio_post`, `lam_converged` arrays; non-converged nodes and `p_accept_min`
  in the JSON). The per-event ledger `imbalance = T1 + m_p + M_t - sum E_lab - M_r - T_r` is the
  signed tally `nuclear_imbalance`; `diagnostics` holds 2e4 seeded events at 18 nodes per target
  (Delta_lab, Delta_CM, P(Delta<0), mean |sum p_CM|, E' ratios, local deposit);
  `transport_path_bound_terms` (`n_max`, `t_lab_max_mev` for p and d) give the capacity bound
  `B_L = 1.25 (mixed_path_bound(E_hi) + sum_s n_max,s mixed_path_bound_s(t_lab_max,s))`;
  `transport_energy_bound_mev` is the largest per-particle `t_lab_max_mev` (<= 500 MeV).
- **D6.** `gate_d6.numbers["150"|"250"]` hold `G`, `D`, `p_event`, the 99.9th-percentile lab alpha
  energy and its ASTAR CSDA range; `tier1_pass`, `tier2_pass`, `ceiling_pass` (D <= 2e-2, 99.9th-percentile range <= 3 g/cm2) are
  recorded in the JSON, but the gate is the `qualification` array of the npz (ceiling flags,
  all nodes converged, tiers), authenticated by the id: the loader derives it from the arrays,
  cross-checks the JSON copies (`NuclearTableStaleError` on any disagreement) and refuses a table
  whose D6 ceiling or convergence gate failed (`NuclearTableUnqualifiedError`; the tiers are
  reported, only the ceiling and the convergence are enforced). The `bounds` array of the npz holds
  the capacity bounds the transport uses; the JSON bounds are a view of it.
- **Loading.** `ionmc.nuclear.tables.NuclearTable.load(cache_dir, id)` uses
  `np.load(allow_pickle=False)`, re-hashes the npz against the JSON `npz_sha256`, recomputes the id from
  the canonical sidecar without its id and requires it to equal the id it was asked for (the
  config pin), checks the source pins against the registry, derives qualification and bounds from
  the npz arrays and freezes them (`NuclearTableMissingError`, `NuclearTableStaleError`,
  `NuclearTablePinError`, `NuclearTableUnqualifiedError`); editing any sidecar field keeping
  the id is stale, resealing the id gives another id that the pin does not name. `material_rows(material, f_e)` gives
  `Sigma_mass`, cumulative target fractions (and the unnormalised cumulative partial Sigma) and two majorant arrays: the window majorant
  `1.02 max Sigma` over `[E_k (1 - 2 f_E - 0.01), E_{k+1}]` and the end-of-range majorant
  (running maximum); use the step lookup.
- **Checks.** `validation/scripts/transport/nuclear_checks.py` runs N1, V1, V1b, V4, V4b and the D6 report.

## Nuclear pipeline end to end (decision 0041)

1. **Sources.** `ionmc data fetch|import` places the hash-pinned datasets in the cache: ENDF/B-VIII.0
   proton sublibrary (`endf-b8.0-protons`, MF3/MT5 and MF6/MT5 per target), AME2020 masses
   (`ame2020-mass`), NIST ASTAR (`nist-astar-water-2005`, alpha ranges for D6) and the EXFOR and
   geant-val entries of the V1b manifest. Roles are assigned per decision 0038; a dataset that is
   missing or whose hash differs stops every later stage.
2. **Parsers.** `ionmc.data.endf6`, `exfor`, `ame` and `geant_val` read the sources into immutable
   records (ENDF TAB1/TAB2/LIST with INT 1-5, LAW=1 LANG=1/2, LCT=3; EXFOR units; AME estimated values).
3. **Event model.** `ionmc.nuclear.events` builds per-target models (AME masses) and the shared
   functions (floor+Bernoulli multiplicities, inverse-CDF E', Kalbach mu, residual acceptance,
   ledger) used by the builder, the python reference and the Warp twins.
4. **Build.** `ionmc data build nuclear-proton` writes `derived/nuclear-proton-<id>.npz` and `.json`
   (id = sha256 of the canonical sidecar without its id, which contains the npz sha256) with sigma, product rows, lambda,
   `p_accept`, the D6 gate, the capacity-bound terms and the diagnostics block.
5. **Load.** `NuclearTable.load` re-hashes the arrays, recomputes the id from the sidecar, checks the registry pins
   and freezes the arrays; `material_rows(material, f_e)` composes `Sigma_mass`, the target
   fractions and the two majorants for a material of the geometry.
6. **Configuration.** `SimulationConfig` with `nuclear=True` and `nuclear_table_id` is validated
   (python backend, proton source, supported elements, E0 + 6 sigma_E <= 250 MeV, table present and
   current); the effective configuration records the table id, the file hashes and the row grids.
7. **Transport.** The python reference runs the thinning, event and stack logic of
   `docs/architecture/transport.md` ("Nuclear interactions") and returns the nuclear tallies and
   counters.
8. **Checks.** `validation/scripts/transport/nuclear_checks.py` (N1, V1, V1b, V4, V4b, D6) and the
   `lv5` suite read the same table; every `lv5` document records the table id and the sha256 of
   the `.json` and `.npz` files.

## Derived elastic table (decision 0041 slice C, V3-005C step C2)

`ionmc data build elastic-proton` writes `derived/elastic-proton-<id>.npz` and `.json` from the
hash-pinned LA150 proton sublibrary, AME2020 and the Geant4 11.4.2 Barashenkov/BGG sources (registry
ids `geant4-*-11.4.2`; the PDG p-p compilation `pdg-rpp2022-pp-elastic` has role "evaluation" and is not a construction input). Git holds only the registry entries, the EXFOR
manifest and `src/ionmc/data/elastic_table_pin.json` (id and hashes); the table is not in Git. This
step builds data only: no transport kernel uses the table yet (C3 to C5).

- **Targets.** `H-1` followed by the seven LA150 targets of the nuclear table (same
  `ELEMENT_TARGET` map and `(A_el/A_ref)^(2/3)` surrogate scaling; hydrogen has its own target).
  Grid: the nuclear union-grid convention (1 MeV, the H-1 ENDF nodes, uniform ln E with at least
  50 points per decade) plus every Barashenkov node, cut at 250 MeV, with midpoint refinement below
  14 MeV (BGG Coulomb rule) and for H-1 above 13 MeV so that lin-lin interpolation of `sigma` in E
  reproduces the direct evaluation to 5e-4 (`nodes_added_*` in the JSON).
- **p + A.** `sigma_el` is the BGG/Barashenkov value (`ionmc.nuclear.bgg`: arrays hand-transcribed
  from `G4BarashenkovData.hh` and compared with the pinned file by the builder, the exact Geant4
  rules including the Coulomb factor below 14 MeV). The angular shape is the black-disk form
  `|2 J1(qR)/(qR)|^2`, `q = 2 p_CM sin(theta_CM/2)`, with `pi (R + lambdabar)^2 = sigma_nonel` of
  the transport's own LA150 MT5 (Tripathi shape above 150 MeV); R, lambdabar, `sigma_nonel` (and
  its sha256) and the inversion residual are stored per node. The inverse CDF of `mu_CM` over
  [-1, 1] is stored as 257 edges of 256 equiprobable bins per node. Nodes without a real radius
  (`R <= 0` or `sigma_nonel = 0`) set `sigma_el = 0` below `e_min_shape` (per target, recorded).
- **H-1.** The LA150 Hale LAW=5 LTP=1 reconstruction (`ionmc.nuclear.law5`, ENDF-102 eqs 6.9, 6.10,
  6.14, b/sr CM) with ratio interpolation `sigma_e/sigma_c` in E; the transported density is
  `NI = sigma_e - sigma_c` over `|mu_CM| <= 0.96` (theta_CM >= 16.26 deg). `sigma` is the
  half-sphere integral (`2 pi`, each event of identical protons once); the edges cover the
  symmetric range [-0.96, 0.96]. `e_min_pp` (12.53 MeV) is the lowest node from which every node and
  midpoint has a non-negative density, `sigma = 0` below it. The frozen row P6 (negativity) FAILS
  for H-1 below `e_min_pp`; this is kept as a recorded failure and not converted into a pass:
  `negative_density_below_e_min_pp: true` with every negative node/midpoint and its value
  (`first_negative_nodes`), while `e_min_pp` and `no_negative_density_in_domain: true` record the
  separate revised-domain check (the qualification gate). Above 150 MeV `sigma = sigma_NI(150)
  S(E)/S(150)` with the Geant4 `G4HadronNucleonXsc` p-p formula only (no PDG data enter the table
  or its identity) and the shape fixed.
- **Domain (declared limitation).** The JSON `elastic_domain` block gives `[e_min, 250]` MeV per target (H-1: `E_min,pp`;
  p+A: `e_min_shape`) and `sigma_bgg_below_domain_max_mb` per target; `ElasticTable.elastic_domain()`
  returns a dict, the loader asserts `sigma = 0` and `valid = 0` below it, the npz carries
  `target_e_min_mev` and the packed `ElasticDevice` the fields `e_min_shape` and `e_min_pp`, so that
  the transport can count `elastic_below_domain` and `pp_below_domain` crossings. The builder also
  fails if `E_min,pp` exceeds 15 MeV. P6-D is evaluated in a pass separate from the one that fixes
  `E_min,pp`. The builder fails (`BuildError`) if any
  `e_min_shape` exceeds 10 MeV. The JSON records the NI cross section above the cut at 15 and 20
  MeV (`sigma_ni_above_cut_mb`), the p+O estimate below `e_min_shape(O-16)` with the BGG `sigma_el`
  (`pa_omitted_events_per_history_150mev`) and the omitted p-p correction below `E_min,pp`
  as a total-variation diagnostic: `pp_omitted_ni_correction_total_variation_per_history_150mev`
  and `pp_omitted_ni_correction_weighted_total_variation_mev`, with the nodes, M(E), W(E), the
  residual range and a `semantics` string in the block `pp_omitted_ni_correction`. The quantity is
  the total variation of the signed nuclear-plus-interference correction to Rutherford scattering
  over the half sphere above the cut (`|rho_NI|` integrated over `0 <= mu_CM <= 0.96`), integrated
  along the residual path in water from `E_min,pp` to 1 MeV; it is an event-equivalent magnitude,
  not an expected number of physical events, and the weighted form (recoil-energy weighted) is not
  energy transferred. Limitations: the signed correction can cancel, and the physical sigma_pp below
  `E_min,pp` is not available from the evaluation. (The S-wave unitarity bound of the C2b build was
  not a bound on the full amplitude and was withdrawn.) The P6 failure record, P6-D and the
  MF6 per-product gating of the O-16 check are described under the findings below.
  `model_revisions` references Amendment 15.
- **Findings and report-only data.** The builder asserts, fail closed, that the LA150 O-16 MT2 is a
  numerical copy of C-12 from 24 MeV (recorded), that the O-16 MF3/MT5 sigma_nonel (the input of the
  disk radius) is not a copy of C-12's (`o16_mt5_copy_check`: relative difference above 1e-3
  at more than 50 % of the shared nodes and a coefficient of variation of the O/C ratio above 1e-3;
  flag `o16_mt5_sigma_not_c12_copy`, recorded as `o16_mt5_finding`). The MF6/MT5 per-product yield rule
  as frozen (all compared energies in the denominator) fails on LA150 (zap 3007 14/30 and zap 5012
  10/30, both-zero nodes 16 and 20) and is kept as `mf6_rule_as_frozen` with `passes: false`; the
  informative-node revision (`mf6_rule_revised`) was written after that outcome and is report-only
  and influences no flag (Amendment 16 item 1). Independence of the O-16 product data is instead
  checked on the energy spectra (`o16_mt5_spectra_check`, row X-MT5-SPEC): for n, p and alpha at
  50, 100 and 150 MeV the normalised angle-integrated MF6/MT5 secondary-energy distributions
  (LAW=1, LEP=1, read by the same `SpeciesTables` reader as the nuclear table builder; exact
  histogram cumulative at a tabulated incident energy, the builder's quantile interpolation
  otherwise) of O-16 and C-12 are compared by the Kolmogorov-Smirnov distance D on the union E'
  grid; `o16_mt5_spectra_independent` is true iff D > 0.02 at at least 8 of the 9 combinations and
  D > 0 at all 9 (D below 1e-12 counts as 0, a floating-point floor), recorded in
  `o16_mt5_spectra_finding`, and the build fails closed otherwise. The qualification flags are
  `no_negative_density_in_domain`, `o16_mt5_sigma_not_c12_copy`, `o16_mt5_spectra_independent`
  and `shape_normalised`. The builder never uses LAW=5 data other than H-1 as a
  construction input. The LA150 C-12/N-14/Ca-40 LTP=12 NI densities (ratio interpolation in mu) at
  20-40 deg and 50/100/150 MeV are stored next to the model values (row X-ENDF, report-only).
  `sigma_el` at 10 MeV per target is recorded for the comparison with the decision-0041 planning
  value.
- **Load and device.** `ElasticTable.load` re-hashes the npz, recomputes the id, checks the registry
  pins and the qualification flags and refuses missing, stale, mis-pinned and unqualified tables
  (all `UnsupportedCombinationError`); `material_rows` gives `Sigma_el` and the majorants;
  `ElasticDevice` / `cached_elastic_device` (`ionmc.transport.elastic_device`) pack the arrays per
  process and key as `NuclearDevice` does. Exact two-body kinematics: `ionmc.nuclear.elastic_kin`.
