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
(`ionmc.nuclear.build.build_nuclear_proton`, builder `ionmc-nuclear-proton-builder-2`, schema
`ionmc-nuclear-proton-table-2`) turns the pinned sources `endf-b8.0-protons`, `ame2020-mass` and
`nist-astar-water-2005` into `<cache>/derived/nuclear-proton-<id>.npz` and `.json`;
`id = sha256(source hashes, builder version, canonical options)`. The build is single-process and
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
  recorded, never enforced.
- **Loading.** `ionmc.nuclear.tables.NuclearTable.load(cache_dir, id)` uses
  `np.load(allow_pickle=False)`, re-hashes the npz against the JSON, re-derives the id, checks the
  source pins against the registry and freezes the arrays (`NuclearTableMissingError`,
  `NuclearTableStaleError`, `NuclearTablePinError`). `material_rows(material, f_e)` gives
  `Sigma_mass`, cumulative target fractions (and the unnormalised cumulative partial Sigma) and two majorant arrays: the window majorant
  `1.02 max Sigma` over `[E_k (1 - 2 f_E - 0.01), E_{k+1}]` and the end-of-range majorant
  (running maximum); use the step lookup.
- **Checks.** `validation/scripts/transport/nuclear_checks.py` runs N1, V1, V1b, V4, V4b and the D6 report.
