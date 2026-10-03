# External data layer

Package `ionmc.data`. External datasets are never committed to Git; they are
downloaded on demand into a user-controlled cache, verified against a pinned SHA-256
and used offline afterwards. Design decision:
[0038](../generated/decisions/0038-electronic-stopping-power-and-data-roles.md).

## Registry

`ionmc.data.registry.DATASETS` maps a dataset id to a frozen `Dataset(id, version, url,
method, post_body, sha256, bytes, license, citation, parser, description)`. Registered:

| id | Content | Licence |
|---|---|---|
| `nist-pstar-water-2005` | NIST PSTAR liquid water (POST to `physics.nist.gov/cgi-bin/Star/apdata.pl`) | NIST SRD 124, copyright all rights reserved: use-only, downloaded by each user, never redistributed |
| `nist-astar-water-2005` | NIST ASTAR liquid water (same endpoint, `prog=ASTAR`) | NIST SRD 124, copyright all rights reserved: use-only, downloaded by each user, never redistributed |
| `geant4-icru90-stopping-11.4.2` | `G4ICRU90StoppingData.cc` of Geant4 v11.4.2 (ICRU 90 proton/alpha arrays) | Geant4 Software License |

A downloaded payload whose SHA-256 or size differs from the registry is an integrity
failure (`IntegrityError`): nothing is stored.

## Cache

Directory precedence: explicit argument (`--cache-dir`), environment variable
`IONMC_CACHE_DIR`, `~/.cache/ionmc`. Layout:

    objects/<sha256>              raw bytes
    manifests/<dataset_id>.json   dataset_id, version, url, method, post_body, sha256,
                                  bytes, retrieved_at (UTC ISO-8601), license, citation

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

## Command line

    ionmc data list
    ionmc data fetch <id> [--cache-dir DIR] [--offline]
    ionmc data verify <id> [--cache-dir DIR]
    ionmc data path <id> [--cache-dir DIR]

`path` never uses the network. Exit code 1 reports a missing, offline-unavailable or
corrupt dataset.
