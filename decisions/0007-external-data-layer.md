# 0007 — External physics-data layer: acquire/use split, versioning, integrity

- Status: accepted
- Date: 2026-09-10
- Task: DEV-003
- Affects: reproducibility, provenance, dependency footprint, maintainability, validation strategy

## Problem

`REQUIREMENTS.md` (*External Data Management*) requires that large physics
datasets are not committed to Git, are downloadable and cached automatically,
carry explicit immutable versions and integrity checks, keep enough provenance
for results to identify the exact data used, and work offline once fetched,
with a configurable cache location. The tabulated stopping-power layer
(decision `0008`) is the first consumer and needs this layer to exist first.

A hard constraint from the experiment environment: neither the orchestrator
sandbox nor the controlled host runner is guaranteed network access (the host
runner has none), while GitHub CI and the Codex worker do. The design must let
a dataset be acquired where a network exists and then used where none does.

## Context and evidence

- The MCsquare repository (Apache-2.0, UCLouvain) redistributes PSTAR-derived
  and Geant4-derived per-material stopping-power tables at immutable Git
  commits; the raw GitHub URL for a commit-pinned file is stable and
  content-addressable. Two independent fetches of the water file agreed
  (DEV-002). Sandbox egress to `*.githubusercontent.com` is pre-authorised in
  `.claude/settings.json`.
- NIST SRD 124 (PSTAR) itself is copyrighted and must not be vendored as a
  bulk dataset; only de-minimis factual values are quoted, cited to the URL
  (decision `0008`).
- The host runner mounts a dedicated cache at `/cache` with
  `XDG_CACHE_HOME=/cache`, so a cache resolved from `XDG_CACHE_HOME` is
  visible to host validation without any new mount.

## Candidate approaches

1. **Commit the tables into the repository.** Rejected for the general case
   by the *No large physics datasets in Git* MUST; a small committed subset is
   kept only as regression reference (`ionmc.reference_data.pstar_water`,
   20 numbers), not as the data layer.
2. **Download on first use, transparently, inside the physics code.** Rejected:
   it couples transport code to network availability and would fail on the
   host runner exactly when validation needs it; it also hides when and from
   where data were obtained.
3. **Explicit acquire/use split with a content-addressed cache (selected).**

## Selected approach

Two phases, implemented in `ionmc.data`:

- **acquire** (`ionmc.data.acquire`, CLI `python -m ionmc.data acquire`):
  needs network. Downloads a `DatasetSpec` (name, immutable `version`,
  commit-pinned `url`, `sha256`, `filename`, `license`, `provenance`) to a
  temporary file in the target directory, verifies the SHA-256 against the
  spec, and only then atomically renames it into place and writes a JSON
  manifest (spec, checksum, size, retrieval time). A checksum mismatch writes
  nothing and raises `IntegrityError`.
- **use** (`ionmc.data.load_path`): offline. Returns the cached file for a
  spec or raises `DatasetNotCached` with the exact acquire command; verifies
  the cached bytes against the spec checksum on load (tamper/corruption
  detection), never touches the network.

Cache location precedence (`ionmc.data.cache_dir`): explicit argument,
`IONMC_CACHE_DIR`, `$XDG_CACHE_HOME/ionmc`, `~/.cache/ionmc`. Datasets are
stored under `datasets/<name>/<version>/`, so different versions coexist.
The registry (`ionmc.data.registry`) pins each dataset to a Git commit and
its file SHA-256; only redistributable data (Apache-2.0 MCsquare tables) are
registered.

## Rationale

Content-addressing by SHA-256 makes acquisition reproducible and
tamper-evident and makes "which data were used" a checksum in the provenance
record rather than a path. The acquire/use split is the only structure that
satisfies both the offline-reuse MUST and the no-network reality of the host
runner: CI or the Codex worker (or a networked orchestrator session, as here)
acquire into a cache the host runner can read. Atomic rename plus
verify-before-place means an interrupted or corrupted download never yields a
half-written "cached" file.

## Expected tradeoffs

- A dataset must be acquired before host validation can use it; the validation
  script fails with exit code 4 and an actionable message if it is missing,
  rather than attempting a download in the no-network sandbox.
- The registry must be updated (new version + checksum) when an upstream file
  changes; this is deliberate, since a silent upstream change must not alter
  results.
- Only redistributable data are registered; NIST-hosted originals are not
  vendored (licensing), so the project depends on MCsquare's redistribution
  for automatic acquisition. The committed 20-point subset is the fallback
  regression reference required by the roadmap's network contingency.

## Validation strategy

`tests/ionmc/test_data_layer.py`: acquire/verify/offline-load round trip with
a `file://` source, checksum-mismatch rejection writing nothing,
tamper detection on load, side-by-side versions, cache-dir precedence, CLI.
A `@pytest.mark.network` test acquires the real registered table when a
network is present and is skipped otherwise. The datasets were acquired into
the host-runner cache (`/cache/ionmc`) for the DEV-003 host validation.

## Later validation outcome

DEV-003: both registered water tables (`mcsquare-pstar-water`,
`mcsquare-g4-water`) were acquired over the network in the orchestrator
session and verified by SHA-256; the offline `use` phase then served them to
the tabulated layer and to the host runner without network access. See the
DEV-003 validation record.
