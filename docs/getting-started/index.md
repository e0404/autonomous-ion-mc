# Getting started

`ionmc` is a Python package. At this stage it installs, reports its exact
source identity and exposes the command-line entry point; transport,
scoring and data acquisition are added by subsequent tasks and documented
here as they become usable.

## Installation from source

Requirements: Python 3.12 or newer and a C/C++-free environment (NVIDIA Warp
ships its own compilers). CUDA execution additionally needs an NVIDIA driver.

```bash
git clone https://github.com/e0404/autonomous-ion-mc.git
cd autonomous-ion-mc
python -m venv .venv && . .venv/bin/activate
pip install .
ionmc version
```

With [uv](https://docs.astral.sh/uv/): `uv sync` creates the locked
environment and `uv run ionmc version` runs the CLI.

The runtime dependencies are `numpy` and `warp-lang` only.

## Source identity

`ionmc version` prints JSON with the package version, the Git commit SHA and
whether the source tree was dirty. From a Git checkout the identity is read
live; an installed wheel reports the identity stamped at build time and says
so (`"source": "wheel:build-time-git"`). If neither is available the fields
are `null` and `"source": "unknown"` rather than a guessed value.

## External data cache

Large external datasets are never stored in Git. `ionmc.data` downloads them
on first use, verifies SHA-256, and reuses the local copy afterwards without
network access.

- Location: `$IONMC_DATA_DIR` if set, else `$XDG_CACHE_HOME/ionmc/data`, else
  `~/.cache/ionmc/data`.
- `IONMC_OFFLINE=1` forbids downloads; a dataset missing from the cache then
  raises `DatasetUnavailableError` naming the dataset, version and URL.
- Currently downloaded: NIST PSTAR/ASTAR stopping-power tables (NIST SRD 124,
  ICRU Report 49) for selected materials, fetched from
  `physics.nist.gov` as a plain-text table (about 9 kB each).
- Layout: exact bytes in `objects/<sha256>`; a JSON provenance record (URL,
  query fields, SHA-256, size, retrieval time UTC, license basis, citation)
  in `datasets/<dataset_id with / as __>/<version>.json`. A corrupted object
  is detected by re-hashing and fetched again.

## Stopping-power tables from the command line

```bash
ionmc data fetch                      # acquire PSTAR/ASTAR for materials with NIST codes
ionmc data fetch --material pmma --offline   # verify from the cache only
ionmc stopping --species c12 --material water --energy 100 290
```

`ionmc stopping` prints the electronic mass stopping power (MeV cm²/g) and
CSDA range (g/cm²) at the given kinetic energies per nucleon together with
the table's provenance (sources, dataset hashes, mean excitation energy,
blend window). Water needs no network: its low-energy data are the shipped
ICRU 90 tables.
