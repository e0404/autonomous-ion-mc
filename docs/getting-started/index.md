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

## Running a simulation

Write a JSON configuration (lengths mm, energies MeV/u):

```json
{
  "source": {"species": "proton", "energy_mev_per_u": 100, "position_mm": [0, 0, -1]},
  "geometry": {"type": "box", "size_mm": [40, 40, 100], "spacing_mm": 1, "material": "water"},
  "scoring": {"spacing_mm": [40, 40, 1]},
  "histories": 500, "batches": 5, "seed": 1, "backend": "python"
}
```

```bash
ionmc run config.json --output results/p100 --offline
```

This writes `results/p100.npz` (dose and energy per primary with standard
errors and voxel masses), `results/p100.json` (requested and effective
configuration, code identity, geometry and scoring coordinates, units, table
provenance, energy accounting, RNG identity) and `results/p100.txt` (a
human-readable summary with the integral depth-dose profile). The `python`
backend is a slow float64 reference (about 4×10⁴ steps/s); heterogeneous
phantoms use `"slabs": [[z0, z1, "bone_cortical", null]]` in the geometry.
`ionmc capabilities` prints which physics, scorers and backends are
implemented; unsupported requests fail before any transport.

Backends: `"backend": "python"` (float64 reference, slow), `"warp-cpu"`
(single host thread, ≈ 1000× faster than the reference for protons) and
`"warp-cuda"` (NVIDIA GPU). Warp backends accept `"precision": "float32"`
(default) or `"float64"`. The first run of a session compiles the kernels
(seconds); later runs reuse the kernel cache.
