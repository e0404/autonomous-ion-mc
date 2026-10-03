# Getting started

IonMC is a Python package (`ionmc`) built on NumPy and NVIDIA Warp. Requirements:
Python 3.12 or newer. The runtime dependencies are `numpy>=2` and
`warp-lang>=1.17,<2`; Warp runs on the CPU and, when an NVIDIA driver is
available, on CUDA devices.

At the current stage the package provides the project scaffold (version and
environment report), the external data layer, materials and the electronic
stopping-power model and tables. Transport is not implemented yet.

## Installation

From a clone of the repository, with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
uv run ionmc version
```

or with pip in a virtual environment:

```bash
python -m venv .venv
. .venv/bin/activate
pip install .
```

## Command line

```bash
ionmc version   # print the package version
ionmc --version # same, in "ionmc <version>" form
ionmc info      # print versions and Warp devices as JSON
```

`ionmc info` reports the ionmc, Python, NumPy and Warp versions and the Warp
devices (alias, whether CUDA, name, architecture). Git SHA and dirty flag are
reported only if a tracked `src/ionmc/_build_info.json` exists; otherwise they
are `null`. Warp may print CUDA driver messages on stderr on hosts without a GPU.

## External data

Reference tables (NIST PSTAR/ASTAR and the ICRU 90 arrays embedded in Geant4) are not
part of the repository. Download them once into the cache (`IONMC_CACHE_DIR` or
`~/.cache/ionmc`); later use is offline and re-verified against pinned SHA-256 hashes:

```bash
ionmc data list
ionmc data fetch nist-pstar-water-2005
ionmc data fetch nist-astar-water-2005
ionmc data fetch geant4-icru90-stopping-11.4.2
ionmc data verify nist-pstar-water-2005
ionmc data fetch nist-pstar-water-2005 --offline   # cache only, never the network
```

See [data layer](../architecture/data-layer.md) and
[stopping power](../physics/stopping-power.md).
