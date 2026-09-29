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
