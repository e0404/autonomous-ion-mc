# Getting started

IonMC is a Python package (`ionmc`) built on NumPy and NVIDIA Warp. Requirements:
Python 3.12 or newer. The runtime dependencies are `numpy>=2` and
`warp-lang>=1.17,<2`; Warp runs on the CPU and, when an NVIDIA driver is
available, on CUDA devices.

At the current stage the package provides only the project scaffold: a version
and an environment report. No transport physics is implemented yet.

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
