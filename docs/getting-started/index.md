# Getting started

IonMC is a Python package (`ionmc`) built on NumPy and NVIDIA Warp. Requirements:
Python 3.12 or newer. The runtime dependencies are `numpy>=2` and
`warp-lang>=1.17,<2`; Warp runs on the CPU and, when an NVIDIA driver is
available, on CUDA devices.

At the current stage the package provides the project scaffold (version and
environment report), the external data layer, materials, the electronic
stopping-power model and tables, and the proton electromagnetic transport engine with
its Python reference backend (the Warp CPU and CUDA backends are not available yet).

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

## A minimal transport run

The reference backend (`backend="python"`, float64) runs a few histories of a 30 MeV proton
beam in a water box with the analytic Bethe stopping power (no downloaded data). It is slow by
design: use small numbers of histories.

```python
from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource

config = SimulationConfig(
    source=PencilBeamSource(
        PROTON, position_mm=(0.0, 0.0, 0.0), direction=(0.0, 0.0, 1.0), kinetic_energy_mev=30.0
    ),
    geometry=BoxPhantom(lower_mm=(-10.0, -10.0, 0.0), size_mm=(20.0, 20.0, 20.0), material=WATER),
    scoring=(ScoringGrid(origin_mm=(-10.0, -10.0, 0.0), spacing_mm=(2.0, 2.0, 1.0), shape=(10, 10, 20)),),
    physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource(), max_step_mm=1.0),
    run=RunOptions(backend="python", precision="float64", seed=1, n_histories=4, n_batches=2),
)
result = Simulation(config).run()
dose = result.grid("dose")
print(result.valid, dose.dose_gy.shape, result.energy_balance.relative_residual)
```

`result.grid("dose")` holds the mean energy and dose per primary (`energy_mev`, `dose_gy`), their
standard errors and the defined-value mask; `result.energy_balance` and `result.counters` report
where the energy went and whether any transport limit was hit (a nonzero counter raises
`TransportLimitError`). `nuclear` has no default and `nuclear=True` is rejected until nuclear
interactions are implemented. See [transport engine](../architecture/transport.md) and
[EM transport](../physics/em-transport.md).
