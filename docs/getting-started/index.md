# Getting started

## Installation (development)

`ionmc` is an installable, `src/`-layout Python package (see
`pyproject.toml`) with `numpy` and `warp-lang` as runtime dependencies.
Development uses [`uv`](https://docs.astral.sh/uv/) to create the virtual
environment and install the package together with its `dev` extra
(`pytest`, `ruff`, `mypy`, `pre-commit`):

```bash
uv sync --extra dev
```

Run the test suite with:

```bash
uv run python -m pytest
```

Optionally install the pre-commit hooks (formatter, linter, type checker):

```bash
uv run pre-commit install
```

## Current functionality

The package is under active, autonomous development. After task `DEV-002`
it provides:

- unit conventions and constants (`ionmc.units`, `ionmc.constants`);
- materials with provenance-tagged mean excitation energies
  (`ionmc.materials`, e.g. `WATER`, `WATER_ICRU90`) and projectiles
  (`ionmc.particles`);
- the analytical electronic stopping power and CSDA range of protons and
  ions (`ionmc.AnalyticStoppingPower`) on three execution paths:

```python
import numpy as np
from ionmc import WATER, PROTON, AnalyticStoppingPower

model = AnalyticStoppingPower(WATER, PROTON, path="numpy")   # or "python", "warp"
energies = np.array([10.0, 100.0, 200.0])                   # MeV
model.mass_stopping_power(energies)                          # MeV cm^2/g
model.csda_range(energies)                                   # g/cm^2
model.provenance()                                           # I value, corrections, path
```

- a pure-Python mirror of Warp's random-number generator (`ionmc.rng`);
- a tabulated stopping-power layer fed by an external, versioned, cached
  dataset (`ionmc.TabulatedStoppingPower`, `ionmc.data`):

```python
from ionmc import WATER, PROTON, TabulatedStoppingPower
from ionmc.data import MCSQUARE_PSTAR_WATER, acquire

acquire(MCSQUARE_PSTAR_WATER)                    # once, with network access
model = TabulatedStoppingPower.from_dataset(     # offline afterwards
    MCSQUARE_PSTAR_WATER, WATER, PROTON, path="numpy"
)
model.mass_stopping_power([10.0, 100.0, 200.0])  # MeV cm^2/g
model.csda_range([100.0, 200.0])                 # g/cm^2 from the 0.5 MeV floor
```

The dataset can also be fetched from the command line:

```bash
python -m ionmc.data acquire --all      # download and cache the registered tables
python -m ionmc.data status             # show what is cached
```

Transport, geometry, scoring and tabulated physics data are planned
(`docs/development/roadmap.md`) but not yet implemented.
