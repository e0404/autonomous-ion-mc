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

- a pure-Python mirror of Warp's random-number generator (`ionmc.rng`).

Transport, geometry, scoring and tabulated physics data are planned
(`docs/development/roadmap.md`) but not yet implemented.
