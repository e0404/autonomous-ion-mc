# Development

## Repository layout

| Path | Contents |
|---|---|
| `src/ionmc/` | The `ionmc` package (Python-first physics, Warp backends, scoring, results, CLI). |
| `tests/ionmc/` | Package tests. Tests marked `local` or `cuda` are excluded from GitHub CI and run under the local exact-SHA validation gate. |
| `tests/infrastructure/` | Tests of the experiment infrastructure. |
| `validation/` | Requirement ledger, release plan and validation methods/results. |
| `decisions/` | Consequential scientific and technical decision records. |
| `infrastructure/` | Autonomous-development tooling (MCP services, host runner, reference-engine execution, release qualification). |

## Environment and tests

```bash
uv sync --extra dev          # locked environment with pytest, ruff, mypy
uv run pytest -q             # all tests discoverable on this host
uv run pytest -q -m "not local and not cuda"   # the GitHub CI subset
pre-commit run --all-files   # formatting, linting, typing, hygiene
python3 infrastructure/docs/generate_docs.py && uvx zensical build --clean --strict
```

`pytest` imports the package from `src/` through the `pythonpath` setting in
`pyproject.toml`, so no installation is needed in CI, in the agent sandbox or
inside the controlled host runner, whose fixed virtual environment provides
only `numpy`, `warp-lang` and `pytest`.

## Execution environments

| Environment | Devices | Use |
|---|---|---|
| GitHub CI | CPU only | Lightweight software-quality checks (`ci.yml`). |
| Agent sandbox | Warp CPU | Development, reference-path and CPU tests. |
| Controlled host runner | Warp CPU + CUDA (RTX A6000) | Exact-SHA local scientific validation and performance measurement. |

## Requirement ledger

`validation/requirement-ledger.json` maps every MUST requirement to its
status, evidence suites and tasks; the rendered table is in the Validation
section.
