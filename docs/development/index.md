# Development

This section documents development tooling and contributor-facing workflows.

## Setup

The package uses a `src/` layout and the Hatchling build backend. Create the
development environment with [uv](https://docs.astral.sh/uv/):

```bash
uv sync --extra dev
```

The `dev` extra installs pytest, pytest-cov, ruff, mypy and pre-commit. The
resolved environment is locked in `uv.lock`; CI uses `uv sync --extra dev --locked`.

## Checks

```bash
uv run pytest                          # tests under tests/
uv run pytest -m "not cuda and not host"  # what CI runs
uv run mypy src/ionmc                  # type checking
uv run pre-commit run --all-files      # ruff, ruff-format, mypy, hygiene hooks
```

Pytest markers: `cuda` (needs a CUDA device), `slow`, and `host` (controlled
host runner only). Selection is left to the caller.

## Controlled host runner

The host runner executes tests from the source tree without installing the
package: `pythonpath = ["src"]` in `pyproject.toml` makes `import ionmc` work
directly from `src/`. Tests use pytest's `importlib` import mode so test
directories do not shadow the package.
