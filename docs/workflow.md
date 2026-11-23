# Workflow

From a clean checkout to a checked result bundle. Every step is a command of the [CLI](cli.md); the
configuration keys are in the [configuration reference](configuration/index.md) and the bundle in the
[results format](results/format.md).

## 1. Install

```bash
uv sync
uv run ionmc version
```

(or `pip install .` in a virtual environment, see [getting started](getting-started/index.md)).

## 2. Fetch the data a configuration needs

The analytic model `bethe-ionmc-v1` needs no download. The tabulated stopping models and nuclear runs need
cached datasets, which are downloaded once into `IONMC_CACHE_DIR` (default `~/.cache/ionmc`) and
re-verified against pinned hashes at every use:

```bash
ionmc data list
ionmc data fetch nist-pstar-water-2005
ionmc data verify nist-pstar-water-2005
```

## 3. Write a configuration

Start from a documented example in `tests/ionmc/data/v3_006/` (`proton_box.toml` is the smallest) and
follow the [configuration reference](configuration/index.md). `schema`, `schema_version`, `[source]`,
`[geometry]`, `[[scoring]]`, `[physics]` (with `nuclear` and `stopping_model`) and `[run]` are required;
every file reference needs its sha256.

## 4. Check

```bash
ionmc check proton_box.toml
```

Loads the file, verifies the datasets and the input hashes, validates and reconciles the configuration
and prints the effective summary. Nothing is run and nothing is written. A refusal is exit code 2 (the
file), 1 (a dataset), 4 (the combination is not supported) or 5 (no such backend); see the
[exit codes](cli.md#exit-codes).

## 5. Run

```bash
ionmc run proton_box.toml --out results/proton_box
```

The run writes `result.json`, `arrays.npz`, `report.md` and `plots/*.svg` to the destination. Use
`--overwrite` to replace an existing bundle; the old one is kept next to the new one as
`<path>.previous-<token>`.

## 6. Read the report and the plots

`report.md` states the identity of the code, the validity of the run, the requested and the effective
configuration, the physics and the stopping model with their evidence records as recorded, the data
identities, the random-number setup, the statistics definition, a per-grid summary, the energy balance and
counters and the list of missing items. `plots/*.svg` hold the depth dose profile with a one-standard-error
band, the depth-integrated energy profile and one lateral profile per grid, and a LET depth profile for
LET tallies. The report is regenerated without rerunning with `ionmc report results/proton_box`.

## 7. Verify

```bash
ionmc verify results/proton_box                 # tier 1, bundle only
ionmc report results/proton_box --check         # derived files equal their regeneration
ionmc verify results/proton_box --rederive --cache-dir "$IONMC_CACHE_DIR"
```

Tier 1 reads only the bundle and may be run anywhere. Tier 2 trusts your installation and cache and
compares the stored effective configuration with a rebuilt one; for a bundle with lookup tallies add
`--inputs lookup:<name>=PATH` to supply the trusted table files. The `binding_limits` printed by tier 1
list what tier 1 does not bind. For an untrusted bundle lower the reader budgets
(`--max-total-bytes` and the others).

## 8. Reproduce from `requested`

The `requested` block of `result.json` is the complete configuration. Writing it back as a TOML file
(with the same input files, which `input_manifest` identifies by sha256) and running it with the same seed, backend
and precision uses the same random streams; compare the two reports and `ionmc verify` output. The recorded `provenance.code` identifies the
code that produced the bundle.

## 9. Offline reuse

After the data have been fetched, `ionmc check`, `ionmc run`, `ionmc verify`, `ionmc report` and
`ionmc data fetch ID --offline` need no network. A bundle can be moved to another machine and read there
with tier 1 alone.

## What the capability contract refuses

`ionmc capabilities` declares what the engine supports; the declaration is derived from the rules the
validation enforces, so a configuration outside it is refused before transport with
`UnsupportedCombinationError` (exit 4) or `BackendUnavailableError` (exit 5), never run with a changed
backend, precision or model. Examples are nuclear interactions with an ion source, nuclear interactions on a
Warp backend, an energy outside the domain of the selected stopping model and a species the model does not
cover. A requested feature that would leave no trace in the result is an error, not a thinner bundle.

Not a clinical claim.
