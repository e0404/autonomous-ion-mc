# Command line

`ionmc` is installed with the package. Global option: `--version` (prints `ionmc <version>`).

| Command | Purpose |
|---|---|
| `ionmc version` | print the package version |
| `ionmc info` | print versions and Warp devices as JSON |
| `ionmc notices` | print the packaged third-party notices |
| `ionmc data ...` | manage cached external datasets |
| `ionmc run CONFIG.toml` | run a configuration and write a result bundle |
| `ionmc check CONFIG.toml` | load, validate and reconcile a configuration; no transport |
| `ionmc capabilities` | print the capability report as JSON |
| `ionmc verify DIR` | verify a result bundle (tier 1; `--rederive` adds tier 2) |
| `ionmc report DIR` | print the report regenerated from a bundle |

Warp prints a banner and kernel messages on stdout. The commands `run`, `check`, `capabilities`,
`verify` and `report` therefore send everything except their result (JSON or report text) to stderr,
so stdout can be piped.

## data

`ionmc data list | fetch ID | verify ID | path ID | import PATH --dataset ID | build nuclear-proton`.
`list`, `fetch`, `verify` and `path` take `--cache-dir` and (except `list`) `--offline`; see
[getting started](getting-started/index.md#external-data) and the [data layer](architecture/data-layer.md).
They exit with 0, 1 (an error message on stderr) or 2 (usage).

## run

```bash
ionmc run CONFIG.toml [--out DIR] [--overwrite]
```

Loads the file, verifies the datasets the configuration needs against their pinned hashes, validates
and reconciles it (all before transport), runs, and writes the bundle to `--out` or to `[output].path`.
Exactly one of the two must exist: `--out` wins. The command prints `bundle <path>` on stdout.

- An existing destination is refused (exit 2) unless `--overwrite`. `--overwrite` never deletes: the
  previous bundle is kept as `<path>.previous-<token>` and its location is printed on stderr
  (`previous bundle kept at ...`). Removing it is your decision.
- A run that completes with `valid = false` (only with `allow_invalid_result = true`) writes the bundle
  with a banner in the report, prints a note on stderr and exits 3.
- If the write fails after its temporary directory exists, that directory is kept as
  `<path>.failed-<token>`, and a `note:` line on stderr names it. The destination never holds a partial
  bundle.

## check

```bash
ionmc check CONFIG.toml
```

Prints the effective summary as canonical JSON. Use it to see the resolved configuration, and any refusal,
before spending run time.

## capabilities

```bash
ionmc capabilities [--nuclear]
```

Prints the capability report as JSON (`--nuclear`: the report for `physics.nuclear = true`).

## verify

```bash
ionmc verify DIR [reader budgets] [--rederive --cache-dir CACHE [--inputs ROLE=PATH ...]]
```

Tier 1 (see [persisted results](results/index.md)): reads `result.json` and `arrays.npz` only and prints
`{"binding_limits": [...], "tier1": "ok"}`. With `--rederive` (needs `--cache-dir`; `--cache-dir` and
`--inputs` are accepted only with it) tier 2 rebuilds the effective configuration under explicit trust and
adds a `tier2` record (`ok`, `status`, `differences`, `exclusions`, `not_rederivable`, `lookups_from`);
a failed re-derivation exits 6. `--inputs ROLE=PATH` names a trusted input file for a role of the manifest
(`material_index`, `density`, `lookup:<name>`); each file is checked against the manifest sha256.

## report

```bash
ionmc report DIR [reader budgets] [--check]
```

Prints `report.md` regenerated from the bundle (tier-1 load first). With `--check` nothing is printed
but `report ok (N files)`: the regenerated `report.md` and every `plots/*.svg` are compared byte for byte
with the files on disk, and a difference exits 6. Nothing is written.

## Reader budgets

`verify` and `report` take every budget of `ionmc.results.load` as a positive-integer option; an absent
option keeps the default of the [results format](results/format.md#reader-budgets):

| Option | Bounds |
|---|---|
| `--max-total-bytes` | `arrays.npz` and all arrays together |
| `--max-array-bytes` | one array |
| `--max-json-bytes` | `result.json` |
| `--max-central-directory-bytes` | the npz central directory |
| `--max-json-tokens` | JSON tokens of `result.json` |
| `--max-json-depth` | JSON nesting depth of `result.json` |
| `--max-derived-bytes` | `report.md` and `plots/*.svg`, each and together |

## Exit codes

The commands `run`, `check`, `capabilities`, `verify` and `report` map a failure to an exit code by the
first matching row of this table (subclasses precede their bases). Failures print
`error: <ExceptionName>: <message>` on stderr; a defect (code 70) also prints the traceback.

| Exception | Exit | Meaning | Bundle |
|---|---|---|---|
| `_UsageError` | 2 | misuse found after parsing (no output path, existing destination without `--overwrite`, `--rederive` without `--cache-dir`); argparse usage errors also exit 2 | none |
| `ConfigError` | 2 | TOML syntax, unknown table or key, wrong type, missing key, unknown `stopping_model`, sha256 mismatch of an input file | none |
| `DatasetError` | 1 | a needed cache object is missing or does not match its pinned hash | none |
| `UnsupportedCombinationError` | 1 | as above, reached through the stopping selection (a missing or mismatched cache object, a dataset parse failure) | none |
| `UnsupportedCombinationError` | 4 | the combination is outside the declared capabilities, or reconciliation failed | none |
| `BackendUnavailableError` | 5 | the backend or device is not available; no fallback is made | none |
| `TransportLimitError` | 3 | the run hit a transport limit and `allow_invalid_result = false` | none |
| `ResultFormatError` | 6 | a bundle defect found by `verify`, `report`, `report --check` or by the writer | none |
| `VerificationError` | 6 | tier 2 could not be carried out or failed (a failed `--rederive` also exits 6) | none |
| any other exception | 70 | a defect | none, traceback printed |

Two outcomes are not exceptions: success exits 0, and a run that completed with `valid = false` exits 3
with the bundle written.

Not a clinical claim.
