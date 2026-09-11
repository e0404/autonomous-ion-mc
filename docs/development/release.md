# Release process and validation

`main` contains releases only; it is updated only from a release-ready state of
`develop` (`EXPERIMENT.md`, *Releases*). A release must have passed the **defined
release validation suite**, with documented validation status, reproducible benchmark
information, a version tag, and identification of the experiment configuration.

## The release validation suite

`validation/release_validation.py` (decision `0033`) **is** the defined suite. It runs
every milestone validation (V0–V5) plus the benchmark physics gates and the
performance-regression checks against the committed baselines (V6), aggregating a
single machine-readable report and exiting 0 **iff every** suite passes. Run it on the
controlled host runner:

```bash
env PYTHONPATH=/workspace/src python validation/release_validation.py \
    --require-cuda --cache-dir /cache/ionmc
```

The report stamps the experiment configuration (`ionmc` / Python / NumPy / Warp
versions, `--require-cuda`) and lists each suite's return code and pass/fail; its
`release_ready` flag is the gate. Composition is discovered by globbing
`validation/v*.py` (excluding the `warp_cuda_smoke` diagnostic), so new milestone
validations are included automatically; each is run with only the flags it supports.

## Release criteria

A `develop` state is release-ready when all of the following hold:

1. the release SHA has **passing required CI** (unit tests, lint/pre-commit,
   type-check, docs build) — the release validation suite is deliberately physics-only
   and does not re-run the unit/packaging checks, so CI health is a separate criterion;
2. the release validation suite passes on the host runner (CPU + CUDA);
3. the result is recorded for the exact SHA (`record_local_validation`);
4. the documented validation status and the reproducible benchmark baselines
   (`benchmarks/baselines/`) are committed;
5. the version (`ionmc.__version__`) is set for the release;
6. the experiment configuration is identified (stamped in the report).

## Cutting the release

Once the criteria hold at a release SHA, the release is tagged from `develop` to
`main` (a version tag), with the validation status and benchmark information
referenced. `EXPERIMENT.md` delegates the release decision to the autonomous system
unless release authorization is reserved to the operator; force-pushes and history
rewrites on `develop`/`main` remain prohibited.

## Deferred

Wiring the suite into GitHub CI (it needs a GPU and the datasets), a packaging /
wheel-build gate, and a persisted release-history record.
