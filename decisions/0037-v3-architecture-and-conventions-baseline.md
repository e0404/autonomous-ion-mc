# 0037 — V3 architecture, dependency and convention baseline

- Status: accepted
- Date: 2026-10-03
- Task: V3-001
- Affects: architecture, dependency footprint, reproducibility, numerical accuracy, maintainability

## Problem

Experiment v3 starts from infrastructure only. Before any physics is written,
the project needs a fixed set of baseline choices that every later task builds
on: how the reference Python path and the Warp CPU/CUDA paths share one
physical model, which runtime dependencies are allowed, and which units,
coordinate and precision conventions API boundaries use.

## Context

- `EXPERIMENT.md` requires a Python-first package with a reference Python
  path and Warp CPU/CUDA acceleration sharing the same physical model.
- The controlled GPU host runner executes committed code from a read-only
  snapshot in a virtual environment that contains only `numpy`, `warp-lang`
  (1.17.0) and `pytest`, without network access. Anything the package needs
  at runtime for validation must therefore come from those packages or the
  standard library.
- Verified on the installed Warp 1.17.0: functions decorated with `@wp.func`
  can be called directly from Python scope with scalar, vector and
  `wp.struct` arguments and execute their Python source with Warp builtins
  (`wp.sqrt`, `wp.exp`, loops). Array item indexing and `wp.rand_init` are
  not available in that mode.

## Candidate approaches

1. **Two separate implementations** (numpy reference, Warp kernels) kept in
   sync by tests. Rejected: this is exactly the "independently maintained
   implementations whose behaviour can silently diverge" the protocol warns
   against; the reference would validate only itself.
2. **Warp-only physics, reference path = Warp CPU.** Rejected: a separately
   compiled backend is not an independent correctness oracle, and the
   protocol requires a native Python execution path.
3. **Shared pure physics functions executed by three drivers (selected).**
   Physics is written once as `@wp.func` functions over scalars, vectors and
   structs. A reference driver loops over histories in Python and calls the
   same functions with Python floats; Warp CPU and CUDA kernels call them
   compiled. Random numbers and table values are *inputs* to these
   functions. Randomness comes from one shared counter-based generator
   (Philox4x32-10, see *Precision and randomness* below) implemented twice
   with bit-identical output: as a `@wp.func` for kernels and over Python
   integers for the reference driver. Table interpolation is split into a
   shared bin-location function, a backend-specific memory read (numpy
   array in the reference driver, `wp.array` in kernels) and a shared
   interpolation function, so only the read differs between backends.

## Selected approach

### Shared physics functions

- `ionmc.physics.*` modules contain `@wp.func` functions that are pure:
  they take sampled uniform random numbers and already-interpolated table
  values as arguments and return the physical result (energy loss, scattering
  angle, interaction outcome). They never index arrays or draw random numbers.
- The reference backend (`backend="python"`) executes these functions in
  CPython with float64 arithmetic, history by history. It is slow by design
  and is the correctness and inspection path.
- Warp backends (`backend="warp-cpu"`, `backend="warp-cuda"`) execute the
  same functions inside kernels. The two Warp devices share kernel source;
  their equivalence is checked by the parity methodology of decision 0001
  for deterministic kernels and by statistical criteria for transport.
- Independence of *physics* is never claimed from backend agreement. It
  comes from analytic, tabulated, measured and independent Monte Carlo
  evidence recorded per suite.

### Runtime dependency footprint

- Runtime dependencies: `numpy>=2` and `warp-lang>=1.17,<2` only.
- Persisted results use NumPy `.npz` plus JSON sidecars; human-readable
  reports are Markdown/plain text; plots are generated as SVG by a small
  internal writer. No plotting or HDF5 library is required.
- Development tooling (pytest, ruff, mypy, pre-commit) is an optional
  dependency group managed with `uv` and a committed `uv.lock`.

### Units (API boundaries and persisted outputs)

| Quantity | Unit |
|---|---|
| Kinetic energy | MeV (per ion); tables for ions may be stored per nucleon and are labelled `MeV/u` |
| Length, position, spacing | mm |
| Mass density | g/cm³ |
| Linear stopping power | MeV/mm; mass stopping power MeV·cm²/g where tabulated |
| Energy deposition | MeV per voxel per primary unless stated |
| Absorbed dose | Gy per primary (1 MeV/g = 1.602176634e-10 Gy) |
| LET | keV/µm |
| Fluence | mm⁻² per primary |
| Angles | rad |

Every persisted array carries a `units` attribute in its sidecar; functions
whose arguments are not in these units state the unit in the parameter name
(for example `energy_per_nucleon_mev`).

### Coordinates and grids

- Right-handed Cartesian world coordinates in mm.
- A regular grid is defined by `origin` (coordinates of the corner of voxel
  `(0,0,0)`, not its centre), `spacing` (3 positive values) and `shape`
  `(nx, ny, nz)`. Voxel `(i,j,k)` spans
  `origin + (i,j,k)*spacing` to `origin + (i+1,j+1,k+1)*spacing`.
- Arrays are stored C-ordered with index order `[ix, iy, iz]`, so
  `array.shape == (nx, ny, nz)`.
- Scoring grids are independent regular grids with their own origin,
  spacing and shape; they need not align with the transport grid.
- A source is defined by a position, a unit direction vector and optional
  lateral/energy spread; no axis is privileged.

### Precision and randomness

- Warp transport state and tables default to float32; the reference
  backend is float64 throughout. Scoring on the Warp backends accumulates
  into per-batch float32 arrays (one accumulator set per statistical batch,
  which also provides the batch variance) and reduces across batches in
  float64. A single float32 accumulator is not used because the archived
  measurement (`precision.py`, step 08) shows a relative accumulation error
  of about 1e-3 at 3e6 deposits per voxel, above the statistical error,
  whereas spreading the same deposits over 40 batch accumulators reduces it
  to below 1e-6. A float64 build of the same kernels (float64 state, tables
  and accumulators) is provided for numerical falsification probes and
  trajectory-level parity with the reference; it is not a production mode.
- Random numbers: all backends use the same counter-based generator,
  Philox4x32-10 (Salmon et al., SC'11; the Random123 construction). Warp's
  built-in `wp.rand_init`/`wp.randf` is **not** used anywhere: it is a
  stateful 32-bit PCG hash whose streams share one 2³² state space. The
  measurement scripts committed under `validation/scripts/rng/` (results
  tabulated in its README and in `docs/research/warp-architecture.md`) found,
  on Warp 1.17.0, 20 % of draws revisiting states already used by other
  histories at 1e6 histories × 2000 draws, and 260 of 1e6 histories
  bit-identical between seeds `s` and `s+1`. The same scripts check the
  Philox `@wp.func` against the Random123 known-answer vectors and against a
  pure-Python integer implementation (20 000 random blocks, bit-identical).
  The package implementation and its unit tests are delivered by the
  transport task (V3-003); until then this section specifies the design.

### Random stream allocation

Philox4x32-10 maps a 128-bit counter `(c0, c1, c2, c3)` and a 64-bit key
`(k0, k1)` to four 32-bit words; each word yields one uniform in (0, 1) as
`(w >> 8 + 0.5) · 2⁻²⁴` in float64 and `(w >> 9 + 0.5) · 2⁻²³` in float32
(never exactly 0 or 1; the 24-bit form rounds to exactly 1.0 in float32 for
the largest word, an error found during V3-003 planning and amended here on
2026-10-03, decision 0039). Streams are disjoint if and
only if no (counter, key) pair is used twice, which the following fixed
encoding guarantees by construction:

| Field | Content | Bound (exceeding it fails closed) |
|---|---|---|
| `k0`, `k1` | low and high 32 bits of the user seed (64-bit integer) | — |
| `c0` | global history index within the run, `0 … N−1`, never reset between batches or beamlets | N ≤ 2³² histories per run |
| `c1` | particle identifier within the history: `0` for the primary; a secondary born as the `b`-th child (`b` from 1) of a particle at generation `g` (primary has `g = 0`) gets `id = parent_id + b · 32^g` | at most 31 children per particle and 6 generations (`id < 2³⁰`); a 32nd child or a 7th generation is **never transported and never deposited as dose**: its kinetic energy is added to a separate `unaccounted_energy` tally and an `overflow` counter, and a nonzero counter makes the result invalid in every mode (the API raises `TransportLimitError` after the batch; a caller may opt in to receive the invalid result object, whose `valid` flag is `False` and whose dose arrays are marked not physically meaningful) |
| `c2` | draw-block index within the particle, incremented at every Philox call (four uniforms per block) | 2³² blocks per particle (unreachable; step limits bind first) |
| `c3` | purpose: `0` transport, `1` source sampling of the primary, `2` reserved | — |

Consequences:

- The same fail-closed rule applies to step-count truncation and queue
  overflow: truncated energy is tallied separately, counted, and invalidates
  the result; it is never folded into dose.
- Batch membership is `history_index mod n_batches` and beamlet membership
  is a range of history indices, so batches and beamlets never share
  counters and the batch estimator sees independent streams.
- A secondary's stream depends only on its genealogy (parent identifier and
  birth order), never on the slot it receives in a global queue, so results
  are independent of thread scheduling and launch partitioning apart from
  floating-point summation order in accumulators.
- The reference backend uses the identical encoding with Python integers, so
  reference and float64 Warp trajectories can be compared history by
  history; production float32 Warp results are compared statistically.
- Different seeds give different keys; there is no relation between the
  streams of seed `s` and `s+1`.

## Expected tradeoffs

- The reference backend will be orders of magnitude slower than Warp; it is
  used for small histories counts and cross-checks.
- Keeping physics functions free of array access and RNG calls adds adapter
  code in each backend; this is the price of one physical model.
- Pinning to numpy + warp restricts output formats; SVG/JSON/NPZ are
  sufficient and dependency-free.

## Validation strategy

- Unit tests call each shared function from Python and inside a Warp CPU
  kernel with identical arguments and compare results within float32
  tolerance (decision 0001 classes).
- The Philox implementation is tested against the Random123 known-answer
  vectors in both the kernel and the Python form, and stream disjointness
  (no repeated counter/key pairs across histories, secondaries and batches)
  is asserted by construction tests.
- Trajectory-level parity: reference backend versus the float64 Warp CPU
  build for the first K histories with identical random streams.
- Transport-level comparisons between backends use statistical criteria
  defined in the numerical parity suite.

## Later validation outcome

To be appended as tasks merge.
