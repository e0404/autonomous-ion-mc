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
   compiled. Randomness and table interpolation are *inputs* to these
   functions, supplied by backend-specific adapters (numpy `Generator` and
   `numpy.interp` in the reference driver; `wp.rand_init`/`wp.randf` and
   array lookups in kernels).

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

- Warp transport state defaults to float32 with float64 accumulators for
  scoring; the reference backend is float64 throughout. A float64 transport
  mode on the Warp backends is provided for numerical falsification probes.
- Random numbers: all backends use the same counter-based generator,
  Philox4x32-10 (Salmon et al., SC'11; the Random123 construction). Warp's
  built-in `wp.rand_init`/`wp.randf` is **not** used: it is a stateful
  32-bit PCG hash whose streams share one 2³² state space, and a measurement
  on the installed Warp 1.17.0 (research report `warp-architecture.md`)
  found 20 % of draws revisiting states already used by other histories at
  1e6 histories × 2000 draws, and 260 of 1e6 histories bit-identical between
  seeds `s` and `s+1`. Philox is implemented once as a `@wp.func` over
  `uint32` words for the Warp kernels and once over Python integers for the
  reference backend; the two were verified bit-identical on the published
  known-answer vectors and 20 000 random blocks. The counter is
  `(history_id, genealogy_id, step, sub-draw)` and the key is `(seed,
  stream)`, so every history, every secondary particle (identified by its
  parent and birth index, never by a queue slot) and every batch has a
  disjoint, reproducible stream independent of thread count or launch
  partitioning. Because the streams are identical across backends, the
  reference backend and the float64 build of the Warp kernels can be
  compared trajectory by trajectory; production float32 Warp results are
  compared statistically.

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
