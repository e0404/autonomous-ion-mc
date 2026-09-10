# 0005 — Shared-source execution model, reference path, RNG and equivalence criteria

- Status: accepted (provisional until Stage 1 re-examines it under real transport)
- Date: 2026-09-10
- Task: DEV-002
- Affects: architecture, reproducibility, numerical accuracy, validation strategy, maintainability

## Problem

`EXPERIMENT.md` requires one physical model executed through a reference
Python path, Warp CPU and Warp CUDA, with the reference path being an
"independent or minimally transformed correctness oracle", and forbids
independently maintained implementations that can silently diverge. The
project must decide *how* a single source can serve all three paths, what
precision each path has, how random numbers are shared between them, and by
what criterion their results are declared equivalent (decision 0001 covers
only deterministic Warp CPU-vs-CUDA parity).

## Context

Facts established by the host-runner probe
`infrastructure/diagnostics/warp_execution_model_probe.py` on the experiment
workstation (Warp 1.17.0, CUDA 12.9, RTX A6000; run
`RUN-20260910T004325Z-472e71df`, task DEV-002, SHA `090cd34`):

- A `@wp.func`-decorated function (fixed `float` or `typing.Any` argument
  types) **is callable from Python scope** and returns a Python float; the
  kernel object's underlying Python function is accessible via `.func`.
- Called from Python scope, Warp's math builtins evaluate in **float32**
  (`2 log(100)/100 + sqrt(100)` came out 1.3 × 10⁻⁹ off the float64 value,
  exactly the float32 rounding of `log`). Python-scope execution of
  Warp-bound source is therefore *not* a float64 oracle.
- `wp.rand_init`, `wp.randf` etc. are **not callable from Python scope**
  (`Couldn't find a function 'rand_init' compatible with the arguments`).
- The kernel-side generator is bit-exactly reproduced by a 3-line Python
  PCG mirror on both CPU and CUDA (states for `rand_init(42, i)`, `i<8`).
- `wp.atomic_add` on 3-D `float32` and `float64` arrays works on CPU and
  CUDA; `wp.DeterministicMode` exists.
- Generic (`Any`-typed) functions instantiate for float32 and float64; the
  float64 kernel reproduced the numpy value exactly, float32 to 7 × 10⁻⁸.
- A data-dependent `while` loop with per-thread state (iterative range
  integration, up to 550 steps) runs on both devices; float32 accumulation
  differed from float64 Python by 1.5 × 10⁻⁶ relative and between CPU and
  CUDA by 9 × 10⁻⁸.
- Warp's documentation: literals inside generic functions are typed
  float32/int32 unless wrapped as `type(x)(literal)`; kernels on the CPU
  device execute serially.

## Candidate approaches

1. **Two implementations** (numpy reference, Warp kernels) kept in sync by
   tests. Rejected: this is the "independently maintained implementations"
   failure mode the constraint forbids; divergence is detected only where a
   test happens to look.
2. **Warp-bound source executed in Python scope as the reference.** Works
   (probe), but the reference then computes transcendentals in float32 and
   requires Warp to be installed; it cannot serve as a float64 oracle nor run
   in a Warp-less CI job.
3. **One source, bound at import time to one of three math namespaces
   (selected).** Physics functions are written against a small namespace
   `m` (`m.sqrt`, `m.log`, `m.where`, …) and decorated with `@func`. Under
   the `warp` binding `m.*` are Warp builtins and `@func` is `wp.func`;
   under the `python` binding `m.*` are `math` functions in float64 and
   `@func` is the identity; under the `numpy` binding `m.*` are ufuncs and
   the same source evaluates on arrays. A loader executes the identical
   source file under a second module name with a different binding.

## Selected approach

Approach 3, implemented in `src/ionmc/backend/mathlib.py` and
`src/ionmc/backend/reference.py`, with these rules:

1. **Shared source** (`ionmc/physics/*`): plain `float` scalar arguments
   (float32 in Warp), small tables as array arguments, `m.where` instead of
   data-dependent `if`, no Python `**`/`%` on values, module-level Python
   floats as constants. The Warp instantiation is **single precision**; a
   float64 Warp variant is deferred until a precision study needs it, because
   Warp requires `type(x)(literal)` wrapping in generic code, which the
   numpy binding cannot express.
2. **Reference path** = the `python` binding: identical source text, CPython,
   float64, no Warp dependency. The `numpy` binding is the fast reference
   for validation scripts and must be bitwise identical to the `python`
   binding (same operation order).
3. **Kernels** (`ionmc/backend/warp_kernels.py`) contain no physics: they
   unpack a parameter struct, index the per-thread inputs and call the shared
   functions. The Python-scope driver of the reference path mirrors exactly
   this thin layer.
4. **Random numbers.** Physics sampling functions are pure maps from uniform
   variates to outcomes; drawing variates happens in the execution layer.
   Streams are counter-based per history: `state = rand_init(seed, history_index)`
   with a globally unique history index, secondaries continue the parent's
   stream. The reference path draws from `ionmc.rng`, a pure-Python mirror
   of Warp's PCG hash that is bit-exact for states, integers and uniforms
   (`randf` uses 24 bits, exactly representable) and float32-accurate for
   `randn`. Identical streams on Python, Warp CPU and Warp CUDA make
   cross-backend comparisons of stochastic output far tighter than
   statistics alone.
5. **Equivalence criteria** (fixed here, before the first comparison):
   - `numpy` vs `python` binding: bitwise equal (`rtol` 1e-14 guard).
   - float64 reference vs float32 Warp kernel (CPU or CUDA): mixed criterion
     `|c − r| ≤ atol + rtol |r|` with `rtol = 1e-5, atol = 0` for
     stopping power and `rtol = 2e-5` for the 200-step range integral. These
     bound float32 rounding on ~10² operations with margin; a genuine
     algorithmic divergence is orders of magnitude larger.
   - Warp CPU vs Warp CUDA: decision 0001's methodology with per-kernel
     classes recorded in decision 0006.
   - Stochastic transport output (Stage 1) needs its own criterion; not
     covered here.
6. **Determinism.** `wp.DeterministicMode` will be exposed as an option when
   atomics are introduced (Stage 1).

## Rationale

The binding mechanism gives literally one source for all paths while keeping
the oracle in float64 and Warp-independent, which is stronger than approach 2
on both counts and avoids approach 1's divergence risk. Approach 2 was in
fact found to be incomplete during the DEV-002 host test run: a Warp-bound
function that uses `wp.where` raises "function is undefined" when called
from Python scope, because that builtin has no Python-scope implementation,
so the Python-scope fallback cannot serve as a general reference path. The probe showed all
Warp features the design relies on (Python-scope fallback, generic functions,
dynamic loops with per-thread state, atomics) work on both devices.

## Expected tradeoffs

- Shared functions carry long explicit parameter lists (no dataclasses in
  Warp); kernels bundle them in a `wp.struct`.
- Vectorisation constraints (`m.where`, no value-dependent `if`) shape the
  physics code style; loops with data-dependent trip counts are allowed but
  then execute per element in the numpy binding only if the trip count is
  data-independent — the CSDA integration uses a fixed step count for this
  reason.
- Warp physics runs in float32 only for now.
- The RNG mirror must be re-verified against a kernel dump whenever Warp is
  upgraded (test `tests/ionmc/test_rng.py` pins the values).
- The stress test so far is deterministic kernels with loops and per-thread
  state; real transport (secondary stacks, atomics, early exits) may expose
  limits — hence **provisional**.

## Validation strategy

- Unit tests of the bindings and loader (`tests/ionmc/test_backend.py`).
- RNG mirror pinned to the kernel dump (`tests/ionmc/test_rng.py`).
- `validation/v0_stopping_power.py` on the host runner: reference vs Warp CPU
  vs Warp CUDA under the criteria above (results in decision 0006 and the
  DEV-002 validation record).

## Later validation outcome

DEV-002 host run `RUN-20260910T072536Z-f9f8c5ff` (SHA `5ef05c7`): the
shared stopping-power source compiled and ran on Warp CPU and CUDA once loop
accumulators were declared as dynamic variables (`float(0.0)`; Warp refuses
to mutate a literal-initialised constant inside a dynamic loop — recorded
as a shared-source rule). float32 Warp versus the float64 Python reference
agreed to 8.7 × 10⁻⁶ relative for the stopping power (rtol 1e-5 criterion
met with only 1.15× margin) and to 3 × 10⁻⁷ for the 200-step range
integral; CPU and CUDA agreed to 9.5 × 10⁻⁷ absolute (transcendental class,
normalized 0.028). The numpy and Python bindings were bitwise identical.

The small margin was traced (decision 0006, *Finding on the
reference-vs-float32 margin*) to float32 cancellation in ``gamma^2 - 1``,
not to the execution model; rewriting the kinematics as ``tau (tau + 2)``
reduced the float32 error of S to 4.2 × 10⁻⁷ on Warp CPU and CUDA alike
(host run `RUN-20260910T080455Z-e666b469`, SHA `d23c78a`; CPU and CUDA
differ by at most one float32 ULP, see decision 0006). A second shared-source rule follows from it: quantities
that lose precision by cancellation in float32 must be written in
cancellation-free form, because the Warp instantiation is single precision
while the reference is float64. To be revisited at Stage 1 under real
transport.
