# Architecture

This section documents the implemented Monte Carlo software architecture and
major internal interfaces. It describes the state after task `DEV-002`
(Stage 0 of the [roadmap](../development/roadmap.md)); transport, geometry
and scoring do not exist yet.

## Package layout

```
src/ionmc/
  constants.py        CODATA 2018 constants (MeV, cm, mol)
  units.py            internal unit system and conversion helpers
  materials.py        Element, Material, mean excitation energy, density-effect parameters
  particles.py        Particle definitions (proton, alpha, carbon-12, oxygen-16)
  rng.py              pure-Python mirror of Warp's PCG random-number generator
  physics/
    stopping.py       shared-source Bethe stopping power and CSDA range
    barkas_table.py   tabulated Ashley-Ritchie-Brandt function and b parameters
  backend/
    mathlib.py        math-namespace bindings (warp / python / numpy) and @func
    reference.py      loader executing a shared-source module under a binding
    warp_kernels.py   Warp kernels and launch helpers (float32)
  stopping_power.py   Python API: AnalyticStoppingPower, Corrections, parameters
  reference_data/     small provenance-tagged reference subsets used by tests
```

## Shared-source execution model

Decision `0005` fixes how one physics source runs on three execution paths.

Physics modules import a *math namespace* ``m`` from
``ionmc.backend.mathlib`` and decorate their functions with ``@func``. The
namespace is chosen when the module is **loaded**:

| binding | ``m.sqrt`` etc. | ``@func`` | precision | role |
|---|---|---|---|---|
| ``warp`` (default if Warp is installed) | Warp builtins | ``wp.func`` | float32 | kernels on ``cpu`` and ``cuda:N`` |
| ``python`` | ``math`` | identity | float64 | **reference path** (oracle) |
| ``numpy`` | numpy ufuncs | identity | float64, vectorised | fast reference / validation |

``ionmc.backend.reference.load_bound_module("ionmc.physics.stopping", "python")``
executes the identical source file under a second module name with the
requested binding. The reference path therefore needs no Warp installation
and no second implementation: the only transformation is the binding of
``m``. Warp kernels (``backend/warp_kernels.py``) contain no physics; they
unpack a ``wp.struct`` of parameters, index the per-thread inputs and call
the shared functions.

Rules that shared source must follow (so that the three bindings agree):

- plain ``float`` scalar arguments (float32 in Warp) and small tables passed
  as arrays; no Python containers or dataclasses;
- data-dependent selection through ``m.where(cond, a, b)`` rather than
  ``if`` on values;
- loops with data-independent trip counts where the numpy binding must
  vectorise (the CSDA integration uses a fixed step count);
- no reliance on Python ``**`` or ``%`` semantics on values.

Warp's Python-scope fallback (``Function.__call__`` executes the original
Python function) also makes the Warp-bound functions callable from Python,
but with float32 builtins; the ``python`` binding is the float64 oracle.

## Random numbers

Streams are counter-based per history: ``rand_init(seed, history_index)``
with a globally unique history index; secondaries continue the parent's
stream. ``ionmc.rng`` mirrors Warp's PCG hash bit-exactly for states,
integers and 24-bit uniforms (verified against a kernel dump on Warp CPU and
CUDA), so the reference path draws the same variates as the kernels. Physics
sampling functions take variates as arguments; drawing them is the execution
layer's job. Transport itself is Stage 1 work.

## Execution backends and devices

The Python API selects the path explicitly:

```python
from ionmc import WATER, PROTON, AnalyticStoppingPower

ref = AnalyticStoppingPower(WATER, PROTON, path="python")       # float64 oracle
fast = AnalyticStoppingPower(WATER, PROTON, path="numpy")       # vectorised
gpu = AnalyticStoppingPower(WATER, PROTON, path="warp", device="cuda:0")
```

Warp CPU kernels execute serially (Warp documentation); the CPU device is a
correctness and portability path, CUDA is the performance path.

## Units and conventions

See decision `0004` and ``ionmc.units``: MeV, mm, g/cm³; mass stopping power
in MeV cm²/g and mass range in g/cm² as in the ICRU tables; right-handed
Cartesian coordinates in mm with no privileged beam axis.
