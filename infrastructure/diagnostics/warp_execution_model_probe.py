"""Probe Warp 1.17 behaviours that decide the shared-source execution model.

Run on the host runner (warp + numpy available):

    python warp_scope_probe.py

Each probe prints one JSON line so the result is machine-readable.
"""

from __future__ import annotations

import json
import traceback
from typing import Any

import numpy as np
import warp as wp


def report(probe_name, **kw):
    print(json.dumps({"probe": probe_name, **kw}, default=str))


wp.init()

# --- 1. @wp.func called from Python scope -------------------------------


@wp.func
def bethe_like(e: float, a: float) -> float:
    return a * wp.log(e) / e + wp.sqrt(e)


try:
    v = bethe_like(100.0, 2.0)
    expected = 2.0 * np.log(100.0) / 100.0 + np.sqrt(100.0)
    report("wp_func_python_scope", ok=True, value=v, expected=expected,
           type=type(v).__name__)
except Exception as exc:  # noqa: BLE001
    report("wp_func_python_scope", ok=False, error=repr(exc))

# --- 2. Kernel underlying python function accessible? --------------------


@wp.kernel
def eval_kernel(e: wp.array(dtype=float), a: float, out: wp.array(dtype=float)):
    i = wp.tid()
    out[i] = bethe_like(e[i], a)


try:
    f = getattr(eval_kernel, "func", None)
    report("kernel_func_attr", present=f is not None, callable=callable(f))
except Exception as exc:  # noqa: BLE001
    report("kernel_func_attr", ok=False, error=repr(exc))

# --- 3. Python-scope RNG builtins ----------------------------------------

for fn_name in ("rand_init", "randf", "randn", "randi", "randu"):
    report("rng_builtin_exported", name=fn_name, present=hasattr(wp, fn_name))

try:
    s = wp.rand_init(42, 7)
    r1 = wp.randf(s)
    r2 = wp.randf(s)
    report("rng_python_scope", ok=True, state=s, state_type=type(s).__name__,
           r1=r1, r2=r2, state_advances=(r1 != r2))
except Exception as exc:  # noqa: BLE001
    report("rng_python_scope", ok=False, error=repr(exc))

# --- 4. Kernel-side RNG dump for bit-exact pure-Python mirror -------------


@wp.kernel
def rng_dump(seed: int, states: wp.array(dtype=wp.uint32),
             floats: wp.array(dtype=float), ints: wp.array(dtype=wp.uint32),
             normals: wp.array(dtype=float)):
    i = wp.tid()
    st = wp.rand_init(seed, i)
    states[i] = st
    floats[i] = wp.randf(st)
    ints[i] = wp.randu(st)
    normals[i] = wp.randn(st)


def pcg(state: int) -> int:
    state = (state * 747796405 + 2891336453) & 0xFFFFFFFF
    word = (((state >> ((state >> 28) + 4)) ^ state) * 277803737) & 0xFFFFFFFF
    return ((word >> 22) ^ word) & 0xFFFFFFFF


def py_rand_init(seed: int, offset: int) -> int:
    return pcg((seed & 0xFFFFFFFF) + pcg(offset & 0xFFFFFFFF) & 0xFFFFFFFF)


for device in [d.alias for d in wp.get_devices()]:
    try:
        n = 8
        states = wp.zeros(n, dtype=wp.uint32, device=device)
        floats = wp.zeros(n, dtype=float, device=device)
        ints = wp.zeros(n, dtype=wp.uint32, device=device)
        normals = wp.zeros(n, dtype=float, device=device)
        wp.launch(rng_dump, dim=n, inputs=[42, states, floats, ints, normals],
                  device=device)
        wp.synchronize_device(device)
        st = states.numpy().tolist()
        mirror = [py_rand_init(42, i) for i in range(n)]
        report("rng_kernel_dump", device=device, states=st,
               py_mirror_states=mirror, mirror_matches=(st == mirror),
               floats=floats.numpy().tolist(), ints=ints.numpy().tolist(),
               normals=normals.numpy().tolist())
    except Exception as exc:  # noqa: BLE001
        report("rng_kernel_dump", device=device, ok=False,
               error=repr(exc), tb=traceback.format_exc()[-800:])

# --- 5. float64 atomics on 3D arrays -------------------------------------


@wp.kernel
def atomic3d_f64(out: wp.array3d(dtype=wp.float64)):
    i = wp.tid()
    wp.atomic_add(out, i % 2, i % 3, i % 5, wp.float64(1.0))


@wp.kernel
def atomic3d_f32(out: wp.array3d(dtype=wp.float32)):
    i = wp.tid()
    wp.atomic_add(out, i % 2, i % 3, i % 5, 1.0)


for device in [d.alias for d in wp.get_devices()]:
    for name, kern, dt in (("f64", atomic3d_f64, wp.float64),
                           ("f32", atomic3d_f32, wp.float32)):
        try:
            out = wp.zeros((2, 3, 5), dtype=dt, device=device)
            wp.launch(kern, dim=3000, inputs=[out], device=device)
            wp.synchronize_device(device)
            total = float(out.numpy().sum())
            report("atomic_add_3d", device=device, dtype=name, ok=True,
                   total=total, expected=3000.0)
        except Exception as exc:  # noqa: BLE001
            report("atomic_add_3d", device=device, dtype=name, ok=False,
                   error=repr(exc)[:400])

# --- 6. Same wp.func executed in kernel vs Python scope ------------------

for device in [d.alias for d in wp.get_devices()]:
    try:
        e = np.linspace(1.0, 250.0, 64, dtype=np.float32)
        ea = wp.array(e, dtype=float, device=device)
        out = wp.zeros(64, dtype=float, device=device)
        wp.launch(eval_kernel, dim=64, inputs=[ea, 2.0, out], device=device)
        wp.synchronize_device(device)
        k = out.numpy()
        py = np.array([bethe_like(float(x), 2.0) for x in e], dtype=np.float64)
        report("func_kernel_vs_python", device=device, ok=True,
               max_abs_diff=float(np.max(np.abs(k.astype(np.float64) - py))),
               max_rel_diff=float(np.max(np.abs(k - py) / np.abs(py))))
    except Exception as exc:  # noqa: BLE001
        report("func_kernel_vs_python", device=device, ok=False,
               error=repr(exc)[:400])

# --- 6b. dtype-generic wp.func (Any) instantiated for float32 and float64 -


@wp.func
def generic_bethe_like(e: Any, a: Any):
    return a * wp.log(e) / e + wp.sqrt(e)


@wp.kernel
def generic_eval_f32(e: wp.array(dtype=wp.float32), out: wp.array(dtype=wp.float32)):
    i = wp.tid()
    out[i] = generic_bethe_like(e[i], wp.float32(2.0))


@wp.kernel
def generic_eval_f64(e: wp.array(dtype=wp.float64), out: wp.array(dtype=wp.float64)):
    i = wp.tid()
    out[i] = generic_bethe_like(e[i], wp.float64(2.0))


for device in [d.alias for d in wp.get_devices()]:
    for name, kern, npdt, wpdt in (("f32", generic_eval_f32, np.float32, wp.float32),
                                   ("f64", generic_eval_f64, np.float64, wp.float64)):
        try:
            e = np.linspace(1.0, 250.0, 64, dtype=npdt)
            ea = wp.array(e, dtype=wpdt, device=device)
            out = wp.zeros(64, dtype=wpdt, device=device)
            wp.launch(kern, dim=64, inputs=[ea, out], device=device)
            wp.synchronize_device(device)
            k = out.numpy().astype(np.float64)
            py = np.array([2.0 * np.log(float(x)) / float(x) + np.sqrt(float(x)) for x in e])
            report("generic_func_dtype", device=device, dtype=name, ok=True,
                   max_rel_diff=float(np.max(np.abs(k - py) / np.abs(py))))
        except Exception as exc:  # noqa: BLE001
            report("generic_func_dtype", device=device, dtype=name, ok=False,
                   error=repr(exc)[:600])

try:
    v = generic_bethe_like(100.0, 2.0)
    report("generic_func_python_scope", ok=True, value=v, type=type(v).__name__)
except Exception as exc:  # noqa: BLE001
    report("generic_func_python_scope", ok=False, error=repr(exc)[:600])

# --- 6c. data-dependent while loop with per-thread state (range integration)


@wp.func
def inv_speed(e: float) -> float:
    return e / (wp.log(e) + 1.0)


@wp.kernel
def iterative_range(e0: wp.array(dtype=float), out: wp.array(dtype=float),
                    steps: wp.array(dtype=wp.int32)):
    i = wp.tid()
    e = e0[i]
    r = float(0.0)
    n = int(0)
    while e > 1.0 and n < 100000:
        de = wp.min(0.01 * e, e - 1.0)
        r += de * inv_speed(e - 0.5 * de)
        e -= de
        n += 1
    out[i] = r
    steps[i] = n


def py_iterative_range(e):
    r = 0.0
    n = 0
    while e > 1.0 and n < 100000:
        de = min(0.01 * e, e - 1.0)
        r += de * inv_speed(e - 0.5 * de)
        e -= de
        n += 1
    return r, n


for device in [d.alias for d in wp.get_devices()]:
    try:
        e = np.array([5.0, 20.0, 100.0, 250.0], dtype=np.float32)
        ea = wp.array(e, dtype=float, device=device)
        out = wp.zeros(4, dtype=float, device=device)
        steps = wp.zeros(4, dtype=wp.int32, device=device)
        wp.launch(iterative_range, dim=4, inputs=[ea, out, steps], device=device)
        wp.synchronize_device(device)
        k = out.numpy().astype(np.float64)
        py = np.array([py_iterative_range(float(x))[0] for x in e])
        report("iterative_while_kernel", device=device, ok=True,
               values=k.tolist(), py_values=py.tolist(),
               steps=steps.numpy().tolist(),
               max_rel_diff=float(np.max(np.abs(k - py) / np.abs(py))))
    except Exception as exc:  # noqa: BLE001
        report("iterative_while_kernel", device=device, ok=False,
               error=repr(exc)[:600])

# --- 7. deterministic mode availability ----------------------------------

report("deterministic_mode", present=hasattr(wp, "DeterministicMode"),
       config_attr=hasattr(wp.config, "deterministic"))
report("warp_version", version=wp.config.version,
       devices=[d.alias for d in wp.get_devices()])
