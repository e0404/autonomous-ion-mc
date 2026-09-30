"""Warp kernels used by the shared-physics tests.

Kept outside ``test_*.py`` modules: pytest rewrites assertions in test modules
by compiling them itself, which breaks Warp's source inspection of kernels.
"""

from typing import Any

import warp as wp

from ionmc.transport.shared import physics

wpv = physics("warp")


@wp.kernel
def two_losses(out: wp.array2d(dtype=Any), mean: Any, var: Any):
    i = wp.tid()
    st = wp.rand_init(7, i)
    a, st = wpv.sample_energy_loss(st, mean, var)
    b, st = wpv.sample_energy_loss(st, mean, var)
    out[i, 0] = a
    out[i, 1] = b


@wp.kernel
def sequential_losses(out: wp.array(dtype=Any), mean: Any, var: Any, n: int):
    st = wp.rand_init(11, wp.tid())
    for k in range(n):
        x, st = wpv.sample_energy_loss(st, mean, var)
        out[k] = x


@wp.kernel
def f_dm_kernel(out: wp.array(dtype=Any), pv: Any, p1v1: Any):
    i = wp.tid()
    out[i] = wpv.f_dm(pv, p1v1)
