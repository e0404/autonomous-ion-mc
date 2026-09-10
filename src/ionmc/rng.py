"""Pure-Python mirror of NVIDIA Warp's random-number generator.

Warp's kernel-side generator (``warp/native/rand.h``, Warp 1.17.0) is the
32-bit PCG-style hash of Jarzynski & Olano, *Hash Functions for GPU
Rendering*, JCGT 9(3), 2020. Its builtins (``wp.rand_init``, ``wp.randf`` ...)
cannot be called from Python scope (verified on the experiment workstation,
host run ``RUN-20260910T004325Z-472e71df``), so the reference Python execution
path uses this module instead. The integer functions reproduce the kernel
values **bit-exactly** (verified against a kernel dump on Warp CPU and CUDA in
the same run); ``randf`` is exact as well because ``(state >> 8) / 2**24`` is
exactly representable in both float32 and float64. ``randn`` uses float32
transcendentals in Warp and float64 here, so it agrees only to float32
precision.

Stream design (decision ``0005``): every history ``h`` of a run with seed ``s``
draws from the state ``rand_init(s, h)``; secondaries continue the parent's
stream. The same mapping is used by the Warp kernels, so reference Python and
Warp draw identical uniforms for the same history.
"""

from __future__ import annotations

import math

_MASK32 = 0xFFFFFFFF
_RANDN_EPSILON = 5.96e-8
_INV_2_24 = 1.0 / 16777216.0


def rand_pcg(state: int) -> int:
    """One step of the PCG-style hash on a 32-bit state."""
    b = (state * 747796405 + 2891336453) & _MASK32
    c = (((b >> ((b >> 28) + 4)) ^ b) * 277803737) & _MASK32
    return ((c >> 22) ^ c) & _MASK32


def rand_init(seed: int, offset: int | None = None) -> int:
    """Initial state for ``seed`` (and stream ``offset``), as ``wp.rand_init``."""
    if offset is None:
        return rand_pcg(seed & _MASK32)
    return rand_pcg(((seed & _MASK32) + rand_pcg(offset & _MASK32)) & _MASK32)


class RandomState:
    """Mutable 32-bit generator state mirroring a kernel-local ``uint32``.

    Warp's builtins take the state by reference; Python integers are
    immutable, so the state lives in this small object.
    """

    __slots__ = ("state",)

    def __init__(self, seed: int, offset: int | None = None) -> None:
        self.state = rand_init(seed, offset)

    @classmethod
    def from_state(cls, state: int) -> RandomState:
        obj = cls.__new__(cls)
        obj.state = state & _MASK32
        return obj

    def randu(self) -> int:
        """Uniform unsigned 32-bit integer (``wp.randu``)."""
        self.state = rand_pcg(self.state)
        return self.state

    def randi(self, low: int | None = None, high: int | None = None) -> int:
        """Signed 32-bit integer (``wp.randi``), optionally in ``[low, high)``."""
        self.state = rand_pcg(self.state)
        if low is None:
            value = self.state
            return value - (1 << 32) if value >= (1 << 31) else value
        assert high is not None
        return self.state % (high - low) + low

    def randf(self, low: float = 0.0, high: float = 1.0) -> float:
        """Uniform float in ``[low, high)`` (``wp.randf``); 24 random bits.

        The unit-interval value is bit-exact with the kernel; the scaled form
        is computed in float64 here and in float32 by Warp, so it agrees to
        float32 rounding only.
        """
        self.state = rand_pcg(self.state)
        u = (self.state >> 8) * _INV_2_24
        if low == 0.0 and high == 1.0:
            return u
        return (high - low) * u + low

    def randn(self) -> float:
        """Standard normal via Box-Muller (``wp.randn``); consumes two uniforms.

        Warp evaluates this in float32; agreement is to float32 precision.

        Warp's ``randn`` (``rand.h``) draws its two uniforms as the two
        operands of one C++ multiplication, ``sqrt(-2 log(randf(state) + eps))
        * cos(2 pi randf(state))``, whose evaluation order the language does
        not fix. The order mirrored here (``u1`` for the radius, ``u2`` for
        the angle) is the one observed on the experiment workstation for both
        the CPU and CUDA back ends; ``tests/ionmc/test_rng.py`` pins it and
        must be re-checked after a Warp or compiler upgrade.
        """
        u1 = self.randf()
        u2 = self.randf()
        return math.sqrt(-2.0 * math.log(u1 + _RANDN_EPSILON)) * math.cos(
            2.0 * math.pi * u2
        )
