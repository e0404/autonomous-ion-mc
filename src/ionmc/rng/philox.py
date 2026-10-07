# mypy: ignore-errors
# (Warp function sources use runtime precision types in annotations; see _wpfunc.py.)
"""Philox4x32-10 counter-based generator, bit-identical in Warp kernels and Python ints.

Construction: Salmon, Moraes, Dror, Shaw, "Parallel random numbers: as easy as 1, 2, 3",
SC'11 (the Random123 library); constants and round function as in Random123
``philox4x32_R(10)``. The Python form below and the ``@wp.func`` are checked against the
Random123 known-answer vectors and against each other in ``tests/ionmc/test_rng.py``.

Counter/key encoding (decision 0037)::

    key     = (seed & 0xFFFFFFFF, seed >> 32)                      seed in [0, 2^64)
    counter = (c0 history, c1 genealogy id, c2 block, c3 purpose)
    c0  global history index, 0 <= c0 < 2^32
    c1  genealogy id: 0 for the primary; child b (1..31) of a particle of generation g
        (primary g = 0, g <= 5) has id = parent_id + b * 32**g  (id < 2^30)
    c2  draw-block index within the particle (four uniforms per block), < 2^32
    c3  purpose: 0 transport (EM), 1 source sampling of the primary, 2 nuclear interactions
        (``PURPOSE_NUCLEAR``; formerly reserved; decision 0041 section 4 amends 0037)

Uniform mapping (decision 0039, amending 0037): ``u = (k + 0.5) * 2^-bits`` with
``k = w >> 8, bits = 24`` for float64 and ``k = w >> 9, bits = 23`` for float32, so that
``u`` is strictly inside (0, 1) in the precision of use (the 24-bit form rounds to 1.0 in
float32 for the largest word).
"""

import functools
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func
from ionmc.errors import CounterOverflowError

wp.set_module_options({"enable_backward": False})

_M0 = 0xD2511F53
_M1 = 0xCD9E8D57
_W0 = 0x9E3779B9
_W1 = 0xBB67AE85
_MASK32 = 0xFFFFFFFF

MAX_HISTORIES = 2**32
MAX_GENEALOGY_ID = 2**30
MAX_CHILDREN = 31
MAX_PARENT_GENERATION = 5  # a child of a generation-5 particle (generation 6) is the last allowed
PURPOSE_TRANSPORT = 0
PURPOSE_SOURCE = 1
PURPOSE_NUCLEAR = 2
PURPOSE_RESERVED = PURPOSE_NUCLEAR  # alias kept for existing imports (decision 0041 section 4)

_K0 = wp.constant(wp.uint64(_M0))
_K1 = wp.constant(wp.uint64(_M1))
_WEYL0 = wp.constant(wp.uint32(_W0))
_WEYL1 = wp.constant(wp.uint32(_W1))
_MASK = wp.constant(wp.uint64(_MASK32))


@wp.func
def philox4x32_10(c: wp.vec4ui, k: wp.vec2ui) -> wp.vec4ui:
    """Philox4x32-10 block function: counter ``c`` and key ``k`` to four 32-bit words."""
    c0 = c[0]
    c1 = c[1]
    c2 = c[2]
    c3 = c[3]
    k0 = k[0]
    k1 = k[1]
    for _r in range(10):
        p0 = _K0 * wp.uint64(c0)
        p1 = _K1 * wp.uint64(c2)
        n0 = wp.uint32(p1 >> wp.uint64(32)) ^ c1 ^ k0
        n2 = wp.uint32(p0 >> wp.uint64(32)) ^ c3 ^ k1
        c0 = n0
        c1 = wp.uint32(p1 & _MASK)
        c2 = n2
        c3 = wp.uint32(p0 & _MASK)
        k0 = k0 + _WEYL0
        k1 = k1 + _WEYL1
    return wp.vec4ui(c0, c1, c2, c3)


@wp.func
def philox_block(
    history: wp.uint32, genealogy: wp.uint32, block: wp.uint32, purpose: wp.uint32, k: wp.vec2ui
) -> wp.vec4ui:
    """Philox block for the counter ``(history, genealogy, block, purpose)``."""
    return philox4x32_10(wp.vec4ui(history, genealogy, block, purpose), k)


@functools.cache
def make_philox(real: type) -> SimpleNamespace:
    """Return ``(philox4x32_10, philox_block, u01)`` for precision ``real``.

    ``u01(w)`` maps a 32-bit word to a uniform in (0, 1) of type ``real`` (module docstring).
    """
    name = check_real(real)
    if name == "float64":
        shift = 8
        scale = 1.0 / 16777216.0  # 2**-24
    else:
        shift = 9
        scale = 1.0 / 8388608.0  # 2**-23
    sh = wp.constant(wp.uint32(shift))
    sc = wp.constant(real(scale))

    @named_func(name)
    def u01(x: wp.uint32) -> real:
        return (real(x >> sh) + real(0.5)) * sc

    return SimpleNamespace(
        philox4x32_10=philox4x32_10, philox_block=philox_block, u01=u01, real=name
    )


def philox4x32_10_py(
    counter: tuple[int, int, int, int], key: tuple[int, int]
) -> tuple[int, int, int, int]:
    """Philox4x32-10 over Python integers (reference backend and tests)."""
    c0, c1, c2, c3 = counter
    k0, k1 = key
    for v in (c0, c1, c2, c3, k0, k1):
        if not 0 <= v <= _MASK32:
            raise ValueError("counter and key words must be 32-bit unsigned integers")
    for _ in range(10):
        p0 = _M0 * c0
        p1 = _M1 * c2
        c0, c1, c2, c3 = (
            ((p1 >> 32) ^ c1 ^ k0) & _MASK32,
            p1 & _MASK32,
            ((p0 >> 32) ^ c3 ^ k1) & _MASK32,
            p0 & _MASK32,
        )
        k0 = (k0 + _W0) & _MASK32
        k1 = (k1 + _W1) & _MASK32
    return c0, c1, c2, c3


def u01_py(word: int, precision: str) -> float:
    """Uniform in (0, 1) from a 32-bit word, identical to the ``u01`` Warp function.

    The returned Python float is exactly representable in the requested precision.
    """
    if not 0 <= word <= _MASK32:
        raise ValueError("word must be a 32-bit unsigned integer")
    if precision == "float64":
        return ((word >> 8) + 0.5) * 2.0**-24
    if precision == "float32":
        return ((word >> 9) + 0.5) * 2.0**-23
    raise ValueError(f"precision must be 'float32' or 'float64', got {precision!r}")


def key_from_seed(seed: int) -> tuple[int, int]:
    """Philox key ``(low, high)`` 32-bit words of a 64-bit seed; raises outside [0, 2^64)."""
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**64:
        raise CounterOverflowError("seed must be an integer in [0, 2**64)")
    return seed & _MASK32, seed >> 32


def child_genealogy_id(parent_id: int, parent_generation: int, birth_order: int) -> int:
    """Genealogy id of the ``birth_order``-th child (1..31) of a particle (decision 0037).

    ``id = parent_id + birth_order * 32**parent_generation``; the primary has generation 0.
    Raises ``CounterOverflowError`` for a 32nd child, a child of a generation-6 particle
    (a 7th generation) or an id of 2^30 or more; callers must then not transport the child.
    """
    if not 0 <= parent_id < MAX_GENEALOGY_ID or parent_generation < 0:
        raise CounterOverflowError("invalid parent genealogy")
    if not 1 <= birth_order <= MAX_CHILDREN:
        raise CounterOverflowError(
            f"birth order {birth_order} outside 1..{MAX_CHILDREN}: genealogy overflow"
        )
    if parent_generation > MAX_PARENT_GENERATION:
        raise CounterOverflowError(
            f"children of generation {parent_generation} exceed the 6-generation bound"
        )
    child = parent_id + birth_order * 32**parent_generation
    if child >= MAX_GENEALOGY_ID:
        raise CounterOverflowError("genealogy id exceeds 2**30")
    return child


def make_counter(
    history: int, genealogy_id: int, block: int, purpose: int
) -> tuple[int, int, int, int]:
    """Philox counter ``(c0, c1, c2, c3)``; raises ``CounterOverflowError`` beyond a bound."""
    if not 0 <= history < MAX_HISTORIES:
        raise CounterOverflowError(f"history index {history} outside [0, 2**32)")
    if not 0 <= genealogy_id < MAX_GENEALOGY_ID:
        raise CounterOverflowError(f"genealogy id {genealogy_id} outside [0, 2**30)")
    if not 0 <= block < 2**32:
        raise CounterOverflowError(f"block index {block} outside [0, 2**32)")
    if purpose not in (PURPOSE_TRANSPORT, PURPOSE_SOURCE, PURPOSE_RESERVED):
        raise CounterOverflowError(f"unknown purpose {purpose}")
    return history, genealogy_id, block, purpose


def draw_block(
    key: tuple[int, int], history: int, genealogy_id: int, block: int, purpose: int
) -> tuple[int, int, int, int]:
    """Four 32-bit words of the stream block ``(history, genealogy_id, block, purpose)``."""
    return philox4x32_10_py(make_counter(history, genealogy_id, block, purpose), key)
