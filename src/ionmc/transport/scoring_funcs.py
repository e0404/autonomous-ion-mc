# mypy: ignore-errors
# (Warp function sources use runtime precision types in annotations; see _wpfunc.py.)
"""Shared scoring functions of the channel hook (decision 0040, sections 1, 2 and 7).

``make_scoring_funcs(real)`` returns precision-generic Warp functions for kernels; the reference
backend uses their pure-Python twins (``ionmc._wpfunc.python_twin``: the same source text without
Warp). Everything scientifically relevant about a scoring piece is here and nowhere else; array
access, quantization and atomics stay glue in the backends.

A transport step with path length ``s_act`` and CSDA mean loss ``dE_m`` has, at the midpoint
energy ``E_mid = E0 - dE_m / 2``, the water stopping power ``S_mid`` and the log-log slope
``gamma = d ln S / d ln E`` of its table interpolant. The water stopping power is taken as a linear
ramp along the step, anchored at the path midpoint::

    k    = dS/dt = -gamma S_mid dE_m / (E_mid s_act)        (MeV/mm^2)
    Edot = dE_m / s_act                                       (MeV/mm)

A piece ``p`` of length ``l`` at ``tau = t_a + l/2 - s_act/2`` (distance of its midpoint from the
step midpoint) has ``S_p = S_mid + k tau``, ``E_p = E_mid - Edot tau`` and the exact moments of a
linear ramp, ``m1 = l S_p`` and ``m2 = l (S_p^2 + k^2 l^2 / 12)``; the sum of ``m1`` over the pieces
of a step is ``s_act S_mid`` because ``sum l tau = 0``. ``channel_value`` is the per-piece
increment of each channel kind (codes ``0 E, 1 L, 2 LS, 3 LS2, 4 ES, 5 FE, 6 FL, 7 N``, the order
of ``ionmc.transport.channels.KINDS``). ``lookup_bin`` and ``spectrum_bin`` locate a value on a
uniform (linear or logarithmic) axis; a spectrum has an underflow bin 0 and an overflow bin
``n_bins + 1`` with bins closed at the lower and open at the upper edge.

Units: energies MeV, lengths mm, stopping power MeV/mm (numerically keV/um).
"""

import functools
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func

wp.set_module_options({"enable_backward": False})

KIND_CODES = {"E": 0, "L": 1, "LS": 2, "LS2": 3, "ES": 4, "FE": 5, "FL": 6, "N": 7}
"""Integer code of every channel kind (the order of ``channels.KINDS``)."""
TINY_X = 1.0e-30


@functools.cache
def make_scoring_funcs(real: type) -> SimpleNamespace:
    """Return the shared scoring functions for precision ``real``."""
    name = check_real(real)
    tiny = wp.constant(real(TINY_X))
    twelve = wp.constant(real(12.0))

    @named_func(name)
    def loglog_slope(ly0: real, ly1: real, inv_dl: real) -> real:
        """Log-log slope ``(ln S_{i+1} - ln S_i) / dln E`` of a table bin: the exact derivative
        ``d ln S / d ln E`` of the interpolant (``inv_dl`` is the inverse step in ``ln E``)."""
        return (ly1 - ly0) * inv_dl

    @named_func(name)
    def let_ramp_slope(s_mid: real, gamma: real, de_mean: real, e_mid: real, s_act: real) -> real:
        """Slope ``k = dS/dt = -gamma S_mid dE_m / (E_mid s_act)`` [MeV/mm^2] of the water
        stopping-power ramp along a step; 0 for a step without length or energy."""
        k = real(0.0)
        if s_act > real(0.0) and e_mid > real(0.0):
            k = -gamma * s_mid * de_mean / (e_mid * s_act)
        return k

    @named_func(name)
    def piece_state(s_mid: real, k: real, e_mid: real, e_dot: real, tau: real) -> tuple[real, real]:
        """Water stopping power ``S_mid + k tau`` and energy ``E_mid - Edot tau`` at the
        midpoint of a piece (``tau``: its distance from the step midpoint along the path)."""
        return s_mid + k * tau, e_mid - e_dot * tau

    @named_func(name)
    def piece_moments(s_bar: real, k: real, length: real) -> tuple[real, real]:
        """``(l S, l (S^2 + k^2 l^2 / 12))``: the exact integrals of ``S`` and ``S^2`` over a
        piece of length ``l`` for the linear ramp with midpoint value ``s_bar`` and slope ``k``."""
        return length * s_bar, length * (s_bar * s_bar + k * k * length * length / twelve)

    @named_func(name)
    def channel_value(
        kind: int, eps: real, length: real, m1: real, m2: real, s_bar: real, f: real
    ) -> real:
        """Increment of a channel kind for one piece: E ``eps``, L ``l``, LS ``m1``, LS2 ``m2``,
        ES ``eps S``, FE ``eps f``, FL ``l f`` and N ``1``; an unknown kind gives 0."""
        out = real(0.0)
        if kind == 0:
            out = eps
        elif kind == 1:
            out = length
        elif kind == 2:
            out = m1
        elif kind == 3:
            out = m2
        elif kind == 4:
            out = eps * s_bar
        elif kind == 5:
            out = eps * f
        elif kind == 6:
            out = length * f
        elif kind == 7:
            out = real(1.0)
        return out

    @named_func(name)
    def lookup_bin(x: real, a0: real, inv_da: real, n: int, log_axis: int) -> tuple[int, real, int]:
        """Bin ``i`` in [0, n-2], fraction ``f`` in [0, 1] and ``in_domain`` (1 if ``x`` lies on
        the axis) of ``x`` on a uniform axis of ``n`` points: ``t = (x - a0) inv_da`` (linear) or
        ``(ln x - a0) inv_da`` (``log_axis`` = 1, ``a0`` is then ``ln`` of the first point). The
        index and fraction are clamped like ``log_bin_index`` (no extrapolation)."""
        v = x
        if log_axis != 0:
            v = wp.log(wp.max(x, tiny))
        t = (v - a0) * inv_da
        inside = 0
        if t >= real(0.0) and t <= real(n - 1):
            inside = 1
        i = int(wp.floor(t))
        i = wp.max(wp.min(i, n - 2), 0)
        f = wp.min(wp.max(t - real(i), real(0.0)), real(1.0))
        return i, f, inside

    @named_func(name)
    def spectrum_bin(x: real, a0: real, inv_da: real, n_bins: int, log_axis: int) -> int:
        """Bin of ``x`` in a spectrum of ``n_bins`` uniform bins: 0 underflow (also for NaN),
        ``1 + floor(t)`` for ``0 <= t < n_bins`` (lower edge included) and ``n_bins + 1``
        (overflow) for ``t >= n_bins``."""
        v = x
        if log_axis != 0:
            v = wp.log(wp.max(x, tiny))
        t = (v - a0) * inv_da
        b = 0
        if t >= real(0.0):
            if t >= real(n_bins):
                b = n_bins + 1
            else:
                b = int(wp.floor(t)) + 1
        return b

    return SimpleNamespace(
        loglog_slope=loglog_slope,
        let_ramp_slope=let_ramp_slope,
        piece_state=piece_state,
        piece_moments=piece_moments,
        channel_value=channel_value,
        lookup_bin=lookup_bin,
        spectrum_bin=spectrum_bin,
        real=name,
    )
