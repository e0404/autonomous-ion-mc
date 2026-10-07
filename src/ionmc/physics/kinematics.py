# mypy: ignore-errors
# (Warp function sources use runtime precision types in annotations; see _wpfunc.py.)
"""Relativistic kinematics as precision-generic Warp functions (decision 0039).

``make_kinematics(R)`` returns functions for ``R`` = ``wp.float32`` or ``wp.float64``.
Every function is for use in kernels; the Python reference backend uses its pure-Python twin
(``ionmc._wpfunc.python_twin(make_kinematics)``, same source text, no Warp). Units:
kinetic energy ``t_mev`` and rest energy ``m_mev`` in MeV, momentum times velocity ``pv`` in
MeV (``p c * beta``).
"""

import functools
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func
from ionmc.physics.projectiles import ELECTRON_MASS_MEV

wp.set_module_options({"enable_backward": False})


@functools.cache
def make_kinematics(real: type) -> SimpleNamespace:
    """Return ``pv_mev, beta2, gamma, tmax_mev`` for precision ``real``."""
    name = check_real(real)
    me = wp.constant(real(ELECTRON_MASS_MEV))

    @named_func(name)
    def pv_mev(t_mev: real, m_mev: real) -> real:
        """p v = T (T + 2 M) / (T + M) [MeV]."""
        return t_mev * (t_mev + real(2.0) * m_mev) / (t_mev + m_mev)

    @named_func(name)
    def beta2(t_mev: real, m_mev: real) -> real:
        """beta^2 = T (T + 2 M) / (T + M)^2."""
        tm = t_mev + m_mev
        return t_mev * (t_mev + real(2.0) * m_mev) / (tm * tm)

    @named_func(name)
    def gamma(t_mev: real, m_mev: real) -> real:
        """Lorentz factor 1 + T / M."""
        return real(1.0) + t_mev / m_mev

    @named_func(name)
    def tmax_mev(t_mev: real, m_mev: real) -> real:
        """Maximum energy transfer to a free electron [MeV], exact two-body kinematics."""
        g = real(1.0) + t_mev / m_mev
        bg2 = t_mev * (t_mev + real(2.0) * m_mev) / (m_mev * m_mev)
        ratio = me / m_mev
        return real(2.0) * me * bg2 / (real(1.0) + real(2.0) * g * ratio + ratio * ratio)

    return SimpleNamespace(pv_mev=pv_mev, beta2=beta2, gamma=gamma, tmax_mev=tmax_mev, real=name)
