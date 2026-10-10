# mypy: ignore-errors
# (Warp function sources use runtime precision types in annotations; see _wpfunc.py.)
"""Hadronic elastic proton scattering: shared sampling and kinematics functions (decision 0041
slice C, V3-005C step C3).

Written as ``@wp.func`` factories like :mod:`ionmc.physics.nuclear`: the Python reference backend
calls the pure-Python twins (``ionmc._wpfunc.python_twin``), the Warp kernels of C4 will call the
same source. Scalar arithmetic only, so that the twin and the kernel agree.

``elastic_two_body`` is the exact relativistic two-body kinematics of
:func:`ionmc.nuclear.elastic_kin.two_body` (that function is the independent numpy oracle used by
the closure tests, row P6 (7)): a projectile of mass ``m1`` and kinetic energy ``t`` on a target of
mass ``m2`` at rest scatters to the CM cosine ``mu``; the two lab four-momenta are boosted from
their own CM momenta, so that the sum can be compared with the incident one. For ``m1 = m2`` (p-p)
the CM cosine of "particle 1" is sampled from the table and the faster of the two protons is the
one that continues as the primary (:func:`faster_is_recoil`).

``sample_mu_edges`` is the inverse-CDF draw of the CM cosine from the quantile table
(``n_q`` equiprobable bins, quantile interpolation linear in E between the two bracketing grid
nodes, piecewise-uniform density inside a bin), the reference being
:meth:`ionmc.nuclear.elastic_tables.ElasticTable.sample_mu`.

Random-number layout of an elastic event (PURPOSE_NUCLEAR stream of the particle): the candidate
block ``nc`` supplies ``u0`` (acceptance), ``u1`` (the redrawn optical depth), ``u2`` (target) and
``u3`` (the channel: elastic iff ``u3 < Sigma_el / Sigma_tot``); an accepted elastic event consumes
exactly one further block ``nc + 1``: slot 0 is the CM-cosine uniform, slot 1 the azimuth uniform,
and the block counter advances by one (the particle goes on). Non-elastic events use the blocks
from ``nc + 1`` as before (the particle ends there), so the non-elastic layout is unchanged.
"""

from __future__ import annotations

import functools
from types import SimpleNamespace

import warp as wp

from ionmc._wpfunc import check_real, named_func

__all__ = ["make_elastic"]


@functools.cache
def make_elastic(real: type) -> SimpleNamespace:
    """Return the shared elastic-scattering functions for precision ``real``."""
    name = check_real(real)

    @named_func(name)
    def sample_mu_edges(
        u: real,
        edges: wp.array(dtype=real),
        base_lo: int,
        base_hi: int,
        t: real,
        n_q: int,
    ) -> real:
        """CM cosine from the quantile table. ``edges[base_lo + i]`` and ``edges[base_hi + i]``
        (``i = 0 .. n_q``) are the ``n_q + 1`` quantile edges at the two grid nodes that bracket
        the energy, ``t`` in [0, 1] the lin-in-E weight of the upper node. ``x = u n_q``,
        ``i = min(floor(x), n_q - 1)``, ``mu = e_i + (x - i)(e_{i+1} - e_i)`` with
        ``e_i = (1 - t) edges_lo[i] + t edges_hi[i]``."""
        x = u * real(n_q)
        i = int(wp.floor(x))
        i = wp.min(i, n_q - 1)
        lo0 = (real(1.0) - t) * edges[base_lo + i] + t * edges[base_hi + i]
        lo1 = (real(1.0) - t) * edges[base_lo + i + 1] + t * edges[base_hi + i + 1]
        return lo0 + (x - real(i)) * (lo1 - lo0)

    @named_func(name)
    def elastic_two_body(
        t_mev: real, m1: real, m2: real, mu: real
    ) -> tuple[real, real, real, real, real, real]:
        """Lab kinetic energies and momentum components of both particles of the elastic
        scattering to the CM cosine ``mu`` (azimuth 0, scattering plane x-z, beam along +z):
        ``(t1, px1, pz1, t2, px2, pz2)`` [MeV]. Particle 1 is the projectile (CM cosine ``mu``),
        particle 2 the target recoil; each is boosted from its own CM momentum (closure of the sum
        with the incident four-momentum is the test, row P6 (7))."""
        e1 = t_mev + m1
        p = wp.sqrt(t_mev * (t_mev + real(2.0) * m1))
        rs = wp.sqrt(m1 * m1 + m2 * m2 + real(2.0) * e1 * m2)
        gam = (e1 + m2) / rs
        bgam = p / rs
        pcm = p * m2 / rs
        sin_t = wp.sqrt(wp.max(real(1.0) - mu * mu, real(0.0)))
        ecm1 = wp.sqrt(pcm * pcm + m1 * m1)
        ecm2 = wp.sqrt(pcm * pcm + m2 * m2)
        e1p = gam * ecm1 + bgam * pcm * mu
        pz1 = gam * pcm * mu + bgam * ecm1
        px1 = pcm * sin_t
        e2p = gam * ecm2 - bgam * pcm * mu
        pz2 = -gam * pcm * mu + bgam * ecm2
        px2 = -pcm * sin_t
        return e1p - m1, px1, pz1, e2p - m2, px2, pz2

    @named_func(name)
    def faster_is_recoil(t1: real, t2: real) -> int:
        """1 iff the second particle is strictly faster (identical-particle rule of p-p: the
        faster proton continues as the primary), else 0."""
        r = int(0)
        if t2 > t1:
            r = 1
        return r

    @named_func(name)
    def channel_is_elastic(u_channel: real, sigma_el: real, sigma_tot: real) -> int:
        """Channel choice of an accepted candidate: elastic iff ``u < Sigma_el / Sigma_tot``
        written as ``u Sigma_tot < Sigma_el`` (no division); 0 for a non-positive total."""
        r = int(0)
        if sigma_tot > real(0.0):
            if u_channel * sigma_tot < sigma_el:
                r = 1
        return r

    return SimpleNamespace(
        sample_mu_edges=sample_mu_edges,
        elastic_two_body=elastic_two_body,
        faster_is_recoil=faster_is_recoil,
        channel_is_elastic=channel_is_elastic,
        real=name,
    )
