"""One hadronic elastic event on recorded inputs (decision 0041 slice C, V3-005C C3).

:func:`sample_elastic_event` is the whole event decision of the elastic channel, free of transport
state: given the material's elastic rows, the proton energy, the three uniforms of the event (the
target uniform of the candidate block and the CM-cosine and azimuth uniforms of the event block,
see :mod:`ionmc.physics.elastic`) and the incident direction, it chooses the target, samples the CM
cosine from the table, applies the exact two-body kinematics and returns the outgoing states. The
Python reference calls it from its history loop; the recorded-input parity hook (row P5-ext) replays
recorded inputs through it (and C4 will replay them on the Warp twins of the same
:mod:`ionmc.physics.elastic` functions).

Outputs: ``primary`` is the particle that continues as the primary (for p-p the faster proton);
``other`` is the slower proton of p-p, to be transported as a secondary, or the p + A recoil
nucleus, deposited locally (``recoil_t_mev``). Energy is conserved exactly per event:
``t1 = t_primary + t_other``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from ionmc.nuclear.tables import interp_row

__all__ = ["ElasticEvent", "cum_fraction_at", "sample_elastic_event"]

TWO_PI = 2.0 * math.pi


@dataclass(frozen=True)
class ElasticEvent:
    target: int  # elastic-table target index (0 = H-1, a p-p event)
    pp: bool
    mu_cm: float
    phi: float
    t_primary_mev: float
    dir_primary: tuple[float, float, float]
    t_other_mev: float  # p-p: kinetic energy of the slower proton; p + A: recoil energy
    dir_other: tuple[float, float, float]


def cum_fraction_at(rows: Any, e_mev: float) -> np.ndarray:
    """Cumulative target fractions of a ``MaterialElastic`` at ``e_mev`` (lin-lin cumulative partial
    Sigma over the lin-lin total, as ``MaterialNuclear.cum_fraction_at``)."""
    tot = rows.sigma_at(e_mev)
    cum = np.array([interp_row(rows.grid_e_mev, c, e_mev) for c in rows.cum_sigma_mass_cm2_g])
    return cum / tot if tot > 0.0 else np.zeros_like(cum)


def _frame_dir(
    px: float, pz: float, phi: float, e1: Any, e2: Any, d: tuple[float, float, float]
) -> tuple[float, float, float]:
    pm = math.sqrt(px * px + pz * pz)
    if pm <= 0.0:
        return d
    cx, cy, cz = px * math.cos(phi) / pm, px * math.sin(phi) / pm, pz / pm
    return (
        cx * float(e1[0]) + cy * float(e2[0]) + cz * d[0],
        cx * float(e1[1]) + cy * float(e2[1]) + cz * d[1],
        cx * float(e1[2]) + cy * float(e2[2]) + cz * d[2],
    )


def sample_elastic_event(
    EL: Any,
    NU: Any,
    arrays: dict[str, np.ndarray],
    rows: Any,
    m_proton_mev: float,
    t1_mev: float,
    u_target: float,
    u_mu: float,
    u_phi: float,
    direction: tuple[float, float, float],
    e1: Any,
    e2: Any,
) -> ElasticEvent:
    """The event (module docstring). ``EL`` / ``NU`` are the python twins of
    :func:`ionmc.physics.elastic.make_elastic` / :func:`ionmc.physics.nuclear.make_nuclear`,
    ``arrays`` the table arrays, ``(e1, e2)`` an orthonormal basis perpendicular to
    ``direction``."""
    cum = cum_fraction_at(rows, t1_mev)
    n_t = len(rows.target_index)
    chosen = -1
    for k in range(n_t):
        chosen = int(NU.select_target(u_target, float(cum[k]), k, chosen, int(k == n_t - 1)))
    tgt = int(rows.target_index[chosen])
    grid = arrays["grid_e_mev"]
    k = int(NU.grid_locate(t1_mev, grid, grid.size))
    t = min(max((t1_mev - grid[k]) / (grid[k + 1] - grid[k]), 0.0), 1.0)
    edges = arrays["edges_mu"]
    n_g, n_e = edges.shape[1], edges.shape[2]
    flat = edges.reshape(-1)
    mu = float(
        EL.sample_mu_edges(u_mu, flat, (tgt * n_g + k) * n_e, (tgt * n_g + k + 1) * n_e, t, n_e - 1)
    )
    phi = TWO_PI * u_phi
    m2 = float(arrays["target_mass_mev"][tgt])
    ta, pxa, pza, tb, pxb, pzb = (
        float(x) for x in EL.elastic_two_body(t1_mev, m_proton_mev, m2, mu)
    )
    pp = tgt == 0
    d = (float(direction[0]), float(direction[1]), float(direction[2]))
    swap = pp and int(EL.faster_is_recoil(ta, tb)) == 1  # identical particles: faster continues
    if swap:
        ta, pxa, pza, tb, pxb, pzb = tb, pxb, pzb, ta, pxa, pza
    return ElasticEvent(
        tgt, pp, mu, phi, ta, _frame_dir(pxa, pza, phi, e1, e2, d), tb,
        _frame_dir(pxb, pzb, phi, e1, e2, d),
    )  # fmt: skip
