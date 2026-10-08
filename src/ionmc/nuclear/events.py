"""Non-elastic proton event sampler: one algorithm, a numpy batch path and a scalar twin path.

Decision 0041 section 3 (amendment 2026-10-07). An event is assembled from the shared functions of
``ionmc.physics.nuclear`` (floor + Bernoulli multiplicities, inverse-CDF E', Kalbach-Mann mu,
CM-to-lab boost); the residual only has to exist, it is not given a four-momentum:

1. up to ``MAX_ATTEMPTS`` attempts; attempt ``a`` draws, for the species ``n, p, d, alpha, gamma``
   (in this order, ``SPECIES``), a multiplicity
   ``n_s = floor(lam_s) + [u_s < lam_s - floor(lam_s)]``
   (``multiplicity_round``, cap 16) from the uniform of slot ``s``;
2. the attempt is accepted iff the residual nuclide ``(Z_c - sum z, A_c - sum a)`` has an AME2020
   mass ``M_r`` (``m_res`` finite): no mass or energy test; an exhausted event (64 attempts) is
   reported ``accepted = False``;
3. every product of the accepted attempt (species-major order, running index ``j``) draws
   ``u_bin, u_frac, u_branch, u_mu, u_phi``: the bin of the 64-bin equiprobable table, the
   position in the bin (E' linear in the bin), the (unused, layout-reserving) branch uniform, the
   Kalbach mu uniform and the azimuth; gamma is isotropic (``mu = 2 u_mu - 1``); each product is
   boosted to the laboratory frame;
4. the residual carries the ENDF mean heavy-recoil energy ``T_r = rows.recoil_t_mev`` and is
   deposited locally with the alphas (local deposit = sum of the alpha lab kinetic energies +
   ``T_r``); the event ledger is ``imbalance = T_1 + m_p + M_t - sum E_lab - M_r - T_r`` (signed,
   the tally ``nuclear_imbalance``), so that ``T_1 = sum T_lab + T_r + binding + imbalance``
   with ``binding = sum m_out + M_r - m_p - M_t``.

``sample_events`` is the vectorised numpy implementation (used by the table builder, which needs
10^5-10^6 events per node); ``sample_event_scalar`` evaluates the identical algorithm with the
pure-Python twins of the Warp functions (``ionmc._wpfunc.python_twin(make_nuclear)``) and is the
path the Python reference backend reuses. Both consume uniforms through the same addressing
``(event, attempt, slot)``; ``slot = 0..4`` are the multiplicity uniforms of the species in
``SPECIES`` order and the product with running index ``j`` uses the slots ``8 + 8 j + k`` for
``k = 0..4`` (``u_bin, u_frac, u_branch, u_mu, u_phi``). ``CounterUniforms`` is a stateless
counter-based generator (splitmix64 finalizer), so a batch and the scalar path see identical
uniforms. ``RngUniforms`` wraps a numpy ``Generator`` for plain sampling.
:func:`exact_post_acceptance` enumerates the at most 32 multiplicity outcomes exactly (the builder
solves ``lam`` with it).

Units: MeV (c = 1), angles in radians. All numpy functions are float64 and mirror the shared
functions' operation order; they agree with the twins to rounding of the transcendental
functions (not bitwise).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Protocol

import numpy as np
from numpy.typing import NDArray

from ionmc._wpfunc import python_twin
from ionmc.data.ame import AmeEntry, nuclear_mass_mev
from ionmc.physics.nuclear import (
    DA_MAX,
    DZ_MAX,
    EVENT_ACCEPTED,
    EVENT_BLOCKS_PER_ATTEMPT,
    EVENT_EXHAUSTED,
    EVF_STRIDE,
    EVI_STRIDE,
    KALBACH_A_MAX,
    KALBACH_A_MIN,
    KALBACH_C1,
    KALBACH_C2,
    KALBACH_C3,
    KALBACH_E_T1,
    KALBACH_E_T3,
    KALBACH_I,
    KALBACH_M_A_PROTON,
    KALBACH_M_B,
    KALBACH_NEWTON_ITERATIONS,
    KALBACH_TOLERANCE,
    MAX_ATTEMPTS,
    MAX_PRODUCTS,
    N_BINS,
    N_PARTICLE_UNIFORMS,
    N_SPECIES,
    POISSON_N_MAX,
    PROD_STRIDE,
    SLOT_PARTICLE_BASE,
    SLOT_PARTICLE_STRIDE,
    SPECIES_A,
    SPECIES_Z,
    TCONST_STRIDE,
    kalbach_separation_energy,
    make_nuclear,
)

SPECIES = ("n", "p", "d", "a", "g")
_KALBACH_KEYS = ("n", "p", "d", "a")  # the charged/neutral nucleons; gamma has no Kalbach a
_M_B = np.array([KALBACH_M_B[k] for k in _KALBACH_KEYS] + [0.0], dtype=np.float64)


__all__ = [
    "DA_MAX",
    "DZ_MAX",
    "MAX_ATTEMPTS",
    "N_BINS",
    "N_PARTICLE_UNIFORMS",
    "N_SPECIES",
    "SLOT_PARTICLE_BASE",
    "SLOT_PARTICLE_STRIDE",
    "SPECIES",
    "SPECIES_A",
    "SPECIES_Z",
]  # fmt: skip  (constants shared with ionmc.physics.nuclear, re-exported)


def particle_slot(j: int, k: int) -> int:
    """Uniform slot ``k`` (0..4) of the product with running index ``j``."""
    return SLOT_PARTICLE_BASE + SLOT_PARTICLE_STRIDE * j + k


# ---------------------------------------------------------------------------------------------
# data classes
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class EventModel:
    """Constants of one target nucleus.

    ``species_mass_mev`` (5,) nuclear masses of n, p, d, alpha, gamma (0); ``m_p_mev``,
    ``m_t_mev`` projectile and target masses; ``s_a_mev`` separation energy of the incident
    proton from the compound nucleus; ``s_b_mev`` (5,) those of the ejectiles (gamma 0);
    ``m_res_mev`` (DZ_MAX + 1, DA_MAX + 1) ground-state mass of the residual after the products
    of total charge ``dz`` and mass number ``da`` (``+inf`` where the residual has no AME2020
    mass or does not exist).
    """

    z_t: int
    a_t: int
    m_p_mev: float
    m_t_mev: float
    species_mass_mev: NDArray[np.float64]
    s_a_mev: float
    s_b_mev: NDArray[np.float64]
    m_res_mev: NDArray[np.float64]


@dataclass(frozen=True)
class EnergyRows:
    """Rows of one (target, energy): ``lam`` (5,) multiplicity means (floor + Bernoulli),
    ``edges_mev`` (5, 65) bin edges of the 64 equiprobable bins of E'_CM, ``r`` (5, 64)
    pre-compound fraction per bin, ``recoil_t_mev`` the ENDF mean heavy-recoil energy
    ``sum_r y_r <E_r>`` [MeV] given to the residual."""

    lam: NDArray[np.float64]
    edges_mev: NDArray[np.float64]
    r: NDArray[np.float64]
    recoil_t_mev: float


def interp_rows(
    grid: NDArray[np.float64],
    lam: NDArray[np.float64],
    edges: NDArray[np.float64],
    r: NDArray[np.float64],
    recoil_t: NDArray[np.float64],
    e_mev: float,
) -> EnergyRows:
    """Rows at ``e_mev`` of one target by lin-lin interpolation in E on the table grid (the
    runtime rule; ``lam`` (5, N), ``edges`` (5, N, 65), ``r`` (5, N, 64), ``recoil_t`` (N,)).
    Outside the grid the end rows are used."""
    n = grid.size
    k = int(np.clip(np.searchsorted(grid, e_mev, side="right") - 1, 0, n - 2))
    t = min(max((e_mev - grid[k]) / (grid[k + 1] - grid[k]), 0.0), 1.0)
    return EnergyRows(
        (1.0 - t) * lam[:, k] + t * lam[:, k + 1],
        (1.0 - t) * edges[:, k] + t * edges[:, k + 1],
        (1.0 - t) * r[:, k] + t * r[:, k + 1],
        float((1.0 - t) * recoil_t[k] + t * recoil_t[k + 1]),
    )


def build_event_model(ame: dict[tuple[int, int], AmeEntry], z_t: int, a_t: int) -> EventModel:
    """The total break-up ``(Z_r, A_r) = (0, 0)`` (for example p + C-12 -> p + 3 alpha) is an
    existing residual of mass 0 (orchestrator decision 2026-10-07); every other residual needs an
    AME2020 mass. Masses, Kalbach separation energies and the residual-mass table of the target
    ``(z_t, a_t)``; masses from AME2020 (``nuclear_mass_mev``), separation energies from the Kalbach
    (1988) systematics formula (not from AME: the systematics were fitted with that formula)."""
    z_c, a_c = z_t + 1, a_t + 1
    masses = np.array(
        [
            nuclear_mass_mev(ame, 0, 1),
            nuclear_mass_mev(ame, 1, 1),
            nuclear_mass_mev(ame, 1, 2),
            nuclear_mass_mev(ame, 2, 4),
            0.0,
        ]
    )
    s_a = kalbach_separation_energy(z_c, a_c, z_t, a_t, KALBACH_I["p"])
    s_b = np.zeros(N_SPECIES)
    for i, key in enumerate(_KALBACH_KEYS):
        s_b[i] = kalbach_separation_energy(
            z_c, a_c, z_c - SPECIES_Z[i], a_c - SPECIES_A[i], KALBACH_I[key]
        )
    m_res = np.full((DZ_MAX + 1, DA_MAX + 1), np.inf)
    for dz in range(DZ_MAX + 1):
        for da in range(DA_MAX + 1):
            z_r, a_r = z_c - dz, a_c - da
            if z_r >= 0 and a_r >= 1 and (z_r, a_r) in ame:
                m_res[dz, da] = nuclear_mass_mev(ame, z_r, a_r)
            elif z_r == 0 and a_r == 0:
                m_res[dz, da] = 0.0
    return EventModel(
        z_t=z_t,
        a_t=a_t,
        m_p_mev=float(masses[1]),
        m_t_mev=nuclear_mass_mev(ame, z_t, a_t),
        species_mass_mev=masses,
        s_a_mev=float(s_a),
        s_b_mev=s_b,
        m_res_mev=m_res,
    )


def event_constants(model: EventModel) -> NDArray[np.float64]:
    """The ``TCONST_STRIDE`` row of ``model`` for ``sample_event``: ``m_p, m_t, s_a, a_t``, the five
    species masses, the five separation energies, ``z_t``, 0."""
    row = np.zeros(TCONST_STRIDE)
    row[0:4] = (model.m_p_mev, model.m_t_mev, model.s_a_mev, float(model.a_t))
    row[4:9] = model.species_mass_mev
    row[9:14] = model.s_b_mev
    row[14] = float(model.z_t)
    return row


def reachable_residuals(model: EventModel) -> NDArray[np.bool_]:
    """Boolean ``(DZ_MAX + 1, DA_MAX + 1)``: ``(dz, da)`` reachable by a multiplicity vector with
    every ``n_s <= 16`` (n, p, d, alpha)."""
    out = np.zeros(model.m_res_mev.shape, dtype=bool)
    n = np.arange(POISSON_N_MAX + 1)
    nn, npr, nd, na = np.meshgrid(n, n, n, n, indexing="ij")
    out[(npr + nd + 2 * na).ravel(), (nn + npr + 2 * nd + 4 * na).ravel()] = True
    return out


def reaction_q_values(model: EventModel) -> NDArray[np.float64]:
    """Q [MeV] of every reachable channel ``(dz, da)`` with the residual in its ground state and
    the lightest product assignment of that ``(dz, da)`` irrelevant (Q depends on the product
    masses, so the maximum over the assignments is returned; ``-inf`` where unreachable or the
    residual has no AME mass)."""
    n = np.arange(POISSON_N_MAX + 1)
    nn, npr, nd, na = (x.ravel() for x in np.meshgrid(n, n, n, n, indexing="ij"))
    m = model.species_mass_mev
    m_out = nn * m[0] + npr * m[1] + nd * m[2] + na * m[3]
    dz = npr + nd + 2 * na
    da = nn + npr + 2 * nd + 4 * na
    q = model.m_p_mev + model.m_t_mev - m_out - model.m_res_mev[dz, da]
    q = np.where(np.isfinite(model.m_res_mev[dz, da]), q, -np.inf)
    out = np.full(model.m_res_mev.shape, -np.inf)
    np.maximum.at(out, (dz, da), q)
    return out


# ---------------------------------------------------------------------------------------------
# uniform sources
# ---------------------------------------------------------------------------------------------
class UniformSource(Protocol):
    """Provider of uniforms addressed by ``(event, attempt, slot)``."""

    def uniform(
        self, events: NDArray[np.int64], attempt: int, slot: NDArray[np.int64] | int
    ) -> NDArray[np.float64]:
        """Uniform in [0, 1) for every entry of ``events`` (and ``slot`` when it is an array)."""
        ...


_M64 = (1 << 64) - 1
_GOLDEN = 0x9E3779B97F4A7C15
_K_ATTEMPT = 0xC2B2AE3D27D4EB4F
_K_SLOT = 0x165667B19E3779F9
_MIX1 = 0xBF58476D1CE4E5B9
_MIX2 = 0x94D049BB133111EB
_TWO53 = float(1 << 53)


def _splitmix_scalar(x: int) -> int:
    x &= _M64
    x = ((x ^ (x >> 30)) * _MIX1) & _M64
    x = ((x ^ (x >> 27)) * _MIX2) & _M64
    return x ^ (x >> 31)


class CounterUniforms:
    """Stateless counter-based uniforms: ``u = (splitmix64(key + c_e e + c_a a + c_s s) >> 11 +
    0.5) / 2^53`` in (0, 1). Identical in the numpy and scalar forms; the same ``key`` always
    gives the same uniforms for the same address (common random numbers)."""

    def __init__(self, key: int) -> None:
        self.key = int(key) & _M64

    def uniform(
        self, events: NDArray[np.int64], attempt: int, slot: NDArray[np.int64] | int
    ) -> NDArray[np.float64]:
        ev = np.asarray(events).astype(np.uint64)
        s = np.asarray(slot).astype(np.uint64)
        with np.errstate(over="ignore"):
            x = (
                np.uint64(self.key)
                + ev * np.uint64(_GOLDEN)
                + np.uint64((attempt * _K_ATTEMPT) & _M64)
                + s * np.uint64(_K_SLOT)
            )
            x = (x ^ (x >> np.uint64(30))) * np.uint64(_MIX1)
            x = (x ^ (x >> np.uint64(27))) * np.uint64(_MIX2)
            x = x ^ (x >> np.uint64(31))
        out: NDArray[np.float64] = ((x >> np.uint64(11)).astype(np.float64) + 0.5) / _TWO53
        return out

    def scalar(self, event: int, attempt: int, slot: int) -> float:
        """The same uniform for one address, in pure Python integers."""
        x = self.key + event * _GOLDEN + attempt * _K_ATTEMPT + slot * _K_SLOT
        return ((_splitmix_scalar(x) >> 11) + 0.5) / _TWO53

    def for_event(self, event: int) -> Callable[[int, int], float]:
        """``uni(attempt, slot)`` of one event for :func:`sample_event_scalar`."""
        return lambda attempt, slot: self.scalar(event, attempt, slot)


class RngUniforms:
    """Uniforms from a numpy ``Generator`` (the addresses are ignored: plain sampling)."""

    def __init__(self, rng: np.random.Generator) -> None:
        self.rng = rng

    def uniform(
        self, events: NDArray[np.int64], attempt: int, slot: NDArray[np.int64] | int
    ) -> NDArray[np.float64]:
        return self.rng.random(np.asarray(events).shape)


# ---------------------------------------------------------------------------------------------
# numpy versions of the shared functions (same operation order)
# ---------------------------------------------------------------------------------------------
def multiplicity_round_np(
    u: NDArray[np.float64], lam: float | NDArray[np.float64]
) -> NDArray[np.int64]:
    """Vectorised ``multiplicity_round`` (``n_max = POISSON_N_MAX``)."""
    lam_a = np.broadcast_to(np.asarray(lam, dtype=np.float64), u.shape)
    base = np.floor(lam_a)
    n = base.astype(np.int64) + (u < lam_a - base)
    return np.minimum(n, POISSON_N_MAX).astype(np.int64)


_OUTCOME_BITS = (np.arange(1 << N_SPECIES)[:, None] >> np.arange(N_SPECIES)) & 1


def exact_post_acceptance(
    model: EventModel, lam: NDArray[np.float64]
) -> tuple[float, NDArray[np.float64]]:
    """``(p_accept, mean_counts)`` of one attempt by exact enumeration of the 2^5 outcomes of the
    floor + Bernoulli multiplicities: ``p_accept`` the probability that the residual exists, and
    ``mean_counts`` (5,) the mean multiplicities of an accepted attempt (the post-acceptance
    yields; zeros if ``p_accept = 0``)."""
    lam_a = np.asarray(lam, dtype=np.float64)
    base = np.floor(lam_a)
    frac = lam_a - base
    bits = _OUTCOME_BITS
    prob = np.prod(np.where(bits == 1, frac, 1.0 - frac), axis=1)
    counts = np.minimum(base.astype(np.int64) + bits, POISSON_N_MAX)
    dz = counts @ np.asarray(SPECIES_Z)
    da = counts @ np.asarray(SPECIES_A)
    inside = (dz <= DZ_MAX) & (da <= DA_MAX)
    big_m = np.where(
        inside, model.m_res_mev[np.minimum(dz, DZ_MAX), np.minimum(da, DA_MAX)], np.inf
    )
    ok = np.isfinite(big_m) & (prob > 0.0)
    p_acc = float(prob[ok].sum())
    if p_acc <= 0.0:
        return 0.0, np.zeros(N_SPECIES)
    mean = (prob[ok, None] * counts[ok]).sum(axis=0) / p_acc
    return p_acc, np.asarray(mean, dtype=np.float64)


def kalbach_a_np(
    e_a: float, e_b: NDArray[np.float64], m_b: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Vectorised ``kalbach_a`` (``e_a`` the incident-channel energy incl. separation energy)."""
    x1 = np.minimum(e_a, KALBACH_E_T1) * e_b / e_a
    x3 = np.minimum(e_a, KALBACH_E_T3) * e_b / e_a
    a: NDArray[np.float64] = (
        KALBACH_C1 * x1
        + KALBACH_C2 * x1 * x1 * x1
        + KALBACH_C3 * KALBACH_M_A_PROTON * m_b * x3 * x3 * x3 * x3
    )
    return a


def _kalbach_cdf_np(
    mu: NDArray[np.float64], a: NDArray[np.float64], r: NDArray[np.float64]
) -> NDArray[np.float64]:
    p = 0.5 * a * (1.0 + mu)
    q = 0.5 * a * (1.0 - mu)
    return np.sinh(p) * (np.cosh(q) - r * np.sinh(q)) / np.sinh(a)


def _kalbach_pdf_np(
    mu: NDArray[np.float64], a: NDArray[np.float64], r: NDArray[np.float64]
) -> NDArray[np.float64]:
    return a * (np.cosh(a * mu) + r * np.sinh(a * mu)) / (2.0 * np.sinh(a))


def kalbach_mu_np(
    u: NDArray[np.float64], a_in: NDArray[np.float64], r: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Vectorised ``kalbach_mu`` (float64): the same guarded Newton iteration with bisection
    fallback, 40 iterations at most, tolerance ``KALBACH_TOLERANCE['float64']``; converged
    entries are dropped from the working set."""
    tol = KALBACH_TOLERANCE["float64"]
    a = np.minimum(np.maximum(a_in, 0.0), KALBACH_A_MAX)
    mu = 2.0 * u - 1.0
    idx = np.nonzero(a >= KALBACH_A_MIN)[0]
    if idx.size:
        ua, aa, rr = u[idx], a[idx], r[idx]
        s = (2.0 * ua - 1.0) * np.sinh(aa)
        t = np.abs(s)
        g = np.log(t + np.sqrt(t * t + 1.0)) / aa
        g = np.where(s < 0.0, -g, g)
        m = np.minimum(np.maximum(g, -1.0), 1.0)
        lo = np.full(idx.size, -1.0)
        hi = np.full(idx.size, 1.0)
        work = np.arange(idx.size)
        for _ in range(KALBACH_NEWTON_ITERATIONS):
            if work.size == 0:
                break
            mw, aw, rw, uw = m[work], aa[work], rr[work], ua[work]
            diff = _kalbach_cdf_np(mw, aw, rw) - uw
            above = diff > 0.0
            hi_w = np.where(above, mw, hi[work])
            lo_w = np.where(above, lo[work], mw)
            pdf = _kalbach_pdf_np(mw, aw, rw)
            nxt = 0.5 * (lo_w + hi_w)
            with np.errstate(divide="ignore", invalid="ignore"):
                cand = mw - diff / pdf
            use = (pdf > 0.0) & (cand >= lo_w) & (cand <= hi_w)
            nxt = np.where(use, cand, nxt)
            converged = np.abs(nxt - mw) <= tol
            hi[work], lo[work], m[work] = hi_w, lo_w, nxt
            work = work[~converged]
        mu[idx] = m
    return np.minimum(np.maximum(mu, -1.0), 1.0)


def cm_boost_np(t_mev: float, m_p: float, m_t: float) -> tuple[float, float, float]:
    """``cm_boost``: ``(beta, gamma, sqrt_s)`` of the p + target centre of mass."""
    e_lab = t_mev + m_p
    p_lab = math.sqrt(t_mev * (t_mev + 2.0 * m_p))
    e_tot = e_lab + m_t
    s = m_p * m_p + m_t * m_t + 2.0 * e_lab * m_t
    sqrt_s = math.sqrt(s)
    return p_lab / e_tot, e_tot / sqrt_s, sqrt_s


# ---------------------------------------------------------------------------------------------
# numpy batch sampler
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class EventBatch:
    """Result of :func:`sample_events` for ``n`` events (arrays of length ``n``; NaN / -1 for an
    event that was not accepted).

    ``accepted``; ``attempts`` (``MAX_ATTEMPTS`` for an exhausted event); ``counts`` (n, 5)
    multiplicities of the accepted attempt; ``z_r``, ``a_r`` the residual nuclide;
    ``recoil_t_mev`` the residual kinetic energy ``T_r``; ``imbalance_mev`` the signed ledger
    ``T_1 + m_p + M_t - sum E_lab - M_r - T_r`` (tally ``nuclear_imbalance``); ``binding_mev``
    ``sum m_out + M_r - m_p - M_t`` (tally ``nuclear_binding``); ``local_deposit_mev`` the sum of
    the alpha lab kinetic energies + ``T_r``; diagnostics ``imbalance_cm_mev`` (``sqrt(s) - sum
    E_cm - M_r - T_r``), ``p_cm_sum_mev`` (|sum of the product CM momenta|), ``sum_e_prime_mev``
    (n, 5) the sum of the sampled E'_CM per species. With ``want_particles``: the accepted events'
    products ``particle_event`` (index into the batch), ``particle_species``, ``e_cm_mev`` (E'),
    ``mu``, ``phi`` and the lab ``(e_lab, px, py, pz)`` [MeV] in ``particle_lab`` (P, 4)."""

    accepted: NDArray[np.bool_]
    attempts: NDArray[np.int64]
    counts: NDArray[np.int64]
    z_r: NDArray[np.int64]
    a_r: NDArray[np.int64]
    recoil_t_mev: NDArray[np.float64]
    imbalance_mev: NDArray[np.float64]
    binding_mev: NDArray[np.float64]
    local_deposit_mev: NDArray[np.float64]
    imbalance_cm_mev: NDArray[np.float64]
    p_cm_sum_mev: NDArray[np.float64]
    sum_e_prime_mev: NDArray[np.float64]
    particle_event: NDArray[np.int64] | None = None
    particle_species: NDArray[np.int64] | None = None
    e_cm_mev: NDArray[np.float64] | None = None
    mu: NDArray[np.float64] | None = None
    phi: NDArray[np.float64] | None = None
    particle_lab: NDArray[np.float64] | None = None


def sample_events(
    model: EventModel,
    rows: EnergyRows,
    t_lab_mev: float,
    n_events: int,
    source: UniformSource,
    *,
    want_particles: bool = False,
    event_offset: int = 0,
) -> EventBatch:
    """Sample ``n_events`` non-elastic events of the proton of kinetic energy ``t_lab_mev``."""
    beta, gamma, sqrt_s = cm_boost_np(t_lab_mev, model.m_p_mev, model.m_t_mev)
    e_a = t_lab_mev * model.a_t / (model.a_t + 1.0) + model.s_a_mev
    masses = model.species_mass_mev
    sp_z = np.asarray(SPECIES_Z, dtype=np.int64)
    sp_a = np.asarray(SPECIES_A, dtype=np.int64)
    z_c, a_c = model.z_t + 1, model.a_t + 1
    nan = np.full(n_events, np.nan)
    accepted = np.zeros(n_events, dtype=bool)
    attempts = np.zeros(n_events, dtype=np.int64)
    counts_out = np.zeros((n_events, N_SPECIES), dtype=np.int64)
    z_r = np.full(n_events, -1, dtype=np.int64)
    a_r = np.full(n_events, -1, dtype=np.int64)
    t_r = nan.copy()
    imb, bind, local, imb_cm, p_sum = (nan.copy() for _ in range(5))
    sum_e_prime = np.full((n_events, N_SPECIES), np.nan)
    keep: dict[str, list[Any]] = {k: [] for k in ("ev", "sp", "ecm", "mu", "phi", "lab")}
    active = np.arange(n_events, dtype=np.int64)
    for attempt in range(1, MAX_ATTEMPTS + 1):
        if active.size == 0:
            break
        ev_id = active + event_offset
        counts = np.empty((active.size, N_SPECIES), dtype=np.int64)
        for s in range(N_SPECIES):
            counts[:, s] = multiplicity_round_np(source.uniform(ev_id, attempt, s), rows.lam[s])
        dz = counts @ sp_z
        da = counts @ sp_a
        inside = (dz <= DZ_MAX) & (da <= DA_MAX)
        big_m_all = np.where(
            inside, model.m_res_mev[np.minimum(dz, DZ_MAX), np.minimum(da, DA_MAX)], np.inf
        )
        ok = np.isfinite(big_m_all)
        attempts[active] = attempt
        if np.any(ok):
            acc = active[ok]
            counts_a = counts[ok]
            big_m = big_m_all[ok]
            m = acc.size
            ev_a = acc + event_offset
            tot = counts_a.sum(axis=1)
            n_part = int(tot.sum())
            sp_flat = np.repeat(np.tile(np.arange(N_SPECIES), m), counts_a.ravel())
            ev_loc = np.repeat(np.arange(m), tot)
            start = np.cumsum(tot) - tot
            j = np.arange(n_part) - np.repeat(start, tot)
            ev_part = ev_a[ev_loc]
            u = [
                source.uniform(ev_part, attempt, SLOT_PARTICLE_BASE + SLOT_PARTICLE_STRIDE * j + k)
                for k in range(N_PARTICLE_UNIFORMS)
            ]
            k_bin = np.minimum(np.maximum(np.floor(u[0] * N_BINS).astype(np.int64), 0), N_BINS - 1)
            e_lo = rows.edges_mev[sp_flat, k_bin]
            e_hi = rows.edges_mev[sp_flat, k_bin + 1]
            e_p = e_lo + u[1] * (e_hi - e_lo)
            r_p = rows.r[sp_flat, k_bin]
            a_k = kalbach_a_np(e_a, e_p + model.s_b_mev[sp_flat], _M_B[sp_flat])
            mu = kalbach_mu_np(u[3], a_k, r_p)
            mu = np.where(sp_flat == 4, 2.0 * u[3] - 1.0, mu)
            phi = 2.0 * math.pi * u[4]
            m_s = masses[sp_flat]
            p_cm = np.sqrt(e_p * (e_p + 2.0 * m_s))
            st = np.sqrt(np.maximum(1.0 - mu * mu, 0.0))
            px = p_cm * st * np.cos(phi)
            py = p_cm * st * np.sin(phi)
            pz = p_cm * mu
            e_cm_tot = e_p + m_s
            e_lab = gamma * (e_cm_tot + beta * pz)
            pz_lab = gamma * (pz + beta * e_cm_tot)
            sum_e_lab = np.bincount(ev_loc, weights=e_lab, minlength=m)
            sum_e_cm = np.bincount(ev_loc, weights=e_cm_tot, minlength=m)
            sum_p = np.sqrt(
                sum(np.bincount(ev_loc, weights=c, minlength=m) ** 2 for c in (px, py, pz))
            )
            alpha = sp_flat == 3
            alpha_t = np.bincount(ev_loc[alpha], weights=(e_lab - m_s)[alpha], minlength=m)
            sum_ep = np.bincount(
                ev_loc * N_SPECIES + sp_flat, weights=e_p, minlength=m * N_SPECIES
            ).reshape(m, N_SPECIES)
            m_out = counts_a @ masses
            recoil = rows.recoil_t_mev
            accepted[acc] = True
            counts_out[acc] = counts_a
            z_r[acc] = z_c - counts_a @ sp_z
            a_r[acc] = a_c - counts_a @ sp_a
            t_r[acc] = recoil
            imb[acc] = t_lab_mev + model.m_p_mev + model.m_t_mev - sum_e_lab - big_m - recoil
            bind[acc] = m_out + big_m - model.m_p_mev - model.m_t_mev
            local[acc] = alpha_t + recoil
            imb_cm[acc] = sqrt_s - sum_e_cm - big_m - recoil
            p_sum[acc] = sum_p
            sum_e_prime[acc] = sum_ep
            if want_particles:
                keep["ev"].append(acc[ev_loc])
                keep["sp"].append(sp_flat)
                keep["ecm"].append(e_p)
                keep["mu"].append(mu)
                keep["phi"].append(phi)
                keep["lab"].append(np.stack([e_lab, px, py, pz_lab], axis=1))
        active = active[~ok]
    common: dict[str, Any] = {
        "accepted": accepted,
        "attempts": attempts,
        "counts": counts_out,
        "z_r": z_r,
        "a_r": a_r,
        "recoil_t_mev": t_r,
        "imbalance_mev": imb,
        "binding_mev": bind,
        "local_deposit_mev": local,
        "imbalance_cm_mev": imb_cm,
        "p_cm_sum_mev": p_sum,
        "sum_e_prime_mev": sum_e_prime,
    }
    if not want_particles:
        return EventBatch(**common)
    cat = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in keep.items()}
    lab = np.concatenate(keep["lab"]) if keep["lab"] else np.zeros((0, 4))
    order = np.argsort(cat["ev"], kind="stable")
    return EventBatch(
        **common,
        particle_event=cat["ev"].astype(np.int64)[order],
        particle_species=cat["sp"].astype(np.int64)[order],
        e_cm_mev=cat["ecm"][order],
        mu=cat["mu"][order],
        phi=cat["phi"][order],
        particle_lab=lab[order],
    )


# ---------------------------------------------------------------------------------------------
# scalar twin path
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ScalarEvent:
    """Result of :func:`sample_event_scalar`: ``accepted``, ``attempts``, ``counts`` (5 ints),
    the residual ``z_r``, ``a_r``, ``recoil_t_mev``, the ledger ``imbalance_mev`` and
    ``binding_mev``, ``local_deposit_mev`` and the products ``particles`` as ``(species, e_cm,
    mu, phi, e_lab, px, py, pz_lab)`` tuples (species index into ``SPECIES``; lab four-momentum in
    MeV) in running order (fields as in :class:`EventBatch`)."""

    accepted: bool
    attempts: int
    counts: tuple[int, ...]
    z_r: int
    a_r: int
    recoil_t_mev: float
    imbalance_mev: float
    binding_mev: float
    local_deposit_mev: float
    particles: tuple[tuple[float, ...], ...]


def sample_event_scalar(
    model: EventModel,
    rows: EnergyRows,
    t_lab_mev: float,
    uni: Callable[[int, int], float],
    nu: Any = None,
) -> ScalarEvent:
    """One event with the pure-Python twins of the shared functions (the algorithm of the module
    docstring). ``uni(attempt, slot)`` returns the uniform of an address; ``nu`` defaults to
    ``python_twin(make_nuclear)``."""
    nu = python_twin(make_nuclear) if nu is None else nu
    t = float(t_lab_mev)
    # rows already interpolated at t: a two-node grid [t, t + 1] gives the weight 0, so that the
    # shared sampler reproduces the rows exactly ((1 - 0) x + 0 x = x)
    two = np.array([t, t + 1.0])
    lam = np.repeat(np.asarray(rows.lam, dtype=np.float64), 2)
    edges = np.repeat(np.asarray(rows.edges_mev, dtype=np.float64)[:, None, :], 2, axis=1).ravel()
    rpre = np.repeat(np.asarray(rows.r, dtype=np.float64)[:, None, :], 2, axis=1).ravel()
    recoil = np.full(2, float(rows.recoil_t_mev))
    evi = np.zeros(EVI_STRIDE, dtype=np.int64)
    evf = np.zeros(EVF_STRIDE)
    prod = np.zeros(MAX_PRODUCTS * PROD_STRIDE)

    def key(_h: int, _gid: int, blk: int, word: int) -> float:
        return uni(blk // EVENT_BLOCKS_PER_ATTEMPT + 1, (blk % EVENT_BLOCKS_PER_ATTEMPT) * 4 + word)

    status = int(
        nu.sample_event(
            0, 0, 0, key, 0, t, two, 2, lam, edges, rpre, recoil,
            event_constants(model), np.ravel(model.m_res_mev), evi, evf, prod, 0,
        )
    )  # fmt: skip
    if status == EVENT_EXHAUSTED:
        return ScalarEvent(
            False,
            MAX_ATTEMPTS,
            (0,) * N_SPECIES,
            -1,
            -1,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            (),
        )
    if status != EVENT_ACCEPTED:
        raise RuntimeError(
            f"nuclear event with more than {MAX_PRODUCTS} products (status {status})"
        )
    n_prod = int(evi[8])
    return ScalarEvent(
        True,
        int(evi[0]),
        tuple(int(x) for x in evi[1:6]),
        int(evi[6]),
        int(evi[7]),
        float(evf[0]),
        float(evf[1]),
        float(evf[2]),
        float(evf[3]),
        tuple(
            tuple(float(x) for x in prod[j * PROD_STRIDE : (j + 1) * PROD_STRIDE])
            for j in range(n_prod)
        ),
    )
