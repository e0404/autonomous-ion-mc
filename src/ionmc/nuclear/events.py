"""Non-elastic proton event sampler: one algorithm, a numpy batch path and a scalar twin path.

Decision 0041 section 3. An event is assembled from the shared functions of
``ionmc.physics.nuclear`` (Poisson multiplicities, inverse-CDF E', Kalbach-Mann mu, residual
four-momentum by difference, AME2020 residual-mass test, at most 64 attempts, CM-to-lab boost):

1. up to ``MAX_ATTEMPTS`` attempts; attempt ``a`` draws, for the species ``n, p, d, alpha, gamma``
   (in this order, ``SPECIES``), a Poisson multiplicity ``n_s`` of mean ``lam_s`` (cap 16);
2. every product (species-major order, running index ``j``) draws ``u_bin, u_frac, u_branch,
   u_mu, u_phi``: the bin of the 64-bin equiprobable table, the position in the bin (E' linear
   in the bin), the (unused, layout-reserving) branch uniform, the Kalbach mu uniform and the
   azimuth; gamma is isotropic (``mu = 2 u_mu - 1``);
3. the residual four-momentum is ``(sqrt(s), 0) - sum`` of the product CM four-momenta; the
   attempt is accepted iff the residual nuclide ``(Z_c - sum z, A_c - sum a)`` has an AME2020
   mass ``M_r`` and its invariant mass ``m_r >= M_r`` (``E* = m_r - M_r``);
4. an accepted event is boosted to the laboratory frame.

``sample_events`` is the vectorised numpy implementation (used by the table builder, which needs
10^5-10^6 events per node); ``sample_event_scalar`` evaluates the identical algorithm with the
pure-Python twins of the Warp functions (``ionmc._wpfunc.python_twin(make_nuclear)``) and is the
path the Python reference backend reuses. Both consume uniforms through the same addressing
``(event, attempt, slot)``; ``slot = 0..4`` are the multiplicity uniforms of the species in
``SPECIES`` order and the product with running index ``j`` uses the slots ``8 + 8 j + k`` for
``k = 0..4`` (``u_bin, u_frac, u_branch, u_mu, u_phi``). ``CounterUniforms`` is a stateless
counter-based generator (splitmix64 finalizer), so a batch and the scalar path see identical
uniforms and repeated batches with the same key are common random numbers (the lambda fixed
point of the builder needs this). ``RngUniforms`` wraps a numpy ``Generator`` for plain sampling.

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
    POISSON_N_MAX,
    kalbach_separation_energy,
    make_nuclear,
)

SPECIES = ("n", "p", "d", "a", "g")
"""Sampled product species: neutron, proton, deuteron, alpha, gamma (decision 0041 section 3)."""
N_SPECIES = 5
N_BINS = 64
MAX_ATTEMPTS = 64
SPECIES_Z = (0, 1, 1, 2, 0)
SPECIES_A = (1, 1, 2, 4, 0)
DZ_MAX = POISSON_N_MAX * (1 + 1 + 2)
"""Largest total charge of the products of the multiplicity bound (16 p + 16 d + 16 alpha)."""
DA_MAX = POISSON_N_MAX * (1 + 1 + 2 + 4)
"""Largest total mass number of the products (16 of each of n, p, d, alpha)."""
_KALBACH_KEYS = ("n", "p", "d", "a")  # the charged/neutral nucleons; gamma has no Kalbach a
_M_B = np.array([KALBACH_M_B[k] for k in _KALBACH_KEYS] + [0.0], dtype=np.float64)
SLOT_PARTICLE_BASE = 8
SLOT_PARTICLE_STRIDE = 8
N_PARTICLE_UNIFORMS = 5


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
    """Rows of one (target, energy): ``lam`` (5,) Poisson means, ``edges_mev`` (5, 65) bin edges
    of the 64 equiprobable bins of E'_CM, ``r`` (5, 64) pre-compound fraction per bin."""

    lam: NDArray[np.float64]
    edges_mev: NDArray[np.float64]
    r: NDArray[np.float64]


def build_event_model(ame: dict[tuple[int, int], AmeEntry], z_t: int, a_t: int) -> EventModel:
    """Masses, Kalbach separation energies and the residual-mass table of the target ``(z_t,
    a_t)``; masses from AME2020 (``nuclear_mass_mev``), separation energies from the Kalbach
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
def poisson_inverse_np(
    u: NDArray[np.float64], lam: float | NDArray[np.float64]
) -> NDArray[np.int64]:
    """Vectorised ``poisson_inverse`` (``n_max = POISSON_N_MAX``)."""
    lam_a = np.broadcast_to(np.asarray(lam, dtype=np.float64), u.shape)
    p = np.exp(-lam_a)
    cdf = p.copy()
    n = np.zeros(u.shape, dtype=np.int64)
    for _ in range(16):
        step = (u > cdf) & (n < POISSON_N_MAX)
        n = n + step
        p = np.where(step, p * lam_a / np.maximum(n, 1), p)
        cdf = np.where(step, cdf + p, cdf)
    return n


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
    """Result of :func:`sample_events` for ``n`` events.

    ``accepted`` (n,) bool; ``attempts`` (n,) number of attempts used (``MAX_ATTEMPTS`` for an
    exhausted event); ``counts`` (n, 5) multiplicities of the accepted attempt (0 for exhausted
    events); ``e_star_mev`` (n,) residual excitation (NaN if not accepted). With
    ``want_particles``: arrays of the accepted events' products ``particle_event`` (index into the
    batch), ``particle_species``, ``e_cm_mev`` (E'), ``mu``, ``phi``, lab ``(e_lab, px, py, pz)``
    [MeV] in ``particle_lab`` (P, 4), and the residual lab four-momentum ``residual_lab`` (n, 4)
    (NaN if not accepted)."""

    accepted: NDArray[np.bool_]
    attempts: NDArray[np.int64]
    counts: NDArray[np.int64]
    e_star_mev: NDArray[np.float64]
    particle_event: NDArray[np.int64] | None = None
    particle_species: NDArray[np.int64] | None = None
    e_cm_mev: NDArray[np.float64] | None = None
    mu: NDArray[np.float64] | None = None
    phi: NDArray[np.float64] | None = None
    particle_lab: NDArray[np.float64] | None = None
    residual_lab: NDArray[np.float64] | None = None
    m_r_mev: NDArray[np.float64] | None = None


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
    accepted = np.zeros(n_events, dtype=bool)
    attempts = np.zeros(n_events, dtype=np.int64)
    counts_out = np.zeros((n_events, N_SPECIES), dtype=np.int64)
    e_star = np.full(n_events, np.nan)
    keep: dict[str, list[Any]] = {k: [] for k in ("ev", "sp", "ecm", "mu", "phi", "lab")}
    res_lab = np.full((n_events, 4), np.nan) if want_particles else None
    m_r_out = np.full(n_events, np.nan) if want_particles else None
    active = np.arange(n_events, dtype=np.int64)
    for attempt in range(1, MAX_ATTEMPTS + 1):
        if active.size == 0:
            break
        m = active.size
        ev_id = active + event_offset
        counts = np.empty((m, N_SPECIES), dtype=np.int64)
        for s in range(N_SPECIES):
            counts[:, s] = poisson_inverse_np(source.uniform(ev_id, attempt, s), rows.lam[s])
        tot = counts.sum(axis=1)
        n_part = int(tot.sum())
        sp_flat = np.repeat(np.tile(np.arange(N_SPECIES), m), counts.ravel())
        ev_loc = np.repeat(np.arange(m), tot)
        start = np.cumsum(tot) - tot
        j = np.arange(n_part) - np.repeat(start, tot)
        ev_part = ev_id[ev_loc]
        u = [
            source.uniform(ev_part, attempt, SLOT_PARTICLE_BASE + SLOT_PARTICLE_STRIDE * j + k)
            for k in range(N_PARTICLE_UNIFORMS)
        ]
        k_bin = np.minimum(np.maximum(np.floor(u[0] * N_BINS).astype(np.int64), 0), N_BINS - 1)
        e_lo = rows.edges_mev[sp_flat, k_bin]
        e_hi = rows.edges_mev[sp_flat, k_bin + 1]
        e_p = e_lo + u[1] * (e_hi - e_lo)
        r_p = rows.r[sp_flat, k_bin]
        is_g = sp_flat == 4
        a_k = kalbach_a_np(e_a, e_p + model.s_b_mev[sp_flat], _M_B[sp_flat])
        mu = kalbach_mu_np(u[3], a_k, r_p)
        mu = np.where(is_g, 2.0 * u[3] - 1.0, mu)
        phi = 2.0 * math.pi * u[4]
        m_s = masses[sp_flat]
        p_cm = np.sqrt(e_p * (e_p + 2.0 * m_s))
        st = np.sqrt(np.maximum(1.0 - mu * mu, 0.0))
        px = p_cm * st * np.cos(phi)
        py = p_cm * st * np.sin(phi)
        pz = p_cm * mu
        e_cm_tot = e_p + m_s
        sum_e = np.bincount(ev_loc, weights=e_cm_tot, minlength=m)
        sum_px = np.bincount(ev_loc, weights=px, minlength=m)
        sum_py = np.bincount(ev_loc, weights=py, minlength=m)
        sum_pz = np.bincount(ev_loc, weights=pz, minlength=m)
        e_r = sqrt_s - sum_e
        pr = -sum_px, -sum_py, -sum_pz
        p_res = np.sqrt(pr[0] * pr[0] + pr[1] * pr[1] + pr[2] * pr[2])
        m2 = (e_r - p_res) * (e_r + p_res)
        m_r = np.sqrt(np.abs(m2))
        m_r = np.where(m2 < 0.0, -m_r, m_r)
        dz = (counts * sp_z).sum(axis=1)
        da = (counts * sp_a).sum(axis=1)
        inside = (dz <= DZ_MAX) & (da <= DA_MAX)
        big_m = np.where(
            inside, model.m_res_mev[np.minimum(dz, DZ_MAX), np.minimum(da, DA_MAX)], np.inf
        )
        ok = m_r >= big_m
        attempts[active] = attempt
        if np.any(ok):
            acc = active[ok]
            accepted[acc] = True
            counts_out[acc] = counts[ok]
            e_star[acc] = m_r[ok] - big_m[ok]
            if want_particles:
                assert res_lab is not None and m_r_out is not None
                sel = ok[ev_loc]
                e_lab = gamma * (e_cm_tot + beta * pz)
                pz_lab = gamma * (pz + beta * e_cm_tot)
                keep["ev"].append(active[ev_loc[sel]])
                keep["sp"].append(sp_flat[sel])
                keep["ecm"].append(e_p[sel])
                keep["mu"].append(mu[sel])
                keep["phi"].append(phi[sel])
                keep["lab"].append(np.stack([e_lab[sel], px[sel], py[sel], pz_lab[sel]], axis=1))
                res_lab[acc, 0] = gamma * (e_r[ok] + beta * pr[2][ok])
                res_lab[acc, 1] = pr[0][ok]
                res_lab[acc, 2] = pr[1][ok]
                res_lab[acc, 3] = gamma * (pr[2][ok] + beta * e_r[ok])
                m_r_out[acc] = m_r[ok]
        active = active[~ok]
    if not want_particles:
        return EventBatch(accepted, attempts, counts_out, e_star)
    cat = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in keep.items()}
    lab = np.concatenate(keep["lab"]) if keep["lab"] else np.zeros((0, 4))
    order = np.argsort(cat["ev"], kind="stable")
    return EventBatch(
        accepted,
        attempts,
        counts_out,
        e_star,
        particle_event=cat["ev"].astype(np.int64)[order],
        particle_species=cat["sp"].astype(np.int64)[order],
        e_cm_mev=cat["ecm"][order],
        mu=cat["mu"][order],
        phi=cat["phi"][order],
        particle_lab=lab[order],
        residual_lab=res_lab,
        m_r_mev=m_r_out,
    )


# ---------------------------------------------------------------------------------------------
# scalar twin path
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ScalarEvent:
    """Result of :func:`sample_event_scalar`: ``accepted``, ``attempts``, ``counts`` (5 ints),
    ``e_star_mev``, ``m_r_mev``, products ``particles`` as ``(species, e_cm, mu, phi, e_lab, px,
    py, pz_lab)`` tuples (species index into ``SPECIES``; lab four-momentum in MeV) in running
    order, and the residual lab four-momentum ``residual_lab`` ``(E, px, py, pz)``."""

    accepted: bool
    attempts: int
    counts: tuple[int, ...]
    e_star_mev: float
    m_r_mev: float
    particles: tuple[tuple[float, ...], ...]
    residual_lab: tuple[float, float, float, float]


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
    beta, gamma, sqrt_s = nu.cm_boost(t_lab_mev, model.m_p_mev, model.m_t_mev)
    e_a = t_lab_mev * model.a_t / (model.a_t + 1.0) + model.s_a_mev
    masses = [float(x) for x in model.species_mass_mev]
    for attempt in range(1, MAX_ATTEMPTS + 1):
        n = [
            nu.poisson_inverse(uni(attempt, s), float(rows.lam[s]), POISSON_N_MAX) for s in range(5)
        ]
        parts: list[tuple[int, float, float, float]] = []
        j = 0
        for s in range(N_SPECIES):
            for _ in range(n[s]):
                u = [uni(attempt, particle_slot(j, k)) for k in range(N_PARTICLE_UNIFORMS)]
                kb = nu.inv_cdf_bin(u[0], N_BINS)
                e_p = nu.inv_cdf_sample(
                    u[1], float(rows.edges_mev[s, kb]), float(rows.edges_mev[s, kb + 1])
                )
                if s < 4:
                    a = nu.kalbach_a(e_a, e_p + float(model.s_b_mev[s]), float(_M_B[s]))
                    mu = nu.kalbach_mu(u[3], a, float(rows.r[s, kb]))
                else:
                    mu = 2.0 * u[3] - 1.0
                parts.append((s, e_p, mu, 2.0 * math.pi * u[4]))
                j += 1
        dz = sum(SPECIES_Z[p[0]] for p in parts)
        da = sum(SPECIES_A[p[0]] for p in parts)
        big_m = float(model.m_res_mev[dz, da]) if dz <= DZ_MAX and da <= DA_MAX else math.inf
        if math.isinf(big_m):
            continue
        cm = []
        for s, e_p, mu, phi in parts:
            p = math.sqrt(e_p * (e_p + 2.0 * masses[s]))
            st = math.sqrt(max(1.0 - mu * mu, 0.0))
            cm.append((e_p + masses[s], p * st * math.cos(phi), p * st * math.sin(phi), p * mu))
        e_r = sqrt_s - sum(c[0] for c in cm)
        pr = [-sum(c[i] for c in cm) for i in (1, 2, 3)]
        m_r = nu.residual_invariant_mass(e_r, pr[0], pr[1], pr[2])
        if not nu.residual_mass_ok(m_r, big_m):
            continue
        out = []
        for (s, e_p, mu, phi), _c in zip(parts, cm, strict=True):
            lab = nu.cm_to_lab(e_p, mu, phi, masses[s], beta, gamma)
            out.append((float(s), e_p, mu, phi, lab[0], lab[1], lab[2], lab[3]))
        e_r_lab, pz_r_lab = nu.boost_z(e_r, pr[2], beta, gamma)
        return ScalarEvent(
            True,
            attempt,
            tuple(int(x) for x in n),
            m_r - big_m,
            m_r,
            tuple(out),
            (e_r_lab, pr[0], pr[1], pz_r_lab),
        )
    return ScalarEvent(
        False, MAX_ATTEMPTS, (0,) * N_SPECIES, math.nan, math.nan, (), (math.nan,) * 4
    )
