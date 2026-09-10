"""Host-side generation of secondary protons from nonelastic reactions (0013).

When a primary proton is removed by a nonelastic nuclear reaction (decision
`0012`), a fraction ``f_p`` of its residual kinetic energy is carried by
secondary protons that travel on and deposit a broad low-level dose. This module
samples those secondary protons from the reaction records emitted by the
transport drivers and returns a batch (a :class:`~ionmc.transport.state.ParticleState`)
ready to be transported by the *same* proton engine, plus the local and escaping
energy adjustments needed to keep the energy budget closed.

The sampling is pure host code (numpy / Python), identical for the reference and
Warp backends: it consumes only the per-reaction records (vertex depth, residual
energy, weight), which DEV-007 already guarantees are the same on both backends,
and it draws from the counter-based RNG (:mod:`ionmc.rng`) keyed by a stable
reaction ordinal, so the secondary set is bit-reproducible across backends
(decision `0001`).

Model (decision `0013`): mean multiplicity ``nu_p(E) = 0.5 + 0.004 E`` (Poisson);
a two-component energy spectrum (evaporation Maxwellian ``E exp(-E/T)``, ``T = 2``
MeV, mixed with a forward cascade component uniform on ``[10 MeV, E]``); energies
renormalised so their per-reaction sum equals ``f_p w E`` exactly (closing the
budget); forward emission along ``+z``; secondaries below a transport cut are
deposited locally at the vertex.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ionmc.rng import RandomState, rand_init

if TYPE_CHECKING:
    from ionmc.transport.state import ParticleState

#: Per-reaction energy partition (decision 0013): local heavy fragments,
#: transported secondary protons; the remainder (neutrons, gammas, binding) is
#: the escaping fraction ``1 - F_HEAVY - F_SECONDARY``.
F_HEAVY: float = 0.12
F_SECONDARY: float = 0.50

#: Mean secondary-proton multiplicity ``nu_p(E) = NU_A + NU_B E`` [protons].
NU_A: float = 0.5
NU_B: float = 0.004

#: Energy-spectrum parameters: evaporation nuclear temperature [MeV], the
#: evaporation vs cascade mixing probability, and the cascade lower edge [MeV].
EVAP_TEMPERATURE_MEV: float = 2.0
P_EVAPORATION: float = 0.5
CASCADE_MIN_MEV: float = 10.0

#: Below this kinetic energy [MeV] a secondary proton is deposited locally at the
#: vertex instead of being transported (short-range evaporation tail).
SECONDARY_TRANSPORT_CUT_MEV: float = 1.0

#: RNG stream offsets so secondary sampling and secondary transport use distinct
#: counter streams that do not collide with the primary streams.
_SAMPLING_STREAM_BASE: int = 0x5EC0_0000
_TRANSPORT_STREAM_BASE: int = 0x5EC1_0000


@dataclass
class SecondaryBatch:
    """Secondary protons to transport plus the energy adjustments they imply.

    ``state`` is ``None`` when no secondary is above the transport cut.
    ``local_deposit`` lists ``(z_mm, weighted_energy_mev)`` deposits made at the
    reaction vertices (sub-cut secondaries). ``escaped_reduction_mev`` is the
    total weighted energy moved *out* of the escaping channel into the secondary
    (transported + local) channel; the caller subtracts it from the pass-1
    escaping accumulator.
    """

    state: ParticleState | None
    local_deposit: list[tuple[float, float]]
    escaped_reduction_mev: float
    n_secondaries: int


def mean_multiplicity(energy_mev: float) -> float:
    """Mean number of transportable secondary protons ``nu_p(E)`` (>= 0)."""
    return max(NU_A + NU_B * energy_mev, 0.0)


def _sample_poisson(rng: RandomState, mean: float) -> int:
    """Knuth's Poisson sampler using the counter RNG (small means here)."""
    if mean <= 0.0:
        return 0
    target = np.exp(-mean)
    k = 0
    product = 1.0
    while True:
        product *= rng.randf()
        if product <= target:
            return k
        k += 1


def _sample_raw_energy(rng: RandomState, energy_mev: float) -> float:
    """One raw secondary energy [MeV] from the evaporation+cascade mixture."""
    if rng.randf() < P_EVAPORATION or energy_mev <= CASCADE_MIN_MEV:
        # evaporation: E ~ Gamma(shape=2, scale=T) = -T (ln u1 + ln u2)
        u1 = rng.randf()
        u2 = rng.randf()
        e = -EVAP_TEMPERATURE_MEV * (np.log(max(u1, 1e-30)) + np.log(max(u2, 1e-30)))
        return min(e, energy_mev)
    # cascade: forward knockout, uniform on [CASCADE_MIN, E]
    return CASCADE_MIN_MEV + rng.randf() * (energy_mev - CASCADE_MIN_MEV)


def generate_secondaries(
    react_z_mm: np.ndarray,
    react_energy_mev: np.ndarray,
    react_weight: np.ndarray,
    seed: int,
    f_heavy: float = F_HEAVY,
    f_secondary: float = F_SECONDARY,
    transport_cut_mev: float = SECONDARY_TRANSPORT_CUT_MEV,
) -> SecondaryBatch:
    """Generate secondary protons from the per-reaction records.

    Parameters
    ----------
    react_z_mm, react_energy_mev, react_weight:
        Per-reacting-history vertex depth [mm], primary residual kinetic energy
        at the vertex [MeV], and history weight. Only entries with
        ``react_energy_mev > 0`` are reactions; the arrays are indexed by the
        primary history index, which is the stable reaction ordinal.
    seed:
        Base seed; secondary sampling and transport streams are derived from it.
    f_heavy, f_secondary:
        The local-fragment and secondary-proton energy fractions (decision 0013).

    Returns
    -------
    SecondaryBatch
        The transportable secondary state (or ``None``), the sub-cut local
        deposits, and the energy to remove from the escaping channel.
    """
    from ionmc.transport.state import ParticleState, Species, Status

    positions: list[float] = []
    energies: list[float] = []
    weights: list[float] = []
    rng_states: list[int] = []
    local_deposit: list[tuple[float, float]] = []
    escaped_reduction = 0.0
    global_index = 0

    n = int(react_energy_mev.shape[0])
    for h in range(n):
        e_primary = float(react_energy_mev[h])
        if e_primary <= 0.0:
            continue  # this history did not react
        z = float(react_z_mm[h])
        w = float(react_weight[h])
        pool = f_secondary * e_primary  # unweighted secondary energy budget [MeV]
        sampler = RandomState.from_state(rand_init(seed ^ _SAMPLING_STREAM_BASE, h))
        count = _sample_poisson(sampler, mean_multiplicity(e_primary))
        if count == 0:
            continue  # no secondary produced: the pool stays in the escaping channel
        raw = np.array([_sample_raw_energy(sampler, e_primary) for _ in range(count)])
        total_raw = float(raw.sum())
        if total_raw <= 0.0:
            continue
        scaled = raw * (pool / total_raw)  # renormalise: sum == f_secondary * E
        # all of this reaction's secondary pool leaves the escaping channel
        escaped_reduction += w * float(scaled.sum())
        for e_sec in scaled:
            e_sec = float(e_sec)
            if e_sec < transport_cut_mev:
                local_deposit.append((z, w * e_sec))  # short-range: deposit local
                continue
            positions.append(z)
            energies.append(e_sec)
            weights.append(w)
            rng_states.append(rand_init(seed ^ _TRANSPORT_STREAM_BASE, global_index))
            global_index += 1

    n_transported = len(energies)
    if n_transported == 0:
        return SecondaryBatch(None, local_deposit, escaped_reduction, 0)

    state = ParticleState.allocate(n_transported)
    state.position_mm[:, 2] = np.asarray(positions, dtype=np.float64)
    state.direction[:] = np.array([0.0, 0.0, 1.0])
    state.energy_mev[:] = np.asarray(energies, dtype=np.float64)
    state.weight[:] = np.asarray(weights, dtype=np.float64)
    state.species[:] = int(Species.PROTON)
    state.status[:] = int(Status.ALIVE)
    state.rng_state[:] = np.asarray(rng_states, dtype=np.uint32)
    return SecondaryBatch(state, local_deposit, escaped_reduction, n_transported)
