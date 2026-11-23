"""Frozen rules of result schema ``ionmc.result/1.0`` (decision 0045 amendment A-2, tier 1).

``ionmc.results.load`` must validate a schema-1.0 bundle the same way at every later library
version, so it takes the constants, built-in definitions and small numerical rules that the
validation needs from THIS data-only module (plus numpy and the stdlib), never from the live
engine registries. A test compares every value here with the live definition at the commit that
introduced the module: a later change of the engine then fails that test and is handled by a new
schema version (a new ``_result_schema_<major>_<minor>`` module), not by silently re-judging old
bundles.

Contents: the projectile table (z, a, mass), the element table and the built-in material
definitions (the Bragg-additivity effective I), the numerical constants of the engine and the
fixed-point channel rules (``fe_exponent``, ``adaptive_exponent``, the rounding-bound interval
arithmetic), the minimum defined batch count of a ratio, the nuclear tally names, the physical
constant of the dose conversion, and the RNG identity record.

Declared limits (tier 1 still reads these LIVE; each is a binding limit recorded in
``LoadedResult.binding_limits`` or covered by tier 2): the syntax validator of ``requested``
(``ionmc.configio.validate_dict_structure``), the channel-plan helpers ``_generation_range``,
``_species_row`` and ``_spectrum_spec``, the batch reductions ``reduce_batches``/``reduce_ratio``,
``overlap_matrix`` and the immutable stopping-model registry used for the identity and pin checks.
The user stopping object of ``requested.physics.stopping`` is bound field by field where a field
is derivable from the bundle; its table and source content are not in the bundle, so tier 2
reports it as not re-derivable.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

SCATTERING_E_S_MEV = 15.0
MAX_REJECTION_ATTEMPTS = 64
MEV_PER_G_TO_GY = 1.602176634e-10
K_MAX = 40
CAPACITY_EXPONENT = 62
E_QUANTUM_EXPONENT = 30
FE_F_MIN_EXPONENT = -60
FE_F_MAX_EXPONENT = 60
GEN_MAX = 2**31 - 1
CLASS_STEP = 1
CLASS_LOCAL = 2
DEUTERON_NAME = "deuteron"
NUCLEAR_TALLY_NAMES = (
    "nuclear_local",
    "nuclear_alpha_local",
    "nuclear_escaped_neutron",
    "nuclear_escaped_gamma",
    "nuclear_binding",
    "nuclear_imbalance",
)
# Destinations of the nuclear energy balance that are genuinely unsigned (>= 0); nuclear_binding
# and nuclear_imbalance are SIGNED (exothermic reactions give a negative binding), and
# nuclear_alpha_local is a part of nuclear_local, not a destination.
NUCLEAR_UNSIGNED_DESTINATIONS = (
    "nuclear_local",
    "nuclear_escaped_neutron",
    "nuclear_escaped_gamma",
)
NUCLEAR_SIGNED_DESTINATIONS = ("nuclear_binding", "nuclear_imbalance")

U01_MAPPING = {
    "float64": "((w >> 8) + 0.5) * 2**-24",
    "float32": "((w >> 9) + 0.5) * 2**-23",
}


def rng_identity(precision: str) -> dict[str, str]:
    """The ``run.rng`` / ``effective.rng`` record of schema 1.0 (Philox4x32-10)."""
    return {
        "generator": "philox4x32-10",
        "key": "seed low/high 32 bits",
        "counter": "(history, genealogy id, block, purpose)",
        "u01": U01_MAPPING[precision],
    }


PROJECTILES: dict[str, tuple[int, int, float]] = {  # name -> (z, a, nuclear mass [MeV/c2])
    "proton": (1, 1, 938.27208816),
    "deuteron": (1, 2, 1875.61294257),
    "triton": (1, 3, 2808.92113298),
    "helium3": (2, 3, 2808.3915264803227),
    "alpha": (2, 4, 3727.3794066),
    "carbon12": (6, 12, 11174.86323534),
    "oxygen16": (8, 16, 14895.080645240894),
}

ELEMENTS: dict[str, tuple[int, float, float]] = {  # symbol -> (Z, A [g/mol], I [eV])
    "H": (1, 1.00794, 19.2),
    "C": (6, 12.0107, 81.0),
    "N": (7, 14.0067, 82.0),
    "O": (8, 15.9994, 95.0),
    "Na": (11, 22.98977, 149.0),
    "Be": (4, 9.012182, 63.7),
    "Mg": (12, 24.305, 156.0),
    "Al": (13, 26.981538, 166.0),
    "P": (15, 30.97376, 173.0),
    "S": (16, 32.065, 180.0),
    "Cl": (17, 35.453, 174.0),
    "Ar": (18, 39.948, 188.0),
    "K": (19, 39.0983, 190.0),
    "Ca": (20, 40.078, 191.0),
    "Cu": (29, 63.546, 322.0),
    "Pb": (82, 207.2, 823.0),
}

# name -> (density [g/cm3], mass fractions, requested I [eV] or None, sternheimer
# (x0, x1, cbar, a, m) or None, source)
MATERIALS: dict[
    str, tuple[float, dict[str, float], float | None, tuple[float, ...] | None, str]
] = {
    "water": (
        1.0,
        {"H": 0.11189834407236524, "O": 0.8881016559276347},
        78.0,
        (0.24, 2.8004, 3.5017, 0.09116, 3.4773),
        "G4NistMaterialBuilder v11.4.2: G4_WATER; I and Sternheimer parameters as cited in module docstring",  # noqa: E501
    ),
    "air": (
        0.00120479,
        {
            "C": 0.000124000124000124,
            "N": 0.7552677552677554,
            "O": 0.23178123178123178,
            "Ar": 0.012827012827012827,
        },
        85.7,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_AIR (dry air, near sea level)",
    ),
    "pmma": (
        1.19,
        {"C": 0.5998402650050711, "H": 0.08054191635248056, "O": 0.3196178186424483},
        74.0,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_PLEXIGLASS (C5H8O2)",
    ),
    "adipose_tissue_icrp": (
        0.95,
        {"H": 0.114, "C": 0.598, "N": 0.007, "O": 0.278, "Na": 0.001, "S": 0.001, "Cl": 0.001},
        63.2,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_ADIPOSE_TISSUE_ICRP (ICRU 46 / ICRP 23)",
    ),
    "muscle_skeletal_icrp": (
        1.05,
        {
            "H": 0.102,
            "C": 0.143,
            "N": 0.034,
            "O": 0.71,
            "Na": 0.001,
            "P": 0.002,
            "S": 0.003,
            "Cl": 0.001,
            "K": 0.004,
        },
        75.3,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_MUSCLE_SKELETAL_ICRP (ICRU 46 / ICRP 23)",
    ),
    "bone_compact_icru": (
        1.85,
        {
            "H": 0.064,
            "C": 0.278,
            "N": 0.027,
            "O": 0.41,
            "Mg": 0.002,
            "P": 0.07,
            "S": 0.002,
            "Ca": 0.147,
        },
        91.9,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_BONE_COMPACT_ICRU (ICRU 44)",
    ),
    "lung_icrp": (
        1.04,
        {
            "H": 0.105,
            "C": 0.083,
            "N": 0.023,
            "O": 0.779,
            "Na": 0.002,
            "P": 0.001,
            "S": 0.002,
            "Cl": 0.003,
            "K": 0.002,
        },
        75.3,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_LUNG_ICRP (ICRU 46 / ICRP 23)",
    ),
    "beryllium": (
        1.848,
        {"Be": 1.0},
        63.7,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_Be",
    ),
    "aluminium": (
        2.699,
        {"Al": 1.0},
        166.0,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_Al",
    ),
    "copper": (
        8.96,
        {"Cu": 1.0},
        322.0,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_Cu",
    ),
    "lead": (
        11.35,
        {"Pb": 1.0},
        823.0,
        None,
        "G4NistMaterialBuilder v11.4.2: G4_Pb",
    ),
}


def material_definition(record: dict[str, Any]) -> dict[str, Any]:
    """The complete definition of a ``requested.materials`` record (``{"name"}`` of a built-in
    material, else the inline definition) as persisted in ``effective_extra.materials``, with
    the effective mean excitation energy (Bragg additivity when none is requested).
    Raises ``ValueError`` for an invalid record."""
    if set(record) == {"name"}:
        if record["name"] not in MATERIALS:
            raise ValueError(f"{record['name']!r} is not a built-in material")
        rho, fractions, i_ev, st, source = MATERIALS[record["name"]]
        name = record["name"]
    else:
        name, rho = record["name"], float(record["density_g_cm3"])
        fractions = {k: float(v) for k, v in record["mass_fractions"].items()}
        i_ev, source = record.get("I_eV"), record.get("source", "")
        s = record.get("sternheimer")
        st = None if s is None else tuple(float(s[k]) for k in ("x0", "x1", "cbar", "a", "m"))
    if not fractions or any(k not in ELEMENTS for k in fractions):
        raise ValueError(f"{name}: empty or unknown element in the mass fractions")
    if not all(math.isfinite(w) and w > 0.0 for w in fractions.values()):
        raise ValueError(f"{name}: mass fractions must be finite and positive")
    total = math.fsum(fractions.values())
    if abs(total - 1.0) > 1e-6:
        raise ValueError(f"{name}: mass fractions sum to {total}, not 1")
    if not (math.isfinite(rho) and rho > 0.0):
        raise ValueError(f"{name}: density must be finite and positive")
    if i_ev is not None and not (math.isfinite(i_ev) and i_ev > 0.0):
        raise ValueError(f"{name}: I_eV must be finite and positive")
    if st is not None and not all(math.isfinite(x) for x in st):
        raise ValueError(f"{name}: Sternheimer parameters must be finite")
    if i_ev is None:
        z_over_a = sum(w / total * ELEMENTS[s][0] / ELEMENTS[s][1] for s, w in fractions.items())
        ln_i = (
            sum(
                w / total * ELEMENTS[s][0] / ELEMENTS[s][1] * math.log(ELEMENTS[s][2])
                for s, w in fractions.items()
            )
            / z_over_a
        )
    else:
        ln_i = math.log(i_ev)
    return {
        "name": name,
        "density_g_cm3": rho,
        "mass_fractions": {k: float(v) for k, v in sorted(fractions.items())},
        "I_eV_requested": i_ev,
        "I_eV_effective": float(math.exp(ln_i)),
        "sternheimer": None if st is None else dict(zip(("x0", "x1", "cbar", "a", "m"), st,
                                                        strict=True)),
        "source": source,
    }  # fmt: skip


def fe_exponent(f_max: float) -> int:
    """Exponent of the FE quantum ``2^-30 2^ceil(log2 f_max)`` (channels.fe_exponent)."""
    return E_QUANTUM_EXPONENT - math.ceil(math.log2(f_max))


def adaptive_exponent(histories_per_batch: int, bound_per_history: float) -> int:
    """``k = min(40, floor(62 - log2(hpb B_c)))`` with ``hpb B_c 2^k < 2^62``
    (channels.adaptive_exponent); the bound must be positive and finite."""
    if not (bound_per_history > 0.0 and math.isfinite(bound_per_history)):
        raise ValueError(
            f"channel bound per history must be positive and finite: {bound_per_history}"
        )
    total = float(histories_per_batch) * bound_per_history
    k = min(K_MAX, math.floor(CAPACITY_EXPONENT - math.log2(total)))
    while math.ldexp(total, k) >= math.ldexp(1.0, CAPACITY_EXPONENT):
        k -= 1
    return int(k)


def min_defined_batches(n_batches: int) -> int:
    """Smallest ``n_nonzero`` of a defined ratio bin: ``max(2, ceil(B / 2))``."""
    return max(2, math.ceil(n_batches / 2))


def ratio_rounding_bound(
    xhat: np.ndarray[Any, Any],
    a: np.ndarray[Any, Any],
    yhat: np.ndarray[Any, Any],
    b: np.ndarray[Any, Any],
    rhat: np.ndarray[Any, Any],
) -> np.ndarray[Any, Any]:
    """Interval-arithmetic bound of a ratio of quantized sums (simulation.ratio_rounding_bound):
    ``+inf`` where ``yhat - b <= 0``."""
    lower = yhat - b
    with np.errstate(divide="ignore", invalid="ignore"):
        hi = np.abs((xhat + a) / lower - rhat)
        lo = np.abs((xhat - a) / (yhat + b) - rhat)
    return np.where(lower > 0.0, np.maximum(hi, lo), np.inf)
