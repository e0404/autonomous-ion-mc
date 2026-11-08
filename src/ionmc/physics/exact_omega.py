# mypy: ignore-errors
# (Warp function sources use runtime precision types in annotations; see _wpfunc.py.)
"""Exact-Omega mapping of the qualified proton ``f_dM`` to light ions (decision 0043 section 3).

The model (decision 0043 section 3, "Selected (r5)")::

    T_ion(pv, p1v1; material) = f_dM(pv_eq, p1v1_eq) (z E_s / pv)^2 / X_S
    p1v1_eq / pv_eq = p1v1 / pv                 (the memory term 1 - (pv/p1v1)^2 is unchanged)
    pv_eq is the root of  b_p(pv_eq, p1v1_eq; material) = b_ion(pv, p1v1; z, M, material)
    b = ln Omega + delta_b_Fano                 (Moliere collision number of the slab traversed
                                                 from p1v1 down to pv along the CSDA path)

with, for a projectile of charge ``z`` and mass ``M`` (kinetic energy ``T`` from ``pv`` by
``T (T + 2M) / (T + M) = pv``), the slab of mass thickness ``t = R(E_1) - R(E)``::

    chi_c^2     = K z^2 sum_i w_i Z_i^2 / A_i  int_0^t dt' / (pv)^2
    <ln chi_a^2> = int ln chi_a^2 dt' / (pv)^2 / int dt' / (pv)^2
    delta_b     = int delta_b(beta) dt' / (pv)^2 / int dt' / (pv)^2     (Fano 1954)
    b           = ln chi_c^2 - <ln chi_a^2> - ln 1.167 + delta_b

Every formula, constant and the quadrature (64 Gauss-Legendre nodes on ``t' = t (1 - (1 - s)^2)``)
is that of the committed P1/P11 reference ``validation/scripts/ion_em/exact_omega.py`` (equation
source; ``electrons = "fano"``, the reference default). The differences are only arithmetic: the
composition sums are folded into per-material constants (:func:`pack_composition`), the logarithm
of ``beta^2 / (1 - beta^2)`` is evaluated as ``(pc / M)^2``, and the CSDA path is read from the
project transport-table layout (exact closed-form range in the log-log bin, the inverse-range
table; ``ionmc.transport.tables``) instead of the interpolating ``StoppingTable.range_at``.

Evaluation. Two forms, selected by ``PhysicsOptions.mcs_mapping_eval`` (decision 0043 section 3,
"Evaluation"): ``"root"``, the per-step root solve below, and ``"table"`` (the default for z >= 2),
the decision's allowed alternative: a per-(material, ion) tabulation of ``pv_eq/pv`` on a uniform
grid in ``(ln pv, ln(p1v1/pv))`` (a shear of the decision's ``(ln pv, ln p1v1)`` plane), whose
nodes are computed by this same root solve (:func:`build_pv_eq_table`) and which is read by
bilinear interpolation (:func:`pv_equivalent_table`). The error of the resulting per-step variance
is certified at build for every produced table by a domain-wide analytic bound
(:func:`certify_pv_eq_table`, threshold :data:`TABLE_VARIANCE_BOUND` = 1e-3;
``tests/ionmc/test_ion_mcs_table.py``). The root solve:
:func:`pv_equivalent` solves ``g(u) = b_p(e^u, r e^u) - b_ion(pv, p1v1) = 0`` for
``u = ln pv_eq``, ``r = p1v1 / pv``, by a bracketed Illinois (modified regula falsi) iteration
with a bisection safeguard. The bracket is the proton path itself: from its lowest energy to the
largest ``pv_eq`` whose ``p1v1_eq = r pv_eq`` is still on the proton path. ``g`` is monotone
increasing in ``pv_eq`` at fixed ``r`` (P11), the bracket stays valid, and the iteration is capped
at :data:`MAX_ITERATIONS`; there is no failure mode other than the returned status.
Status codes (:data:`STATUS_NAMES`): 0 converged, 1 invalid ion state (outside the ion path or
``pv >= p1v1``), 2 no sign change in the bracket (``b_ion`` outside the proton ``b`` range at this
``r``), 3 the proton path is too short for this ``r``, 4 not converged within the iteration cap,
5 (table form only) the state lies outside the table or in a cell with an invalid node.
A caller must treat any nonzero status as an invalid state (fail closed); ``f_dM`` is then 0.

The proton path of the mapping is the Bethe proton table of the same material to
:data:`PROTON_PATH_EMAX_MEV` (:func:`proton_path`), independent of the transport stopping source,
exactly as in the reference. The default 500 MeV proton transport tables are not touched.

z = 1 and M = m_p: ``ion_is_proton = 1`` returns ``pv_eq = pv`` and ``p1v1_eq = p1v1`` without
evaluating anything (exact identity: ``f_dM`` is then ``em.fdm(pv, p1v1)`` bit for bit).

Function arguments (all ``float`` arrays are 1-D, in the precision of the kernel):

* ``pth``: a packed CSDA path (:func:`pack_path`, energy axis = TOTAL kinetic energy ``A e`` as in
  the transport tables). Header ``[ln_e0, inv_dln_e, ln_r0, inv_dln_r, n_e, n_r, e_min, e_max]``,
  then ``r_mass[n_e]``, ``f_mass[n_e]``, ``d_f[n_e]``, ``ln_e_of_r[n_r]`` (the arrays of
  ``TransportTables`` for one material);
* ``comp``: packed material composition (:func:`pack_composition`);
* ``quad``: the Gauss-Legendre rule (:func:`pack_quadrature`).

Units: energies and ``pv`` in MeV, ranges in g/cm2.
"""

import functools
import math
from collections.abc import Sequence
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any

import numpy as np
import warp as wp

from ionmc._validate import fail
from ionmc._wpfunc import check_real, named_func
from ionmc.materials import ELEMENTS, N_A, Material
from ionmc.physics.em import FDM_COEFFICIENTS, make_em
from ionmc.physics.kinematics import make_kinematics
from ionmc.physics.projectiles import PROTON, Projectile
from ionmc.physics.stopping import (
    BetheStoppingSource,
    StoppingTable,
    exact_loglog_range_increments,
)
from ionmc.transport.funcs import make_transport_funcs

wp.set_module_options({"enable_backward": False})

# --- constants of the Moliere quantities (identical to validation/scripts/ion_em/moliere_hanson.py)
ALPHA = 1.0 / 137.035999084
R_E_CM = 2.8179403262e-13
ME_MEV = 0.51099895000
EULER_GAMMA = 0.5772156649015329
K_CHI_C = 4.0 * math.pi * N_A * R_E_CM**2 * ME_MEV**2
CHI_0_MEV = ALPHA * ME_MEV / 0.885
SCREEN_A, SCREEN_B = 1.13, 3.76
LN_OMEGA_FACTOR = 2.0 * EULER_GAMMA - 1.0  # ln 1.167
FANO_CONST = 1130.0
U_IN_NEG = {1: 3.6}
U_IN_NEG_DEFAULT = 5.0

N_QUAD = 64
"""Gauss-Legendre nodes of the slab integrals inside b (the reference ``N_SLAB``)."""
MAX_ITERATIONS = 100
"""Iteration cap of the root solve (a failure is status 4, never a silent stop)."""
RESIDUAL_TOL = {"float64": 1.0e-13, "float32": 5.0e-5}
"""Convergence criterion |b_p - b_ion| of the root solve (the I11-map bound is 1e-10, float64)."""
STEP_TOL = {"float64": 1.0e-14, "float32": 1.0e-6}
"""A step of the iteration (or the bracket) smaller than this in ``ln pv_eq`` (relative change of
``pv_eq``) ends the solve as converged: for very thin slabs the rounding noise of ``b`` (the
cancellation in ``t = R(E_1) - R(E)``) exceeds the residual tolerance."""
EDGE_TOL = {"float64": 1.0e-12, "float32": 1.0e-5}
"""Relative tolerance of the path-table edges (the reference uses 1e-12)."""
BRACKET_HI_MARGIN = {"float64": 1.0e-7, "float32": 1.0e-5}
"""Relative offset below the proton path's top energy of the upper bracket end."""
BRACKET_LO_MARGIN = 1.0e-6
"""Relative offset above the proton path's lowest energy of the lower bracket end."""

PROTON_PATH_EMAX_MEV = 4000.0
"""Upper edge [MeV] of the proton path used by the mapping. The proton-equivalent states of the ion
domain reach 0.8 GeV (He 250 MeV/u) and 2.5 GeV (C, O at 430 MeV/u) in water (P11) and, at the
domain edge, up to 3.4 GeV in beryllium (C, O); the edge is 4 GeV (at least the 3 GeV of row I0)
so that every ``materials.py`` material is covered at 430 MeV/u. Grid: ``log_grid(1, 4000, 200
per decade)`` = 722 nodes, the same Bethe model and options as the 500 MeV transport table."""
PROTON_PATH_POINTS_PER_DECADE = 200
FDM_DOMAIN_EDGE_MEV = 300.0
"""Proton kinetic energy up to which f_dM is qualified (decision 0043, context and section 5)."""
DEPTH_FRACTIONS = (0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97)
"""Depth fractions x/R of the source state scan (:func:`fdm_domain_status`)."""

STATUS_NAMES = {
    0: "converged",
    1: "invalid ion state",
    2: "no sign change (b_ion outside the proton b range)",
    3: "proton path too short",
    4: "not converged",
    5: "state outside the pv_eq/pv table",
}

TABLE_N_PV = 128
"""Default nodes of the ``ln pv`` axis of the ``pv_eq/pv`` table (the variance error is certified
at build, :data:`TABLE_VARIANCE_BOUND`)."""
TABLE_N_R = 128
"""Default nodes of the ``ln(p1v1/pv)`` axis of the table."""
TABLE_INVALID_NODE = -1.0
"""Value stored at a node whose root solve failed (outside the paths); a cell touching it returns
status 5."""
TABLE_EDGE_TOL = {"float64": 1.0e-9, "float32": 1.0e-4}
"""Tolerance (in grid-index units) of the table edges."""
TABLE_PV_MARGIN = 1.0e-3
"""Relative margin of the table's largest ``p1v1`` above ``pv(E0 + 6 sigma_E)``."""
TABLE_HEADER = 8
TABLE_VARIANCE_BOUND = 1.0e-3
"""Decision 0043 bound of the table (0.1 %): the per-step MCS variance (proportional to T and
theta^2) from the interpolated ``pv_eq`` against the root solve, relative. Certified at build for
every produced table (:func:`certify_pv_eq_table`); a table above it fails closed."""
TABLE_SAFETY_FACTOR = 2.0
"""Safety factor on the second-difference estimate of the bilinear error bound
(:func:`certify_pv_eq_table`); checked empirically by the sampled cross-check (2026-11-08: with
factor 1 the cross-check passes in all 70 tables of He, C, O at both domain ends in the eleven
built-in materials and a user material, with 0.75 it fails in every one; 2 doubles that margin)."""
TABLE_SAMPLE_N = 5
"""Interior sample points per cell and axis of the cross-check (``N x N`` per cell)."""
TABLE_ROOT_TOL_REL = 1.0e-10
"""Relative allowance for the root solve's convergence at the nodes and the stencil (float64)."""
TABLE_MIN_LOSS = 1.0e-4
"""Scope of the relative contract (:func:`certify_pv_eq_table`): steps of a particle that has lost
at least this fraction of its initial ``pv`` at the state that enters the memory term, i.e.
``pv_mid <= (1 - L) p1v1`` for a later step and ``pv_end <= (1 - L) p1v1`` for the first step
(a few micrometres of path for 50 MeV/u helium). Before that the memory term
``1 - (pv/p1v1)^2 < 2e-4`` drives the f_dM fit towards its zero (clamped to 0 for both paths;
e.g. near ``1 - (pv/p1v1)^2 ~ 1e-5`` at the lowest helium energies in high-Z materials), where a
relative error is not defined; the absolute error there is still bounded by
``|c2 + c3 l1| |dl2|`` times the common prefactor."""
TABLE_CERT_CONTRACT = (
    "certified at build for every produced table: analytic bilinear error bound of pv_eq/pv from "
    "second differences on a 2x-refined root-solve grid (3 x 3 stencil per solved cell, maximum "
    "over the cell and its neighbours) with safety factor 2, converted exactly into the relative "
    "per-step MCS variance error (T, theta^2) for later steps and for the first step over all "
    "physical end points pv_end in "
    "(0, min(pv_mid, (1 - 1e-4) p1v1)]; cross-checked by 5 x 5 interior samples per cell (later "
    "steps; first step at the upper end point and at the table's lowest pv); fails closed when "
    "the analytic bound exceeds the threshold or a sample exceeds its cell's bound; states of a "
    "particle that has lost less than 1e-4 of its initial pv (TABLE_MIN_LOSS: the f_dM memory "
    "term near the fit's zero) are outside the relative contract"
)
"""The certification contract, as recorded in the summary and the capability report."""
LOG10_E = 0.4342944819032518
"""log10(e): the shift of the first step's analytic memory average (``ionmc.physics.ion_mcs``)."""

_HEADER = 8


# --- host-side packing -------------------------------------------------------------------------
def pack_quadrature(n: int = N_QUAD) -> np.ndarray:
    """``[n, qn_0 .. qn_{n-1}, qw_0 .. qw_{n-1}]`` of the reference rule: the node ``t' = t qn_k``
    and the weight ``t qw_k`` of a slab of thickness ``t`` (Gauss-Legendre ``s`` mapped by
    ``t' = t (1 - (1 - s)^2)``, which clusters nodes where ``1/(pv)^2`` grows fastest)."""
    gx, gw = np.polynomial.legendre.leggauss(n)
    s = 0.5 * (gx + 1.0)
    qn = 1.0 - (1.0 - s) ** 2
    qw = 2.0 * (1.0 - s) * 0.5 * gw
    return np.concatenate(([float(n)], qn, qw))


def pack_composition(material: Material) -> np.ndarray:
    """Packed composition of ``material`` for the Fano form of b (see the module docstring)::

        [n_el, zsum, ca0, d0, d1, Z_0, W_0, Z_1, W_1, ...]

    ``zsum = sum w Z^2/A``; ``W_i = w_i Z_i^2/A_i / zsum``; ``ca0 = sum W_i 2 ln(chi_0 Z_i^(1/3))``
    (``ln chi_a^2 = ca0 - ln (pc)^2 + sum_i W_i ln(1.13 + 3.76 (z Z_i alpha)^2 / beta^2)``);
    ``d0 + d1 (ln((pc/M)^2) - beta^2/2)`` is Fano's increase ``delta_b`` of b with
    ``d1 = sum w Z/A / zsum`` and ``d0 = sum (w Z/A) (ln(1130 Z^(-4/3)) - u_in) / zsum``."""
    els = [(ELEMENTS[sym], w) for sym, w in material.mass_fractions.items()]
    zsum = math.fsum(w * e.Z * e.Z / e.A_g_mol for e, w in els)
    weights = [w * e.Z * e.Z / e.A_g_mol / zsum for e, w in els]
    ca0 = math.fsum(
        wi * 2.0 * math.log(CHI_0_MEV * e.Z ** (1.0 / 3.0))
        for (e, _w), wi in zip(els, weights, strict=False)
    )
    d1 = math.fsum(w * e.Z / e.A_g_mol for e, w in els) / zsum
    d0 = (
        math.fsum(
            w
            * e.Z
            / e.A_g_mol
            * (math.log(FANO_CONST * e.Z ** (-4.0 / 3.0)) + U_IN_NEG.get(e.Z, U_IN_NEG_DEFAULT))
            for e, w in els
        )
        / zsum
    )
    out = [float(len(els)), zsum, ca0, d0, d1]
    for (e, _w), wi in zip(els, weights, strict=False):
        out += [float(e.Z), wi]
    return np.array(out, dtype=np.float64)


@dataclass(frozen=True, eq=False)
class MappingPath:
    """A packed CSDA path (float64) of one projectile in one material (module docstring)."""

    projectile: Projectile
    pack: np.ndarray
    e_min_mev: float
    e_max_mev: float


def pack_path(table: StoppingTable) -> MappingPath:
    """Pack ``table`` into the transport-table layout of one material, on the TOTAL kinetic energy
    axis ``E = A e`` (``ionmc.transport.tables.TransportTables.from_stopping_tables`` for a
    proton gives bit-identical arrays, tested; that class is restricted to z = 1 projectiles)."""
    from ionmc.transport.tables import (  # noqa: PLC0415
        DEFAULT_RANGE_POINTS_PER_DECADE,
        MIN_POINTS_PER_DECADE,
        _count,
        _exact_inverse,
    )

    proj = table.projectile
    a_proj = float(proj.a)
    ln_a = math.log(a_proj)
    n_e = _count(table.energy_per_u[-1] / table.energy_per_u[0], MIN_POINTS_PER_DECADE)
    n_r = _count(
        table.csda_range_g_cm2[-1] / table.csda_range_g_cm2[0], DEFAULT_RANGE_POINTS_PER_DECADE
    )
    lo, hi = math.log(table.energy_per_u[0]), math.log(table.energy_per_u[-1])
    grid_u = np.linspace(lo, hi, n_e)
    lns = np.interp(grid_u, np.log(table.energy_per_u), np.log(table.s_el_mass))
    grid = grid_u + ln_a
    ln_f = grid - lns
    inc = exact_loglog_range_increments(grid, lns, 1.0)
    r = float(table.csda_range_g_cm2[0]) + np.concatenate(([0.0], np.cumsum(inc)))
    f_mass = np.exp(ln_f)
    d_f = np.zeros(n_e)
    d_f[:-1] = np.diff(ln_f)
    inv_dln_e = (n_e - 1) / (hi - lo)
    rlo, rhi = float(np.log(r[0])), float(np.log(r[-1]))
    rgrid = np.linspace(rlo, rhi, n_r)
    inv_dln_r = (n_r - 1) / (rhi - rlo)
    ln_e_of_r = _exact_inverse(np.exp(rgrid), r, f_mass, d_f, 1.0 / inv_dln_e, grid)
    ln_e_of_r[0], ln_e_of_r[-1] = grid[0], grid[-1]
    e_min = float(table.energy_per_u[0] * a_proj)
    e_max = float(table.energy_per_u[-1] * a_proj)
    header = [grid[0], inv_dln_e, rgrid[0], inv_dln_r, float(n_e), float(n_r), e_min, e_max]
    pack = np.concatenate((header, r, f_mass, d_f, ln_e_of_r)).astype(np.float64)
    return MappingPath(proj, pack, e_min, e_max)


def pack_transport_row(tables: Any, material: int = 0) -> np.ndarray:
    """The same packing from the arrays of a ``TransportTables`` row (an ion's or the proton's);
    equals :func:`pack_path` of the table the tables were built from."""
    n_e, n_r = tables.n_e, tables.n_r
    header = [
        tables.ln_e0[material],
        tables.inv_dln_e[material],
        tables.ln_r0[material],
        tables.inv_dln_r[material],
        float(n_e),
        float(n_r),
        tables.e_min_mev[material],
        tables.e_max_mev[material],
    ]
    return np.concatenate(
        (
            header,
            tables.r_mass[material],
            tables.f_mass[material],
            tables.d_f[material],
            tables.ln_e_of_r[material],
        )
    ).astype(np.float64)


_PROTON_PATHS: dict[tuple[str, float], MappingPath] = {}


def proton_path(material: Material, e_max_mev: float = PROTON_PATH_EMAX_MEV) -> MappingPath:
    """The proton path of the mapping: the Bethe proton table of ``material`` from 1 MeV to
    ``e_max_mev`` (default :data:`PROTON_PATH_EMAX_MEV`), 200 points per decade, default
    :class:`~ionmc.physics.stopping.BetheOptions`. Cached per material definition."""
    from ionmc.transport.tables import material_fingerprint  # noqa: PLC0415

    key = (material_fingerprint(material), float(e_max_mev))
    if key not in _PROTON_PATHS:
        src = BetheStoppingSource(
            e_min_per_u=1.0,
            e_max_per_u=float(e_max_mev),
            points_per_decade=PROTON_PATH_POINTS_PER_DECADE,
        )
        _PROTON_PATHS[key] = pack_path(src.table(material, PROTON))
    return _PROTON_PATHS[key]


# --- Warp functions ----------------------------------------------------------------------------
@functools.cache
def make_exact_omega(real: type) -> SimpleNamespace:
    """Return the exact-Omega mapping functions for precision ``real``."""
    name = check_real(real)
    kin = make_kinematics(real)
    em = make_em(real)
    tf = make_transport_funcs(real)
    k_chi_c = wp.constant(real(K_CHI_C))
    alpha = wp.constant(real(ALPHA))
    scr_a = wp.constant(real(SCREEN_A))
    scr_b = wp.constant(real(SCREEN_B))
    ln_omega = wp.constant(real(LN_OMEGA_FACTOR))
    m_p = wp.constant(real(PROTON.mass_mev))
    edge_tol = wp.constant(real(EDGE_TOL[name]))
    res_tol = wp.constant(real(RESIDUAL_TOL[name]))
    step_tol = wp.constant(real(STEP_TOL[name]))
    hi_margin = wp.constant(real(BRACKET_HI_MARGIN[name]))
    lo_margin = wp.constant(real(BRACKET_LO_MARGIN))

    @named_func(name)
    def kinetic_from_pv(pv: real, m_mev: real) -> real:
        """Kinetic energy T [MeV] with ``T (T + 2M) / (T + M) = pv``: the positive root of
        ``T^2 + (2M - pv) T - pv M = 0``, in the form without cancellation for either sign of
        ``q = 2M - pv``."""
        q = real(2.0) * m_mev - pv
        root = wp.sqrt(q * q + real(4.0) * pv * m_mev)
        t = real(0.5) * (root - q)
        if q >= real(0.0):
            t = real(2.0) * pv * m_mev / (q + root)
        return t

    @named_func(name)
    def path_range(pth: wp.array(dtype=real), e_mev: real) -> real:
        """CSDA range [g/cm2] at total kinetic energy ``e_mev`` on the packed path: the exact
        closed form in the log-log bin (shared ``range_in_bin``), clamped to the path."""
        n_e = int(pth[4])
        i, f = tf.log_bin_index(e_mev, pth[0], pth[1], n_e)
        return tf.range_in_bin(
            pth[8 + i], pth[8 + n_e + i], pth[8 + 2 * n_e + i], real(1.0) / pth[1], f
        )

    @named_func(name)
    def path_energy(pth: wp.array(dtype=real), r_g_cm2: real) -> real:
        """Total kinetic energy [MeV] with CSDA range ``r_g_cm2`` on the packed path (the
        inverse-range table, shared bin location and interpolation), clamped to the path."""
        n_e = int(pth[4])
        n_r = int(pth[5])
        i, f = tf.log_bin_index(r_g_cm2, pth[2], pth[3], n_r)
        base = 8 + 3 * n_e
        return tf.interp_exp(pth[base + i], pth[base + i + 1], f)

    @named_func(name)
    def ln_chi_a2(comp: wp.array(dtype=real), z: real, b2: real, pc2: real) -> real:
        """Compound ``ln chi_a^2`` (Moliere screening angle with the Born screening factor) at
        ``beta^2 = b2`` and ``(pc)^2 = pc2`` of a projectile of charge ``z``."""
        n_el = int(comp[0])
        scr = real(0.0)
        for k in range(n_el):
            za = z * comp[5 + 2 * k] * alpha
            scr = scr + comp[6 + 2 * k] * wp.log(scr_a + scr_b * za * za / b2)
        return comp[2] - wp.log(pc2) + scr

    @named_func(name)
    def fano_delta_b(comp: wp.array(dtype=real), b2: real, bg2: real) -> real:
        """Fano's increase of Moliere's b at ``beta^2 = b2`` with ``bg2 = (pc/M)^2 =
        beta^2/(1 - beta^2)``: ``d0 + d1 (ln bg2 - beta^2/2)`` (compound weighting of the
        reference ``fano_delta_b``)."""
        return comp[3] + comp[4] * (wp.log(bg2) - real(0.5) * b2)

    @named_func(name)
    def b_slab(
        pth: wp.array(dtype=real),
        comp: wp.array(dtype=real),
        quad: wp.array(dtype=real),
        z: real,
        m_mev: real,
        pv: real,
        p1v1: real,
    ) -> tuple[real, int]:
        """``b = ln Omega + delta_b_Fano`` of the slab from ``p1v1`` (initial) down to ``pv``
        (current) along the CSDA path ``pth`` of the projectile (charge ``z``, mass ``m_mev``)
        and the status: 0 ok, 1 ``pv`` not in ``(0, p1v1)`` or an energy outside the path
        table, 2 slab thickness not positive (``b`` is 0 for a nonzero status)."""
        b = real(0.0)
        status = int(0)
        e = real(0.0)
        e1 = real(0.0)
        if pv > real(0.0) and pv < p1v1:
            e = kinetic_from_pv(pv, m_mev)
            e1 = kinetic_from_pv(p1v1, m_mev)
            if e < pth[6] * (real(1.0) - edge_tol) or e1 > pth[7] * (real(1.0) + edge_tol):
                status = 1
            if not e < e1:
                status = 1
        else:
            status = 1
        if status == 0:
            r_e = path_range(pth, e)
            r_1 = path_range(pth, e1)
            t = r_1 - r_e
            if not t > real(0.0):
                status = 2
            else:
                n_q = int(quad[0])
                integral = real(0.0)
                lna_acc = real(0.0)
                db_acc = real(0.0)
                for k in range(n_q):
                    # depth t' = t qn_k from the p1v1 end; the range left at t' is
                    # R(E) + (t - t'), formed without cancellation
                    rem = r_e + (t - t * quad[1 + k])
                    en = path_energy(pth, rem)
                    pc2 = en * (en + real(2.0) * m_mev)
                    pvn = kin.pv_mev(en, m_mev)
                    b2 = kin.beta2(en, m_mev)
                    w = t * quad[1 + n_q + k] / (pvn * pvn)
                    integral = integral + w
                    lna_acc = lna_acc + w * ln_chi_a2(comp, z, b2, pc2)
                    db_acc = db_acc + w * fano_delta_b(comp, b2, pc2 / (m_mev * m_mev))
                chi_c2 = k_chi_c * z * z * comp[1] * integral
                b = wp.log(chi_c2) - lna_acc / integral - ln_omega + db_acc / integral
        return b, status

    @named_func(name)
    def b_residual(
        pro_pth: wp.array(dtype=real),
        comp: wp.array(dtype=real),
        quad: wp.array(dtype=real),
        u: real,
        r: real,
        target: real,
    ) -> tuple[real, int]:
        """Root function ``g(u) = b_p(pv_eq, r pv_eq) - target`` at ``pv_eq = exp(u)`` on the
        proton path (``r = p1v1/pv > 1``), and the status of the slab."""
        p = wp.exp(u)
        bp, st = b_slab(pro_pth, comp, quad, real(1.0), m_p, p, r * p)
        return bp - target, st

    @named_func(name)
    def pv_equivalent(
        ion_pth: wp.array(dtype=real),
        pro_pth: wp.array(dtype=real),
        comp: wp.array(dtype=real),
        quad: wp.array(dtype=real),
        z: real,
        m_ion_mev: real,
        pv: real,
        p1v1: real,
        ion_is_proton: int,
    ) -> tuple[real, int, int]:
        """Proton-equivalent ``pv_eq`` of the ion state ``(pv, p1v1)``, a status (see the module
        docstring: 0 converged) and the number of root iterations (``MAX_ITERATIONS`` bound).

        ``ion_is_proton = 1`` (z = 1, M = m_p, ``ion_pth`` the proton path) returns ``pv``
        exactly. Otherwise ``target = b_ion(pv, p1v1)`` is solved against the proton slab
        ``b_p(pv_eq, r pv_eq)``, ``r = p1v1/pv``, in ``u = ln pv_eq`` by the Illinois iteration
        on a valid bracket (bisection step whenever the secant point leaves the bracket)."""
        out = pv
        status = int(0)
        iters = int(0)
        if ion_is_proton == 0:
            target, st_ion = b_slab(ion_pth, comp, quad, z, m_ion_mev, pv, p1v1)
            if st_ion != 0:
                status = 1
            else:
                r = p1v1 / pv
                pv_lo = kin.pv_mev(pro_pth[6] * (real(1.0) + lo_margin), m_p)
                pv_hi = kin.pv_mev(pro_pth[7] * (real(1.0) - hi_margin), m_p) / r
                if not pv_hi > pv_lo:
                    status = 3
                else:
                    ua = wp.log(pv_lo)
                    ub = wp.log(pv_hi)
                    ga, sta = b_residual(pro_pth, comp, quad, ua, r, target)
                    gb, stb = b_residual(pro_pth, comp, quad, ub, r, target)
                    if sta != 0 or stb != 0:
                        status = 3
                    elif ga > real(0.0) or gb < real(0.0):
                        status = 2
                    else:
                        done = int(0)
                        side = int(0)
                        uc = ub
                        up = ua
                        if wp.abs(ga) <= res_tol:
                            uc = ua
                            done = 1
                        elif wp.abs(gb) <= res_tol:
                            done = 1
                        for it in range(MAX_ITERATIONS):
                            if done == 0:
                                iters = it + 1
                                uc = real(0.5) * (ua + ub)
                                if gb - ga > real(0.0):
                                    us = (gb * ua - ga * ub) / (gb - ga)
                                    if us > ua and us < ub:
                                        uc = us
                                gc, stc = b_residual(pro_pth, comp, quad, uc, r, target)
                                step = wp.abs(uc - up)
                                up = uc
                                if stc != 0:
                                    status = 3
                                    done = 1
                                elif (
                                    wp.abs(gc) <= res_tol or ub - ua <= step_tol or step <= step_tol
                                ):
                                    done = 1
                                elif gc > real(0.0):
                                    ub = uc
                                    gb = gc
                                    if side == 1:
                                        ga = real(0.5) * ga
                                    side = 1
                                else:
                                    ua = uc
                                    ga = gc
                                    if side == -1:
                                        gb = real(0.5) * gb
                                    side = -1
                        if done == 0:
                            status = 4
                        out = wp.exp(uc)
        return out, status, iters

    @named_func(name)
    def fdm_equivalent(
        ion_pth: wp.array(dtype=real),
        pro_pth: wp.array(dtype=real),
        comp: wp.array(dtype=real),
        quad: wp.array(dtype=real),
        z: real,
        m_ion_mev: real,
        pv: real,
        p1v1: real,
        ion_is_proton: int,
    ) -> tuple[real, real, real, int]:
        """``(f_dM(pv_eq, p1v1_eq), pv_eq, p1v1_eq, status)`` with ``p1v1_eq = pv_eq p1v1/pv``
        (``p1v1_eq = p1v1`` exactly for ``ion_is_proton = 1``). ``f_dM`` is the project proton
        function ``em.fdm`` (the qualified fit, zero where it diverges); it is 0 for a nonzero
        status."""
        peq, status, iters = pv_equivalent(
            ion_pth, pro_pth, comp, quad, z, m_ion_mev, pv, p1v1, ion_is_proton
        )
        p1eq = p1v1
        if ion_is_proton == 0:
            p1eq = peq * (p1v1 / pv)
        f = real(0.0)
        if status == 0:
            f = em.fdm(peq, p1eq)
        return f, peq, p1eq, status

    tab_tol = wp.constant(real(TABLE_EDGE_TOL[name]))

    @named_func(name)
    def pv_equivalent_table(tab: wp.array(dtype=real), pv: real, p1v1: real) -> tuple[real, int]:
        """Tabulated ``(pv_eq, status)`` of the ion state ``(pv, p1v1)``: bilinear interpolation of
        ``pv_eq/pv`` in ``(x, y) = (ln pv, ln(p1v1/pv))`` on the packed table row ``tab``
        (:func:`build_pv_eq_table`). Status 5 (``pv_eq`` = ``pv``) when the state lies outside the
        grid (tolerance :data:`TABLE_EDGE_TOL` in index units) or a corner node of its cell is
        invalid; 0 otherwise. The caller applies the thin-slab guard (``p1v1/pv >= 1 + 1e-6``,
        the first node of the ``y`` axis)."""
        out = pv
        status = int(5)
        n_x = int(tab[2])
        n_y = int(tab[5])
        fx = (wp.log(pv) - tab[0]) * tab[1]
        fy = wp.log(p1v1 / pv) * tab[4] - tab[3] * tab[4]
        if (
            fx >= -tab_tol
            and fx <= real(n_x - 1) + tab_tol
            and fy >= -tab_tol
            and fy <= real(n_y - 1) + tab_tol
        ):
            fx = wp.min(wp.max(fx, real(0.0)), real(n_x - 1))
            fy = wp.min(wp.max(fy, real(0.0)), real(n_y - 1))
            i = wp.min(int(wp.floor(fx)), n_x - 2)
            j = wp.min(int(wp.floor(fy)), n_y - 2)
            tx = fx - real(i)
            ty = fy - real(j)
            base = 8 + i * n_y + j
            q00 = tab[base]
            q01 = tab[base + 1]
            q10 = tab[base + n_y]
            q11 = tab[base + n_y + 1]
            if q00 > real(0.0) and q01 > real(0.0) and q10 > real(0.0) and q11 > real(0.0):
                q = (real(1.0) - tx) * ((real(1.0) - ty) * q00 + ty * q01) + tx * (
                    (real(1.0) - ty) * q10 + ty * q11
                )
                out = q * pv
                status = 0
        return out, status

    return SimpleNamespace(
        pv_equivalent_table=pv_equivalent_table,
        kinetic_from_pv=kinetic_from_pv,
        path_range=path_range,
        path_energy=path_energy,
        ln_chi_a2=ln_chi_a2,
        fano_delta_b=fano_delta_b,
        b_slab=b_slab,
        b_residual=b_residual,
        pv_equivalent=pv_equivalent,
        fdm_equivalent=fdm_equivalent,
        real=name,
    )


# --- host-side evaluation (python twin, float64) ------------------------------------------------
@dataclass(frozen=True, eq=False)
class ExactOmegaMapping:
    """The mapping of one projectile in one material, evaluated by the float64 python twin.

    ``ion_path`` is the projectile's own path (the proton path itself for ``PROTON``); ``comp``
    and ``quad`` as in the module docstring. Lists are used for the twin's array arguments (plain
    indexing is the fastest array access in Python scope)."""

    projectile: Projectile
    material: Material
    ion_path: MappingPath
    proton: MappingPath
    comp: np.ndarray
    quad: np.ndarray
    ion_is_proton: int

    @classmethod
    def build(
        cls,
        projectile: Projectile,
        material: Material,
        ion_table: StoppingTable | None = None,
        proton: MappingPath | None = None,
    ) -> "ExactOmegaMapping":
        """Mapping of ``projectile`` in ``material``. ``ion_table`` is the projectile's stopping
        table (default: the Bethe table 1 to 500 MeV/u with default options, as the reference);
        ``proton`` the proton path (default :func:`proton_path`)."""
        pro = proton_path(material) if proton is None else proton
        if projectile == PROTON:
            ion = pro
            flag = 1
        else:
            tab = (
                BetheStoppingSource().table(material, projectile)
                if ion_table is None
                else ion_table
            )
            if tab.projectile != projectile:
                raise fail("ion_table does not describe the requested projectile")
            ion = pack_path(tab)
            flag = 0
        return cls(
            projectile,
            material,
            ion,
            pro,
            pack_composition(material),
            pack_quadrature(),
            flag,
        )

    @functools.cached_property
    def _twin(self) -> SimpleNamespace:
        from ionmc._wpfunc import python_twin  # noqa: PLC0415

        return python_twin(make_exact_omega)

    @functools.cached_property
    def _lists(self) -> tuple[list[float], list[float], list[float], list[float]]:
        return (
            self.ion_path.pack.tolist(),
            self.proton.pack.tolist(),
            self.comp.tolist(),
            self.quad.tolist(),
        )

    def pv_of_energy(self, e_mev: float) -> float:
        """``pv`` [MeV] of the projectile at total kinetic energy ``e_mev``."""
        m = self.projectile.mass_mev
        return e_mev * (e_mev + 2.0 * m) / (e_mev + m)

    def b_ion(self, pv: float, p1v1: float) -> tuple[float, int]:
        """``(b_ion(pv, p1v1), status)`` on the projectile's own path."""
        ion, _pro, comp, quad = self._lists
        return self._twin.b_slab(
            ion, comp, quad, float(self.projectile.z), self.projectile.mass_mev, pv, p1v1
        )

    def b_proton(self, pv: float, p1v1: float) -> tuple[float, int]:
        """``(b_p(pv, p1v1), status)`` on the proton path."""
        _ion, pro, comp, quad = self._lists
        return self._twin.b_slab(pro, comp, quad, 1.0, PROTON.mass_mev, pv, p1v1)

    def pv_equivalent(self, pv: float, p1v1: float) -> tuple[float, int, int]:
        """``(pv_eq, status, iterations)`` of the state (module docstring)."""
        ion, pro, comp, quad = self._lists
        return self._twin.pv_equivalent(
            ion,
            pro,
            comp,
            quad,
            float(self.projectile.z),
            self.projectile.mass_mev,
            pv,
            p1v1,
            self.ion_is_proton,
        )

    def fdm_equivalent(self, pv: float, p1v1: float) -> tuple[float, float, float, int]:
        """``(f_dM(pv_eq, p1v1_eq), pv_eq, p1v1_eq, status)``."""
        ion, pro, comp, quad = self._lists
        return self._twin.fdm_equivalent(
            ion,
            pro,
            comp,
            quad,
            float(self.projectile.z),
            self.projectile.mass_mev,
            pv,
            p1v1,
            self.ion_is_proton,
        )

    def residual(self, pv: float, p1v1: float, pv_eq: float) -> float:
        """``b_p(pv_eq, r pv_eq) - b_ion(pv, p1v1)``, ``r = p1v1/pv`` (the I11-map residual)."""
        r = p1v1 / pv
        return self.b_proton(pv_eq, r * pv_eq)[0] - self.b_ion(pv, p1v1)[0]

    def state_at_depth(self, e0_mev: float, fraction: float) -> tuple[float, float]:
        """``(pv, p1v1)`` of the projectile at depth ``fraction`` of its CSDA range from an initial
        total kinetic energy ``e0_mev`` (the ion path of the packed table)."""
        ion = self._lists[0]
        r1 = float(self._twin.path_range(ion, e0_mev))
        e = float(self._twin.path_energy(ion, (1.0 - fraction) * r1))
        return self.pv_of_energy(e), self.pv_of_energy(e0_mev)


def proton_kinetic_from_pv(pv: float) -> float:
    """Proton kinetic energy [MeV] of ``pv`` (float64)."""
    m = PROTON.mass_mev
    q = 2.0 * m - pv
    root = math.sqrt(q * q + 4.0 * pv * m)
    return 2.0 * pv * m / (q + root) if q >= 0.0 else 0.5 * (root - q)


# --- source status (decision 0043 section 3 r6, acceptance rows I0 and P11b) ---------------------
@dataclass(frozen=True)
class FdmDomainStatus:
    """The f_dM domain status of a source (row I0), persisted in the effective configuration.

    ``p1v1_eq_kinetic_energy_mev`` is the proton kinetic energy of ``p1v1_eq`` of the source:
    the largest over the scanned depth fractions ``x/R`` (:data:`DEPTH_FRACTIONS`) of the ion
    slowing down from the nominal energy, and over the materials; ``..._hi`` the same for the
    energy ``E0 + 6 sigma_E`` (the fail-closed bound). ``fdm_extrapolated`` is True when the
    nominal value exceeds :data:`FDM_DOMAIN_EDGE_MEV` (declared, not fail-closed: decision 0043
    section 5). None for a proton, whose f_dM is evaluated directly."""

    mapping: str
    fdm_extrapolated: bool
    p1v1_eq_kinetic_energy_mev: float | None
    p1v1_eq_kinetic_energy_hi_mev: float | None
    fdm_domain_edge_mev: float
    proton_table_edge_mev: float
    depth_fractions: tuple[float, ...]

    def to_dict(self) -> dict[str, Any]:
        """JSON-serialisable form."""
        return {
            "mapping": self.mapping,
            "fdm_extrapolated": self.fdm_extrapolated,
            "p1v1_eq_kinetic_energy_mev": self.p1v1_eq_kinetic_energy_mev,
            "p1v1_eq_kinetic_energy_hi_mev": self.p1v1_eq_kinetic_energy_hi_mev,
            "fdm_domain_edge_mev": self.fdm_domain_edge_mev,
            "proton_table_edge_mev": self.proton_table_edge_mev,
            "depth_fractions": list(self.depth_fractions),
        }


PROTON_FDM_STATUS = FdmDomainStatus(
    mapping="none (proton f_dM)",
    fdm_extrapolated=False,
    p1v1_eq_kinetic_energy_mev=None,
    p1v1_eq_kinetic_energy_hi_mev=None,
    fdm_domain_edge_mev=FDM_DOMAIN_EDGE_MEV,
    proton_table_edge_mev=PROTON_PATH_EMAX_MEV,
    depth_fractions=(),
)
"""The status recorded for a proton source (no mapping: ``fdm_extrapolated`` is False)."""


def _max_p1v1_eq_energy(
    mapping: ExactOmegaMapping, e_total_mev: float, fractions: Sequence[float]
) -> float:
    best = 0.0
    for frac in fractions:
        pv, p1v1 = mapping.state_at_depth(e_total_mev, frac)
        peq, status, _it = mapping.pv_equivalent(pv, p1v1)
        if status != 0:
            raise fail(
                f"exact-Omega mapping of the source state failed for {mapping.projectile.name} "
                f"in {mapping.material.name!r} at E = {e_total_mev:g} MeV, x/R = {frac}: "
                f"{STATUS_NAMES[status]}"
            )
        best = max(best, proton_kinetic_from_pv(peq * p1v1 / pv))
    return best


def check_proton_table_edge(p1v1_eq_hi_mev: float, edge_mev: float = PROTON_PATH_EMAX_MEV) -> None:
    """Fail closed (``UnsupportedCombinationError``) when the proton-equivalent energy of the
    source (at ``E0 + 6 sigma_E``) exceeds the proton table edge ``edge_mev`` (row I0)."""
    if not p1v1_eq_hi_mev <= edge_mev:
        raise fail(
            f"the proton-equivalent energy of the source state ({p1v1_eq_hi_mev:g} MeV at "
            f"E0 + 6 sigma_E) exceeds the proton b-table edge ({edge_mev:g} MeV): the exact-Omega "
            "mapping cannot be evaluated"
        )


def fdm_domain_status(
    projectile: Projectile,
    materials: Sequence[Material],
    e0_mev: float,
    energy_sigma_mev: float = 0.0,
    ion_tables: dict[str, StoppingTable] | None = None,
    fractions: Sequence[float] = DEPTH_FRACTIONS,
    proton_edge_mev: float = PROTON_PATH_EMAX_MEV,
) -> FdmDomainStatus:
    """The row-I0 status of an ion source (``PROTON_FDM_STATUS`` for a proton).

    Evaluates the mapping for the ion slowing down in each of ``materials`` from the total
    kinetic energy ``e0_mev`` and from ``e0_mev + 6 sigma_E``, at the depth fractions
    ``fractions``, takes the largest proton-equivalent ``p1v1_eq`` kinetic energy and raises
    ``UnsupportedCombinationError`` if the ``E0 + 6 sigma_E`` value exceeds ``proton_edge_mev``
    (or the mapping of any scanned state fails). ``ion_tables`` maps the full material identity
    (``ionmc.transport.tables.material_fingerprint``: name, density, composition, I, density-effect
    parameters; not the name alone, since distinct materials may share a name) to the ion's
    stopping table (default: Bethe tables, see :meth:`ExactOmegaMapping.build`)."""
    if projectile == PROTON:
        return PROTON_FDM_STATUS
    if not materials:
        raise fail("at least one material is required")
    from ionmc.transport.tables import material_fingerprint  # noqa: PLC0415

    e_hi = e0_mev + 6.0 * energy_sigma_mev
    nominal = 0.0
    upper = 0.0
    for mat in materials:
        tab = None if ion_tables is None else ion_tables.get(material_fingerprint(mat))
        mp = ExactOmegaMapping.build(projectile, mat, tab, proton_path(mat, proton_edge_mev))
        nominal = max(nominal, _max_p1v1_eq_energy(mp, e0_mev, fractions))
        upper = max(upper, _max_p1v1_eq_energy(mp, e_hi, fractions) if e_hi > e0_mev else nominal)
    check_proton_table_edge(upper, proton_edge_mev)
    return FdmDomainStatus(
        mapping="exact-Omega b-matching (decision 0043 section 3)",
        fdm_extrapolated=nominal > FDM_DOMAIN_EDGE_MEV,
        p1v1_eq_kinetic_energy_mev=nominal,
        p1v1_eq_kinetic_energy_hi_mev=upper,
        fdm_domain_edge_mev=FDM_DOMAIN_EDGE_MEV,
        proton_table_edge_mev=proton_edge_mev,
        depth_fractions=tuple(float(f) for f in fractions),
    )


# --- tabulated evaluation (decision 0043 section 3, "Evaluation", the allowed alternative) -------
@functools.cache
def _table_node_kernel() -> Any:
    """Warp CPU kernel (float64): :func:`pv_equivalent` (the root solve) at the table nodes."""
    xo = make_exact_omega(wp.float64)
    d = wp.float64

    @wp.kernel(module="unique")
    def nodes(
        ion: wp.array(dtype=d),
        pro: wp.array(dtype=d),
        comp: wp.array(dtype=d),
        quad: wp.array(dtype=d),
        z: d,
        mass: d,
        pv: wp.array(dtype=d),
        p1v1: wp.array(dtype=d),
        out: wp.array(dtype=d),
        status: wp.array(dtype=int),
    ):  # type: ignore[no-untyped-def]
        k = wp.tid()
        peq, st, _it = xo.pv_equivalent(ion, pro, comp, quad, z, mass, pv[k], p1v1[k], 0)
        out[k] = peq
        status[k] = st

    return nodes


def root_pv_equivalent_batch(
    ion_pack: np.ndarray,
    proton_pack: np.ndarray,
    comp: np.ndarray,
    quad: np.ndarray,
    z: float,
    mass_mev: float,
    pv: np.ndarray,
    p1v1: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """``(pv_eq, status)`` of many states by the root solve :func:`pv_equivalent` (float64, run as a
    Warp CPU kernel; the same function source as the python twin and the transport kernels)."""
    pv = np.ascontiguousarray(pv, dtype=np.float64)
    p1v1 = np.ascontiguousarray(p1v1, dtype=np.float64)
    n = pv.size
    dev = "cpu"

    def arr(a: np.ndarray) -> Any:
        return wp.array(np.ascontiguousarray(a, dtype=np.float64), dtype=wp.float64, device=dev)

    out = wp.zeros(n, dtype=wp.float64, device=dev)
    st = wp.zeros(n, dtype=int, device=dev)
    wp.launch(
        _table_node_kernel(),
        dim=n,
        inputs=[
            arr(ion_pack), arr(proton_pack), arr(comp), arr(quad), wp.float64(z),
            wp.float64(mass_mev), arr(pv), arr(p1v1), out, st,
        ],
        device=dev,
    )  # fmt: skip
    return out.numpy().copy(), st.numpy().copy()


def table_domain(
    projectile: Projectile, e_mid_min_mev: float, e_lo_mev: float, e_hi_mev: float
) -> tuple[float, float, float]:
    """``(pv_lo, p1v1_min, p1v1_max)`` [MeV] of the table of a run. ``pv_lo = pv(e_mid_min_mev)``
    bounds the step midpoint states from below (the caller passes ``(1 - f) E_cut`` with ``f`` the
    step's energy-loss fraction limit: a step starts above ``E_cut`` and its midpoint loses less
    than ``f E``). The initial states ``p1v1`` of the run's particles (the source: no secondaries
    for z >= 2) lie in ``[pv(e_lo_mev) / (1 + m), pv(e_hi_mev) (1 + m)]`` with
    ``m = TABLE_PV_MARGIN``, ``e_lo_mev = max(E0 - 6 sigma_E, E_cut)`` and
    ``e_hi_mev = min(E0 + 6 sigma_E, table maximum)`` (the source sampler's largest Gaussian
    deviate is below 6: ``sqrt(-2 ln u)`` with ``u >= 2^-25``)."""
    m = projectile.mass_mev

    def pv(e: float) -> float:
        return e * (e + 2.0 * m) / (e + m)

    return (
        pv(e_mid_min_mev),
        pv(e_lo_mev) / (1.0 + TABLE_PV_MARGIN),
        pv(e_hi_mev) * (1.0 + TABLE_PV_MARGIN),
    )


@dataclass(frozen=True)
class TableCertificate:
    """Certificate of one ``pv_eq/pv`` table (:func:`certify_pv_eq_table`).

    ``certified_bound_analytic`` is the certified value: the largest per-cell analytic bound of
    ``|var_table / var_root - 1|`` (``bound_later_steps`` and ``bound_birth_step`` per step kind)
    over the ``n_cells`` solved cells, from ``n_refined_points`` root solves of the 2x-refined
    grid and the safety factor ``safety_factor``. ``observed_max_sampled`` is the largest error
    observed at the ``n_sampled_points`` interior samples (the empirical cross-check; never above
    the analytic bound of its cell). ``bound`` is the threshold it was checked against."""

    certified_bound_analytic: float
    bound_later_steps: float
    bound_birth_step: float
    observed_max_sampled: float
    n_cells: int
    n_refined_points: int
    n_sampled_points: int
    safety_factor: float
    bound: float

    def to_dict(self) -> dict[str, Any]:
        return dict(self.__dict__)


@functools.cache
def _table_sample_kernel() -> Any:
    """Warp CPU kernel (float64): the per-step variance of a state from the root solve
    (``mapped_state``) and from the table (``mapped_state_table``), for a later step and for a
    first step with ``pv_end = min(pv_mid, keep p1v1)`` (``keep = 1 - TABLE_MIN_LOSS``) and
    with ``pv_end = min(pv_end_lo, that)``; the transport functions of
    ``ionmc.physics.ion_mcs``."""
    from ionmc.physics.ion_mcs import make_ion_mcs  # noqa: PLC0415 (ion_mcs imports this module)

    im = make_ion_mcs(wp.float64)
    d = wp.float64

    @wp.kernel(module="unique")
    def sample(
        ion: wp.array2d(dtype=d),
        pro: wp.array2d(dtype=d),
        comp: wp.array2d(dtype=d),
        quad: wp.array(dtype=d),
        tab: wp.array2d(dtype=d),
        z: d,
        mass: d,
        pv_end_lo: d,
        keep: d,
        pv: wp.array(dtype=d),
        p1v1: wp.array(dtype=d),
        out: wp.array2d(dtype=d),
        status: wp.array2d(dtype=int),
    ):  # type: ignore[no-untyped-def]
        k = wp.tid()
        one = d(1.0)
        pr, p1r, sr = im.mapped_state(ion, pro, comp, quad, 0, z, mass, pv[k], p1v1[k])
        pt, p1t, st = im.mapped_state_table(tab, 0, pv[k], p1v1[k])
        for c in range(3):
            birth = int(0)
            pv_end = pv[k]
            if c > 0:
                birth = 1
                pv_end = wp.min(pv[k], keep * p1v1[k])
            if c == 2:
                pv_end = wp.min(pv_end_lo, pv_end)
            out[k, 2 * c] = im.variance_of_state(
                pr, p1r, sr, z, pv[k], pv_end, p1v1[k], birth, one, one, one
            )
            out[k, 2 * c + 1] = im.variance_of_state(
                pt, p1t, st, z, pv[k], pv_end, p1v1[k], birth, one, one, one
            )
        status[k, 0] = sr
        status[k, 1] = st

    return sample


def _band_cells(xs: np.ndarray, ys: np.ndarray, p1v1_min: float, p1v1_max: float) -> np.ndarray:
    """Mask ``[n_x - 1, n_y - 1]`` of the cells meeting the band
    ``ln p1v1_min <= x + y <= ln p1v1_max`` (``x + y`` is extreme at the corners)."""
    s_lo = xs[:-1, None] + ys[None, :-1]
    s_hi = xs[1:, None] + ys[None, 1:]
    return (s_hi >= math.log(p1v1_min) - 1e-12) & (s_lo <= math.log(p1v1_max) + 1e-12)


def _table_axes(row: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = row[:TABLE_HEADER]
    return np.linspace(h[0], h[6], int(h[2])), np.linspace(h[3], h[7], int(h[5]))


def _neighbour_max(v: np.ndarray, cells: np.ndarray) -> np.ndarray:
    """Per cell, the maximum of ``v`` over the cell and its (up to 8) solved neighbours."""
    ci, cj = cells[:, 0], cells[:, 1]
    g = np.zeros((ci.max() + 3, cj.max() + 3))
    g[ci + 1, cj + 1] = v
    out = np.zeros_like(v)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            out = np.maximum(out, g[ci + 1 + di, cj + 1 + dj])
    return out


def fdm_sensitivity_max(
    l1_lo: np.ndarray, l1_hi: np.ndarray, l2_lo: np.ndarray, l2_hi: np.ndarray
) -> np.ndarray:
    """Exact maximum of ``S = |c2 + c3 l1| / f(l1, l2)``, ``f = c0 + c1 l1 + c2 l2 + c3 l1 l2``
    (``FDM_COEFFICIENTS``), over the rectangles ``[l1_lo, l1_hi] x [l2_lo, l2_hi]`` (elementwise);
    ``inf`` where ``f <= 0`` at a corner or ``c2 + c3 l1`` changes sign. At fixed ``l2``, ``S`` is a
    ratio of two linear functions of ``l1`` (monotone where ``f > 0``); at fixed ``l1`` it is a
    constant over a linear function of ``l2`` (monotone); ``f`` is bilinear, so ``f > 0`` at the
    four corners implies ``f > 0`` on the rectangle, and the maximum is at a corner."""
    c0, c1, c2, c3 = FDM_COEFFICIENTS
    out = np.zeros(np.broadcast(l1_lo, l2_lo).shape)
    ok = np.ones(out.shape, dtype=bool)
    n_lo, n_hi = c2 + c3 * l1_lo, c2 + c3 * l1_hi
    ok &= np.sign(n_lo) == np.sign(n_hi)
    for l1, n in ((l1_lo, n_lo), (l1_hi, n_hi)):
        for l2 in (l2_lo, l2_hi):
            f = c0 + c1 * l1 + c2 * l2 + c3 * l1 * l2
            ok &= f > 0.0
            out = np.maximum(out, np.abs(n) / np.where(f > 0.0, f, 1.0))
    return np.where(ok, out, np.inf)


def certify_pv_eq_table(
    row: np.ndarray,
    ion_pack: np.ndarray,
    proton_pack: np.ndarray,
    comp: np.ndarray,
    quad: np.ndarray,
    z: float,
    mass_mev: float,
    p1v1_min: float,
    p1v1_max: float,
) -> TableCertificate:
    """Certificate of the table ``row`` (:func:`build_pv_eq_table`): a domain-wide analytic bound
    of the table's per-step variance error, per solved cell (the cells meeting the band of the
    initial states), with a sampled cross-check. Deterministic (float64). Checked against
    :data:`TABLE_VARIANCE_BOUND` by :func:`check_table_certificate`.

    **Interpolation error of q = pv_eq/pv** (the interpolated quantity). On a cell with spacings
    ``hx`` (``ln pv``) and ``hy`` (``ln r``, ``r = p1v1/pv``) the tensor-product bilinear
    interpolant satisfies ``|q - Pi q| <= hx^2/8 sup|q_xx| + hy^2/8 sup|q_yy|``. The suprema are
    estimated from second differences of the exact root solve on the 2x-refined grid of the cell
    (its 3 x 3 stencil: corners, edge midpoints and centre, all root-solved): ``hx^2/8 |q_xx|``
    ~ ``|q(0) - 2 q(1/2) + q(1)| / 2`` along each of the three x-lines (likewise in y), the maximum
    over the lines and over the stencils of the cell and its (up to 8) solved neighbours (a cell
    where a second derivative changes sign has a vanishing own difference), times the safety
    factor :data:`TABLE_SAFETY_FACTOR` = 2 (``q`` is a smooth, slowly varying function of
    ``ln pv`` and ``ln r``, P11; the factor covers the variation of the second derivative within
    a cell, and the sampled cross-check below falsifies it if it does not), plus
    :data:`TABLE_ROOT_TOL_REL` relative for the root solve's own convergence. This is ``dq``;
    ``q`` lies in ``[q_lo, q_hi]`` = stencil range widened by ``dq``, and
    ``|ln q_table - ln q_root| <= dq / q_lo``.

    **Variance error.** The table changes only ``pv_eq``; the prefactor ``(z E_s / pv_mid)^2 rho
    / X_S s`` is the same for both paths, so ``var_table / var_root = f_t / f_r`` with ``f`` the
    f_dM fit ``c0 + c1 l1 + c2 l2 + c3 l1 l2``, ``l2 = log10 pv_eq``:
    later steps ``l1 = log10(1 - (pv_eq / p1v1_eq)^2) = log10(1 - 1/r^2)`` (``p1v1_eq = r pv_eq``,
    independent of ``pv_eq``); the first step (birth, analytic memory average)
    ``l1 = log10(1 - (pv_end / p1v1)^2) - log10(e)``, independent of ``pv_eq``. ``f`` is linear in
    ``l2`` at fixed ``l1``, hence exactly ``f_t / f_r - 1 = (c2 + c3 l1) (l2_t - l2_r) / f_r`` and
    ``|var_t / var_r - 1| <= S_max |l2_t - l2_r|`` with ``S_max`` the maximum of
    ``|c2 + c3 l1| / f`` over the cell's ``(l1, l2)`` rectangle (:func:`fdm_sensitivity_max`,
    exact at its corners). The ratio depends on the birth step's end point through ``l1``
    (``c3 != 0``), so it is not endpoint-independent; the rectangle covers every physical end
    point: ``pv_end`` in ``(0, min(pv_mid, (1 - L) p1v1)]`` with ``L`` =
    :data:`TABLE_MIN_LOSS`, i.e. ``1 - (pv_end/p1v1)^2`` in
    ``[max(1 - 1/r^2, 1 - (1 - L)^2), 1)``, which includes ``pv_end`` at the cutoff and below it
    (straggling); later steps cover ``1 - 1/r^2`` in ``[max(1 - 1/r_lo^2, 1 - (1 - L)^2),
    max(1 - 1/r_hi^2, 1 - (1 - L)^2)]``. ``l2`` covers ``[q_lo pv_cell_min, q_hi pv_cell_max]``.
    Contract scope: states of the table (``r >= 1 + THIN_SLAB_R``) of a particle that has lost at
    least ``L`` of its initial ``pv`` (:data:`TABLE_MIN_LOSS`: before that the memory term drives
    the f_dM fit to its zero, where a relative error is not defined).

    **Sampled cross-check.** At :data:`TABLE_SAMPLE_N` x :data:`TABLE_SAMPLE_N` interior points
    per cell the variance from the table is compared with the root solve for a later step and a
    first step with ``pv_end = min(pv_mid, (1 - L) p1v1)`` and ``pv_end`` = the table's lowest
    ``pv`` (the cutoff side); later-step samples with ``pv > (1 - L) p1v1`` are outside the
    contract and skipped. Fails closed (``UnsupportedCombinationError``) when an observed error
    exceeds the analytic bound of its cell (the smoothness assumption is falsified), or when a
    root solve or a table lookup is invalid at a stencil or sample point."""
    row = np.ascontiguousarray(row, dtype=np.float64)
    xs, ys = _table_axes(row)
    n_y = len(ys)
    cells = np.argwhere(_band_cells(xs, ys, p1v1_min, p1v1_max))
    n_c = len(cells)
    if n_c == 0:
        raise fail("the pv_eq/pv table has no solved cell to certify")
    hx, hy = xs[1] - xs[0], ys[1] - ys[0]
    ci, cj = cells[:, 0], cells[:, 1]

    # (1) 3 x 3 stencil of every cell on the 2x-refined grid, all root-solved (shared points once)
    a = np.arange(3)
    ri, rj = np.broadcast_arrays(
        2 * ci[:, None, None] + a[None, :, None], 2 * cj[:, None, None] + a[None, None, :]
    )
    n_ry = 2 * n_y - 1
    uniq, inv = np.unique((ri * n_ry + rj).reshape(-1), return_inverse=True)
    pv_s = np.exp(xs[0] + (uniq // n_ry) * (0.5 * hx))
    p1_s = pv_s * np.exp(ys[0] + (uniq % n_ry) * (0.5 * hy))
    peq, st = root_pv_equivalent_batch(ion_pack, proton_pack, comp, quad, z, mass_mev, pv_s, p1_s)
    if np.any(st != 0):
        k = int(np.argwhere(st != 0)[0, 0])
        raise fail(
            "the pv_eq/pv table cannot be certified: the root solve fails at the stencil point "
            f"pv = {pv_s[k]:.6g} MeV, p1v1 = {p1_s[k]:.6g} MeV (status {st[k]}); use "
            "mcs_mapping_eval='root'"
        )
    q = (peq / pv_s)[inv].reshape(n_c, 3, 3)
    d_xx = _neighbour_max(np.abs(q[:, 0, :] - 2.0 * q[:, 1, :] + q[:, 2, :]).max(axis=1), cells)
    d_yy = _neighbour_max(np.abs(q[:, :, 0] - 2.0 * q[:, :, 1] + q[:, :, 2]).max(axis=1), cells)
    q_max = q.max(axis=(1, 2))
    dq = TABLE_SAFETY_FACTOR * 0.5 * (d_xx + d_yy) + TABLE_ROOT_TOL_REL * q_max
    q_lo = q.min(axis=(1, 2)) - dq
    q_hi = q_max + dq
    if np.any(q_lo <= 0.0):
        raise fail("the pv_eq/pv table cannot be certified: pv_eq/pv not bounded away from 0")
    dl2 = dq / (q_lo * math.log(10.0))

    # (2) sensitivity of the variance to log10 pv_eq over the cell (both step kinds)
    l2_lo = np.log10(q_lo) + xs[ci] / math.log(10.0)
    l2_hi = np.log10(q_hi) + xs[ci + 1] / math.log(10.0)
    om_min = -math.expm1(2.0 * math.log1p(-TABLE_MIN_LOSS))  # 1 - (1 - L)^2
    om_lo = np.maximum(-np.expm1(-2.0 * ys[cj]), om_min)
    om_hi = np.maximum(-np.expm1(-2.0 * ys[cj + 1]), om_min)
    s_later = fdm_sensitivity_max(np.log10(om_lo), np.log10(om_hi), l2_lo, l2_hi)
    s_birth = fdm_sensitivity_max(np.log10(om_lo) - LOG10_E, np.full(n_c, -LOG10_E), l2_lo, l2_hi)
    b_later = s_later * dl2
    b_birth = s_birth * dl2

    # (3) sampled cross-check: N x N interior points per cell
    off = (np.arange(TABLE_SAMPLE_N) + 0.5) / TABLE_SAMPLE_N
    ox, oy = (g.reshape(-1) for g in np.meshgrid(off, off, indexing="ij"))
    x = (xs[ci][:, None] + ox[None, :] * hx).reshape(-1)
    y = (ys[cj][:, None] + oy[None, :] * hy).reshape(-1)
    pv = np.exp(x)
    p1 = pv * np.exp(y)
    n = pv.size
    dev = "cpu"

    def a1(v: np.ndarray) -> Any:
        return wp.array(np.ascontiguousarray(v, dtype=np.float64), dtype=wp.float64, device=dev)

    def a2(v: np.ndarray) -> Any:
        return wp.array(
            np.ascontiguousarray(v, dtype=np.float64).reshape(1, -1), dtype=wp.float64, device=dev
        )

    out = wp.zeros((n, 6), dtype=wp.float64, device=dev)
    sst = wp.zeros((n, 2), dtype=int, device=dev)
    wp.launch(
        _table_sample_kernel(),
        dim=n,
        inputs=[
            a2(ion_pack), a2(proton_pack), a2(comp), a1(quad), a2(row), wp.float64(z),
            wp.float64(mass_mev), wp.float64(math.exp(xs[0])),
            wp.float64(1.0 - TABLE_MIN_LOSS), a1(pv), a1(p1), out, sst,
        ],
        device=dev,
    )  # fmt: skip
    v = out.numpy()
    s = sst.numpy()
    exact = v[:, 0::2]
    tab = v[:, 1::2]
    incl = np.ones(exact.shape, dtype=bool)
    # later steps: states that lost less than L of pv are outside the relative contract
    incl[:, 0] = p1 * (1.0 - TABLE_MIN_LOSS) >= pv
    bad = np.any(s != 0, axis=1) | np.any(incl & ~((exact > 0.0) & (tab > 0.0)), axis=1)
    if np.any(bad):
        k = int(np.argwhere(bad)[0, 0])
        raise fail(
            "the pv_eq/pv table cannot be certified: invalid mapping at the sample point "
            f"pv = {pv[k]:.6g} MeV, p1v1 = {p1[k]:.6g} MeV (root status {s[k, 0]}, table status "
            f"{s[k, 1]}); use mcs_mapping_eval='root'"
        )
    rel = np.where(incl, np.abs(tab / np.where(incl, exact, 1.0) - 1.0), 0.0)
    rel = rel.reshape(n_c, TABLE_SAMPLE_N**2, 3)
    o_later = rel[:, :, 0].max(axis=1)
    o_birth = rel[:, :, 1:].max(axis=(1, 2))
    over = (o_later > b_later) | (o_birth > b_birth)
    if np.any(over):
        k = int(np.argwhere(over)[0, 0])
        raise fail(
            "the pv_eq/pv table cannot be certified: the sampled error exceeds the analytic "
            f"bound in the cell at pv = {math.exp(xs[ci[k]]):.6g} MeV, p1v1/pv = "
            f"{math.exp(ys[cj[k]]):.6g} (observed {o_later[k]:.3e} / {o_birth[k]:.3e}, bound "
            f"{b_later[k]:.3e} / {b_birth[k]:.3e} for later / first steps): the smoothness "
            "assumption of the bound is falsified (a corrupted node?); use mcs_mapping_eval='root'"
        )
    return TableCertificate(
        certified_bound_analytic=float(max(b_later.max(), b_birth.max())),
        bound_later_steps=float(b_later.max()),
        bound_birth_step=float(b_birth.max()),
        observed_max_sampled=float(rel.max()),
        n_cells=n_c,
        n_refined_points=int(uniq.size),
        n_sampled_points=n,
        safety_factor=TABLE_SAFETY_FACTOR,
        bound=TABLE_VARIANCE_BOUND,
    )


def check_table_certificate(cert: TableCertificate) -> None:
    """Fail closed (``UnsupportedCombinationError``) when the certified analytic bound exceeds
    :data:`TABLE_VARIANCE_BOUND` (read at call time)."""
    if not cert.certified_bound_analytic <= TABLE_VARIANCE_BOUND:
        raise fail(
            "the pv_eq/pv table exceeds the decision 0043 interpolation bound: certified analytic "
            f"bound of |var_table / var_root - 1| = {cert.certified_bound_analytic:.3e} > "
            f"{TABLE_VARIANCE_BOUND:.1e} over {cert.n_cells} cells; use mcs_mapping_eval='root'"
        )


_TABLES: dict[str, tuple[np.ndarray, TableCertificate]] = {}


def build_pv_eq_table(
    ion_pack: np.ndarray,
    proton_pack: np.ndarray,
    comp: np.ndarray,
    quad: np.ndarray,
    z: float,
    mass_mev: float,
    pv_lo: float,
    p1v1_min: float,
    p1v1_max: float,
    n_pv: int = TABLE_N_PV,
    n_r: int = TABLE_N_R,
) -> np.ndarray:
    """The packed table row of :func:`build_certified_pv_eq_table` (certified, fails closed above
    the bound)."""
    return build_certified_pv_eq_table(
        ion_pack, proton_pack, comp, quad, z, mass_mev, pv_lo, p1v1_min, p1v1_max, n_pv, n_r
    )[0]


def build_certified_pv_eq_table(
    ion_pack: np.ndarray,
    proton_pack: np.ndarray,
    comp: np.ndarray,
    quad: np.ndarray,
    z: float,
    mass_mev: float,
    pv_lo: float,
    p1v1_min: float,
    p1v1_max: float,
    n_pv: int = TABLE_N_PV,
    n_r: int = TABLE_N_R,
) -> tuple[np.ndarray, TableCertificate]:
    """``(row, certificate)``: the packed ``pv_eq/pv`` table row of one ion in one material
    (float64)::

        [x0, inv_dx, n_x, y0, inv_dy, n_y, x_last, y_last, q[0, 0], q[0, 1], ..., q[n_x-1, n_y-1]]

    ``x = ln pv`` from ``ln pv_lo`` to ``ln(p1v1_max / (1 + THIN_SLAB_R))`` (``n_pv`` nodes) and
    ``y = ln(p1v1/pv)`` from ``ln(1 + THIN_SLAB_R)`` (the thin-slab guard of ``ion_mcs``) to
    ``ln(p1v1_max / pv_lo)`` (``n_r`` nodes), ``q = pv_eq/pv`` of the root solve at the node. Only
    the nodes of cells that meet the band of the run's initial states, ``p1v1_min <= p1v1 <=
    p1v1_max`` (:func:`table_domain`), are solved; every other node, and a node whose solve fails,
    holds :data:`TABLE_INVALID_NODE` (a state there returns status 5, fail closed). Fails closed
    (``UnsupportedCombinationError``) unless every cell meeting the band has four valid nodes.
    Every produced table is certified (:func:`certify_pv_eq_table`: analytic bilinear error
    bound from second differences on a 2x-refined root-solve grid with safety factor 2,
    cross-checked by 5 x 5 interior samples per cell) and fails closed when the certified bound
    exceeds :data:`TABLE_VARIANCE_BOUND` (checked on every call, also for a cached table).
    Cached per input."""
    from ionmc.physics.ion_mcs import THIN_SLAB_R  # noqa: PLC0415

    if not (0.0 < pv_lo < p1v1_min <= p1v1_max) or n_pv < 2 or n_r < 2:
        raise fail(
            f"invalid pv_eq table domain: pv_lo={pv_lo}, p1v1 in [{p1v1_min}, {p1v1_max}], "
            f"n_pv={n_pv}, n_r={n_r}"
        )
    import hashlib  # noqa: PLC0415

    h = hashlib.sha256()
    for a in (ion_pack, proton_pack, comp, quad):
        h.update(np.ascontiguousarray(a, dtype=np.float64).tobytes())
    h.update(
        np.array([z, mass_mev, pv_lo, p1v1_min, p1v1_max, n_pv, n_r], dtype=np.float64).tobytes()
    )
    key = h.hexdigest()
    if key in _TABLES:
        check_table_certificate(_TABLES[key][1])
        return _TABLES[key]
    y0 = math.log1p(THIN_SLAB_R)
    x0 = math.log(pv_lo)
    x1 = math.log(p1v1_max) - y0
    y1 = math.log(p1v1_max / pv_lo)
    xs = np.linspace(x0, x1, n_pv)
    ys = np.linspace(y0, y1, n_r)
    gx, gy = np.meshgrid(xs, ys, indexing="ij")
    pv = np.exp(gx).reshape(-1)
    p1 = pv * np.exp(gy).reshape(-1)
    q = np.full(pv.size, TABLE_INVALID_NODE)
    band = _band_cells(xs, ys, p1v1_min, p1v1_max)
    need = np.zeros((n_pv, n_r), dtype=bool)
    need[:-1, :-1] |= band
    need[1:, :-1] |= band
    need[:-1, 1:] |= band
    need[1:, 1:] |= band
    reach = need.reshape(-1)
    if np.any(reach):
        peq, st = root_pv_equivalent_batch(
            ion_pack, proton_pack, comp, quad, z, mass_mev, pv[reach], p1[reach]
        )
        q[reach] = np.where(st == 0, peq / pv[reach], TABLE_INVALID_NODE)
    grid = q.reshape(n_pv, n_r)
    cell_ok = (
        (grid[:-1, :-1] > 0.0)
        & (grid[1:, :-1] > 0.0)
        & (grid[:-1, 1:] > 0.0)
        & (grid[1:, 1:] > 0.0)
    )
    if np.any(band & ~cell_ok):
        i, j = np.argwhere(band & ~cell_ok)[0]
        raise fail(
            "the pv_eq/pv table cannot be built over the source domain: the root solve fails at a "
            f"node of the cell at pv = {math.exp(xs[i]):.6g} MeV, p1v1/pv = {math.exp(ys[j]):.6g} "
            "(the ion or proton path does not reach one cell beyond pv(E0 + 6 sigma_E)); use "
            "mcs_mapping_eval='root' or extend the stopping table"
        )
    header = [
        x0,
        (n_pv - 1) / (x1 - x0),
        float(n_pv),
        y0,
        (n_r - 1) / (y1 - y0),
        float(n_r),
        x1,
        y1,
    ]
    row = np.concatenate((header, q)).astype(np.float64)
    row.setflags(write=False)
    cert = certify_pv_eq_table(
        row, ion_pack, proton_pack, comp, quad, z, mass_mev, p1v1_min, p1v1_max
    )
    check_table_certificate(cert)
    _TABLES[key] = (row, cert)
    return row, cert
