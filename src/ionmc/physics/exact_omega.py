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
nodes are computed by this same root solve (:func:`build_pv_eq_table`). The table is an
accelerator only (Codex REVIEW-ee53867c): its bilinear value (:func:`pv_equivalent_table`) is the
starting guess of the same root solve, started from a tight bracket around it
(:func:`pv_equivalent_seeded`; the full bracket when the tight one shows no sign change), so every
step returns the exact root to the root path's tolerance (:data:`TABLE_EVAL_CONTRACT`). The table's
interpolation error is a performance diagnostic (:func:`estimate_pv_eq_table`). The root solve:
:func:`pv_equivalent` solves ``g(u) = b_p(e^u, r e^u) - b_ion(pv, p1v1) = 0`` for
``u = ln pv_eq``, ``r = p1v1 / pv``, by a bracketed Illinois (modified regula falsi) iteration
with a bisection safeguard. The bracket is the proton path itself: from its lowest energy to the
largest ``pv_eq`` whose ``p1v1_eq = r pv_eq`` is still on the proton path. ``g`` is monotone
increasing in ``pv_eq`` at fixed ``r`` (P11), the bracket stays valid, and the iteration is capped
at :data:`MAX_ITERATIONS`; there is no failure mode other than the returned status.
Status codes (:data:`STATUS_NAMES`): 0 converged, 1 invalid ion state (outside the ion path or
``pv >= p1v1``), 2 no sign change in the bracket (``b_ion`` outside the proton ``b`` range at this
``r``), 3 the proton path is too short for this ``r``, 4 not converged within the iteration cap,
5 (:func:`pv_equivalent_table` only) the state lies outside the table or in a cell with an
invalid node (no starting guess; the seeded solve then uses the full bracket).
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
from ionmc.physics.em import make_em
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
"""Default nodes of the ``ln pv`` axis of the ``pv_eq/pv`` table."""
TABLE_N_R = 128
"""Default nodes of the ``ln(p1v1/pv)`` axis of the table."""
TABLE_INVALID_NODE = -1.0
"""Value stored at a node whose root solve failed (outside the paths); a cell touching it returns
status 5."""
TABLE_EDGE_TOL = {"float64": 1.0e-9, "float32": 1.0e-4}
"""Tolerance (in grid-index units) of the table edges."""
TABLE_PV_MARGIN = 1.0e-3
"""Relative margin of the table's largest ``p1v1`` above ``pv(E0 + 6 sigma_E)``."""
TABLE_HEADER = 9
"""Header length of a packed table row (:func:`build_pv_eq_table`); entry 8 is the row's starting
bracket half-width."""
TABLE_BRACKET_FACTOR = 10.0
"""Starting bracket half-width of the table-accelerated root, in units of the table's build-time
interpolation estimate (:func:`estimate_pv_eq_table`)."""
TABLE_BRACKET_MIN = 1.0e-5
"""Smallest starting bracket half-width (relative in ``pv_eq``, i.e. in ``ln pv_eq``)."""
TABLE_SAFETY_FACTOR = 2.0
"""Safety factor of the build-time interpolation estimate (:func:`estimate_pv_eq_table`; it
only sets the starting bracket, a performance quantity)."""
TABLE_ROOT_TOL_REL = 1.0e-10
"""Relative allowance for the root solve's convergence in the interpolation estimate."""
TABLE_EVAL_CONTRACT = (
    "mapping evaluation = exact root solve on every step (Illinois iteration of b_p(pv_eq, "
    "pv_eq r) = b_ion(pv, p1v1) to |residual| <= 1e-13 or a step/bracket <= 1e-14 in ln pv_eq, "
    "float64, the root path's tolerance and status codes); the pv_eq/pv table only supplies the "
    "starting bracket (half-width max(10 x the table's interpolation estimate, 1e-5) around the "
    "bilinear guess, widened once x10, else the full proton-path bracket, counted as "
    "mcs_table_bracket_miss); the interpolation error of the table is a performance diagnostic "
    "(build-time estimate recorded), not a correctness claim"
)
"""The evaluation contract of ``mcs_mapping_eval = "table"``, as recorded in the summary and the
capability report."""

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
    def illinois(
        pro_pth: wp.array(dtype=real),
        comp: wp.array(dtype=real),
        quad: wp.array(dtype=real),
        r: real,
        target: real,
        ua_in: real,
        ub_in: real,
        ga_in: real,
        gb_in: real,
    ) -> tuple[real, int, int]:
        """``(u, status, iterations)``: the Illinois iteration (modified regula falsi, bisection
        step whenever the secant point leaves the bracket) for ``g(u) = 0`` on the valid bracket
        ``[ua_in, ub_in]`` with ``g(ua_in) = ga_in <= 0 <= g(ub_in) = gb_in``. Converged (status
        0) when ``|g| <= RESIDUAL_TOL``, the bracket or the last step is below ``STEP_TOL``;
        status 3 when a slab fails, 4 when not converged within ``MAX_ITERATIONS``."""
        ua = ua_in
        ub = ub_in
        ga = ga_in
        gb = gb_in
        status = int(0)
        iters = int(0)
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
                elif wp.abs(gc) <= res_tol or ub - ua <= step_tol or step <= step_tol:
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
        return uc, status, iters

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
        (:func:`illinois`) on the full bracket of the proton path."""
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
                        uc, status, iters = illinois(pro_pth, comp, quad, r, target, ua, ub, ga, gb)
                        out = wp.exp(uc)
        return out, status, iters

    @named_func(name)
    def pv_equivalent_seeded(
        ion_pth: wp.array(dtype=real),
        pro_pth: wp.array(dtype=real),
        comp: wp.array(dtype=real),
        quad: wp.array(dtype=real),
        z: real,
        m_ion_mev: real,
        pv: real,
        p1v1: real,
        q0: real,
        w: real,
    ) -> tuple[real, int, int, int]:
        """``(pv_eq, status, iterations, bracket_miss)``: the same root as :func:`pv_equivalent`
        (the same equation, iteration, tolerances and status codes) started from a tight bracket
        around a starting guess ``pv_eq ~ q0 pv`` (the table). One end is the guess
        ``u0 = ln(q0 pv)``, the other ``u0 + w`` or ``u0 - w`` on the side where the residual
        changes sign (``w`` a relative half-width), clipped to the proton path's full bracket.
        Without a sign change the other end is moved once to ``u0 -+ 10 w``; failing that (or for
        ``q0 <= 0``, no guess) the full bracket of :func:`pv_equivalent` is used.
        ``bracket_miss`` = 1 whenever the INITIAL tight bracket had no sign change (the widened
        bracket or the full bracket was needed); it is a performance diagnostic, not a failure.
        Ion states only (z >= 2)."""
        out = pv
        status = int(0)
        iters = int(0)
        miss = int(1)
        target, st_ion = b_slab(ion_pth, comp, quad, z, m_ion_mev, pv, p1v1)
        if st_ion != 0:
            status = 1
            miss = 0
        else:
            r = p1v1 / pv
            pv_lo = kin.pv_mev(pro_pth[6] * (real(1.0) + lo_margin), m_p)
            pv_hi = kin.pv_mev(pro_pth[7] * (real(1.0) - hi_margin), m_p) / r
            if not pv_hi > pv_lo:
                status = 3
                miss = 0
            else:
                u_lo = wp.log(pv_lo)
                u_hi = wp.log(pv_hi)
                ua = u_lo
                ub = u_hi
                ga = real(0.0)
                gb = real(0.0)
                found = int(0)
                if q0 > real(0.0):
                    # one end at the guess u0, the other on the side of the sign change
                    u0 = wp.min(wp.max(wp.log(q0 * pv), u_lo), u_hi)
                    g0, s0 = b_residual(pro_pth, comp, quad, u0, r, target)
                    if s0 == 0:
                        hw = w
                        for k in range(2):
                            if found == 0:
                                if g0 <= real(0.0):
                                    b = wp.min(u0 + hw, u_hi)
                                    if b > u0:
                                        g2, s2 = b_residual(pro_pth, comp, quad, b, r, target)
                                        if s2 == 0 and g2 >= real(0.0):
                                            found = 1
                                            miss = k
                                            ua = u0
                                            ub = b
                                            ga = g0
                                            gb = g2
                                else:
                                    a = wp.max(u0 - hw, u_lo)
                                    if a < u0:
                                        g1, s1 = b_residual(pro_pth, comp, quad, a, r, target)
                                        if s1 == 0 and g1 <= real(0.0):
                                            found = 1
                                            miss = k
                                            ua = a
                                            ub = u0
                                            ga = g1
                                            gb = g0
                                hw = real(10.0) * hw
                if found == 0:
                    ga, sta = b_residual(pro_pth, comp, quad, ua, r, target)
                    gb, stb = b_residual(pro_pth, comp, quad, ub, r, target)
                    if sta != 0 or stb != 0:
                        status = 3
                    elif ga > real(0.0) or gb < real(0.0):
                        status = 2
                    else:
                        found = 1
                if found == 1:
                    uc, status, iters = illinois(pro_pth, comp, quad, r, target, ua, ub, ga, gb)
                    out = wp.exp(uc)
        return out, status, iters, miss

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
    table_header = wp.constant(TABLE_HEADER)

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
            base = table_header + i * n_y + j
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
        illinois=illinois,
        pv_equivalent=pv_equivalent,
        pv_equivalent_seeded=pv_equivalent_seeded,
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
class TableEstimate:
    """Build-time interpolation estimate of one ``pv_eq/pv`` table (:func:`estimate_pv_eq_table`):
    ``interp_estimate`` the largest per-cell estimate of ``|ln pv_eq(table) - ln pv_eq(root)|``
    over the ``n_cells`` solved cells (``n_refined_points`` stencil root solves, safety factor
    ``safety_factor``) and ``bracket_half_width`` the starting bracket half-width it sets. A
    performance diagnostic: every step is solved exactly (module docstring)."""

    interp_estimate: float
    bracket_half_width: float
    n_cells: int
    n_refined_points: int
    safety_factor: float

    def to_dict(self) -> dict[str, Any]:
        return dict(self.__dict__)


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


def estimate_pv_eq_table(
    row: np.ndarray,
    ion_pack: np.ndarray,
    proton_pack: np.ndarray,
    comp: np.ndarray,
    quad: np.ndarray,
    z: float,
    mass_mev: float,
    p1v1_min: float,
    p1v1_max: float,
) -> TableEstimate:
    """Build-time interpolation estimate of the table ``row`` (a performance diagnostic that sets
    the starting bracket of the table-accelerated root; not a correctness claim).

    On a cell with spacings ``hx`` (``ln pv``) and ``hy`` (``ln r``) the bilinear interpolant of
    ``q = pv_eq/pv`` deviates by ``<= hx^2/8 sup|q_xx| + hy^2/8 sup|q_yy|``. The second
    derivatives are estimated from second differences of the root solve on the 2x-refined grid of
    the cell (its 3 x 3 stencil, all root-solved): ``|q(0) - 2 q(1/2) + q(1)| / 2`` along each of
    the three lines per axis, the maximum over the lines and over the cell and its solved
    neighbours, times :data:`TABLE_SAFETY_FACTOR`, plus :data:`TABLE_ROOT_TOL_REL` relative:
    ``dq``. The estimate of ``|ln q_table - ln q_root|`` is ``dq / q_lo`` (``q_lo`` the stencil
    minimum less ``dq``), its maximum over the solved cells is ``interp_estimate``, and the
    starting bracket half-width ``max(TABLE_BRACKET_FACTOR x interp_estimate,
    TABLE_BRACKET_MIN)``. Fails closed (``UnsupportedCombinationError``) when a root solve fails
    at a stencil point."""
    row = np.ascontiguousarray(row, dtype=np.float64)
    xs, ys = _table_axes(row)
    n_y = len(ys)
    cells = np.argwhere(_band_cells(xs, ys, p1v1_min, p1v1_max))
    n_c = len(cells)
    if n_c == 0:
        raise fail("the pv_eq/pv table has no solved cell")
    hx, hy = xs[1] - xs[0], ys[1] - ys[0]
    ci, cj = cells[:, 0], cells[:, 1]
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
            "the pv_eq/pv table cannot be built: the root solve fails at the stencil point "
            f"pv = {pv_s[k]:.6g} MeV, p1v1 = {p1_s[k]:.6g} MeV (status {st[k]}); use "
            "mcs_mapping_eval='root'"
        )
    q = (peq / pv_s)[inv].reshape(n_c, 3, 3)
    d_xx = _neighbour_max(np.abs(q[:, 0, :] - 2.0 * q[:, 1, :] + q[:, 2, :]).max(axis=1), cells)
    d_yy = _neighbour_max(np.abs(q[:, :, 0] - 2.0 * q[:, :, 1] + q[:, :, 2]).max(axis=1), cells)
    q_max = q.max(axis=(1, 2))
    dq = TABLE_SAFETY_FACTOR * 0.5 * (d_xx + d_yy) + TABLE_ROOT_TOL_REL * q_max
    q_lo = q.min(axis=(1, 2)) - dq
    if np.any(q_lo <= 0.0):
        raise fail("the pv_eq/pv table cannot be built: pv_eq/pv not bounded away from 0")
    est = float(np.max(dq / q_lo))
    return TableEstimate(
        interp_estimate=est,
        bracket_half_width=max(TABLE_BRACKET_FACTOR * est, TABLE_BRACKET_MIN),
        n_cells=n_c,
        n_refined_points=int(uniq.size),
        safety_factor=TABLE_SAFETY_FACTOR,
    )


_TABLES: dict[str, tuple[np.ndarray, TableEstimate]] = {}


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
    """The packed table row of :func:`build_pv_eq_table_and_estimate`."""
    return build_pv_eq_table_and_estimate(
        ion_pack, proton_pack, comp, quad, z, mass_mev, pv_lo, p1v1_min, p1v1_max, n_pv, n_r
    )[0]


def build_pv_eq_table_and_estimate(
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
) -> tuple[np.ndarray, TableEstimate]:
    """``(row, estimate)``: the packed ``pv_eq/pv`` table row of one ion in one material
    (float64)::

        [x0, inv_dx, n_x, y0, inv_dy, n_y, x_last, y_last, w, q[0, 0], ..., q[n_x-1, n_y-1]]

    ``x = ln pv`` from ``ln pv_lo`` to ``ln(p1v1_max / (1 + THIN_SLAB_R))`` (``n_pv`` nodes) and
    ``y = ln(p1v1/pv)`` from ``ln(1 + THIN_SLAB_R)`` (the thin-slab guard of ``ion_mcs``) to
    ``ln(p1v1_max / pv_lo)`` (``n_r`` nodes), ``q = pv_eq/pv`` of the root solve at the node and
    ``w`` the starting bracket half-width (:func:`estimate_pv_eq_table`). Only the nodes of cells
    that meet the band of the run's initial states, ``p1v1_min <= p1v1 <= p1v1_max``
    (:func:`table_domain`), are solved; every other node, and a node whose solve fails, holds
    :data:`TABLE_INVALID_NODE` (no starting guess there: the full bracket is used). Fails closed
    (``UnsupportedCombinationError``) unless every cell meeting the band has four valid nodes.
    The table only accelerates the root (:data:`TABLE_EVAL_CONTRACT`). Cached per input."""
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
        0.0,
    ]
    row = np.concatenate((header, q)).astype(np.float64)
    est = estimate_pv_eq_table(
        row, ion_pack, proton_pack, comp, quad, z, mass_mev, p1v1_min, p1v1_max
    )
    row[8] = est.bracket_half_width
    row.setflags(write=False)
    _TABLES[key] = (row, est)
    return row, est
