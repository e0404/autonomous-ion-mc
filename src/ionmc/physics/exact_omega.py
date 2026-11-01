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

Evaluation (the P11 record declares the per-step root solve; it is the choice implemented, no
tabulation): :func:`pv_equivalent` solves ``g(u) = b_p(e^u, r e^u) - b_ion(pv, p1v1) = 0`` for
``u = ln pv_eq``, ``r = p1v1 / pv``, by a bracketed Illinois (modified regula falsi) iteration
with a bisection safeguard. The bracket is the proton path itself: from its lowest energy to the
largest ``pv_eq`` whose ``p1v1_eq = r pv_eq`` is still on the proton path. ``g`` is monotone
increasing in ``pv_eq`` at fixed ``r`` (P11), the bracket stays valid, and the iteration is
capped at :data:`MAX_ITERATIONS`; there is no failure mode other than the returned status.
Status codes (:data:`STATUS_NAMES`): 0 converged, 1 invalid ion state (outside the ion path or
``pv >= p1v1``), 2 no sign change in the bracket (``b_ion`` outside the proton ``b`` range at this
``r``), 3 the proton path is too short for this ``r``, 4 not converged within the iteration cap.
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
}

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

    return SimpleNamespace(
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
    (or the mapping of any scanned state fails). ``ion_tables`` maps a material name to the ion's
    stopping table (default: Bethe tables, see :meth:`ExactOmegaMapping.build`)."""
    if projectile == PROTON:
        return PROTON_FDM_STATUS
    if not materials:
        raise fail("at least one material is required")
    e_hi = e0_mev + 6.0 * energy_sigma_mev
    nominal = 0.0
    upper = 0.0
    for mat in materials:
        tab = None if ion_tables is None else ion_tables.get(mat.name)
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
