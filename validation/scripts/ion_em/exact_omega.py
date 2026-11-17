"""Exact-Omega mapping of the proton MCS model to light ions (V3-008, check P11; decision 0043 §3).

Pure functions, numpy only (the report driver imports ``ionmc`` lazily for the stopping tables).
Every formula is that of the committed P1 reference ``moliere_hanson.py`` (Moliere chi_c, chi_a with
the Born screening factor, Fano's delta_b, the 1/(pv)^2 thick-slab average, project range tables).

Model (decision 0043 §3, "Selected (r5)", quoted)::

    T_ion(pv, p1v1; material) = f_dM(pv_eq, p1v1_eq) (z E_s / pv)^2 / X_S
    p1v1_eq / pv_eq = p1v1 / pv          (the memory term 1 - (pv/p1v1)^2 is unchanged)
    pv_eq is the root of  b_p(pv_eq, p1v1_eq; material) = b_ion(pv, p1v1; z, M, material)
    b = ln Omega + delta_b_Fano   (Moliere collision number of the slab traversed from p1v1 down to
                                   pv along the CSDA path, with the Fano term evaluated at each
                                   projectile's own velocities)

Hence ``p1v1_eq`` is NOT "the proton with the ion's initial velocity": it is fixed by the decision
as ``p1v1_eq = pv_eq * (p1v1 / pv)``, so the memory argument ``1 - (pv/p1v1)^2`` of f_dM is the
ion's own. The proton slab is the one from p1v1_eq down to pv_eq on the proton CSDA path.

b of a slab (kinetic energies from ``pv``: T = [-(2M - pv) + sqrt((2M - pv)^2 + 4 M pv)]/2)::

    t       = R(E_1) - R(E)                       (range table of the projectile, per-u energies)
    chi_c^2 = K z^2 sum_i w_i Z_i^2/A_i  int_0^t dt'/(pv)^2          (``electrons="fano"``)
    <ln chi_a^2> = int ln chi_a^2 dt'/(pv)^2 / int dt'/(pv)^2
    delta_b = int delta_b(beta) dt'/(pv)^2 / int dt'/(pv)^2        (Fano; 0 for ``"bethe"``)
    b = ln chi_c^2 - <ln chi_a^2> - ln 1.167 + delta_b

Evaluation choice: the default is the per-step root solve (:func:`pv_eq`, Brent on ln pv_eq, the
bracket being the proton table; fail-closed without a sign change). The allowed alternative, a
(pv, p1v1) tabulation of pv_eq/pv with bilinear interpolation, is :class:`PvEqTable`; the test
``test_tabulation_*`` bounds its f_dM (equivalently T) error; the root solve is what P9 binds unless
that test passes and the tabulation is chosen explicitly.

Note on the domain of the proton equivalent: matching b of a heavy ion needs a proton slab of far
larger collision number than the equal-velocity proton (b grows by about ln(A*/F_scr) = 1 to 3), so
the proton-equivalent state has pv_eq/pv = 0.13 to 0.69 but a proton p1v1_eq of up to 0.8 GeV
(He), 2.5 GeV (C, O) kinetic energy on the I11 grid (depth fractions 0.02 to 0.97 included). The
proton path used for the root is therefore built on an extended Bethe table
(``PROTON_EMAX_PER_U`` = 8000 MeV), identical formulas to the project table.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

# --- the committed P1 reference (loaded by path; the tests do the same) --------------------------
_MH_PATH = Path(__file__).resolve().parent / "moliere_hanson.py"
if "moliere_hanson" in sys.modules:
    MH: Any = sys.modules["moliere_hanson"]
else:
    _spec = importlib.util.spec_from_file_location("moliere_hanson", _MH_PATH)
    assert _spec is not None and _spec.loader is not None
    MH = importlib.util.module_from_spec(_spec)
    sys.modules["moliere_hanson"] = MH
    _spec.loader.exec_module(MH)

Element = tuple[int, float, float]

PROTON_EMAX_PER_U = 8000.0
"""Upper end [MeV/u] of the extended proton table used for the proton-equivalent state."""
N_SLAB = 64
"""Gauss-Legendre nodes of the slab integrals inside b (``MH`` uses 400; see the P11 note)."""
N_THETA = 400
"""Gauss-Legendre nodes of the theta_0 integral along the path (the value of ``MH``)."""
ROOT_XTOL = 1e-15
"""Relative tolerance of the Brent iteration in ln(pv_eq)."""
RESIDUAL_TOL = 1e-10
"""Plan I11-map: |b_p(pv_eq) - b_ion| <= 1e-10 (absolute)."""


# --- projectile path ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Path_:
    """A projectile (charge ``z``, rest energy ``mass_mev``, mass number ``a``) with its CSDA range
    relation: ``range_of(E_per_u) -> g/cm^2`` and ``energy_of(R) -> E_per_u`` (exact inverses)."""

    z: float
    mass_mev: float
    a: float
    range_of: Callable[[Any], Any]
    energy_of: Callable[[Any], Any]
    e_min_per_u: float
    e_max_per_u: float


def path_from_table(table: Any) -> Path_:
    """:class:`Path_` of an ``ionmc.physics.stopping.StoppingTable`` (MeV/u, log-log)."""
    proj = table.projectile
    return Path_(
        z=float(proj.z),
        mass_mev=float(proj.mass_mev),
        a=float(proj.a),
        range_of=lambda e: np.asarray(table.range_at(e)),
        energy_of=lambda r: np.asarray(table.energy_from_range(r)),
        e_min_per_u=float(table.energy_per_u[0]),
        e_max_per_u=float(table.energy_per_u[-1]),
    )


def bethe_path(name: str, e_max_per_u: float | None = None, material: Any = None) -> Path_:
    """Path of a project projectile in ``material`` (default water) from the Bethe table."""
    from ionmc.materials import WATER as W  # noqa: PLC0415
    from ionmc.physics.projectiles import PROJECTILES  # noqa: PLC0415
    from ionmc.physics.stopping import BetheStoppingSource  # noqa: PLC0415

    src = (
        BetheStoppingSource()
        if e_max_per_u is None
        else BetheStoppingSource(e_max_per_u=e_max_per_u)
    )
    return path_from_table(src.table(W if material is None else material, PROJECTILES[name]))


def scaled_path(base: Path_, z: float, a_star: float) -> Path_:
    """Synthetic ion on a proton path with the exact scaling R_ion(E_u) = (A*/z^2) R_p(E_u) and
    M = A* m_p (``base`` must be the proton). Used only to test the exact Omega identity."""
    return Path_(
        z=z,
        mass_mev=a_star * MH.M_P_MEV,
        a=a_star,
        range_of=lambda e: a_star / (z * z) * base.range_of(e),
        energy_of=lambda r: base.energy_of(r * z * z / a_star),
        e_min_per_u=base.e_min_per_u,
        e_max_per_u=base.e_max_per_u,
    )


# --- kinematics --------------------------------------------------------------------------------
def kinetic_from_pv(pv: float, m_mev: float) -> float:
    """T [MeV] with T (T + 2M)/(T + M) = pv  (positive root of T^2 + (2M - pv) T - pv M = 0)."""
    q = 2.0 * m_mev - pv
    root = math.sqrt(q * q + 4.0 * pv * m_mev)
    if q >= 0.0:
        return 2.0 * pv * m_mev / (q + root)  # no cancellation for q >= 0
    return 0.5 * (root - q)  # q < 0: -q + root has no cancellation


def pv_of_path(path: Path_, e_per_u: float) -> float:
    """p v [MeV] of the projectile at kinetic energy ``e_per_u`` * a."""
    return float(MH.pv_mev(e_per_u * path.a, path.mass_mev))


# --- b of a slab ---------------------------------------------------------------------------------
def b_slab(
    path: Path_,
    pv: float,
    p1v1: float,
    elements: Sequence[Element] = MH.WATER,
    electrons: str = "fano",
    n: int = N_SLAB,
) -> float:
    """b = ln Omega + delta_b_Fano of the slab from p1v1 (initial) down to pv (current) along the
    CSDA path of ``path`` (module docstring). Requires ``pv < p1v1``; the reference definitions."""
    if not 0.0 < pv < p1v1:
        raise ValueError(f"need 0 < pv < p1v1, got pv={pv}, p1v1={p1v1}")
    e = kinetic_from_pv(pv, path.mass_mev) / path.a
    e1 = kinetic_from_pv(p1v1, path.mass_mev) / path.a
    if not (path.e_min_per_u * (1 - 1e-12) <= e < e1 <= path.e_max_per_u * (1 + 1e-12)):
        raise ValueError(f"energies {e}, {e1} MeV/u outside the path table or not ordered")
    r_e = float(path.range_of(np.float64(e)))
    r_1 = float(path.range_of(np.float64(e1)))
    t = r_1 - r_e
    if not t > 0.0:
        raise ValueError("slab thickness not positive")
    nodes, weights = MH._gauss_nodes(t, n)
    # depth t' from the p1v1 end; the range remaining at t' is R(E) + (t - t') (no cancellation)
    rem = r_e + (t - nodes)
    en = np.asarray(path.energy_of(rem)) * path.a
    inv_pv2 = 1.0 / MH.pv_mev(en, path.mass_mev) ** 2
    integral = float(np.sum(weights * inv_pv2))
    zsum = MH.z2_sum(elements) if electrons == "fano" else MH.zz1_sum(elements)
    chi_c2 = MH.K_CHI_C * path.z**2 * zsum * integral
    lna = (
        float(
            np.sum(weights * inv_pv2 * MH.ln_chi_a2(path.z, elements, en, path.mass_mev, electrons))
        )
        / integral
    )
    delta_b = 0.0
    if electrons == "fano":
        db = MH.fano_delta_b(elements, MH.beta2(en, path.mass_mev))
        delta_b = float(np.sum(weights * inv_pv2 * db)) / integral
    return math.log(chi_c2) - lna - math.log(MH.OMEGA_FACTOR) + delta_b


def b_proton(
    proton: Path_,
    pv: float,
    p1v1: float,
    elements: Sequence[Element] = MH.WATER,
    electrons: str = "fano",
    n: int = N_SLAB,
) -> float:
    """b_p(pv, p1v1; material): :func:`b_slab` of the proton path (z = 1)."""
    return b_slab(proton, pv, p1v1, elements, electrons, n)


def b_ion(
    ion: Path_,
    pv: float,
    p1v1: float,
    elements: Sequence[Element] = MH.WATER,
    electrons: str = "fano",
    n: int = N_SLAB,
) -> float:
    """b_ion(pv, p1v1; z, M, material): :func:`b_slab` of the ion path (charge z, mass M)."""
    return b_slab(ion, pv, p1v1, elements, electrons, n)


# --- the root ------------------------------------------------------------------------------------
def brent(
    f: Callable[[float], float], lo: float, hi: float, xtol: float = ROOT_XTOL, maxit: int = 200
) -> float:
    """Brent's root finder on [lo, hi]; fail-closed (ValueError) without a sign change."""
    fa, fb = f(lo), f(hi)
    if fa == 0.0:
        return lo
    if fb == 0.0:
        return hi
    if fa * fb > 0.0:
        raise ValueError(f"no bracket: f({lo}) = {fa}, f({hi}) = {fb}")
    a, b, c = lo, hi, lo
    fc = fa
    d = e = b - a
    for _ in range(maxit):
        if fb * fc > 0.0:
            c, fc = a, fa
            d = e = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol = 2.0 * 2.220446049250313e-16 * abs(b) + 0.5 * xtol * max(1.0, abs(b))
        m = 0.5 * (c - b)
        if abs(m) <= tol or fb == 0.0:
            return b
        if abs(e) >= tol and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p, q = 2.0 * m * s, 1.0 - s
            else:
                q0, r = fa / fc, fb / fc
                p = s * (2.0 * m * q0 * (q0 - r) - (b - a) * (r - 1.0))
                q = (q0 - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0.0:
                q = -q
            p = abs(p)
            if 2.0 * p < min(3.0 * m * q - abs(tol * q), abs(e * q)):
                e, d = d, p / q
            else:
                d = e = m
        else:
            d = e = m
        a, fa = b, fb
        b += d if abs(d) > tol else math.copysign(tol, m)
        fb = f(b)
    raise RuntimeError("Brent iteration did not converge")


def pv_eq(
    ion: Path_,
    proton: Path_,
    pv: float,
    p1v1: float,
    elements: Sequence[Element] = MH.WATER,
    electrons: str = "fano",
    n: int = N_SLAB,
) -> float:
    """Root pv_eq of  b_p(pv_eq, r pv_eq) = b_ion(pv, p1v1),  r = p1v1/pv  (decision 0043 §3).

    Brent on u = ln pv_eq. Bracket: from the proton table's lowest energy (offset by 1e-9, thin
    slab, b_p small) to the largest pv_eq whose p1v1_eq = r pv_eq is still in the table. Raises
    ValueError if b_ion is not inside the bracket (fail-closed). The ion's own slab must lie in
    the ion's table.
    """
    r = p1v1 / pv
    if not r > 1.0:
        raise ValueError("pv must be below p1v1")
    target = b_ion(ion, pv, p1v1, elements, electrons, n)
    pv_lo = pv_of_path(proton, proton.e_min_per_u * (1.0 + 1e-6))
    pv_hi = pv_of_path(proton, proton.e_max_per_u * (1.0 - 1e-7)) / r
    if not pv_hi > pv_lo:
        raise ValueError("proton table too small for this p1v1/pv")

    def g(u: float) -> float:
        p = math.exp(u)
        return b_slab(proton, p, r * p, elements, electrons, n) - target

    if (
        p1v1 <= pv_of_path(proton, proton.e_max_per_u)
        and pv >= pv_lo
        and b_slab(proton, pv, p1v1, elements, electrons, n) == target
    ):
        return pv  # the ion is the proton of this table (z = 1, M = m_p): identity, bit-exact

    return math.exp(brent(g, math.log(pv_lo), math.log(pv_hi)))


def residual(
    ion: Path_,
    proton: Path_,
    pv: float,
    p1v1: float,
    pv_eq_value: float,
    elements: Sequence[Element] = MH.WATER,
    electrons: str = "fano",
    n: int = N_SLAB,
) -> float:
    """b_p(pv_eq, r pv_eq) - b_ion(pv, p1v1) (I11-map root residual)."""
    r = p1v1 / pv
    return b_slab(proton, pv_eq_value, r * pv_eq_value, elements, electrons, n) - b_ion(
        ion, pv, p1v1, elements, electrons, n
    )


# --- tabulation (allowed alternative) ------------------------------------------------------------
class PvEqTable:
    """Tabulation of pv_eq/pv on a (ln pv, ln(p1v1/pv)) grid, bilinear interpolation.

    The decision allows a (pv, p1v1) grid read by bilinear interpolation in (ln pv, ln p1v1); the
    second axis here is ln(p1v1/pv) = ln p1v1 - ln pv, a shear of the same plane, chosen so that the
    unphysical region p1v1 < pv is excluded and the diagonal p1v1 -> pv is a grid edge.
    """

    def __init__(
        self,
        ion: Path_,
        proton: Path_,
        pv_range: tuple[float, float],
        ln_r_range: tuple[float, float],
        n_pv: int = 40,
        n_r: int = 40,
        elements: Sequence[Element] = MH.WATER,
        electrons: str = "fano",
    ) -> None:
        self.x = np.linspace(math.log(pv_range[0]), math.log(pv_range[1]), n_pv)
        self.y = np.linspace(ln_r_range[0], ln_r_range[1], n_r)
        self.ratio = np.empty((n_pv, n_r))
        for i, lx in enumerate(self.x):
            for j, ly in enumerate(self.y):
                pv = math.exp(lx)
                self.ratio[i, j] = (
                    pv_eq(ion, proton, pv, pv * math.exp(ly), elements, electrons) / pv
                )

    def __call__(self, pv: float, p1v1: float) -> float:
        """Interpolated pv_eq."""
        lx, ly = math.log(pv), math.log(p1v1 / pv)
        i = int(np.clip(np.searchsorted(self.x, lx) - 1, 0, len(self.x) - 2))
        j = int(np.clip(np.searchsorted(self.y, ly) - 1, 0, len(self.y) - 2))
        tx = (lx - self.x[i]) / (self.x[i + 1] - self.x[i])
        ty = (ly - self.y[j]) / (self.y[j + 1] - self.y[j])
        rr = self.ratio
        val = (
            (1 - tx) * (1 - ty) * rr[i, j]
            + tx * (1 - ty) * rr[i + 1, j]
            + (1 - tx) * ty * rr[i, j + 1]
            + tx * ty * rr[i + 1, j + 1]
        )
        return float(val) * pv


# --- theta_0 of the model -----------------------------------------------------------------------
def f_dm_equivalent(
    ion: Path_,
    proton: Path_,
    pv: float,
    p1v1: float,
    elements: Sequence[Element] = MH.WATER,
    electrons: str = "fano",
    n: int = N_SLAB,
    mapping: Callable[[float, float], float] | None = None,
) -> float:
    """f_dM(pv_eq, p1v1_eq) with p1v1_eq = pv_eq p1v1/pv: the proton model at the equivalent state.

    ``mapping(pv, p1v1) -> pv_eq`` replaces the root solve (tabulation variant)."""
    if mapping is None:
        peq = pv_eq(ion, proton, pv, p1v1, elements, electrons, n)
    else:
        peq = mapping(pv, p1v1)
    return float(MH.fdm(np.array([peq]), peq * p1v1 / pv)[0])


def theta0_from_power(
    path: Path_,
    elements: Sequence[Element],
    e1_per_u: float,
    x_over_r: float,
    f_of_state: Callable[[float, float], float],
    n: int = N_THETA,
) -> float:
    """theta_0 = sqrt(int T dt) with T = f(pv, p1v1) (z E_s/pv)^2 / X_S along the CSDA path of
    ``path`` from E_1 = ``e1_per_u`` to the depth ``x_over_r`` R(E_1) (the Gauss nodes of
    ``moliere_hanson.theta0_candidate``)."""
    r1 = float(path.range_of(np.float64(e1_per_u)))
    t = x_over_r * r1
    nodes, weights = MH._gauss_nodes(t, n)
    en = np.asarray(path.energy_of(r1 - nodes)) * path.a
    pvs = MH.pv_mev(en, path.mass_mev)
    p1v1 = float(MH.pv_mev(e1_per_u * path.a, path.mass_mev))
    inv_xs = MH.inverse_scattering_length(elements)
    total = 0.0
    for w, p in zip(weights, pvs, strict=True):
        total += w * f_of_state(float(p), p1v1) * (path.z * MH.E_S_MEV / float(p)) ** 2 * inv_xs
    return math.sqrt(total)


def theta0_ion(
    ion: Path_,
    proton: Path_,
    e_u: float,
    x_over_r: float,
    elements: Sequence[Element] = MH.WATER,
    electrons: str = "fano",
    n: int = N_THETA,
    n_slab: int = N_SLAB,
    mapping: Callable[[float, float], float] | None = None,
) -> float:
    """Projected theta_0 [rad] of the exact-Omega model: ion prefactor z/(pv) times the proton f_dM
    at (pv_eq, p1v1_eq), integrated along the ion CSDA path (E_1 = ``e_u`` MeV/u)."""

    def f(pv: float, p1v1: float) -> float:
        if not pv < p1v1 * (1.0 - 1e-12):
            return 0.0  # f_dM = 0 at pv = p1v1 (clamped fit); the node never sits there
        return f_dm_equivalent(ion, proton, pv, p1v1, elements, electrons, n_slab, mapping)

    return theta0_from_power(ion, elements, e_u, x_over_r, f, n)


def theta0_proton_model(
    proton: Path_,
    e_u: float,
    x_over_r: float,
    elements: Sequence[Element] = MH.WATER,
    n: int = N_THETA,
) -> float:
    """theta_0 of the unchanged proton model (f_dM(pv, p1v1), z = 1) on the proton path."""

    def f(pv: float, p1v1: float) -> float:
        return float(MH.fdm(np.array([pv]), p1v1)[0])

    return theta0_from_power(proton, elements, e_u, x_over_r, f, n)


def moliere_hanson_theta0(
    path: Path_,
    e_u: float,
    x_over_r: float,
    elements: Sequence[Element] = MH.WATER,
    electrons: str = "fano",
) -> float:
    """Molière-Hanson theta_0 of the slab to x/R (reference module, same path as the model)."""
    r1 = float(path.range_of(np.float64(e_u)))

    def energy_at(t: NDArray[np.float64]) -> NDArray[np.float64]:
        rem = np.maximum(r1 - t, 0.0)
        return np.asarray(path.energy_of(rem)) * path.a

    return float(
        MH.theta0_hanson(
            path.z, path.mass_mev, elements, x_over_r * r1, energy_at, electrons=electrons
        )
    )


# --- I11 report driver ----------------------------------------
def equal_velocity_proton_energy(ion: Path_, e_u: float) -> float:
    """Proton kinetic energy [MeV] with the velocity of the ion at ``e_u`` MeV/u.

    beta equal  <=>  T/M equal:  E_p = E_u A m_p / M  (A mass number, M the ion mass as used in
    ``ion``). With the proton path's per-nucleon axis (a = 1) this number is also its ``e_per_u``.
    """
    return e_u * ion.a * MH.M_P_MEV / ion.mass_mev


def beta_of_kinetic(t_mev: float, m_mev: float) -> float:
    """beta = sqrt(1 - (M/(T+M))^2) of kinetic energy T [MeV] and mass M [MeV]."""
    return math.sqrt(1.0 - (m_mev / (t_mev + m_mev)) ** 2)


def i11_tables(
    fractions: Sequence[float] = MH.I11_FRACTIONS,
    grid: dict[str, tuple[float, ...]] | None = None,
    electrons: str = "fano",
    n: int = N_THETA,
    n_slab: int = N_SLAB,
) -> list[dict[str, float | str]]:
    """Report-only I11-map / I11-rel / I11-abs ratios on the exact I11 grid (water).

    Per (ion, E [MeV/u], x/R) row:
    ``map_residual``: |b_p(pv_eq) - b_ion| at the grid point's (pv, p1v1) = (the ion's pv at the
    depth, p1v1) of the slab end;
    ``pv_eq_over_pv``: pv_eq/pv at that state;
    ``rel``: (theta_ion,model/theta_p,model) / (theta_MH,ion/theta_MH,p) - 1 at EQUAL VELOCITY (the
    proton at E_p = E_u A m_p / M, see :func:`equal_velocity_proton_energy`) and equal x/R;
    ``abs``: theta_ion,model/theta_MH,ion - 1 (informative).
    """
    grid = MH.I11_GRID_MEV_U if grid is None else grid
    proton = bethe_path("proton", PROTON_EMAX_PER_U)
    rows: list[dict[str, float | str]] = []
    for name, energies in grid.items():
        ion = bethe_path(name)
        for e_u in energies:
            for frac in fractions:
                th_ion = theta0_ion(ion, proton, e_u, frac, MH.WATER, electrons, n, n_slab)
                e_p = equal_velocity_proton_energy(ion, e_u)
                th_p = theta0_proton_model(proton, e_p, frac, MH.WATER, n)
                mh_ion = moliere_hanson_theta0(ion, e_u, frac, MH.WATER, electrons)
                mh_p = moliere_hanson_theta0(proton, e_p, frac, MH.WATER, electrons)
                # I11-map state: end of the slab (x/R of the grid point)
                r1 = float(ion.range_of(np.float64(e_u)))
                e_end = float(ion.energy_of(np.float64((1.0 - frac) * r1)))
                pv_end = pv_of_path(ion, e_end)
                p1v1 = pv_of_path(ion, e_u)
                peq = pv_eq(ion, proton, pv_end, p1v1, MH.WATER, electrons, n_slab)
                res = residual(ion, proton, pv_end, p1v1, peq, MH.WATER, electrons, n_slab)
                rows.append(
                    {
                        "ion": name,
                        "e_u": e_u,
                        "x_over_r": frac,
                        "pv_eq_over_pv": peq / pv_end,
                        "map_residual": abs(res),
                        "theta_model": th_ion,
                        "theta_p_model": th_p,
                        "theta_mh": mh_ion,
                        "theta_mh_p": mh_p,
                        "rel": (th_ion / th_p) / (mh_ion / mh_p) - 1.0,
                        "abs": th_ion / mh_ion - 1.0,
                    }
                )
    return rows


def main() -> None:
    rows = i11_tables()
    for r in rows:
        print(
            f"{r['ion']:9s} {float(r['e_u']):5.0f} {float(r['x_over_r']):4.2f} "
            f"pveq/pv={float(r['pv_eq_over_pv']):7.3f} res={float(r['map_residual']):.1e} "
            f"rel={100 * float(r['rel']):+6.2f}% abs={100 * float(r['abs']):+6.2f}%"
        )
    for name in MH.I11_GRID_MEV_U:
        sel = [abs(float(r["rel"])) for r in rows if r["ion"] == name]
        print(f"max |rel| {name}: {100 * max(sel):.3f} %")


if __name__ == "__main__":
    main()
