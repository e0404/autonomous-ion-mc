# mypy: ignore-errors
# (Warp kernel annotations use runtime types; see _wpfunc.py.)
"""Ion counterpart of :class:`ionmc.transport.mcs_checks.ProtonPath`, shared by the V3-008 CI
rows (``tests/ionmc/test_v3_008_rows_ci.py``: I10, I11-map) and the verification steps
(``validation/scripts/transport/steps_v8.py``: I9, I12, I15 and the lateral grid of ``lv8``).

:class:`IonPath` is the CSDA path of an ion in one material, with the projected scattering
power ``T_ion`` [rad^2/mm] evaluated by the engine's own shared Warp functions
(``make_ion_mcs``: the exact-Omega mapped state and ``mcs_variance``), so the reference cannot
drift from the physics code. The scattering-power inputs are the ones the transport uses:
``1/X_S`` is the material value ``inverse_scattering_length_cm2_per_g(material)`` (the value the
transport loads into ``TransportTables.inv_rho_xs_cm2_g``) and ``rho`` is the material density
``material.density_g_cm3`` (the voxel density of a homogeneous phantom). Both stay explicit
keyword parameters so that ratio-only checks can pass any value; their defaults are the material
values, so the absolute scale of ``power``, ``theta2_quadrature`` and the Fermi-Eyges A2 built on
it is the implemented one.

Units: mm, MeV, g/cm2, rad.
"""

from __future__ import annotations

import math

import numpy as np
import warp as wp

from ionmc.materials import WATER, Material
from ionmc.physics import exact_omega as XO
from ionmc.physics.ion_mcs import THIN_SLAB_R, make_ion_mcs
from ionmc.physics.projectiles import Projectile
from ionmc.physics.scattering import inverse_scattering_length_cm2_per_g
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.transport.tables import TransportTables

D = wp.float64


def ion_tables(proj: Projectile, material: Material = WATER) -> TransportTables:
    """Transport tables of ``proj`` in ``material`` with the default Bethe stopping source."""
    return TransportTables.from_stopping_tables(
        [BetheStoppingSource().table(material, proj)], projectile=proj
    )


def mapping_arrays(
    proj: Projectile, tables: TransportTables | None = None, material: Material = WATER
) -> tuple:
    """Packed exact-Omega arrays ``(ion, proton, composition, quadrature)`` of one material."""
    t = tables or ion_tables(proj, material)
    ion = XO.pack_transport_row(t, 0)[None, :]
    pro = XO.proton_path(material).pack[None, :]
    comp = XO.pack_composition(material)[None, :]
    return ion, pro, comp, XO.pack_quadrature()


def mcs_kernel():
    """Warp CPU kernel: per state the step variance (``mcs_variance``), the exact-Omega residual
    ``b_p(pv_eq, r pv_eq) - b_ion(pv, p1v1)`` and ``pv_eq``; ``inv_rho_xs`` [cm2/g] and ``rho``
    [g/cm3] are kernel arguments (no literal)."""
    IM = make_ion_mcs(D)
    XOF = XO.make_exact_omega(D)

    @wp.kernel(module="unique")
    def k(
        ion: wp.array2d(dtype=D),
        pro: wp.array2d(dtype=D),
        comp: wp.array2d(dtype=D),
        quad: wp.array(dtype=D),
        z: D,
        mass: D,
        pv: wp.array(dtype=D),
        pv_end: wp.array(dtype=D),
        p1v1: D,
        birth: int,
        inv_rho_xs: D,
        rho: D,
        s_mm: wp.array(dtype=D),
        out: wp.array2d(dtype=D),
        status: wp.array(dtype=int),
    ):
        i = wp.tid()
        v, st = IM.mcs_variance(
            ion, pro, comp, quad, 0, z, mass, pv[i], pv_end[i], p1v1, birth, inv_rho_xs, rho,
            s_mm[i],
        )  # fmt: skip
        out[i, 0] = v
        peq, p1eq, st2 = IM.mapped_state(ion, pro, comp, quad, 0, z, mass, pv[i], p1v1)
        pv_map = wp.min(pv[i], p1v1 / (D(1.0) + D(THIN_SLAB_R)))
        target, sb = XOF.b_slab(ion[0], comp[0], quad, z, mass, pv_map, p1v1)
        g, sg = XOF.b_residual(pro[0], comp[0], quad, wp.log(peq), p1v1 / pv_map, target)
        out[i, 1] = g
        out[i, 2] = peq
        status[i] = st + st2 + sb + sg

    return k


class IonPath:
    """CSDA path of an ion in one material (default water): the implemented projected
    scattering power ``T_ion`` [rad^2/mm] (``power``), the stepped and the quadrature theta^2
    with the shared functions.

    ``inv_rho_xs_cm2_g`` defaults to ``inverse_scattering_length_cm2_per_g(material)`` and
    ``rho_g_cm3`` to ``material.density_g_cm3`` (the transport's values)."""

    def __init__(
        self,
        proj: Projectile,
        e_u: float,
        *,
        material: Material = WATER,
        inv_rho_xs_cm2_g: float | None = None,
        rho_g_cm3: float | None = None,
    ) -> None:
        self.proj = proj
        self.material = material
        self.inv_rho_xs_cm2_g = (
            inverse_scattering_length_cm2_per_g(material)
            if inv_rho_xs_cm2_g is None
            else float(inv_rho_xs_cm2_g)
        )
        self.rho_g_cm3 = material.density_g_cm3 if rho_g_cm3 is None else float(rho_g_cm3)
        self.t = ion_tables(proj, material)
        self.arrs = mapping_arrays(proj, self.t, material)
        self.e0 = e_u * proj.a
        self.r1 = self.t.range_g_cm2(0, self.e0)
        self.r1_mm = self.r1 * 10.0 / self.rho_g_cm3
        self.p1v1 = self._pv(np.array([self.e0]))[0]
        self.kern = mcs_kernel()
        self.last_status = np.zeros(0, dtype=int)

    def _pv(self, e: np.ndarray) -> np.ndarray:
        m = self.proj.mass_mev
        return e * (e + 2.0 * m) / (e + m)

    def energy(self, x_mm: np.ndarray) -> np.ndarray:
        t = self.t
        lnr = np.log(np.maximum(self.r1 - self.rho_g_cm3 * x_mm / 10.0, 1e-30))
        grid = t.ln_r0[0] + np.arange(t.n_r) / t.inv_dln_r[0]
        return np.exp(np.interp(lnr, grid, t.ln_e_of_r[0]))

    def _run(self, pv: np.ndarray, pv_end: np.ndarray, s: np.ndarray, birth: int) -> np.ndarray:
        n = len(pv)
        out = wp.zeros((n, 3), dtype=D, device="cpu")
        st = wp.zeros(n, dtype=int, device="cpu")
        ion, pro, comp, quad = self.arrs
        wp.launch(
            self.kern, dim=n,
            inputs=[wp.array(ion, dtype=D, device="cpu"), wp.array(pro, dtype=D, device="cpu"),
                    wp.array(comp, dtype=D, device="cpu"), wp.array(quad, dtype=D, device="cpu"),
                    D(float(self.proj.z)), D(self.proj.mass_mev),
                    wp.array(pv, dtype=D, device="cpu"), wp.array(pv_end, dtype=D, device="cpu"),
                    D(float(self.p1v1)), birth, D(self.inv_rho_xs_cm2_g), D(self.rho_g_cm3),
                    wp.array(s, dtype=D, device="cpu"), out, st],
            device="cpu",
        )  # fmt: skip
        self.last_status = st.numpy()
        res: np.ndarray = out.numpy()
        return res

    def power(self, x_mm: np.ndarray) -> np.ndarray:
        """``T_ion`` [rad^2/mm] at depths ``x_mm`` (the variance of a unit step, non-birth)."""
        pv = self._pv(self.energy(x_mm))
        return self._run(pv, pv, np.ones_like(pv), 0)[:, 0]

    def birth_variance(self, s_mm: float) -> float:
        pv = self._pv(self.energy(np.array([0.5 * s_mm])))
        pe = self._pv(self.energy(np.array([s_mm])))
        return float(self._run(pv, pe, np.array([s_mm]), 1)[0, 0])

    def theta2_quadrature(self, x_mm: float) -> float:
        """Integral of ``T_ion`` over [0, x] by Gauss-Legendre in s with x' = x s^4."""
        gx, gw = np.polynomial.legendre.leggauss(200)
        s = 0.5 * (gx + 1.0)
        pts = x_mm * s**4
        jac = 4.0 * x_mm * s**3 * 0.5 * gw
        return float(np.sum(jac * self.power(pts)))

    def theta2_steps(self, x_mm: float, step_mm: float) -> float:
        """Sum over uniform steps (last one shortened) as in the transport: the birth-step
        variance on the first step, ``s T(E_mid)`` on the others."""
        n = int(math.floor(x_mm / step_mm))
        lens = np.full(n, step_mm)
        starts = np.arange(n) * step_mm
        rem = x_mm - n * step_mm
        if rem > 0.0:
            lens = np.append(lens, rem)
            starts = np.append(starts, n * step_mm)
        total = self.birth_variance(float(lens[0]))
        if len(lens) > 1:
            total += float(np.sum(lens[1:] * self.power(starts[1:] + 0.5 * lens[1:])))
        return total
