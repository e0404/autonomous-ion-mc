# mypy: ignore-errors
# (Warp kernel annotations use runtime types; see _wpfunc.py.)
"""Scattering-power checks shared by the tests and ``validation/scripts/transport/mcs_checks.py``.

Everything here uses the engine's own shared Warp functions (``scattering_power_dm``,
``scattering_variance_birth``) so that the checks cannot drift from the physics code:

* U4: the quadrature of ``T_dM`` along the CSDA path against ``theta_Hanson (1 + dM %/100)`` of
  the table of Gottschalk, Med. Phys. 37 (2010) 352, parsed at run time from the LaTeX source of
  arXiv:0908.1413 (a hash-verified file in the git-ignored ``.ionmc-cache``; never committed);
* the comparison of ``X_S`` with the table ``tbl:LS`` of the same paper;
* U5: step independence of the stepped variance sum against the quadrature;
* a generalised-Highland cross-check (Gottschalk et al. 1993, Tsai radiation length).

Units: mm, MeV, g/cm2, rad.
"""

import hashlib
import io
import math
import re
import tarfile
from pathlib import Path

import numpy as np
import warp as wp

from ionmc.materials import N_A, Material
from ionmc.materials import WATER as WATER_MATERIAL
from ionmc.physics.em import make_em
from ionmc.physics.projectiles import PROTON
from ionmc.physics.scattering import (
    inverse_scattering_length_cm2_per_g,
    scattering_length_g_cm2,
)
from ionmc.transport.tables import TransportTables

F64 = wp.float64
EM = make_em(F64)
M_P = PROTON.mass_mev

GOTTSCHALK_SHA256 = "67fb1478e51534064f0e5363b1ee4160da2a37a600f51897b56dd74ea22ba5ed"
GOTTSCHALK_CACHE = Path(__file__).resolve().parents[3] / ".ionmc-cache" / "reference"
U4_FRACTIONS = (0.01, 0.1, 0.5, 0.9)
U4_TOLERANCE = 0.015
U4B_TOLERANCE = 0.045
"""Frozen U4b tolerance: the paper's own T_dM is within 2.74 % of theta_Hanson over the frozen
points (largest |dM %| of theta0Single for x/R1 >= 0.01) and U4 bounds our reproduction of T_dM
at 1.5 %, so (1 + 0.0274)(1 + 0.015) - 1 = 4.3 % -> 4.5 %."""
STEPS_MM = (0.01, 0.1, 0.5, 1.0, 2.0, 5.0)
DEPTH_FRACTIONS = (0.05, 0.25, 0.5, 0.9)
_ROW = re.compile(r"^\s*(" + r"\s*&\s*".join([r"(-?\d+\.\d+)"] * 11) + r")\s*(?:\\\\)?\s*$")


@wp.kernel(module="unique")
def _power_kernel(
    pv: wp.array(dtype=wp.float64),  # type: ignore[valid-type]
    p1v1: wp.float64,
    inv_xs: wp.float64,
    rho: wp.float64,
    out: wp.array(dtype=wp.float64),  # type: ignore[valid-type]
) -> None:
    i = wp.tid()
    out[i] = EM.scattering_power_dm(pv[i], p1v1, wp.float64(1.0), inv_xs, rho)


def pv_np(e: np.ndarray) -> np.ndarray:
    return e * (e + 2.0 * M_P) / (e + M_P)


class ProtonPath:
    """CSDA path of a proton in one material: power T_dM [rad2/mm] at depth x [mm]."""

    def __init__(self, tables: TransportTables, material: Material, e0: float) -> None:
        self.t = tables
        self.rho = material.density_g_cm3
        self.inv_xs = inverse_scattering_length_cm2_per_g(material)
        self.e0 = e0
        self.r1_g = tables.range_g_cm2(0, e0)
        self.r1_mm = self.r1_g * 10.0 / self.rho
        self.p1v1 = float(pv_np(np.array([e0]))[0])

    def energy(self, x_mm: np.ndarray) -> np.ndarray:
        t = self.t
        lnr = np.log(np.maximum(self.r1_g - self.rho * x_mm / 10.0, 1e-30))
        grid = t.ln_r0[0] + np.arange(t.n_r) / t.inv_dln_r[0]
        return np.exp(np.interp(lnr, grid, t.ln_e_of_r[0]))

    def power(self, x_mm: np.ndarray) -> np.ndarray:
        pv = pv_np(self.energy(x_mm))
        out = wp.zeros(len(pv), dtype=wp.float64, device="cpu")
        wp.launch(
            _power_kernel,
            dim=len(pv),
            inputs=[
                wp.array(pv, dtype=wp.float64, device="cpu"),
                F64(self.p1v1),
                F64(self.inv_xs),
                F64(self.rho),
                out,
            ],
            device="cpu",
        )
        res: np.ndarray = out.numpy()
        return res

    def theta2_quadrature(self, x_mm: float) -> float:
        """Integral of T over [0, x] by Gauss-Legendre in s with x' = x s^4 (log singularity)."""
        gx, gw = np.polynomial.legendre.leggauss(200)
        s = 0.5 * (gx + 1.0)
        pts = x_mm * s**4
        jac = 4.0 * x_mm * s**3 * 0.5 * gw
        return float(np.sum(jac * self.power(pts)))

    def birth_variance(self, s_mm: float) -> float:
        """Variance of the birth step from the shared ``scattering_variance_birth``."""
        e_mid = self.energy(np.array([0.5 * s_mm]))
        e_end = self.energy(np.array([s_mm]))
        return float(
            EM.scattering_variance_birth(
                F64(float(pv_np(e_mid)[0])),
                F64(float(pv_np(e_end)[0])),
                F64(self.p1v1),
                F64(1.0),
                F64(self.inv_xs),
                F64(self.rho),
                F64(s_mm),
            )
        )

    def theta2_steps(self, x_mm: float, step_mm: float) -> float:
        """Sum over uniform steps (last one shortened) as in the transport: the linearized
        analytic log-average variance (residual ~1e-3) on the first step,
        ``s T(E_mid)`` on the others."""
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


_ALPHA = 1.0 / 137.035999084
_R_E = 2.8179403262e-13
_TSAI_L = {1: (5.31, 6.144), 2: (4.79, 5.621), 3: (4.74, 5.805), 4: (4.71, 5.924)}


def radiation_length_g_cm2(material: Material) -> float:
    """Radiation length from Tsai's formula as given by the PDG (independent of the engine)."""
    inv = 0.0
    for el, w in material._fractions:
        a = _ALPHA * el.Z
        coulomb = a * a * (1 / (1 + a * a) + 0.20206 - 0.0369 * a**2 + 0.0083 * a**4 - 0.002 * a**6)
        if el.Z in _TSAI_L:
            lrad, lrad_p = _TSAI_L[el.Z]
        else:
            lrad = math.log(184.15 * el.Z ** (-1 / 3))
            lrad_p = math.log(1194.0 * el.Z ** (-2 / 3))
        inv += (
            w
            * 4
            * _ALPHA
            * _R_E**2
            * N_A
            / el.A_g_mol
            * (el.Z**2 * (lrad - coulomb) + el.Z * lrad_p)
        )
    return 1.0 / inv


def generalised_highland_theta2(path: ProtonPath, x_mm: float, x0_g_cm2: float) -> float:
    """Generalised Highland (Gottschalk et al. 1993): projected theta_0^2 with the pv of the
    CSDA path, ``(14.1 MeV)^2 [1 + log10(x/X0)/9]^2 int dx / (X0 pv^2)`` (x in g/cm2)."""
    gx, gw = np.polynomial.legendre.leggauss(200)
    s = 0.5 * (gx + 1.0)
    pts = x_mm * s**4
    jac = 4.0 * x_mm * s**3 * 0.5 * gw
    pv = pv_np(path.energy(pts))
    x_g = x_mm * path.rho / 10.0
    integral = float(np.sum(jac * (path.rho / 10.0) / (x0_g_cm2 * pv**2)))
    return (14.1 * (1.0 + math.log10(x_g / x0_g_cm2) / 9.0)) ** 2 * integral


# ---------------------------------------------------------------------------------------------
def gottschalk_tex(cache_dir: Path | None = None) -> str:
    """LaTeX source of arXiv:0908.1413 from the cache; raises ``FileNotFoundError`` if the file
    is missing and ``ValueError`` if its SHA-256 differs from the pinned value."""
    path = (cache_dir or GOTTSCHALK_CACHE) / GOTTSCHALK_SHA256
    if not path.is_file():
        raise FileNotFoundError(f"Gottschalk source not in the cache: {path}")
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != GOTTSCHALK_SHA256:
        raise ValueError(f"cached file {path} does not have sha256 {GOTTSCHALK_SHA256}")
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as tar:
        member = next(m for m in tar.getmembers() if m.name.endswith("ScatPowerV2.tex"))
        handle = tar.extractfile(member)
        if handle is None:
            raise ValueError("ScatPowerV2.tex is not a regular file")
        return handle.read().decode("utf-8", errors="replace")


def parse_tables(
    tex: str,
) -> tuple[dict[str, dict[float, tuple[float, float]]], dict[str, float], dict[str, float]]:
    """Parse ``tbl:theta0Single`` (rows ``x/R1 & rho x & E & theta_Hanson & ... & dM %``) and
    ``tbl:LS`` (rho X_S). Returns ({material: {x/R1: (theta_Hanson [mrad], dM %)}},
    {material: rho R1}, {material: rho X_S})."""
    start = tex.index("%%% Theta0 table:")
    end = tex.index("tbl:theta0Single", start)
    rows: dict[str, dict[float, tuple[float, float]]] = {}
    rho_r1: dict[str, float] = {}
    material = ""
    for line in tex[start:end].splitlines():
        head = re.search(r"158\.60\\,MeV p on (\w+), \$\\rho R_1=([\d.]+)\$", line)
        if head:
            material = head.group(1)
            rho_r1[material] = float(head.group(2))
            rows[material] = {}
            continue
        m = _ROW.match(line)
        if m and material:
            v = [float(x) for x in re.split(r"\s*&\s*", m.group(1))]
            rows[material][v[0]] = (v[3], v[10])
    ls = re.search(r"\\rho X_S\$\\quad\(g/cm\$\^2\$\)&([^\\]*)\\\\", tex)
    if ls is None:
        raise ValueError("table tbl:LS not found")
    values = [float(x) for x in ls.group(1).split("&")]
    x_s = dict(zip(["Be", "Lexan", "H2O", "Al", "Cu", "Pb"], values, strict=True))
    return rows, rho_r1, x_s


def u4_deviations(
    tables: TransportTables,
    material: Material,
    rows: dict[float, tuple[float, float]],
) -> tuple[dict[float, float], float]:
    """theta_dM(quadrature) / (theta_Hanson (1 + dM %/100)) - 1 per x/R1, and our rho R1."""
    path = ProtonPath(tables, material, 158.6)
    out = {}
    for frac in U4_FRACTIONS:
        theta_hanson, dm_percent = rows[frac]
        reference = theta_hanson * (1.0 + dm_percent / 100.0)
        out[frac] = 1e3 * math.sqrt(path.theta2_quadrature(frac * path.r1_mm)) / reference - 1.0
    return out, path.r1_g


def u4b_deviations(
    tables: TransportTables,
    material: Material,
    rows: dict[float, tuple[float, float]],
) -> tuple[dict[float, float], float]:
    """theta_dM(quadrature) / theta_Hanson - 1 per x/R1 (the measured/Moliere-theory column,
    independent of the paper's dM fit) and the paper's max |dM %| over the same points."""
    path = ProtonPath(tables, material, 158.6)
    out = {}
    for frac in U4_FRACTIONS:
        theta_hanson, _ = rows[frac]
        out[frac] = 1e3 * math.sqrt(path.theta2_quadrature(frac * path.r1_mm)) / theta_hanson - 1.0
    return out, max(abs(rows[f][1]) for f in U4_FRACTIONS)


def xs_deviations(x_s_table: dict[str, float], materials: dict[str, Material]) -> dict[str, float]:
    """Our X_S over the tabulated one minus 1, per material key of the table."""
    return {k: scattering_length_g_cm2(m) / x_s_table[k] - 1.0 for k, m in materials.items()}


def u5_deviations(tables: TransportTables) -> dict[tuple[float, float], float]:
    """Relative deviation of the stepped sum from the quadrature, keyed by (x/R1, s [mm]);
    water, 150 MeV."""
    path = ProtonPath(tables, WATER_MATERIAL, 150.0)
    out = {}
    for frac in DEPTH_FRACTIONS:
        x = frac * path.r1_mm
        exact = path.theta2_quadrature(x)
        for s in STEPS_MM:
            out[(frac, s)] = path.theta2_steps(x, s) / exact - 1.0
    return out


def u5_negative_control(tables: TransportTables) -> dict[str, float]:
    """Relative spread (max - min) / min over the steps of the summed theta^2 at 0.5 R1 for the
    per-step Highland formula (all steps) and for the engine's stepped T_dM (steps <= 1 mm)."""
    path = ProtonPath(tables, WATER_MATERIAL, 150.0)
    x0 = radiation_length_g_cm2(WATER_MATERIAL)
    x = 0.5 * path.r1_mm
    totals = []
    for s in STEPS_MM:
        n = int(round(x / s))
        step = x / n
        mids = (np.arange(n) + 0.5) * step
        pv = pv_np(path.energy(mids))
        step_g = step * path.rho / 10.0
        theta2 = (14.1 / pv) ** 2 * (step_g / x0) * (1.0 + np.log10(step_g / x0) / 9.0) ** 2
        totals.append(float(np.sum(theta2)))
    tdm = [path.theta2_steps(x, s) for s in (0.01, 0.1, 0.5, 1.0)]
    return {
        "per_step_highland_relative_spread": (max(totals) - min(totals)) / min(totals),
        "stepped_tdm_relative_spread_up_to_1mm": (max(tdm) - min(tdm)) / min(tdm),
    }


def highland_cross_check(tables: TransportTables, material: Material) -> dict[float, float]:
    """theta_dM(quadrature) / theta_generalised-Highland - 1 at 158.6 MeV per frozen x/R1."""
    path = ProtonPath(tables, material, 158.6)
    x0 = radiation_length_g_cm2(material)
    out = {}
    for frac in U4_FRACTIONS:
        x = frac * path.r1_mm
        out[frac] = (
            math.sqrt(path.theta2_quadrature(x) / generalised_highland_theta2(path, x, x0)) - 1.0
        )
    return out
