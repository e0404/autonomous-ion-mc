"""Row V6 / V6-sens (Gottschalk 2015, arXiv:1409.1938): data parser, beam construction, E0 and the
comparator (plan Amendment 14 (g), Amendment 17 (c)).

This is the ONLY module that reads the measured dose table (``tbl:dmlg``): the beam builder, the E0
procedure and every fit never see it (isolation test). The e-print is read only from the
hash-verified data cache (registry id ``arxiv-1409-1938-source``); the tarball hash is verified by the
data layer and the member ``halo.tex`` against :data:`HALO_TEX_SHA256`; a mismatch is recorded in the
document and fails the row. Roles (plan (g)): ``tbl:dmlg`` evaluation (measured); ``tbl:parms``
(p1, p2, p5-p7) construction for the beam only, never for the dose comparison; the row is
"evaluation with shared beam lineage (calibration-type inputs)".

Observable (Amendment 17 (c)1): at every data point (r_i, z_j) the dose is the overlap-volume-weighted
mean over the 0.4 cm diameter sphere centred at the point of the dose [MeV/g/primary] of a 0.5 mm x
0.5 mm axisymmetric (r, z) annular tally. The overlap weights are computed by 16 x 16 x 16 midpoint
sub-samples per cell in (rho, phi, z) and recorded (shape and sha256).
"""

from __future__ import annotations

import hashlib
import io
import math
import re
import tarfile
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

DATASET_ID = "arxiv-1409-1938-source"
HALO_TEX_SHA256 = "e5a270e79c5a6d0e91f6c06936892a34d691c22b72f28ebfe6f0af7f8a5023d6"
TARBALL_SHA256 = "92cdaffbc4237e5178e8a83b3ceee6e15f708c2f4ebbd16aff63f249dc990ac2"
LINEAGE = "evaluation with shared beam lineage (calibration-type inputs)"
BEAM_ROLE = "construction for the beam only (paper's 25-parameter MD fit to the same tbl:dmlg data)"
DATA_ROLE = "evaluation (measured)"

RADII_CM = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10)
DEPTHS_CM = (
    1.533, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
    18.5, 19.0, 19.2, 19.4, 19.5, 19.6, 19.8, 20.0, 20.2, 20.4, 20.6, 20.8, 21.0, 21.2, 21.4, 21.5, 21.6,
    21.8, 22.0, 22.2, 22.4, 22.5, 22.6, 22.8, 23.0, 23.5, 24.0, 25.0,
)  # fmt: skip
"""Measurement depths [cm] of ``tbl:dmlg`` (coordinates, not dose values; the producers sample the
ionmc dose at them without reading the table; the comparator asserts equality with the parsed depths)."""
N_DEPTHS = 46
Z_FIRST_CM, Z_LAST_CM = 1.533, 25.000
PARMS_FROZEN = {"p1": 21.321, "p2": 0.01372, "p5": 5.948, "p6": 31.55, "p7": -0.1462}
RHO_PLAN = -0.521  # the plan's rounding of the correlation coefficient
SIGMA_R_FRACTION = 0.008  # beam range spread, fraction of R0 (plan (g))
SPHERE_DIAMETER_MM = 4.0
SPHERE_SUBSAMPLES = 16
GRID_DR_MM = 0.5
GRID_DZ_MM = 0.5
GRID_NR = 210  # r to 105 mm
GRID_NZ = 520  # z to 260 mm
R0_VARIATION_CM = 0.05  # R0 +- 0.5 mm
VARIATION = 0.2
E0_BAND_RECORD = (176.0, 178.0)  # recording band (plan (g)): outside = range-model finding
E0_BAND_EXPECT = (176.5, 177.5)  # "177 +- 0.5 MeV" expectation
E0_BISECTION_TOL = 1.0e-6
E0_BRACKET_MEV = (150.0, 200.0)

TOL_R0 = 0.02
TOL_HALO_POINT = 0.15
HALO_FRAC = 0.90
HALO_MEDIAN = 0.10
SIGN_P = 0.01
SENS_TOL = 0.15
SENS_NSIGMA = 2.0
MASK_FRACTION = 1.0e-4
Z90 = 1.6448536269514722
BACKEND_Z_FINDING = 4.0
VARIANTS = ("sx-lo", "sx-hi", "st-lo", "st-hi", "rho-lo", "rho-hi")
R0_VARIANTS = ("r0-lo", "r0-hi")
SENS_NAMES = ("nominal", *VARIANTS, *R0_VARIANTS)

PLAN_QUOTE_STATUS = (
    "Amendment 14 (g)/row V6-sens: the r >= 1 cm halo gates only if, for every gating r >= 1 cm point and "
    "each of the six beam variants, |D_var/D_nom - 1| + 2 sigma_Delta < 0.15; otherwise V6 as a whole becomes "
    "report-only (Amendment 17 (c)7: the r = 0 verdict is still computed and reported)"
)


TRANSPORT_FINDING_NOTE = (
    "the comparison of the effective transport table range_g_cm2(E0) with the target at 1e-6 was an "
    "orchestrator-added consistency check, not a plan gate (Amendment 14 (g) fixes one table, the V3-003D "
    "exact-quadrature water table, and the 1e-6 relative tolerance of the bisection); it is reclassified "
    "as a recorded report-only finding after the value was observed (orchestrator decision 2026-10-11, "
    "to be recorded in decision 0041 by the orchestrator)"
)
V3_003D_TOLERANCES = (
    "tests/ionmc/test_range_quadrature.py D1b round trip: max |Rinv(R(E))/E - 1| <= 1e-6 and max |R(Rinv(r))/r - 1| "
    "<= 2e-6; U2 as frozen: both <= 1e-5 (validation/plans/v3-003d-acceptance.md); D1: the transport table equals "
    "the independent exact quadrature to 1e-12. NOTE: StoppingTable.range_at interpolates ln R log-log between "
    "the exact-quadrature nodes, so it differs from the transport table's exact in-bin quadrature by the "
    "interpolation residual recorded here"
)


class V6DataError(RuntimeError):
    """The e-print is missing or does not match its pin (the document records ``record``)."""

    def __init__(self, message: str, record: dict[str, Any]) -> None:
        super().__init__(message)
        self.record = record


# ---- data ---------------------------------------------------------------------------------------
def load_halo_tex(cache_dir: str | Path | None = None) -> tuple[str, dict[str, Any]]:
    """``(halo.tex text, record)`` from the hash-verified cache object of the registry entry.

    ``record`` = ``{dataset, sha256_tarball, sha256_halo_tex, tarball_ok, halo_tex_ok}``; any
    mismatch raises :class:`V6DataError` carrying the record (so the document can record it)."""
    from ionmc.data import cache
    from ionmc.data.registry import DATASETS

    record: dict[str, Any] = {"dataset": DATASET_ID, "pinned_sha256_tarball": TARBALL_SHA256,
                              "pinned_sha256_halo_tex": HALO_TEX_SHA256, "tarball_ok": False,
                              "halo_tex_ok": False, "sha256_tarball": None, "sha256_halo_tex": None}  # fmt: skip
    if DATASETS[DATASET_ID].sha256 != TARBALL_SHA256:
        raise V6DataError("registry pin differs from the plan's tarball hash", record)
    try:
        path = cache.verify(DATASET_ID, Path(cache_dir) if cache_dir else None)
        raw = path.read_bytes()
    except (FileNotFoundError, cache.IntegrityError) as exc:
        raise V6DataError(f"{type(exc).__name__}: {exc}", record) from exc
    record["sha256_tarball"] = hashlib.sha256(raw).hexdigest()
    record["tarball_ok"] = record["sha256_tarball"] == TARBALL_SHA256
    if not record["tarball_ok"]:
        raise V6DataError("tarball sha256 mismatch", record)
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tar:
        handle = tar.extractfile("halo.tex")
        if handle is None:
            raise V6DataError("halo.tex missing from the e-print", record)
        tex_bytes = handle.read()
    record["sha256_halo_tex"] = hashlib.sha256(tex_bytes).hexdigest()
    record["halo_tex_ok"] = record["sha256_halo_tex"] == HALO_TEX_SHA256
    if not record["halo_tex_ok"]:
        raise V6DataError("halo.tex sha256 mismatch (recorded)", record)
    return tex_bytes.decode("utf-8", errors="replace"), record


_ROW = re.compile(r"^\s*(-?\d+\.\d+)\s*&(.*?)\\\\")


def parse_dmlg(tex: str) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], int]:
    """``(z_cm[46], r_cm[10], log10_dose[46, 10] (NaN for 'M'), n_missing)`` of ``tbl:dmlg``.

    The rows are those of the ``rrrrrrrrrrr`` tabular that precedes ``\\label{tbl:dmlg}``; the
    radii are the ten header cells. Asserts 46 depths from 1.533 to 25.000 cm and 10 radii."""
    label = tex.index("\\label{tbl:dmlg}")
    begin = tex.rindex("\\begin{tabular}{rrrrrrrrrrr}", 0, label)
    end = tex.index("\\end{tabular}", begin)
    body = tex[begin:end]
    radii = [float(x) for x in re.findall(r"\\multicolumn\{1\}\{c\}\{([\d.]+)\}", body)]
    if radii[: len(RADII_CM)] != [float(r) for r in RADII_CM] or len(radii) < len(RADII_CM):
        raise ValueError(f"tbl:dmlg radii {radii} differ from {RADII_CM}")
    z: list[float] = []
    vals: list[list[float]] = []
    for line in body.splitlines():
        m = _ROW.match(line)
        if m is None:
            continue
        cells = [c.strip() for c in m.group(2).split("&")]
        if len(cells) != len(RADII_CM):
            raise ValueError(f"tbl:dmlg row {m.group(1)}: {len(cells)} values, expected {len(RADII_CM)}")
        z.append(float(m.group(1)))
        vals.append([math.nan if c == "M" else float(c) for c in cells])
    if len(z) != N_DEPTHS:
        raise ValueError(f"tbl:dmlg has {len(z)} depths, expected {N_DEPTHS}")
    if abs(z[0] - Z_FIRST_CM) > 1e-9 or abs(z[-1] - Z_LAST_CM) > 1e-9 or sorted(z) != z:
        raise ValueError("tbl:dmlg depth range differs from 1.533 .. 25.000 cm")
    if tuple(z) != DEPTHS_CM:
        raise ValueError("tbl:dmlg depths differ from DEPTHS_CM")
    arr = np.array(vals, dtype=np.float64)
    return (np.array(z), np.array(RADII_CM, dtype=np.float64), arr, int(np.isnan(arr).sum()))


def parse_parms(tex: str) -> dict[str, float]:
    """``{p1..p25: value}`` of ``tbl:parms``; asserts the five frozen values exactly."""
    label = tex.index("\\label{tbl:parms}")
    begin = tex.rindex("\\begin{tabular}{lrll}", 0, label)
    end = tex.index("\\end{tabular}", begin)
    out: dict[str, float] = {}
    for line in tex[begin:end].splitlines():
        m = re.match(r"^\$p_\{?(\d+)\}?\$&\s*\$?(-?[\d.]+)\$?&", line)
        if m:
            out[f"p{m.group(1)}"] = float(m.group(2))
    for k, v in PARMS_FROZEN.items():
        if out.get(k) != v:
            raise ValueError(f"tbl:parms {k} = {out.get(k)!r}, the plan freezes {v}")
    return out


# ---- beam ---------------------------------------------------------------------------------------
def beam_numbers(p: dict[str, float]) -> dict[str, Any]:
    """Derived beam numbers (plan (g)): f = -(1 + p7), rho = -sqrt(1 - f^2), A1 = rho sigma_x sigma_theta
    [cm rad], emittance pi sigma_x sigma_theta sqrt(1 - rho^2) [mm mrad]."""
    sx_mm, st_mrad, p7 = p["p5"], p["p6"], p["p7"]
    f = -(1.0 + p7)
    rho = -math.sqrt(1.0 - f * f)
    sx_cm, st_rad = sx_mm / 10.0, st_mrad * 1.0e-3
    a2, a0 = sx_cm**2, st_rad**2
    a1 = -math.sqrt(1.0 - f * f) * math.sqrt(a0 * a2)
    return {
        "sigma_x_mm": sx_mm, "sigma_theta_mrad": st_mrad, "p7": p7, "f": f,
        "rho_xtheta": rho, "rho_xtheta_plan_rounded": round(rho, 3), "A0_rad2": a0, "A2_cm2": a2,
        "A1_cm_rad": a1, "emittance_pi_mm_mrad": sx_mm * st_mrad * math.sqrt(1.0 - rho * rho),
        "sigma_range_fraction": SIGMA_R_FRACTION, "convention": "(tan thx, tan thy, 1) normalised",
    }  # fmt: skip


def solve_e0(range_fn: Callable[[float], float], target_g_cm2: float,
             bracket: tuple[float, float] = E0_BRACKET_MEV,
             tol: float = E0_BISECTION_TOL) -> tuple[float, int, list[float]]:
    """Bisection of the strictly increasing ``range_fn(E)`` [g/cm2] to ``range_fn(E) = target`` with
    relative tolerance ``tol`` on E; returns ``(E0, iterations, trace)``."""
    lo, hi = bracket
    if not range_fn(lo) < target_g_cm2 < range_fn(hi):
        raise ValueError("the target range is not bracketed")
    trace: list[float] = []
    e = 0.5 * (lo + hi)
    for it in range(1, 200):
        e = 0.5 * (lo + hi)
        trace.append(e)
        if range_fn(e) < target_g_cm2:
            lo = e
        else:
            hi = e
        if (hi - lo) <= tol * e:
            e = 0.5 * (lo + hi)
            return e, it, trace
    raise ValueError("bisection did not converge")


def e0_record(
    r0_cm: float,
    density_g_cm3: float,
    range_fn: Callable[[float], float],
    exact_inverse_fn: Callable[[float], float],
    stopping_mass_fn: Callable[[float], float],
    transport_range_fn: Callable[[float], float],
    table_identity: dict[str, Any],
    code_sha: str | None = None,
) -> dict[str, Any]:
    """The deterministic E0 record (plan (g)): bisection of the water range table to R0 * rho, the two
    cross-checks (exact inverse, transport table), sigma_E from sigma_R,beam = 0.008 R0 and the
    stopping power at E0 (dR/dE = 1/S of the same table). Nothing here reads the dose table."""
    target = r0_cm * density_g_cm3
    e0, it, trace = solve_e0(range_fn, target)
    inv = float(exact_inverse_fn(target))
    tr_range = float(transport_range_fn(e0))
    tr_rel = abs(tr_range - target) / target
    inv_rel = abs(inv - e0) / e0
    s_mass = float(stopping_mass_fn(e0))
    sigma_r = SIGMA_R_FRACTION * r0_cm
    sigma_e = sigma_r * density_g_cm3 * s_mass
    h = 1.0e-3 * e0  # finite-difference dR/dE of the interpolated range (recorded, not used)
    drde_fd = (range_fn(e0 + h) - range_fn(e0 - h)) / (2.0 * h)
    sigma_e_fd = sigma_r * density_g_cm3 / drde_fd
    rec: dict[str, Any] = {
        "row": "V6", "R0_cm": r0_cm, "density_g_cm3": density_g_cm3, "target_range_g_cm2": target,
        "E0_MeV": e0, "bisection": {"tolerance_relative": E0_BISECTION_TOL, "iterations": it,
                                    "bracket_MeV": list(E0_BRACKET_MEV), "trace_length": len(trace)},
        "range_at_E0_g_cm2": float(range_fn(e0)),
        "cross_check_exact_inverse": {"E_MeV": inv, "relative_difference": inv_rel,
                                      "ok": inv_rel <= 1e-6},
        "transport_table_finding": {
            "transport_table_range_g_cm2": tr_range, "target_range_g_cm2": target, "rel_diff": tr_rel,
            "energy_equivalent_rel": tr_rel * target * s_mass / e0,
            "within_1e-6": tr_rel <= 1e-6,
            "status": "finding (report-only): resampling residual of the transport table against the "
                      "V3-003D source table",
            "reclassified_after_value_observed": TRANSPORT_FINDING_NOTE,
            "v3_003d_tolerances": V3_003D_TOLERANCES,
            "within_v3_003d_u2": tr_rel <= 1e-5,
        },
        "sigma_R_beam_cm": sigma_r, "S_mass_at_E0_MeV_cm2_g": s_mass,
        "dR_dE_g_cm2_per_MeV": 1.0 / s_mass, "sigma_E_MeV": sigma_e,
        "sigma_E_finite_difference_MeV": sigma_e_fd,
        "inside_recording_band": E0_BAND_RECORD[0] <= e0 <= E0_BAND_RECORD[1],
        "inside_expectation_band": E0_BAND_EXPECT[0] <= e0 <= E0_BAND_EXPECT[1],
        "finding_outside_176_178": not (E0_BAND_RECORD[0] <= e0 <= E0_BAND_RECORD[1]),
        "table_identity": table_identity, "code_sha": code_sha,
        "no_detour_factor": True, "dose_table_read": False,
    }  # fmt: skip
    rec["pass"] = bool(rec["cross_check_exact_inverse"]["ok"])  # the bisection itself converged to 1e-6
    return rec


def beam_parameters(p: dict[str, float], e0_for_r0: Callable[[float], tuple[float, float]]) -> dict[str, dict[str, float]]:
    """Source parameters of the nominal run and the eight V6-sens variants.

    ``e0_for_r0(R0_cm) -> (E0, sigma_E)`` is the E0 procedure. The six beam variants scale sigma_x,
    sigma_theta and the correlation coefficient by 0.8 and 1.2 at fixed E0 (Amendment 17 (c)5: rho_xtheta
    itself is varied, to -0.417 and -0.625, with sigma_x, sigma_theta fixed); the two R0 variants
    recompute E0 and sigma_E (0.8 % of the varied R0) at R0 +- 0.5 mm."""
    b = beam_numbers(p)
    e0, se = e0_for_r0(p["p1"])
    nominal = {"E0_MeV": e0, "sigma_E_MeV": se, "sigma_x_mm": b["sigma_x_mm"],
               "sigma_theta_mrad": b["sigma_theta_mrad"], "rho_xtheta": b["rho_xtheta"]}  # fmt: skip
    out = {"nominal": dict(nominal)}
    lo, hi = 1.0 - VARIATION, 1.0 + VARIATION
    for name, key, s in (("sx-lo", "sigma_x_mm", lo), ("sx-hi", "sigma_x_mm", hi),
                         ("st-lo", "sigma_theta_mrad", lo), ("st-hi", "sigma_theta_mrad", hi),
                         ("rho-lo", "rho_xtheta", lo), ("rho-hi", "rho_xtheta", hi)):  # fmt: skip
        out[name] = {**nominal, key: nominal[key] * s}
    for name, dr in (("r0-lo", -R0_VARIATION_CM), ("r0-hi", R0_VARIATION_CM)):
        e, s_e = e0_for_r0(p["p1"] + dr)
        out[name] = {**nominal, "E0_MeV": e, "sigma_E_MeV": s_e}
    return out


# ---- sphere averaging ---------------------------------------------------------------------------
def sphere_weights(r_mm: float, z_mm: float, dr: float = GRID_DR_MM, dz: float = GRID_DZ_MM,
                   n_r: int = GRID_NR, n_z: int = GRID_NZ, diameter_mm: float = SPHERE_DIAMETER_MM,
                   n_sub: int = SPHERE_SUBSAMPLES) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.float64]]:
    """Overlap volumes [mm3] of the sphere of ``diameter_mm`` centred at ``(r_mm, 0, z_mm)`` with every
    annular cell it touches: ``(ir, iz, w)``. Each cell is sampled by ``n_sub`` midpoints in each of
    (rho, phi, z) with the volume element rho drho dphi dz; the phi range is the angular half-width
    ``asin(R / r)`` of the sphere's footprint (the full circle for r <= R)."""
    rad = 0.5 * diameter_mm
    ir0 = max(int(math.floor((r_mm - rad) / dr)), 0)
    ir1 = min(int(math.floor((r_mm + rad) / dr)), n_r - 1)
    iz0 = max(int(math.floor((z_mm - rad) / dz)), 0)
    iz1 = min(int(math.floor((z_mm + rad) / dz)), n_z - 1)
    phi_half = math.pi if r_mm <= rad else math.asin(rad / r_mm)
    mid = (np.arange(n_sub) + 0.5) / n_sub
    phi = (2.0 * mid - 1.0) * phi_half
    dphi = 2.0 * phi_half / n_sub
    cphi = np.cos(phi)
    ir_l: list[int] = []
    iz_l: list[int] = []
    w_l: list[float] = []
    for ir in range(ir0, ir1 + 1):
        rho = (ir + mid) * dr
        for iz in range(iz0, iz1 + 1):
            zz = (iz + mid) * dz - z_mm
            d2 = (rho[:, None, None] ** 2 + r_mm**2 - 2.0 * rho[:, None, None] * r_mm * cphi[None, :, None]
                  + zz[None, None, :] ** 2)  # fmt: skip
            inside = d2 <= rad * rad
            vol = float(np.sum(inside * rho[:, None, None])) * (dr / n_sub) * dphi * (dz / n_sub)
            if vol > 0.0:
                ir_l.append(ir)
                iz_l.append(iz)
                w_l.append(vol)
    return np.array(ir_l, dtype=np.int64), np.array(iz_l, dtype=np.int64), np.array(w_l)


def all_sphere_weights(z_cm: Sequence[float], r_cm: Sequence[float]) -> dict[tuple[int, int], tuple[Any, Any, Any]]:
    """Weights of every data point ``(j, i)`` (depth index, radius index)."""
    return {(j, i): sphere_weights(10.0 * float(r_cm[i]), 10.0 * float(z_cm[j]))
            for j in range(len(z_cm)) for i in range(len(r_cm))}  # fmt: skip


def weights_digest(weights: dict[tuple[int, int], tuple[Any, Any, Any]]) -> str:
    h = hashlib.sha256()
    for key in sorted(weights):
        ir, iz, w = weights[key]
        h.update(repr(key).encode())
        h.update(np.ascontiguousarray(ir).tobytes())
        h.update(np.ascontiguousarray(iz).tobytes())
        h.update(np.ascontiguousarray(w).tobytes())
    return h.hexdigest()


def sphere_average_batches(
    dose_rz: NDArray[np.float64], weights: dict[tuple[int, int], tuple[Any, Any, Any]], n_z: int, n_r: int
) -> NDArray[np.float64]:
    """Per-batch sphere-averaged dose ``[B, n_z, n_r]`` from the annular dose ``dose_rz[B, n_rcells,
    n_zcells]`` (MeV/g per primary): ``sum_c w_c D_c / sum_c w_c`` for every data point."""
    b = dose_rz.shape[0]
    out = np.full((b, n_z, n_r), np.nan)
    for (j, i), (ir, iz, w) in weights.items():
        out[:, j, i] = (dose_rz[:, ir, iz] * w[None, :]).sum(axis=1) / w.sum()
    return out


# ---- statistics ---------------------------------------------------------------------------------
def binom_two_sided_p(n: int, k: int) -> float:
    """Two-sided exact binomial (p = 1/2) p-value for ``k`` successes in ``n`` trials."""
    if n == 0:
        return 1.0
    m = min(k, n - k)
    tail = sum(math.comb(n, i) for i in range(m + 1)) / 2.0**n
    return float(min(1.0, 2.0 * tail))


def batch_mean_se(batches: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Mean and standard error of the mean over axis 0 (equal-size batches)."""
    nb = batches.shape[0]
    return batches.mean(axis=0), batches.std(axis=0, ddof=1) / math.sqrt(nb)


def ratio_paired(var: NDArray[np.float64], nom: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Ratio estimator ``mean(var)/mean(nom)`` and its paired batch standard error (linearisation,
    ``se^2 = sum (x_b - R y_b)^2 / (B (B - 1) ybar^2)``) per element; batches on axis 0 pair up."""
    nb = var.shape[0]
    ym = nom.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = var.mean(axis=0) / ym
        se = np.sqrt(((var - r[None] * nom) ** 2).sum(axis=0) / (nb * (nb - 1))) / np.abs(ym)
    return r, se


# ---- comparator ---------------------------------------------------------------------------------
def measured_dose(log10_dose: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.power(10.0, log10_dose)


def _range_shift_mm(z_cm: NDArray[np.float64], d_meas: NDArray[np.float64], j: int, residual: float) -> float | None:
    """Range-shift equivalent ``residual / (d ln D / dz)`` [mm] at depth index ``j`` of the r = 0 column
    from the measured points (central difference over the neighbouring non-missing depths)."""
    lo = j - 1 if j > 0 and np.isfinite(d_meas[j - 1]) else j
    hi = j + 1 if j + 1 < len(z_cm) and np.isfinite(d_meas[j + 1]) else j
    if lo == hi:
        return None
    g = (math.log(d_meas[hi]) - math.log(d_meas[lo])) / (10.0 * (z_cm[hi] - z_cm[lo]))  # per mm
    if abs(g) < 1e-6:
        return None
    return float(residual / g)


def evaluate(
    z_cm: NDArray[np.float64], r_cm: NDArray[np.float64], d_meas: NDArray[np.float64],
    batches: NDArray[np.float64], sens: dict[str, NDArray[np.float64]] | None = None,
    sens_nominal: NDArray[np.float64] | None = None,
) -> dict[str, Any]:
    """The V6 verdict of one backend (Amendment 14 (g), Amendment 17 (c)).

    ``d_meas[46, 10]`` measured dose (NaN missing), ``batches[B, 46, 10]`` the sphere-averaged dose
    per primary of the V6 sample, ``sens`` the six beam-variant batches ``name -> [B, 46, 10]`` and
    ``sens_nominal`` the CRN nominal (``None``: V6-sens not available, the status is "pending")."""
    mean, se = batch_mean_se(batches)
    r0_col = int(np.where(r_cm == 0)[0][0])
    present = np.isfinite(d_meas)
    r0_pts: list[dict[str, Any]] = []
    for j in range(len(z_cm)):
        if not present[j, r0_col]:
            continue
        dm, di, s = float(d_meas[j, r0_col]), float(mean[j, r0_col]), float(se[j, r0_col])
        res = di / dm - 1.0
        half = Z90 * s / dm
        r0_pts.append({
            "z_cm": float(z_cm[j]), "D_meas": dm, "D_ionmc": di, "rel_se": s / di if di > 0 else None,
            "residual": res, "interval90": [res - half, res + half],
            "range_shift_mm": _range_shift_mm(z_cm, d_meas[:, r0_col], j, res),
            "pass": bool(abs(res) <= TOL_R0),
        })  # fmt: skip
    r0_pass = bool(r0_pts) and all(p["pass"] for p in r0_pts)
    r0_max = float(np.nanmax(d_meas[:, r0_col]))
    halo_pts: list[dict[str, Any]] = []
    for j in range(len(z_cm)):
        for i in range(len(r_cm)):
            if r_cm[i] < 1 or not present[j, i]:
                continue
            dm, di, s = float(d_meas[j, i]), float(mean[j, i]), float(se[j, i])
            res = di / dm - 1.0
            halo_pts.append({
                "z_cm": float(z_cm[j]), "r_cm": float(r_cm[i]), "D_meas": dm, "D_ionmc": di,
                "rel_se": s / di if di > 0 else None, "residual": res,
                "measured_over_r0_max": dm / r0_max, "in_sizing_set": bool(dm >= MASK_FRACTION * r0_max),
                "pass": bool(abs(res) <= TOL_HALO_POINT), "range_shift_mm": None,
            })  # fmt: skip
    absres = np.array([abs(p["residual"]) for p in halo_pts])
    frac = float(np.mean(absres <= TOL_HALO_POINT)) if len(absres) else 0.0
    med = float(np.median(absres)) if len(absres) else math.inf
    # no r-trend: per depth with points in both groups, sign of median(4-10 cm) - median(1-3 cm)
    n_pos = n_neg = n_tie = 0
    for j in range(len(z_cm)):
        a = [p["residual"] for p in halo_pts if p["z_cm"] == float(z_cm[j]) and 1 <= p["r_cm"] <= 3]
        c = [p["residual"] for p in halo_pts if p["z_cm"] == float(z_cm[j]) and p["r_cm"] >= 4]
        if not a or not c:
            continue
        diff = float(np.median(c) - np.median(a))
        if diff > 0:
            n_pos += 1
        elif diff < 0:
            n_neg += 1
        else:
            n_tie += 1
    n_used = n_pos + n_neg
    p_sign = binom_two_sided_p(n_used, n_pos)
    sizing = [p for p in halo_pts if p["in_sizing_set"]]
    halo_pass = bool(len(halo_pts) > 0 and frac >= HALO_FRAC and med <= HALO_MEDIAN and p_sign > SIGN_P)
    sens_doc = _sens_status(halo_pts, z_cm, r_cm, sens, sens_nominal)
    gating = sens_doc["condition_holds"]
    status = "gating" if gating else ("pending" if gating is None else "report-only")
    gate_pass: bool | None
    if gating:
        gate_pass = bool(r0_pass and halo_pass)
    elif gating is None:
        gate_pass = None
    else:
        gate_pass = None  # report-only: recorded, not gating
    return {
        "r0": {"points": r0_pts, "pass": r0_pass, "n_points": len(r0_pts), "tolerance": TOL_R0,
               "rule": "|D_ionmc/D_meas - 1| <= 0.02 at every non-missing r = 0 depth"},
        "halo": {"points": halo_pts, "n_points": len(halo_pts), "frac_pass": frac, "median_abs": med,
                 "sign_test": {"n_depths": n_used + n_tie, "n_pos": n_pos, "n_neg": n_neg, "n_ties_excluded": n_tie,
                               "p": p_sign, "pass": bool(p_sign > SIGN_P)},
                 "sizing_set": {"n_points": len(sizing),
                                "frac_pass": float(np.mean([p["pass"] for p in sizing])) if sizing else None},
                 "pass": halo_pass, "gating": gating,
                 "rule": "per point 0.15; >= 90 % of points pass; median |residual| <= 0.10; sign-test p > 0.01"},
        "v6_sens": sens_doc, "status": status, "gate_pass": gate_pass,
        "r0_max_measured": r0_max, "plan_rule_status": PLAN_QUOTE_STATUS,
    }  # fmt: skip


def _sens_status(
    halo_pts: list[dict[str, Any]], z_cm: NDArray[np.float64], r_cm: NDArray[np.float64],
    sens: dict[str, NDArray[np.float64]] | None, nominal: NDArray[np.float64] | None,
) -> dict[str, Any]:
    if sens is None or nominal is None:
        return {"condition_holds": None, "variants": {}, "note": "V6-sens samples not supplied"}
    zi = {float(z): j for j, z in enumerate(z_cm)}
    ri = {float(r): i for i, r in enumerate(r_cm)}
    variants: dict[str, Any] = {}
    holds = True
    for name in VARIANTS:
        if name not in sens:
            return {"condition_holds": None, "variants": variants, "note": f"variant {name} missing"}
        ratio, se = ratio_paired(sens[name], nominal)
        worst = -math.inf
        pts: list[dict[str, Any]] = []
        for p in halo_pts:
            j, i = zi[p["z_cm"]], ri[p["r_cm"]]
            m = abs(float(ratio[j, i]) - 1.0) + SENS_NSIGMA * float(se[j, i])
            if not math.isfinite(m):
                m = math.inf
            worst = max(worst, m)
            pts.append({"z_cm": p["z_cm"], "r_cm": p["r_cm"], "ratio": float(ratio[j, i]),
                        "sigma_delta": float(se[j, i]), "metric": m, "ok": bool(m < SENS_TOL)})  # fmt: skip
        ok = bool(worst < SENS_TOL)
        holds = holds and ok
        variants[name] = {"max_metric": worst, "ok": ok, "n_points": len(pts), "per_point": pts}
    return {"condition_holds": bool(holds), "variants": variants, "tolerance": SENS_TOL, "n_sigma": SENS_NSIGMA,
            "rule": PLAN_QUOTE_STATUS}  # fmt: skip


def r0_variants_report(
    z_cm: NDArray[np.float64], r_cm: NDArray[np.float64], sens: dict[str, NDArray[np.float64]] | None,
    nominal: NDArray[np.float64] | None, halo_points: list[dict[str, Any]],
) -> dict[str, Any]:
    """The R0 +- 0.5 mm variants (report-only): the same metric as the beam variants."""
    if sens is None or nominal is None:
        return {}
    zi = {float(z): j for j, z in enumerate(z_cm)}
    ri = {float(r): i for i, r in enumerate(r_cm)}
    out: dict[str, Any] = {}
    for name in R0_VARIANTS:
        if name not in sens:
            continue
        ratio, se = ratio_paired(sens[name], nominal)
        worst = max(abs(float(ratio[zi[p["z_cm"]], ri[p["r_cm"]]]) - 1.0)
                    + SENS_NSIGMA * float(se[zi[p["z_cm"]], ri[p["r_cm"]]]) for p in halo_points)
        out[name] = {"max_metric": worst, "report_only": True}
    return out


def combine_backends(cpu: dict[str, Any], cuda: dict[str, Any]) -> dict[str, Any]:
    """V6 over the two backends (Amendment 17 (c)8): V6 fails if the lv5c and hr5c verdicts differ on any
    gate component (r = 0, halo, overall); the per-point paired z between the backends is reported with
    any |z| > 4 listed as a finding. V6 passes only if both backends pass (and the status is gating)."""
    comps = {"r0": (cpu["r0"]["pass"], cuda["r0"]["pass"]), "halo": (cpu["halo"]["pass"], cuda["halo"]["pass"]),
             "overall": (cpu["gate_pass"], cuda["gate_pass"])}  # fmt: skip
    disagree = sorted(k for k, (a, b) in comps.items() if a != b)
    zs: list[dict[str, Any]] = []
    findings: list[dict[str, Any]] = []
    for grp in ("r0", "halo"):
        for pc, pg in zip(cpu[grp]["points"], cuda[grp]["points"], strict=True):
            sc, sg = (pc["rel_se"] or 0.0) * pc["D_ionmc"], (pg["rel_se"] or 0.0) * pg["D_ionmc"]
            den = math.hypot(sc, sg)
            z = (pc["D_ionmc"] - pg["D_ionmc"]) / den if den > 0 else math.nan
            rec = {"group": grp, "z_cm": pc["z_cm"], "r_cm": pc.get("r_cm", 0.0), "z": z}
            zs.append(rec)
            if math.isfinite(z) and abs(z) > BACKEND_Z_FINDING:
                findings.append(rec)
    status = cpu["status"]
    both = bool(cpu["gate_pass"]) and bool(cuda["gate_pass"])
    return {
        "disagreement_components": disagree, "paired_z": zs, "findings_abs_z_gt_4": findings,
        "status": status, "pass": (None if status != "gating" else bool(both and not disagree)),
        "r0_pass_both": bool(cpu["r0"]["pass"] and cuda["r0"]["pass"]),
        "rule": "V6 fails if the verdicts of lv5c and hr5c differ on any gate component",
    }  # fmt: skip
