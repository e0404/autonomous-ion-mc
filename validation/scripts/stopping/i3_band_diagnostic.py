"""Row I3 diagnostic: band-wise ionmc Bethe vs libdEdx ICRU73_NEW water for Z = 1-10 (no gate).

Usage::

    uv run python validation/scripts/stopping/i3_band_diagnostic.py --archive PATH \
        --output-dir DIR --code-sha "$(git rev-parse HEAD)" \
        [--expect-i3 validation/results/stopping/v3-008-i3/icru73_comparison.json] [--repo DIR]

Purpose (amendment I3-FU-1): attribute the row I3 FAIL (``compare_icru73.py``) by looking at
how the ionmc-vs-table stopping deviation is distributed over energy bands and over the ion
charge Z. This is a diagnostic only: it has no gate, and row I3 is not re-graded by it.

Reference and handling are those of ``compare_icru73.py``, which this script loads from the same
directory and whose fail-closed loader it reuses unchanged (archive sha256, the two member
sha256, ``tarfile`` filter ``"data"``, strict block parsing, 53-point grid, finite positive
values). The archive is untrusted data and nothing from it is executed.

Content of ``icru73_band_diagnostic.json``:

* ``ions``: for every ion Z = 1-10 in liquid water (target 276) and every band
  [2, 5), [5, 10), [10, 20), [20, 50), [50, 1000] MeV/u, only ``n_points``, ``band_per_u_mev`` and
  ``mean_relative_deviation`` = mean over the band's table grid points of S_model/S_table - 1,
  with the default :class:`BetheOptions` at ``water(78.0)``, ROUNDED to 0.1 percentage point
  (``round(d, 3)`` in relative units; see "Rounding" below). No band max/min/RMS is written (with
  2-4 points per band they would, with the mean, expose per-energy deviations); a band with fewer
  than two grid points raises ``ValueError`` (its mean would be a per-energy value). Carbon and
  oxygen additionally carry ``max_abs_relative_deviation_2_to_10`` (2 <= E < 10 MeV/u) and
  ``max_abs_relative_deviation_ge_10`` (E >= 10 MeV/u to the grid end) at full precision, the
  two band maxima the committed row I3 file already publishes; they are checked against that
  file (``--expect-i3``, see "Provenance"; rtol 1e-9; ``RuntimeError`` and no output on mismatch).
* ``hypothesis_check`` (amendment I3-FU-1 rule 1a): the recorded C/O hypothesis (``HYPOTHESIS``,
  ``HYPOTHESIS_HIGH_ABS``) against the rounded band means. CONFIRMED iff, for both ions, every
  rounded band mean below 50 MeV/u equals the quoted value (equality at 0.1 pp, i.e.
  |delta| <= 0.05 pp before rounding) and |rounded mean of [50, 1000]| <= the quoted bound;
  otherwise WITHDRAWN.
* ``z_scan``: for the bands below 50 MeV/u, single-term least squares without intercept of the
  rounded band means d_z over Z = 3-10 (the ICRU 73 ion set) against four shapes: ``charge``
  f = band mean of 1 - q^2/z^2 (q the model's effective charge, ``pierce_blann_charge``),
  ``barkas_like`` f = z, ``bloch_like`` f = z^2 and the z-independent control ``constant``
  f = 1. Coefficient, residual RMS and the uncentred
  r^2 = 1 - sum(res^2)/sum(d^2) (the no-intercept convention), a ranking by residual RMS, the
  ratios d_z/d_6 and the rule-1b ``verdict``: with RMS_1 <= RMS_2 the two smallest residual RMS
  and r^2_1 that of the best shape, ``attributed_to`` = the best shape iff RMS_1 <= 0.5 RMS_2,
  r^2_1 >= 0.90 and the rounded d_z over Z = 3-10 have one common sign (no zero, no sign change);
  otherwise ``attributed_to`` = null with the failed conditions in ``reason`` (``rms_ratio``,
  ``r_squared``, ``sign_change``). Z = 1 and 2 are reported separately (``informative_z1_z2``):
  their table lineage is OPEN.
* ``i3b`` (informative, carbon12 and oxygen16): model variants with one correction switched off
  or the Lindhard-Sorensen modes. Only the ``default`` variant carries the table-relative band
  mean (a duplicate of the ``ions`` entry); every variant carries model-only quantities only
  (CSDA-range shift vs the default model at the row I3 standard energies, band mean of the
  stopping ratio r_v to the default model). Reason: eight variant-vs-table band means at full
  precision are eight linear equations in the inverse table values of the band, so the
  per-energy table values of a band with <= 8 grid points could be solved for from the output
  and the public model. The omitted variant-vs-table mean is approximately
  (1 + d_default) r_v - 1.

Ions other than the ionmc projectiles (Li-7, Be-9, B-11, N-14, F-19, Ne-20) use the nuclear mass
approximated as A u - Z m_e (``mass_note``; the recorded model-only sensitivity shows the
approximation is immaterial here).

Rounding (``rounding_note``): every table-relative band mean in the output (all ions, all bands,
``informative_z1_z2``, the ``i3b`` default duplicate) is rounded to 0.1 percentage point before
it is written and before any derived quantity is computed; the z-scan fits, ``shape_ratio_to_z6``,
``prediction_from_z3_10_fit``, ``hypothesis_check`` and the z-scan verdicts are computed from the
rounded means only, so nothing in the output depends on an unrounded mean. Reason: on the real
53-point grid the [10, 20) band holds two points, and the row I3 file already publishes the
>= 10 MeV/u maximum, which for a monotone deficit sits at the 10 MeV/u point; a full-precision
[10, 20) mean would then give d(15) = 2 mean - (+-max) and, with the public model, the table value
at 15 MeV/u, and the z-scan coefficients and d_z/d_6 ratios would expose the C/O means again. At
0.1 pp the same solve fixes d(15) only to +-0.1 pp, which is the precision already present in Git
(plan results note of 2026-11-15); the fits are functions of the rounded means and add nothing.
Cost: the quantisation limits the z-scan in the [20, 50) band, where the deviations are a few
0.1 pp.

Provenance (fail closed; ``RuntimeError`` and nothing written on any failure; the checkout,
environment and ``--expect-i3`` checks run before the computation, the source sealing after it
and before the output is written):

* ``--code-sha`` is required and must be 40 lowercase hex digits; in ``--repo`` (default: the
  repository containing this script) ``git rev-parse HEAD`` must equal it and
  ``git status --porcelain --untracked-files=no`` must be empty (no modified tracked file).
* Environment (``provenance.environment``, checked against the ``uv.lock`` committed at
  ``code_sha``): the running ``numpy.__version__`` must equal the single numpy ``[[package]]``
  version of the lock, and the running Python must satisfy the lock's ``requires-python``
  (recorded as is when the lock has none; an unsupported specifier fails closed). The
  environment seal (``env_seal.py``, Codex REVIEW-c4f3fdd1 finding 3; ``--environment-seal``,
  default ``validation/results/stopping/v3-008-i3fu1/environment-seal.json``, read as committed
  at ``code_sha``): the installed numpy distribution is verified file by file against its wheel
  RECORD, and the numpy version, the RECORD sha256, the sha256 of numpy's
  ``_multiarray_umath`` extension, the sha256 of the resolved Python executable,
  ``sys.version`` and ``platform.platform()`` must equal the seal (recorded as
  ``provenance.environment.seal``). ``--write-environment-seal`` writes the seal of the running
  environment to ``--repo/--environment-seal`` and exits (run once by the orchestrator and
  committed with the review pin). Runtime controls (Codex REVIEW-183057ce finding 3), after
  the computation: every loaded numpy module resolves under the verified distribution root and
  is hashed in its RECORD, ``sys.dont_write_bytecode`` is true (PYTHONDONTWRITEBYTECODE=1) and
  no sealed or numpy module has its ``__cached__`` bytecode file on disk (the ``__pycache__``
  directories are removed before the run); recorded as ``provenance.environment.runtime``.
* Sealed sources after the computation (the import closure is then complete, including the
  lazily loaded ``compare_icru73`` and anything ``BetheStoppingSource`` imports): every module
  in ``sys.modules`` of the ``ionmc`` package with a ``__file__`` (path relative to the root
  that contains ``src/ionmc``), this script and ``compare_icru73.py`` (loaded by path). For each,
  the file actually loaded by this process, the file in ``--repo`` and
  ``git show <code_sha>:<path>`` must have identical bytes; a module loaded from outside a
  ``src/ionmc`` tree of a repository (e.g. a site-packages install) has no committed path and
  fails. ``uv.lock`` and ``pyproject.toml`` are sealed too (disk == committed bytes). The
  ``{path: sha256}`` map is ``provenance.sealed_sources``; ``provenance.ionmc_loaded_from_repo``
  says whether the loaded ionmc files lie under ``--repo`` itself (byte identity, not location,
  is the requirement). ``env_seal.py`` (loaded by path) and the environment seal are sealed
  too. Declared limitation: shared libraries outside the numpy distribution and the
  interpreter binary (libc, system libraries) are not sealed beyond ``platform.platform()``.
* ``--expect-i3`` is a repo-relative path (default ``DEFAULT_EXPECT_I3``) read as committed
  bytes at ``code_sha`` (``git show``), not from the working copy; the working copy must be
  byte-identical. Recorded as ``provenance.expect_i3`` = ``{path, sha256, committed_at, matched}``.
* ``provenance.code_sha`` and ``provenance.ionmc_version`` (``ionmc.__version__``).

Git is run as a subprocess without a shell and with inherited ``GIT_*`` variables removed.

The JSON output contains aggregates and model-only quantities only; no libdEdx/ICRU table
values, per-energy deviations or ratios are written (fields ``redistribution_note`` and
``rounding_note``).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
import tomllib
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.materials import Material, water
from ionmc.physics.projectiles import (
    ALPHA,
    AMU_MEV,
    CARBON12,
    ELECTRON_MASS_MEV,
    OXYGEN16,
    PROTON,
    Projectile,
)
from ionmc.physics.stopping import (
    BetheOptions,
    BetheStoppingSource,
    bethe_mass_stopping,
    kinematics,
    pierce_blann_charge,
)


def _load_env_seal() -> ModuleType:
    """Load the shared environment seal (``env_seal.py``) from this directory."""
    name = "env_seal"
    if name in sys.modules:
        return sys.modules[name]
    path = Path(__file__).resolve().parent / "env_seal.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _load_i3() -> ModuleType:
    """Load the row I3 script (``compare_icru73.py``) from this directory."""
    path = Path(__file__).resolve().with_name("compare_icru73.py")
    spec = importlib.util.spec_from_file_location("_i3_compare_icru73", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["_i3_compare_icru73"] = module
    spec.loader.exec_module(module)
    return module


I3 = _load_i3()
ENV = _load_env_seal()
ARCHIVE_SHA256: str = I3.ARCHIVE_SHA256
MEMBER_SHA256: dict[str, str] = I3.MEMBER_SHA256
TARGET_WATER: int = I3.TARGET_WATER
I_EV_MODEL: float = I3.I_EV_MODEL
REDISTRIBUTION_NOTE: str = I3.REDISTRIBUTION_NOTE
I_VALUE_NOTE: str = I3.I_VALUE_NOTE

E_MIN = 2.0  # MeV/u
BANDS: tuple[tuple[float, float], ...] = (
    (2.0, 5.0),
    (5.0, 10.0),
    (10.0, 20.0),
    (20.0, 50.0),
    (50.0, 1000.0),  # closed at the grid end
)
Z_SCAN_BANDS = BANDS[:4]
Z_SCAN_IONS = tuple(range(3, 11))
Z_REFERENCE = 6
MIN_BAND_POINTS = 2
MASS_SENSITIVITY_FACTOR = 1.001


def _approx_nucleus(name: str, symbol: str, z: int, a: int) -> Projectile:
    """Nucleus with mass A u - Z m_e (atomic-mass approximation; see ``MASS_NOTE``)."""
    return Projectile(name, symbol, z, a, a * AMU_MEV - z * ELECTRON_MASS_MEV)


PROJECTILES: dict[int, Projectile] = {
    1: PROTON,
    2: ALPHA,
    3: _approx_nucleus("lithium7", "Li7", 3, 7),
    4: _approx_nucleus("beryllium9", "Be9", 4, 9),
    5: _approx_nucleus("boron11", "B11", 5, 11),
    6: CARBON12,
    7: _approx_nucleus("nitrogen14", "N14", 7, 14),
    8: OXYGEN16,
    9: _approx_nucleus("fluorine19", "F19", 9, 19),
    10: _approx_nucleus("neon20", "Ne20", 10, 20),
}
APPROXIMATE_MASS_Z = (3, 4, 5, 7, 9, 10)
MASS_SENSITIVITY_Z = 7
I3_MAX_Z = (6, 8)  # the ions whose band maxima the row I3 file already publishes
E_I3_SPLIT = 10.0  # MeV/u, the row I3 band split
EXPECT_I3_RTOL = 1e-9
ROUND_DIGITS = 3  # table-relative band means: 0.1 percentage point (relative units)

# Rule 1a (amendment I3-FU-1): the recorded hypothesis (plan results note of 2026-11-15), band
# means below 50 MeV/u at 0.1 pp and the bound on |mean| of [50, 1000] MeV/u (relative units).
HYPOTHESIS: dict[str, dict[str, float]] = {
    "carbon12": {"2_to_5": -0.086, "5_to_10": -0.043, "10_to_20": -0.018, "20_to_50": -0.003},
    "oxygen16": {"2_to_5": -0.122, "5_to_10": -0.075, "10_to_20": -0.039, "20_to_50": -0.008},
}
HYPOTHESIS_HIGH_ABS: dict[str, float] = {"carbon12": 0.0003, "oxygen16": 0.0006}
# Rule 1b (amendment I3-FU-1): attribution thresholds of the z-scan.
ATTRIBUTION_RMS_RATIO = 0.5
ATTRIBUTION_R_SQUARED = 0.90

I3B_IONS: dict[str, tuple[Projectile, tuple[float, ...]]] = {
    "carbon12": (CARBON12, (80.0, 290.0, 430.0)),
    "oxygen16": (OXYGEN16, (100.0, 350.0, 430.0)),
}
VARIANTS: dict[str, BetheOptions] = {
    "default": BetheOptions(),
    "shell_off": BetheOptions(shell=False),
    "barkas_off": BetheOptions(barkas=False),
    "bloch_off": BetheOptions(bloch=False),
    "mott_off": BetheOptions(mott=False),
    "effective_charge_off": BetheOptions(effective_charge=False),
    "ls_eff": BetheOptions(ls="eff"),  # LS replaces Bloch + Mott (both on)
    "ls_bare_mixed": BetheOptions(ls="bare_mixed"),
}

# Provenance (Codex REVIEW-22f384ed finding 6): exact clean commit, sealed sources, committed I3.
CODE_SHA_RE = re.compile(r"[0-9a-f]{40}")
DEFAULT_REPO = Path(__file__).resolve().parents[3]
DEFAULT_EXPECT_I3 = "validation/results/stopping/v3-008-i3/icru73_comparison.json"
SCRIPT_REL = "validation/scripts/stopping/i3_band_diagnostic.py"
I3_SCRIPT_REL = "validation/scripts/stopping/compare_icru73.py"
ENV_SCRIPT_REL = "validation/scripts/stopping/env_seal.py"
LOCK_FILES = ("uv.lock", "pyproject.toml")  # sealed (disk == committed bytes)
ENVIRONMENT_NOTE = (
    "numpy version and Python version checked against the uv.lock committed at code_sha; the "
    "numpy distribution (every RECORD-listed file), its _multiarray_umath extension and the "
    "Python executable are sealed by content against the committed environment seal; shared "
    "libraries outside them (libc, system libraries) are NOT sealed (declared limitation)"
)
_SPECIFIER_RE = re.compile(r"\s*(>=|<=|==|!=|>|<)\s*(\d+(?:\.\d+)*)(\.\*)?\s*")

STATUS_NOTE = (
    "diagnostic of the row I3 FAIL under amendment I3-FU-1; no gate; row I3 is not re-graded"
)
ROUNDING_NOTE = (
    "every table-relative band mean (ions, informative_z1_z2, i3b default) is rounded to 0.1 "
    "percentage point (round(d, 3) in relative units) before it is written and before any "
    "derived quantity is computed; z-scan fits, shape_ratio_to_z6, prediction_from_z3_10_fit, "
    "hypothesis_check and the z-scan verdicts use the rounded means only, so nothing here "
    "depends on an unrounded mean. 0.1 pp is the precision already present in Git (plan results "
    "note of 2026-11-15). Reason: the [10, 20) band holds two grid points and row I3 publishes "
    "the >= 10 MeV/u maximum (at 10 MeV/u for a monotone deficit), so a full-precision mean would "
    "give the deviation at the other point; at 0.1 pp it is fixed only to +-0.1 pp. The "
    "quantisation limits the z-scan in the [20, 50) band (deviations of a few 0.1 pp). The C/O "
    "band maxima (2-10, >= 10 MeV/u) are at full precision: row I3 already publishes them"
)
LINEAGE_NOTE = (
    "table lineage for Z = 1, 2 OPEN: excluded from the z-scan fits (Z = 3-10); band means and "
    "the predictions of the Z = 3-10 fits are informative only"
)


def _band_key(band: tuple[float, float]) -> str:
    return f"{band[0]:g}_to_{band[1]:g}_MeV_u"


def _band_mask(e_u: NDArray[np.float64], band: tuple[float, float]) -> NDArray[np.bool_]:
    lo, hi = band
    upper = e_u <= hi if band == BANDS[-1] else e_u < hi
    mask: NDArray[np.bool_] = (e_u >= lo) & upper
    n = int(np.count_nonzero(mask))
    if n < MIN_BAND_POINTS:
        raise ValueError(
            f"band {band} has {n} grid points (< {MIN_BAND_POINTS}): its mean would be a "
            "per-energy value"
        )
    return mask


def round_deviation(value: float) -> float:
    """A table-relative band mean rounded to 0.1 percentage point (``-0.0`` written as ``0.0``)."""
    return round(float(value), ROUND_DIGITS) + 0.0


def model_stopping(
    e_u: NDArray[np.float64], proj: Projectile, material: Material, o: BetheOptions
) -> NDArray[np.float64]:
    """ionmc Bethe mass stopping [MeV cm2/g] with the options ``o``."""
    return bethe_mass_stopping(
        e_u,
        proj,
        material,
        shell=o.shell,
        barkas=o.barkas,
        bloch=o.bloch,
        mott=o.mott,
        effective_charge=o.effective_charge,
        ls=o.ls,
    )


def _options_dict(o: BetheOptions) -> dict[str, Any]:
    return {
        "shell": o.shell,
        "barkas": o.barkas,
        "bloch": o.bloch,
        "mott": o.mott,
        "effective_charge": o.effective_charge,
        "ls": o.ls,
    }


def charge_shape(e_u: NDArray[np.float64], proj: Projectile) -> NDArray[np.float64]:
    """1 - q^2/z^2 with the model's effective charge (q = z for z = 1, as in the model)."""
    if proj.z < 2:
        return np.zeros_like(e_u)
    beta = np.sqrt(kinematics(e_u, proj)[0])
    q = pierce_blann_charge(beta, proj)
    out: NDArray[np.float64] = 1.0 - (q / proj.z) ** 2
    return out


def band_means(
    e_u: NDArray[np.float64], values: NDArray[np.float64]
) -> dict[str, tuple[int, float]]:
    """``{band key: (n_points, mean of values over the band)}`` for all bands."""
    out = {}
    for band in BANDS:
        mask = _band_mask(e_u, band)
        out[_band_key(band)] = (int(np.count_nonzero(mask)), float(np.mean(values[mask])))
    return out


def ion_bands(
    e_u: NDArray[np.float64], s_model: NDArray[np.float64], s_table: NDArray[np.float64]
) -> dict[str, Any]:
    """Rounded band means of S_model/S_table - 1 (no per-energy values, no band max/min/RMS)."""
    means = band_means(e_u, s_model / s_table - 1.0)
    return {
        _band_key(band): {
            "n_points": means[_band_key(band)][0],
            "band_per_u_mev": list(band),
            "mean_relative_deviation": round_deviation(means[_band_key(band)][1]),
        }
        for band in BANDS
    }


def fit_shape(d: NDArray[np.float64], f: NDArray[np.float64]) -> dict[str, float | None]:
    """Least squares d = c f without intercept: coefficient, residual RMS, uncentred r^2."""
    ff = float(np.dot(f, f))
    c = float(np.dot(f, d)) / ff if ff > 0.0 else 0.0
    res = d - c * f
    dd = float(np.dot(d, d))
    return {
        "coefficient": c,
        "residual_rms": float(np.sqrt(np.mean(res**2))),
        "r_squared": 1.0 - float(np.dot(res, res)) / dd if dd > 0.0 else None,
    }


def attribution_verdict(
    d: NDArray[np.float64], fits: dict[str, dict[str, float | None]], ranking: list[str]
) -> dict[str, Any]:
    """Rule 1b of amendment I3-FU-1 for one band (``d`` the rounded band means over Z = 3-10).

    With RMS_1 <= RMS_2 the two smallest residual RMS and r^2_1 that of the best shape: attributed
    to the best shape iff RMS_1 <= 0.5 RMS_2, r^2_1 >= 0.90 and all d_z share one strict sign.
    """
    rms_1 = float(fits[ranking[0]]["residual_rms"] or 0.0)
    rms_2 = float(fits[ranking[1]]["residual_rms"] or 0.0)
    r2_1 = fits[ranking[0]]["r_squared"]
    common_sign = bool(np.all(d > 0.0) or np.all(d < 0.0))
    reason = []
    if not rms_1 <= ATTRIBUTION_RMS_RATIO * rms_2:
        reason.append("rms_ratio")
    if r2_1 is None or not r2_1 >= ATTRIBUTION_R_SQUARED:
        reason.append("r_squared")
    if not common_sign:
        reason.append("sign_change")
    return {
        "attributed_to": None if reason else ranking[0],
        "reason": reason or None,
        "rms_ratio_1_to_2": rms_1 / rms_2 if rms_2 > 0.0 else None,
        "common_sign": common_sign,
    }


def z_scan(
    means: dict[int, dict[str, float]], f_charge: dict[int, dict[str, float]]
) -> tuple[dict[str, Any], dict[str, Any]]:
    """z-scan fits over Z = 3-10 per band and the informative Z = 1, 2 entries.

    ``means`` must hold the ROUNDED band means (the values written to the output).
    """
    scan: dict[str, Any] = {}
    informative: dict[str, Any] = {"lineage_note": LINEAGE_NOTE, "bands": {}}
    zs = np.asarray(Z_SCAN_IONS, dtype=np.float64)
    for band in Z_SCAN_BANDS:
        key = _band_key(band)
        d = np.asarray([means[z][key] for z in Z_SCAN_IONS])
        shapes = {
            "charge": np.asarray([f_charge[z][key] for z in Z_SCAN_IONS]),
            "barkas_like": zs,
            "bloch_like": zs**2,
            "constant": np.ones_like(zs),
        }
        fits = {name: fit_shape(d, f) for name, f in shapes.items()}
        coef = {name: float(fit["coefficient"] or 0.0) for name, fit in fits.items()}
        ranking = sorted(fits, key=lambda name: float(fits[name]["residual_rms"] or 0.0))
        d6 = means[Z_REFERENCE][key]
        scan[key] = {
            "band_per_u_mev": list(band),
            "hypotheses": fits,
            "ranking": ranking,
            "shape_ratio_to_z6": {
                str(z): (means[z][key] / d6 if d6 != 0.0 else None) for z in Z_SCAN_IONS
            },
            "f_charge_band_mean": {str(z): f_charge[z][key] for z in Z_SCAN_IONS},
            "verdict": attribution_verdict(d, fits, ranking),
        }
        informative["bands"][key] = {
            "band_per_u_mev": list(band),
            **{
                f"z{z}": {
                    "mean_relative_deviation": means[z][key],
                    "f_charge_band_mean": f_charge[z][key],
                    "prediction_from_z3_10_fit": {
                        "charge": coef["charge"] * f_charge[z][key],
                        "barkas_like": coef["barkas_like"] * z,
                        "bloch_like": coef["bloch_like"] * z**2,
                        "constant": coef["constant"],
                    },
                }
                for z in (1, 2)
            },
        }
    return scan, informative


def hypothesis_check(ions: dict[str, Any]) -> dict[str, Any]:
    """Rule 1a of amendment I3-FU-1 on the rounded C/O band means of the ``ions`` output."""
    per_ion: dict[str, Any] = {}
    confirmed = True
    for name, quoted in HYPOTHESIS.items():
        bands = ions[name]["bands"]
        entry: dict[str, Any] = {}
        for band, value in quoted.items():
            committed = bands[f"{band}_MeV_u"]["mean_relative_deviation"]
            match = committed == value
            entry[f"{band}_MeV_u"] = {"quoted": value, "committed": committed, "match": match}
            confirmed &= match
        high = _band_key(BANDS[-1])
        committed = bands[high]["mean_relative_deviation"]
        match = abs(committed) <= HYPOTHESIS_HIGH_ABS[name]
        entry[high] = {
            "quoted_abs_bound": HYPOTHESIS_HIGH_ABS[name],
            "committed": committed,
            "match": match,
        }
        confirmed &= match
        per_ion[name] = entry
    return {
        "rule": (
            "amendment I3-FU-1 rule 1a: CONFIRMED iff, for carbon12 and oxygen16, every band "
            "mean of [2,5), [5,10), [10,20), [20,50) MeV/u rounded to 0.1 pp equals the quoted "
            "value (|delta| <= 0.05 pp before rounding) and |rounded mean of [50,1000]| <= the "
            "quoted bound (0.03 % C, 0.06 % O); otherwise WITHDRAWN. The bounds are below the "
            "0.1 pp quantum, so the last condition means a rounded mean of 0 (|mean| < 0.05 pp). "
            "The quoted >= 50 MeV/u figure was a max-type statement; a mean checks it only partly"
        ),
        "per_ion": per_ion,
        "verdict": "CONFIRMED" if confirmed else "WITHDRAWN",
    }


def i3b_ion(
    proj: Projectile,
    e_u: NDArray[np.float64],
    s_table: NDArray[np.float64],
    range_energies: tuple[float, ...],
    material: Material,
) -> dict[str, Any]:
    """Informative row I3b of one ion: model-only variant quantities, table mean for default only.

    Non-default variants carry no table-relative value (see the module docstring: several
    variant-vs-table band means would let the per-energy table values be solved for).
    """
    s_default = model_stopping(e_u, proj, material, VARIANTS["default"])
    mm_per_g_cm2 = 10.0 / material.density_g_cm3
    r_default = BetheStoppingSource(options=VARIANTS["default"]).table(material, proj)
    out: dict[str, Any] = {}
    for name, o in VARIANTS.items():
        s_v = model_stopping(e_u, proj, material, o)
        t_v = BetheStoppingSource(options=o).table(material, proj)
        ratio = band_means(e_u, s_v / s_default)
        bands: dict[str, dict[str, Any]] = {
            _band_key(band): {
                "n_points": ratio[_band_key(band)][0],
                "band_per_u_mev": list(band),
                "model_only_stopping_ratio_to_default_band_mean": ratio[_band_key(band)][1],
            }
            for band in BANDS
        }
        if name == "default":
            dev = band_means(e_u, s_v / s_table - 1.0)
            for key, band_out in bands.items():
                band_out["mean_relative_deviation"] = round_deviation(dev[key][1])
        out[name] = {
            "options": _options_dict(o),
            "bands": bands,
            "csda_range_shift_mm_vs_default": {
                f"{e:g}": float(t_v.range_at(e) - r_default.range_at(e)) * mm_per_g_cm2
                for e in range_energies
            },
        }
    return out


def mass_sensitivity(e_u: NDArray[np.float64], material: Material) -> float:
    """Model-only max |S(1.001 m)/S(m) - 1| over ``e_u`` for the ion ``MASS_SENSITIVITY_Z``."""
    p = PROJECTILES[MASS_SENSITIVITY_Z]
    heavier = Projectile(p.name, p.symbol, p.z, p.a, p.mass_mev * MASS_SENSITIVITY_FACTOR)
    o = BetheOptions()
    ratio = model_stopping(e_u, heavier, material, o) / model_stopping(e_u, p, material, o)
    return float(np.max(np.abs(ratio - 1.0)))


def run(
    archive: Path,
    code_sha: str,
    *,
    archive_sha256: str = ARCHIVE_SHA256,
    member_sha256: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Compute the band diagnostic as a JSON-serialisable dictionary."""
    from ionmc import __version__

    member_sha256 = MEMBER_SHA256 if member_sha256 is None else member_sha256
    tables = I3.load_icru73_water(
        archive,
        archive_sha256=archive_sha256,
        member_sha256=member_sha256,
        ions=tuple(PROJECTILES),
    )
    material = water(I_EV_MODEL)
    o = BetheOptions()
    ions: dict[str, Any] = {}
    means: dict[int, dict[str, float]] = {}
    f_charge: dict[int, dict[str, float]] = {}
    e_sel: NDArray[np.float64] | None = None
    for z, proj in PROJECTILES.items():
        table = tables[z]
        if table.z != proj.z:
            raise ValueError(f"table Z = {table.z} does not match projectile z = {proj.z}")
        sel = table.energy_per_u_mev >= E_MIN
        e_u = table.energy_per_u_mev[sel]
        s_table = table.stopping_mev_cm2_g[sel]
        e_sel = e_u
        s_model = model_stopping(e_u, proj, material, o)
        bands = ion_bands(e_u, s_model, s_table)
        entry: dict[str, Any] = {
            "z": proj.z,
            "a": proj.a,
            "mass_mev": proj.mass_mev,
            "bands": bands,
        }
        if z in I3_MAX_Z:
            abs_dev = np.abs(s_model / s_table - 1.0)
            entry["max_abs_relative_deviation_2_to_10"] = float(np.max(abs_dev[e_u < E_I3_SPLIT]))
            entry["max_abs_relative_deviation_ge_10"] = float(np.max(abs_dev[e_u >= E_I3_SPLIT]))
        ions[proj.name] = entry
        means[z] = {k: b["mean_relative_deviation"] for k, b in bands.items()}
        f_charge[z] = {k: v[1] for k, v in band_means(e_u, charge_shape(e_u, proj)).items()}
    assert e_sel is not None
    scan, informative = z_scan(means, f_charge)

    i3b: dict[str, Any] = {
        "note": (
            "informative; table-relative band means for the default variant only (duplicate of "
            "the rounded ions entry); all variants carry model-only quantities (CSDA-range shift "
            "vs the default model from BetheStoppingSource tables, 1-500 MeV/u, 200 points/decade; "
            "band-mean stopping ratio r_v to the default model). Variant-vs-table band means are "
            "not written: eight of them at full precision would let the per-energy table values "
            "of a band with <= 8 grid points be solved for; the variant-vs-table mean is "
            "approximately (1 + d_default) r_v - 1"
        ),
        "ions": {},
    }
    for name, (proj, energies) in I3B_IONS.items():
        t = tables[proj.z]
        sel = t.energy_per_u_mev >= E_MIN
        i3b["ions"][name] = {
            "standard_energies_per_u_mev": list(energies),
            "variants": i3b_ion(
                proj, t.energy_per_u_mev[sel], t.stopping_mev_cm2_g[sel], energies, material
            ),
        }

    sensitivity = mass_sensitivity(e_sel, material)
    mass_note = (
        "Li-7, Be-9, B-11, N-14, F-19 and Ne-20 use the nuclear mass A u - Z m_e (mass number "
        "times the atomic mass unit minus Z electron masses; mass excess and electron binding "
        "neglected). Model-only sensitivity: a 0.1 % mass increase changes the "
        f"{PROJECTILES[MASS_SENSITIVITY_Z].name} stopping by at most {sensitivity:.3g} (relative) "
        "over the grid at E >= 2 MeV/u."
    )
    return {
        "row": "I3-diagnostic",
        "status_note": STATUS_NOTE,
        "reference": "libdEdx ICRU73_NEW, target 276 (liquid water), Z = 1-10",
        "redistribution_note": REDISTRIBUTION_NOTE,
        "rounding_note": ROUNDING_NOTE,
        "i_value_note": I_VALUE_NOTE,
        "bands": {
            _band_key(b): {"band_per_u_mev": list(b), "upper_closed": b == BANDS[-1]} for b in BANDS
        },
        "ions": ions,
        "hypothesis_check": hypothesis_check(ions),
        "z_scan": {
            "ions_z": [Z_SCAN_IONS[0], Z_SCAN_IONS[-1]],
            "method": (
                "per band, least squares d_z = c f(z) without intercept over Z = 3-10; d_z the "
                "band mean of S_model/S_table - 1 rounded to 0.1 pp; charge: f = band mean of "
                "1 - q^2/z^2 (model effective charge); barkas_like: f = z; bloch_like: f = z^2; "
                "constant: f = 1 (z-independent control); "
                "r_squared uncentred (1 - sum res^2 / sum d^2); ranking by residual RMS"
            ),
            "verdict_rule": (
                "amendment I3-FU-1 rule 1b: with RMS_1 <= RMS_2 the two smallest residual RMS and "
                "r^2_1 that of the best shape, attributed_to = best iff RMS_1 <= 0.5 RMS_2, "
                "r^2_1 >= 0.90 and the rounded d_z over Z = 3-10 have one common sign (no zero, "
                "no sign change); otherwise attributed_to = null (no single-term scaling "
                "explains the z-dependence in this band) with reason in {rms_ratio, r_squared, "
                "sign_change}. Interpretive, not a grade"
            ),
            "bands": scan,
        },
        "informative_z1_z2": informative,
        "i3b": i3b,
        "provenance": {
            "archive_sha256": archive_sha256,
            "member_sha256": dict(member_sha256),
            "target": TARGET_WATER,
            "code_sha": code_sha,
            "ionmc_version": __version__,
            "numpy_version": np.__version__,
            "model": {
                "source": "bethe",
                "I_eV": I_EV_MODEL,
                "material": "water",
                "options": _options_dict(o),
            },
            "projectiles": {
                p.name: {
                    "z": p.z,
                    "a": p.a,
                    "mass_mev": p.mass_mev,
                    "mass": "A u - Z m_e" if z in APPROXIMATE_MASS_Z else "ionmc projectile",
                }
                for z, p in PROJECTILES.items()
            },
            "mass_note": mass_note,
            "mass_sensitivity_max_rel": sensitivity,
            "sealed_sources": None,
            "ionmc_loaded_from_repo": None,
            "environment": None,
            "expect_i3": None,
        },
    }


def _git(repo: Path, *args: str) -> bytes:
    """Run ``git -C repo args`` (no shell, ``GIT_*`` environment removed); stdout bytes."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    proc = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, env=env)
    if proc.returncode != 0:
        detail = proc.stderr.decode(errors="replace").strip()
        raise RuntimeError(f"git {' '.join(args)} failed in {repo}: {detail}")
    return proc.stdout


def committed_bytes(repo: Path, code_sha: str, rel: str) -> bytes:
    """Bytes of ``rel`` committed at ``code_sha``; the working copy in ``repo`` must match."""
    raw = _git(repo, "show", f"{code_sha}:{rel}")
    disk = Path(repo, rel)
    if not disk.is_file() or disk.read_bytes() != raw:
        raise RuntimeError(f"{rel} in {repo} differs from its committed bytes at {code_sha}")
    return raw


def verify_checkout(repo: Path, code_sha: str) -> None:
    """Fail closed unless ``repo`` is a clean checkout at HEAD == ``code_sha`` (40 lowercase hex).

    Raises ``RuntimeError`` on a malformed SHA, HEAD != ``code_sha`` or a modified tracked file.
    """
    if not isinstance(code_sha, str) or not CODE_SHA_RE.fullmatch(code_sha):
        raise RuntimeError(f"--code-sha must be 40 lowercase hex digits, got {code_sha!r}")
    head = _git(repo, "rev-parse", "HEAD").decode().strip()
    if head != code_sha:
        raise RuntimeError(f"HEAD of {repo} is {head}, not --code-sha {code_sha}")
    dirty = _git(repo, "status", "--porcelain", "--untracked-files=no").decode().strip()
    if dirty:
        raise RuntimeError(f"{repo} has modified tracked files:\n{dirty}")


def _python_satisfies(version: tuple[int, ...], spec: str) -> bool:
    """Whether ``version`` satisfies a ``requires-python`` specifier (>=, <=, ==, !=, >, <;
    ``==X.Y.*``/``!=X.Y.*`` prefix forms). Any other form raises ``RuntimeError``."""
    ok = True
    for part in spec.split(","):
        m = _SPECIFIER_RE.fullmatch(part)
        if m is None:
            raise RuntimeError(f"unsupported requires-python specifier {spec!r} in uv.lock")
        op, num, star = m.group(1), tuple(int(x) for x in m.group(2).split(".")), m.group(3)
        if star and op not in ("==", "!="):
            raise RuntimeError(f"unsupported requires-python specifier {spec!r} in uv.lock")
        width = len(num) if star else max(len(num), len(version))
        a = (version + (0,) * width)[:width]
        b = (num + (0,) * width)[:width]
        ok &= {
            ">=": a >= b,
            "<=": a <= b,
            ">": a > b,
            "<": a < b,
            "==": a == b,
            "!=": a != b,
        }[op]
    return ok


def check_environment(lock_raw: bytes) -> dict[str, Any]:
    """Running numpy and Python against the committed ``uv.lock``; ``RuntimeError`` on mismatch."""
    lock = tomllib.loads(lock_raw.decode("utf-8"))
    numpy_entries = [p for p in lock.get("package", []) if p.get("name") == "numpy"]
    if len(numpy_entries) != 1:
        raise RuntimeError(f"uv.lock has {len(numpy_entries)} numpy entries (need exactly 1)")
    locked = str(numpy_entries[0].get("version"))
    running = str(np.__version__)
    if running != locked:
        raise RuntimeError(f"running numpy {running} != numpy {locked} pinned in uv.lock")
    py = tuple(sys.version_info[:3])
    spec = lock.get("requires-python")
    python: dict[str, Any] = {
        "running": ".".join(map(str, py)),
        "major_minor": f"{py[0]}.{py[1]}",
        "requires_python": spec,
    }
    if spec is None:
        python["check"] = "recorded only (no requires-python in uv.lock)"
    else:
        if not _python_satisfies(py, str(spec)):
            raise RuntimeError(f"running Python {python['running']} violates uv.lock {spec!r}")
        python["check"] = "satisfies uv.lock requires-python"
    return {
        "numpy": {"running": running, "locked": locked, "match": True},
        "python": python,
        "note": ENVIRONMENT_NOTE,
    }


def ionmc_closure() -> dict[str, Path]:
    """``{repo-relative path: loaded file}`` of every ``ionmc`` module in ``sys.modules``.

    The path is relative to the directory that contains ``src/ionmc`` (two levels above the
    package directory); modules without ``__file__`` (namespace packages) are skipped.
    """
    import ionmc

    base = Path(str(ionmc.__file__)).resolve().parents[2]
    out: dict[str, Path] = {}
    for name, module in list(sys.modules.items()):
        if name != "ionmc" and not name.startswith("ionmc."):
            continue
        file = getattr(module, "__file__", None)
        if file is None:
            continue
        path = Path(file).resolve()
        out[path.relative_to(base).as_posix()] = path
    return out


def seal_sources(
    repo: Path, code_sha: str, extra: tuple[str, ...] = ()
) -> tuple[dict[str, str], bool]:
    """Seal the executable closure and the lock files (see the module docstring).

    Returns ``({path: sha256}, whether every loaded ionmc file lies under repo)``; raises
    ``RuntimeError`` if a loaded, on-disk or committed copy differs or a path is not committed.
    """
    closure = ionmc_closure()
    loaded: dict[str, Path] = {
        SCRIPT_REL: Path(__file__).resolve(),
        I3_SCRIPT_REL: Path(str(I3.__file__)).resolve(),
        ENV_SCRIPT_REL: Path(str(ENV.__file__)).resolve(),
        **closure,
    }
    sealed: dict[str, str] = {}
    for rel, path in sorted(loaded.items()):
        raw = committed_bytes(repo, code_sha, rel)
        if path.read_bytes() != raw:
            raise RuntimeError(f"loaded {rel} ({path}) differs from its committed bytes")
        sealed[rel] = hashlib.sha256(raw).hexdigest()
    for rel in (*LOCK_FILES, *extra):
        sealed[rel] = hashlib.sha256(committed_bytes(repo, code_sha, rel)).hexdigest()
    root = Path(repo).resolve()
    return sealed, all(path.is_relative_to(root) for path in closure.values())


def check_expect_i3(result: dict[str, Any], raw: bytes, rel: str, code_sha: str) -> dict[str, Any]:
    """Compare the C/O band maxima with the committed row I3 file (``compare_icru73.py`` output).

    ``raw`` are the committed bytes of ``rel`` at ``code_sha`` (:func:`committed_bytes`). Raises
    ``RuntimeError`` on any mismatch (rtol ``EXPECT_I3_RTOL``) or missing field; returns the
    provenance record ``{"path", "sha256", "committed_at", "matched": True}``.
    """
    i3 = json.loads(raw)
    fields = {
        "max_abs_relative_deviation_2_to_10": "2_to_10_MeV_u",
        "max_abs_relative_deviation_ge_10": "ge_10_MeV_u",
    }
    for z in I3_MAX_Z:
        name = PROJECTILES[z].name
        for ours, theirs in fields.items():
            try:
                expected = float(i3["ions"][name]["stopping"][theirs]["max_abs_relative_deviation"])
            except (KeyError, TypeError, ValueError) as exc:
                raise RuntimeError(f"row I3 file {rel} lacks {name} {theirs} maximum") from exc
            got = float(result["ions"][name][ours])
            if not np.isclose(got, expected, rtol=EXPECT_I3_RTOL, atol=0.0):
                raise RuntimeError(
                    f"{name} {ours} does not match the row I3 file {rel} (rtol {EXPECT_I3_RTOL})"
                )
    return {
        "path": rel,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "committed_at": code_sha,
        "matched": True,
    }


def main() -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    parser.add_argument("--archive", help="libdEdx archive (sha256-verified); required")
    parser.add_argument("--output-dir", help="required")
    parser.add_argument(
        "--code-sha",
        help="40-hex SHA that must equal HEAD of the clean --repo checkout; required",
    )
    parser.add_argument(
        "--environment-seal",
        default=ENV.SEAL_REL,
        help="repo-relative committed environment seal (read at --code-sha); the running "
        "environment must match every field",
    )
    parser.add_argument(
        "--write-environment-seal",
        action="store_true",
        help="write the seal of the running environment to --repo/--environment-seal and exit "
        "(nothing else is computed)",
    )
    parser.add_argument(
        "--repo",
        default=str(DEFAULT_REPO),
        help="repository root (default: the one containing this script)",
    )
    parser.add_argument(
        "--expect-i3",
        default=DEFAULT_EXPECT_I3,
        help="repo-relative row I3 JSON, read as committed at --code-sha; the C/O band maxima "
        "must match it (rtol 1e-9)",
    )
    args = parser.parse_args()
    repo = Path(args.repo)
    if args.write_environment_seal:
        target = repo / args.environment_seal
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(ENV.seal_text(ENV.current_environment()), encoding="utf-8")
        print(target)
        return 0
    missing = [f"--{k.replace('_', '-')}" for k in ("archive", "output_dir", "code_sha")
               if getattr(args, k) is None]  # fmt: skip
    if missing:
        parser.error(f"the following arguments are required: {', '.join(missing)}")
    verify_checkout(repo, args.code_sha)
    environment = check_environment(committed_bytes(repo, args.code_sha, "uv.lock"))
    environment["seal"] = ENV.verify_environment(
        committed_bytes(repo, args.code_sha, args.environment_seal)
    )
    i3_raw = committed_bytes(repo, args.code_sha, args.expect_i3)
    result = run(Path(args.archive), args.code_sha)
    # after the computation: the import closure is complete (nothing written before this passes)
    sealed, from_repo = seal_sources(repo, args.code_sha, (args.environment_seal,))
    sealed_modules = {
        name: mod
        for name, mod in sys.modules.items()
        if getattr(mod, "__file__", None)
        and (
            name == "ionmc" or name.startswith("ionmc.") or mod in (sys.modules[__name__], I3, ENV)
        )
    }
    environment["runtime"] = {
        **ENV.numpy_modules_check(),
        **ENV.bytecode_check({**sealed_modules, **ENV.numpy_modules()}),
    }
    result["provenance"]["sealed_sources"] = sealed
    result["provenance"]["ionmc_loaded_from_repo"] = from_repo
    result["provenance"]["environment"] = environment
    result["provenance"]["expect_i3"] = check_expect_i3(
        result, i3_raw, args.expect_i3, args.code_sha
    )
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "icru73_band_diagnostic.json"
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
