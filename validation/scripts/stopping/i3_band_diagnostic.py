"""Row I3 diagnostic: band-wise ionmc Bethe vs libdEdx ICRU73_NEW water for Z = 1-10 (no gate).

Usage::

    uv run python validation/scripts/stopping/i3_band_diagnostic.py --archive PATH \
        --output-dir DIR --code-sha SHA

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
  with the default :class:`BetheOptions` at ``water(78.0)``. No band max/min/RMS is written (with
  2-4 points per band they would, with the mean, expose per-energy deviations); a band with fewer
  than two grid points raises ``ValueError`` (its mean would be a per-energy value). Carbon and
  oxygen additionally carry ``max_abs_relative_deviation_2_to_10`` (2 <= E < 10 MeV/u) and
  ``max_abs_relative_deviation_ge_10`` (E >= 10 MeV/u to the grid end), the two band maxima the
  committed row I3 file already publishes; ``--expect-i3 PATH`` checks them against that file
  (rtol 1e-9; ``RuntimeError`` and no output on mismatch) and records its path and sha256.
* ``z_scan``: for the bands below 50 MeV/u, single-term least squares without intercept of the
  band means d_z over Z = 3-10 (the ICRU 73 ion set) against four shapes: ``charge``
  f = band mean of 1 - q^2/z^2 (q the model's effective charge, ``pierce_blann_charge``),
  ``barkas_like`` f = z, ``bloch_like`` f = z^2 and the z-independent control ``constant``
  f = 1. Coefficient, residual RMS and the uncentred
  r^2 = 1 - sum(res^2)/sum(d^2) (the no-intercept convention), a ranking by residual RMS and the
  ratios d_z/d_6. Z = 1 and 2 are reported separately (``informative_z1_z2``): their table
  lineage is OPEN.
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

The JSON output contains aggregates and model-only quantities only; no libdEdx/ICRU table
values, per-energy deviations or ratios are written (field ``redistribution_note``).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
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

STATUS_NOTE = (
    "diagnostic of the row I3 FAIL under amendment I3-FU-1; no gate; row I3 is not re-graded"
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
    """Band means of S_model/S_table - 1 (no per-energy values, no band max/min/RMS)."""
    means = band_means(e_u, s_model / s_table - 1.0)
    return {
        _band_key(band): {
            "n_points": means[_band_key(band)][0],
            "band_per_u_mev": list(band),
            "mean_relative_deviation": means[_band_key(band)][1],
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


def z_scan(
    means: dict[int, dict[str, float]], f_charge: dict[int, dict[str, float]]
) -> tuple[dict[str, Any], dict[str, Any]]:
    """z-scan fits over Z = 3-10 per band and the informative Z = 1, 2 entries."""
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
                band_out["mean_relative_deviation"] = dev[key][1]
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
    code_sha: str | None,
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
            "the ions entry); all variants carry model-only quantities (CSDA-range shift vs the "
            "default model from BetheStoppingSource tables, 1-500 MeV/u, 200 points/decade; "
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
        "i_value_note": I_VALUE_NOTE,
        "bands": {
            _band_key(b): {"band_per_u_mev": list(b), "upper_closed": b == BANDS[-1]} for b in BANDS
        },
        "ions": ions,
        "z_scan": {
            "ions_z": [Z_SCAN_IONS[0], Z_SCAN_IONS[-1]],
            "method": (
                "per band, least squares d_z = c f(z) without intercept over Z = 3-10; d_z the "
                "band mean of S_model/S_table - 1; charge: f = band mean of 1 - q^2/z^2 "
                "(model effective charge); barkas_like: f = z; bloch_like: f = z^2; "
                "constant: f = 1 (z-independent control); "
                "r_squared uncentred (1 - sum res^2 / sum d^2); ranking by residual RMS"
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
            "expect_i3": None,
        },
    }


def check_expect_i3(result: dict[str, Any], path: Path) -> dict[str, Any]:
    """Compare the C/O band maxima with a row I3 file (``compare_icru73.py`` output).

    Raises ``RuntimeError`` on any mismatch (rtol ``EXPECT_I3_RTOL``) or missing field; returns
    the provenance record ``{"path", "sha256", "matched": True}``.
    """
    raw = Path(path).read_bytes()
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
                raise RuntimeError(f"row I3 file {path} lacks {name} {theirs} maximum") from exc
            got = float(result["ions"][name][ours])
            if not np.isclose(got, expected, rtol=EXPECT_I3_RTOL, atol=0.0):
                raise RuntimeError(
                    f"{name} {ours} does not match the row I3 file {path} (rtol {EXPECT_I3_RTOL})"
                )
    return {"path": str(path), "sha256": hashlib.sha256(raw).hexdigest(), "matched": True}


def main() -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    parser.add_argument("--archive", required=True, help="libdEdx archive (sha256-verified)")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--code-sha",
        default=None,
        help="git SHA of the clean checkout that produced this diagnostic (recorded verbatim)",
    )
    parser.add_argument(
        "--expect-i3",
        default=None,
        help="row I3 JSON (compare_icru73.py); the C/O band maxima must match it (rtol 1e-9)",
    )
    args = parser.parse_args()
    result = run(Path(args.archive), args.code_sha)
    if args.expect_i3 is not None:
        result["provenance"]["expect_i3"] = check_expect_i3(result, Path(args.expect_i3))
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "icru73_band_diagnostic.json"
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
