"""Compare the analytic stopping tables with cached NIST PSTAR/ASTAR and ICRU 90 data.

Usage::

    uv run python validation/scripts/stopping/compare_nist.py --output-dir DIR \
        [--cache-dir CACHE] [--online]

Requires the datasets ``nist-pstar-water-2005``, ``nist-astar-water-2005`` and
``geant4-icru90-stopping-11.4.2`` in the cache (``ionmc data fetch <id>``); no network
access is used unless ``--online`` is given. Writes ``stopping_comparison.json`` to
DIR. The file contains aggregates only (NIST SRD 124 content must not be reproducible from
it): for each comparison the number of points, the energy range, and the maximum and RMS
relative deviation (Bethe/reference - 1) of the electronic stopping power for energies per
nucleon >= 2 and >= 10 MeV/u (Bethe at I = 75 eV against PSTAR/ASTAR, Bethe at I = 78 eV
against the ICRU 90 arrays, ICRU 90 against NIST on shared energies); for CSDA ranges at
100/150/200/250 MeV protons and 100/150/200 MeV/u alpha particles the maximum absolute
relative deviation of Bethe (75 eV) from the NIST range and of Bethe (78 eV) from the ICRU 90
integral (200 points per decade, from 1 MeV/u, started from the NIST range at 1 MeV/u), and
the model-only Bethe 78 eV - 75 eV range shift in mm per energy.

Row I2 of ``validation/plans/v3-008-acceptance.md`` (He table reproduction; a regression
guard on SEEN values, not evidence) is evaluated under the top-level key ``"i2"``: the alpha
stopping-power aggregates (Bethe 75 eV vs ASTAR, Bethe 78 eV vs ICRU 90) in the separate bands
2 <= E < 10 MeV/u and 10 <= E <= 500 MeV/u, the maximum absolute relative CSDA-range deviation
over 50/100/150/200/250 MeV/u per reference together with the energy attaining it (no
per-energy deviation is written), the pre-registered gates (S <= 1.0 % for 2-10 MeV/u and
<= 0.30 % for >= 10 MeV/u per reference; range <= 0.20 % vs ASTAR and <= 0.25 % vs the
integrated ICRU 90 table), per-criterion booleans and an overall ``pass``, the statement of the
ASTAR (75 eV) versus model-default (78 eV) I-value difference and the model-only Bethe
78 eV - 75 eV range shift at the five energies.

NIST tables are SRD 124 data (attribution required) and are not redistributed here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from ionmc.data.acquire import fetch
from ionmc.data.icru90 import load_icru90_water
from ionmc.data.nist_star import StarTable, load_star_table
from ionmc.materials import water
from ionmc.physics.projectiles import ALPHA, PROTON, Projectile
from ionmc.physics.stopping import BetheStoppingSource, StoppingTable, build_table, log_grid

E_COMPARE_MIN = 2.0  # MeV/u
E_MAX = 500.0  # MeV/u
PROTON_ENERGIES = (100.0, 150.0, 200.0, 250.0)  # MeV
ALPHA_ENERGIES_PER_U = (100.0, 150.0, 200.0)  # MeV/u
# Row I2 (validation/plans/v3-008-acceptance.md): CSDA-range energies and pre-registered gates
# (bounds = 1.5x the SEEN decision-0038 maxima, rounded up; a regression guard, not evidence).
I2_ENERGIES_PER_U = (50.0, 100.0, 150.0, 200.0, 250.0)  # MeV/u
I2_BAND_LOW = "2_to_10_MeV_u"  # 2 <= E < 10 MeV/u
I2_BAND_HIGH = "ge_10_MeV_u"  # 10 <= E <= E_MAX MeV/u
I2_S_MAX_ABS = {I2_BAND_LOW: 0.010, I2_BAND_HIGH: 0.0030}
I2_RANGE_MAX_ABS = {"vs_astar_75eV": 0.0020, "vs_icru90_integral_78eV": 0.0025}
I2_I_VALUE_STATEMENT = (
    "NIST ASTAR (SRD 124) liquid water uses I = 75 eV (ICRU Report 49); the ionmc default "
    "water uses I = 78 eV (ICRU Report 90). Bethe is therefore compared at I = 75 eV with "
    "ASTAR and at the default I = 78 eV with ICRU 90 (stopping and integrated CSDA range); "
    "the model-only Bethe 78 eV - 75 eV CSDA-range shift at the I2 energies is reported in "
    "bethe_78eV_minus_75eV_range_shift_mm_by_energy_per_u_mev."
)


def _nist_range(star: StarTable, e_total_mev: float) -> float:
    """NIST CSDA range [g/cm2] at total kinetic energy [MeV] (log-log interpolation)."""
    return float(
        np.exp(np.interp(np.log(e_total_mev), np.log(star.energy_mev), np.log(star.csda_range)))
    )


def window_mask(e_u: np.ndarray) -> np.ndarray:
    """Boolean mask of energies per nucleon inside the comparison window [2, E_MAX] MeV/u."""
    return (e_u >= E_COMPARE_MIN) & (e_u <= E_MAX)


def _aggregate(e_u: np.ndarray, dev: np.ndarray) -> dict[str, Any]:
    """Non-reversible summary of relative deviations ``dev`` at energies per nucleon ``e_u``."""
    out: dict[str, Any] = {
        "n_points": int(e_u.size),
        "energy_per_u_range_mev": [float(e_u.min()), float(e_u.max())],
    }
    for label, cut in (("ge_2_MeV_u", 2.0), ("ge_10_MeV_u", 10.0)):
        m = e_u >= cut
        out[f"max_abs_relative_deviation_{label}"] = float(np.max(np.abs(dev[m])))
        out[f"rms_relative_deviation_{label}"] = float(np.sqrt(np.mean(dev[m] ** 2)))
        out[f"n_points_{label}"] = int(np.count_nonzero(m))
    return out


def band_aggregate(
    e_u: np.ndarray, dev: np.ndarray, lo: float, hi: float, *, include_hi: bool
) -> dict[str, Any]:
    """Max abs and RMS of ``dev`` for ``lo <= e_u < hi`` (``<= hi`` if ``include_hi``).

    An empty band gives ``n_points = 0`` and ``None`` statistics (gates on it fail closed).
    """
    m = (e_u >= lo) & ((e_u <= hi) if include_hi else (e_u < hi))
    n = int(np.count_nonzero(m))
    out: dict[str, Any] = {
        "band_per_u_mev": [lo, hi],
        "upper_bound_inclusive": include_hi,
        "n_points": n,
        "energy_per_u_range_mev": None,
        "max_abs_relative_deviation": None,
        "rms_relative_deviation": None,
    }
    if n:
        out["energy_per_u_range_mev"] = [float(e_u[m].min()), float(e_u[m].max())]
        out["max_abs_relative_deviation"] = float(np.max(np.abs(dev[m])))
        out["rms_relative_deviation"] = float(np.sqrt(np.mean(dev[m] ** 2)))
    return out


def _deviation(
    proj: Projectile, e_total: np.ndarray, s_ref: np.ndarray, i_ev: float
) -> tuple[np.ndarray, np.ndarray]:
    """Energies per nucleon in the window and relative deviations Bethe(i_ev)/reference - 1."""
    e_u = e_total / proj.a
    sel = window_mask(e_u)
    s = BetheStoppingSource().table(water(i_ev), proj)
    return e_u[sel], s.stopping_at(e_u[sel]) / s_ref[sel] - 1.0


def _stopping_deviation(
    proj: Projectile, e_total: np.ndarray, s_ref: np.ndarray, i_ev: float
) -> dict[str, Any]:
    """Aggregate relative deviation of the Bethe model from reference stopping powers."""
    e_u, dev = _deviation(proj, e_total, s_ref, i_ev)
    return {"I_eV": i_ev, **_aggregate(e_u, dev)}


def _icru_range_table(
    proj: Projectile, star: StarTable, e_icru: np.ndarray, s_icru: np.ndarray
) -> StoppingTable:
    """ICRU 90 range table: log-log interpolation of the coarse grid onto 200 points/decade from
    1 MeV/u, integrated like the analytic tables (`build_table`), started at the NIST range."""
    e_u_icru = e_icru / proj.a
    e_fine = log_grid(1.0, min(E_MAX, float(e_u_icru[-1])), 200)
    s_fine = np.exp(np.interp(np.log(e_fine), np.log(e_u_icru), np.log(s_icru)))
    r_start = _nist_range(star, float(e_fine[0] * proj.a))
    return build_table(proj, water(78.0), e_fine, s_fine, r_start, {"source": "icru90"})


def _max_abs_at(energies: tuple[float, ...], dev: list[float]) -> dict[str, float]:
    """Maximum absolute deviation and the energy attaining it (no per-energy values)."""
    k = int(np.argmax(np.abs(dev)))
    return {
        "max_abs_relative_deviation": float(abs(dev[k])),
        "energy_per_u_mev_at_max": float(energies[k]),
    }


def i2_gates(
    stopping: dict[str, dict[str, Any]], range_dev: dict[str, dict[str, float]]
) -> dict[str, Any]:
    """Evaluate the pre-registered I2 gates; a missing or empty aggregate fails closed.

    ``stopping[ref][band]`` are :func:`band_aggregate` results for ``ref`` in
    ``("vs_astar_75eV", "vs_icru90_78eV")``; ``range_dev[ref]`` are :func:`_max_abs_at` results
    for the keys of ``I2_RANGE_MAX_ABS``.
    """
    criteria: dict[str, bool] = {}
    for ref in ("vs_astar_75eV", "vs_icru90_78eV"):
        for band, bound in I2_S_MAX_ABS.items():
            value = stopping.get(ref, {}).get(band, {}).get("max_abs_relative_deviation")
            criteria[f"stopping_{ref}_{band}_max_abs_le_{bound:g}"] = (
                value is not None and bool(np.isfinite(value)) and value <= bound
            )
    for ref, bound in I2_RANGE_MAX_ABS.items():
        value = range_dev.get(ref, {}).get("max_abs_relative_deviation")
        criteria[f"csda_range_{ref}_max_abs_le_{bound:g}"] = (
            value is not None and bool(np.isfinite(value)) and value <= bound
        )
    return {"criteria": criteria, "pass": all(criteria.values())}


def i2_block(astar: StarTable, e_icru: np.ndarray, s_icru: np.ndarray) -> dict[str, Any]:
    """Row I2 (He table reproduction, regression guard): aggregates, gates and I-value note.

    ``astar`` is the ASTAR water table (total alpha energy); ``e_icru``/``s_icru`` the ICRU 90
    alpha arrays (total alpha energy [MeV], MeV cm2/g). Only aggregates are returned.
    """
    bands = {I2_BAND_LOW: (2.0, 10.0, False), I2_BAND_HIGH: (10.0, E_MAX, True)}
    stopping: dict[str, dict[str, Any]] = {}
    for ref, (e_tot, s_ref, i_ev) in {
        "vs_astar_75eV": (astar.energy_mev, astar.s_electronic, 75.0),
        "vs_icru90_78eV": (e_icru, s_icru, 78.0),
    }.items():
        e_u, dev = _deviation(ALPHA, e_tot, s_ref, i_ev)
        stopping[ref] = {"I_eV": i_ev}
        for band, (lo, hi, inc) in bands.items():
            stopping[ref][band] = band_aggregate(e_u, dev, lo, hi, include_hi=inc)
    t75 = BetheStoppingSource().table(water(75.0), ALPHA)
    t78 = BetheStoppingSource().table(water(78.0), ALPHA)
    t_icru = _icru_range_table(ALPHA, astar, e_icru, s_icru)
    dev_astar, dev_icru, shift = [], [], {}
    for e_u in I2_ENERGIES_PER_U:
        r75 = float(t75.range_at(e_u))
        r78 = float(t78.range_at(e_u))
        dev_astar.append(r75 / _nist_range(astar, e_u * ALPHA.a) - 1.0)
        dev_icru.append(r78 / float(t_icru.range_at(e_u)) - 1.0)
        shift[f"{e_u:g}"] = (r78 - r75) * 10.0
    range_dev = {
        "vs_astar_75eV": _max_abs_at(I2_ENERGIES_PER_U, dev_astar),
        "vs_icru90_integral_78eV": _max_abs_at(I2_ENERGIES_PER_U, dev_icru),
    }
    gates = i2_gates(stopping, range_dev)
    return {
        "row": "I2",
        "role": "regression guard on SEEN values (decision 0038), not evidence",
        "stopping": stopping,
        "csda_range": {
            "energies_per_u_mev": list(I2_ENERGIES_PER_U),
            "construction": (
                "Bethe tables 200 points/decade from 1 MeV/u (constant-S start); ICRU 90 "
                "log-log interpolated to 200 points/decade from 1 MeV/u, started at the ASTAR "
                "range at 1 MeV/u; ASTAR range log-log interpolated"
            ),
            **range_dev,
            "bethe_78eV_minus_75eV_range_shift_mm_by_energy_per_u_mev": shift,
        },
        "thresholds": {
            "stopping_max_abs": dict(I2_S_MAX_ABS),
            "csda_range_max_abs": dict(I2_RANGE_MAX_ABS),
        },
        "gates": gates["criteria"],
        "pass": gates["pass"],
        "i_value_statement": I2_I_VALUE_STATEMENT,
    }


def run(cache_dir: str | None, offline: bool) -> dict[str, Any]:
    """Compute the aggregate comparison as a JSON-serialisable dictionary."""
    pstar = load_star_table(fetch("nist-pstar-water-2005", cache_dir, offline=offline))
    astar = load_star_table(fetch("nist-astar-water-2005", cache_dir, offline=offline))
    icru = load_icru90_water(fetch("geant4-icru90-stopping-11.4.2", cache_dir, offline=offline))
    result: dict[str, Any] = {
        "redistribution_note": (
            "This file contains aggregate statistics only; no NIST SRD 124 or ICRU 90 "
            "table values, per-energy deviations or ratios are reproduced."
        ),
        "stopping": {},
        "range": {},
    }
    result["i2"] = i2_block(astar, icru.alpha_energy_mev, icru.alpha_stopping)
    refs = {
        "proton": (PROTON, pstar, icru.proton_energy_mev, icru.proton_stopping),
        "alpha": (ALPHA, astar, icru.alpha_energy_mev, icru.alpha_stopping),
    }
    for name, (proj, star, e_icru, s_icru) in refs.items():
        result["stopping"][name] = {
            "vs_nist_75eV": _stopping_deviation(proj, star.energy_mev, star.s_electronic, 75.0),
            "vs_icru90_78eV": _stopping_deviation(proj, e_icru, s_icru, 78.0),
        }
        # ICRU 90 against NIST on the shared energies (checks the energy-grid reading).
        common = np.intersect1d(star.energy_mev, e_icru)
        common = common[window_mask(common / proj.a)]
        ratio = np.array(
            [
                s_icru[np.searchsorted(e_icru, e)]
                / star.s_electronic[np.searchsorted(star.energy_mev, e)]
                for e in common
            ]
        )
        result["stopping"][name]["icru90_vs_nist_shared_energies"] = _aggregate(
            common / proj.a, ratio - 1.0
        )
        t75 = BetheStoppingSource().table(water(75.0), proj)
        t78 = BetheStoppingSource().table(water(78.0), proj)
        # ICRU 90 range: log-log interpolation of its coarse grid onto a fine log grid,
        # then the same exact log-log quadrature as for the analytic tables (`build_table`).
        t_icru = _icru_range_table(proj, star, e_icru, s_icru)
        energies_u = PROTON_ENERGIES if proj is PROTON else ALPHA_ENERGIES_PER_U
        dev75, dev78, shift = [], [], {}
        for e_u in energies_u:
            r75 = float(t75.range_at(e_u))
            r78 = float(t78.range_at(e_u))
            dev75.append(r75 / _nist_range(star, e_u * proj.a) - 1.0)
            dev78.append(r78 / float(t_icru.range_at(e_u)) - 1.0)
            shift[f"{e_u:g}"] = (r78 - r75) * 10.0
        result["range"][name] = {
            "energies_per_u_mev": list(energies_u),
            "max_abs_relative_deviation_bethe75_vs_nist": float(np.max(np.abs(dev75))),
            "max_abs_relative_deviation_bethe78_vs_icru90_integral": float(np.max(np.abs(dev78))),
            "bethe_78eV_minus_75eV_shift_mm_by_energy_per_u_mev": shift,
        }
    return result


def provenance(cache_dir: str | None, code_sha: str | None) -> dict[str, Any]:
    """Dataset identities (id, version, pinned sha256, verified cached bytes) and versions."""
    from ionmc import __version__
    from ionmc.data.cache import read_manifest, resolve_cache_dir
    from ionmc.data.registry import DATASETS

    root = resolve_cache_dir(Path(cache_dir) if cache_dir else None)
    datasets = {}
    ids = ("nist-pstar-water-2005", "nist-astar-water-2005", "geant4-icru90-stopping-11.4.2")
    for dataset_id in ids:
        dataset = DATASETS[dataset_id]
        manifest = read_manifest(dataset_id, root)
        datasets[dataset_id] = {
            "version": dataset.version,
            "pinned_sha256": dataset.sha256,
            "cached_sha256": manifest.get("sha256") if manifest else None,
            "retrieved_at": manifest.get("retrieved_at") if manifest else None,
        }
    return {
        "ionmc_version": __version__,
        "numpy_version": np.__version__,
        "code_sha": code_sha,
        "datasets": datasets,
        "i_values_eV": {"bethe_vs_nist": 75.0, "bethe_vs_icru90": 78.0},
    }


def main() -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--cache-dir", default=None)
    parser.add_argument("--online", action="store_true", help="allow downloading missing data")
    parser.add_argument(
        "--code-sha",
        default=None,
        help="git SHA of the clean checkout that produced this comparison (recorded verbatim)",
    )
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    result = run(args.cache_dir, offline=not args.online)
    result["provenance"] = provenance(args.cache_dir, args.code_sha)
    path = out / "stopping_comparison.json"
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
