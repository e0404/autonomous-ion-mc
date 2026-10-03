"""Compare the analytic stopping tables with cached NIST PSTAR/ASTAR and ICRU 90 data.

Usage::

    uv run python validation/scripts/stopping/compare_nist.py --output-dir DIR \
        [--cache-dir CACHE] [--online]

Requires the datasets ``nist-pstar-water-2005``, ``nist-astar-water-2005`` and
``geant4-icru90-stopping-11.4.2`` in the cache (``ionmc data fetch <id>``); no network
access is used unless ``--online`` is given. Writes ``stopping_comparison.json`` to
DIR, containing

* relative deviations (Bethe/reference - 1) of the electronic stopping power on the
  reference energy grid for energies per nucleon >= 2 MeV/u: Bethe at I = 75 eV against
  PSTAR/ASTAR and Bethe at I = 78 eV against the ICRU 90 arrays;
* CSDA range differences at 100/150/200/250 MeV protons and 100/150/200 MeV/u alpha
  particles: Bethe (I = 75 eV) against the NIST CSDA range; Bethe (I = 78 eV) against the
  range obtained by integrating the ICRU 90 stopping powers (log-log interpolated onto a
  grid of 200 points per decade) from 1 MeV/u with the same trapezoid rule, starting from
  the NIST CSDA range at 1 MeV/u; and the Bethe 78 eV - 75 eV range shift in mm.

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
from ionmc.physics.stopping import BetheStoppingSource, build_table, log_grid

E_COMPARE_MIN = 2.0  # MeV/u
E_MAX = 500.0  # MeV/u
PROTON_ENERGIES = (100.0, 150.0, 200.0, 250.0)  # MeV
ALPHA_ENERGIES_PER_U = (100.0, 150.0, 200.0)  # MeV/u


def _nist_range(star: StarTable, e_total_mev: float) -> float:
    """NIST CSDA range [g/cm2] at total kinetic energy [MeV] (log-log interpolation)."""
    return float(
        np.exp(np.interp(np.log(e_total_mev), np.log(star.energy_mev), np.log(star.csda_range)))
    )


def _stopping_deviation(
    proj: Projectile, e_total: np.ndarray, s_ref: np.ndarray, i_ev: float
) -> dict[str, Any]:
    """Relative deviation of the Bethe model from reference stopping powers."""
    e_u = e_total / proj.a
    sel = (e_u >= E_COMPARE_MIN) & (e_u <= E_MAX)
    s = BetheStoppingSource().table(water(i_ev), proj)
    dev = s.stopping_at(e_u[sel]) / s_ref[sel] - 1.0
    above10 = e_u[sel] >= 10.0
    return {
        "I_eV": i_ev,
        "energy_per_u_mev": e_u[sel].tolist(),
        "relative_deviation": dev.tolist(),
        "max_abs_deviation_all_ge_2_MeV_u": float(np.max(np.abs(dev))),
        "max_abs_deviation_ge_10_MeV_u": float(np.max(np.abs(dev[above10]))),
    }


def run(cache_dir: str | None, offline: bool) -> dict[str, Any]:
    """Compute the comparison and return it as a JSON-serialisable dictionary."""
    pstar = load_star_table(fetch("nist-pstar-water-2005", cache_dir, offline=offline))
    astar = load_star_table(fetch("nist-astar-water-2005", cache_dir, offline=offline))
    icru = load_icru90_water(fetch("geant4-icru90-stopping-11.4.2", cache_dir, offline=offline))
    result: dict[str, Any] = {"stopping": {}, "range": {}}
    refs = {
        "proton": (PROTON, pstar, icru.proton_energy_mev, icru.proton_stopping),
        "alpha": (ALPHA, astar, icru.alpha_energy_mev, icru.alpha_stopping),
    }
    for name, (proj, star, e_icru, s_icru) in refs.items():
        result["stopping"][name] = {
            "vs_nist_75eV": _stopping_deviation(proj, star.energy_mev, star.s_electronic, 75.0),
            "vs_icru90_78eV": _stopping_deviation(proj, e_icru, s_icru, 78.0),
        }
        # NIST ASTAR vs ICRU 90 alpha / PSTAR vs ICRU 90 proton on common energies (data check).
        common = np.intersect1d(star.energy_mev, e_icru)
        ratio = [
            float(
                s_icru[np.searchsorted(e_icru, e)]
                / star.s_electronic[np.searchsorted(star.energy_mev, e)]
            )
            for e in common
            if e / proj.a >= E_COMPARE_MIN
        ]
        result["stopping"][name]["icru90_over_nist_on_common_energies"] = {
            "energy_per_u_mev": [float(e / proj.a) for e in common if e / proj.a >= E_COMPARE_MIN],
            "ratio": ratio,
        }
        # Range comparison.
        t75 = BetheStoppingSource().table(water(75.0), proj)
        t78 = BetheStoppingSource().table(water(78.0), proj)
        # ICRU 90 range: log-log interpolation of its coarse grid onto a fine log grid,
        # then the same trapezoid integral as for the analytic tables.
        e_u_icru = e_icru / proj.a
        e_fine = log_grid(1.0, min(E_MAX, float(e_u_icru[-1])), 200)
        s_fine = np.exp(np.interp(np.log(e_fine), np.log(e_u_icru), np.log(s_icru)))
        r_start = _nist_range(star, float(e_fine[0] * proj.a))
        t_icru = build_table(proj, water(78.0), e_fine, s_fine, r_start, {"source": "icru90"})
        energies_u = PROTON_ENERGIES if proj is PROTON else ALPHA_ENERGIES_PER_U
        rows = []
        for e_u in energies_u:
            r_nist = _nist_range(star, e_u * proj.a)
            r75 = float(t75.range_at(e_u))
            r78 = float(t78.range_at(e_u))
            r_icru = float(t_icru.range_at(e_u))
            rows.append(
                {
                    "energy_per_u_mev": e_u,
                    "nist_csda_g_cm2": r_nist,
                    "bethe_75eV_g_cm2": r75,
                    "bethe_75eV_over_nist_minus_1": r75 / r_nist - 1.0,
                    "icru90_integrated_g_cm2": r_icru,
                    "bethe_78eV_g_cm2": r78,
                    "bethe_78eV_over_icru90_minus_1": r78 / r_icru - 1.0,
                    "bethe_78eV_minus_75eV_mm": (r78 - r75) * 10.0,
                    "icru90_minus_nist_mm": (r_icru - r_nist) * 10.0,
                }
            )
        result["range"][name] = rows
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
