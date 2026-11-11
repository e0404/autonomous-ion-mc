"""D-3 range attribution (V3-008 amendment R2-STOP3, item 1(b)): report-only diagnostic.

Compares ionmc's Bethe CSDA range with (He) the integrated ICRU 90 alpha array and (C, O) a
TOPAS-like rebuild of the Geant4 GenericIon stopping power, and estimates the lengthening of a
range from pre-step-energy stopping at 0.1 mm steps. It then splits the committed TOPAS-EM
R80 excess over ionmc's R80 into these pieces. Exploratory attribution only: never a pass
criterion, never construction, nothing is fitted.

Conventions (all ranges from E_CUT = 2 MeV/u, water rho = 1 g/cm3, ionmc's exact log-log
increments on a 400 points/decade grid, as in ``p3_driver.range_shift``):

* ``R_ionmc``: integral of ionmc's Bethe stopping (all options on, I = 78 eV, ls "off").
  ``R_ionmc_table_mm`` is the committed table range (``BetheStoppingSource`` default grid, with
  its constant-S start range at 1 MeV/u), reported for orientation.
* He: ICRU 90 alpha water array, log-log interpolated, total energy = A E/u.
* C, O ("TOPAS-like", an APPROXIMATION of the Geant4 construction, not a transcription):
  ``S = S_p(E_p) q^2 (1 + (dL_LS + dL_Bark)/L_p)`` with ``E_p`` the proton energy of equal
  velocity, ``S_p`` the ICRU 90 proton array (log-log), ``q`` ionmc's effective charge,
  ``L_p = S_p beta^2 / (K Z/A)``,
  ``dL_LS = dL_LS(Z bare, beta) - Bloch(beta^2, 1) - Mott1(beta, 1)``
  (the ICRU 90 proton value already contains the proton's Bloch and Mott terms; ionmc's
  ``bare_mixed`` indexing: LS at the nominal Z), ``dL_Bark = barkas_term(beta^2, q) -
  barkas_term(beta^2, 1)`` (ionmc's Barkas form is linear in q: the extra (q - 1) L1).
  Further sensitivity variants: "barkas_only" and "halfls_barkas" (dL_LS x 0.5, the exploratory
  half-magnitude reading of the Geant4 LS table, p3-lindhard-sorensen.md; unverified).
  The LS correction is evaluated on a 60 points/decade grid and interpolated linearly in ln E.
  Variant "plain": the same without the (dL_LS + dL_Bark) factor.
* Pre-step estimate: forward Euler with step 0.1 mm, S at the step's entry energy, vs the exact
  integral of the same S, as a fraction of the range. Primary ``prestep_euler_0p1mm_pct``:
  steps use pre-step S while the full-step loss is <= 2 % of the energy (Geant4's limit), the
  rest of the range (E at the first violation down to E_CUT) is added exactly;
  ``prestep_pure_euler_0p1mm_pct`` is unconditional Euler (upper bound). Planning estimate only.

Outputs aggregates only (ranges in mm, ratios in percent, hashes); no stopping values.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
WORKTREE = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(WORKTREE / "src"))
import lindhard_sorensen as lsm  # noqa: E402

from ionmc.data.icru90 import load_icru90_water  # noqa: E402
from ionmc.materials import water  # noqa: E402
from ionmc.physics import projectiles as pj  # noqa: E402
from ionmc.physics import stopping as st  # noqa: E402

E_CUT = 2.0
POINTS_PER_DECADE = 400
LS_POINTS_PER_DECADE = 60
STEP_MM = 0.1
ICRU_ID = "geant4-icru90-stopping-11.4.2"
CASES = (  # (ion, MeV/u, TOPAS-EM report case number or None for held-out energies)
    ("he", 50.0, 1),
    ("he", 100.0, None),
    ("he", 150.0, 2),
    ("he", 250.0, 3),
    ("carbon", 80.0, 4),
    ("carbon", 150.0, None),
    ("carbon", 290.0, 5),
    ("carbon", 430.0, 6),
    ("oxygen", 100.0, 7),
    ("oxygen", 200.0, None),
    ("oxygen", 350.0, 8),
    ("oxygen", 430.0, 9),
)
PROJ = {"he": pj.ALPHA, "carbon": pj.CARBON12, "oxygen": pj.OXYGEN16}
REPORT_DIR = "validation/results/transport/v8r2-topas/reports"
SOURCES = (
    "src/ionmc/physics/stopping.py",
    "src/ionmc/data/icru90.py",
    "validation/scripts/ion_em/lindhard_sorensen.py",
)


ICRU_SHA = "afd1479497a605d9675c1fe2a3870ddb3f305110fb2d76e13bb600fa3ab281cd"


def find_icru(cache_dir: Path) -> Path:
    """Registered ICRU 90 source object (content-addressed; the loader re-verifies the hash)."""
    for base in (cache_dir, WORKTREE.parent / "v3-005b" / ".ionmc-cache"):
        cand = base / "ionmc-data" / "objects" / ICRU_SHA
        if cand.is_file():
            return cand
    raise FileNotFoundError(f"object {ICRU_SHA} not found; run 'ionmc data fetch {ICRU_ID}'")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def loglog(x: np.ndarray, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    if np.any(x < xs[0] * (1 - 1e-12)) or np.any(x > xs[-1] * (1 + 1e-12)):
        raise ValueError("log-log interpolation outside the registered array")
    return np.exp(np.interp(np.log(x), np.log(xs), np.log(ys)))


def cumulative_range_mm(e: np.ndarray, s: np.ndarray, a: float) -> np.ndarray:
    inc = st.exact_loglog_range_increments(np.log(e), np.log(s), float(a))
    return np.concatenate([[0.0], np.cumsum(inc)]) * 10.0


def euler_range_mm(
    e: np.ndarray, s: np.ndarray, a: int, e_start: float, *, switch_fraction: float | None = None
) -> float:
    """Range from e_start down to E_CUT with S frozen at each step's entry energy.

    ``switch_fraction``: once the pre-step loss of a step exceeds this fraction of the energy
    (Geant4's 2 % limit, where the mean loss is taken from the range table) the rest of the
    range is added exactly; ``None``: pure forward Euler all the way down.
    """
    ln_e, ln_s = np.log(e), np.log(s)
    ds_cm = STEP_MM / 10.0
    cum = cumulative_range_mm(e, s, a)
    en, n_steps = e_start, 0
    while True:
        sn = math.exp(float(np.interp(math.log(en), ln_e, ln_s)))
        de = sn * ds_cm / a  # MeV/u lost in one full step
        if switch_fraction is not None and de > switch_fraction * en:
            r_rest = float(np.interp(math.log(en), ln_e, cum))  # exact range en -> E_CUT
            return n_steps * STEP_MM + r_rest
        if en - de <= E_CUT:
            return (n_steps + (en - E_CUT) / de) * STEP_MM
        en -= de
        n_steps += 1


def proton_energy_equal_velocity(beta2: np.ndarray) -> np.ndarray:
    gamma = 1.0 / np.sqrt(1.0 - beta2)
    return (gamma - 1.0) * pj.PROTON.mass_mev


def ls_fraction(
    z: int, e_u: np.ndarray, beta2: np.ndarray, q: np.ndarray, l_p: np.ndarray, mat
) -> tuple[np.ndarray, np.ndarray]:
    """(dL_LS / L_p, dL_Barkas / L_p), LS evaluated on a coarse grid and interpolated in ln E."""
    n = int(math.ceil(LS_POINTS_PER_DECADE * math.log10(e_u[-1] / e_u[0]))) + 1
    ec = np.geomspace(e_u[0], e_u[-1], n)
    b2c = np.interp(np.log(ec), np.log(e_u), beta2)
    dl = np.array(
        [
            lsm.delta_l_ls(lsm.ALPHA_FS * z, math.sqrt(b2)).delta_l
            - lsm.bloch_ionmc(b2, 1.0)
            - lsm.mott1_ionmc(math.sqrt(b2), 1.0)
            for b2 in b2c
        ]
    )
    dl_ls = np.interp(np.log(e_u), np.log(ec), dl)
    dl_bark = st.barkas_term(beta2, q, mat) - st.barkas_term(beta2, np.ones_like(q), mat)
    return dl_ls / l_p, dl_bark / l_p


def compute(worktree: Path, cache_dir: Path) -> dict:
    mat = water(78.0)
    icru_path = find_icru(cache_dir)
    icru = load_icru90_water(icru_path)
    e_full = np.unique(
        np.concatenate(
            [
                np.geomspace(
                    E_CUT, 430.0, int(math.ceil(POINTS_PER_DECADE * math.log10(430.0 / E_CUT))) + 1
                ),
                np.array([c[1] for c in CASES]),
            ]
        )
    )
    inputs = {ICRU_ID + " (cache object)": sha256(icru_path)}
    reports = {}
    for ion, mev_u, num in CASES:
        if num is not None:
            p = worktree / REPORT_DIR / f"TD-EM-case{num}.json"
            inputs[f"{REPORT_DIR}/TD-EM-case{num}.json"] = sha256(p)
            crit = [c for c in json.loads(p.read_text())["criteria"] if c["criterion"] == "(a)"][0]
            reports[(ion, mev_u)] = (crit["x"], crit["r"])
    for rel in SOURCES:
        inputs[rel] = sha256(worktree / rel)

    table_src = st.BetheStoppingSource()
    rows = []
    per_ion = {}
    for ion in PROJ:
        p = PROJ[ion]
        # the ICRU 90 alpha array ends at 1000 MeV total = 250 MeV/u
        e_all = e_full[e_full <= 250.0] if ion == "he" else e_full
        s_ion = st.bethe_mass_stopping(e_all, p, mat)  # all options on, ls off: as transported
        r_ion = cumulative_range_mm(e_all, s_ion, p.a)
        tab = table_src.table(mat, p)
        beta2, _, _ = st.kinematics(e_all, p)
        d = {"e": e_all, "s_ion": s_ion, "r_ion": r_ion, "tab": tab}
        if ion == "he":
            s_dir = loglog(e_all * p.a, icru.alpha_energy_mev, icru.alpha_stopping)
            d["variants"] = {"icru90_alpha_direct": s_dir}
        else:
            q = st.pierce_blann_charge(np.sqrt(beta2), p)
            s_p = loglog(
                proton_energy_equal_velocity(beta2), icru.proton_energy_mev, icru.proton_stopping
            )
            l_p = s_p * beta2 / (st.K_MEV_CM2_MOL * mat.z_over_a)
            f_ls, f_bark = ls_fraction(p.z, e_all, beta2, q, l_p, mat)
            base = s_p * q**2
            d["variants"] = {
                "topas_like_plain": base,
                "topas_like_barkas_only": base * (1.0 + f_bark),
                "topas_like_halfls_barkas": base * (1.0 + 0.5 * f_ls + f_bark),
                "topas_like_ls_barkas": base * (1.0 + f_ls + f_bark),
            }
        per_ion[ion] = d
    for ion, mev_u, _num in CASES:
        d = per_ion[ion]
        p = PROJ[ion]
        e_all = d["e"]
        i = int(np.searchsorted(e_all, mev_u))
        r0 = float(d["r_ion"][i])
        row = {
            "ion": ion,
            "mev_u": mev_u,
            "R_ionmc_mm": r0,
            "R_ionmc_table_mm": float(d["tab"].range_at(mev_u)) * 10.0,
        }
        for name, s_var in d["variants"].items():
            r_var = float(cumulative_range_mm(e_all, s_var, p.a)[i])
            row[f"{name}_vs_ionmc_pct"] = (r_var / r0 - 1) * 100
        r_euler = euler_range_mm(e_all, d["s_ion"], p.a, mev_u, switch_fraction=0.02)
        row["prestep_euler_0p1mm_pct"] = (r_euler / r0 - 1) * 100
        r_naive = euler_range_mm(e_all, d["s_ion"], p.a, mev_u)
        row["prestep_pure_euler_0p1mm_pct"] = (r_naive / r0 - 1) * 100
        if (ion, mev_u) in reports:
            x, r = reports[(ion, mev_u)]
            excess = (r / x - 1) * 100
            primary = (
                row["icru90_alpha_direct_vs_ionmc_pct"]
                if ion == "he"
                else row["topas_like_ls_barkas_vs_ionmc_pct"]
            )
            row["topas_r80_over_ionmc_r80_pct"] = excess
            row["construction_diff_pct"] = primary
            row["residual_pct"] = excess - primary - row["prestep_euler_0p1mm_pct"]
            if ion != "he":
                for name in d["variants"]:
                    row[f"residual_{name}_pct"] = (
                        excess - row[f"{name}_vs_ionmc_pct"] - row["prestep_euler_0p1mm_pct"]
                    )
        rows.append(row)
    return {
        "note": (
            "D-3 report-only attribution; aggregates only (ranges in mm of water, ratios in "
            "percent); no NIST/ICRU/Geant4 stopping values are reproduced. C/O "
            "'topas_like_*' are approximations of the Geant4 construction (see script "
            "docstring); pre-step estimate is a planning estimate."
        ),
        "e_cut_mev_per_u": E_CUT,
        "points_per_decade": POINTS_PER_DECADE,
        "ls_points_per_decade": LS_POINTS_PER_DECADE,
        "step_mm": STEP_MM,
        "residual_definition": (
            "topas_r80_over_ionmc_r80_pct - construction_diff_pct - prestep_euler_0p1mm_pct"
        ),
        "cases": rows,
        "input_sha256": inputs,
        "script_sha256": sha256(Path(__file__).resolve()),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worktree", type=Path, default=WORKTREE)
    ap.add_argument("--cache-dir", type=Path, default=WORKTREE / ".ionmc-cache")
    ap.add_argument(
        "--out",
        type=Path,
        default=WORKTREE / "validation/results/stopping/v3-008-d3-range-attribution.json",
    )
    args = ap.parse_args()
    res = compute(args.worktree, args.cache_dir)
    args.out.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
