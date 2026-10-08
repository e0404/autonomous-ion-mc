"""Acceptance row V5 (V3-005): absolute integral depth dose (r = 20 cm) of ionmc against TOPAS and MCsquare.

Usage::

    compare_idd_v5.py --ionmc-dir DIR --topas-runs RUN_DIR... --mcsquare-runs RUN_DIR... --output OUT.json

``DIR`` holds the four partials ``v5-{150,200}-{on,off}.json`` written by
``validation/scripts/transport/steps_v5b.py v5-ionmc`` (per-batch IDD in MeV/(g/cm^2)/primary,
20 batches, 0.5 mm bins, 400 x 400 mm water box). The engine runs are materialised reference runs of
the committed cases ``topas/proton-water-{150,200}mev-idd-r20[-emonly]-seed{1,2,3}`` and
``mcsquare/proton-water-{150,200}mev-idd-r20-{on,off}-seed{1,2,3}`` (each case.json carries a ``v5``
block: energy, nuclear flag). Engine runs are grouped by (engine, energy, nuclear flag); every group
needs >= 3 runs with distinct seeds (read from the manifested case.json and cross-checked against the
native input by ``compare_batches.run_seed``) and 1e5 histories.

Units. TOPAS: dose [Gy] x bin mass [kg] / 1.602176634e-13 J/MeV = MeV in the bin, divided by
histories x rho x dz [g/cm^2] (bin mass from the CSV header widths, rho = 1 g/cm^3).
MCsquare: IDD = sum over x,z voxels of dose [MeV/g/primary] x dx x dz [cm^2] (density cancels), with
dose [MeV/g] = file value x 1e-6 (the eV/g unit of Dose.mhd is INFERRED, see the case.json; the in-grid
total over the beam energy must lie in [0.85, 1.02] or the engine is flagged and excluded).
All curves are cut to the depth range of the ionmc grid (1.1 R(CSDA)) so that "total deposit" is the
same quantity for every code.

Criteria (frozen row V5, equivalence by TOST, alpha = 0.05 per side, i.e. the 90 % interval of the
difference must lie inside the tolerance): plateau (mean over 20-60 mm) within 2 % relative; R80
within 0.5 mm; peak/plateau within 2 % relative; plateau-integrated Delta IDD (on - off; TOPAS
EM-only is the "off" run) within 10 % of the reference Delta; total deposit within 1 %. Estimates:
ionmc = metric of the pooled batch-mean curve with a delete-one jackknife standard error over its 20
batches (df 19); engines = mean of the per-seed metrics, SE = SD/sqrt(n_seeds) (df n-1). The difference
has the Welch-Satterthwaite df (floored) and the one-sided 95 % t quantile; relative differences use
the delta-method SE. TOPAS gates; MCsquare (shared sigma lineage) is report-only (Amendment 7 (b)).

Kappa rule (Amendment 7 (c)), evaluated on TOPAS only, D = (dep_ionmc - dep_TOPAS) / dep_TOPAS:
both total-deposit criteria pass -> kappa not used; a failing energy with D > 0 and D <= 1.5 % ->
build kappa (if every failing energy qualifies); a failing energy with D < 0 or |D| > 1.5 % ->
kappa cannot be the remedy, stop (this also applies when failures are mixed).
The analysis code SHA is git HEAD of this repository with a dirty flag.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

from ionmc.reference import metrics as rm
from ionmc.reference.parsers import ParseError, parse_metaimage, parse_topas_csv
from ionmc.reference.runs import ReferenceRun, RunError, file_hashes, load_run, read_verified

MEV_J = 1.602176634e-13
ENERGIES = (150, 200)
BIN_MM = 0.5
N_BATCHES = 20
HISTORIES = 100_000
MIN_SEEDS = 3
PLATEAU_MM = (20.0, 60.0)
UNIT = "MeV/(g/cm^2)/primary"
DEP_CONTAINED = (0.85, 1.02)  # in-grid total / beam energy; catches a wrong unit
MCSQUARE_UNIT_MEV = 1e-6  # Dose.mhd value -> MeV/g (inferred eV/g)
KAPPA_MAX_D = 0.015
# frozen V5 tolerances: (name, kind, tolerance)
TOLERANCES = {
    "plateau": ("relative", 0.02),
    "r80_mm": ("absolute_mm", 0.5),
    "peak_over_plateau": ("relative", 0.02),
    "delta_idd_plateau_integral": ("relative", 0.10),
    "total_deposit": ("relative", 0.01),
}

# one-sided 95 % Student-t quantiles, df = 1..30; beyond 30 the df=30 value (conservative)
_T95 = (6.3138, 2.9200, 2.3534, 2.1318, 2.0150, 1.9432, 1.8946, 1.8595, 1.8331, 1.8125,
        1.7959, 1.7823, 1.7709, 1.7613, 1.7531, 1.7459, 1.7396, 1.7341, 1.7291, 1.7247,
        1.7207, 1.7171, 1.7139, 1.7109, 1.7081, 1.7056, 1.7033, 1.7011, 1.6991, 1.6973)  # fmt: skip


def t95(df: float) -> float:
    return _T95[min(max(int(math.floor(df)), 1), 30) - 1]


class IddError(ValueError):
    """Raised when inputs are inconsistent (fail closed)."""


# -- estimates ----------------------------------------------------------------------------------
class Est:
    """A value with standard error and degrees of freedom."""

    def __init__(self, value: float, se: float, df: float) -> None:
        self.value, self.se, self.df = float(value), float(se), float(df)

    def as_dict(self) -> dict[str, float]:
        return {"value": self.value, "se": self.se, "df": self.df}


def welch_df(sa: float, dfa: float, sb: float, dfb: float) -> float:
    num = (sa * sa + sb * sb) ** 2
    den = (sa**4 / dfa if sa > 0 else 0.0) + (sb**4 / dfb if sb > 0 else 0.0)
    return num / den if den > 0 else max(dfa, dfb)


def difference(a: Est, b: Est, relative: bool) -> Est:
    """``a - b`` (absolute) or ``(a - b) / b`` (relative, delta method) with Welch df."""
    df = welch_df(a.se, a.df, b.se, b.df)
    if not relative:
        return Est(a.value - b.value, math.hypot(a.se, b.se), df)
    if b.value == 0.0:
        raise IddError("relative difference to a zero reference")
    ratio = a.value / b.value
    return Est(ratio - 1.0, math.hypot(a.se, ratio * b.se) / abs(b.value), df)


def tost(d: Est, tol: float) -> dict[str, Any]:
    """Two one-sided tests: equivalent iff the 90 % interval of ``d`` lies inside ``(-tol, tol)``."""
    q = t95(d.df)
    lo, hi = d.value - q * d.se, d.value + q * d.se
    return {"diff": d.value, "se": d.se, "df": d.df, "t95": q, "ci90": [lo, hi],
            "tolerance": tol, "pass": bool(lo > -tol and hi < tol)}  # fmt: skip


# -- curve metrics ------------------------------------------------------------------------------
def depth_centres(n: int) -> np.ndarray:
    return (np.arange(n) + 0.5) * BIN_MM


def curve_metrics(idd: np.ndarray, rho: float = 1.0) -> dict[str, float]:
    """Plateau mean (20-60 mm), R80, peak/plateau, plateau integral and total deposit of one IDD curve."""
    depth = depth_centres(idd.size)
    sel = (depth >= PLATEAU_MM[0]) & (depth <= PLATEAU_MM[1])
    plateau = float(idd[sel].mean())
    bin_g_cm2 = rho * BIN_MM / 10.0
    return {
        "plateau": plateau,
        "r80_mm": rm.r80(depth, rm.normalize_to_peak(idd)),
        "peak_over_plateau": float(idd.max()) / plateau,
        "plateau_integral_mev": float(idd[sel].sum()) * bin_g_cm2,  # MeV per primary
        "total_deposit": float(idd.sum()) * bin_g_cm2,  # MeV per primary in the grid
    }


METRIC_KEYS = ("plateau", "r80_mm", "peak_over_plateau", "plateau_integral_mev", "total_deposit")


def jackknife_metrics(batches: np.ndarray, rho: float = 1.0) -> dict[str, Est]:
    """Metrics of the pooled batch-mean curve; delete-one jackknife SE over the batches."""
    nb = batches.shape[0]
    full = curve_metrics(batches.mean(axis=0), rho)
    total = batches.sum(axis=0)
    loo = [curve_metrics((total - batches[i]) / (nb - 1), rho) for i in range(nb)]
    out = {}
    for k in METRIC_KEYS:
        vals = np.array([m[k] for m in loo])
        se = math.sqrt((nb - 1) / nb * float(((vals - vals.mean()) ** 2).sum()))
        out[k] = Est(full[k], se, nb - 1)
    return out


def seed_metrics(curves: np.ndarray, rho: float = 1.0) -> dict[str, Est]:
    """Mean of the per-seed metrics; SE = SD / sqrt(n)."""
    n = curves.shape[0]
    if n < 2:
        raise IddError("at least 2 seeds are needed for a standard error")
    per = [curve_metrics(c, rho) for c in curves]
    out = {}
    for k in METRIC_KEYS:
        vals = np.array([m[k] for m in per])
        out[k] = Est(vals.mean(), vals.std(ddof=1) / math.sqrt(n), n - 1)
    return out


def delta_integral(on: dict[str, Est], off: dict[str, Est]) -> Est:
    a, b = on["plateau_integral_mev"], off["plateau_integral_mev"]
    return Est(a.value - b.value, math.hypot(a.se, b.se), welch_df(a.se, a.df, b.se, b.df))


# -- criteria -----------------------------------------------------------------------------------
def criteria(ion_on: dict[str, Est], ion_off: dict[str, Est],
             ref_on: dict[str, Est], ref_off: dict[str, Est]) -> dict[str, Any]:  # fmt: skip
    """The five V5 criteria of ionmc against one reference (nuclear-on for four, on-off for one)."""
    rows: dict[str, Any] = {}
    for name, key in (("plateau", "plateau"), ("r80_mm", "r80_mm"),
                      ("peak_over_plateau", "peak_over_plateau"), ("total_deposit", "total_deposit")):  # fmt: skip
        kind, tol = TOLERANCES[name]
        d = difference(ion_on[key], ref_on[key], relative=kind == "relative")
        rows[name] = {"ionmc": ion_on[key].as_dict(), "reference": ref_on[key].as_dict(),
                      "kind": kind, **tost(d, tol)}  # fmt: skip
    di, dr = delta_integral(ion_on, ion_off), delta_integral(ref_on, ref_off)
    d = difference(di, dr, relative=True)
    name = "delta_idd_plateau_integral"
    rows[name] = {"ionmc": di.as_dict(), "reference": dr.as_dict(), "kind": "relative",
                  "reference_delta_significant": bool(abs(dr.value) > 2.0 * dr.se),
                  **tost(d, TOLERANCES[name][1])}  # fmt: skip
    rows["_all_pass"] = bool(all(rows[k]["pass"] for k in TOLERANCES))
    return rows


def kappa_rule(total_rows: dict[int, dict[str, Any]]) -> dict[str, Any]:
    """Amendment 7 (c) on the TOPAS total-deposit rows per energy (``D`` = relative difference)."""
    d = {e: r["diff"] for e, r in total_rows.items()}
    failing = [e for e, r in total_rows.items() if not r["pass"]]
    if not failing:
        return {"decision": "kappa-not-used", "D": d, "failing_energies": []}
    blocked = [e for e in failing if d[e] < 0.0 or abs(d[e]) > KAPPA_MAX_D]
    decision = "stop-kappa-cannot-be-the-remedy" if blocked else "build-kappa"
    return {"decision": decision, "D": d, "failing_energies": failing, "blocked_energies": blocked}


# -- inputs -------------------------------------------------------------------------------------
def content_digest(doc: dict[str, Any]) -> str:
    """sha256 of the canonical JSON of ``doc`` without ``content_sha256`` (as steps_v5.content_digest)."""
    norm = json.loads(json.dumps(doc))
    norm.pop("content_sha256", None)
    return hashlib.sha256(json.dumps(norm, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_ionmc(path: Path, energy: int, nuclear: bool) -> tuple[np.ndarray, dict[str, Any]]:
    blob = path.read_bytes()
    doc = json.loads(blob)
    if doc.get("content_sha256") != content_digest(doc):
        raise IddError(f"{path.name}: content_sha256 does not match the content")
    expect = {"energy": float(energy), "nuclear": nuclear, "unit": UNIT, "n_batches": N_BATCHES,
              "bin_mm": BIN_MM, "valid": True, "reduced": False}  # fmt: skip
    for k, v in expect.items():
        if doc.get(k) != v:
            raise IddError(f"{path.name}: {k} = {doc.get(k)!r}, expected {v!r}")
    if doc.get("n") != HISTORIES:
        raise IddError(f"{path.name}: n = {doc.get('n')!r}, expected {HISTORIES}")
    batches = np.asarray(doc["idd_batches"], dtype=np.float64)
    if batches.ndim != 2 or batches.shape[0] != N_BATCHES or not np.all(np.isfinite(batches)):
        raise IddError(f"{path.name}: idd_batches must be finite [{N_BATCHES}, nz]")
    return batches, {"file": path.name, "file_sha256": hashlib.sha256(blob).hexdigest(),
                     "content_sha256": doc["content_sha256"], "seed": doc.get("seed"),
                     "density_g_cm3": doc.get("density_g_cm3")}  # fmt: skip


def _batches_module() -> Any:
    path = Path(__file__).resolve().parent / "compare_batches.py"
    spec = importlib.util.spec_from_file_location("compare_batches", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("compare_batches", mod)
    spec.loader.exec_module(mod)
    return mod


def _topas_idd(run: ReferenceRun) -> tuple[np.ndarray, list[str]]:
    rel = "work/idd_dose.csv"
    s = parse_topas_csv(read_verified(run, rel).decode("utf-8"), rel)
    if s.bins[0] != 1 or s.bins[1] != 1 or "Sum" not in s.values:
        raise IddError(f"{run.run_id}: expected a 1x1xN IDD with a Sum column, got {s.bins}")
    if not s.unit.strip().startswith("Gy"):
        raise IddError(f"{run.run_id}: TOPAS unit {s.unit!r} is not Gy")
    to_cm = {"mm": 0.1, "cm": 1.0, "m": 100.0}
    try:
        wx, wy, wz = (w * to_cm[u] for w, u in zip(s.bin_width, s.bin_unit, strict=True))
    except KeyError as exc:
        raise IddError(f"{run.run_id}: unknown bin unit {exc}") from exc
    if abs(wz - BIN_MM / 10.0) > 1e-9 or wx < 40.0 - 1e-9 or wy < 40.0 - 1e-9:
        raise IddError(f"{run.run_id}: bins {wx:g} x {wy:g} x {wz:g} cm are not 0.5 mm deep and >= 400 x 400 mm")
    dose = s.values["Sum"][0, 0, :]
    if not np.all(np.isfinite(dose)):
        raise IddError(f"{run.run_id}: missing bins")
    mass_kg = 1.0 * wx * wy * wz / 1000.0  # 1 g/cm3
    mev = dose * mass_kg / MEV_J
    return mev / (run.histories * 1.0 * wz), [rel]


def _mcsquare_idd(run: ReferenceRun) -> tuple[np.ndarray, list[str]]:
    rel = "work/Outputs/Dose.mhd"

    def loader(name: str) -> bytes:
        if "/" in name or name in ("", ".", ".."):
            raise IddError(f"{run.run_id}: unsafe ElementDataFile {name!r}")
        return read_verified(run, f"work/Outputs/{name}")

    img = parse_metaimage(read_verified(run, rel), rel, loader)
    if img.header.get("ElementDataFile") != "Dose.raw":
        raise IddError(f"{run.run_id}: ElementDataFile must be Dose.raw")
    sx, sy, sz = (s / 10.0 for s in img.spacing)  # cm
    if abs(sy - BIN_MM / 10.0) > 1e-9 or img.dims[0] * sx < 40.0 - 1e-9 or img.dims[2] * sz < 40.0 - 1e-9:
        raise IddError(f"{run.run_id}: grid {img.dims} x {img.spacing} mm is not 0.5 mm deep and >= 400 x 400 mm")
    dose_y = img.data.astype(np.float64).sum(axis=(0, 2)) * MCSQUARE_UNIT_MEV  # MeV/g per y, lateral-summed
    return dose_y[::-1] * sx * sz, [rel, "work/Outputs/Dose.raw"]  # y reversed -> depth; cm^2


_ENGINE_IDD: dict[str, Callable[[ReferenceRun], tuple[np.ndarray, list[str]]]] = {
    "topas": _topas_idd, "mcsquare": _mcsquare_idd}  # fmt: skip


def load_engine_groups(run_dirs: list[Path], engine: str, nz: dict[int, int]
                       ) -> dict[tuple[int, bool], dict[str, Any]]:  # fmt: skip
    """Per (energy, nuclear): curves ``[seeds, nz]`` cut to the ionmc grid, with seeds and hashes."""
    cb = _batches_module()
    groups: dict[tuple[int, bool], dict[str, Any]] = {}
    for rd in run_dirs:
        try:
            run = load_run(rd)
            if run.engine != engine:
                raise IddError(f"{rd}: engine {run.engine!r}, expected {engine!r}")
            v5 = run.case.get("v5")
            if not isinstance(v5, dict) or v5.get("row") != "V5":
                raise IddError(f"{run.run_id}: case.json has no V5 block")
            e, nuc = int(v5["energy_mev"]), bool(v5["nuclear"])
            if e not in nz:
                raise IddError(f"{run.run_id}: energy {e} not in {sorted(nz)}")
            if run.histories != HISTORIES:
                raise IddError(f"{run.run_id}: histories {run.histories} != {HISTORIES}")
            seed = cb.run_seed(run)
            idd, files = _ENGINE_IDD[engine](run)
        except (RunError, ParseError, cb.BatchError) as exc:
            raise IddError(str(exc)) from exc
        if idd.size < nz[e]:
            raise IddError(f"{run.run_id}: {idd.size} depth bins < ionmc grid {nz[e]}")
        cut = idd[: nz[e]]
        total = float(cut.sum()) * BIN_MM / 10.0
        if not DEP_CONTAINED[0] <= total / e <= DEP_CONTAINED[1]:
            raise IddError(
                f"{run.run_id}: in-grid total {total:.3f} MeV is {total / e:.3f} of the beam energy, "
                f"outside {DEP_CONTAINED}: unit or geometry error")  # fmt: skip
        g = groups.setdefault((e, nuc), {"curves": [], "seeds": [], "runs": [], "hashes": {}})
        if seed in g["seeds"]:
            raise IddError(f"{run.run_id}: duplicate seed {seed} in group {(e, nuc)}")
        g["curves"].append(cut)
        g["seeds"].append(seed)
        g["runs"].append(run.run_id)
        g["hashes"][run.run_id] = file_hashes(run, [*files, "inputs/case.json"])
    for key, g in groups.items():
        if len(g["curves"]) < MIN_SEEDS:
            raise IddError(f"{engine} group {key}: {len(g['curves'])} runs, need >= {MIN_SEEDS}")
        g["curves"] = np.stack(g["curves"])
    return groups


# -- driver -------------------------------------------------------------------------------------
def _git(*args: str) -> str:
    repo = Path(__file__).resolve().parent
    return subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True,
                          text=True).stdout.strip()  # fmt: skip


def evaluate(ion: dict[tuple[int, bool], np.ndarray],
             engines: dict[str, dict[tuple[int, bool], np.ndarray]]) -> dict[str, Any]:  # fmt: skip
    """Verdict document body from ionmc batch curves and engine per-seed curves (all cut to one grid)."""
    ion_m = {k: jackknife_metrics(v) for k, v in ion.items()}
    out: dict[str, Any] = {"engines": {}}
    for engine, groups in engines.items():
        gating = engine == "topas"
        per_e: dict[str, Any] = {}
        totals: dict[int, dict[str, Any]] = {}
        for e in ENERGIES:
            ref = {k: seed_metrics(groups[(e, k)]) for k in (True, False)}
            rows = criteria(ion_m[(e, True)], ion_m[(e, False)], ref[True], ref[False])
            per_e[str(e)] = rows
            totals[e] = rows["total_deposit"]
        out["engines"][engine] = {
            "role": "gating" if gating else "report-only",
            "energies": per_e,
            "all_pass": bool(all(per_e[str(e)]["_all_pass"] for e in ENERGIES)),
        }
        if gating:
            out["kappa_rule"] = kappa_rule(totals)
    out["gating_engine"] = "topas"
    out["pass"] = out["engines"]["topas"]["all_pass"]
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ionmc-dir", required=True, type=Path)
    ap.add_argument("--topas-runs", nargs="+", required=True, type=Path)
    ap.add_argument("--mcsquare-runs", nargs="+", type=Path, default=[])
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args(argv)

    ion: dict[tuple[int, bool], np.ndarray] = {}
    ion_src: dict[str, Any] = {}
    for e in ENERGIES:
        for nuc in (True, False):
            name = f"v5-{e}-{'on' if nuc else 'off'}.json"
            ion[(e, nuc)], ion_src[name] = load_ionmc(args.ionmc_dir / name, e, nuc)
    nz = {e: ion[(e, True)].shape[1] for e in ENERGIES}
    if any(ion[(e, False)].shape[1] != nz[e] for e in ENERGIES):
        raise IddError("ionmc on/off grids differ")
    sources = {"topas": args.topas_runs, "mcsquare": args.mcsquare_runs}
    groups = {eng: load_engine_groups(dirs, eng, nz) for eng, dirs in sources.items() if dirs}
    missing = [(eng, e, n) for eng, g in groups.items() for e in ENERGIES for n in (True, False)
               if (e, n) not in g]  # fmt: skip
    if missing:
        raise IddError(f"missing engine groups (engine, energy, nuclear): {missing}")
    doc = evaluate(ion, {eng: {k: g["curves"] for k, g in gr.items()} for eng, gr in groups.items()})
    code_sha, dirty = _git("rev-parse", "HEAD"), bool(_git("status", "--porcelain"))
    doc.update({
        "row": "V5", "unit": UNIT, "analysis_code_sha": code_sha, "analysis_code_dirty": dirty,
        "tolerances": {k: {"kind": v[0], "tolerance": v[1]} for k, v in TOLERANCES.items()},
        "ionmc_inputs": ion_src,
        "engine_inputs": {eng: {f"{e}-{'on' if n else 'off'}":
                                {"seeds": g["seeds"], "runs": g["runs"], "output_sha256": g["hashes"]}
                                for (e, n), g in gr.items()} for eng, gr in groups.items()},
    })  # fmt: skip
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=2) + "\n")
    for eng, ev in doc["engines"].items():
        for e, rows in ev["energies"].items():
            line = "  ".join(f"{k}:{'ok' if rows[k]['pass'] else 'FAIL'}" for k in TOLERANCES)
            print(f"{eng:9s} {e} MeV  {line}")
    print(f"kappa rule: {doc['kappa_rule']['decision']}; V5 (TOPAS gating): {'PASS' if doc['pass'] else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
