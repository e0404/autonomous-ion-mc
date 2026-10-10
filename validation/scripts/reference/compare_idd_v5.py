"""Acceptance row V5 (V3-005): absolute integral depth dose (r = 20 cm) of ionmc against TOPAS and MCsquare.

Usage::

    compare_idd_v5.py --ionmc-dir DIR --topas-runs RUN_DIR... --mcsquare-runs RUN_DIR... --output OUT.json
                      [--cases-dir validation/reference_cases]

``DIR`` holds the four partials ``v5-{150,200}-{on,off}.json`` written by
``validation/scripts/transport/steps_v5b.py v5-ionmc`` (per-batch IDD in MeV/(g/cm^2)/primary,
20 batches, 0.5 mm bins, 400 x 400 mm water box). The engine runs are materialised reference runs of
the committed cases ``topas/proton-water-{150,200}mev-idd-r20[-emonly]-seed{1,2,3}`` and
``mcsquare/proton-water-{150,200}mev-idd-r20-{on,off}-seed{1,2,3}`` (each case.json carries a ``v5``
block, a LABEL only). Lineage (C19 F2): every engine run passes ``compare_batches.run_fingerprint``
(clean-commit request.json, manifested inputs, histories == native input == engine summary); the
beam energy and the configuration (TOPAS full vs EM-only from the ``Ph/Default/Modules`` line,
MCsquare nuclear on/off from ``Simulate_Nuclear_Interactions``) are DERIVED from the verified native
input (MCsquare energy: the ``####Energy (MeV)`` entry of the manifested Plan.txt), the case.json
``v5`` label must agree or the run is refused. Engine runs are grouped by (engine, derived energy,
derived nuclear flag); the replicates of a group must have the same fingerprint configuration and
engine identity, and every group needs >= 3 runs with distinct seeds (read from the manifested
case.json and cross-checked against the native input by ``compare_batches.run_seed``) and 1e5
histories. The fingerprints are recorded in the verdict.

Binding to the frozen cases (C20 G1, review 0f1aa5d9): the fingerprint only compares replicates with
each other, so every run must also be BYTE-IDENTICAL (all of ``inputs/``, case.json included, equal
file-name sets) to exactly one committed case directory ``--cases-dir/<engine>/<name>``; that case must
carry the V5 block, may be bound by one run only, and supplies the (energy, nuclear) group. The native
input is cross-checked against the frozen physics (TOPAS: Modules exactly the full QGSP_BIC_HP set or
exactly ``g4em-standard_opt4``, BeamEnergy integral and equal, BeamEnergySpread 0, a
CutForAllParticles line; MCsquare: the four Simulate_* flags all True or all False, Plan.txt energy
integral and equal). The bound case names of a group must be the seed variants ``<family>-seed<k>`` of
one family; case, file hashes and family are recorded in the verdict. The CLI exits 1 when the gating verdict
(TOPAS) fails, after writing the verdict JSON.

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

V3-005C (C7) adds, without changing any of the above, the TOPAS modes ``emelastic`` (Modules exactly
``g4em-standard_opt4`` + ``g4h-elastic_HP``; cases ``...-r20-emelastic-seed{1,2,3}``, row V11, plan rows
10/11) and ``noelastic`` (the six frozen modules minus ``g4h-elastic_HP``; ``...-r20-noelastic-seed{1,2,3}``,
``v5.row`` X-elastic-factor, rows 12/13, report-only) with :func:`build_v11_verdict` (row V11), the rung-2
test :func:`f_ne_test` and :func:`elastic_factor_report`; see the V11 section below and the plan quotes
in ``V11_PLAN_RULE``.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import re
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
DEFAULT_CASES_DIR = Path(__file__).resolve().parents[3] / "validation" / "reference_cases"
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


TOPAS_MODULES_FULL = frozenset({"g4em-standard_opt4", "g4h-phy_QGSP_BIC_HP", "g4h-elastic_HP",
                                "g4stopping", "g4ion-binarycascade", "g4decay"})
TOPAS_MODULES_EMONLY = frozenset({"g4em-standard_opt4"})
MC_NUCLEAR_FLAGS = ("Simulate_Nuclear_Interactions", "Simulate_Secondary_Protons",
                    "Simulate_Secondary_Deuterons", "Simulate_Secondary_Alphas")
_TOPAS_MODULES = re.compile(r"^[ \t]*sv:Ph/Default/Modules[ \t]*=[ \t]*(\d+)((?:[ \t]+\"[^\"\n]*\")*)[ \t]*$", re.M)
_TOPAS_SPREAD = re.compile(r"^[ \t]*[ud]:So/Beam/BeamEnergySpread[ \t]*=[ \t]*([0-9.eE+-]+)[ \t]*(?:%)?[ \t]*$", re.M)
_TOPAS_CUT = re.compile(r"^[ \t]*d:Ph/Default/CutForAllParticles[ \t]*=[ \t]*[0-9.eE+-]+[ \t]*\w+[ \t]*$", re.M)
_SEED_FAMILY = re.compile(r"^(?P<family>.+)-seed(?P<k>\d+)$")
_TOPAS_ENERGY = re.compile(r"^[ \t]*d:So/Beam/BeamEnergy[ \t]*=[ \t]*([0-9.eE+-]+)[ \t]*MeV[ \t]*$", re.M)
_MC_PLAN_ENERGY = re.compile(r"^####Energy \(MeV\)[ \t]*\n[ \t]*([0-9.eE+-]+)[ \t]*$", re.M)


def _mc_flag(native: str, name: str, run_id: str) -> bool:
    hits = re.findall(rf"^[ \t]*{name}[ \t]+(True|False)[ \t]*$", native, re.M)
    if len(hits) != 1:
        raise IddError(f"{run_id}: expected exactly one {name} line, found {len(hits)}")
    return hits[0] == "True"


def derive_config(run: ReferenceRun, fp: dict[str, Any]) -> tuple[int, bool]:
    """(energy in MeV, nuclear flag) from the verified native input of the run (fail closed).

    The frozen physics is enforced: TOPAS Modules exactly the full set (nuclear on) or exactly
    ``g4em-standard_opt4`` (EM-only), BeamEnergySpread 0, CutForAllParticles present; MCsquare the four
    Simulate_* flags all True or all False."""
    native = read_verified(run, f"inputs/{run.case['input']}").decode("utf-8", errors="replace")
    if run.engine == "topas":
        hits = _TOPAS_ENERGY.findall(native)
        mods = _TOPAS_MODULES.findall(native)
        if len(mods) != 1:
            raise IddError(f"{run.run_id}: expected exactly one parsable Ph/Default/Modules line")
        names = re.findall(r'"([^"\n]*)"', mods[0][1])
        if int(mods[0][0]) != len(names) or len(set(names)) != len(names):
            raise IddError(f"{run.run_id}: Modules count/duplicates inconsistent: {names}")
        if set(names) == TOPAS_MODULES_FULL:
            nuclear = True
        elif set(names) == TOPAS_MODULES_EMONLY:
            nuclear = False
        else:
            raise IddError(
                f"{run.run_id}: TOPAS Modules {sorted(names)} are neither the frozen full set "
                f"{sorted(TOPAS_MODULES_FULL)} nor {sorted(TOPAS_MODULES_EMONLY)}")  # fmt: skip
        if nuclear != (fp["group"] != "topas-emonly"):
            raise IddError(f"{run.run_id}: fingerprint group {fp['group']!r} contradicts the Modules line")
        spread = _TOPAS_SPREAD.findall(native)
        if len(spread) != 1 or float(spread[0]) != 0.0:
            raise IddError(f"{run.run_id}: BeamEnergySpread must be exactly one line equal to 0: {spread}")
        if len(_TOPAS_CUT.findall(native)) != 1:
            raise IddError(f"{run.run_id}: expected exactly one Ph/Default/CutForAllParticles line")
    else:
        hits = _MC_PLAN_ENERGY.findall(read_verified(run, "inputs/Plan.txt").decode("utf-8"))
        flags = {n: _mc_flag(native, n, run.run_id) for n in MC_NUCLEAR_FLAGS}
        if len(set(flags.values())) != 1:
            raise IddError(f"{run.run_id}: nuclear/secondary flags are not all equal: {flags}")
        nuclear = flags[MC_NUCLEAR_FLAGS[0]]
    if len(hits) != 1 or float(hits[0]) != int(float(hits[0])):
        raise IddError(f"{run.run_id}: beam energy not unique/integral in the native input: {hits}")
    return int(float(hits[0])), nuclear


# -- binding to the frozen committed cases (C20 G1) -----------------------------------------------
def _committed_cases(cases_dir: Path, engine: str) -> dict[str, dict[str, bytes]]:
    """``{case name: {relative file name: bytes}}`` of every committed case of ``engine``."""
    root = cases_dir / engine
    if not root.is_dir():
        raise IddError(f"no committed cases directory {root}")
    out: dict[str, dict[str, bytes]] = {}
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        out[d.name] = {f.relative_to(d).as_posix(): f.read_bytes()
                       for f in sorted(d.rglob("*")) if f.is_file()}  # fmt: skip
    return out


def bind_run(run: ReferenceRun, committed: dict[str, dict[str, bytes]]) -> tuple[str, dict[str, str]]:
    """The one committed case whose files equal ``inputs/`` of the run byte for byte (names included)."""
    inputs = {r[len("inputs/"):]: read_verified(run, r)
              for r in sorted(run.files) if r.startswith("inputs/")}  # fmt: skip
    hits = [n for n, files in committed.items() if files == inputs]
    if len(hits) != 1:
        raise IddError(
            f"{run.run_id}: inputs/ ({sorted(inputs)}) match {len(hits)} committed {run.engine} cases "
            f"{hits}, exactly one is required (run not frozen in the repository?)")  # fmt: skip
    return hits[0], {f: hashlib.sha256(b).hexdigest() for f, b in sorted(committed[hits[0]].items())}


def load_engine_groups(run_dirs: list[Path], engine: str, nz: dict[int, int], cases_dir: Path
                       ) -> dict[tuple[int, bool], dict[str, Any]]:  # fmt: skip
    """Per (energy, nuclear): curves ``[seeds, nz]`` cut to the ionmc grid, with seeds and hashes.
    Every run is bound to its frozen committed case under ``cases_dir/<engine>``."""

    def classify(run: ReferenceRun, case_name: str, case_doc: dict[str, Any], cb: Any
                 ) -> tuple[tuple[int, bool], int, dict[str, Any]]:  # fmt: skip
        v5 = case_doc.get("v5")
        if not isinstance(v5, dict) or v5.get("row") != "V5":
            raise IddError(f"{run.run_id}: committed case {case_name} has no V5 block")
        fp = cb.run_fingerprint(run)  # fail closed: lineage of the run
        e, nuc = derive_config(run, fp)
        if (int(v5["energy_mev"]), bool(v5["nuclear"])) != (e, nuc):
            raise IddError(
                f"{run.run_id}: committed case {case_name} v5 label ({v5['energy_mev']}, "
                f"{v5['nuclear']}) != native input ({e}, {nuc}): mislabeled run")  # fmt: skip
        return (e, nuc), e, fp

    return _load_groups(run_dirs, engine, nz, cases_dir, classify)  # type: ignore[return-value]


def _load_groups(run_dirs: list[Path], engine: str, nz: dict[int, int], cases_dir: Path,
                 classify: Callable[..., Any]) -> dict[Any, dict[str, Any]]:  # fmt: skip
    """Shared loader: ``classify(run, case_name, case_doc, cb)`` returns ``(group key, energy,
    fingerprint)`` after the label and native-input checks of its mode (fail closed)."""
    cb = _batches_module()
    committed = _committed_cases(cases_dir, engine)
    bound: dict[str, str] = {}
    groups: dict[Any, dict[str, Any]] = {}
    for rd in run_dirs:
        try:
            run = load_run(rd)
            if run.engine != engine:
                raise IddError(f"{rd}: engine {run.engine!r}, expected {engine!r}")
            case_name, case_hashes = bind_run(run, committed)  # byte-identical to a frozen case
            if case_name in bound:
                raise IddError(f"{run.run_id}: committed case {case_name} is already bound by "
                               f"{bound[case_name]}")  # fmt: skip
            bound[case_name] = run.run_id
            key, e, fp = classify(run, case_name, json.loads(committed[case_name]["case.json"]), cb)
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
        g = groups.setdefault(
            key, {"curves": [], "seeds": [], "runs": [], "hashes": {}, "fingerprints": {},
                       "bound_cases": {}})  # fmt: skip
        if seed in g["seeds"]:
            raise IddError(f"{run.run_id}: duplicate seed {seed} in group {key}")
        g["curves"].append(cut)
        g["seeds"].append(seed)
        g["runs"].append(run.run_id)
        g["hashes"][run.run_id] = file_hashes(run, [*files, "inputs/case.json"])
        g["fingerprints"][run.run_id] = {k: fp[k] for k in ("group", "config_sha256", "identity_sha256")}
        g["bound_cases"][run.run_id] = {"bound_case": f"{engine}/{case_name}",
                                        "committed_sha256": case_hashes}  # fmt: skip
    for key, g in groups.items():
        if len(g["curves"]) < MIN_SEEDS:
            raise IddError(f"{engine} group {key}: {len(g['curves'])} runs, need >= {MIN_SEEDS}")
        for field in ("config_sha256", "identity_sha256"):
            if len({f[field] for f in g["fingerprints"].values()}) != 1:
                raise IddError(f"{engine} group {key}: replicates differ in {field} "
                               f"(configuration or engine identity): {g['fingerprints']}")  # fmt: skip
        names = [c["bound_case"].split("/", 1)[1] for c in g["bound_cases"].values()]
        parsed = [_SEED_FAMILY.match(n) for n in names]
        if not all(parsed):
            raise IddError(f"{engine} group {key}: bound cases {names} are not '<family>-seed<k>' variants")
        families = {m.group("family") for m in parsed if m}
        if len(families) != 1:
            raise IddError(f"{engine} group {key}: bound cases belong to several families {sorted(families)}")
        g["family"] = families.pop()
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


def build_verdict(ionmc_dir: Path, topas_runs: list[Path], mcsquare_runs: list[Path],
                  cases_dir: Path | None = None) -> dict[str, Any]:  # fmt: skip
    """The V5 verdict document (``pass`` = the gating TOPAS verdict) from the four ionmc partials
    of ``ionmc_dir`` and the engine run directories; raises :class:`IddError` (fail closed)."""
    ion: dict[tuple[int, bool], np.ndarray] = {}
    ion_src: dict[str, Any] = {}
    for e in ENERGIES:
        for nuc in (True, False):
            name = f"v5-{e}-{'on' if nuc else 'off'}.json"
            ion[(e, nuc)], ion_src[name] = load_ionmc(ionmc_dir / name, e, nuc)
    nz = {e: ion[(e, True)].shape[1] for e in ENERGIES}
    if any(ion[(e, False)].shape[1] != nz[e] for e in ENERGIES):
        raise IddError("ionmc on/off grids differ")
    sources = {"topas": topas_runs, "mcsquare": mcsquare_runs}
    absent = [eng for eng, dirs in sources.items() if not dirs]
    if absent:
        raise IddError(f"no {absent[0]} reference runs: V5 needs both engines (missing: {absent})")
    cases = DEFAULT_CASES_DIR if cases_dir is None else Path(cases_dir)
    groups = {eng: load_engine_groups(dirs, eng, nz, cases) for eng, dirs in sources.items()}
    missing = [(eng, e, n) for eng, g in groups.items() for e in ENERGIES for n in (True, False)
               if (e, n) not in g]  # fmt: skip
    if missing:
        raise IddError(f"missing engine groups (engine, energy, nuclear): {missing}")
    doc = evaluate(ion, {eng: {k: g["curves"] for k, g in gr.items()} for eng, gr in groups.items()})
    code_sha, dirty = _git("rev-parse", "HEAD"), bool(_git("status", "--porcelain"))
    doc.update({
        "row": "V5", "unit": UNIT, "analysis_code_sha": code_sha, "analysis_code_dirty": dirty,
        "tolerances": {k: {"kind": v[0], "tolerance": v[1]} for k, v in TOLERANCES.items()},
        "ionmc_inputs": ion_src, "cases_dir": str(cases),
        "engine_inputs": {eng: {f"{e}-{'on' if n else 'off'}":
                                {"seeds": g["seeds"], "runs": g["runs"], "output_sha256": g["hashes"],
                                 "fingerprints": g["fingerprints"], "family": g["family"],
                                 "binding": g["bound_cases"]}
                                for (e, n), g in gr.items()} for eng, gr in groups.items()},
    })  # fmt: skip
    return doc


# -- V11 (plan Amendment 14 (c) row V11, (f) rung 2 and rows 10-13) ------------------------------------
V11_TOL_F = 0.02  # |F_ionmc / F_TOPAS - 1|, frozen (plan row V11); never a CLI option
V11_TOL_GAIN = 0.005  # absolute, plateau gain; frozen
V11_NE_TOL = 0.02  # rung-2 threshold on F_ne,ionmc / F_ne,TOPAS - 1
V11_WINDOW_MM = 10.0  # depth window before R80 of the F_ne depth-resolved report
V11_PLAN_RULE = {
    "source": "validation/plans/v3-005-acceptance.md, Amendment 14 (c) row V11 (l.204), (f) (l.253-277)",
    "V11": ("Elastic-only effect, ionmc vs TOPAS, at 150 and 200 MeV. F = (peak/plateau)_EM+el / "
            "(peak/plateau)_EM | |F_ionmc/F_TOPAS - 1| within 0.02, and the plateau gain "
            "(plateau_EM+el/plateau_EM - 1) within 0.005 absolute. Both use the V5 TOST on the 90 % "
            "interval with Welch df (Amendment 9)"),
    "F_ne": ("F_ne = (peak/plateau)_nuclear on / (peak/plateau)_elastic only, ionmc (V5 \"on\" vs V11) "
             "against TOPAS (full vs rows 10/11). Attributed iff the 90 % interval of "
             "F_ne,ionmc/F_ne,TOPAS - 1 lies beyond 0.02 on the side of the residual. The "
             "depth-resolved difference in the last 10 mm before R80 is reported."),
    "rows_12_13": "Rows 12/13 give the full/X1 elastic factor, report-only.",
}  # fmt: skip
MODES = ("full", "emonly", "emelastic", "noelastic")
TOPAS_MODULES_ELASTIC_ONLY = frozenset({"g4em-standard_opt4", "g4h-elastic_HP"})
TOPAS_MODULES_NOELASTIC = TOPAS_MODULES_FULL - {"g4h-elastic_HP"}
MODE_MODULES = {"full": TOPAS_MODULES_FULL, "emonly": TOPAS_MODULES_EMONLY,
                "emelastic": TOPAS_MODULES_ELASTIC_ONLY, "noelastic": TOPAS_MODULES_NOELASTIC}  # fmt: skip
MODE_FP_GROUP = {"full": "topas", "emonly": "topas-emonly", "emelastic": "topas", "noelastic": "topas"}
"""``compare_batches.run_fingerprint`` labels every Modules line with a ``g4h-``/``g4ion`` module
``topas``: emelastic and noelastic share that label with full, and ``config_sha256`` plus the exact
module-set check separate them."""
MODE_ROW = {"full": "V5", "emonly": "V5", "emelastic": "V11", "noelastic": "X-elastic-factor"}
MODE_NUCLEAR = {"full": True, "emonly": False, "emelastic": False, "noelastic": True}


def _topas_module_names(run: ReferenceRun, native: str) -> list[str]:
    mods = _TOPAS_MODULES.findall(native)
    if len(mods) != 1:
        raise IddError(f"{run.run_id}: expected exactly one parsable Ph/Default/Modules line")
    names = re.findall(r'"([^"\n]*)"', mods[0][1])
    if int(mods[0][0]) != len(names) or len(set(names)) != len(names):
        raise IddError(f"{run.run_id}: Modules count/duplicates inconsistent: {names}")
    return names


def derive_mode(run: ReferenceRun, fp: dict[str, Any]) -> tuple[int, str]:
    """(energy in MeV, mode) of a TOPAS run from its verified native input: the Modules set must be
    EXACTLY one of ``MODE_MODULES`` (full, emonly, emelastic, noelastic), the fingerprint group the
    one of the mode, BeamEnergySpread 0, one CutForAllParticles line, an integral unique energy."""
    if run.engine != "topas":
        raise IddError(f"{run.run_id}: V11 modes are TOPAS only, engine {run.engine!r}")
    native = read_verified(run, f"inputs/{run.case['input']}").decode("utf-8", errors="replace")
    names = _topas_module_names(run, native)
    modes = [m for m in MODES if set(names) == MODE_MODULES[m]]
    if len(modes) != 1:
        raise IddError(f"{run.run_id}: TOPAS Modules {sorted(names)} are none of the frozen sets "
                       f"{ {m: sorted(v) for m, v in MODE_MODULES.items()} }")  # fmt: skip
    mode = modes[0]
    if fp["group"] != MODE_FP_GROUP[mode]:
        raise IddError(f"{run.run_id}: fingerprint group {fp['group']!r} contradicts the {mode} Modules line")
    spread = _TOPAS_SPREAD.findall(native)
    if len(spread) != 1 or float(spread[0]) != 0.0:
        raise IddError(f"{run.run_id}: BeamEnergySpread must be exactly one line equal to 0: {spread}")
    if len(_TOPAS_CUT.findall(native)) != 1:
        raise IddError(f"{run.run_id}: expected exactly one Ph/Default/CutForAllParticles line")
    hits = _TOPAS_ENERGY.findall(native)
    if len(hits) != 1 or float(hits[0]) != int(float(hits[0])):
        raise IddError(f"{run.run_id}: beam energy not unique/integral in the native input: {hits}")
    return int(float(hits[0])), mode


def load_topas_mode_groups(run_dirs: list[Path], nz: dict[int, int], cases_dir: Path,
                           modes: tuple[str, ...] = MODES) -> dict[tuple[int, str], dict[str, Any]]:  # fmt: skip
    """TOPAS groups keyed ``(energy, mode)`` for ``modes`` with the V5 lineage and binding checks; a run
    of another mode, a label (row, mode, energy, nuclear flag, modules) that disagrees with the verified
    native input, or an unbound run is refused (fail closed)."""

    def classify(run: ReferenceRun, case_name: str, case_doc: dict[str, Any], cb: Any
                 ) -> tuple[tuple[int, str], int, dict[str, Any]]:  # fmt: skip
        v5 = case_doc.get("v5")
        if not isinstance(v5, dict):
            raise IddError(f"{run.run_id}: committed case {case_name} has no v5 block")
        fp = cb.run_fingerprint(run)
        e, mode = derive_mode(run, fp)
        if mode not in modes:
            raise IddError(f"{run.run_id}: derived mode {mode!r} is not accepted here {modes}")
        label = (v5.get("row"), v5.get("mode", "full" if v5.get("nuclear") else "emonly"),
                 int(v5["energy_mev"]), bool(v5["nuclear"]))  # fmt: skip
        if label != (MODE_ROW[mode], mode, e, MODE_NUCLEAR[mode]):
            raise IddError(
                f"{run.run_id}: committed case {case_name} v5 label {label} != native input "
                f"{(MODE_ROW[mode], mode, e, MODE_NUCLEAR[mode])}: mislabeled run")  # fmt: skip
        if "modules" in v5 and set(v5["modules"]) != MODE_MODULES[mode]:
            raise IddError(f"{run.run_id}: committed case {case_name} lists modules {v5['modules']}")
        return (e, mode), e, fp

    return _load_groups(run_dirs, "topas", nz, cases_dir, classify)


def ratio_est(a: Est, b: Est) -> Est:
    """``a / b`` for independent estimates (delta method, Welch df)."""
    if b.value == 0.0:
        raise IddError("ratio to a zero reference")
    r = a.value / b.value
    return Est(r, math.hypot(a.se, r * b.se) / abs(b.value), welch_df(a.se, a.df, b.se, b.df))


def paired_ratio_jackknife(num: np.ndarray, den: np.ndarray) -> dict[str, Est]:
    """``peak_over_plateau`` and ``plateau`` ratios (num / den) of the pooled batch-mean curves of two
    common-random-number runs, with the delete-one jackknife that removes the SAME batch from both
    (df = batches - 1); keys ``pp_ratio`` and ``plateau_ratio``."""
    if num.shape != den.shape or num.ndim != 2 or num.shape[0] < 2:
        raise IddError("paired jackknife needs two equal [batches >= 2, nz] arrays")
    nb = num.shape[0]

    def ratios(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
        ma, mb = curve_metrics(a), curve_metrics(b)
        return {"pp_ratio": ma["peak_over_plateau"] / mb["peak_over_plateau"],
                "plateau_ratio": ma["plateau"] / mb["plateau"]}  # fmt: skip

    full = ratios(num.mean(axis=0), den.mean(axis=0))
    tn, td = num.sum(axis=0), den.sum(axis=0)
    loo = [ratios((tn - num[i]) / (nb - 1), (td - den[i]) / (nb - 1)) for i in range(nb)]
    out = {}
    for k, v in full.items():
        vals = np.array([r[k] for r in loo])
        out[k] = Est(v, math.sqrt((nb - 1) / nb * float(((vals - vals.mean()) ** 2).sum())), nb - 1)
    return out


def _ci90(d: Est) -> dict[str, Any]:
    q = t95(d.df)
    return {"value": d.value, "se": d.se, "df": d.df, "t95": q,
            "ci90": [d.value - q * d.se, d.value + q * d.se]}  # fmt: skip


def v11_energy(ion_emel: np.ndarray, ion_emonly: np.ndarray, top_emel: np.ndarray,
               top_emonly: np.ndarray) -> dict[str, Any]:  # fmt: skip
    """Row V11 at one energy: F and the plateau gain, ionmc (paired jackknife) vs TOPAS (seed means)."""
    ion = paired_ratio_jackknife(ion_emel, ion_emonly)
    te, tn = seed_metrics(top_emel), seed_metrics(top_emonly)
    f_top = ratio_est(te["peak_over_plateau"], tn["peak_over_plateau"])
    g_top_ratio = ratio_est(te["plateau"], tn["plateau"])
    g_top = Est(g_top_ratio.value - 1.0, g_top_ratio.se, g_top_ratio.df)
    g_ion = Est(ion["plateau_ratio"].value - 1.0, ion["plateau_ratio"].se, ion["plateau_ratio"].df)
    f_row = {"ionmc": ion["pp_ratio"].as_dict(), "reference": f_top.as_dict(), "kind": "relative",
             **tost(difference(ion["pp_ratio"], f_top, relative=True), V11_TOL_F)}  # fmt: skip
    g_row = {"ionmc": g_ion.as_dict(), "reference": g_top.as_dict(), "kind": "absolute",
             **tost(difference(g_ion, g_top, relative=False), V11_TOL_GAIN)}  # fmt: skip
    return {"F": f_row, "plateau_gain": g_row, "pass": bool(f_row["pass"] and g_row["pass"]),
            "F_diff_sign": "ionmc above TOPAS" if f_row["diff"] > 0 else "ionmc below TOPAS"}  # fmt: skip


def f_ne_test(ion_on: np.ndarray, ion_emel: np.ndarray, top_full: np.ndarray,
              top_emel: np.ndarray) -> dict[str, Any]:  # fmt: skip
    """Rung-2 non-elastic factor test (plan (f)): F_ne = (peak/plateau)_on / (peak/plateau)_elastic-only,
    ionmc (V5 "on" vs V11 emelastic; independent samples) against TOPAS (full vs row 10/11). The 90 %
    interval of ``F_ne,ionmc / F_ne,TOPAS - 1`` is attributed iff it lies beyond 0.02 on the side of the
    residual (both sides reported; the residual side belongs to the V5 verdict). Also reports the
    depth-resolved relative differences in the last 10 mm before the TOPAS-full R80. Never gates."""
    on_m, el_m = jackknife_metrics(ion_on), jackknife_metrics(ion_emel)
    f_ion = ratio_est(on_m["peak_over_plateau"], el_m["peak_over_plateau"])
    tf, te = seed_metrics(top_full), seed_metrics(top_emel)
    f_top = ratio_est(tf["peak_over_plateau"], te["peak_over_plateau"])
    d = difference(f_ion, f_top, relative=True)
    ci = _ci90(d)
    lo, hi = ci["ci90"]
    top_full_mean, top_el_mean = top_full.mean(axis=0), top_emel.mean(axis=0)
    depth = depth_centres(top_full_mean.size)
    r80 = rm.r80(depth, rm.normalize_to_peak(top_full_mean))
    sel = (depth >= r80 - V11_WINDOW_MM) & (depth <= r80)
    on_mean, el_mean = ion_on.mean(axis=0), ion_emel.mean(axis=0)
    idd_diff = (on_mean[sel] / top_full_mean[sel] - 1.0)
    ratio_diff = (on_mean[sel] / el_mean[sel]) / (top_full_mean[sel] / top_el_mean[sel]) - 1.0
    return {
        "F_ne_ionmc": f_ion.as_dict(), "F_ne_TOPAS": f_top.as_dict(), "ratio_minus_1": ci,
        "threshold": V11_NE_TOL,
        "attributed_if_residual_above": bool(lo > V11_NE_TOL),
        "attributed_if_residual_below": bool(hi < -V11_NE_TOL),
        "window_mm": [float(r80 - V11_WINDOW_MM), float(r80)], "r80_topas_full_mm": float(r80),
        "window_depth_mm": depth[sel].tolist(),
        "window_idd_on_over_topas_full_minus_1": idd_diff.tolist(),
        "window_idd_on_over_topas_full_mean": float(idd_diff.mean()),
        "window_Fne_depth_ratio_minus_1": ratio_diff.tolist(),
        "window_Fne_depth_ratio_mean": float(ratio_diff.mean()),
        "role": "report-only for V11; read with the V5 residual side (rung 2)",
    }  # fmt: skip


def elastic_factor_report(top_full: np.ndarray, top_noel: np.ndarray) -> dict[str, Any]:
    """Rows 12/13, report-only: TOPAS full / no-hadronic-elastic factors of peak/plateau and plateau."""
    a, b = seed_metrics(top_full), seed_metrics(top_noel)
    return {"F_el_pp": _ci90(ratio_est(a["peak_over_plateau"], b["peak_over_plateau"])),
            "plateau_ratio": _ci90(ratio_est(a["plateau"], b["plateau"])), "role": "report-only"}  # fmt: skip


def load_ionmc_v11(path: Path, energy: int, mode: str) -> tuple[np.ndarray, dict[str, Any]]:
    """A ``v11-{energy}-{emel|emonly}`` partial: ``load_ionmc`` checks plus the V11 labels."""
    nuclear = mode == "emel"
    batches, src = load_ionmc(path, energy, nuclear)
    doc = json.loads(path.read_bytes())
    expect = {"row": f"v11-{energy}-{mode}", "mode": mode, "elastic": True, "elastic_only": nuclear}
    for k, v in expect.items():
        if doc.get(k) != v:
            raise IddError(f"{path.name}: {k} = {doc.get(k)!r}, expected {v!r}")
    return batches, src


def build_v11_verdict(ionmc_dir: Path, topas_runs: list[Path],
                      cases_dir: Path | None = None) -> dict[str, Any]:  # fmt: skip
    """Row V11 verdict (``pass``: both TOSTs at both energies) from the four ``v11-*`` partials of
    ``ionmc_dir`` and the TOPAS run directories (emelastic and emonly required; full and noelastic
    optional). F_ne is computed iff ``v5-{e}-on.json`` partials exist in ``ionmc_dir`` and the full
    group is present; the elastic factor (rows 12/13) iff full and noelastic are present. Both are
    report-only. Raises :class:`IddError`."""
    ion: dict[tuple[int, str], np.ndarray] = {}
    ion_src: dict[str, Any] = {}
    for e in ENERGIES:
        for mode in ("emel", "emonly"):
            name = f"v11-{e}-{mode}.json"
            ion[(e, mode)], ion_src[name] = load_ionmc_v11(ionmc_dir / name, e, mode)
    nz = {e: ion[(e, "emel")].shape[1] for e in ENERGIES}
    if any(ion[(e, "emonly")].shape[1] != nz[e] for e in ENERGIES):
        raise IddError("ionmc emel/emonly grids differ")
    for e in ENERGIES:
        if ion_src[f"v11-{e}-emel.json"]["seed"] != ion_src[f"v11-{e}-emonly.json"]["seed"]:
            raise IddError(f"{e} MeV: emel and emonly must share the seed (common random numbers)")
    if not topas_runs:
        raise IddError("no TOPAS reference runs: V11 needs the emelastic and the frozen emonly runs")
    cases = DEFAULT_CASES_DIR if cases_dir is None else Path(cases_dir)
    groups = load_topas_mode_groups(topas_runs, nz, cases)
    missing = [(e, m) for e in ENERGIES for m in ("emelastic", "emonly") if (e, m) not in groups]
    if missing:
        raise IddError(f"missing TOPAS groups (energy, mode): {missing}")
    per_e: dict[str, Any] = {}
    f_ne: dict[str, Any] = {}
    el: dict[str, Any] = {}
    for e in ENERGIES:
        c = {m: groups[(e, m)]["curves"] for m in MODES if (e, m) in groups}
        per_e[str(e)] = v11_energy(ion[(e, "emel")], ion[(e, "emonly")], c["emelastic"], c["emonly"])
        on_path = ionmc_dir / f"v5-{e}-on.json"
        if "full" not in c:
            f_ne[str(e)] = {"computed": False, "reason": "TOPAS full runs not supplied"}
        elif not on_path.is_file():
            f_ne[str(e)] = {"computed": False, "reason": "v5-on partial absent; computed in v5-attribution"}
        else:
            on, on_src = load_ionmc(on_path, e, True)
            if on.shape[1] != nz[e]:
                raise IddError(f"v5-{e}-on grid differs from the V11 grid")
            if on_src["seed"] == ion_src[f"v11-{e}-emel.json"]["seed"]:
                raise IddError(f"{e} MeV: V5-on and V11 share a seed (F_ne needs independent samples)")
            f_ne[str(e)] = {"computed": True, **f_ne_test(on, ion[(e, "emel")], c["full"], c["emelastic"]),
                            "ionmc_on_input": on_src}  # fmt: skip
        if "full" in c and "noelastic" in c:
            el[str(e)] = {"computed": True, **elastic_factor_report(c["full"], c["noelastic"])}
        else:
            el[str(e)] = {"computed": False, "reason": "TOPAS full and noelastic groups not both supplied"}
    code_sha, dirty = _git("rev-parse", "HEAD"), bool(_git("status", "--porcelain"))
    return {
        "row": "V11", "pass": bool(all(per_e[str(e)]["pass"] for e in ENERGIES)),
        "plan_rule": V11_PLAN_RULE,
        "tolerances": {"F_relative": V11_TOL_F, "plateau_gain_absolute": V11_TOL_GAIN},
        "energies": per_e, "f_ne": f_ne, "elastic_factor_rows_12_13": el,
        "analysis_code_sha": code_sha, "analysis_code_dirty": dirty, "cases_dir": str(cases),
        "ionmc_inputs": ion_src,
        "engine_inputs": {f"{e}-{m}": {"seeds": g["seeds"], "runs": g["runs"], "output_sha256": g["hashes"],
                                       "fingerprints": g["fingerprints"], "family": g["family"],
                                       "binding": g["bound_cases"]}
                          for (e, m), g in sorted(groups.items())},
    }  # fmt: skip


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ionmc-dir", required=True, type=Path)
    ap.add_argument("--topas-runs", nargs="+", required=True, type=Path)
    ap.add_argument("--mcsquare-runs", nargs="+", type=Path, default=[])
    ap.add_argument("--cases-dir", type=Path, default=DEFAULT_CASES_DIR,
                    help="committed reference cases the runs must be byte-identical to")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args(argv)
    doc = build_verdict(args.ionmc_dir, args.topas_runs, args.mcsquare_runs, args.cases_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=2) + "\n")
    for eng, ev in doc["engines"].items():
        for e, rows in ev["energies"].items():
            line = "  ".join(f"{k}:{'ok' if rows[k]['pass'] else 'FAIL'}" for k in TOLERANCES)
            print(f"{eng:9s} {e} MeV  {line}")
    print(f"kappa rule: {doc['kappa_rule']['decision']}; V5 (TOPAS gating): {'PASS' if doc['pass'] else 'FAIL'}")
    return 0 if doc["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
