"""Amendment I3-FU-1 step 2: thin-slab mean energy loss in TOPAS vs ionmc Bethe and ICRU73_NEW.

Usage::

    python validation/scripts/stopping/i3fu2_evaluate.py --runs DIR --archive PATH \
        --manifest validation/reference_cases/topas/i3fu2-frozen-manifest.json \
        --output-dir validation/results/stopping/v3-008-i3fu2 --code-sha SHA

Diagnostic (plan ``validation/plans/v3-008-acceptance.md``, amendment I3-FU-1 step 2, decision
rule 2a); no gate, no validation claim.

**Run layout (``--runs``).** One subdirectory per bundle name (``gen_i3fu2_bundles.bundle_name``,
e.g. ``i3fu2-carbon-water-2p4mevu-case13-part0``) holding the materialized TOPAS outputs
``ss_exit.header``, ``ss_exit.phsp`` and ``edep.csv`` (the ``work/`` files of the run) and,
optionally, the run's ``record.json`` (its ``run_id`` is copied into the provenance) and its
``case.json`` / ``input.txt`` (if present they must equal the frozen manifest entry). Bundle
parameters (histories, seed, ion, E0, slab thickness) are read from the committed bundles next to
the manifest (``--bundles``, default the manifest's directory); the manifest must equal the one
regenerated from those bundle files, and its bundle set must be the generator's.

**Per run.** Phase-space header parsed with ``ss_evaluate.parse_header`` (the column list
:data:`EXPECTED_COLUMNS` verified on the V5-SS scorer smokes); rows parsed here because
``ss_evaluate.parse_phasespace`` fixes the exit plane at the 1 mm V5-SS slab and returns angles.
Validity (fail closed, ``RecordError``): every primary scored exactly once (deficit 0), only the
beam ion's PDG code (C 1000060120, O 1000080160, He 1000020040 as ``ss_evaluate`` case 22), Z at
the far face (world coordinates, z = t), finite values, first-scored flag 1, forward direction,
0 < <Delta E> < 20 % of the entrance energy (guards the total-kinetic-energy convention of the
energy column), and the energy-deposit cross check |EDep Sum / N / <Delta E>_T - 1| <= 0.1 %.
<Delta E>_T = A E0 - mean(E_exit), SE = std(E_exit, ddof=1) / sqrt(N).

**Predictions** (fixed by this committed code before any run): dE/dz = -rho S(E/A) with rho =
1 g/cm3, integrated over the slab by RK4 (:data:`N_RK4` steps) from E_in = A E0; S from (M) the
ionmc ``bethe_mass_stopping`` default options at ``water(78.0)`` and (I) the libdEdx ICRU73_NEW
liquid-water table (``compare_icru73.load_icru73_water``, pinned hashes), log-log interpolated
on its grid. The slab must stay above 2.0 MeV/u. r_M = Delta E_T / Delta E_M - 1, r_I =
Delta E_T / Delta E_I - 1, s = Delta E_I / Delta E_M - 1 (kept in memory only).

**Aggregates** (``i3fu2_result.json``; per ion C/O and band [2,5), [5,10), [10,20), [20,50] over
the slab-t runs): n, mean r_M, mean r_I, mean s, SE_b = sqrt(sum SE_rel^2) / n, sigma_b =
sqrt(SE_b^2 + 0.0025^2), X_b = min(0.01, |s_b| / 3), class (U checked first: X_b <= 1.645
sigma_b; I: |r_I| + 1.645 sigma_b <= X_b; M: |r_M| + 1.645 sigma_b <= X_b; L: r_M + 1.645
sigma_b < -X_b; else N); pattern ICRU / MODEL / MIXED. Controls: thickness (|r_M(t) -
r_M(t/2)| <= 0.002 at 2.4/10/50 MeV/u per ion), He (|r_M| <= 0.005 at 10 and 50 MeV/u against
ionmc's alpha Bethe; no r_I for He), fluctuations off (C, |r_M(off) - r_M(on)| <= 0.002 at 2.4
and 50 MeV/u). Nothing per point except the control differences and the two He r_M values
(TOPAS vs ionmc only). A missing evidence run (parts 0, 2) fails closed; a missing control run
makes ``conclusive`` false.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import importlib.util
import io as _io
import json
import math
import platform
import sys
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.materials import water
from ionmc.physics.stopping import bethe_mass_stopping
from ionmc.reference.parsers import ParseError, parse_topas_csv

HERE = Path(__file__).resolve().parent


def _load(name: str, path: Path) -> ModuleType:
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


gen = _load("gen_i3fu2_bundles", HERE / "gen_i3fu2_bundles.py")
icru = _load("compare_icru73", HERE / "compare_icru73.py")
ss = _load("ss_evaluate", HERE.parent / "ion_em" / "ss_evaluate.py")

RecordError = ss.RecordError
EXPECTED_COLUMNS: tuple[str, ...] = ss.EXPECTED_COLUMNS
COL_Z = EXPECTED_COLUMNS.index("Position Z [cm]")
COL_U = EXPECTED_COLUMNS.index("Direction Cosine X")
COL_V = EXPECTED_COLUMNS.index("Direction Cosine Y")
COL_E = EXPECTED_COLUMNS.index("Energy [MeV]")
COL_PDG = EXPECTED_COLUMNS.index("Particle Type (in PDG Format)")
COL_NEG = EXPECTED_COLUMNS.index(
    "Flag to tell if Third Direction Cosine is Negative (1 means true)"
)
COL_FIRST = EXPECTED_COLUMNS.index(
    "Flag to tell if this is the First Scored Particle from this History (1 means true)"
)

RESULT_SCHEMA = "ionmc-i3fu2-result-1"
RESULT_NAME = "i3fu2_result.json"
PDG = {"carbon": 1000060120, "oxygen": 1000080160, "he": 1000020040}
RHO_G_CM3 = 1.0
N_RK4 = 200  # >= 200 (amendment I3-FU-1 step 2)
E_MIN_PER_U = 2.0  # the slab must stay above this [MeV/u]
EDEP_RTOL = 0.001
MAX_LOSS_FRACTION = 0.2  # sanity bound on <Delta E> / E_in (design 0.03)
Z_RTOL = 1e-4
Z_ATOL_CM = 1e-9
METHOD_SYSTEMATIC = 0.0025
X_CAP = 0.01
K_ONE_SIDED = 1.645
THICKNESS_TOL = 0.002
HE_TOL = 0.005
FLUCT_TOL = 0.002
BANDS: tuple[tuple[str, float, float, bool], ...] = (
    ("[2,5)", 2.0, 5.0, False),
    ("[5,10)", 5.0, 10.0, False),
    ("[10,20)", 10.0, 20.0, False),
    ("[20,50]", 20.0, 50.0, True),
)
EVIDENCE_PARTS = {"carbon": 0, "oxygen": 2}
HALF_PARTS = {"carbon": 1, "oxygen": 3}
HE_PART, FLUCT_PART = 4, 5
THICKNESS_CASES = (13, 17, 21)
HE_CASES = (17, 21)
FLUCT_CASES = (13, 21)
REDISTRIBUTION_NOTE = (
    "Aggregates only: no libdEdx/ICRU 73 or Geant4 table values, no per-energy deviations, "
    "ratios or TOPAS energy-loss values (decision 0038); the control differences and the two "
    "He r_M values (TOPAS vs ionmc only) are the declared exceptions."
)

Stopping = Callable[[float], float]


# =============================================================================================
# predictions
# =============================================================================================


def model_stopping(proj: Any) -> Stopping:
    """ionmc default Bethe (``bethe_mass_stopping`` defaults, I = 78 eV) [MeV cm2/g] vs MeV/u."""
    mat = water(78.0)

    def s(e_u: float) -> float:
        return float(bethe_mass_stopping(e_u, proj, mat)[0])

    return s


def table_stopping(e_grid: NDArray[np.float64], s_grid: NDArray[np.float64]) -> Stopping:
    """Log-log interpolation of a tabulated S(E/A) on its own grid (no extrapolation)."""
    le, ls = np.log(e_grid), np.log(s_grid)

    def s(e_u: float) -> float:
        if not e_grid[0] <= e_u <= e_grid[-1]:
            raise ValueError(f"{e_u} MeV/u outside the table grid")
        return float(np.exp(np.interp(math.log(e_u), le, ls)))

    return s


def slab_energy_loss(
    e_in: float, a: int, slab_mm: float, stopping: Stopping, n_steps: int = N_RK4
) -> float:
    """Energy loss [MeV] of an ion of mass number ``a`` entering with total kinetic energy
    ``e_in`` [MeV] over ``slab_mm`` of water: RK4 on dE/dz = -rho S(E/a) [MeV/cm]. Raises
    ValueError if E/a falls to :data:`E_MIN_PER_U` or below anywhere along the slab."""
    if n_steps < 200:
        raise ValueError("at least 200 RK4 steps")
    h = slab_mm / 10.0 / n_steps  # cm

    def f(e: float) -> float:
        if e / a <= E_MIN_PER_U:
            raise ValueError(f"slab leaves the stopping domain: {e / a} MeV/u <= {E_MIN_PER_U}")
        return -RHO_G_CM3 * stopping(e / a)

    e = e_in
    for _ in range(n_steps):
        k1 = f(e)
        k2 = f(e + 0.5 * h * k1)
        k3 = f(e + 0.5 * h * k2)
        k4 = f(e + h * k3)
        e += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    f(e)  # end point inside the domain
    return e_in - e


# =============================================================================================
# per-run parsing
# =============================================================================================


def parse_exit_energies(
    header: str, rows: str, histories: int, pdg: int, slab_mm: float
) -> NDArray[np.float64]:
    """Kinetic exit energies [MeV] of the primaries of one ASCII phase space; every primary must
    be scored exactly once at the far face (z = slab thickness, world frame)."""
    h = ss.parse_header(header)
    if h["original_histories"] != histories:
        raise RecordError(
            f"phase space: original histories {h['original_histories']} != {histories}"
        )
    if h["reached"] != h["scored"] or h["scored"] != histories:
        raise RecordError(
            f"phase space: reached {h['reached']}, scored {h['scored']} of {histories} "
            "(deficit must be 0)"
        )
    try:
        data = np.loadtxt(_io.StringIO(rows), dtype=np.float64, ndmin=2)
    except ValueError as exc:
        raise RecordError(f"phase space rows are not numeric: {exc}") from exc
    if data.shape != (histories, len(EXPECTED_COLUMNS)):
        raise RecordError(f"phase space has shape {data.shape}, not ({histories}, 10)")
    if not np.all(np.isfinite(data)):
        raise RecordError("phase space holds non-finite values")
    if not np.all(data[:, COL_PDG] == pdg):
        raise RecordError(f"phase space holds particles other than PDG {pdg}")
    if not np.all(data[:, COL_NEG] == 0):
        raise RecordError("phase space: a particle with negative third direction cosine")
    if not np.all(data[:, COL_FIRST] == 1):
        raise RecordError("phase space: a row is not the first scored particle of its history")
    z_cm = slab_mm / 10.0
    if not np.all(np.abs(data[:, COL_Z] - z_cm) <= Z_RTOL * z_cm + Z_ATOL_CM):
        raise RecordError(f"phase space: Z differs from the far face {z_cm} cm")
    if not np.all(data[:, COL_U] ** 2 + data[:, COL_V] ** 2 < 1.0):
        raise RecordError("phase space: direction cosines with u^2 + v^2 >= 1")
    e = data[:, COL_E]
    if not np.all(e > 0.0):
        raise RecordError("phase space: non-positive kinetic energy")
    return np.asarray(e, dtype=np.float64)


def edep_per_history(text: str, histories: int) -> float:
    """Mean energy deposit per history [MeV] of the 1 x 1 x 1 EnergyDeposit CSV scorer."""
    try:
        s = parse_topas_csv(text, "edep.csv")
    except ParseError as exc:
        raise RecordError(str(exc)) from exc
    if s.bins != (1, 1, 1) or "Sum" not in s.values:
        raise RecordError("edep.csv must be a single bin with a Sum")
    if s.unit != "MeV":
        raise RecordError(f"edep.csv unit {s.unit!r} is not MeV")
    v = float(s.values["Sum"][0, 0, 0])
    if not math.isfinite(v) or v <= 0.0:
        raise RecordError("edep.csv Sum is missing or not positive")
    return v / histories


def measure_run(run_dir: Path, case_doc: dict[str, Any]) -> dict[str, Any]:
    """<Delta E>_T and its SE of one run directory (validity checks fail closed)."""
    info = case_doc["i3fu2"]
    n = int(case_doc["histories"])
    a, e0 = int(info["a"]), float(info["e0_per_u_mev"])
    e_in = a * e0
    try:
        header = (run_dir / "ss_exit.header").read_text()
        rows = (run_dir / "ss_exit.phsp").read_text()
        edep_text = (run_dir / "edep.csv").read_text()
    except OSError as exc:
        raise RecordError(f"{run_dir.name}: missing output ({exc})") from exc
    try:
        e = parse_exit_energies(header, rows, n, PDG[info["ion"]], float(info["slab_mm"]))
    except RecordError as exc:
        raise RecordError(f"{run_dir.name}: {exc}") from exc
    de = e_in - float(np.mean(e))
    if not 0.0 < de < MAX_LOSS_FRACTION * e_in:
        raise RecordError(
            f"{run_dir.name}: mean energy loss outside (0, {MAX_LOSS_FRACTION} E_in); exit "
            "energy not the total kinetic energy?"
        )
    se = float(np.std(e, ddof=1)) / math.sqrt(n)
    ed = edep_per_history(edep_text, n)
    if abs(ed / de - 1.0) > EDEP_RTOL:
        raise RecordError(
            f"{run_dir.name}: energy-deposit cross check fails (|EDep/N / <dE> - 1| > {EDEP_RTOL})"
        )
    run_id = None
    rec = run_dir / "record.json"
    if rec.is_file():
        run_id = json.loads(rec.read_text()).get("run_id")
    return {"de": de, "se": se, "scored": n, "run_id": run_id}


# =============================================================================================
# bundles and manifest
# =============================================================================================


def load_bundles(bundles_dir: Path, manifest_path: Path) -> tuple[dict[str, Any], str]:
    """Verify the manifest against the committed bundle files; return ({name: case.json},
    manifest text)."""
    try:
        text = manifest_path.read_text()
        doc = json.loads(text)
    except (OSError, ValueError) as exc:
        raise RecordError(f"manifest unreadable: {exc}") from exc
    if doc.get("schema") != gen.MANIFEST_SCHEMA:
        raise RecordError(f"manifest schema {doc.get('schema')!r} != {gen.MANIFEST_SCHEMA}")
    if sorted(doc.get("bundles", {})) != gen.bundle_names():
        raise RecordError("manifest bundle set differs from the generator's")
    try:
        fresh = gen.manifest_text(bundles_dir)
    except OSError as exc:
        raise RecordError(f"committed bundle missing: {exc}") from exc
    if fresh != text:
        raise RecordError("manifest does not match the committed bundle files")
    cases = {
        name: json.loads((bundles_dir / name / "case.json").read_text())
        for name in gen.bundle_names()
        if name != gen.SMOKE_NAME
    }
    return cases, text


def check_run_inputs(run_dir: Path, entry: dict[str, str]) -> None:
    for f in ("case.json", "input.txt"):
        p = run_dir / f
        if p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest() != entry[f]:
            raise RecordError(f"{run_dir.name}/{f} differs from the frozen manifest entry")


# =============================================================================================
# aggregation
# =============================================================================================


def classify(r_m: float, r_i: float, s: float, sigma: float) -> str:
    """Decision rule 2a for one band (U checked first)."""
    x = min(X_CAP, abs(s) / 3.0)
    k = K_ONE_SIDED * sigma
    if x <= k:
        return "U"
    if abs(r_i) + k <= x:
        return "I"
    if abs(r_m) + k <= x:
        return "M"
    if r_m + k < -x:
        return "L"
    return "N"


def band_of(e0: float) -> str:
    for name, lo, hi, closed in BANDS:
        if lo <= e0 < hi or (closed and e0 == hi):
            return name
    raise ValueError(f"E0 {e0} outside the bands")


def band_aggregates(points: list[dict[str, float]]) -> dict[str, dict[str, Any]]:
    """Per band: n, means of r_M, r_I, s, SE_b, sigma_b, X_b, class."""
    out: dict[str, dict[str, Any]] = {}
    for name, *_ in BANDS:
        pts = [p for p in points if band_of(p["e0"]) == name]
        n = len(pts)
        if n == 0:
            raise RecordError(f"band {name} has no run")
        r_m = float(np.mean([p["r_M"] for p in pts]))
        r_i = float(np.mean([p["r_I"] for p in pts]))
        s = float(np.mean([p["s"] for p in pts]))
        se_b = math.sqrt(sum(p["se_rel"] ** 2 for p in pts)) / n
        sigma = math.sqrt(se_b**2 + METHOD_SYSTEMATIC**2)
        out[name] = {
            "n": n,
            "mean_r_M": r_m,
            "mean_r_I": r_i,
            "mean_s": s,
            "se_b": se_b,
            "sigma_b": sigma,
            "x_b": min(X_CAP, abs(s) / 3.0),
            "class": classify(r_m, r_i, s, sigma),
        }
    return out


def pattern(bands: dict[str, dict[str, dict[str, Any]]]) -> str:
    def cls(ion: str, b: str) -> str:
        return str(bands[ion][b]["class"])

    ions = ("carbon", "oxygen")
    if all(
        cls(i, "[2,5)") == "I" and cls(i, "[5,10)") == "I" and cls(i, "[10,20)") in ("I", "U")
        for i in ions
    ):
        return "ICRU"
    if all(cls(i, b) in ("M", "L") for i in ions for b in ("[2,5)", "[5,10)")):
        return "MODEL"
    return "MIXED"


# =============================================================================================
# evaluation
# =============================================================================================


def evaluate(
    runs_dir: Path,
    archive: Path,
    manifest_path: Path,
    code_sha: str,
    bundles_dir: Path | None = None,
    *,
    archive_sha256: str | None = None,
    member_sha256: dict[str, str] | None = None,
) -> dict[str, Any]:
    """The result document (aggregates only). ``archive_sha256``/``member_sha256`` override the
    pinned libdEdx hashes of ``compare_icru73`` (tests with a synthetic archive only)."""
    bundles_dir = manifest_path.parent if bundles_dir is None else bundles_dir
    archive_sha = icru.ARCHIVE_SHA256 if archive_sha256 is None else archive_sha256
    members = dict(icru.MEMBER_SHA256 if member_sha256 is None else member_sha256)
    cases, manifest = load_bundles(bundles_dir, manifest_path)
    entries = json.loads(manifest)["bundles"]
    tables = icru.load_icru73_water(
        archive, archive_sha256=archive_sha, member_sha256=members, ions=(6, 8)
    )

    def stopping_m(ion: str) -> Stopping:
        return model_stopping(gen.IONS[ion][0])

    s_model = {ion: stopping_m(ion) for ion in gen.IONS}
    s_table = {
        ion: table_stopping(tables[z].energy_per_u_mev, tables[z].stopping_mev_cm2_g)
        for ion, z in (("carbon", 6), ("oxygen", 8))
    }

    results: dict[tuple[int, int], dict[str, float]] = {}
    runs: dict[str, Any] = {}
    missing: list[str] = []
    for case, part in gen.evidence_keys():
        name = gen.bundle_name(case, part)
        run_dir = runs_dir / name
        if not run_dir.is_dir():
            if part in EVIDENCE_PARTS.values():
                raise RecordError(f"missing evidence run {name}")
            missing.append(name)
            continue
        check_run_inputs(run_dir, entries[name])
        doc = cases[name]
        info = doc["i3fu2"]
        m = measure_run(run_dir, doc)
        ion, a, e0, t = info["ion"], int(info["a"]), float(info["e0_per_u_mev"]), info["slab_mm"]
        try:
            de_m = slab_energy_loss(a * e0, a, float(t), s_model[ion])
            de_i = slab_energy_loss(a * e0, a, float(t), s_table[ion]) if ion in s_table else None
        except ValueError as exc:
            raise RecordError(f"{name}: prediction failed: {exc}") from exc
        rec: dict[str, float] = {
            "e0": e0,
            "r_M": m["de"] / de_m - 1.0,
            "se_rel": m["se"] / m["de"],
        }
        if de_i is not None:
            rec["r_I"] = m["de"] / de_i - 1.0
            rec["s"] = de_i / de_m - 1.0
        results[(case, part)] = rec
        runs[name] = {
            "seed": int(doc["seeds"][0]),
            "histories": int(doc["histories"]),
            "scored": m["scored"],
            "run_id": m["run_id"],
        }

    bands = {
        ion: band_aggregates([results[(c, p)] for (c, p) in results if p == part])
        for ion, part in EVIDENCE_PARTS.items()
    }

    def diffs(pairs: list[tuple[tuple[int, int], tuple[int, int]]], tol: float) -> dict:
        d = [
            results[x]["r_M"] - results[y]["r_M"] for x, y in pairs if x in results and y in results
        ]
        complete = len(d) == len(pairs)
        return {
            "differences": d,
            "tolerance": tol,
            "complete": complete,
            "pass": complete and all(abs(v) <= tol for v in d),
        }

    thickness = {
        ion: {
            "e0_per_u_mev": [float(gen.CASES[c]) for c in THICKNESS_CASES],
            **diffs(
                [((c, EVIDENCE_PARTS[ion]), (c, HALF_PARTS[ion])) for c in THICKNESS_CASES],
                THICKNESS_TOL,
            ),
        }
        for ion in ("carbon", "oxygen")
    }
    he_r = [results[(c, HE_PART)]["r_M"] for c in HE_CASES if (c, HE_PART) in results]
    he_complete = len(he_r) == len(HE_CASES)
    helium = {
        "e0_per_u_mev": [float(gen.CASES[c]) for c in HE_CASES],
        "r_M": he_r,
        "tolerance": HE_TOL,
        "complete": he_complete,
        "pass": he_complete and all(abs(v) <= HE_TOL for v in he_r),
    }
    fluct = {
        "e0_per_u_mev": [float(gen.CASES[c]) for c in FLUCT_CASES],
        **diffs([((c, FLUCT_PART), (c, EVIDENCE_PARTS["carbon"])) for c in FLUCT_CASES], FLUCT_TOL),
    }
    controls_pass = (
        all(thickness[i]["pass"] for i in thickness) and helium["pass"] and fluct["pass"]
    )
    low_u = any(bands[i][b]["class"] == "U" for i in bands for b in ("[2,5)", "[5,10)"))
    try:
        ionmc_version: str | None = importlib.metadata.version("ionmc")
    except importlib.metadata.PackageNotFoundError:
        ionmc_version = None
    return {
        "schema": RESULT_SCHEMA,
        "task_id": "V3-008",
        "amendment": "I3-FU-1 step 2",
        "claim": "diagnostic; no gate; no validation claim",
        "code_sha": code_sha,
        "method": {
            "quantity": "<Delta E>_T = A E0 - mean(E_exit); predictions RK4 over the slab",
            "rk4_steps": N_RK4,
            "model": "ionmc bethe_mass_stopping default options, water(78.0)",
            "table": "libdEdx ICRU73_NEW liquid water (target 276), log-log interpolated",
            "method_systematic": METHOD_SYSTEMATIC,
            "x_cap": X_CAP,
            "k_one_sided": K_ONE_SIDED,
        },
        "bands": bands,
        "pattern": pattern(bands),
        "controls": {
            "thickness": thickness,
            "helium": helium,
            "fluctuations_off": fluct,
            "missing_runs": missing,
            "controls_pass": controls_pass,
        },
        "conclusive": controls_pass and not low_u,
        "provenance": {
            "archive_sha256": archive_sha,
            "member_sha256": members,
            "manifest_bundle_set_digest": json.loads(manifest)["bundle_set_digest"],
            "code_sha": code_sha,
            "runs": runs,
            "versions": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "ionmc": ionmc_version,
            },
        },
        "redistribution_note": REDISTRIBUTION_NOTE,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, default=gen.DEFAULT_MANIFEST)
    ap.add_argument("--bundles", type=Path, default=None)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--code-sha", required=True)
    a = ap.parse_args(argv)
    doc = evaluate(a.runs, a.archive, a.manifest, a.code_sha, a.bundles)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    out = a.output_dir / RESULT_NAME
    out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    print("wrote", out, "pattern", doc["pattern"], "conclusive", doc["conclusive"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
