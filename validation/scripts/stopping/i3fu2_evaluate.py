"""Amendment I3-FU-1 step 2: thin-slab mean energy loss in TOPAS vs ionmc Bethe and ICRU73_NEW.

Usage::

    python validation/scripts/stopping/i3fu2_evaluate.py --runs DIR --archive PATH \
        --code-sha SHA --smokes DIR \
        [--records validation/results/stopping/v3-008-i3fu2/topas-runs] \
        [--review-pin validation/results/stopping/v3-008-i3fu2/i3fu2-review.json] \
        [--markers-pin validation/results/stopping/v3-008-i3fu2/i3fu2-markers-review.json] \
        [--markers validation/results/stopping/v3-008-i3fu2/log-markers.json] \
        [--manifest validation/reference_cases/topas/i3fu2-frozen-manifest.json] \
        [--repo .] [--output-dir validation/results/stopping/v3-008-i3fu2]

Diagnostic (plan ``validation/plans/v3-008-acceptance.md``, amendment I3-FU-1 step 2 and decision
rule 2a, revision 2: the dispositions of Codex REVIEW-22f384ed and REVIEW-aa4e7b1b and of the
independent reviewer's passes 2 and 3); no gate, no validation claim.

**Inputs.** ``--repo`` (default: this checkout) is read only as committed at its HEAD
(``r2_io.CommittedReader``: a working copy that differs fails closed); ``--records``,
``--smokes``, both review pins, ``--markers`` and ``--manifest`` are paths in it. Records
(committed): ``<dir>/index.json`` = ``{"runs": [{"bundle", "run_id", ...}, ...]}`` (``--records``:
one entry per evidence or offered repeat bundle with its ``case`` and ``part``; ``--smokes``:
exactly the four smokes of :data:`gen.SMOKES`) and per run ``<run_id>/record.json`` (the signed
reference-runner record, verbatim) and ``<run_id>/artifacts.json`` (the
``list_reference_artifacts`` result). Materialized outputs (``--runs``, not committed, evidence,
repeats and smokes alike): ``<run_id>/`` holding the listed files at their listed paths
(``request.json``, ``stdout.txt``, ``stderr.txt``, ``work/ss_exit.header``,
``work/ss_exit.phsp``, ``work/edep.csv``).

**Sealed executable code (Codex REVIEW-22f384ed finding 1, REVIEW-aa4e7b1b finding 5).** The
evaluator refuses to run unless the ``--repo`` worktree is clean (``git status --porcelain
--untracked-files=no`` empty) and every file of the executable closure is byte-identical on
disk, in the file the interpreter loaded, and as committed at HEAD, at ``--code-sha``, at the
evaluator review's ``head_sha`` (pin A) and at every smoke's code_sha. The closure
(:func:`executable_closure`) is computed from ``sys.modules`` after all imports: every ``ionmc``
module and every validation script loaded by path (:data:`SCRIPT_MODULES`, with the
attestation module), each with a ``__file__`` under this checkout, plus :data:`SEALED_FILES`
(``uv.lock``, ``pyproject.toml``). A module imported after the seal is a ``RecordError``. The
running numpy version must equal the version pinned for numpy in ``uv.lock`` and Python
major.minor the lower bound of its ``requires-python``. Declared limitation: the bytes of
third-party binaries (numpy, the Python interpreter) are not sealed; the lock file is. Written:
``provenance.sealed_sources`` {path: sha256} and ``provenance.environment``.

**Review pins and the plan note.** Two verbatim codex-worker review records: pin A
``i3fu2-review.json`` (the review of the commit that seals the evaluator) and pin B
``i3fu2-markers-review.json`` (the review of the commit that adds the log-marker file, derived
from the smoke logs, and pin A). Each must hold :data:`rio.REVIEW_FIELDS` and ``report.verdict``,
``status`` "passed", ``report.verdict`` "pass", ``task_id`` "V3-008", a 40-hex ``head_sha`` (equal
to ``report.head_sha`` where present), a 64-hex ``report_sha256`` and ISO ``started_at`` <
``finished_at``. The review service does not sign its records (a signed review record is the
future item V3-INFRA-005): as for the V5-SS chain, the binding is to the reviewed, committed
bytes. The acceptance plan at HEAD and at ``code_sha`` must contain exactly one line matching
:data:`PLAN_NOTE_RE`, identical at both::

    **Dated note (<YYYY-MM-DD>, orchestrator; I3-FU-1 step-2 review pins).** evaluator review
    <REVIEW-...> at <40-hex>; markers review <REVIEW-...> at <40-hex>

(one line in the plan), whose ids and SHAs must equal pin A's and pin B's ``review_id`` /
``head_sha``. Pin A's ``head_sha`` is an ancestor of pin B's, pin B's an ancestor of
``code_sha``, and ``code_sha`` an ancestor of (or equal to) HEAD (HEAD cannot equal ``code_sha``:
the signed records are committed after the runs). Pin A and the marker file are byte-identical
at pin B's ``head_sha``, at ``code_sha`` and at HEAD; pin B at ``code_sha`` and HEAD. The frozen
manifest is byte-identical at HEAD, at ``code_sha``, at pin A's ``head_sha`` and at every smoke's
code_sha, has the generator's bundle set and reproduces its own ``bundle_set_digest``. The
trust map is ``infrastructure/keys/trusted_attestation_keys.json`` at HEAD.

**Chronology.** Smokes (all four; missing or late: ``RecordError``, nothing graded): the
record's code_sha carries pin A (pin A byte-identical there, pin A's ``head_sha`` an ancestor of
it) and is an ancestor of (or equal to) pin B's ``head_sha``; ``started_at`` after pin A's
``finished_at`` and ``finished_at`` before pin B's ``started_at`` (the marker file is derived from
the smoke logs, so its review begins after the smokes end). Evidence and repeat runs: code_sha
equals ``--code-sha`` and ``started_at`` is after pin B's (and pin A's) ``finished_at``. Every
run: ``started_at`` <= ``finished_at``.

**Log markers.** ``--markers`` is the committed JSON ``{"fluct_on_marker", "fluct_off_marker",
"step_line_re", "track_header_re"}``; every value must be a non-empty regular expression
(compiled with ``re.M``), else ``RecordError``. Neither fluctuation marker may match a run's
native input text (they must identify the Geant4 ``G4EmParameters`` printout, not the TOPAS
parameter echo): ``RecordError``. The fluctuation markers must discriminate over all
authenticated runs and smokes: ``fluct_off_marker`` matches every fluctuation-off stdout and no
fluctuation-on stdout, ``fluct_on_marker`` the converse; otherwise the fluctuation control is
FAILED (never passed).

**Per-run authentication and binding (any failure: ``RecordError`` naming the run; Codex
REVIEW-cfdff398 finding 1, mirroring ``ss_evaluate.load_run`` / ``r2_io.load_signed_run``).**
Ed25519 attestation of ``record.json`` by a trusted ``reference-runner`` key; schema 3, not timed
out, engine topas, ``succeeded`` true, ``exit_code`` 0, no missing outputs; run id derived from
the configuration id; the chronology above; the listing reproduces the signed
``output_files_sha256`` / ``output_count`` and the signed stdout/stderr/request digests; every
file read hashes to its listed sha256; the request's ``input_files`` equal the bundle committed
at the run's code_sha, whose ``case.json`` / ``input.txt`` equal the manifest entry; the request
case equals that ``case.json``; the seed is 20619000 + 100 part + case (smokes: the generator's
rehearsal seed); the ``i3fu2`` block matches the index entry and the kind (evidence / repeat /
smoke); configuration markers in stdout: the ``Ph/ListProcesses`` block, the opt4 module marker
and the EM-only check of ``r2_io.module_check`` (configuration "idd-em"), the frozen physics
lines of the bundle and an ``ionIoni:  for <GenericIon | alpha>`` block.

**Run validity (an INVALID run is recorded with its reason and triggers the repeat rule).** stderr
empty, no stdout line outside "Begin processing" mentions killed or stuck tracks, every ``Particle
source Beam: Total number of histories: N`` equals the bundle's histories; phase space: every
primary scored exactly once (deficit 0), only the beam ion's PDG code, Z at the far face, finite
values, first-scored flag 1, forward direction, 0 < <Delta E> < 20 % of the entrance energy; the
energy-deposit cross check |EDep Sum / N / <Delta E>_T - 1| <= 0.1 %.

**One-step premise (step-count smokes).** For both step-count smokes: every one of the 50
primaries has exactly one step whose pre-step point lies in the Water slab (0 <= z_pre < slab
thickness, with the printed-precision tolerance of :func:`count_slab_steps`), that step ends at
or beyond the far face (its limiting process is expected to be Transportation), and the smoke's
far-face phase space holds exactly one primary row per history (the evidence-run checks:
deficit 0, PDG, far-face z). Otherwise ``conclusive`` is false with the reason
``one_step_premise``. Written: ``one_step_premise = {checked, smokes: {name: {primaries,
max_steps_in_slab, traversing_primaries, phase_space_complete}}, pass}`` (counts only).

**Predictions** (fixed by this committed code before any run): dE/dz = -rho S(E/A), rho = 1 g/cm3,
RK4 over the slab (:data:`N_RK4` steps) from E_in = A E0; S from (M) ``bethe_mass_stopping``
default options at ``water(78.0)`` and (I) the libdEdx ICRU73_NEW liquid-water table
(``compare_icru73.load_icru73_water``, pinned hashes), log-log interpolated on its grid. Per run
r_M = Delta E_T / Delta E_M - 1, r_I = Delta E_T / Delta E_I - 1, s = Delta E_I / Delta E_M - 1,
SE_rel = SE / Delta E_T.

**Stepping model.** Geant4's range-table branch freezes the effective charge at the pre-step
energy for a one-step slab. Emulated with ionmc's own model: S(E) = q^2(E) S_0(E) with q from
``pierce_blann_charge`` and S_0 = S / q^2 numerically; the algorithm's loss integrates
S(E) q^2(E_pre) / q^2(E) over the slab with the same RK4 integrator as the exact loss. g(f) =
Delta E_alg / Delta E_exact - 1 at f1 = 0.03 (slab t) and f2 = 0.048 (slab 1.6 t) per point (He
with z = 2 through the same callable). Per point: Delta_pred = g(f2) - g(f1); Delta_meas =
r_M(f2) - r_M(f1); consistency |Delta_meas - Delta_pred| <= :data:`CONSISTENCY_TOL` (a failure
makes the point INVALID: repeat rule); model residual of the extrapolation R_model = (f2 g(f1) -
f1 g(f2)) / (f2 - f1) with |R_model| <= :data:`MODEL_RESIDUAL_LIMIT` at every C/O/He point (else
``RecordError``: the design premise fails). The 0.25 % method systematic in sigma_b stays (its
0.2 % stepping term is twice the residual limit).

**Points.** Per ion (C parts 0/1, O parts 2/3) and case: the zero-thickness (Richardson)
extrapolation r_0 = (f2 r(f1) - f1 r(f2)) / (f2 - f1) of r_M, r_I and s, f1, f2 read from the
bundles' ``thickness_factor`` and asserted (:data:`FACTORS`). SE_rel,0 = sqrt((f2 SE_rel(f1))^2
+ (f1 SE_rel(f2))^2) / (f2 - f1).

**Repeat rule (deterministic).** Repeat bundle = part + 8 (:data:`gen.REPEAT_OFFSET`). A validity
failure of a run requires the repeat of exactly that run; a consistency failure of a point (both
runs valid) requires the repeat of BOTH slabs, and the repeat pair replaces the pair. The offered
repeats of a point must be exactly the required set: a repeat offered without a recorded
failure, a repeat of a valid run, or a partial set is a ``RecordError``. Required repeats not
offered: the point is invalid ("no repeat"). A second failure makes the point invalid; an invalid
evidence point makes the result INCONCLUSIVE.

**Aggregates** (``i3fu2_result.json``; per ion C/O and band [2,5), [5,10), [10,20), [20,50] over
the valid points): n, mean r_M, mean r_I, mean s, SE_b = sqrt(sum SE_rel,0^2) / n, sigma_b =
sqrt(SE_b^2 + 0.0025^2), X_b = min(0.01, |s_b| / 3), the class (relative to the sign of s_b;
U first: X_b <= 1.645 sigma_b; I: |r_I| + 1.645 sigma_b <= X_b; H: sign(s_b) r_I - 1.645
sigma_b > X_b; M: |r_M| + 1.645 sigma_b <= X_b; L: -sign(s_b) r_M - 1.645 sigma_b > X_b; else N)
and the band maxima of |r_M(f2) - r_M(f1)|, |Delta_meas - Delta_pred|, |Delta_pred| and
|R_model|. The table-relative aggregates (mean r_M, mean r_I, mean s, X_b and the He
``max_abs_r_m``) are written rounded to 0.1 percentage point; classes and the pattern are
computed from the full-precision values. Pattern: ICRU iff both ions are I/H in [2,5) and [5,10)
and I/H/U in [10,20); MODEL iff both are M/L in [2,5) and [5,10) and M/L/U in [10,20); else
MIXED (``pattern_before_controls``). Robustness to the declared systematic (pre-declared): the
classes and the pattern are computed a second time with the method systematic doubled
(``classes_sigma_doubled``, ``pattern_sigma_doubled``); ``robust`` = (``pattern_before_controls``
== ``pattern_sigma_doubled``). ``pattern`` is INCONCLUSIVE whenever ``conclusive`` is false: a
control failed or is incomplete, an evidence point is invalid, a band of [2,5) / [5,10) is U,
the one-step premise fails, or the pattern is not robust (``not_robust_to_systematic``).
Controls: He (part 4, cases 17/21): pass iff both points valid and |r_M| <= 0.005; fluctuations
off (part 5 vs the C slab-t run of the same case, cases 13/21): pass iff complete, the marker
set discriminates and |r_M(off) - r_M(on)| <= 0.002. Nothing per point is written; reasons carry
no values (decision 0038).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import importlib.util
import io as _io
import json
import math
import platform
import re
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.materials import water
from ionmc.physics.projectiles import Projectile
from ionmc.physics.stopping import bethe_mass_stopping, kinematics, pierce_blann_charge
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
rio = ss.io  # validation/scripts/ion_em/r2_io.py (record verification helpers)

RecordError = rio.RecordError
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

RESULT_SCHEMA = "ionmc-i3fu2-result-4"
RESULT_NAME = "i3fu2_result.json"
RESULTS_DIR = "validation/results/stopping/v3-008-i3fu2"
REVIEW_PIN = f"{RESULTS_DIR}/i3fu2-review.json"  # pin A: the evaluator seal review
MARKERS_PIN = f"{RESULTS_DIR}/i3fu2-markers-review.json"  # pin B: markers + pin A review
MARKERS = f"{RESULTS_DIR}/log-markers.json"
RECORDS = f"{RESULTS_DIR}/topas-runs"
MANIFEST = f"{gen.REL_PREFIX}/i3fu2-frozen-manifest.json"
PLAN_PATH: str = rio.PLAN_PATH
TASK_ID = "V3-008"
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
#: The dated plan note naming both review pins (one line; parsed at HEAD and at code_sha).
PLAN_NOTE_RE = re.compile(
    r"^\*\*Dated note \((?P<date>\d{4}-\d{2}-\d{2}), orchestrator; I3-FU-1 step-2 review "
    r"pins\)\.\*\* evaluator review (?P<eval_id>REVIEW-[0-9A-Za-z_-]+) at "
    r"(?P<eval_sha>[0-9a-f]{40}); markers review (?P<markers_id>REVIEW-[0-9A-Za-z_-]+) at "
    r"(?P<markers_sha>[0-9a-f]{40})(?![0-9A-Za-z])",
    re.M,
)
MARKER_KEYS = ("fluct_on_marker", "fluct_off_marker", "step_line_re", "track_header_re")
#: Validation scripts loaded by path (module names) whose files belong to the executable
#: closure, besides this module and every ``ionmc`` module (:func:`executable_closure`).
SCRIPT_MODULES = (
    "gen_i3fu2_bundles",
    "compare_icru73",
    "ss_evaluate",
    "r2_io",
    "r2_evaluate",
    "ss_budget",
    "ss_wrapper_template",
    "infrastructure.experiment_v3.attest",
)
#: Non-module files sealed with the closure (the dependency lock and the project metadata).
SEALED_FILES = ("uv.lock", "pyproject.toml")
REPO_ROOT: Path = gen.REPO_ROOT
ZEPS_MM = 1e-12  # float guard of the printed-precision comparisons [mm]
POST_STEP_MIN_TOL_MM = 1e-6  # the counted step must end at >= slab - max(tol, this) [mm]
PDG = {"carbon": 1000060120, "oxygen": 1000080160, "he": 1000020040}
BEAM_PARTICLE = {"carbon": "GenericIon", "oxygen": "GenericIon", "he": "alpha"}
EM_CONFIG = "idd-em"  # r2 configuration whose frozen module set / physics lines these bundles use
RHO_G_CM3 = 1.0
N_RK4 = 200  # >= 200 (amendment I3-FU-1 step 2)
E_MIN_PER_U = 2.0  # the slab must stay above this [MeV/u]
EDEP_RTOL = 0.001
MAX_LOSS_FRACTION = 0.2  # sanity bound on <Delta E> / E_in (design 0.03 / 0.048)
Z_RTOL = 1e-4
Z_ATOL_CM = 1e-9
METHOD_SYSTEMATIC = 0.0025
X_CAP = 0.01
K_ONE_SIDED = 1.645
CONSISTENCY_TOL = 0.003  # |Delta_meas - Delta_pred| per point (Codex REVIEW-22f384ed finding 4)
MODEL_RESIDUAL_LIMIT = 0.001  # |R_model| per point (design premise)
FACTORS = (0.03, 0.048)  # thickness factors f1 (slab t), f2 (slab 1.6 t) of every pair
HE_TOL = 0.005
FLUCT_TOL = 0.002
ROUND_DECIMALS = 3  # 0.1 percentage point in relative units
BANDS: tuple[tuple[str, float, float, bool], ...] = (
    ("[2,5)", 2.0, 5.0, False),
    ("[5,10)", 5.0, 10.0, False),
    ("[10,20)", 10.0, 20.0, False),
    ("[20,50]", 20.0, 50.0, True),
)
LOW_BANDS = ("[2,5)", "[5,10)")
PAIR_PARTS = {"carbon": (0, 1), "oxygen": (2, 3)}  # (slab t, slab 1.6 t)
HE_PART, FLUCT_PART = 4, 5
HE_CASES = (17, 21)
FLUCT_CASES = (13, 21)
STEP_SMOKES = ("i3fu2-smoke-carbon-water-2p4mevu-steps", "i3fu2-smoke-carbon-water-50mevu-steps")
ALL_SMOKES: tuple[str, ...] = tuple(gen.SMOKES)  # the four smokes (two execution, two step-count)
ICRU_LIKE, MODEL_LIKE = ("I", "H"), ("M", "L")
_IONIONI = "^ionIoni:\\s+for\\s+{}\\s"
#: printed length unit -> mm (Geant4 G4BestUnit "Length" symbols)
LENGTH_UNITS_MM = {
    "fm": 1e-12,
    "Ang": 1e-7,
    "nm": 1e-6,
    "um": 1e-3,
    "mm": 1.0,
    "cm": 10.0,
    "m": 1e3,
    "km": 1e6,
    "pc": 3.0856775814913673e19,
}
CONSISTENCY_REASON = f"thickness consistency |Delta_meas - Delta_pred| > {CONSISTENCY_TOL}"
PREMISE_REASON = "one_step_premise"
ROBUST_REASON = "not_robust_to_systematic"
REDISTRIBUTION_NOTE = (
    "Aggregates only: no libdEdx/ICRU 73 or Geant4 table values, no per-energy deviations, "
    "ratios or TOPAS energy-loss values (decision 0038); the controls are written as pass/fail "
    "and maxima, the thickness and stepping checks as maxima per ion and band, reasons without "
    "values."
)
ROUNDING_NOTE = (
    "mean_r_M, mean_r_I, mean_s, x_b and controls.helium.max_abs_r_m are written rounded to "
    "0.1 percentage point (3 decimals in relative units); the classes and the pattern were "
    "computed from the full-precision values (pre-declared)."
)

Stopping = Callable[[float], float]


class RunInvalid(Exception):
    """A run failed a validity criterion (recorded with its reason; the repeat rule applies)."""


def _round(x: float | None) -> float | None:
    return None if x is None else round(float(x), ROUND_DECIMALS)


# =============================================================================================
# predictions and the stepping model
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


_PREDICTIONS: dict[tuple[str, int, float, float], float] = {}


def predicted_loss(key: str, stopping: Stopping, a: int, e0: float, slab_mm: float) -> float:
    """:func:`slab_energy_loss` memoised by (stopping identity ``key``, a, E0, slab): the
    predictions depend on nothing else."""
    k = (key, a, e0, slab_mm)
    if k not in _PREDICTIONS:
        _PREDICTIONS[k] = slab_energy_loss(a * e0, a, slab_mm, stopping)
    return _PREDICTIONS[k]


def effective_charge_sq(proj: Projectile, e_u: float) -> float:
    """q^2 of ionmc's default model (``pierce_blann_charge``, Barkas-125 form) at E/A = ``e_u``."""
    beta2 = kinematics(np.atleast_1d(np.float64(e_u)), proj)[0]
    q = pierce_blann_charge(np.sqrt(beta2), proj)
    return float(q[0]) ** 2


def frozen_charge_stopping(proj: Projectile, stopping: Stopping, e_pre_u: float) -> Stopping:
    """The one-step algorithm's stopping: q^2(E_pre) S_0(E) with S_0 = S / q^2 of the default
    model, i.e. S(E) q^2(E_pre) / q^2(E) (the effective charge frozen at the pre-step energy)."""
    q2_pre = effective_charge_sq(proj, e_pre_u)

    def s(e_u: float) -> float:
        return stopping(e_u) * q2_pre / effective_charge_sq(proj, e_u)

    return s


_STEPPING: dict[tuple[str, int], dict[str, float]] = {}


def stepping_bias(ion: str, case: int) -> dict[str, float]:
    """Model-side emulation of the single-step algorithm for one point (memoised): g(f1),
    g(f2) with g = Delta E_alg / Delta E_exact - 1 over the generator's slabs t and 1.6 t,
    ``d_pred`` = g(f2) - g(f1) and ``r_model`` = (f2 g(f1) - f1 g(f2)) / (f2 - f1)."""
    key = (ion, case)
    if key not in _STEPPING:
        proj, _z, a, _p = gen.IONS[ion]
        e0 = float(gen.CASES[case])
        exact = model_stopping(proj)
        alg = frozen_charge_stopping(proj, exact, e0)
        g = []
        for mult in (1, 2):
            t = gen.slab_mm(ion, case, mult)
            de_exact = predicted_loss(f"model:{ion}", exact, a, e0, t)
            de_alg = slab_energy_loss(a * e0, a, t, alg)
            g.append(de_alg / de_exact - 1.0)
        _STEPPING[key] = {
            "g1": g[0],
            "g2": g[1],
            "d_pred": g[1] - g[0],
            "r_model": richardson(g[0], g[1]),
        }
    return _STEPPING[key]


def stepping_model_check() -> dict[str, dict[int, dict[str, float]]]:
    """:func:`stepping_bias` for every C/O point and both He points; any |R_model| above
    :data:`MODEL_RESIDUAL_LIMIT` is a ``RecordError`` (the design premise fails)."""
    out: dict[str, dict[int, dict[str, float]]] = {}
    for ion, cases in (("carbon", gen.ALL_CASES), ("oxygen", gen.ALL_CASES), ("he", HE_CASES)):
        out[ion] = {}
        for case in cases:
            b = stepping_bias(ion, case)
            if not abs(b["r_model"]) <= MODEL_RESIDUAL_LIMIT:
                raise RecordError(
                    f"stepping model: |R_model| > {MODEL_RESIDUAL_LIMIT} for {ion} case {case}; "
                    "the zero-thickness design premise fails"
                )
            out[ion][case] = b
    return out


# =============================================================================================
# per-run parsing (validity)
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


def measure_outputs(header: str, rows: str, edep_text: str, case_doc: dict[str, Any]) -> dict:
    """<Delta E>_T [MeV] and its SE from one run's outputs; a validity failure raises
    :class:`RunInvalid` (reason without values)."""
    info = case_doc["i3fu2"]
    n = int(case_doc["histories"])
    a, e0 = int(info["a"]), float(info["e0_per_u_mev"])
    e_in = a * e0
    try:
        e = parse_exit_energies(header, rows, n, PDG[info["ion"]], float(info["slab_mm"]))
        ed = edep_per_history(edep_text, n)
    except RecordError as exc:
        raise RunInvalid(str(exc)) from exc
    de = e_in - float(np.mean(e))
    if not 0.0 < de < MAX_LOSS_FRACTION * e_in:
        raise RunInvalid(
            f"mean energy loss outside (0, {MAX_LOSS_FRACTION} E_in); exit energy not the total "
            "kinetic energy?"
        )
    if abs(ed / de - 1.0) > EDEP_RTOL:
        raise RunInvalid(f"energy-deposit cross check fails (|EDep/N / <dE> - 1| > {EDEP_RTOL})")
    return {"de": de, "se": float(np.std(e, ddof=1)) / math.sqrt(n), "scored": n}


def validity_problems(stdout: str, stderr: bytes, histories: int) -> list[str]:
    """Log-level run-validity problems (as ``ss_evaluate.log_checks`` / ``r2_io.run_checks``)."""
    problems = []
    if stderr.strip():
        problems.append("stderr is not empty")
    bad = [ln for ln in stdout.splitlines() if rio._KILLED.search(ln) and not rio._BEGIN.search(ln)]
    if bad:
        problems.append(f"stdout reports killed/stuck tracks ({len(bad)} line(s))")
    found = [int(v) for v in rio._HISTORIES.findall(stdout)]
    if not found or any(v != histories for v in found):
        problems.append(f"histories completed {found} != {histories}")
    return problems


def config_problems(native: str, stdout: str, ion: str) -> list[str]:
    """Configuration markers (a problem means the run did not execute the frozen configuration:
    ``RecordError``, not a validity failure). The fluctuation state is checked separately by the
    discriminating marker set (:func:`marker_hits`)."""
    problems = []
    if rio.registered_processes(stdout) is None:
        problems.append("the Ph/ListProcesses block (Registered Physics Processes) is absent")
    problems += rio.module_check(native, stdout, EM_CONFIG)[0]
    problems += rio.frozen_lines_check(native, EM_CONFIG)
    particle = BEAM_PARTICLE[ion]
    if not re.search(_IONIONI.format(re.escape(particle)), stdout, re.M):
        problems.append(f"no ionIoni block for {particle} in the EM printout")
    return problems


# =============================================================================================
# log markers and the step-count parser
# =============================================================================================


def parse_markers(blob: bytes, where: str) -> dict[str, re.Pattern[str]]:
    """The committed marker file: exactly :data:`MARKER_KEYS`, each a non-empty regular
    expression (compiled with ``re.M``); the step parser's named groups are required."""
    try:
        doc = json.loads(blob.decode("utf-8"))
    except ValueError as exc:
        raise RecordError(f"log-marker file {where} is not JSON: {exc}") from exc
    if not isinstance(doc, dict) or sorted(doc) != sorted(MARKER_KEYS):
        raise RecordError(f"log-marker file {where} must hold exactly the keys {MARKER_KEYS}")
    out: dict[str, re.Pattern[str]] = {}
    for key in MARKER_KEYS:
        v = doc[key]
        if not isinstance(v, str) or not v.strip():
            raise RecordError(f"log-marker file {where}: {key} is not a non-empty regex")
        try:
            out[key] = re.compile(v, re.M)
        except re.error as exc:
            raise RecordError(f"log-marker file {where}: {key} does not compile: {exc}") from exc
    groups = {
        "step_line_re": {"step", "z"},
        "track_header_re": {"track", "parent"},
    }
    for key, need in groups.items():
        if not need <= set(out[key].groupindex):
            raise RecordError(f"log-marker file {where}: {key} lacks the named groups {need}")
    return out


_NUM_RE = re.compile(r"^[+-]?(\d*)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$")


def printed_tolerance(text: str, unit_mm: float) -> float:
    """Half a unit of the last printed digit of a number printed as ``text`` [mm]: e.g. "65.7"
    um -> 0.05 um; "1.25e-02" mm -> 0.005e-2 mm. A printed zero is exact (tolerance 0: Geant4
    prints a vanishing coordinate as a bare 0, whose digit count carries no precision)."""
    m = _NUM_RE.match(text.strip())
    if not m or not (m.group(1) or m.group(2)):
        raise RecordError(f"step line: z {text!r} is not a number")
    if float(text) == 0.0:
        return 0.0
    decimals = len(m.group(2) or "")
    exponent = int(m.group(3) or 0)
    return 0.5 * 10.0 ** (exponent - decimals) * unit_mm


def count_slab_steps(
    stdout: str,
    slab_mm: float,
    step_line_re: re.Pattern[str],
    track_header_re: re.Pattern[str],
) -> dict[int, dict[str, Any]]:
    """Steps per PRIMARY track whose pre-step point lies in the Water slab, from a Geant4
    SteppingVerbose printout (``Ts/TrackingVerbosity = 1``), and whether the (last) counted step
    ends at or beyond the far face.

    Parser contract (the two regexes come from the committed marker file, pinned from the
    step-count smoke logs): ``track_header_re`` matches the header line of every track with the
    named groups ``track`` (track ID) and ``parent`` (parent ID; 0 = primary);
    ``step_line_re`` matches every step line of the step table and nothing else inside a
    primary's block (e.g. not the secondary-list lines), with the named groups ``step`` (step
    number) and ``z`` (the printed z position of the step's post-step point) and optionally
    ``zunit`` (a :data:`LENGTH_UNITS_MM` symbol; absent: mm). Lines are matched one at a time.
    Step 0 is the initial point; the pre-step point of step n >= 1 is the position printed for
    step n - 1. With tol = half a unit of the last digit printed for that z
    (:func:`printed_tolerance`), step n counts iff -tol <= z_pre < slab - tol (world frame, the
    slab spans z = 0 ... slab; a pre-step z within tol of the far face is outside the slab).
    The counted step traverses iff its own printed post-step z >= slab - max(tol_post, 1e-6
    mm). Step numbers of a primary must run 0, 1, 2, ... (else ``RecordError``). Returns
    {ordinal of the primary in the printout: {"steps": n, "traverses": bool}}; Geant4 restarts
    the track IDs in every event, so the track ID itself is not a key. Expected for a one-step
    slab: every one of the primaries has exactly one step in the slab, ending at the far face
    (limiting process Transportation)."""
    out: dict[int, dict[str, Any]] = {}
    current: int | None = None
    prev: tuple[int, float, float] | None = None  # (step, z [mm], tol [mm])
    for line in stdout.splitlines():
        h = track_header_re.search(line)
        if h:
            if int(h.group("parent")) == 0:
                current = len(out)
                out[current] = {"steps": 0, "traverses": False}
            else:
                current = None
            prev = None
            continue
        if current is None:
            continue
        m = step_line_re.search(line)
        if not m:
            continue
        step = int(m.group("step"))
        unit = (m.groupdict().get("zunit") or "mm").strip()
        if unit not in LENGTH_UNITS_MM:
            raise RecordError(f"step line: unknown length unit {unit!r}")
        z_text = m.group("z")
        z = float(z_text) * LENGTH_UNITS_MM[unit]
        tol = printed_tolerance(z_text, LENGTH_UNITS_MM[unit])
        if prev is None:
            if step != 0:
                raise RecordError(f"primary {current}: the first step line is not step 0")
        else:
            if step != prev[0] + 1:
                raise RecordError(f"primary {current}: step numbers are not consecutive")
            z_pre, tol_pre = prev[1], prev[2]
            if -tol_pre - ZEPS_MM <= z_pre < slab_mm - tol_pre - ZEPS_MM:
                out[current]["steps"] += 1
                limit = slab_mm - max(tol, POST_STEP_MIN_TOL_MM) - ZEPS_MM
                out[current]["traverses"] = z >= limit
        prev = (step, z, tol)
    return out


def native_marker_problems(markers: dict[str, re.Pattern[str]], native: str) -> list[str]:
    """Neither fluctuation marker may match the run's native input text: they must identify the
    Geant4 ``G4EmParameters`` printout, not the TOPAS parameter echo."""
    return [
        f"{key} matches the bundle's native input text"
        for key in ("fluct_on_marker", "fluct_off_marker")
        if markers[key].search(native)
    ]


def marker_hits(markers: dict[str, re.Pattern[str]], stdout: str) -> tuple[bool, bool]:
    """(fluct_on_marker matches, fluct_off_marker matches) for one run's stdout."""
    return (
        bool(markers["fluct_on_marker"].search(stdout)),
        bool(markers["fluct_off_marker"].search(stdout)),
    )


# =============================================================================================
# authentication and binding
# =============================================================================================


@dataclass
class Context:
    """Everything bound once per evaluation (see the module docstring)."""

    repo: Path
    reader: Any
    records: str
    smoke_records: str
    outputs: Path
    code_sha: str
    pins: dict[str, dict[str, Any]]
    plan_note: dict[str, str]
    manifest_rel: str
    manifest: dict[str, Any]
    manifest_sha256: str
    markers_rel: str
    markers_sha256: str
    markers: dict[str, re.Pattern[str]]
    sealed: dict[str, str]
    trust: dict[str, Any]
    trust_info: dict[str, Any]
    index: dict[str, str]
    smoke_index: dict[str, str]
    environment: dict[str, str]
    bundles: dict[str, tuple[int, int]] = field(default_factory=dict)
    _committed: dict[tuple[str, str], dict[str, Any]] = field(default_factory=dict)

    @property
    def runs_after(self) -> Any:
        """The later of both pins' ``finished_at`` (evidence and repeat runs start after it)."""
        return max(p["finished"] for p in self.pins.values())


def run_bundles() -> dict[str, tuple[int, int]]:
    """Bundle name -> (case, part) of the 40 evidence and 40 repeat bundles."""
    return {gen.bundle_name(c, p): (c, p) for c, p in gen.evidence_keys() + gen.repeat_keys()}


def manifest_digest(bundles: dict[str, Any]) -> str:
    """``bundle_set_digest`` of the generator (sha256 of the sorted sha256sum-style listing)."""
    lines = [
        f"{entry[f]}  {gen.REL_PREFIX}/{name}/{f}"
        for name, entry in bundles.items()
        for f in ("case.json", "input.txt")
    ]
    listing = "".join(ln + "\n" for ln in sorted(lines, key=lambda s: s.split("  ", 1)[1]))
    return hashlib.sha256(listing.encode()).hexdigest()


def require_clean(repo: Path) -> None:
    """The evaluator worktree has no modified tracked file."""
    proc = subprocess.run(  # noqa: S603 - fixed argv
        ["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=no"],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RecordError(f"git status failed: {proc.stderr.strip()}")
    if proc.stdout.strip():
        raise RecordError("the evaluator worktree is not clean (modified tracked files)")


def executable_closure() -> dict[str, Path]:
    """{path relative to this checkout: loaded file} of every module of the running interpreter
    that belongs to the decision code: this module, every ``ionmc`` module and every script of
    :data:`SCRIPT_MODULES` (namespace packages without a file are skipped). A closure module
    loaded from outside this checkout is a ``RecordError``."""
    rio.attest_module()
    root = REPO_ROOT.resolve()
    out: dict[str, Path] = {}
    for name, mod in sorted(sys.modules.items()):
        if not (
            name == __name__
            or name == "ionmc"
            or name.startswith("ionmc.")
            or name in SCRIPT_MODULES
        ):
            continue
        f = getattr(mod, "__file__", None)
        if not f:
            continue
        path = Path(f).resolve()
        try:
            rel = path.relative_to(root).as_posix()
        except ValueError as exc:
            raise RecordError(f"closure module {name} is loaded from outside {root}") from exc
        out[rel] = path
    return out


def _blob_id(blob: bytes) -> str:
    """git's blob object id of ``blob``."""
    return hashlib.sha1(b"blob %d\0" % len(blob) + blob, usedforsecurity=False).hexdigest()


def tree_blob_ids(repo: Path, sha: str, rels: list[str]) -> dict[str, str]:
    """{path: blob id} of ``rels`` in the tree of commit ``sha`` (one ``git ls-tree`` call)."""
    proc = subprocess.run(  # noqa: S603 - fixed argv
        ["git", "-C", str(repo), "ls-tree", "-r", "--full-tree", sha, "--", *rels],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RecordError(f"git ls-tree {sha}: {proc.stderr.strip()}")
    out = {}
    for ln in proc.stdout.splitlines():
        meta, path = ln.split("\t", 1)
        out[path] = meta.split()[2]
    return out


def sealed_sources(reader: Any, shas: dict[str, str]) -> dict[str, str]:
    """Every file of the executable closure and of :data:`SEALED_FILES` byte-identical on disk
    in ``--repo``, in the loaded file (modules) and as committed at HEAD and at every commit of
    ``shas`` {what: sha}; returns {path: sha256}."""
    closure = executable_closure()
    rels = sorted(set(closure) | set(SEALED_FILES))
    disk: dict[str, bytes] = {}
    for rel in rels:
        path = reader.repo / rel
        if path.is_symlink() or not path.is_file():
            raise RecordError(f"sealed source {rel} is not a regular file in {reader.repo}")
        disk[rel] = path.read_bytes()
        if rel in closure and closure[rel].read_bytes() != disk[rel]:
            raise RecordError(f"sealed source {rel}: the loaded module differs from the worktree")
    for what, sha in {"HEAD": reader.head, **shas}.items():
        ids = tree_blob_ids(reader.repo, sha, rels)
        for rel in rels:
            if ids.get(rel) != _blob_id(disk[rel]):
                raise RecordError(f"sealed source {rel} differs between {what} {sha} and disk")
    return {rel: hashlib.sha256(disk[rel]).hexdigest() for rel in rels}


_LOCK_NUMPY_RE = re.compile(r'^\[\[package\]\]\s*\nname = "numpy"\s*\nversion = "([^"]+)"', re.M)
_LOCK_PYTHON_RE = re.compile(r'^requires-python = ">=(\d+)\.(\d+)', re.M)


def environment_check(lock_text: str) -> dict[str, str]:
    """The running numpy equals the version locked in ``uv.lock`` and Python major.minor equals
    the lower bound of its ``requires-python``; else ``RecordError``."""
    found = _LOCK_NUMPY_RE.findall(lock_text)
    py = _LOCK_PYTHON_RE.search(lock_text)
    if len(found) != 1 or py is None:
        raise RecordError("uv.lock: numpy version or requires-python not found exactly once")
    locked_py = f"{py.group(1)}.{py.group(2)}"
    running_py = f"{sys.version_info[0]}.{sys.version_info[1]}"
    if np.__version__ != found[0]:
        raise RecordError(f"numpy {np.__version__} differs from the uv.lock pin {found[0]}")
    if running_py != locked_py:
        raise RecordError(f"Python {running_py} differs from uv.lock requires-python {locked_py}")
    return {
        "numpy": np.__version__,
        "numpy_locked": found[0],
        "python": running_py,
        "python_locked": locked_py,
        "note": "third-party binary bytes (numpy, the interpreter) are not sealed; uv.lock is",
    }


def check_review_record(doc: Any, what: str) -> dict[str, Any]:
    """A verbatim passed codex-worker review record (schema of ``ss-seal-review.json``)."""
    for key in rio.REVIEW_FIELDS:
        rio.require(doc, key, what)
    rio.require(doc["report"], "verdict", f"{what} report")
    if doc["status"] != "passed":
        raise RecordError(f"{what}: status {doc['status']!r} is not 'passed'")
    if doc["report"]["verdict"] != "pass":
        raise RecordError(f"{what}: report verdict {doc['report']['verdict']!r} is not 'pass'")
    if doc.get("task_id") != TASK_ID:
        raise RecordError(f"{what}: task_id {doc.get('task_id')!r} is not {TASK_ID}")
    if not isinstance(doc["head_sha"], str) or not SHA_RE.match(doc["head_sha"]):
        raise RecordError(f"{what}: head_sha is not a 40-hex commit SHA")
    if "head_sha" in doc["report"] and doc["report"]["head_sha"] != doc["head_sha"]:
        raise RecordError(f"{what}: report head_sha differs from head_sha")
    digest = rio.require(doc, "report_sha256", what)
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise RecordError(f"{what}: report_sha256 is not a 64-hex digest")
    started = rio._parse_time(rio.require(doc, "started_at", what), f"{what} started_at")
    if not started < rio._parse_time(doc["finished_at"], f"{what} finished_at"):
        raise RecordError(f"{what}: started_at is not before finished_at")
    return dict(doc)


def parse_plan_note(text: str, where: str) -> dict[str, str]:
    """The one :data:`PLAN_NOTE_RE` line of the acceptance plan (zero or several: fail closed)."""
    found = list(PLAN_NOTE_RE.finditer(text))
    if len(found) != 1:
        raise RecordError(
            f"{PLAN_PATH} at {where}: expected exactly one I3-FU-1 step-2 review-pin note, "
            f"found {len(found)}"
        )
    return {**found[0].groupdict(), "line": found[0].group(0)}


def load_pins(reader: Any, rel_a: str, rel_b: str, code_sha: str) -> tuple[dict, dict[str, str]]:
    """Both review pins and the plan note naming them (see the module docstring)."""
    note = parse_plan_note(reader.read(PLAN_PATH).decode("utf-8"), "HEAD")
    at_code = parse_plan_note(rio.git_show(reader.repo, code_sha, PLAN_PATH).decode(), code_sha)
    if at_code != note:
        raise RecordError(f"the review-pin plan note differs between code_sha {code_sha} and HEAD")
    pins: dict[str, dict[str, Any]] = {}
    for name, rel, id_key, sha_key in (
        ("evaluator", rel_a, "eval_id", "eval_sha"),
        ("markers", rel_b, "markers_id", "markers_sha"),
    ):
        blob = reader.read(rel)
        if rio.git_show(reader.repo, code_sha, rel) != blob:
            raise RecordError(f"review pin {rel} differs between code_sha {code_sha} and HEAD")
        try:
            doc = json.loads(blob.decode("utf-8"))
        except ValueError as exc:
            raise RecordError(f"review pin {rel} is not JSON: {exc}") from exc
        doc = check_review_record(doc, f"review pin {rel}")
        if doc["review_id"] != note[id_key] or doc["head_sha"] != note[sha_key]:
            raise RecordError(f"review pin {rel} is not the {name} review named by the plan note")
        pins[name] = {
            "path": rel,
            "doc": doc,
            "sha256": hashlib.sha256(blob).hexdigest(),
            "blob": blob,
            "finished": rio._parse_time(doc["finished_at"], f"{rel} finished_at"),
            "started": rio._parse_time(doc["started_at"], f"{rel} started_at"),
        }
    a, b = pins["evaluator"]["doc"]["head_sha"], pins["markers"]["doc"]["head_sha"]
    if not rio.is_ancestor(reader.repo, a, b):
        raise RecordError(f"evaluator review head_sha {a} is not an ancestor of the markers one")
    if not rio.is_ancestor(reader.repo, b, code_sha):
        raise RecordError(f"markers review head_sha {b} is not an ancestor of code_sha {code_sha}")
    if rio.git_show(reader.repo, b, rel_a) != pins["evaluator"]["blob"]:
        raise RecordError(f"review pin {rel_a} differs between the markers review {b} and HEAD")
    return pins, note


def load_manifest(reader: Any, rel: str, shas: dict[str, str]) -> tuple[dict, str]:
    """The frozen manifest, byte-identical at HEAD and at every commit of ``shas``."""
    at_head = reader.read(rel)
    for what, sha in shas.items():
        if rio.git_show(reader.repo, sha, rel) != at_head:
            raise RecordError(f"frozen manifest {rel} differs between {what} {sha} and HEAD")
    try:
        doc = json.loads(at_head.decode("utf-8"))
    except ValueError as exc:
        raise RecordError(f"frozen manifest {rel} is not JSON: {exc}") from exc
    if not isinstance(doc, dict) or doc.get("schema") != gen.MANIFEST_SCHEMA:
        raise RecordError(f"frozen manifest schema is not {gen.MANIFEST_SCHEMA}")
    bundles = doc.get("bundles")
    if not isinstance(bundles, dict) or sorted(bundles) != gen.bundle_names():
        raise RecordError("frozen manifest bundle set differs from the generator's")
    if manifest_digest(bundles) != doc.get("bundle_set_digest"):
        raise RecordError("frozen manifest entries do not reproduce its bundle_set_digest")
    return doc, hashlib.sha256(at_head).hexdigest()


def load_index(reader: Any, records: str) -> dict[str, str]:
    """bundle -> run id from ``<records>/index.json`` (evidence and repeat bundles only)."""
    idx = reader.json(f"{records}/index.json", "i3fu2 run index")
    runs = rio.require(idx, "runs", "i3fu2 run index")
    if not isinstance(runs, list):
        raise RecordError("i3fu2 run index: runs must be a list")
    known = run_bundles()
    out: dict[str, str] = {}
    for e in runs:
        for key in ("case", "part", "run_id", "bundle"):
            rio.require(e, key, "i3fu2 run index entry")
        bundle, run_id = str(e["bundle"]), str(e["run_id"])
        if bundle not in known or known[bundle] != (e["case"], e["part"]):
            raise RecordError(f"index entry {bundle} ({e['case']}, {e['part']}) is not i3fu2")
        if bundle in out:
            raise RecordError(f"index lists bundle {bundle} twice")
        if run_id in out.values():
            raise RecordError(f"index lists run {run_id} twice")
        out[bundle] = run_id
    return out


def load_smoke_index(reader: Any, records: str) -> dict[str, str]:
    """smoke bundle -> run id from ``<smokes>/index.json``: exactly the four smokes."""
    idx = reader.json(f"{records}/index.json", "i3fu2 smoke index")
    runs = rio.require(idx, "runs", "i3fu2 smoke index")
    if not isinstance(runs, list):
        raise RecordError("i3fu2 smoke index: runs must be a list")
    out: dict[str, str] = {}
    for e in runs:
        bundle = str(rio.require(e, "bundle", "smoke index entry"))
        run_id = str(rio.require(e, "run_id", "smoke index entry"))
        if bundle not in ALL_SMOKES:
            raise RecordError(f"smoke index entry {bundle} is not an i3fu2 smoke")
        if bundle in out:
            raise RecordError(f"smoke index lists {bundle} twice")
        out[bundle] = run_id
    missing = sorted(set(ALL_SMOKES) - set(out))
    if missing:
        raise RecordError(f"smoke index lacks {missing} (all four smokes are required)")
    return out


def context(
    repo: Path | str,
    records: str,
    outputs: Path | str,
    review_pin: str,
    manifest: str,
    code_sha: str,
    *,
    smokes: str,
    markers_pin: str = MARKERS_PIN,
    markers: str = MARKERS,
) -> Context:
    reader = rio.CommittedReader(repo)
    require_clean(reader.repo)
    if not SHA_RE.match(code_sha):
        raise RecordError(f"--code-sha {code_sha!r} is not a 40-hex commit SHA")
    if not rio.is_ancestor(reader.repo, code_sha, reader.head):
        raise RecordError(f"code_sha {code_sha} is not an ancestor of HEAD {reader.head}")
    rels = [reader.rel(p) for p in (records, smokes, review_pin, markers_pin, manifest)]
    records_rel, smoke_rel, pin_a, pin_b, man_rel = rels
    pins, note = load_pins(reader, pin_a, pin_b, code_sha)
    reviewed = pins["evaluator"]["doc"]["head_sha"]
    markers_reviewed = pins["markers"]["doc"]["head_sha"]
    seal_shas = {"code_sha": code_sha, "the evaluator review head_sha": reviewed}
    sealed = sealed_sources(reader, seal_shas)
    man, man_sha = load_manifest(reader, man_rel, seal_shas)
    markers_rel = reader.rel(markers)
    m_blob = reader.read(markers_rel)
    for what, sha in (("code_sha", code_sha), ("the markers review head_sha", markers_reviewed)):
        if rio.git_show(reader.repo, sha, markers_rel) != m_blob:
            raise RecordError(
                f"log-marker file {markers_rel} differs between {what} {sha} and HEAD"
            )
    trust, trust_info = rio.committed_trust(reader)
    environment = environment_check(reader.read("uv.lock").decode("utf-8"))
    return Context(
        repo=reader.repo,
        reader=reader,
        records=records_rel,
        smoke_records=smoke_rel,
        outputs=Path(outputs),
        code_sha=code_sha,
        pins=pins,
        plan_note=note,
        manifest_rel=man_rel,
        manifest=man,
        manifest_sha256=man_sha,
        markers_rel=markers_rel,
        markers_sha256=hashlib.sha256(m_blob).hexdigest(),
        markers=parse_markers(m_blob, markers_rel),
        sealed=sealed,
        trust=trust,
        trust_info=trust_info,
        index=load_index(reader, records_rel),
        smoke_index=load_smoke_index(reader, smoke_rel),
        environment=environment,
        bundles=run_bundles(),
    )


def committed_bundle(ctx: Context, bundle: str, sha: str) -> dict[str, Any]:
    """{file: {"sha256", "bytes"}} of the bundle committed at ``sha`` (cached)."""
    key = (sha, bundle)
    if key not in ctx._committed:
        bdir = f"{gen.REL_PREFIX}/{bundle}"
        files: dict[str, Any] = {}
        blobs: dict[str, bytes] = {}
        for rel in rio.git_ls(ctx.repo, sha, bdir):
            blob = rio.git_show(ctx.repo, sha, rel)
            name = rel[len(bdir) + 1 :]
            files[name] = {"sha256": hashlib.sha256(blob).hexdigest(), "bytes": len(blob)}
            blobs[name] = blob
        if not files:
            raise RecordError(f"bundle {bundle} is not committed at {sha}")
        ctx._committed[key] = {"files": files, "blobs": blobs}
    return ctx._committed[key]


@dataclass(frozen=True)
class Expect:
    """What a run's bundle must be: (case, part, kind, seed) and whether it is a step smoke."""

    case: int
    part: int
    kind: str
    seed: int
    smoke: bool = False


def expect_for(ctx: Context, bundle: str) -> Expect:
    if bundle in ALL_SMOKES:
        seed, part, case, _hist, _tracking = gen.SMOKES[bundle]
        return Expect(case, part, "smoke", seed, smoke=True)
    case, part = ctx.bundles[bundle]
    kind = "repeat" if part in gen.REPEAT_PARTS else "evidence"
    return Expect(case, part, kind, gen.topas_seed(case, part))


def check_smoke_sha(ctx: Context, sha: str) -> None:
    """A smoke ran a commit that carries pin A (pin A byte-identical there and pin A's
    ``head_sha`` an ancestor of it), that is an ancestor of (or equal to) pin B's ``head_sha``,
    and at which the sealed closure and the manifest are byte-identical."""
    pin_a = ctx.pins["evaluator"]
    a_head, b_head = pin_a["doc"]["head_sha"], ctx.pins["markers"]["doc"]["head_sha"]
    if not SHA_RE.match(sha) or not rio.is_ancestor(ctx.repo, a_head, sha):
        raise RecordError(f"smoke code_sha {sha} is not a descendant of the evaluator review")
    if rio.git_show(ctx.repo, sha, pin_a["path"]) != pin_a["blob"]:
        raise RecordError(f"smoke code_sha {sha} does not carry the evaluator review pin")
    if not rio.is_ancestor(ctx.repo, sha, b_head):
        raise RecordError(f"smoke code_sha {sha} is not an ancestor of the markers review")
    rels = sorted(ctx.sealed)
    ids = tree_blob_ids(ctx.repo, sha, rels)
    for rel in rels:
        blob = rio.git_show(ctx.repo, sha, rel) if rel in ids else b""
        if rel not in ids or hashlib.sha256(blob).hexdigest() != ctx.sealed[rel]:
            raise RecordError(f"sealed source {rel} differs at the smoke code_sha {sha}")
    if hashlib.sha256(rio.git_show(ctx.repo, sha, ctx.manifest_rel)).hexdigest() != (
        ctx.manifest_sha256
    ):
        raise RecordError(f"frozen manifest differs at the smoke code_sha {sha}")


def load_run(ctx: Context, bundle: str, run_id: str, records: str | None = None) -> Any:
    """Authenticate and bind one signed reference-runner run of ``bundle`` (``r2_io.TopasRun``);
    every failure raises RecordError."""
    exp = expect_for(ctx, bundle)
    rdir = f"{records or ctx.records}/{run_id}"
    record_path = f"{rdir}/record.json"
    result = ctx.reader.json(record_path, f"signed run record {run_id}")
    status, why = rio.attest_module().classify_record(result, ctx.trust, rio.REFERENCE_SIGNER)
    if status != "verified":
        raise RecordError(f"record.json attestation {status}: {why}")
    for key in rio.RESULT_FIELDS:
        rio.require(result, key, f"result record {run_id}")
    if result["schema_version"] != 3 or result.get("timed_out") is not False:
        raise RecordError("schema_version/timed_out not 3/false")
    if (
        not re.fullmatch(
            r"REF-" + re.escape(str(result["configuration_id"])[:20]) + r"-[0-9a-f]{8}", run_id
        )
        or result["run_id"] != run_id
    ):
        raise RecordError("run id not bound to the record/configuration")
    if result["engine"] != "topas" or not (
        result["succeeded"] is True and result["exit_code"] == 0 and result["missing_outputs"] == []
    ):
        raise RecordError(
            f"not a successful topas run (succeeded={result['succeeded']!r}, exit_code="
            f"{result['exit_code']!r}, missing_outputs={result['missing_outputs']!r})"
        )
    code_sha = str(result["code_sha"])
    if exp.smoke:
        check_smoke_sha(ctx, code_sha)
        after, after_what = ctx.pins["evaluator"]["finished"], "the evaluator review pin's"
    else:
        if code_sha != ctx.code_sha:
            raise RecordError(f"record code_sha {code_sha} != --code-sha {ctx.code_sha}")
        after, after_what = ctx.runs_after, "both review pins'"
    started = rio._parse_time(result["started_at"], f"{run_id} started_at")
    finished = rio._parse_time(result["finished_at"], f"{run_id} finished_at")
    if started > finished:
        raise RecordError("started_at is later than finished_at")
    if not started > after:
        raise RecordError(
            f"reference_started {result['started_at']} is not after {after_what} finished_at"
        )
    if exp.smoke and not finished < ctx.pins["markers"]["started"]:
        raise RecordError(
            f"smoke finished_at {result['finished_at']} is not before the markers review's "
            "started_at (the marker file is derived from the smoke logs)"
        )
    sizes: dict[str, Any] = {}
    listing = rio._listing(
        ctx.reader.json(f"{rdir}/artifacts.json", "artifact listing"), run_id, code_sha, sizes
    )
    work = {}
    for p, d in listing.items():
        if p.startswith("work/"):
            b = sizes.get(p)
            if not isinstance(b, int) or isinstance(b, bool) or b < 0:
                raise RecordError(f"artifact listing lacks the byte size of {p}")
            work[p[len("work/") :]] = {"sha256": d, "bytes": b}
    if rio.canonical_digest(work) != result["output_files_sha256"] or result["output_count"] != len(
        work
    ):
        raise RecordError("listed work/ outputs do not reproduce the signed digest")
    for name in ("request.json", "stdout.txt", "stderr.txt", "record.json"):
        if name not in listing:
            raise RecordError(f"artifact listing lacks {name}")
    if ctx.reader.reads[record_path] != listing["record.json"]:
        raise RecordError("record.json differs from the listed record.json")
    for name, key in (("stdout.txt", "stdout_sha256"), ("stderr.txt", "stderr_sha256")):
        if result[key] != listing[name]:
            raise RecordError(f"signed {key} differs from the listing")
    run = rio.TopasRun(
        run_id=run_id,
        case=exp.case,
        config="i3fu2",
        part=exp.part,
        bundle=bundle,
        configuration_id=str(result["configuration_id"]),
        code_sha=code_sha,
        finished_at=str(result["finished_at"]),
        seconds=float(result["seconds"]),
        reference_started=str(result["started_at"]),
        histories=0,
        expected_outputs=[],
        listing=listing,
        outputs_dir=ctx.outputs / run_id,
    )
    request_blob = rio.read_verified(run, "request.json")
    if hashlib.sha256(request_blob).hexdigest() != result["request_sha256"]:
        raise RecordError("request.json does not hash to the signed request_sha256")
    try:
        request = json.loads(request_blob.decode("utf-8"))
    except ValueError as exc:
        raise RecordError(f"request.json is not JSON: {exc}") from exc
    if rio.canonical_digest(request) != run.configuration_id or request.get("code_sha") != code_sha:
        raise RecordError("request.json not bound to the record")
    committed = committed_bundle(ctx, bundle, code_sha)
    if request.get("input_files") != committed["files"]:
        raise RecordError(f"request input_files differ from the bundle {bundle} at code_sha")
    entry = ctx.manifest["bundles"].get(bundle)
    for name in ("case.json", "input.txt"):
        if not isinstance(entry, dict) or committed["files"].get(name, {}).get(
            "sha256"
        ) != entry.get(name):
            raise RecordError(f"{bundle}/{name} differs from the frozen manifest")
    case_doc = json.loads(committed["blobs"]["case.json"].decode("utf-8"))
    if request.get("case") != case_doc:
        raise RecordError("request case differs from the committed case.json")
    info = case_doc.get("i3fu2", {})
    if (info.get("case"), info.get("part"), info.get("kind")) != (exp.case, exp.part, exp.kind):
        raise RecordError(f"case.json i3fu2 block does not describe {exp.kind} ({exp.case}, "
                          f"{exp.part})")  # fmt: skip
    if case_doc.get("seeds") != [exp.seed]:
        raise RecordError(f"seed {case_doc.get('seeds')} breaks the formula")
    hist = case_doc.get("histories")
    if not isinstance(hist, int) or isinstance(hist, bool) or hist <= 0:
        raise RecordError(f"histories {hist!r} is not a positive integer")
    run.histories = hist
    run.case_doc, run.request = case_doc, request
    run.native = committed["blobs"][case_doc["input"]].decode("utf-8", "replace")
    return run


def check_run(ctx: Context, bundle: str, run_id: str) -> dict[str, Any]:
    """Authenticate one run (RecordError naming the run), then measure it; a validity failure
    gives ``valid`` False with its ``reason``. ``markers`` = (fluct_on hit, fluct_off hit)."""
    try:
        run = load_run(ctx, bundle, run_id)
        info = run.case_doc["i3fu2"]
        stdout = rio.read_verified(run, "stdout.txt").decode("utf-8", "replace")
        stderr = rio.read_verified(run, "stderr.txt")
        problems = config_problems(run.native, stdout, info["ion"])
        problems += native_marker_problems(ctx.markers, run.native)
        if problems:
            raise RecordError("configuration markers: " + "; ".join(problems))
        files = {n: rio.read_verified(run, f"work/{n}") for n in gen.EXPECTED_OUTPUTS}
    except RecordError as exc:
        raise RecordError(f"run {bundle} ({run_id}): {exc}") from exc
    out = {
        "bundle": bundle,
        "run_id": run_id,
        "info": info,
        "seed": run.case_doc["seeds"][0],
        "histories": run.histories,
        "valid": True,
        "reason": None,
        "scored": None,
        "markers": marker_hits(ctx.markers, stdout),
    }
    try:
        bad = validity_problems(stdout, stderr, run.histories)
        if bad:
            raise RunInvalid("; ".join(bad))
        header, rows, edep = (files[n].decode("utf-8", "replace") for n in gen.EXPECTED_OUTPUTS)
        m = measure_outputs(header, rows, edep, run.case_doc)
    except RunInvalid as exc:
        out.update(valid=False, reason=str(exc))
        return out
    out.update(m)
    return out


def check_smoke(ctx: Context, bundle: str, run_id: str) -> dict[str, Any]:
    """Authenticate one smoke (RecordError naming it); for a step-count smoke also the one-step
    premise: :func:`count_slab_steps` and the far-face phase space (the evidence-run checks;
    a failure there fails the premise, it is not a ``RecordError``)."""
    try:
        run = load_run(ctx, bundle, run_id, ctx.smoke_records)
        info = run.case_doc["i3fu2"]
        stdout = rio.read_verified(run, "stdout.txt").decode("utf-8", "replace")
        problems = config_problems(run.native, stdout, info["ion"])
        problems += native_marker_problems(ctx.markers, run.native)
        if problems:
            raise RecordError("configuration markers: " + "; ".join(problems))
        out: dict[str, Any] = {
            "fluctuations": bool(info["fluctuations"]),
            "markers": marker_hits(ctx.markers, stdout),
            "run_id": run_id,
            "seed": run.case_doc["seeds"][0],
            "code_sha": run.code_sha,
        }
        if bundle not in STEP_SMOKES:
            return out
        slab = float(info["slab_mm"])
        counts = count_slab_steps(
            stdout, slab, ctx.markers["step_line_re"], ctx.markers["track_header_re"]
        )
        header, rows = (
            rio.read_verified(run, f"work/{n}").decode("utf-8", "replace")
            for n in ("ss_exit.header", "ss_exit.phsp")
        )
    except RecordError as exc:
        raise RecordError(f"smoke {bundle} ({run_id}): {exc}") from exc
    try:
        parse_exit_energies(header, rows, run.histories, PDG[info["ion"]], slab)
        complete = True
    except RecordError:
        complete = False
    traversing = sum(1 for c in counts.values() if c["steps"] == 1 and c["traverses"])
    out.update(
        primaries=len(counts),
        max_steps_in_slab=max((c["steps"] for c in counts.values()), default=None),
        traversing_primaries=traversing,
        phase_space_complete=complete,
        ok=len(counts) == run.histories and traversing == run.histories and complete,
    )
    return out


# =============================================================================================
# points, repeats and aggregation
# =============================================================================================


def classify(r_m: float, r_i: float, s: float, sigma: float) -> str:
    """Decision rule 2a for one band, relative to the sign of s_b; U checked first (X_b > 0
    implies s_b != 0 for the other classes)."""
    x = min(X_CAP, abs(s) / 3.0)
    k = K_ONE_SIDED * sigma
    if x <= k:
        return "U"
    sign = 1.0 if s > 0 else -1.0
    if abs(r_i) + k <= x:
        return "I"
    if sign * r_i - k > x:
        return "H"
    if abs(r_m) + k <= x:
        return "M"
    if -sign * r_m - k > x:
        return "L"
    return "N"


def band_of(e0: float) -> str:
    for name, lo, hi, closed in BANDS:
        if lo <= e0 < hi or (closed and e0 == hi):
            return name
    raise ValueError(f"E0 {e0} outside the bands")


def richardson(r1: float, r2: float, f1: float = FACTORS[0], f2: float = FACTORS[1]) -> float:
    """Zero-thickness extrapolation of a quantity linear in the thickness factor f, from its
    values ``r1`` at ``f1`` and ``r2`` at ``f2``: (f2 r1 - f1 r2) / (f2 - f1)."""
    return (f2 * r1 - f1 * r2) / (f2 - f1)


def richardson_se(se1: float, se2: float, f1: float = FACTORS[0], f2: float = FACTORS[1]) -> float:
    """SE of :func:`richardson` for independent runs."""
    return math.hypot(f2 * se1, f1 * se2) / (f2 - f1)


def run_ratios(run: dict[str, Any], s_model: dict, s_table: dict, keys: dict) -> dict[str, Any]:
    """r_M (and r_I, s for C/O), SE_rel of one valid run (predictions memoised)."""
    info = run["info"]
    ion, a, e0, t = info["ion"], int(info["a"]), float(info["e0_per_u_mev"]), info["slab_mm"]
    try:
        de_m = predicted_loss(f"model:{ion}", s_model[ion], a, e0, float(t))
        de_i = predicted_loss(keys[ion], s_table[ion], a, e0, float(t)) if ion in s_table else None
    except ValueError as exc:
        raise RecordError(f"run {run['bundle']}: prediction failed: {exc}") from exc
    out = {"e0": e0, "r_M": run["de"] / de_m - 1.0, "se_rel": run["se"] / run["de"]}
    if de_i is not None:
        out["r_I"] = run["de"] / de_i - 1.0
        out["s"] = de_i / de_m - 1.0
    return out


@dataclass
class Point:
    """Outcome of the repeat rule for one point: ``runs`` the final attempt (with ratios where
    valid; None if an original is missing), ``valid``, ``reason``, the repeats used and the
    predicted thickness difference (pairs)."""

    runs: list[dict[str, Any]] | None
    valid: bool
    reason: str | None
    used: dict[str, str]
    d_pred: float | None = None

    def thickness_difference(self) -> float | None:
        """Delta_meas = r_M(f2) - r_M(f1) of the final attempt if both of its runs are valid."""
        if self.runs is None or len(self.runs) != 2 or not all(r["valid"] for r in self.runs):
            return None
        return float(self.runs[1]["r_M"] - self.runs[0]["r_M"])

    def consistency_residual(self) -> float | None:
        """Delta_meas - Delta_pred of the final attempt (pairs with both runs valid)."""
        d = self.thickness_difference()
        return None if d is None or self.d_pred is None else d - self.d_pred


def resolve_point(
    runs: dict[str, dict[str, Any]],
    case: int,
    parts: tuple[int, ...],
    ratios: Callable,
    d_pred: float | None = None,
) -> Point:
    """Deterministic repeat rule for one point (one slab part, or the t / 1.6 t pair with its
    predicted difference ``d_pred``); see the module docstring."""
    if len(parts) == 2 and d_pred is None:
        raise ValueError("a pair needs its predicted thickness difference")
    names = [gen.bundle_name(case, p) for p in parts]
    reps = [gen.bundle_name(case, p + gen.REPEAT_OFFSET) for p in parts]

    def consistent(rs: list[dict[str, Any]]) -> bool:
        if len(rs) == 1:
            return True
        d_meas = float(rs[1]["r_M"]) - float(rs[0]["r_M"])
        return abs(d_meas - float(d_pred or 0.0)) <= CONSISTENCY_TOL

    def with_ratios(r: dict[str, Any]) -> dict[str, Any]:
        return {**r, **ratios(r)} if r["valid"] else r

    if any(n not in runs for n in names):
        for n, rep in zip(names, reps, strict=True):
            if n not in runs and rep in runs:
                raise RecordError(f"repeat {rep} offered but the original {n} is missing")
        return Point(None, False, "missing run", {}, d_pred)
    first = [with_ratios(runs[n]) for n in names]
    failed = {i: str(r["reason"]) for i, r in enumerate(first) if not r["valid"]}
    by_consistency = not failed and not consistent(first)
    if by_consistency:
        failed = dict.fromkeys(range(len(parts)), CONSISTENCY_REASON)
    offered = {i for i, rep in enumerate(reps) if rep in runs}
    if not failed:
        if offered:
            raise RecordError(
                f"repeat {reps[min(offered)]} offered without a recorded failure of the point"
            )
        return Point(first, True, None, {}, d_pred)
    for i in sorted(offered - set(failed)):
        raise RecordError(f"repeat {reps[i]} offered but the original {names[i]} is valid")
    why = CONSISTENCY_REASON if by_consistency else "invalid run"
    if not offered:
        return Point(first, False, f"{why}; no repeat", {}, d_pred)
    if offered != set(failed):
        missing = [reps[i] for i in sorted(set(failed) - offered)]
        raise RecordError(f"partial repeat set: {missing} required with the offered repeats")
    chosen, used = list(first), {}
    for i in sorted(offered):
        chosen[i] = with_ratios(runs[reps[i]])
        used[reps[i]] = f"original {names[i]} invalid: {failed[i]}"
    if any(not r["valid"] for r in chosen):
        return Point(
            chosen, False, "repeat run invalid (second failure of the point)", used, d_pred
        )
    if not consistent(chosen):
        return Point(
            chosen, False, f"{CONSISTENCY_REASON} (second failure of the point)", used, d_pred
        )
    return Point(chosen, True, None, used, d_pred)


def _max_abs(values: list[float | None]) -> float | None:
    vals = [abs(v) for v in values if v is not None]
    return max(vals) if vals else None


def band_classes(points: list[dict[str, Any]], systematic: float) -> dict[str, str]:
    """Class per band of the valid points with the method systematic ``systematic`` (the
    robustness check uses twice :data:`METHOD_SYSTEMATIC`); "none" for an empty band."""
    out = {}
    for name, *_ in BANDS:
        valid = [p for p in points if p["valid"] and band_of(p["e0"]) == name]
        if not valid:
            out[name] = "none"
            continue
        r_m = float(np.mean([p["r_M"] for p in valid]))
        r_i = float(np.mean([p["r_I"] for p in valid]))
        s = float(np.mean([p["s"] for p in valid]))
        se_b = math.sqrt(sum(p["se_rel"] ** 2 for p in valid)) / len(valid)
        out[name] = classify(r_m, r_i, s, math.sqrt(se_b**2 + systematic**2))
    return out


def band_aggregates(points: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Per band: n, means of r_M, r_I, s (extrapolated points; written rounded, classified at
    full precision), SE_b, sigma_b, X_b, class and the band maxima of the thickness difference,
    the consistency residual, the predicted difference and the model residual."""
    out: dict[str, dict[str, Any]] = {}
    for name, *_ in BANDS:
        pts = [p for p in points if band_of(p["e0"]) == name]
        valid = [p for p in pts if p["valid"]]
        n = len(valid)
        band: dict[str, Any] = {
            "n": n,
            "mean_r_M": None,
            "mean_r_I": None,
            "mean_s": None,
            "se_b": None,
            "sigma_b": None,
            "x_b": None,
            "class": "none",
            "max_abs_thickness_difference": _max_abs([p["d_thick"] for p in pts]),
            "max_abs_consistency_residual": _max_abs([p["d_resid"] for p in pts]),
            "max_abs_predicted_difference": _max_abs([p["d_pred"] for p in pts]),
            "max_abs_model_residual": _max_abs([p["r_model"] for p in pts]),
        }
        if n:
            r_m = float(np.mean([p["r_M"] for p in valid]))
            r_i = float(np.mean([p["r_I"] for p in valid]))
            s = float(np.mean([p["s"] for p in valid]))
            se_b = math.sqrt(sum(p["se_rel"] ** 2 for p in valid)) / n
            sigma = math.sqrt(se_b**2 + METHOD_SYSTEMATIC**2)
            band.update(
                mean_r_M=_round(r_m),
                mean_r_I=_round(r_i),
                mean_s=_round(s),
                se_b=se_b,
                sigma_b=sigma,
                x_b=_round(min(X_CAP, abs(s) / 3.0)),
                **{"class": classify(r_m, r_i, s, sigma)},
            )
        out[name] = band
    return out


def pattern(bands: dict[str, dict[str, dict[str, Any]]]) -> str:
    """ICRU / MODEL / MIXED over both ions (revision 1, symmetric)."""

    def like(ion: str, allowed: tuple[str, ...]) -> bool:
        cls = {b: str(bands[ion][b]["class"]) for b in ("[2,5)", "[5,10)", "[10,20)")}
        return all(cls[b] in allowed for b in LOW_BANDS) and cls["[10,20)"] in (*allowed, "U")

    ions = ("carbon", "oxygen")
    if all(like(i, ICRU_LIKE) for i in ions):
        return "ICRU"
    if all(like(i, MODEL_LIKE) for i in ions):
        return "MODEL"
    return "MIXED"


def markers_discriminate(hits: list[tuple[bool, bool, bool]]) -> bool:
    """[(fluctuations on, on-marker hit, off-marker hit)] over all authenticated runs: each
    run's own marker matches and the other does not; at least one run of each state."""
    states = {f for f, _on, _off in hits}
    ok = all((on and not off) if f else (off and not on) for f, on, off in hits)
    return ok and states == {True, False}


# =============================================================================================
# evaluation
# =============================================================================================


def evaluate(
    repo: Path | str,
    records: str,
    runs_dir: Path | str,
    review_pin: str,
    manifest: str,
    archive: Path,
    code_sha: str,
    *,
    smokes: str,
    markers_pin: str = MARKERS_PIN,
    markers: str = MARKERS,
    archive_sha256: str | None = None,
    member_sha256: dict[str, str] | None = None,
) -> dict[str, Any]:
    """The result document (aggregates only). ``archive_sha256``/``member_sha256`` override the
    pinned libdEdx hashes of ``compare_icru73`` (tests with a synthetic archive only)."""
    ctx = context(
        repo,
        records,
        runs_dir,
        review_pin,
        manifest,
        code_sha,
        smokes=smokes,
        markers_pin=markers_pin,
        markers=markers,
    )
    stepping = stepping_model_check()
    archive_sha = icru.ARCHIVE_SHA256 if archive_sha256 is None else archive_sha256
    members = dict(icru.MEMBER_SHA256 if member_sha256 is None else member_sha256)
    tables = icru.load_icru73_water(
        archive, archive_sha256=archive_sha, member_sha256=members, ions=(6, 8)
    )
    s_model = {ion: model_stopping(gen.IONS[ion][0]) for ion in gen.IONS}
    s_table = {
        ion: table_stopping(tables[z].energy_per_u_mev, tables[z].stopping_mev_cm2_g)
        for ion, z in (("carbon", 6), ("oxygen", 8))
    }
    keys = {ion: f"table:{archive_sha}:{ion}" for ion in s_table}

    runs = {b: check_run(ctx, b, rid) for b, rid in sorted(ctx.index.items())}
    for name in runs:
        c, p = ctx.bundles[name]
        if p in gen.REPEAT_PARTS and gen.bundle_name(c, p - gen.REPEAT_OFFSET) not in runs:
            raise RecordError(f"repeat {name} offered but its original run is not indexed")
        info = runs[name]["info"]
        mult = 2 if gen.PARTS[gen.base_part(p)][1] == 2 else 1
        if float(info["slab_mm"]) != gen.slab_mm(info["ion"], c, mult):
            raise RecordError(f"run {name}: slab thickness differs from the generator's")

    smoke_runs = {b: check_smoke(ctx, b, rid) for b, rid in sorted(ctx.smoke_index.items())}
    premise_keys = ("primaries", "max_steps_in_slab", "traversing_primaries")
    premise = {
        "checked": True,  # the smoke index holds all four smokes (else RecordError)
        "smokes": {
            b: {
                **{k: smoke_runs[b][k] for k in premise_keys},
                "phase_space_complete": smoke_runs[b]["phase_space_complete"],
            }
            for b in STEP_SMOKES
        },  # fmt: skip
        "pass": all(smoke_runs[b]["ok"] for b in STEP_SMOKES),
    }

    def ratios(r: dict[str, Any]) -> dict[str, Any]:
        return run_ratios(r, s_model, s_table, keys)

    repeats_used: dict[str, str] = {}
    invalid_points: dict[str, str] = {}
    points: dict[str, list[dict[str, Any]]] = {}
    chosen_t: dict[tuple[str, int], dict[str, Any]] = {}
    for ion, parts in PAIR_PARTS.items():
        points[ion] = []
        for case in gen.PARTS[parts[0]][3]:
            for p in parts:
                if gen.bundle_name(case, p) not in runs:
                    raise RecordError(f"missing evidence run {gen.bundle_name(case, p)}")
            model = stepping[ion][case]
            pt = resolve_point(runs, case, parts, ratios, model["d_pred"])
            repeats_used.update(pt.used)
            e0 = float(gen.CASES[case])
            base = {
                "e0": e0,
                "d_thick": pt.thickness_difference(),
                "d_resid": pt.consistency_residual(),
                "d_pred": model["d_pred"],
                "r_model": model["r_model"],
            }
            if not pt.valid or pt.runs is None:
                invalid_points[f"{ion} case {case}"] = str(pt.reason)
                points[ion].append({**base, "valid": False})
                continue
            t, t2 = pt.runs
            f1, f2 = (float(r["info"]["thickness_factor"]) for r in (t, t2))
            if (f1, f2) != FACTORS:
                raise RecordError(f"{ion} case {case}: thickness factors {(f1, f2)} != {FACTORS}")
            chosen_t[(ion, case)] = t
            points[ion].append(
                {
                    **base,
                    "valid": True,
                    "r_M": richardson(t["r_M"], t2["r_M"], f1, f2),
                    "r_I": richardson(t["r_I"], t2["r_I"], f1, f2),
                    "s": richardson(t["s"], t2["s"], f1, f2),
                    "se_rel": richardson_se(t["se_rel"], t2["se_rel"], f1, f2),
                }
            )
    bands = {ion: band_aggregates(points[ion]) for ion in PAIR_PARTS}

    control_invalid: dict[str, str] = {}
    he_r: list[float] = []
    for case in HE_CASES:
        pt = resolve_point(runs, case, (HE_PART,), ratios)
        repeats_used.update(pt.used)
        if not pt.valid or pt.runs is None:
            control_invalid[f"he case {case}"] = str(pt.reason)
        else:
            he_r.append(pt.runs[0]["r_M"])
    max_he = max((abs(v) for v in he_r), default=None)
    helium = {
        "pass": len(he_r) == len(HE_CASES) and all(abs(v) <= HE_TOL for v in he_r),
        "max_abs_r_m": _round(max_he),
    }
    fl_d: list[float] = []
    for case in FLUCT_CASES:
        pt = resolve_point(runs, case, (FLUCT_PART,), ratios)
        repeats_used.update(pt.used)
        if not pt.valid or pt.runs is None:
            control_invalid[f"carbon fluctuations-off case {case}"] = str(pt.reason)
        elif ("carbon", case) not in chosen_t:
            control_invalid[f"carbon fluctuations-off case {case}"] = "no valid slab-t C point"
        else:
            fl_d.append(pt.runs[0]["r_M"] - chosen_t[("carbon", case)]["r_M"])
    hits = [(bool(r["info"]["fluctuations"]), *r["markers"]) for r in runs.values()]
    hits += [(bool(v["fluctuations"]), *v["markers"]) for v in smoke_runs.values()]
    discriminate = markers_discriminate(hits)
    fluct = {
        "pass": discriminate
        and len(fl_d) == len(FLUCT_CASES)
        and all(abs(v) <= FLUCT_TOL for v in fl_d),
        "max_abs_difference": max((abs(v) for v in fl_d), default=None),
        "markers_discriminate": discriminate,
    }
    controls_pass = bool(helium["pass"] and fluct["pass"])
    low_u = any(bands[i][b]["class"] == "U" for i in bands for b in LOW_BANDS)
    reasons = []
    if not controls_pass:
        reasons.append("a control failed or is incomplete")
    if invalid_points:
        reasons.append("an evidence point is invalid without a valid repeat")
    if low_u:
        reasons.append("class U in [2,5) or [5,10)")
    if not premise["pass"]:
        reasons.append(PREMISE_REASON)
    before = pattern(bands)
    classes2 = {ion: band_classes(points[ion], 2.0 * METHOD_SYSTEMATIC) for ion in PAIR_PARTS}
    pattern2 = pattern(
        {ion: {b: {"class": c} for b, c in cl.items()} for ion, cl in classes2.items()}
    )
    robust = before == pattern2
    if not robust:
        reasons.append(ROBUST_REASON)
    conclusive = not reasons
    late = sorted(set(executable_closure()) - set(ctx.sealed))
    if late:
        raise RecordError(f"modules imported after the seal: {late}")
    try:
        ionmc_version: str | None = importlib.metadata.version("ionmc")
    except importlib.metadata.PackageNotFoundError:
        ionmc_version = None
    he_model = _max_abs([stepping["he"][c]["r_model"] for c in HE_CASES])
    return {
        "schema": RESULT_SCHEMA,
        "task_id": TASK_ID,
        "amendment": "I3-FU-1 step 2 (revision 2)",
        "claim": "diagnostic; no gate; no validation claim",
        "code_sha": code_sha,
        "method": {
            "quantity": "<Delta E>_T = A E0 - mean(E_exit); predictions RK4 over the slab",
            "point": "zero-thickness extrapolation (f2 r1 - f1 r2) / (f2 - f1) of r_M, r_I, s",
            "thickness_factors": list(FACTORS),
            "rk4_steps": N_RK4,
            "model": "ionmc bethe_mass_stopping default options, water(78.0)",
            "table": "libdEdx ICRU73_NEW liquid water (target 276), log-log interpolated",
            "stepping_model": "q^2 frozen at the pre-step energy over the one-step slab, "
            "S_0 = S / q^2 of the default model; g(f) = Delta E_alg / Delta E_exact - 1",
            "method_systematic": METHOD_SYSTEMATIC,
            "method_systematic_robustness": 2.0 * METHOD_SYSTEMATIC,
            "x_cap": X_CAP,
            "k_one_sided": K_ONE_SIDED,
            "consistency_tolerance": CONSISTENCY_TOL,
            "model_residual_limit": MODEL_RESIDUAL_LIMIT,
            "he_model_max_abs_r_model": he_model,
            "he_tolerance": HE_TOL,
            "fluctuation_tolerance": FLUCT_TOL,
        },
        "bands": bands,
        "pattern_before_controls": before,
        "classes_sigma_doubled": classes2,
        "pattern_sigma_doubled": pattern2,
        "robust": robust,
        "pattern": before if conclusive else "INCONCLUSIVE",
        "conclusive": conclusive,
        "inconclusive_reasons": reasons,
        "one_step_premise": premise,
        "controls": {
            "helium": helium,
            "fluctuations_off": fluct,
            "controls_pass": controls_pass,
        },
        "invalid_runs": {b: str(r["reason"]) for b, r in sorted(runs.items()) if not r["valid"]},
        "invalid_points": dict(sorted({**invalid_points, **control_invalid}.items())),
        "repeats_used": dict(sorted(repeats_used.items())),
        "provenance": {
            "archive_sha256": archive_sha,
            "member_sha256": members,
            "manifest": {
                "path": ctx.manifest_rel,
                "sha256": ctx.manifest_sha256,
                "bundle_set_digest": ctx.manifest["bundle_set_digest"],
            },
            "review_pins": {
                name: {
                    "path": p["path"],
                    "sha256": p["sha256"],
                    **{k: p["doc"][k] for k in ("review_id", "head_sha", "finished_at")},
                }
                for name, p in ctx.pins.items()
            },
            "plan_note": {"path": PLAN_PATH, "line": ctx.plan_note["line"]},
            "log_markers": {"path": ctx.markers_rel, "sha256": ctx.markers_sha256},
            "sealed_sources": ctx.sealed,
            "environment": ctx.environment,
            "records": ctx.records,
            "smoke_records": ctx.smoke_records,
            "code_sha": code_sha,
            "evaluator_head_sha": ctx.reader.head,
            "trust_file": ctx.trust_info,
            "committed_inputs": dict(sorted(ctx.reader.reads.items())),
            "runs": {
                b: {k: r[k] for k in ("run_id", "seed", "histories", "scored")}
                for b, r in sorted(runs.items())
            },
            "smokes": {
                b: {k: v[k] for k in ("run_id", "seed", "code_sha")}
                for b, v in sorted(smoke_runs.items())
            },
            "versions": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "ionmc": ionmc_version,
            },
        },
        "redistribution_note": REDISTRIBUTION_NOTE,
        "rounding_note": ROUNDING_NOTE,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--records", default=RECORDS)
    ap.add_argument("--smokes", required=True)
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--review-pin", default=REVIEW_PIN)
    ap.add_argument("--markers-pin", default=MARKERS_PIN)
    ap.add_argument("--markers", default=MARKERS)
    ap.add_argument("--repo", type=Path, default=gen.REPO_ROOT)
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--manifest", default=MANIFEST)
    ap.add_argument("--output-dir", type=Path, default=gen.REPO_ROOT / RESULTS_DIR)
    ap.add_argument("--code-sha", required=True)
    a = ap.parse_args(argv)
    doc = evaluate(
        a.repo,
        a.records,
        a.runs,
        a.review_pin,
        a.manifest,
        a.archive,
        a.code_sha,
        smokes=a.smokes,
        markers_pin=a.markers_pin,
        markers=a.markers,
    )
    a.output_dir.mkdir(parents=True, exist_ok=True)
    out = a.output_dir / RESULT_NAME
    out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    print("wrote", out, "pattern", doc["pattern"], "conclusive", doc["conclusive"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
