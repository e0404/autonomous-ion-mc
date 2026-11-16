"""Amendment I3-FU-1 step 2: thin-slab mean energy loss in TOPAS vs ionmc Bethe and ICRU73_NEW.

Usage::

    python validation/scripts/stopping/i3fu2_evaluate.py --runs DIR --archive PATH \
        --code-sha SHA [--records validation/results/stopping/v3-008-i3fu2/topas-runs] \
        [--review-pin validation/results/stopping/v3-008-i3fu2/i3fu2-review.json] \
        [--manifest validation/reference_cases/topas/i3fu2-frozen-manifest.json] \
        [--repo .] [--output-dir validation/results/stopping/v3-008-i3fu2]

Diagnostic (plan ``validation/plans/v3-008-acceptance.md``, amendment I3-FU-1 step 2 and decision
rule 2a, revision 1); no gate, no validation claim.

**Inputs.** ``--repo`` (default: this checkout) is read only as committed at its HEAD
(``r2_io.CommittedReader``: a working copy that differs fails closed); ``--records``,
``--review-pin`` and ``--manifest`` are paths in it. Records (``--records``, committed):
``index.json`` = ``{"runs": [{"case", "part", "run_id", "bundle"}, ...]}`` (one entry per
evidence or used repeat bundle, as ``validation/results/transport/v8-ss/topas-runs/index.json``)
and per run ``<run_id>/record.json`` (the signed reference-runner record, verbatim) and
``<run_id>/artifacts.json`` (the ``list_reference_artifacts`` result). Materialized outputs
(``--runs``, not committed): ``<run_id>/`` holding the listed files at their listed paths,
i.e. ``request.json``, ``stdout.txt``, ``stderr.txt``, ``work/ss_exit.header``,
``work/ss_exit.phsp`` and ``work/edep.csv`` (the layout ``ss_evaluate.load_run`` reads). The
review pin is the committed JSON ``{"review_id", "head_sha", "finished_at", ...}`` of the Codex
review that approved the bundles.

**Authentication and binding (any failure: ``RecordError`` naming the run; Codex
REVIEW-cfdff398 finding 1, mirroring ``ss_evaluate.load_run`` / ``r2_io.load_signed_run``).** Once:
``--code-sha`` is a 40-hex ancestor of HEAD; the review pin is byte-identical at ``code_sha`` and
at HEAD, its ``head_sha`` is an ancestor of (or equal to) ``code_sha``; the frozen manifest is
byte-identical at HEAD, at ``code_sha`` and at the pin's ``head_sha`` (the reviewed bundles),
has the generator's bundle set and reproduces its own ``bundle_set_digest``; the trust map is
``infrastructure/keys/trusted_attestation_keys.json`` at HEAD. Per run: Ed25519 attestation of
``record.json`` by a trusted ``reference-runner`` key (``attest.classify_record``); schema 3, not
timed out, engine topas, ``succeeded`` true, ``exit_code`` 0, no missing outputs; run id derived
from the configuration id; the record's ``code_sha`` equals ``--code-sha``; its ``started_at``
(the controlled ``reference_started``) is not before ``finished_at`` and later than the pin's
``finished_at``; the listing reproduces the signed ``output_files_sha256`` / ``output_count``
and the signed stdout/stderr/request digests; every file read hashes to its listed sha256
(``r2_io.read_verified``); the request's ``input_files`` equal the bundle committed at
``code_sha``, whose ``case.json`` / ``input.txt`` equal the manifest entry; the request case equals
that ``case.json``; seed = 20619000 + 100 part + case; the case's ``i3fu2`` block matches the
index entry and the bundle kind (evidence / repeat); configuration markers in stdout:
the ``Ph/ListProcesses`` block (``Registered Physics Processes:``), the opt4 module marker and
the EM-only check of ``r2_io.module_check`` (configuration "idd-em": ``Use ICRU90 data 1`` and the
GoudsmitSaunderson e- msc model; no hadronic marker; the physics-list name itself is not printed
by OpenTOPAS 4.3), the frozen physics lines of the bundle (``r2_io.frozen_lines_check``), an
``ionIoni:  for <GenericIon | alpha>`` block, and the fluctuation state: ``Fluctuations of
dE/dx are enabled 1`` (format verified on a materialized opt4 log) for the standard runs and
:data:`FLUCT_OFF_MARKER` for the fluctuation-off runs (to be pinned by the orchestrator from the
fluctuation-off smoke log; empty = every fluctuation-off run fails closed).

**Run validity (an INVALID run is recorded with its reason and triggers the repeat rule; plan
"Validity per run").** stderr empty, no stdout line outside "Begin processing" mentions killed or
stuck tracks, every ``Particle source Beam: Total number of histories: N`` equals the bundle's
histories; phase space (header by ``ss_evaluate.parse_header``, :data:`EXPECTED_COLUMNS`):
every primary scored exactly once (deficit 0), only the beam ion's PDG code, Z at the far face
(world frame, z = slab), finite values, first-scored flag 1, forward direction, 0 < <Delta E> <
20 % of the entrance energy; the energy-deposit cross check |EDep Sum / N / <Delta E>_T - 1| <=
0.1 %. <Delta E>_T = A E0 - mean(E_exit), SE = std(E_exit, ddof=1) / sqrt(N).

**Predictions** (fixed by this committed code before any run): dE/dz = -rho S(E/A), rho = 1 g/cm3,
RK4 over the slab (:data:`N_RK4` steps) from E_in = A E0; S from (M) ``bethe_mass_stopping``
default options at ``water(78.0)`` and (I) the libdEdx ICRU73_NEW liquid-water table
(``compare_icru73.load_icru73_water``, pinned hashes), log-log interpolated on its grid. Per run
r_M = Delta E_T / Delta E_M - 1, r_I = Delta E_T / Delta E_I - 1, s = Delta E_I / Delta E_M - 1,
SE_rel = SE / Delta E_T.

**Points (revision 1).** Per ion (C parts 0/1, O parts 2/3) and case: the point value is the
zero-thickness (Richardson) extrapolation r_0 = (f2 r(f1) - f1 r(f2)) / (f2 - f1) of r_M, r_I
and s, with f1 = 0.03 (slab t) and f2 = 0.048 (slab 1.6 t) read from the bundles'
``thickness_factor`` and asserted (:data:`FACTORS`): a bias linear in the thickness factor
cancels; s, a model-only number, is extrapolated the same way so that r_M, r_I and s refer to the
same limit. SE_rel,0 = sqrt((f2 SE_rel(f1))^2 + (f1 SE_rel(f2))^2) / (f2 - f1) (independent
runs). Consistency: |r_M(f2) - r_M(f1)| <=
:data:`CONSISTENCY_TOL`; a larger difference makes the point INVALID. **Repeat rule:** a repeat
bundle (part + 6) is used for a point only if that point's original run is present and INVALID
(validity failure; or the point's thickness-consistency failure, after which a repeat of either
slab may be used); it is recorded in ``repeats_used`` with the original's reason. A repeat offered
without an invalid original (or for a missing original) is a ``RecordError``. A point that stays
invalid (no repeat, an invalid repeat, or a second consistency failure) is excluded from the bands
and makes the result not conclusive.

**Aggregates** (``i3fu2_result.json``; per ion C/O and band [2,5), [5,10), [10,20), [20,50] over
the valid points): n, mean r_M, mean r_I, mean s, SE_b = sqrt(sum SE_rel,0^2) / n, sigma_b =
sqrt(SE_b^2 + 0.0025^2), X_b = min(0.01, |s_b| / 3), the class (U first: X_b <= 1.645 sigma_b; I:
|r_I| + 1.645 sigma_b <= X_b; H: r_I - 1.645 sigma_b > X_b; M: |r_M| + 1.645 sigma_b <= X_b; L:
r_M + 1.645 sigma_b < -X_b; else N) and the band maximum of |r_M(f2) - r_M(f1)|. Pattern:
ICRU iff both ions are I/H in [2,5) and [5,10) and I/H/U in [10,20); MODEL iff both are M/L in
[2,5) and [5,10) and M/L/U in [10,20); else MIXED (``pattern_before_controls``). ``pattern`` is
INCONCLUSIVE whenever ``conclusive`` is false: a control failed or is incomplete, an evidence
point is invalid, or a band of [2,5) / [5,10) is U. Controls: He (parts 4, cases 17/21): pass iff
both points valid and |r_M| <= 0.005 (written: ``pass``, ``max_abs_r_m``); fluctuations off (part
5 vs the C slab-t run of the same case, cases 13/21): pass iff complete and |r_M(off) - r_M(on)|
<= 0.002 (written: ``pass``, ``max_abs_difference``). Nothing per point is written; reasons carry
no values (decision 0038).
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
import re
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
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

RESULT_SCHEMA = "ionmc-i3fu2-result-2"
RESULT_NAME = "i3fu2_result.json"
RESULTS_DIR = "validation/results/stopping/v3-008-i3fu2"
REVIEW_PIN = f"{RESULTS_DIR}/i3fu2-review.json"
RECORDS = f"{RESULTS_DIR}/topas-runs"
MANIFEST = f"{gen.REL_PREFIX}/i3fu2-frozen-manifest.json"
PIN_KEYS = ("review_id", "head_sha", "finished_at")
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
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
CONSISTENCY_TOL = 0.008  # |r_M(f2) - r_M(f1)| per point (revision 1)
FACTORS = (0.03, 0.048)  # thickness factors f1 (slab t), f2 (slab 1.6 t) of every pair
HE_TOL = 0.005
FLUCT_TOL = 0.002
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
ICRU_LIKE, MODEL_LIKE = ("I", "H"), ("M", "L")
#: stdout line of the Geant4 EM-parameter dump with fluctuations enabled (format verified on a
#: materialized opt4 TOPAS log of this cache, ``Fluctuations of dE/dx are enabled  1``).
FLUCT_ON_RE = re.compile(r"^Fluctuations of dE/dx are enabled\s+1\s*$", re.M)
#: Regular expression (re.M) that a fluctuation-off run's stdout must match. TO BE PINNED BY THE
#: ORCHESTRATOR from the fluctuation-off smoke log (seed 20600441) before the evidence runs: the
#: region parameter ``EnableEnergyLossFluctuations`` may show in the global dump line above
#: (value 0) or only in the per-process lines (``fluct: 0``) -- unknown until the smoke is read.
#: Empty: every fluctuation-off run fails closed.
FLUCT_OFF_MARKER = ""
_IONIONI = "^ionIoni:\\s+for\\s+{}\\s"
CONSISTENCY_REASON = f"thickness consistency |r_M(f2) - r_M(f1)| > {CONSISTENCY_TOL}"
REDISTRIBUTION_NOTE = (
    "Aggregates only: no libdEdx/ICRU 73 or Geant4 table values, no per-energy deviations, "
    "ratios or TOPAS energy-loss values (decision 0038); the controls are written as pass/fail "
    "and maxima, the thickness check as one maximum per ion and band, reasons without values."
)

Stopping = Callable[[float], float]


class RunInvalid(Exception):
    """A run failed a validity criterion (recorded with its reason; the repeat rule applies)."""


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


_PREDICTIONS: dict[tuple[str, int, float, float], float] = {}


def predicted_loss(key: str, stopping: Stopping, a: int, e0: float, slab_mm: float) -> float:
    """:func:`slab_energy_loss` memoised by (stopping identity ``key``, a, E0, slab): the
    predictions depend on nothing else."""
    k = (key, a, e0, slab_mm)
    if k not in _PREDICTIONS:
        _PREDICTIONS[k] = slab_energy_loss(a * e0, a, slab_mm, stopping)
    return _PREDICTIONS[k]


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


def config_problems(native: str, stdout: str, ion: str, fluctuations: bool) -> list[str]:
    """Configuration markers (a problem means the run did not execute the frozen configuration:
    ``RecordError``, not a validity failure)."""
    problems = []
    if rio.registered_processes(stdout) is None:
        problems.append("the Ph/ListProcesses block (Registered Physics Processes) is absent")
    problems += rio.module_check(native, stdout, EM_CONFIG)[0]
    problems += rio.frozen_lines_check(native, EM_CONFIG)
    particle = BEAM_PARTICLE[ion]
    if not re.search(_IONIONI.format(re.escape(particle)), stdout, re.M):
        problems.append(f"no ionIoni block for {particle} in the EM printout")
    if fluctuations:
        if not FLUCT_ON_RE.search(stdout):
            problems.append("EM parameter dump does not show fluctuations enabled (1)")
    elif not FLUCT_OFF_MARKER:
        problems.append("FLUCT_OFF_MARKER is not pinned (fluctuation-off runs fail closed)")
    elif not re.search(FLUCT_OFF_MARKER, stdout, re.M):
        problems.append("the pinned fluctuation-off marker is missing from stdout")
    return problems


# =============================================================================================
# authentication and binding
# =============================================================================================


@dataclass
class Context:
    """Everything bound once per evaluation (see the module docstring)."""

    repo: Path
    reader: Any
    records: str
    pin_rel: str
    manifest_rel: str
    outputs: Path
    code_sha: str
    pin: dict[str, Any]
    pin_finished: Any
    pin_sha256: str
    manifest: dict[str, Any]
    manifest_sha256: str
    trust: dict[str, Any]
    trust_info: dict[str, Any]
    index: dict[str, str]
    bundles: dict[str, tuple[int, int]] = field(default_factory=dict)
    _committed: dict[str, dict[str, Any]] = field(default_factory=dict)


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


def load_pin(reader: Any, rel: str, code_sha: str) -> tuple[dict[str, Any], Any, str]:
    """The review pin, byte-identical at HEAD and at ``code_sha``; its ``head_sha`` must be an
    ancestor of (or equal to) ``code_sha``."""
    at_head = reader.read(rel)
    if rio.git_show(reader.repo, code_sha, rel) != at_head:
        raise RecordError(f"review pin {rel} differs between code_sha {code_sha} and HEAD")
    try:
        pin = json.loads(at_head.decode("utf-8"))
    except ValueError as exc:
        raise RecordError(f"review pin {rel} is not JSON: {exc}") from exc
    for key in PIN_KEYS:
        rio.require(pin, key, "review pin")
    if not isinstance(pin["head_sha"], str) or not SHA_RE.match(pin["head_sha"]):
        raise RecordError("review pin: head_sha is not a 40-hex commit SHA")
    finished = rio._parse_time(pin["finished_at"], "review pin finished_at")
    if not rio.is_ancestor(reader.repo, pin["head_sha"], code_sha):
        raise RecordError(
            f"review pin head_sha {pin['head_sha']} is not an ancestor of code_sha {code_sha}"
        )
    return pin, finished, hashlib.sha256(at_head).hexdigest()


def load_manifest(reader: Any, rel: str, code_sha: str, reviewed: str) -> tuple[dict, str]:
    """The frozen manifest, byte-identical at HEAD, ``code_sha`` and the reviewed commit."""
    at_head = reader.read(rel)
    for sha, what in ((code_sha, "the runs' code_sha"), (reviewed, "the reviewed head_sha")):
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


def context(
    repo: Path | str,
    records: str,
    outputs: Path | str,
    review_pin: str,
    manifest: str,
    code_sha: str,
) -> Context:
    reader = rio.CommittedReader(repo)
    if not SHA_RE.match(code_sha):
        raise RecordError(f"--code-sha {code_sha!r} is not a 40-hex commit SHA")
    if not rio.is_ancestor(reader.repo, code_sha, reader.head):
        raise RecordError(f"code_sha {code_sha} is not an ancestor of HEAD {reader.head}")
    records_rel, pin_rel, man_rel = (reader.rel(p) for p in (records, review_pin, manifest))
    pin, finished, pin_sha = load_pin(reader, pin_rel, code_sha)
    man, man_sha = load_manifest(reader, man_rel, code_sha, pin["head_sha"])
    trust, trust_info = rio.committed_trust(reader)
    return Context(
        repo=reader.repo,
        reader=reader,
        records=records_rel,
        pin_rel=pin_rel,
        manifest_rel=man_rel,
        outputs=Path(outputs),
        code_sha=code_sha,
        pin=pin,
        pin_finished=finished,
        pin_sha256=pin_sha,
        manifest=man,
        manifest_sha256=man_sha,
        trust=trust,
        trust_info=trust_info,
        index=load_index(reader, records_rel),
        bundles=run_bundles(),
    )


def committed_bundle(ctx: Context, bundle: str) -> dict[str, Any]:
    """{file: {"sha256", "bytes"}} of the bundle committed at ``code_sha`` (cached)."""
    if bundle not in ctx._committed:
        bdir = f"{gen.REL_PREFIX}/{bundle}"
        files: dict[str, Any] = {}
        blobs: dict[str, bytes] = {}
        for rel in rio.git_ls(ctx.repo, ctx.code_sha, bdir):
            blob = rio.git_show(ctx.repo, ctx.code_sha, rel)
            name = rel[len(bdir) + 1 :]
            files[name] = {"sha256": hashlib.sha256(blob).hexdigest(), "bytes": len(blob)}
            blobs[name] = blob
        if not files:
            raise RecordError(f"bundle {bundle} is not committed at {ctx.code_sha}")
        ctx._committed[bundle] = {"files": files, "blobs": blobs}
    return ctx._committed[bundle]


def load_run(ctx: Context, bundle: str, run_id: str) -> Any:
    """Authenticate and bind one signed reference-runner run of ``bundle`` (``r2_io.TopasRun``);
    every failure raises RecordError."""
    case, part = ctx.bundles[bundle]
    rdir = f"{ctx.records}/{run_id}"
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
    if code_sha != ctx.code_sha:
        raise RecordError(f"record code_sha {code_sha} != --code-sha {ctx.code_sha}")
    started = rio._parse_time(result["started_at"], f"{run_id} started_at")
    if started > rio._parse_time(result["finished_at"], f"{run_id} finished_at"):
        raise RecordError("started_at is later than finished_at")
    if not started > ctx.pin_finished:
        raise RecordError(
            f"reference_started {result['started_at']} is not after the review pin's "
            f"finished_at {ctx.pin['finished_at']} (review {ctx.pin['review_id']})"
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
        case=case,
        config="i3fu2",
        part=part,
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
    committed = committed_bundle(ctx, bundle)
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
    kind = "repeat" if part in gen.REPEAT_PARTS else "evidence"
    if (info.get("case"), info.get("part"), info.get("kind")) != (case, part, kind):
        raise RecordError(f"case.json i3fu2 block does not describe {kind} ({case}, {part})")
    if case_doc.get("seeds") != [gen.topas_seed(case, part)]:
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
    gives ``valid`` False with its ``reason``."""
    try:
        run = load_run(ctx, bundle, run_id)
        info = run.case_doc["i3fu2"]
        stdout = rio.read_verified(run, "stdout.txt").decode("utf-8", "replace")
        stderr = rio.read_verified(run, "stderr.txt")
        problems = config_problems(run.native, stdout, info["ion"], bool(info["fluctuations"]))
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


# =============================================================================================
# points, repeats and aggregation
# =============================================================================================


def classify(r_m: float, r_i: float, s: float, sigma: float) -> str:
    """Decision rule 2a (revision 1) for one band; U checked first."""
    x = min(X_CAP, abs(s) / 3.0)
    k = K_ONE_SIDED * sigma
    if x <= k:
        return "U"
    if abs(r_i) + k <= x:
        return "I"
    if r_i - k > x:
        return "H"
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
    valid; None if an original is missing), ``valid``, ``reason`` and the repeats used."""

    runs: list[dict[str, Any]] | None
    valid: bool
    reason: str | None
    used: dict[str, str]

    def thickness_difference(self) -> float | None:
        """r_M(f2) - r_M(f1) of the final attempt if both of its runs are valid."""
        if self.runs is None or len(self.runs) != 2 or not all(r["valid"] for r in self.runs):
            return None
        return float(self.runs[1]["r_M"] - self.runs[0]["r_M"])


def resolve_point(
    runs: dict[str, dict[str, Any]], case: int, parts: tuple[int, ...], ratios: Callable
) -> Point:
    """Repeat rule for one point (one slab part, or the t / 1.6 t pair)."""
    names = [gen.bundle_name(case, p) for p in parts]
    reps = [gen.bundle_name(case, p + gen.REPEAT_OFFSET) for p in parts]

    def consistent(rs: list[dict[str, Any]]) -> bool:
        return len(rs) == 1 or abs(rs[1]["r_M"] - rs[0]["r_M"]) <= CONSISTENCY_TOL

    def with_ratios(r: dict[str, Any]) -> dict[str, Any]:
        return {**r, **ratios(r)} if r["valid"] else r

    if any(n not in runs for n in names):
        for n, rep in zip(names, reps, strict=True):
            if n not in runs and rep in runs:
                raise RecordError(f"repeat {rep} offered but the original {n} is missing")
        return Point(None, False, "missing run", {})
    first = [with_ratios(runs[n]) for n in names]
    failed = {i: str(r["reason"]) for i, r in enumerate(first) if not r["valid"]}
    by_consistency = not failed and not consistent(first)
    if by_consistency:
        failed = dict.fromkeys(range(len(parts)), CONSISTENCY_REASON)
    if not failed:
        for rep in reps:
            if rep in runs:
                raise RecordError(f"repeat {rep} offered without an invalid original")
        return Point(first, True, None, {})
    chosen, used = list(first), {}
    for i, rep in enumerate(reps):
        if rep in runs:
            if i not in failed:
                raise RecordError(f"repeat {rep} offered but the original {names[i]} is valid")
            chosen[i] = with_ratios(runs[rep])
            used[rep] = f"original {names[i]} invalid: {failed[i]}"
    if not used or (not by_consistency and any(reps[i] not in runs for i in failed)):
        why = CONSISTENCY_REASON if by_consistency else "invalid run"
        return Point(chosen, False, f"{why}; no valid repeat", used)
    if any(not r["valid"] for r in chosen):
        return Point(chosen, False, "repeat run invalid (second invalid of the point)", used)
    if not consistent(chosen):
        return Point(chosen, False, f"{CONSISTENCY_REASON} (second invalid of the point)", used)
    return Point(chosen, True, None, used)


def band_aggregates(points: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Per band: n, means of r_M, r_I, s (extrapolated points), SE_b, sigma_b, X_b, class and the
    maximum |r_M(f2) - r_M(f1)| over the band's points with both runs valid."""
    out: dict[str, dict[str, Any]] = {}
    for name, *_ in BANDS:
        pts = [p for p in points if band_of(p["e0"]) == name]
        valid = [p for p in pts if p["valid"]]
        diffs = [abs(p["d_thick"]) for p in pts if p["d_thick"] is not None]
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
            "max_abs_thickness_difference": max(diffs) if diffs else None,
        }
        if n:
            r_m = float(np.mean([p["r_M"] for p in valid]))
            r_i = float(np.mean([p["r_I"] for p in valid]))
            s = float(np.mean([p["s"] for p in valid]))
            se_b = math.sqrt(sum(p["se_rel"] ** 2 for p in valid)) / n
            sigma = math.sqrt(se_b**2 + METHOD_SYSTEMATIC**2)
            band.update(
                mean_r_M=r_m,
                mean_r_I=r_i,
                mean_s=s,
                se_b=se_b,
                sigma_b=sigma,
                x_b=min(X_CAP, abs(s) / 3.0),
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
    archive_sha256: str | None = None,
    member_sha256: dict[str, str] | None = None,
) -> dict[str, Any]:
    """The result document (aggregates only). ``archive_sha256``/``member_sha256`` override the
    pinned libdEdx hashes of ``compare_icru73`` (tests with a synthetic archive only)."""
    ctx = context(repo, records, runs_dir, review_pin, manifest, code_sha)
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
            pt = resolve_point(runs, case, parts, ratios)
            repeats_used.update(pt.used)
            e0 = float(gen.CASES[case])
            d = pt.thickness_difference()
            if not pt.valid or pt.runs is None:
                invalid_points[f"{ion} case {case}"] = str(pt.reason)
                points[ion].append({"e0": e0, "valid": False, "d_thick": d})
                continue
            t, t2 = pt.runs
            f1, f2 = (float(r["info"]["thickness_factor"]) for r in (t, t2))
            if (f1, f2) != FACTORS:
                raise RecordError(f"{ion} case {case}: thickness factors {(f1, f2)} != {FACTORS}")
            chosen_t[(ion, case)] = t
            points[ion].append(
                {
                    "e0": e0,
                    "valid": True,
                    "d_thick": d,
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
    helium = {
        "pass": len(he_r) == len(HE_CASES) and all(abs(v) <= HE_TOL for v in he_r),
        "max_abs_r_m": max((abs(v) for v in he_r), default=None),
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
    fluct = {
        "pass": len(fl_d) == len(FLUCT_CASES) and all(abs(v) <= FLUCT_TOL for v in fl_d),
        "max_abs_difference": max((abs(v) for v in fl_d), default=None),
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
    conclusive = not reasons
    before = pattern(bands)
    try:
        ionmc_version: str | None = importlib.metadata.version("ionmc")
    except importlib.metadata.PackageNotFoundError:
        ionmc_version = None
    return {
        "schema": RESULT_SCHEMA,
        "task_id": "V3-008",
        "amendment": "I3-FU-1 step 2 (revision 1)",
        "claim": "diagnostic; no gate; no validation claim",
        "code_sha": code_sha,
        "method": {
            "quantity": "<Delta E>_T = A E0 - mean(E_exit); predictions RK4 over the slab",
            "point": "zero-thickness extrapolation (f2 r1 - f1 r2) / (f2 - f1) of r_M, r_I, s",
            "thickness_factors": list(FACTORS),
            "rk4_steps": N_RK4,
            "model": "ionmc bethe_mass_stopping default options, water(78.0)",
            "table": "libdEdx ICRU73_NEW liquid water (target 276), log-log interpolated",
            "method_systematic": METHOD_SYSTEMATIC,
            "x_cap": X_CAP,
            "k_one_sided": K_ONE_SIDED,
            "consistency_tolerance": CONSISTENCY_TOL,
            "he_tolerance": HE_TOL,
            "fluctuation_tolerance": FLUCT_TOL,
        },
        "bands": bands,
        "pattern_before_controls": before,
        "pattern": before if conclusive else "INCONCLUSIVE",
        "conclusive": conclusive,
        "inconclusive_reasons": reasons,
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
            "review_pin": {
                "path": ctx.pin_rel,
                "sha256": ctx.pin_sha256,
                **{k: ctx.pin[k] for k in PIN_KEYS},
            },
            "records": ctx.records,
            "code_sha": code_sha,
            "evaluator_head_sha": ctx.reader.head,
            "trust_file": ctx.trust_info,
            "committed_inputs": dict(sorted(ctx.reader.reads.items())),
            "runs": {
                b: {k: r[k] for k in ("run_id", "seed", "histories", "scored")}
                for b, r in sorted(runs.items())
            },
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
    ap.add_argument("--records", default=RECORDS)
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--review-pin", default=REVIEW_PIN)
    ap.add_argument("--repo", type=Path, default=gen.REPO_ROOT)
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--manifest", default=MANIFEST)
    ap.add_argument("--output-dir", type=Path, default=gen.REPO_ROOT / RESULTS_DIR)
    ap.add_argument("--code-sha", required=True)
    a = ap.parse_args(argv)
    doc = evaluate(a.repo, a.records, a.runs, a.review_pin, a.manifest, a.archive, a.code_sha)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    out = a.output_dir / RESULT_NAME
    out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    print("wrote", out, "pattern", doc["pattern"], "conclusive", doc["conclusive"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
