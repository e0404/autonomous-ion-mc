"""Validation steps of the V3-004 acceptance plan (rows A7, A8, A9, A11, A13, A15, A16).

Usage::

    python validation/scripts/transport/steps_v4.py <step> [--workers N] [--scale F] ...

The steps run in the suites ``lv4`` and ``hr4`` of ``run_suite.py`` and follow the conventions of
``steps.py`` (the V3-003 steps, whose helpers are reused): each prints one JSON document between
``#JSON-BEGIN`` and ``#JSON-END`` and exits 0 iff every frozen criterion of
``validation/plans/v3-004-acceptance.md`` passes; the history counts are the frozen ones unless
``--scale`` (or ``--seeds`` for A9) reduces them, in which case the document says
``"reduced": true`` and a pass is not a conformant result. Everything uses the offline analytic
Bethe stopping source (I = 78 eV) and the synthetic lookup fixture ``tests/data/synthetic_lookup.json``.

Seeds derive deterministically from one ``--seed-base`` (default 20361004, the V3-004 qualification
base; 20351004 is the rehearsal base, see the plan): A7 ``base + i`` (i-th step length), A8
``base + 8``, A9 ``base + 10000 + i`` (i-th seed), A11 LV ``base + 11``, A11 HR ``base + 1000 k``
(k = 1 python, 2 cpu32, 3 cpu64, 4 cuda32, as T12), A13 ``base + 13``, A15 and A16 fixed small
configurations at ``base`` (A16 digests use the fixed seed documented in ``a16_digest.py``).

Definitions made here where the frozen rows leave them open (all recorded in the documents):

* plateau: the depth bins whose centre lies in ``[0.1 R, 0.7 R]``, ``R`` the CSDA range;
* A7 "LET_d vs LET_d^eps" (plan footnote 2): ``D = LS2/LS - ES/E_step`` per bin must satisfy
  ``|D| <= 4 sigma_D + delta_v``. ``sigma_D`` is the delta-method standard error from the batch
  covariance of the four channels (the two ratios are strongly correlated). ``delta_v`` is the
  analytic bound of the expected-share shift of the straggled end-energy ramp (plan delta item 6):
  the expected deposit of a piece differs from the mean-loss apportionment by
  ``l_p g kappa (t_mid/s - 1/2) / S`` with ``g = -dS/dE``, ``kappa = Var(eps)/s`` the Bohr variance
  per path length (``EM.bohr_variance`` at 1 mm) and ``|t_mid/s - 1/2| <= 1/2``, i.e. a relative
  shift of at most ``r = g kappa / (2 S^2)``. It sums to zero over a whole step and survives only
  for step starts or ends without their partners in the voxel. Both estimators are weighted sums
  of the same ``S`` ramp values, so the shift moves ``LET_d^eps`` by at most ``r`` of its value
  (the sensitivity ``(S_p - R)/E_step`` of the ratio to a piece shift is bounded by the value
  itself): ``delta_v = r_v LET_d,v`` with ``S`` and ``E`` of the voxel taken at ``S_v = LET_d,v``
  (``E_v`` from the inverse of the table ``S_w``, ``g`` from its numerical derivative). No free
  constant enters.
* A9 reference: the leave-one-out pooled ratio of the other seeds, ``z_i = (R_i - R_ref,i) /
  sqrt(V_i + V_ref,i)`` with ``V_ref,i`` the mean variance over the other seeds divided by their
  number; the three depths are the plateau mid point (``0.5 R``), the maximum and the distal 80 %
  point of the pooled dose profile.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

import steps as base
from ionmc.config import DiagnosticsOptions, PhysicsOptions, RunOptions, SimulationConfig
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.let_offline import let_from_spectrum, log_edges
from ionmc.lookup import LookupTable
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.scoring import ScoringGrid, TallyRequest, reduce_ratio
from ionmc.simulation import Result, Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport import parity
from ionmc.transport.parity_channels import (
    channel_labels,
    compare_channel_partition,
    compare_channel_runs,
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SYNTHETIC = REPO / "tests" / "data" / "synthetic_lookup.json"
LK_LET = LookupTable.from_file(SYNTHETIC)
LK_E = LookupTable(
    "e_lin", "synthetic energy lookup (test fixture, no biological meaning)", "1",
    "energy_per_nucleon_mev", "log", np.geomspace(0.5, 200.0, 50),
    {"proton": np.linspace(1.0, 3.0, 50)}, "none: mathematical test function", "CC0-1.0",
    "validation/scripts/transport/steps_v4.py (f = 1 .. 3 linear in the sample index)", True,
)  # fmt: skip
LOOKUPS = (LK_LET, LK_E)
QUALIFICATION_SEED_BASE = 20361004
REHEARSAL_SEED_BASE = 20351004

ENERGY_MEV = 150.0
PLATEAU = (0.1, 0.7)  # fractions of the CSDA range
A7_STEPS_MM = (0.1, 0.5, 1.0)
A7_HISTORIES = 100_000
A7_BATCHES = 20
A7_TRACE_HISTORIES = 400
A8_HISTORIES = 100_000
A8_BATCHES = 20
A8_BINS_PER_DECADE = 200
A8_E_MIN, A8_DECADES = 1.0, 2.4  # 1 .. 251 MeV per nucleon
A9_SEEDS = 200
A9_HISTORIES = 100_000
A9_BATCHES = 100
A13_HISTORIES = 400_000
A13_BATCHES = 20
A15_HISTORIES = {"warp-cpu": 200_000, "warp-cuda": 1_000_000}
A15_BATCHES = 20
A11_LV_K = 256
A16_BASELINE = "a524f209"
A16_SPECS = (
    "t1:python:float64",
    "t1:warp-cpu:float64",
    "t1:warp-cpu:float32",
    "t13:python:float64",
    "t13:warp-cpu:float32",
    "t13:warp-cpu:float64",
)
GRAZ_NOTE = (
    "Grassberger and Paganetti 2011 (Phys. Med. Biol. 56:6677): distal maximum of the primary "
    "LET_d about 12 keV/um (another beam; related model only, no pass/fail)"
)


# -- helpers -------------------------------------------------------------------------------------
def scaled(n: int, scale: float, floor: int, multiple: int = 1) -> int:
    """``n`` times ``scale`` (a multiple of ``multiple``, at least ``floor``); exactly ``n`` for
    ``scale == 1``."""
    if scale == 1.0:
        return n
    return max(floor, int(n * scale) // multiple * multiple)


def finish4(doc: dict[str, Any], frozen: Any, used: Any, reduced: bool) -> int:
    doc["reduced"] = bool(reduced)
    doc["frozen_histories"] = frozen
    doc["histories"] = used
    return base.emit(doc)


def depth_setup(
    energy: float, dz: float = 1.0, half_mm: float = 40.0
) -> tuple[BoxPhantom, ScoringGrid, int, float]:
    """Water box of depth ``1.1 R`` and a ``(1, 1, nz)`` depth grid of ``dz`` bins named
    ``dose`` (the lateral extent contains the beam; deposits outside it are tallied separately)."""
    r = base.r_csda_mm(energy)
    nz = int(math.ceil(1.1 * r / dz))
    geo = BoxPhantom((-half_mm, -half_mm, 0.0), (2 * half_mm, 2 * half_mm, nz * dz), WATER)
    grid = ScoringGrid((-half_mm, -half_mm, 0.0), (2 * half_mm, 2 * half_mm, dz), (1, 1, nz),
                       name="dose")  # fmt: skip
    return geo, grid, nz, r


def all_tallies(grid: str, edges: tuple[float, ...]) -> tuple[TallyRequest, ...]:
    """One request of every quantity kind (and a species/generation-restricted one)."""
    return (
        TallyRequest("edep", grid, "edep"),
        TallyRequest("lt", grid, "let_t"),
        TallyRequest("ld", grid, "let_d"),
        TallyRequest("le", grid, "let_d_eps"),
        TallyRequest("fl", grid, "fluence"),
        TallyRequest("fe", grid, "lookup_sum", lookup=LK_LET.name),
        TallyRequest("fa", grid, "lookup_dose_avg", lookup=LK_E.name),
        TallyRequest("sp", grid, "fluence_spectrum", energy_edges_mev_per_u=edges),
        TallyRequest("lt_p", grid, "let_t", species=("proton",), generation="primary"),
    )


def batch_values(res: Result, ci: int) -> NDArray[np.float64]:
    """Per-batch, per-primary values ``[B, size]`` of channel ``ci`` in channel units."""
    plan = res.effective_config.channels
    assert plan is not None
    hpb = res.n_histories // res.n_batches
    return res.channel_batches(ci).astype(np.float64) * plan.channels[ci].quantum / hpb


def ratio_difference(
    x1: NDArray[Any], y1: NDArray[Any], x2: NDArray[Any], y2: NDArray[Any]
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """``D = X1/Y1 - X2/Y2`` of the batch means (arrays ``[B, m]``), its delta-method standard
    error from the batch covariance of the four channels, and the mask of bins where both
    denominators are positive in at least half of the batches."""
    b = x1.shape[0]
    v = np.stack([x1, y1, x2, y2])  # (4, B, m)
    m = v.mean(axis=1)  # (4, m)
    d = v - m[:, None, :]
    cov = np.einsum("ibm,jbm->ijm", d, d) / (b * (b - 1))
    ok = (m[1] > 0.0) & (m[3] > 0.0) & ((y1 > 0).sum(0) >= max(2, math.ceil(b / 2)))
    ok &= (y2 > 0).sum(0) >= max(2, math.ceil(b / 2))
    y1m, y2m = np.where(ok, m[1], 1.0), np.where(ok, m[3], 1.0)
    g = np.stack([1.0 / y1m, -m[0] / y1m**2, -1.0 / y2m, m[2] / y2m**2])
    var = np.einsum("im,ijm,jm->m", g, cov, g)
    diff = m[0] / y1m - m[2] / y2m
    return np.where(ok, diff, np.nan), np.where(ok, np.sqrt(np.maximum(var, 0.0)), np.nan), ok


def profile(res: Result, grid: str = "dose") -> dict[str, Any]:
    """LET_t, LET_d, LET_d^eps profiles (delta-method) and the dose profile of one run."""
    lab = channel_labels(res, grid)
    v = {k: batch_values(res, i) for k, i in lab.items()}
    lt = reduce_ratio(v["LS"], v["L"])
    ld = reduce_ratio(v["LS2"], v["LS"])
    le = reduce_ratio(v["ES"], v["E_step"])
    dose = res.grid(grid).energy_mev.reshape(-1)
    return {"v": v, "let_t": lt, "let_d": ld, "let_d_eps": le, "dose": dose}


def expected_share_bound(res: Result, let_d: NDArray[np.float64]) -> NDArray[np.float64]:
    """``delta_v = r_v LET_d,v`` of the A7 criterion (module docstring); NaN where undefined."""
    from ionmc._wpfunc import python_twin
    from ionmc.physics.em import make_em

    em = python_twin(make_em)
    tab = res.effective_config.tables
    mass = float(res.requested_config.source.projectile.mass_mev)
    zoa = float(WATER.z_over_a)
    rho = float(WATER.density_g_cm3)
    out = np.full(let_d.shape, np.nan)
    for i, sv in enumerate(let_d):
        if not np.isfinite(sv):
            continue
        lo, hi = 1.0, 300.0  # S_w decreases monotonically with E above the stopping maximum
        for _ in range(80):
            mid = math.sqrt(lo * hi)
            if tab.s_water(mid) > sv:
                lo = mid
            else:
                hi = mid
        e = math.sqrt(lo * hi)
        h = 1e-4
        g = -(tab.s_water(e * (1 + h)) - tab.s_water(e * (1 - h))) / (2 * h * e)
        kappa = float(em.bohr_variance(e, mass, 1.0, zoa, rho, 1.0))
        out[i] = 0.5 * g * kappa / (sv * sv) * sv
    return out


def plateau_mask(nz: int, dz: float, r: float) -> NDArray[np.bool_]:
    z = (np.arange(nz) + 0.5) * dz
    return (z >= PLATEAU[0] * r) & (z <= PLATEAU[1] * r)


def counters_clean(res: Result) -> bool:
    return bool(res.valid and not res.counters.any_nonzero)


# -- A7 ------------------------------------------------------------------------------------------
def step_a7(a: argparse.Namespace) -> int:
    n = scaled(A7_HISTORIES, a.scale, 2000, A7_BATCHES)
    geo, grid, nz, r = depth_setup(ENERGY_MEV)
    pm = plateau_mask(nz, 1.0, r)
    reqs = (
        TallyRequest("edep", "dose", "edep"),
        TallyRequest("lt", "dose", "let_t"),
        TallyRequest("ld", "dose", "let_d"),
        TallyRequest("le", "dose", "let_d_eps"),
    )
    runs: dict[float, dict[str, Any]] = {}
    res_by_s: dict[float, Result] = {}
    est: dict[float, float] = {}
    k_trace = min(a.trace_histories, n)
    for i, s in enumerate(A7_STEPS_MM):
        res = base.run_cfg(
            energy=ENERGY_MEV, geometry=geo, scoring=(grid,), seed=a.seed + i, n=n,
            n_batches=A7_BATCHES, workers=a.workers, timeout=a.timeout, max_step=s,
            tallies=reqs,
        )  # fmt: skip
        p = profile(res)
        p["valid"] = counters_clean(res)
        runs[s] = p
        res_by_s[s] = res
        # negative control: the trace-based eps/l estimator over the plateau steps of the first
        # ``k_trace`` histories (the per-step trace exists in float64 only: a separate small run)
        tres = base.run_cfg(
            energy=ENERGY_MEV, geometry=geo, scoring=(grid,), seed=a.seed + i, n=k_trace,
            n_batches=2, workers=1, timeout=a.timeout, max_step=s, precision="float64",
            diag=DiagnosticsOptions(trace_histories=k_trace),
        )  # fmt: skip
        tr = tres.diagnostics["trace"]
        eps, ln = np.asarray(tr["deposit_mev"]), np.asarray(tr["step_mm"])
        zmid = np.asarray(tr["z_mm"]) - 0.5 * ln * np.asarray(tr["uz"])
        sel = (zmid >= PLATEAU[0] * r) & (zmid <= PLATEAU[1] * r) & (ln > 0.0) & (eps > 0.0)
        est[s] = float((eps[sel] ** 2 / ln[sel]).sum() / eps[sel].sum())
        runs[s]["trace_steps_in_plateau"] = int(sel.sum())
    # plateau LET_d, pairwise
    pairs = {}
    ok_pairs = True
    for i, s1 in enumerate(A7_STEPS_MM):
        for s2 in A7_STEPS_MM[i + 1 :]:
            m1, m2 = runs[s1]["let_d"], runs[s2]["let_d"]
            sel = pm & m1.defined_mask & m2.defined_mask
            rel = np.abs(m1.mean - m2.mean)[sel] / (0.5 * (m1.mean + m2.mean)[sel])
            worst = float(rel.max()) if rel.size else math.nan
            ok = bool(rel.size == pm.sum() and worst <= 0.01)
            pairs[f"{s1}_vs_{s2}"] = {"max_rel_diff": worst, "bins": int(rel.size), "pass": ok}
            ok_pairs &= ok
    # LET_d vs LET_d^eps within 4 sigma where dose > 10 % of max
    eps_cmp = {}
    ok_eps = True
    for s in A7_STEPS_MM:
        p = runs[s]
        lab_v = p["v"]
        d, sd, defined = ratio_difference(lab_v["LS2"], lab_v["LS"], lab_v["ES"], lab_v["E_step"])
        sel = (p["dose"] > 0.1 * p["dose"].max()) & defined
        delta = expected_share_bound(res_by_s[s], p["let_d"].mean)
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.where(sd > 0.0, d / np.where(sd > 0.0, sd, 1.0), np.nan)
        zs = z[sel]
        allow = 4.0 * sd + delta
        ratio_allow = np.abs(d[sel]) / allow[sel]
        ratio_delta = np.abs(d[sel]) / delta[sel]
        n_sel = int(((p["dose"] > 0.1 * p["dose"].max())).sum())
        worst = float(np.nanmax(np.abs(zs))) if zs.size else math.nan
        ok = bool(
            zs.size == n_sel and np.all(np.isfinite(ratio_allow)) and ratio_allow.max() <= 1.0
            and p["valid"]
        )  # fmt: skip
        z_pl = pm[sel]
        z_di = (np.arange(nz) + 0.5 > r * 0.9)[sel]
        rel = np.abs(d[sel]) / p["let_d"].mean[sel]
        eps_cmp[str(s)] = {
            "max_abs_z": worst,
            "bins": int(zs.size),
            "bins_dose_above_10pct": n_sel,
            "max_rel_diff": float(rel.max()) if rel.size else math.nan,
            "rms_z": float(np.sqrt(np.nanmean(zs**2))) if zs.size else math.nan,
            "max_abs_d_over_4sigma_plus_delta": float(ratio_allow.max()),
            "max_abs_d_over_delta": float(ratio_delta.max()),
            "max_delta_rel": float(np.nanmax(delta[sel] / p["let_d"].mean[sel])),
            "plateau_max_d_over_4sigma_plus_delta": float(ratio_allow[z_pl].max())
            if z_pl.any() else None,
            "distal_max_d_over_4sigma_plus_delta": float(ratio_allow[z_di].max())
            if z_di.any() else None,
            "pass": ok,
        }
        ok_eps &= ok
    change = {
        f"{s}_vs_1.0": abs(est[s] - est[1.0]) / est[1.0] for s in A7_STEPS_MM if s != 1.0
    }
    ctrl_ok = bool(change["0.1_vs_1.0"] >= 0.10)
    zc = (np.arange(nz) + 0.5)[pm]
    doc = {
        "step": "a7",
        "energy_mev": ENERGY_MEV,
        "n_batches": A7_BATCHES,
        "plateau_mm": [float(zc[0]), float(zc[-1])],
        "range_mm": r,
        "plateau_let_d_kev_um": {
            str(s): {
                "min": float(runs[s]["let_d"].mean[pm].min()),
                "max": float(runs[s]["let_d"].mean[pm].max()),
            }
            for s in A7_STEPS_MM
        },
        "step_independence": {"pairs": pairs, "bound": 0.01, "pass": ok_pairs},
        "let_d_vs_let_d_eps": {"by_step": eps_cmp, "criterion": "|D| <= 4 sigma_D + delta_v (footnote 2)", "pass": ok_eps},
        "negative_control": {
            "estimator": "sum(eps^2/l) / sum(eps) over the traced plateau steps",
            "value_kev_um": {str(s): est[s] for s in A7_STEPS_MM},
            "relative_change": change,
            "traced_histories": k_trace,
            "traced_plateau_steps": {str(s): runs[s]["trace_steps_in_plateau"] for s in runs},
            "required_change": 0.10,
            "pass": ctrl_ok,
        },
        "pass": bool(ok_pairs and ok_eps and ctrl_ok),
    }
    return finish4(doc, A7_HISTORIES, n, n < A7_HISTORIES)


# -- A8 ------------------------------------------------------------------------------------------
def step_a8(a: argparse.Namespace) -> int:
    n = scaled(A8_HISTORIES, a.scale, 2000, A8_BATCHES)
    geo, grid, nz, r = depth_setup(ENERGY_MEV)
    edges = tuple(float(e) for e in log_edges(A8_E_MIN, A8_BINS_PER_DECADE, A8_DECADES))
    reqs = (
        TallyRequest("lt", "dose", "let_t"),
        TallyRequest("ld", "dose", "let_d"),
        TallyRequest("sp", "dose", "fluence_spectrum", energy_edges_mev_per_u=edges),
    )
    res = base.run_cfg(
        energy=ENERGY_MEV, geometry=geo, scoring=(grid,), seed=a.seed + 8, n=n,
        n_batches=A8_BATCHES, workers=a.workers, timeout=a.timeout, max_step=1.0, tallies=reqs,
    )  # fmt: skip
    q = res.grid("dose").quantities
    off = let_from_spectrum(
        q["sp"].mean.reshape(nz, -1), edges, res.effective_config.tables.s_water
    )
    dose = res.grid("dose").energy_mev.reshape(-1)
    out: dict[str, Any] = {}
    ok = counters_clean(res)
    for name, offline in (("let_t", off.let_t), ("let_d", off.let_d)):
        on = q["lt" if name == "let_t" else "ld"]
        defined = on.defined_mask.reshape(-1) & np.isfinite(offline)
        n_def = int(on.defined_mask.sum())
        rel = np.abs(offline - on.mean.reshape(-1))[defined] / on.mean.reshape(-1)[defined]
        hi = defined & (dose > 0.1 * dose.max())
        rel_hi = np.abs(offline - on.mean.reshape(-1))[hi] / on.mean.reshape(-1)[hi]
        worst = float(rel.max()) if rel.size else math.nan
        passed = bool(rel.size == n_def and n_def > 0 and worst <= 0.005)
        out[name] = {
            "defined_bins": n_def,
            "compared_bins": int(rel.size),
            "max_rel_diff": worst,
            "worst_depth_mm": float(np.nonzero(defined)[0][np.argmax(rel)] + 0.5) if rel.size
            else None,
            "max_rel_diff_dose_above_10pct": float(rel_hi.max()) if rel_hi.size else math.nan,
            "bound": 0.005,
            "pass": passed,
        }
        ok &= passed
    outside = float(np.nanmax(off.outside_fraction)) if np.isfinite(off.outside_fraction).any() \
        else 0.0  # fmt: skip
    ok &= outside == 0.0
    doc = {
        "step": "a8",
        "energy_mev": ENERGY_MEV,
        "spectrum": {
            "bins_per_decade": A8_BINS_PER_DECADE,
            "edges_mev_per_u": [edges[0], edges[-1]],
            "n_bins": len(edges) - 1,
            "max_weight_outside_edges": outside,
        },
        "a8": out,
        "n_batches": A8_BATCHES,
        "pass": bool(ok),
    }
    return finish4(doc, A8_HISTORIES, n, n < A8_HISTORIES)


# -- A9 (two parts that a9-compare merges) --------------------------------------------------------
A9_PARTS = 2
A9_BAND = (0.91, 0.99)  # 99 % binomial interval for p = 0.95 at n = 200 (plan footnote 3)


def a9_split(part: tuple[int, int], seeds: int) -> range:
    i, n = part
    return range((i - 1) * seeds // n, i * seeds // n)


def a9_file(i: int, n: int) -> str:
    return f"a9-part-{i}-of-{n}.npz"


def step_a9_part(a: argparse.Namespace) -> int:
    i, nparts = (int(x) for x in a.part.split("/"))
    if not 1 <= i <= nparts:
        raise SystemExit("--part must be i/n with 1 <= i <= n")
    n = scaled(A9_HISTORIES, a.scale, A9_BATCHES * 10, A9_BATCHES)
    geo, grid, nz, r = depth_setup(ENERGY_MEV)
    reqs = (TallyRequest("ld", "dose", "let_d"),)
    idx = list(a9_split((i, nparts), a.seeds))
    xs, ys, doses = [], [], []
    ok = True
    t_start = time.perf_counter()
    for c, k in enumerate(idx):
        res = base.run_cfg(
            energy=ENERGY_MEV, geometry=geo, scoring=(grid,), seed=a.seed + 10_000 + k, n=n,
            n_batches=A9_BATCHES, workers=a.workers, timeout=a.timeout, max_step=1.0,
            tallies=reqs,
        )  # fmt: skip
        ok &= counters_clean(res)
        lab = channel_labels(res, "dose")
        xs.append(batch_values(res, lab["LS2"]))
        ys.append(batch_values(res, lab["LS"]))
        doses.append(res.grid("dose").energy_mev.reshape(-1))
        if (c + 1) % 10 == 0:
            el = time.perf_counter() - t_start
            print(f"# a9 seed {c + 1}/{len(idx)} elapsed {el:.0f} s", file=sys.stderr, flush=True)
    elapsed = time.perf_counter() - t_start
    meta = {
        "format": 1, "part": i, "n_parts": nparts, "seed_base": a.seed, "seeds_total": a.seeds,
        "seed_indices": idx, "histories_per_seed": n, "batches": A9_BATCHES,
        "git_sha": os.environ.get("IONMC_RUN_SHA", "unknown"),
    }  # fmt: skip
    path = Path(a.out_dir) / a9_file(i, nparts)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, meta=np.array(json.dumps(meta, sort_keys=True)), x=np.stack(xs),
             y=np.stack(ys), dose=np.stack(doses))  # fmt: skip
    doc = {
        "step": "a9-part",
        "file": path.name,
        "sample_file": f"samples/{path.name}",
        "sha256": base.file_sha256(path),
        "meta": meta,
        "timing": {
            "elapsed_s": elapsed,
            "seconds_per_seed": elapsed / max(1, len(idx)),
            "histories_per_second": len(idx) * n / elapsed,
        },
        "pass": bool(ok),
    }
    return finish4(doc, A9_HISTORIES, n, n < A9_HISTORIES or a.seeds < A9_SEEDS)


def _a9_load(dirs: list[Path], expected_sha: str, a: argparse.Namespace) -> tuple[Any, ...]:
    base.verify_archives(dirs, expected_sha)
    found: dict[str, tuple[dict[str, Any], str]] = {}
    for d in dirs:
        for txt in sorted(d.glob("[0-9][0-9]-*.txt")):
            text = txt.read_text()
            if "#JSON-BEGIN" not in text:
                continue
            try:
                doc = json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])
            except json.JSONDecodeError:
                continue
            if doc.get("step") == "a9-part":
                sha = next((ln.split(": ", 1)[1] for ln in text.splitlines()
                            if ln.startswith("# git_sha: ")), "")  # fmt: skip
                if doc.get("pass") is not True:
                    raise SystemExit(f"{txt.name}: producer step did not pass")
                found[doc["file"]] = (doc, sha)
    xs, ys, ds = {}, {}, {}
    n_hist = None
    for i in range(1, A9_PARTS + 1):
        fname = a9_file(i, A9_PARTS)
        if fname not in found:
            raise SystemExit(f"missing producer of {fname}")
        doc, sha = found[fname]
        path = next((d / "samples" / fname for d in dirs if (d / "samples" / fname).exists()), None)
        if path is None:
            raise SystemExit(f"{fname} recorded but not found in {dirs}")
        if base.file_sha256(path) != doc["sha256"]:
            raise SystemExit(f"sha256 mismatch for {fname}")
        data = np.load(path)
        meta = json.loads(str(data["meta"]))
        want = {"part": i, "n_parts": A9_PARTS, "seed_base": a.seed, "seeds_total": a.seeds,
                "batches": A9_BATCHES, "git_sha": expected_sha}  # fmt: skip
        for k, v in want.items():
            if meta.get(k) != v or doc["meta"].get(k) != v:
                raise SystemExit(f"{fname}: metadata {k}={meta.get(k)!r}, expected {v!r}")
        if sha != expected_sha:
            raise SystemExit(f"{fname} produced at SHA {sha}")
        if meta["seed_indices"] != list(a9_split((i, A9_PARTS), a.seeds)):
            raise SystemExit(f"{fname}: seed indices differ from the deterministic split")
        if n_hist is not None and meta["histories_per_seed"] != n_hist:
            raise SystemExit("parts differ in histories per seed")
        n_hist = meta["histories_per_seed"]
        for j, k in enumerate(meta["seed_indices"]):
            xs[k], ys[k], ds[k] = data["x"][j], data["y"][j], data["dose"][j]
    if sorted(xs) != list(range(a.seeds)):
        raise SystemExit("the parts do not tile the seeds")
    order = range(a.seeds)
    return (np.stack([xs[k] for k in order]), np.stack([ys[k] for k in order]),
            np.stack([ds[k] for k in order]), n_hist)  # fmt: skip


def step_a9_compare(a: argparse.Namespace) -> int:
    expected = os.environ.get("IONMC_RUN_SHA", "")
    x, y, doses, n = _a9_load([Path(d) for d in a.dirs], expected, a)
    seeds = x.shape[0]
    _, _, nz, r = depth_setup(ENERGY_MEV)
    dose = doses.mean(axis=0)
    ok = True
    st = [reduce_ratio(x[i], y[i]) for i in range(seeds)]
    mean = np.stack([s.mean for s in st])  # (S, nz)
    var = np.stack([s.variance_of_mean for s in st])
    defined = np.stack([s.defined_mask for s in st])
    tot_x, tot_y = x.sum(axis=1), y.sum(axis=1)  # (S, nz) sums over batches
    loo_x, loo_y = tot_x.sum(0) - tot_x, tot_y.sum(0) - tot_y
    with np.errstate(divide="ignore", invalid="ignore"):
        ref = loo_x / loo_y  # leave-one-out pooled ratio per seed
        v_ref = np.nanmean(var, axis=0) / (seeds - 1)
        z = (mean - ref) / np.sqrt(var + v_ref)
    i_peak = int(np.argmax(dose))
    r80 = parity.r80_of(dose, 1.0)
    depths = {
        "plateau_0.5R": int(0.5 * r),
        "bragg_peak": i_peak,
        "distal_80pct": int(r80) if math.isfinite(r80) else i_peak,
    }
    cov: dict[str, Any] = {}
    for name, b in depths.items():
        zb = z[:, b]
        good = bool(np.all(defined[:, b]) and np.all(np.isfinite(zb)))
        frac = float(np.mean(np.abs(zb) < 1.96)) if good else math.nan
        passed = bool(good and A9_BAND[0] <= frac <= A9_BAND[1])
        cov[name] = {
            "bin": b, "depth_mm": b + 0.5,
            "dose_fraction_of_max": float(dose[b] / dose.max()),
            "coverage": frac,
            "z_mean": float(np.mean(zb)) if good else math.nan,
            "z_std": float(np.std(zb, ddof=1)) if good else math.nan,
            "let_d_mean_kev_um": float(np.mean(mean[:, b])),
            "all_seeds_defined": bool(np.all(defined[:, b])),
            "pass": passed,
        }  # fmt: skip
        ok &= passed
    sel = (dose > 0.1 * dose.max()) & np.all(defined, axis=0)
    cover_all = np.mean(np.abs(z[:, sel]) < 1.96, axis=0)
    doc = {
        "step": "a9-compare",
        "energy_mev": ENERGY_MEV,
        "seeds": seeds,
        "frozen_seeds": A9_SEEDS,
        "histories_per_seed": n,
        "batches": A9_BATCHES,
        "coverage": cov,
        "band": list(A9_BAND),
        "r80_mm": r80,
        "informative_all_bins_dose_above_10pct": {
            "bins": int(sel.sum()),
            "coverage_min": float(cover_all.min()) if cover_all.size else math.nan,
            "coverage_mean": float(cover_all.mean()) if cover_all.size else math.nan,
            "coverage_max": float(cover_all.max()) if cover_all.size else math.nan,
        },
        "reference": "leave-one-out pooled ratio of the other seeds (see module docstring)",
        "pass": bool(ok),
    }
    return finish4(doc, A9_HISTORIES, n, n < A9_HISTORIES or seeds < A9_SEEDS)


# -- A13 -----------------------------------------------------------------------------------------
def step_a13(a: argparse.Namespace) -> int:
    n = scaled(A13_HISTORIES, a.scale, 2000, A13_BATCHES)
    dz = 0.5
    geo, grid, nz, r = depth_setup(ENERGY_MEV, dz)
    reqs = (
        TallyRequest("ld", "dose", "let_d", species=("proton",), generation="primary"),
        TallyRequest("lt", "dose", "let_t", species=("proton",), generation="primary"),
    )
    res = base.run_cfg(
        energy=ENERGY_MEV, geometry=geo, scoring=(grid,), seed=a.seed + 13, n=n,
        n_batches=A13_BATCHES, workers=a.workers, timeout=a.timeout, max_step=dz, tallies=reqs,
    )  # fmt: skip
    q = res.grid("dose").quantities
    dose = res.grid("dose").energy_mev.reshape(-1)
    ld, lt = q["ld"], q["lt"]
    m = ld.mean.reshape(-1)
    dfn = ld.defined_mask.reshape(-1)
    z = (np.arange(nz) + 0.5) * dz
    i_max = int(np.nanargmax(np.where(dfn, m, np.nan)))
    sample = [10.0, 0.5 * r, 0.9 * r, r - 10.0, r - 5.0, r - 2.0, r, r + 2.0]
    rows = []
    for zz in sample:
        i = int(zz / dz)
        if 0 <= i < nz:
            rows.append(
                {
                    "depth_mm": float(z[i]),
                    "let_d_kev_um": float(m[i]) if dfn[i] else None,
                    "let_d_std": float(ld.std.reshape(-1)[i]) if dfn[i] else None,
                    "let_t_kev_um": float(lt.mean.reshape(-1)[i])
                    if lt.defined_mask.reshape(-1)[i] else None,
                    "dose_fraction_of_max": float(dose[i] / dose.max()),
                }
            )
    doc = {
        "step": "a13",
        "gating": False,
        "energy_mev": ENERGY_MEV,
        "range_mm": r,
        "bin_mm": dz,
        "let_d_primaries": {
            "distal_maximum_kev_um": float(m[i_max]),
            "distal_maximum_depth_mm": float(z[i_max]),
            "distal_maximum_defined_bins": int(dfn.sum()),
            "profile_samples": rows,
        },
        "literature": GRAZ_NOTE,
        "note": "exploratory and non-gating: the gating LET comparison is TOPAS ProtonLET (V3-010)",
        "pass": bool(counters_clean(res) and np.isfinite(m[i_max])),  # the run itself is valid
    }
    return finish4(doc, A13_HISTORIES, n, n < A13_HISTORIES)


# -- A15 -----------------------------------------------------------------------------------------
def step_a15(a: argparse.Namespace) -> int:
    frozen = A15_HISTORIES[a.backend]
    n = scaled(frozen, a.scale, 4000, A15_BATCHES)
    geo, grid, nz, r = depth_setup(ENERGY_MEV)
    edges = tuple(float(e) for e in np.geomspace(2.0, 160.0, 12))
    kw: dict[str, Any] = {
        "energy": ENERGY_MEV, "geometry": geo, "scoring": (grid,), "seed": a.seed, "n": n,
        "n_batches": A15_BATCHES, "max_step": 1.0, "tallies": all_tallies("dose", edges),
        "lookups": LOOKUPS,
    }  # fmt: skip
    out: dict[str, Any] = {}
    ok = True
    if a.mode == "workers":
        if a.workers < 2:
            raise SystemExit("a15 workers mode needs --workers >= 2")
        for prec in ("float32", "float64"):
            one = base.run_cfg(backend="warp-cpu", precision=prec, workers=1, **kw)
            many = base.run_cfg(
                backend="warp-cpu", precision=prec, workers=a.workers, timeout=a.timeout, **kw
            )
            v = compare_channel_partition(one, many)
            v["workers_used"] = many.transport_report.get("workers")
            out[f"workers_1_vs_{a.workers}_{prec}"] = v
            ok &= bool(v["pass"] and counters_clean(one) and counters_clean(many))
    else:
        for prec in ("float32", "float64"):
            small = base.run_cfg(
                backend=a.backend, precision=prec, chunk=2**10, workers=1, **kw
            )
            large = base.run_cfg(
                backend=a.backend, precision=prec, chunk=2**18, workers=1, **kw
            )
            v = compare_channel_partition(small, large)
            ns = small.transport_report["partials"][0]["n_chunks"]
            nl = large.transport_report["partials"][0]["n_chunks"]
            assert small.channel_raw is not None
            v.update(
                n_chunks_small=ns,
                n_chunks_large=nl,
                different_partitions=bool(ns > nl),
                accumulator_nonzero_entries=int(np.count_nonzero(small.channel_raw.acc)),
                channels=len(small.effective_config.channels.channels),  # type: ignore[union-attr]
            )
            out[f"chunks_1024_vs_262144_{a.backend}_{prec}"] = v
            ok &= bool(
                v["pass"] and v["different_partitions"] and v["accumulator_nonzero_entries"] > 0
                and counters_clean(small) and counters_clean(large)
            )  # fmt: skip
    doc = {"step": "a15", "mode": a.mode, "backend": a.backend, "a15": out, "pass": bool(ok)}
    return finish4(doc, frozen, n, n < frozen)


# -- A11 (LV: python vs warp-cpu float64; HR: cpu32 vs cuda32) -----------------------------------
def compare_channels_ci(a: Result, b: Result) -> dict[str, Any]:
    """The A11 CI/LV rule: per batch and voxel ``|a - b| <= 1e-10 max(|a|, |b|) + n_v q_c``; the N
    channels, ``n_nonzero`` and ``defined_mask`` identical; equal out-of-domain counts."""
    plan = a.effective_config.channels
    assert plan is not None and a.channel_raw is not None and b.channel_raw is not None
    n_ci = {c.grid: i for i, c in enumerate(plan.channels) if c.kind == "N"}
    worst: dict[str, float] = {}
    ok = True
    for ci, c in enumerate(plan.channels):
        x = a.channel_batches(ci) * c.quantum
        y = b.channel_batches(ci) * c.quantum
        n_per = a.channel_batches(n_ci[c.grid]) * plan.channels[n_ci[c.grid]].quantum
        nv = np.repeat(n_per, c.size // n_per.shape[1], axis=1)
        bound = 1e-10 * np.maximum(np.abs(x), np.abs(y)) + nv * c.quantum
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(bound > 0.0, np.abs(x - y) / np.where(bound > 0, bound, 1.0), 0.0)
        ratio = np.where((bound == 0.0) & (x != y), np.inf, ratio)
        key = f"{ci}:{c.kind}:grid{c.grid}"
        worst[key] = float(ratio.max())
        ok &= bool(ratio.max() <= 1.0)
        if c.kind == "N":
            ok &= bool(np.array_equal(a.channel_batches(ci), b.channel_batches(ci)))
    qa, qb = a.grids[0].quantities, b.grids[0].quantities
    same_q = set(qa) == set(qb) and all(
        np.array_equal(qa[k].n_nonzero, qb[k].n_nonzero)
        and np.array_equal(qa[k].defined_mask, qb[k].defined_mask)
        for k in qa
    )
    ood = (a.channel_raw.lookup_out_of_domain, b.channel_raw.lookup_out_of_domain)
    resid = bool(
        np.allclose(a.channel_raw.residual, b.channel_raw.residual, rtol=1e-9, atol=1e-18)
    )
    return {
        "worst_over_bound_by_channel": worst,
        "max_over_bound": max(worst.values()),
        "n_nonzero_and_defined_identical": same_q,
        "lookup_out_of_domain": list(ood),
        "residual_close": resid,
        "pass": bool(ok and same_q and ood == (0, 0) and resid),
    }


def step_a11_lv(a: argparse.Namespace) -> int:
    k = scaled(A11_LV_K, a.scale, 8, 2)
    shape = (12, 12, 40)
    geo = VoxelGeometry(
        (-30.0, -30.0, 0.0), (5.0, 5.0, 5.0), shape, (WATER,), np.zeros(shape, dtype=np.int32)
    )
    grid = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 100), name="dose")
    edges = tuple(float(e) for e in np.geomspace(2.0, 160.0, 12))
    common: dict[str, Any] = {
        "energy": ENERGY_MEV, "geometry": geo, "scoring": (grid,), "seed": a.seed + 11, "n": k,
        "n_batches": 2, "max_step": 2.0, "precision": "float64",
        "tallies": all_tallies("dose", edges), "lookups": LOOKUPS,
    }  # fmt: skip
    t0 = time.perf_counter()
    ref = base.run_cfg(backend="python", workers=1, timeout=a.timeout, **common)
    t_py = time.perf_counter() - t0
    t0 = time.perf_counter()
    wrp = base.run_cfg(backend="warp-cpu", **common)
    t_wp = time.perf_counter() - t0
    v = compare_channels_ci(ref, wrp)
    v["python_seconds"], v["warp_cpu_seconds"] = t_py, t_wp
    v["counters_python"], v["counters_warp"] = ref.counters.as_dict(), wrp.counters.as_dict()
    v["pass"] = bool(v["pass"] and counters_clean(ref) and counters_clean(wrp))
    doc = {"step": "a11", "tier": "lv", "k": k, "energy_mev": ENERGY_MEV, "a11": v,
           "pass": v["pass"]}  # fmt: skip
    return finish4(doc, A11_LV_K, k, k < A11_LV_K)


HR_EDGES = tuple(float(e) for e in np.geomspace(2.0, 160.0, 16))


def hr_sample_config(name: str, a: argparse.Namespace) -> SimulationConfig:
    """Sample ``name`` (``base.SAMPLES``: backend, precision, frozen count, batches, seed index) in
    the T12 configuration (0.2 mm steps, all physics on) with the depth grid ``idd`` only and the
    linear channels of A11 on it."""
    backend, prec, _, batches, kidx = base.SAMPLES[name]
    n = base.sample_histories(name, a.scale)
    cfg, _ = parity.t12_config(
        energy_mev=a.energy, backend=backend, precision=prec, seed=a.seed + 1000 * kidx,
        n_histories=n, n_batches=batches, lateral_bin_mm=0.2, half_width_mm=20.0, workers=1,
        timeout_s=a.timeout,
    )  # fmt: skip
    tallies = (
        TallyRequest("lt", "idd", "let_t"),
        TallyRequest("ld", "idd", "let_d"),
        TallyRequest("le", "idd", "let_d_eps"),
        TallyRequest("fe", "idd", "lookup_sum", lookup=LK_LET.name),
        TallyRequest("sp", "idd", "fluence_spectrum", energy_edges_mev_per_u=HR_EDGES),
    )
    return replace(cfg, scoring=(cfg.scoring[0],), tallies=tallies, lookups=(LK_LET,))


def step_a11_hr(a: argparse.Namespace) -> int:
    pairs = [p.split(":") for p in a.pairs.split(",")]
    names = sorted({x for p in pairs for x in p})
    res: dict[str, Result] = {}
    info: dict[str, Any] = {}
    ok = True
    for name in names:
        t0 = time.perf_counter()
        r = Simulation(hr_sample_config(name, a)).run()
        res[name] = r
        info[name] = {
            "backend": r.backend,
            "precision": r.precision,
            "n_histories": r.n_histories,
            "frozen_histories": base.SAMPLES[name][2],
            "n_batches": r.n_batches,
            "seed": r.seed,
            "device": r.device,
            "seconds": time.perf_counter() - t0,
            "counters": r.counters.as_dict(),
            "valid": counters_clean(r),
        }
        ok &= info[name]["valid"]
    out = {}
    for x, y in pairs:
        v = compare_channel_runs(res[x], res[y], grid="idd")
        out[f"{x}_vs_{y}"] = v
        ok &= bool(v["pass"])
    reduced = any(i["n_histories"] < i["frozen_histories"] for i in info.values())
    doc = {
        "step": "a11",
        "tier": "hr",
        "energy_mev": a.energy,
        "samples": info,
        "a11_hr": out,
        "lookup_file_sha256": LK_LET.file_sha256,
        "pass": bool(ok),
        "frozen_histories": {n: i["frozen_histories"] for n, i in info.items()},
        "histories": {n: i["n_histories"] for n, i in info.items()},
        "reduced": reduced,
    }
    return base.emit(doc)


# -- A16 -----------------------------------------------------------------------------------------
def _git(*args: str) -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                           timeout=300)  # fmt: skip
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout.strip() if r.returncode == 0 else None


def _digests(src: Path, tallies: str, specs: str, env_extra: dict[str, str]) -> dict[str, Any]:
    env = dict(os.environ, PYTHONPATH=str(src), PYTHONDONTWRITEBYTECODE="1", **env_extra)
    cmd = [sys.executable, str(HERE / "a16_digest.py"), "--specs", specs,
           "--tallies", tallies]  # fmt: skip
    p = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=3000)
    if p.returncode != 0 or "#DIGEST-BEGIN" not in p.stdout:
        raise SystemExit(f"digest run failed ({src}, tallies={tallies}):\n{p.stderr[-3000:]}")
    doc: dict[str, Any] = json.loads(
        p.stdout.split("#DIGEST-BEGIN", 1)[1].split("#DIGEST-END", 1)[0]
    )
    if not Path(doc["source"]).is_relative_to(src):
        raise SystemExit(f"digest run imported {doc['source']}, not the intended tree {src}")
    return doc


def step_a16(a: argparse.Namespace) -> int:
    full = _git("rev-parse", "--verify", f"{A16_BASELINE}^{{commit}}")
    if full is None:
        raise SystemExit(
            f"A16 needs the baseline commit {A16_BASELINE} in the git repository (fail closed)"
        )
    specs = ",".join(A16_SPECS)
    env_extra = {k: v for k, v in os.environ.items() if k.startswith("WARP")}
    with tempfile.TemporaryDirectory(prefix="a16-baseline-") as tmp:
        tar = subprocess.run(["git", "-C", str(REPO), "archive", "--format=tar", full, "src"],
                             capture_output=True, timeout=300)  # fmt: skip
        if tar.returncode != 0:
            raise SystemExit(f"git archive of {full} failed: {tar.stderr.decode()[-500:]}")
        subprocess.run(["tar", "-x", "-C", tmp], input=tar.stdout, check=True, timeout=300)
        baseline = _digests(Path(tmp) / "src", "none", specs, env_extra)
    current = {
        "no_tallies": _digests(REPO / "src", "none", specs, env_extra),
        "with_tallies": _digests(REPO / "src", "all", specs, env_extra),
    }
    out: dict[str, Any] = {}
    ok = True
    for spec in A16_SPECS:
        b = baseline["digests"][spec]
        entry: dict[str, Any] = {"fields": len(b)}
        for label, doc in current.items():
            c = doc["digests"][spec]
            diff = sorted(k for k in set(b) | set(c) if b.get(k) != c.get(k))
            entry[label] = {"identical": not diff, "differing_fields": diff[:10]}
            ok &= not diff
        out[spec] = entry
    doc2 = {
        "step": "a16",
        "baseline_commit": full,
        "baseline_ref": A16_BASELINE,
        "fields": "per-grid batch energy, energy balance, counters, end state and trace digests",
        "specs": out,
        "pass": bool(ok),
    }
    return finish4(doc2, None, None, False)


STEPS = {
    "a7": step_a7,
    "a8": step_a8,
    "a9-part": step_a9_part,
    "a9-compare": step_a9_compare,
    "a11-lv": step_a11_lv,
    "a11-hr": step_a11_hr,
    "a13": step_a13,
    "a15": step_a15,
    "a16": step_a16,
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--workers", default="1", help="integer or 'auto' (all cores)")
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--part", default="1/2", help="a9-part: seed range i/n")
    ap.add_argument("--out-dir", default="samples", help="a9-part: directory of the part files")
    ap.add_argument("--dirs", nargs="+", default=["."], help="a9-compare: archive directories")
    ap.add_argument("--seeds", type=int, default=A9_SEEDS, help="a9: number of seeds")
    ap.add_argument("--trace-histories", type=int, default=A7_TRACE_HISTORIES)
    ap.add_argument("--energy", type=float, default=ENERGY_MEV, help="a11-hr: beam energy")
    ap.add_argument("--backend", default="warp-cpu", help="a15: warp-cpu or warp-cuda")
    ap.add_argument("--mode", choices=("workers", "chunks"), default="chunks", help="a15")
    ap.add_argument("--pairs", default="cpu32:cuda32,python:cpu64", help="a11-hr: sample pairs")
    ap.add_argument(
        "--seed-base",
        "--seed",
        dest="seed",
        type=int,
        default=QUALIFICATION_SEED_BASE,
        help="base of every seed of the step (derivation: see the module docstring)",
    )
    ap.add_argument("--timeout", type=float, default=None)
    args = ap.parse_args(argv)
    if not 0.0 < args.scale <= 1.0 or args.seeds < 3:
        raise SystemExit("need 0 < --scale <= 1 and --seeds >= 3")
    base.SEED_BASE = args.seed
    args.workers = (os.cpu_count() or 1) if args.workers == "auto" else int(args.workers)
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
