"""Builder of the derived proton non-elastic nuclear table (decision 0041 sections 1, 3, 6, 7).

``build_nuclear_proton(cache_dir, options)`` reads the hash-pinned sources from the cache
(ENDF/B-VIII.0 proton sublibrary = LA150 evaluations, AME2020 masses, NIST ASTAR water for the D6
gate) and writes ``<cache>/derived/nuclear-proton-<id>.npz`` and ``.json``. Everything is numpy
and stdlib, single process, deterministic: the table id is ``sha256`` of the source hashes, the
builder version and the canonical JSON of the options, and two builds with equal options write
identical bytes.

Cross sections
--------------
Per target nucleus (C-12, N-14, O-16, Al-27, Si-28, P-31, Ca-40) the non-elastic cross section is
ENDF MF3/MT5 on its native TAB1 interpolation up to 150 MeV (zero below the first positive
tabulated value, i.e. below threshold) and, from 150 to 250 MeV, ``sigma(150) * sigma_TL(E) /
sigma_TL(150)`` with the Tripathi light-system shape (``ionmc.physics.tripathi``). The common grid
is uniform in ln E from 1 MeV with 150 MeV as a grid node (spacing ``h = ln 150 / n``,
``n = ceil(points_per_decade log10 150)``) and extends to the first node >= 250 MeV. Elements
without evaluation use a surrogate target and the factor ``(A_el / A_ref)^(2/3)`` (Na, Mg -> Al-27;
S, Cl -> P-31; K, Ar -> Ca-40; the factor of each element is recorded); hydrogen contributes 0.

Product rows (MF6/MT5)
----------------------
For the species n, p, d, alpha (Kalbach-Mann, LANG=2) and gamma (LANG=1) the tables give, at every
grid node: the yield ``y_s(E)`` (TAB1), 64 equiprobable bin edges of E'_CM [MeV] and the
pre-compound fraction ``r`` at the bin midpoints. Each ENDF distribution (histogram, LEP=1) is
converted to its quantile function at ``q = k/64``; between incident energies the 65 quantiles
(and ``r`` per bin) are interpolated linearly in E (quantile interpolation, the inverse of
interpolating the normalised CDF at fixed probability). Above 150 MeV the 150 MeV rows are used
with the edges multiplied by ``E_avail(E) / E_avail(150)``; ``E_avail = sqrt(s) - m_p - M_t`` is the
kinetic energy available in the proton + target centre of mass; yields and ``r`` are held.
Residual recoils (LANG=1 heavy products) are not sampled; the y-weighted sum of their mean energies
is stored per node (information for the V4 energy check). Products t and 3He (present in some
evaluations) are not sampled (recorded in the JSON).

Multiplicities
--------------
Independent Poisson per species (cap 16); the means ``lam_s(E)`` are fixed-point adjusted at
energy nodes so that the post-acceptance mean yields (the sampler of ``ionmc.nuclear.events``, the
numpy implementation of the shared algorithm, with the AME2020 residual-mass test and at most 64
attempts) equal the ENDF yields; ``lam`` is interpolated linearly in ln E between the nodes.
The iteration uses common random numbers (counter-based uniforms, one key per node) so that it
solves a deterministic equation; the residuals are then checked on an independent sample.

D6 gate and Q+
--------------
See :func:`gate_d6` and :func:`q_plus`.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import time
import zipfile
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.data import cache, endf6
from ionmc.data.ame import AmeEntry, load_ame2020
from ionmc.data.endf6 import EndfMaterial, Product, Tab1
from ionmc.data.nist_star import load_star_table
from ionmc.data.registry import DATASETS
from ionmc.materials import ELEMENTS, N_A, WATER, Material
from ionmc.nuclear import events as ev
from ionmc.physics.tripathi import extension_factor

BUILDER_VERSION = "ionmc-nuclear-proton-builder-1"
SCHEMA = "ionmc-nuclear-proton-table-1"
SOURCE_IDS = ("endf-b8.0-protons", "ame2020-mass", "nist-astar-water-2005")
E_MIN_MEV = 1.0
E_ANCHOR_MEV = 150.0
E_MAX_MEV = 250.0
LAMBDA_CAP = 20.0
ENDF_ZAP = {0: 1, 1: 1001, 2: 1002, 3: 2004, 4: 0}
"""ENDF ZAP of the sampled species (n, p, d, alpha, gamma)."""
UNSAMPLED_LIGHT_ZAP = {1003: "t", 2003: "He-3"}
RANGE_THRESHOLD_G_CM2 = 0.01
"""0.1 mm of water (density 1 g/cm3), the D6 range threshold."""
RANGE_TIER2_G_CM2 = 0.2
"""2 mm of water, the tier-2 limit on the 99.9th-percentile alpha range."""
PERCENTILE = 99.9
PILOT_EVENTS = 5000
BISECTION_STEPS = 4


@dataclass(frozen=True)
class TargetSpec:
    """A target nucleus of the table: name, Z, A and the ENDF member file stem."""

    name: str
    z: int
    a: int
    member: str


TARGETS: tuple[TargetSpec, ...] = (
    TargetSpec("C-12", 6, 12, "p-006_C_012"),
    TargetSpec("N-14", 7, 14, "p-007_N_014"),
    TargetSpec("O-16", 8, 16, "p-008_O_016"),
    TargetSpec("Al-27", 13, 27, "p-013_Al_027"),
    TargetSpec("Si-28", 14, 28, "p-014_Si_028"),
    TargetSpec("P-31", 15, 31, "p-015_P_031"),
    TargetSpec("Ca-40", 20, 40, "p-020_Ca_040"),
)
ELEMENT_TARGET: dict[str, tuple[str, bool]] = {
    "C": ("C-12", False),
    "N": ("N-14", False),
    "O": ("O-16", False),
    "Al": ("Al-27", False),
    "P": ("P-31", False),
    "Ca": ("Ca-40", False),
    "Na": ("Al-27", True),
    "Mg": ("Al-27", True),
    "S": ("P-31", True),
    "Cl": ("P-31", True),
    "K": ("Ca-40", True),
    "Ar": ("Ca-40", True),
}
"""Element -> (target name, surrogate flag). Hydrogen contributes 0 and is not listed; every other
element is absent from the table (the runtime raises ``UnsupportedCombinationError``)."""
DEFAULT_NODES_MEV = (2.0, 3.0, 5.0, 7.0, 10.0, 14.0, 20.0, 30.0, 45.0, 60.0, 80.0, 100.0, 125.0,
                     150.0, 175.0, 200.0, 225.0, 250.0)  # fmt: skip


@dataclass(frozen=True)
class BuildOptions:
    """Builder options (their canonical JSON enters the table id)."""

    points_per_decade: int = 100
    lambda_events: int = 200_000
    lambda_seed: int = 20450731
    lambda_max_iterations: int = 8
    lambda_tolerance: float = 0.01
    max_exhausted_fraction: float = 0.01
    lambda_nodes_mev: tuple[float, ...] | None = None
    d6_events: int = 200_000
    d6_seed: int = 20450801
    strict: bool = True

    def canonical(self) -> str:
        """Canonical JSON (sorted keys, no spaces) of the options."""
        return json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))


@dataclass(frozen=True)
class BuildResult:
    """Paths and id of a built table, with the JSON content (``info``)."""

    npz_path: Path
    json_path: Path
    table_id: str
    info: dict[str, Any] = field(default_factory=dict, compare=False)


class BuildError(RuntimeError):
    """The build failed (a lambda node did not converge, or the sources are unusable)."""


def table_id(source_hashes: dict[str, str], options: BuildOptions) -> str:
    """``sha256`` over the source hashes, the builder version and the canonical options."""
    payload = json.dumps(
        {"sources": dict(sorted(source_hashes.items())), "builder": BUILDER_VERSION},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256((payload + "\n" + options.canonical()).encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------------------------
# grid
# ---------------------------------------------------------------------------------------------
def build_grid(points_per_decade: int) -> tuple[NDArray[np.float64], int]:
    """Uniform ln E grid ``(energies [MeV], index of the 150 MeV node)``: ``E_k = exp(k h)``,
    ``h = ln 150 / n``, ``n = ceil(points_per_decade log10 150)``, up to the first node >= 250."""
    if points_per_decade < 50:
        raise ValueError("points_per_decade must be >= 50")
    n = math.ceil(points_per_decade * math.log10(E_ANCHOR_MEV))
    h = math.log(E_ANCHOR_MEV) / n
    n_above = math.ceil(math.log(E_MAX_MEV / E_ANCHOR_MEV) / h)
    e = np.exp(h * np.arange(n + n_above + 1))
    e[0] = E_MIN_MEV
    e[n] = E_ANCHOR_MEV
    return e, n


# ---------------------------------------------------------------------------------------------
# ENDF distributions -> equiprobable bins
# ---------------------------------------------------------------------------------------------
def quantile_edges(
    x_mev: NDArray[np.float64], f: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """Histogram pdf (``f[j]`` on ``[x[j], x[j+1])``, the last ``f`` ignored): the 65 edges of the
    64 equiprobable bins, the segment weights (normalised) and the total unnormalised weight."""
    dx = np.diff(x_mev)
    w = f[:-1] * dx
    total = float(w.sum())
    if not total > 0.0:
        raise BuildError("an ENDF energy distribution has no positive weight")
    cum = np.concatenate(([0.0], np.cumsum(w) / total))
    q = np.arange(N_BINS + 1) / N_BINS
    j = np.clip(np.searchsorted(cum, q, side="right") - 1, 0, w.size - 1)
    pos = np.nonzero(w > 0.0)[0]
    j = np.where(w[j] > 0.0, j, pos[np.clip(np.searchsorted(pos, j), 0, pos.size - 1)])
    seg = (cum[j + 1] - cum[j]) / np.where(dx[j] > 0.0, dx[j], 1.0)
    edges = x_mev[j] + (q - cum[j]) / np.where(seg > 0.0, seg, 1.0)
    edges[0] = x_mev[pos[0]]
    edges[-1] = x_mev[pos[-1] + 1]
    return np.maximum.accumulate(edges), w / total, total


N_BINS = ev.N_BINS


@dataclass(frozen=True)
class _Dist:
    energy_mev: float
    edges: NDArray[np.float64]
    r: NDArray[np.float64]
    mean_mev: float


class SpeciesTables:
    """MF6 data of one product: yield TAB1 and, per incident energy, edges, r and the mean E'."""

    def __init__(self, product: Product) -> None:
        if product.lep != 1:
            raise BuildError(f"LEP={product.lep} is not supported (histogram LEP=1 only)")
        if np.any(product.energy_interp != 2):
            raise BuildError("incident-energy interpolation other than lin-lin is not supported")
        self.zap = product.zap
        self.yield_tab: Tab1 = product.yield_
        dists: list[_Dist] = []
        self.max_norm_dev = 0.0
        for d in product.distributions:
            x = d.rows[:, 0] * 1.0e-6
            f = d.rows[:, 1]
            edges, w, total_ev = quantile_edges(x, f)
            self.max_norm_dev = max(self.max_norm_dev, abs(total_ev * 1.0e6 - 1.0))
            mid = 0.5 * (edges[:-1] + edges[1:])
            if product.lang == 2:
                seg = np.clip(np.searchsorted(x, mid, side="right") - 1, 0, x.size - 2)
                r = np.clip(d.rows[seg, 2], 0.0, 1.0)
            else:
                r = np.zeros(N_BINS)
            mean = float(np.sum(w * 0.5 * (x[:-1] + x[1:])))
            dists.append(_Dist(d.energy * 1.0e-6, edges, r, mean))
        self.dists = dists
        self.energies = np.array([d.energy_mev for d in dists])
        self.edges_all = np.array([d.edges for d in dists])
        self.r_all = np.array([d.r for d in dists])
        self.mean_all = np.array([d.mean_mev for d in dists])

    def yield_at(self, e_mev: float) -> float:
        """Yield at ``e_mev`` <= 150 (0 below the first tabulated energy)."""
        x = self.yield_tab.x * 1.0e-6
        if e_mev < x[0]:
            return 0.0
        return float(self.yield_tab.interpolate(min(e_mev, x[-1]) * 1.0e6)[0])

    def rows_at(self, e_mev: float) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
        """``(edges (65,), r (64,), mean E')`` at ``e_mev`` <= 150 by quantile interpolation (the
        first/last distribution outside the tabulated range)."""
        en = self.energies
        if e_mev <= en[0]:
            return self.edges_all[0], self.r_all[0], float(self.mean_all[0])
        if e_mev >= en[-1]:
            return self.edges_all[-1], self.r_all[-1], float(self.mean_all[-1])
        i = int(np.searchsorted(en, e_mev, side="right") - 1)
        t = (e_mev - en[i]) / (en[i + 1] - en[i])
        return (
            (1.0 - t) * self.edges_all[i] + t * self.edges_all[i + 1],
            (1.0 - t) * self.r_all[i] + t * self.r_all[i + 1],
            float((1.0 - t) * self.mean_all[i] + t * self.mean_all[i + 1]),
        )


def e_avail_mev(model: ev.EventModel, t_lab_mev: float) -> float:
    """Kinetic energy available in the p + target centre of mass, ``sqrt(s) - m_p - M_t``."""
    return (
        ev.cm_boost_np(t_lab_mev, model.m_p_mev, model.m_t_mev)[2] - model.m_p_mev - model.m_t_mev
    )


class TargetTables:
    """All build-time data of one target: cross section, species tables, recoils, event model."""

    def __init__(
        self,
        spec: TargetSpec,
        material: EndfMaterial,
        model: ev.EventModel,
    ) -> None:
        self.spec = spec
        self.model = model
        self.sigma_tab = material.cross_section(5)
        sec = material.products(5)
        if sec.lct != 3:
            raise BuildError(f"{spec.name}: LCT={sec.lct}, only LCT=3 is supported")
        self.species: list[SpeciesTables] = []
        by_zap = {p.zap: p for p in sec.products}
        for s in range(ev.N_SPECIES):
            prod = by_zap.get(ENDF_ZAP[s])
            if prod is None or prod.law != 1:
                raise BuildError(f"{spec.name}: no LAW=1 product with ZAP={ENDF_ZAP[s]}")
            if (s == 4) != (prod.lang == 1) or (s < 4 and prod.lang != 2):
                raise BuildError(f"{spec.name}: unexpected LANG of ZAP={ENDF_ZAP[s]}")
            self.species.append(SpeciesTables(prod))
        self.unsampled_light = {
            UNSAMPLED_LIGHT_ZAP[p.zap]: float(p.yield_.interpolate(150.0e6)[0])
            for p in sec.products
            if p.zap in UNSAMPLED_LIGHT_ZAP
        }
        self.recoils = [
            SpeciesTables(p)
            for p in sec.products
            if p.lang == 1 and p.zap != 0 and p.zap not in UNSAMPLED_LIGHT_ZAP
        ]
        x, y = self.sigma_tab.x * 1.0e-6, self.sigma_tab.y
        i0 = int(np.argmax(y > 0.0))
        self.threshold_mev = float(x[max(i0 - 1, 0)])
        self.e_avail_150 = e_avail_mev(model, E_ANCHOR_MEV)
        self.max_norm_dev = max(s.max_norm_dev for s in self.species)

    def sigma_native_barn(self, e_mev: NDArray[np.float64]) -> NDArray[np.float64]:
        """MT5 on its native interpolation, ``e_mev`` <= 150 (0 below the first tabulated
        energy)."""
        x = self.sigma_tab.x * 1.0e-6
        e = np.asarray(e_mev, dtype=np.float64)
        out = np.zeros(e.shape)
        ok = e >= x[0]
        out[ok] = self.sigma_tab.interpolate(
            np.clip(e[ok] * 1.0e6, self.sigma_tab.x[0], self.sigma_tab.x[-1])
        )
        return out

    def sigma_barn(self, e_mev: NDArray[np.float64]) -> NDArray[np.float64]:
        """Cross section [b] at any energy up to 250 MeV: native up to 150, ``sigma(150)
        * sigma_TL(E) / sigma_TL(150)`` above."""
        e = np.asarray(e_mev, dtype=np.float64)
        s150 = float(self.sigma_native_barn(np.array([E_ANCHOR_MEV]))[0])
        out = self.sigma_native_barn(np.minimum(e, E_ANCHOR_MEV))
        hi = e > E_ANCHOR_MEV
        if np.any(hi):
            out[hi] = s150 * extension_factor(e[hi], self.spec.z, self.spec.a)
        return out

    def rows_at(self, e_mev: float) -> dict[str, Any]:
        """Rows at ``e_mev`` (<= 250): ``yield`` (5,), ``edges`` (5, 65), ``r`` (5, 64),
        ``mean_ecm`` (5,) and ``recoil_energy`` (sum over recoils of yield x mean energy)."""
        e_use = min(e_mev, E_ANCHOR_MEV)
        stretch = 1.0
        if e_mev > E_ANCHOR_MEV:
            stretch = e_avail_mev(self.model, e_mev) / self.e_avail_150
        y = np.zeros(ev.N_SPECIES)
        edges = np.zeros((ev.N_SPECIES, N_BINS + 1))
        r = np.zeros((ev.N_SPECIES, N_BINS))
        mean = np.zeros(ev.N_SPECIES)
        for s, sp in enumerate(self.species):
            y[s] = sp.yield_at(e_use)
            edges[s], r[s], mean[s] = sp.rows_at(e_use)
        recoil = sum((rc.yield_at(e_use) * rc.rows_at(e_use)[2] for rc in self.recoils), 0.0)
        return {
            "yield": y,
            "edges": edges * stretch,
            "r": r,
            "mean_ecm": mean * stretch,
            "recoil_energy": recoil * stretch,
        }


# ---------------------------------------------------------------------------------------------
# lambda adjustment
# ---------------------------------------------------------------------------------------------
def _node_key(seed: int, target: int, node: int, phase: int) -> int:
    ss = np.random.SeedSequence([seed, target, node, phase])
    return int(np.random.PCG64(ss).random_raw())


def batch_means(
    model: ev.EventModel,
    rows: ev.EnergyRows,
    e_mev: float,
    n_events: int,
    key: int,
) -> dict[str, Any]:
    """Post-acceptance mean multiplicities of ``n_events`` events with their standard errors,
    the accepted fraction and the attempt statistics."""
    b = ev.sample_events(model, rows, e_mev, n_events, ev.CounterUniforms(key))
    acc = b.accepted
    n_acc = int(acc.sum())
    c = b.counts[acc].astype(np.float64)
    mean = c.mean(axis=0) if n_acc else np.zeros(ev.N_SPECIES)
    sem = c.std(axis=0) / math.sqrt(max(n_acc, 1)) if n_acc else np.zeros(ev.N_SPECIES)
    total_attempts = int(b.attempts.sum())
    return {
        "mean": mean,
        "sem": sem,
        "n_accepted": n_acc,
        "accepted_fraction": n_acc / n_events,
        "mean_attempts": total_attempts / n_events,
        "rejection_fraction": 1.0 - n_acc / max(total_attempts, 1),
    }


def adjust_lambda_node(
    model: ev.EventModel,
    rows0: dict[str, Any],
    e_mev: float,
    options: BuildOptions,
    keys: tuple[int, int],
) -> dict[str, Any]:
    """Fixed-point adjustment of the Poisson means at one (target, energy) node.

    Iteration ``k`` evaluates the post-acceptance means ``<n_s>_k`` with the common random
    numbers of ``keys[0]`` and updates ``ln lam_s <- ln lam_s - (ln <n_s>_k - ln y_s) / b_s``
    with ``b_s`` the secant slope of ``ln <n_s>`` against ``ln lam_s`` (1 for the first step,
    clipped to [0.25, 1.5]). A species whose expected count ``N y_s`` is below 100 is not resolved
    by the sample and is exempt from the in-sample criterion (recorded).

    A solution is usable only if at most ``options.max_exhausted_fraction`` of the events
    exhaust the 64 attempts (an exhausted event invalidates a transport run, decision 0041
    section 3). A pilot of ``PILOT_EVENTS`` events (the first events of the same stream) checks
    this before each full evaluation. When an iterate leaves the usable region the search stops,
    the last usable iterate is refined by ``BISECTION_STEPS`` bisections towards the unusable one
    and that result is kept with ``converged = False`` if its residuals exceed the tolerance. The
    kept means are finally checked on an independent sample (``keys[1]``):
    ``|<n>/y - 1| <= tol + 3 sem/y``."""
    y = rows0["yield"]
    n = options.lambda_events
    n_pilot = min(n, PILOT_EVENTS)
    tol = options.lambda_tolerance
    min_acc = 1.0 - options.max_exhausted_fraction
    want = y > 0.0
    resolved = want & (n * y >= 100.0)

    def rows_for(lam_v: NDArray[np.float64]) -> ev.EnergyRows:
        return ev.EnergyRows(lam_v.copy(), rows0["edges"], rows0["r"])

    def residual(mean: NDArray[np.float64]) -> NDArray[np.float64]:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(want, mean / np.where(want, y, 1.0) - 1.0, 0.0)

    def pilot_ok(lam_v: NDArray[np.float64]) -> bool:
        return bool(
            batch_means(model, rows_for(lam_v), e_mev, n_pilot, keys[0])["accepted_fraction"]
            >= min_acc
        )

    lam = np.clip(y.copy(), 0.0, LAMBDA_CAP)
    hist: list[dict[str, Any]] = []
    good: dict[str, Any] | None = None
    bad: NDArray[np.float64] | None = None
    prev: tuple[NDArray[np.float64], NDArray[np.float64]] | None = None
    converged = False

    def evaluate(lam_v: NDArray[np.float64], it: int) -> dict[str, Any]:
        full = batch_means(model, rows_for(lam_v), e_mev, n, keys[0])
        res = residual(full["mean"])
        entry = {
            "iteration": it,
            "lam": lam_v.tolist(),
            "residual": res.tolist(),
            "accepted_fraction": full["accepted_fraction"],
            "mean": full["mean"],
        }
        hist.append({k: v for k, v in entry.items() if k != "mean"})
        return entry

    for it in range(1, options.lambda_max_iterations + 1):
        if not pilot_ok(lam):
            bad = lam.copy()
            if good is None:  # even the ENDF yields as means exhaust: record them as they are
                good = evaluate(lam, it)
            break
        cur = evaluate(lam, it)
        if cur["accepted_fraction"] < min_acc:
            bad = lam.copy()
            if good is None:
                good = cur
            break
        good = cur
        res = np.array(cur["residual"])
        if np.all(np.abs(res[resolved]) <= tol):
            converged = True
            break
        if np.any((lam >= LAMBDA_CAP) & resolved & (res < -tol)):
            break  # the Poisson mean is saturated: the yield cannot be reached
        mean = cur["mean"]
        slope = np.ones(ev.N_SPECIES)
        if prev is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                dl = np.log(np.where(lam > 0, lam, 1.0)) - np.log(
                    np.where(prev[0] > 0, prev[0], 1.0)
                )
                dm = np.log(np.where(mean > 0, mean, 1.0)) - np.log(
                    np.where(prev[1] > 0, prev[1], 1.0)
                )
                ok = (np.abs(dl) > 1e-6) & (np.abs(dm) > 1e-9)
                slope = np.where(ok, np.clip(dm / np.where(ok, dl, 1.0), 0.25, 1.5), 1.0)
        prev = (lam.copy(), mean.copy())
        new = lam.copy()
        for s in range(ev.N_SPECIES):
            if not want[s]:
                new[s] = 0.0
            elif mean[s] > 0.0:
                new[s] = lam[s] * math.exp(-math.log(mean[s] / y[s]) / slope[s])
            else:
                new[s] = lam[s] * 2.0
        lam = np.clip(new, 0.0, LAMBDA_CAP)
    assert good is not None
    refined = False
    if not converged and bad is not None:
        lo = np.array(good["lam"])
        hi = bad
        for _ in range(BISECTION_STEPS):
            mid = 0.5 * (lo + hi)
            if pilot_ok(mid):
                lo = mid
            else:
                hi = mid
        if not np.array_equal(lo, np.array(good["lam"])):
            cur = evaluate(lo, len(hist) + 1)
            if cur["accepted_fraction"] >= min_acc:
                good = cur
                refined = True
                converged = bool(np.all(np.abs(np.array(cur["residual"])[resolved]) <= tol))
    lam = np.array(good["lam"])
    res_in = np.array(good["residual"])
    oos = batch_means(model, rows_for(lam), e_mev, n, keys[1])
    res_oos = residual(oos["mean"])
    with np.errstate(divide="ignore", invalid="ignore"):
        sem_rel = np.where(want, oos["sem"] / np.where(want, y, 1.0), 0.0)
    oos_ok = bool(np.all(np.abs(res_oos[want]) <= tol + 3.0 * sem_rel[want]))
    return {
        "e_mev": e_mev,
        "lam": lam.tolist(),
        "yield_endf": y.tolist(),
        "iterations": len(hist),
        "converged": bool(converged),
        "usable_search_limited": bad is not None,
        "bisection_refined": refined,
        "residual_in_sample": res_in.tolist(),
        "unresolved_species": [
            ev.SPECIES[s] for s in range(ev.N_SPECIES) if want[s] and not resolved[s]
        ],
        "residual_independent": res_oos.tolist(),
        "sem_rel_independent": sem_rel.tolist(),
        "independent_ok": oos_ok,
        "lambda_capped": bool(np.any(lam >= LAMBDA_CAP)),
        "max_abs_residual_in_sample": float(np.max(np.abs(res_in[resolved]), initial=0.0)),
        "max_abs_residual_independent": float(np.max(np.abs(res_oos[want]), initial=0.0)),
        "accepted_fraction": oos["accepted_fraction"],
        "mean_attempts": oos["mean_attempts"],
        "rejection_fraction": oos["rejection_fraction"],
        "history": hist,
    }


# ---------------------------------------------------------------------------------------------
# Q+, element mapping, D6 gate
# ---------------------------------------------------------------------------------------------
def q_plus(model: ev.EventModel) -> tuple[float, int]:
    """``(Q+, unreachable)``: ``max(0, max Q)`` over the reachable channels (multiplicity bound
    16 per species, residual in its ground state) and the number of reachable ``(Z_r, A_r)``
    residuals without an AME2020 mass (rejected by the sampler)."""
    q = ev.reaction_q_values(model)
    reach = ev.reachable_residuals(model)
    unreachable = int(np.sum(reach & ~np.isfinite(model.m_res_mev)))
    return max(0.0, float(q.max())), unreachable


def element_rows(targets: tuple[TargetSpec, ...]) -> dict[str, dict[str, Any]]:
    """Per element of the table: target name, surrogate flag, ``A_el`` [g/mol] and the cross
    section factor ``(A_el / A_ref)^(2/3)`` of a surrogate (1 for a direct target)."""
    a_ref = {t.name: t.a for t in targets}
    out: dict[str, dict[str, Any]] = {}
    for sym, (tname, surrogate) in ELEMENT_TARGET.items():
        a_el = ELEMENTS[sym].A_g_mol
        scale = (a_el / a_ref[tname]) ** (2.0 / 3.0) if surrogate else 1.0
        out[sym] = {
            "target": tname,
            "surrogate": surrogate,
            "a_g_mol": a_el,
            "a_ref": a_ref[tname],
            "sigma_scale": scale,
        }
    return out


def alpha_range_g_cm2(
    t_mev: NDArray[np.float64], energy: NDArray[np.float64], csda: NDArray[np.float64]
) -> NDArray[np.float64]:
    """CSDA range [g/cm2] of alphas of total kinetic energy ``t_mev`` from a table (log-log
    interpolation; linear in T below the first node)."""
    t = np.asarray(t_mev, dtype=np.float64)
    out = np.exp(np.interp(np.log(np.maximum(t, 1e-300)), np.log(energy), np.log(csda)))
    return np.where(t < energy[0], csda[0] * t / energy[0], out)


def gate_numbers(
    t_alpha_mev: NDArray[np.float64],
    n_events: int,
    p_event: float,
    e_beam_mev: float,
    range_fn: Callable[[NDArray[np.float64]], NDArray[np.float64]],
) -> dict[str, float]:
    """The D6 numbers of decision 0041 section 7 from the lab alpha energies of ``n_events``
    events: ``G`` = share of the alpha energy carried by alphas with CSDA range > 0.1 mm,
    ``D = (such alpha energy per event) P_event / E_beam``, the 99.9th percentile of the alpha
    energies (count-weighted; the energy-weighted percentile is recorded too) and its range."""
    t = np.asarray(t_alpha_mev, dtype=np.float64)
    r = range_fn(t) if t.size else t
    long = r > RANGE_THRESHOLD_G_CM2
    e_all = float(t.sum())
    e_long = float(t[long].sum())
    g = e_long / e_all if e_all > 0.0 else 0.0
    d = (e_long / n_events) * p_event / e_beam_mev
    out = {
        "G": g,
        "D": d,
        "p_event": p_event,
        "alpha_energy_per_event_mev": e_all / n_events,
        "long_range_alpha_energy_per_event_mev": e_long / n_events,
        "alphas_per_event": t.size / n_events,
        "fraction_alphas_long_range": float(long.mean()) if t.size else 0.0,
    }
    if t.size:
        p = float(np.percentile(t, PERCENTILE))
        order = np.argsort(t)
        cw = np.cumsum(t[order]) / e_all
        pw = float(t[order][min(int(np.searchsorted(cw, PERCENTILE / 100.0)), t.size - 1)])
        out.update(
            {
                "percentile_energy_mev": p,
                "percentile_range_g_cm2": float(range_fn(np.array([p]))[0]),
                "energy_weighted_percentile_energy_mev": pw,
                "energy_weighted_percentile_range_g_cm2": float(range_fn(np.array([pw]))[0]),
            }
        )
    else:
        out.update({"percentile_energy_mev": 0.0, "percentile_range_g_cm2": 0.0,
                    "energy_weighted_percentile_energy_mev": 0.0,
                    "energy_weighted_percentile_range_g_cm2": 0.0})  # fmt: skip
    return out


def tiers(gates: dict[str, dict[str, float]]) -> tuple[bool, bool]:
    """``(tier1_pass, tier2_pass)`` of row D6 from ``{"150": numbers, "250": numbers}``."""
    t1 = all(gates[k]["G"] <= 0.05 for k in ("150", "250"))
    t2 = all(
        gates[k]["D"] <= 1.0e-3 and gates[k]["percentile_range_g_cm2"] <= RANGE_TIER2_G_CM2
        for k in ("150", "250")
    )
    return bool(t1), bool(t2)


def material_sigma_mass(
    material: Material,
    sigma_by_target: dict[str, NDArray[np.float64]],
    elements: dict[str, dict[str, Any]],
) -> NDArray[np.float64]:
    """``Sigma_mass = N_A sum_el w_el sigma_el / A_el`` [cm2/g] from per-target cross sections in
    barn (``1 b = 1e-24 cm2``); hydrogen contributes 0, any other missing element raises."""
    total: NDArray[np.float64] | None = None
    for sym, w in material.mass_fractions.items():
        if sym == "H":
            continue
        row = elements.get(sym)
        if row is None:
            raise KeyError(f"element {sym} is not in the nuclear table")
        term = w * row["sigma_scale"] * sigma_by_target[row["target"]] / row["a_g_mol"]
        total = term if total is None else total + term
    if total is None:
        return np.zeros_like(next(iter(sigma_by_target.values())))
    return np.asarray(N_A * 1.0e-24 * total)


# ---------------------------------------------------------------------------------------------
# deterministic npz
# ---------------------------------------------------------------------------------------------
def write_npz_deterministic(path: Path, arrays: dict[str, NDArray[Any]]) -> str:
    """Write an uncompressed ``.npz`` (sorted members, fixed timestamps, no pickle) and return
    the SHA-256 of the file."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_STORED) as zf:
        for name in sorted(arrays):
            member = io.BytesIO()
            np.lib.format.write_array(
                member, np.ascontiguousarray(arrays[name]), allow_pickle=False
            )
            info = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED
            info.external_attr = 0o600 << 16
            zf.writestr(info, member.getvalue())
    data = buf.getvalue()
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------------------------
def lambda_nodes_for(target: TargetTables, options: BuildOptions) -> list[float]:
    """Energy nodes [MeV] of the lambda adjustment of a target: the default list (or the option)
    restricted to ``[threshold, 250]``, with 150 and 250 MeV always included."""
    base = options.lambda_nodes_mev if options.lambda_nodes_mev is not None else DEFAULT_NODES_MEV
    nodes = {float(e) for e in base if target.threshold_mev <= e <= E_MAX_MEV}
    nodes |= {E_ANCHOR_MEV, E_MAX_MEV}
    return sorted(nodes)


def _noop(_: str) -> None:
    return None


def build_nuclear_proton(
    cache_dir: str | Path | None = None,
    options: BuildOptions | None = None,
    log: Callable[[str], None] = _noop,
) -> BuildResult:
    """Build the table (module docstring); returns the paths and the id. Raises
    :class:`BuildError` when ``options.strict`` and a lambda node does not converge (the files
    are then not written). The D6 gate never raises: its numbers and tier booleans go to the
    JSON."""
    t_start = time.perf_counter()
    opt = options or BuildOptions()
    cdir = cache.resolve_cache_dir(cache_dir)
    zip_path = cache.verify("endf-b8.0-protons", cdir)
    ame_path = cache.verify("ame2020-mass", cdir)
    astar_path = cache.verify("nist-astar-water-2005", cdir)
    sources = {sid: DATASETS[sid].sha256 for sid in SOURCE_IDS}
    tid = table_id(sources, opt)
    ame_tab: dict[tuple[int, int], AmeEntry] = load_ame2020(ame_path.read_text(encoding="ascii"))
    astar = load_star_table(astar_path)

    grid, i150 = build_grid(opt.points_per_decade)
    n_grid = grid.size
    ln_grid = np.log(grid)
    n_t = len(TARGETS)
    elements = element_rows(TARGETS)

    arr: dict[str, NDArray[Any]] = {
        "grid_e_mev": grid,
        "sigma_barn": np.zeros((n_t, n_grid)),
        "yield_endf": np.zeros((n_t, ev.N_SPECIES, n_grid)),
        "lam": np.zeros((n_t, ev.N_SPECIES, n_grid)),
        "edges_mev": np.zeros((n_t, ev.N_SPECIES, n_grid, N_BINS + 1)),
        "r_pre": np.zeros((n_t, ev.N_SPECIES, n_grid, N_BINS)),
        "mean_ecm_mev": np.zeros((n_t, ev.N_SPECIES, n_grid)),
        "recoil_energy_mev": np.zeros((n_t, n_grid)),
        "target_z": np.array([t.z for t in TARGETS], dtype=np.int64),
        "target_a": np.array([t.a for t in TARGETS], dtype=np.int64),
        "target_mass_mev": np.zeros(n_t),
        "s_a_mev": np.zeros(n_t),
        "s_b_mev": np.zeros((n_t, ev.N_SPECIES)),
        "m_res_mev": np.zeros((n_t, ev.DZ_MAX + 1, ev.DA_MAX + 1)),
    }
    models: list[ev.EventModel] = []
    tts: list[TargetTables] = []
    target_info: list[dict[str, Any]] = []
    lam_nodes: list[tuple[NDArray[np.float64], NDArray[np.float64]]] = []
    lam_info: dict[str, Any] = {}
    q_table = 0.0
    for it, spec in enumerate(TARGETS):
        mat = endf6.parse_endf(
            endf6.read_member(zip_path, f"ENDF-B-VIII.0_protons/{spec.member}.endf")
        )
        model = ev.build_event_model(ame_tab, spec.z, spec.a)
        tt = TargetTables(spec, mat, model)
        models.append(model)
        tts.append(tt)
        qp, unreach = q_plus(model)
        q_table = max(q_table, qp)
        arr["target_mass_mev"][it] = model.m_t_mev
        arr["s_a_mev"][it] = model.s_a_mev
        arr["s_b_mev"][it] = model.s_b_mev
        arr["m_res_mev"][it] = model.m_res_mev
        arr["sigma_barn"][it] = tt.sigma_barn(grid)
        nodes = lambda_nodes_for(tt, opt)
        node_lam = np.zeros((len(nodes), ev.N_SPECIES))
        node_results = []
        for ie, e in enumerate(nodes):
            t0 = time.perf_counter()
            rows0 = tt.rows_at(e)
            res = adjust_lambda_node(
                model,
                rows0,
                e,
                opt,
                (_node_key(opt.lambda_seed, it, ie, 0), _node_key(opt.lambda_seed, it, ie, 1)),
            )
            node_lam[ie] = res["lam"]
            node_results.append(res)
            log(
                f"{spec.name} E={e:g} MeV: converged={res['converged']} "
                f"max|res| in={res['max_abs_residual_in_sample']:.4f} "
                f"indep={res['max_abs_residual_independent']:.4f} "
                f"acc={res['accepted_fraction']:.3f} "
                f"[{time.perf_counter() - t0:.1f} s]"
            )
        lam_nodes.append((np.array(nodes), node_lam))
        lam_info[spec.name] = {
            "nodes_mev": nodes,
            "nodes": node_results,
            "all_converged": all(r["converged"] for r in node_results),
        }
        for k in range(n_grid):
            rows = tt.rows_at(float(grid[k]))
            arr["yield_endf"][it, :, k] = rows["yield"]
            arr["edges_mev"][it, :, k] = rows["edges"]
            arr["r_pre"][it, :, k] = rows["r"]
            arr["mean_ecm_mev"][it, :, k] = rows["mean_ecm"]
            arr["recoil_energy_mev"][it, k] = rows["recoil_energy"]
        for s in range(ev.N_SPECIES):
            arr["lam"][it, s] = np.interp(ln_grid, np.log(nodes), node_lam[:, s])
        target_info.append(
            {
                "name": spec.name,
                "z": spec.z,
                "a": spec.a,
                "endf_member": f"ENDF-B-VIII.0_protons/{spec.member}.endf",
                "mass_mev": model.m_t_mev,
                "threshold_mev": tt.threshold_mev,
                "q_plus_mev": qp,
                "unreachable_residual_nuclides": unreach,
                "kalbach_s_a_mev": model.s_a_mev,
                "kalbach_s_b_mev": {ev.SPECIES[s]: float(model.s_b_mev[s]) for s in range(4)},
                "sigma_150_barn": float(tt.sigma_native_barn(np.array([E_ANCHOR_MEV]))[0]),
                "unsampled_light_product_yield_150_mev": tt.unsampled_light,
                "max_distribution_normalisation_deviation": tt.max_norm_dev,
            }
        )
    arr["species_mass_mev"] = models[0].species_mass_mev.copy()
    arr["m_p_mev"] = np.array(models[0].m_p_mev)

    # ---- D6 gate ---------------------------------------------------------------------------
    from ionmc.physics.projectiles import PROTON
    from ionmc.physics.stopping import BetheStoppingSource

    stop = BetheStoppingSource().table(WATER, PROTON)
    sig_by_target = {t.name: arr["sigma_barn"][i] for i, t in enumerate(TARGETS)}
    sigma_water = material_sigma_mass(WATER, sig_by_target, elements)

    def p_event(e_beam: float) -> float:
        lnm = np.linspace(0.0, math.log(e_beam), 4001)
        em = np.exp(lnm)
        sig = np.interp(lnm, ln_grid, sigma_water)
        f = sig / stop.stopping_at(em) * em
        return float(1.0 - math.exp(-np.trapezoid(f, lnm)))

    def range_fn(t: NDArray[np.float64]) -> NDArray[np.float64]:
        return alpha_range_g_cm2(t, astar.energy_mev, astar.csda_range)

    gates: dict[str, dict[str, float]] = {}
    for ie, e_beam in enumerate((150.0, 250.0)):
        weights: dict[int, float] = {}
        for sym, w in WATER.mass_fractions.items():
            if sym == "H":
                continue
            row = elements[sym]
            ti = next(i for i, t in enumerate(TARGETS) if t.name == row["target"])
            sg = float(tts[ti].sigma_barn(np.array([e_beam]))[0])
            weights[ti] = weights.get(ti, 0.0) + w * row["sigma_scale"] * sg / row["a_g_mol"]
        tot_w = sum(weights.values())
        t_alpha: list[NDArray[np.float64]] = []
        n_acc = 0
        for ti, wt in sorted(weights.items()):
            n_ev = int(round(opt.d6_events * wt / tot_w))
            nodes_e, nodes_l = lam_nodes[ti]
            lam = np.array(
                [
                    np.interp(math.log(e_beam), np.log(nodes_e), nodes_l[:, s])
                    for s in range(ev.N_SPECIES)
                ]
            )
            rows_d6 = tts[ti].rows_at(e_beam)
            b = ev.sample_events(
                models[ti],
                ev.EnergyRows(lam, rows_d6["edges"], rows_d6["r"]),
                e_beam,
                n_ev,
                ev.CounterUniforms(_node_key(opt.d6_seed, ti, ie, 2)),
                want_particles=True,
            )
            sel = b.particle_species == 3
            assert b.particle_lab is not None and b.particle_species is not None
            t_alpha.append(b.particle_lab[sel, 0] - models[ti].species_mass_mev[3])
            n_acc += int(b.accepted.sum())
        g = gate_numbers(np.concatenate(t_alpha), max(n_acc, 1), p_event(e_beam), e_beam, range_fn)
        g["accepted_events"] = float(n_acc)
        gates[f"{e_beam:g}"] = g
        log(
            f"D6 {e_beam:g} MeV: G={g['G']:.4f} D={g['D']:.3e} "
            f"R99.9={g['percentile_range_g_cm2']:.4f} g/cm2"
        )
    tier1, tier2 = tiers(gates)

    # ---- JSON and files --------------------------------------------------------------------
    out_dir = cdir / "derived"
    out_dir.mkdir(parents=True, exist_ok=True)
    npz_path = out_dir / f"nuclear-proton-{tid}.npz"
    json_path = out_dir / f"nuclear-proton-{tid}.json"
    all_ok = all(v["all_converged"] for v in lam_info.values())
    if opt.strict and not all_ok:
        bad = [
            f"{n} E={r['e_mev']:g}"
            for n, v in lam_info.items()
            for r in v["nodes"]
            if not r["converged"]
        ]
        raise BuildError(f"lambda nodes did not converge to {opt.lambda_tolerance:.0%}: {bad}")
    npz_sha = write_npz_deterministic(npz_path, arr)
    info: dict[str, Any] = {
        "schema": SCHEMA,
        "table_id": tid,
        "builder_version": BUILDER_VERSION,
        "options": asdict(opt),
        "sources": {
            sid: {"sha256": DATASETS[sid].sha256, "version": DATASETS[sid].version}
            for sid in SOURCE_IDS
        },
        "npz_sha256": npz_sha,
        "npz_arrays": {k: list(v.shape) for k, v in sorted(arr.items())},
        "units": {
            "grid_e_mev": "MeV (proton kinetic energy), uniform in ln E",
            "sigma_barn": "barn per target nucleus [target, node]",
            "yield_endf": "ENDF mean multiplicity [target, species n p d a g, node]",
            "lam": "Poisson mean after the fixed-point adjustment [target, species, node]",
            "edges_mev": "MeV, 65 edges of the 64 equiprobable E'_CM bins [target, species, node]",
            "r_pre": "Kalbach pre-compound fraction per bin [target, species, node, 64]",
            "mean_ecm_mev": "MeV, ENDF mean E'_CM of the product [target, species, node]",
            "recoil_energy_mev": "MeV, sum over recoils of yield x mean energy (information)",
            "m_res_mev": "MeV, residual ground-state mass by (dz, da), inf if no AME2020 mass",
            "s_a_mev, s_b_mev": "MeV, Kalbach separation energies (systematics formula)",
        },
        "species": list(ev.SPECIES),
        "grid": {
            "e_min_mev": E_MIN_MEV,
            "e_max_mev": float(grid[-1]),
            "anchor_node_mev": E_ANCHOR_MEV,
            "anchor_index": i150,
            "n_points": n_grid,
            "ln_step": float(ln_grid[1] - ln_grid[0]) if n_grid > 1 else 0.0,
            "points_per_decade_requested": opt.points_per_decade,
            "points_per_decade_actual": float(1.0 / (math.log10(grid[1]) - math.log10(grid[0]))),
        },
        "targets": target_info,
        "elements": elements,
        "hydrogen": "non-elastic cross section 0 (p-p elastic belongs to V3-005B)",
        "extension": {
            "method": "sigma(150) x sigma_TL(E)/sigma_TL(150) (Tripathi light system, 150-250 MeV)",
            "stretched_fraction": (
                "above 150 MeV the 150 MeV rows are used with E' multiplied by E_avail(E)/"
                "E_avail(150); E_avail = sqrt(s) - m_p - M_t is the kinetic energy available in "
                "the p + target centre of mass; yields and r are held at their 150 MeV values"
            ),
        },
        "product_interpolation": (
            "ENDF histogram distributions (LEP=1) converted to 64 equiprobable bins; between "
            "incident energies the 65 quantiles and r per bin are interpolated linearly in E "
            "(quantile interpolation)"
        ),
        "kalbach_separation": "systematics-formula",
        "multiplicity": {
            "model": "independent Poisson per species n p d a g, cap 16, <= 64 attempts",
            "lambda_cap": LAMBDA_CAP,
            "seed": opt.lambda_seed,
            "events_per_node_per_iteration": opt.lambda_events,
            "max_iterations": opt.lambda_max_iterations,
            "tolerance": opt.lambda_tolerance,
            "random_numbers": "counter-based (splitmix64), one common key per node and phase "
            "derived from numpy PCG64(SeedSequence([seed, target, node, phase]))",
            "interpolation": "linear in ln E between nodes, constant outside",
            "all_nodes_converged": all_ok,
            "targets": lam_info,
        },
        "q_plus_table_mev": q_table,
        "gate_d6": {
            "stopping_source": "Bethe (BetheStoppingSource defaults, I = 78 eV), water, protons",
            "alpha_range_source": "nist-astar-water-2005 CSDA range, log-log interpolation",
            "range_threshold_g_cm2": RANGE_THRESHOLD_G_CM2,
            "tier2_range_limit_g_cm2": RANGE_TIER2_G_CM2,
            "percentile": PERCENTILE,
            "events": opt.d6_events,
            "seed": opt.d6_seed,
            "numbers": gates,
            "tier1_pass": tier1,
            "tier2_pass": tier2,
        },
    }
    json_path.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    info_out = dict(info)
    info_out["timing_s"] = time.perf_counter() - t_start
    return BuildResult(npz_path, json_path, tid, info_out)
