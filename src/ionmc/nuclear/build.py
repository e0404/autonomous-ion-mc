"""Builder of the derived proton non-elastic nuclear table (decision 0041 sections 1, 3, 6, 7).

``build_nuclear_proton(cache_dir, options)`` reads the hash-pinned sources from the cache
(ENDF/B-VIII.0 proton sublibrary = LA150 evaluations, AME2020 masses, NIST ASTAR water for the D6
gate) and writes ``<cache>/derived/nuclear-proton-<id>.npz`` and ``.json``. Everything is numpy
and stdlib, single process, deterministic: the table id is ``sha256`` of the source hashes, the
builder version and the canonical JSON of the options, and two builds with equal options write
identical bytes.

Grid (amendment 2026-10-07, acceptance Amendment 3)
---------------------------------------------------
The runtime grid is the union of 1 MeV, every ENDF MF3/MT5 node, every MF6/MT5 incident energy
and every MF6/MT5 yield TAB1 node of all targets inside [1, 150] MeV, and a uniform ln E grid
with ``points_per_decade >= 50`` (spacing ``h = ln 150 / n``, ``n = ceil(ppd log10 150)``, 150 MeV
a node) continued to the first node >= 250 MeV; uniform nodes within 1e-12 (relative) of an ENDF
node are replaced by it. All these ENDF functions are lin-lin (INT=2) in [1, 150] MeV (the build
fails otherwise), so the runtime interpolation, lin-lin in E between union nodes, reproduces
sigma, the yields and the quantile rows exactly. Lookup: ``ionmc.physics.nuclear.grid_locate``.

Cross sections
--------------
Per target nucleus (C-12, N-14, O-16, Al-27, Si-28, P-31, Ca-40) the non-elastic cross section is
ENDF MF3/MT5 up to 150 MeV (zero below the first positive tabulated value, i.e. below threshold)
and, from 150 to 250 MeV, ``sigma(150) * sigma_TL(E) / sigma_TL(150)`` with the Tripathi
light-system shape (``ionmc.physics.tripathi``). Elements without evaluation use a surrogate
target and the factor ``(A_el / A_ref)^(2/3)`` (Na, Mg -> Al-27; S, Cl -> P-31; K, Ar -> Ca-40;
the factor of each element is recorded); hydrogen contributes 0.

Product rows (MF6/MT5)
----------------------
For the species n, p, d, alpha (Kalbach-Mann, LANG=2) and gamma (LANG=1) the tables give, at every
grid node: the yield ``y_s(E)`` (TAB1), 64 equiprobable bin edges of E'_CM [MeV] and the
pre-compound fraction ``r`` at the bin midpoints. Each ENDF distribution (histogram, LEP=1) is
converted to its quantile function at ``q = k/64``; between incident energies the 65 quantiles
(and ``r`` per bin) are interpolated linearly in E (quantile interpolation). Above 150 MeV the
150 MeV rows are used with the edges multiplied by ``E_avail(E) / E_avail(150)``;
``E_avail = sqrt(s) - m_p - M_t``; yields and ``r`` are held. The residual receives the ENDF mean
heavy-recoil energy ``sum_r y_r <E_r>`` (stored per node, stretched above 150 MeV).

Multiplicities
--------------
Floor + Bernoulli per species (``ionmc.nuclear.events``). The means ``lam_s(E)`` are solved at
every grid node by the fixed point ``lam <- clip(lam y / ybar_post(lam), 0, 16)`` on the exact
enumeration of the 2^5 outcomes (:func:`solve_lambda`; no Monte Carlo), so that the
post-acceptance mean yields equal ENDF; ``lam`` is interpolated lin-lin in E between nodes.

Fail-closed rules
-----------------
``BuildError`` for any non-INT=2 law in [1, 150] MeV, ``P_accept < 0.5`` at any node or node
midpoint (so the runtime exhaustion probability is below 2^-64), or ``sigma > 0`` with
``sum_s y_s = 0`` at a node.

Diagnostics, D6 and the capacity bound
--------------------------------------
See :func:`diagnostics_block`, :func:`gate_d6` numbers (``gate_numbers``) and
:func:`transport_path_terms`.
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

BUILDER_VERSION = "ionmc-nuclear-proton-builder-4"
SCHEMA = "ionmc-nuclear-proton-table-3"
SOURCE_IDS = ("endf-b8.0-protons", "ame2020-mass", "nist-astar-water-2005")
E_MIN_MEV = 1.0
E_ANCHOR_MEV = 150.0
E_MAX_MEV = 250.0
STOPPING_TABLE_MAX_MEV = 500.0
LAMBDA_CAP = 16.0
LAMBDA_SOLVER_ITERATIONS = 50
LAMBDA_CONVERGED_TOLERANCE = 1.0e-3
P_ACCEPT_MIN = 0.5
CEILING_D = 2.0e-2
CEILING_RANGE_G_CM2 = 3.0
ENDF_ZAP = {0: 1, 1: 1001, 2: 1002, 3: 2004, 4: 0}
"""ENDF ZAP of the sampled species (n, p, d, alpha, gamma)."""
UNSAMPLED_LIGHT_ZAP = {1003: "t", 2003: "He-3"}
RANGE_THRESHOLD_G_CM2 = 0.01
"""0.1 mm of water (density 1 g/cm3), the D6 range threshold."""
RANGE_TIER2_G_CM2 = 0.2
"""2 mm of water, the tier-2 limit on the 99.9th-percentile alpha range."""
PERCENTILE = 99.9


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
    """Builder options (their canonical JSON enters the table id). ``strict`` raises
    :class:`BuildError` when a lambda node does not reach ``LAMBDA_CONVERGED_TOLERANCE``
    (otherwise it is recorded ``converged = False``); ``diagnostic_nodes_mev`` restricts the
    diagnostic nodes (default ``DEFAULT_NODES_MEV``).
    """

    points_per_decade: int = 50
    lambda_tolerance: float = 1.0e-6
    multiplicity_model: str = "floor-bernoulli"
    diagnostic_events: int = 20_000
    diagnostic_seed: int = 20450731
    diagnostic_nodes_mev: tuple[float, ...] | None = None
    d6_events: int = 200_000
    d6_seed: int = 20450801
    strict: bool = False

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
    """The build failed (a fail-closed rule, or the sources are unusable)."""


def table_id(source_hashes: dict[str, str], options: BuildOptions, npz_sha256: str) -> str:
    """``sha256`` over the canonical JSON of the source hashes, the builder version, the options
    and the SHA-256 of the npz bytes. The id therefore authenticates the arrays, including the
    ``qualification`` and ``bounds`` arrays that gate transport (decision 0041 section 5)."""
    payload = json.dumps(
        {
            "sources": dict(sorted(source_hashes.items())),
            "builder": BUILDER_VERSION,
            "options": asdict(options),
            "npz_sha256": npz_sha256,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


QUALIFICATION_FIELDS = (
    "ceiling_pass",
    "ceiling_pass_energy_weighted_range",
    "all_nodes_converged",
    "tier1_pass",
    "tier2_pass",
)
"""Order of the ``qualification`` int8 array of the npz (1 = true)."""

BOUND_SPECIES = ("n", "p", "d", "a", "g")
BOUND_FIELDS = (
    "history_energy_bound_mev",
    "transport_energy_bound_mev",
    "recoil_t_max_mev",
    *(f"n_max_{k}" for k in BOUND_SPECIES),
    *(f"t_lab_max_mev_{k}" for k in BOUND_SPECIES),
)
"""Order of the ``bounds`` float64 array of the npz."""


def bounds_array(
    history_bound: float,
    transport_bound: float,
    recoil_max: float,
    path_terms: dict[str, dict[str, float]],
) -> NDArray[np.float64]:
    """The ``bounds`` array of the npz (:data:`BOUND_FIELDS`)."""
    return np.array(
        [history_bound, transport_bound, recoil_max]
        + [path_terms[k]["n_max"] for k in BOUND_SPECIES]
        + [path_terms[k]["t_lab_max_mev"] for k in BOUND_SPECIES],
        dtype=np.float64,
    )


def bounds_from_array(b: NDArray[np.float64]) -> dict[str, Any]:
    """Inverse of :func:`bounds_array`: the scalar bounds and ``transport_path_bound_terms``."""
    if b.shape != (len(BOUND_FIELDS),):
        raise ValueError(f"bounds array has shape {b.shape}, expected ({len(BOUND_FIELDS)},)")
    v = [float(x) for x in b]
    return {
        "history_energy_bound_mev": v[0],
        "transport_energy_bound_mev": v[1],
        "recoil_t_max_mev": v[2],
        "transport_path_bound_terms": {
            k: {"n_max": v[3 + i], "t_lab_max_mev": v[8 + i]} for i, k in enumerate(BOUND_SPECIES)
        },
    }


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


def union_grid(
    endf_nodes_mev: NDArray[np.float64], points_per_decade: int
) -> tuple[NDArray[np.float64], int]:
    """``(grid, index of 150 MeV)``: 1 MeV, the ENDF nodes inside [1, 150] MeV and the uniform ln E
    grid of :func:`build_grid`; a uniform node within 1e-12 (relative) of an ENDF node is dropped
    in favour of the ENDF node (150 MeV stays exactly 150)."""
    uni, _ = build_grid(points_per_decade)
    endf = np.unique(
        np.concatenate(([E_MIN_MEV, E_ANCHOR_MEV], endf_nodes_mev[
            (endf_nodes_mev >= E_MIN_MEV) & (endf_nodes_mev <= E_ANCHOR_MEV)
        ]))
    )  # fmt: skip
    keep = np.ones(endf.size, dtype=bool)
    keep[1:] = np.diff(endf) > 1.0e-12 * endf[1:]  # collapse ENDF nodes closer than 1e-12
    endf = endf[keep]
    pos = np.searchsorted(endf, uni)
    lo = endf[np.clip(pos - 1, 0, endf.size - 1)]
    hi = endf[np.clip(pos, 0, endf.size - 1)]
    near = (np.abs(uni - lo) <= 1.0e-12 * uni) | (np.abs(uni - hi) <= 1.0e-12 * uni)
    grid = np.unique(np.concatenate((endf, uni[~near])))
    return grid, int(np.searchsorted(grid, E_ANCHOR_MEV))


def check_lin_lin(tab: Tab1, label: str) -> None:
    """``BuildError`` unless every interval of ``tab`` overlapping [1, 150] MeV (x in eV) with a
    positive width is interpolated lin-lin (INT=2)."""
    x = tab.x * 1.0e-6
    n = x.size
    region = np.searchsorted(tab.nbt, np.arange(n - 1) + 2, side="left")
    law = tab.interp[region]
    inside = (x[1:] > E_MIN_MEV) & (x[:-1] < E_ANCHOR_MEV) & (np.diff(x) > 0.0)
    if np.any(inside & (law != 2)):
        bad = sorted({int(v) for v in law[inside & (law != 2)]})
        raise BuildError(f"{label}: interpolation law(s) {bad} other than INT=2 in [1, 150] MeV")


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
        self.extend_below_mev = math.inf
        self.extend_to_mev = -math.inf
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
        """Yield at ``e_mev`` <= 150: 0 below the first tabulated energy. Between
        ``extend_below_mev`` (the MT5 threshold, set by :class:`TargetTables`) and
        ``extend_to_mev`` (the first energy at which the target has a positive total yield) the
        yield of ``extend_to_mev`` is used (constant extrapolation downward)."""
        if self.extend_below_mev <= e_mev < self.extend_to_mev:
            e_mev = self.extend_to_mev
        x = self.yield_tab.x * 1.0e-6
        if e_mev < x[0]:
            return 0.0
        e_ev = float(np.clip(e_mev * 1.0e6, self.yield_tab.x[0], self.yield_tab.x[-1]))
        return float(self.yield_tab.interpolate(e_ev)[0])

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
        cand = np.unique(np.concatenate([sp.yield_tab.x * 1.0e-6 for sp in self.species]))
        cand = cand[cand >= self.threshold_mev]
        self.first_mf6_mev = float(min(sp.yield_tab.x[0] for sp in self.species) * 1.0e-6)
        self.first_yield_mev = next(
            (float(e) for e in cand if sum(sp.yield_at(float(e)) for sp in self.species) > 0.0),
            float(cand[-1]),
        )
        for sp in (*self.species, *self.recoils):
            sp.extend_below_mev = self.threshold_mev
            sp.extend_to_mev = self.first_yield_mev
        self.e_avail_150 = e_avail_mev(model, E_ANCHOR_MEV)
        self.max_norm_dev = max(s.max_norm_dev for s in self.species)
        parts = [self.sigma_tab.x * 1.0e-6]
        check_lin_lin(self.sigma_tab, f"{spec.name} MF3/MT5")
        for sp in (*self.species, *self.recoils):
            check_lin_lin(sp.yield_tab, f"{spec.name} MF6/MT5 yield ZAP={sp.zap}")
            parts += [sp.yield_tab.x * 1.0e-6, sp.energies]
        self.endf_nodes_mev = np.unique(np.concatenate(parts))

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
# lambda solver (exact enumeration)
# ---------------------------------------------------------------------------------------------
def _node_key(seed: int, target: int, node: int, phase: int) -> int:
    ss = np.random.SeedSequence([seed, target, node, phase])
    return int(np.random.PCG64(ss).random_raw())


def solve_lambda(model: ev.EventModel, y: NDArray[np.float64], tolerance: float) -> dict[str, Any]:
    """Multiplicity means ``lam`` with ``ybar_post(lam) = y`` for the post-acceptance mean yields
    ``ybar_post`` of the exact enumeration (:func:`ev.exact_post_acceptance`).

    Fixed point ``lam <- clip(lam y / ybar_post(lam), 0, LAMBDA_CAP)`` from ``lam = y``, at most
    ``LAMBDA_SOLVER_ITERATIONS`` iterations, stopped when ``max |ybar/y - 1| <= tolerance`` over
    the species with ``y > 0`` (a species with ``y = 0`` has ``lam = 0``). The returned ``lam`` is
    the iterate of smallest residual among those with ``P_accept >= P_ACCEPT_MIN``: when the yield
    is not reachable (C-12 near 14-25 MeV: with one proton, three alphas leave no residual, so
    the post-acceptance alpha mean cannot exceed about 2) the unguarded fixed point runs to
    ``P_accept ~ 0.003``; the iteration then stops at the last feasible iterate and the node is
    recorded ``converged = False`` with its residual. Returns ``lam``,
    ``p_accept`` and the post-acceptance yield ratios at the returned ``lam``, the maximum
    absolute residual, the iterations and ``converged`` (residual <= ``LAMBDA_CONVERGED_
    TOLERANCE``)."""
    want = y > 0.0
    lam = np.where(want, np.minimum(y, LAMBDA_CAP), 0.0)
    best = (math.inf, lam.copy())
    p_acc, mean = ev.exact_post_acceptance(model, lam)
    it = 0
    for it in range(LAMBDA_SOLVER_ITERATIONS + 1):
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(want, mean / np.where(want, y, 1.0), 1.0)
        res = float(np.max(np.abs(ratio[want] - 1.0), initial=0.0))
        if p_acc < P_ACCEPT_MIN and it > 0:
            break  # infeasible iterate: the yield is not reachable inside the usable region
        if res < best[0]:
            best = (res, lam.copy())
        if res <= tolerance or it == LAMBDA_SOLVER_ITERATIONS:
            break
        new = np.where(
            want,
            np.where(
                ratio > 0.0, lam / np.where(ratio > 0.0, ratio, 1.0), np.maximum(lam, 1e-3) * 2
            ),
            0.0,
        )
        lam = np.clip(new, 0.0, LAMBDA_CAP)
        p_acc, mean = ev.exact_post_acceptance(model, lam)
    lam = best[1]
    p_acc, mean = ev.exact_post_acceptance(model, lam)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(want, mean / np.where(want, y, 1.0), 1.0)
    return {
        "lam": lam,
        "p_accept": p_acc,
        "yield_ratio": ratio,
        "max_residual": float(np.max(np.abs(ratio[want] - 1.0), initial=0.0)),
        "iterations": it,
        "converged": bool(
            np.max(np.abs(ratio[want] - 1.0), initial=0.0) <= LAMBDA_CONVERGED_TOLERANCE
        ),
    }


PATH_SPECIES = (("n", 0), ("p", 1), ("d", 2), ("a", 3), ("g", 4))
"""Species keys of ``transport_path_bound_terms`` and their index into ``ev.SPECIES``."""


def transport_path_terms(
    grid: NDArray[np.float64],
    lam: NDArray[np.float64],
    edges: NDArray[np.float64],
    models: list[ev.EventModel],
    active: NDArray[np.bool_],
) -> dict[str, dict[str, float]]:
    """Per product species ``s`` in {n, p, d, a, g} (decision 0041 section 5, amended 2026-10-08):
    ``n_max`` = the maximum over targets and grid nodes with ``sigma > 0`` (``active`` (targets,
    nodes)) of ``ceil(lam_s)`` (cap 16) and ``t_lab_max_mev`` = the maximum of the lab kinetic
    energy of a product of species ``s`` with the top bin edge as CM energy emitted along the
    incident direction (``mu = 1``). The path bound is ``B_L = 1.25 (mixed_path_bound(E_hi) +
    sum_s n_max,s mixed_path_bound_s(t_lab_max,s))`` over the transported species p and d; the
    per-history energy bound sums all five species (:func:`history_energy_bound`)."""
    out = {key: {"n_max": 0.0, "t_lab_max_mev": 0.0} for key, _ in PATH_SPECIES}
    for it, model in enumerate(models):
        for k in np.nonzero(active[it])[0]:
            beta, gamma, _ = ev.cm_boost_np(float(grid[k]), model.m_p_mev, model.m_t_mev)
            for key, s in PATH_SPECIES:
                m = float(model.species_mass_mev[s])
                t_cm = float(edges[it, s, k, -1])
                p = math.sqrt(t_cm * (t_cm + 2.0 * m))
                t_lab = gamma * (t_cm + m + beta * p) - m
                out[key]["t_lab_max_mev"] = max(out[key]["t_lab_max_mev"], t_lab)
                n = min(math.ceil(float(lam[it, s, k]) - 1e-12), 16)
                out[key]["n_max"] = max(out[key]["n_max"], float(n))
    return out


def recoil_t_max(recoil_t_cm: NDArray[np.float64], active: NDArray[np.bool_]) -> float:
    """``T_r,max`` [MeV]: the largest heavy-recoil energy ``T_r`` (the ENDF mean, applied as is
    by the event sampler) over the active (target, node) pairs."""
    return float(np.max(np.where(active, recoil_t_cm, 0.0), initial=0.0))


def history_energy_bound(
    terms: dict[str, dict[str, float]], recoil_max_mev: float, e_hi_mev: float
) -> float:
    """Per-history energy bound of the table (decision 0041 section 5, amended 2026-10-08):
    ``max(E_hi, sum over ALL species s of N_s,max T_lab,max,s + T_r,max)``; ``e_hi_mev`` is the
    largest tabulated proton energy of the table (a history never needs less)."""
    total = sum(v["n_max"] * v["t_lab_max_mev"] for v in terms.values()) + recoil_max_mev
    return max(e_hi_mev, total)


def diagnostics_block(
    model: ev.EventModel,
    rows: ev.EnergyRows,
    endf: dict[str, Any],
    e_mev: float,
    n_events: int,
    key: int,
) -> dict[str, Any]:
    """Diagnostics of ``n_events`` seeded events at one (target, energy) with the table rows:
    exact ``p_accept``; post-acceptance yield ratios sampled/ENDF with 3 sem; per-species mean
    ``E'`` ratios sampled/ENDF; the ledger ``imbalance`` (``Delta_lab``: mean, sem, sd,
    ``P(Delta < 0)``), ``Delta_CM = sqrt(s) - sum E_cm - M_r - T_r`` mean, mean ``|sum p_CM|``,
    mean local deposit; ``E_avail``."""
    b = ev.sample_events(model, rows, e_mev, n_events, ev.CounterUniforms(key))
    acc = b.accepted
    n_acc = int(acc.sum())
    c = b.counts[acc].astype(np.float64)
    mean_c = c.mean(axis=0)
    sem_c = c.std(axis=0) / math.sqrt(n_acc)
    y = endf["yield"]
    sum_ep = b.sum_e_prime_mev[acc].sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_y = np.where(y > 0, mean_c / np.where(y > 0, y, 1.0), np.nan)
        sem_y = np.where(y > 0, 3.0 * sem_c / np.where(y > 0, y, 1.0), np.nan)
        e_mean = np.where(
            c.sum(axis=0) > 0, sum_ep / np.where(c.sum(axis=0) > 0, c.sum(axis=0), 1.0), np.nan
        )
        ratio_e = np.where(
            endf["mean_ecm"] > 0,
            e_mean / np.where(endf["mean_ecm"] > 0, endf["mean_ecm"], 1.0),
            np.nan,
        )
    d = b.imbalance_mev[acc]
    p_acc, _ = ev.exact_post_acceptance(model, rows.lam)
    sum_ye_sampled = float(sum_ep.sum() / n_acc + rows.recoil_t_mev)
    sum_ye_endf = float(np.sum(y * endf["mean_ecm"]) + endf["recoil_energy"])
    return {
        "e_mev": e_mev,
        "events": n_events,
        "accepted_fraction_sampled": n_acc / n_events,
        "p_accept_exact": p_acc,
        "e_avail_mev": e_avail_mev(model, e_mev),
        "yield_ratio": ratio_y.tolist(),
        "yield_ratio_3sem": sem_y.tolist(),
        "mean_e_prime_ratio": ratio_e.tolist(),
        "sum_y_e_prime_plus_recoil_mev": sum_ye_sampled,
        "sum_y_e_prime_plus_recoil_endf_mev": sum_ye_endf,
        "sum_y_e_prime_ratio": sum_ye_sampled / sum_ye_endf,
        "delta_lab_mean_mev": float(d.mean()),
        "delta_lab_sem_mev": float(d.std() / math.sqrt(n_acc)),
        "delta_lab_sd_mev": float(d.std()),
        "delta_cm_mean_mev": float(b.imbalance_cm_mev[acc].mean()),
        "p_delta_negative": float((d < 0.0).mean()),
        "mean_abs_sum_p_cm_mev": float(b.p_cm_sum_mev[acc].mean()),
        "mean_local_deposit_mev": float(b.local_deposit_mev[acc].mean()),
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


def ceilings(gates: dict[str, dict[str, float]]) -> tuple[bool, bool]:
    """``(count_weighted, energy_weighted)``: the B1 validity ceiling of the alpha local-deposition
    decision (acceptance Amendment 2): ``D <= CEILING_D`` and the 99.9th-percentile alpha range
    ``<= CEILING_RANGE_G_CM2`` at both energies (the range of the count-weighted, respectively
    the energy-weighted, 99.9th-percentile energy)."""
    d_ok = all(gates[k]["D"] <= CEILING_D for k in ("150", "250"))
    r_c = all(gates[k]["percentile_range_g_cm2"] <= CEILING_RANGE_G_CM2 for k in ("150", "250"))
    r_e = all(
        gates[k]["energy_weighted_percentile_range_g_cm2"] <= CEILING_RANGE_G_CM2
        for k in ("150", "250")
    )
    return bool(d_ok and r_c), bool(d_ok and r_e)


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
def _noop(_: str) -> None:
    return None


def build_nuclear_proton(
    cache_dir: str | Path | None = None,
    options: BuildOptions | None = None,
    log: Callable[[str], None] = _noop,
) -> BuildResult:
    """Build the table (module docstring); returns the paths and the id. Raises
    :class:`BuildError` on a fail-closed rule (the files are then not written) and, with
    ``options.strict``, on a lambda node that did not converge to ``LAMBDA_CONVERGED_TOLERANCE``.
    The D6 gate never raises: its numbers and the tier and ceiling booleans go to the JSON."""
    t_start = time.perf_counter()
    opt = options or BuildOptions()
    if opt.multiplicity_model != "floor-bernoulli":
        raise BuildError(f"unknown multiplicity model {opt.multiplicity_model!r}")
    cdir = cache.resolve_cache_dir(cache_dir)
    zip_path = cache.verify("endf-b8.0-protons", cdir)
    ame_path = cache.verify("ame2020-mass", cdir)
    astar_path = cache.verify("nist-astar-water-2005", cdir)
    sources = {sid: DATASETS[sid].sha256 for sid in SOURCE_IDS}
    ame_tab: dict[tuple[int, int], AmeEntry] = load_ame2020(ame_path.read_text(encoding="ascii"))
    astar = load_star_table(astar_path)

    models: list[ev.EventModel] = []
    tts: list[TargetTables] = []
    for spec in TARGETS:
        mat = endf6.parse_endf(
            endf6.read_member(zip_path, f"ENDF-B-VIII.0_protons/{spec.member}.endf")
        )
        model = ev.build_event_model(ame_tab, spec.z, spec.a)
        models.append(model)
        tts.append(TargetTables(spec, mat, model))
    grid, i150 = union_grid(
        np.unique(np.concatenate([tt.endf_nodes_mev for tt in tts])), opt.points_per_decade
    )
    n_grid = grid.size
    n_t = len(TARGETS)
    elements = element_rows(TARGETS)
    n_endf = int(np.unique(np.concatenate([tt.endf_nodes_mev for tt in tts])).size)
    log(f"grid: {n_grid} nodes ({n_endf} distinct ENDF energies over all targets)")

    arr: dict[str, NDArray[Any]] = {
        "grid_e_mev": grid,
        "sigma_barn": np.zeros((n_t, n_grid)),
        "yield_endf": np.zeros((n_t, ev.N_SPECIES, n_grid)),
        "lam": np.zeros((n_t, ev.N_SPECIES, n_grid)),
        "p_accept": np.zeros((n_t, n_grid)),
        "yield_ratio_post": np.zeros((n_t, ev.N_SPECIES, n_grid)),
        "lam_converged": np.zeros((n_t, n_grid), dtype=np.int8),
        "edges_mev": np.zeros((n_t, ev.N_SPECIES, n_grid, N_BINS + 1)),
        "r_pre": np.zeros((n_t, ev.N_SPECIES, n_grid, N_BINS)),
        "mean_ecm_mev": np.zeros((n_t, ev.N_SPECIES, n_grid)),
        "recoil_t_cm_mev": np.zeros((n_t, n_grid)),
        "target_z": np.array([t.z for t in TARGETS], dtype=np.int64),
        "target_a": np.array([t.a for t in TARGETS], dtype=np.int64),
        "target_mass_mev": np.zeros(n_t),
        "s_a_mev": np.zeros(n_t),
        "s_b_mev": np.zeros((n_t, ev.N_SPECIES)),
        "m_res_mev": np.zeros((n_t, ev.DZ_MAX + 1, ev.DA_MAX + 1)),
    }
    target_info: list[dict[str, Any]] = []
    lam_info: dict[str, Any] = {}
    nonconverged: list[dict[str, Any]] = []
    p_min: dict[str, Any] = {"p_accept": math.inf}
    p_fail: list[str] = []
    no_prod: list[str] = []
    q_table = 0.0
    for it, (spec, model, tt) in enumerate(zip(TARGETS, models, tts, strict=True)):
        qp, unreach = q_plus(model)
        q_table = max(q_table, qp)
        arr["target_mass_mev"][it] = model.m_t_mev
        arr["s_a_mev"][it] = model.s_a_mev
        arr["s_b_mev"][it] = model.s_b_mev
        arr["m_res_mev"][it] = model.m_res_mev
        arr["sigma_barn"][it] = tt.sigma_barn(grid)
        n_iter = 0
        for k in range(n_grid):
            rows = tt.rows_at(float(grid[k]))
            arr["yield_endf"][it, :, k] = rows["yield"]
            arr["edges_mev"][it, :, k] = rows["edges"]
            arr["r_pre"][it, :, k] = rows["r"]
            arr["mean_ecm_mev"][it, :, k] = rows["mean_ecm"]
            arr["recoil_t_cm_mev"][it, k] = rows["recoil_energy"]
            if arr["sigma_barn"][it, k] > 0.0 and not np.sum(rows["yield"]) > 0.0:
                no_prod.append(
                    f"{spec.name} E={grid[k]:.6g} MeV sigma={arr['sigma_barn'][it, k]:.4g} b"
                )
            sol = solve_lambda(model, rows["yield"], opt.lambda_tolerance)
            arr["lam"][it, :, k] = sol["lam"]
            arr["p_accept"][it, k] = sol["p_accept"]
            arr["yield_ratio_post"][it, :, k] = sol["yield_ratio"]
            arr["lam_converged"][it, k] = int(sol["converged"])
            n_iter += sol["iterations"]
            if not sol["converged"]:
                nonconverged.append(
                    {
                        "target": spec.name,
                        "e_mev": float(grid[k]),
                        "max_residual": sol["max_residual"],
                    }
                )
            if sol["p_accept"] < p_min["p_accept"]:
                p_min = {"p_accept": sol["p_accept"], "target": spec.name,
                         "e_mev": float(grid[k]), "kind": "node"}  # fmt: skip
        pa_nodes = arr["p_accept"][it]
        lam_mid = 0.5 * (arr["lam"][it, :, :-1] + arr["lam"][it, :, 1:])
        pa_mid = np.array(
            [ev.exact_post_acceptance(model, lam_mid[:, k])[0] for k in range(n_grid - 1)]
        )
        k_mid = int(np.argmin(pa_mid))
        if pa_mid[k_mid] < p_min["p_accept"]:
            p_min = {"p_accept": float(pa_mid[k_mid]), "target": spec.name,
                     "e_mev": float(0.5 * (grid[k_mid] + grid[k_mid + 1])),
                     "kind": "midpoint"}  # fmt: skip
        if min(float(pa_nodes.min()), float(pa_mid.min())) < P_ACCEPT_MIN:
            kn = int(np.argmin(pa_nodes))
            p_fail.append(
                f"{spec.name}: min node {pa_nodes.min():.3f} at {grid[kn]:.4g} MeV, min midpoint "
                f"{pa_mid.min():.3f} at {0.5 * (grid[k_mid] + grid[k_mid + 1]):.4g} MeV, "
                f"{int(np.sum(pa_nodes < P_ACCEPT_MIN))} nodes and "
                f"{int(np.sum(pa_mid < P_ACCEPT_MIN))} midpoints below {P_ACCEPT_MIN}"
            )
        lam_info[spec.name] = {
            "p_accept_min_nodes": float(pa_nodes.min()),
            "p_accept_min_midpoints": float(pa_mid.min()),
            "solver_iterations_total": n_iter,
            "non_converged_nodes": int(np.sum(arr["lam_converged"][it] == 0)),
            "max_lambda": float(arr["lam"][it].max()),
        }
        log(f"{spec.name}: P_accept min nodes {pa_nodes.min():.4f} midpoints {pa_mid.min():.4f}")
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
                "yield_extended_below_mev": [tt.first_yield_mev, tt.threshold_mev],
                "first_mf6_energy_mev": tt.first_mf6_mev,
                "yield_extended": bool(tt.threshold_mev + 1e-9 < tt.first_yield_mev),
                "endf_nodes_in_range": int(
                    np.sum((tt.endf_nodes_mev >= E_MIN_MEV) & (tt.endf_nodes_mev <= E_ANCHOR_MEV))
                ),
            }
        )
    if p_fail or no_prod:
        raise BuildError(
            "fail-closed: exact P_accept below 0.5: "
            + ("; ".join(p_fail) or "none")
            + " | sigma > 0 without product yields: "
            + (", ".join(no_prod) or "none")
        )
    arr["species_mass_mev"] = models[0].species_mass_mev.copy()
    arr["m_p_mev"] = np.array(models[0].m_p_mev)
    if opt.strict and nonconverged:
        raise BuildError(f"lambda nodes did not converge to 1e-3: {nonconverged[:10]} ...")
    path_terms = transport_path_terms(
        grid, arr["lam"], arr["edges_mev"], models, arr["sigma_barn"] > 0.0
    )
    particle_bound = max(path_terms[k]["t_lab_max_mev"] for k in ("p", "d"))
    if particle_bound > STOPPING_TABLE_MAX_MEV:
        raise BuildError(
            f"per-particle energy bound exceeds {STOPPING_TABLE_MAX_MEV} MeV: {path_terms}"
        )
    recoil_max = recoil_t_max(arr["recoil_t_cm_mev"], arr["sigma_barn"] > 0.0)
    hist_bound = history_energy_bound(path_terms, recoil_max, float(grid[-1]))
    arr["bounds"] = bounds_array(hist_bound, particle_bound, recoil_max, path_terms)

    def rows_of(it: int, e: float) -> ev.EnergyRows:
        return ev.interp_rows(
            grid, arr["lam"][it], arr["edges_mev"][it], arr["r_pre"][it],
            arr["recoil_t_cm_mev"][it], e,
        )  # fmt: skip

    # ---- diagnostics -----------------------------------------------------------------------
    diag_nodes = (
        opt.diagnostic_nodes_mev if opt.diagnostic_nodes_mev is not None else DEFAULT_NODES_MEV
    )
    diagnostics: dict[str, Any] = {}
    for it, (spec, model, tt) in enumerate(zip(TARGETS, models, tts, strict=True)):
        entries = []
        for ie, e in enumerate(diag_nodes):
            if not (tt.threshold_mev <= e <= E_MAX_MEV):
                continue
            entries.append(
                diagnostics_block(
                    model, rows_of(it, float(e)), tt.rows_at(float(e)), float(e),
                    opt.diagnostic_events, _node_key(opt.diagnostic_seed, it, ie, 0),
                )
            )  # fmt: skip
        diagnostics[spec.name] = entries
        log(f"{spec.name}: {len(entries)} diagnostic nodes")

    # ---- D6 gate ---------------------------------------------------------------------------
    from ionmc.physics.projectiles import PROTON
    from ionmc.physics.stopping import BetheStoppingSource

    stop = BetheStoppingSource().table(WATER, PROTON)
    sig_by_target = {t.name: arr["sigma_barn"][i] for i, t in enumerate(TARGETS)}
    sigma_water = material_sigma_mass(WATER, sig_by_target, elements)

    def p_event(e_beam: float) -> float:
        lnm = np.linspace(0.0, math.log(e_beam), 4001)
        em = np.exp(lnm)
        sig = np.interp(em, grid, sigma_water)
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
            b = ev.sample_events(
                models[ti],
                rows_of(ti, e_beam),
                e_beam,
                n_ev,
                ev.CounterUniforms(_node_key(opt.d6_seed, ti, ie, 2)),
                want_particles=True,
            )
            assert b.particle_lab is not None and b.particle_species is not None
            sel = b.particle_species == 3
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
    ceiling, ceiling_energy = ceilings(gates)
    arr["qualification"] = np.array(
        [ceiling, ceiling_energy, bool(np.all(arr["lam_converged"] == 1)), tier1, tier2],
        dtype=np.int8,
    )

    # ---- JSON and files --------------------------------------------------------------------
    out_dir = cdir / "derived"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = out_dir / "nuclear-proton-building.npz"
    npz_sha = write_npz_deterministic(tmp_path, arr)
    tid = table_id(sources, opt, npz_sha)  # the id covers the npz bytes (decision 0041 section 5)
    npz_path = out_dir / f"nuclear-proton-{tid}.npz"
    json_path = out_dir / f"nuclear-proton-{tid}.json"
    tmp_path.replace(npz_path)
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
            "grid_e_mev": "MeV (proton kinetic energy), union grid, interpolated lin-lin in E",
            "sigma_barn": "barn per target nucleus [target, node]",
            "yield_endf": "ENDF mean multiplicity [target, species n p d a g, node]",
            "lam": "floor+Bernoulli mean after the exact solve [target, species, node]",
            "p_accept": "exact probability that the residual exists in one attempt [target, node]",
            "yield_ratio_post": "exact post-acceptance yield / ENDF yield [target, species, node]",
            "lam_converged": "1 iff max |ratio - 1| <= 1e-3 [target, node]",
            "qualification": f"int8 flags {list(QUALIFICATION_FIELDS)} (1 = true); the loader's "
            "gate, the JSON copies are cross-checked",
            "bounds": f"float64 {list(BOUND_FIELDS)}; the loader's capacity bounds",
            "edges_mev": "MeV, 65 edges of the 64 equiprobable E'_CM bins [target, species, node]",
            "r_pre": "Kalbach pre-compound fraction per bin [target, species, node, 64]",
            "mean_ecm_mev": "MeV, ENDF mean E'_CM of the product [target, species, node]",
            "recoil_t_cm_mev": "MeV, ENDF mean heavy-recoil energy [target, node]",
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
            "n_distinct_endf_energies": n_endf,
            "points_per_decade_requested": opt.points_per_decade,
            "interpolation": "lin-lin in E between union nodes; lookup grid_locate",
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
            "model": "floor+Bernoulli per species n p d a g, cap 16, <= 64 attempts, residual-"
            "existence acceptance",
            "lambda_cap": LAMBDA_CAP,
            "solver": "fixed point lam <- clip(lam y / ybar_post(lam), 0, 16), exact enumeration",
            "tolerance": opt.lambda_tolerance,
            "converged_tolerance": LAMBDA_CONVERGED_TOLERANCE,
            "p_accept_min": p_min,
            "non_converged_nodes": nonconverged,
            "all_nodes_converged": not nonconverged,
            "targets": lam_info,
        },
        "empty_residual_allowed": True,
        "transport_energy_bound_mev": particle_bound,
        "history_energy_bound_mev": hist_bound,
        "recoil_t_max_mev": recoil_max,
        "transport_path_bound_terms": path_terms,
        "q_plus_table_mev": q_table,
        "diagnostics": {
            "events_per_node": opt.diagnostic_events,
            "seed": opt.diagnostic_seed,
            "nodes_mev": [float(e) for e in diag_nodes],
            "targets": diagnostics,
        },
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
            "ceiling_d_max": CEILING_D,
            "ceiling_range_g_cm2": CEILING_RANGE_G_CM2,
            "ceiling_pass": ceiling,
            "ceiling_pass_energy_weighted_range": ceiling_energy,
        },
    }
    json_path.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    info_out = dict(info)
    info_out["timing_s"] = time.perf_counter() - t_start
    return BuildResult(npz_path, json_path, tid, info_out)
