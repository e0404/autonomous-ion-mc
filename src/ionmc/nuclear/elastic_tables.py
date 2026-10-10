"""Loader of the derived elastic table and its per-material rows (decision 0041 slice C).

``ElasticTable.load(cache_dir, table_id)`` mirrors :meth:`ionmc.nuclear.tables.NuclearTable.load`:
it reads ``<cache>/derived/elastic-proton-<id>.npz`` and ``.json`` with ``np.load(allow_pickle=
False)``, re-hashes the npz bytes against the JSON, checks that the id is the sha256 of the
canonical
sidecar without its id, that every source pin equals the registry pin, derives the qualification
flags from the npz arrays and freezes them. Failures are fail-closed with specific exception types
(all subclass :class:`ElasticTableError`, itself an ``UnsupportedCombinationError``): missing, stale
(re-hash or id mismatch, schema/builder), source pin mismatch, unqualified (a qualification flag is
false). The per-target lower domain bound ``e_min`` (``elastic_domain``) is the declared limitation
of the table: below it ``sigma = 0`` (H-1: ``E_min,pp``; p + A: ``e_min_shape``); the transport
counts such crossings as ``elastic_below_domain``. An element without elastic data and an energy
outside [1, 250] MeV raise
``UnsupportedCombinationError``.

``material_rows`` composes ``Sigma_mass(E) = N_A sum_el w_el scale_el sigma_el / A_el`` [cm2/g]
(hydrogen contributes through the p-p target) and the cumulative partial ``Sigma`` per target with
the majorants of :func:`ionmc.nuclear.tables.majorant`, on the table grid, lin-lin in E.
"""

from __future__ import annotations

import hashlib
import io
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc._frozen import freeze_array, freeze_json
from ionmc.data import cache
from ionmc.data.registry import DATASETS
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import N_A, Material
from ionmc.nuclear.build import table_id as compute_table_id
from ionmc.nuclear.elastic_build import BUILDER_VERSION, QUALIFICATION_FIELDS, SCHEMA
from ionmc.nuclear.tables import DEFAULT_F_E, interp_row, locate, majorant

E_MIN_MEV = 1.0
E_MAX_MEV = 250.0


class ElasticTableError(UnsupportedCombinationError, RuntimeError):
    """Base class of elastic-table load failures."""


class ElasticTableMissingError(ElasticTableError, FileNotFoundError):
    """The table files are not in the cache."""


class ElasticTableStaleError(ElasticTableError):
    """The npz bytes, the id, the schema or the builder do not match the JSON record."""


class ElasticTablePinError(ElasticTableError):
    """A recorded source hash differs from the registry pin."""


class ElasticTableUnqualifiedError(ElasticTableError):
    """A qualification flag of the npz is false (no negative density in the domain, O-16 MT2 copy
    asserted, O-16 MT5 not a copy, shape normalised); the builder writes such tables as evidence,
    the loader refuses them."""


def flags_from_array(q: NDArray[Any]) -> dict[str, bool]:
    if q.shape != (len(QUALIFICATION_FIELDS),):
        raise ElasticTableStaleError("elastic table npz: malformed qualification array (stale)")
    return {k: bool(q[i] == 1) for i, k in enumerate(QUALIFICATION_FIELDS)}


@dataclass(frozen=True)
class MaterialElastic:
    """Per-material rows on the table grid (read-only), as ``MaterialNuclear``:
    ``sigma_mass_cm2_g``,
    ``sigma_hat_window``, ``sigma_hat_end``, ``target_index`` (into the table targets, ``H-1`` = 0),
    ``cum_fraction`` and ``cum_sigma_mass_cm2_g`` (K, N)."""

    material: str
    grid_e_mev: NDArray[np.float64]
    sigma_mass_cm2_g: NDArray[np.float64]
    sigma_hat_window: NDArray[np.float64]
    sigma_hat_end: NDArray[np.float64]
    target_index: tuple[int, ...]
    cum_fraction: NDArray[np.float64]
    cum_sigma_mass_cm2_g: NDArray[np.float64]
    f_e: float

    def sigma_at(self, e_mev: float) -> float:
        return interp_row(self.grid_e_mev, self.sigma_mass_cm2_g, e_mev)


class ElasticTable:
    """A loaded, verified, frozen elastic table."""

    def __init__(self, info: dict[str, Any], arrays: dict[str, NDArray[Any]], npz_path: Path):
        self.info = freeze_json(info)
        self.arrays = arrays
        self.npz_path = npz_path
        self.table_id: str = info["table_id"]
        self.target_names: tuple[str, ...] = tuple(info["target_names"])

    @classmethod
    def load(cls, cache_dir: str | Path | None, table_id: str) -> ElasticTable:
        cdir = cache.resolve_cache_dir(cache_dir)
        npz_path = cdir / "derived" / f"elastic-proton-{table_id}.npz"
        json_path = npz_path.with_suffix(".json")
        if not npz_path.is_file() or not json_path.is_file():
            raise ElasticTableMissingError(f"elastic table {table_id} is not in {cdir / 'derived'}")
        raw = npz_path.read_bytes()
        info = json.loads(json_path.read_text(encoding="utf-8"))
        if hashlib.sha256(raw).hexdigest() != info.get("npz_sha256"):
            raise ElasticTableStaleError(f"elastic table {table_id}: npz bytes changed (stale)")
        for sid, rec in info["sources"].items():
            pin = DATASETS.get(sid)
            if pin is None or pin.sha256 != rec["sha256"]:
                raise ElasticTablePinError(f"elastic table {table_id}: source {sid} pin mismatch")
        if info.get("schema") != SCHEMA or info.get("builder_version") != BUILDER_VERSION:
            raise ElasticTableStaleError(
                f"elastic table {table_id}: schema/builder {info.get('schema')!r}/"
                f"{info.get('builder_version')!r}, expected {SCHEMA!r}/{BUILDER_VERSION!r}"
            )
        if info.get("table_id") != table_id or compute_table_id(info) != table_id:
            raise ElasticTableStaleError(f"elastic table {table_id}: id not reproduced (stale)")
        with np.load(io.BytesIO(raw), allow_pickle=False) as z:
            arrays = {k: freeze_array(z[k], z[k].dtype, k) for k in z.files}
        if "qualification" not in arrays:
            raise ElasticTableStaleError("elastic table npz has no qualification array (stale)")
        if "target_e_min_mev" not in arrays:
            raise ElasticTableStaleError("elastic table npz has no target_e_min_mev (stale)")
        grid, e_min = arrays["grid_e_mev"], arrays["target_e_min_mev"]
        for i in range(len(e_min)):
            below = grid < e_min[i]
            if np.any(arrays["sigma_barn"][i][below] != 0.0) or np.any(arrays["valid"][i][below]):
                raise ElasticTableUnqualifiedError(
                    f"elastic table {table_id}: target {i} has sigma != 0 or valid != 0 below "
                    "its domain (stale)"
                )
        flags = flags_from_array(arrays["qualification"])
        failed = [k for k, v in flags.items() if not v]
        if failed:
            raise ElasticTableUnqualifiedError(
                f"elastic table {table_id} is not qualified for transport: {failed}"
            )
        return cls(info, arrays, npz_path)

    # -- domain --------------------------------------------------------------------------------
    def check_energy(self, e_mev: float) -> None:
        """``UnsupportedCombinationError`` outside [1, 250] MeV (the table domain)."""
        if not E_MIN_MEV <= e_mev <= E_MAX_MEV:
            raise UnsupportedCombinationError(
                f"elastic scattering: energy {e_mev} MeV outside the table domain "
                f"[{E_MIN_MEV}, {E_MAX_MEV}] MeV"
            )

    def elastic_domain(self) -> dict[str, tuple[float, float]]:
        """``{target: (e_min, e_max)}`` [MeV] over which each table target is defined (H-1:
        ``E_min,pp``; p + A: ``e_min_shape``). Below ``e_min`` the elastic channel is omitted
        (declared limitation; ``sigma = 0`` and ``valid = 0``, asserted at load)."""
        e_min = self.arrays["target_e_min_mev"]
        return {n: (float(e_min[i]), E_MAX_MEV) for i, n in enumerate(self.target_names)}

    def sigma_barn(self, target: int, e_mev: float) -> float:
        """sigma_el of table target ``target`` at ``e_mev`` [b], lin-lin in E."""
        self.check_energy(e_mev)
        return interp_row(self.arrays["grid_e_mev"], self.arrays["sigma_barn"][target], e_mev)

    def edges_at(self, target: int, e_mev: float) -> NDArray[np.float64]:
        """The ``n_q + 1`` mu_CM edges at ``e_mev`` (quantile interpolation lin in E)."""
        self.check_energy(e_mev)
        grid = self.arrays["grid_e_mev"]
        k = locate(grid, e_mev)
        t = min(max((e_mev - grid[k]) / (grid[k + 1] - grid[k]), 0.0), 1.0)
        e = self.arrays["edges_mu"][target]
        return np.asarray((1.0 - t) * e[k] + t * e[k + 1])

    def sample_mu(self, target: int, e_mev: float, u: NDArray[np.float64]) -> NDArray[np.float64]:
        """Reference inverse-CDF sample of mu_CM for uniforms ``u`` in [0, 1): bin
        ``i = floor(n_q u)``, linear in the bin (the table's piecewise-uniform density)."""
        edges = self.edges_at(target, e_mev)
        n_q = edges.size - 1
        x = np.asarray(u, dtype=np.float64) * n_q
        i = np.minimum(x.astype(np.int64), n_q - 1)
        return np.asarray(edges[i] + (x - i) * (edges[i + 1] - edges[i]))

    # -- materials -----------------------------------------------------------------------------
    def material_rows(self, material: Material, f_e: float = DEFAULT_F_E) -> MaterialElastic:
        """Rows of ``material``; an element absent from the table raises
        ``UnsupportedCombinationError``."""
        elements = self.info["elements"]
        grid = self.arrays["grid_e_mev"]
        sig = self.arrays["sigma_barn"]
        per_target: dict[int, NDArray[np.float64]] = {}
        for sym, w in material.mass_fractions.items():
            row = elements.get(sym)
            if row is None:
                raise UnsupportedCombinationError(
                    f"elastic=True: element {sym} of material {material.name!r} has no elastic "
                    "data (evaluated, BGG or surrogate) in the elastic table"
                )
            t = self.target_names.index(row["target"])
            term = N_A * 1.0e-24 * w * row["sigma_scale"] * sig[t] / row["a_g_mol"]
            per_target[t] = per_target.get(t, 0.0) + term
        order = sorted(per_target)
        total = sum(per_target.values())
        cum_sigma = np.cumsum([per_target[t] for t in order], axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            cum = np.where(total > 0.0, cum_sigma / np.where(total > 0.0, total, 1.0), 0.0)
        cum[-1] = 1.0
        win, end = majorant(grid, np.asarray(total), f_e)
        return MaterialElastic(
            material.name,
            grid,
            freeze_array(np.asarray(total), np.float64),
            freeze_array(win, np.float64),
            freeze_array(end, np.float64),
            tuple(order),
            freeze_array(cum, np.float64),
            freeze_array(cum_sigma, np.float64),
            f_e,
        )


def chi2_equiprobable(
    mu: NDArray[np.float64], edges: NDArray[np.float64], n_bins: int = 64
) -> tuple[float, float, int]:
    """P7 harness: chi-square of the sampled ``mu`` histogram in ``n_bins`` equiprobable bins of
    the table CDF (the table edges grouped ``n_q / n_bins`` at a time) against the uniform
    expectation. Returns ``(chi2, p_value, dof)``; the p-value uses the regularised upper incomplete
    gamma function by series / continued fraction (no scipy)."""
    n_q = edges.size - 1
    if n_q % n_bins:
        raise ValueError("n_q must be a multiple of n_bins")
    sel = edges[:: n_q // n_bins]
    counts = np.histogram(np.asarray(mu), bins=sel)[0].astype(float)
    exp = counts.sum() / n_bins
    chi2 = float(np.sum((counts - exp) ** 2 / exp))
    dof = n_bins - 1
    return chi2, gamma_q(dof / 2.0, chi2 / 2.0), dof


def gamma_q(a: float, x: float) -> float:
    """Regularised upper incomplete gamma function Q(a, x) (Numerical Recipes series/CF)."""
    import math

    if x <= 0.0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(1000):
            ap += 1.0
            d *= x / ap
            s += d
            if abs(d) < abs(s) * 1e-16:
                break
        return 1.0 - s * math.exp(-x + a * math.log(x) - gln)
    b = x + 1.0 - a
    c, dd = 1e300, 1.0 / b
    h = dd
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        dd = an * dd + b
        dd = 1e-300 if abs(dd) < 1e-300 else dd
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        dd = 1.0 / dd
        delta = dd * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return math.exp(-x + a * math.log(x) - gln) * h
