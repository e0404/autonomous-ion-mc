"""Loader of the derived nuclear table and its per-material rows (decision 0041 sections 1, 6).

``NuclearTable.load(cache_dir, table_id)`` reads ``<cache>/derived/nuclear-proton-<id>.npz`` and
``.json`` with ``np.load(allow_pickle=False)``, re-hashes the npz bytes against the JSON, checks
that the recorded table id is reproduced by the recorded sources/builder/options and that every
source pin equals the registry pin, and freezes the arrays. Failures are fail-closed with
specific exception types (all subclass :class:`NuclearTableError`).

``material_rows`` composes ``Sigma_mass(E) = N_A sum_el w_el sigma_el / A_el`` [cm2/g] on the
table's union grid, the cumulative partial ``Sigma`` per target and the majorants
(:func:`majorant`). Every runtime row is interpolated lin-lin in E between the union nodes; the
node is found by ``grid_locate`` (fixed-step bisection, shared function with a bitwise python
twin; :func:`locate` is the python entry point).
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
from ionmc._wpfunc import python_twin
from ionmc.data import cache
from ionmc.data.registry import DATASETS
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import N_A, Material
from ionmc.nuclear.build import (
    BUILDER_VERSION,
    QUALIFICATION_FIELDS,
    SCHEMA,
    BuildOptions,
    bounds_from_array,
)
from ionmc.nuclear.build import table_id as compute_table_id
from ionmc.physics.nuclear import make_nuclear

MAJORANT_FACTOR = 1.02
DEFAULT_F_E = 0.02


class NuclearTableError(UnsupportedCombinationError, RuntimeError):
    """Base class of nuclear-table load failures (also an ``UnsupportedCombinationError``, so
    that ``validate()`` lets them propagate as fail-closed configuration errors)."""


class NuclearTableMissingError(NuclearTableError, FileNotFoundError):
    """The table files are not in the cache."""


class NuclearTableStaleError(NuclearTableError):
    """The npz bytes or the table id do not match the JSON record."""


class NuclearTableUnqualifiedError(NuclearTableError):
    """The table was built but is not qualified for transport: the D6 ceiling failed (count- or
    energy-weighted range) or a multiplicity node did not converge. The builder writes such tables
    as evidence; the loader refuses them (decision 0041 section 5)."""


class NuclearTablePinError(NuclearTableError):
    """A recorded source hash differs from the registry pin."""


def qualification_failures(info: dict[str, Any]) -> list[str]:
    """Why the table JSON ``info`` is not qualified for transport (empty list: qualified): the D6
    ceiling must pass in both weightings (``gate_d6.ceiling_pass`` and
    ``ceiling_pass_energy_weighted_range``) and every lambda node must have converged. A missing
    field counts as a failure (fail closed)."""
    out: list[str] = []
    gate = info.get("gate_d6", {})
    if gate.get("ceiling_pass") is not True:
        out.append("gate_d6.ceiling_pass is not true (count-weighted range ceiling)")
    if gate.get("ceiling_pass_energy_weighted_range") is not True:
        out.append("gate_d6.ceiling_pass_energy_weighted_range is not true")
    mult = info.get("multiplicity", {})
    nodes = mult.get("non_converged_nodes")
    if nodes is None or len(nodes) > 0 or mult.get("all_nodes_converged", True) is not True:
        out.append("multiplicity has non-converged lambda nodes")
    return out


def derive_gates(arrays: dict[str, NDArray[Any]]) -> dict[str, Any]:
    """The gating quantities read from the npz arrays alone: ``qualification`` (flags of
    :data:`QUALIFICATION_FIELDS`, with ``all_nodes_converged`` also recomputed from
    ``lam_converged``), the bounds of ``bounds`` and the number of non-converged lambda nodes. A
    missing array raises :class:`NuclearTableStaleError`."""
    for name in ("qualification", "bounds", "lam_converged"):
        if name not in arrays:
            raise NuclearTableStaleError(f"nuclear table npz has no {name!r} array (stale)")
    q = arrays["qualification"]
    if q.shape != (len(QUALIFICATION_FIELDS),):
        raise NuclearTableStaleError("nuclear table npz: malformed qualification array (stale)")
    flags = {k: bool(q[i] == 1) for i, k in enumerate(QUALIFICATION_FIELDS)}
    n_bad = int(np.sum(arrays["lam_converged"] != 1))
    out: dict[str, Any] = {"flags": flags, "non_converged_count": n_bad}
    try:
        out.update(bounds_from_array(arrays["bounds"]))
    except ValueError as exc:
        raise NuclearTableStaleError(f"nuclear table npz: {exc} (stale)") from exc
    return out


def qualification_failures_npz(derived: dict[str, Any]) -> list[str]:
    """Why the npz-derived gates (:func:`derive_gates`) do not qualify the table (empty list:
    qualified): both D6 ceilings must pass and every lambda node must have converged."""
    f = derived["flags"]
    out: list[str] = []
    if not f["ceiling_pass"]:
        out.append("ceiling_pass is not true (count-weighted range ceiling)")
    if not f["ceiling_pass_energy_weighted_range"]:
        out.append("ceiling_pass_energy_weighted_range is not true")
    if derived["non_converged_count"] > 0 or not f["all_nodes_converged"]:
        out.append("multiplicity has non-converged lambda nodes")
    return out


def sidecar_disagreements(info: dict[str, Any], derived: dict[str, Any]) -> list[str]:
    """Fields of the JSON sidecar ``info`` that differ from the npz-derived ``derived`` gates."""
    out: list[str] = []
    gate = info.get("gate_d6", {})
    mult = info.get("multiplicity", {})
    f = derived["flags"]
    json_flags = {
        "ceiling_pass": gate.get("ceiling_pass"),
        "ceiling_pass_energy_weighted_range": gate.get("ceiling_pass_energy_weighted_range"),
        "all_nodes_converged": mult.get("all_nodes_converged"),
        "tier1_pass": gate.get("tier1_pass"),
        "tier2_pass": gate.get("tier2_pass"),
    }
    for k, v in json_flags.items():
        if v is not f[k]:
            out.append(f"{k}: JSON {v!r}, npz {f[k]!r}")
    n_json = mult.get("non_converged_nodes")
    if not isinstance(n_json, list) or (len(n_json) > 0) != (derived["non_converged_count"] > 0):
        out.append("multiplicity.non_converged_nodes")
    for k in (
        "history_energy_bound_mev",
        "transport_energy_bound_mev",
        "recoil_t_max_mev",
        "transport_path_bound_terms",
    ):
        if info.get(k) != derived[k]:
            out.append(f"{k}: JSON {info.get(k)!r}, npz {derived[k]!r}")
    return out


def locate(grid: NDArray[np.float64], e_mev: float) -> int:
    """Interval ``k`` with ``grid[k] <= e < grid[k+1]`` (clamped to ``[0, n-2]``) by the shared
    ``grid_locate`` (python twin): the lookup of every runtime row."""
    nu = python_twin(make_nuclear)
    return int(nu.grid_locate(float(e_mev), grid, int(grid.size)))


def interp_row(grid: NDArray[np.float64], row: NDArray[np.float64], e_mev: float) -> float:
    """``row`` (values on the grid nodes) at ``e_mev``, lin-lin in E, constant outside the grid."""
    k = locate(grid, e_mev)
    t = min(max((e_mev - grid[k]) / (grid[k + 1] - grid[k]), 0.0), 1.0)
    return float((1.0 - t) * row[k] + t * row[k + 1])


def majorant(
    grid_e: NDArray[np.float64],
    sigma: NDArray[np.float64],
    f_e: float = DEFAULT_F_E,
    factor: float = MAJORANT_FACTOR,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """``(window, end_of_range)`` majorant arrays on the grid nodes.

    ``sigma`` is interpolated lin-lin in E between the (non-uniform) nodes (the runtime rule); the
    maximum of that piecewise-linear interpolant over an interval is attained at a node or at the
    interval ends. The window array at node ``k`` is ``factor`` times the maximum over
    ``[E_k (1 - 2 f_E - 0.01), E_{k+1}]`` (the frozen window ``[E(1-2f_E-0.01), E]`` extended to
    the next node, so that the STEP lookup ``win[k]`` (node at or below ``E0``, no interpolation)
    is a majorant for every ``E0`` in ``[E_k, E_{k+1}]``; interpolating between nodes would not
    be); the end-of-range array is ``factor`` times the running maximum
    over ``[0, E_{k+1}]``."""
    n = grid_e.size
    window = np.empty(n)
    for k in range(n):
        lo = grid_e[k] * (1.0 - 2.0 * f_e - 0.01)
        hi_i = min(k + 1, n - 1)
        inside = (grid_e >= lo) & (np.arange(n) <= hi_i)
        cand = [float(sigma[inside].max(initial=0.0))]
        if lo > grid_e[0]:
            cand.append(float(np.interp(lo, grid_e, sigma)))
        window[k] = factor * max(cand)
    run = np.maximum.accumulate(sigma)
    end = factor * np.concatenate((run[1:], run[-1:]))
    return window, end


@dataclass(frozen=True)
class MaterialNuclear:
    """Per-material rows on the table grid (read-only): ``grid_e_mev``; ``sigma_mass_cm2_g``;
    ``sigma_hat_window`` and ``sigma_hat_end`` majorants [cm2/g]; ``target_index`` (K,) indices
    into the table targets and ``cum_fraction`` (K, N) cumulative probability of each target per
    node (the last row is 1; all 0 except the last where Sigma = 0); ``cum_sigma_mass_cm2_g``
    (K, N) the unnormalised cumulative partial Sigma per node (runtime: ``cum_k(E) / Sigma(E)``
    from lin-lin interpolation of both, so ``select_target`` is unchanged); ``f_e`` of the
    window."""

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
        """``Sigma_mass(E)`` [cm2/g], lin-lin in E on the union grid."""
        return interp_row(self.grid_e_mev, self.sigma_mass_cm2_g, e_mev)

    def cum_fraction_at(self, e_mev: float) -> NDArray[np.float64]:
        """Cumulative target fractions at ``E``: the lin-lin cumulative partial Sigma over the
        lin-lin total (zeros where ``Sigma = 0``)."""
        tot = self.sigma_at(e_mev)
        cum = np.array([interp_row(self.grid_e_mev, c, e_mev) for c in self.cum_sigma_mass_cm2_g])
        return cum / tot if tot > 0.0 else np.zeros_like(cum)


class NuclearTable:
    """A loaded, verified, frozen nuclear table."""

    def __init__(self, info: dict[str, Any], arrays: dict[str, NDArray[Any]], npz_path: Path):
        self.info = freeze_json(info)
        self.arrays = arrays
        self.npz_path = npz_path
        self.table_id: str = info["table_id"]
        self.target_names = tuple(t["name"] for t in info["targets"])

    @classmethod
    def load(cls, cache_dir: str | Path | None, table_id: str) -> NuclearTable:
        """Load and verify table ``table_id`` (see the module docstring)."""
        cdir = cache.resolve_cache_dir(cache_dir)
        npz_path = cdir / "derived" / f"nuclear-proton-{table_id}.npz"
        json_path = npz_path.with_suffix(".json")
        if not npz_path.is_file() or not json_path.is_file():
            raise NuclearTableMissingError(f"nuclear table {table_id} is not in {cdir / 'derived'}")
        raw = npz_path.read_bytes()
        info = json.loads(json_path.read_text(encoding="utf-8"))
        npz_sha = hashlib.sha256(raw).hexdigest()
        if npz_sha != info.get("npz_sha256"):
            raise NuclearTableStaleError(f"nuclear table {table_id}: npz bytes changed (stale)")
        for sid, rec in info["sources"].items():
            pin = DATASETS.get(sid)
            if pin is None or pin.sha256 != rec["sha256"]:
                raise NuclearTablePinError(f"nuclear table {table_id}: source {sid} pin mismatch")
        if info.get("schema") != SCHEMA or info.get("builder_version") != BUILDER_VERSION:
            raise NuclearTableStaleError(
                f"nuclear table {table_id}: schema/builder {info.get('schema')!r}/"
                f"{info.get('builder_version')!r}, expected {SCHEMA!r}/{BUILDER_VERSION!r}"
            )
        opts = dict(info["options"])
        opts["diagnostic_nodes_mev"] = (
            None
            if opts.get("diagnostic_nodes_mev") is None
            else tuple(opts["diagnostic_nodes_mev"])
        )
        try:
            build_options = BuildOptions(**opts)
        except TypeError as exc:  # options of another builder revision: the id cannot be checked
            raise NuclearTableStaleError(
                f"nuclear table {table_id}: unknown build options ({exc}); stale"
            ) from exc
        # the id covers the npz bytes recomputed here (not the sidecar's digest), so a resealed
        # sidecar cannot reproduce the pinned id of altered arrays
        expected = compute_table_id(
            {sid: rec["sha256"] for sid, rec in info["sources"].items()}, build_options, npz_sha
        )
        if info["table_id"] != table_id or expected != table_id:
            raise NuclearTableStaleError(f"nuclear table {table_id}: id not reproduced (stale)")
        with np.load(io.BytesIO(raw), allow_pickle=False) as z:
            arrays = {k: freeze_array(z[k], z[k].dtype, k) for k in z.files}
        # the gating quantities are those of the authenticated arrays; the JSON is a view of them
        derived = derive_gates(arrays)
        disagree = sidecar_disagreements(info, derived)
        if disagree:
            raise NuclearTableStaleError(
                f"nuclear table {table_id}: the JSON sidecar disagrees with the npz ("
                + "; ".join(disagree)
                + "; stale)"
            )
        reasons = qualification_failures_npz(derived)
        if reasons:
            raise NuclearTableUnqualifiedError(
                f"nuclear table {table_id} is not qualified for transport: " + "; ".join(reasons)
            )
        return cls(info, arrays, npz_path)

    def material_rows(self, material: Material, f_e: float = DEFAULT_F_E) -> MaterialNuclear:
        """Rows of ``material``; an element absent from the table raises
        ``UnsupportedCombinationError`` (hydrogen contributes 0)."""
        elements = self.info["elements"]
        grid = self.arrays["grid_e_mev"]
        sig = self.arrays["sigma_barn"]
        per_target: dict[int, NDArray[np.float64]] = {}
        for sym, w in material.mass_fractions.items():
            if sym == "H":
                continue
            row = elements.get(sym)
            if row is None:
                raise UnsupportedCombinationError(
                    f"nuclear=True: element {sym} of material {material.name!r} has no "
                    "evaluated or surrogate entry in the nuclear table"
                )
            t = self.target_names.index(row["target"])
            term = N_A * 1.0e-24 * w * row["sigma_scale"] * sig[t] / row["a_g_mol"]
            per_target[t] = per_target.get(t, 0.0) + term
        if per_target:
            total = sum(per_target.values())
            order = sorted(per_target)
            cum_sigma = np.cumsum([per_target[t] for t in order], axis=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                cum = np.where(total > 0.0, cum_sigma / np.where(total > 0.0, total, 1.0), 0.0)
            cum[-1] = 1.0
        else:
            total, order = np.zeros_like(grid), []
            cum, cum_sigma = np.zeros((0, grid.size)), np.zeros((0, grid.size))
        win, end = majorant(grid, np.asarray(total), f_e)
        return MaterialNuclear(
            material.name,
            grid,
            freeze_array(total, np.float64),
            freeze_array(win, np.float64),
            freeze_array(end, np.float64),
            tuple(order),
            freeze_array(cum, np.float64),
            freeze_array(cum_sigma, np.float64),
            f_e,
        )

    def product_rows(self, target: int) -> dict[str, NDArray[np.float64]]:
        """Per-target product rows (read-only views): ``lam`` (5, N), ``edges_mev`` (5, N, 65),
        ``r_pre`` (5, N, 64), ``yield_endf`` (5, N), ``mean_ecm_mev`` (5, N), ``recoil_t_cm_mev``
        (N,), ``p_accept`` (N,)."""
        keys = (
            "lam", "edges_mev", "r_pre", "yield_endf", "mean_ecm_mev", "recoil_t_cm_mev",
            "p_accept",
        )  # fmt: skip
        return {k: self.arrays[k][target] for k in keys}
