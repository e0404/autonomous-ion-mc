"""Loader of the derived nuclear table and its per-material rows (decision 0041 sections 1, 6).

``NuclearTable.load(cache_dir, table_id)`` reads ``<cache>/derived/nuclear-proton-<id>.npz`` and
``.json`` with ``np.load(allow_pickle=False)``, re-hashes the npz bytes against the JSON, checks
that the recorded table id is reproduced by the recorded sources/builder/options and that every
source pin equals the registry pin, and freezes the arrays. Failures are fail-closed with
specific exception types (all subclass :class:`NuclearTableError`).

``material_rows`` composes ``Sigma_mass(E) = N_A sum_el w_el sigma_el / A_el`` [cm2/g] on the
table's ln E grid, the cumulative target fractions per node and the majorants (:func:`majorant`).
"""

from __future__ import annotations

import hashlib
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
from ionmc.nuclear.build import BuildOptions
from ionmc.nuclear.build import table_id as compute_table_id

MAJORANT_FACTOR = 1.02
DEFAULT_F_E = 0.02


class NuclearTableError(RuntimeError):
    """Base class of nuclear-table load failures."""


class NuclearTableMissingError(NuclearTableError, FileNotFoundError):
    """The table files are not in the cache."""


class NuclearTableStaleError(NuclearTableError):
    """The npz bytes or the table id do not match the JSON record."""


class NuclearTablePinError(NuclearTableError):
    """A recorded source hash differs from the registry pin."""


def majorant(
    grid_e: NDArray[np.float64],
    sigma: NDArray[np.float64],
    f_e: float = DEFAULT_F_E,
    factor: float = MAJORANT_FACTOR,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """``(window, end_of_range)`` majorant arrays on the grid nodes.

    ``sigma`` is interpolated linearly in ln E between nodes (the runtime rule). The window array
    at node ``k`` is ``factor`` times the maximum of that interpolant over
    ``[E_k (1 - 2 f_E - 0.01), E_{k+1}]`` (the frozen window ``[E(1-2f_E-0.01), E]`` extended to
    the next node, so that the STEP lookup ``win[k]`` (node at or below ``E0``, no interpolation)
    is a majorant for every ``E0`` in ``[E_k, E_{k+1}]``; interpolating between nodes would not
    be); the end-of-range array is ``factor`` times the running maximum
    over ``[0, E_{k+1}]``."""
    ln_e = np.log(grid_e)
    n = grid_e.size
    window = np.empty(n)
    for k in range(n):
        lo = grid_e[k] * (1.0 - 2.0 * f_e - 0.01)
        hi_i = min(k + 1, n - 1)
        inside = (grid_e >= lo) & (np.arange(n) <= hi_i)
        cand = [float(sigma[inside].max(initial=0.0))]
        if lo > grid_e[0]:
            cand.append(float(np.interp(np.log(lo), ln_e, sigma)))
        window[k] = factor * max(cand)
    run = np.maximum.accumulate(sigma)
    end = factor * np.concatenate((run[1:], run[-1:]))
    return window, end


@dataclass(frozen=True)
class MaterialNuclear:
    """Per-material rows on the table grid (read-only): ``grid_e_mev``; ``sigma_mass_cm2_g``;
    ``sigma_hat_window`` and ``sigma_hat_end`` majorants [cm2/g]; ``target_index`` (K,) indices
    into the table targets and ``cum_fraction`` (K, N) cumulative probability of each target per
    node (the last row is 1; all 0 except the last where Sigma = 0); ``f_e`` of the window."""

    material: str
    grid_e_mev: NDArray[np.float64]
    sigma_mass_cm2_g: NDArray[np.float64]
    sigma_hat_window: NDArray[np.float64]
    sigma_hat_end: NDArray[np.float64]
    target_index: tuple[int, ...]
    cum_fraction: NDArray[np.float64]
    f_e: float


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
        if hashlib.sha256(raw).hexdigest() != info.get("npz_sha256"):
            raise NuclearTableStaleError(f"nuclear table {table_id}: npz bytes changed (stale)")
        for sid, rec in info["sources"].items():
            pin = DATASETS.get(sid)
            if pin is None or pin.sha256 != rec["sha256"]:
                raise NuclearTablePinError(f"nuclear table {table_id}: source {sid} pin mismatch")
        opts = info["options"]
        opts["lambda_nodes_mev"] = (
            None if opts.get("lambda_nodes_mev") is None else tuple(opts["lambda_nodes_mev"])
        )
        expected = compute_table_id(
            {sid: rec["sha256"] for sid, rec in info["sources"].items()}, BuildOptions(**opts)
        )
        if info["table_id"] != table_id or expected != table_id:
            raise NuclearTableStaleError(f"nuclear table {table_id}: id not reproduced (stale)")
        with np.load(npz_path, allow_pickle=False) as z:
            arrays = {k: freeze_array(z[k], z[k].dtype, k) for k in z.files}
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
            cum = np.cumsum([per_target[t] for t in order], axis=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                cum = np.where(total > 0.0, cum / np.where(total > 0.0, total, 1.0), 0.0)
            cum[-1] = 1.0
        else:
            total, order, cum = np.zeros_like(grid), [], np.zeros((0, grid.size))
        win, end = majorant(grid, np.asarray(total), f_e)
        return MaterialNuclear(
            material.name,
            grid,
            freeze_array(total, np.float64),
            freeze_array(win, np.float64),
            freeze_array(end, np.float64),
            tuple(order),
            freeze_array(cum, np.float64),
            f_e,
        )

    def product_rows(self, target: int) -> dict[str, NDArray[np.float64]]:
        """Per-target product rows (read-only views): ``lam`` (5, N), ``edges_mev`` (5, N, 65),
        ``r_pre`` (5, N, 64), ``yield_endf`` (5, N), ``mean_ecm_mev`` (5, N)."""
        keys = ("lam", "edges_mev", "r_pre", "yield_endf", "mean_ecm_mev")
        return {k: self.arrays[k][target] for k in keys}
