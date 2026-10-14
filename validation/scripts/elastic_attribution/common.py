"""Shared loaders and report helpers of the D-E steps (exploratory, report-only; Amendment 19).

The built elastic table is loaded through the repository loader ``ionmc.nuclear.elastic_tables.
ElasticTable.load`` (the pinned, qualified table; ``allow_exploratory_table`` is not needed). The
archived inputs are read from ``validation/generated`` and ``.ionmc-cache``.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

REPO = Path(__file__).resolve().parents[3]
CACHE_DIR = Path(os.environ.get("IONMC_CACHE_DIR") or REPO / ".ionmc-cache" / "ionmc-data")
os.environ.setdefault(
    "IONMC_CACHE_DIR", str(CACHE_DIR)
)  # the transport helpers resolve the cache from it
TABLE_ID = "1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591"
GEN = REPO / "validation" / "generated"
RUN_V10A = GEN / "RUN-20261012T142343Z-081123c3" / "transport" / "lv5c-dd976c5-03" / "03-v10-a.txt"
RUN_V11 = (
    GEN / "RUN-20261012T142957Z-a8e41216" / "transport" / "lv5c-dd976c5-16" / "16-v11-compare.txt"
)
SAMPLES = REPO / ".ionmc-cache" / "transport-samples" / "dd976c5"
V11_PARTIALS = {
    (150, "emel"): SAMPLES / "12" / "samples" / "v11-150-emel.json",
    (150, "emonly"): SAMPLES / "13" / "samples" / "v11-150-emonly.json",
    (200, "emel"): SAMPLES / "14" / "samples" / "v11-200-emel.json",
    (200, "emonly"): SAMPLES / "15" / "samples" / "v11-200-emonly.json",
}
OUT = GEN / "elastic_attribution"
BANNER = "Exploratory, report-only (Amendment 19); not evidence; changes no row."
N_AVOGADRO = 6.02214076e23


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO))
    except ValueError:
        return str(path)


def load_table() -> Any:
    """The pinned elastic table through the repository loader (verified pins, no exploratory flag)."""
    from ionmc.nuclear.elastic_tables import ElasticTable

    return ElasticTable.load(CACHE_DIR, TABLE_ID)


def table_files() -> list[Path]:
    d = CACHE_DIR / "derived"
    return [d / f"elastic-proton-{TABLE_ID}.npz", d / f"elastic-proton-{TABLE_ID}.json"]


def json_payload(path: Path) -> dict[str, Any]:
    """The JSON document of a step ``.txt`` (between ``#JSON-BEGIN`` and ``#JSON-END`` if present)."""
    t = Path(path).read_text()
    if "#JSON-BEGIN" in t:
        t = t.split("#JSON-BEGIN", 1)[1].split("#JSON-END")[0]
    return json.loads(t)  # type: ignore[no-any-return]


def v11_reference() -> dict[str, Any]:
    """Archived V11 numbers (step 16): ionmc and TOPAS F and plateau gain per energy."""
    v = json_payload(RUN_V11)["verdict"]["energies"]
    return {
        int(e): {
            "F_ionmc": d["F"]["ionmc"]["value"],
            "F_ionmc_se": d["F"]["ionmc"]["se"],
            "F_topas": d["F"]["reference"]["value"],
            "F_topas_se": d["F"]["reference"]["se"],
            "gain_ionmc": d["plateau_gain"]["ionmc"]["value"],
            "gain_topas": d["plateau_gain"]["reference"]["value"],
        }  # fmt: skip
        for e, d in v.items()
    }


def idd_mean(
    energy: int, mode: str
) -> tuple[NDArray[np.float64], NDArray[np.float64], dict[str, Any]]:
    """``(depth bin centres [mm], mean IDD over batches [MeV/(g/cm^2)/primary], metrics)`` of a V11 partial."""
    j = json.loads(V11_PARTIALS[(energy, mode)].read_text())
    idd = np.asarray(j["idd_batches"], dtype=float).mean(axis=0)
    z = (np.arange(idd.size) + 0.5) * float(j["bin_mm"])
    return z, idd, j["metrics"]


class Water:
    """Liquid water of the V5 geometry: Bethe CSDA range table (the source the V5/V11 runs use) and the
    atom densities of the elastic table targets (``H-1``, ``O-16``) with their element ``sigma_scale``."""

    def __init__(self, table: Any) -> None:
        from ionmc.materials import WATER
        from ionmc.physics.projectiles import PROTON
        from ionmc.physics.stopping import BetheStoppingSource

        self.stop = BetheStoppingSource().table(WATER, PROTON)
        self.rho = float(WATER.density_g_cm3)
        self.targets: dict[str, dict[str, float]] = {}
        for sym, w in WATER.mass_fractions.items():
            el = table.info["elements"][sym]
            self.targets[str(el["target"])] = {
                "n_per_cm3": N_AVOGADRO * self.rho * float(w) / float(el["a_g_mol"]),
                "scale": float(el["sigma_scale"]),
            }

    def range_mm(self, e_mev: Any) -> Any:
        return 10.0 * np.asarray(self.stop.range_at(np.asarray(e_mev, dtype=float))) / self.rho

    def energy_at_depth(self, e0: float, z_mm: NDArray[np.float64]) -> NDArray[np.float64]:
        """CSDA energy [MeV] at depth ``z_mm`` of a proton entering at ``e0``."""
        r = float(self.range_mm(e0)) - np.asarray(z_mm, dtype=float)
        return np.asarray(self.stop.energy_from_range(r * self.rho / 10.0))

    def stopping_mev_per_cm(self, e_mev: Any) -> Any:
        return np.asarray(self.stop.stopping_at(np.asarray(e_mev, dtype=float))) * self.rho


def header(title: str, command: str, inputs: list[Path], extra: list[str] | None = None) -> str:
    lines = [BANNER, "", f"# {title}", "", f"Command (from the worktree root, `.venv` python, single thread): `{command}`", "",
             "Input files (sha256):", ""]  # fmt: skip
    lines += [f"- `{rel(p)}`  {sha256_file(p)}" for p in inputs]
    lines += ["", *(extra or [])]
    return "\n".join(lines) + "\n"


# Git projection (decision 0041: Git holds no raw or derived nuclear data, only aggregate statistics).
# Keys that carry experimental or compilation values, or per-point arrays from which they are losslessly
# recoverable (model value / (1 + ratio)), are removed from the copy committed under validation/results/.
GIT_FORBIDDEN_KEYS: frozenset[str] = frozenset(
    {"data_mb_sr", "rel_err", "model_total_mb_sr", "model_ratio_minus_1", "rutherford_mb_sr",
     "family_total_mb_sr", "pdg_lowest_points"}
)  # fmt: skip
GIT_MAX_NUMERIC_LIST = 8  # longer numeric lists are per-point arrays, never aggregate statistics
GIT_NOTE = (
    "Git projection (decision 0041): experimental / compilation values and every per-point array are removed; "
    "only aggregate statistics (window integrals, ratios, minima diagnostics, fits) remain. The complete record is "
    "regenerated into validation/generated/ by the committed module."
)


def _is_numeric_list(x: Any) -> bool:
    return isinstance(x, list) and len(x) > 0 and all(v is None or isinstance(v, (int, float)) for v in x)


def git_projection(doc: Any) -> Any:
    """Return a deep copy of ``doc`` without the forbidden keys and without per-point numeric arrays."""
    if isinstance(doc, dict):
        return {k: git_projection(v) for k, v in doc.items()
                if k not in GIT_FORBIDDEN_KEYS and not (_is_numeric_list(v) and len(v) > GIT_MAX_NUMERIC_LIST)}  # fmt: skip
    if isinstance(doc, list):
        return [git_projection(v) for v in doc]
    return doc


def git_safety_violations(doc: Any, path: str = "") -> list[str]:
    """Paths inside ``doc`` that a committed record must not contain (empty when the record is Git-safe)."""
    bad: list[str] = []
    if isinstance(doc, dict):
        for k, v in doc.items():
            here = f"{path}/{k}"
            if k in GIT_FORBIDDEN_KEYS:
                bad.append(here)
            elif _is_numeric_list(v) and len(v) > GIT_MAX_NUMERIC_LIST:
                bad.append(f"{here}[len={len(v)}]")
            else:
                bad += git_safety_violations(v, here)
    elif isinstance(doc, list):
        for i, v in enumerate(doc):
            bad += git_safety_violations(v, f"{path}[{i}]")
    return bad


def write_outputs(outdir: Path, json_name: str, doc: dict[str, Any], md_name: str, md: str) -> None:
    """Write the complete record to ``outdir`` and its Git projection to ``outdir / "git"``.

    Only the ``git`` copy may be committed under ``validation/results/elastic_attribution/``.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / json_name).write_text(json.dumps(doc, indent=1, sort_keys=True, default=float) + "\n")
    (outdir / md_name).write_text(md)
    proj = git_projection(json.loads(json.dumps(doc, default=float)))
    proj.setdefault("note_git", GIT_NOTE)
    assert not git_safety_violations(proj), git_safety_violations(proj)
    (outdir / "git").mkdir(exist_ok=True)
    (outdir / "git" / json_name).write_text(json.dumps(proj, indent=1, sort_keys=True) + "\n")


def import_transport_module(name: str) -> Any:
    """Import a module of ``validation/scripts/transport`` (read only; no bytecode is written)."""
    sys.dont_write_bytecode = True
    p = str(REPO / "validation" / "scripts" / "transport")
    if p not in sys.path:
        sys.path.insert(0, p)
    import importlib

    return importlib.import_module(name)
