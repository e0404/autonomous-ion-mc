"""Shared loaders and report helpers of the D-E steps (exploratory, report-only; Amendment 19).

The built elastic table is loaded through the repository loader ``ionmc.nuclear.elastic_tables.
ElasticTable.load`` (the pinned, qualified table; ``allow_exploratory_table`` is not needed). The
archived inputs are read from ``validation/generated`` and ``.ionmc-cache``.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from collections.abc import Callable, Iterator
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


def v11_archived_block(path: Path | None = None) -> dict[str, dict[str, float]]:
    """Archived V11 numbers of the step-16 record per energy (the block recorded in the committed H5 record as ``archived_v11``): F and its
    SE for ionmc and TOPAS, the plateau gains and the SE of the relative F difference."""
    p = Path(path or RUN_V11)
    v = json_payload(p)["verdict"]["energies"]
    return {
        str(int(e)): {"F_ionmc": d["F"]["ionmc"]["value"], "F_ionmc_se": d["F"]["ionmc"]["se"], "F_topas": d["F"]["reference"]["value"],
                      "F_topas_se": d["F"]["reference"]["se"], "gain_ionmc": d["plateau_gain"]["ionmc"]["value"],
                      "gain_topas": d["plateau_gain"]["reference"]["value"], "F_relative_diff_se": d["F"]["se"]}
        for e, d in v.items()
    }  # fmt: skip


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
# Policy for committed records under validation/results/elastic_attribution/: ONLY dimensionless ratios, per-primary
# probabilities and rates, kinematic / geometric statistics (angles, medians, exceedance probabilities, tail moments, mean
# energies of secondaries), kernel integrals, counts and flags. NO absolute cross-section values (mb, mb/sr: table, formula or
# data), NO sampled or tabulated distributions, histograms or CDFs, NO per-point arrays.
#
# The projection is an explicit ALLOWLIST: every module declares ``GIT_ALLOW``, a tuple of path patterns, and only the leaves
# whose path matches a pattern are kept. A path joins dict keys and list indices with "/". Pattern syntax: a segment is
# literal, or ``*`` (exactly one segment: any dict key or list index), and a final ``**`` matches one or more remaining segments
# (keys may contain "/", e.g. file paths). ``banner``, ``note`` and ``note_git`` strings are always kept. A numeric list is one
# leaf (its path is the path of the list) and is kept only if it has at most 2 entries (a range or an interval) or its path
# matches an entry of the module's ``GIT_ALLOW_LONG``. Containers left empty are dropped.
GIT_FORBIDDEN_KEYS: frozenset[str] = frozenset(
    {"data_mb_sr", "rel_err", "model_total_mb_sr", "model_ratio_minus_1", "rutherford_mb_sr",
     "family_total_mb_sr", "theta_cm_deg_points", "neighbour_ratio",
     "pdg_lowest_points", "pdg_lowest_energies_mev", "pdg_points_within_15pct_in_T"}
)  # fmt: skip
XSEC_KEY = re.compile(
    r"(^|_)(mb(_sr|_mev)?|barn)$"
)  # absolute cross-section names (sigma_el_mb, sigma_mb, data_mb_sr, ...)
ALWAYS_KEYS = ("banner", "note", "note_git")
BODY_MARKER = "<!-- body below is rendered from the Git projection of the JSON record -->"
GIT_MAX_NUMERIC_LIST = 2  # a numeric list longer than this is a distribution / per-point array unless explicitly allowed
GIT_NOTE = (
    "Git projection (decision 0041): explicit allowlist of dimensionless ratios, per-primary rates and probabilities, kinematic "
    "statistics, kernel integrals, counts and flags; no absolute cross sections, no distributions or per-point arrays. "
    "The complete record is regenerated into validation/generated/ by the committed module."
)


def _is_numeric_list(x: Any) -> bool:
    return (
        isinstance(x, list)
        and len(x) > 0
        and all(v is None or isinstance(v, (int, float)) for v in x)
    )


def path_matches(path: str, patterns: tuple[str, ...] | list[str]) -> bool:
    """True if the "/"-joined ``path`` matches a pattern (``*`` = one segment, trailing ``**`` = the rest)."""
    parts = path.split("/") if path else []
    for pat in patterns:
        pp = pat.split("/")
        if pp[-1] == "**":
            head = pp[:-1]
            if len(parts) > len(head) and all(
                h in ("*", p) for h, p in zip(head, parts, strict=False)
            ):
                return True
        elif len(pp) == len(parts) and all(h in ("*", p) for h, p in zip(pp, parts, strict=True)):
            return True
    return False


def _list_ok(path: str, value: list[Any], allow_long: tuple[str, ...]) -> bool:
    return len(value) <= GIT_MAX_NUMERIC_LIST or path_matches(path, allow_long)


def git_projection(
    doc: Any, allow: tuple[str, ...], allow_long: tuple[str, ...] = (), _path: str = ""
) -> Any:
    """Deep copy of ``doc`` keeping only the leaves whose path matches ``allow`` (see the policy above)."""
    keep = _project(doc, allow, allow_long, _path)
    return {} if keep is _DROP else keep


_DROP = object()


def _project(doc: Any, allow: tuple[str, ...], allow_long: tuple[str, ...], path: str) -> Any:
    if isinstance(doc, dict):
        out = {}
        for k, v in doc.items():
            here = f"{path}/{k}" if path else str(k)
            if not path and k in ALWAYS_KEYS and isinstance(v, str):
                out[k] = v
                continue
            if k in GIT_FORBIDDEN_KEYS or XSEC_KEY.search(str(k)):
                continue
            r = _project(v, allow, allow_long, here)
            if r is not _DROP:
                out[k] = r
        return out if out else _DROP
    if isinstance(doc, list) and not _is_numeric_list(doc):
        outl = [_project(v, allow, allow_long, f"{path}/{i}") for i, v in enumerate(doc)]
        outl = [v for v in outl if v is not _DROP]
        return outl if outl else _DROP
    if not path_matches(path, allow):
        return _DROP
    if _is_numeric_list(doc) and not _list_ok(path, doc, allow_long):
        return _DROP
    return doc


def git_safety_violations(
    doc: Any, allow: tuple[str, ...] | None = None, allow_long: tuple[str, ...] = (), path: str = ""
) -> list[str]:
    """Paths inside ``doc`` that a committed record must not contain: forbidden keys, absolute cross-section key names, numeric
    lists longer than 2 that are not in ``allow_long``, and (if ``allow`` is given) leaves that match no allowlist pattern."""
    bad: list[str] = []
    if isinstance(doc, dict):
        for k, v in doc.items():
            here = f"{path}/{k}" if path else str(k)
            if k in GIT_FORBIDDEN_KEYS or XSEC_KEY.search(str(k)):
                bad.append(here)
            else:
                bad += git_safety_violations(v, allow, allow_long, here)
    elif isinstance(doc, list) and not _is_numeric_list(doc):
        for i, v in enumerate(doc):
            bad += git_safety_violations(v, allow, allow_long, f"{path}/{i}")
    else:
        if _is_numeric_list(doc) and not _list_ok(path, doc, allow_long):
            bad.append(f"{path}[len={len(doc)}]")
        if allow is not None and path not in ALWAYS_KEYS and not path_matches(path, allow):
            bad.append(f"{path} (not in the allowlist)")
    return bad


def write_outputs(
    outdir: Path,
    json_name: str,
    doc: dict[str, Any],
    md_name: str,
    header_text: str,
    render: Callable[[dict[str, Any]], str],
    allow: tuple[str, ...],
    allow_long: tuple[str, ...] = (),
) -> None:
    """Write the complete JSON record to ``outdir`` and the Git-safe pair to ``outdir / "git"``.

    The Git copy keeps only the leaves of the module allowlist ``allow`` (policy above). The Markdown report is ``header_text``
    (command and input hashes) + ``BODY_MARKER`` + ``render(projection)``: its body is a pure function of the Git projection, so it
    cannot carry anything the JSON projection removed. Only the ``git`` copies may be committed under
    ``validation/results/elastic_attribution/`` (``tests/ionmc/test_elastic_attribution_de0.py`` checks every committed leaf against
    the owning module's allowlist and re-renders the committed Markdown from the committed JSON).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / json_name).write_text(json.dumps(doc, indent=1, sort_keys=True, default=float) + "\n")
    proj = git_projection(json.loads(json.dumps(doc, default=float)), allow, allow_long)
    proj.setdefault("note_git", GIT_NOTE)
    assert not git_safety_violations(proj, allow, allow_long), git_safety_violations(
        proj, allow, allow_long
    )
    (outdir / "git").mkdir(exist_ok=True)
    text = json.dumps(proj, indent=1, sort_keys=True) + "\n"
    (outdir / "git" / json_name).write_text(text)
    # render from the committed form (sorted keys), so a re-render of the committed JSON reproduces the body
    md = header_text + BODY_MARKER + "\n" + render(json.loads(text))
    (outdir / "git" / md_name).write_text(
        md.rstrip("\n") + "\n"
    )  # exactly one final newline (end-of-file-fixer)


@contextlib.contextmanager
def scoped_cache_env(path: Path) -> Iterator[Path]:
    """Set ``IONMC_CACHE_DIR`` to ``path`` for the scope only (the previous value is restored) and say so."""
    old = os.environ.get("IONMC_CACHE_DIR")
    print(
        f"IONMC_CACHE_DIR = {path} for this run (was {old if old is not None else 'unset'})",
        flush=True,
    )
    os.environ["IONMC_CACHE_DIR"] = str(path)
    try:
        yield Path(path)
    finally:
        if old is None:
            os.environ.pop("IONMC_CACHE_DIR", None)
        else:
            os.environ["IONMC_CACHE_DIR"] = old


def import_transport_module(name: str) -> Any:
    """Import a module of ``validation/scripts/transport`` (read only; no bytecode is written)."""
    sys.dont_write_bytecode = True
    p = str(REPO / "validation" / "scripts" / "transport")
    if p not in sys.path:
        sys.path.insert(0, p)
    import importlib

    return importlib.import_module(name)
