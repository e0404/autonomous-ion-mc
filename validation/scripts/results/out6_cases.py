"""The eight ``out6`` cases of plan v3-006 (step S8): tiny, synthetic, not a clinical claim.

``seed = OUT6_BASE + 100 * r_index + k`` is the only source of seeds here (plan v3-006 "Seeds":
``k <= 99``, ``r_index <= 49``). k = 99 is reserved for tests, k = 0 and k = 1 for the LV/HR
evidence runs. The TOML templates in ``out6/`` carry the k = 0 seed; :func:`write_config` rewrites
the seed (and, for an override, the backend and precision), resolves the nuclear table id of the
cache for r_index 5 and copies the referenced files of r_index 3 beside the written TOML.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

OUT6_BASE = 20700004
OUT6_R_MAX = 49
OUT6_K_MAX = 99
HERE = Path(__file__).resolve().parent
TEMPLATES = HERE / "out6"
REPO = HERE.parents[2]
SLAB_FILES = (
    REPO / "tests" / "ionmc" / "data" / "v3_006" / "slab_material_index.npy",
    REPO / "tests" / "ionmc" / "data" / "v3_006" / "slab_density.npy",
    REPO / "tests" / "data" / "synthetic_lookup.json",
)

# r_index -> (template, backend, precision, notes)
_CASES: dict[int, tuple[str, str, str, str]] = {
    0: ("r0_water_box_python.toml", "python", "float64", "proton water box, CLI water box"),
    1: ("r1_water_box_warp_cpu.toml", "warp-cpu", "float32", "the r0 box on warp-cpu"),
    2: ("r2_water_box_warp_cuda.toml", "warp-cuda", "float32", "the r0 box on CUDA (HR only)"),
    3: ("r3_voxel_slab.toml", "warp-cpu", "float64", "water/PMMA voxel slab, LET, lookup tally"),
    4: ("r4_carbon_box.toml", "warp-cpu", "float64", "carbon-12 box, species identity"),
    5: ("r5_nuclear_proton.toml", "python", "float64", "nuclear proton (cached derived table)"),
    6: ("r6_reproduction_box.toml", "python", "float64", "OUT-3 reproduction originals (k 0, 1)"),
    7: ("r7_icru90_alpha.toml", "warp-cpu", "float64", "icru90-water-v1 alpha box (cache)"),
}


def seed_of(r_index: int, k: int) -> int:
    if not (0 <= r_index <= OUT6_R_MAX and 0 <= k <= OUT6_K_MAX):
        raise ValueError(f"out6 index out of range: r_index={r_index}, k={k}")
    return OUT6_BASE + 100 * r_index + k


def cases() -> list[tuple[Path, int, str, str, str]]:
    """The eight cases as ``(toml path, seed at k = 0, backend, precision, notes)``."""
    return [
        (TEMPLATES / name, seed_of(r, 0), backend, precision, notes)
        for r, (name, backend, precision, notes) in sorted(_CASES.items())
    ]


def nuclear_table_id(cache_dir: str | Path) -> str:
    """The id of a built, qualified nuclear table of ``cache_dir``; fails when there is none."""
    from ionmc.nuclear.tables import NuclearTable, NuclearTableError

    for p in sorted((Path(cache_dir) / "derived").glob("nuclear-proton-*.json")):
        tid = p.stem.removeprefix("nuclear-proton-")
        try:
            if "transport_path_bound_terms" in NuclearTable.load(Path(cache_dir), tid).info:
                return tid
        except NuclearTableError:
            continue
    raise RuntimeError("no built nuclear table in the cache")


def _set_key(text: str, key: str, literal: str) -> str:
    out, n = re.subn(rf"^{key} = .*$", f"{key} = {literal}", text, flags=re.MULTILINE)
    if n != 1:
        raise ValueError(f"template has {n} lines for {key!r}")
    return out


def write_config(
    r_index: int,
    k: int,
    out_dir: str | Path,
    cache_dir: str | Path,
    *,
    backend: str | None = None,
    precision: str | None = None,
) -> Path:
    """Write ``out_dir/r<r>k<k>.toml`` from the template of ``r_index`` and return its path."""
    name = _CASES[r_index][0]
    text = (TEMPLATES / name).read_text(encoding="utf-8")
    text = _set_key(text, "seed", str(seed_of(r_index, k)))
    if backend is not None:
        text = _set_key(text, "backend", f'"{backend}"')
    if precision is not None:
        text = _set_key(text, "precision", f'"{precision}"')
    if "@TABLE_ID@" in text:
        text = text.replace("@TABLE_ID@", nuclear_table_id(cache_dir))
    dest_dir = Path(out_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    if r_index == 3:
        for src in SLAB_FILES:
            shutil.copyfile(src, dest_dir / src.name)
    dest = dest_dir / f"r{r_index}k{k}.toml"
    dest.write_text(text, encoding="utf-8")
    return dest
