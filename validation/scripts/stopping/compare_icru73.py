"""Row I3: ionmc Bethe (I = 78 eV) for carbon and oxygen vs libdEdx ICRU73_NEW liquid water.

Usage::

    uv run python validation/scripts/stopping/compare_icru73.py --archive PATH \
        --output-dir DIR --code-sha SHA

Reference (record ``validation/plans/v3-008/p5b-libdedx-icru73.md``): the libdEdx source
archive (gzip of a tar, sha256 ``ARCHIVE_SHA256``; GPL-3 data, used as a cached evaluation
reference outside Git). The archive is untrusted data: its sha256 is verified, only the two
members below are extracted (``tarfile`` filter ``"data"``) into a fresh temporary directory,
their sha256 are verified before parsing, and nothing from the archive is executed.

* ``data/raw/icru73_newEng.dat``: first non-blank line the number of energies (53), then one
  energy per line, kinetic energy per nucleon [MeV/u] (``tools/dat2c.py`` registers
  ``icru73new`` with energy scale 1.0 and the generated ``dedx_icru73new.h`` states "Energy
  values are in MeV/nucl").
* ``data/raw/ICRU73_NEW.dat``: blocks ``#<target>:<ion Z>:<n>`` followed by ``n`` values, one
  per line, blank lines between blocks; targets 104 (air) and 276 (liquid water), ions Z = 1-18.
  The values are mass electronic stopping powers in MeV cm2/g: libdEdx ``README.md`` line 41
  ("stp in MeV cm² / g at 100 MeV/nucl") for the library output, which for ``DEDX_ICRU73`` is
  the table generated from this file without value scaling (``tools/dat2c.py`` ``read_dat``;
  ``src/dedx_embedded_data.c`` maps ``dedx_embedded_icru73`` to ``dedx_icru73new_*``).

Comparison (row I3 of ``validation/plans/v3-008-acceptance.md``): the ionmc Bethe model with
the default :class:`BetheOptions` (``ls = "off"``) in ``water(78.0)`` for ``CARBON12`` and
``OXYGEN16``, evaluated analytically on the table's own grid points with E >= 2 MeV/u (the
full 53-point grid above 2 MeV/u, up to its 1000 MeV/u end; an informative 10-500 MeV/u
aggregate is added because the Bethe tables used for transport end at 500 MeV/u). Stopping
aggregates (max abs and RMS relative deviation Bethe/table - 1, number of points, energy range)
for 2 <= E < 10 MeV/u and E >= 10 MeV/u; gates: E >= 10 MeV/u max <= 1.5 % and RMS <= 0.8 %,
2-10 MeV/u max <= 6 %. CSDA range: the table gives no range, so the table is log-log
interpolated onto a 200-points-per-decade grid from 2 MeV/u and integrated with the same exact
log-log quadrature as the analytic tables (``build_table``), **started at the ionmc model's own
CSDA range at 2 MeV/u** (Bethe 78 eV table, which starts at 1 MeV/u with the constant-S
approximation). It is compared with the Bethe (78 eV) table range at the ends and standard
points (C: 80, 290, 430 MeV/u; O: 100, 350, 430 MeV/u); gate max abs rel dev <= 0.5 %; only the
maximum and the energy attaining it are written.

The JSON output contains aggregates only; no libdEdx/ICRU table values, per-energy deviations or
ratios are written (field ``redistribution_note``).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.data.cache import IntegrityError
from ionmc.materials import water
from ionmc.physics.projectiles import CARBON12, OXYGEN16, Projectile
from ionmc.physics.stopping import (
    BetheOptions,
    BetheStoppingSource,
    bethe_mass_stopping,
    build_table,
    log_grid,
)

ARCHIVE_SHA256 = "b941a0c30793daddc4d119a2cc664bfbdd430554c3a4788be04e97ff6fcc0863"
STOPPING_MEMBER = "data/raw/ICRU73_NEW.dat"
ENERGY_MEMBER = "data/raw/icru73_newEng.dat"
MEMBER_SHA256 = {
    STOPPING_MEMBER: "749768c1c2de475c928cc2c50a354903577139c11debd4b5948b035ca3a3d87e",
    ENERGY_MEMBER: "8b5b809f5c0cc4f6018f9faa89e9f73c7b9b81a3c588f98bef556712d4c77dd1",
}
TARGET_WATER = 276  # libdEdx / NIST material id of liquid water
N_GRID = 53
I_EV_MODEL = 78.0
LIBDEDX_I_EV_WATER = 75.0  # libdEdx dedx_metadata.h, target 276 (material constant)

E_COMPARE_MIN = 2.0  # MeV/u
E_BAND_SPLIT = 10.0  # MeV/u
E_INFORMATIVE_MAX = 500.0  # MeV/u, end of the transport Bethe tables
GATE_S_HIGH_MAX = 0.015
GATE_S_HIGH_RMS = 0.008
GATE_S_LOW_MAX = 0.06
GATE_RANGE_MAX = 0.005
RANGE_POINTS_PER_DECADE = 200
IONS: dict[str, tuple[Projectile, tuple[float, ...]]] = {
    "carbon12": (CARBON12, (80.0, 290.0, 430.0)),
    "oxygen16": (OXYGEN16, (100.0, 350.0, 430.0)),
}

I_VALUE_NOTE = (
    "The ionmc model is evaluated at I = 78 eV (water(78.0), not tuned). libdEdx lists "
    "I = 75.0 eV for target 276 (dedx_metadata.h material constant). Whether the ICRU73_NEW "
    "water values were computed with 75 eV or 78 eV: TO-VERIFY (P5b); the archive documentation "
    "(README.md, data/README.md, tools/dat2c.py, the dedx_icru73new.h header) states no I-value "
    "for this table. If 75 eV, the known bounded contribution of the difference is about "
    "0.4-0.6 % in S above 10 MeV/u and 0.4-0.5 % in CSDA range (plan row I3, L0 ~ -ln I; model S "
    "lower and range longer at 78 eV); a range-clause failure within that bound is recorded as "
    "FAIL with this attribution, never re-graded."
)
REDISTRIBUTION_NOTE = (
    "This file contains aggregate statistics only; no libdEdx/ICRU 73 table values, per-energy "
    "deviations or ratios are reproduced."
)


@dataclass(frozen=True)
class Icru73Water:
    """ICRU73_NEW liquid-water electronic stopping of one ion: energies [MeV/u], S [MeV cm2/g]."""

    z: int
    energy_per_u_mev: NDArray[np.float64]
    stopping_mev_cm2_g: NDArray[np.float64]


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_energy_file(text: str) -> NDArray[np.float64]:
    """Parse ``icru73_newEng.dat``: count line, then that many positive increasing energies."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("empty energy file")
    n = int(lines[0])
    values = np.asarray([float(v) for v in lines[1:]], dtype=np.float64)
    if values.size != n:
        raise ValueError(f"energy file declares {n} values, found {values.size}")
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0) or np.any(np.diff(values) <= 0.0):
        raise ValueError("energies must be finite, positive and strictly increasing")
    return values


def parse_stopping_file(text: str) -> dict[tuple[int, int], NDArray[np.float64]]:
    """Parse ``ICRU73_NEW.dat`` into ``{(target, ion Z): values}``.

    Each block is a header ``#<target>:<ion>:<n>`` followed by exactly ``n`` numeric lines;
    a short block, values outside a block, a duplicate (target, ion) or a malformed header
    raise ``ValueError``.
    """
    blocks: dict[tuple[int, int], NDArray[np.float64]] = {}
    key: tuple[int, int] | None = None
    n_expected = 0
    buf: list[float] = []

    def close() -> None:
        if key is not None:
            if len(buf) != n_expected:
                raise ValueError(f"block {key}: declared {n_expected} values, found {len(buf)}")
            blocks[key] = np.asarray(buf, dtype=np.float64)

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            close()
            parts = line[1:].split(":")
            if len(parts) != 3:
                raise ValueError(f"malformed block header {line!r}")
            target, ion, n_expected = (int(p) for p in parts)
            key = (target, ion)
            if key in blocks:
                raise ValueError(f"duplicate block {key}")
            buf = []
        else:
            if key is None:
                raise ValueError("value before the first block header")
            if len(buf) >= n_expected:
                raise ValueError(f"block {key}: more than {n_expected} values")
            buf.append(float(line))
    close()
    return blocks


def _member(tar: tarfile.TarFile, path: str) -> tarfile.TarInfo:
    """The unique regular-file member whose path, without the top-level directory, is ``path``."""
    hits = [
        m
        for m in tar.getmembers()
        if m.isfile() and PurePosixPath(*PurePosixPath(m.name).parts[1:]) == PurePosixPath(path)
    ]
    if len(hits) != 1:
        raise IntegrityError(f"archive must contain exactly one member */{path}, found {len(hits)}")
    return hits[0]


def read_members(
    archive: Path,
    *,
    archive_sha256: str = ARCHIVE_SHA256,
    member_sha256: dict[str, str] | None = None,
) -> dict[str, bytes]:
    """Verify the archive hash, extract the two members safely and verify their hashes."""
    member_sha256 = MEMBER_SHA256 if member_sha256 is None else member_sha256
    raw = Path(archive).read_bytes()
    if _sha256(raw) != archive_sha256:
        raise IntegrityError(f"archive sha256 {_sha256(raw)} != expected {archive_sha256}")
    out: dict[str, bytes] = {}
    with tempfile.TemporaryDirectory(prefix="ionmc-libdedx-") as tmp:
        with tarfile.open(archive, "r:gz") as tar:
            members = {path: _member(tar, path) for path in member_sha256}
            tar.extractall(tmp, members=list(members.values()), filter="data")
        for path, info in members.items():
            data = (Path(tmp) / info.name).read_bytes()
            if _sha256(data) != member_sha256[path]:
                raise IntegrityError(f"member {path} sha256 {_sha256(data)} != expected")
            out[path] = data
    return out


def load_icru73_water(
    archive: Path,
    *,
    archive_sha256: str = ARCHIVE_SHA256,
    member_sha256: dict[str, str] | None = None,
    ions: tuple[int, ...] = (6, 8),
) -> dict[int, Icru73Water]:
    """Liquid-water (target 276) ICRU73_NEW tables of the ions ``ions`` from the archive."""
    data = read_members(archive, archive_sha256=archive_sha256, member_sha256=member_sha256)
    energy = parse_energy_file(data[ENERGY_MEMBER].decode("ascii"))
    if energy.size != N_GRID:
        raise ValueError(f"expected a {N_GRID}-point energy grid, found {energy.size}")
    blocks = parse_stopping_file(data[STOPPING_MEMBER].decode("ascii"))
    out = {}
    for z in ions:
        s = blocks.get((TARGET_WATER, z))
        if s is None:
            raise ValueError(f"no block for target {TARGET_WATER}, ion Z = {z}")
        if s.size != energy.size:
            raise ValueError(f"block ({TARGET_WATER}, {z}) has {s.size} values, grid {energy.size}")
        if not np.all(np.isfinite(s)) or np.any(s <= 0.0):
            raise ValueError(f"block ({TARGET_WATER}, {z}) has non-positive or non-finite values")
        out[z] = Icru73Water(z, energy.copy(), s.copy())
    return out


def _band(e_u: NDArray[np.float64], dev: NDArray[np.float64], mask: NDArray[np.bool_]) -> dict:
    n = int(np.count_nonzero(mask))
    if n == 0:
        return {
            "n_points": 0,
            "energy_per_u_range_mev": None,
            "max_abs_relative_deviation": None,
            "rms_relative_deviation": None,
        }
    return {
        "n_points": n,
        "energy_per_u_range_mev": [float(e_u[mask].min()), float(e_u[mask].max())],
        "max_abs_relative_deviation": float(np.max(np.abs(dev[mask]))),
        "rms_relative_deviation": float(np.sqrt(np.mean(dev[mask] ** 2))),
    }


def _le(value: float | None, bound: float) -> bool:
    return value is not None and bool(np.isfinite(value)) and value <= bound


def compare_ion(proj: Projectile, table: Icru73Water, range_energies: tuple[float, ...]) -> dict:
    """Aggregates and gate booleans of one ion (no per-energy values)."""
    if table.z != proj.z:
        raise ValueError(f"table Z = {table.z} does not match projectile z = {proj.z}")
    o = BetheOptions()
    material = water(I_EV_MODEL)
    sel = table.energy_per_u_mev >= E_COMPARE_MIN
    e_u = table.energy_per_u_mev[sel]
    s_ref = table.stopping_mev_cm2_g[sel]
    s_model = bethe_mass_stopping(
        e_u,
        proj,
        material,
        shell=o.shell,
        barkas=o.barkas,
        bloch=o.bloch,
        mott=o.mott,
        effective_charge=o.effective_charge,
        ls=o.ls,
    )
    dev = s_model / s_ref - 1.0
    low = _band(e_u, dev, e_u < E_BAND_SPLIT)
    high = _band(e_u, dev, e_u >= E_BAND_SPLIT)
    informative = _band(e_u, dev, (e_u >= E_BAND_SPLIT) & (e_u <= E_INFORMATIVE_MAX))

    t_model = BetheStoppingSource(options=o).table(material, proj)
    full = table.energy_per_u_mev
    e_fine = log_grid(E_COMPARE_MIN, float(full[-1]), RANGE_POINTS_PER_DECADE)
    s_fine = np.exp(np.interp(np.log(e_fine), np.log(full), np.log(table.stopping_mev_cm2_g)))
    r_start = float(t_model.range_at(E_COMPARE_MIN))
    t_ref = build_table(proj, material, e_fine, s_fine, r_start, {"source": "libdedx-icru73new"})
    rdev = [float(t_model.range_at(e)) / float(t_ref.range_at(e)) - 1.0 for e in range_energies]
    k = int(np.argmax(np.abs(rdev)))
    rmax = float(abs(rdev[k]))

    gates = {
        f"stopping_ge_10_MeV_u_max_abs_le_{GATE_S_HIGH_MAX:g}": _le(
            high["max_abs_relative_deviation"], GATE_S_HIGH_MAX
        ),
        f"stopping_ge_10_MeV_u_rms_le_{GATE_S_HIGH_RMS:g}": _le(
            high["rms_relative_deviation"], GATE_S_HIGH_RMS
        ),
        f"stopping_2_to_10_MeV_u_max_abs_le_{GATE_S_LOW_MAX:g}": _le(
            low["max_abs_relative_deviation"], GATE_S_LOW_MAX
        ),
        f"csda_range_max_abs_le_{GATE_RANGE_MAX:g}": _le(rmax, GATE_RANGE_MAX),
    }
    return {
        "z": proj.z,
        "a": proj.a,
        "stopping": {
            "2_to_10_MeV_u": {"band_per_u_mev": [E_COMPARE_MIN, E_BAND_SPLIT], **low},
            "ge_10_MeV_u": {"band_per_u_mev": [E_BAND_SPLIT, float(full[-1])], **high},
            "informative_10_to_500_MeV_u": {
                "band_per_u_mev": [E_BAND_SPLIT, E_INFORMATIVE_MAX],
                **informative,
            },
        },
        "csda_range": {
            "energies_per_u_mev": list(range_energies),
            "max_abs_relative_deviation": rmax,
            "energy_per_u_mev_at_max": float(range_energies[k]),
        },
        "gates": gates,
        "pass": all(gates.values()),
    }


def run(
    archive: Path,
    code_sha: str | None,
    *,
    archive_sha256: str = ARCHIVE_SHA256,
    member_sha256: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Compute the aggregate I3 comparison as a JSON-serialisable dictionary."""
    from ionmc import __version__

    member_sha256 = MEMBER_SHA256 if member_sha256 is None else member_sha256
    tables = load_icru73_water(
        archive,
        archive_sha256=archive_sha256,
        member_sha256=member_sha256,
        ions=tuple(p.z for p, _ in IONS.values()),
    )
    per_ion = {name: compare_ion(p, tables[p.z], energies) for name, (p, energies) in IONS.items()}
    o = BetheOptions()
    return {
        "row": "I3",
        "reference": "libdEdx ICRU73_NEW, target 276 (liquid water), Z = 6 and 8",
        "redistribution_note": REDISTRIBUTION_NOTE,
        "i_value_note": I_VALUE_NOTE,
        "thresholds": {
            "stopping_ge_10_MeV_u_max_abs": GATE_S_HIGH_MAX,
            "stopping_ge_10_MeV_u_rms": GATE_S_HIGH_RMS,
            "stopping_2_to_10_MeV_u_max_abs": GATE_S_LOW_MAX,
            "csda_range_max_abs": GATE_RANGE_MAX,
        },
        "range_construction": (
            "table log-log interpolated to 200 points/decade from 2 MeV/u, exact log-log "
            "quadrature (build_table), started at the ionmc Bethe 78 eV CSDA range at 2 MeV/u; "
            "model: BetheStoppingSource table (200 points/decade from 1 MeV/u, constant-S start)"
        ),
        "ions": per_ion,
        "pass": all(v["pass"] for v in per_ion.values()),
        "provenance": {
            "archive_sha256": archive_sha256,
            "member_sha256": dict(member_sha256),
            "target": TARGET_WATER,
            "code_sha": code_sha,
            "ionmc_version": __version__,
            "numpy_version": np.__version__,
            "model": {
                "source": "bethe",
                "I_eV": I_EV_MODEL,
                "material": "water",
                "options": {
                    "shell": o.shell,
                    "barkas": o.barkas,
                    "bloch": o.bloch,
                    "mott": o.mott,
                    "effective_charge": o.effective_charge,
                    "ls": o.ls,
                },
            },
        },
    }


def main() -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--archive", required=True, help="libdEdx archive (sha256-verified)")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--code-sha",
        default=None,
        help="git SHA of the clean checkout that produced this comparison (recorded verbatim)",
    )
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    result = run(Path(args.archive), args.code_sha)
    path = out / "icru73_comparison.json"
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
