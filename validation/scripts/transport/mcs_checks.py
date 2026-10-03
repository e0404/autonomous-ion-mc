"""Multiple-scattering checks of the transport engine (V3-003A), aggregates only.

Usage::

    uv run python validation/scripts/transport/mcs_checks.py --output PATH [--code-sha SHA]

Requires the hash-verified LaTeX source of Gottschalk, arXiv:0908.1413, in
``.ionmc-cache/reference/<sha256>`` (materialized by the research tool; git-ignored). The
script fails closed: a missing or wrongly hashed file, or a ``--code-sha`` that differs from
``git rev-parse HEAD`` of this repository, exits non-zero and writes nothing.

Writes one JSON file with: the sha256 of the consumed source; U4 maximum absolute deviation of
the quadrature theta_dM from theta_Hanson (1 + dM %/100) per material (and per frozen x/R1,
deviations only, no table values), our rho R1 against the paper's; X_S deviations from the
table ``tbl:LS``; the full U5 table (relative deviation of the stepped variance sum from the
quadrature per step length and x/R1) and the negative control (relative spread of the per-step
Highland formula over the step lengths); the generalised-Highland cross-check for water; the
U4b comparison with theta_Hanson itself and the paper's maximum |dM %| per material. No table rows
of the paper are written. All physics comes from ``ionmc.transport.mcs_checks`` (the same
functions as the tests).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from ionmc.materials import ALUMINIUM, BERYLLIUM, COPPER, LEAD, WATER
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.transport import mcs_checks as mc
from ionmc.transport.tables import TransportTables

SOURCE_MODEL = "source-model reproduction (Gottschalk 2010 formulae and tables; not independent)"
REPO = Path(__file__).resolve().parents[3]


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(REPO), *args], check=True, capture_output=True, text=True
    ).stdout


def run() -> dict[str, Any]:
    tex = mc.gottschalk_tex()  # raises if missing or wrongly hashed
    rows, rho_r1, x_s = mc.parse_tables(tex)
    bethe = BetheStoppingSource()

    def tables(material: Any) -> TransportTables:
        return TransportTables.from_stopping_tables([bethe.table(material, PROTON)])

    u4: dict[str, Any] = {}
    u4b: dict[str, Any] = {}
    for key, material in {
        "Be": BERYLLIUM,
        "Al": ALUMINIUM,
        "Cu": COPPER,
        "Pb": LEAD,
    }.items():
        dev, our_r1 = mc.u4_deviations(tables(material), material, rows[key])
        dev_b, paper_dm = mc.u4b_deviations(tables(material), material, rows[key])
        u4b[key] = {
            "deviation_by_x_over_R1": {str(f): d for f, d in dev_b.items()},
            "max_abs_deviation": max(abs(d) for d in dev_b.values()),
            "paper_max_abs_dM_percent": paper_dm,
        }
        u4[key] = {
            "deviation_by_x_over_R1": {str(f): d for f, d in dev.items()},
            "max_abs_deviation": max(abs(d) for d in dev.values()),
            "rho_R1_ours_g_cm2": our_r1,
            "rho_R1_paper_g_cm2": rho_r1[key],
        }
    xs = mc.xs_deviations(
        x_s,
        {"Be": BERYLLIUM, "H2O": WATER, "Al": ALUMINIUM, "Cu": COPPER, "Pb": LEAD},
    )
    water = tables(WATER)
    u5 = mc.u5_deviations(water)
    return {
        "source": {
            "arxiv": "0908.1413",
            "sha256": mc.GOTTSCHALK_SHA256,
            "note": "aggregates only; no table rows of the paper are reproduced",
        },
        "evidence_class": {
            "u4": SOURCE_MODEL,
            "u4b": SOURCE_MODEL,
            "x_s": SOURCE_MODEL,
            "u5": "deterministic self-consistency",
            "highland_cross_check": "related model",
        },
        "u4": {
            "energy_mev": 158.6,
            "x_over_R1": list(mc.U4_FRACTIONS),
            "tolerance": mc.U4_TOLERANCE,
            "reference": "theta_Hanson * (1 + dM %/100)",
            "materials": u4,
            "max_abs_deviation_all": max(m["max_abs_deviation"] for m in u4.values()),
        },
        "u4b": {
            "reference": "theta_Hanson column (independent of the dM fit)",
            "tolerance": mc.U4B_TOLERANCE,
            "materials": u4b,
            "max_abs_deviation_all": max(m["max_abs_deviation"] for m in u4b.values()),
        },
        "x_s_relative_deviation_vs_tbl_LS": xs,
        "u5": {
            "material": "water",
            "energy_mev": 150.0,
            "criterion": 2e-3,
            "steps_mm": list(mc.STEPS_MM),
            "x_over_R1": list(mc.DEPTH_FRACTIONS),
            "relative_deviation": {
                str(frac): {str(s): u5[(frac, s)] for s in mc.STEPS_MM}
                for frac in mc.DEPTH_FRACTIONS
            },
            "max_abs_deviation": max(abs(v) for v in u5.values()),
            "negative_control": mc.u5_negative_control(water),
        },
        "generalised_highland_cross_check_water": {
            str(f): d for f, d in mc.highland_cross_check(water, WATER).items()
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--code-sha", default=None, help="must equal `git rev-parse HEAD`")
    args = ap.parse_args()
    head = _git("rev-parse", "HEAD").strip()
    if args.code_sha is not None and args.code_sha != head:
        print(f"--code-sha {args.code_sha} != HEAD {head}", file=sys.stderr)
        return 2
    dirty = bool(_git("status", "--porcelain").strip())
    try:
        result = run()
    except (FileNotFoundError, ValueError) as exc:
        print(f"fail closed: {exc}", file=sys.stderr)
        return 1
    result["provenance"] = {"analysis_code_sha": head, "analysis_code_dirty": dirty}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
