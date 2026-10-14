"""Re-evaluate an archived V5 comparator verdict (``v5-compare`` step output) under the corrected Welch df.

Usage::

    v5_df_recheck.py ARCHIVED_STEP_OUTPUT.txt [--json]

The archived step output (``NN-v5-compare.txt`` of an lv5b archive) holds, per comparison, the ionmc and
reference estimates (value, se, df), the difference with its df and the TOST outcome. Nothing is read
from or written to the archive or ``validation/results``; the table goes to stdout.

Correction (C7 review, 2026-10-11): for relative differences ``(a - b) / b`` the standard error is
``hypot(a.se, r b.se) / |b|`` with ``r = a / b``, and the Welch-Satterthwaite df must use the same
components ``a.se`` and ``|r| b.se``; the recorded code used the unscaled ``b.se`` in the df. This script
recomputes the df and the 90 % interval with the corrected formula from the recorded estimates (the
difference SE and the floored t table are unchanged) and reports whether any pass/fail component
changes. Absolute differences (R80) are unaffected.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent


def _cmp() -> Any:
    spec = importlib.util.spec_from_file_location("compare_idd_v5", HERE / "compare_idd_v5.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("compare_idd_v5", mod)
    spec.loader.exec_module(mod)
    return mod


def load_document(path: Path) -> dict[str, Any]:
    text = path.read_text()
    a, b = text.index("#JSON-BEGIN"), text.index("#JSON-END")
    doc = json.loads(text[a + len("#JSON-BEGIN"):b])
    if not isinstance(doc.get("verdict"), dict):
        raise SystemExit(f"{path}: no verdict in the archived document (error: {doc.get('error')!r})")
    return doc


def recheck(verdict: dict[str, Any], m: Any) -> list[dict[str, Any]]:
    rows = []
    for engine, ev in verdict["engines"].items():
        for e, crit in ev["energies"].items():
            for name, r in crit.items():
                if name.startswith("_"):
                    continue
                ion, ref = m.Est(**r["ionmc"]), m.Est(**r["reference"])
                relative = r["kind"] == "relative"
                # the recorded df, recomputed with the recorded (old) formula as a consistency check
                old_df = m.welch_df(ion.se, ion.df, ref.se, ref.df)
                new_df = old_df
                if relative:
                    ratio = ion.value / ref.value
                    new_df = m.welch_df(ion.se, ion.df, abs(ratio) * ref.se, ref.df)
                tol = r["tolerance"]
                lo0, hi0 = r["ci90"]
                q = m.t95(new_df)
                lo1, hi1 = r["diff"] - q * r["se"], r["diff"] + q * r["se"]
                p1 = bool(lo1 > -tol and hi1 < tol)
                rows.append({"engine": engine, "energy": e, "criterion": name, "kind": r["kind"],
                             "recorded_df": r["df"], "old_df_recomputed": old_df, "new_df": new_df,
                             "recorded_t95": r["t95"], "new_t95": q, "old_ci90": [lo0, hi0],
                             "new_ci90": [lo1, hi1], "tolerance": tol, "old_pass": r["pass"],
                             "new_pass": p1, "changed": bool(p1 != r["pass"])})  # fmt: skip
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archive", type=Path)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)
    m = _cmp()
    doc = load_document(args.archive)
    rows = recheck(doc["verdict"], m)
    if args.json:
        print(json.dumps(rows, indent=1))
        return 0
    print(f"{'engine':8s} {'E':>3s} {'criterion':28s} {'df old':>7s} {'df new':>7s}  "
          f"{'old 90% CI':>24s}  {'new 90% CI':>24s}  old new chg")  # fmt: skip
    for r in rows:
        tag = "CHANGED" if r["changed"] else "-"
        print(f"{r['engine']:8s} {r['energy']:>3s} {r['criterion']:28s} {r['recorded_df']:7.2f} {r['new_df']:7.2f}  "
              f"[{r['old_ci90'][0]:+.5f},{r['old_ci90'][1]:+.5f}]  [{r['new_ci90'][0]:+.5f},{r['new_ci90'][1]:+.5f}]  "
              f"{'P' if r['old_pass'] else 'F'}   {'P' if r['new_pass'] else 'F'}   {tag}")  # fmt: skip
    bad = [r for r in rows if abs(r["old_df_recomputed"] - r["recorded_df"]) > 1e-9 * max(1.0, r["recorded_df"])]
    if bad:
        print(f"NOTE: {len(bad)} recorded df values differ from the old formula recomputed from the recorded estimates")
    print(f"components changed: {sum(r['changed'] for r in rows)} of {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
