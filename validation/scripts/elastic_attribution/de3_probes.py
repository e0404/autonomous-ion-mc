"""D-E3 probes (exploratory, report-only; Amendment 19 (b), (f)): the V11 harness on the D-E3 variant tables.

One invocation = one transport run (warp-cpu float64, single process, V5 geometry r = 20 cm, 20 batches,
0.5 mm depth bins) of one variant (``--variant``; ``emonly`` = EM only, ``pinned`` = the pinned table) at one energy,
with the common seed of the energy: **20990100 (150 MeV) and 20990101 (200 MeV)**, shared by every variant and by
EM-only (CRN). ``--summary`` combines the run files of an output directory into F, gain, recoil, R80, lateral rms and the K*
prediction. ``--check-base`` is the pre-D-E3 check (plan 3.1): the base copy ``all`` must reproduce the pinned
table bit for bit at the same seed.

Reused definitions (existing harness, not re-derived):

- configuration and run: ``steps_v5c.v11_config`` (l.567-577: warp-cpu float64, ``nuc_config_c_elastic_only`` for EM+elastic,
  ``nuc_config_c(nuclear=False)`` for EM-only), ``steps_v5b.v5_geometry`` (l.587), ``steps_v5b.V5_BATCHES / V5_DZ_MM``;
- depth-dose curves, R80, plateau: ``steps_v5c.step_v11_ionmc`` (l.580-607; ``idd_b = batch_energy / bin_g_cm2``, plateau 20-60 mm);
- F = (peak/plateau)_EM+el / (peak/plateau)_EM and gain = plateau ratio - 1: ``compare_idd_v5.paired_ratio_jackknife``
  (l.680-700) and ``v11_energy`` (l.708-716) with ``curve_metrics`` (l.167);
- recoil channel: ``steps_v5d.REC`` (l.84) and ``steps_v5b.tally_batches`` (steps_v5d.step_v5_ionmc l.660-673);
- lateral rms: ``steps_v5d.r1d2_measure`` (l.342-376) grid (+-50 mm, 1 mm, slab 0.45-0.55 R80) and ``lateral_rms_batches`` (l.317).

Only the elastic table differs from the V11 producer: the variant is loaded through
``PhysicsConfig(elastic_table_id=<variant id>, allow_exploratory_table=True)`` (the pinned id keeps the evidence flags).
The tables come from the scratch cache of ``make_variant_tables`` (``--scratch-cache DIR``, else the environment variable
``IONMC_DE3_SCRATCH``, else ``validation/generated/elastic_attribution/de3-cache/``, git-ignored). ``IONMC_CACHE_DIR`` is set to it
only for the scope of the run (printed, then restored); an already-set value is never overwritten silently, and provenance
records only the placeholder ``<scratch-cache>``.

Run: ``python -m validation.scripts.elastic_attribution.de3_probes --variant pA-s0.5 --energy 150 [--histories 100000]``.
"""

from __future__ import annotations

import os
import argparse  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402

from . import common as c  # noqa: E402

SEEDS = {150: 20990100, 200: 20990101}  # Amendment 19 (f)
GIT_ALLOW = (
    "bit_for_bit", "energy", "n", "n_batches", "seed", "pinned_id", "all_id", "detail/*", "counters/*",
    "sha256_idd_batches/*", "sha256_recoil_idd_batches/*",
    "rows/*/*/*", "runs/*/*/*", "derived/**",
    "provenance/git_head", "provenance/git_status_clean", "provenance/git_status_elastic_attribution/*", "provenance/run_logs/**", "provenance/run_groups/**", "provenance/scripts_sha256/*",
    "provenance/table_files_sha256/*", "provenance/backend", "provenance/precision",
    "provenance/single_process", "provenance/cache_dir", "provenance/seeds/*",
)  # fmt: skip
"""Git allowlist of ``de3-summary.json`` and ``check-base-*.json`` (syntax: :mod:`common`): F, gain, peak removal, R80, recoil
and lateral statistics per variant (simulation output, no nuclear data), run bookkeeping and provenance (code hashes; host
paths are never recorded, the scratch cache appears as ``<scratch-cache>``)."""
PINNED = c.TABLE_ID
GIT_FILES = {"de3-summary.json": (GIT_ALLOW, ()), "check-base-*.json": (GIT_ALLOW, ())}
OUT_DEFAULT = c.OUT / "de3"
R1D2_HALF_MM, R1D2_BIN_MM, R1D2_SLAB = 50.0, 1.0, (0.45, 0.55)  # steps_v5d.py l.88


def _t() -> Any:
    """The transport step modules (read only)."""
    v5c = c.import_transport_module("steps_v5c")
    v5d = c.import_transport_module("steps_v5d")
    return v5c, v5d, c.import_transport_module("steps_v5b"), c.import_transport_module("steps_v5")


def _cmp() -> Any:
    path = c.REPO / "validation" / "scripts" / "reference" / "compare_idd_v5.py"
    spec = importlib.util.spec_from_file_location("compare_idd_v5", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def variant_table_id(variant: str) -> str:
    if variant in ("pinned", "emonly"):
        return PINNED
    from .make_variant_tables import variant_ids

    return variant_ids(Path(os.environ["IONMC_CACHE_DIR"]))[variant]


def build_config(
    energy: int,
    mode: str,
    n: int,
    seed: int,
    table_id: str,
    batches: int = 20,
    grid_geo: Any = None,
    tallies: tuple[Any, ...] | None = None,
) -> Any:
    """V11 configuration (``steps_v5c.v11_config``) with the elastic table ``table_id``. ``mode``: ``emel`` or ``emonly``."""
    from dataclasses import replace as rep

    v5c, v5d, v5b, _ = _t()
    geo, grid, _nz, _ = v5b.v5_geometry(float(energy))
    if grid_geo is not None:
        geo, grid = grid_geo
    if mode == "emonly":
        cfg = v5c.nuc_config_c(energy=float(energy), n=n, seed=seed, geometry=geo, grid=grid, nuclear=False,
                               n_batches=batches, tallies=tallies or ())  # fmt: skip
    else:
        cfg = v5c.nuc_config_c_elastic_only(energy=float(energy), n=n, seed=seed, geometry=geo, grid=grid, nuclear=True,
                                            n_batches=batches, tallies=tallies if tallies is not None else (v5d.REC,))  # fmt: skip
        if table_id != PINNED:
            cfg = rep(
                cfg,
                physics=rep(cfg.physics, elastic_table_id=table_id, allow_exploratory_table=True),
            )
    return rep(cfg, run=rep(cfg.run, backend="warp-cpu", precision="float64"))


def run_one(
    variant: str, energy: int, n: int, batches: int, lateral: bool, out: Path
) -> dict[str, Any]:
    """One run -> ``out/de3-<energy>-<variant>.json``."""
    v5c, v5d, v5b, v5 = _t()
    from ionmc import materials as M
    from ionmc.reference import metrics as rm
    from ionmc.scoring import ScoringGrid

    mode = "emonly" if variant == "emonly" else "emel"
    tid = variant_table_id(variant)
    seed = SEEDS[energy]
    geo, grid, nz, r_mm = v5b.v5_geometry(float(energy))
    cfg = build_config(energy, mode, n, seed, tid, batches)
    t0 = time.perf_counter()
    res = v5b.Simulation(cfg).run()
    wall = time.perf_counter() - t0
    rho = M.WATER.density_g_cm3
    bin_g = rho * v5b.V5_DZ_MM / 10.0
    idd_b = res.grid("dose").batch_energy_mev.reshape(batches, nz) / bin_g
    idd = idd_b.mean(axis=0)
    depth = (np.arange(nz) + 0.5) * v5b.V5_DZ_MM
    curve = rm.normalize_to_peak(idd)
    plateau = (depth >= 20.0) & (depth <= 60.0)
    doc: dict[str, Any] = {
        "banner": c.BANNER, "variant": variant, "mode": mode, "energy": energy, "seed": seed, "n": n, "n_batches": batches,
        "elastic_table_id": tid, "backend": cfg.run.backend, "precision": cfg.run.precision, "scoring": {"origin_mm": [float(x) for x in grid.origin_mm], "spacing_mm": [float(x) for x in grid.spacing_mm],
                    "shape": [int(x) for x in grid.shape], "name": grid.name}, "tallies": [q.name for q in cfg.tallies], "bin_mm": v5b.V5_DZ_MM, "idd_batches": idd_b.tolist(), "unit": "MeV/(g/cm^2)/primary",
        "metrics": {"peak_depth_mm": rm.peak_depth(depth, curve), "r80_mm": rm.r80(depth, curve),
                    "plateau_mean_20_60_mm": float(idd[plateau].mean()), "peak_over_plateau": float(idd.max() / idd[plateau].mean())},
        "valid": bool(v5.clean(res) and res.energy_balance.relative_residual <= 1e-12),
        "relative_residual": float(res.energy_balance.relative_residual),
        "counters": res.counters.as_dict(), "wall_s": wall, "hist_per_s": n / wall,
    }  # fmt: skip
    if mode == "emel":
        rec_b = v5b.tally_batches(res, "recoil") / bin_g
        doc["recoil_idd_batches"] = rec_b.tolist()
        doc["recoil_total_mev_per_primary"] = float((rec_b.sum(axis=1) * bin_g).mean())
    if lateral:
        r80_off = float(
            c.idd_mean(energy, "emonly")[2]["r80_mm"]
        )  # slab from the archived V11 EM-only R80
        z0, dz = R1D2_SLAB[0] * r80_off, (R1D2_SLAB[1] - R1D2_SLAB[0]) * r80_off
        nl = int(round(2 * R1D2_HALF_MM / R1D2_BIN_MM))
        lgrid = ScoringGrid(
            (-R1D2_HALF_MM, -R1D2_HALF_MM, z0),
            (R1D2_BIN_MM, R1D2_BIN_MM, dz),
            (nl, nl, 1),
            name="dose",
        )
        lcfg = build_config(energy, mode, n, seed, tid, batches, grid_geo=(geo, lgrid), tallies=())
        t1 = time.perf_counter()
        lres = v5b.Simulation(lcfg).run()
        doc["lateral"] = {"slab_mm": [z0, z0 + dz], "r80_ref_mm": r80_off, "wall_s": time.perf_counter() - t1,
                          "rms_mm_batches": v5d.lateral_rms_batches(lres.grid("dose").batch_energy_mev, R1D2_HALF_MM, R1D2_BIN_MM).tolist(),
                          "relative_residual": float(lres.energy_balance.relative_residual),
                          "valid": bool(v5.clean(lres) and lres.energy_balance.relative_residual <= 1e-12)}  # fmt: skip
    out.mkdir(parents=True, exist_ok=True)
    (out / f"de3-{energy}-{variant}.json").write_text(json.dumps(doc, sort_keys=True) + "\n")
    return doc


def _array_digest(x: Any) -> str:
    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(x, dtype=np.float64)).tobytes()
    ).hexdigest()


def render_check_base(doc: dict[str, Any]) -> str:
    d = doc["detail"]
    lines = ["", f"Energy {doc['energy']} MeV, n = {doc['n']}, {doc['n_batches']} batches, seed {doc['seed']}; pinned table {doc['pinned_id'][:8]}, base copy {doc['all_id'][:8]}.", "",
             f"**Bit for bit: {doc['bit_for_bit']}** (depth-dose batches {d['idd_batches']}, recoil channel {d['recoil_idd_batches']}, counters {d['counters']}).", "",
             "| array | sha256 pinned | sha256 base copy |", "|---|---|---|"]  # fmt: skip
    for k in ("sha256_idd_batches", "sha256_recoil_idd_batches"):
        lines.append(f"| {k[7:]} | {doc[k]['pinned']} | {doc[k]['all']} |")
    return "\n".join(lines) + "\n"


def check_base(
    energy: int,
    n: int,
    batches: int,
    out: Path,
    reuse_runs: bool = False,
    logs: tuple[str, ...] = (),
) -> dict[str, Any]:
    """Pre-D-E3 check (plan 3.1): the base copy ``all`` equals the pinned table bit for bit at the same seed.

    ``reuse_runs`` compares the run files already in ``out`` (no transport) and records the array digests.
    """
    if reuse_runs:
        a = json.loads((out / f"de3-{energy}-pinned.json").read_text())
        b = json.loads((out / f"de3-{energy}-all.json").read_text())
    else:
        a = run_one("pinned", energy, n, batches, False, out)
        b = run_one("all", energy, n, batches, False, out)
    same = {
        k: bool(np.array_equal(np.asarray(a[k]), np.asarray(b[k])))
        for k in ("idd_batches", "recoil_idd_batches")
    }
    same["counters"] = a["counters"] == b["counters"]
    res = {
        "energy": energy,
        "n": n,
        "n_batches": a["n_batches"],
        "seed": a["seed"],
        "bit_for_bit": bool(all(same.values())),
        "detail": same,
        "pinned_id": a["elastic_table_id"],
        "all_id": b["elastic_table_id"],
        "sha256_idd_batches": {"pinned": _array_digest(a["idd_batches"]), "all": _array_digest(b["idd_batches"])},
        "sha256_recoil_idd_batches": {"pinned": _array_digest(a["recoil_idd_batches"]),
                                      "all": _array_digest(b["recoil_idd_batches"])},
        "counters": a["counters"],
        "provenance": provenance(logs),
    }  # fmt: skip
    inputs = [out / f"de3-{energy}-pinned.json", out / f"de3-{energy}-all.json"]
    c.write_outputs(out, f"check-base-{energy}.json", res, f"check-base-{energy}-report.md",
                    c.header(f"D-E3 base check at {energy} MeV: pinned table against the base copy",
                             f"python -m validation.scripts.elastic_attribution.de3_probes --check-base --reuse-runs --energy {energy}",
                             inputs), render_check_base, GIT_ALLOW)  # fmt: skip
    return res


# -- provenance ----------------------------------------------------------------------------------------
SCRIPTS = ("make_variant_tables.py", "de3_probes.py", "de0_kernels.py", "common.py")


LOG_ROOT_DEFAULT = "~/.local/share/ionmc-experiment/experiment-v3/telemetry"
SCRIPT_REL = "validation/scripts/elastic_attribution/de3_probes.py"


def read_logs(names: list[str], root: Path | None = None) -> dict[str, Any]:
    """Provenance of logged runs from their ``meta.json`` (``<root>/test-logs/<name>/meta.json``): per run the label (variant and energy parsed
    from ``de3-`` / ``de6b-`` / ``de6c-<variant>-<E>``), commit, git status (clean / dirty with the lines), start / end times, the
    ``scripts_sha256`` of the logged runner if present (the strongest binding of the executed bytes) and the binding of the executed
    ``de3_probes.py``: **clean** commit -> the executed bytes ARE the blob at the recorded commit (``git rev-parse <commit>:<path>``, "bound: clean
    commit"); **dirty** working tree -> "unbound (dirty working tree)", no blob is derived. Groups: runs sharing commit and status."""
    import re
    import subprocess

    base = (
        Path(
            root or os.path.expanduser(os.environ.get("IONMC_TELEMETRY_DIR", LOG_ROOT_DEFAULT))
        ).expanduser()
        / "test-logs"
    )
    runs = []
    for n in names:
        meta = json.loads((base / n / "meta.json").read_text())
        m = re.fullmatch(r"(?:de3|de6b|de6c)-(.+)-(\d+)", meta["label"])
        status = [x for x in (meta.get("git_status") or "").splitlines() if x.strip()]
        sha = meta.get("scripts_sha256")
        blob = None
        if not status:
            out = subprocess.run(
                ["git", "rev-parse", f"{meta['commit']}:{SCRIPT_REL}"],
                cwd=c.REPO,
                capture_output=True,
                text=True,
            )
            blob = out.stdout.strip() if out.returncode == 0 else None
        binding = ("bound: scripts_sha256 recorded by the runner" if sha else
                   "bound: clean commit (executed bytes = blob at the recorded commit)" if (not status and blob) else
                   "unbound (dirty working tree)" if status else "unbound (blob not found)")  # fmt: skip
        runs.append({"log": n, "label": meta["label"], "variant": m.group(1) if m else None, "energy": int(m.group(2)) if m else None,
                     "commit": meta["commit"], "status_clean": not status, "status_lines": status, "script_blob": blob,
                     "scripts_sha256": sha, "binding": binding, "start_utc": meta.get("start_utc"), "end_utc": meta.get("end_utc"),
                     "exit_code": meta.get("exit_code")})  # fmt: skip
    groups: dict[tuple[str, bool, tuple[str, ...], str], list[dict[str, Any]]] = {}
    for r in runs:
        groups.setdefault(
            (r["commit"], r["status_clean"], tuple(r["status_lines"]), r["binding"]), []
        ).append(r)
    grp = [{"commit": k[0], "status_clean": k[1], "status_lines": list(k[2]), "binding": k[3], "runs": len(v),
            "first_start_utc": min(x["start_utc"] for x in v), "last_end_utc": max(x["end_utc"] for x in v),
            "script_blob": v[0]["script_blob"]} for k, v in sorted(groups.items())]  # fmt: skip
    return {"run_logs": runs, "run_groups": grp}


def provenance(logs: list[str] | None = None, logs_root: Path | None = None) -> dict[str, Any]:
    """Code state of the record: commit, dirty/untracked state, script and table hashes, backend."""
    import subprocess

    here = Path(__file__).resolve().parent
    git = lambda *a: subprocess.run(
        ["git", *a], cwd=c.REPO, capture_output=True, text=True
    ).stdout.strip()  # noqa: E731
    status = git(
        "status", "--short", "--untracked-files=all", "--", "validation/scripts/elastic_attribution"
    )
    return {
        "git_head": git("rev-parse", "HEAD"),
        "git_status_elastic_attribution": status.splitlines(),
        "git_status_clean": not status.strip(),
        "scripts_sha256": {k: c.sha256_file(here / k) for k in SCRIPTS},
        "table_files_sha256": {f.name: c.sha256_file(f) for f in c.table_files()},
        "backend": "warp-cpu", "precision": "float64", "single_process": True,
        "cache_dir": "<scratch-cache>", "seeds": SEEDS, **(read_logs(list(logs), logs_root) if logs else {}),
    }  # fmt: skip


BALANCE_TOL = 1.0e-12


def run_valid(d: dict[str, Any]) -> bool:
    """Fail-closed validity of one run file: the run is clean and its energy balance closes (``relative_residual`` is required,
    an absent value is invalid), and, if a lateral run is present, the same holds for it."""
    ok = bool(d.get("valid")) and d.get("relative_residual", float("inf")) <= BALANCE_TOL
    lat = d.get("lateral")
    if lat is not None:
        ok = (
            ok
            and bool(lat.get("valid"))
            and lat.get("relative_residual", float("inf")) <= BALANCE_TOL
        )
    return bool(ok)


def _jackknife(stat: Any, arrays: list[np.ndarray], den: np.ndarray) -> tuple[float, float]:
    """Delete-one-batch jackknife of ``stat(list of pooled numerator curves, pooled denominator curve)``; the same batch is
    removed from every array (common random numbers: the batches are shared across variants)."""
    nb = den.shape[0]
    full = stat([a.mean(axis=0) for a in arrays], den.mean(axis=0))
    tot, td = [a.sum(axis=0) for a in arrays], den.sum(axis=0)
    loo = np.array(
        [
            stat(
                [(t - a[i]) / (nb - 1) for t, a in zip(tot, arrays, strict=True)],
                (td - den[i]) / (nb - 1),
            )
            for i in range(nb)
        ]
    )
    return float(full), float(np.sqrt((nb - 1) / nb * ((loo - loo.mean()) ** 2).sum()))


def contrasts(
    cmp: Any, num_all: np.ndarray, num_pa: np.ndarray, num_pp: np.ndarray, den: np.ndarray
) -> dict[str, Any]:
    """Paired-jackknife value and SE of ``F_all - F_pA F_pp`` and of ``removal_pA + removal_pp - removal_all``."""
    pop = lambda x: cmp.curve_metrics(x)["peak_over_plateau"]  # noqa: E731
    f = lambda n, d: pop(n) / pop(d)  # noqa: E731
    rem = lambda n, d: 1.0 - n.max() / d.max()  # noqa: E731
    f_c, f_se = _jackknife(
        lambda n, d: f(n[0], d) - f(n[1], d) * f(n[2], d), [num_all, num_pa, num_pp], den
    )
    r_c, r_se = _jackknife(
        lambda n, d: rem(n[1], d) + rem(n[2], d) - rem(n[0], d), [num_all, num_pa, num_pp], den
    )
    return {"F_all_minus_product": f_c, "F_all_minus_product_se": f_se, "F_all_minus_product_z": f_c / f_se if f_se else None,
            "removal_additivity": r_c, "removal_additivity_se": r_se, "removal_additivity_z": r_c / r_se if r_se else None,
            "consistent_with_multiplicative": bool(f_se > 0 and abs(f_c) < 2.0 * f_se),
            "removal_consistent_with_additive": bool(r_se > 0 and abs(r_c) < 2.0 * r_se)}  # fmt: skip


def composition_and_gap(
    rows: dict[str, Any], out: Path, curves: dict[str, np.ndarray] | None = None, cmp: Any = None
) -> dict[str, Any]:
    """Derived rows requested by the D-E3 review: multiplicative composition F_pA x F_pp vs F_all, and the
    fraction of the archived V11 gap (F_ionmc / F_TOPAS - 1, from the H5 record) closed by each p+A shape variant.
    F_TOPAS is reconstructed as F_all / (1 + Delta_pk) with this run's F_all (equal to the archived F within 0.6 sigma)."""
    h5 = out.parent / "de0h5" / "h5.json"
    delta = {}
    if h5.is_file():
        delta = {
            int(k): v["delta_pk_F_ionmc_over_F_topas_minus_1"]
            for k, v in json.loads(h5.read_text())["per_energy"].items()
        }
    res: dict[str, Any] = {}
    for e in (150, 200):
        if f"{e}/all" not in rows:
            continue
        f_all, f_pa, f_pp = (
            rows.get(f"{e}/{v}", {}).get("F") for v in ("all", "pA-only", "pp-only")
        )
        block: dict[str, Any] = {}
        if f_pa is not None and f_pp is not None:
            block["F_pA_times_F_pp"] = f_pa * f_pp
            block["F_all"] = f_all
            block["F_all_minus_product"] = f_all - f_pa * f_pp
            block["removal_pA_plus_removal_pp_minus_removal_all"] = (
                rows[f"{e}/pA-only"]["peak_removal"]
                + rows[f"{e}/pp-only"]["peak_removal"]
                - rows[f"{e}/all"]["peak_removal"]
            )
            if (
                curves is not None
                and cmp is not None
                and all(f"{e}/{v}" in curves for v in ("all", "pA-only", "pp-only", "emonly"))
            ):
                block["paired_jackknife"] = contrasts(
                    cmp,
                    curves[f"{e}/all"],
                    curves[f"{e}/pA-only"],
                    curves[f"{e}/pp-only"],
                    curves[f"{e}/emonly"],
                )
        if e in delta:
            f_topas = f_all / (1.0 + delta[e])
            block["F_topas_reconstructed"] = f_topas
            block["delta_pk_archived"] = delta[e]
            block["gap_closed_fraction"] = {
                v: (rows[f"{e}/{v}"]["F"] - f_all) / (f_topas - f_all)
                for v in ("pA-sigma0.8", "pA-Rx1.2", "pA-s0.5", "pA-s0.9", "pA-gauss")
                if f"{e}/{v}" in rows
            }
            block["pp_only_gain_over_topas_total_gain"] = (
                None  # filled by the reviewer from the archived V11 gain; not computed here
            )
        res[str(e)] = block
    return res


# -- K* prediction and summary --------------------------------------------------------------------------
PRED_KIND = {"all": "table", "pp-only": "table", "pA-only": "table", "pA-sigma0.8": "table", "pA-Rx1.2": "Rx1.2",
             "pA-s0.5": "s0.5", "pA-s0.9": "s0.9", "pA-gauss": "gauss", "pinned": "table"}  # fmt: skip
PRED_WEIGHT = {
    "pp-only": (0.0, 1.0),
    "pA-only": (1.0, 0.0),
    "pA-sigma0.8": (0.8, 1.0),
}  # (O, H-1) cross-section factors


def kstar_prediction(variant: str, energy: int) -> dict[str, float]:
    """Peak removal predicted by D-E0 ``kernels``: ``K*`` shifted (physical) and as frozen (sign error; no conclusion),
    O + p-p, with the variant's cross-section factors (K* is linear in sigma)."""
    from . import de0_kernels as de0
    from .de1_shapes import Shapes

    tab = c.load_table()
    wo, wh = PRED_WEIGHT.get(variant, (1.0, 1.0))
    k = de0.kernels(Shapes(tab), c.Water(tab), energy, PRED_KIND[variant], ("O-16", "H-1"), True)
    return {"K_star_shifted": wo * k["O-16"]["K_star_shifted"] + wh * k["H-1"]["K_star_shifted"],
            "K_star_as_frozen_sign_error": wo * k["O-16"]["K_star_literal"] + wh * k["H-1"]["K_star_literal"]}  # fmt: skip


def _composition_sentence(derived: dict[str, Any]) -> str:
    parts = []
    for e, b in derived.items():
        pj = b.get("paired_jackknife")
        if pj:
            parts.append(
                f"{e} MeV: F_all - F_pA F_pp = {pj['F_all_minus_product']:+.5f} +- {pj['F_all_minus_product_se']:.5f} "
                f"({'consistent' if pj['consistent_with_multiplicative'] else 'NOT consistent'} with multiplicative within the paired uncertainty); "
                f"removal additivity {pj['removal_additivity']:+.5f} +- {pj['removal_additivity_se']:.5f} "
                f"({'consistent' if pj['removal_consistent_with_additive'] else 'NOT consistent'} with additive)"
            )
    return (
        "; ".join(parts) + "."
        if parts
        else "not evaluated (the all, pA-only, pp-only and EM-only run files are required)."
    )


def _run_provenance_lines(pv: dict[str, Any]) -> list[str]:
    """Derived from the recorded log metadata: which commit and git status produced which group of runs, and how the executed bytes are bound."""
    groups = pv.get("run_groups")
    if not groups:
        return [
            "- the logged run metadata (commit, git status) of the runs is not recorded in this record"
        ]
    out = []
    for g in groups:
        status = (
            "clean" if g.get("status_clean") else "dirty: " + "; ".join(g.get("status_lines", []))
        )
        blob = (f"; de3_probes.py blob {g['script_blob']} (`git rev-parse <commit>:validation/scripts/elastic_attribution/de3_probes.py`)"
                if g.get("script_blob") else "")  # fmt: skip
        out.append(
            f"- {g['runs']} logged run(s) at commit {g['commit']} with git status {status} (first start {g['first_start_utc']}, "
            f"last end {g['last_end_utc']}); binding of the executed bytes: {g.get('binding')}{blob}"
        )
    out.append(
        "- runs: " + ", ".join(f"{r['variant']}/{r['energy']}" for r in pv.get("run_logs", []))
    )
    return out


def render_summary(doc: dict[str, Any]) -> str:
    """Markdown body of ``de3-summary.json``, a pure function of the (projected) record."""
    rows, derived, pv = doc.get("rows", {}), doc.get("derived", {}), doc.get("provenance", {})
    lines = [
        "",
        "| E | variant | F | gain | peak removal | K* shifted (pred) | pred/meas removal | F pred (meas. gain) | recoil MeV | R80 mm | R80 EM-only | lateral rms mm (EM+el / EM) | valid |",
        "|" + "---|" * 13,
    ]
    for k, r in rows.items():
        e, v = k.split("/")
        lat = (
            f"{r['lateral_rms_mm']:.3f} / {r['lateral_rms_emonly_mm']:.3f}"
            if "lateral_rms_mm" in r
            else "-"
        )
        lines.append(
            f"| {e} | {v} | {r['F']:.4f} +- {r['F_se']:.4f} | {r['gain']:.4f} +- {r['gain_se']:.4f} | {r['peak_removal']:.4f} | "
            f"{r.get('K_star_shifted', float('nan')):.4f} | {r.get('peak_removal_pred_over_measured') or float('nan'):.3f} | "
            f"{r.get('F_pred_shifted_with_measured_gain', float('nan')):.4f} | {r['recoil_mev_per_primary'] or 0:.4f} | {r['r80_mm']:.3f} | {r['r80_emonly_mm']:.3f} | {lat} | {r['valid']} |"
        )
    lines += ["", "## Derived (D-E3 review): composition and gap closure", "",
              "| E | F_all | F_pA x F_pp | F_all - F_pA F_pp (paired jackknife) | removal additivity pA + pp - all (paired jackknife) | F_TOPAS (reconstructed) | gap closed: sigma0.8 / Rx1.2 / s0.5 / s0.9 / gauss |",
              "|---|---|---|---|---|---|---|"]  # fmt: skip
    for e, b in derived.items():
        gc = b.get("gap_closed_fraction", {})
        pj = b.get("paired_jackknife")
        dif = (
            f"{pj['F_all_minus_product']:+.5f} +- {pj['F_all_minus_product_se']:.5f}"
            if pj
            else f"{b.get('F_all_minus_product', float('nan')):+.5f}"
        )
        add = (
            f"{pj['removal_additivity']:+.5f} +- {pj['removal_additivity_se']:.5f}"
            if pj
            else f"{b.get('removal_pA_plus_removal_pp_minus_removal_all', float('nan')):+.5f}"
        )
        lines.append(
            f"| {e} | {b.get('F_all', float('nan')):.4f} | {b.get('F_pA_times_F_pp', float('nan')):.4f} | {dif} | "
            f"{add} | {b.get('F_topas_reconstructed', float('nan')):.4f} | "
            + " / ".join(
                f"{gc.get(v, float('nan')):.2f}"
                for v in ("pA-sigma0.8", "pA-Rx1.2", "pA-s0.5", "pA-s0.9", "pA-gauss")
            )
            + " |"
        )
    lines += ["", "## Findings recorded with the independent D-E3 review (exploratory; nothing enters the section 5 rule)", "",
              "- **Composition:** " + _composition_sentence(derived) + " Where it is consistent, D-E6 effects are to be read on -ln F; on F itself the product term would consume part of the 30 % interaction budget of section 5.",
              "- **Gap closure toward TOPAS F** (F_TOPAS reconstructed from the archived V11 Delta_pk): pA-s0.9 closes most of the gap, pA-s0.5 and pA-Rx1.2 about half; the shape (core width) is the lever on F.",
              "- **Contrary to a pure tail story:** pA-gauss removes the large-angle tail yet leaves F unchanged while lowering the gain; the peak deficit is driven by the core width, the plateau gain by the tail.",
              "- **Contrary evidence on p-p (H3 stays live):** the pp-only plateau gain alone is of the order of TOPAS's total elastic gain (84 % / 108 % per the review, from the archived V11 gain); no p+A shape variant closes the gain excess, so the swap-pp arm of D-E6 is essential.",
              "- **K* model:** the shifted K* overpredicts the measured peak removal by 8-14 % for every p+A variant and by 1-3 % for pp-only, the bias anticipated in the D-E0 dated note (adequate for ranking, not for a 0.010 closure).",
              "", "## Provenance", "",
              f"- recorded git HEAD {pv.get('git_head')}; recorded status of validation/scripts/elastic_attribution: "
              + ("clean" if pv.get("git_status_clean") else "dirty: " + "; ".join(pv.get("git_status_elastic_attribution") or ["(not recorded)"])),
              "- script sha256 recorded at summary time: " + ", ".join(f"{k} {v[:16]}" for k, v in pv.get("scripts_sha256", {}).items()),
              "- script sha256 recorded for the run (staged copies, if any): " + (", ".join(f"{k} {v[:16]}" for k, v in pv.get("runtime_scripts_sha256", {}).items()) or "none recorded"),
              "- table files sha256: " + ", ".join(f"{k} {v[:16]}" for k, v in pv.get("table_files_sha256", {}).items()),
              f"- backend {pv.get('backend')} {pv.get('precision')}, single process {pv.get('single_process')}; seeds {pv.get('seeds')}; runs (seed, n, batches, table id, wall s) in `runs` of de3-summary.json",
              *_run_provenance_lines(pv)]  # fmt: skip
    return "\n".join(lines) + "\n"


def summarize(out: Path, logs: list[str] | None = None) -> dict[str, Any]:
    cmp = _cmp()
    rows: dict[str, Any] = {}
    by_key: dict[str, dict[str, Any]] = {}
    for e in (150, 200):
        ref_p = out / f"de3-{e}-emonly.json"
        if not ref_p.is_file():
            continue
        ref = json.loads(ref_p.read_text())
        by_key[f"{e}/emonly"] = ref
        den = np.asarray(ref["idd_batches"])
        for p in sorted(out.glob(f"de3-{e}-*.json")):
            d = json.loads(p.read_text())
            if d["variant"] == "emonly":
                continue
            by_key[f"{e}/{d['variant']}"] = d
            num = np.asarray(d["idd_batches"])
            jk = cmp.paired_ratio_jackknife(num, den)
            f, g = jk["pp_ratio"], jk["plateau_ratio"]
            peak_removal = 1.0 - (num.mean(axis=0).max() / den.mean(axis=0).max())
            row: dict[str, Any] = {
                "F": f.value, "F_se": f.se, "gain": g.value - 1.0, "gain_se": g.se, "peak_removal": float(peak_removal),
                "r80_mm": d["metrics"]["r80_mm"], "r80_emonly_mm": ref["metrics"]["r80_mm"],
                "recoil_mev_per_primary": d.get("recoil_total_mev_per_primary"), "valid": run_valid(d) and run_valid(ref),
                "n": d["n"], "wall_s": d["wall_s"],
            }  # fmt: skip
            if "lateral" in d and "lateral" in ref:
                r_on, r_off = (
                    np.asarray(d["lateral"]["rms_mm_batches"]),
                    np.asarray(ref["lateral"]["rms_mm_batches"]),
                )
                row["lateral_rms_mm"], row["lateral_rms_emonly_mm"] = (
                    float(r_on.mean()),
                    float(r_off.mean()),
                )
                row["lateral_rms_diff_se"] = float((r_on - r_off).std(ddof=1) / np.sqrt(r_on.size))
            if d["variant"] in PRED_KIND:
                kp = kstar_prediction(d["variant"], e)
                row.update(
                    kp,
                    F_pred_shifted_with_measured_gain=(1.0 - kp["K_star_shifted"]) / g.value,
                    peak_removal_pred_over_measured=kp["K_star_shifted"] / peak_removal
                    if peak_removal
                    else None,
                )
            rows[f"{e}/{d['variant']}"] = row
    runs_prov = {k: {"seed": d["seed"], "n": d["n"], "n_batches": d["n_batches"], "elastic_table_id": d["elastic_table_id"],
                     "mode": d["mode"], "wall_s": d["wall_s"], "lateral_wall_s": (d.get("lateral") or {}).get("wall_s")}
                 for k, d in sorted(by_key.items())}  # fmt: skip
    curves = {k: np.asarray(d["idd_batches"]) for k, d in by_key.items()}
    derived = composition_and_gap(rows, out, curves, cmp)
    doc = {"banner": c.BANNER, "rows": rows, "runs": runs_prov, "derived": derived, "provenance": provenance(logs),
           "note": "F = (peak/plateau)_EM+el/(peak/plateau)_EM by the paired jackknife of compare_idd_v5; K* shifted is the physical "
                   "kernel, the as-frozen K* carries a sign error (D-E0 dated note); probes only, never used in the section 5 rule"}  # fmt: skip
    run_files = sorted(p for e in (150, 200) for p in out.glob(f"de3-{e}-*.json"))
    c.write_outputs(out, "de3-summary.json", doc, "de3-summary.md",
                    c.header("D-E3 probes: summary", "python -m validation.scripts.elastic_attribution.de3_probes --summary",
                             run_files), render_summary, GIT_ALLOW)  # fmt: skip
    return doc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument(
        "--variant",
        default="all",
        help="all, pp-only, pA-only, pA-sigma0.8, pA-Rx1.2, pA-s0.5, pA-s0.9, pA-gauss, emonly, pinned",
    )
    ap.add_argument("--energy", type=int, choices=(150, 200), default=150)
    ap.add_argument("--histories", type=int, default=100_000)
    ap.add_argument("--batches", type=int, default=20)
    ap.add_argument(
        "--no-lateral",
        action="store_true",
        help="skip the lateral-rms run (a second transport run of the same size)",
    )
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument(
        "--check-base", action="store_true", help="pinned vs base copy 'all' bit for bit (two runs)"
    )
    ap.add_argument(
        "--logs",
        nargs="+",
        default=None,
        help="--summary: telemetry log directory names of the runs (recorded in the provenance)",
    )
    ap.add_argument("--summary", action="store_true", help="combine the run files of --out")
    ap.add_argument(
        "--scratch-cache",
        type=Path,
        default=None,
        help="scratch cache of make_variant_tables (env IONMC_DE3_SCRATCH; default validation/generated/elastic_attribution/de3-cache)",
    )
    ap.add_argument(
        "--reuse-runs",
        action="store_true",
        help="--check-base from the existing run files (no transport)",
    )
    a = ap.parse_args(argv)
    if a.histories % a.batches:
        raise SystemExit("--histories must be a multiple of --batches")
    from .make_variant_tables import resolve_scratch

    scratch = resolve_scratch(a.scratch_cache)
    if a.summary:
        print(json.dumps(summarize(a.out, a.logs)["rows"], indent=1))
    elif a.check_base:
        with c.scoped_cache_env(scratch):
            print(
                json.dumps(
                    check_base(
                        a.energy,
                        a.histories,
                        a.batches,
                        a.out,
                        reuse_runs=a.reuse_runs,
                        logs=tuple(a.logs or ()),
                    ),
                    indent=1,
                )
            )
    else:
        with c.scoped_cache_env(scratch):
            d = run_one(a.variant, a.energy, a.histories, a.batches, not a.no_lateral, a.out)
        print(
            json.dumps(
                {
                    k: d[k]
                    for k in (
                        "variant",
                        "energy",
                        "seed",
                        "n",
                        "metrics",
                        "valid",
                        "wall_s",
                        "hist_per_s",
                    )
                },
                indent=1,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
