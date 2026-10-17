"""D-E6 swap tables and summary (exploratory, report-only; plan row D-E6, Amendment 19 (e), addendum (b)4; plan section 5).

``--build``: writes the variants ``swap-pA`` (the O-16 mu_CM inverse CDF replaced by the Geant4 empirical CDF of the water rows 14-18),
``swap-pp`` (the H-1 CDF above 16.26 deg CM replaced by the Geant4 p-p CDF of rows 14-18), ``swap-both``, and the post-hoc arms of
addendum (c)3 ``swap-pp-sigma`` (ionmc shapes; sigma_pp scaled per energy by the Geant4/table p-p rate ratio after the corrected rule,
linear in E between the D-E4 energies, factor 1 below 30 MeV) and ``swap-all`` (both shapes + scaled sigma_pp) into the D-E3 scratch cache
(``make_variant_tables``; run that first; it rebuilds the cache and therefore removes the swap tables) and records them in its
``variants.json``, so that ``de3_probes --variant swap-pA --energy 150`` runs them with the D-E3 seeds 20990100 / 20990101 (CRN with
the base, ``all``, and EM-only runs). sigma is unchanged. The sources are the D-E4 CDF files ``de4-case0..4-cdf.npz`` (never
committed; sidecar sha256 recorded in the ``exploratory`` key).

Construction of the swapped nodes (declared): for a table node of energy E with 30 MeV <= E <= E_last (E_last = the first grid node
>= 200 MeV; above 200 MeV the 200 MeV q-edges are held) the 257 edges of the 256 equiprobable mu bins are obtained at fixed CDF level j
from ``q_ij = 2 p_CM sin(theta_CM / 2)`` of the D-E4 energies E_i, **linear in E** (``np.interp``), converted with the node's own
CM momentum, ``mu = 1 - 2 (q / 2 p_CM)^2``; the end edges are the sample minimum and maximum of mu_CM (addendum (c)4, interpolated like the
others) and the edges are made strictly ascending. sigma is unchanged in the shape arms. Below 30 MeV the ionmc shape is kept (the transport interpolates the edges linearly in E between the last
ionmc node and the first swapped node).

Finite-sample CDF error (plan 3.2, addendum (c)4): ``finite_sample_error`` = quadrature of (a) a 200-resample bootstrap of the empirical CDFs
(seed 20990200) propagated through the D-E0 kernels K* (shifted) and G at 150 and 200 MeV, (b) the tail-edge systematic (end edges
moved from the sample extremes to the middle of the extreme probability bins), (c) a leave-one-energy-out check at 150 MeV (the CDF predicted from
100 and 200 MeV by the q-scaled interpolation; used for both energies, flagged); expressed on the predicted F = (1 - K*)/(1 + G) and entering SE_r.

``--de6-summary``: from the run files of the D-E3 directory: Delta_pk(X), Delta_g(X) against ``all`` (paired batch jackknife, CRN) for
X in pA, pp, both; the interaction I = Delta(both) - Delta(pA) - Delta(pp); the closure residual r = Delta_obs - Delta(both) with
SE_r = sqrt(SE_obs^2 + SE_swap^2 + SE_cdf^2), on 1 - F (rule of record: Delta = F_X - F_base, positive closes) and on -ln F (addendum (b)4
sensitivity); the extended set of addendum (c)3 (Delta(all), main effects pA shape / pp shape / pp sigma, I_all, r_all); the H4 recoil-channel
share of the gain excess; clauses 1-3 of section 5 as computable now; the section 2.2 signature inputs (from the D-E4 ratios record) printed, not decided.

Run: ``python -m validation.scripts.elastic_attribution.de6_swap --build`` and ``... --de6-summary``.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from . import common as c
from . import make_variant_tables as mvt
from .de1_shapes import Shapes

SWAP_VARIANTS = ("swap-pA", "swap-pp", "swap-both", "swap-pp-sigma", "swap-all")
SWAP_CLASSES = {
    "swap-pA": ("O-16",),
    "swap-pp": ("H-1",),
    "swap-both": ("O-16", "H-1"),
    "swap-pp-sigma": (),
    "swap-all": ("O-16", "H-1"),
}
SIGMA_ARMS = ("swap-pp-sigma", "swap-all")  # sigma_pp scaled by the Geant4/table rate ratio
RATIOS_PATH = c.OUT / "de4" / "de4-ratios.json"
SWAP_E_MIN = 30.0
CDF_SEED = 20990200  # inside the investigation family 2099xxxx; distinct from 20990100/1 (D-E3/D-E6 runs) and 20991000+ (D-E5)
N_Q = 256
CDF_DIR = c.OUT / "de4" / "cdf"
OUT = c.OUT / "de6"
WATER_CASES = (0, 1, 2, 3, 4)
MU_RANGE = {"O-16": (-1.0, 1.0), "H-1": (-0.96, 0.96)}
ARCHIVED_SE_REL = {}  # filled from the V11 step-16 record at run time


@dataclass
class ClassSource:
    """Geant4 empirical mu_CM CDFs of one class at the D-E4 energies."""

    energies: np.ndarray  # (nE,)
    mu: list[np.ndarray]  # sorted samples
    p_cm: np.ndarray  # (nE,) CM momentum [MeV/c]
    lo: float  # physical range of mu (clip)
    hi: float
    edges_q: np.ndarray  # (nE, 257) q at the 257 equiprobable mu edges [MeV/c]
    symmetrised: bool = (
        False  # p-p: the stored sample holds the (+mu, -mu) pair of every physical event
    )


def _edges_q(mu: np.ndarray, p_cm: float, tail_variant: bool = False) -> np.ndarray:
    """q at the 257 equiprobable edges of the sample (end edges = sample min / max; ``tail_variant``: the end edges at the middle of the
    extreme probability bins, the tail-edge systematic of addendum (c)4)."""
    srt = np.sort(mu)
    e = np.quantile(srt, np.arange(N_Q + 1) / N_Q)
    if tail_variant:
        e[0], e[-1] = np.quantile(srt, 0.5 / N_Q), np.quantile(srt, 1.0 - 0.5 / N_Q)
    return 2.0 * p_cm * np.sin(np.arccos(np.clip(e, -1.0, 1.0)) / 2.0)


def load_sources(cdf_dir: Path = CDF_DIR) -> tuple[dict[str, ClassSource], list[dict[str, Any]]]:
    """Sources of the water rows 14-18 (cases 0-4): classes O-16 and H-1, and their provenance (sidecar / npz sha256)."""
    per: dict[str, dict[str, list[Any]]] = {"O-16": {}, "H-1": {}}
    prov = []
    for n in WATER_CASES:
        side_p = Path(cdf_dir) / f"de4-case{n}-cdf.json"
        side = json.loads(side_p.read_text())
        npz_p = Path(cdf_dir) / side["npz"]
        assert c.sha256_file(npz_p) == side["npz_sha256"], (
            f"{npz_p}: sha256 differs from its sidecar"
        )
        z = np.load(npz_p)
        assert side["material"] == "water" and side["case"] == n
        prov.append({"case": n, "row": side["row"], "energy_mev": side["energy_mev"], "primaries": side["primaries"],
                     "npz_sha256": side["npz_sha256"], "sidecar_sha256": c.sha256_file(side_p)})  # fmt: skip
        for cls, d in per.items():
            d.setdefault("E", []).append(float(side["energy_mev"]))
            d.setdefault("mu", []).append(np.sort(z[f"{cls}_mu_sorted"]))
            d.setdefault("p", []).append(float(side["classes"][cls]["p_cm_mev"]))
    out = {}
    for cls, d in per.items():
        order = np.argsort(d["E"])
        e = np.asarray(d["E"])[order]
        mus = [d["mu"][i] for i in order]
        p = np.asarray(d["p"])[order]
        lo, hi = MU_RANGE[cls]
        out[cls] = ClassSource(
            e, mus, p, lo, hi, np.array([_edges_q(m, pc) for m, pc in zip(mus, p, strict=True)])
        )
    return out, prov


def swap_edges(e_node: float, src: ClassSource, p_cm_node: float) -> np.ndarray:
    """The 257 mu edges of the swapped node (see the module docstring): linear in E at fixed CDF level in q."""
    q = np.array([np.interp(e_node, src.energies, src.edges_q[:, j]) for j in range(N_Q + 1)])
    s = np.clip(q / (2.0 * p_cm_node), 0.0, 1.0)
    mu = np.clip(1.0 - 2.0 * s * s, src.lo, src.hi)
    mu = np.maximum.accumulate(mu)
    for j in range(1, mu.size):  # strictly ascending
        if mu[j] <= mu[j - 1]:
            mu[j] = np.nextafter(mu[j - 1], np.inf)
    return mu


def sigma_pp_factor(e: float, energies: np.ndarray, factors: np.ndarray) -> float:
    """sigma_pp scale at ``e``: 1 below 30 MeV, linear in E between the D-E4 energies, held above the last."""
    return 1.0 if e < SWAP_E_MIN else float(np.interp(e, energies, factors))


def load_sigma_factors(ratios_path: Path = RATIOS_PATH) -> dict[str, list[float]]:
    """Geant4 / table p-p rate ratio (1 / ratio of the D-E4 ratios record, rate per primary after the rule) at the D-E4 energies."""
    doc = json.loads(Path(ratios_path).read_text())
    ratios = [
        doc["cases"][str(n)]["H-1"]["signatures"]["rate_per_primary_after_rule"]["ratio"]
        for n in WATER_CASES
    ]
    return {
        "energies_mev": [30.0, 50.0, 100.0, 150.0, 200.0],
        "table_over_geant4_rate": ratios,
        "factors": [1.0 / r for r in ratios],
    }


def swap_arrays(base: dict[str, np.ndarray], names: tuple[str, ...], sources: dict[str, ClassSource],
                classes: tuple[str, ...], sigma: dict[str, list[float]] | None = None) -> tuple[dict[str, np.ndarray], dict[str, Any]]:  # fmt: skip
    """Copy of the table arrays with the mu edges of ``classes`` swapped and, if ``sigma`` is given, sigma_pp (H-1) scaled by the
    per-energy factor (everything else unchanged)."""
    a = {k: v.copy() for k, v in base.items()}
    grid = np.asarray(a["grid_e_mev"], dtype=float)
    stats: dict[str, Any] = {}
    for cls in classes:
        it = names.index(cls)
        src = sources[cls]
        e_last = float(grid[np.nonzero(grid >= src.energies.max())[0][0]])
        done = 0
        for k, e in enumerate(grid):
            if not a["valid"][it, k] or e < SWAP_E_MIN or e > e_last:
                continue
            a["edges_mu"][it, k] = swap_edges(float(e), src, float(a["p_cm_mev"][it, k]))
            done += 1
        stats[cls] = {"nodes_swapped": done, "e_min_mev": SWAP_E_MIN, "e_last_mev": e_last}
    if sigma is not None:
        ih = names.index("H-1")
        en, fa = np.asarray(sigma["energies_mev"]), np.asarray(sigma["factors"])
        for k, e in enumerate(grid):
            a["sigma_barn"][ih, k] *= sigma_pp_factor(float(e), en, fa)
        stats["sigma_pp_scaled"] = True
    return a, stats


def build_swap_tables(
    scratch: Path = mvt.SCRATCH, cdf_dir: Path = CDF_DIR, ratios_path: Path = RATIOS_PATH
) -> dict[str, Any]:
    """Write the swap tables into the D-E3 scratch cache and register them in its ``variants.json``."""
    from ionmc.nuclear.build import table_id

    mvt.check_scratch(scratch)
    vj = Path(scratch) / "variants.json"
    if not vj.is_file():
        raise SystemExit(f"{vj} not found: run make_variant_tables first")
    doc = json.loads(vj.read_text())
    der = mvt._derived(scratch)
    for v in SWAP_VARIANTS:  # remove an earlier build of the swap tables
        old = doc["variants"].pop(v, None)
        if old:
            for ext in (".json", ".npz"):
                (der / f"elastic-proton-{old['table_id']}{ext}").unlink(missing_ok=True)
    src_cache = mvt.REAL_CACHE / "derived"
    pinned = json.loads((src_cache / f"elastic-proton-{c.TABLE_ID}.json").read_text())
    base = dict(np.load(src_cache / f"elastic-proton-{c.TABLE_ID}.npz", allow_pickle=False).items())
    names = tuple(pinned["target_names"])
    sources, prov = load_sources(cdf_dir)
    sigma = load_sigma_factors(ratios_path)
    for v in SWAP_VARIANTS:
        arrays, stats = swap_arrays(
            base, names, sources, SWAP_CLASSES[v], sigma if v in SIGMA_ARMS else None
        )
        mvt._assert_table(arrays, names, check_ends=not SWAP_CLASSES[v])
        tmp = der / "tmp-swap.npz"
        np.savez(tmp, **arrays)
        sha = c.sha256_file(tmp)
        side = json.loads(json.dumps(pinned))
        side["npz_sha256"] = sha
        side["exploratory"] = {
            "d_e6_variant": v, "source_table_id": c.TABLE_ID, "purpose": mvt.PURPOSE.replace("D-E3 sensitivity probe (Amendment 19 (b))", "D-E6 swap test (Amendment 19 (e))"),
            "parameters": {"classes": list(SWAP_CLASSES[v]), "interpolation": "linear in E at fixed CDF level of q = 2 p_CM sin(theta_CM/2)",
                           "below_mev_ionmc_shape": SWAP_E_MIN,
                           "sigma": sigma if v in SIGMA_ARMS else "ionmc (unchanged)", "end_edges": "sample min / max of mu_CM"},
            "normalisation": stats, "cdf_sources": prov,
        }  # fmt: skip
        side["table_id"] = table_id(side)
        tid = side["table_id"]
        tmp.rename(der / f"elastic-proton-{tid}.npz")
        (der / f"elastic-proton-{tid}.json").write_text(
            json.dumps(side, indent=1, sort_keys=True) + "\n"
        )
        doc["variants"][v] = {
            "variant": v,
            "table_id": tid,
            "npz_sha256": sha,
            "normalisation": stats,
        }
    vj.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    return {v: doc["variants"][v]["table_id"] for v in SWAP_VARIANTS}


# -- finite-sample CDF error ------------------------------------------------------------------------------
class SwapShapes(Shapes):
    """``Shapes`` whose O-16 / H-1 mu quantiles come from (bootstrapped) Geant4 sources, for the K* propagation."""

    def __init__(
        self, table: Any, sources: dict[str, ClassSource], classes: tuple[str, ...]
    ) -> None:
        super().__init__(table)
        self.sources, self.classes = sources, classes

    def mu_quantiles(self, kind: str, name: str, e_mev: float, u: np.ndarray) -> np.ndarray:
        if name in self.classes and e_mev >= SWAP_E_MIN:
            pcm = float(
                np.interp(e_mev, self.a["grid_e_mev"], self.a["p_cm_mev"][self.index(name)])
            )
            edges = swap_edges(e_mev, self.sources[name], pcm)
            x = np.asarray(u, dtype=float) * N_Q
            i = np.minimum(x.astype(np.int64), N_Q - 1)
            return np.asarray(edges[i] + (x - i) * (edges[i + 1] - edges[i]))
        return super().mu_quantiles(kind, name, e_mev, u)


def _resample(m: np.ndarray, symmetrised: bool, rng: np.random.Generator) -> np.ndarray:
    """One bootstrap replicate of a stored sample. p-p: the sample holds the (+mu, -mu) pair of each physical event, so the N physical
    events (the non-negative half) are resampled with replacement and each value is mirrored again (effective sample size N, not 2N)."""
    srt = np.sort(m)
    if symmetrised:
        half = srt[srt.size // 2 :]
        x = rng.choice(half, size=half.size, replace=True)
        return np.sort(np.concatenate([x, -x]))
    return np.sort(rng.choice(srt, size=srt.size, replace=True))


def _bootstrap(
    sources: dict[str, ClassSource], rng: np.random.Generator, classes: tuple[str, ...]
) -> dict[str, ClassSource]:
    out = dict(sources)
    for cls in classes:
        s = sources[cls]
        mus = [_resample(m, s.symmetrised, rng) for m in s.mu]
        out[cls] = ClassSource(s.energies, mus, s.p_cm, s.lo, s.hi, np.array([_edges_q(m, pc) for m, pc in zip(mus, s.p_cm, strict=True)]), s.symmetrised)  # fmt: skip
    return out


def _tail_variant(sources: dict[str, ClassSource]) -> dict[str, ClassSource]:
    out = {}
    for cls, s in sources.items():
        out[cls] = ClassSource(s.energies, s.mu, s.p_cm, s.lo, s.hi, np.array([_edges_q(m, pc, True) for m, pc in zip(s.mu, s.p_cm, strict=True)]), s.symmetrised)  # fmt: skip
    return out


def _drop_energy(sources: dict[str, ClassSource], e_drop: float) -> dict[str, ClassSource]:
    out = {}
    for cls, s in sources.items():
        keep = np.abs(s.energies - e_drop) > 1e-9
        out[cls] = ClassSource(
            s.energies[keep],
            [m for m, k in zip(s.mu, keep, strict=True) if k],
            s.p_cm[keep],
            s.lo,
            s.hi,
            s.edges_q[keep],
        )
    return out


def _kg(
    tab: Any, water: Any, sources: dict[str, ClassSource], e0: int, n_z: int, dz_mm: float
) -> tuple[float, float, float]:
    """(K* shifted [O + p-p], G total, predicted F = (1 - K*) / (1 + G)) of the swapped shape (both classes) at ``e0``."""
    from . import de0_kernels as de0

    shapes = SwapShapes(tab, sources, ("O-16", "H-1"))
    k = sum(
        v["K_star_shifted"]
        for v in de0.kernels(shapes, water, e0, "table", ("O-16", "H-1"), True, n_z=n_z).values()
    )
    g = de0.plateau_gain(shapes, water, e0, "table", dz_mm=dz_mm)["G_total"]
    return float(k), float(g), float((1.0 - k) / (1.0 + g))


def finite_sample_error(
    sources: dict[str, ClassSource],
    n_boot: int = 200,
    n_z: int = 100,
    dz_mm: float = 2.0,
    prov: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """The finite-sample CDF error of plan 3.2 / addendum (c)4, on the predicted F = (1 - K*)/(1 + G) of the swapped shape at 150 and 200 MeV:
    (a) sd over ``n_boot`` bootstrap CDFs (both classes resampled), (b) the tail-edge systematic, (c) the leave-one-energy-out check at 150 MeV
    (CDF predicted from 100 and 200 MeV); total = quadrature (the 150 MeV value of (c) is used at 200 MeV too, flagged)."""
    tab = c.load_table()
    from . import common as _c

    water = _c.Water(tab)
    res: dict[str, Any] = {"n_boot": n_boot, "seed": CDF_SEED, "n_z": n_z, "gain_dz_mm": dz_mm,
                           "events_per_energy": {cls: [int(m.size) for m in s.mu] for cls, s in sources.items()},
                           "source_energies_mev": [float(x) for x in sources["O-16"].energies],
                           "lineage": {"cdf_sources": [{"case": q["case"], "energy_mev": q["energy_mev"], "npz_sha256": q["npz_sha256"]} for q in (prov or [])],
                                       "kernels": "D-E0 kernels: K* shifted (O-16 + H-1) and G total, swapped shapes, 'table' kind",
                                       "bootstrap": "N physical events resampled with replacement per class and energy; p-p mirrored",
                                       "n_boot": n_boot, "seed": CDF_SEED, "n_z": n_z, "gain_dz_mm": dz_mm}}  # fmt: skip
    loo = None
    k0, g0, f0 = _kg(tab, water, sources, 150, n_z, dz_mm)
    k1, g1, f1 = _kg(tab, water, _drop_energy(sources, 150.0), 150, n_z, dz_mm)
    loo = {"dK_star": k1 - k0, "dG": g1 - g0, "dF": f1 - f0}
    for e0 in (150, 200):
        kn, gn, fn = _kg(tab, water, sources, e0, n_z, dz_mm)
        kt, gt, ft = _kg(tab, water, _tail_variant(sources), e0, n_z, dz_mm)
        rng = np.random.default_rng(np.random.SeedSequence(CDF_SEED + e0))
        ks, gs, fs = [], [], []
        for _ in range(n_boot):
            k, g, f = _kg(tab, water, _bootstrap(sources, rng, ("O-16", "H-1")), e0, n_z, dz_mm)
            ks.append(k), gs.append(g), fs.append(f)
        boot = {
            "K_star_sd": float(np.std(ks, ddof=1)),
            "G_sd": float(np.std(gs, ddof=1)),
            "F_sd": float(np.std(fs, ddof=1)),
        }
        tail = {"dK_star": kt - kn, "dG": gt - gn, "dF": ft - fn}
        total = math.sqrt(boot["F_sd"] ** 2 + tail["dF"] ** 2 + loo["dF"] ** 2)
        res[str(e0)] = {"nominal": {"K_star_shifted": kn, "G_total": gn, "F_pred": fn}, "bootstrap": boot, "tail_edge": tail,
                        "leave_one_energy_out_at_150": loo, "loo_is_proxy_at_200": e0 == 200, "F_total_error": total}  # fmt: skip
    return res


# -- summary ---------------------------------------------------------------------------------------------
def signature_inputs(
    ratios: dict[str, Any], effects: dict[str, dict[str, float]]
) -> dict[str, Any]:
    """The section 2.2 signature inputs of H1 / H2 / H3 from the D-E4 ratios record (``ratios``) and the swap closures (``effects[E][X]`` =
    Delta(X) / Delta_pk,obs). Printed as inputs; no decision is made here (clauses 4 and 5 of section 5)."""
    water = {30.0: "0", 50.0: "1", 100.0: "2", 150.0: "3", 200.0: "4"}
    graphite = {30.0: "5", 50.0: "6", 100.0: "7", 150.0: "8", 200.0: "9"}

    def series(cases: dict[float, str], cls: str, sig: str) -> list[float | None]:
        return [ratios["cases"][n][cls]["signatures"][sig]["ratio"] for n in cases.values()]

    def outside(xs: list[float | None], tol: float) -> list[bool]:
        return [x is not None and abs(x - 1.0) > tol for x in xs]

    rate_o, rate_c = (
        series(water, "O-16", "rate_per_primary_after_rule"),
        series(graphite, "C-12", "rate_per_primary_after_rule"),
    )
    p2, p3 = (
        series(water, "O-16", "P_theta_lab_gt_theta_2mm"),
        series(water, "O-16", "P_theta_lab_gt_theta_3mm"),
    )
    med, tail = (
        series(water, "O-16", "theta_lab_median_deg"),
        series(water, "O-16", "T_tail_moment"),
    )
    pp_rate, pp_t = (
        series(water, "H-1", "rate_per_primary_after_rule"),
        series(water, "H-1", "mean_t_other_mev"),
    )
    pp_exc = series(water, "H-1", "P_theta_lab_gt_theta_2mm")
    n_ge = lambda xs, thr: sum(1 for x in xs if x is not None and x >= thr)  # noqa: E731
    return {
        "energies_mev": [30, 50, 100, 150, 200],
        "H1": {
            "rate_ratio_O": rate_o,
            "rate_ratio_C": rate_c,
            "n_energies_O_outside_10pct": sum(outside(rate_o, 0.10)),
            "all_within_5pct_O_C": not any(outside(rate_o, 0.05) + outside(rate_c, 0.05)),
            "note": "swap-sigma is conditional on the identified signature; not decided here",
        },  # fmt: skip
        "H2": {
            "exceedance_ratio_2mm": p2,
            "exceedance_ratio_3mm": p3,
            "theta_med_ratio": med,
            "T_ratio": tail,
            "n_energies_exceedance_ge_1p15": max(n_ge(p2, 1.15), n_ge(p3, 1.15)),
            "T_ratio_ge_1p5_energies": n_ge(tail, 1.5),
            "all_within_10pct": not any(
                outside(p2, 0.10) + outside(p3, 0.10) + outside(med, 0.10) + outside(tail, 0.10)
            ),
            "swap_pA_closure_fraction": {e: v.get("pA") for e, v in effects.items()},
        },  # fmt: skip
        "H3": {
            "pp_rate_ratio": pp_rate,
            "pp_T_slow_ratio": pp_t,
            "pp_exceedance_ratio_2mm": pp_exc,
            "any_outside_10pct": any(outside(pp_rate, 0.10) + outside(pp_t, 0.10)),
            "all_within_5pct": not any(outside(pp_rate, 0.05) + outside(pp_t, 0.05)),
            "swap_pp_closure_fraction": {e: v.get("pp") for e, v in effects.items()},
        },  # fmt: skip
        "decision": "inputs only; clauses 4 and 5 of section 5 are decided by the synthesis, not here",
    }


def _clause3(d_pa: float, d_pp: float, d_both: float, inter: float) -> dict[str, Any]:
    fa, fp = d_pa / d_both, d_pp / d_both
    named = "pA" if fa >= 0.70 else "pp" if fp >= 0.70 else None
    return {"main_effect_pA_fraction": fa, "main_effect_pp_fraction": fp, "named_component": named,
            "interaction_fraction": abs(inter) / abs(d_both), "holds": bool(named is not None and abs(inter) <= 0.30 * abs(d_both))}  # fmt: skip


def arithmetic(delta: dict[str, float], se_both: float, se_i: float, obs: float, se_obs: float, se_cdf: float,
               threshold_floor: float = 0.010) -> dict[str, Any]:  # fmt: skip
    """Interaction, closure residual and clauses 2-3 from the effects ``delta`` (pA, pp, both) on one scale."""
    inter = delta["both"] - delta["pA"] - delta["pp"]
    r = obs - delta["both"]
    se_r = math.sqrt(se_obs**2 + se_both**2 + se_cdf**2)
    thr = max(threshold_floor, 2.0 * se_r)
    return {"delta": dict(delta), "interaction": inter, "interaction_se": se_i, "observed": obs, "observed_se": se_obs,
            "closure_residual": r, "se_swap_both": se_both, "se_cdf": se_cdf, "se_r": se_r, "closure_threshold": thr,
            "clause2_holds": bool(abs(r) <= thr), "clause3": _clause3(delta["pA"], delta["pp"], delta["both"], inter)}  # fmt: skip


N_FROZEN, BATCHES_FROZEN = 100_000, 20


def expected_table_ids(scratch: Path) -> dict[str, Any]:
    """The expected elastic table id of every run variant (pinned for ``emonly``; otherwise the id recorded in the scratch cache's
    ``variants.json``), and the CDF sources recorded in the swap tables (sha256 per case)."""
    ids = mvt.variant_ids(scratch)
    ids = {"emonly": c.TABLE_ID, **ids}
    src: dict[str, list[tuple[int, str]]] = {}
    for v in SWAP_VARIANTS:
        j = json.loads((mvt._derived(scratch) / f"elastic-proton-{ids[v]}.json").read_text())
        src[v] = sorted((q["case"], q["npz_sha256"]) for q in j["exploratory"]["cdf_sources"])
    return {"ids": ids, "swap_cdf_sources": src}


def bind_run_files(files: dict[str, dict[str, Any]], energy: int, ids: dict[str, str]) -> None:
    """Fail closed (``SystemExit``) unless every run file is bound to the frozen D-E6 configuration: variant, energy, n = 1e5, 20 batches, the
    common seed of the energy, warp-cpu float64, the expected table id, equal 20 x bins array dimensions, and valid (clean, closed balance)."""
    from . import de3_probes as d3

    problems: list[str] = []
    dims = set()
    for v, d in files.items():
        want = {"variant": v, "energy": energy, "n": N_FROZEN, "n_batches": BATCHES_FROZEN, "seed": d3.SEEDS[energy],
                "backend": "warp-cpu", "precision": "float64", "elastic_table_id": ids[v],
                "mode": "emonly" if v == "emonly" else "emel"}  # fmt: skip
        got = {k: d.get(k) for k in want}
        if got != want:
            problems.append(f"{v}: {{found: {got}}} != frozen {want}")
        arr = np.asarray(d.get("idd_batches", []), dtype=float)
        if arr.ndim != 2 or arr.shape[0] != BATCHES_FROZEN or arr.shape[1] < 121:
            problems.append(f"{v}: idd_batches shape {arr.shape} (need {BATCHES_FROZEN} x >= 121)")
        dims.add(arr.shape)
        if v != "emonly" and np.asarray(d.get("recoil_idd_batches", [])).shape != arr.shape:
            problems.append(f"{v}: recoil_idd_batches shape differs from idd_batches")
        if not d3.run_valid(d):
            problems.append(
                f"{v}: run not valid (clean run, relative residual <= 1e-12, lateral part)"
            )
    if len(dims) > 1:
        problems.append(f"array dimensions differ between the runs: {sorted(dims)}")
    if problems:
        raise SystemExit(
            f"{energy} MeV: run files are not bound to the frozen D-E6 configuration (fail closed):\n  "
            + "\n  ".join(problems)
        )


def verify_cdf_error_lineage(
    err: dict[str, Any], swap_sources: dict[str, list[tuple[int, str]]]
) -> None:
    """Fail closed unless the CDF error record carries the frozen settings and the sha256 of exactly the CDF sources recorded in the swap tables."""
    problems = []
    lin = err.get("lineage", {})
    if (err.get("n_boot"), err.get("seed"), err.get("n_z"), err.get("gain_dz_mm")) != (
        200,
        CDF_SEED,
        100,
        2.0,
    ):
        problems.append(
            f"settings {(err.get('n_boot'), err.get('seed'), err.get('n_z'), err.get('gain_dz_mm'))} != (200, {CDF_SEED}, 100, 2.0)"
        )
    if err.get("source_energies_mev") != [30.0, 50.0, 100.0, 150.0, 200.0]:
        problems.append(f"source energies {err.get('source_energies_mev')}")
    if not all(len(v) == 5 for v in err.get("events_per_energy", {}).values()) or set(
        err.get("events_per_energy", {})
    ) != {"O-16", "H-1"}:
        problems.append("events per energy missing or incomplete")
    rec = sorted((q["case"], q["npz_sha256"]) for q in lin.get("cdf_sources", []))
    if not rec:
        problems.append("no CDF source lineage recorded")
    for v, src in swap_sources.items():
        if src != rec:
            problems.append(
                f"CDF sources of the swap table {v} differ from the lineage of the error record"
            )
    if problems:
        raise SystemExit("CDF error record not bound (fail closed):\n  " + "\n  ".join(problems))


def extended_arithmetic(delta: dict[str, float], se_all: float, se_i: float, obs: float, se_obs: float, se_cdf: float,
                        threshold_floor: float = 0.010) -> dict[str, Any]:  # fmt: skip
    """Addendum (c)3: ``delta`` has pA, pp, ppS (pp sigma arm) and all. I_all = Delta(all) - sum of the main effects; the closure residual
    r_all = obs - Delta(all) with SE_r = sqrt(SE_obs^2 + SE_swap^2 + SE_cdf^2); thresholds as the frozen rule."""
    main = delta["pA"] + delta["pp"] + delta["ppS"]
    inter = delta["all"] - main
    r = obs - delta["all"]
    se_r = math.sqrt(se_obs**2 + se_all**2 + se_cdf**2)
    thr = max(threshold_floor, 2.0 * se_r)
    return {"delta": dict(delta), "sum_main_effects": main, "interaction_all": inter, "interaction_all_se": se_i,
            "closure_residual_all": r, "se_swap_all": se_all, "se_cdf": se_cdf, "se_r": se_r, "closure_threshold": thr,
            "clause2_holds_extended": bool(abs(r) <= thr),
            "main_effect_fractions": {k: delta[k] / delta["all"] for k in ("pA", "pp", "ppS")},
            "interaction_fraction": abs(inter) / abs(delta["all"])}  # fmt: skip


def outcome_label(frozen_clause2: bool, extended_clause2: bool) -> str:
    """Addendum (c)3: the label of the closure outcome (clauses 4 and 5 are decided elsewhere)."""
    if frozen_clause2:
        return "closure holds on the frozen 2x2 (clause 2)"
    if extended_clause2:
        return "closure reached only with the post-hoc p-p sigma arm (declared after D-E4)"
    return "closure not reached (clause 2 fails on the frozen 2x2 and on the extended set)"


def clause3_status(holds_one_minus_f: bool, holds_minus_ln_f: bool) -> str:
    """Addendum (b)4: metric-dependent clause 3 is recorded as inconclusive, never as satisfied."""
    if holds_one_minus_f and holds_minus_ln_f:
        return "holds under both metrics"
    if not holds_one_minus_f and not holds_minus_ln_f:
        return "fails under both metrics"
    return "inconclusive (metric-dependent)"


def recoil_h4(
    recoil_batches: dict[str, np.ndarray], den_plateau: float, gain_excess: float
) -> dict[str, Any]:
    """H4: recoil-channel contribution to the plateau (plateau mean of the recoil depth channel over the EM-only plateau) per arm, its
    change against ``all`` and the share of the observed gain excess (the section 2.2 exclusion bound is 20 %)."""
    plateau = slice(40, 121)  # depth bins 20-60 mm (0.5 mm bins, centres (i + 0.5) * 0.5 mm)
    gain = {
        k: float(v.mean(axis=0)[plateau].mean() / den_plateau) for k, v in recoil_batches.items()
    }
    out = {}
    for k, g in gain.items():
        d = gain["all"] - g
        out[k] = {"recoil_plateau_gain_contribution": g, "reduction_vs_base": d, "share_of_gain_excess": d / gain_excess if gain_excess else None,
                  "ge_20pct_bound": bool(gain_excess and abs(d / gain_excess) >= 0.20)}  # fmt: skip
    return out


def summarize(
    runs: Path,
    out: Path,
    cdf_error: Path,
    ratios_path: Path | None = None,
    scratch: Path | None = None,
) -> dict[str, Any]:
    from . import de3_probes as d3

    cmp = d3._cmp()
    ref = c.v11_reference()
    arch = c.json_payload(c.RUN_V11)["verdict"]["energies"]
    err = json.loads(Path(cdf_error).read_text())
    expected = expected_table_ids(mvt.resolve_scratch(scratch))
    verify_cdf_error_lineage(err, expected["swap_cdf_sources"])
    pop = lambda x: cmp.curve_metrics(x)["peak_over_plateau"]  # noqa: E731
    plat = lambda x: cmp.curve_metrics(x)["plateau"]  # noqa: E731
    f_of = lambda n, d: pop(n) / pop(d)  # noqa: E731
    g_of = lambda n, d: plat(n) / plat(d) - 1.0  # noqa: E731
    fns = {
        "one_minus_F": lambda fx, fb_: fx - fb_,
        "minus_ln_F": lambda fx, fb_: math.log(fx / fb_),
    }
    energies: dict[str, Any] = {}
    effects: dict[str, dict[str, float]] = {}
    for e in (150, 200):
        names = ("emonly", "all", *SWAP_VARIANTS)
        files = {v: json.loads((runs / f"de3-{e}-{v}.json").read_text()) for v in names}
        bind_run_files(files, e, expected["ids"])
        a = {v: np.asarray(d["idd_batches"]) for v, d in files.items()}
        den, base = a["emonly"], a["all"]
        x = {
            "pA": a["swap-pA"],
            "pp": a["swap-pp"],
            "both": a["swap-both"],
            "ppS": a["swap-pp-sigma"],
            "all": a["swap-all"],
        }
        fb, fb_se = d3._jackknife(lambda n, d: f_of(n[0], d), [base], den)
        gb, gb_se = d3._jackknife(lambda n, d: g_of(n[0], d), [base], den)
        r = ref[e]
        se_rel = float(arch[str(e)]["F"]["se"])
        se_obs = se_rel * r["F_topas"]
        chk_p = runs / f"check-base-{e}.json"
        chk = json.loads(chk_p.read_text()) if chk_p.is_file() else {}
        diff, se_d = fb - r["F_ionmc"], math.hypot(fb_se, r["F_ionmc_se"])
        block: dict[str, Any] = {
            "base": {
                "F": fb,
                "F_se": fb_se,
                "gain": gb,
                "gain_se": gb_se,
                "F_ionmc_archived": r["F_ionmc"],
                "F_topas_archived": r["F_topas"],
                "gain_ionmc_archived": r["gain_ionmc"],
                "gain_topas_archived": r["gain_topas"],
            },  # fmt: skip
            "clause1": {
                "base_minus_archived_F": diff,
                "combined_se": se_d,
                "within_2se": bool(abs(diff) <= 2.0 * se_d),
                "table_identity_bit_for_bit": chk.get("bit_for_bit"),
                "holds": bool(abs(diff) <= 2.0 * se_d and chk.get("bit_for_bit") is True),
            },  # fmt: skip
        }
        obs = {"one_minus_F": (r["F_topas"] - r["F_ionmc"], se_obs, err[str(e)]["F_total_error"]),
               "minus_ln_F": (math.log(r["F_topas"] / r["F_ionmc"]), se_rel, err[str(e)]["F_total_error"] / fb)}  # fmt: skip
        for sc, fn in fns.items():
            o, se_o, se_c = obs[sc]
            delta, ses = {}, {}
            for k in ("pA", "pp", "both", "ppS", "all"):
                delta[k], ses[k] = d3._jackknife(
                    lambda n, d, fn=fn: fn(f_of(n[1], d), f_of(n[0], d)), [base, x[k]], den
                )
            ival, i_se = d3._jackknife(
                lambda n, d, fn=fn: fn(f_of(n[1], d), f_of(n[0], d)) - fn(f_of(n[2], d), f_of(n[0], d)) - fn(f_of(n[3], d), f_of(n[0], d)),
                [base, x["both"], x["pA"], x["pp"]], den)  # fmt: skip
            frozen = arithmetic(
                {k: delta[k] for k in ("pA", "pp", "both")}, ses["both"], i_se, o, se_o, se_c
            )
            frozen["delta_se"] = {k: ses[k] for k in ("pA", "pp", "both")}
            frozen["interaction_jackknife"] = ival
            i_all, i_all_se = d3._jackknife(
                lambda n, d, fn=fn: fn(f_of(n[1], d), f_of(n[0], d)) - sum(fn(f_of(n[j], d), f_of(n[0], d)) for j in (2, 3, 4)),
                [base, x["all"], x["pA"], x["pp"], x["ppS"]], den)  # fmt: skip
            ext = extended_arithmetic(
                {k: delta[k] for k in ("pA", "pp", "ppS", "all")},
                ses["all"],
                i_all_se,
                o,
                se_o,
                se_c,
            )
            ext["delta_se"] = {k: ses[k] for k in ("pA", "pp", "ppS", "all")}
            ext["interaction_all_jackknife"] = i_all
            block[sc] = frozen
            block[sc + "_extended"] = ext
        block["clause3_status"] = clause3_status(
            block["one_minus_F"]["clause3"]["holds"], block["minus_ln_F"]["clause3"]["holds"]
        )
        block["outcome_label"] = outcome_label(
            block["one_minus_F"]["clause2_holds"],
            block["one_minus_F_extended"]["clause2_holds_extended"],
        )
        # gain (descriptive): Delta_g(X) = gain_base - gain_X (positive = reduction, towards TOPAS)
        gain = {}
        for k in ("pA", "pp", "both", "ppS", "all"):
            dg, dg_se = d3._jackknife(lambda n, d: g_of(n[0], d) - g_of(n[1], d), [base, x[k]], den)
            rp, rp_se = d3._jackknife(
                lambda n, d: (f_of(n[1], d) - f_of(n[0], d)) / (1.0 - f_of(n[0], d)),
                [base, x[k]],
                den,
            )
            rg, rg_se = d3._jackknife(
                lambda n, d: (g_of(n[0], d) - g_of(n[1], d)) / g_of(n[0], d), [base, x[k]], den
            )
            dr, dr_se = d3._jackknife(lambda n, d: (g_of(n[0], d) - g_of(n[1], d)) / g_of(n[0], d) - (f_of(n[1], d) - f_of(n[0], d)) / (1.0 - f_of(n[0], d)), [base, x[k]], den)  # fmt: skip
            order = (
                "indeterminate"
                if abs(dr) <= 2.0 * dr_se
                else ("tail-like" if dr > 0 else "core-like")
            )
            gain[k] = {"delta_gain": dg, "delta_gain_se": dg_se, "relative_peak_change": rp, "relative_peak_change_se": rp_se,
                       "relative_gain_change": rg, "relative_gain_change_se": rg_se, "gain_minus_peak_relative": dr, "gain_minus_peak_se": dr_se,
                       "ordering": order}  # fmt: skip
        gain["observed_delta_gain"] = r["gain_ionmc"] - r["gain_topas"]
        block["gain"] = gain
        rec = {"all": np.asarray(files["all"]["recoil_idd_batches"]), **{k: np.asarray(files[v]["recoil_idd_batches"]) for k, v in
                                                                         (("pA", "swap-pA"), ("pp", "swap-pp"), ("both", "swap-both"), ("ppS", "swap-pp-sigma"), ("swap_all", "swap-all"))}}  # fmt: skip
        plat_den = float(den.mean(axis=0)[40:121].mean())
        h4 = recoil_h4(rec, plat_den, gain["observed_delta_gain"])
        for k, v in (
            ("all", "all"),
            ("pA", "swap-pA"),
            ("pp", "swap-pp"),
            ("both", "swap-both"),
            ("ppS", "swap-pp-sigma"),
            ("swap_all", "swap-all"),
        ):
            h4[k]["recoil_mev_per_primary"] = files[v].get("recoil_total_mev_per_primary")
        block["h4"] = h4
        energies[str(e)] = block
        effects[str(e)] = {
            k: block["one_minus_F"]["delta"][k] / block["one_minus_F"]["observed"]
            for k in ("pA", "pp")
        }
    rpath = ratios_path or RATIOS_PATH
    sig_in = (
        signature_inputs(json.loads(Path(rpath).read_text()), effects)
        if Path(rpath).is_file()
        else {"decision": "ratios record not found"}
    )
    doc = {"banner": c.BANNER, "energies": energies, "signature_inputs": sig_in, "cdf_finite_sample": err,
           "clauses_1_to_3_hold_at_both_energies": {sc: bool(all(energies[e]["clause1"]["holds"] and energies[e][sc]["clause2_holds"] and energies[e][sc]["clause3"]["holds"] for e in energies)) for sc in ("one_minus_F", "minus_ln_F")},  # fmt: skip
           "outcome_labels": {e: energies[e]["outcome_label"] for e in energies},
           "note": "diagnostics only, never a model; clauses 4 and 5 of plan section 5 and the outcome are decided by the synthesis"}  # fmt: skip
    inputs = sorted(
        p
        for p in Path(runs).glob("de3-*-*.json")
        if any(v in p.name for v in ("emonly", "-all", "swap-"))
    )
    inputs += [Path(cdf_error)] + ([Path(rpath)] if Path(rpath).is_file() else [])
    c.write_outputs(out, "de6-summary.json", doc, "de6-report.md",
                    c.header("D-E6 swap test: summary", "python -m validation.scripts.elastic_attribution.de6_swap --de6-summary", inputs),
                    render_report, GIT_ALLOW, GIT_ALLOW_LONG)  # fmt: skip
    return doc


GIT_ALLOW = (
    "energies/**",
    "signature_inputs/**",
    "cdf_finite_sample/**",
    "clauses_1_to_3_hold_at_both_energies/*",
    "outcome_labels/*",
)
GIT_ALLOW_LONG = (
    "signature_inputs/*/*",
    "cdf_finite_sample/source_energies_mev",
    "cdf_finite_sample/events_per_energy/*",
    "signature_inputs/energies_mev",
)
GIT_FILES = {
    "de6-summary.json": (GIT_ALLOW, GIT_ALLOW_LONG),
    "cdf-error.json": (("**",), ("source_energies_mev", "events_per_energy/*")),
}


def _f(x: Any, fmt: str = "+.4f") -> str:
    return "-" if x is None else format(x, fmt)


def render_report(doc: dict[str, Any]) -> str:
    out = [
        "",
        f"Clauses 1-3 of plan section 5 hold at both energies (rule of record, 1 - F): **{doc['clauses_1_to_3_hold_at_both_energies']['one_minus_F']}**; "
        f"(-ln F sensitivity): {doc['clauses_1_to_3_hold_at_both_energies']['minus_ln_F']}. Clauses 4 and 5 are not decided here.",
        "",
    ]
    for e, b in doc["energies"].items():
        c1, base = b["clause1"], b["base"]
        out += [
            f"## {e} MeV",
            "",
            f"Clause 1: base F {base['F']:.4f} +- {base['F_se']:.4f} against archived {base['F_ionmc_archived']:.4f}: "
            f"difference {c1['base_minus_archived_F']:+.4f} (2 SE = {2 * c1['combined_se']:.4f}), table identity bit for bit: {c1['table_identity_bit_for_bit']}, holds: {c1['holds']}.",
            "",
        ]
        for sc, label in (
            ("one_minus_F", "1 - F (rule of record; Delta = F_X - F_base, positive closes)"),
            ("minus_ln_F", "-ln F (sensitivity)"),
        ):
            s = b[sc]
            d, se = s["delta"], s["delta_se"]
            out += [f"### {label}", "", "| quantity | value | SE |", "|---|---|---|"]
            out += [
                f"| Delta({k}) | {_f(d[k])} | {_f(se[k], '.4f')} |" for k in ("pA", "pp", "both")
            ]
            out += [
                f"| interaction I | {_f(s['interaction'])} | {_f(s['interaction_se'], '.4f')} |",
                f"| observed (archived V11) | {_f(s['observed'])} | {_f(s['observed_se'], '.4f')} |",
                f"| closure residual r | {_f(s['closure_residual'])} | SE_r {_f(s['se_r'], '.4f')} (swap {_f(s['se_swap_both'], '.4f')}, CDF {_f(s['se_cdf'], '.4f')}) |",
                f"| threshold max(0.010, 2 SE_r) | {_f(s['closure_threshold'], '.4f')} | clause 2 holds: {s['clause2_holds']} |",
                f"| main effects pA / pp (fraction of Delta(both)) | {_f(s['clause3']['main_effect_pA_fraction'], '.2f')} / {_f(s['clause3']['main_effect_pp_fraction'], '.2f')} | "
                f"named {s['clause3']['named_component']}, |I|/|Delta(both)| {_f(s['clause3']['interaction_fraction'], '.2f')}, clause 3 holds: {s['clause3']['holds']} |",
                "",
            ]
            x_ = b[sc + "_extended"]
            out += [
                "Extended set (addendum (c)3; post-hoc arms declared after D-E4):",
                "",
                "| quantity | value | SE |",
                "|---|---|---|",
            ]
            out += [
                f"| Delta({k}) | {_f(x_['delta'][k])} | {_f(x_['delta_se'][k], '.4f')} |"
                for k in ("pA", "pp", "ppS", "all")
            ]
            out += [
                f"| I_all = Delta(all) - sum of main effects | {_f(x_['interaction_all'])} | {_f(x_['interaction_all_se'], '.4f')} |",
                f"| closure residual r_all | {_f(x_['closure_residual_all'])} | SE_r {_f(x_['se_r'], '.4f')} (swap {_f(x_['se_swap_all'], '.4f')}, CDF {_f(x_['se_cdf'], '.4f')}) |",
                f"| threshold max(0.010, 2 SE_r) | {_f(x_['closure_threshold'], '.4f')} | clause 2 (extended) holds: {x_['clause2_holds_extended']} |",
                f"| main effects pA shape / pp shape / pp sigma (fraction of Delta(all)) | "
                + " / ".join(_f(x_["main_effect_fractions"][k], ".2f") for k in ("pA", "pp", "ppS"))
                + f" | |I_all|/|Delta(all)| {_f(x_['interaction_fraction'], '.2f')} |",
                "",
            ]
        out += [
            f"Clause 3 of the frozen 2x2: **{b['clause3_status']}**. Outcome label: **{b['outcome_label']}**.",
            "",
        ]
        out += [
            "### Gain and peak-versus-gain ordering (descriptive)",
            "",
            "| X | Delta_g | SE | relative peak change | relative gain change | gain - peak (SE) | ordering |",
            "|---|---|---|---|---|---|---|",
        ]
        for k in ("pA", "pp", "both", "ppS", "all"):
            g = b["gain"][k]
            out.append(
                f"| {k} | {_f(g['delta_gain'], '+.5f')} | {_f(g['delta_gain_se'], '.5f')} | {_f(g['relative_peak_change'], '+.3f')} | {_f(g['relative_gain_change'], '+.3f')} | "
                f"{_f(g['gain_minus_peak_relative'], '+.3f')} ({_f(g['gain_minus_peak_se'], '.3f')}) | {g['ordering']} |"
            )
        out += [
            "",
            f"Observed gain reduction (ionmc - TOPAS, archived): {_f(b['gain']['observed_delta_gain'], '+.5f')}.",
            "",
            "### H4: recoil channel (share of the observed gain excess; the section 2.2 exclusion bound is 20 %)",
            "",
            "| arm | recoil MeV/primary | recoil plateau contribution | reduction vs base | share of gain excess | >= 20 % |",
            "|---|---|---|---|---|---|",
        ]
        for k, v in b["h4"].items():
            out.append(
                f"| {k} | {_f(v.get('recoil_mev_per_primary'), '.4f')} | {_f(v['recoil_plateau_gain_contribution'], '.5f')} | {_f(v['reduction_vs_base'], '+.5f')} | {_f(v['share_of_gain_excess'], '+.2f')} | {v['ge_20pct_bound']} |"
            )
        out.append("")
    si = doc["signature_inputs"]
    out += [
        "## Section 2.2 signature inputs (from the D-E4 ratios record; inputs, not decisions)",
        "",
    ]
    if "H1" in si:
        out += ["| quantity (30 / 50 / 100 / 150 / 200 MeV) | values |", "|---|---|"]
        for h, keys in (
            ("H1", ("rate_ratio_O", "rate_ratio_C")),
            ("H2", ("exceedance_ratio_2mm", "exceedance_ratio_3mm", "theta_med_ratio", "T_ratio")),
            ("H3", ("pp_rate_ratio", "pp_T_slow_ratio", "pp_exceedance_ratio_2mm")),
        ):
            out += [
                f"| {h} {k} | " + " / ".join(_f(x, ".3f") for x in si[h][k]) + " |" for k in keys
            ]
        out += [
            "",
            f"H1: O rate ratio outside 1 +- 0.10 at {si['H1']['n_energies_O_outside_10pct']} energies; all O and C within 1 +- 0.05: {si['H1']['all_within_5pct_O_C']}. "
            f"H2: exceedance ratio >= 1.15 at {si['H2']['n_energies_exceedance_ge_1p15']} energies, T ratio >= 1.5 at {si['H2']['T_ratio_ge_1p5_energies']}, all within 1 +- 0.10: {si['H2']['all_within_10pct']}. "
            f"H3: any p-p rate or T_slow ratio outside 1 +- 0.10: {si['H3']['any_outside_10pct']}; all within 1 +- 0.05: {si['H3']['all_within_5pct']}.",
            "",
        ]
    ce = doc["cdf_finite_sample"]
    out += [
        "## Finite-sample CDF error on the predicted F (addendum (c)4): bootstrap sd, tail-edge systematic, leave-one-energy-out (150 MeV value, proxy at 200), quadrature",
        "",
        "| E | (a) bootstrap sd | (b) tail edge | (c) leave-one-out | total |",
        "|---|---|---|---|---|",
    ]
    for e in ("150", "200"):
        out.append(
            f"| {e} | {ce[e]['bootstrap']['F_sd']:.5f} | {ce[e]['tail_edge']['dF']:+.5f} | {ce[e]['leave_one_energy_out_at_150']['dF']:+.5f} | {ce[e]['F_total_error']:.5f} |"
        )
    return "\n".join(out) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument(
        "--build", action="store_true", help="write the swap tables and the finite-sample CDF error"
    )
    ap.add_argument("--de6-summary", action="store_true")
    ap.add_argument("--scratch-cache", type=Path, default=None)
    ap.add_argument("--cdf-dir", type=Path, default=CDF_DIR)
    ap.add_argument("--runs", type=Path, default=c.OUT / "de3")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--bootstrap", type=int, default=200)
    ap.add_argument("--skip-error", action="store_true", help="--build: tables only")
    ap.add_argument(
        "--cdf-error-only",
        action="store_true",
        help="recompute only the finite-sample CDF error record (no table build)",
    )
    a = ap.parse_args(argv)
    if a.cdf_error_only:
        sources, prov = load_sources(a.cdf_dir)
        err = finite_sample_error(sources, a.bootstrap, prov=prov)
        a.out.mkdir(parents=True, exist_ok=True)
        (a.out / "cdf-error.json").write_text(json.dumps(err, indent=1, sort_keys=True) + "\n")
        print(json.dumps({e: err[e]["F_total_error"] for e in ("150", "200")}))
    if a.build:
        scratch = mvt.resolve_scratch(a.scratch_cache)
        ids = build_swap_tables(scratch, a.cdf_dir)
        print(json.dumps(ids, indent=1))
        print(json.dumps(mvt.preflight(scratch)["per_variant"].keys().__repr__()))
        if not a.skip_error:
            sources, prov = load_sources(a.cdf_dir)
            err = finite_sample_error(sources, a.bootstrap, prov=prov)
            a.out.mkdir(parents=True, exist_ok=True)
            (a.out / "cdf-error.json").write_text(json.dumps(err, indent=1, sort_keys=True) + "\n")
            print(json.dumps(err, indent=1))
    if a.de6_summary:
        summarize(a.runs, a.out, a.out / "cdf-error.json", scratch=a.scratch_cache)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
