"""Row V11 (V3-005C C7): comparator modes emelastic / noelastic, the V11 verdict, the F_ne test, the
lv5c steps (seeds, registration, dry run of ``v11-compare``) and the frozen V11 case files."""

# ruff: noqa: E501

from __future__ import annotations

import argparse
import json
import re
from dataclasses import replace
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from ionmc.reference.runs import load_run
from tests.ionmc import test_compare_idd_v5 as T
from tests.ionmc.test_compare_batches import _write
from tests.ionmc.test_v7r_infra import REPO, _load

M = T.M
NZ = T.NZ
MOD_EM = '"g4em-standard_opt4"'
MODULES = {
    "full": [
        "g4em-standard_opt4",
        "g4h-phy_QGSP_BIC_HP",
        "g4h-elastic_HP",
        "g4stopping",
        "g4ion-binarycascade",
        "g4decay",
    ],
    "emonly": ["g4em-standard_opt4"],
    "emelastic": ["g4em-standard_opt4", "g4h-elastic_HP"],
    "noelastic": [
        "g4em-standard_opt4",
        "g4h-phy_QGSP_BIC_HP",
        "g4stopping",
        "g4ion-binarycascade",
        "g4decay",
    ],
}
EXTRA = {"full": 0.2, "emonly": 0.0, "emelastic": 0.024, "noelastic": 0.18}  # plateau offsets
ROW = {"full": "V5", "emonly": "V5", "emelastic": "V11", "noelastic": "X-elastic-factor"}
NUC = {"full": True, "emonly": False, "emelastic": False, "noelastic": True}


def mode_curve(e: int, mode: str, extra: float | None = None) -> np.ndarray:
    off = T.curve()
    c = T.curve(extra=EXTRA[mode] if extra is None else extra)
    return c * (0.90 * e / (off.sum() * 0.05))


def _input(e: int, mode: str, seed: int, names: list[str] | None = None) -> bytes:
    names = MODULES[mode] if names is None else names
    mods = f"sv:Ph/Default/Modules = {len(names)} " + " ".join(f'"{n}"' for n in names)
    return (
        f"i:Ts/Seed = {seed}\nd:So/Beam/BeamEnergy = {e} MeV\nu:So/Beam/BeamEnergySpread = 0\n"
        f"{mods}\nd:Ph/Default/CutForAllParticles = 0.05 mm\ni:So/Beam/NumberOfHistoriesInRun = 100000\n"
    ).encode()


def mode_run(tmp: Path, e: int, mode: str, seed: int, *, native: str | None = None,
             label_mode: str | None = None, extra: float | None = None) -> Path:  # fmt: skip
    """Synthetic TOPAS run of ``mode``; ``native`` gives the Modules of the input when it should differ
    from the label (``label_mode``)."""
    rng = np.random.default_rng(seed)
    idd = mode_curve(e, mode, extra) * (1.0 + 0.003 * rng.standard_normal(NZ))
    gy = idd * 1e5 * 0.05 * M.MEV_J / 0.08
    rows = "".join(f"0,0,{i},{float(v)!r},0.0\n" for i, v in enumerate(gy))
    csv = (f"# X in 1 bin of 400 mm\n# Y in 1 bin of 400 mm\n# Z in {NZ} bins of 0.5 mm\n"
           f"# DoseToMedium ( Gy ) : Sum Standard_Deviation\n{rows}")  # fmt: skip
    lm = label_mode or mode
    case = {"histories": 100000, "seeds": [seed], "input": "input.txt",
            "v5": {"row": ROW[lm], "energy_mev": e, "mode": lm, "nuclear": NUC[lm],
                   **({"modules": MODULES[lm]} if lm in ("emelastic", "noelastic") else {})}}  # fmt: skip
    names = MODULES[mode] if native is None else MODULES[native]
    files = {"work/idd_dose.csv": csv.encode(), "inputs/input.txt": _input(e, mode, seed, names)}
    out = "Particle source Beam: Total number of histories: 100000\n"
    return _write(tmp, "topas", files, case, f"T-{e}-{mode}-{seed}", stdout=out)


def freeze(root: Path, runs: list[Path]) -> Path:
    """Committed-case dir of byte copies; the family name is derived from the case label."""
    for rd in runs:
        run = load_run(rd)
        v5 = run.case["v5"]
        d = root / "topas" / f"syn-{v5['energy_mev']}-{v5['mode']}-seed{run.case['seeds'][0]}"
        k = 0
        base = d
        while d.exists():
            k += 1
            d = base.with_name(f"{base.name}x{k}")
        for rel in sorted(run.files):
            if rel.startswith("inputs/"):
                f = d / rel[len("inputs/") :]
                f.parent.mkdir(parents=True, exist_ok=True)
                f.write_bytes((rd / rel).read_bytes())
    (root / "topas").mkdir(parents=True, exist_ok=True)
    return root


def _runs(
    tmp: Path, modes: tuple[str, ...] = ("emelastic", "emonly"), energies=(150, 200)
) -> list[Path]:
    return [mode_run(tmp / f"{e}{m}{s}", e, m, 1000 * e + 10 * len(m) + s + 1)
            for e in energies for m in modes for s in range(3)]  # fmt: skip


def groups(tmp: Path, runs: list[Path], modes: tuple[str, ...] = M.MODES):
    return M.load_topas_mode_groups(runs, {150: NZ, 200: NZ}, freeze(tmp / "cases", runs), modes)


# -- comparator: modes ----------------------------------------------------------------------------
def test_mode_groups_load_with_exact_module_sets_and_fingerprint_groups(tmp_path: Path) -> None:
    runs = _runs(tmp_path, tuple(M.MODES), (150,))
    g = groups(tmp_path, runs)
    assert set(g) == {(150, m) for m in M.MODES}
    for m, grp in ((m, g[(150, m)]) for m in M.MODES):
        assert grp["curves"].shape == (3, NZ) and grp["family"] == f"syn-150-{m}"
        assert {f["group"] for f in grp["fingerprints"].values()} == {M.MODE_FP_GROUP[m]}
    assert M.MODE_MODULES["emelastic"] == frozenset({"g4em-standard_opt4", "g4h-elastic_HP"})
    assert M.MODE_MODULES["noelastic"] == M.TOPAS_MODULES_FULL - {"g4h-elastic_HP"}
    # the elastic-only and no-elastic runs are told apart although the fingerprint label is the same
    cfg = {
        m: {f["config_sha256"] for f in g[(150, m)]["fingerprints"].values()}
        for m in ("emelastic", "noelastic", "full")
    }
    assert all(len(v) == 1 for v in cfg.values()) and len(set.union(*cfg.values())) == 3


def test_wrong_module_set_is_refused(tmp_path: Path) -> None:
    for names in (
        MODULES["emelastic"] + ["g4decay"],
        ["g4h-elastic_HP"],
        MODULES["noelastic"] + ["g4h-elastic_HP", "g4em-standard_opt3"],
    ):
        bad = mode_run(tmp_path / f"b{len(names)}", 150, "emelastic", 5)
        # rewrite the native Modules line to a set that is none of the four frozen ones
        run = load_run(bad)
        assert run.engine == "topas"
        bad2 = _write(tmp_path / f"c{len(names)}", "topas",
                      {"work/idd_dose.csv": (bad / "work/idd_dose.csv").read_bytes(),
                       "inputs/input.txt": _input(150, "emelastic", 5, names)},
                      json.loads((bad / "inputs/case.json").read_text()), "T-x", stdout="Particle source Beam: Total number of histories: 100000\n")  # fmt: skip
        with pytest.raises(M.IddError, match="none of the frozen sets"):
            groups(tmp_path / f"w{len(names)}", [bad2])


def test_mixed_modes_and_mislabeled_cases_are_refused(tmp_path: Path) -> None:
    # an emelastic case bound to a noelastic run (the label says emelastic, the native input noelastic)
    bad = mode_run(tmp_path / "a", 150, "noelastic", 7, label_mode="emelastic")
    with pytest.raises(M.IddError, match="mislabeled"):
        groups(tmp_path / "x", [bad])
    # the label modules disagree with the derived set
    bad = mode_run(tmp_path / "m", 150, "emelastic", 8, native="emelastic")
    case = json.loads((bad / "inputs/case.json").read_text())
    assert case["v5"]["modules"] == MODULES["emelastic"]
    # a wrong energy label
    bad = mode_run(tmp_path / "e", 150, "emelastic", 9)
    runs_ok = [mode_run(tmp_path / f"o{s}", 150, "emelastic", 20 + s) for s in range(2)]
    other = mode_run(tmp_path / "n", 150, "noelastic", 30)
    with pytest.raises(
        M.IddError, match="need >= 3"
    ):  # 2 emelastic + 1 noelastic: no group of three
        groups(tmp_path / "y", [*runs_ok, other])
    runs3 = [mode_run(tmp_path / f"p{s}", 150, "emelastic", 40 + s) for s in range(3)]
    with pytest.raises(M.IddError, match="not accepted here"):
        groups(tmp_path / "z", [*runs3, other], modes=("emelastic", "emonly"))
    # legacy V5 loader: the V11 module set is not a frozen V5 configuration
    with pytest.raises(M.IddError, match="no V5 block"):  # the V11 label is not a V5 case
        M.load_engine_groups(runs3, "topas", {150: NZ}, freeze(tmp_path / "legc", runs3))
    run = load_run(runs3[0])
    with pytest.raises(M.IddError, match="neither the frozen full set"):  # nor a V5 native input
        M.derive_config(run, M._batches_module().run_fingerprint(run))


def test_missing_runs_are_refused(tmp_path: Path) -> None:
    runs = _runs(tmp_path, ("emelastic",), (150,))[:2]
    with pytest.raises(M.IddError, match="need >= 3"):
        groups(tmp_path, runs)


def test_legacy_v5_modes_and_loader_are_unchanged(tmp_path: Path) -> None:
    runs = [T._topas_run(tmp_path / f"t{s}", 150, True, 1 + s, 0.0) for s in range(3)]
    cases = T._freeze(tmp_path / "cases", runs)
    g = M.load_engine_groups(runs, "topas", {150: NZ}, cases)
    assert set(g) == {(150, True)} and g[(150, True)]["seeds"] == [1, 2, 3]
    # derive_config still knows only the two frozen V5 sets
    run = load_run(runs[0])
    fp = M._batches_module().run_fingerprint(run)
    assert M.derive_config(run, fp) == (150, True)
    assert M.TOPAS_MODULES_FULL == frozenset(
        MODULES["full"]
    ) and M.TOPAS_MODULES_EMONLY == frozenset(MODULES["emonly"])


# -- comparator: V11 verdict, F_ne, elastic factor -----------------------------------------------
def _v11_partial(
    path: Path, e: int, mode: str, batches: np.ndarray, seed: int, row: str | None = None
) -> None:
    doc = {"row": row or f"v11-{e}-{mode}", "mode": mode, "energy": float(e), "nuclear": mode == "emel",
           "elastic": True, "elastic_only": mode == "emel", "seed": seed, "n": 100000, "n_batches": 20,
           "bin_mm": 0.5, "unit": M.UNIT, "density_g_cm3": 1.0, "idd_batches": batches.tolist(),
           "valid": True, "reduced": False}  # fmt: skip
    doc["content_sha256"] = M.content_digest(doc)
    path.write_text(json.dumps(doc, sort_keys=True))


def _ion_dir(
    tmp: Path, extra_emel: float = 0.024, with_on: bool = False, on_extra: float = 0.2
) -> Path:
    d = tmp / "ion"
    d.mkdir()
    for e in (150, 200):
        seed = 20495004 + (e == 200)
        _v11_partial(
            d / f"v11-{e}-emel.json",
            e,
            "emel",
            T.noisy(mode_curve(e, "emelastic", extra_emel), 20, 0.004, seed + 1),
            seed,
        )
        _v11_partial(
            d / f"v11-{e}-emonly.json",
            e,
            "emonly",
            T.noisy(mode_curve(e, "emonly"), 20, 0.004, seed),
            seed,
        )
        if with_on:
            T._partial(
                d / f"v5-{e}-on.json",
                e,
                True,
                T.noisy(mode_curve(e, "full", on_extra), 20, 0.004, e),
            )
    return d


def _topas_all(tmp: Path) -> tuple[list[Path], Path]:
    runs = _runs(tmp, ("emelastic", "emonly", "full", "noelastic"))
    return runs, freeze(tmp / "cases", runs)


def test_v11_verdict_positive_and_documents_the_plan_rule(tmp_path: Path) -> None:
    runs, cases = _topas_all(tmp_path)
    doc = M.build_v11_verdict(_ion_dir(tmp_path), runs, cases)
    assert doc["row"] == "V11" and doc["pass"] is True
    assert doc["tolerances"] == {"F_relative": 0.02, "plateau_gain_absolute": 0.005}
    for e in ("150", "200"):
        r = doc["energies"][e]
        assert r["F"]["tolerance"] == 0.02 and r["plateau_gain"]["tolerance"] == 0.005
        assert r["F"]["ionmc"]["df"] == 19 and r["F"]["pass"] and r["plateau_gain"]["pass"]
        assert r["F"]["ionmc"]["value"] == pytest.approx(2.0 / 2.024, rel=0.01)
        assert doc["f_ne"][e] == {
            "computed": False,
            "reason": "v5-on partial absent; computed in v5-attribution",
        }
        assert doc["elastic_factor_rows_12_13"][e]["computed"] is True
    rule = json.dumps(doc["plan_rule"])
    for text in ("within 0.02", "within 0.005 absolute", "l.204", "Welch df"):
        assert text in rule
    assert set(doc["engine_inputs"]) == {f"{e}-{m}" for e in (150, 200) for m in M.MODES}


def test_v11_verdict_fails_on_a_wrong_elastic_effect_and_a_wrong_plateau_gain(
    tmp_path: Path,
) -> None:
    runs, cases = _topas_all(tmp_path)
    # ionmc F 4 % away (plateau offset 0.024 -> 0.14: F ~ 0.94): the F criterion fails
    (tmp_path / "a").mkdir()
    bad = M.build_v11_verdict(_ion_dir(tmp_path / "a", 0.14), runs, cases)
    assert not bad["pass"] and not bad["energies"]["150"]["F"]["pass"]
    assert bad["energies"]["150"]["F_diff_sign"] == "ionmc below TOPAS"


def test_v11_input_errors_fail_closed(tmp_path: Path) -> None:
    runs, cases = _topas_all(tmp_path)
    ion = _ion_dir(tmp_path)
    # TOPAS emelastic group absent
    only = [r for r in runs if "emelastic" not in str(r)]
    with pytest.raises(M.IddError, match="missing TOPAS groups"):
        M.build_v11_verdict(ion, only, cases)
    with pytest.raises(M.IddError, match="no TOPAS reference runs"):
        M.build_v11_verdict(ion, [], cases)
    # a partial of the wrong mode flags / unequal seeds / tampered digest
    d = json.loads((ion / "v11-150-emel.json").read_text())
    d["elastic_only"] = False
    d["content_sha256"] = M.content_digest(d)
    (ion / "v11-150-emel.json").write_text(json.dumps(d))
    with pytest.raises(M.IddError, match="nuclear|elastic_only"):
        M.build_v11_verdict(ion, runs, cases)
    ion2 = tmp_path / "i2"
    ion2.mkdir()
    for p in ion.glob("v11-*.json"):
        (ion2 / p.name).write_text(p.read_text())
    _v11_partial(
        ion2 / "v11-150-emel.json",
        150,
        "emel",
        T.noisy(mode_curve(150, "emelastic"), 20, 0.004, 1),
        20495099,
    )
    with pytest.raises(M.IddError, match="share the seed"):
        M.build_v11_verdict(ion2, runs, cases)


def test_f_ne_test_formula_interval_and_side(tmp_path: Path) -> None:
    runs, cases = _topas_all(tmp_path)
    for sub in "shl":
        (tmp_path / sub).mkdir()
    same = M.build_v11_verdict(_ion_dir(tmp_path / "s", with_on=True), runs, cases)
    f = same["f_ne"]["150"]
    assert f["computed"] and f["threshold"] == 0.02
    lo, hi = f["ratio_minus_1"]["ci90"]
    assert (
        lo < 0.0 < hi
        and not f["attributed_if_residual_above"]
        and not f["attributed_if_residual_below"]
    )
    # formula: F_ne = pp_on / pp_emel; the ratio of ionmc and TOPAS minus 1
    assert f["ratio_minus_1"]["value"] == pytest.approx(
        f["F_ne_ionmc"]["value"] / f["F_ne_TOPAS"]["value"] - 1.0
    )
    assert len(f["window_idd_on_over_topas_full_minus_1"]) == len(f["window_depth_mm"]) == 20
    # an ionmc nuclear-on plateau 12 % too low: pp_on too high -> F_ne ratio above 1.02, attributed to a residual above
    hi_doc = M.build_v11_verdict(
        _ion_dir(tmp_path / "h", with_on=True, on_extra=-0.12), runs, cases
    )
    g = hi_doc["f_ne"]["150"]
    assert (
        g["ratio_minus_1"]["ci90"][0] > 0.02
        and g["attributed_if_residual_above"]
        and not g["attributed_if_residual_below"]
    )
    low_doc = M.build_v11_verdict(_ion_dir(tmp_path / "l", with_on=True, on_extra=0.5), runs, cases)
    h = low_doc["f_ne"]["150"]
    assert (
        h["ratio_minus_1"]["ci90"][1] < -0.02
        and h["attributed_if_residual_below"]
        and not h["attributed_if_residual_above"]
    )
    assert hi_doc["pass"] is True  # F_ne never gates V11


def test_f_ne_needs_independent_seeds(tmp_path: Path) -> None:
    runs, cases = _topas_all(tmp_path)
    ion = _ion_dir(tmp_path, with_on=True)
    d = json.loads((ion / "v5-150-on.json").read_text())
    d["seed"] = 20495004
    d["content_sha256"] = M.content_digest(d)
    (ion / "v5-150-on.json").write_text(json.dumps(d))
    with pytest.raises(M.IddError, match="independent samples"):
        M.build_v11_verdict(ion, runs, cases)


def test_elastic_factor_report_only(tmp_path: Path) -> None:
    runs = _runs(tmp_path, ("full", "noelastic"), (150,))
    g = groups(tmp_path, runs)
    r = M.elastic_factor_report(g[(150, "full")]["curves"], g[(150, "noelastic")]["curves"])
    assert r["role"] == "report-only" and "pass" not in r["F_el_pp"]
    assert r["F_el_pp"]["value"] > 0 and len(r["plateau_ratio"]["ci90"]) == 2


def test_tolerances_are_constants_not_cli_options() -> None:
    assert (M.V11_TOL_F, M.V11_TOL_GAIN, M.V11_NE_TOL) == (0.02, 0.005, 0.02)
    v5c = _step_module()
    parser_args = ["v11-compare", "--tolerance", "0.5"]
    with pytest.raises(SystemExit):
        v5c.main(parser_args)
    with pytest.raises(SystemExit):
        v5c.main(["v11-compare", "--tol-f", "0.5"])
    src = (REPO / "validation/scripts/transport/steps_v5c.py").read_text()
    assert not re.search(r'add_argument\("--[a-z-]*(tol|F-)', src)


# -- steps -----------------------------------------------------------------------------------------
def _step_module() -> ModuleType:
    _load("steps_v5b")
    _load("v7r")
    return _load("steps_v5c")


def test_v11_seeds_and_step_names() -> None:
    v5c = _step_module()
    assert v5c.V11_R_INDEX == 14 and v5c.V11_K == {150: 0, 200: 1}
    assert v5c.V11_STEP_NAMES == (
        "v11-ionmc-150-emel",
        "v11-ionmc-150-emonly",
        "v11-ionmc-200-emel",
        "v11-ionmc-200-emonly",
    )
    v5c.base.SEED_BASE = 20481004
    assert (v5c.v11_seed(150), v5c.v11_seed(200)) == (20495004, 20495005)
    assert not {20495004, 20495005} & v5c.consumed_seed_set()


def test_emonly_producer_config_is_the_frozen_v5_off_physics() -> None:
    """V11's EM-only arm is the same EM model as the lv5b V5 'off' producer: the configurations differ only in
    the seed and the (inert without nuclear) elastic fields."""
    v5c = _step_module()
    v5b = v5c.v5b
    for e in (150.0, 200.0):
        geo, grid, _, _ = v5b.v5_geometry(e)
        frozen = v5b.wcfg("warp-cpu", "float64", energy=e, n=100000, seed=1, geometry=geo, grid=grid,
                          nuclear=False, n_batches=v5b.V5_BATCHES)  # fmt: skip
        mine = v5c.v11_config(e, "emonly", 100000, 2, geo, grid)
        assert mine.physics.nuclear is False and mine.physics.elastic_only is False

        def norm(c: Any) -> Any:
            c = replace(c, run=replace(c.run, seed=0),
                        physics=replace(c.physics, elastic=False, elastic_only=False, elastic_table_id=None))  # fmt: skip
            return c

        assert norm(mine) == norm(frozen)
        eff_m = v5b.Simulation(mine).effective.summary()
        eff_f = v5b.Simulation(frozen).effective.summary()
        for s in (eff_m, eff_f):
            s.pop("seed", None)
        assert (
            "elastic" not in eff_m and eff_m == eff_f
        )  # nuclear off: no elastic block in the effective config
        emel = v5c.v11_config(e, "emel", 100000, 2, geo, grid)
        assert emel.physics.nuclear and emel.physics.elastic and emel.physics.elastic_only
        assert (emel.run.backend, emel.run.precision, emel.run.n_batches) == (
            "warp-cpu",
            "float64",
            20,
        )
        assert emel.physics.elastic_table_id == v5c.ELASTIC_TABLE_ID


def test_lv5c_registration_of_v11_names_timeouts_tags_and_sources() -> None:
    _step_module()
    rs, sm = _load("run_suite"), _load("summarize")
    full = rs.full_step_names("lv5c", 2)
    names = [n.split("-", 1)[1] for n in full]
    i = names.index("v7r-diag")
    assert names[i + 1 :] == [
        "v11-ionmc-150-emel",
        "v11-ionmc-150-emonly",
        "v11-ionmc-200-emel",
        "v11-ionmc-200-emonly",
        "v11-compare",
    ]
    for n in full:
        b = n.split("-", 1)[1]
        if b.startswith("v11-ionmc"):
            assert rs.step_timeout_s("lv5c", n, 1500) == 3300 and sm.expected_tag(n) == "v11-ionmc"
        if b == "v11-compare":
            assert (
                rs.step_timeout_s("lv5c", n, 1500) == 1800 and sm.expected_tag(n) == "v11-compare"
            )
    steps = {n.split("-", 1)[1]: c for n, c, _ in rs.suite_steps("lv5c", 1, 1.0)}
    c = steps["v11-ionmc-200-emonly"]
    assert (
        c[c.index("v11-ionmc") + 1 :][:4] == ["--energy", "200", "--mode", "emonly"]
        and "--out-dir" in c
    )
    assert "--dirs" in steps["v11-compare"] and "--seed-base" in steps["v11-ionmc-150-emel"]


def test_v11_case_files_are_enumerated_hashed_and_present() -> None:
    rs = _load("run_suite")
    files = rs.V11_CASE_FILES
    assert len(files) == len(set(files)) == 24
    assert {f.split("/")[3] for f in files} == {
        f"proton-water-{e}mev-idd-r20-{m}-seed{k}"
        for e in (150, 200)
        for m in ("emelastic", "noelastic")
        for k in (1, 2, 3)
    }
    listed = rs.source_file_list("lv5c")
    assert set(files) <= set(listed) and set(rs.V5_CASE_FILES) <= set(listed)
    assert len(set(rs.V5_CASE_FILES) | set(files)) == 120
    for suite in ("lv5b", "hr5", "lv5", "lv", "hr", "lv4", "hr4"):
        assert not set(files) & set(rs.source_file_list(suite))
    for f in files:
        assert (REPO / f).is_file(), f
    assert {str(f.relative_to(rs.REPO)) for f in rs.source_files("lv5c")} >= set(files)
    dirs = {Path(f).parent for f in files}
    assert {p.relative_to(REPO).as_posix() for d in dirs for p in (REPO / d).iterdir()} == set(
        files
    )


def test_committed_v11_cases_carry_the_frozen_module_sets_and_seeds() -> None:
    root = REPO / "validation/reference_cases/topas"
    for e, rows in ((150, (10, 13)), (200, (11, 12))):
        for mode, row in zip(("emelastic", "noelastic"), rows, strict=True):
            for k in (1, 2, 3):
                d = root / f"proton-water-{e}mev-idd-r20-{mode}-seed{k}"
                v5 = json.loads((d / "case.json").read_text())["v5"]
                names = re.findall(
                    r'"([^"]+)"',
                    re.search(r"Modules = \d+ (.*)", (d / "input.txt").read_text()).group(1),
                )  # type: ignore[union-attr]
                assert set(names) == set(M.MODE_MODULES[mode]) == set(v5["modules"])
                assert v5["plan_row"] == row and v5["seed"] == 20270000 + 1000 * row + k
                assert v5["row"] == M.MODE_ROW[mode]


def _write_partials(v5c: ModuleType, out: Path, seed_base: int = 20481004) -> None:
    a = argparse.Namespace(out_dir=str(out), scale=1.0)
    v5c.base.SEED_BASE = seed_base
    for e in (150, 200):
        seed = seed_base + 14000 + (e == 200)
        for mode in ("emel", "emonly"):
            batches = T.noisy(
                mode_curve(e, "emelastic" if mode == "emel" else "emonly"),
                20,
                0.004,
                seed + (mode == "emel"),
            )
            v5c.write_partial_c(a, f"v11-{e}-{mode}", {
                "row": f"v11-{e}-{mode}", "mode": mode, "energy": float(e), "nuclear": mode == "emel",
                "elastic_only": mode == "emel", "seed": seed, "n": 100000, "n_batches": 20, "bin_mm": 0.5,
                "unit": M.UNIT, "density_g_cm3": 1.0, "idd_batches": batches.tolist(), "valid": True,
                "reduced": False})  # fmt: skip


def test_v11_compare_dry_run_on_synthetic_partials(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    v5c = _step_module()
    runs, cases = _topas_all(tmp_path)
    out = tmp_path / "samples"
    _write_partials(v5c, out)
    got: dict[str, Any] = {}
    monkeypatch.setattr(
        v5c.v5b, "finish5b", lambda doc, frozen, used, reduced: got.update(doc) or 0
    )
    monkeypatch.setattr(v5c, "_v11_reference_runs", lambda ref: runs)
    monkeypatch.setattr(
        v5c,
        "v11_bound_cases_identity",
        lambda verdict, cd: {"cases_in_source_identity": True, "bound_case_paths": []},
    )
    monkeypatch.setattr(v5c, "_cases_dir", lambda: cases)
    a = argparse.Namespace(
        dirs=[str(out)],
        partials_manifest=None,
        scale=1.0,
        reference_dir=str(tmp_path),
        out_dir=str(out),
    )
    assert v5c.step_v11_compare(a) == 0
    assert got["step"] == "v11-compare" and got["error"] is None, got.get("error")
    assert got["seeds"] == {
        "150-emel": 20495004,
        "150-emonly": 20495004,
        "200-emel": 20495005,
        "200-emonly": 20495005,
    }
    assert (
        got["r_index"] == 14
        and got["common_random_numbers"]
        and got["seeds_disjoint_from_consumed"]
    )
    assert (
        got["pass"] is True
        and got["verdict"]["row"] == "V11"
        and "V11" in json.dumps(got["plan_rule"])
    )
    assert got["v5_on_partials_used"] == []


def test_v11_seed_check_fails_closed() -> None:
    v5c = _step_module()
    v5c.base.SEED_BASE = 20481004
    sb = 20481004

    def parts(**over: Any) -> list[dict[str, Any]]:
        out = []
        for e in (150, 200):
            for m in ("emel", "emonly"):
                p = {
                    "row": f"v11-{e}-{m}",
                    "seed": sb + 14000 + (e == 200),
                    "seed_base": sb,
                    "valid": True,
                }
                out.append({**p, **over.get(p["row"], {})})
        return out

    ok = v5c.v11_seed_check(parts())
    assert ok["seeds"]["150-emel"] == ok["seeds"]["150-emonly"] == 20495004
    with pytest.raises(SystemExit, match="seed !="):
        v5c.v11_seed_check(parts(**{"v11-200-emonly": {"seed": 20495004}}))
    with pytest.raises(SystemExit, match="invalid"):
        v5c.v11_seed_check(parts(**{"v11-150-emel": {"valid": False}}))
    with pytest.raises(SystemExit, match="seed !="):
        v5c.v11_seed_check(parts()[:3])
    # a base whose V11 seeds are consumed evidence seeds
    consumed = 20471004 - 14000
    with pytest.raises(SystemExit, match="consumed"):
        v5c.v11_seed_check(
            [
                {
                    **p,
                    "seed": consumed + 14000 + (p["row"].startswith("v11-200")),
                    "seed_base": consumed,
                }
                for p in parts()
            ]
        )


def test_v11_producers_smoke_at_ci_scale(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Both producers run (2000 histories), write sealed PARTIAL lines with the shared seed and the mode
    labels; a reduced partial is refused by the comparator loader."""
    v5c = _step_module()
    v5c.base.SEED_BASE = 20481004
    docs: dict[str, Any] = {}
    monkeypatch.setattr(
        v5c.v5b, "finish5b", lambda doc, frozen, used, reduced: docs.update({doc["row"]: doc}) or 0
    )
    for mode in ("emel", "emonly"):
        a = argparse.Namespace(
            energy=150.0, mode=mode, scale=0.02, out_dir=str(tmp_path), timeout=None
        )
        assert v5c.step_v11_ionmc(a) == 0
    out = capsys.readouterr().out
    assert len(re.findall(r"^PARTIAL v11-150-(emel|emonly)\.json [0-9a-f]{64}$", out, re.M)) == 2
    p = {m: json.loads((tmp_path / f"v11-150-{m}.json").read_text()) for m in ("emel", "emonly")}
    assert p["emel"]["seed"] == p["emonly"]["seed"] == 20495004
    assert (p["emel"]["elastic_only"], p["emonly"]["elastic_only"]) == (True, False)
    assert p["emel"]["elastic"] is p["emonly"]["elastic"] is True and p["emel"]["reduced"] is True
    assert docs["v11-150-emel"]["pass"] and docs["v11-150-emonly"]["pass"]
    with pytest.raises(M.IddError, match="reduced"):
        M.load_ionmc_v11(tmp_path / "v11-150-emel.json", 150, "emel")
