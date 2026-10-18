"""D-E8 go/no-go evaluator (Amendments 20 (c), 21 (e), reviews A21-3 / A21-4) on synthetic inputs: CONTINUE and STOP outcomes, band recomputation
against a hand-computed Welch-t value, (iv) as report-only, the operational definitions, and every refusal (dirty tree, review id, table identity,
registered seeds, reused-document hashes, missing documents, invalid run = STOP "not established")."""

# ruff: noqa: E501
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from validation.scripts.elastic_attribution import common as c
from validation.scripts.elastic_attribution import de8_diagnostic as d8

NZ = 400
STATE = {"git_head": "a" * 40, "git_clean": True}
SCRIPTS = {n: hashlib.sha256(n.encode()).hexdigest() for n in c.CODE_FILES}
ARCH = {  # the archived values of Amendment 20 (c) (committed single-process record; synthetic file, same numbers)
    150: {
        "F": (0.8567840337079885, 0.891584942125469, 0.002885429972284705, 6.484105642092071),
        "g": (0.0220987766133649, 0.016252044281034106, 0.0008090603470251417, 19.82534536304928),
        "PP": (4.2406039099585335, 4.401920390776041, 0.004130309269421311, 16.07329835966843),
        "I": (1.9332, 1.7520, 0.03093, 16.89),
    },
    200: {
        "F": (0.8440660646539855, 0.8873703022321767, 0.0038376917410936247, 8.92),
        "g": (0.02672217441583613, 0.018419142339446415, 0.000783, 20.08),
        "PP": (3.8286, 4.0255, 0.005482, 7.41),
        "I": (2.2847836591541615, 2.140906708859845, 0.018813288632313716, 21.236523327510504),
    },
}
F_BASE = {150: 0.8580, 200: 0.8455}
GAIN_BASE = {150: 0.0221, 200: 0.0267}
PLATEAU_EM = 6.5
TABLES = {
    "m2prime": ("addd" * 16, "11" * 32),
    "m2": ("ad0d" * 16, "22" * 32),
    "rb": ("beef" * 16, "33" * 32),
}


def _curve(plateau: float, peak: float, jitter: float) -> np.ndarray:
    base = np.full(NZ, plateau)
    base[300] = peak
    rows = []
    for i in range(20):
        r = base.copy()
        r[300] *= 1.0 + 1e-3 * math.sin(3 * i + jitter)
        r[40:121] *= 1.0 + 5e-4 * math.cos(i + jitter)
        rows.append(r)
    return np.array(rows)


def _doc(
    energy: int,
    mode: str,
    variant: str,
    arr: np.ndarray,
    *,
    seed: int | None = None,
    valid: bool = True,
    table: str | None = None,
    sha: str | None = None,
) -> dict[str, Any]:
    tid, nsha = TABLES.get(variant, (None, None))
    doc: dict[str, Any] = {"energy": energy, "mode": mode, "variant": variant, "n": 100000, "n_batches": 20, "seed": seed or d8.CRN_SEEDS[energy], "backend": "warp-cpu",
                           "precision": "float64", "valid": valid, "relative_residual": 1e-16, "idd_batches": arr.tolist(), "bin_mm": 0.5}  # fmt: skip
    if variant not in ("all", "emonly"):
        doc.update(configuration="elastic_only", elastic_table_id=table or tid, elastic_table_npz_sha256=sha or nsha, recoil_idd_batches=(arr * 0.01).tolist(),
                   code={**STATE, "scripts_sha256": SCRIPTS})  # fmt: skip
    if variant == "all":
        doc["recoil_idd_batches"] = (arr * 0.012).tolist()
    return doc


def _arm_array(e: int, closure: float, dgain: float) -> np.ndarray:
    pp_em = 5.3
    delta_obs = ARCH[e]["F"][1] - ARCH[e]["F"][0]
    f_arm = F_BASE[e] + closure * delta_obs
    plateau = PLATEAU_EM * (1 + GAIN_BASE[e]) + dgain * PLATEAU_EM
    return _curve(plateau, plateau * pp_em * f_arm, 0.0)


def _inputs(
    tmp: Path,
    closure: float = 1.0,
    dgain: tuple[float, float] = (-0.0058, -0.0083),
    with_rb: bool = False,
    arm_valid: bool = True,
    p_arm: float = 0.8,
) -> dict[str, Any]:
    de3 = tmp / "de3"
    de3.mkdir()
    pp_em = 5.3
    files = {}
    for e in d8.ENERGIES:
        em = _curve(PLATEAU_EM, PLATEAU_EM * pp_em, 1.0)
        pl_b = PLATEAU_EM * (1 + GAIN_BASE[e])
        arrs = {"emonly": em, "all": _curve(pl_b, pl_b * pp_em * F_BASE[e], 0.0)}
        arrs["m2prime"] = _arm_array(e, closure, dgain[d8.ENERGIES.index(e)])
        arrs["m2"] = _arm_array(e, 0.5 * closure, dgain[d8.ENERGIES.index(e)])
        if with_rb:
            arrs["rb"] = _arm_array(e, 0.9 * closure, dgain[d8.ENERGIES.index(e)])
        for v, a in arrs.items():
            doc = _doc(
                e,
                "emonly" if v == "emonly" else "emel",
                v,
                a,
                valid=arm_valid if v == "m2prime" else True,
            )
            p = de3 / f"de3-{e}-{v}.json"
            p.write_text(json.dumps(doc))
            files[(e, v)] = p
    report = tmp / "de6-report.md"
    report.write_text(
        "# r\n\nInput files (sha256):\n\n"
        + "".join(
            f"- `validation/generated/x/{files[(e, v)].name}`  {hashlib.sha256(files[(e, v)].read_bytes()).hexdigest()}\n"
            for e in d8.ENERGIES
            for v in ("all", "emonly")
        )
    )
    steps = {"16-v11-compare": {"verdict": {"energies": {str(e): {"F": {"ionmc": {"value": ARCH[e]["F"][0]}, "reference": {"value": ARCH[e]["F"][1]}, "se": ARCH[e]["F"][2], "df": ARCH[e]["F"][3]},
                                                         "plateau_gain": {"ionmc": {"value": ARCH[e]["g"][0]}, "reference": {"value": ARCH[e]["g"][1]}, "se": ARCH[e]["g"][2], "df": ARCH[e]["g"][3]}} for e in d8.ENERGIES}}},
             "21-v5-compare": {"verdict": {"engines": {"topas": {"energies": {str(e): {"peak_over_plateau": {"ionmc": {"value": ARCH[e]["PP"][0]}, "reference": {"value": ARCH[e]["PP"][1]}, "se": ARCH[e]["PP"][2], "df": ARCH[e]["PP"][3]},
                                                                                   "delta_idd_plateau_integral": {"ionmc": {"value": ARCH[e]["I"][0]}, "reference": {"value": ARCH[e]["I"][1]}, "se": ARCH[e]["I"][2], "df": ARCH[e]["I"][3]}} for e in d8.ENERGIES}}}}}}  # fmt: skip
    archive = tmp / "archive.json"
    archive.write_text(json.dumps({"steps": steps}))
    h5 = tmp / "h5.json"
    h5.write_text(
        json.dumps(
            {
                "archived_v11": {
                    "energies": {
                        str(e): {
                            "F_ionmc": ARCH[e]["F"][0],
                            "F_topas": ARCH[e]["F"][1],
                            "gain_ionmc": ARCH[e]["g"][0],
                            "gain_topas": ARCH[e]["g"][1],
                            "F_relative_diff_se": ARCH[e]["F"][2],
                        }
                        for e in d8.ENERGIES
                    }
                }
            }
        )
    )
    de7 = tmp / "de7.json"
    de7.write_text(
        json.dumps(
            {
                "variants": {
                    "M2prime": {
                        "identity": {
                            "table_id": TABLES["m2prime"][0],
                            "npz_sha256": TABLES["m2prime"][1],
                        },
                        "trigger2": {"passed": True},
                        "trigger3": {"passes": True},
                        "trigger4": {"passes": True},
                    },
                    "M2": {
                        "identity": {"table_id": TABLES["m2"][0], "npz_sha256": TABLES["m2"][1]}
                    },
                },
                "reproduction": {
                    "rule_2pct_verdict": "all six cells within 2 % in sigma_tot, sigma_R, sigma_el"
                },
            }
        )
    )
    ratios = {
        "cases": {
            str(k): {"O-16": {"signatures": {"P_theta_lab_gt_theta_2mm": {"ratio": 1.4}}}}
            for k in range(5)
        }
    }
    (tmp / "de4.json").write_text(json.dumps(ratios))

    def d5(base: int, tid: str, p_scale: float, ident: bool = False) -> dict[str, Any]:
        cs = {
            str(k): {
                "seed": base + k,
                "classes": {
                    "O-16": {
                        "P_theta_lab_gt_theta_2mm": 0.7 * (1.0 if ident and k < 2 else p_scale),
                        "T_tail_moment": 0.1 * p_scale,
                        "theta_lab_median_deg": 3.0 * p_scale,
                    }
                },
            }
            for k in range(10)
        }
        return {
            "seed_base": base,
            "table_id": tid,
            "cases": cs,
            "code": {**STATE, "scripts_sha256": SCRIPTS},
        }

    table5 = {
        "cases": {
            str(k): {
                "classes": {
                    "O-16": {
                        "P_theta_lab_gt_theta_2mm": 0.7,
                        "T_tail_moment": 0.1,
                        "theta_lab_median_deg": 3.0,
                    }
                }
            }
            for k in range(10)
        }
    }
    (tmp / "de5-table.json").write_text(json.dumps(table5))
    for arm, base in d8.D4_SEED_BASES.items():
        (tmp / f"de5-{arm}").mkdir()
        (tmp / f"de5-{arm}" / "de5.json").write_text(
            json.dumps(d5(base, TABLES[arm][0], p_arm, ident=arm == "m2prime"))
        )
    expected = {
        a: {"table_id": TABLES[a][0], "npz_sha256": TABLES[a][1]} for a in ("m2prime", "m2")
    }
    if with_rb:
        expected["rb"] = {"table_id": TABLES["rb"][0], "npz_sha256": TABLES["rb"][1]}
    return {
        "de3": de3,
        "de5": tmp,
        "archive": archive,
        "h5": h5,
        "de7": de7,
        "de4": tmp / "de4.json",
        "de5_table": tmp / "de5-table.json",
        "report": report,
        "expected": expected,
        "files": files,
    }


def _eval(inp: dict[str, Any], **over: Any) -> dict[str, Any]:
    kw = {
        "archive": inp["archive"],
        "de7_path": inp["de7"],
        "de4_path": inp["de4"],
        "de5_table_path": inp["de5_table"],
        "de6_report": inp["report"],
        "h5_path": inp["h5"],
        "expected_tables": inp["expected"],
    }
    kw.update(over)
    return d8.evaluate(inp["de3"], inp["de5"], None, STATE, SCRIPTS, **kw)


# ---------------------------------------------------------------------------------------------------------------------------------
def test_bands_are_recomputed_and_match_a_hand_computed_welch_t_value() -> None:
    cmp = d8._cmp()
    a = {e: {"F_ion": ARCH[e]["F"][0], "F_top": ARCH[e]["F"][1], "F_se": ARCH[e]["F"][2], "F_df": ARCH[e]["F"][3], "g_se": ARCH[e]["g"][2], "g_df": ARCH[e]["g"][3],
             "PP_ion": ARCH[e]["PP"][0], "PP_top": ARCH[e]["PP"][1], "PP_se": ARCH[e]["PP"][2], "PP_df": ARCH[e]["PP"][3], "I_se": ARCH[e]["I"][2], "I_df": ARCH[e]["I"][3]} for e in d8.ENERGIES}  # fmt: skip
    b = d8.compute_bands(cmp, a, F_BASE)
    # hand computation at 150 MeV: Welch-Satterthwaite df of the relative F difference from the archived components, floored, one-sided 95 % t = 1.9432 (6 df)
    sa, sb = 0.0015196783488711665, (0.8567840337079885 / 0.891584942125469) * 0.0021600982500448217
    df = (sa**2 + sb**2) ** 2 / (sa**4 / 19.0 + sb**4 / 2.8675951621584)
    assert math.floor(df) == 6 and df == pytest.approx(6.4841, abs=1e-3)
    allowed = 0.02 - 1.9432 * 0.002885429972284705
    assert (
        b[150]["t95_F"] == 1.9432
        and b[150]["a_F"] == pytest.approx(allowed, abs=1e-12)
        and allowed == pytest.approx(0.014393, abs=1e-6)
    )
    d = a[150]["F_top"] - a[150]["F_ion"]
    assert b[150]["F_closure_lo"] == pytest.approx(
        (a[150]["F_top"] * (1 - allowed) - a[150]["F_ion"]) / d, abs=1e-12
    )
    # the printed (rounded) bands of Amendment 20 (c)3
    assert (round(100 * b[150]["F_closure_lo"], 1), round(100 * b[150]["F_closure_hi"], 1)) == (
        63.1,
        136.9,
    )
    assert (round(100 * b[200]["F_closure_lo"], 1), round(100 * b[200]["F_closure_hi"], 1)) == (
        73.6,
        126.4,
    )
    assert (round(100 * b[150]["PP_closure_lo"], 1), round(100 * b[150]["PP_closure_hi"], 1)) == (
        61.1,
        126.5,
    )
    assert (
        round(b[150]["a_gain"], 5) == 0.00360
        and round(b[200]["a_gain"], 5) == 0.00365
        and round(b[150]["a_I"], 4) == 0.0460
        and round(b[200]["a_I"], 4) == 0.0676
    )


def test_continue_outcome_and_operational_definitions(tmp_path: Path) -> None:
    inp = _inputs(tmp_path)
    doc = _eval(inp)
    assert (
        doc["verdict"] == "CONTINUE"
        and doc["failing_conditions"] == []
        and doc["established"] is True
    )
    g = doc["arms"]["m2prime"]["per_energy"]["150"]
    assert (
        g["D1_closure"] == pytest.approx(1.0, abs=2e-3)
        and abs(g["D1b_q"]) < 0.01
        and abs(g["D2_dg"]) < 0.0036
        and abs(g["D3_p"]) < 0.046
    )
    assert all(
        v["holds"]
        for e in ("150", "200")
        for v in doc["arms"]["m2prime"]["per_energy"][e]["conditions"].values()
    )
    assert g["D1_se"] >= 0.0 and g["D1b_se"] >= 0.0  # reported, not used
    # F, gain and the plateau integral by the 20-60 mm definitions
    cmp = d8._cmp()
    arm = np.asarray(json.loads((inp["de3"] / "de3-150-m2prime.json").read_text())["idd_batches"])
    base = np.asarray(json.loads((inp["de3"] / "de3-150-all.json").read_text())["idd_batches"])
    dl = (
        cmp.curve_metrics(arm.mean(axis=0))["plateau_integral_mev"]
        - cmp.curve_metrics(base.mean(axis=0))["plateau_integral_mev"]
    )
    assert g["delta_plateau_integral"] == pytest.approx(dl, rel=1e-9)
    assert dl == pytest.approx(
        arm.mean(axis=0)[40:120].sum() * 0.05 - base.mean(axis=0)[40:120].sum() * 0.05, rel=1e-9
    )
    assert g["D3_p"] == pytest.approx((ARCH[150]["I"][0] + dl) / ARCH[150]["I"][1] - 1.0, abs=1e-9)
    assert g["gain_arm"] - g["gain_base"] == pytest.approx(-0.0058, abs=2e-4)


def test_stop_outcome_lists_every_failing_condition(tmp_path: Path) -> None:
    doc = _eval(_inputs(tmp_path, closure=0.3, dgain=(0.0, 0.0)))
    assert doc["verdict"] == "STOP" and doc["established"] is True
    assert (
        "(i) at 150 MeV" in doc["failing_conditions"]
        and "(ii) at 200 MeV" in doc["failing_conditions"]
        and "(iii) at 150 MeV" in doc["failing_conditions"]
    )


def test_condition_iv_is_report_only_and_flags_identical_by_construction(tmp_path: Path) -> None:
    doc = _eval(
        _inputs(tmp_path, p_arm=1.0)
    )  # arm = table sampling: D4 never moves toward 1 in the strict sense
    d4 = doc["arms"]["m2prime"]["D4"]
    assert d4["report_only"] and d4["iv_toward_count_5"] == 0 and d4["iv_in_domain_3_of_3"] is False
    assert [
        d4["rows"][e]["identical_by_construction"] for e in ("30", "50", "100", "150", "200")
    ] == [True, True, False, False, False]
    assert not doc["arms"]["m2"]["D4"]["rows"]["30"]["identical_by_construction"]
    assert doc["verdict"] == "CONTINUE"  # (iv) is not part of the verdict
    better = _eval(
        _inputs(tmp_path / "b" if (tmp_path / "b").mkdir() is None else tmp_path, p_arm=0.8)
    )["arms"]["m2prime"]["D4"]
    assert better["rows"]["150"]["moves_toward_1"] is True and better["iv_in_domain_3_of_3"] is True


def test_construction_record_v_refers_to_the_committed_de7_record(tmp_path: Path) -> None:
    inp = _inputs(tmp_path)
    v = json.loads(inp["de7"].read_text())
    v["variants"]["M2prime"]["trigger3"]["passes"] = False
    inp["de7"].write_text(json.dumps(v))
    doc = _eval(inp)
    assert doc["verdict"] == "STOP" and "(v) construction record" in doc["failing_conditions"]
    assert (
        doc["conditions"]["v"]["construction_record"]["de7_sha256"]
        == hashlib.sha256(inp["de7"].read_bytes()).hexdigest()
    )


def test_an_invalid_published_run_is_stop_not_established_and_there_is_no_rerun_path(
    tmp_path: Path,
) -> None:
    doc = _eval(_inputs(tmp_path, arm_valid=False))
    assert (
        doc["verdict"] == "STOP"
        and doc["established"] is False
        and doc["failing_conditions"] == ["run invalid: not established"]
    )
    assert (
        doc["runs"]["arms"]["m2prime"]["150"]["executed"] is True
        and doc["runs"]["arms"]["m2prime"]["150"]["valid"] is False
    )
    assert not hasattr(d8, "rerun") and "--rerun" not in (d8.__doc__ or "")


def test_missing_documents_and_reused_hash_mismatch_refuse(tmp_path: Path) -> None:
    inp = _inputs(tmp_path)
    (inp["de3"] / "de3-200-m2prime.json").unlink()
    with pytest.raises(SystemExit, match="documents are missing"):
        _eval(inp)
    inp2 = _inputs(tmp_path / "x" if (tmp_path / "x").mkdir() is None else tmp_path)
    p = inp2["de3"] / "de3-150-all.json"
    doc = json.loads(p.read_text())
    doc["energy"] = 150
    p.write_text(json.dumps(doc) + " ")
    with pytest.raises(SystemExit, match="sha256 differs"):
        _eval(inp2)


def test_table_identity_refusals(tmp_path: Path) -> None:
    inp = _inputs(tmp_path)
    p = inp["de3"] / "de3-150-m2prime.json"
    doc = json.loads(p.read_text())
    doc["elastic_table_id"] = "f" * 64
    p.write_text(json.dumps(doc))
    with pytest.raises(SystemExit, match="table id / npz sha256 differ"):
        _eval(inp)
    doc["elastic_table_id"] = TABLES["m2prime"][0]
    doc["elastic_table_npz_sha256"] = "0" * 64
    p.write_text(json.dumps(doc))
    with pytest.raises(SystemExit, match="table id / npz sha256 differ"):
        _eval(inp)
    inp2 = _inputs(tmp_path / "e" if (tmp_path / "e").mkdir() is None else tmp_path)
    bad = copy.deepcopy(inp2["expected"])
    bad["m2prime"]["table_id"] = "c" * 64  # the manifest input must agree with de7.json
    with pytest.raises(SystemExit, match="differs from the committed de7 record"):
        _eval(inp2, expected_tables=bad)
    with pytest.raises(SystemExit, match="no expected table"):
        _eval(inp2, expected_tables={"m2": inp2["expected"]["m2"]})


def test_seed_refusals(tmp_path: Path) -> None:
    inp = _inputs(tmp_path)
    p = inp["de3"] / "de3-200-m2prime.json"
    doc = json.loads(p.read_text())
    doc["seed"] = 20990102
    p.write_text(json.dumps(doc))
    with pytest.raises(SystemExit, match="not a registered harness seed"):
        _eval(inp)
    inp2 = _inputs(tmp_path / "s" if (tmp_path / "s").mkdir() is None else tmp_path)
    p5 = inp2["de5"] / "de5-m2prime" / "de5.json"
    d5 = json.loads(p5.read_text())
    d5["cases"]["3"]["seed"] = 20991200  # the M2 block, not the M2' block
    p5.write_text(json.dumps(d5))
    with pytest.raises(SystemExit, match="registered block 20991300-20991309"):
        _eval(inp2)


def test_every_generated_document_read_is_recorded_by_sha256(tmp_path: Path) -> None:
    inp = _inputs(tmp_path)
    doc = _eval(inp)
    r = doc["runs"]
    assert (
        r["reused"]["de3-150-all.json"]
        == hashlib.sha256((inp["de3"] / "de3-150-all.json").read_bytes()).hexdigest()
    )
    assert (
        r["arms"]["m2prime"]["200"]["sha256"]
        == hashlib.sha256((inp["de3"] / "de3-200-m2prime.json").read_bytes()).hexdigest()
    )
    assert (
        r["arms"]["m2"]["de5_sha256"]
        == hashlib.sha256((inp["de5"] / "de5-m2" / "de5.json").read_bytes()).hexdigest()
    )
    assert len(r["reused"]) == 4 and doc["provenance"]["git_head"] == STATE["git_head"]


def test_code_provenance_of_documents_must_match(tmp_path: Path) -> None:
    inp = _inputs(tmp_path)
    p = inp["de3"] / "de3-150-m2prime.json"
    doc = json.loads(p.read_text())
    doc["code"]["scripts_sha256"]["de3_probes.py"] = "0" * 64
    p.write_text(json.dumps(doc))
    with pytest.raises(SystemExit, match="sha256 of de3_probes.py"):
        _eval(inp)
    doc["code"]["scripts_sha256"] = SCRIPTS
    doc["code"]["git_clean"] = False
    p.write_text(json.dumps(doc))
    with pytest.raises(SystemExit, match="clean tree required"):
        _eval(inp)


def test_review_and_clean_tree_refusals(tmp_path: Path) -> None:
    man = tmp_path / "de8-manifest.json"
    with pytest.raises(SystemExit, match="not clean"):
        d8.check_review("REVIEW-0123abcd", {**STATE, "git_clean": False}, man)
    for bad in (None, "", "review-1", "REVIEW-xyz"):
        with pytest.raises(SystemExit, match="review-id"):
            d8.check_review(bad, STATE, man)
    assert d8.check_review("REVIEW-0123abcd", STATE, man) == {}
    man.write_text(json.dumps({"reviews": {STATE["git_head"]: "REVIEW-0123abcd"}}))
    with pytest.raises(SystemExit, match="binds REVIEW-0123abcd"):
        d8.check_review("REVIEW-89abcdef", STATE, man)
    assert d8.check_review("REVIEW-0123abcd", STATE, man)["reviews"]


def test_run_writes_projection_manifest_and_a_second_run_must_reproduce(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inp = _inputs(tmp_path / "in" if (tmp_path / "in").mkdir() is None else tmp_path)
    for name, val in (
        ("ARCHIVE_JSON", inp["archive"]),
        ("DE7_JSON", inp["de7"]),
        ("DE4_RATIOS", inp["de4"]),
        ("DE5_COMMITTED", inp["de5_table"]),
        ("DE6_REPORT", inp["report"]),
        ("H5_JSON", inp["h5"]),
    ):
        monkeypatch.setattr(d8, name, val)
    monkeypatch.setattr(c, "repo_state", lambda repo=None: dict(STATE))
    monkeypatch.setattr(
        c, "code_provenance", lambda repo=None: {**STATE, "scripts_sha256": dict(SCRIPTS)}
    )
    exp = tmp_path / "expected.json"
    exp.write_text(json.dumps(inp["expected"]))
    out = tmp_path / "out"
    doc = d8.run(inp["de3"], inp["de5"], None, out, "REVIEW-0123abcd", exp)
    assert doc["verdict"] == "CONTINUE"
    proj = json.loads((out / "git" / "de8.json").read_text())
    assert (
        c.git_safety_violations(proj, d8.GIT_ALLOW) == []
        and proj["verdict"] == "CONTINUE"
        and proj["disclosure"].startswith("Go/no-go screen")
    )
    md = (out / "git" / "de8-report.md").read_text()
    head, body = md.split(c.BODY_MARKER + "\n", 1)
    assert body.rstrip("\n") == d8.render_report(proj).rstrip("\n") and d8.DISCLOSURE in body
    man = json.loads((out / "de8-manifest.json").read_text())
    assert (
        man["reviews"][STATE["git_head"]] == "REVIEW-0123abcd"
        and man["evaluation"]["verdict"] == "CONTINUE"
    )
    d8.run(
        inp["de3"], inp["de5"], None, out, "REVIEW-0123abcd", exp
    )  # the same evaluation reproduces
    p = inp["de3"] / "de3-150-m2prime.json"
    d = json.loads(p.read_text())
    d["idd_batches"] = (np.asarray(d["idd_batches"]) * 1.01).tolist()
    p.write_text(json.dumps(d))
    with pytest.raises(SystemExit, match="executed once"):
        d8.run(inp["de3"], inp["de5"], None, out, "REVIEW-0123abcd", exp)


# ---------------------------------------------------------------------------------------------------------------------------------
# variant loader, seed-base and table-id parameters
# ---------------------------------------------------------------------------------------------------------------------------------
def _builder_pair(d: Path, model: str, pauli: str | None, payload: bytes = b"npz-bytes") -> str:
    from ionmc.nuclear.build import table_id

    d.mkdir(parents=True, exist_ok=True)
    side: dict[str, Any] = {
        "npz_sha256": hashlib.sha256(payload).hexdigest(),
        "exploratory": {"model": model, "pauli_variant": pauli},
        "target_names": ["x"],
    }
    tid = table_id(side)
    side["table_id"] = tid
    (d / f"elastic-proton-{tid}.npz").write_bytes(payload)
    (d / f"elastic-proton-{tid}.json").write_text(json.dumps(side))
    return tid


def _scratch(tmp: Path) -> Path:
    from validation.scripts.elastic_attribution import make_variant_tables as mv

    s = tmp / "scratch"
    (s / "derived").mkdir(parents=True)
    (s / mv.MARKER).write_text("marker\n")
    (s / "variants.json").write_text(json.dumps({"variants": {"all": {"table_id": "a" * 64}}}))
    return s


def test_builder_table_is_installed_as_a_named_variant_with_its_binding(tmp_path: Path) -> None:
    from validation.scripts.elastic_attribution import make_variant_tables as mv

    s = _scratch(tmp_path)
    tid = _builder_pair(tmp_path / "m2p", "m2prime-pauli-blocked", "clementel-villi")
    e = mv.install_builder_table(s, "m2prime", tmp_path / "m2p")
    assert e["table_id"] == tid and e["npz_sha256"] == hashlib.sha256(b"npz-bytes").hexdigest()
    assert mv.variant_ids(s)["m2prime"] == tid and mv.verify_variant(s, "m2prime") == {
        "table_id": tid,
        "npz_sha256": e["npz_sha256"],
    }
    assert mv.install_builder_table(s, "m2prime", tmp_path / "m2p")["table_id"] == tid  # idempotent
    other = _builder_pair(
        tmp_path / "other", "m2prime-pauli-blocked", "clementel-villi", payload=b"different"
    )
    with pytest.raises(SystemExit, match="already installed"):
        mv.install_builder_table(s, "m2prime", tmp_path / "other")
    assert other != tid
    (s / "derived" / f"elastic-proton-{tid}.npz").write_bytes(b"tampered")
    with pytest.raises(SystemExit, match="differ from the recorded binding"):
        mv.verify_variant(s, "m2prime")


def test_builder_install_refuses_wrong_model_unknown_variant_and_inconsistent_ids(
    tmp_path: Path,
) -> None:
    from validation.scripts.elastic_attribution import make_variant_tables as mv

    s = _scratch(tmp_path)
    _builder_pair(tmp_path / "m2", "m2-optical-glauber", None)
    with pytest.raises(SystemExit, match="not the 'm2prime' construction"):
        mv.install_builder_table(s, "m2prime", tmp_path / "m2")
    with pytest.raises(SystemExit, match="unknown builder variant"):
        mv.install_builder_table(s, "pinned", tmp_path / "m2")
    tid = _builder_pair(tmp_path / "bad", "m2prime-pauli-blocked", "bdc-pn")
    side = json.loads((tmp_path / "bad" / f"elastic-proton-{tid}.json").read_text())
    side["exploratory"]["beta"] = (
        "changed"  # the recomputed table id no longer equals the file name
    )
    (tmp_path / "bad" / f"elastic-proton-{tid}.json").write_text(json.dumps(side))
    with pytest.raises(SystemExit, match="recomputed table id differs"):
        mv.install_builder_table(s, "rb", tmp_path / "bad")


def test_d4_seed_bases_and_table_id_arguments() -> None:
    from validation.scripts.elastic_attribution import de0_kernels as de0
    from validation.scripts.elastic_attribution import de5_sampling as de5

    assert de5.D4_SEED_BASES == d8.D4_SEED_BASES == {"m2prime": 20991300, "m2": 20991200}
    assert [cs["seed"] for cs in de5.case_table(20991300)] == list(range(20991300, 20991310)) and [
        cs["seed"] for cs in de5.case_table()
    ] == list(range(20991000, 20991010))
    assert [{k: v for k, v in a.items() if k != "seed"} for a in de5.case_table(20991200)] == [
        {k: v for k, v in b.items() if k != "seed"} for b in de5.case_table()
    ]
    with pytest.raises(SystemExit, match="needs --tag"):
        de5.main(["--table-id", "x" * 64, "--seed-base", "20991300"])
    with pytest.raises(SystemExit, match="own seed base"):
        de5.main(["--table-id", "x" * 64, "--tag", "m2prime"])
    with pytest.raises(SystemExit, match="needs --tag"):
        de0.main(["--table-id", "x" * 64])
