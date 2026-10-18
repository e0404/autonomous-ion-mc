"""Row V6-F / V6-sens-F (``validation/scripts/transport/steps_v6f.py``; plan Amendment 18 (b), Amendment 20 (e)-(g),
``validation/plans/v3-005c-v6f-followup.md``): seeds and the checks (a)-(e), the point sets, the aura and beam
screens, the coverage minimum, the CPU/CUDA agreement rule, the shard plans and the producers. Synthetic inputs
except the marked producer test; sandbox results here are never evidence."""

# ruff: noqa: E501

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from tests.ionmc.test_v7r_infra import _load

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"


@pytest.fixture(scope="module")
def f6() -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    for m in ("steps_v6f", "steps_v6"):
        sys.modules.pop(m, None)
    return _load("steps_v6f")


@pytest.fixture(scope="module")
def g6(f6: ModuleType) -> ModuleType:
    return f6.s6.gv6()


PARMS = {"p1": 21.321, "p2": 0.01372, "p23": 0.0008981, "p24": 0.8092, "p25": -0.1431}


# -- seeds ---------------------------------------------------------------------------------------
def test_constants_names_and_seed_ranges(f6: ModuleType) -> None:
    s = f6.seeds_of(20551004)
    assert s["cpu"] == list(range(20558004, 20558020)) and s["sens"] == list(
        range(20559004, 20559012)
    )
    assert s["cuda"] == list(range(20560004, 20560020))
    assert (f6.V6F_R_CPU, f6.V6F_R_SENS, f6.V6F_R_CUDA, f6.V6F_BASE) == (7, 8, 9, 20551004)
    assert f6.V6F_TIMING_SEED == 20505000 + 1000 * 9 + 0 == 20514000
    assert (
        len(f6.cpu_names()) == 16 and len(f6.cuda_names()) == 16 and len(f6.sens_names()) == 7 * 8
    )
    assert f6.V6F_SENS_RUNS == ("nominal", "sx-lo", "sx-hi", "st-lo", "st-hi", "rho-lo", "rho-hi")
    assert 20551004 + 7000 > f6.V6F_OLDEST_REHEARSAL_MAX


def test_checks_a_to_e_pass_for_the_registered_seeds_and_fail_closed(f6: ModuleType) -> None:
    s = f6.seeds_of(20551004)
    assert all(f6.check_v6f_seeds(s["cpu"], s["sens"], s["cuda"]).values())
    with pytest.raises(SystemExit, match=r"\(a\)"):
        f6.check_v6f_seeds(s["cpu"], s["sens"], s["cuda"], consumed={s["cpu"][3]})
    with pytest.raises(SystemExit, match=r"\(b\)"):
        f6.check_v6f_seeds([20525999, *s["cpu"][1:]], s["sens"], s["cuda"])
    with pytest.raises(SystemExit, match=r"\(c\)"):
        f6.check_v6f_seeds(
            [20493004 + 5, *s["cpu"][1:]], s["sens"], s["cuda"], rehearsal_used_max=0
        )
    with pytest.raises(SystemExit, match=r"\(d\)"):
        f6.check_v6f_seeds(s["cpu"], s["sens"], s["cuda"], other_rows={s["cuda"][0]})
    with pytest.raises(SystemExit, match=r"\(e\)"):
        f6.check_v6f_seeds(s["cpu"], s["sens"], [*s["cuda"][:-1], s["cpu"][0]])
    with pytest.raises(SystemExit, match=r"\(e\)"):
        f6.check_v6f_seeds([*s["cpu"][:-1], s["cpu"][0]], s["sens"], s["cuda"])
    with pytest.raises(
        SystemExit
    ):  # the lv5d base shares seeds with other rows of the executing task
        t = f6.seeds_of(20531004)
        f6.check_v6f_seeds(t["cpu"], t["sens"], t["cuda"])


# -- aura, point sets ----------------------------------------------------------------------------
def test_aura_term_values(f6: ModuleType) -> None:
    vals = [f6.aura_dose(float(r), PARMS) * 1e4 for r in range(1, 7)]
    assert [round(v, 2) for v in vals] == [1.90, 1.65, 1.42, 1.21, 1.03, 0.86]  # v6f-followup 2


def _grid(g: ModuleType) -> tuple[Any, Any]:
    return np.array(g.DEPTHS_CM, float), np.array(g.RADII_CM, float)


def _meas(z: Any, r: Any) -> Any:
    """Smooth synthetic measured dose, large enough that the aura share and the floor never exclude (r0 max 10)."""
    d = 0.4 * np.exp(-0.25 * r)[None, :] * (1.0 + 0.02 * np.arange(len(z)))[:, None]
    d[:, 0] = 10.0 * np.exp(-(((z - 21.0) / 12.0) ** 2))
    return d


def test_point_sets_depth_cut_floor_aura_and_missing(f6: ModuleType, g6: ModuleType) -> None:
    z, r = _grid(g6)
    d = _meas(z, r)
    d[0, 0] = math.nan  # a missing r = 0 point
    d[10, 2] = 3.0e-4  # below the 1e-4 floor of the r = 0 maximum (10 -> 1e-3)
    d[10, 6] = 2.0e-3  # r = 6 cm: aura share ~ 0.86e-4 / 2e-3 = 0.043 -> kept
    d[11, 6] = 1.0e-3  # aura share 0.086 > 0.05 -> excluded (floor 1e-3 exactly kept)
    sets = f6.point_sets(z, r, d, PARMS)
    assert sets["z_f"] == 21.6  # deepest data depth <= 21.321 (1 + 0.01372) = 21.6135
    j_f = list(z).index(21.6)
    assert len(sets["p_f0"]) == j_f + 1 - 1 and all(
        z[j] <= 21.6 for j, _ in sets["p_f0"]
    )  # one missing
    assert len(sets["r0_report_only"]) == len(z) - (j_f + 1)  # the five deeper depths and beyond
    assert (
        (10, 2) in sets["excluded"]["below_floor"]
        and (10, 6) in sets["h_f"]
        and (11, 6) in sets["excluded"]["aura"]
    )
    assert all(j <= j_f and r[i] >= 1 for j, i in sets["h_f"])
    assert len(sets["excluded"]["beyond_z_f"]) == (len(z) - j_f - 1) * 9
    assert sets["aura_share"][(10, 6)] == pytest.approx(f6.aura_dose(6.0, PARMS) / 2.0e-3)
    assert f6.point_sets(z, r, d, PARMS, aura_max=None)["h_f"].__contains__((11, 6))
    assert len(f6.point_sets(z, r, d, PARMS, z_f=21.8)["p_f0"]) == len(sets["p_f0"]) + 1


# -- beam screen ---------------------------------------------------------------------------------
def _batches(rng: Any, mean: Any, rel: float = 0.01, nb: int = 20) -> Any:
    return mean[None] * (1.0 + rel * rng.standard_normal((nb, *mean.shape)))


def test_beam_screen_per_point_metric_and_measured_denominator(
    f6: ModuleType, g6: ModuleType
) -> None:
    z, r = _grid(g6)
    d = _meas(z, r)
    rng = np.random.default_rng(1)
    nominal = _batches(rng, d, 0.002)
    sets = f6.point_sets(z, r, d, PARMS)
    sens = {v: nominal * (1.0 + 0.001) for v in f6.V6F_BEAM_VARIANTS}
    scr = f6.beam_screen(sets["h_f"], d, nominal, sens)
    assert len(scr["h_g"]) == len(sets["h_f"]) and scr["n_screened_out"] == 0
    p = sets["h_f"][5]
    sens["st-hi"] = nominal.copy()
    sens["st-hi"][:, p[0], p[1]] *= 1.30  # a 30 % shift at one point only
    scr = f6.beam_screen(sets["h_f"], d, nominal, sens)
    assert p not in scr["h_g"] and len(scr["h_g"]) == len(sets["h_f"]) - 1
    assert scr["metric_max"][p] == pytest.approx(0.30, abs=0.02)
    sens["st-hi"] = nominal.copy()
    sens["st-hi"][:, p[0], p[1]] *= 1.10  # 10 % + 2 sigma' < 0.15 stays
    assert p in f6.beam_screen(sets["h_f"], d, nominal, sens)["h_g"]
    # contrary readings (c) 1 sigma and (d) the model-normalised metric run on the same inputs
    assert p in f6.beam_screen(sets["h_f"], d, nominal, sens, n_sigma=1.0)["h_g"]
    assert len(f6.beam_screen(sets["h_f"], d, nominal, sens, model_normalised=True)["h_g"]) <= len(
        sets["h_f"]
    )
    with pytest.raises(ValueError):
        f6.beam_screen(sets["h_f"], d, nominal, {"sx-lo": nominal})
    bad = {v: nominal for v in f6.V6F_BEAM_VARIANTS}
    bad["sx-lo"] = nominal[:5]
    with pytest.raises(ValueError):
        f6.beam_screen(sets["h_f"], d, nominal, bad)


# -- sign test and coverage minimum ----------------------------------------------------------------
def test_sign_test_ties_and_coverage_minimum(f6: ModuleType, g6: ModuleType) -> None:
    z, r = _grid(g6)
    n_z, n_r = len(z), len(r)
    resid = np.zeros((n_z, n_r))
    # depths 0..11: r = 1 cm and r = 4 cm points; the 4 cm residual is larger -> sign +
    h_g = [(j, 1) for j in range(14)] + [(j, 4) for j in range(14)]
    resid[:, 4] = 0.05
    resid[12:14, 4] = 0.0  # two ties
    st = f6.sign_test(h_g, z, r, resid)
    assert (st["n_pos"], st["n_neg"], st["n_ties_excluded"], st["n_depths_after_ties"]) == (
        12,
        0,
        2,
        12,
    )
    assert st["p"] == pytest.approx(2 / 2**12) and st["pass"] is False  # 12/12 same sign fails
    cov = f6.coverage_minimum(h_g, r, st)
    assert (
        cov["ok"] is False and cov["n_points"] == 28 and cov["n_radii"] == 2
    )  # < 40 points, < 3 radii
    more = [*h_g, *[(j, 3) for j in range(14)]]
    assert f6.coverage_minimum(more, r, st)["ok"] is True  # 42 points, 3 radii, 12 depths
    resid2 = resid.copy()
    resid2[:4, 4] = -0.05  # 8 +, 4 - : p = 0.39
    assert f6.sign_test(h_g, z, r, resid2)["pass"] is True
    assert f6.sign_test([], z, r, resid)["p"] == 1.0
    one = {"n_depths_after_ties": 11}
    assert f6.coverage_minimum(more, r, one)["ok"] is False  # 11 depths < 12


# -- the verdict ---------------------------------------------------------------------------------
def _scenario(
    f6: ModuleType,
    g6: ModuleType,
    *,
    cpu_bias: Any = None,
    cuda_bias: Any = None,
    sens_scale: float = 1.001,
    seed: int = 3,
) -> dict[str, Any]:
    z, r = _grid(g6)
    d = _meas(z, r)
    rng = np.random.default_rng(seed)

    def model(bias: Any) -> Any:
        m = d.copy()
        if bias is not None:
            m = bias(m)
        return _batches(rng, m, 0.001)

    nominal = _batches(rng, d, 0.002)
    sens = {v: nominal * sens_scale for v in f6.V6F_BEAM_VARIANTS}
    return f6.evaluate_v6f(z, r, d, PARMS, model(cpu_bias), model(cuda_bias), nominal, sens)


def test_perfect_model_passes_and_reports_contrary_readings(f6: ModuleType, g6: ModuleType) -> None:
    v = _scenario(f6, g6)
    assert v["verdict"] == "pass" and v["pass"] is True
    assert v["sets"]["n_p_f0"] == 35 and v["sets"]["n_h_g"] == v["sets"]["n_h_f"] > 40
    assert v["cpu"]["coverage"]["ok"] and v["backend_agreement"]["agree"] is True
    assert set(v["contrary_readings"]) >= {"a_no_aura_exclusion", "b_z_f_next_depth", "c_one_sigma_screen",
                                           "d_model_normalised_screen", "e_screened_out_points"}  # fmt: skip
    assert v["contrary_readings"]["e_screened_out_points"] == {"n_points": 0}
    assert len(v["sets"]["aura_share_points"]) == 35 * 9 - 0 or v["sets"]["aura_share_points"]
    json.dumps(v, default=float)  # the document is serialisable


def test_r0_miss_fails_and_halo_clause_failure_fails(f6: ModuleType, g6: ModuleType) -> None:
    def r0_off(m: Any) -> Any:
        m[5, 0] *= 1.03
        return m

    v = _scenario(f6, g6, cpu_bias=r0_off, cuda_bias=r0_off)
    assert v["cpu"]["r0"]["pass"] is False and v["verdict"] == "fail" and v["pass"] is False

    def halo_off(m: Any) -> Any:
        m[:, 1:] *= 1.20  # every halo point 20 % off: 0 % pass the 15 % tolerance
        return m

    v = _scenario(f6, g6, cpu_bias=halo_off, cuda_bias=halo_off)
    assert v["cpu"]["halo"]["status"] == "fail" and v["verdict"] == "fail"
    assert v["cpu"]["halo"]["frac_pass"] == 0.0


def test_beam_sensitive_halo_is_not_established_never_a_pass(
    f6: ModuleType, g6: ModuleType
) -> None:
    v = _scenario(f6, g6, sens_scale=1.30)  # every point moves 30 % under every variant
    assert v["sets"]["n_h_g"] == 0 and v["beam_screen"]["n_screened_out"] == v["sets"]["n_h_f"]
    assert v["cpu"]["halo"]["status"] == f6.NOT_ESTABLISHED and v["cuda"]["coverage"]["ok"] is False
    assert v["verdict"] == f6.NOT_ESTABLISHED and v["pass"] is False


def test_backend_disagreement_fails_and_large_z_is_listed(f6: ModuleType, g6: ModuleType) -> None:
    def halo_off(m: Any) -> Any:
        m[:, 1:] *= 1.20
        return m

    v = _scenario(f6, g6, cuda_bias=halo_off)  # CPU passes, CUDA halo fails
    ag = v["backend_agreement"]
    assert ag["agree"] is False and {"halo", "overall"} <= set(ag["disagreement"])
    assert v["verdict"] == "fail" and v["pass"] is False
    assert len(ag["findings_abs_z_gt_4"]) > 0 and all(
        abs(q["z"]) > 4 for q in ag["findings_abs_z_gt_4"]
    )

    def r0_off(m: Any) -> Any:
        m[3, 0] *= 1.05
        return m

    v = _scenario(f6, g6, cuda_bias=r0_off)
    assert "r0" in v["backend_agreement"]["disagreement"] and v["pass"] is False
    # tiny backend difference inside the noise: no finding, agreement holds
    assert _scenario(f6, g6, seed=9)["backend_agreement"]["findings_abs_z_gt_4"] == []


# -- shard plans, timing -------------------------------------------------------------------------
def test_cpu_plans_fit_and_cuda_plan_follows_the_measured_rate(f6: ModuleType) -> None:
    v5d = _load("steps_v5d")
    rate = 1e6 / (184 * 1.672)  # planning: ~308 s per 1e6-history shard
    table = v5d.shard_plan_c(f6.v6f_plans(), {"v5": rate}, 1500.0, 1.25)
    assert all(t["ok"] for t in table.values()) and len(table) == 2
    assert [t["declared_shards"] for t in table.values()] == [16, 56]
    assert all(
        t["planned_s_per_declared_shard"] == pytest.approx(184 * 1.672, rel=1e-6)
        for t in table.values()
    )
    slow = v5d.shard_plan_c(f6.v6f_plans(), {"v5": 1e6 / 1300}, 1500.0, 1.25)
    assert not any(t["ok"] for t in slow.values())  # 1300 s per shard > 1200 s effective
    ok = f6.cuda_shard_plan(1e6 / 37.0)  # planning ~ 37 s per 1e6 CUDA histories
    assert ok["ok"] is True and ok["planned_s_per_declared_shard"] == pytest.approx(37.0)
    bad = f6.cuda_shard_plan(1e6 / 1300.0)
    assert bad["ok"] is False and bad["shards_needed_for_planned_limit"] == math.ceil(
        16 * 1300 / 1200
    )
    with pytest.raises(ValueError):
        f6.cuda_shard_plan(0.0)


# -- producers -----------------------------------------------------------------------------------
@pytest.mark.requires_tables("nuclear", "elastic")
def test_producers_seeds_labels_and_reduced_compare(
    f6: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    s6 = f6.s6
    monkeypatch.setattr(s6, "_CUDA_TARGET", ("warp-cpu", "float64"))  # no GPU in the sandbox
    out = tmp_path / "samples"
    common = ["--scale", "0.002", "--out-dir", str(out), "--seed-base", "20551004"]
    assert f6.main(["v6f-shard", "--shard", "1", *common]) == 0
    assert f6.main(["v6f-sens-shard", "--variant", "st-hi", "--shard", "2", *common]) == 0
    assert f6.main(["v6f-cuda-shard", "--shard", "3", *common]) == 0
    text = capsys.readouterr().out
    assert text.count("PARTIAL ") == 3
    cpu = json.loads((out / "v6f-s1.json").read_text())
    sens = json.loads((out / "v6f-sens-st-hi-s2.json").read_text())
    cuda = json.loads((out / "v6f-cuda-s3.json").read_text())
    assert (cpu["seed"], cpu["r_index"], cpu["row"], cpu["kind"]) == (20558005, 7, "V6-F", "v6f")
    assert (sens["seed"], sens["r_index"], sens["row"], sens["variant"]) == (
        20559006,
        8,
        "V6-sens-F",
        "st-hi",
    )
    assert (cuda["seed"], cuda["r_index"], cuda["row"]) == (20560007, 9, "V6-F")
    assert cpu["beam"]["E0_MeV"] == s6.beam_table()["nominal"]["E0_MeV"] and cpu["valid"] is True
    assert sens["beam"] == s6.beam_table()["st-hi"]
    with pytest.raises(SystemExit, match="--variant"):
        f6.main(["v6f-sens-shard", "--variant", "r0-lo", *common])
    with pytest.raises(SystemExit, match="outside"):
        f6.main(["v6f-shard", "--shard", "16", *common])
    with pytest.raises(SystemExit, match="consumed"):
        f6.main(
            [
                "v6f-shard",
                "--shard",
                "0",
                "--scale",
                "0.002",
                "--seed-base",
                "20481004",
                "--out-dir",
                str(out),
            ]
        )
    # the existing V6 steps keep their labels
    assert s6.step_v6_shard is not None and s6.STEPS["v6-shard"] is s6.step_v6_shard


@pytest.mark.requires_tables("nuclear", "elastic")
def test_timing_step_is_timing_only(f6: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                    capsys: pytest.CaptureFixture[str]) -> None:  # fmt: skip
    monkeypatch.setattr(f6.s6, "_CUDA_TARGET", ("warp-cpu", "float64"))
    rc = f6.main(
        ["v6f-timing", "--scale", "0.002", "--out-dir", str(tmp_path), "--seed-base", "20505000"]
    )
    text = capsys.readouterr().out
    doc = json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])
    assert (
        rc == 0
        and doc["seed"] == 20514000
        and doc["reduced"] is True
        and "timing only" in doc["row"]
    )
    assert doc["plan"]["declared_shards"] == 16 and doc["hist_per_s"] > 0
    assert "PARTIAL" not in text and not list(tmp_path.glob("*.json"))


# -- the consumed-seed guard of the producers (deterministic, no tables, no cache) ----------------------
def test_producers_refuse_consumed_rehearsal_and_original_v6_seeds(
    f6: ModuleType, tmp_path: Path
) -> None:
    ok = f6.seeds_of(20551004)
    for s in (*ok["cpu"], *ok["sens"], *ok["cuda"]):
        f6.check_producer_seed(s)  # the registered seeds pass
    for base in (20481004, 20491004, 20471004):  # lv5c / hr5c / lv5b evidence bases: seeds observed
        with pytest.raises(SystemExit, match="consumed"):
            f6.check_producer_seed(base + 1000 * 7 + 3)
    with pytest.raises(SystemExit, match=r"\(b\)"):
        f6.check_producer_seed(20505000 + 7000, consumed=set())  # rehearsal family
    with pytest.raises(SystemExit, match=r"\(c\)"):
        f6.check_producer_seed(
            20493004 + 5, consumed=set()
        )  # original V6 seeds, never run, not reused
    with pytest.raises(SystemExit, match="consumed"):
        f6.check_producer_seed(1, consumed={1})
    common = ["--scale", "0.002", "--out-dir", str(tmp_path)]
    for step, extra in (
        ("v6f-shard", []),
        ("v6f-sens-shard", ["--variant", "st-hi"]),
        ("v6f-cuda-shard", []),
    ):
        with pytest.raises(SystemExit, match="consumed"):  # refused before any transport runs
            f6.main([step, "--shard", "0", *extra, *common, "--seed-base", "20481004"])
    assert not list(tmp_path.glob("*"))
