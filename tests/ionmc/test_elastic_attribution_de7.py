"""D-E7 (Amendments 20 (b''), 21 (f)): projection of the M2 / M2' builder records on a small synthetic record, the withheld sigma_R of M2,
the 2 % rule on relative differences, and the leak check (no absolute cross section in the Git projection). Synthetic numbers only."""

# ruff: noqa: E501
from __future__ import annotations

import copy
import json
import math
from pathlib import Path

import pytest

from validation.scripts.elastic_attribution import common as c
from validation.scripts.elastic_attribution import de7_m2prime_record as d7

ENERGIES = (60.0, 80.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0)
EL_RATIO = {"C-12": {100: 0.60, 150: 0.98, 200: 1.13}, "O-16": {100: 0.60, 150: 1.06, 200: 1.29}}


def _nodes(target: str, scale: float = 1.0, pauli: bool = False) -> list[dict]:
    rows = []
    for e in ENERGIES:
        r_el = float(math.fsum([0.0]) + _lin(EL_RATIO[target], e)) * scale
        rows.append({
            "e_mev": e, "m2_over_bgg": r_el, "sigma_el_bgg_mb": 211.0 + e, "sigma_el_m2_mb": (211.0 + e) * r_el,
            "sigma_tot_m2_mb": 1234.5 + 3 * e, "sigma_react_m2_mb": 345.6 + e * 0.5 * scale,
            "q_first_min_fm_inv": (1.38 + 0.0005 * e) if not (target == "C-12" and e < 160) else None,
            **({"pauli_factor_peak": 0.5 + e / 1000.0} if pauli else {}),
        })  # fmt: skip
    return rows


def _lin(tab: dict[int, float], e: float) -> float:
    xs = sorted(tab)
    if e <= xs[0]:
        return tab[xs[0]]
    if e >= xs[-1]:
        return tab[xs[-1]]
    for a, b in zip(xs, xs[1:], strict=False):
        if a <= e <= b:
            return tab[a] + (tab[b] - tab[a]) * (e - a) / (b - a)
    raise AssertionError


def _info(pauli: bool, scale: float = 1.0) -> dict:
    comp = {t: {f"{e}_mev": {"measured_mb": 400.0 + e, "model_mb": 500.0 + e, "model_over_measured": f} for e, f in zip((100, 150, 200), (0.8, 1.02, 1.08), strict=True)}
            for t in ("C-12", "O-16")}  # fmt: skip
    nodes = {t: _nodes(t, scale, pauli) for t in ("C-12", "O-16")}
    level = {
        t: {f"{e}_mev": d7.interp(nodes[t], "m2_over_bgg", e) for e in (100, 150, 200)}
        for t in nodes
    }
    return {
        "table_id": "ab" * 32, "npz_sha256": "cd" * 32,
        "exploratory": {"model": "m2prime-pauli-blocked" if pauli else "m2-optical-glauber", "variant": "v", "pauli_variant": "clementel-villi" if pauli else None,
                        "pauli_emin_mev": 80.0 if pauli else None, "beta": "free-unitarity", "pp_level": "none", "unfold": "rms-matched", "inputs_unpinned": False},
        "m2": {
            "targets": {t: {"nodes": v} for t, v in nodes.items()},
            "checks": {"m2_over_bgg_sigma_el_at_100_150_200": level},
            "trigger_2_convergence": {"tol": 1e-3, "max": 2e-4, "passed": True, "n_evaluated": 9,
                                      "component_max": {"b_step": 1e-4, "s_grid": 2e-4, "radial_step": None, "q_step": None, "q_cutoff": 5e-5},
                                      "candidate_energies": [{"target": "C-12", "e_mev": 150.0, "components": {"b_step": 1e-4}, "max": 1e-4}]},
            "radius_sanity_trigger_4": {"band": [0.6, 1.0],
                                        **{t: {"helm": {"q_h_fm_inv": 1.8, "r_h_fm": 2.5}, "rows": {
                                            "150_mev": {"q_min_m2_fm_inv": "not evaluable (no true minimum)" if t == "C-12" else 1.5, "ratio_to_q_h": None if t == "C-12" else 0.83,
                                                        "in_band": None if t == "C-12" else True},
                                            "200_mev": {"q_min_m2_fm_inv": 1.6, "ratio_to_q_h": 0.89, "in_band": True}}} for t in ("C-12", "O-16")}},
            "sigma_tot_comparison_report_only": comp,
            "provenance": {"said": {"read": {"PP_SGT": {"sha256": ["11" * 32], "n_rows": 5}}}, "inputs_json_sha256": "22" * 32, "horiuchi_json_sha256": None, "inputs_unpinned": False},
            "pp_reanchor": {"scan_rows": [{"files_sha256": ["33" * 32, "44" * 32]}]},
            "nodes_below_said_range_black_disk_kept": {"count_per_target": {"C-12": 0, "O-16": 0}},
            **({"pauli": {"pauli_variant": "clementel-villi", "pauli_scale": 1.0, "pauli_emin_mev": 80.0, "rho0_fm3": 0.17, "nodes_applied_per_target": {"C-12": 43, "O-16": 43},
                          "nodes_below_emin_black_disk_per_target": {"C-12": 120, "O-16": 120}, "pauli_factor_at_peak_density": {"C-12": {"150_mev": 0.64}}}} if pauli else {}),
        },
    }  # fmt: skip


def _scratch(info: dict, bump: float = 0.0) -> str:
    lines = []
    for t in ("C-12", "O-16"):
        for e in (100, 150, 200):
            nd = info["m2"]["targets"][t]["nodes"]
            f = 1.0 + (bump if (t, e) == ("O-16", 150) else 0.0)
            tot, sr, sel = (
                d7.interp(nd, k, e) / f
                for k in ("sigma_tot_m2_mb", "sigma_react_m2_mb", "sigma_el_m2_mb")
            )
            lines.append(
                f"{t} {e} {d7.SCRATCH_LABEL:28s} sTot {tot:6.1f} (+1.0%) sR {sr:6.1f} sEl {sel:6.1f} el/Bar {d7.interp(nd, 'm2_over_bgg', e):.3f} qmin {d7.nearest(nd, 'q_first_min_fm_inv', e) or 1.0:.3f} P2 0.5000"
            )
        lines.append("   factor at peak density 0.700: R-B 0.1  R-C 0.2  R-A 0.3")
    return "\n".join(lines) + "\n"


def _doc(bump: float = 0.0) -> tuple[dict, dict, dict, str]:
    m2, m2p = _info(False), _info(True)
    text = _scratch(m2p, bump)
    return d7.build_document(m2, m2p, text), m2, m2p, text


def test_projection_has_no_absolute_cross_section_and_passes_the_leak_check() -> None:
    doc, m2, m2p, text = _doc()
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), d7.GIT_ALLOW)
    assert c.git_safety_violations(proj, d7.GIT_ALLOW) == []
    ab = d7.collect_absolutes([m2, m2p], d7.parse_scratch(text))
    assert len(ab) > 100 and d7.leak_violations(proj, ab) == []
    assert doc["leak_check"]["passed"] is True
    md = d7.render_report(proj)
    assert "withheld for M2" in md and "Trigger 3" in md and "relative differences" in md.lower()


def test_leak_check_catches_an_absolute_cross_section_and_forbidden_keys() -> None:
    doc, m2, m2p, text = _doc()
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), d7.GIT_ALLOW)
    ab = d7.collect_absolutes([m2, m2p], d7.parse_scratch(text))
    leaked = copy.deepcopy(proj)
    leaked["variants"]["M2prime"]["cells"]["C-12"]["100_mev"]["el_over_barashenkov"] = (
        345.6 + 100 * 0.5
    )  # sigma_R of a node
    assert any("equals an absolute" in b for b in d7.leak_violations(leaked, ab))
    leaked2 = copy.deepcopy(proj)
    leaked2["variants"]["M2prime"]["cells"]["C-12"]["sigma_tot_mb"] = 1.0
    assert d7.leak_violations(leaked2, ab)
    leaked3 = copy.deepcopy(proj)
    leaked3["variants"]["M2prime"]["sigma_tot_band"]["C-12"]["measured_value"] = 1.0
    assert any("measured" in b for b in d7.leak_violations(leaked3, ab))


def test_sigma_r_is_withheld_for_m2_and_projected_for_m2prime() -> None:
    doc, *_ = _doc()
    assert "react_over_barashenkov" not in doc["variants"]["M2"]["cells"]["O-16"]["150_mev"]
    assert doc["variants"]["M2"]["withheld"]
    r = doc["variants"]["M2prime"]["cells"]["O-16"]["150_mev"]["react_over_barashenkov"]
    assert r is not None and r > 0.0


def test_triggers_band_and_margins() -> None:
    doc, *_ = _doc()
    v = doc["variants"]["M2prime"]
    assert v["trigger3"]["passes"] is True  # in band at 150 and 200 for both targets
    row = v["trigger3"]["targets"]["O-16"]["rows"][2]
    assert (
        row["in_band"] is True
        and row["margin_to_nearest_bound"] == pytest.approx(0.11, abs=1e-12)
        and row["excess_over_upper_bound_rel"] < 0
    )
    assert all(
        r["matches_builder_record"] for t in v["trigger3"]["targets"].values() for r in t["rows"]
    )
    assert (
        v["trigger4"]["targets"]["C-12"]["verdict_cell"] == "200_mev"
        and v["trigger4"]["passes"] is True
    )
    b = v["sigma_tot_band"]["O-16"]
    assert (
        b["pass_per_cell"] == {"100_mev": False, "150_mev": True, "200_mev": True}
        and b["pass_at_least_2_of_3"]
        and not b["pass_all_3"]
    )
    assert v["trigger2"]["passed"] is True and v["trigger2"]["component_max"]["b_step"] == 1e-4
    assert v["pauli"]["nodes_applied_per_target"] == {"C-12": 43, "O-16": 43}
    assert v["pins"]["dsg_files_sha256"] == ["33" * 32, "44" * 32] and v["pins"][
        "said_files_sha256"
    ]["PP_SGT"] == ["11" * 32]
    # a ratio above the band: 1.407-style margin
    m2p = _info(True)
    for r in m2p["m2"]["targets"]["O-16"]["nodes"]:
        r["m2_over_bgg"] *= 1.0
    t3 = d7._trigger3({"C-12": {f"{e}_mev": {"el_over_barashenkov": x} for e, x in zip((100, 150, 200), (0.85, 1.407, 1.622), strict=True)},
                       "O-16": {f"{e}_mev": {"el_over_barashenkov": x} for e, x in zip((100, 150, 200), (0.851, 1.407, 1.622), strict=True)}}, {})  # fmt: skip
    assert t3["passes"] is False and t3["targets"]["O-16"]["rows"][1][
        "excess_over_upper_bound_rel"
    ] == pytest.approx(0.005, abs=1e-12)


def test_reproduction_record_is_relative_and_applies_the_2_percent_rule() -> None:
    doc, *_ = _doc()
    r = doc["reproduction"]
    assert r["complete"] and r["rule_2pct_verdict"].startswith("all six cells within 2 %")
    cell = r["cells"]["O-16"]["150_mev"]
    assert (
        abs(cell["rel_diff_sigma_tot"]) < 1e-3 and cell["within_2pct"] is True
    )  # print rounding only
    bad, *_ = _doc(bump=0.05)
    rb = bad["reproduction"]
    assert (
        rb["cells"]["O-16"]["150_mev"]["within_2pct"] is False
        and "defect review" in rb["rule_2pct_verdict"]
    )
    assert rb["cells"]["O-16"]["150_mev"]["rel_diff_sigma_tot"] == pytest.approx(0.05, abs=2e-3)
    missing = d7.build_document(_info(False), _info(True), None)["reproduction"]
    assert missing["complete"] is False and missing["rule_2pct_verdict"].startswith("INCOMPLETE")


def test_run_writes_projection_and_a_body_rendered_from_it(tmp_path: Path) -> None:
    m2, m2p = _info(False), _info(True)
    (tmp_path / "m2").mkdir()
    (tmp_path / "m2p").mkdir()
    (tmp_path / "m2" / "elastic-proton-x.json").write_text(json.dumps(m2))
    (tmp_path / "m2p" / "elastic-proton-y.json").write_text(json.dumps(m2p))
    (tmp_path / "scratch.txt").write_text(_scratch(m2p))
    d7.run(tmp_path / "m2", tmp_path / "m2p", tmp_path / "scratch.txt", tmp_path / "out")
    proj = json.loads((tmp_path / "out" / "git" / "de7.json").read_text())
    assert c.git_safety_violations(proj, d7.GIT_ALLOW) == []
    md = (tmp_path / "out" / "git" / "de7-report.md").read_text()
    head, body = md.split(c.BODY_MARKER + "\n", 1)
    assert (
        body.rstrip("\n") == d7.render_report(proj).rstrip("\n")
        and md.endswith("\n")
        and not md.endswith("\n\n")
    )
    assert "Input files (sha256):" in head
    full = json.loads((tmp_path / "out" / "de7.json").read_text())
    ab = d7.collect_absolutes([m2, m2p], d7.parse_scratch(_scratch(m2p)))
    assert (
        d7.leak_violations(full, ab) == []
    )  # not even the complete record carries an absolute value
    with pytest.raises(SystemExit):
        d7.load_builder_record(tmp_path)  # a directory without exactly one builder record


def test_leak_check_ignores_declared_constants_and_energy_labels() -> None:
    """False positives of the real run: the +-10 % band constant and an energy label equal to an absolute value must not fire; other floats do."""
    doc, m2, m2p, text = _doc()
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), d7.GIT_ALLOW)
    ab = d7.collect_absolutes([m2, m2p], d7.parse_scratch(text)) + [0.1, 100.0, 150.0, 80.0, 1.0]
    assert d7.leak_violations(proj, ab) == []
    assert proj["variants"]["M2"]["sigma_tot_band"]["C-12"]["band"] == 0.1
    assert proj["variants"]["M2"]["trigger2"]["candidate_energies"][0]["e_mev"] == 150.0
    leaked = copy.deepcopy(proj)
    leaked["variants"]["M2prime"]["cells"]["C-12"]["100_mev"]["q_min_fm_inv"] = 0.1234567
    assert d7.leak_violations(leaked, [*ab, 0.1234567])


def test_report_renders_when_empty_blocks_are_dropped_by_the_projection() -> None:
    """M2 has no DSG files and no Pauli block: empty lists and dicts vanish in the projection and the renderer must cope."""
    m2 = _info(False)
    m2["m2"]["pp_reanchor"] = None
    doc = d7.build_document(m2, _info(True), _scratch(_info(True)))
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), d7.GIT_ALLOW)
    assert "dsg_files_sha256" not in proj["variants"]["M2"]["pins"]
    assert "DSG files: none" in d7.render_report(proj)


def test_rb_variant_is_projected_with_both_ratios_and_checked_by_the_leak_check() -> None:
    m2, m2p, rb = _info(False), _info(True), _info(True)
    rb["table_id"], rb["exploratory"]["pauli_variant"] = "cd" * 32, "bdc-pn"
    doc = d7.build_document(m2, m2p, _scratch(m2p), rb)
    v = doc["variants"]["RB"]
    assert v["identity"]["table_id"] == "cd" * 32 and v["exploratory"]["pauli_variant"] == "bdc-pn"
    assert v["cells"]["O-16"]["150_mev"]["react_over_barashenkov"] is not None and not v.get(
        "withheld"
    )  # R-B keeps both ratios (Amendment 21 (f))
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), d7.GIT_ALLOW)
    assert (
        d7.leak_violations(
            proj, d7.collect_absolutes([m2, m2p, rb], d7.parse_scratch(_scratch(m2p)))
        )
        == []
    )
    assert "## RB" in d7.render_report(proj)
