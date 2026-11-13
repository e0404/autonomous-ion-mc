"""V3-008 P1c (FC-1): committed attribution of the Pb residual and the compound/water check.

Reproduces the recorded aggregates of ``validation/scripts/ion_em/p1c_attribution.py`` from the
hash-pinned inputs: the arXiv:0908.1413 LaTeX object, the arXiv:1610.01279 PDF (text via
``mutool``/``pdftotext``) and the cached NIST PSTAR responses (never copied into Git; only
ratios are asserted). Inputs are looked up in ``$IONMC_CACHE_DIR[/objects]`` and the repository
cache (``.ionmc-cache/reference``, ``.ionmc-cache/ionmc-data/objects``). When absent the tests
skip, and FAIL under ``IONMC_REQUIRE_DATA=1`` (repository cached-data convention).

The Pb clearing is an ATTRIBUTION (input difference in the range table), NOT a pass of the
committed Pb test, which keeps its strict xfail on ionmc Bethe ranges.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc.materials import ALUMINIUM, BERYLLIUM, COPPER, LEAD
from ionmc.physics.stopping import BetheStoppingSource

_PATH = Path(__file__).resolve().parents[2] / "validation/scripts/ion_em/p1c_attribution.py"
_SPEC = importlib.util.spec_from_file_location("p1c_attribution", _PATH)
assert _SPEC is not None and _SPEC.loader is not None
P1C = importlib.util.module_from_spec(_SPEC)
sys.modules["p1c_attribution"] = P1C
_SPEC.loader.exec_module(P1C)

MATS = {"Be": BERYLLIUM, "Al": ALUMINIUM, "Cu": COPPER}


def _need(fn: Any, *args: Any) -> Any:
    try:
        return fn(*args)
    except P1C.InputUnavailable as exc:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail(f"IONMC_REQUIRE_DATA=1 but {exc}")
        pytest.skip(str(exc))


@pytest.fixture(scope="module")
def tbl() -> dict[str, Any]:
    return P1C.theta_table(_need(P1C.gottschalk_tex))


@pytest.fixture(scope="module")
def makarova() -> dict[str, Any]:
    """Table 1 of arXiv:1610.01279 (published Go93 / LOOKUP theta values), parsed from the PDF."""
    return P1C.parse_makarova_table1(_need(P1C.makarova_pdf_text))


@pytest.fixture(scope="module")
def pb(tbl: dict[str, Any]) -> dict[str, Any]:
    return P1C.pb_attribution("Pb", LEAD, tbl, quad_ns=(400, 3200), dense_table=False)


@pytest.fixture(scope="module")
def pb_pstar(tbl: dict[str, Any]) -> dict[str, Any]:
    return P1C.pstar_element_check("Pb", LEAD, tbl)


def _mx(d: dict[Any, Any]) -> float:
    return max(abs(v) for v in d.values() if v is not None)


# ------------------------------------------------------------------------------- the parsers
@pytest.mark.parametrize("matno", sorted(P1C.PSTAR_SHA256))
def test_pstar_parser(matno: str) -> None:
    name, e, r = _need(P1C.pstar_data, matno)
    assert name == P1C.PSTAR_NAMES[matno]
    assert len(e) == len(r) > 100
    assert np.all(np.diff(e) > 0) and np.all(np.diff(r) > 0)
    assert 1e-3 <= e[0] < 0.01 and e[-1] >= 1e4  # the default PSTAR grid, 1 keV to 10 GeV


def test_gottschalk_table_shape(tbl: dict[str, Any]) -> None:
    assert set(tbl) == {"Be", "Al", "Cu", "Pb"}
    assert all(len(v["rows"]) == 7 for v in tbl.values())
    assert [f for f, *_ in tbl["Pb"]["rows"]] == [0.001, 0.01, 0.1, 0.2, 0.5, 0.9, 0.97]


def test_makarova_table_parsed_from_pdf(makarova: dict[str, Any]) -> None:
    """Transcription check against published values of arXiv:1610.01279 Table 1 (Hanson column
    of the LOOKUP; independent of any NIST/ICRU data)."""
    assert list(makarova) == list(P1C.MAKAROVA_MATERIALS)
    assert makarova["Polystyrene"][0][:3] == (0.347, 3.346, 3.289)
    assert makarova["Polystyrene"][1][:3] == (15.751, 42.031, 41.973)
    assert makarova["Lead"][2][:3] == (31.566, 175.421, 186.292)
    assert makarova["Teflon"][0][:3] == (0.055, 1.626, 1.353)
    assert [len(v) for v in makarova.values()] == [3, 2, 2, 2, 2, 2, 3, 3, 3, 2, 2, 2, 3, 2]


def test_slab_theta_equals_committed_reference_at_defaults(tbl: dict[str, Any]) -> None:
    """The local variant function is the committed reference when no variant is requested."""
    el = P1C.elements(LEAD)
    _t, r1, path = P1C.ours_path(LEAD)
    for f in (0.1, 0.9):
        ref = P1C.MH.moliere_slab(1.0, P1C.MP, el, f * r1, path)
        mine = P1C.slab_theta(el, f * r1, path)
        for k, v in ref.items():
            assert mine[k] == pytest.approx(v, rel=1e-13, abs=0.0)
    assert P1C.theta(el, 1.0, path, neg_u_in_default=P1C.MH.U_IN_NEG_DEFAULT) == P1C.theta(
        el, 1.0, path
    )


# --------------------------------------------------------------------- (1) Pb attribution
def test_pb_baseline_on_ionmc_ranges_fails_as_recorded(pb: dict[str, Any]) -> None:
    """The committed Pb test (ionmc Bethe ranges) exceeds 1 %: +1.10 / +1.32 % at x/R1 0.9 / 0.97
    (it keeps its strict xfail; this is the failure being attributed)."""
    base = pb["baseline"]
    assert base[0.9] == pytest.approx(0.0110, abs=2e-4)
    assert base[0.97] == pytest.approx(0.0132, abs=2e-4)
    assert _mx(base) > 0.01


def test_pb_cleared_by_attribution_on_icru49_node_path(pb: dict[str, Any]) -> None:
    """Attribution, not a pass of the committed test: on Gottschalk's own ln R(ln T) spline through
    his (E, residual range) nodes the reference reproduces his Pb theta_Hanson to <= 0.05 %."""
    spl = {k: v for k, v in pb["icru49_path"].items() if k.startswith("spline") and "max" not in k}
    assert set(spl) == {"spline_lnR_lnT_nodes_ge_0.1", "spline_lnR_lnT_all_nodes"}
    for v in spl.values():
        assert _mx(v) <= 5e-4
    kap = pb["icru49_path"]["kappa_nodes_ge_0.1"]
    assert _mx(kap) <= 5e-4  # the kappa(E)-scaled ionmc table gives the same answer


def test_pb_cleared_on_pstar_path(pb_pstar: dict[str, Any]) -> None:
    assert _mx(pb_pstar["theta_on_pstar_path_over_thetaH_minus_1"]) <= 3e-4
    assert abs(pb_pstar["R1_pstar_over_G_minus_1"]) <= 7e-4  # PSTAR vs Gottschalk's rho R1


def test_pb_range_input_difference(pb: dict[str, Any]) -> None:
    assert pb["rhoR1_ours_over_G_minus_1"] == pytest.approx(0.00846, abs=1e-4)
    # up to +3 % in residual range at 20-41 MeV: the cause of the baseline deviation
    assert max(pb["residual_range_ours_at_EG_over_G_minus_1"].values()) > 0.03


def test_pb_fano_and_quadrature_controls(pb: dict[str, Any]) -> None:
    """Quadrature/table density do not matter: change vs baseline <= 1.2e-6; the Fano -u_in
    variants move theta by far less than the baseline residual explained by the ranges."""
    base = pb["baseline"]
    for key, vals in pb["numerics"].items():
        assert max(abs(vals[f] - base[f]) for f in base) <= 1.2e-6, key
    assert "uniform_GL_4000" in pb["numerics"]
    # -u_in = 5 (the committed default) reproduces the baseline; 3.6 and 6.0 differ
    d5 = pb["electrons"]["fano_neg_u_in_5.0"]
    assert max(abs(d5[f] - base[f]) for f in base) < 1e-12


@pytest.mark.parametrize("name", sorted(MATS))
def test_controls_be_al_cu(tbl: dict[str, Any], name: str) -> None:
    out = P1C.pb_attribution(name, MATS[name], tbl, quad_ns=(400,), dense_table=False)
    assert max(_mx(v) for k, v in out["icru49_path"].items() if not k.endswith("_max")) <= 1.1e-3
    ps = P1C.pstar_element_check(name, MATS[name], tbl)
    assert _mx(ps["theta_on_pstar_path_over_thetaH_minus_1"]) <= 1.1e-3


# ------------------------------------------------------------ (2) compound / water check
@pytest.fixture(scope="module")
def compound(makarova: dict[str, Any]) -> dict[str, Any]:
    return P1C.compound_check(makarova, P1C.gottschalk_tex(), BetheStoppingSource())


def _lookup(c: dict[str, Any], key: str) -> list[float]:
    return [e["ours_over_theta_hanson_LOOKUP_minus_1"] for e in c["makarova_tab_results"][key]]


def test_compound_lookup_ratios(compound: dict[str, Any]) -> None:
    """ours/LOOKUP - 1 (theory against theory; bound 1 %), tolerance +-0.03 pp of the record."""
    tol = 3e-4
    poly = _lookup(compound, "Polystyrene")
    assert poly[0] == pytest.approx(-0.0016, abs=tol)
    assert poly[1] == pytest.approx(-0.0020, abs=tol)
    lex = _lookup(compound, "Lexan_C16H14O3")  # stoichiometric polycarbonate
    assert lex[0] == pytest.approx(-0.0015, abs=tol)
    assert lex[1] == pytest.approx(-0.0010, abs=tol)
    luc = _lookup(compound, "Lucite")
    assert luc[0] == pytest.approx(-0.0018, abs=tol)
    assert luc[1] == pytest.approx(-0.0015, abs=tol)
    ny = _lookup(compound, "Nylon")
    assert ny[0] == pytest.approx(-0.0025, abs=tol)
    assert ny[1] == pytest.approx(-0.0016, abs=tol)
    assert _lookup(compound, "Teflon")[0] == pytest.approx(-0.0005, abs=tol)
    # the Go93-listed Lexan composition (H 0.074) is worse (-0.57 / -0.54 %) but still < 1 %
    assert max(abs(v) for v in _lookup(compound, "Lexan")) < 0.01
    for key, entries in compound["makarova_tab_results"].items():
        for e in entries:
            if key in ("Polystyrene", "Lexan_C16H14O3", "Lucite", "Nylon", "Teflon", "Carbon"):
                assert abs(e["ours_over_theta_hanson_LOOKUP_minus_1"]) < 0.0030


def test_lexan_230_mev_modstack(compound: dict[str, Any]) -> None:
    row = compound["gottschalk2010_modStack"]["row12_Lexan_C16H14O3"]
    assert row["ours_over_thetaH_minus_1"] == pytest.approx(-0.0004, abs=3e-4)
    assert abs(row["E_out_ours_minus_table_MeV"]) < 0.06  # exit energy within 0.04 MeV


def test_compound_rules_are_discriminated(makarova: dict[str, Any]) -> None:
    """The comparison excludes the electron-/mass-weighted Fano alternatives
    (+4.1 % / +1.8 % for polystyrene)."""
    alt = P1C.alt_compound_rules(makarova)
    d = {r: alt["Polystyrene_0.347"][r] - alt["Polystyrene_0.347"]["committed"]
         for r in ("electron_weighted", "mass_weighted")}  # fmt: skip
    assert d["electron_weighted"] == pytest.approx(0.041, abs=3e-3)
    assert d["mass_weighted"] == pytest.approx(0.018, abs=3e-3)
    # thin-slab constant-energy approximation: the committed rule itself is within 1 %
    assert abs(alt["Polystyrene_0.347"]["committed"]) < 0.01
    luc = alt["Lucite_0.366"]
    assert luc["electron_weighted"] - luc["committed"] == pytest.approx(0.046, abs=3e-3)
    assert luc["mass_weighted"] - luc["committed"] == pytest.approx(0.022, abs=3e-3)


def test_water_b_check_against_tbl_b(compound: dict[str, Any]) -> None:
    """Only a two-digit water theta_Hanson exists (arXiv:0908.1413 tbl:B): a b check, not a
    1 % validation. Our b is low by 0.040 +- 0.002 at B = 4.5; Be/Cu/Pb within 0.007."""
    b = compound["gottschalk2010_tblB"]
    assert b["water"]["b_ours_minus_b_G"] == pytest.approx(-0.040, abs=2e-3)
    for k in ("Be", "Cu", "Pb"):
        assert abs(b[k]["b_ours_minus_b_G"]) < 0.007
    eff = P1C.water_b_shift_effect()
    for v in eff.values():
        assert -0.0020 < v["theta_rel_change"] < -0.0010
