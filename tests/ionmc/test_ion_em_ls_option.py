"""Tests of the Lindhard-Sorensen option ``BetheOptions.ls`` (V3-008, P3 rule triggered).

The src port ``ionmc.physics.lindhard_sorensen`` is compared with the validation reference module
(loaded by file path); the option's effect is checked against the P3 record
``validation/plans/v3-008/p3-output.json`` (SEEN values; tolerances from the P3 record).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from ionmc.config import validate
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import WATER
from ionmc.physics.lindhard_sorensen import delta_l_ls
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, PROTON
from ionmc.physics.stopping import (
    BetheOptions,
    BetheStoppingSource,
    bethe_mass_stopping,
)

_ROOT = Path(__file__).resolve().parents[2]
_PATH = _ROOT / "validation/scripts/ion_em/lindhard_sorensen.py"
_SPEC = importlib.util.spec_from_file_location("lindhard_sorensen_ref_for_option", _PATH)
assert _SPEC is not None and _SPEC.loader is not None
ref = importlib.util.module_from_spec(_SPEC)
sys.modules["lindhard_sorensen_ref_for_option"] = ref
_SPEC.loader.exec_module(ref)

_P3 = json.loads((_ROOT / "validation/plans/v3-008/p3-output.json").read_text())["summary"]
E_GRID = np.geomspace(2.0, 430.0, 60)
IONS = {"12C": CARBON12, "16O": OXYGEN16}


def _s(proj, e, ls, **kw):
    return bethe_mass_stopping(e, proj, WATER, ls=ls, **kw)


def test_src_port_equals_validation_reference() -> None:
    rng = np.random.default_rng(3)
    zs = rng.choice([2, 6, 8], size=20)
    betas = rng.uniform(0.04, 0.8, size=20)
    betas[:3] = (0.0457, 0.2, 0.75)  # near the 1 MeV/u floor and relativistic end
    got = np.array([delta_l_ls(z, b)[0] for z, b in zip(zs, betas, strict=True)])
    want = np.array(
        [ref.delta_l_ls(ref.ALPHA_FS * z, b).delta_l for z, b in zip(zs, betas, strict=True)]
    )
    assert np.max(np.abs(got - want)) <= 1e-10
    # batched call with per-point (non-integer) charge equals point-wise calls
    q = np.array([5.2, 7.9, 1.5])
    b = np.array([0.1, 0.3, 0.6])
    batch = delta_l_ls(q, b)
    one = np.array([delta_l_ls(qq, bb)[0] for qq, bb in zip(q, b, strict=True)])
    assert np.max(np.abs(batch - one)) <= 1e-12


@pytest.mark.parametrize("proj", [ALPHA, CARBON12, OXYGEN16, PROTON])
def test_off_is_the_current_model_bit_for_bit(proj) -> None:
    default = bethe_mass_stopping(E_GRID, proj, WATER)
    assert np.array_equal(_s(proj, E_GRID, "off"), default)
    assert BetheOptions().ls == "off"
    t_def = BetheStoppingSource().table(WATER, proj)
    t_off = BetheStoppingSource(options=BetheOptions(ls="off")).table(WATER, proj)
    assert np.array_equal(t_def.s_el_mass, t_off.s_el_mass)
    assert t_def.metadata == t_off.metadata  # default metadata (and table hash) unchanged


@pytest.mark.parametrize("ion", ["12C", "16O"])
def test_eff_differs_from_off_by_p3_bound_and_is_positive(ion: str) -> None:
    proj = IONS[ion]
    off = _s(proj, E_GRID, "off")
    eff = _s(proj, E_GRID, "eff")
    rel = eff / off - 1.0
    assert np.all(rel > 0.0)  # ionmc S increases
    assert np.max(rel) <= 0.0011  # P3 record: <= 0.0965 % (O), 0.0598 % (C); tolerance 0.11 %


@pytest.mark.parametrize("ion", ["12C", "16O"])
def test_bare_mixed_at_2_mev_per_u_and_coincides_with_eff_above_30(ion: str) -> None:
    proj = IONS[ion]
    e = np.array([2.0])
    off = _s(proj, e, "off")
    mixed = _s(proj, e, "bare_mixed")
    pct = 100.0 * (mixed / off - 1.0)[0]
    expected = _P3[ion]["mixed_charge_convention"]["r_over_S_pct_at_2"]
    assert pct == pytest.approx(expected, abs=0.1)
    assert pct < -1.0  # sign: ionmc S decreases in the mixed convention
    hi = np.geomspace(30.0, 430.0, 20)
    assert np.max(np.abs(_s(proj, hi, "bare_mixed") / _s(proj, hi, "eff") - 1.0)) <= 1e-4


@pytest.mark.parametrize("ls", ["off", "eff", "bare_mixed"])
def test_protons_are_bit_identical_for_every_value(ls: str) -> None:
    base = bethe_mass_stopping(E_GRID, PROTON, WATER)
    assert np.array_equal(_s(PROTON, E_GRID, ls), base)
    t = BetheStoppingSource(options=BetheOptions(ls=ls)).table(WATER, PROTON)
    assert np.array_equal(t.s_el_mass, BetheStoppingSource().table(WATER, PROTON).s_el_mass)


@pytest.mark.parametrize("bad", ["on", "", "EFF", None, 1, True])
def test_invalid_value_is_rejected(bad) -> None:
    with pytest.raises(ValueError):
        BetheOptions(ls=bad)
    with pytest.raises(ValueError):
        bethe_mass_stopping(10.0, CARBON12, WATER, ls=bad)


def test_ls_needs_bloch_and_mott_for_ions() -> None:
    with pytest.raises(ValueError):
        bethe_mass_stopping(10.0, CARBON12, WATER, ls="eff", bloch=False)
    bethe_mass_stopping(10.0, PROTON, WATER, ls="eff", bloch=False)  # ignored for protons


def test_config_validation_rejects_corrupted_option(make_config) -> None:
    src = BetheStoppingSource(options=BetheOptions(ls="eff"))
    object.__setattr__(src.options, "ls", "bogus")  # bypass the frozen dataclass
    cfg = make_config()
    with pytest.raises(UnsupportedCombinationError):
        replace(cfg, physics=replace(cfg.physics, stopping=src))


def test_option_is_persisted_in_the_effective_configuration(make_config) -> None:
    """A proton run is bit-identical whatever the (ignored) LS request: same table hash, no LS in
    the table metadata; the configuration states requested and effective mode (new tests in
    ``test_ion_em_mapping_src.py``)."""
    cfg = make_config()
    default = validate(cfg)
    assert "ls" not in default.summary()["tables"]["materials"][0]["metadata"]["options"]
    src = BetheStoppingSource(options=BetheOptions(ls="eff"))
    eff = validate(replace(cfg, physics=replace(cfg.physics, stopping=src)))
    meta = eff.summary()["tables"]["materials"][0]["metadata"]
    assert "ls" not in meta["options"]
    assert eff.tables.sha256 == default.tables.sha256  # protons: the request is not an identity
    phys = eff.summary()["physics"]
    assert (phys["ls_requested"], phys["ls_effective"]) == ("eff", "off")


def test_ion_table_build_uses_ls_in_table() -> None:
    t_off = BetheStoppingSource(e_max_per_u=20.0, points_per_decade=20).table(WATER, CARBON12)
    src = BetheStoppingSource(
        options=BetheOptions(ls="eff"), e_max_per_u=20.0, points_per_decade=20
    )
    t_eff = src.table(WATER, CARBON12)
    assert t_eff.metadata["options"]["ls"] == "eff"
    assert np.all(t_eff.s_el_mass > t_off.s_el_mass)
    assert np.all(t_eff.range_mm < t_off.range_mm)
