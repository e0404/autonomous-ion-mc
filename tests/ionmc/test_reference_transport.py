"""Reference-backend transport validation (electromagnetic physics only).

Small cases run in CI; larger statistical cases are marked ``local`` and run
under the exact-SHA validation gate. Ranges are compared with the ICRU 90
CSDA range of water (same mean-excitation-energy lineage, decision 0040).
"""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pytest

from ionmc.config import (
    PhysicsConfig,
    SimulationConfig,
    UnsupportedConfigurationError,
    capabilities,
)
from ionmc.data import icru90
from ionmc.geometry import homogeneous_box
from ionmc.physics.scattering import (
    fermi_eyges_lateral_sigma_mm,
    scattering_length_g_cm2,
)
from ionmc.results import SimulationResult
from ionmc.scoring import ScoringGrid
from ionmc.simulation import build_tables, run
from ionmc.sources import PencilBeam
from ionmc.species import get_species


def r80_mm(result: SimulationResult) -> float:
    z, prof = result.depth_dose(2)
    peak = prof.max()
    i = int(np.where(prof >= 0.8 * peak)[0][-1])
    if i + 1 < len(prof) and prof[i + 1] != prof[i]:
        return float(
            z[i] + (0.8 * peak - prof[i]) / (prof[i + 1] - prof[i]) * (z[i + 1] - z[i])
        )
    return float(z[i])


def water_box(size=(40.0, 40.0, 100.0), spacing=1.0):
    return homogeneous_box(size, spacing, "water")


@pytest.fixture(scope="module")
def proton_tables():
    geo = water_box()
    cfg = SimulationConfig(PencilBeam("proton", 100.0), geo, histories=1, batches=1)
    return build_tables(cfg, offline=True)


def test_energy_is_conserved_and_accounted(proton_tables):
    geo = water_box()
    cfg = SimulationConfig(
        PencilBeam("proton", 100.0, (0.0, 0.0, -0.5)),
        geo,
        histories=40,
        batches=4,
        seed=1,
        scoring=ScoringGrid.coarse(geo, (40.0, 40.0, 1.0)),
    )
    res = run(cfg, tables=proton_tables, offline=True)
    acc = res.metadata["energy_accounting_per_primary"]
    total = (
        acc["deposited_continuous_mev"]
        + acc["deposited_cutoff_mev"]
        + acc["escaped_mev"]
        + acc["truncated_max_steps_mev"]
    )
    assert total == pytest.approx(acc["initial_mev"], rel=1e-9)
    assert res.energy_mean.sum() == pytest.approx(
        acc["deposited_continuous_mev"] + acc["deposited_cutoff_mev"], rel=1e-9
    )
    assert acc["escaped_mev"] == 0.0
    assert acc["steps_per_history"] > 100
    assert np.isfinite(res.dose_mean).all() and res.dose_mean.max() > 0


def test_range_matches_icru90_csda_and_is_step_size_independent(proton_tables):
    geo = water_box()
    ref = float(icru90.water_table("proton").csda_range_at(100.0)) * 10.0
    results = {}
    for label, ph in (
        ("default", PhysicsConfig()),
        ("fine", PhysicsConfig(max_step_mm=0.5, energy_step_fraction=0.05)),
    ):
        cfg = SimulationConfig(
            PencilBeam("proton", 100.0, (0.0, 0.0, -0.5)),
            geo,
            histories=120,
            batches=4,
            seed=5,
            physics=ph,
            scoring=ScoringGrid.coarse(geo, (40.0, 40.0, 1.0)),
        )
        results[label] = r80_mm(run(cfg, tables=proton_tables, offline=True))
    for label, r80 in results.items():
        assert abs(r80 - ref) < max(0.5, 0.005 * ref), (label, r80, ref)
    assert abs(results["default"] - results["fine"]) < 0.3


def test_scoring_grid_shift_and_coarsening_conserve_energy(proton_tables):
    geo = water_box()
    base = ScoringGrid.coarse(geo, 2.0)
    totals = []
    for grid in (base, base.shifted((1.0, 0.0, 0.5)), ScoringGrid.coarse(geo, 4.0)):
        cfg = SimulationConfig(
            PencilBeam("proton", 100.0, (0.0, 0.0, -0.5)),
            geo,
            histories=30,
            batches=3,
            seed=9,
            scoring=grid,
        )
        res = run(cfg, tables=proton_tables, offline=True)
        totals.append(res.energy_mean.sum())
    # the shifted grid loses only the energy deposited in the half-voxel strip pushed outside (none here: beam on axis)
    assert totals[0] == pytest.approx(totals[1], rel=0.02)
    assert totals[0] == pytest.approx(totals[2], rel=1e-9)


def test_capability_contract_fails_closed():
    geo = water_box()
    beam = PencilBeam("proton", 100.0)
    with pytest.raises(UnsupportedConfigurationError, match="nuclear"):
        run(
            SimulationConfig(
                beam, geo, histories=2, batches=1, physics=PhysicsConfig(nuclear=True)
            ),
            offline=True,
        )
    with pytest.raises(UnsupportedConfigurationError, match="backend"):
        run(
            SimulationConfig(beam, geo, histories=2, batches=1, backend="opencl"),
            offline=True,
        )
    with pytest.raises(UnsupportedConfigurationError, match="precision"):
        run(
            SimulationConfig(beam, geo, histories=2, batches=1, precision="float32"),
            offline=True,
        )
    with pytest.raises(UnsupportedConfigurationError, match="scorer"):
        run(
            SimulationConfig(beam, geo, histories=2, batches=1, scorers=("let",)),
            offline=True,
        )
    with pytest.raises(UnsupportedConfigurationError, match="model"):
        PhysicsConfig(mcs_model="highland").validate()
    assert "python" in capabilities()


def test_result_reopens_in_a_separate_process(tmp_path, proton_tables):
    geo = water_box()
    cfg = SimulationConfig(
        PencilBeam("proton", 100.0, (0.0, 0.0, -0.5)),
        geo,
        histories=12,
        batches=3,
        seed=2,
        scoring=ScoringGrid.coarse(geo, 4.0),
    )
    res = run(cfg, tables=proton_tables, offline=True)
    npz, js = res.save(tmp_path / "result")
    code = (
        "import json,sys,numpy as np;"
        f"m=json.load(open({str(js)!r}));d=np.load({str(npz)!r});"
        "print(json.dumps({'sha':m['code']['git_sha'],'dirty':m['code']['git_dirty'],'hist':m['effective']['histories'],"
        "'batches':m['effective']['batches'],'seed':m['effective']['seed'],'phys':m['effective']['physics']['straggling'],"
        "'units':m['arrays']['dose_mean']['units'],'shape':list(d['dose_mean'].shape),'origin':m['scoring']['origin_mm'],"
        "'tables':list(m['effective']['tables']['tables'].keys()),'acc':m['energy_accounting_per_primary']['initial_mev']}))"
    )
    out = json.loads(
        subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, check=True
        ).stdout
    )
    assert (
        out["hist"] == 12
        and out["batches"] == 3
        and out["seed"] == 2
        and out["phys"] is True
    )
    assert out["units"] == "Gy per primary" and out["shape"] == [10, 10, 25]
    assert out["tables"] == ["proton/water"] and out["acc"] == 100.0
    assert out["sha"] is None or len(out["sha"]) == 40
    reloaded = SimulationResult.load(tmp_path / "result")
    assert np.array_equal(reloaded.dose_mean, res.dose_mean)


@pytest.mark.local
def test_lateral_spread_at_30mm_matches_fermi_eyges(proton_tables):
    """Gaussian core of MCS: sigma_x in a 1 mm slab at 30 mm depth vs Fermi–Eyges with T_dM."""
    geo = water_box((60.0, 60.0, 40.0), 1.0)
    slab = ScoringGrid((-30.0, -30.0, 30.0), (0.5, 60.0, 1.0), (120, 1, 1))
    cfg = SimulationConfig(
        PencilBeam("proton", 100.0, (0.0, 0.0, -0.5)),
        geo,
        histories=600,
        batches=4,
        seed=11,
        physics=PhysicsConfig(straggling=False),
        scoring=slab,
    )
    res = run(cfg, tables=proton_tables, offline=True)
    x = slab.centers_mm(0)
    w = res.energy_mean[:, 0, 0]
    sigma_mc = float(np.sqrt(np.sum(w * x * x) / np.sum(w)))
    sp = get_species("proton")
    depth = np.linspace(0.0, 30.5, 306)  # mm, to the slab centre
    t = np.array([proton_tables_energy(proton_tables, 100.0, d / 10.0) for d in depth])
    sigma_an = fermi_eyges_lateral_sigma_mm(
        sp, t, depth / 10.0, scattering_length_g_cm2(geo.materials[0]), 1.0
    )
    assert sigma_mc == pytest.approx(sigma_an, rel=0.10), (sigma_mc, sigma_an)


def proton_tables_energy(tables, t0: float, path_g_cm2: float) -> float:
    from ionmc.transport.shared import physics

    py = physics("python")
    return py.energy_after_path(
        tables.log_range,
        0,
        0,
        t0,
        path_g_cm2,
        tables.t_log_min,
        tables.t_inv_dlog,
        tables.t_grid.size,
    )


@pytest.mark.local
def test_carbon_range_in_water(proton_tables):
    geo = homogeneous_box((40.0, 40.0, 200.0), 1.0, "water")
    cfg = SimulationConfig(
        PencilBeam("c12", 290.0, (0.0, 0.0, -0.5)),
        geo,
        histories=40,
        batches=4,
        seed=21,
        scoring=ScoringGrid.coarse(geo, (40.0, 40.0, 1.0)),
    )
    res = run(cfg, offline=True)
    ref = float(icru90.water_table("c12").csda_range_at(290.0)) * 10.0
    assert abs(r80_mm(res) - ref) < max(0.7, 0.007 * ref)
