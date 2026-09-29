import json

from ionmc import cli
from ionmc.runconfig import config_from_dict


def test_config_from_dict_builds_slab_phantom_and_scoring():
    cfg = config_from_dict(
        {
            "source": {
                "species": "c12",
                "energy_mev_per_u": 290,
                "position_mm": [0, 0, -2],
            },
            "geometry": {
                "type": "box",
                "size_mm": [40, 40, 200],
                "spacing_mm": 2,
                "slabs": [[20, 40, "bone_cortical", None]],
            },
            "scoring": {"spacing_mm": 4, "shift_mm": [1, 0, 0]},
            "histories": 100,
            "batches": 5,
            "physics": {"straggling": False},
        }
    )
    assert cfg.source.species.name == "c12" and cfg.geometry.shape == (20, 20, 100)
    assert [m.name for m in cfg.geometry.materials] == ["water", "bone_cortical"]
    assert (
        cfg.scoring is not None
        and cfg.scoring.shape == (10, 10, 50)
        and cfg.scoring.origin_mm[0] == -19.0
    )
    assert cfg.physics.straggling is False and cfg.batches == 5


def test_cli_run_writes_arrays_metadata_and_summary(tmp_path, capsys):
    spec = {
        "source": {
            "species": "proton",
            "energy_mev_per_u": 100,
            "position_mm": [0, 0, -0.5],
        },
        "geometry": {"size_mm": [20, 20, 90], "spacing_mm": 1},
        "scoring": {"spacing_mm": [20, 20, 1]},
        "histories": 8,
        "batches": 2,
        "seed": 4,
    }
    (tmp_path / "cfg.json").write_text(json.dumps(spec))
    assert (
        cli.main(
            [
                "run",
                str(tmp_path / "cfg.json"),
                "--output",
                str(tmp_path / "out" / "r"),
                "--offline",
            ]
        )
        == 0
    )
    out = json.loads(capsys.readouterr().out)
    assert (tmp_path / "out" / "r.npz").exists() and (
        tmp_path / "out" / "r.json"
    ).exists()
    summary = (tmp_path / "out" / "r.txt").read_text()
    assert "integral depth-energy profile" in summary and "histories: 8" in summary
    assert out["wall_seconds"] > 0
    assert cli.main(["capabilities"]) == 0
    assert "python" in json.loads(capsys.readouterr().out)
