"""Tests for the external data layer: cache, registry, acquisition, parsers (offline)."""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pytest

from ionmc.cli import main
from ionmc.data import acquire, cache
from ionmc.data.icru90 import parse_icru90_source
from ionmc.data.nist_star import parse_star_text
from ionmc.data.registry import DATASETS, Dataset

# Excerpt (header and first 12 data rows) of the NIST PSTAR liquid-water table as returned
# for the request in DATASETS["nist-pstar-water-2005"]; NIST SRD 124, attribution:
# Berger, Coursey, Zucker, Chang (2005), https://doi.org/10.18434/T4NC7P.
PSTAR_EXCERPT = """PSTAR: Stopping Powers and Range Tables for Protons

WATER, LIQUID

Kinetic   Electron. Nuclear   Total     CSDA      Projected Detour
Energy    Stp. Pow. Stp. Pow. Stp. Pow. Range     Range     Factor
MeV       MeV cm2/g MeV cm2/g MeV cm2/g g/cm2     g/cm2

1.000E-03 1.337E+02 4.315E+01 1.769E+02 6.319E-06 2.878E-06 0.4555
1.500E-03 1.638E+02 3.460E+01 1.984E+02 8.969E-06 4.400E-06 0.4906
2.000E-03 1.891E+02 2.927E+01 2.184E+02 1.137E-05 5.909E-06 0.5197
2.500E-03 2.114E+02 2.557E+01 2.370E+02 1.357E-05 7.380E-06 0.5440
3.000E-03 2.316E+02 2.281E+01 2.544E+02 1.560E-05 8.811E-06 0.5647
4.000E-03 2.675E+02 1.894E+01 2.864E+02 1.930E-05 1.155E-05 0.5986
5.000E-03 2.990E+02 1.631E+01 3.153E+02 2.262E-05 1.415E-05 0.6254
6.000E-03 3.276E+02 1.439E+01 3.420E+02 2.567E-05 1.661E-05 0.6473
7.000E-03 3.538E+02 1.292E+01 3.667E+02 2.849E-05 1.896E-05 0.6656
8.000E-03 3.782E+02 1.175E+01 3.900E+02 3.113E-05 2.121E-05 0.6813
9.000E-03 4.012E+02 1.080E+01 4.120E+02 3.363E-05 2.337E-05 0.6950
1.000E-02 4.229E+02 1.000E+01 4.329E+02 3.599E-05 2.545E-05 0.7070
"""

# Excerpt: first five entries of each water-relevant array of G4ICRU90StoppingData.cc
# (Geant4 v11.4.2, Geant4 Software License) in the original syntax. Index 1 is G4_WATER.
ICRU90_EXCERPT = """
  static const G4float T0_proton[57] = {  0.0010f, 0.001500f, 0.0020f, 0.0030f, 0.0040f };
  static const G4float T0_alpha[49] = {  0.0010f, 0.001500f, 0.0020f, 0.0030f, 0.0040f };
  static const G4float e0_proton[57] = {  119.70f, 146.70f, 169.30f, 207.40f, 239.50f  };
  static const G4float e0_alpha[49] = {  87.50f, 108.60f, 126.70f, 157.30f, 183.50f  };
  static const G4float e1_proton[57] = {  133.70f, 163.80f, 189.10f, 231.60f, 267.50f  };
  static const G4float e1_alpha[49] = {  98.910f, 122.80f, 143.10f, 177.60f, 206.90f  };
  static const G4float e2_proton[57] = {  118.50f, 145.10f, 167.60f, 205.30f, 237.00f  };
  static const G4float e2_alpha[49] = {  192.30f, 228.90f, 259.00f, 308.30f, 348.90f  };
"""

PAYLOAD = b"test dataset payload\n"
PAYLOAD_SHA = cache.sha256_bytes(PAYLOAD)


def _fake_dataset(sha: str = PAYLOAD_SHA, size: int = len(PAYLOAD)) -> Dataset:
    return Dataset(
        id="fake",
        version="1",
        url="https://example.invalid/data",
        method="GET",
        post_body=None,
        sha256=sha,
        bytes=size,
        license="test",
        citation="test",
        parser="none",
        description="fake",
    )


def test_registry_has_exactly_the_three_datasets() -> None:
    assert set(DATASETS) == {
        "nist-pstar-water-2005",
        "nist-astar-water-2005",
        "geant4-icru90-stopping-11.4.2",
    }
    for ds in DATASETS.values():
        assert ds.url.startswith("https://")
        assert len(ds.sha256) == 64
        assert ds.method in {"GET", "POST"}
        assert (ds.post_body is not None) == (ds.method == "POST")
    assert DATASETS["nist-pstar-water-2005"].bytes == 9308
    assert "prog=ASTAR" in (DATASETS["nist-astar-water-2005"].post_body or "")


def test_cache_dir_precedence(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.delenv("IONMC_CACHE_DIR", raising=False)
    assert cache.resolve_cache_dir() == Path.home() / ".cache" / "ionmc"
    monkeypatch.setenv("IONMC_CACHE_DIR", str(tmp_path / "env"))
    assert cache.resolve_cache_dir() == tmp_path / "env"
    assert cache.resolve_cache_dir(tmp_path / "explicit") == tmp_path / "explicit"


def test_fetch_stores_object_and_manifest(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setitem(DATASETS, "fake", _fake_dataset())
    monkeypatch.setattr(acquire, "_download", lambda ds: PAYLOAD)
    path = acquire.fetch("fake", tmp_path)
    assert path == tmp_path / "objects" / PAYLOAD_SHA
    assert path.read_bytes() == PAYLOAD
    manifest = json.loads((tmp_path / "manifests" / "fake.json").read_text())
    assert set(cache.MANIFEST_FIELDS) <= set(manifest)
    assert manifest["sha256"] == PAYLOAD_SHA and manifest["bytes"] == len(PAYLOAD)
    assert manifest["retrieved_at"].endswith("+00:00")
    assert cache.verify("fake", tmp_path) == path


def test_cached_fetch_needs_no_network(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setitem(DATASETS, "fake", _fake_dataset())
    monkeypatch.setattr(acquire, "_download", lambda ds: PAYLOAD)
    acquire.fetch("fake", tmp_path)

    def boom(ds: Dataset) -> bytes:
        raise AssertionError("network used")

    monkeypatch.setattr(acquire, "_download", boom)
    assert acquire.fetch("fake", tmp_path, offline=True).read_bytes() == PAYLOAD


def test_hash_mismatch_is_rejected_and_nothing_stored(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setitem(DATASETS, "fake", _fake_dataset())
    monkeypatch.setattr(acquire, "_download", lambda ds: b"tampered")
    with pytest.raises(cache.IntegrityError):
        acquire.fetch("fake", tmp_path)
    assert not (tmp_path / "objects").exists()
    assert not (tmp_path / "manifests").exists()


def test_offline_raises_when_not_cached(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setitem(DATASETS, "fake", _fake_dataset())
    with pytest.raises(acquire.OfflineError):
        acquire.fetch("fake", tmp_path, offline=True)
    monkeypatch.setenv("IONMC_OFFLINE", "1")
    with pytest.raises(acquire.OfflineError):
        acquire.fetch("fake", tmp_path)


def test_unknown_dataset(tmp_path: Path) -> None:
    with pytest.raises(KeyError):
        acquire.fetch("nope", tmp_path)


def test_verify_detects_corruption(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setitem(DATASETS, "fake", _fake_dataset())
    monkeypatch.setattr(acquire, "_download", lambda ds: PAYLOAD)
    path = acquire.fetch("fake", tmp_path)
    path.write_bytes(b"corrupted")
    with pytest.raises(cache.IntegrityError):
        cache.verify("fake", tmp_path)
    with pytest.raises(cache.IntegrityError):
        acquire.fetch("fake", tmp_path)


def test_https_only_redirects() -> None:
    handler = acquire._HttpsOnlyRedirect()
    req = urllib.request.Request("https://example.invalid/a")
    with pytest.raises(urllib.error.URLError):
        handler.redirect_request(req, None, 302, "Found", {}, "http://example.invalid/b")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        acquire._download(
            _fake_dataset().__class__(**{**_fake_dataset().__dict__, "url": "http://x"})
        )


def test_cli_data_commands(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setitem(DATASETS, "fake", _fake_dataset())
    monkeypatch.setattr(acquire, "_download", lambda ds: PAYLOAD)
    cdir = str(tmp_path)
    assert main(["data", "list", "--cache-dir", cdir]) == 0
    assert "fake\tmissing" in capsys.readouterr().out
    assert main(["data", "path", "fake", "--cache-dir", cdir]) == 1
    assert main(["data", "fetch", "fake", "--cache-dir", cdir]) == 0
    assert main(["data", "verify", "fake", "--cache-dir", cdir]) == 0
    assert main(["data", "path", "fake", "--cache-dir", cdir]) == 0
    assert main(["data", "fetch", "fake", "--cache-dir", str(tmp_path / "other"), "--offline"]) == 1
    assert main(["data"]) == 2


def test_star_parser_excerpt() -> None:
    table = parse_star_text(PSTAR_EXCERPT)
    assert table.program == "PSTAR"
    assert table.material.startswith("WATER, LIQUID")
    assert table.energy_mev.shape == (12,)
    assert table.energy_mev[0] == 1.0e-3 and table.energy_mev[-1] == 1.0e-2
    assert table.s_electronic[0] == 133.7 and table.csda_range[-1] == 3.599e-5
    assert table.detour[2] == 0.5197
    assert table.energy_mev.dtype == np.float64


def test_star_parser_rejects_bad_input() -> None:
    lines = PSTAR_EXCERPT.splitlines()
    swapped = "\n".join(lines[:8] + [lines[9], lines[8]] + lines[10:])
    with pytest.raises(ValueError, match="increasing"):
        parse_star_text(swapped)
    short = "\n".join(lines[:8] + ["1.0E-03 1.0 2.0"])
    with pytest.raises(ValueError, match="columns"):
        parse_star_text(short)
    with pytest.raises(ValueError):
        parse_star_text("\n".join(lines[:7]))


def test_icru90_parser_excerpt() -> None:
    data = parse_icru90_source(ICRU90_EXCERPT)
    assert data.proton_energy_mev.tolist() == [0.001, 0.0015, 0.002, 0.003, 0.004]
    assert data.proton_stopping.tolist() == [133.7, 163.8, 189.1, 231.6, 267.5]
    assert data.alpha_stopping[0] == 98.91
    graphite = parse_icru90_source(ICRU90_EXCERPT, "G4_GRAPHITE")
    assert graphite.proton_stopping[0] == 118.5
    with pytest.raises(ValueError):
        parse_icru90_source(ICRU90_EXCERPT.replace("T0_alpha", "X"))
