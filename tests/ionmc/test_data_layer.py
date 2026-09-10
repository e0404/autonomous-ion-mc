"""Tests of the external physics-data layer (cache, integrity, parsing)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from ionmc.data import (
    DATASETS,
    MCSQUARE_PSTAR_WATER,
    DatasetNotCached,
    DatasetSpec,
    IntegrityError,
    acquire,
    cache,
    cache_dir,
    is_cached,
    load_path,
    read_manifest,
)
from ionmc.data.__main__ import main as data_cli
from ionmc.data.stopping_tables import (
    TableFormatError,
    hermite_eval,
    load_stopping_table,
    parse_two_column_table,
    pchip_slopes,
    prepare_stopping_table,
    segment_range_increments,
)

SAMPLE_ROWS = "0.0\t129.7\n0.5\t413.2\n1.0\t260.8\n1.5\t195.7\n2.0\t158.6\n"


def _local_spec(tmp_path: Path, content: str = SAMPLE_ROWS, **overrides) -> DatasetSpec:
    version = overrides.get("version", "v1")
    source = tmp_path / f"source-{version}.dat"
    source.write_text(content)
    fields = dict(
        name="test-table",
        version=version,
        url=source.as_uri(),
        sha256=hashlib.sha256(content.encode()).hexdigest(),
        filename="table.dat",
        license="test",
        provenance={"dataset": "synthetic"},
    )
    fields.update(overrides)
    return DatasetSpec(**fields)


def test_cache_dir_precedence(tmp_path, monkeypatch) -> None:
    monkeypatch.delenv(cache.CACHE_ENV_VAR, raising=False)
    monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
    assert cache_dir() == Path.home() / ".cache" / "ionmc"
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    assert cache_dir() == tmp_path / "xdg" / "ionmc"
    monkeypatch.setenv(cache.CACHE_ENV_VAR, str(tmp_path / "env"))
    assert cache_dir() == tmp_path / "env"
    assert cache_dir(tmp_path / "explicit") == tmp_path / "explicit"


def test_spec_validation() -> None:
    with pytest.raises(ValueError):
        DatasetSpec("bad name", "v", "u", "0" * 64, "f", "l", {})
    with pytest.raises(ValueError):
        DatasetSpec("ok", "v", "u", "xyz", "f", "l", {})
    with pytest.raises(ValueError):
        DatasetSpec("ok", "v", "u", "0" * 64, "../f", "l", {})


def test_acquire_then_offline_load(tmp_path) -> None:
    root = tmp_path / "cache"
    spec = _local_spec(tmp_path)
    with pytest.raises(DatasetNotCached) as excinfo:
        load_path(spec, root)
    assert "python -m ionmc.data acquire test-table" in str(excinfo.value)
    assert not is_cached(spec, root)
    path = acquire(spec, root)
    assert path.read_text() == SAMPLE_ROWS
    assert is_cached(spec, root)
    manifest = read_manifest(spec, root)
    assert manifest["sha256"] == spec.sha256
    assert manifest["spec"]["version"] == "v1"
    assert manifest["size_bytes"] == len(SAMPLE_ROWS)
    # offline use returns the same file and verifies its checksum
    assert load_path(spec, root) == path
    # a second acquire is a no-op without network (the URL may be gone)
    (tmp_path / "source-v1.dat").unlink()
    assert acquire(spec, root) == path


def test_acquire_rejects_checksum_mismatch_and_writes_nothing(tmp_path) -> None:
    root = tmp_path / "cache"
    spec = _local_spec(tmp_path, sha256="0" * 64)
    with pytest.raises(IntegrityError):
        acquire(spec, root)
    assert not cache.cached_file_path(spec, root).exists()
    assert not list(root.rglob(".partial-*")) if root.exists() else True


def test_load_detects_tampered_cache(tmp_path) -> None:
    root = tmp_path / "cache"
    spec = _local_spec(tmp_path)
    path = acquire(spec, root)
    path.write_text(SAMPLE_ROWS + "2.5\t134.4\n")
    with pytest.raises(IntegrityError):
        load_path(spec, root)
    assert load_path(spec, root, verify=False) == path


def test_versions_are_cached_side_by_side(tmp_path) -> None:
    root = tmp_path / "cache"
    a = _local_spec(tmp_path, version="v1")
    b_content = SAMPLE_ROWS + "2.5\t134.4\n"
    b = _local_spec(tmp_path, content=b_content, version="v2")
    pa, pb = acquire(a, root), acquire(b, root)
    assert pa != pb and pa.parent.parent == pb.parent.parent
    assert pa.read_text() == SAMPLE_ROWS and pb.read_text() == b_content


def test_cli_status_and_cache_dir(tmp_path, capsys) -> None:
    assert data_cli(["--cache-dir", str(tmp_path), "cache-dir"]) == 0
    assert capsys.readouterr().out.strip() == str(tmp_path)
    assert data_cli(["--cache-dir", str(tmp_path), "status"]) == 0
    out = capsys.readouterr().out
    for name in DATASETS:
        assert f"{name}\t" in out and "missing" in out
    with pytest.raises(SystemExit):
        data_cli(["--cache-dir", str(tmp_path), "acquire", "no-such-dataset"])


def test_registry_specs_are_commit_pinned() -> None:
    for spec in DATASETS.values():
        assert spec.version in spec.url, "URL must be pinned to the immutable version"
        assert "githubusercontent.com" in spec.url
        assert len(spec.sha256) == 64
        assert "Apache-2.0" in spec.license
    assert MCSQUARE_PSTAR_WATER.filename == "PSTAR_Stop_Pow.dat"


# --- parsing and preparation -------------------------------------------------


def test_parse_two_column_table(tmp_path) -> None:
    p = tmp_path / "t.dat"
    p.write_text("# comment\n1.0 2.0\n\n2.0\t3.0\n")
    e, s = parse_two_column_table(p)
    np.testing.assert_array_equal(e, [1.0, 2.0])
    np.testing.assert_array_equal(s, [2.0, 3.0])
    p.write_text("1.0 2.0 3.0\n")
    with pytest.raises(TableFormatError):
        parse_two_column_table(p)
    p.write_text("1.0 x\n")
    with pytest.raises(TableFormatError):
        parse_two_column_table(p)
    p.write_text("")
    with pytest.raises(TableFormatError):
        parse_two_column_table(p)


def test_prepare_drops_zero_energy_row_and_validates(tmp_path) -> None:
    p = tmp_path / "t.dat"
    p.write_text(SAMPLE_ROWS)
    table = load_stopping_table(p, {"dataset": "sample"})
    assert table.size == 4 and table.range_floor_energy_mev == 0.5
    assert table.bisection_steps == 2
    assert table.csda_range_g_per_cm2[0] == 0.0
    assert np.all(np.diff(table.csda_range_g_per_cm2) > 0.0)
    assert np.all(table.slope < 0.0)  # S decreases with E above 0.5 MeV
    with pytest.raises(TableFormatError):
        prepare_stopping_table(np.array([1.0, 1.0]), np.array([1.0, 2.0]))
    with pytest.raises(TableFormatError):
        prepare_stopping_table(np.array([1.0, 2.0]), np.array([1.0, -2.0]))
    with pytest.raises(TableFormatError):
        prepare_stopping_table(np.array([0.0]), np.array([1.0]))


def test_pchip_slopes_reproduce_cubic_data_and_preserve_monotonicity() -> None:
    # Fritsch-Carlson slopes interpolate a monotone quadratic well and are
    # exactly zero where the data are flat, never overshooting.
    x = np.linspace(0.0, 4.0, 9)
    y = x**2 + 1.0
    d = pchip_slopes(x, y)
    assert np.all(d >= 0.0)
    assert d[0] == pytest.approx(0.0, abs=0.06)  # true slope 0 at x = 0
    assert d[-1] == pytest.approx(8.0, rel=0.02)
    flat = pchip_slopes(np.array([0.0, 1.0, 2.0, 3.0]), np.array([1.0, 2.0, 2.0, 3.0]))
    assert flat[1] == 0.0 and flat[2] == 0.0
    # a sign change in secants gives a zero slope at the extremum
    peak = pchip_slopes(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 0.0]))
    assert peak[1] == 0.0


def test_hermite_eval_reproduces_end_values_and_a_cubic() -> None:
    h, y0, y1, d0, d1 = 2.0, 1.0, 5.0, 0.5, -3.0
    assert hermite_eval(np.array(0.0), h, y0, y1, d0, d1) == y0
    assert hermite_eval(np.array(1.0), h, y0, y1, d0, d1) == y1
    # exact for cubics: p(x) = 1 + 2 x - x^2 + 0.5 x^3 on [0, 2]
    p = np.polynomial.Polynomial([1.0, 2.0, -1.0, 0.5])
    dp = p.deriv()
    t = np.linspace(0.0, 1.0, 11)
    got = hermite_eval(t, 2.0, p(0.0), p(2.0), dp(0.0), dp(2.0))
    np.testing.assert_allclose(got, p(2.0 * t), rtol=1e-13)


def test_segment_range_increments_match_fine_trapezoid() -> None:
    x = np.linspace(1.0, 5.0, 9)
    y = 100.0 / x
    d = pchip_slopes(x, y)
    inc = segment_range_increments(x, y, d)  # 4-point Gauss-Legendre
    fine = []
    for k in range(x.size - 1):
        t = np.linspace(0.0, 1.0, 20001)
        s = hermite_eval(t, x[k + 1] - x[k], y[k], y[k + 1], d[k], d[k + 1])
        fine.append(np.trapezoid(1.0 / s, t) * (x[k + 1] - x[k]))
    np.testing.assert_allclose(inc, fine, rtol=1e-6)
    # the 8-point rule is a check of the 4-point one
    np.testing.assert_allclose(inc, segment_range_increments(x, y, d, 8), rtol=1e-6)


@pytest.mark.network
def test_acquire_registered_pstar_table_from_network(tmp_path) -> None:
    """Real acquisition (skipped when the network is unreachable)."""
    import urllib.error

    try:
        path = acquire(MCSQUARE_PSTAR_WATER, tmp_path)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        pytest.skip(f"network unavailable: {exc}")
    table = load_stopping_table(path)
    assert (
        table.size == 800
        and table.energy_mev[0] == 0.5
        and table.energy_mev[-1] == 400.0
    )
    assert float(table.stopping_mev_cm2_per_g[table.energy_mev == 100.0][0]) == 7.289
    manifest = json.loads(
        cache.manifest_path(MCSQUARE_PSTAR_WATER, tmp_path).read_text()
    )
    assert manifest["spec"]["sha256"] == MCSQUARE_PSTAR_WATER.sha256
