"""Packaging test: third-party notices are shipped in the wheel and the sdist."""

from __future__ import annotations

import shutil
import subprocess
import tarfile
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
ACK = "This product includes software developed by Members of the Geant4 Collaboration"


def test_wheel_and_sdist_contain_third_party_notices(tmp_path: Path) -> None:
    uv = shutil.which("uv")
    assert uv is not None, "uv must be on PATH for the packaging test"
    result = subprocess.run(
        [uv, "build", "--wheel", "--sdist", "--out-dir", str(tmp_path)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert result.returncode == 0, result.stderr
    wheels = list(tmp_path.glob("*.whl"))
    sdists = list(tmp_path.glob("*.tar.gz"))
    assert len(wheels) == 1 and len(sdists) == 1
    with zipfile.ZipFile(wheels[0]) as wheel:
        names = wheel.namelist()
        assert "ionmc/THIRD_PARTY_NOTICES.md" in names
        # hatchling keeps the source-relative path below dist-info/licenses
        licenses = [
            n
            for n in names
            if ".dist-info/licenses/" in n and n.endswith("/THIRD_PARTY_NOTICES.md")
        ]
        assert len(licenses) == 1
        text = wheel.read(licenses[0]).decode("utf-8")
        assert ACK in text and "Geant4 Software License" in text
        assert ACK in wheel.read("ionmc/THIRD_PARTY_NOTICES.md").decode("utf-8")
    with tarfile.open(sdists[0]) as sdist:
        assert any(n.endswith("/src/ionmc/THIRD_PARTY_NOTICES.md") for n in sdist.getnames())
