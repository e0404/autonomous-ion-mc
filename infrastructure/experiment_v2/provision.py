"""Maintainer-only registration of dedicated reference runtime trees."""

import argparse
import json
from pathlib import Path

from infrastructure.experiment_v2.common import STATE, write_json
from infrastructure.experiment_v2.reference import runtime_inventory


def provision(engine, root, version, source, *, geant4_data=None, registry=None):
    root = Path(root).resolve()
    dest = "/opt/reference/" + engine
    if engine == "fred":
        dest += "/curr"
    mounts = [
        {"source": str(root), "destination": dest, "sha256": runtime_inventory(root)}
    ]
    env = {}
    if engine == "topas":
        if not geant4_data:
            raise ValueError("TOPAS requires a dedicated Geant4 data directory")
        data = Path(geant4_data).resolve()
        mounts.append(
            {
                "source": str(data),
                "destination": "/opt/reference/g4data",
                "sha256": runtime_inventory(data),
            }
        )
        executable = dest + "/bin/topas"
        env = {
            "TOPAS_G4_DATA_DIR": "/opt/reference/g4data",
            "LD_LIBRARY_PATH": dest + "/lib",
            "QT_QPA_PLATFORM": "offscreen",
        }
    elif engine == "mcsquare":
        executable = dest + "/MCsquare_linux"
        env = {"MCsquare_Materials_Dir": dest + "/Materials"}
    elif engine == "fred":
        executable = dest + "/scripts/fred.py"
        env = {
            "FREDDIR": "/opt/reference/fred",
            "LD_LIBRARY_PATH": dest + ":/usr/lib/wsl/lib",
            "FRED_CACHE_DIR": "/tmp/fred-cache",
            "FRED_GPU_AVAILABLE": "0",
            "FRED_MAX_THREADS": "1",
            "FRED_FORCE_RESOURCE_SCAN": "0",
        }
    else:
        raise ValueError("Unknown engine")
    local_exe = root / executable.removeprefix(dest + "/")
    if not local_exe.is_file():
        raise ValueError(f"Executable missing: {local_exe}")
    path = Path(registry or STATE / "engines.json")
    registry_data = (
        json.loads(path.read_text())
        if path.exists()
        else {"schema_version": 2, "engines": {}}
    )
    registry_data["engines"][engine] = {
        "engine": engine,
        "version": version,
        "source": source,
        "executable": executable,
        "mounts": mounts,
        "environment": env,
    }
    write_json(path, registry_data)
    path.chmod(0o600)
    return registry_data["engines"][engine]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("engine", choices=["topas", "mcsquare", "fred"])
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--version", required=True)
    p.add_argument("--source", required=True)
    p.add_argument("--geant4-data", type=Path)
    p.add_argument("--registry", type=Path)
    a = p.parse_args()
    print(
        json.dumps(
            provision(
                a.engine,
                a.root,
                a.version,
                a.source,
                geant4_data=a.geant4_data,
                registry=a.registry,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
