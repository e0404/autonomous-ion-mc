"""Admission checks for protected host-run storage, including WSL backing disks."""

import json
import os
import platform
import shutil
from pathlib import Path

GIB = 1024**3


class StorageGuard:
    """Keep headroom at admission; this is not a reservation or runtime quota."""

    def __init__(self, state):
        self.state = Path(state)
        config_path = self.state / "host-storage.json"
        config = json.loads(config_path.read_text()) if config_path.exists() else {}
        reserve = config.get("reserve_bytes", 10 * GIB)
        if type(reserve) is not int or reserve <= 0:
            raise ValueError(
                "host-storage.json reserve_bytes must be a positive integer"
            )
        self.reserve = reserve
        volume = os.environ.get("IONMC_HOST_BACKING_VOLUME") or config.get(
            "backing_volume"
        )
        if not volume and "microsoft" in platform.release().lower():
            raise RuntimeError(
                "Host validation paused: WSL backing-volume check is unconfigured. "
                "Set backing_volume in protected state/host-storage.json (for example "
                '"/mnt/c", only if that drive holds this distro VHD), or set '
                "IONMC_HOST_BACKING_VOLUME in the MCP service environment. "
                "Request operator intervention; do not retry expensive validations."
            )
        self.paths = [self.state]
        if volume:
            path = Path(volume)
            if not path.is_absolute() or not path.is_dir():
                raise ValueError(
                    "Host backing_volume must be an existing absolute directory"
                )
            self.paths.append(path)

    def require(self, additional_bytes=0):
        for path in self.paths:
            free = shutil.disk_usage(path).free
            required = self.reserve + additional_bytes
            if free < required:
                raise RuntimeError(
                    f"Host validation paused: insufficient free space at {path}: "
                    f"{free} bytes free; need {required} bytes "
                    f"({additional_bytes} new bytes plus {self.reserve} reserve). "
                    "Request operator storage intervention; do not retry expensive "
                    "validations or delete archived evidence."
                )
