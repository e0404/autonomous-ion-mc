"""Immutable, content-addressed host inputs, never linked to live task files.

The store and snapshots belong to the trusted service and must remain outside
agent-writable mounts. Read-only file modes alone are not a security boundary.
"""

import fcntl
import hashlib
import os
import stat
import tempfile
from pathlib import Path

DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
FILE_FLAGS = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK
CHUNK = 1024 * 1024


def identity(info):
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def digest_stream(stream, outgoing=None, *, size):
    digest = hashlib.sha256()
    count = 0
    while block := stream.read(CHUNK):
        count += len(block)
        if count > size:
            raise ValueError("Input grew during capture")
        digest.update(block)
        if outgoing is not None:
            outgoing.write(block)
    if count != size:
        raise ValueError("Input shrank during capture")
    return digest.hexdigest()


def sync_dir(path):
    fd = os.open(path, DIRECTORY_FLAGS)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


class InputStore:
    def __init__(self, path, guard):
        self.path = Path(path)
        self.guard = guard
        # Ancestors and this directory are protected service storage, not inputs.
        if self.path.absolute() != self.path.resolve():
            raise ValueError("Host input store must not contain symlink components")
        self.path.mkdir(mode=0o700, exist_ok=True)
        info = self.path.stat(follow_symlinks=False)
        if (
            not stat.S_ISDIR(info.st_mode)
            or info.st_uid != os.getuid()
            or info.st_mode & 0o077
        ):
            raise ValueError(
                "Host input store must be a private service-owned directory"
            )
        self.verified = set()

    def __enter__(self):
        self.lock = os.open(
            self.path / "store.lock", os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600
        )
        try:
            if not stat.S_ISREG(os.fstat(self.lock).st_mode):
                raise ValueError("Host input store lock is not a regular file")
            fcntl.flock(self.lock, fcntl.LOCK_EX)
            # This lock is held across every staging write. Leftovers therefore
            # belong to interrupted captures; never remove named content objects.
            for staging in self.path.glob("staging-*"):
                info = staging.stat(follow_symlinks=False)
                if not stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid():
                    raise ValueError("Unsafe interrupted input staging file")
                staging.unlink()
            sync_dir(self.path)
        except BaseException:
            os.close(self.lock)
            raise
        return self

    def __exit__(self, *args):
        os.close(self.lock)

    def verify(self, path, digest, size, mode):
        with os.fdopen(os.open(path, FILE_FLAGS), "rb") as incoming:
            before = os.fstat(incoming.fileno())
            if (
                not stat.S_ISREG(before.st_mode)
                or before.st_size != size
                or stat.S_IMODE(before.st_mode) != mode
                or before.st_uid != os.getuid()
            ):
                raise ValueError(f"Invalid protected input object: {path.name}")
            if path not in self.verified:
                if digest_stream(incoming, size=size) != digest:
                    raise ValueError(f"Corrupt protected input object: {path.name}")
                if identity(before) != identity(os.fstat(incoming.fileno())):
                    raise ValueError(
                        "Protected input object changed during verification"
                    )
                self.verified.add(path)

    def capture(self, source_fd, name, destination):
        fd = os.open(name, FILE_FLAGS, dir_fd=source_fd)
        with os.fdopen(fd, "rb") as incoming:
            before = os.fstat(fd)
            if not stat.S_ISREG(before.st_mode):
                raise ValueError("Snapshot input changed to a special file")
            digest = digest_stream(incoming, size=before.st_size)
            if identity(before) != identity(os.fstat(fd)):
                raise ValueError("Snapshot input changed during hashing")
            mode = 0o444 | (stat.S_IMODE(before.st_mode) & 0o111)
            object_name = f"{digest}-{mode:o}"
            target = self.path / object_name
            # lexists also detects dangling symlinks, which verify must reject.
            if not os.path.lexists(target):
                self.guard.require(before.st_size)
                temp_fd, temp_name = tempfile.mkstemp(prefix="staging-", dir=self.path)
                try:
                    with os.fdopen(temp_fd, "wb") as outgoing:
                        incoming.seek(0)
                        copied = digest_stream(incoming, outgoing, size=before.st_size)
                        if copied != digest or identity(before) != identity(
                            os.fstat(fd)
                        ):
                            raise ValueError("Snapshot input changed during copying")
                        outgoing.flush()
                        os.fchmod(outgoing.fileno(), mode)
                        os.fsync(outgoing.fileno())
                    # No overwrite of an existing object, even on an unexpected race.
                    os.link(temp_name, target)
                    sync_dir(self.path)
                finally:
                    os.unlink(temp_name)
            self.verify(target, digest, before.st_size, mode)
            if identity(before) != identity(
                os.stat(name, dir_fd=source_fd, follow_symlinks=False)
            ):
                raise ValueError("Snapshot input replaced during capture")
            os.link(target, destination, follow_symlinks=False)
        return {
            "sha256": digest,
            "size": before.st_size,
            "mode": mode,
            "object": object_name,
        }

    def capture_tree(self, source_fd, destination, *, relative=".ionmc-cache"):
        """Capture an anchored tree, rejecting links, special files and mutations."""
        before = os.fstat(source_fd)
        destination.mkdir(mode=0o700)
        entries = [{"path": relative, "type": "directory"}]
        for name in sorted(os.listdir(source_fd)):
            info = os.stat(name, dir_fd=source_fd, follow_symlinks=False)
            child_path = f"{relative}/{name}"
            if stat.S_ISDIR(info.st_mode):
                child = os.open(name, DIRECTORY_FLAGS, dir_fd=source_fd)
                try:
                    if identity(info) != identity(os.fstat(child)):
                        raise ValueError("Snapshot input directory replaced")
                    entries.extend(
                        self.capture_tree(
                            child, destination / name, relative=child_path
                        )
                    )
                finally:
                    os.close(child)
            elif stat.S_ISREG(info.st_mode):
                entry = self.capture(source_fd, name, destination / name)
                entries.append({"path": child_path, "type": "file", **entry})
            else:
                raise ValueError(
                    "Snapshot inputs must not contain links or special files"
                )
        if identity(before) != identity(os.fstat(source_fd)):
            raise ValueError("Snapshot input directory changed during capture")
        os.chmod(destination, 0o555)
        sync_dir(destination)
        return entries
