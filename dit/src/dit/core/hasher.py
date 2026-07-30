from __future__ import annotations

from pathlib import Path

import blake3

HASH_PREFIX = "blake3:"
HASH_CHUNK_SIZE = 1024 * 1024


def hash_file(path: Path) -> str:
    hasher = blake3.blake3()
    try:
        with path.open("rb") as f:
            while True:
                chunk = f.read(HASH_CHUNK_SIZE)
                if not chunk:
                    break
                hasher.update(chunk)
    except OSError as exc:
        raise OSError(f"failed to hash {path}: {exc}") from exc
    return HASH_PREFIX + hasher.hexdigest()


def strip_hash_prefix(value: str) -> str:
    if value.startswith(HASH_PREFIX):
        return value[len(HASH_PREFIX) :]
    return value
