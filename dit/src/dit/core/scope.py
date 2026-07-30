from __future__ import annotations

import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

import tomli_w
from natsort import natsorted

from dit.core.repo import Repo


class Scope:
    def __init__(self, repo: Repo) -> None:
        self.repo = repo
        self._directories: set[str] = set()
        self._load()

    @property
    def directories(self) -> set[str]:
        return set(self._directories)

    def add(self, dir_path: Path) -> str:
        if not dir_path.is_dir():
            raise NotADirectoryError(f"{dir_path} is not a directory")
        rel = self.repo.rel(dir_path)
        self._directories.add(rel)
        self._save()
        return rel

    def remove(self, dir_path: Path | str) -> str:
        rel = dir_path if isinstance(dir_path, str) else self.repo.rel(Path(dir_path))
        if rel not in self._directories:
            raise KeyError(f"{rel} is not in scope: {natsorted(self._directories)}")
        self._directories.remove(rel)
        self._save()
        return rel

    def list(self) -> list[str]:
        return natsorted(self._directories)

    def contains(self, rel_path: str) -> bool:
        if not self._directories:
            return False
        for directory in self._directories:
            if rel_path == directory or rel_path.startswith(directory.rstrip("/") + "/"):
                return True
        return False

    def _load(self) -> None:
        path = self.repo.scope_toml
        if not path.is_file():
            self._directories = set()
            return
        with path.open("rb") as f:
            data = tomllib.load(f)
        directories = data.get("directories") or []
        if not isinstance(directories, list):
            raise ValueError(f"invalid scope file: {path}")
        self._directories = {str(d) for d in directories}

    def _save(self) -> None:
        self.repo.dit_dir.mkdir(parents=True, exist_ok=True)
        with self.repo.scope_toml.open("wb") as f:
            tomli_w.dump({"directories": natsorted(self._directories)}, f)
