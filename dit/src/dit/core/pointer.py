from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

import tomli_w

from dit.core.repo import POINTER_SUFFIX


@dataclass(frozen=True)
class Pointer:
    path: str
    hash: str
    size: int

    @property
    def pointer_relpath(self) -> str:
        return self.path + POINTER_SUFFIX

    def to_dict(self) -> dict[str, str | int]:
        return {"path": self.path, "hash": self.hash, "size": self.size}


def pointer_path_for(data_path: Path) -> Path:
    return Path(str(data_path) + POINTER_SUFFIX)


def data_path_for_pointer(pointer_path: Path) -> Path:
    text = str(pointer_path)
    if not text.endswith(POINTER_SUFFIX):
        raise ValueError(f"not a pointer path: {pointer_path}")
    return Path(text[: -len(POINTER_SUFFIX)])


def read_pointer(path: Path) -> Pointer:
    try:
        with path.open("rb") as f:
            data = tomllib.load(f)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"pointer not found: {path}") from exc
    except tomllib.TOMLDecodeError as exc:
        raise ValueError(f"invalid pointer TOML {path}: {exc}") from exc
    try:
        return Pointer(path=str(data["path"]), hash=str(data["hash"]), size=int(data["size"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid pointer fields in {path}: {exc}") from exc


def write_pointer(repo_root: Path, pointer: Pointer) -> Path:
    out = repo_root / pointer.pointer_relpath
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as f:
        tomli_w.dump(pointer.to_dict(), f)
    return out


def ensure_gitignore_entry(data_file: Path, rel_name: str) -> None:
    gitignore = data_file.parent / ".gitignore"
    line = rel_name
    if gitignore.is_file():
        existing = gitignore.read_text(encoding="utf-8").splitlines()
        if line in existing or f"/{line}" in existing:
            return
        with gitignore.open("a", encoding="utf-8") as f:
            f.write(f"\n{line}\n")
        return
    gitignore.write_text(f"{line}\n", encoding="utf-8")
