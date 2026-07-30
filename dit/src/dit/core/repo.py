from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

DIT_DIR_NAME = ".dit"
DIT_TOML_NAME = "dit.toml"
POINTER_SUFFIX = ".dit"
INDEX_DB_NAME = "index.db"
SCOPE_TOML_NAME = "scope.toml"


@dataclass(frozen=True)
class Repo:
    root: Path

    @property
    def dit_dir(self) -> Path:
        return self.root / DIT_DIR_NAME

    @property
    def dit_toml(self) -> Path:
        return self.root / DIT_TOML_NAME

    @property
    def index_db(self) -> Path:
        return self.dit_dir / INDEX_DB_NAME

    @property
    def scope_toml(self) -> Path:
        return self.dit_dir / SCOPE_TOML_NAME

    def rel(self, path: Path) -> str:
        return path.resolve().relative_to(self.root.resolve()).as_posix()

    def abs(self, rel_path: str) -> Path:
        return (self.root / rel_path).resolve()


def find_repo(start: Path | None = None) -> Repo:
    origin = (start or Path.cwd()).resolve()
    current = origin
    while True:
        if (current / ".git").exists() and (current / DIT_TOML_NAME).is_file():
            return Repo(root=current)
        if (current / ".git").exists():
            # git repo without dit.toml yet (e.g. during init)
            return Repo(root=current)
        if current.parent == current:
            raise FileNotFoundError(f"git repository not found from {origin}")
        current = current.parent


def require_initialized(start: Path | None = None) -> Repo:
    repo = find_repo(start)
    if not repo.dit_toml.is_file():
        raise FileNotFoundError(
            f"{DIT_TOML_NAME} not found in {repo.root}. Run `dit init` first."
        )
    if not repo.dit_dir.is_dir():
        raise FileNotFoundError(
            f"{DIT_DIR_NAME}/ not found in {repo.root}. Run `dit init` first."
        )
    return repo
