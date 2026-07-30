from __future__ import annotations

import os
from pathlib import Path

import pathspec

from dit.core.config import DitConfig
from dit.core.repo import POINTER_SUFFIX, Repo

SKIP_DIR_NAMES = {".git", ".dit", "__pycache__", ".venv", "venv"}


def build_track_spec(patterns: list[str]) -> pathspec.PathSpec:
    return pathspec.PathSpec.from_lines("gitwildmatch", patterns)


def load_gitignore_spec(repo: Repo) -> pathspec.PathSpec:
    patterns: list[str] = []
    gitignore = repo.root / ".gitignore"
    if gitignore.is_file():
        patterns.extend(gitignore.read_text(encoding="utf-8").splitlines())
    return pathspec.PathSpec.from_lines("gitwildmatch", patterns)


def is_tracked_path(rel_posix: str, is_dir: bool, track_spec: pathspec.PathSpec) -> bool:
    candidate = rel_posix + ("/" if is_dir and not rel_posix.endswith("/") else "")
    return track_spec.match_file(candidate)


def iter_tracked_files(repo: Repo, config: DitConfig) -> list[Path]:
    track_spec = build_track_spec(config.track_patterns)
    ignore_spec = load_gitignore_spec(repo)
    results: list[Path] = []
    for path in _walk(repo.root):
        rel = path.relative_to(repo.root).as_posix()
        if rel.endswith(POINTER_SUFFIX):
            continue
        if ignore_spec.match_file(rel):
            # still allow tracked large files that are gitignored
            if not is_tracked_path(rel, False, track_spec):
                continue
        if is_tracked_path(rel, False, track_spec):
            results.append(path)
    return sorted(results)


def iter_pointer_files(repo: Repo) -> list[Path]:
    results: list[Path] = []
    for path in _walk(repo.root):
        if path.name.endswith(POINTER_SUFFIX) and path.is_file():
            results.append(path)
    return sorted(results)


def _walk(root: Path) -> list[Path]:
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIR_NAMES]
        for name in filenames:
            found.append(Path(dirpath) / name)
    return found
