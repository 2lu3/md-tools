from __future__ import annotations

import subprocess
from pathlib import Path

from dit.core.config import load_config
from dit.core.content import resolve_content_hash, write_pointer_for_file
from dit.core.index import StatIndex
from dit.core.pointer import read_pointer
from dit.core.repo import Repo
from dit.core.tracker import iter_pointer_files, iter_tracked_files


def run_add(repo: Repo, quiet: bool = False, prune: bool = False) -> int:
    config = load_config(repo)
    changed = 0
    with StatIndex(repo.index_db) as index:
        tracked = iter_tracked_files(repo, config)
        tracked_rels = {repo.rel(p) for p in tracked}
        for data_path in tracked:
            digest = resolve_content_hash(repo, index, data_path)
            pointer_path = Path(str(data_path) + ".dit")
            needs_write = True
            if pointer_path.is_file():
                existing = read_pointer(pointer_path)
                if existing.hash == digest and existing.size == data_path.stat().st_size:
                    needs_write = False
            if needs_write:
                pointer = write_pointer_for_file(repo, index, data_path, digest)
                _git_add(repo, repo.root / pointer.pointer_relpath)
                changed += 1
                if not quiet:
                    print(f"add {pointer.path}")
        if prune:
            changed += _prune_orphaned_pointers(repo, tracked_rels, quiet)
    return changed


def _prune_orphaned_pointers(repo: Repo, tracked_rels: set[str], quiet: bool) -> int:
    removed = 0
    for pointer_path in iter_pointer_files(repo):
        try:
            pointer = read_pointer(pointer_path)
        except ValueError:
            continue
        if pointer.path in tracked_rels:
            continue
        pointer_path.unlink(missing_ok=True)
        _git_add(repo, pointer_path)
        removed += 1
        if not quiet:
            print(f"prune {pointer.path}")
    return removed


def _git_add(repo: Repo, path: Path) -> None:
    try:
        subprocess.run(
            ["git", "add", "--", str(path)],
            cwd=repo.root,
            check=False,
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"git add failed for {path}: {exc}") from exc
