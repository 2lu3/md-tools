from __future__ import annotations

import subprocess
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from dit.core.config import DitConfig, load_config
from dit.core.content import resolve_content_hash, utc_now_iso, write_pointer_for_file
from dit.core.index import StatIndex
from dit.core.pointer import Pointer, read_pointer
from dit.core.remote.base import Remote
from dit.core.remote.s3 import open_remote
from dit.core.repo import Repo
from dit.core.scope import Scope
from dit.core.tracker import iter_pointer_files, iter_tracked_files


class SyncAction(str, Enum):
    OK = "ok"
    PUSH = "push"
    PULL = "pull"
    UPDATE_POINTER = "update_pointer"
    DELETE_LOCAL = "delete_local"
    DELETE_REMOTE = "delete_remote"
    WARNING = "warning"
    ERROR = "error"


@dataclass
class SyncResult:
    path: str
    action: SyncAction
    message: str


def require_remote(config: DitConfig) -> Remote:
    if config.remote is None:
        raise ValueError("remote is not configured in dit.toml")
    return open_remote(config.remote)


def run_push(repo: Repo, dry_run: bool = False) -> list[SyncResult]:
    config = load_config(repo)
    remote = require_remote(config)
    results: list[SyncResult] = []
    with StatIndex(repo.index_db) as index:
        for pointer_path in iter_pointer_files(repo):
            pointer = read_pointer(pointer_path)
            data_path = repo.abs(pointer.path)
            if not data_path.is_file():
                continue
            if remote.exists(pointer.hash):
                continue
            results.append(SyncResult(pointer.path, SyncAction.PUSH, "upload"))
            if not dry_run:
                remote.upload(data_path, pointer.hash)
                index.mark_pushed(pointer.path, utc_now_iso())
    return results


def run_pull(repo: Repo, dry_run: bool = False) -> list[SyncResult]:
    config = load_config(repo)
    remote = require_remote(config)
    scope = Scope(repo)
    results: list[SyncResult] = []
    for pointer_path in iter_pointer_files(repo):
        pointer = read_pointer(pointer_path)
        if not scope.contains(pointer.path):
            continue
        data_path = repo.abs(pointer.path)
        if data_path.is_file():
            continue
        results.append(SyncResult(pointer.path, SyncAction.PULL, "download"))
        if not dry_run:
            remote.download(pointer.hash, data_path)
    return results


def run_sync(
    repo: Repo,
    dry_run: bool = False,
    prune_remote: bool = False,
) -> list[SyncResult]:
    config = load_config(repo)
    remote = require_remote(config)
    scope = Scope(repo)
    results: list[SyncResult] = []
    with StatIndex(repo.index_db) as index:
        tracked = {repo.rel(path): path for path in iter_tracked_files(repo, config)}
        pointers: dict[str, Pointer] = {}
        for pointer_path in iter_pointer_files(repo):
            pointer = read_pointer(pointer_path)
            pointers[pointer.path] = pointer

        for rel in sorted(set(tracked) | set(pointers)):
            data_path = tracked.get(rel) or repo.abs(rel)
            results.extend(
                _sync_one(
                    repo=repo,
                    index=index,
                    remote=remote,
                    scope=scope,
                    rel=rel,
                    data_path=data_path,
                    pointer=pointers.get(rel),
                    dry_run=dry_run,
                )
            )
        if prune_remote:
            results.extend(_prune_remote_orphans(repo, remote, dry_run))
    return results


def _sync_one(
    repo: Repo,
    index: StatIndex,
    remote: Remote,
    scope: Scope,
    rel: str,
    data_path: Path,
    pointer: Pointer | None,
    dry_run: bool,
) -> list[SyncResult]:
    has_data = data_path.is_file()
    if pointer is not None and has_data:
        return _sync_both(repo, index, remote, scope, rel, data_path, pointer, dry_run)
    if pointer is not None and not has_data:
        return _sync_pointer_only(remote, scope, rel, data_path, pointer, dry_run)
    if pointer is None and has_data:
        return [SyncResult(rel, SyncAction.WARNING, "missing pointer or untracked")]
    return []


def _sync_both(
    repo: Repo,
    index: StatIndex,
    remote: Remote,
    scope: Scope,
    rel: str,
    data_path: Path,
    pointer: Pointer,
    dry_run: bool,
) -> list[SyncResult]:
    local_hash = resolve_content_hash(repo, index, data_path)
    active = pointer
    results: list[SyncResult] = []

    if local_hash != pointer.hash:
        results.extend(
            _resolve_hash_mismatch(
                repo, index, remote, scope, rel, data_path, pointer, local_hash, dry_run
            )
        )
        if results and results[-1].action == SyncAction.ERROR:
            return results
        if not dry_run:
            active = read_pointer(repo.root / pointer.pointer_relpath)
        else:
            data_mtime = data_path.stat().st_mtime_ns
            pointer_mtime = (repo.root / pointer.pointer_relpath).stat().st_mtime_ns
            if data_mtime >= pointer_mtime:
                active = Pointer(path=rel, hash=local_hash, size=data_path.stat().st_size)

    remote_has = remote.exists(active.hash)
    if has_local_needing_push(data_path, remote_has):
        results.append(SyncResult(rel, SyncAction.PUSH, "upload"))
        if not dry_run:
            remote.upload(data_path, active.hash)
            index.mark_pushed(rel, utc_now_iso())
            remote_has = True

    if not scope.contains(rel):
        return results + _delete_local_if_safe(remote, rel, data_path, active.hash, dry_run)

    if not results:
        results.append(SyncResult(rel, SyncAction.OK, "in sync"))
    return results


def has_local_needing_push(data_path: Path, remote_has: bool) -> bool:
    return data_path.is_file() and not remote_has


def _resolve_hash_mismatch(
    repo: Repo,
    index: StatIndex,
    remote: Remote,
    scope: Scope,
    rel: str,
    data_path: Path,
    pointer: Pointer,
    local_hash: str,
    dry_run: bool,
) -> list[SyncResult]:
    data_mtime = data_path.stat().st_mtime_ns
    pointer_mtime = (repo.root / pointer.pointer_relpath).stat().st_mtime_ns
    if data_mtime >= pointer_mtime:
        results = [SyncResult(rel, SyncAction.UPDATE_POINTER, "local newer")]
        if not dry_run:
            write_pointer_for_file(repo, index, data_path, local_hash)
        return results

    if not remote.exists(pointer.hash):
        return [SyncResult(rel, SyncAction.ERROR, "pointer newer but remote missing")]
    if scope.contains(rel):
        if not dry_run:
            remote.download(pointer.hash, data_path)
        return [SyncResult(rel, SyncAction.PULL, "pointer newer")]
    return _delete_local_if_safe(remote, rel, data_path, pointer.hash, dry_run)


def _sync_pointer_only(
    remote: Remote,
    scope: Scope,
    rel: str,
    data_path: Path,
    pointer: Pointer,
    dry_run: bool,
) -> list[SyncResult]:
    if not remote.exists(pointer.hash):
        return [SyncResult(rel, SyncAction.ERROR, "file not found locally or remotely")]
    if not scope.contains(rel):
        return [SyncResult(rel, SyncAction.OK, "out of scope; no local copy")]
    if not dry_run:
        remote.download(pointer.hash, data_path)
    return [SyncResult(rel, SyncAction.PULL, "download")]


def _delete_local_if_safe(
    remote: Remote,
    rel: str,
    data_path: Path,
    content_hash: str,
    dry_run: bool,
) -> list[SyncResult]:
    if not remote.exists(content_hash):
        return [SyncResult(rel, SyncAction.ERROR, "cannot delete local; remote missing")]
    if not dry_run and data_path.is_file():
        data_path.unlink()
    return [SyncResult(rel, SyncAction.DELETE_LOCAL, "out of scope")]


def _prune_remote_orphans(repo: Repo, remote: Remote, dry_run: bool) -> list[SyncResult]:
    _git_fetch_all(repo)
    referenced = _hashes_from_all_refs(repo)
    results: list[SyncResult] = []
    for content_hash in remote.list_hashes():
        if content_hash in referenced:
            continue
        results.append(SyncResult(content_hash, SyncAction.DELETE_REMOTE, "orphan"))
        if not dry_run:
            remote.delete(content_hash)
    return results


def _git_fetch_all(repo: Repo) -> None:
    try:
        subprocess.run(
            ["git", "fetch", "--all", "--prune"],
            cwd=repo.root,
            check=True,
            capture_output=True,
            text=True,
            timeout=600,
        )
    except (OSError, subprocess.TimeoutExpired, subprocess.CalledProcessError) as exc:
        raise RuntimeError(f"git fetch --all --prune failed: {exc}") from exc


def _hashes_from_all_refs(repo: Repo) -> set[str]:
    try:
        commits = subprocess.run(
            ["git", "rev-list", "--all"],
            cwd=repo.root,
            check=True,
            capture_output=True,
            text=True,
            timeout=600,
        ).stdout.splitlines()
    except (OSError, subprocess.TimeoutExpired, subprocess.CalledProcessError) as exc:
        raise RuntimeError(f"git rev-list --all failed: {exc}") from exc

    hashes: set[str] = set()
    for commit in commits:
        names = _ls_tree_names(repo, commit)
        for name in names:
            if not name.endswith(".dit"):
                continue
            digest = _hash_from_blob(repo, commit, name)
            if digest:
                hashes.add(digest)
    for pointer_path in iter_pointer_files(repo):
        hashes.add(read_pointer(pointer_path).hash)
    return hashes


def _ls_tree_names(repo: Repo, commit: str) -> list[str]:
    try:
        result = subprocess.run(
            ["git", "ls-tree", "-r", "--name-only", commit],
            cwd=repo.root,
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except (OSError, subprocess.TimeoutExpired, subprocess.CalledProcessError):
        return []
    return result.stdout.splitlines()


def _hash_from_blob(repo: Repo, commit: str, name: str) -> str | None:
    try:
        blob = subprocess.run(
            ["git", "show", f"{commit}:{name}"],
            cwd=repo.root,
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        ).stdout
    except (OSError, subprocess.TimeoutExpired, subprocess.CalledProcessError):
        return None
    for line in blob.splitlines():
        if line.startswith("hash = "):
            return line.split("=", 1)[1].strip().strip('"')
    return None
