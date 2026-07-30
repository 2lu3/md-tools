from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

import tomli_w

from dit.core.repo import Repo

DEFAULT_TRACK_PATTERNS: list[str] = [
    "*.dcd",
    "*.dvl",
    "*.rst",
    "*.npy",
    "*.pkl",
    "*.tar",
]

S3_REQUEST_TIMEOUT_SEC = 60


@dataclass(frozen=True)
class RemoteConfig:
    url: str
    endpoint_url: str | None = None

    @property
    def bucket(self) -> str:
        parsed = urlparse(self.url)
        if parsed.scheme != "s3":
            raise ValueError(f"unsupported remote scheme: {parsed.scheme!r} in {self.url}")
        if not parsed.netloc:
            raise ValueError(f"missing bucket in remote url: {self.url}")
        return parsed.netloc

    @property
    def prefix(self) -> str:
        parsed = urlparse(self.url)
        return parsed.path.lstrip("/")


@dataclass(frozen=True)
class DitConfig:
    remote: RemoteConfig | None = None
    track_patterns: list[str] = field(default_factory=lambda: list(DEFAULT_TRACK_PATTERNS))

    @classmethod
    def load(cls, path: Path) -> DitConfig:
        try:
            with path.open("rb") as f:
                data = tomllib.load(f)
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"config not found: {path}") from exc
        except tomllib.TOMLDecodeError as exc:
            raise ValueError(f"invalid TOML in {path}: {exc}") from exc
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DitConfig:
        remote_data = data.get("remote") or {}
        remote: RemoteConfig | None = None
        if remote_data.get("url"):
            remote = RemoteConfig(
                url=str(remote_data["url"]),
                endpoint_url=remote_data.get("endpoint_url"),
            )
        track = data.get("track") or {}
        patterns = track.get("patterns") or list(DEFAULT_TRACK_PATTERNS)
        if not isinstance(patterns, list):
            raise ValueError("[track].patterns must be a list of strings")
        return cls(remote=remote, track_patterns=[str(p) for p in patterns])

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {"track": {"patterns": list(self.track_patterns)}}
        if self.remote is not None:
            remote: dict[str, str] = {"url": self.remote.url}
            if self.remote.endpoint_url:
                remote["endpoint_url"] = self.remote.endpoint_url
            result["remote"] = remote
        return result

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as f:
            tomli_w.dump(self.to_dict(), f)


def load_config(repo: Repo) -> DitConfig:
    return DitConfig.load(repo.dit_toml)


def default_init_config(
    remote_url: str | None = None,
    endpoint_url: str | None = None,
) -> DitConfig:
    remote = None
    if remote_url:
        remote = RemoteConfig(url=remote_url, endpoint_url=endpoint_url)
    return DitConfig(remote=remote, track_patterns=list(DEFAULT_TRACK_PATTERNS))
