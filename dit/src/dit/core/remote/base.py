from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path


class Remote(ABC):
    @abstractmethod
    def exists(self, content_hash: str) -> bool:
        pass

    @abstractmethod
    def upload(self, local_path: Path, content_hash: str) -> None:
        pass

    @abstractmethod
    def download(self, content_hash: str, local_path: Path) -> None:
        pass

    @abstractmethod
    def delete(self, content_hash: str) -> None:
        pass

    @abstractmethod
    def list_hashes(self) -> list[str]:
        pass
