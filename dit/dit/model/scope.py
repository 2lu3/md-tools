import os
from glob import glob
from .config import Configuration
from .git import Git
from loguru import logger


class Scope:
    """ditが管理するディレクトリを示す"""

    directories: set[str] = set()

    def __init__(self):
        self.git = Git()
        self._load_config()

    def find_patterns(self, file_patterns: list[str]):
        files = []
        for directory in Scope.directories:
            for pattern in file_patterns:
                files.extend(glob(os.path.join(directory, f"**/{pattern}")))

        files = set(files)
        logger.debug(f"patterns: {file_patterns}: {len(files)} files found")
        return files

    def add(self, dir_path: str):
        Scope.directories.add(self._norm_path(dir_path))
        self._save_config()
        logger.debug(f"added {dir_path} to scope")

    def remove(self, dir_path: str):
        assert os.path.isdir(dir_path), f"{dir_path} is not a directory"
        assert (
            self._norm_path(dir_path) in Scope.directories
        ), f"{self._norm_path(dir_path)} is not in scope. {Scope.directories}"
        Scope.directories.remove(self._norm_path(dir_path))
        self._save_config()
        logger.debug(f"removed {dir_path} from scope")

    def _load_config(self):
        config = Configuration()

        directories: list[str] = config.load_config().get("directories", [])
        assert type(directories) == list

        Scope.directories = set(directories)
        logger.debug(f"loaded {len(directories)} directories")

    def _save_config(self):
        config = Configuration()
        config.update({"directories": list(Scope.directories)})

    def _norm_path(self, path: str):
        return os.path.relpath(path, start=self.git.root_dir())
