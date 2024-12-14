import os
from glob import glob
from .config import Configuration
from .git import Git



class Directory:
    """ditが管理するディレクトリを示す"""

    def __init__(self):
        self.directories = self._load_config()

    def find_patterns(self, patterns: list[str]):
        files = []
        for directory in self.directories:
            for pattern in patterns:
                files.extend(glob(os.path.join(directory, pattern)))

        files = set(files)
        return files

    def add_directory(self, dir_path: str):
        self.directories.append(os.path.normpath(dir_path))
        self._save_config()

    def remove_directory(self, dir_path: str):
        self.directories.remove(os.path.normpath(dir_path))
        self._save_config()

    def _load_config(self):
        git = Git()
        config = Configuration(git.root_dir())

        return config.load_config().get("directories", [])

    def _save_config(self):
        git = Git()
        config = Configuration(git.root_dir())
        config.update({"directories": self.directories})

