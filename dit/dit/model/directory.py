import os
from glob import glob

class Directory:
    """ditが管理するディレクトリを示す"""

    def __init__(self, dir_path: str):
        self.dir_path = os.path.normpath(dir_path)

    def find_patterns(self, patterns: list[str]):
        files = []
        for pattern in patterns:
            files.extend(glob(os.path.join(self.dir_path, pattern)))

        files = set(files)

        return files

    @classmethod
    def load_config(cls):
        pass
