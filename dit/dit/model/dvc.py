import glob
import os

from dit.model.extension import Extension
from dit.model.git import Git
from dit.model.scope import Scope


class DVC:
    def __init__(self):
        self.git = Git()
        self.scope = Scope()
        self.extension = Extension()

    def find_dvc_files(self):
        files = []
        for directory in self.scope.directories:
            relpath = os.path.relpath(os.path.join(self.git.root_dir(), directory))
            files.extend(
                glob.glob(os.path.join(relpath, "**", "*.dvc"), recursive=True)
            )

        files = set(files)
        return files
