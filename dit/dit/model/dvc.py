import glob
import os

from dit.model.extension import Extension
from dit.model.scope import Scope


class DVC:
    def __init__(self):
        self.scope = Scope()
        self.extension = Extension()

    def find_dvc_files(self):
        files = []
        for directory in self.scope.directories:
            relpath = os.path.relpath(directory, start=os.getcwd())
            print(os.path.join(relpath, "**", "*.dvc"))
            print(glob.glob(os.path.join(relpath, "**", "*.dvc"), recursive=True))
            files.extend(
                glob.glob(os.path.join(relpath, "**", "*.dvc"), recursive=True)
            )

        files = set(files)
        return files

