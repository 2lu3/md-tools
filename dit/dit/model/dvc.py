import glob
import os
from natsort import natsorted
import subprocess
import shutil
from alive_progress import alive_bar

from dit.model.extension import Extension
from dit.model.git import Git
from dit.model.scope import Scope


class DVC:
    def __init__(self):
        self.git = Git()
        self.scope = Scope()
        self.extension = Extension()

    def find_dvc_files_in_scope(self):
        files = []
        for directory in self.scope.directories:
            relpath = os.path.relpath(os.path.join(self.git.root_dir(), directory))
            files.extend(
                glob.glob(os.path.join(relpath, "**", "*.dvc"), recursive=True)
            )

        files = set(files)
        return natsorted(list(files))

    def find_dvc_files_in_project(self):
        files = []
        relpath = os.path.relpath(self.git.root_dir())
        files.extend(
            glob.glob(os.path.join(relpath, "**", "*.dvc"), recursive=True)
        )
        return natsorted(files)


    def clean(self):
        """cacheを全て削除し、scope以外のdvc管理下の実態ファイルを削除する"""

        # delete cache
        process = subprocess.Popen(
            ["dvc", "gc", "-w", "--not-in-remote", "-f"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True  # Python3.7以降なら、universal_newlines=TrueでもOK
        )
        
        # 出力を1行ずつ読みながらprint
        while True:
            output = process.stdout.readline() # type: ignore
            if output == "" and process.poll() is not None:
                break
            if output:
                print(output.strip())

        # delete files
        all_files = self.find_dvc_files_in_project()
        scope_files = self.find_dvc_files_in_scope()
        files = set(all_files) - set(scope_files)

        with alive_bar(len(files)) as bar:
            for file in files:
                large_file = os.path.splitext(file)[0]
                if os.path.exists(large_file):
                    os.remove(large_file)
                bar()
