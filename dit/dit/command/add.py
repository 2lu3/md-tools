import argparse
import os
from dit.model.extension import Extension
from dit.model.scope import Scope
import glob
import subprocess

def search_files(path_list: list[str]):
    files_to_consider = []
    for path in path_list:
        path = os.path.normpath(path)

        if os.path.isfile(path):
            files_to_consider.append(path)
        else:
            for file in glob.glob(os.path.join(path, "**", "*"), recursive=True):
                files_to_consider.append(file)
    return files_to_consider


def add(path_list: list[str], is_all: bool):
    files_to_consider = search_files(path_list)
    if is_all:
        scope = Scope()
        files_to_consider.extend(
                search_files(scope.directories)
                )

    extension = Extension()
    files_by_git = []
    files_by_dvc = []

    for file in files_to_consider:
        if extension.is_match(file):
            files_by_dvc.append(file)
        else:
            files_by_git.append(file)

    git_result = subprocess.run(["git", "add", " ".join(files_by_git)], capture_output=True, text=True)
    print(git_result)
    dvc_result = subprocess.run(["dvc", "add", " ".join(files_by_dvc)], capture_output=True, text=True)
    print(dvc_result)


def register_subparser(subparser):
    parser = subparser.add_parser("add", help="git add/dvc addを行う。")
    parser.add_argument("path", help="git add/dvc addする対象のディレクトリ/ファイルを選ぶ。それ以外は実行されない。Scopeには影響しない", type=str, nargs='*')
    parser.add_argument("-A", "--all", help="Scopeの中の全てのファイルに対してgit add/dvc addを行う。", action="store_true")


def handle(args: argparse.Namespace):
    if args.command != "add":
        return
    add(args.path, args.all)

