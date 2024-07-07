import os
import glob
import shutil
from typing import Optional
from loguru import logger

from .util import copy_file_safe


def copy_toppar_files(
    toppar_dir: str,
    output_dir: str,
    files_to_copy: list[tuple[str, str]] = [
        ("rtf", "top_all36_prot.rtf"),
        ("prm", "par_all36m_prot.prm"),
        ("str", "toppar_water_ions.str"),
    ],
):
    """topparからMDに必要なファイルをコピーする

    Args:
        toppar_dir (str): topparディレクトリへのパス
        output_dir (str): 出力先ディレクトリへのパス
        files_tocopy (list[tuple[str, str]], optional): コピーするファイルのリスト. Defaults to [("rtf", "top_all36_prot.rtf"), ("prm", "par_all36m_prot.prm"), ("str", "toppar_water_ions.str")].

    Examples:
        >>> copy_toppar_files("../01_data/toppar", "production_1")
    """

    def find_toppar_file(toppar_dir: str, file_name: str):
        assert os.path.exists(toppar_dir), f"directory {toppar_dir} not found"
        file_paths = glob.glob(os.path.join(toppar_dir, file_name))
        if len(file_paths) == 0:
            raise FileNotFoundError(f"{file_name} not found in {toppar_dir}")
        elif len(file_paths) > 1:
            raise RuntimeError(f"Multiple {file_name} found in {toppar_dir}")
        else:
            logger.debug(f"Found {file_name} in {toppar_dir}")
            return file_paths[0]

    for subdir, filename in files_to_copy:
        file_path = find_toppar_file(toppar_dir, filename)
        copy_file_safe(file_path, output_dir, subdir, filename)


def copy_structure_files(
    output_dir: str, pdb_file: Optional[str] = None, psf_file: Optional[str] = None, rst_file: Optional[str] = None
):
    """MDに必要な構造ファイル(PDB, PSF, RST)をコピーする

    Args:
        output_dir (str): 出力先ディレクトリへのパス
        pdb_file (str): pdb_file
        psf_file (str): psf_file
        rst_file (Optional[str]): rst_file
    """
    if pdb_file is not None:
        copy_file_safe(pdb_file, output_dir, "pdb", "input.pdb")
    if psf_file is not None:
        copy_file_safe(psf_file, output_dir, "psf", "input.psf")
    if rst_file is not None:
        copy_file_safe(rst_file, output_dir, "rst", "input.rst")


def clean_directory(output_dir: str, is_delete_output: bool = False):
    """output_dirを新しいプロジェクトファイルを生成できるように初期化する

    - .dvcファイルは削除しない
    - .gitignoreファイルは削除しない

    Args:
        is_delete_output (str): out/を削除するか(Trueなら削除)
    """

    files = glob.glob(os.path.join(output_dir, "**", "*"), recursive=True)
    for file in files:
        if not os.path.isfile(file):
            continue
        if file.endswith(".dvc") or file.endswith(".gitignore"):
            continue

        if file.startswith(os.path.join(output_dir, "out")):
            if not is_delete_output:
                continue

        os.remove(file)

    # 空のディレクトリを削除
    for root, dirs, files in os.walk(output_dir, topdown=False):
        if len(dirs) == 0 and len(files) == 0:
            os.rmdir(root)

    _create_directory(output_dir)


def init_directory(output_dir: str):
    """output_dirを初期化する
    すべてのファイルが削除される

    Args:
        output_dir (str): output_dir
    """
    raise NotImplementedError("This function is deprecated. Use clean_directory instead.")

def _create_directory(output_dir: str):
    directories = ["inp", "out", "rtf", "prm", "str", "pdb", "psf", "rst"]
    for directory in directories:
        os.makedirs(os.path.join(output_dir, directory), exist_ok=True)
        with open(os.path.join(output_dir, directory, ".gitkeep"), "w"):
            pass
        logger.debug(f"Create {directory} directory")

