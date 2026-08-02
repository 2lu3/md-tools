"""Copy input files into an MD simulation directory layout."""

from pathlib import Path

from loguru import logger

from .util import copy_file_safe

DEFAULT_TOPPAR_FILES = [
    ("rtf", "top_all36_prot.rtf"),
    ("prm", "par_all36m_prot.prm"),
    ("str", "toppar_water_ions.str"),
]


def copy_toppar_files(
    toppar_dir: str,
    output_dir: str,
    files_to_copy: list[tuple[str, str]] | None = None,
) -> None:
    """topparからMDに必要なファイルをコピーする.

    Args:
        toppar_dir: topparディレクトリへのパス.
        output_dir: 出力先ディレクトリへのパス.
        files_to_copy: コピーするファイルのリスト.

    Examples:
        >>> copy_toppar_files("../01_data/toppar", "production_1")
    """
    if files_to_copy is None:
        files_to_copy = DEFAULT_TOPPAR_FILES

    def find_toppar_file(toppar_dir: str, file_name: str) -> Path:
        toppar_path = Path(toppar_dir)
        if not toppar_path.exists():
            msg = f"directory {toppar_dir} not found"
            raise FileNotFoundError(msg)
        file_paths = list(toppar_path.glob(file_name))
        if len(file_paths) == 0:
            msg = f"{file_name} not found in {toppar_dir}"
            raise FileNotFoundError(msg)
        if len(file_paths) > 1:
            msg = f"Multiple {file_name} found in {toppar_dir}"
            raise RuntimeError(msg)
        logger.debug(f"Found {file_name} in {toppar_dir}")
        return file_paths[0]

    for subdir, filename in files_to_copy:
        file_path = find_toppar_file(toppar_dir, filename)
        copy_file_safe(file_path, output_dir, subdir, filename)


def copy_structure_files(
    output_dir: str,
    pdb_file: str | None = None,
    psf_file: str | None = None,
    rst_file: str | None = None,
) -> None:
    """MDに必要な構造ファイル(PDB, PSF, RST)をコピーする.

    Args:
        output_dir: 出力先ディレクトリへのパス.
        pdb_file: PDBファイルへのパス.
        psf_file: PSFファイルへのパス.
        rst_file: RSTファイルへのパス.
    """
    if pdb_file is not None:
        copy_file_safe(pdb_file, output_dir, "pdb", "input.pdb")
    if psf_file is not None:
        copy_file_safe(psf_file, output_dir, "psf", "input.psf")
    if rst_file is not None:
        copy_file_safe(rst_file, output_dir, "rst", "input.rst")


def clean_directory(
    output_dir: str,
    is_delete_output: bool = False,  # noqa: FBT001, FBT002
    additional_dirs: list[str] | None = None,
) -> None:
    """output_dirを新しいプロジェクトファイルを生成できるように初期化する.

    - .dvcファイルは削除しない
    - .gitignoreファイルは削除しない

    Args:
        output_dir: プロジェクトのディレクトリ.
        is_delete_output: out/を削除するか(Trueなら削除).
        additional_dirs: 追加で生成するディレクトリ.
    """
    if additional_dirs is None:
        additional_dirs = []
    output_path = Path(output_dir)
    out_dir = output_path / "out"
    for file_path in output_path.rglob("*"):
        if not file_path.is_file():
            continue
        if file_path.suffix == ".dvc" or file_path.name == ".gitignore":
            continue

        if file_path.is_relative_to(out_dir) and not is_delete_output:
            continue

        file_path.unlink()

    # 空のディレクトリを削除
    directories = sorted(
        (path for path in output_path.rglob("*") if path.is_dir()),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    for directory in [*directories, output_path]:
        if directory.is_dir() and not any(directory.iterdir()):
            directory.rmdir()

    _create_directory(output_dir, additional_dirs)


def init_directory(output_dir: str) -> None:
    """output_dirを初期化する.

    すべてのファイルが削除される.

    Args:
        output_dir: output_dir.
    """
    msg = "This function is deprecated. Use clean_directory instead."
    raise NotImplementedError(msg)


def _create_directory(
    output_dir: str,
    additional_dirs: list[str] | None = None,
) -> None:
    if additional_dirs is None:
        additional_dirs = []
    directories = ["inp", "out", "rtf", "prm", "str", "pdb", "psf", "rst"]
    directories.extend(additional_dirs)

    for directory in directories:
        directory_path = Path(output_dir) / directory
        directory_path.mkdir(parents=True, exist_ok=True)
        (directory_path / ".gitkeep").touch()
        logger.debug(f"Create {directory} directory")

