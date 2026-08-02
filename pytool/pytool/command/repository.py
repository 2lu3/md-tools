"""Create a project repository layout for MD simulation."""

import shutil
import subprocess
from argparse import ArgumentParser
from pathlib import Path

import requests
from loguru import logger

REQUEST_TIMEOUT_SECONDS = 30


def _copy_from_template(file_name: str, out_path: str) -> None:
    root_dir = Path(__file__).parent / "template"
    shutil.copy(root_dir / file_name, out_path)


def _create_dirs() -> None:
    dirs = ["software", "tool", "data"]
    for d in dirs:
        Path(d).mkdir(exist_ok=True)


def _download_latest_genesis() -> None:
    github_url = "https://api.github.com/repos/genesis-release-r-ccs/genesis/releases/latest"

    response = requests.get(github_url, timeout=REQUEST_TIMEOUT_SECONDS)
    response.raise_for_status()
    data = response.json()

    url = data["tarball_url"]
    response = requests.get(url, timeout=REQUEST_TIMEOUT_SECONDS)
    response.raise_for_status()
    with Path("software/genesis.tar.gz").open("wb") as f:
        f.write(response.content)


def _resolve_executable(executable: str) -> str:
    resolved = shutil.which(executable)
    if resolved is None:
        msg = f"{executable} executable not found"
        raise FileNotFoundError(msg)
    return resolved


def _install_genesis() -> None:
    genesis_dir = Path("software/genesis").resolve()
    genesis_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(  # noqa: S603
        [
            _resolve_executable("tar"),
            "xvf",
            "software/genesis.tar.gz",
            "-C",
            "software/genesis",
            "--strip-components=1",
        ],
        check=True,
    )

    subprocess.run(  # noqa: S603
        [str(genesis_dir / "configure")],
        check=True,
        cwd=genesis_dir,
    )
    make = _resolve_executable("make")
    subprocess.run([make, "-j"], check=True, cwd=genesis_dir)  # noqa: S603
    subprocess.run([make, "install"], check=True, cwd=genesis_dir)  # noqa: S603


def _install_requirements() -> None:
    _copy_from_template("install_requirements.sh", "tool")
    subprocess.run(  # noqa: S603
        [_resolve_executable("bash"), "tool/install_requirements.sh"],
        check=True,
    )


def _create_env() -> None:
    Path(".env").touch()

    with Path(".envrc").open("w") as f:
        f.write("dotenv\n")
        f.write('export PATH="$PATH:$(pwd)/tool"\n')
        f.write('export PATH="$PATH:$(pwd)/software/genesis/bin"\n')


def create_repository() -> None:
    """Create the standard MD simulation repository layout."""
    _create_dirs()

    _copy_from_template("switch.py", "software")

    _install_requirements()
    _download_latest_genesis()
    _install_genesis()

    _create_env()

    logger.warning("AWSの環境変数やFUGAKU_USER_IDを設定してください")


def to_command() -> None:
    """Run the create-repository command-line interface."""
    parser = ArgumentParser("MD simulation用のレポジトリを自動作成する")

    _ = parser.parse_args()

    create_repository()


if __name__ == "__main__":
    to_command()
