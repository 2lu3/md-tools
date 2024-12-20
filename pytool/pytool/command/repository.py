from argparse import ArgumentParser
import os
import shutil
import requests
import subprocess
from loguru import logger

def _copy_from_template(file_name: str, out_path: str):
    root_dir = os.path.join(os.path.dirname(__file__), "template")
    shutil.copy(os.path.join(root_dir, file_name), out_path)

def _create_dirs():
    dirs = ["software", "tool", "data"]
    for d in dirs:
        os.makedirs(d, exist_ok=True)

def _download_latest_genesis():
    github_url = "https://api.github.com/repos/genesis-release-r-ccs/genesis/releases/latest"

    response = requests.get(github_url)
    response.raise_for_status()
    data = response.json()

    url = data["tarball_url"]
    response = requests.get(url)
    response.raise_for_status()
    with open("software/genesis.tar.gz", "wb") as f:
        f.write(response.content)

def _install_genesis():
    os.makedirs("software/genesis", exist_ok=True)
    subprocess.run(["tar", "xvf", "software/genesis.tar.gz", "-C", "software/genesis", "--strip-components=1"])

    os.chdir("software/genesis")
    subprocess.run(["./configure"])
    subprocess.run(["make", "-j"])
    subprocess.run(["make", "install"])
    os.chdir("../..")

def _install_requirements():
    _copy_from_template("install_requirements.sh", "tool")
    subprocess.run(["bash", "tool/install_requirements.sh"])

def _create_env():
    with open(".env", "w") as f:
        pass

    with open(".envrc", "w") as f:
        f.write("dotenv\n")
        f.write("export PATH=\"$PATH:$(pwd)/tool\"\n")
        f.write("export PATH=\"$PATH:$(pwd)/software/genesis/bin\"\n")

def create_repository():
    _create_dirs()

    _copy_from_template("switch.py", "software")

    _install_requirements()
    _download_latest_genesis()
    _install_genesis()

    _create_env()

    logger.warning("AWSの環境変数やFUGAKU_USER_IDを設定してください")

def to_command():
    parser = ArgumentParser("MD simulation用のレポジトリを自動作成する")

    _ = parser.parse_args()

    create_repository()

if __name__ == "__main__":
    to_command()
