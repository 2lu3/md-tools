import subprocess
import os
import click


@click.command()
def install_genesis_to_command():
    script_path = os.path.join(os.path.dirname(__file__), "install_genesis.sh")

    subprocess.call(["bash", script_path])

@click.command()
def install_genesis_requirements_to_command():
    script_path = os.path.join(os.path.dirname(__file__), "install_genesis_requirements.sh")

    subprocess.call(["bash", script_path])

