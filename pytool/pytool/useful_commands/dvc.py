import subprocess
import click

def dvc_add(additional_globs=[]):
    subprocess.run(["dvc", "add", "--glob", "**/*.dcd", "**/*.rst", "**/*.dvl", "**/*.npy", "**/*.pkl", *additional_globs])

@click.command()
@click.option("--additional_globs", "-a", multiple=True, help="additional globs to add")
def dvc_add_to_command(additional_globs):
    dvc_add(additional_globs)
