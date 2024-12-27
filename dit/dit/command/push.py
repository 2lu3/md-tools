import click
import subprocess

@click.command()
def push():
    dvc_result= subprocess.run(["dvc", "push"], capture_output=True, text=True)
    print(dvc_result)
    git_result = subprocess.run(["git", "push"], capture_output=True, text=True)
    print(git_result)

