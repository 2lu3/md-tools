import click
import subprocess

@click.command()
@click.argument("message")
def commit(message: str):
    result = subprocess.run(["git", "commit", "-m", message], capture_output=True, text=True)
    print(result.stdout)
