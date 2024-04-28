import os
import shutil
import click
from loguru import logger

def copy_template(output_dir: str):
    """テンプレートのすべてのファイルをコピーする

    mode: min, eq, pr
    """

    source_dir = os.path.join(_template_dir())
    target_dir = os.path.join(output_dir, "template")

    shutil.copytree(source_dir, target_dir, dirs_exist_ok=True)
    logger.info("copy template files to {}".format(target_dir))

    if not os.path.exists(os.path.join(output_dir, "builder.py")):
        shutil.copy(os.path.join(target_dir, "builder.py"), output_dir)
        logger.info("copy builder.py to root direcotry")
    else:
        logger.warning("builder.py did not copy to root direcotry")


@click.command()
@click.option("--output_dir", "-o", default=".", help="output directory")
def copy_template_to_command(output_dir: str):
    copy_template(output_dir)


def _template_dir():
    return os.path.normpath(
        os.path.join(os.path.dirname(__file__), "..", "config", "template")
    )
