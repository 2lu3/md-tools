import os
import shutil
from pathlib import Path

from pytool.input_builder import write_template


def test_input_builder():

    sandbox = Path("sandbox")
    shutil.rmtree("sandbox", ignore_errors=True)
    sandbox.mkdir()
    os.chdir(sandbox)

    write_template("min")

    os.chdir("..")
    shutil.rmtree("sandbox", ignore_errors=True)
