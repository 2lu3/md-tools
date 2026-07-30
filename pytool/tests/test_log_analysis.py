import os

from pytool.log_analyzer import analyze_minimization


def test_png_existance():
    analyze_minimization("tests/data/min.log", "min", "Minimization", 10, False)

    if os.path.isfile("energy_min.png"):
        os.remove("energy_min.png")
    else:
        msg = "File not found"
        raise Exception(msg)


