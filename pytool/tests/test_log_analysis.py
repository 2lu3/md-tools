from pathlib import Path

from pytool.log_analyzer import analyze_minimization


def test_png_existance():
    analyze_minimization("tests/data/min.log", "min", "Minimization", 10, popup=False)

    output_path = Path("energy_min.png")
    if output_path.is_file():
        output_path.unlink()
    else:
        msg = "File not found"
        raise FileNotFoundError(msg)


