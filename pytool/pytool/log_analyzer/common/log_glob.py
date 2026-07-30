"""Find log files for analysis."""

from pathlib import Path

import natsort


def glob_log_files(log_path: str) -> list[str]:
    """Return naturally sorted log file paths.

    Args:
        log_path: Directory containing ``*.log*`` files or a log file path.

    Returns:
        Naturally sorted log file paths.
    """
    path = Path(log_path)
    if path.is_file():
        return [log_path]
    if path.is_dir():
        return natsort.natsorted(
            str(log_file) for log_file in path.rglob("*.log*")
        )
    msg = f"Invalid path: {log_path}"
    raise ValueError(msg)
