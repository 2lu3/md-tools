"""Read and normalize log files."""

import itertools
from pathlib import Path

import pandas as pd


def read_column_by_index(log_paths: list[str], column: int) -> list[float]:
    """Read one column from log files by index.

    Args:
        log_paths: Log file paths to read.
        column: Zero-based column index.

    Returns:
        Values from the requested column.
    """
    df = read_log(log_paths)
    return list(df.iloc[:, column])


def read_column_by_name(log_paths: list[str], column: str) -> list[float]:
    """Read one column from log files by name.

    Args:
        log_paths: Log file paths to read.
        column: Column name to read.

    Returns:
        Values from the requested column.
    """
    df = read_log(log_paths)
    return list(df[column])


def read_log(log_paths: list[str]) -> pd.DataFrame:
    """Read log files into a data frame.

    Args:
        log_paths: Log file paths to read.

    Returns:
        Parsed numeric log rows.
    """
    df = pd.DataFrame(columns=_read_column_names(log_paths))

    for log_path in log_paths:
        with Path(log_path).open() as f:
            for line in f:
                if "INFO:" not in line:
                    continue

                if not _is_numeric_row(line):
                    continue

                elements = [x for x in line.split(" ") if x != ""]

                df = pd.concat(
                    [
                        df,
                        pd.DataFrame(
                            [[float(element) for element in elements[1:]]],
                            columns=df.columns,
                        ),
                    ],
                    ignore_index=True,
                )

    if "TIME" in df.columns:
        df["TIME"] = _normalize_time(df["TIME"])
    if "STEP" in df.columns:
        df["STEP"] = _normalize_time(df["STEP"])

    return df


def _is_numeric_row(row: str) -> bool:
    elements = [x for x in row.split(" ") if x != ""]
    return any(x.isnumeric() for x in elements)


def _read_column_names(log_paths: list[str]) -> list[str]:
    column_names: list[str] | None = None
    for lines in _read_info_lines(log_paths):
        for line in lines:
            if _is_numeric_row(line):
                continue

            elements = [x for x in line.split(" ") if x != ""]

            names = elements[1:]

            if column_names is None:
                column_names = names
            elif column_names != names:
                msg = "Column names are not consistent"
                raise ValueError(msg)
    return list(column_names or [])


def _read_info_lines(log_paths: list[str]) -> list[list[str]]:
    lines = []
    for log_path in log_paths:
        with Path(log_path).open() as f:
            lines.append([line for line in f if "INFO:" in line])

    return lines


def _normalize_time(time_list: list[float]) -> list[float]:
    """Normalize time values that reset between split MD runs.

    Args:
        time_list: Time values.

    Returns:
        Corrected time values.

    Examples:
        >>> normalize_time([1, 2, 3, 4, 5])
        [1, 2, 3, 4, 5]
        >>> normalize_time([10, 20, 30, 10, 20, 30])
        [10, 20, 30, 40, 50, 60]
    """
    # A single time value does not need correction.
    if len(time_list) <= 1:
        return list(time_list)

    # Strictly increasing time values have not reset between segments.
    if all(a < b for a, b in itertools.pairwise(time_list)):
        return list(time_list)

    result = []

    offset = 0
    prev = time_list[0]
    for t in time_list:
        # A lower value starts a new segment.
        if t < prev:
            offset += prev
        result.append(t + offset)
        prev = t

    return result
