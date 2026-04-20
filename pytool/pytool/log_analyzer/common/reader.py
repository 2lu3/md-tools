from tarfile import LENGTH_NAME
from loguru import logger
import pandas as pd


def read_column_by_index(log_paths: list[str], column: int) -> list[float]:
    df = read_log(log_paths)
    return list(df.iloc[:, column])


def read_column_by_name(log_paths: list[str], column: str) -> list[float]:
    df = read_log(log_paths)
    return list(df[column])


def read_log(log_paths: list[str]) -> pd.DataFrame:
    df = pd.DataFrame(columns=_read_column_names(log_paths))

    for log_path in log_paths:
        with open(log_path, "r") as f:
            for line in f.readlines():
                if "INFO:" not in line:
                    continue

                if not _is_numeric_row(line):
                    continue

                elements = [x for x in line.split(" ") if x != ""]

                df = df.append(
                    pd.Series(
                        [float(element) for element in elements[1:]], index=df.columns
                    ),
                    ignore_index=True,
                )

    if "TIME" in df.columns:
        df["TIME"] = _normalize_time(df["TIME"])
    if "STEP" in df.columns:
        df["STEP"] = _normalize_time(df["STEP"])

    return df


def _is_numeric_row(row: str):
    elements = [x for x in row.split(" ") if x != ""]
    return any([x.isnumeric() for x in elements])


def _read_column_names(log_paths: list[str]):
    column_names: set[str] = set()
    for lines in _read_info_lines(log_paths):
        for line in lines:
            if _is_numeric_row(line):
                continue

            elements = [x for x in line.split(" ") if x != ""]

            names = elements[1:]

            if len(column_names) == 0:
                column_names = names
            else:
                assert column_names == names, "Column names are not consistent"
    return list(column_names)


def _read_info_lines(log_paths: list[str]) -> list[list[str]]:
    lines = []
    for log_path in log_paths:
        with open(log_path, "r") as f:
            lines.append([line for line in f.readlines() if "INFO:" in line])

    return lines


def _normalize_time(time_list: list[float]) -> list[float]:
    """
    MD計算を分割して行った場合、分割された計算ごとに時刻が0スタートになるので修正する

    Args:
        time_list (list[float]): 時刻のリスト

    Returns:
        list[float]: 修正された時刻のリスト

    Examples:
        >>> normalize_time([1, 2, 3, 4, 5])
        [1, 2, 3, 4, 5]
        >>> normalize_time([10, 20, 30, 10, 20, 30])
        [10, 20, 30, 40, 50, 60]
    """

    # 時刻が1つなら修正しなくて良い
    if len(time_list) <= 1:
        return list(time_list)

    # 常に t[i] < t[i+1] ならtimeがリセットされていないので、修正しなくていい
    if all(a < b for a, b in zip(time_list, time_list[1:])):
        return list(time_list)

    result = []

    offset = 0
    prev = time_list[0]
    for t in time_list:
        # リセット検出（前より小さくなったら新しいセグメント）
        if t < prev:
            offset += prev
        result.append(t + offset)
        prev = t

    return result
