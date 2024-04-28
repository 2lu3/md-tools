from loguru import logger

def normalize_time(time_list: list[float]):
    result: list[float] = [time_list[0]]
    for i in range(len(time_list) - 1):
        if time_list[i + 1] < time_list[i]:
            result.append(time_list[i + 1] + result[i])
        else:
            result.append(time_list[i + 1] - time_list[i] + result[i])
    return result



def read_info_lines(log_paths: list[str]):
    for log_path in log_paths:
        with open(log_path, "r") as f:
            for line in f.readlines():
                if "INFO:" in line:
                    yield line


def read_column_names(log_paths: list[str]):
    def is_numeric_row(row: str):
        elements = [x for x in row.split(" ") if x != ""]
        return any([x.isnumeric() for x in elements])

    result: list[str] = []
    for line in read_info_lines(log_paths):
        if is_numeric_row(line):
            continue


        elements = [x for x in line.split(" ") if x != ""]

        if len(result) == 0:
            result = elements[1:]
            continue

        if result != elements[1:]:
            raise ValueError("Column names are not consistent")
    return result


def read_column_by_index(log_paths: list[str], column: int):
    result: list[float] = []
    for line in read_info_lines(log_paths):
        elements = [x for x in line.split(" ") if x != ""]
        try:
            result.append(float(elements[column+1]))
        except IndexError:
            logger.warning(f"Column {column} does not exist")
            continue
        except ValueError:
            continue

    return result


def read_column_by_name(log_paths: list[str], column: str):
    logger.debug(f"Reading column {column}")
    logger.debug(f"Found columns: {read_column_names(log_paths)}")
    target_column = read_column_names(log_paths).index(column)

    if column in ["TIME", "STEP"]:
        logger.debug("Normalizing time")
        return normalize_time(read_column_by_index(log_paths, target_column))

    return read_column_by_index(log_paths, target_column)
