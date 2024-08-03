import os
import glob
import natsort


def glob_log_files(log_path: str):
    """ログファイルのパスを全て取得する
    取得したログファイルはnatsortでソートされる

    Args:
        log_path (str): *.log*ファイルが存在するディレクトリのパスもしくはログファイルそのもの
    """
    if os.path.isfile(log_path):
        return [log_path]
    if os.path.isdir(log_path):
        return natsort.natsorted(
            glob.glob(os.path.join(log_path, "**", "*.log*"), recursive=True)
        )
    raise Exception(f"Invalid path: {log_path}")
