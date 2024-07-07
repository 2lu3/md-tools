import os
import shutil
from loguru import logger

def copy_file_safe(source_file: str, output_dir: str, subdir: str, dest_filename: str):
    """Copy file after checking source file exists, creating output directory

    Args:
        source_file (str): source_file
        output_dir (str): output_dir
        subdir (str): subdir
        dest_filename (str): dest_filename
    """
    if not os.path.exists(source_file):
        raise FileNotFoundError(f"{source_file} not found")

    dest_dir = os.path.join(output_dir, subdir)
    os.makedirs(dest_dir, exist_ok=True)

    dest_file = os.path.join(dest_dir, dest_filename)
    shutil.copy(source_file, dest_file)
    logger.debug(f"Copy {source_file} to {dest_file}")


