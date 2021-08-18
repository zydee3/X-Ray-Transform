from datetime import datetime
from pathlib import Path

log_folder_path = ""
log_folder_name = ""


def logger_init():
    global log_folder_path
    global log_folder_name

    log_folder_name = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")

    # go back in paths 1 by 1
    for path in reversed(Path(__file__).parents):
        if "python" in str(path):
            log_folder_path = path
