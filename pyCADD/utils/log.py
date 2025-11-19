import logging
import os
from datetime import datetime

from concurrent_log_handler import ConcurrentRotatingFileHandler
from rich.logging import RichHandler


class ModifiedFileHandler(ConcurrentRotatingFileHandler):
    def emit(self, record):
        if not os.path.exists(os.path.dirname(self.baseFilename)):
            os.makedirs(os.path.dirname(self.baseFilename))
        return super().emit(record)


def get_logfile_name():
    """
    Get log file name. Generate log file name based on current date.

    Returns:
        str: log file name
    """
    log_dir = os.path.join(os.getcwd(), "logs")

    date = datetime.now()
    year = str(date.year)
    month = str(date.month)
    day = str(date.day)
    now = year + month.rjust(2, "0") + day.rjust(2, "0")

    logfile_name = os.path.join(log_dir, f"{now}.log")
    return logfile_name


def _init_log(logname):
    """
    Initialize log file.

    Args:
        logname (str): logger name

    Returns:
        logging.Logger: logger object
    """
    logger = logging.getLogger(logname)
    logger.setLevel(logging.DEBUG)
    file_fmt = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    console_fmt = logging.Formatter("%(message)s")
    logfile = get_logfile_name()

    file_handler = ModifiedFileHandler(logfile, "a")
    file_handler.setFormatter(file_fmt)
    file_handler.setLevel(logging.DEBUG)

    console_handler = RichHandler(
        show_path=False, rich_tracebacks=True, log_time_format="[%Y-%m-%d %H:%M:%S]"
    )
    console_handler.setFormatter(console_fmt)
    console_handler.setLevel(logging.INFO)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger
