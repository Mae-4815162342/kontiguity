"""Global logging configuration for the kontiguity Python side.

Every top-level command (load, retrieve, map, describe, pipeline) gets a
logger via get_logger(outpath), which writes to:
    <outpath>/logs/kontiguity.log     (persistent, shared across a whole run)
and, only if verbose=True, to the console. File logging is always on;
console output is opt-in (via each command's --verbose flag) since a large
dataset can produce a lot of INFO lines that are more useful in the log
file than scrolling past in the terminal.

Calling get_logger() again with the same outpath returns the same logger
(and does not duplicate handlers) - the verbose value from the *first* call
for a given outpath is the one that sticks, so pass it consistently.
"""

import logging
import os

_LOGGERS = {}


def get_logger(outpath, name="kontiguity", verbose=False):
    """Returns a logger writing to <outpath>/logs/kontiguity.log, and to the
    console only if verbose is True.

    Safe to call multiple times for the same outpath: returns the cached logger
    instead of re-adding handlers.
    """
    log_dir = f"{outpath}/logs"
    os.makedirs(log_dir, exist_ok=True)
    log_file = f"{log_dir}/kontiguity.log"

    key = os.path.abspath(log_file)
    if key in _LOGGERS:
        return _LOGGERS[key]

    logger = logging.getLogger(f"{name}.{key}")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    if verbose:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    _LOGGERS[key] = logger
    return logger


def global_log_path(outpath):
    """Returns the path of the shared global log file for a given outpath,
    creating the logs/ folder if needed. Used to export GLOBAL_LOG to
    generated bash scripts so they log to the same file as the Python side."""
    log_dir = f"{outpath}/logs"
    os.makedirs(log_dir, exist_ok=True)
    return f"{log_dir}/kontiguity.log"