import os
import logging
from .utils import ADDON_DIR


def start_logging(logfile_name: str = "side-packages-install") -> logging.Logger:
    """
    Configure and start logging to a file.

    Parameters
    ----------
    logfile_name : str, optional
        The name of the log file. Defaults to 'side-packages-install'.

    Returns
    -------
    logging.Logger
        A Logger object that can be used to write log messages.

    This function sets up a logging configuration with a specified log file name and logging level.
    The log file will be created in the `ADDON_DIR/logs` directory. If the directory
    does not exist, it will be created. The function returns a Logger object that can be used to
    write log messages.

    """
    # Create the logs directory if it doesn't exist
    logs_dir = os.path.join(os.path.abspath(ADDON_DIR), "logs")
    os.makedirs(logs_dir, exist_ok=True)

    # Set up logging configuration
    logfile_path = os.path.join(logs_dir, f"{logfile_name}.log")
    logging.basicConfig(filename=logfile_path, level=logging.INFO)

    # Return logger object
    return logging.getLogger()
