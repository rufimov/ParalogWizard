import gzip
import logging
import logging.handlers
import multiprocessing
import os
import shutil
import sys
import threading
from functools import wraps


# =============================================================================
# LOGGER SETUP
# =============================================================================
def listener_configurer(log_file: str):
    root = logging.getLogger()
    handler = logging.FileHandler(log_file, mode="a")
    formatter = logging.Formatter(
        "%(asctime)s - %(processName)s - %(name)s - %(levelname)s - %(message)s",
        "%d-%b-%y %H:%M:%S",
    )
    handler.setFormatter(formatter)
    root.addHandler(handler)
    root.setLevel(logging.INFO)


def listener_process(log_queue: multiprocessing.Queue, log_file: str):
    listener_configurer(log_file)
    while True:
        record = log_queue.get()
        if record is None:  # Sentinel value to stop.
            break
        logger = logging.getLogger(record.name)
        logger.handle(record)


def worker_initializer(log_queue: multiprocessing.Queue):
    """
    Initializer for each worker process in a multiprocessing pool.
    Attaches a QueueHandler to the worker's root logger so that all log records
    are sent to the shared logging queue.
    """
    queue_handler = logging.handlers.QueueHandler(log_queue)
    root = logging.getLogger()
    # Remove any pre-existing handlers
    root.handlers = []
    root.addHandler(queue_handler)
    root.setLevel(logging.INFO)


def setup_logging():
    import logging
    import os
    import logging.handlers

    logger = logging.getLogger("ParalogWizard")

    # If a QueueHandler is already attached, assume this is a worker process.
    if any(isinstance(h, logging.handlers.QueueHandler) for h in logger.handlers):
        return logger

    if not logger.handlers:
        log_file = os.environ.get("PARALOGWIZARD_LOGFILE", "ParalogWizard.log")
        file_handler = logging.FileHandler(log_file)
        stream_handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "[%(asctime)s] [%(processName)s:%(process)d] [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler.setFormatter(formatter)
        stream_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.addHandler(stream_handler)
        logger.setLevel(logging.INFO)
    return logger


# =============================================================================
# Logging Decorator
# =============================================================================
def log_exceptions(func):
    """
    Decorator that logs function entry, exit, and any exceptions.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger = setup_logging()
        try:
            result = func(*args, **kwargs)
            return result
        except Exception as e:
            logger.exception(f"Exception in {func.__name__}: {e}")
            raise

    return wrapper


# =============================================================================
# DECOMPRESS AND COMPRESS FASTQ FILES
# =============================================================================
def decompress_fastq_files(files):
    """
    For each file in the list, if it ends with '.gz',
    decompress it to a file without the .gz extension.
    Returns a list of file paths that are now uncompressed.
    """
    logger = logging.getLogger("ParalogWizard")
    new_files = []
    for f in files:
        if f.endswith(".gz"):
            uncompressed = f[:-3]
            logger.info("Decompressing %s to %s", f, uncompressed)
            try:
                # Explicitly cast the returned objects to BinaryIO.
                fin = gzip.open(f, "rb")
                fout = open(uncompressed, "wb")
                with fin, fout:
                    shutil.copyfileobj(fin, fout)
            except Exception as e:
                logger.error("Error decompressing %s: %s", f, e)
                sys.exit(1)
            new_files.append(uncompressed)
        else:
            new_files.append(f)
    return new_files


def compress_fastq_files(files):
    """
    For each file in the list, if it is not gzipped (i.e. does not end with '.gz'),
    compress it back to a .gz file and remove the uncompressed version.
    """
    logger = logging.getLogger("ParalogWizard")
    for f in files:
        if not f.endswith(".gz"):
            gz_file = f + ".gz"
            logger.info("Compressing %s to %s", f, gz_file)
            try:
                # Explicitly cast to BinaryIO.
                fin = open(f, "rb")
                fout = gzip.open(gz_file, "wb")
                with fin, fout:
                    shutil.copyfileobj(fin, fout)
                os.remove(f)
            except Exception as e:
                logger.error("Error compressing %s: %s", f, e)
