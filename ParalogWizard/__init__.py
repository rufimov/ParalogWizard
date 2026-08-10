import gzip
import logging
import logging.handlers
import multiprocessing
import os
import shutil
import sys
import threading
from contextlib import contextmanager
from functools import wraps


# =============================================================================
# LOGGER SETUP
# =============================================================================
class LevelRangeFilter(logging.Filter):
    """Pass only records with min_level <= levelno <= max_level."""

    def __init__(self, min_level, max_level=logging.CRITICAL):
        super(LevelRangeFilter, self).__init__()
        self.min_level = min_level
        self.max_level = max_level

    def filter(self, record):
        return self.min_level <= record.levelno <= self.max_level


def _log_formatter():
    return logging.Formatter(
        "%(asctime)s - %(processName)s - %(name)s - %(levelname)s - %(message)s",
        "%d-%b-%y %H:%M:%S",
    )


def _make_file_handler(path, min_level, max_level=logging.CRITICAL):
    handler = logging.FileHandler(path, mode="a")
    handler.setLevel(min_level)
    handler.addFilter(LevelRangeFilter(min_level, max_level))
    handler.setFormatter(_log_formatter())
    return handler


def log_paths_from_base(log_file, debug=False):
    """
    Derive split log paths from the run base name.

    Example base ParalogWizard_cast_call_10.Aug.26_12:00.log yields:
      .log, .errors.log, and optionally .debug.log
    """
    base, ext = os.path.splitext(log_file)
    if ext.lower() != ".log":
        base = log_file
    paths = {
        "info": f"{base}.log",
        "errors": f"{base}.errors.log",
    }
    if debug:
        paths["debug"] = f"{base}.debug.log"
    return paths


def listener_configurer(log_file: str, debug: bool = False):
    """
    Configure the logging listener process.

    Always writes:
      <base>.log         — INFO and WARNING (normal run log)
      <base>.errors.log  — ERROR and CRITICAL
    With debug=True also writes:
      <base>.debug.log   — DEBUG and above (full detail)
    """
    root = logging.getLogger()
    root.handlers = []
    root.setLevel(logging.DEBUG if debug else logging.INFO)

    paths = log_paths_from_base(log_file, debug=debug)
    root.addHandler(
        _make_file_handler(paths["info"], logging.INFO, logging.WARNING)
    )
    root.addHandler(
        _make_file_handler(paths["errors"], logging.ERROR, logging.CRITICAL)
    )
    if debug:
        root.addHandler(
            _make_file_handler(paths["debug"], logging.DEBUG, logging.CRITICAL)
        )


def listener_process(
    log_queue: multiprocessing.Queue, log_file: str, debug: bool = False
):
    listener_configurer(log_file, debug=debug)
    while True:
        try:
            record = log_queue.get()
        except (EOFError, OSError, KeyboardInterrupt):
            break
        if record is None:  # Sentinel value to stop.
            break
        logger = logging.getLogger(record.name)
        logger.handle(record)


def _has_queue_handler(logger):
    return any(isinstance(h, logging.handlers.QueueHandler) for h in logger.handlers)


def _debug_enabled():
    return os.environ.get("PARALOGWIZARD_DEBUG", "").strip() in (
        "1",
        "true",
        "True",
        "yes",
    )


def worker_initializer(log_queue: multiprocessing.Queue):
    """
    Attach a QueueHandler on the worker root logger and clear the ParalogWizard
    logger so records only go through the shared queue (no per-worker file I/O).
    """
    level = logging.DEBUG if _debug_enabled() else logging.INFO
    queue_handler = logging.handlers.QueueHandler(log_queue)
    root = logging.getLogger()
    root.handlers = []
    root.addHandler(queue_handler)
    root.setLevel(level)

    pw = logging.getLogger("ParalogWizard")
    pw.handlers = []
    pw.propagate = True
    pw.setLevel(level)


def setup_logging():
    logger = logging.getLogger("ParalogWizard")
    level = logging.DEBUG if _debug_enabled() else logging.INFO

    # Worker / queue mode: never attach FileHandlers here (deadlocks / hangs).
    if _has_queue_handler(logger) or _has_queue_handler(logging.getLogger()):
        logger.handlers = [
            h for h in logger.handlers if isinstance(h, logging.handlers.QueueHandler)
        ]
        logger.propagate = True
        logger.setLevel(level)
        return logger

    if not logger.handlers:
        log_file = os.environ.get("PARALOGWIZARD_LOGFILE", "ParalogWizard.log")
        paths = log_paths_from_base(log_file, debug=_debug_enabled())
        logger.addHandler(
            _make_file_handler(paths["info"], logging.INFO, logging.WARNING)
        )
        logger.addHandler(
            _make_file_handler(paths["errors"], logging.ERROR, logging.CRITICAL)
        )
        if "debug" in paths:
            logger.addHandler(
                _make_file_handler(paths["debug"], logging.DEBUG, logging.CRITICAL)
            )

        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(level)
        stream_handler.setFormatter(
            logging.Formatter(
                "[%(asctime)s] [%(processName)s:%(process)d] [%(levelname)s] %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        logger.addHandler(stream_handler)
        logger.setLevel(level)
    else:
        logger.setLevel(level)
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


@contextmanager
def managed_pool(processes, log_queue):
    """
    multiprocessing.Pool that always terminates/joins, including on errors.

    Avoids hangs from close()+join() while workers are stuck, and from leaving
    workers alive after an exception (which can block the logging listener).
    """
    logger = logging.getLogger("ParalogWizard")
    pool = multiprocessing.Pool(
        processes=processes,
        initializer=worker_initializer,
        initargs=(log_queue,),
    )
    abort = False
    try:
        yield pool
    except BaseException:
        abort = True
        raise
    finally:
        if abort:
            logger.warning(
                "Terminating multiprocessing pool (%d worker(s)) due to error",
                processes,
            )
            try:
                pool.terminate()
            except Exception as exc:
                logger.debug("pool.terminate() failed: %s", exc)
        else:
            try:
                pool.close()
            except Exception as exc:
                logger.debug("pool.close() failed: %s", exc)
                try:
                    pool.terminate()
                except Exception:
                    pass
        try:
            pool.join()
            logger.debug("Multiprocessing pool shut down cleanly")
        except Exception as exc:
            logger.warning("pool.join() failed (%s); forcing terminate", exc)
            try:
                pool.terminate()
                pool.join()
            except Exception:
                pass


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
                fin = open(f, "rb")
                fout = gzip.open(gz_file, "wb")
                with fin, fout:
                    shutil.copyfileobj(fin, fout)
                os.remove(f)
            except Exception as e:
                logger.error("Error compressing %s: %s", f, e)
