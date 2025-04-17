# -*- coding: utf-8 -*-
import logging


class Logging:
    """
    Control level of detail being printed.
    """

    plain_logger = logging.StreamHandler()
    plain_logger.setFormatter(logging.Formatter("%(message)s"))

    detail_logger = logging.StreamHandler()
    detail_logger.setFormatter(logging.Formatter("[%(name)8s] %(message)s"))

    loggers = (
        "batch",
        "bdd",
        "bins",
        "bkw",
        "dual",
        "gb",
        "guess",
        "repeat",
        "sweep",
        "usvp",
    )

    QUIET = logging.ERROR
    CRITICAL = logging.CRITICAL
    ERROR = logging.ERROR
    WARNING = logging.WARNING
    INFO = logging.INFO
    LEVEL0 = logging.INFO
    LEVEL1 = logging.INFO - 2
    LEVEL2 = logging.INFO - 4
    LEVEL3 = logging.INFO - 6
    LEVEL4 = logging.INFO - 8
    LEVEL5 = logging.DEBUG
    DEBUG = logging.DEBUG
    NOTSET = logging.NOTSET

    logging.getLogger("estimator").handlers = [plain_logger]
    logging.getLogger("estimator").setLevel(logging.INFO)

    for logger in loggers:
        logging.getLogger(logger).handlers = [detail_logger]
        logging.getLogger(logger).setLevel(logging.INFO)

    @staticmethod
    def set_level(lvl, loggers=None):
        """Set logging level

        The supported levels are:
        - `QUIET`,
        - `CRITICAL`,
        - `ERROR`,
        - `WARNING`,
        - `INFO`,
        - `LEVELX` with `X` ∈ [0,5]
        - `DEBUG` or
        - `NOTSET`

        :param lvl: see above
        :param loggers: one of `Logging.loggers`, if `None` all loggers are used.

        """
        if loggers is None:
            loggers = Logging.loggers

        for logger in loggers:
            logging.getLogger(logger).setLevel(lvl)

    @classmethod
    def log(cls, logger, level, msg, *args, **kwds):
        """
        Log `msg` indicating level.

        :param logger: logger to use
        :param level: level to log at
        :param msg: message to log

        """
        level = int(level)
        return logging.getLogger(logger).log(
            cls.INFO - 2 * level, f"{{{level}}} {msg}", *args, **kwds
        )

    @classmethod
    def print(cls, logger, level, msg, *args, **kwds):
        """
        Print `msg` maybe.

        :param logger: logger to use
        :param level: level to log at
        :param msg: message to log

        """
        level = int(level)
        if logging.getLogger(logger).level <= cls.INFO - 2 * level:
            print(msg, *args, **kwds)
