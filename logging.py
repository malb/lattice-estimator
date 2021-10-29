# -*- coding: utf-8 -*-
import logging


class Logging:
    """
    Control level of detail being printed.
    """

    plain_logger = logging.StreamHandler()
    plain_logger.setFormatter(logging.Formatter("%(message)s"))

    detail_logger = logging.StreamHandler()
    detail_logger.setFormatter(logging.Formatter("%(levelname)s:%(name)s: %(message)s"))

    logging.getLogger("estimator").handlers = [plain_logger]
    logging.getLogger("estimator").setLevel(logging.INFO)

    loggers = ("binsearch", "repeat", "guess")

    for logger in loggers:
        logging.getLogger(logger).handlers = [detail_logger]
        logging.getLogger(logger).setLevel(logging.INFO)

    CRITICAL = logging.CRITICAL
    ERROR = logging.ERROR
    WARNING = logging.WARNING
    INFO = logging.INFO
    DEBUG = logging.DEBUG
    NOTSET = logging.NOTSET

    @staticmethod
    def set_level(lvl, loggers=None):
        """Set logging level

        :param lvl: one of `CRITICAL`, `ERROR`, `WARNING`, `INFO`, `DEBUG`, `NOTSET`
        :param loggers: one of `Logging.loggers`, if `None` all loggers are used.

        """
        if loggers is None:
            loggers = Logging.loggers

        for logger in loggers:
            logging.getLogger(logger).setLevel(lvl)
