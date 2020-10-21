#!/usr/bin/env python

"""Make logging fancy again"""

from typing import Dict, List, Optional

__all__: List[str] = [
    'StrippedFormatter',
    'ColoredFormatter',
    'fmttime'
]

import os
import sys
import time
import logging
import warnings

LOG_FORMAT = '%(asctime)s %(levelname)s:\t%(message)s' # type: str
DATE_FORMAT = '%Y-%m-%d %H:%M:%S' # type: str
LOG_LEVELS: Dict[str, int] = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

class StrippedFormatter(logging.Formatter):
    """A formatter where all ANSI formatting is removed"""

    def __init__(self, *args, **kwargs) -> None:
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record: logging.LogRecord) -> str:
        """Strip ANSI formatting from log messages"""
        message: str = logging.Formatter.format(self, record)
        while True:
            #   In Python, '\x1b' == '\033', so both codes for ANSI formatting are covered
            start: int = message.find('\x1b')
            #   If no ASI formatting is found break
            if start == -1:
                break
            #   Find the first 'm' after the ANSI code start
            #   and remove everything between and including
            #   the ANSI code start and the 'm'
            m_pos: int = message.find('m', start)
            message: str = message[:start] + message[m_pos + 1:]
        return message


class ColoredFormatter(logging.Formatter):
    """A colorized formatter for logging"""

    _colors: Dict[int, str] = {
        50: '\x1b[1m\x1b[31m', # CRITICAL: bold red
        40: '\x1b[31m', # ERROR: red
        30: '\x1b[33m', # WARNING: yellow
        20: '\x1b[32m', # INFO: green
        10: '\x1b[36m' # DEBUG: cyan
    }

    _default: str = '\x1b[0m' # Anything else: reset

    def __init__(self, *args, **kwargs) -> None:
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record: logging.LogRecord) -> str:
        """Colorize log messages"""
        message: str = logging.Formatter.format(self, record)
        if sys.platform not in ('win32', 'cygwin'):
            color_level: int = min(self._colors.keys(), key=lambda level: abs(level - record.levelno))
            color_level: int = min((color_level, record.levelno))
            color: str = self._colors.get(color_level, self._default)
            message: str = color + message + self._default
        return message


def configure_logging(level: str='info', logstream: bool=True, logfile: Optional[str]=None) -> None:
    """Configure a logger"""
    logging.getLogger().handlers: List = list()
    verbosity: int = LOG_LEVELS.get(level, 'info')
    logging.basicConfig(
        stream=open(os.devnull, 'w'),
        level=verbosity
    )
    if verbosity == LOG_LEVELS['debug']:
        logging.captureWarnings(True)
    else:
        warnings.filterwarnings('ignore')
    if logstream:
        streamhandler: logging.StreamHandler = logging.StreamHandler()
        streamhandler.setFormatter(ColoredFormatter(fmt=LOG_FORMAT, datefmt=DATE_FORMAT))
        logging.getLogger().addHandler(streamhandler)
    if isinstance(logfile, str):
        if os.path.splitext(logfile)[-1] != '.log':
            logfile += '.log'
        filehandler: logging.FileHandler = logging.FileHandler(filename=logfile, mode='w')
        filehandler.setFormatter(StrippedFormatter(fmt=LOG_FORMAT, datefmt=DATE_FORMAT))
        logging.getLogger().addHandler(filehandler)
    logging.debug("Log level set to '%s'", level)


def fmttime(start: float, n: int=3) -> float:
    return round(time.time() - start, 3)
