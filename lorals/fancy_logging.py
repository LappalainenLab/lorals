#!/usr/bin/env python

"""Make logging fancy again"""

from __future__ import division
from __future__ import print_function

__all__ = [ # type: List[str, ...]
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
LOG_LEVELS = { # type: Dict[str, int]
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

class StrippedFormatter(logging.Formatter):
    """A formatter where all ANSI formatting is removed"""

    def __init__(self, *args, **kwargs):
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record): # type: (logging.LogRecord) -> str
        """Strip ANSI formatting from log messages"""
        message = logging.Formatter.format(self, record) # type: str
        while True:
            #   In Python, '\x1b' == '\033', so both codes for ANSI formatting are covered
            start = message.find('\x1b') # type: int
            #   If no ASI formatting is found break
            if start == -1:
                break
            #   Find the first 'm' after the ANSI code start
            #   and remove everything between and including
            #   the ANSI code start and the 'm'
            m_pos = message.find('m', start) # type: int
            message = message[:start] + message[m_pos + 1:]
        return message


class ColoredFormatter(logging.Formatter):
    """A colorized formatter for logging"""

    _colors = { # type: Dict[int, str]
        50: '\x1b[1m\x1b[31m', # CRITICAL: bold red
        40: '\x1b[31m', # ERROR: red
        30: '\x1b[33m', # WARNING: yellow
        20: '\x1b[32m', # INFO: green
        10: '\x1b[36m' # DEBUG: cyan
    }

    _default = '\x1b[0m' # Anything else: reset

    def __init__(self, *args, **kwargs):
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record): # type: (logging.LogRecord) -> str
        """Colorize log messages"""
        message = logging.Formatter.format(self, record) # type: str
        if sys.platform not in ('win32', 'cygwin'):
            color_level = min(self._colors.keys(), key=lambda level: abs(level - record.levelno)) # type: int
            color_level = min((color_level, record.levelno)) # type: int
            color = self._colors.get(color_level, self._default) # type: str
            message = color + message + self._default # type: str
        return message


def configure_logging(level='info', logstream=True, logfile=None): # type: (str, bool, Optional[str]) -> None
    """Configure a logger"""
    logging.getLogger().handlers = list() # type: List
    verbosity = LOG_LEVELS.get(level, 'info') # type: int
    logging.basicConfig(
        stream=open(os.devnull, 'w'),
        level=verbosity
    )
    if verbosity == LOG_LEVELS['debug']:
        logging.captureWarnings(True)
    else:
        warnings.filterwarnings('ignore')
    if logstream:
        streamhandler = logging.StreamHandler() # type: logging.StreamHandler
        streamhandler.setFormatter(ColoredFormatter(fmt=LOG_FORMAT, datefmt=DATE_FORMAT))
        logging.getLogger().addHandler(streamhandler)
    if isinstance(logfile, str):
        if os.path.splitext(logfile)[-1] != '.log':
            logfile += '.log'
        filehandler = logging.FileHandler(filename=logfile, mode='w') # type: logging.FileHandler
        filehandler.setFormatter(StrippedFormatter(fmt=LOG_FORMAT, datefmt=DATE_FORMAT))
        logging.getLogger().addHandler(filehandler)
    logging.debug("Log level set to '%s'", level)


def fmttime(start, n=3): # type: (float, int) -> float
    return round(time.time() - start, 3)
