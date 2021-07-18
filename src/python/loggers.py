#!/bin/python3
# author: Jan Hybs
import enum
import sys
import os
import time
import re
from loguru import logger as _logger

__flow__ = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../' * 2)
)


class LogLevels(object):
    NORMAL = 'NORMAL'
    IMPORTANT = 'IMPORTANT'
    SUCCESS = 'SUCCESS'
    FAILED = 'FAILED'
    WARN = 'WARN'


class GlobalTimer(object):
    def __init__(self):
        self.start_time = time.time()

    @property
    def elapsed(self):
        return time.time() - self.start_time

    @property
    def elapsed5(self):
        return str(time.time() - self.start_time)[:5]


class PrintFormat(object):
    SEP = '-' * 60
    NOTHING = ' ' * 60
    _verbosity = None
    _tty = None
    logger = _logger

    class OutputVerbosity(enum.Enum):
        # will display entire output
        FULL = 'full'

        # will display only minimal info about the run, but when upon error
        # displays entire output
        SMART = 'smart'

        # will display minimal info
        MINIMAL = 'minimal'

        def __repr__(self):
            return self.value

    @classmethod
    def init_logger(cls, verbosity=None, logfile=None, prefix=''):
        if verbosity:
            cls.set_verbosity(verbosity)

        if logfile is None:
            logfile = os.path.join(
                __flow__, '.runtest.log'
            )

        _logger.configure(
            handlers=[
                dict(sink=sys.stdout, format='<level>{extra[prefix]}{message}</level>'),
                dict(sink=logfile),
            ],
            levels=[
                dict(name=LogLevels.NORMAL, no=60,  color=''),
                dict(name=LogLevels.IMPORTANT, no=60, color='<b>'),
                dict(name=LogLevels.SUCCESS, color='<b><g>'),
                dict(name=LogLevels.FAILED, no=60, color='<b><r>'),
                dict(name=LogLevels.WARN, no=60, color='<b><y>'),
            ]
        )
        cls.logger = _logger.bind(timer=GlobalTimer(), prefix=prefix)
        # increase depth for methods out and _write
        cls.logger._depth = 2
        return cls.logger

    @classmethod
    def verbosity(cls):
        return cls._verbosity if cls._verbosity else cls.OutputVerbosity.FULL

    @classmethod
    def tty(cls):
        if cls._tty is None:
            try:
                if sys.stdout.isatty():
                    cls._tty = True
                else:
                    cls._tty = False
            except:
                cls._tty = False
        return cls._tty

    @classmethod
    def set_verbosity(cls, verbosity):
        cls._verbosity = verbosity

    def __init__(self, raw=False, level=0, indent='  ', default=LogLevels.NORMAL, ansi=False):
        self.raw = raw
        self.level = level
        self.indent = indent
        self.ansi = ansi,
        self.default = default

    def indent_str(self):
        return self.level * self.indent

    def opt(self, **kwargs):
        opts = dict(
            raw=self.raw, level=self.level, ansi=self.ansi,
            indent=self.indent, default=self.default
        )
        opts.update(**kwargs)

        return PrintFormat(**opts)

    def sep(self):
        return self._write(self.SEP, self.default)

    def _write(self, msg, lvl, *args, **kwargs):
        indent_msg = self._get_msg(msg, *args, **kwargs)
        if self.ansi:
            self.logger.opt(ansi=self.ansi).log(lvl, indent_msg)
        else:
            self.logger.log(lvl, indent_msg)
        return self

    def _get_msg(self, msg, *args, **kwargs):
        fmt_msg = str(msg)
        if not self.raw:
            fmt_msg = fmt_msg.format(*args, **kwargs)
        fmt_msg = '' + self.indent_str() + fmt_msg
        # escape brackets as loguru tries to use them for formatting.
        fmt_msg = re.sub(r"(\\?</?((?:[fb]g\s)?[^<>\s]*)>)", r"\\\1", fmt_msg)
        return fmt_msg

    def __enter__(self):
        self.level += 1
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.level -= 1
        return False

    def out(self, msg, *args, **kwargs):
        return self._write(msg, LogLevels.NORMAL, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        return self._write(msg, LogLevels.FAILED, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        return self._write(msg, LogLevels.WARN, *args, **kwargs)

    def success(self, msg, *args, **kwargs):
        return self._write(msg, LogLevels.SUCCESS, *args, **kwargs)

    def important(self, msg, *args, **kwargs):
        return self._write(msg, LogLevels.IMPORTANT, *args, **kwargs)

    def stream(self, iterable):
        for l in iterable:
            self._write(l, self.default)
        return self

    def rewrite(self, msg='\n', *args, **kwargs):
        self.logger.opt(raw=True).log(self.default, '\r' + self.NOTHING + '\r')
        self.logger.opt(raw=True).log(self.default, '%-60s' % self._get_msg(msg, *args, **kwargs) + '\r')
        return self

    def finish_rewrite(self):
        self._write('', self.default)
        return self


printf = PrintFormat()
