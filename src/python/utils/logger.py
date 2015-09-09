# encoding: utf-8
# author:   Jan Hybs
import logging, traceback

logging.basicConfig(level=logging.INFO)


class Logger(object):
    def __init__(self, name):
        self.logger = logging.getLogger(name)

    def _log_traceback(self, method):
        tb = [line.strip() for line in traceback.format_stack()]
        method('Traceback:\n' + '\n'.join(tb))

    def error(self, msg, *args, **kwargs):
        self.logger.error(msg, *args, **kwargs)
        self._log_traceback(self.logger.error)

    def exception(self, msg, exception, *args, **kwargs):
        self.logger.exception(msg, exc_info=exception, *args, **kwargs)