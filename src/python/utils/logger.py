# encoding: utf-8
# author:   Jan Hybs
import logging, traceback
import datetime

logging.basicConfig(
    filename='python.log',
    level=logging.INFO,
    format='%(asctime)s %(name)s %(levelname)-4s: %(message)s'
)


class Logger(object):
    def __init__(self, name):
        self.logger = logging.getLogger(name)

        # add console log
        if __debug__:
            stream = logging.StreamHandler()
            stream.setLevel(logging.INFO)
            self.logger.addHandler(stream)

        with open('python.log', 'a+') as fp:
            fp.write('-' * 110)
            fp.write('\n')
        self.info("{:%d-%m-%Y %H:%M:%S}".format(datetime.datetime.now()))

    def _log_traceback(self, method):
        tb = [line.strip() for line in traceback.format_stack()]
        method('Traceback:\n' + '\n'.join(tb))

    def error(self, msg, *args, **kwargs):
        self.logger.error(msg, *args, **kwargs)
        self._log_traceback(self.logger.error)

    def exception(self, msg, exception, *args, **kwargs):
        self.logger.exception(msg, exc_info=exception, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self.logger.info(msg, *args, **kwargs)

    def debug(self, msg, *args, **kwargs):
        self.logger.debug(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self.logger.warning(msg, *args, **kwargs)