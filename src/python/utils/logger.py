#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from __future__ import absolute_import
import logging, traceback
import datetime
import sys, os


class Logger(object):
    _global_logger = None

    @staticmethod
    def instance():
        """
        :rtype : Logger
        """
        if Logger._global_logger is None:
            log_level = 'warning'
            try:
                for index in range(len(sys.argv)):
                    arg = str(sys.argv[index])
                    if arg.startswith('--log'):
                        if arg == '--log':
                            log_level = sys.argv[index + 1].strip()
                        else:
                            log_level = arg.replace('--log=', '').strip()
            except Exception as e:
                pass
            log_level = getattr(logging, log_level.upper())

            # set global log level
            fmt = logging.Formatter('%(asctime)s %(name)-15s %(levelname)s: %(message)s')
            logging.root.setLevel(log_level)
            Logger._global_logger = Logger('ROOT', log_level, fmt)

        return Logger._global_logger

    def __init__(self, name, level=logging.INFO, fmt=None):
        self.logger = logging.getLogger(name)

        f = os.path.join(os.getcwd(), 'python.log')
        stream_logger = logging.StreamHandler()
        stream_logger.setLevel(level)
        stream_logger.setFormatter(fmt)

        file_logger = logging.FileHandler(f)
        file_logger.setLevel(level)
        file_logger.setFormatter(fmt)

        self.logger.addHandler(stream_logger)
        self.logger.addHandler(file_logger)

        # add empty lines if file contains some previous logs
        if os.stat(f).st_size != 0:
            with open(f, 'a+') as fp:
                fp.write('\n' * 3)
        # start logging
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