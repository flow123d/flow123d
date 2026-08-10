#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import logging
import traceback
import datetime
import sys
import os
import io


class Logger(object):
    """
    Class Logger is basic logger class for logging messages (console and file)
    """

    _global_logger = None
    level = 'warning'

    @staticmethod
    def instance():
        """
        :rtype : Logger
        """
        if Logger._global_logger is None:
            log_level = Logger.level
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
            fmt = logging.Formatter('%(asctime)s %(levelname)s: %(message)s', '%H:%M:%S')
            logging.root.setLevel(log_level)
            Logger._global_logger = Logger('ROOT', log_level, fmt)

        return Logger._global_logger

    @classmethod
    def __loc(cls):
        frame = sys._getframe(2)
        loc = frame.f_code.co_filename
        line = frame.f_lineno
        # return 'File "%s", line %d ' % (loc, line)
        loc = loc.split('/')
        if 'python' in loc:
            i = loc.index('python')
            return '/'.join(loc[i-1:]) + ':' + str(frame.f_lineno)
        return '/'.join(loc) + ':' + str(frame.f_lineno)


    def __init__(self, name, level=logging.INFO, fmt=None):
        self.logger = logging.getLogger(name)

        stream_logger = logging.StreamHandler()
        stream_logger.setLevel(level)
        stream_logger.setFormatter(fmt)
        self.logger.addHandler(stream_logger)

        self.memory_stream = None

        # f = os.path.join(os.getcwd(), 'python.log')
        # file_logger = logging.FileHandler(f)
        # file_logger.setLevel(level)
        # file_logger.setFormatter(fmt)
        # self.logger.addHandler(file_logger)

        # add empty lines if file contains some previous logs
        # if os.stat(f).st_size != 0:
        #     with open(f, 'a+') as fp:
        #         fp.write('\n' * 3)
        # start logging
        self.info("{:%d-%m-%Y %H:%M:%S}".format(datetime.datetime.now()))

    def _log_traceback(self, method):
        tb = [line.strip() for line in traceback.format_stack()]
        method('Traceback:\n' + '\n'.join(tb))
    
    def _get_msg(self, loc, msg):
        loc, msg = str(loc), str(msg)
        lloc, lmsg = len(loc), len(msg)
        msg = msg.replace('\n', '\n' + ' ' * 61)
        return ('{loc:45s} {msg:s}').format(loc=loc, msg=msg)

    def error(self, msg, *args, **kwargs):
        msg = self._get_msg(self.__loc(), msg)
        self.logger.error(msg, *args, **kwargs)
        self._log_traceback(self.logger.error)

    def exception(self, msg, exception, *args, **kwargs):
        self.logger.exception(msg, exc_info=exception, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        msg = self._get_msg(self.__loc(), msg)
        self.logger.info(msg, *args, **kwargs)

    def debug(self, msg, *args, **kwargs):
        msg = self._get_msg(self.__loc(), msg)
        self.logger.debug(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        msg = self._get_msg(self.__loc(), msg)
        self.logger.warning(msg, *args, **kwargs)
    
    def add_memory_handler(self):
        # Create an in-memory stream
        self.memory_stream = io.StringIO()

        # Create a StreamHandler that logs to the in-memory stream
        stream_handler = logging.StreamHandler(log_stream)
        stream_logger.setLevel(level)
        stream_logger.setFormatter(fmt)
        logger.addHandler(stream_handler)
    
    def get_memory_stream(self):
        """
        If the momory stream was added, read the log and return it.
        """
        if self.memory_stream is not None:
            return self.memory_stream
