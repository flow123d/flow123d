"""Module for GeoMop root context logging

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""
import logging
import os
import sys
import traceback


LOGGER_PREFIX = 'GeoMop-'


if 'APPDATA' in os.environ:
    __config_dir__ = os.path.join(os.environ['APPDATA'], 'GeoMop')
else:
    __config_dir__ = os.path.join(os.environ['HOME'], '.geomop')


def log_unhandled_exceptions(context_name, callback=None):
    """Initialize logging for unhandled exceptions.

    :param str context_name: name of the context (for example, ModelEditor)
    :param function callback: a function callback with three parameters (type, exception, tback)
       that is called after the exception is logged
    """
    logger = _get_root_context_logger(context_name)

    def log_excepthook(type_, exception, tback):
        """Set exception logging hook."""
        logger.critical('{0}: {1}\n  Traceback:\n{2}'
                        .format(type_.__name__, exception, ''.join(traceback.format_tb(tback))))

        # call the default handler
        sys.__excepthook__(type_, exception, tback)

        # enable customizable behaviour (i.e. show error dialog)
        if callable(callback):
            callback(type_, exception, tback)

    sys.excepthook = log_excepthook


def _get_root_context_logger(context_name):
    """Create a root logger for a context.

    :param str context_name: name of the context, ie. ModelEditor
    """""
    name = LOGGER_PREFIX + context_name
    format_ = '%(asctime)s %(levelname)s %(message)s'
    filename = os.path.join(__config_dir__, name + '.log.txt')
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(filename)
    formatter = logging.Formatter(format_)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger
