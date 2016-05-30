#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from __future__ import absolute_import
from scripts.core.base import Printer, Paths
from scripts.core.threads import BinExecutor, PyPy


def do_work(parser):
    """
    :type parser: utils.argparser.ArgParser
    """

    # parse arguments
    options, others, rest = parser.parse()
    printer = Printer(Printer.LEVEL_KEY)

    # check commands
    if not rest:
        parser.exit_usage('No command specified!')

    # check limits (at least one limit must be set)
    if (options.time_limit, options.memory_limit) == (None, None):
        parser.exit_usage('No limits specified!')

    # prepare executor
    executor = BinExecutor(rest)
    pypy = PyPy(executor, progress=not options.batch)

    # set limits
    pypy.error_monitor.message = None
    pypy.limit_monitor.time_limit = options.time_limit
    pypy.limit_monitor.memory_limit = options.memory_limit

    # turn on output
    if options.batch:
        pypy.info_monitor.stdout_stderr = None
    else:
        pypy.info_monitor.stdout_stderr = Paths.temp_file('exec-limit.log')

    # start process
    printer.line()
    pypy.start()