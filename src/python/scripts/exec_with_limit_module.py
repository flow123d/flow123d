#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from __future__ import absolute_import
from scripts.core.base import Printer, Paths

from scripts.execs.monitor import PyProcess
from scripts.execs.test_executor import BinExecutor


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
    process_monitor = PyProcess(executor, batch_mode=options.batch)

    # set limits
    process_monitor.limit_monitor.time_limit = options.time_limit
    process_monitor.limit_monitor.memory_limit = options.memory_limit

    # turn on output
    if options.batch:
        process_monitor.info_monitor.stdout_stderr = None
    else:
        process_monitor.info_monitor.stdout_stderr = Paths.temp_file('exec-limit.log')

    # start process
    printer.key('-' * 60)
    process_monitor.start()