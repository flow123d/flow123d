#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from scripts.core.base import Paths, IO
from scripts.core.threads import PyPy
from scripts.core.execution import BinExecutor, OutputMode
from scripts.script_module import ScriptModule
# ----------------------------------------------


class ModuleExecWithLimit(ScriptModule):
    """
    Class ModuleExecWithLimit is backend for script exec_with_limit.py
    """

    def _check_arguments(self):
        """
        Arguments additional check
        """

        # check commands
        if not self.rest:
            self.parser.exit_usage('No command specified!', exit_code=1)

        # check limits (at least one limit must be set)
        if self.arg_options.time_limit is None and self.arg_options.memory_limit is None:
            self.parser.exit_usage('No limits specified!', exit_code=2)

    def _run(self):
        """
        Run method for this module
        """

        # prepare executor
        progress = not self.arg_options.batch
        executor = BinExecutor(self.rest)
        pypy = PyPy(executor, progress=progress)
        n_lines = 0 if self.arg_options.batch else 10

        # set up streams
        log_file = Paths.temp_file('exec-limit-{date}-{time}-{rnd}.log')
        pypy.executor.output = OutputMode.variable_output()
        pypy.full_output = log_file

        # set limits
        pypy.limit_monitor.time_limit = self.arg_options.time_limit
        pypy.limit_monitor.memory_limit = self.arg_options.memory_limit

        # save output to file
        pypy.output_monitor.log_file = log_file

        # start process
        pypy.start()
        pypy.join()

        return pypy


def do_work(parser, args=None):
    """
    Main method which invokes ModuleExecWithLimit
    :rtype: scripts.core.threads.PyPy
    :type args: list
    :type parser: utils.argparser.ArgParser
    """
    module = ModuleExecWithLimit()
    return module.run(parser, args, False)
