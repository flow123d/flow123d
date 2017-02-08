#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from scripts.core.base import Printer, Paths, PathFormat
# ----------------------------------------------


class ScriptModule(object):
    """
    ScriptModule is abstract class for all script modules.
    Children must implement or override method such as _check_arguments and _run

    :type arg_options : utils.argparser.RuntestArgs or utils.argparser.ExecWithLimitArgs or utils.argparser.ExecParallelArgs
    """

    def __init__(self, arg_options):
        """
        :type arg_options: utils.argparser.RuntestArgs or utils.argparser.ExecWithLimitArgs or utils.argparser.ExecParallelArgs
        """
        self.arg_options = arg_options
        self.debug = False
        self.progress, self.batch = None, None

    def _prepare(self):
        # configure printer
        Printer.batch_output = self.arg_options.batch
        Printer.dynamic_output = not self.arg_options.batch

        self.progress = Printer.dynamic_output
        self.batch = Printer.batch_output

        # configure path
        Paths.format = PathFormat.ABSOLUTE
        if self.arg_options.root:
            Paths.init(self.arg_options.root)

    def _check_arguments(self):
        pass

    def _run(self):
        pass

    def run(self, debug=False):
        """
        :type parser: utils.argparser.RuntestArgs
        """
        self.debug = debug

        self._prepare()
        self._check_arguments()
        return self._run()
