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
    """

    def __init__(self):
        self.parser = None
        self.arg_options, self.others, self.rest = None, None, None
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

    def run(self, parser, args=None, debug=False):
        """
        :type parser: utils.argparser.ArgParser
        """
        self.parser = parser
        self.debug = debug
        self.arg_options, self.others, self.rest = self.parser.parse(args)

        self._prepare()
        self._check_arguments()
        return self._run()