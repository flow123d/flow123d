#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from py123d.loggers import printf
from py123d.scripts.core.base import Printer, Paths, PathFormat
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

    def _prepare(self):
        # configure path
        Paths.format = PathFormat.ABSOLUTE
        if self.arg_options.root:
            Paths.init(self.arg_options.root)
            printf.out(
                '[INFO]    | Forcing Flow123d root location to "{}"\n'
                '            Flow123d binary is:               "{}"',
                self.arg_options.root,
                Paths.flow123d()
            ).sep()

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
