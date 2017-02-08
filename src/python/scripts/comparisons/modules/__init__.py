#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from scripts.core.threads import ExtendedThread, BrokenProcess


# ----------------------------------------------


class ExecComparison(object):
    """
    Class ExecComparison is interface for comparison which
     compares two files using external binary file
    """

    def get_command(self, reference_filepath, other_filepath, **kwargs):
        """
        Method will create argument list which will be passed to process.Popen

        :param reference_filepath: reference file
        :param other_filepath: other file which was produces using flow123d binary
        :param kwargs: additional arguments which were specified in yaml config
        :rtype: list
        :return: list of strings where first item is binary filepath (or command)
         other elements in list are arguments passed to command
        """
        raise NotImplementedError("Method must be implemented")


class InPlaceComparison(ExtendedThread):
    """
    Class InPlaceComparison is interface for comparison which
    will compare files in this python execution.
    """

    def __init__(self):
        super(InPlaceComparison, self).__init__('in-place-comparison')
        self.reference_filepath = None
        self.other_filepath = None
        self.kwargs = None
        self.command = []
        self.process = None
        self.broken = False
        self.output = None
        self.exception = ""
        self.escaped_command = "<InPlaceComparison> {}".format(self.__class__)

    def prepare(self, reference_filepath, other_filepath, **kwargs):
        """
        Method will just save argument so later there will be passed to
        compare method
        :param reference_filepath:
        :param other_filepath:
        :param kwargs:
        :return:
        """
        self.reference_filepath = reference_filepath
        self.other_filepath = other_filepath
        self.kwargs = kwargs

    def compare(self, reference_filepath, other_filepath, **kwargs):
        """
        Method will compare two given files.

        :param reference_filepath: reference file
        :param other_filepath: other file which was produces using flow123d binary
        :param kwargs: additional arguments which were specified in yaml config
        :rtype: int
        :return: return code value (int) where 0 means success and any other value
        means failure (even None).
        """
        raise NotImplementedError("Method must be implemented")

    def _run(self):
        self.process = BrokenProcess()
        self.returncode = self.compare(
            self.reference_filepath,
            self.other_filepath,
            **self.kwargs
        )