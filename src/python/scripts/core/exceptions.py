#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs


class ArgumentException(SystemExit):
    """
    Class ArgumentException for throwing invalid arguments exceptions
    """

    def __init__(self, exit_code, message=""):
        super(ArgumentException, self).__init__(exit_code)
        self.message = message