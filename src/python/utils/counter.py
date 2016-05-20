#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from scripts.core.base import Printer


class ProgressCounter(object):
    def __init__(self, fmt):
        self.i = 0
        self.fmt = fmt
        self.printer = Printer(Printer.LEVEL_DBG)

    def reset(self):
        self.i = 0

    def next(self, attributes):
        self.i += 1
        self.printer.dbg(self.fmt.format(
            self.i, **attributes
        ))
