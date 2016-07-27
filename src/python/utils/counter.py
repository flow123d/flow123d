#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import threading
import time
import datetime

from scripts.core.base import Printer


class ProgressCounter(object):
    def __init__(self, fmt):
        self.i = 0
        self.fmt = fmt

    def reset(self):
        self.i = 0

    def next(self, attributes):
        self.i += 1
        Printer.out(self.fmt.format(
            self.i, **attributes
        ))


class ProgressTime(object):
    """
    :type thread : threading.Thread
    """
    def __init__(self, format='{}', period=0.1, dynamic=True):
        self.format = format
        self.period = period

        self.thread = None
        self.start_time = None
        self.running = False
        self.active = True
        self.format_args = dict(s=self)
        self.dynamic = dynamic

    @property
    def elapsed(self):
        t = time.time() - self.start_time
        d = datetime.timedelta(seconds=int(t))
        ms = '{:1.3f}'.format(t - int(t))[2:]
        return '{}:{}'.format(d, ms)

    def update(self):
        self.start_time = self.start_time or time.time()
        Printer.dyn(self.format, self.elapsed, **self.format_args)
        if not self.dynamic:
            Printer.out()

    def __enter__(self):
        if not self.active:
            return self

        def target():
            while self.running:
                self.update()
                time.sleep(self.period)

        self.start_time = time.time()
        self.running = True
        self.thread = threading.Thread(target=target)
        self.thread.start()
        return self

    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        if not self.active:
            return False

        self.running = False
        self.update()

        # print \n if in dynamic mode (ending)
        if self.dynamic:
            Printer.out()
        return False

    stop = __exit__
