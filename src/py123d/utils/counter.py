#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import threading
import time
import datetime

from py123d.loggers import printf


class ProgressCounter(object):
    """
    Class ProgressCounter is simple printer-like class which count to specific target
    """

    def __init__(self, fmt='{:02d}', printer=printf):
        self.i = 0
        self.fmt = fmt
        self.printer = printer

    def reset(self):
        self.i = 0

    def next(self, attributes):
        self.i += 1
        if self.printer:
            self.printer.out(self.fmt.format(
                self.i, **attributes
            ))


class ProgressTime(object):
    """
    Class ProgressTime will measure time for specific scope
    and prints elapsed time
    :type thread : threading.Thread
    """

    def __init__(self, period=1.0, handler=None):
        self.period = period

        self.thread = None
        self.start_time = None
        self.running = False
        self.active = True
        self.format_args = dict(s=self)
        self.handler = handler

    @property
    def elapsed(self):
        t = time.time() - self.start_time
        d = datetime.timedelta(seconds=int(t))
        ms = '{:1.3f}'.format(t - int(t))[2:]
        return '{}:{}'.format(d, ms)

    def update(self):
        self.start_time = self.start_time or time.time()

        if self.handler:
            self.handler(self.elapsed)

    def __enter__(self):
        if not self.active:
            return self

        def target():
            while self.running:
                self.update()
                time.sleep(self.period)

        self.start_time = time.time()
        self.running = True
        # self.thread = threading.Thread(target=target)
        # self.thread.start()
        return self

    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        if not self.active:
            return False

        self.running = False
        self.update()
        return False

    stop = __exit__
