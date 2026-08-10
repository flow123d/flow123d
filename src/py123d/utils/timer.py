#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import time


class Timer(object):
    """
    Class Timer measures elapsed time between tick and tock
    :type app_timer: Timer
    """
    app_timer = None

    def __init__(self, name=None):
        self.time = 0
        self.name = name
        self.duration = 0

    def tick(self):
        self.time = time.time()

    def tock(self):
        self.duration = time.time() - self.time

    def __enter__(self):
        self.tick()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.tock()
        return False

    def __repr__(self):
        if self.name is None:
            return "{:1.6f}".format(self.duration)
        return "{:s}: {:1.6f}".format(self.name, self.duration)


Timer.app_timer = Timer("app")