#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs


class Event(object):
    def __init__(self):
        self.handlers = set()

    def on(self, handler):
        self.handlers.add(handler)
        return self

    def off(self, handler):
        try:
            self.handlers.remove(handler)
        except:
            raise ValueError("Handler is not handling this event, so cannot unhandle it.")
        return self

    def fire(self, *args, **kargs):
        for handler in self.handlers:
            handler(*args, **kargs)

    def count(self):
        return len(self.handlers)

    __iadd__ = on
    __isub__ = off
    __call__ = fire
    __len__ = count