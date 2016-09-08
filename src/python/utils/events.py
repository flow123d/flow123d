#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import operator


class Event(object):
    """
    Class Event maintains a list of its dependents and notifies
    them automatically of any state changes
    """
    
    def __init__(self):
        self.handlers = dict()

    def on(self, handler):
        if type(handler) not in (tuple, list):
            handler = (handler, 1)

        self.handlers[handler[0]] = handler[1]
        return self

    def off(self, handler):
        try:
            del self.handlers[handler]
        except:
            raise ValueError("Handler is not handling this event, so cannot unhandle it.")
        return self

    def fire(self, *args, **kwargs):
        items = self.handlers.items()
        sorted_x = sorted(items, key=operator.itemgetter(1), reverse=True)
        for handler in sorted_x:
            handler[0](*args, **kwargs)

    def set_priority(self, handle, priority):
        self.handlers[handle] = priority

    def count(self):
        return len(self.handlers)

    __iadd__ = on
    __isub__ = off
    __call__ = fire
    __len__ = count
