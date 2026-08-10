#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import sys


class RC(object):
    """
    Class RC represents return code value. Class implements
    standard six rich comparison operators.
    """

    MAX_INT = +sys.maxsize
    MIN_INT = -sys.maxsize

    def __init__(self, value=None):
        if value is None:
            self.rc = None
            self.value = None
            self.reversed = False
        elif isinstance(value, self.__class__):
            self.rc = value.rc
            self.value = value.value
            self.reversed = value.reversed
        else:
            self.rc = value
            self.value = abs(int(value))
            self.reversed = False

    def get(self, default=0):
        if self.value is None:
            return default
        return self.value

    @property
    def failed(self):
        return self not in (0, None)

    @property
    def succeeded(self):
        return self in (0, None)

    def __call__(self, *args, **kwargs):
        """
        :rtype: int
        """
        return self.value

    def __gt__(self, other):
        a = self.get(self.MIN_INT)
        if isinstance(other, self.__class__):
            b = other.get(self.MAX_INT)
            return a > b
        if other is None:
            return False
        return a > abs(int(other))

    def __ge__(self, other):
        a = self.get(self.MIN_INT)
        if isinstance(other, self.__class__):
            b = other.get(self.MAX_INT)
            return a >= b
        if other is None:
            return False
        return a >= abs(int(other))

    def __lt__(self, other):
        a = self.get(self.MAX_INT)
        if isinstance(other, self.__class__):
            b = other.get(self.MIN_INT)
            return a < b
        if other is None:
            return False
        return a > abs(int(other))

    def __le__(self, other):
        a = self.get(self.MAX_INT)
        if isinstance(other, self.__class__):
            b = other.get(self.MIN_INT)
            return a <= b
        if other is None:
            return False
        return a >= abs(int(other))

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.value == other.value
        return self.value == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return '{0}'.format(self.value)

RC_NONE = RC()
RC_OK = RC(0)
RC_BROKEN = RC(666)