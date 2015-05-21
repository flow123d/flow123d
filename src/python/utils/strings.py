# encoding: utf-8
# author:   Jan Hybs

"""
Simple module which provides string/list utilities
"""


def join_iterable (iterable, prefix="", suffix="", separator=",", padding=None, extra_space=2):
    """Joins given iterable object with extra surrounding options"""
    result = ""
    result += prefix

    size = len (iterable)
    for i in range (size):
        if padding is None:
            result += str (iterable[i])
        else:
            result += str (iterable[i]).center (padding[i] + extra_space)
        if i < size - 1:
            result += (suffix + separator + prefix)
    result += suffix

    return result
