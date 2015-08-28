# encoding: utf-8
# author:   Jan Hybs

"""
Simple module which provides string/list utilities
"""


def join_iterable (iterable, prefix="", suffix="", separator=",", padding=None, extra_space=2):
    """
    Joins given iterable object with extra surrounding options
    :param iterable: iterable object
    :param prefix: string part which will be added before iterable item value
    :param suffix: string part which will be added after iterable item value
    :param separator: string which separates items
    :param padding: if list is specified, item string representation
     will be padded with white space. Each entry must contain number - desired item length
    :param extra_space: number of chars to be extra added to padding
    :return:
    """
    size = len (iterable)
    if not size:
        return ""

    result = ""
    result += prefix

    for i in range (size):
        if padding is None:
            result += str (iterable[i])
        else:
            result += str (iterable[i]).center (padding[i] + extra_space)
        if i < size - 1:
            result += (suffix + separator + prefix)
    result += suffix

    return result
