#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

"""
Simple module which provides string/list utilities
"""

from __future__ import absolute_import
from scripts.core.base import Printer

format_n_lines_successful = dict(
    line_prefix='== ',
    line_suffix='',
    first_line='=' * 60,
    last_line='=' * 60,
)

format_n_lines_error = dict(
    line_prefix='## ',
    line_suffix='',
    first_line='#' * 60,
    last_line='#' * 60,
)


def format_n_lines(text, success=True):
    if success:
        return format_n_lines_(text, **format_n_lines_successful)
    return format_n_lines_(text, **format_n_lines_error)


def format_n_lines_(text, line_prefix='## ', line_suffix='',
                   first_line='#' * 60, last_line='#' * 60,
                   empty="<file is empty>"):
    """
    Format given lines and adds prefix to them
    :param text:
    :param line_prefix:
    :param line_suffix:
    :param first_line:
    :param last_line:
    :param empty:
    :return:
    """

    n_lines = 0 if not Printer.batched.is_muted() else 20
    indent = Printer.indent()

    # empty output
    if text is None or not text:
        text = '{:-^54s}'.format(empty)
        line_suffix = line_prefix[::-1]

    # ensure we have list or iterable
    text = text.splitlines() if type(text) is str else text

    # positive n_lines (first n lines)
    if n_lines > 0:
        text = text[:n_lines]

    # negative n_lines (last n lines)
    elif n_lines < 0:
        text = text[n_lines:]

    # otherwise all lines (0)

    result = list()
    if first_line:
        result.append(indent + first_line)

    for line in text:
        result.append(indent + line_prefix + line + line_suffix)

    if last_line:
        result.append(indent + last_line)

    return '\n'.join(result)


def join_iterable(iterable, prefix="", suffix="", separator=",", padding=None, extra_space=2):
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
    size = len(iterable)
    if not size:
        return ""

    result = ""
    result += prefix

    for i in range(size):
        if padding is None:
            result += str(iterable[i])
        else:
            result += str(iterable[i]).center(padding[i] + extra_space)
        if i < size - 1:
            result += (suffix + separator + prefix)
    result += suffix

    return result
