#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

"""
Simple module which provides string/list utilities
"""


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


def format_n_lines(text, success=True, n_lines=0):
    if success:
        return format_n_lines_(text, **format_n_lines_successful, n_lines=n_lines)
    return format_n_lines_(text, **format_n_lines_error, n_lines=n_lines)


def format_n_lines_(text, line_prefix='## ', line_suffix='',
                   first_line='#' * 60, last_line='#' * 60,
                   empty="<file is empty>", n_lines=0):
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

    if first_line:
        yield first_line

    for line in text:
        yield line_prefix + line + line_suffix

    if last_line:
        yield last_line


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


def replace_placeholders(template, **kwargs):
    """
    Method will replace placeholders with actual values
    :param template:
    :param kwargs:
    :return:
    """
    result = str(template)
    _format_ = kwargs.get('_format_', '$${}$$')

    for key, value in list(kwargs.items()):
        result = result.replace(_format_.format(key), value)
    return result


def format_dict(d, align_keys=True, indent=0, eq=' = ', sort=True):
    """
    :type sort: bool
    :type eq: str
    :type indent: int
    :type align_keys: bool
    :type d: dict
    """
    p = '    '
    result = list()
    max_width = max(len(k) for k in list(d.keys()))
    keys = sorted(d.keys()) if sort else list(d.keys())
    for k in keys:
        v = d.get(k)
        i = indent * p
        if align_keys:
            pad = max_width - len(str(k))
            start = i + str(k) + ' '*pad + eq
        else:
            start = i + str(k) + eq

        if type(v) is list:
            if not v:
                result.append(start + '[]')
            else:
                if len(v) == 1 and len(str(v[0])) < 16:
                    result.append(start + '[' + str(v[0]) + ']')
                else:
                    result.append(start)
                    for lv in v:
                        result.append(i + p + '- ' + str(lv))
        else:
            result.append(start + str(v))

    return '\n'.join(result)
