#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import time
import sys


def check_modules(*modules):
    """
    Check whether all modules can be imported, return True if everything is ok
    :param modules:
    :return: bool
    """
    import imp

    log = list()
    error = False
    for module in modules:
        try:
            imp.find_module(module)
            log.append(' - found module "{module:s}"'.format(**locals()))
        except ImportError:
            log.append(' - missing module "{module:s}"'.format(**locals()))
            error = True

    # on error exit application
    if error:
        print('Could not find one of the modules!')
        print('\n'.join(log))
        print('\n')
        print('Sys path:\n  ' + '\n  '.join(sys.path))
        return False

    return True


def ensure_iterable(o):
    """
    Method ensure that given object is iterable(list or tuple)
    :param o: tested object
    :return: list or tuple
    """
    return [o] if type(o) not in (list, tuple, set) else o


def apply_to_all(lst, mtd, *args, **kwargs):
    """
    Method will call mtd on every object in lst
    :param lst: list of objects
    :param mtd: string name of callable method
    :param args:
    :param kwargs:
    """
    return [getattr(x, mtd)(*args, **kwargs) for x in lst]


def wait_for(obj, property, period=0.1, max_wait=5):
    """
    Method will wait until property prop on object o exists
    and is not None
    Careful use since if value exists and value is None, can cause thread block
    :param obj:
    :param property:
    :param period:
    :param max_wait:
    :return:
    """
    wait = 0
    while True:
        attr = getattr(obj, property, None)
        if attr is not None:
            return attr

        time.sleep(period)
        wait += period
        if wait > max_wait:
            print 'done waiting on ', property
            return None


def justify(text, max_length=60, max_spaces=2):
    """
    Method will justify text inserting maximum of max_spaces spaces
    :param text:
    :param max_length:
    :param max_spaces:
    :return:
    """
    spaces = len([c for c in text if c.isspace()])
    needed = max_length - (len(text) - spaces)
    space_lst = [needed] + spaces * [0]

    if len(text) < 10 or len(text) > max_length:
        return text

    if spaces < 3 or needed < 0:
        return text

    while True:
        top = space_lst[0]
        if top > spaces:
            space_lst[0] -= spaces
            space_lst[1:] = [space_lst[i + 1] + 1 for i in range(spaces)]
            continue
        space_lst[0] -= top
        space_lst[1:top + 1] = [space_lst[i + 1] + 1 for i in range(top)]
        space_lst = space_lst[1:]
        space_lst = sorted(space_lst, reverse=False)
        break
    position = -1
    if max(space_lst) > max_spaces:
        return text

    for i in range(0, spaces):
        position = text.find(' ', position + 1)
        text = text[:position] + ' ' * space_lst[i] + text[position + 1:]
        position += space_lst[i]
    return text
