#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

"""
This Module provides method for fixing module search path.
Situation is following:
 Scripts located in bin/python requires
 Modules located in src/python

So to sys.path is appended ../../src/python path
"""

import getpass
import sys
import os


def print_debug():
    """Prints debug information about python"""
    print ("Python {version}, {executable}".format(
        version=str(sys.version).replace("\n", ""),
        executable=sys.executable
    ))
    print ("CWD: {cwd}, USER: {whoami}".format(
        cwd=os.getcwd(),
        whoami=getpass.getuser())
    )
    print('-' * 80)


def add_path(*args):
    """
    Adds path to sys.path
    :param args:
    :return:
    """
    root = os.path.dirname(os.path.realpath(__file__))
    if not args:
        return root
    path = os.path.abspath(
        os.path.join(
            root,
            *args
        )
    )
    sys.path.append(path)
    return path


def append_to_path():
    """Performs path fix"""

    # for now, always print debug info
    print_debug()

    # path to lib
    add_path('lib')
    add_path('..', 'lib')
    # path to src/python if COPY_PYTHON is disabled
    add_path('..', '..', 'src', 'python')
    # path to lib/flow123d after COPY_PYTHON
    add_path('..', '..', 'build_tree', 'lib', 'flow123d')
    # path to lib/flow123d/site-packages additional libs after COPY_PYTHON
    add_path('..', '..', 'build_tree', 'lib', 'flow123d', 'site-packages')


# alias
init = append_to_path
