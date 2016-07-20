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

import sys, os


def print_debug ():
    """Prints debug information about python"""
    print ("Python:    " + str (sys.version).replace ("\n", "") + ", " + str (sys.executable))


def append_to_path ():
    """Performs path fix"""

    # for now, always print debug info
    print_debug ()

    # add src/python into module path
    path = os.path.join ('..', '..', 'src', 'python')
    sys.path.append (path)
    path = os.path.join ('src', 'python')
    sys.path.append (path)
    sys.path.append (os.getcwd())
    # print 'paths: \n{:s}'.format ('\n'.join(sys.path))
