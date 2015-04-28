__author__ = 'Jan Hybs'

import sys, os


def print_debug ():
    # print debug information about python
    print ("Python:    " + str (sys.version).replace("\n", "") + ", " + str (sys.executable))


def append_to_path ():
    # for now, always print debug info
    print_debug()

    # add src/python into module path
    path = os.path.join('..', '..', 'src', 'python')
    sys.path.append (path)
