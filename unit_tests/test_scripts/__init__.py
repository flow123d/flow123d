#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from unittest import TestCase
import os

# to run tests for individual modules use:
#   python -m unittest test_scripts.test_exec_with_limit_module
#
# to run all tests, simple run discover (python 2.7 and higher)
#   python -m unittest discover
# from unit_tests dir


# in console mode try to get size in console window
try:
    rows, cols = [int(x) for x in os.popen('stty size', 'r').read().split()]
except Exception as e:
    rows, cols = 80, 80


def print_test(f):
    import sys

    def wrapper(*args, **kwargs):
        # sys.stdout.write('\n')
        # sys.stdout.write('=' * cols + '\n')
        # msg = ' EXECUTING {:=<40s}'.format(str(f.func_name) + ' ')
        # sys.stdout.write(('{:=^'+str(cols)+'}').format(msg) + '\n')
        # sys.stdout.write('=' * cols + '\n')
        return f(*args, **kwargs)
    return wrapper


class UnitTest(TestCase):

    @classmethod
    def setUpClass(cls):
        super(UnitTest, cls).setUpClass()

        cls.__dir__ = current_dir()
        cls.start_state = os.listdir(cls.__dir__)

    @classmethod
    def tearDownClass(cls):
        super(UnitTest, cls).tearDownClass()

        from scripts.core.base import IO
        cls.end_state = os.listdir(cls.__dir__)

        for f in cls.end_state:
            if f not in cls.start_state:
                full_path = os.path.join(cls.__dir__, f)
                if os.path.isdir(full_path):
                    IO.delete_all(full_path)
                else:
                    IO.delete(full_path)


def fix_paths():
    """
    Add path for src/python and bin/python modules
    """
    import sys
    import os

    # we assume that this file is located in unit_tests/test_scripts
    __dir__ = test_dir()

    # add path to src/python folder
    sys.path.append(
        os.path.join(
            __dir__, '..', '..', 'src', 'python'
        )
    )

    # also add path to scripts in bin/python
    sys.path.append(
        os.path.join(
            __dir__, '..', '..', 'bin', 'python'
        )
    )


def test_dir():
    import os
    # we assume that this file is located in unit_tests/test_scripts
    return os.path.dirname(os.path.realpath(__file__))


def current_dir():
    import os
    return os.getcwd()


def get_consumer():
    __dir__ = current_dir()
    return os.path.join(__dir__, 'consumer_cc')


def prepare_consumer():
    import os
    __dir__ = current_dir()
    extras = os.path.join(__dir__, 'extras')
    consumer_bin = os.path.join(__dir__, 'consumer_cc')
    consumer_src = os.path.join(extras, 'consumer.cc')

    # try to compile consumer.cc if not exists
    if not os.path.exists(consumer_bin):
        import subprocess as sub

        process = sub.Popen([
            'g++', '-std=c++11', consumer_src,
            '-o', consumer_bin
        ])
        rc = process.wait()
        if rc != 0:
            print ('#' * 70)
            print ('# Could not compile consumer.cc file, cannot fully run unittests')
            print ('#' * 70)
        return rc
    return 0
