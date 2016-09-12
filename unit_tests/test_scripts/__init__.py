#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import platform
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


def ensure_iterable(o):
    """
    Method ensure that given object is iterable(list or tuple)
    :param o: tested object
    :return: list or tuple
    """
    return [o] if type(o) not in (list, tuple, set) else o


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


class limit_test(object):
    """
    Class limit_test is conditional wrapper for unit_tests.
    Using this wrapper can turn off specific tests which are supposed to run on
    certain platform
    """

    verbose = False

    def __init__(self, hostname=None):
        self.limiters = list()
        if hostname:
            self.limiters.append(HostnameLimiter(hostname))

    def __call__(self, f):
        if self.limiters:
            if self.verbose:
                print('=' * 80)
                print('{:-^80}'.format(' checking test compatibility '))

            for limiter in self.limiters:
                if not limiter.test():
                    print('-' * 80)
                    print('{:-^80s}'.format(' skipping test "{f.func_name}" '.format(**locals())))
                    print('=' * 80)
                    return
            if self.verbose:
                print('=' * 80)

        def wrapper(other, *args, **kwargs):
            return f(other, *args, **kwargs)
        return wrapper


class HostnameLimiter(object):
    """
    Class HostnameLimiter is Limiter which will give green for test to run
    if current platform hostname matches one of the intended hostnames
    """

    node = platform.node()

    def __init__(self, hostname):
        self.hostname = ensure_iterable(hostname)

    def test(self):
        if not self.hostname or '*' in self.hostname:
            if limit_test.verbose:
                print("Test allowed: found wildcard * in hostname list")
            return True
        if self.node in self.hostname:
            if limit_test.verbose:
                print("Test allowed: found '{self.node}' in {self.hostname}".format(self=self))
            return True
        if limit_test.verbose:
            print("Test denied: '{self.node}' not found in {self.hostname}".format(self=self))
        return False
