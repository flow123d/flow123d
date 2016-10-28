#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
# ----------------------------------------------
import test_scripts
test_scripts.fix_paths()
# ----------------------------------------------
from exec_with_limit import parser
from scripts.exec_with_limit_module import do_work
# ----------------------------------------------
consumer = test_scripts.get_consumer()


class TestDoWork(test_scripts.UnitTest):
    """
    Class TestDoWork tests interface for script exec_with_limit.py and backend
    exec_with_limit_module
    """

    @classmethod
    def setUpClass(cls):
        super(TestDoWork, cls).setUpClass()
        test_scripts.prepare_consumer()

    def test_help(self):
        try:
            do_work(parser, ['--help'])
            self.fail()
        except SystemExit as e:
            self.assertEqual(e.code, 0)

    def test_empty(self):
        try:
            do_work(parser, [])
            self.fail()
            pass
        except SystemExit as e:
            self.assertEqual(e.code, 1)

    def test_no_limits(self):
        try:
            do_work(parser, ['--', 'sleep', '0.1'])
            self.fail()
        except SystemExit as e:
            self.assertEqual(e.code, 2)

    def test_time_limit_ok(self):
        pypy = do_work(parser, ['-t', '1', '--', 'sleep', '0.01'])
        self.assertEqual(pypy.returncode, 0)

        pypy = do_work(parser, ['--limit-time', '1', '--', 'sleep', '0.01'])
        self.assertEqual(pypy.returncode, 0)

    def test_time_limit_over(self):
        pypy = do_work(parser, ['-t', '0.1', '--', consumer, '-t', '2'])
        self.assertNotEqual(pypy.returncode, 0)

    def test_memory_limit_ok(self):
        pypy = do_work(parser, ['-m', '100', '--', consumer, '-m', '10', '-t', '1'])
        self.assertEqual(pypy.returncode, 0)

        pypy = do_work(parser, ['--limit-memory', '100', '--', consumer, '-m', '10', '-t', '1'])
        self.assertEqual(pypy.returncode, 0)

    def test_memory_limit_over(self):
        pypy = do_work(parser, ['-m', '100', '--', consumer, '-m', '200', '-t', '1'])
        self.assertNotEqual(pypy.returncode, 0)

    def test_limits_ok(self):
        pypy = do_work(parser, ['-t', '2', '-m', '100', '--', consumer, '-m', '10', '-t', '1'])
        self.assertEqual(pypy.returncode, 0)

    def test_limits_over(self):
        pypy = do_work(parser, ['-t', '.5', '-m', '100', '--', consumer, '-m', '200', '-t', '1'])
        self.assertNotEqual(pypy.returncode, 0)

        pypy = do_work(parser, ['-t', '.1', '-m', '100', '--', consumer, '-m', '200', '-t', '1'])
        self.assertNotEqual(pypy.returncode, 0)
