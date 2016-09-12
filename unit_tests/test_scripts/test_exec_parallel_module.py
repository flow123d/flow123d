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
from exec_parallel import parser
from scripts.exec_parallel_module import do_work
# ----------------------------------------------
consumer = test_scripts.get_consumer()


class TestDoWork(test_scripts.UnitTest):
    """
    Class TestDoWork tests interface for script exec_parallel.py and backend exec_parallel_module
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
        pypy = do_work(parser, ['--', 'mpirun', 'sleep', '0.1'])
        self.assertEqual(pypy.returncode, 0)

    def test_time_limit_ok(self):
        pypy = do_work(parser, ['-t', '1', '--', 'mpirun', 'sleep', '0.1'])
        self.assertEqual(pypy.returncode, 0)

        np = '2'
        pypy = do_work(parser, ['-n', np, '-t', '1', '--', 'mpirun', 'sleep', '0.1'])
        self.assertEqual(pypy.returncode, 0)
        self.assertEqual(int(pypy.executor.command[2]), int(np))

    def test_time_limit_over(self):
        pypy = do_work(parser, ['-t', '0.1', '--', 'mpirun', 'sleep', '2'])
        self.assertNotEqual(pypy.returncode, 0)

        np = '2'
        pypy = do_work(parser, ['-n', np, '-t', '0.1', '--', 'mpirun', 'sleep', '2'])
        self.assertNotEqual(pypy.returncode, 0)
        self.assertEqual(int(pypy.executor.command[2]), int(np))

    def test_memory_limit_ok(self):
        pypy = do_work(parser, ['-n', '1', '-m', '100', '--', 'mpirun', consumer, '-m', '10', '-t', '1'])
        self.assertEqual(pypy.returncode, 0)

        pypy = do_work(parser, ['-n', '2', '-m', '200', '--', 'mpirun', consumer, '-m', '10', '-t', '1'])
        self.assertEqual(pypy.returncode, 0)

        np = '2'
        pypy = do_work(parser, ['-n', np, '-m', '200', '--', 'mpirun', consumer, '-m', '10', '-t', '1'])
        self.assertEqual(int(pypy.returncode), 0)
        self.assertEqual(int(pypy.executor.command[2]), int(np))

    def test_memory_limit_over(self):
        pypy = do_work(parser, ['-n', '1', '-m', '100', '--', 'mpirun', consumer, '-m', '200', '-t', '5'])
        self.assertNotEqual(pypy.returncode, 0)

        pypy = do_work(parser, ['-n', '2', '-m', '200', '--', 'mpirun', consumer, '-m', '200', '-t', '1'])
        self.assertNotEqual(pypy.returncode, 0)

        np = '2'
        pypy = do_work(parser, ['-n', np, '-m', '200', '--', 'mpirun', consumer, '-m', '200', '-t', '1'])
        self.assertNotEqual(int(pypy.returncode), 0)
        self.assertEqual(int(pypy.executor.command[2]), int(np))