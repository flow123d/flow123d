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
from scripts.core.base import IO
from scripts.exec_parallel_module import do_work
# ----------------------------------------------
consumer = test_scripts.get_consumer()


class TestDoWork(test_scripts.UnitTest):

    @classmethod
    def setUpClass(cls):
        super(TestDoWork, cls).setUpClass()
        test_scripts.prepare_consumer()

    @test_scripts.print_test
    def test_help(self):
        try:
            do_work(parser, ['--help'])
            self.fail()
        except SystemExit as e:
            self.assertEqual(e.code, 0)

    @test_scripts.print_test
    def test_empty(self):
        try:
            do_work(parser, [])
            self.fail()
            pass
        except SystemExit as e:
            self.assertEqual(e.code, 1)

    @test_scripts.print_test
    def test_no_limits(self):
        returncode = do_work(parser, ['--', 'mpirun', 'sleep', '0.1'])
        self.assertEqual(returncode, 0)

    @test_scripts.print_test
    def test_time_limit_ok(self):
        returncode = do_work(parser, ['-t', '1', '--', 'mpirun', 'sleep', '0.1'])
        self.assertEqual(returncode, 0)

        np = '2'
        pypy = do_work(parser, ['-n', np, '-t', '1', '--', 'mpirun', 'sleep', '0.1'], debug=True)
        self.assertEqual(pypy.returncode, 0)
        self.assertEqual(int(pypy.executor.command[2]), int(np))

    @test_scripts.print_test
    def test_time_limit_over(self):
        returncode = do_work(parser, ['-t', '0.1', '--', 'mpirun', 'sleep', '2'])
        self.assertNotEqual(returncode, 0)

        np = '2'
        pypy = do_work(parser, ['-n', np, '-t', '0.1', '--', 'mpirun', 'sleep', '2'], debug=True)
        self.assertNotEqual(pypy.returncode, 0)
        self.assertEqual(int(pypy.executor.command[2]), int(np))

    @test_scripts.print_test
    def test_memory_limit_ok(self):
        returncode = do_work(parser, ['-n', '1', '-m', '100', '--', 'mpirun', consumer, '-m', '10', '-t', '1'])
        self.assertEqual(returncode, 0)

        returncode = do_work(parser, ['-n', '2', '-m', '200', '--', 'mpirun', consumer, '-m', '10', '-t', '1'])
        self.assertEqual(returncode, 0)

        np = '2'
        pypy = do_work(parser, ['-n', np, '-m', '200', '--', 'mpirun', consumer, '-m', '10', '-t', '1'], debug=True)
        self.assertEqual(int(pypy.returncode), 0)
        self.assertEqual(int(pypy.executor.command[2]), int(np))

    @test_scripts.print_test
    def test_memory_limit_over(self):
        returncode = do_work(parser, ['-n', '1', '-m', '100', '--', 'mpirun', consumer, '-m', '200', '-t', '5'])
        self.assertNotEqual(returncode, 0)

        returncode = do_work(parser, ['-n', '2', '-m', '200', '--', 'mpirun', consumer, '-m', '200', '-t', '1'])
        self.assertNotEqual(returncode, 0)

        np = '2'
        pypy = do_work(parser, ['-n', np, '-m', '200', '--', 'mpirun', consumer, '-m', '200', '-t', '1'], debug=True)
        self.assertNotEqual(int(pypy.returncode), 0)
        self.assertEqual(int(pypy.executor.command[2]), int(np))

    @test_scripts.print_test
    def test_pbs_mode(self):
        # this test can be only tested on PBS server
        # on all server we can only test whether pbs script are prepared,
        # assuming local_pbs (debug-like module) is assigned to current platform in host_table.json

        # current_dir = test_scripts.current_dir()
        # do_work(parser, ['-q', '--', 'sleep', '1'], debug=True)
        # mpi_dirs = [f for f in os.listdir(current_dir) if f.startswith('exec-parallel-')]
        # mpi_dirs = [f for f in mpi_dirs if os.path.isdir(os.path.join(current_dir, f))]
        #
        # self.assertGreaterEqual(len(mpi_dirs), 1)
        #
        # # test that script were actually created and contain at least one file
        # # also remove them ...
        # for mpi in mpi_dirs:
        #     self.assertGreaterEqual(len(os.listdir(os.path.join(current_dir, mpi))), 1)
        #     IO.delete_all(os.path.join(current_dir, mpi))

        # no tests for now
        pass
