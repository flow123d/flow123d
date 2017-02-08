#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import os
import sys
import shutil
# ----------------------------------------------
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
# ----------------------------------------------
import test_scripts
test_scripts.fix_paths()
# ----------------------------------------------
import utils.argparser as argparser
from scripts.core.threads import ResultHolder
from scripts.pbs.modules.local_pbs import Module
# ----------------------------------------------


def runtest_call(*args, **kwargs):
    import runtest
    from scripts.runtest_module import do_work as runtest_do_work

    parser = runtest.create_parser()
    arg_options = argparser.Parser.parse_runtest(parser, list(args))
    return runtest_do_work(arg_options, kwargs.get('debug', False))


def exec_call(*args, **kwargs):
    import exec_parallel
    from scripts.exec_parallel_module import do_work as exec_do_work

    parser = exec_parallel.create_parser()
    arg_options = argparser.Parser.parse_exec_parallel(parser, list(args))
    return exec_do_work(arg_options, kwargs.get('debug', False))


# grab current dir and set path to local mock
__dir__ = test_scripts.current_dir()
extras = os.path.join(__dir__, 'extras')
Module.mock = os.path.join(extras, 'pbs_mock.py')

# set exit codes
EXIT_OK = 0
EXIT_ERROR = 1
EXIT_COMPARE_ERROR = 13


class TestDoWork(test_scripts.UnitTest):
    """
    Class TestDoWork tests interface for script runtest.py and backend runtest_module
    """

    @test_scripts.limit_test(hostname='*')
    def test_single_job(self):
        """
        Run exec parallel with CPU = 0 meaning no MPI should be used
        """
        pypy = exec_call(*'-n 0 -q -- mpirun uname -a'.split())
        self.assertEqual(pypy.returncode, EXIT_OK)
        self.assertEqual(pypy.command[0], 'uname')

    @test_scripts.limit_test(hostname='*')
    def test_single_job(self):
        """
        Run exec parallel with CPU = 1
        """
        pypy = exec_call(*'-n1 -q -- mpirun uname -a'.split())
        self.assertEqual(pypy.returncode, EXIT_OK)
        self.assertEqual(pypy.command[0], 'mpirun')
        self.assertEqual(pypy.command[2], '1')


    @test_scripts.limit_test(hostname='*')
    def test_parallel_multijob(self):
        """
        Run exec parallel when CPU specification is range from 2 to 3 CPU
        """
        result = exec_call(*'-n 2:3 -q -- mpirun uname -a'.split())
        self.assertIsInstance(result, ResultHolder)
        self.assertEqual(result.returncode, EXIT_OK)
        self.assertEqual(len(result.items), 2)

    @test_scripts.limit_test(hostname='*')
    def test_parallel_multijob(self):
        """
        Run exec parallel with 2 cpu, first cpu should delete file but
        seconds cpu should not be able to delete tmp file since it does
        not exists anymore.
        """
        dummy_file = os.path.abspath('foobar')
        with open(dummy_file, 'w') as fp:
            fp.write('delete me')

        result = exec_call('-n', '2', '-q', '--', 'mpirun', 'rm', dummy_file)
        self.assertNotEqual(result.returncode, EXIT_OK)