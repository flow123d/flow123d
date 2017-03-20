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
import exec_with_limit as script
import utils.argparser as argparser
from scripts.exec_with_limit_module import do_work
# ----------------------------------------------
consumer = test_scripts.get_consumer()


def parse(command, debug=False):
    arg_options = argparser.Parser.parse_exec_with_limit(
        script.create_parser(),
        command.split()
    )
    return do_work(arg_options, debug)


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
            command = '--help'
            pypy = parse(command)
            self.fail()
        except SystemExit as e:
            self.assertEqual(e.code, 0)

    def test_empty(self):
        try:
            command = ''
            pypy = parse(command)
            self.fail()
            pass
        except SystemExit as e:
            self.assertEqual(e.code, 2)

    def test_no_limits(self):
        try:
            command = '-- sleep 0.1'
            pypy = parse(command)
            self.fail()
        except SystemExit as e:
            self.assertEqual(e.code, 2)

    def test_time_limit_ok(self):
        command = '-t 1 -- sleep 0.01'
        pypy = parse(command)
        self.assertEqual(pypy.returncode, 0)

        command = '--limit-time 1 -- sleep 0.01'
        pypy = parse(command)
        self.assertEqual(pypy.returncode, 0)

    def test_time_limit_over(self):
        command = '-t 0.1 -- ' + consumer + ' -t 2'
        pypy = parse(command)
        self.assertNotEqual(pypy.returncode, 0)

    def test_memory_limit_ok(self):
        command = '-m 100 -- ' + consumer + ' -m 10 -t 1'
        pypy = parse(command)
        self.assertEqual(pypy.returncode, 0)

        command = '--limit-memory 100 -- ' + consumer + ' -m 10 -t 1'
        pypy = parse(command)
        self.assertEqual(pypy.returncode, 0)

    def test_memory_limit_over(self):
        command = '-m 100 -- ' + consumer + ' -m 200 -t 1'
        pypy = parse(command)
        self.assertNotEqual(pypy.returncode, 0)

    def test_limits_ok(self):
        command = '-t 2 -m 100 -- ' + consumer + ' -m 10 -t 1'
        pypy = parse(command)
        self.assertEqual(pypy.returncode, 0)

    def test_limits_over(self):
        command = '-t 0.5 -m 100 -- ' + consumer + ' -m 200 -t 1'
        pypy = parse(command)
        self.assertNotEqual(pypy.returncode, 0)

        command = '-t 0.1 -m 100 -- ' + consumer + ' -m 200 -t 1'
        pypy = parse(command)
        self.assertNotEqual(pypy.returncode, 0)
