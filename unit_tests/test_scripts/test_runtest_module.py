#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import os
import shutil
import sys
import time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
# ----------------------------------------------
import test_scripts
test_scripts.fix_paths()
# ----------------------------------------------
from runtest import parser
from scripts.runtest_module import do_work
# ----------------------------------------------
__dir__ = test_scripts.current_dir()
extras = os.path.join(__dir__, 'extras')

root = r'/home/jan-hybs/Dokumenty/projects/Flow123d/flow123d'
EXIT_OK = 0
EXIT_ERROR = 1
EXIT_COMPARE_ERROR = 13


class TestDoWork(test_scripts.UnitTest):
    flow_path = os.path.join(root, 'bin', 'flow123d')
    flow_path_copy = os.path.join(root, 'bin', 'flow123d_copy')
    backup = os.path.exists(flow_path) or os.path.islink(flow_path)

    @classmethod
    def setUpClass(cls):
        super(TestDoWork, cls).setUpClass()
        if cls.backup:
            os.rename(cls.flow_path, cls.flow_path_copy)

        # copy bash
        shutil.copy(
            os.path.join(extras, 'flow123d.sh'),
            cls.flow_path
        )
        # copy python
        shutil.copy(
            os.path.join(extras, 'flow123d_mock.py'),
            cls.flow_path + '.py'
        )

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.flow_path + '.py')
        os.unlink(cls.flow_path)

        # restore prev file
        if cls.backup:
            os.rename(cls.flow_path_copy, cls.flow_path)
        pass

    @test_scripts.print_test
    def test_help(self):
        try:
            do_work(parser, ['--help'])
            self.fail()
        except SystemExit as e:
            self.assertEqual(e.code, EXIT_OK)

    @test_scripts.print_test
    def test_empty(self):
        try:
            do_work(parser, [])
            self.fail()
            pass
        except SystemExit as e:
            self.assertEqual(e.code, 1)

    @test_scripts.print_test
    def test_compare_correct_output(self):
        # all test must pass since we are copying reference outputs
        try:
            do_work(parser,[
                        '--root', root,
                        os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit.yaml'),
                        '--',
                        '-t', '0.1', '--copy', '--clean'
            ])
        except SystemExit as se:
            self.assertEqual(se.code, EXIT_OK)

    @test_scripts.print_test
    def test_compare_no_output(self):
        # all test must fail since test_outputs are empty
        try:
            do_work(parser,[
                        '--root', root,
                        os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit.yaml'),
                        '--',
                        '-t', '0.1', '--clean'
            ])
        except SystemExit as se:
            self.assertEqual(se.code, EXIT_COMPARE_ERROR)

    @test_scripts.print_test
    def test_compare_wrong_output_error(self):
        # all test must fail since test_outputs are empty
        try:
            do_work(parser,[
                        '--root', root,
                        os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit.yaml'),
                        '--',
                        '-t', '0.1', '--clean', '--random'
            ])
        except SystemExit as se:
            self.assertEqual(se.code, EXIT_COMPARE_ERROR)

    @test_scripts.print_test
    def test_missing_yaml(self):
        # test fail if we pass non_existent yaml file
        try:
            do_work(parser,[
                        '--root', root,
                        os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit_non_existent.yaml'),
                        os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit.yaml'),
                        '--',
                        '-t', '0.01', '--copy', '--clean'
            ])
        except SystemExit as se:
            self.assertNotEqual(se.code, EXIT_OK)

    @test_scripts.print_test
    def test_argument_pass(self):
        # we test that arguments are passed to flow123d by specifying sleep duration for flow123d_mock scripts
        # if passing works, duration for this test must be greater than sleep time passed
        start_time = time.time()
        sleep_time = 2
        try:
            do_work(parser,[
                        '--root', root,
                        os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit.yaml'),
                        '--',
                        '-t', str(sleep_time), '-e'
            ])
        except SystemExit as se:
            self.assertNotEqual(se.code, EXIT_OK)

        end_time = time.time()
        diff = end_time - start_time
        self.assertGreater(diff, sleep_time)
