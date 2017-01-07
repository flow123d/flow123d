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
import runtest as script
import utils.argparser as argparser
from scripts.runtest_module import do_work
# ----------------------------------------------
__dir__ = test_scripts.current_dir()
extras = os.path.join(__dir__, 'extras')


root = os.path.abspath(os.path.join(__dir__, '..', '..'))
EXIT_OK = 0
EXIT_ERROR = 1
EXIT_COMPARE_ERROR = 13


def parse(command, debug=False):
    """

    :type command: list or str
    """
    arg_options = argparser.Parser.parse_runtest(
        script.create_parser(),
        command.split() if type(command) is str else command
    )
    return do_work(arg_options, debug)


class TestDoWork(test_scripts.UnitTest):
    """
    Class TestDoWork tests interface for script runtest.py and backend runtest_module
    """

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

    def test_help(self):
        try:
            command = '--help'
            result = parse(command)
            self.fail()
        except SystemExit as e:
            self.assertEqual(e.code, 0)

    def test_empty(self):
        try:
            command = ''
            result = parse(command)
            self.fail()
            pass
        except SystemExit as e:
            self.assertEqual(e.code, 2)

    def test_compare_correct_output(self):
        # all test must pass since we are copying reference outputs
        try:
            command = [
                '--root', root,
                os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit.yaml'),
                '--',
                '-t', '0.1', '--copy', '--clean'
            ]
        except SystemExit as se:
            self.assertEqual(se.code, EXIT_OK)

    def test_compare_no_output(self):
        # all test must fail since test_outputs are empty
        try:
            command = [
                '--root', root,
                os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit.yaml'),
                '--',
                '-t', '0.1', '--clean'
            ]
        except SystemExit as se:
            self.assertEqual(se.code, EXIT_COMPARE_ERROR)

    def test_compare_wrong_output_error(self):
        # all test must fail since test_outputs are empty
        try:
            command = [
                '--root', root,
                os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit.yaml'),
                '--',
                '-t', '0.1', '--clean', '--random'
            ]
        except SystemExit as se:
            self.assertEqual(se.code, EXIT_COMPARE_ERROR)

    def test_missing_yaml(self):
        # test fail if we pass non_existent yaml file
        try:
            command = [
                '--root', root,
                os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit_non_existent.yaml'),
                os.path.join(root, 'tests', '03_transport_small_12d', 'flow_implicit.yaml'),
                '--',
                '-t', '0.01', '--copy', '--clean'
            ]
        except SystemExit as se:
            self.assertNotEqual(se.code, EXIT_OK)

    def test_argument_pass(self):
        # we test that arguments are passed to flow123d by specifying sleep duration for flow123d_mock scripts
        # if passing works, duration for this test must be greater than sleep time passed
        start_time = time.time()
        sleep_time = 2
        try:
            command = [
                '--root', root,
                'extras',
                '--',
                '-t', str(sleep_time), '-e'
            ]
            result = parse(command)
        except SystemExit as se:
            self.assertNotEqual(se.code, EXIT_OK)

        end_time = time.time()
        diff = end_time - start_time
        self.assertGreater(diff, sleep_time)

    def test_config_load(self):
        config="""
common_config:
  proc: {common_proc}
  memory_limit: 1000

test_cases:
- files: test-02.yaml
  proc: {test02_proc}
""".strip()

        # ------------------------------------------------------------------------
        common_proc = [1, 2, 3]
        test02_proc = [1, 2, 3, 4]
        with open('extras/foo/bar/config.yaml', 'w') as fp:
            fp.write(config.format(**locals()))

        # command = 'extras/foo'
        # self.assertEqual(len(common_proc)+len(test02_proc), len(result.threads))
        #
        # command = '-p', '7', 'extras'
        # self.assertEqual(len(common_proc)+len(test02_proc), len(result.threads))

        command = 'extras/foo/bar/test-01.yaml'
        result = parse(command)
        self.assertEqual(len(result.threads), len(common_proc))
        self.assertEqual(result.returncode, 0)

        command = 'extras/foo/bar/test-02.yaml'
        result = parse(command)
        self.assertEqual(len(result.threads), len(test02_proc))
        self.assertEqual(result.returncode, 0)

        # ------------------------------------------------------------------------
        common_proc = [1]
        test02_proc = []
        with open('extras/foo/bar/config.yaml', 'w') as fp:
            fp.write(config.format(**locals()))

        command = 'extras/foo'
        result = parse(command)
        self.assertEqual(len(common_proc)+len(test02_proc), len(result.threads))
        self.assertEqual(result.returncode, 0)

        command = '-p 2 extras/foo/bar'
        result = parse(command)
        self.assertEqual(len(common_proc)+len(test02_proc), len(result.threads))
        self.assertEqual(result.returncode, 0)

        command = 'extras/foo/bar/test-01.yaml'
        result = parse(command)
        self.assertEqual(len(result.threads), len(common_proc))
        self.assertEqual(result.returncode, 0)

        command = 'extras/foo/bar/test-02.yaml'
        result = parse(command)
        self.assertEqual(len(result.threads), len(test02_proc))
        self.assertEqual(result.returncode, 0)