#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from __future__ import absolute_import
import pathfix
pathfix.init()
# ----------------------------------------------
import sys
# ----------------------------------------------
from scripts.core.base import Paths
from utils.argparser import ArgParser
from utils.duration import Duration
from utils.timer import Timer
from scripts.artifacts.artifacts import ArtifactProcessor
# ----------------------------------------------

parser = ArgParser("runtest.py [<parametes>] [<test set>]  [-- <test arguments>]")
# ----------------------------------------------
parser.add_section('General arguments')
parser.add('-a', '--keep-going', type=True, name='keep_going', docs=[
    'Run all tests, do not stop on the first error.',
    'In PBS mode this arguments is ignored.',
])
parser.add('-v', '--valgrind', type=[True, str], name='valgrind', placeholder='[<VALGRIND ARGS>]', docs=[
    'Run tests under valgrind, with python suppression with optional argument',
    ' <valgrind args> passed to the valgrind. (In PBS mode this arguments is ignored.)?',
])
parser.add('-p', '--parallel', type=int, name='parallel', default=1, placeholder='<N>', docs=[
    'Run at most N jobs in parallel.',
    'In PBS mode this arguments is ignored.',
])
parser.add('', '--batch', type=True, name='batch', docs=[
    'Make output of this script more for an off-line reading',
    'In batch mode, stdout and stderr from executed processes will be printed, not saved'
])
parser.add('', '--include', type=list, subtype=str, name='include', docs=[
    'By default all tags are processed but if --include is set, ',
    'only those matching tags will be processed.',
    'If combined with --exclude, only tags matching --include,',
    'filtered with --exclude will be processed.'
])
parser.add('', '--exclude', type=list, subtype=str, name='exclude', docs=[
    'Filter tags which should be processed.',
    'See --include for more information. By default tag "disabled" is set if no',
    'other --exclude flag is set'
])
# ----------------------------------------------
parser.add_section('Passable arguments to run_parallel.py')
parser.add('-n', '--cpu', type=list, subtype=int, name='cpu', placeholder='<proc set>', docs=[
    'Run for every number of processes in the <proc set>',
    '  The <proc set> can be set as:',
    '     - single number (can be defined multiple times)',
    '     - set             "[1,3,4]" or [1 2 4]',
    '     - range           "1:4"   = "[1,2,3,4]"',
    '     - range with step "1:7:2" = "[1,3,5,7]"',
])
parser.add('-q', '--queue', type=[True, str], name='queue', placeholder='[<queue>]', docs=[
    'Optional PBS queue name to use. If the parameter is not used, ',
    'the application is executed in the same process and without PBS.',
    '',
    'If used without <queue> argument it is executed in the ',
    'background preferably under PBS with the queue selected ',
    'automatically for the given wall clock time limit and number of processes.'
])
parser.add('', '--host', type=str, name='host', placeholder='<host>', docs=[
    'Name of the running host that is used to select system ',
    'specific setup script. Default value of this parameter ',
    'is obtained by first getting the hostname ',
    '(using platform.node() or socket.gethostname()) and then search',
    'it in the table "host_table.json" which assign logical hostname',
    'possibly to multiple different real hostnames. ',
    '',
    'If the real host name is not found in the table ',
    'it is used directly otherwise the logical ',
    'hostname is used to select the setup script.',
])
# ----------------------------------------------
parser.add_section('Passable arguments to exec_with_limit.py')
parser.add('-t', '--limit-time', type=Duration.parse, name='time_limit', placeholder='<time>', docs=[
    'Obligatory wall clock time limit for execution in seconds or in HH:MM:SS format.',
    'For precision use float value'
])
parser.add('-m', '--limit-memory', type=float, name='memory_limit', placeholder='<memory>', docs=[
    'Optional memory limit per node in MB',
    'For precision use float value'
])
parser.add_section('Special options')
parser.add('', '--root', hidden=True, type=str, name='root', placeholder='<ROOT>', docs=[
    'Path to base dir of flow123d'
])
parser.add('', '--json', hidden=True, type=str, name='json', placeholder='<JSON>', docs=[
    'Output result to json file'
])
parser.add('', '--list', type=True, name='list', docs=[
    'List tests structure'
])
parser.add('', '--dump', hidden=True, type=str, name='dump', placeholder='<FILE>', docs=[
    'If set will pickle result to given file'
])
parser.add('', '--log', type=str, name='log', placeholder='<FILE>', docs=[
    'Will also redirect output to file'
])
parser.add('-x', '--export', type=True, name='artifacts', default=False, docs=[
    'If set, will process artifacts yaml file in order to save/copy file',
    'or export them to database.'
])

# ----------------------------------------------

if __name__ == '__main__':
    from utils.globals import check_modules

    required = ('psutil', 'yaml', 'shutil', 'importlib', 'platform')
    if not check_modules(*required):
        sys.exit(1)

    from scripts.core.execution import BinExecutor
    from scripts.runtest_module import do_work

    # determine batched mode after parsing
    from scripts.core.base import Printer
    parser.on_parse += Printer.setup_printer

    # import os
    # Paths.init(os.getcwd())

    # run work
    with Timer.app_timer:
        BinExecutor.register_sigint()
        returncode = do_work(parser)

    # collect artifact if not set otherwise
    # parser.parse()
    if parser.simple_options.artifacts:
        if parser.simple_options.artifacts is True:
            artifact_yml = Paths.artifact_yaml()
        else:
            artifact_yml = parser.simple_options.artifacts

        try:
            ap = ArtifactProcessor(artifact_yml)
            ap.run()
        except Exception as e:
            # we catch all error coming from artifact system
            # so it does not affect regular tests
            pass

    if type(returncode) is int:
        sys.exit(returncode)
    else:
        sys.exit(returncode.returncode)
