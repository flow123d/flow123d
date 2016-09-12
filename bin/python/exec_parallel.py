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
from scripts.core.base import Paths, GlobalResult
from utils.argparser import ArgParser
from utils.duration import Duration
# ----------------------------------------------


parser = ArgParser("exec_parallel.py <parameters>  -- <executable> <executable arguments>")
# ----------------------------------------------
parser.add_section('General arguments')
parser.add('-n', '--cpu', type=list, name='cpu', default=[1], placeholder='<cpu>', docs=[
    'Run executable in <cpu> processes',
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
parser.add('-v', '--valgrind', type=[True, str], name='valgrind', placeholder='[<VALGRIND ARGS>]', docs=[
    'Run command under valgrind, with python suppression with optional argument',
    ' <valgrind args> passed to the valgrind. (In PBS mode this arguments is ignored.)?',
])
parser.add('', '--batch', type=True, name='batch', docs=[
    'Make output of this script more suitable for an off-line logging. Features such as ',
    'progress bars and other output will be disabled.'
    'In batch mode, stdout and stderr from executed processes will be printed, not saved into file.'
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
parser.add('', '--json', hidden=True, type=str, name='json', placeholder='<JSON>', docs=[
    'Output result to json file'
])
parser.add('', '--dump', hidden=True, type=str, name='dump', placeholder='<FILE>', docs=[
    'If set will pickle result to given file'
])
# ----------------------------------------------

if __name__ == '__main__':
    from utils.globals import check_modules

    required = ('psutil', 'importlib', 'platform')
    if not check_modules(*required):
        sys.exit(1)

    from scripts.exec_parallel_module import do_work

    # run work
    returncode = do_work(parser)
    if type(returncode) is int:
        sys.exit(returncode)
