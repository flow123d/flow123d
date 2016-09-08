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
from utils.argparser import ArgParser
from utils.duration import Duration
# ----------------------------------------------

parser = ArgParser("exec_with_limit.py [-t <time>] [-m <memory>] -- <executable> <arguments>")
# ----------------------------------------------
parser.add('-t', '--limit-time', type=Duration.parse, name='time_limit', placeholder='<time>', docs=[
    'Obligatory wall clock time limit for execution in seconds or in HH:MM:SS format.',
    'For precision use float value'
])
parser.add('-m', '--limit-memory', type=float, name='memory_limit', placeholder='<memory>', docs=[
    'Optional memory limit per node in MB',
    'For precision use float value'
])
parser.add('', '--batch', type=True, name='batch', docs=[
    'Make output of this script more for an off-line reading',
    'In batch mode, stdout and stderr from executed processes will be printed, not saved'
])
# ----------------------------------------------

if __name__ == '__main__':
    from scripts.exec_with_limit_module import do_work

    # run work
    returncode = do_work(parser)
    sys.exit(returncode)
