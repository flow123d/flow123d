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
import utils.parsers as parsers
import utils.argparser as argparser
# ----------------------------------------------


def create_parser():
    parser = argparser.Parser.create("exec_with_limit")
    argparser.Parser.add(parser, '--limit-time, -t', dest='time_limit', type=parsers.parse_float, help="""R|
        Obligatory wall clock time limit for execution in seconds or in HH:MM:SS format
        For precision use float value
    """)
    argparser.Parser.add(parser, '--limit-memory, -m', dest='memory_limit', type=float, help="""R|
        Optional memory limit per node in MB
        For precision use float value
    """)
    argparser.Parser.add(parser, '--batch', action=argparser.Parser.STORE_TRUE, help="""R|
        Optional memory limit per node in MB
        For precision use float value
    """)

    group = parser.add_argument_group('Special options', 'Options are debug features or machine specific options')
    argparser.Parser.add(group, '--death', action=argparser.Parser.STORE_TRUE, help="""R|
        Reverse evaluation behaviour. Program will exit with 0 (success)
        if and only if command fails (ends with non-zero).
    """)
    return parser

if __name__ == '__main__':
    from scripts.exec_with_limit_module import do_work

    # determine batched mode after parsing
    from scripts.core.base import Printer
    argparser.Parser.on_parse += Printer.setup_printer
    parser = create_parser()

    arg_options = argparser.Parser.parse_exec_with_limit(parser)

    # run work
    returncode = do_work(arg_options)
    returncode = returncode.returncode if type(returncode) is not int else returncode

    if arg_options.death:
        if returncode == 0:
            Printer.all.err('Command did exit with 0 but should not (--death flag was set)!')
            sys.exit(1)
        else:
            Printer.all.suc('Command did not with 0 (--death flag was set)')
            sys.exit(0)
    else:
        sys.exit(returncode)
