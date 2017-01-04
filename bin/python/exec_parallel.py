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
import utils.parsers as parsers
import utils.argparser as argparser
from utils.timer import Timer
from scripts.artifacts.artifacts import ArtifactProcessor
# ----------------------------------------------


def create_parser():
    parser = argparser.Parser.create("runtest")

    group = parser.add_argument_group('General arguments')
    argparser.Parser.add(group, '--valgrind, -v', nargs='?', default=False, help="""R|
        Run tests under valgrind, with python suppression with optional argument
        <valgrind args> passed to the valgrind. (In PBS mode this arguments is ignored.)?
    """)
    argparser.Parser.add(group, '--parallel, -p', default=1, type=int, help="""R|
        Run at most N jobs in parallel.
        In PBS mode this arguments is ignored.
    """)
    argparser.Parser.add(group, '--batch', action=argparser.Parser.STORE_TRUE, help="""R|
        Make output of this script more for an off-line reading
        In batch mode, stdout and stderr from executed processes will be printed, not saved
    """)
    argparser.Parser.add(group, '--cpu, -n', type=parsers.parse_int_list, action=argparser.Parser.APPEND, help="""R|
        Run for every number of processes in the <proc set>
          The <proc set> can be set as:
             - single number (can be defined multiple times)
             - set             "[1,3,4]" or "[1 2 4]"
             - range           "1:4"   = "[1,2,3,4]"
             - range with step "1:7:2" = "[1,3,5,7]"
    """)
    argparser.Parser.add(group, '--queue, -q', nargs='?', default=False, help="""R|
        Optional PBS queue name to use. If the parameter is not used,
        the application is executed in the same process and without PBS.

        If used without <queue> argument it is executed in the
        background preferably under PBS with the queue selected
        automatically for the given wall clock time limit and number of processes
    """)
    argparser.Parser.add(group, '--host', help="""R|
        Name of the running host that is used to select system
        specific setup script. Default value of this parameter
        is obtained by first getting the hostname
        (using platform.node() or socket.gethostname()) and then search
        it in the table "host_table.json" which assign logical hostname
        possibly to multiple different real hostnames.

        If the real host name is not found in the table
        it is used directly otherwise the logical
        hostname is used to select the setup script.
    """)

    group = parser.add_argument_group('exec_with_limit arguments', 'Passable arguments to exec_with_limit.py')
    argparser.Parser.add(group, '--limit-time, -t', dest='time_limit', type=parsers.parse_float, help="""R|
        Obligatory wall clock time limit for execution in seconds or in HH:MM:SS format.
        For precision use float value
    """)
    argparser.Parser.add(group, '--limit-memory, -m', dest='memory_limit', type=float, help="""R|
        Optional memory limit per node in MB
        For precision use float value
    """)

    group = parser.add_argument_group('Special options', 'Options are debug features or machine specific options')
    argparser.Parser.add(group, '--root', help="""R|
        Path to base directory of flow123d.
    """)
    argparser.Parser.add(group, '--json', help="""R|
        Output result to json file.
    """)
    argparser.Parser.add(group, '--list', action=argparser.Parser.STORE_TRUE, help="""R|
        List tests structure.
    """)
    argparser.Parser.add(group, '--dump', help="""R|
        If set will pickle result to given file.
    """)
    argparser.Parser.add(group, '--log', help="""R|
        Will also redirect output to file.
    """)
    argparser.Parser.add(group, '--death', action=argparser.Parser.STORE_TRUE, help="""R|
        Reverse evaluation behaviour. Program will exit with 0 (success)
        if and only if command fails (ends with non-zero).
    """)
    return parser


if __name__ == '__main__':
    from utils.globals import check_modules

    required = ('psutil', 'importlib', 'platform')
    if not check_modules(*required):
        sys.exit(1)

    from scripts.exec_parallel_module import do_work

    # determine batched mode after parsing
    from scripts.core.base import Printer
    argparser.Parser.on_parse += Printer.setup_printer
    parser = create_parser()
    arg_options = argparser.Parser.parse_exec_parallel(parser)
    #
    # import os
    # Paths.init(os.getcwd())

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
