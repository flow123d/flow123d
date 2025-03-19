#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from py123d.bin import pathfix
import sys
pathfix.init()
# ----------------------------------------------
from py123d.utils.globals import check_modules
if not check_modules('psutil', 'yaml', 'shutil', 'importlib', 'platform', 'loguru'):
    sys.exit(1)
# ----------------------------------------------
from py123d.loggers import printf
import py123d.utils.parsers as parsers
import py123d.utils.argparser as argparser
from py123d.utils.timer import Timer
from py123d.scripts.core.execution import BinExecutor
from py123d.scripts.runtest_module import do_work
# ----------------------------------------------


def create_parser():
    parser = argparser.Parser.create("runtest")

    group = parser.add_argument_group('General arguments')
    argparser.Parser.add(group, '--keep-going, -a', action=argparser.Parser.STORE_TRUE, help="""R|
        Run all tests, do not stop on the first error.
        In PBS mode this arguments is ignored.
    """)
    argparser.Parser.add(group, '--valgrind', nargs='?', default=False, help="""R|
        Run tests under valgrind, with python suppression with optional argument
        <valgrind args> passed to the valgrind. (In PBS mode this arguments is ignored.)?
    """)
    argparser.Parser.add(group, '--massif', action='store_true', default=False, help="""R|
        If set will print valgrind massif memory allocation stats.
    """)
    argparser.Parser.add(group, '--parallel, -p', default=1, type=int, help="""R|
        Run at most N jobs in parallel.
        In PBS mode this arguments is ignored.
    """)
    argparser.Parser.add(group, '--batch', action=argparser.Parser.STORE_TRUE, help="""R|
        Make output of this script more for an off-line reading
        In batch mode, stdout and stderr from executed processes will be printed, not saved
        [DEPRECATED] use --verbosity=full
    """)
    argparser.Parser.add(group, '--verbosity, -v',default=printf.OutputVerbosity.SMART, type=printf.OutputVerbosity, help="""R|
        A verbosity of the output, one of the {smart, minimal, full}
    """)
    argparser.Parser.add(group, '--include', nargs='*', metavar='TAG', help="""R|
        By default all tags are processed but if --include is set
        only those matching tags will be processed
        If combined with --exclude, only tags matching --include
        filtered with --exclude will be processed.
    """)
    argparser.Parser.add(group, '--exclude', nargs='*', metavar='TAG', help="""R|
        Filter tags which should be processed.
        See --include for more information. By default tag "disabled" is set if no
        other --exclude flag is set
    """)

    group = parser.add_argument_group('run_parallel arguments', 'Passable arguments to run_parallel.py')
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
        Path to base directory of flow123d. If not set value will be determined by python scripts
        location.
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
        Will also redirect output to a file.
    """)
    argparser.Parser.add(group, '--random-output-dir', default=False, nargs='?', help="""R|
        If set, output directories will have random suffix appended.
        If value is set, it will be used instead of random hash
    """)
    argparser.Parser.add(group, '--no-clean', action=argparser.Parser.STORE_TRUE, help="""R|
        If set, output directories will not be cleaned beforehand.
    """)
    argparser.Parser.add(group, '--no-compare', action=argparser.Parser.STORE_TRUE, help="""R|
        If set, results will not be compared.
    """)
    argparser.Parser.add(group, '--death-test', action=argparser.Parser.STORE_TRUE, help="""R|
        Perform death test, instead of regular one.
        Death test will succeed if and only if tests fails.
    """)
    argparser.Parser.add(group, '--save-to-db', action=argparser.Parser.STORE_TRUE, help="""R|
        If set, will process artifacts yaml file in order to save/copy file
        or export them to database.
    """)
    argparser.Parser.add(group, '--status-file', action=argparser.Parser.STORE_TRUE, help="""R|
        If set, will also generate status file in test_results directory.
        Status file will have name "runtest.status.json" and will contain
        additional information which cannot be obtained from profiler file
    """)
    return parser


if __name__ == '__main__':
    parser = create_parser()
    arg_options = argparser.Parser.parse_runtest(parser)
    printf.init_logger(verbosity=arg_options.verbosity, logfile=arg_options.log)

    # run work
    with Timer.app_timer:
        BinExecutor.register_sigint()
        returncode, debug = do_work(arg_options)

    # run work
    if arg_options.death:
        if returncode == 0:
            printf.error('Command did exit with 0 but should not (--death flag was set)!')
            sys.exit(1)
        else:
            printf.success('Command did not with 0 (--death flag was set)')
            sys.exit(0)
    else:
        sys.exit(returncode())
