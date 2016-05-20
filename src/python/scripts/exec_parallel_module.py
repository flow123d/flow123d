#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from scripts.core.base import Paths, Printer, CommandEscapee, IO
from scripts.core.base import PathFormat
from scripts.core.prescriptions import PBSModule
from scripts.execs.monitor import PyProcess
from scripts.execs.test_executor import BinExecutor
from scripts.pbs.common import get_pbs_module
import subprocess, time, datetime

# global arguments
from scripts.pbs.job import JobState
from utils.dotdict import Map

arg_options = None
arg_others = None
arg_rest = None
printer = Printer(Printer.LEVEL_KEY)


def run_local_mode():
    # build command
    mpi_binary = 'mpirun' if arg_options.mpirun else Paths.mpiexec()
    command = [
        mpi_binary,
        '-np', str(arg_options.get('cpu', 1))
    ]
    if arg_options.valgrind:
        valgrind = ['valgrind']
        if type(arg_options.valgrind) is str:
            valgrind.extend(arg_options.valgrind.split())
        # append to command
        command = command + valgrind
    # append rest arguments
    command.extend(arg_rest)

    # prepare executor
    executor = BinExecutor(command)
    process_monitor = PyProcess(executor, printer, batch_mode=arg_options.batch)

    # set limits
    process_monitor.limit_monitor.time_limit = arg_options.time_limit
    process_monitor.limit_monitor.memory_limit = arg_options.memory_limit

    # turn on output
    if arg_options.batch:
        process_monitor.info_monitor.stdout_stderr = None
    else:
        process_monitor.info_monitor.stdout_stderr = Paths.temp_file('exec-paral.log')

    # start process
    process_monitor.start()


def run_pbs_mode():
    # build command
    mpi_binary = 'mpirun' if arg_options.mpirun else Paths.mpiexec()
    command = [
        mpi_binary,
        '-np', str(arg_options.get('cpu', 1))
    ]
    # append rest arguments
    command.extend(arg_rest)

    # get module
    pbs_module = get_pbs_module(arg_options.host)

    # create pbs command
    test_case = Map(
        memory_limit=arg_options.get('memory_limit', None) or 400,
        time_limit=arg_options.get('time_limit', None) or 30
    )
    module = pbs_module.Module(test_case, arg_options.cpu, None)
    temp_file = Paths.temp_file('exec-temp.qsub')
    pbs_command = module.get_pbs_command(arg_options, temp_file)

    # create regular command for execution
    escaped_command = ' '.join(CommandEscapee.escape_command(command))

    # create pbs script
    pbs_content = PBSModule.format(
        pbs_module.template,
        command=escaped_command,
        root=arg_options.root
    )

    # print debug info
    printer.dbg('Command : {}', escaped_command)
    printer.dbg('PBS     : {}', ' '.join(pbs_command))
    printer.dbg('script  : {}', temp_file)
    printer.dbg('')

    # save pbs script
    IO.write(temp_file, pbs_content)

    # run qsub command
    output = subprocess.check_output(pbs_command)
    start_time = time.time()
    job = pbs_module.ModuleJob.create(output)
    job.update_status()
    printer.key('Job submitted: {}', job)

    # wait for job to end
    while job.state != JobState.COMPLETED:
        for j in range(6):
            elapsed_str = str(datetime.timedelta(seconds=int(time.time() - start_time)))
            printer.out_rr(' ' * 60)
            printer.out_rr('Job #{job.id} status: {job.state} ({t})', job=job, t=elapsed_str)

            # test job state
            if job.state == JobState.COMPLETED:
                break

            # sleep for a bit
            time.sleep(0.5)

        # update status every 6 * 0.5 seconds (3 sec update)
        job.update_status()
    printer.key('\nJob ended')

    # delete tmp file
    IO.delete(temp_file)


def do_work(parser):
    """
    :type parser: utils.argparser.ArgParser
    """

    # parse arguments
    global arg_options, arg_others, arg_rest
    arg_options, arg_others, arg_rest = parser.parse()
    Paths.format = PathFormat.ABSOLUTE

    # check commands
    if not arg_rest:
        parser.exit_usage('No command specified!')

    # run local or pbs mode
    if arg_options.queue:
        printer.dbg('Running in PBS mode')
        printer.key('-' * 60)
        run_pbs_mode()
    else:
        printer.dbg('Running in LOCAL mode')
        printer.key('-' * 60)
        run_local_mode()
