#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import subprocess
import datetime
import time
# ----------------------------------------------
from scripts.core.base import Paths, Printer, Command, IO
from scripts.core.base import PathFormat
from scripts.core.prescriptions import PBSModule
from scripts.core.threads import BinExecutor, PyPy
from scripts.pbs.common import get_pbs_module, job_ok_string
from scripts.pbs.job import JobState, finish_pbs_job, MultiJob
from utils.dotdict import Map
# ----------------------------------------------

# global arguments
arg_options = None
arg_others = None
arg_rest = None
debug_mode = False


def run_local_mode():
    """
    :rtype: scripts.core.threads.PyPy or int
    """
    global arg_options, arg_others, arg_rest, debug_mode
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
    pypy = PyPy(executor, progress=not arg_options.batch)

    # set limits
    pypy.limit_monitor.time_limit = arg_options.time_limit
    pypy.limit_monitor.memory_limit = arg_options.memory_limit

    # turn on output
    if arg_options.batch:
        pypy.info_monitor.stdout_stderr = None
    else:
        pypy.info_monitor.stdout_stderr = Paths.temp_file('exec-paral.log')

    # start process
    pypy.start()
    pypy.join()

    if debug_mode:
        return pypy
    return pypy.returncode


def run_pbs_mode():
    """
    :rtype: scripts.core.threads.PyPy or int
    """
    global arg_options, arg_others, arg_rest, debug_mode
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
    case = pbs_module.Module(test_case, arg_options.cpu, None)
    pbs_command = case.get_pbs_command(arg_options, case.pbs_script)

    # create regular command for execution
    escaped_command = ' '.join(Command.escape_command(command))

    # create pbs script
    pbs_content = PBSModule.format(
        pbs_module.template,
        command=escaped_command,
        root=case.output_dir,
        output=case.job_output,
        status_ok=job_ok_string,
    )

    # save pbs script
    IO.write(case.pbs_script, pbs_content)

    # run qsub command
    output = subprocess.check_output(pbs_command)
    start_time = time.time()
    job = pbs_module.ModuleJob.create(output, case)
    job.full_name = "MPI exec job"

    multijob = MultiJob(pbs_module.ModuleJob)
    multijob.add(job)

    Printer.dyn('Updating job status...')
    multijob.update()
    Printer.out('Job submitted: {}', job)

    # wait for job to end
    while job.status != JobState.COMPLETED:
        for j in range(6):
            elapsed_str = str(datetime.timedelta(seconds=int(time.time() - start_time)))
            Printer.dyn('Job #{job.id} status: {job.state} ({t})', job=job, t=elapsed_str)

            # test job status
            if job.status == JobState.COMPLETED:
                break

            # sleep for a bit
            time.sleep(0.5)

        # update status every 6 * 0.5 seconds (3 sec update)
        multijob.update()

    returncode = finish_pbs_job(job, arg_options.batch)
    if debug_mode:
        return job
    return returncode


def do_work(parser, args=None, debug=False):
    """
    :type args: list
    :type parser: utils.argparser.ArgParser
    """

    # parse arguments
    global arg_options, arg_others, arg_rest, debug_mode
    arg_options, arg_others, arg_rest = parser.parse(args)
    debug_mode = debug

    # configure path
    Paths.format = PathFormat.ABSOLUTE
    if arg_options.root:
        Paths.base_dir(arg_options.root)

    # check commands
    if not arg_rest:
        parser.exit_usage('No command specified!', exit_code=1)

    # we need flow123d, mpiexec and ndiff to exists in LOCAL mode
    if not arg_options.queue and not Paths.test_paths('mpiexec'):
        Printer.err('Missing obligatory files! Exiting')
        exit(1)

    # run local or pbs mode
    if arg_options.queue:
        Printer.out('Running in PBS mode')
        Printer.separator()
        return run_pbs_mode()
    else:
        Printer.out('Running in LOCAL mode')
        Printer.separator()
        return run_local_mode()
