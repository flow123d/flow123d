#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import os
import subprocess
import time
import sys
# ----------------------------------------------
from scripts.core import prescriptions
from scripts.core.base import Paths, PathFormat, PathFilters, Printer, Command, IO
from scripts.config.yaml_config import YamlConfig
from scripts.core.prescriptions import PBSModule
from scripts.core.threads import BinExecutor, ParallelThreads, SequentialThreads, PyPy
from scripts.pbs.common import get_pbs_module, job_ok_string
from scripts.pbs.job import JobState, MultiJob, finish_pbs_job
# ----------------------------------------------


# global arguments
arg_options = None
arg_others = None
arg_rest = None


def create_process(command, limits=None):
    """
    :type command: list[str]
    :type limits: scripts.config.yaml_config.YamlConfigCase
    """
    test_executor = BinExecutor(command)
    process_monitor = PyPy(test_executor)
    process_monitor.limit_monitor.set_limits(limits)
    return process_monitor


def format_case(case):
    return '{} x {}'.format(
        case.proc_value,
        Paths.path_end(case.test_case.files[0])
    )


def create_process_from_case(case):
    """
    :type case: scripts.core.prescriptions.TestPrescription
    """
    pypy = create_process(case.get_command(), case.test_case)
    pypy.case = case
    pypy.info_monitor.end_fmt = ''
    pypy.info_monitor.start_fmt = 'Running: {}'.format(format_case(case))

    # turn on output
    pypy.progress = not arg_options.batch
    pypy.stdout_stderr = Paths.temp_file('run-test-{datetime}.log')

    seq = SequentialThreads('test-case', progress=False)
    seq.add(case.create_clean_thread())
    seq.add(pypy)

    # if clean-up fails do not run other
    seq.stop_on_error = True
    return seq


def create_pbs_job_content(module, command, case):
    """
    :type module: scripts.pbs.modules.pbs_tarkil_cesnet_cz
    :rtype : str
    """
    escaped_command = ' '.join(Command.escape_command(command))
    template = PBSModule.format(
        module.template,
        command=escaped_command,
        root=Paths.base_dir(),
        output=case.job_output,
        status_ok=job_ok_string,
    )

    return template


def run_pbs_mode(all_yamls):
    # create parallel runner instance
    """
    :type all_yamls: dict[str, scripts.config.yaml_config.YamlConfig]
    """
    global arg_options, arg_others, arg_rest
    pbs_module = get_pbs_module(arg_options.host)
    Printer.dyn_enabled = not arg_options.batch
    Printer.dyn('Parsing yaml files')

    jobs = list()
    """ :type: list[(str, scripts.core.prescriptions.PBSModule)] """
    for yaml_file, config in all_yamls.items():
        # parse config.yaml

        # extract all test cases (product of cpu x files)
        for case in config.get_cases_for_file(pbs_module.Module, yaml_file):
            # create run command
            test_command = case.get_command()
            test_command.extend(arg_rest)
            pbs_content = create_pbs_job_content(pbs_module, test_command, case)
            IO.write(case.pbs_script, pbs_content)

            # create pbs file
            qsub_command = case.get_pbs_command(arg_options, case.pbs_script)
            jobs.append((qsub_command, case))

    # start jobs
    Printer.dyn('Starting jobs')

    total = len(jobs)
    job_id = 0
    multijob = MultiJob(pbs_module.ModuleJob)
    for qsub_command, case in jobs:
        job_id += 1

        Printer.dyn('Starting jobs {:02d} of {:02d}', job_id, total)

        output = subprocess.check_output(qsub_command)
        job = pbs_module.ModuleJob.create(output, case)
        job.full_name = "Case {}".format(format_case(case))
        multijob.add(job)

    Printer.out('{} job/s inserted into queue', total)

    # first update to get more info about multijob jobs
    Printer.out()
    Printer.separator()
    Printer.dyn('Updating job status')
    multijob.update()

    # print jobs statuses
    Printer.out()
    if not arg_options.batch:
        multijob.print_status()

    Printer.separator()
    Printer.dyn(multijob.get_status_line())
    returncodes = dict()

    # wait for finish
    while multijob.is_running():
        Printer.dyn('Updating job status')
        multijob.update()
        Printer.dyn(multijob.get_status_line())

        # if some jobs changed status add new line to dynamic output remains
        jobs_changed = multijob.get_all(status=JobState.COMPLETED)
        if jobs_changed:
            Printer.out()
            Printer.separator()

        # get all jobs where was status update to COMPLETE state
        for job in jobs_changed:
            returncodes[job] = finish_pbs_job(job, arg_options.batch)

        if jobs_changed:
            Printer.separator()
            Printer.out()

        # after printing update status lets sleep for a bit
        if multijob.is_running():
            time.sleep(5)
            
    Printer.out(multijob.get_status_line())
    Printer.out('All jobs finished')

    # get max return code or number 2 if there are no returncodes
    returncode = max(returncodes.values()) if returncodes else 2
    sys.exit(returncode)


def run_local_mode(all_yamls):
    # create parallel runner instance
    """
    :type all_yamls: dict[str, scripts.config.yaml_config.YamlConfig]
    """
    global arg_options, arg_others, arg_rest
    runner = ParallelThreads(arg_options.parallel)
    runner.stop_on_error = not arg_options.keep_going

    # turn on/off MPI mode
    if set(arg_options.cpu) == {1}:
        cls = prescriptions.TestPrescription
        Printer.out('Running WITHOUT MPI')
    else:
        cls = prescriptions.MPIPrescription

    # go through each yaml file
    for yaml_file, config in all_yamls.items():
        # extract all test cases (product of cpu x files)
        config.parse()
        config.update(proc=set(arg_options.cpu) if arg_options.cpu else None)
        for case in config.get_cases_for_file(cls, yaml_file):
            # create main process which first clean output dir
            # and then execute test
            multi_process = create_process_from_case(case)
            # get all comparisons threads and add them to main runner
            multi_process.add(case.create_comparison_threads())
            runner.add(multi_process)

    # run!
    runner.start()
    while runner.is_running():
        time.sleep(1)

    Printer.separator()
    Printer.out('Summary: ')
    Printer.open()
    for thread in runner.threads:
        clean, pypy, comp = getattr(thread, 'threads', [None] * 3)
        if pypy.returncode is None:
            returncode = -1
        else:
            returncode = 1 if thread.returncode is None else thread.returncode
        returncode = str(returncode)
        if pypy:
            Printer.out('[{:^3s}] {:7s}: {}', returncode,
                        PyPy.returncode_map.get(returncode, 'ERROR'),
                        format_case(pypy.case))

    Printer.close()
    # exit with runner's exit code
    sys.exit(runner.returncode)


def read_configs(all_yamls):
    all_configs = list(set([Paths.rename(y, 'config.yaml') for y in all_yamls]))
    all_configs = {c: YamlConfig(c) for c in all_configs}
    configs = {y: all_configs[Paths.rename(y, 'config.yaml')] for y in set(all_yamls)}
    for k, v in configs.items():
        v.include = arg_options.include
        v.exclude = arg_options.exclude
        # load all configs a then limit single config instances while traversing
        v.parse()

    return configs


def do_work(parser):
    """
    :type parser: utils.argparser.ArgParser
    """

    # parse arguments
    global arg_options, arg_others, arg_rest
    arg_options, arg_others, arg_rest = parser.parse()
    Paths.format = PathFormat.ABSOLUTE
    if arg_options.root:
        Paths.base_dir(arg_options.root)

    # configure printer
    Printer.batch_output = arg_options.batch
    Printer.dynamic_output = not arg_options.batch

    # we need flow123d, mpiexec and ndiff to exists in LOCAL mode
    if not arg_options.queue and not Paths.test_paths('flow123d', 'mpiexec', 'ndiff'):
        Printer.err('Missing obligatory files! Exiting')
        exit(1)

    # test yaml args
    if not arg_others:
        parser.exit_usage('Error: No yaml files or folder given')
        exit(1)

    all_yamls = list()
    for path in arg_others:
        if not Paths.exists(path):
            Printer.wrn('Warning! given path does not exists, ignoring path "{}"', path)
            continue

        if Paths.is_dir(path):
            all_yamls.extend(Paths.walk(path, filters=[
                PathFilters.filter_type_is_file(),
                PathFilters.filter_ext('.yaml'),
                PathFilters.filter_not(PathFilters.filter_name('config.yaml'))
            ]))
        else:
            all_yamls.append(path)

    Printer.out("Found {} .yaml file/s", len(all_yamls))
    if not all_yamls:
        Printer.wrn('Warning! No yaml files found in locations: \n  {}', '\n  '.join(arg_others))
        exit(1)

    all_configs = read_configs(all_yamls)

    if arg_options.queue:
        Printer.out('Running in PBS mode')
        run_pbs_mode(all_configs)
    else:
        Printer.out('Running in LOCAL mode')
        run_local_mode(all_configs)