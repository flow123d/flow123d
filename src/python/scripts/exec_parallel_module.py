#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import subprocess
import datetime
import time
# ----------------------------------------------
import sys

from scripts.config.yaml_config import ConfigCase, DEFAULTS
from scripts.core.base import Paths, Printer, Command, IO, GlobalResult
from scripts.core.base import PathFormat
from scripts.core.threads import PyPy
from scripts.core.execution import BinExecutor, OutputMode
from scripts.pbs.common import get_pbs_module, job_ok_string
from scripts.pbs.job import JobState, finish_pbs_job, MultiJob
from scripts.prescriptions.remote_run import exec_parallel_command, PBSModule
from scripts.script_module import ScriptModule
from utils.counter import ProgressCounter
# ----------------------------------------------

# global arguments
from utils.globals import ensure_iterable
from utils.strings import format_n_lines


class ModuleExecParallel(ScriptModule):
    """
    Class ModuleExecParallel is backend for script exec_parallel.py
    """

    def __init__(self):
        super(ModuleExecParallel, self).__init__()
        self.proc, self.time_limit, self.memory_limit = None, None, None

    def _check_arguments(self):
        # no args
        if len(self.rest) == 0:
            self.parser.exit_usage('no MPI executable provided', exit_code=1)
        # only one arg
        if len(self.rest) == 1:
            self.parser.exit_usage('no command provided', exit_code=2)

    def _run(self):
        self.proc = ensure_iterable(self.arg_options.get('cpu'))
        self.time_limit = self.arg_options.get('time_limit')
        self.memory_limit = self.arg_options.get('memory_limit')

        # set default values if not set
        self.proc = self.proc if self.proc else DEFAULTS.get('proc')
        self.time_limit = self.time_limit if self.time_limit else DEFAULTS.get('time_limit')
        self.memory_limit = self.memory_limit if self.memory_limit else DEFAULTS.get('memory_limit')

        # run local or pbs mode
        if self.arg_options.queue:
            return self.run_pbs_mode()
        else:
            return self.run_local_mode()

    def run_pbs_mode(self):
        pbs_module = get_pbs_module()
        jobs = self.prepare_pbs_files(pbs_module)

        if self.debug:
            return 0

        # start jobs
        Printer.dyn('Starting jobs')

        total = len(jobs)
        job_id = 0
        multijob = MultiJob(pbs_module.ModuleJob)
        for qsub_command, pbs_run in jobs:
            job_id += 1

            Printer.dyn('Starting jobs {:02d} of {:02d}', job_id, total)

            output = subprocess.check_output(qsub_command)
            job = pbs_module.ModuleJob.create(output, pbs_run.case)
            job.full_name = "Case {}".format(pbs_run.case)
            multijob.add(job)

        Printer.out()
        Printer.out('{} job/s inserted into queue', total)

        # # first update to get more info about multijob jobs
        Printer.out()
        Printer.separator()
        Printer.dyn('Updating job status')
        multijob.update()

        # print jobs statuses
        Printer.out()
        if self.progress:
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
                returncodes[job] = finish_pbs_job(job, self.batch)

            if jobs_changed:
                Printer.separator()
                Printer.out()

            # after printing update status lets sleep for a bit
            if multijob.is_running():
                time.sleep(5)

        Printer.out(multijob.get_status_line())
        Printer.out('All jobs finished')

        # get max return code or number 2 if there are no returncodes
        return max(returncodes.values()) if returncodes else 2

    def prepare_pbs_files(self, pbs_module):
        """
        :type pbs_module: scripts.pbs.modules.pbs_tarkil_cesnet_cz
        :rtype: list[(str, PBSModule)]
        """

        jobs = list()

        for p in self.proc:
            case = ConfigCase(dict(
                proc=p,
                time_limit=self.time_limit,
                memory_limit=self.memory_limit,
                tmp='exec-parallel'
            ), None)

            pbs_run = pbs_module.Module(case)
            pbs_run.queue = self.arg_options.get('queue', True)
            pbs_run.ppn = self.arg_options.get('ppn', 1)

            pbs_content = self.create_pbs_job_content(pbs_module, case)
            IO.write(case.fs.pbs_script, pbs_content)

            qsub_command = pbs_run.get_pbs_command(case.fs.pbs_script)
            jobs.append((qsub_command, pbs_run))
        return jobs

    def run_local_mode(self):
        total = len(self.proc)

        if total == 1:
            pypy = self.run_local_mode_one(self.proc[0])
            GlobalResult.returncode = pypy.returncode
        else:
            # optionally we use counter
            progress = ProgressCounter('Running {:02d} of {total:02d}')
            for p in self.proc:
                Printer.separator()
                progress.next(locals())
                Printer.separator()

                Printer.open()
                pypy = self.run_local_mode_one(p)
                Printer.close()
                GlobalResult.returncode = max(GlobalResult.returncode, pypy.returncode)

        return GlobalResult.returncode if not self.debug else pypy

    def create_pbs_job_content(self, module, case):
        """
        :type case: scripts.config.yaml_config.ConfigCase
        :type module: scripts.pbs.modules.pbs_tarkil_cesnet_cz
        :rtype : str
        """

        import pkgutil

        command = PBSModule.format(
            exec_parallel_command,

            python=sys.executable,
            script=pkgutil.get_loader('exec_parallel').filename,
            limits="-n {case.proc} -m {case.memory_limit} -t {case.time_limit}".format(case=case),
            args="" if not self.rest else Command.to_string(self.rest),
            json_output=case.fs.json_output
        )

        template = PBSModule.format(
            module.template,
            command=command,
            json_output=case.fs.json_output  # TODO remove
        )

        return template

    def run_local_mode_one(self, proc):
        if self.proc == 0:
            command = self.rest[1:]
        else:
            command = [self.rest[0], '-np', proc] + self.rest[1:]

        n_lines = 0 if self.arg_options.batch else 10
        pypy = PyPy(BinExecutor(command))

        # set limits
        pypy.limit_monitor.time_limit = self.time_limit
        pypy.limit_monitor.memory_limit = self.memory_limit
        pypy.progress = self.progress
        pypy.info_monitor.deactivate()
        pypy.error_monitor.deactivate()

        # catch output to variable
        # in batch mode we will keep the files
        # otherwise we will keep logs only on error
        log_file = Paths.temp_file('exec-parallel-{date}-{time}-{rnd}.log')
        pypy.executor.output = OutputMode.variable_output()
        pypy.full_output = log_file

        # start and wait for exit
        pypy.start()
        pypy.join()

        # add result to global json result
        GlobalResult.add(pypy)

        # in batch mode or on error
        if not pypy.with_success() or self.batch:
            content = pypy.executor.output.read()
            IO.write(log_file, content)
            Printer.close()
            Printer.out(format_n_lines(content, indent='    ', n_lines=-n_lines))
            Printer.open()
        return pypy


def do_work(parser, args=None, debug=False):
    """
    Main method which invokes ModuleExecParallel
    :type debug: bool
    :type args: list
    :type parser: utils.argparser.ArgParser
    """
    module = ModuleExecParallel()
    return module.run(parser, args, debug)
