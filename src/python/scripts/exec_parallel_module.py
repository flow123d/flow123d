#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import subprocess
import datetime
import time
# ----------------------------------------------
import sys

from scripts.pbs import pbs_control
from scripts.pbs.job import finish_pbs_exec
from scripts.yamlc.yaml_config import ConfigCase
from scripts.yamlc import DEFAULTS
from scripts.core.base import Paths, Printer, Command, IO
from scripts.core.threads import PyPy, ResultHolder
from scripts.core.execution import BinExecutor, OutputMode
from scripts.pbs.common import get_pbs_module
from scripts.prescriptions.remote_run import exec_parallel_command, PBSModule
from scripts.script_module import ScriptModule
from utils.counter import ProgressCounter
from utils import strings
# ----------------------------------------------

# global arguments
from utils.globals import ensure_iterable


class ModuleExecParallel(ScriptModule):
    """
    Class ModuleExecParallel is backend for script exec_parallel.py
    """

    def __init__(self, arg_options):
        super(ModuleExecParallel, self).__init__(arg_options)
        self.proc, self.time_limit, self.memory_limit = None, None, None

    def _check_arguments(self):
        """
        Arguments additional check
        """
        return

    def _run(self):
        """
        Run method for this module
        """

        self.proc = ensure_iterable(self.arg_options.cpu)
        self.time_limit = self.arg_options.time_limit
        self.memory_limit = self.arg_options.memory_limit

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
        """
        Runs this module in local mode.
        At this point arguments from cmd where parsed and converted.
        We create pbs start script/s and monitor their status/es. All jobs
        are inserted into queue and after finished evaluated.
        """
        pbs_module = get_pbs_module(self.arg_options.host)
        jobs = self.prepare_pbs_files(pbs_module)
        result, multijob = pbs_control.do_work(jobs, pbs_module, finish_pbs_exec)

        return result.singlify()

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
        """
        Runs this module in local mode.
        At this point arguments from cmd where parsed and converted.
        We create command arguments and prepare threads.

        If we need to run multiple jobs (for example cpu specification required
        to run this command using 1 AND 2 CPU) we call method run_local_mode_one
        repeatedly otherwise only once
        :rtype: PyPy
        """
        total = len(self.proc)
        result = ResultHolder()

        if total == 1:
            pypy = self.run_local_mode_one(self.proc[0])
            result.add(pypy)
        else:
            # optionally we use counter
            progress = ProgressCounter('Running {:02d} of {total:02d}')
            for p in self.proc:
                progress.next(locals())
                Printer.all.sep()

                with Printer.all.with_level():
                    pypy = self.run_local_mode_one(p)
                result.add(pypy)
                Printer.all.sep()

        return result.singlify()

    def create_pbs_job_content(self, module, case):
        """
        :type case: scripts.config.yaml_config.ConfigCase
        :type module: scripts.pbs.modules.pbs_tarkil_cesnet_cz
        :rtype : str
        """

        import pkgutil

        command = strings.replace_placeholders(
            exec_parallel_command,

            python=sys.executable,
            script=pkgutil.get_loader('exec_parallel').filename,
            limits="-n {case.proc} -m {case.memory_limit} -t {case.time_limit}".format(case=case),
            args="" if not self.arg_options.rest else Command.to_string(self.arg_options.rest),
            dump_output=case.fs.dump_output,
            log_file=case.fs.job_output
        )

        template = strings.replace_placeholders(
            module.template,
            command=command,
            dump_output=case.fs.dump_output  # TODO remove
        )

        return template

    def run_local_mode_one(self, proc):
        """
        Method runs single job with specified number of CPU
        :param proc:
        """
        if int(proc) == 0:
            command = self.arg_options.rest[1:]
        else:
            command = [self.arg_options.rest[0], '-np', proc] + self.arg_options.rest[1:]

        n_lines = 0 if self.arg_options.batch else 10
        pypy = PyPy(BinExecutor(command))

        # set limits
        pypy.limit_monitor.time_limit = self.time_limit
        pypy.limit_monitor.memory_limit = self.memory_limit

        # catch output to variable
        # in batched mode we will keep the files
        # otherwise we will keep logs only on error
        log_file = Paths.temp_file('exec-parallel-{date}-{time}-{rnd}.log')
        pypy.executor.output = OutputMode.variable_output()
        pypy.full_output = log_file

        # save output to file
        pypy.output_monitor.log_file = log_file

        # start and wait for exit
        pypy.start()
        pypy.join()

        return pypy


def do_work(arg_options=None, debug=False):
    """
    Main method which invokes ModuleExecParallel
    :rtype: scripts.core.threads.ResultHolder or
            scripts.serialization.ResultHolderResult or
            scripts.serialization.PyPyResult or
            scripts.core.threads.PyPy
    :type debug: bool
    :type arg_options: utils.argparser.ExecParallelArgs
    """
    module = ModuleExecParallel(arg_options)
    result = module.run(debug)

    # pickle out result on demand
    if arg_options.dump:
        try:
            import pickle
            pickle.dump(result.dump(), open(arg_options.dump, 'wb'))
        except:
            pass # TODO implement dump in pbs mode
    return result
