#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import subprocess
import time
import sys
# ----------------------------------------------
from scripts.config.yaml_config import ConfigPool
from scripts.core.base import Paths, PathFilters, Printer, Command, IO, GlobalResult
from scripts.core.threads import ParallelThreads, RuntestMultiThread
from scripts.pbs.common import get_pbs_module
from scripts.pbs.job import JobState, MultiJob, finish_pbs_job
from scripts.prescriptions.local_run import LocalRun
from scripts.prescriptions.remote_run import runtest_command, PBSModule
from scripts.script_module import ScriptModule
# ----------------------------------------------


class ModuleRuntest(ScriptModule):
    @staticmethod
    def read_configs(all_yamls):
        """
        :rtype: scripts.config.yaml_config.ConfigPool
        """
        configs = ConfigPool()
        for y in all_yamls:
            configs += y
        configs.parse()
        return configs

    def create_process_from_case(self, case):
        """
        :type case: scripts.config.yaml_config.ConfigCase
        """
        local_run = LocalRun(case)
        local_run.mpi = case.proc > 1
        local_run.progress = self.progress

        seq = RuntestMultiThread(
            local_run.create_clean_thread(),
            local_run.create_pypy(self.rest),
            local_run.create_comparisons()
        )

        seq.stop_on_error = True
        return seq

    def create_pbs_job_content(self, module, case):
        """
        :type case: scripts.config.yaml_config.ConfigCase
        :type module: scripts.pbs.modules.pbs_tarkil_cesnet_cz
        :rtype : str
        """

        import pkgutil

        command = PBSModule.format(
            runtest_command,

            python=sys.executable,
            script=pkgutil.get_loader('runtest').filename,
            yaml=case.file,
            limits="-n {case.proc} -m {case.memory_limit} -t {case.time_limit}".format(case=case),
            args="" if not self.rest else Command.to_string(self.rest),
            json_output=case.fs.json_output
        )

        template = PBSModule.format(
            module.template,
            command=command,
            json_output=case.fs.json_output
        )

        return template

    def run_pbs_mode(self):
        """
        :type debug: bool
        :type configs: scripts.config.yaml_config.ConfigPool
        """
        pbs_module = get_pbs_module(self.arg_options.host)
        Printer.dynamic_output = not self.arg_options.batch
        Printer.dyn('Parsing yaml files')

        jobs = list()
        """ :type: list[(str, PBSModule)] """

        for yaml_file, yaml_config in self.configs.files.items():
            for case in yaml_config.get_one(yaml_file):
                pbs_run = pbs_module.Module(case)
                pbs_run.queue = self.arg_options.get('queue', True)
                pbs_run.ppn = self.arg_options.get('ppn', 1)

                pbs_content = self.create_pbs_job_content(pbs_module, case)
                IO.write(case.fs.pbs_script, pbs_content)

                qsub_command = pbs_run.get_pbs_command(case.fs.pbs_script)
                jobs.append((qsub_command, pbs_run))

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
        if not self.arg_options.batch:
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
                returncodes[job] = finish_pbs_job(job, self.arg_options.batch)

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

    def run_local_mode(self):
        """
        :type debug: bool
        :type configs: scripts.config.yaml_config.ConfigPool
        """
        runner = ParallelThreads(self.arg_options.parallel)
        runner.stop_on_error = not self.arg_options.keep_going

        for yaml_file, yaml_config in self.configs.files.items():
            for case in yaml_config.get_one(yaml_file):
                # create main process which first clean output dir
                # and then execute test following with comparisons
                multi_process = self.create_process_from_case(case)
                runner.add(multi_process)

        # run!
        runner.start()
        while runner.is_running():
            time.sleep(1)

        Printer.separator()
        Printer.out('Summary: ')
        Printer.open()

        for thread in runner.threads:
            multithread = thread
            """ :type: RuntestMultiThread """

            returncode = multithread.returncode
            GlobalResult.add(multithread)

            if multithread.clean.with_error():
                Printer.out("[{:^6}]:{:3} | Could not clean directory '{}': {}", 'ERROR',
                            multithread.clean.returncode,
                            multithread.clean.dir,
                            multithread.clean.error)
                continue

            if not multithread.pypy.with_success():
                Printer.out("[{:^6}]:{:3} | Run error, case: {}",
                            multithread.pypy.returncode_map.get(str(multithread.pypy.returncode), 'ERROR'),
                            multithread.pypy.returncode, multithread.pypy.case.to_string())
                continue

            if multithread.comp.with_error():
                Printer.out("[{:^6}]:{:3} | Compare error, case: {}, Details: ",
                            'FAILED', multithread.comp.returncode, multithread.pypy.case.to_string())
                Printer.open(2)
                for t in multithread.comp.threads:
                    if t:
                        Printer.out('[{:^6}]: {}', 'OK', t.name)
                    else:
                        Printer.out('[{:^6}]: {}', 'FAILED', t.name)
                Printer.close(2)
                continue

            Printer.out("[{:^6}]:{:3} | Test passed: {}",
                        'PASSED', multithread.pypy.returncode, multithread.pypy.case.to_string())
        Printer.close()

        # exit with runner's exit code
        GlobalResult.returncode = runner.returncode
        return runner if self.debug else runner.returncode

    def __init__(self):
        super(ModuleRuntest, self).__init__()
        self.all_yamls = None
        self.configs = None

    def _check_arguments(self):
        # we need flow123d, mpiexec and ndiff to exists in LOCAL mode
        if not self.arg_options.queue and not Paths.test_paths('flow123d', 'mpiexec', 'ndiff'):
            Printer.err('Missing obligatory files! Exiting')
            GlobalResult.error = "missing obligatory files"
            sys.exit(1)

        # test yaml args
        if not self.others:
            self.parser.exit_usage('Error: No yaml files or folder given')
            GlobalResult.error = "no yaml files or folder given"
            sys.exit(2)

    def _run(self):

        self.all_yamls = list()
        for path in self.others:
            if not Paths.exists(path):
                Printer.err('Error! given path does not exists, ignoring path "{}"', path)
                GlobalResult.error = "path does not exist"
                sys.exit(3)

            if Paths.is_dir(path):
                self.all_yamls.extend(Paths.walk(path, filters=[
                    PathFilters.filter_type_is_file(),
                    PathFilters.filter_ext('.yaml'),
                    PathFilters.filter_not(PathFilters.filter_name('config.yaml'))
                ]))
            else:
                self.all_yamls.append(path)

        Printer.out("Found {} .yaml file/s", len(self.all_yamls))
        if not self.all_yamls:
            Printer.wrn('Warning! No yaml files found in locations: \n  {}', '\n  '.join(self.others))
            GlobalResult.error = "no yaml files or folders given"
            sys.exit(3)

        self.configs = self.read_configs(self.all_yamls)
        self.configs.update(
            proc=self.arg_options.cpu,
            time_limit=self.arg_options.time_limit,
            memory_limit=self.arg_options.memory_limit,
        )

        if self.arg_options.queue:
            Printer.out('Running in PBS mode')
            return self.run_pbs_mode()
        else:
            Printer.out('Running in LOCAL mode')
            return self.run_local_mode()


def do_work(parser, args=None, debug=False):
    """
    :type debug: bool
    :type args: list
    :type parser: utils.argparser.ArgParser
    """
    module = ModuleRuntest()
    return module.run(parser, args, debug)
