#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import time
import sys
# ----------------------------------------------
from scripts.pbs import pbs_control
from scripts.yamlc.yaml_config import ConfigPool
from scripts.core.base import Paths, PathFilters, Printer, Command, IO, GlobalResult, DynamicSleep, StatusPrinter
from scripts.core.threads import ParallelThreads, RuntestMultiThread
from scripts.pbs.common import get_pbs_module
from scripts.pbs.job import JobState, MultiJob, finish_pbs_runtest
from scripts.prescriptions.local_run import LocalRun
from scripts.prescriptions.remote_run import runtest_command, PBSModule
from scripts.script_module import ScriptModule
from utils import strings
# ----------------------------------------------


class ModuleRuntest(ScriptModule):
    """
    Class ModuleRuntest is backend for script runtest.py
    """

    @staticmethod
    def list_tests():
        test_dir = Paths.join(Paths.flow123d_root(), 'tests')
        tests = Paths.walk(test_dir, [
            PathFilters.filter_type_is_file(),
            PathFilters.filter_endswith('.yaml'),
            PathFilters.filter_not(PathFilters.filter_name('config.yaml')),
        ])
        result = dict()
        for r in tests:
            dirname = Paths.dirname(r)
            basename = Paths.basename(r)
            if Paths.dirname(dirname) != test_dir:
                continue

            if dirname not in result:
                result[dirname] = list()
            result[dirname].append(basename)
        keys = sorted(result.keys())

        for dirname in keys:
            Printer.all.out(Paths.relpath(dirname, test_dir))
            with Printer.all.with_level(1):
                for basename in result[dirname]:
                    Printer.all.out('{: >4s} {: <40s} {}', '', basename, Paths.relpath(Paths.join(dirname, basename), test_dir))
            Printer.all.newline()

    @staticmethod
    def read_configs(all_yamls):
        """
        Add yamls to ConfigPool and parse configs
        :rtype: ConfigPool
        """
        configs = ConfigPool()
        for y in all_yamls:
            configs += y
        configs.parse()
        return configs

    @property
    def include(self):
        """
        return all exclude tags in set
        """
        return set(self.arg_options.include) if self.arg_options.include else set()

    @property
    def exclude(self):
        """
        return all exclude tags in set
        """
        return set(self.arg_options.exclude) if self.arg_options.exclude else {'disabled'}

    def create_process_from_case(self, case):
        """
        Method creates main thread where clean-up, pypy and comparison is stored
        :type case: scripts.yamlc.yaml_config.ConfigCase
        """
        local_run = LocalRun(case)
        local_run.mpi = case.proc > 1
        local_run.progress = self.progress

        # on demand we do not clean dirs or run comparisons
        no_clean = self.arg_options.no_clean
        no_compare = self.arg_options.no_compare

        clean = local_run.create_dummy_clean_thread() if no_clean else local_run.create_clean_thread()
        compare = local_run.create_dummy_comparisons() if no_compare else local_run.create_comparisons()
        pypy = local_run.create_pypy(self.arg_options.rest)

        # turn on status file creation on demand
        if self.arg_options.status_file:
            pypy.status_file = case.fs.status_file

        seq = RuntestMultiThread(clean, pypy, compare)

        seq.stop_on_error = True
        return seq

    def create_pbs_job_content(self, module, case):
        """
        Method creates pbs start script which will be passed to
        some qsub command

        :type case: scripts.config.yaml_config.ConfigCase
        :type module: scripts.pbs.modules.pbs_tarkil_cesnet_cz
        :rtype : str
        """

        import pkgutil

        command = strings.replace_placeholders(
            runtest_command,

            python=sys.executable,
            script=pkgutil.get_loader('runtest').filename,
            yaml=case.file,
            limits="-n {case.proc} -m {case.memory_limit} -t {case.time_limit}".format(case=case),
            args="" if not self.arg_options.rest else Command.to_string(self.arg_options.rest),
            dump_output=case.fs.dump_output,
            log_file=case.fs.job_output
        )

        template = strings.replace_placeholders(
            module.template,
            command=command,
            dump_output=case.fs.dump_output
        )

        return template

    def prepare_pbs_files(self, pbs_module):
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
        return jobs

    def run_pbs_mode(self):
        """
        Runs this module in local mode.
        At this point all configuration files has been loaded what is left
        to do is to create pbs scripts and put them to queue (qsub).
        After them we monitor all jobs (qstat) and if some job exits we parse
        result json file and determine ok/error status for the job
        """
        pbs_module = get_pbs_module(self.arg_options.host)
        jobs = self.prepare_pbs_files(pbs_module)
        result, multijob = pbs_control.do_work(jobs, pbs_module, finish_pbs_runtest)

        return result.singlify()

    def run_local_mode(self):
        """
        Runs this module in local mode.
        At this point all configuration files has been loaded what is left
        to do is to prepare execution arguments start whole process
        """
        runner = ParallelThreads(self.arg_options.parallel)
        runner.stop_on_error = not self.arg_options.keep_going

        for yaml_file, yaml_config in self.configs.files.items():
            for case in yaml_config.get_one(yaml_file):
                # create main process which first clean output dir
                # and then execute test following with comparisons
                multi_process = self.create_process_from_case(case)
                runner.add(multi_process)

        if self.include or self.exclude:
            Printer.all.out('Running {} cases ({}{})'.format(
                runner.total,
                'including only tags in set {} '.format(list(self.include)) if self.include else '',
                'excluding all tags in set {}'.format(list(self.exclude)) if self.exclude else ''))
        else:
            Printer.all.out('Running {} cases', runner.total)

        # run!
        runner.start()
        while runner.is_running():
            time.sleep(1)

        Printer.all.sep()
        Printer.all.out('Summary: ')

        with Printer.all.with_level(1):
            for thread in runner.threads:
                multithread = thread
                """ :type: RuntestMultiThread """
                StatusPrinter.print_test_result(multithread)
            Printer.all.sep()
            StatusPrinter.print_runner_stat(runner)

        # exit with runner's exit code
        GlobalResult.returncode = runner.returncode
        return runner

    def __init__(self, arg_options):
        super(ModuleRuntest, self).__init__(arg_options)
        self.all_yamls = None
        self.configs = None

    def _check_arguments(self):
        """
        Arguments additional check
        """

        if self.arg_options.list:
            self.list_tests()
            sys.exit(0)

        # we need flow123d, mpiexec and ndiff to exists in LOCAL mode
        if not self.arg_options.queue and not Paths.test_paths('flow123d', 'mpiexec', 'ndiff'):
            Printer.all.wrn('Missing obligatory files!')

    def _run(self):
        """
        Run method for this module
        """

        if self.arg_options.random_output_dir:
            import scripts.yamlc as yamlc
            yamlc.TEST_RESULTS = 'test_results-{}'.format(self.arg_options.random_output_dir)

        self.all_yamls = list()
        for path in self.arg_options.args:
            if not Paths.exists(path):
                Printer.all.err('given path does not exists, path "{}"', path)
                sys.exit(3)

            # append files to all_yamls
            if Paths.is_dir(path):
                self.all_yamls.extend(Paths.walk(path, ConfigPool.yaml_filters))
            else:
                self.all_yamls.append(path)

        Printer.all.out("Found {} yaml file/s", len(self.all_yamls))
        if not self.all_yamls:
            Printer.all.wrn('No yaml files found in locations: \n  {}', '\n  '.join(self.arg_options.args))
            sys.exit(0)

        self.configs = self.read_configs(self.all_yamls)
        self.configs.update(
            proc=self.arg_options.cpu,
            time_limit=self.arg_options.time_limit,
            memory_limit=self.arg_options.memory_limit,
        )

        # filter tags for includes and excludes
        self.configs.filter_tags(
            include=self.include,
            exclude=self.exclude
        )

        if self.arg_options.queue:
            Printer.all.out('Running in PBS mode')
            return self.run_pbs_mode()
        else:
            Printer.all.out('Running in LOCAL mode')
            return self.run_local_mode()


def do_work(arg_options=None, debug=False):
    """
    Main method which invokes ModuleRuntest
    :rtype: ParallelThreads
    :type debug: bool
    :type arg_options: utils.argparser.RuntestArgs
    """
    module = ModuleRuntest(arg_options)
    result = module.run(debug)

    # pickle out result on demand
    if arg_options.dump:
        try:
            import pickle
            pickle.dump(result.dump(), open(arg_options.dump, 'wb'))
        except: pass
    return result
