#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import shutil
# ----------------------------------------------
from scripts.core.base import Paths, Printer
from scripts.core.threads import PyPy, ExtendedThread, ComparisonMultiThread
from scripts.core.execution import BinExecutor, OutputMode
from scripts.prescriptions import AbstractRun
from scripts.comparisons import file_comparison
from scripts.yamlc import REF_OUTPUT_DIR
# ----------------------------------------------


class LocalRun(AbstractRun):
    """
    Class LocalRun creates PyPy object and creates comparison threads
    """

    def __init__(self, case):
        super(LocalRun, self).__init__(case)
        self.progress = False

    def create_pypy(self, arg_rest):
        executor = BinExecutor(self.get_command(arg_rest))
        pypy = PyPy(executor)
        pypy.case = self.case

        pypy.limit_monitor.set_limits(self.case)
        pypy.end_monitor.deactivate()
        pypy.start_monitor.format = 'Running: {}'.format(self.case)

        pypy.progress = self.progress
        pypy.executor.output = OutputMode.file_write(self.case.fs.job_output)
        pypy.full_output = pypy.executor.output.filename
        return pypy

    def create_comparisons(self):
        comparisons = ComparisonMultiThread(self.case.fs.ndiff_log)
        comparisons.thread_name_property = True

        for check_rule in self.case.check_rules:
            method = str(check_rule.keys()[0])
            module = getattr(file_comparison, 'Compare{}'.format(method.capitalize()), None)
            comp_data = check_rule[method]
            if not module:
                Printer.all.err('Warning! No module for check_rule method "{}"', method)
                continue

            pairs = self._get_ref_output_files(comp_data)
            if pairs:
                for pair in pairs:
                    command = module.get_command(*pair, **comp_data)
                    pm = PyPy(BinExecutor(command), progress=True)

                    # if we fail, set error to 13
                    pm.custom_error = 13
                    pm.start_monitor.deactivate()
                    pm.end_monitor.deactivate()
                    pm.progress_monitor.deactivate()
                    pm.limit_monitor.deactivate() # TODO: maybe some time limit would be useful
                    pm.output_monitor.policy = pm.output_monitor.POLICY_ERROR_ONLY

                    pm.error_monitor.message = 'Comparison using method {} failed!'.format(method)
                    pm.error_monitor.indent = 1

                    # catch output
                    pm.executor.output = OutputMode.variable_output()
                    pm.full_output = self.case.fs.ndiff_log

                    path = Paths.path_end_until(pair[0], REF_OUTPUT_DIR)
                    test_name = Paths.basename(Paths.dirname(Paths.dirname(self.case.fs.ref_output)))
                    size = Paths.filesize(pair[0], True)
                    pm.name = '{}: {} ({})'.format(test_name, path, size)
                    comparisons.add(pm)

        return comparisons

    def create_clean_thread(self):
        return CleanThread("cleaner", self.case.fs.output)


class CleanThread(ExtendedThread):
    """
    Class CleanThread clean directory where results will be stored
    """

    def __init__(self, name, dir):
        super(CleanThread, self).__init__(name)
        self.dir = dir
        self.error = None
        self.returncode = 0

    def _run(self):
        if Paths.exists(self.dir):
            try:
                shutil.rmtree(self.dir)
                self.returncode = 0
            except OSError as e:
                self.returncode = 4
                self.error = str(e)

    def to_json(self):
        json = super(CleanThread, self).to_json()
        json['dir'] = self.dir
        json['error'] = self.error
        return json
