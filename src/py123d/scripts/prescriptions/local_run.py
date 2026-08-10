#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import importlib
import shutil

import py123d.scripts.comparisons.modules as modules
# ----------------------------------------------
from py123d.loggers import printf
from py123d.scripts.core.base import Paths
from py123d.scripts.core.execution import BinExecutor, OutputMode
from py123d.scripts.core.pypy import PyPy
from py123d.scripts.core.returncode import RC, RC_OK
from py123d.scripts.core.threads import ComparisonMultiThread, ExtendedThread
from py123d.scripts.prescriptions import AbstractRun
from py123d.scripts.yamlc import REF_OUTPUT_DIR


# ----------------------------------------------

# this format is used if we are not running in TTY (meaning Jenkins or tee)
# this line will be printed once
_start_format = (
    'exec [{monitor.pypy.extra[index]:02d} of {monitor.pypy.extra[runner].total:02d}] '
    '{monitor.pypy} '
    '[{monitor.process_runtime_formatted}] '
    '[{monitor.process_memory_usage_formatted}]'
)

# this format is used if we are running in TTY (meaning user, or better terminal)
# this line will be updated over and over, until process ends (or is terminated)
_update_format = (
    '{monitor.animation} [{monitor.pypy.extra[index]:02d} of {monitor.pypy.extra[runner].total:02d}] '
    '{monitor.pypy} '
    '[{monitor.process_runtime_formatted}] '
    '[{monitor.process_memory_usage_formatted}]'
)

# this format will be used in both TTY and non TTY mode
# it'll be printed always
#   in TTY mode it will rewrite last _update_format line
#   in non TTY mode this line will be appended
_complete_format = (
    'done [{monitor.pypy.extra[index]:02d} of {monitor.pypy.extra[runner].total:02d}] '
    '{monitor.pypy} '
    '[{monitor.process_runtime_formatted}] '
    '[{monitor.process_memory_used_formatted}] exitcode {monitor.pypy.returncode}'
)


class LocalRun(AbstractRun):
    """
    Class LocalRun creates PyPy object and creates comparison threads
    """

    module_path = 'py123d.scripts.comparisons.modules'

    def __init__(self, case):
        super(LocalRun, self).__init__(case)

    @classmethod
    def get_module(cls, compare_method):
        """
        :rtype: scripts.comparisons.modules.ExecComparison | scripts.comparisons.modules.InPlaceComparison
        """
        try:
            module_path = f'{cls.module_path}.{compare_method}'
            package = importlib.import_module(module_path)
            return getattr(package, compare_method)()
        except Exception as e:
            print("comparison module: ", module_path)
            print(e)
            return None

    def create_pypy(self, arg_rest):
        executor = BinExecutor(self.get_command(arg_rest))
        pypy = PyPy(executor)

        if printf.tty():
            pypy.monitor.update_format = _update_format
        else:
            pypy.monitor.start_format = _start_format

        pypy.monitor.color_complete_format = _complete_format
        pypy.case = self.case
        pypy.monitor.set_limits(self.case)
        pypy.executor.output = OutputMode.file_write(self.case.fs.job_output)
        pypy.full_output = pypy.executor.output.filename

        if self.massif:
            import scripts.prescriptions.modules.valgrind as valgrind
            pypy.on_process_complete += valgrind.massif_hook

        return pypy

    def create_comparisons(self):
        comparisons = ComparisonMultiThread(
            self.case.fs.ndiff_log,
            progress=printf.verbosity() is printf.OutputVerbosity.FULL
        )

        all_ref_files = set(self._get_all_ref_files())
        for check_rule in self.case.check_rules:
            method = str(list(check_rule.keys())[0])
            module = self.get_module(method)
            check_rule_config = check_rule[method]
            if not module:
                printf.error('Warning! No module for check_rule method "{}"', method)
                continue

            comparison_pairs = self._get_ref_output_files(check_rule_config) # [(ref_file, output_file), ...]
            if comparison_pairs:
                for pair in comparison_pairs:
                    all_ref_files.discard(pair[0])

                    # load module and determine whether we are dealing with
                    # exec comparison or inplace comparison
                    if issubclass(module.__class__, modules.ExecComparison):
                        command = module.get_command(*pair, **check_rule_config)
                        pm = PyPy(BinExecutor(command))
                        pm.executor.output = OutputMode.variable_output()
                    else:
                        module = self.get_module(method)
                        module.prepare(*pair, **check_rule_config)
                        pm = PyPy(module)
                        pm.executor.output = OutputMode.dummy_output()
                        # pm.error_monitor.deactivate()

                    # if we fail, set error to 13
                    pm.custom_error = 13
                    # TODO: maybe some time limit would be useful
                    pm.full_output = self.case.fs.ndiff_log

                    path = Paths.path_end_until(pair[0], REF_OUTPUT_DIR)
                    test_name = Paths.basename(Paths.dirname(Paths.dirname(self.case.fs.ref_output)))
                    size = Paths.filesize(pair[0], True)
                    pm.name = '{}: {} ({})'.format(test_name, path, size)

                    if printf.verbosity() is printf.OutputVerbosity.FULL:
                        pm.monitor.color_complete_format = '{}: {} ({})'.format(test_name, path, size)
                    else:
                        pm.monitor.error_complete_format = '{}: {} ({})'.format(test_name, path, size)

                    comparisons.add(pm)
        if all_ref_files:
            message = ["Missing comparison rule for reference files:"]
            message.extend([f"    {f}" for f in all_ref_files])
            raise Exception("\n".join(message))

        return comparisons

    def create_clean_thread(self):
        return CleanThread("cleaner", self.case.fs.output)

    def create_dummy_clean_thread(self):
        return DummyCleanThread("cleaner", self.case.fs.output)

    def create_dummy_comparisons(self):
        return DummyComparisonThread(self.case.fs.ndiff_log)


class CleanThread(ExtendedThread):
    """
    Class CleanThread clean directory where results will be stored
    """

    def __init__(self, name, dir):
        super(CleanThread, self).__init__(name)
        self.dir = dir
        self.error = None
        self.returncode = RC_OK

    def _run(self):
        if Paths.exists(self.dir):
            try:
                shutil.rmtree(self.dir)
                self.returncode = RC_OK
            except OSError as e:
                self.returncode = RC(4)
                self.error = str(e)

    def to_json(self):
        json = super(CleanThread, self).to_json()
        json['dir'] = self.dir
        json['error'] = self.error
        return json


class DummyCleanThread(CleanThread):
    def _run(self):
        self.returncode = RC_OK


class DummyComparisonThread(ComparisonMultiThread):
    def _run(self):
        return

    @property
    def returncode(self):
        return RC_OK
