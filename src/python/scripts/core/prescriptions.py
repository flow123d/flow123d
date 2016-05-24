#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import shutil
from scripts.core.base import Paths, Printer, PathFilters
from scripts.core.threads import BinExecutor, SequentialThreads, ExtendedThread, PyPy
from scripts.comparisons import file_comparison


class TestPrescription(object):
    def __init__(self, test_case, proc_value, filename):
        """
        :type filename: str
        :type proc_value: int
        :type test_case: scripts.config.yaml_config.YamlConfigCase
        """
        self.test_case = test_case
        self.proc_value = proc_value
        self.filename = filename
        self.printer = Printer(Printer.LEVEL_KEY)

        if not self.filename:
            return

        self.shortname = Paths.basename(Paths.without_ext(self.filename))
        self.ref_output = Paths.join(test_case.config.ref_output, self.shortname)
        self.output_name = '_{}.{}'.format(self.shortname, self.proc_value)
        self.output_dir = Paths.join(test_case.config.test_results, self.output_name)
        self.ndiff_log = Paths.join(self.output_dir, 'ndiff.log')
        self.pbs_script = Paths.join(self.output_dir, 'pbs_script.qsub')
        self.output_log = Paths.join(self.output_dir, 'pbs_job.log')

    def _get_command(self):
        return [
            Paths.flow123d(),
            '-s', self.filename,
            '-i', self.test_case.config.input,
            '-o', self.output_dir
        ]

    def get_command(self):
        return [str(x) for x in self._get_command()]

    def __repr__(self):
        return '<{self.__class__.__name__} {self.output_name}>'.format(self=self)

    def _get_mirror_files(self, paths):
        return [
            Paths.join(self.output_dir, Paths.relpath(p, self.ref_output))
            for p in paths
        ]

    def get_ref_output_files(self, comp_data):
        """
        :type comp_data: dict
        """
        # parse filters
        filters = [PathFilters.filter_wildcards(x) for x in comp_data.get('files', [])]

        # browse files and make them relative to ref output so filters works properly
        files = Paths.walk(self.ref_output, [PathFilters.filter_type_is_file()])
        files = [Paths.relpath(f, self.ref_output) for f in files]

        # filter files and make them absolute again
        files = Paths.match(files, filters)
        files = [Paths.join(self.ref_output, f) for f in files]
        return zip(files, self._get_mirror_files(files))

    def create_clean_thread(self):
        def target():
            if Paths.exists(self.output_dir):
                self.printer.dbg('Cleaning output dir {}'.format(self.output_dir))
                shutil.rmtree(self.output_dir)
        return ExtendedThread(name='clean', target=target)

    def create_comparison_threads(self):
        compares = SequentialThreads(name='Comparison', progress=True, indent=True)
        compares.thread_name_property = True

        for check_rule in self.test_case.check_rules:

            method = str(check_rule.keys()[0])
            module = getattr(file_comparison, 'Compare{}'.format(method.capitalize()), None)
            comp_data = check_rule[method]
            if not module:
                self.printer.dbg('Warning! No module for check_rule method "{}"', method)
                continue

            pairs = self.get_ref_output_files(comp_data)
            if pairs:
                for pair in pairs:
                    command = module.get_command(*pair, **comp_data)
                    pm = PyPy(BinExecutor(command), progress=True)

                    pm.info_monitor.active = False
                    pm.limit_monitor.active = False
                    pm.progress_monitor.active = False
                    pm.error_monitor.message = 'Error! Comparison using method {} failed!'.format(method)
                    pm.stdout_stderr = self.ndiff_log

                    path = Paths.path_end_until(pair[0], 'ref_output')
                    test_name = Paths.basename(Paths.dirname(Paths.dirname(self.ref_output)))
                    size = Paths.filesize(pair[0], True)
                    pm.name = '{}: {} ({})'.format(test_name, path, size)
                    compares.add(pm)
        return compares


class MPIPrescription(TestPrescription):
    def __init__(self, test_case, proc_value, filename):
        super(MPIPrescription, self).__init__(test_case, proc_value, filename)

    def _get_command(self):
        return [
            Paths.mpiexec(),
            '-np', self.proc_value
        ] + super(MPIPrescription, self)._get_command()


class PBSModule(TestPrescription):
    def _get_command(self):
        return [
            'mpirun',
            '-n', self.proc_value
        ] + super(PBSModule, self)._get_command()

    def get_pbs_command(self, options, pbs_script_filename):
        """
        Method will generate all command which will then create PBS job
        :type options: utils.argparser.ArgOptions
        """
        raise NotImplementedError('Method must be implemented in sub classes')

    @staticmethod
    def format(template, **kwargs):
        t = str(template)
        for k, v in kwargs.items():
            t = t.replace('$${}$$'.format(k), v)
        return t
