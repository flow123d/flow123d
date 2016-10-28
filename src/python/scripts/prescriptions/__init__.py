#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from scripts.core.base import Paths, PathFilters
# ----------------------------------------------


class AbstractRun(object):
    """
    Abstract class which prepares command list for PyPy class
    :type case : scripts.config.yaml_config.ConfigCase
    """

    def __init__(self, case):
        self.case = case
        self.mpi = False
        self.valgrind = False

    def get_command(self, args=None):
        result = self._get_command() if not args else self._get_command() + args
        return [str(x) for x in result]

    def _get_mpi(self):
        return [
            Paths.mpiexec(),
            '-np', self.case.proc
        ]

    def _get_valgrind(self):
        return [
            'valgrind',
        ]

    def _get_flow123d(self):
        return [
            Paths.flow123d(),
            '-s', self.case.file,
            '-i', self.case.fs.input,
            '-o', self.case.fs.output
        ]

    def _get_command(self):
        result = []

        if self.mpi:
            result.extend(self._get_mpi())

        if self.valgrind:
            result.extend(self._get_valgrind())

        result.extend(self._get_flow123d())
        return result

    def _get_ref_output_files(self, comp_data):
        """
        :type comp_data: dict
        """
        # parse filters
        filters = [PathFilters.filter_wildcards(x) for x in comp_data.get('files', [])]

        # browse files and make them relative to ref output so filters works properly
        files = Paths.walk(self.case.fs.ref_output, [PathFilters.filter_type_is_file()])
        files = [Paths.relpath(f, self.case.fs.ref_output) for f in files]

        # filter files and make them absolute again
        files = Paths.match(files, filters)
        files = [Paths.join(self.case.fs.ref_output, f) for f in files]
        return zip(files, self._get_mirror_files(files))

    def _get_mirror_files(self, paths):
        return [
            Paths.join(self.case.fs.output, Paths.relpath(p, self.case.fs.ref_output))
            for p in paths
        ]
