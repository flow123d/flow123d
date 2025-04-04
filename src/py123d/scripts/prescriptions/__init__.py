#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from py123d.scripts.core.base import Paths, PathFilters
# ----------------------------------------------


class AbstractRun(object):
    """
    Abstract class which prepares command list for PyPy class
    :type case : scripts.yamlc.yaml_config.ConfigCase
    """

    def __init__(self, case):
        self.case = case
        self.mpi = False
        self.valgrind = False
        self.massif = False

    def configure_arguments(self, args):
        """
        Method will replace given arguments placeholders
        available placeholders are:
            - FLOW123D_DIR              - path to repository root  (such as /opt/flow123d)
            - CURRENT_TEST_DIR          - path to current test dir (such as /opt/flow123d/tests/01_cmd_line)
            - CURRENT_OUTPUT_DIR        - path to current test output dir
                                            (such as /opt/flow123d/tests/01_cmd_line/test_results/02_input_format.1)
            - CURRENT_REF_OUTPUT_DIR    - path to current test dir
                                            (such as /opt/flow123d/tests/01_cmd_line/ref_out/02_input_format)
            - TESTS_DIR                 - path to tests dir (such as /opt/flow123d/tests)
                                            NOTE:
                                                This value may not be precise and works
                                                only when running tests in standard flow123d
                                                structure.
                                                This value is essentially CURRENT_TEST_DIR/..
        :type args: list[str]
        """

        # build replacements map for the current case
        print("DBG replacment dict:")
        print("  FLOW123D_DIR:           " + str(Paths.flow123d_dir()))
        print("  CURRENT_TEST_DIR:       " + self.case.fs.root)
        print("  CURRENT_OUTPUT_DIR:     " + self.case.fs.root)
        print("  CURRENT_REF_OUTPUT_DIR: " + self.case.fs.output)
        print("  TESTS_DIR:              " + self.case.fs.root)
        replacements = dict(
            FLOW123D_DIR=str(Paths.flow123d_dir()),
            CURRENT_TEST_DIR=self.case.fs.root,
            CURRENT_OUTPUT_DIR=self.case.fs.output,
            CURRENT_REF_OUTPUT_DIR=self.case.fs.ref_output,
            TESTS_DIR=Paths.dirname(self.case.fs.root),
        )

        for i in range(len(args)):
            for repl, val in replacements.items():
                args[i] = args[i].replace('$%s$' % repl, val)  # works for $VALUE$
                args[i] = args[i].replace('<%s>' % repl, val)  # works for <VALUE>
                args[i] = args[i].replace('{%s}' % repl, val)  # works for {VALUE}
        return args

    def get_command(self, args=None):
        result = self._get_command() if not args else self._get_command() + args
        # append all arguments specified in config.yaml
        if self.case is not None and self.case.args:
            result.extend(self.case.args)

        result = [str(x) for x in result]
        return self.configure_arguments(result)

    def _get_mpi(self):
        return [
            Paths.mpiexec(),
            '-np', self.case.proc
        ]

    def _get_valgrind(self):
        if self.massif:
            return [
                'valgrind', '--tool=massif', '--massif-out-file={}'.format(self.case.fs.valgrind_out)
            ]
        if self.valgrind:
            return [
                'valgrind',
            ]
        return []

    def _get_flow123d(self):
        return [
            Paths.flow123d(),
            '-s', self.case.file,
            '-o', self.case.fs.output,
        ]

    def _get_command(self):
        result = []

        if self.mpi:
            result.extend(self._get_mpi())

        if self.valgrind or self.massif:
            result.extend(self._get_valgrind())

        result.extend(self._get_flow123d())
        return result

    def _get_all_ref_files(self):
        return Paths.walk(self.case.fs.ref_output, [PathFilters.filter_type_is_file()])

    def _get_ref_output_files(self, comp_data):
        """
        :type comp_data: dict
        Return: comparison_pairs, filtered

        comparison_pairs is list of filename pairs to compare:
        [(reference_file, output_file)]

        filtered is list of reference files not matching the the checking rule pattern.

        """
        # parse filters
        filters = [PathFilters.filter_wildcards(x) for x in comp_data.get('files', [])]

        # browse files and make them relative to ref output so filters works properly
        files = self._get_all_ref_files()
        files = [Paths.relpath(f, self.case.fs.ref_output) for f in files]

        # filter files and make them absolute again
        files = Paths.match(files, filters)
        files = [Paths.join(self.case.fs.ref_output, f) for f in files]
        return list(zip(files, self._get_mirror_files(files)))

    def _get_mirror_files(self, paths):
        return [
            Paths.join(self.case.fs.output, Paths.relpath(p, self.case.fs.ref_output))
            for p in paths
        ]
