#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from py123d.scripts.core.base import Paths
# ----------------------------------------------

"""
{
    proc=[1],                   # array of int, where each int represents
                                # number of cores which will be used
    time_limit=30,              # number of seconds (give or take resolution) allowed for test to run
    memory_limit=400,           # number of MB allowed for test to use (include all children processes)
    tags=[],                    # array of string, representing tags for case
    args=[],                    # additional arguments passed to binary
    check_rules=[               # array of objects, representing check rules
        {
            'ndiff': {          # key is tool which will be used
                'files': ['*']  # regexp matching for files which will be compared using this tool
            }
        }
    ]
}


"""


class YamlDeathTest(object):
    TRUE = True
    ANY = 1
    ALL = 2
    FALSE = False

    def __init__(self, value):

        if value == self.TRUE:
            self.value = self.TRUE
        elif value == self.FALSE:
            self.value = self.FALSE
        elif str(value).lower() == 'any':
            self.value = self.ANY
        elif str(value).lower() == 'all':
            self.value = self.ALL

    def reverse_return_code(self, returncode):
        """
        :type returncode: scripts.core.returncode.RC
        """
        if self.value == self.TRUE:
            # do not change None return code (test was skipped)
            if returncode is None:
                return None

            # on success return code 100 indicating manual death_error value
            if returncode == 0:
                return 100

            # otherwise return 0
            return 0

        # no reversal
        if self.value == self.FALSE:
            return returncode

        raise NotImplementedError("Value %s is not implemented yet" % str(self.value))

TAG_FILES = 'files'
TAG_PROC = 'proc'
TAG_TIME_LIMIT = 'time_limit'
TAG_MEMORY_LIMIT = 'memory_limit'
TAG_TEST_CASES = 'test_cases'
TAG_CHECK_RULES = 'check_rules'
TAG_TAGS = 'tags'
TAG_ARGS = 'args'
TAG_DEATH_TEST = 'death_test'
REF_OUTPUT_DIR = 'ref_out'
TEST_RESULTS = 'test_results'

YAML = '.yaml'
CONFIG_YAML = 'config.yaml'

"""
TODO: Remove defaults. As hot fix the check_rules default was removed. This
cause failure if any reference file is present as there is no rule for it. 
"""
DEFAULTS = {
    TAG_PROC:           [1],
    TAG_TIME_LIMIT:     30,
    TAG_MEMORY_LIMIT:   400,
    TAG_DEATH_TEST:     YamlDeathTest.FALSE,
    TAG_TAGS:           [],
    TAG_ARGS:           [],
    TAG_CHECK_RULES: []
}


class ConfigCaseFiles(object):
    """
    Class ConfigCaseFiles is helper class for defining path for ConfigCase
    """

    def __init__(self, root, ref_output, output):
        """
        :type ref_output: str
        :type output: str
        :type root: str
        """
        self.root = root
        self.output = output
        self.ndiff_log = self.in_output('ndiff.log')

        self.pbs_script = self.in_output('pbs_script.qsub')
        self.pbs_output = self.in_output('pbs_output.log')

        self.job_output = self.in_output('job_output.log')
        self.json_output = self.in_output('result.json')
        self.dump_output = self.in_output('result.p')
        self.status_file = self.in_output('runtest.status.json')
        self.valgrind_out = self.in_output('valgrind.out')

        self.ref_output = ref_output

    def in_root(self, *names):
        """
        Will return path for file located in root
        :rtype: str
        """
        return Paths.join(self.root, *names)

    def in_output(self, *names):
        """
        Will return path for file located in output
        :rtype: str
        """
        return Paths.join(self.output, *names)

    def __repr__(self):
        vars = ('root', 'output', 'ndiff_log', 'pbs_script', 'pbs_output',
                'job_output', 'json_output', 'dump_output', 'status_file',
                'valgrind_out', 'ref_output')
        return ' - ' + '\n - '.join(['{:>15s}: {:s}'.format(x, getattr(self, x)) for x in vars])
