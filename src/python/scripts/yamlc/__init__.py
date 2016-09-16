#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from scripts.core.base import PathFilters, Paths
# ----------------------------------------------

DEFAULTS = dict(
    proc=[1],
    time_limit=30,
    memory_limit=400,
    tags=[],
    check_rules=[
        {
            'ndiff': {
                'files': ['*']
            }
        }
    ]
)

TAG_FILES = 'files'
TAG_PROC = 'proc'
TAG_TIME_LIMIT = 'time_limit'
TAG_MEMORY_LIMIT = 'memory_limit'
TAG_TEST_CASES = 'test_cases'
TAG_CHECK_RULES = 'check_rules'
TAG_TAGS = 'tags'
REF_OUTPUT_DIR = 'ref_out'

YAML = '.yaml'
CONFIG_YAML = 'config.yaml'


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

        self.input = self.in_root('input')
        self.ref_output = ref_output

    def in_root(self, *names):
        """
        :rtype: str
        """
        return Paths.join(self.root, *names)

    def in_output(self, *names):
        """
        :rtype: str
        """
        return Paths.join(self.output, *names)