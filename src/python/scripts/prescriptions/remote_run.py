#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from scripts.prescriptions import AbstractRun


runtest_command = """
"$$python$$" "$$script$$" "$$yaml$$" $$limits$$ --json "$$json_output$$" -- $$args$$
""".strip()

exec_parallel_command = """
"$$python$$" "$$script$$" $$limits$$ --json "$$json_output$$" -- $$args$$
""".strip()


class PBSModule(AbstractRun):
    def __init__(self, case):
        super(PBSModule, self).__init__(case)
        self.queue = 'default'
        self.ppn = 1

    def get_pbs_command(self, pbs_script_filename):
        """
        :rtype: list[str]
        """
        pass

    @staticmethod
    def format(template, **kwargs):
        t = str(template)
        for k, v in kwargs.items():
            t = t.replace('$${}$$'.format(k), v)
        return t