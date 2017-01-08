#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from scripts.prescriptions import AbstractRun


runtest_command = """
"$$python$$" "$$script$$" "$$yaml$$" $$limits$$ --dump "$$dump_output$$" --log $$log_file$$ --batch -- $$args$$
""".strip()

exec_parallel_command = """
"$$python$$" "$$script$$" $$limits$$ --dump "$$dump_output$$" --log $$log_file$$ --batch -- $$args$$
""".strip()


class PBSModule(AbstractRun):
    """
    Class PBSModule server for executing PBS jobs.
    Class offers static method which formats given string with placeholders and creates dummy
    get_pbs_command method which MUST be implemented in children
    """

    def __init__(self, case):
        super(PBSModule, self).__init__(case)
        self.queue = 'default'
        self.ppn = 1

    def get_pbs_command(self, pbs_script_filename):
        """
        :rtype: list[str]
        """
        pass