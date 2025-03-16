#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from py123d.scripts.prescriptions import AbstractRun


runtest_command = """
"$$python$$" "$$script$$" "$$yaml$$" $$limits$$ --dump "$$dump_output$$" $$save_to_db$$ $$status_file$$ --log $$log_file$$ $$random_output_dir$$ --batch -- $$args$$
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
        self.limits = list()
        self._queue = False

    @property
    def queue(self):
        return self._queue

    @queue.setter
    def queue(self, value):
        if type(value) is str:
            self._queue = value


    def get_pbs_command(self, pbs_script_filename):
        """
        :rtype: list[str]
        """
        pass


class PBSLimit(object):
    def __init__(self, value, type='-l'):
        self.value = value
        self.type = type

    def __repr__(self):
        return self.type + ' ' + self.value

    def __call__(self, *args, **kwargs):
        return [self.type, self.value]
