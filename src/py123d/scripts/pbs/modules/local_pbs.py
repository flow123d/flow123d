#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import re
import sys
import pkgutil

from py123d.scripts.pbs.job import Job, JobState
from py123d.scripts.prescriptions.remote_run import PBSModule


def determine_mock_location():
    return pkgutil.find_loader('local_pbs_server').path


class Module(PBSModule):
    """
    Class Module dummy local pbs module
    """

    # points to pbs mock python file in unit tests
    mock = determine_mock_location()

    def __init__(self, case):
        super(Module, self).__init__(case)

    def get_pbs_command(self, pbs_script_filename):
        return [
            sys.executable,
            Module.mock,
            'qsub',
            '-o', self.case.fs.pbs_output,
            pbs_script_filename]


template = """
#!/bin/bash
$$command$$

""".lstrip()


class ModuleJob(Job):
    """
    Class ModuleJob dummy local pbs module
    """

    def __init__(self, job_id, case):
        super(ModuleJob, self).__init__(job_id, case)
        self.parser = self.parser_builder(
            self, 1, JobState.UNKNOWN,
            queue=2,
        )

    @classmethod
    def update_command(cls):
        return ["python", Module.mock, 'qstat']

    @classmethod
    def create(cls, command_output, case):
        pid = re.findall(r'(\d+)', command_output)[0]
        return ModuleJob(pid, case)
