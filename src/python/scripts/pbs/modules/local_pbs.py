#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import re
import random

from scripts.pbs.job import Job, JobState
from scripts.prescriptions.remote_run import PBSModule


class Module(PBSModule):
    # must be set
    mock = None

    def __init__(self, case):
        super(Module, self).__init__(case)

    def get_pbs_command(self, pbs_script_filename):
        return [
            "python",
            Module.mock,
            'qsub',
            '-o', self.case.fs.pbs_output,
            pbs_script_filename]


template = """
#!/bin/bash
uname -a
echo JOB START: `date`
pwd

$$command$$

""".lstrip()


class ModuleJob(Job):
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