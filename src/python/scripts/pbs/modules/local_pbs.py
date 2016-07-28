#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import re
import random

from scripts.pbs.job import Job
from scripts.prescriptions.remote_run import PBSModule


class Module(PBSModule):
    def __init__(self, case):
        super(Module, self).__init__(case)

    def get_pbs_command(self, pbs_script_filename):
        return ["bash", pbs_script_filename]


template = """
#!/bin/bash
# $$command$$
cp /home/jan-hybs/Dokumenty/projects/Flow123d-python-utils/src/foo.json $$json_output$$ > /dev/null 2>&1
echo 25309
""".lstrip()


class ModuleJob(Job):
    def __init__(self, job_id, case):
        super(ModuleJob, self).__init__(job_id, case)

    @classmethod
    def update_command(cls):
        return ['ps', '-a']

    def parse_status(self, output=""):
        # if random.random() > 0.5:
        #     return 'Q'
        if output.find(str(self.id)) == -1:
            return 'C'
        return 'R'

    @classmethod
    def create(cls, command_output, case):
        pid = re.findall(r'(\d+)', command_output)[0]
        return ModuleJob(pid, case)