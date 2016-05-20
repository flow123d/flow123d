#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from scripts.core.prescriptions import PBSModule
from scripts.pbs.job import Job, JobState

import re


class Module(PBSModule):
    def get_pbs_command(self, options, pbs_script_filename):
        # total parallel process
        np = self.proc_value
        # processes per node, default value 2 (Master-slave)
        ppn = options.get('ppn', 2)

        # command
        command = [
            'qsub',
            '-pe', 'orte', '{np}'.format(**locals()),
            '-l', 'num_proc={ppn}'.format(**locals()),
            '-o', self.output_log,
            pbs_script_filename
        ]

        return command


class ModuleJob(Job):
    instances = list()

    def __init__(self, job_id, case):
        super(ModuleJob, self).__init__(job_id, case)
        self.parser = self.parser_builder(
            self, 4, JobState.UNKNOWN,
            queue=7,
            name=2,
        )

    @classmethod
    def update_command(cls):
        return ['qstat']

    @classmethod
    def create(cls, command_output, case):
        return ModuleJob(re.findall(r'(\d+)', command_output)[0], case)

template = """
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash

# Disable system rsh / ssh only
# export OMPI_MCA_plm_rsh_disable_qrsh=1

#################

ROOT="$$root$$"

echo "$$command$$"
# $$command$$

""".lstrip()