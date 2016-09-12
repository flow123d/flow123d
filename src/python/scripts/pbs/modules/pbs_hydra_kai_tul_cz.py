#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import re
# ----------------------------------------------
from scripts.pbs.job import Job, JobState
from scripts.prescriptions.remote_run import PBSModule
# ----------------------------------------------


class Module(PBSModule):
    """
    Class Module for metacentrum PBS system
    """

    def get_pbs_command(self, pbs_script_filename):
        # total parallel process
        np = self.case.proc

        # processes per node, default value 2 (Master-slave)
        ppn = self.ppn

        # command
        command = [
            'qsub',
            '-pe', 'orte', '{np}'.format(**locals()),
            '-l', 'num_proc={ppn}'.format(**locals()),
            '-o', self.case.fs.pbs_output,
            pbs_script_filename
        ]

        return command


class ModuleJob(Job):
    """
    Class ModuleJob for metacentrum PBS system job
    """

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


# header - info about task
uname -a
echo JOB START: `date`
pwd

$$command$$

""".lstrip()
