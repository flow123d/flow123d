#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import math
import re
import datetime
import getpass
# ----------------------------------------------
from scripts.core.base import Printer
from scripts.core.prescriptions import PBSModule
from scripts.pbs.job import Job, JobState
# ----------------------------------------------


class Module(PBSModule):
    def get_pbs_command(self, options, pbs_script_filename):
        # total parallel process
        np = self.proc_value
        # processes per node, default value 1
        ppn = options.get('ppn', 1)
        # number of physical machines to be reserved
        nodes = float(np) / ppn
        if int(nodes) != nodes:
            Printer.wrn('Warning: NP is not divisible by PPN')
        nodes = int(math.ceil(nodes))

        # memory limit
        mem = int(self.test_case.memory_limit * ppn)

        # time_limit
        walltime = datetime.timedelta(seconds=int(self.test_case.time_limit))

        # get queue, if only -q is set, 'default' queue will be set
        # otherwise given string value will be used
        queue = options.get('queue', True)
        queue = 'default' if type(queue) is not str else queue

        # command
        command = [
            'qsub',
            '-l', 'nodes={nodes}:ppn={ppn}'.format(**locals()), # :nfs4 option may be set
            '-l', 'mem={mem}mb'.format(**locals()),
            '-l', 'walltime={walltime}'.format(**locals()),
            '-l', 'place=infiniband',
            '-q', '{queue}'.format(**locals()),
            '-o', self.pbs_output,
            pbs_script_filename
        ]

        return command


class ModuleJob(Job):
    username = None
    instances = list()

    def __init__(self, job_id, case):
        super(ModuleJob, self).__init__(job_id, case)
        self.parser = self.parser_builder(
            self, 9, JobState.UNKNOWN,
            queue=2,
            name=3,
        )

    @classmethod
    def update_command(cls):
        cls.username = cls.username or getpass.getuser()
        return ['qstat', '-u', cls.username]

    @classmethod
    def create(cls, command_output, case):
        return ModuleJob(re.findall(r'(\d+)', command_output)[0], case)


template = """
#!/bin/bash
#
# Specific PBS setting
#
#PBS -S /bin/bash
#PBS -N flow123d
#PBS -j oe

# load modules
#################
# $$modules$$
#################
module purge
module add /software/modules/current/metabase
module add svn-1.7.6
module add intelcdk-12
module add boost-1.49
module add mpich-p4-intel
module add cmake-2.8
module add python-2.6.2
module unload mpiexec-0.84
module unload mpich-p4-intel
module add openmpi-1.6-intel
module add gcc-4.7.0
module add perl-5.10.1

module add python26-modules-gcc
module add numpy-py2.6
module add python-2.7.6-gcc
module add
#################


# header - info about task
uname -a
echo JOB START: `date`
pwd

echo "$$command$$" > "$$output$$"
$$command$$ >> "$$output$$" 2>&1

echo "$$status_ok$$" >> "$$output$$"
""".lstrip()