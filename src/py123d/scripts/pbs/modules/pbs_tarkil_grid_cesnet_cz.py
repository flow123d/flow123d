#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import math
import re
import datetime
import getpass
# ----------------------------------------------
from py123d.scripts.core.base import Printer
from py123d.scripts.pbs.job import Job, JobState
from py123d.scripts.prescriptions.remote_run import PBSModule
from py123d.scripts.prescriptions.remote_run import PBSLimit as L
# ----------------------------------------------


def ceil(n):
    return int(math.ceil(n))


class Module(PBSModule):
    """
    Class Module for TUL hydra cluster PBS system
    """

    def get_pbs_command(self, pbs_script_filename):

        # we try to use single node if number of processes is below 4
        per_node = 4
        requested_cpus = self.case.proc
        requested_memory_per_cpu = self.case.memory_limit
        requested_walltime = datetime.timedelta(seconds=int(self.case.time_limit))
        requested_queue = self.queue

        # prepare actual resources
        total_nodes = ceil(requested_cpus / per_node)
        cpus_per_node = ceil(requested_cpus / total_nodes)
        total_cpus = total_nodes * cpus_per_node

        # # this method will allocate precise amount of mem in total, if nodes and cpus are not dividable
        # # memory per single process will be less but spread across more cpus
        # # example: we want 5 cpu each with 600mb mem
        # # there will be 2 nodes each having 3 cpu, each cpu having 500mb (1500mb per node)
        # # in total we will allocate not a single byte of memory extra (5 * 600mb) = 3000mb
        # memory_per_node = ceil((requested_cpus * requested_memory_per_cpu) / total_nodes)

        # this method will allocate precise amount of mem for each cpu
        # upon non dividable cpus, greater amount in total will be allocated
        # example: we want 5 cpu each with 600mb mem
        # there will be 2 nodes each having 3 cpu, each cpu having 600mb (1800mb per node)
        # in total we will allocate more memory (6 * 600mb) = 3600mb
        memory_per_node = ceil(cpus_per_node * requested_memory_per_cpu)

        if total_cpus != requested_cpus:
            Printer.all.wrn('Requsted {} cpus but will request {} cpus ({}n x {}cpu)',
                            requested_cpus, total_cpus, total_nodes, cpus_per_node)

        self.limits = list()
        self.limits.extend([
            L('select={total_nodes}:ncpus={cpus_per_node}:mem={memory_per_node}mb'),
            L('walltime={requested_walltime}'),
            L(self.case.fs.pbs_output, type='-o')
        ])
        if self.queue:
            self.limits.append(
                L('{requested_queue}', type='-q')
            )

        command = ['qsub']
        for limit in self.limits:
            limit.value = limit.value.format(**locals())
            command.extend(limit())
        command.append(pbs_script_filename)

        return command


class ModuleJob(Job):
    """
    Class ModuleJob  for TUL hydra cluster PBS system job
    """

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
        return ['qstat', '-xu', cls.username]

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
module purge
module load /software/modules/current/metabase
module load cmake-2.8.12
module load gcc-4.9.2
module load boost-1.56-gcc
module load perl-5.20.1-gcc
module load openmpi

module load python34-modules-gcc
module load python-3.4.1-gcc

module unload gcc-4.8.1
module unload openmpi-1.8.2-gcc
#################


# header - info about task
echo '---------------------------------'
lscpu
uname -a
echo JOB START: `date`
cat $PBS_NODEFILE
echo '---------------------------------'
pwd

$$command$$

""".lstrip()
#PBS qsub -l select=2:ncpus=3 -l mem=1gb -I
#PBS -l select=2:ncpus=4:mem=8000mb:cl_manegrot=True