#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import subprocess
# ----------------------------------------------
from scripts.core.base import Printer, DynamicSleep
from scripts.core.threads import ResultHolder
from scripts.pbs.job import MultiJob, JobState
# ----------------------------------------------


def insert_jobs(jobs, pbs_module):
    # start jobs
    total = len(jobs)
    Printer.all.out('Starting {} job/s', total)

    job_id = 0
    multijob = MultiJob(pbs_module.ModuleJob)
    with Printer.all.with_level():
        for qsub_command, pbs_run in jobs:
            job_id += 1

            Printer.console.dyn('Starting jobs {:02d} of {:02d}', job_id, total)
            Printer.batched.out('Starting jobs {:02d} of {:02d}', job_id, total)

            output = subprocess.check_output(qsub_command)
            job = pbs_module.ModuleJob.create(output, pbs_run.case)
            job.full_name = "Case {}".format(pbs_run.case)
            multijob.add(job)

    # inform that job were inserted
    Printer.console.newline()
    Printer.all.out('{} job/s inserted into queue', total)

    return multijob


def finish_pbs(multijob, finish_method):
    # first update to get more info about multijob jobs
    Printer.all.sep()
    Printer.all.out('Determining job status')
    multijob.update()

    # print jobs statuses
    with Printer.all.with_level():
        multijob.print_status()

    Printer.console.sep()
    Printer.console.dyn(multijob.get_status_line())
    result = ResultHolder()

    # use dynamic sleeper
    sleeper = DynamicSleep(min=300, max=5000, steps=5)

    # wait for finish
    while multijob.is_running():
        Printer.console.dyn('Updating job status')
        multijob.update()
        Printer.console.dyn(multijob.get_status_line())

        # if some jobs changed status add new line to dynamic output remains
        jobs_changed = multijob.get_all(status=JobState.COMPLETED)
        if jobs_changed:
            Printer.console.newline()
            Printer.all.sep()

        # get all jobs where was status update to COMPLETE state
        for job in jobs_changed:
            pypy = finish_method(job, not Printer.batched.is_muted())
            if pypy:
                result.add(pypy)

        if jobs_changed:
            Printer.console.newline()

        # after printing update status lets sleep for a bit
        if multijob.is_running():
            sleeper.sleep()

    Printer.all.sep()
    # print final result
    Printer.all.out(multijob.get_status_line())
    Printer.all.out('All jobs finished')

    return result, multijob


def do_work(jobs, pbs_module, finish_method):
    """
    :rtype: (ResultHolder, MultiJob)
    """
    multijob = insert_jobs(jobs, pbs_module)
    result, multijob = finish_pbs(multijob, finish_method)
    return result, multijob