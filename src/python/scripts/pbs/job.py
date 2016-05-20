#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import subprocess
import time
import datetime


class JobState(object):
    COMPLETED = 'C'
    EXITING = 'E'
    HELD = 'H'
    QUEUED = 'Q'
    RUNNING = 'R'
    TRANSFER = 'T'
    WAITING = 'W'
    SUSPENDED = 'S'
    UNKNOWN = 'U'
    EXIT_OK = 'K'
    EXIT_ERROR = 'X'
    _map = {}

    def __init__(self, value='U'):
        if not JobState._map:
            JobState._map = {getattr(JobState, v): v for v in dir(JobState) if v.upper() == v}
        self.value = str(value).upper()

    def __eq__(self, other):
        if type(other) is JobState:
            return self.value == other.value
        return self.value == str(other).upper()

    def __ne__(self, other):
        if type(other) is JobState:
            return self.value != other.value
        return self.value != str(other).upper()

    def __bool__(self):
        return self.value == self.COMPLETED
    __nonzero__=__bool__

    def __hash__(self):
        return hash(self.value)

    def __repr__(self):
        return "{enum}".format(cls=self.__class__.__name__, enum=JobState._map.get(self.value))

    def enum(self):
        return self.value


class Job(object):
    """
    :type case : scripts.core.prescriptions.PBSModule
    """
    def __init__(self, job_id, case):
        self.id = job_id
        self.case = case

        self.name = None
        self.queue = None
        self.status_changed = False
        self.parser = lambda x: None
        self.is_active = True

        self.last_status = JobState(JobState.UNKNOWN)
        self._status = JobState(JobState.UNKNOWN)

    @property
    def status(self):
        """
        :rtype : scripts.pbs.job.JobState
        """
        return self._status
    
    @status.setter
    def status(self, value):
        """
        :type value: scripts.pbs.job.JobState or str
        """
        if type(value) is not JobState:
            value = JobState(value)

        # set state
        self.last_status = self._status
        self._status = value
        self.status_changed = self.last_status != self._status

    def update_status(self, output):
        if self.is_active:
            self.status = self.parse_status(output)

    def raise_not_found(self):
        raise Exception('Job with id {self.id} does not exists'.format(self=self))

    def __repr__(self):
        return "<Job #{s.id}, status {s.status} in queue {s.queue}>".format(s=self)

    # --------------------------------------------------------

    @classmethod
    def update_command(cls):
        """
        :rtype : list[str]
        """
        raise NotImplementedError('Method not implemented!')

    def parse_status(self, output=""):
        """
        :rtype : str
        """
        return self.parser(output)

    @classmethod
    def create(cls, output, case):
        """
        Creates instance of Job from qsub output
        :param output: output of qsub command
        """
        raise NotImplementedError('Method not implemented!')

    @classmethod
    def parser_builder(cls, o, index, default=JobState.UNKNOWN, **kwargs):
        def parse(output):
            for line in output.splitlines():
                if line.find(o.id) != -1:
                    info = line.split()
                    for name, i in kwargs.items():
                        setattr(o, name, info[i])
                    return info[index]
            return default
        return parse


class MultiJob(object):
    """
    :type items : list[scripts.pbs.job.Job]
    :type cls   : class
    """
    def __init__(self, cls):
        self.items = list()
        self.cls = cls
        self._iter_index = 0
        self.start_time = None

    def __iter__(self):
        self.iter_index = 0
        return self

    def next(self):
        """
        :rtype : scripts.pbs.job.Job
        """
        if self._iter_index >= len(self.items):
            raise StopIteration
        else:
            self._iter_index += 1
            return self.items[self._iter_index - 1]

    __next__ = next

    def add(self, *items):
        self.items.extend(items)

    def status(self):
        return {item: item.status for item in self.items}

    def update(self):
        self.start_time = self.start_time or time.time()
        output = subprocess.check_output(self.cls.update_command())
        return [item.update_status(output) for item in self.items]

    def is_running(self):
        status = set(self.status().values())
        status = status - {JobState.EXIT_OK, JobState.EXIT_ERROR}
        return bool(status)

    def print_status(self, printer):
        """
        :type printer: scripts.core.base.Printer
        """
        for item in self.items:
            printer.key(str(item))

    def status_changed(self, desired=JobState.COMPLETED):
        """
        :rtype : list[scripts.pbs.job.Job]
        """
        if desired is None:
            return [item for item in self.items if item.status_changed]
        if type(desired) is not set:
            desired = set(desired)

        return [item for item in self.items if item.status_changed and item.status in desired]


    def get_status_line(self):
        status = self.status().values()
        result = dict()

        for s in set(status):
            result[s] = status.count(s)

        return 'Time elapsed: {delta} | {status}'.format(
            delta=datetime.timedelta(seconds=int(time.time() - self.start_time)),
            status=', '.join(['{}: {:d}'.format(k, v) for k, v in result.items()])
        )
