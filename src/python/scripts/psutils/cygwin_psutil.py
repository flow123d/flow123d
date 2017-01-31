#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import os
import signal

import time
import subprocess
import csv


from utils.dotdict import Map

KB = 1000.0
MB = KB ** 2
GB = KB ** 3

KiB = 1024.0
MiB = KiB ** 2
GiB = KiB ** 3

# duration used in process kill
_reasonable_amount_of_time = 1
_just_a_sec = 0.1


class use_cache(object):
    """
    Decorator which uses cache for certain amount of time
    """

    def __init__(self, cache_duration=5):
        self.cache_duration = cache_duration

    def __call__(self, f):
        name = f.func_name

        # call wrapped function or use caches value
        def wrapper(other, *args, **kwargs):
            last_update = set_if_not_exists(other, name + '_update', 0)
            last_result = set_if_not_exists(other, name + '_result', None)

            # was is 5 sec since last update?
            since = time.time() - last_update
            if last_result is not None and since < self.cache_duration:
                return last_result

            last_result = f(other, *args, **kwargs)
            setattr(other, name + '_update', time.time())
            setattr(other, name + '_result', last_result)
            return last_result

        return wrapper


def set_if_not_exists(obj, prop, default):
    """
    creates property on object if it does not already exists
    :param obj:
    :param prop:
    :param default:
    :return:
    """
    if hasattr(obj, prop):
        return getattr(obj, prop, default)
    setattr(obj, prop, default)
    return default


def parse_csv(output):
    """
    Parses csv output into list of list of str
    :type output: str
    """
    output = output.strip()
    return [x for x in csv.reader(output.splitlines()) if x]


class Process(object):
    """
    Implementation of Process under Cygwin
    :type process : subprocess.Popen
    """
    platform = 'cygwin'
    memory_info_command = 'wmic process where ProcessId="{}" get WorkingSetSize /format:csv'
    children_command = 'wmic process where ParentProcessId="{}" get ProcessId /format:csv'
    is_running_command = 'wmic process where ProcessId="{}" get ProcessId /format:csv'

    @classmethod
    def popen(cls, *args, **kwargs):
        process = subprocess.Popen(*args, **kwargs)
        return Process(process.pid, process)

    def __init__(self, pid, process=None):
        self.pid = pid
        self.process = process
        self.active = False
        self.start_time = time.time()

    @use_cache()
    def memory_usage(self, prop='vms', units=MiB):
        return sum([x._memory_usage(prop, units) for x in self.children])

    def runtime(self):
        return time.time() - self.start_time

    @use_cache()
    def _memory_usage(self, prop='vms', units=MiB):
        # call wmic
        output = subprocess.check_output(self.memory_info_command.format(self.pid), shell=True)
        # get second line second column value
        csv_output = parse_csv(output)
        return 0.0 if len(csv_output) != 2 else int(csv_output[1][1])

    @use_cache()
    def children(self, recursive=True):
        # call wmic
        output = subprocess.check_output(self.children_command.format(self.pid), shell=True)

        # parse csv output and cut first line (header info)
        children = list()
        csv_output = parse_csv(output)[1:]
        for node, pid in csv_output:
            children.extend(Process(pid).children())
        children.append(self)
        return children

    def terminate(self):
        try:
            if self.process:
                return self.process.terminate()

            os.kill(self.pid, getattr(signal, 'CTRL_C_EVENT', signal.SIGTERM))
        except OSError as e:
            pass

    def kill(self):
        try:
            if self.process:
                return self.process.kill()

            os.kill(self.pid, getattr(signal, 'CTRL_C_EVENT', getattr(signal, 'SIGKILL', signal.SIGTERM)))
        except OSError as e:
            pass

    @use_cache(cache_duration=1.0)
    def is_running(self):
        # call wmic
        output = subprocess.check_output(self.is_running_command.format(self.pid, shell=True))
        csv_output = parse_csv(output)
        return len(csv_output) == 2

    def wait(self, timeout=None):
        return self.process.wait(timeout)

    def secure_kill(self):
        # first, lets be reasonable and terminate all processes (SIGTERM)
        children = self.children()
        self.apply(children, 'terminate')

        # wait jus a sec to let it all blend in
        time.sleep(_just_a_sec)

        # if some process are still running
        if True in self.apply(children, 'is_running'):
            # wait a bit to they can safely exit...
            time.sleep(_reasonable_amount_of_time)
            # check status again, maybe perhaps they finish what
            # they have started and are done?
            if True in self.apply(children, 'is_running'):
                # looks like they are still running to SIGKILL 'em
                self.apply(children, 'kill')
            else:
                # all processes finish they job in reasonable period since
                # SIGTERM was sent
                return True
        else:
            # all processes finish right after* SIGTERM was announced
            return True

        # return max value either True or False
        # False meaning some process is still running

        return not max(self.apply(children, 'is_running'))

    @classmethod
    def apply(cls, children, prop, *args, **kwargs):
        result = []
        for child in children:
            try:
                result.append(getattr(child, prop)(*args, **kwargs))
            except Exception as e:
                pass
        return result

    def __repr__(self):
        if self.process:
            return '<Process ({self.pid}), {self.process}>'.format(self=self)
        return '<Process ({self.pid})>'.format(self=self)
