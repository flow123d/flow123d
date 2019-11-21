#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import psutil
import time
# ----------------------------------------------
from psutil import NoSuchProcess


class try_catch(object):
    """
    Decorator which uses cache for certain amount of time
    """

    def __init__(self, default=0):
        self.default = default

    def __call__(self, f):
        # call wrapped function or use caches value
        def wrapper(other, *args, **kwargs):
            try:
                return f(other, *args, **kwargs)
            except NoSuchProcess:
                return self.default

        return wrapper


KB = 1000.0
MB = KB ** 2
GB = KB ** 3

KiB = 1024.0
MiB = KiB ** 2
GiB = KiB ** 3

# duration used in process kill
_reasonable_amount_of_time = 1
_just_a_sec = 0.1


class Process(psutil.Process):
    """
    Implementation of Process under windows and linux
    """

    platform = 'windows'

    @classmethod
    def popen(cls, *args, **kwargs):
        process = psutil.Popen(*args, **kwargs)
        return Process(process.pid, process)

    # @classmethod
    # def _popen(cls, *args, **kwargs):
    #     process = psutil.Popen(*args, **kwargs)
    #     return Process(process.pid)

    def __init__(self, pid, process=None):
        """
        :type process: psutil.Popen
        """
        super(Process, self).__init__(pid)
        self.process = process
        self.terminated = False

    def wait(self, timeout=None):
        if self.process:
            return self.process.wait()
        return super(Process, self).wait(timeout)

    @property
    def returncode(self):
        if self.terminated:
            return 1
        if self.process:
            return self.process.returncode

    def children(self, recursive=True):
        children = super(Process, self).children(recursive)
        children = [Process(x.pid) for x in children]
        children.append(self)
        return children

    @try_catch(default=0)
    def memory_usage(self, prop='rss', units=MiB):
        """
        rss: aka “Resident Set Size”, this is the non-swapped physical memory a process has used. On UNIX it matches “top“‘s RES column). On Windows this is an alias for wset field and it matches “Mem Usage” column of taskmgr.exe.
        vms: aka “Virtual Memory Size”, this is the total amount of virtual memory used by the process. On UNIX it matches “top“‘s VIRT column. On Windows this is an alias for pagefile field and it matches “Mem Usage” “VM Size” column of taskmgr.exe.
        shared: (Linux) memory that could be potentially shared with other processes. This matches “top“‘s SHR column).
        text (Linux, BSD): aka TRS (text resident set) the amount of memory devoted to executable code. This matches “top“‘s CODE column).
        data (Linux, BSD): aka DRS (data resident set) the amount of physical memory devoted to other than executable code. It matches “top“‘s DATA column).
        lib (Linux): the memory used by shared libraries.
        dirty (Linux): the number of dirty pages.
        pfaults (macOS): number of page faults.
        pageins (macOS): number of actual pageins.
        """
        total = getattr(self.memory_info(), prop)
        for child in super(Process, self).children(recursive=True):
            total += getattr(child.memory_info(), prop)
        return total / MiB

    @try_catch(default=0)
    def runtime(self):
        return time.time() - self.create_time()

    @try_catch(default=True)
    def secure_kill(self):
        # first, lets be reasonable and terminate all processes (SIGTERM)
        self.terminated = True
        children = self.children()
        self.apply(children, 'terminate')

        # wait just a sec to let it all blend in
        time.sleep(_just_a_sec)

        # if some process are still running
        if True in self.apply(children, 'is_running'):
            # wait a bit to they can safely exit...
            time.sleep(_reasonable_amount_of_time)
            # check status again, maybe perhaps they finish what
            # they have started and are done?
            if True in self.apply(children, 'is_running'):
                # looks like they are still running so SIGKILL 'em
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
            except psutil.NoSuchProcess as e:
                pass
        return result
