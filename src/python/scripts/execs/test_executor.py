#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from __future__ import absolute_import

from subprocess import PIPE
import threading
import time

import sys

import psutil, shutil
from psutil import NoSuchProcess

from scripts.core.base import Paths, PathFilters, Printer
from utils.counter import ProgressCounter
from utils.globals import ensure_iterable


class ProcessUtils(object):
    KB = 1000.0
    MB = KB ** 2
    GB = KB ** 3

    KiB = 1024.0
    MiB = KiB ** 2
    GiB = KiB ** 3

    _reasonable_amount_of_time = 1
    _just_a_sec = 0.1


    @classmethod
    def list_children(cls, process):
        """
        :type process: psutil.Process
        """
        result = []
        for child in process.children():
            result.extend(cls.list_children(child))
        result.append(process)
        return result

    @classmethod
    def get_memory_info(cls, process, prop='vms', units=MiB):
        """
        :type process: psutil.Process
        """
        children = cls.list_children(process)
        total = 0
        for child in children:
            try:
                total += getattr(child.memory_info(), prop)
            except NoSuchProcess:
                continue
        return total / units

    @classmethod
    def apply(cls, children, prop, *args, **kwargs):
        result = []
        for child in children:
            try:
                result.append(getattr(child, prop)(*args, **kwargs))
            except NoSuchProcess as e: pass
        return result

    @classmethod
    def terminate(cls, process):
        cls.terminate_all(cls.list_children(process))

    @classmethod
    def terminate_all(cls, children):
        cls.apply(children, 'terminate')

    @classmethod
    def kill(cls, process):
        cls.kill_all(cls.list_children(process))

    @classmethod
    def kill_all(cls, children):
        cls.apply(children, 'kill')

    @classmethod
    def secure_kill(cls, process):
        # first, lets be reasonable and terminate all processes (SIGTERM)
        children = cls.list_children(process)
        cls.terminate_all(children)

        # wait jus a sec to let it all blend in
        time.sleep(cls._just_a_sec)

        # if some process are still running
        if True in cls.apply(children, 'is_running'):
            # wait a bit to they can safely exit...
            time.sleep(cls._reasonable_amount_of_time)
            # check status again, maybe perhaps they finish what
            # they have started and are done?
            if True in cls.apply(children, 'is_running'):
                # looks like they are still running to SIGKILL 'em
                cls.kill_all(children)
            else:
                # all processes finish they job in reasonable period since
                # SIGTERM was sent
                return True
        else:
            # all processes finish right after* SIGTERM was announced
            return True

        # return max value either True or False
        # False meaning some process is still running

        return not max(cls.apply(children, 'is_running'))


class ExtendedThread(threading.Thread):
    def __init__(self, name, target=None):
        super(ExtendedThread, self).__init__(name=name, target=target)
        self._is_over = True
        self.name = name
        assert type(name) is str

    def _run(self):
        pass

    def run(self):
        self._is_over = False
        self._run()
        self._is_over = True

    def is_over(self):
        return self._is_over

    def is_running(self):
        return not self._is_over


class BinExecutor(ExtendedThread):
    """
    :type process: psutil.Popen
    :type threads: list[scripts.execs.test_executor.BinExecutor]
    """
    threads = list()

    @staticmethod
    def register_SIGINT():
        import signal
        signal.signal(signal.SIGINT, BinExecutor.signal_handler)

    @staticmethod
    def signal_handler(signal, frame):
        sys.stderr.write("\nError: Caught SIGINT! Terminating application in peaceful manner...\n")
        # try to kill all running processes
        for executor in BinExecutor.threads:
            try:
                if executor.process.is_running():
                    sys.stderr.write('\nTerminating process {}...\n'.format(executor.process.pid))
                    ProcessUtils.secure_kill(executor.process)
            except Exception as e:
                pass
        sys.exit(1)
        raise Exception('You pressed Ctrl+C!')

    def __init__(self, command=list(), name='bin'):
        self.threads.append(self)
        self.command = [str(x) for x in ensure_iterable(command)]
        self.process = None
        self.running = False
        self.stdout = PIPE
        self.stderr = PIPE
        self.returncode = None
        super(BinExecutor, self).__init__(name)

    def _run(self):
        # run command and block current thread
        try:
            self.process = psutil.Popen(self.command, stdout=self.stdout, stderr=self.stderr)
            self.process.wait()
        except Exception as e:
            # broken process
            self.process = BrokenProcess(e)
        self.returncode = getattr(self.process, 'returncode', None)


class ParallelProcesses(ExtendedThread):
    """
    :type threads: list[threading.Thread]
    """
    def __init__(self, name, *args):
        super(ParallelProcesses, self).__init__(name)
        self.threads = list()
        self.threads.extend(args)

    def add(self, thread):
        self.threads.append(thread)

    def _run(self):
        for t in self.threads:
            t.start()

        for t in self.threads:
            t.join()


class SequentialProcesses(ExtendedThread):
    """
    :type threads: list[threading.Thread]
    """
    def __init__(self, name, pbar=True, indent=False, *args):
        super(SequentialProcesses, self).__init__(name)
        self.pbar = pbar
        self.threads = list()
        self.threads.extend(args)
        self.thread_name_property = False
        self.stop_on_error = False
        self.returncode = None
        self.indent = indent

    def add(self, thread):
        self.threads.append(thread)

    def _run(self):
        total = len(self.threads)
        rcs = [None]
        pc = None

        if self.indent:
            Printer.open()

        if self.pbar:
            if self.thread_name_property:
                pc = ProgressCounter('{self.name}: {:02d} of {total:02d} | {t.name}')
            else:
                pc = ProgressCounter('{self.name}: {:02d} of {total:02d}')

        for t in self.threads:
            if self.pbar:
                pc.next(locals())

            t.start()
            t.join()
            rc = getattr(t, 'returncode', None)
            rcs.append(rc)

            if self.stop_on_error and rc > 0:
                # Printer.out('Aborted next operations due to error')
                break

        self.returncode = max(rcs)

        if self.indent:
            Printer.close()


class ParallelRunner(object):
    """
    :type threads: list[scripts.execs.test_executor.ExtendedThread]
    """
    def __init__(self, n=4):
        self.n = n if type(n) is int else 1
        self.i = 0
        self.threads = list()
        self.printer = Printer(Printer.LEVEL_KEY)

    def add(self, process):
        self.threads.append(process)

    def active_count(self):
        total = 0
        for process in self.threads:
            total += 1 if process.is_alive() else 0
        return total

    def complete_count(self):
        total = 0
        for process in self.threads:
            total += 1 if process.is_over() else 0
        return total

    def run(self):
        self.i = 0
        total = len(self.threads)

        pc = ProgressCounter('Case {:02d} of {total:02d}')
        while self.i < total:
            while self.active_count() < self.n:
                self.printer.key('-' * 60)
                pc.next(locals())
                self.threads[self.i].start()
                self.i += 1

                # no more processes
                if self.i >= total:
                    break
                # sleep a bit to other threads can be active again
                time.sleep(0.1)
            time.sleep(0.1)

    def __repr__(self):
        return '<ParallelRunner x {self.n}>'.format(self=self)


class FD(object):
    def __init__(self):
        self.data = ''

    def write(self, data=''):
        self.data = data[:-1] if data.endswith('\r') else data

class BrokenProcess(object):
    def __init__(self, exception=None):
        self.exception = exception
        self.pid = -1
        self.returncode = 666

    def is_running(self):
        return False