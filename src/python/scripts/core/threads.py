#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import threading
import subprocess
import sys
import time
import math
# ----------------------------------------------
from scripts import psutils
from scripts.core import monitors
from scripts.core.base import Printer, Paths, Command, DynamicSleep
from utils.counter import ProgressCounter
from utils.events import Event
from utils.globals import ensure_iterable, wait_for
# ----------------------------------------------


class ExtendedThread(threading.Thread):
    def __init__(self, name, target=None):
        super(ExtendedThread, self).__init__(name=name)
        self.started = None
        self.target = target

        self._is_over = True
        self._returncode = None

        # create event objects
        self.on_start = Event()
        self.on_complete = Event()
        self.on_update = Event()

    def _run(self):
        if self.target:
            self.target()

    @property
    def returncode(self):
        return self._returncode

    @returncode.setter
    def returncode(self, value):
        self._returncode = value

    def start(self):
        self._is_over = False
        super(ExtendedThread, self).start()

    def run(self):
        self._is_over = False
        self.on_start(self)
        self._run()
        self._is_over = True
        self.on_complete(self)

    def is_over(self):
        return self._is_over

    def is_running(self):
        return not self._is_over

    def __repr__(self):
        return "<{self.__class__.__name__}:{running} E:{self.returncode}>".format(
            self=self, running="RUNNING" if self.is_alive() else "STOPPED")


class BinExecutor(ExtendedThread):
    """
    :type process: scripts.psutils.Process
    :type threads: list[scripts.core.threads.BinExecutor]
    """
    threads = list()
    stopped = False

    @staticmethod
    def register_sigint():
        import signal
        signal.signal(signal.SIGINT, BinExecutor.signal_handler)

    @staticmethod
    def signal_handler(signal, frame):
        BinExecutor.stopped = True

        if signal:
            sys.stderr.write("\nError: Caught SIGINT! Terminating application in peaceful manner...\n")
        else:
            sys.stderr.write("\nError: Terminating application threads\n")
        # try to kill all running processes
        for executor in BinExecutor.threads:
            try:
                if executor.process.is_running():
                    sys.stderr.write('\nTerminating process {}...\n'.format(executor.process.pid))
                    executor.process.secure_kill()
            except Exception as e:
                pass
        sys.exit(1)

    def __init__(self, command, name='exec-thread'):
        super(BinExecutor, self).__init__(name)
        BinExecutor.threads.append(self)
        self.command = [str(x) for x in ensure_iterable(command)]
        self.process = None
        self.broken = False
        self.stdout = subprocess.PIPE
        self.stderr = subprocess.PIPE

    def _run(self):
        if self.stopped:
            process = BrokenProcess(Exception('Application terminating'))
            self.returncode = process.returncode
            self.broken = True
            self.process = process
            return

        # run command and block current thread
        try:
            self.process = psutils.Process.popen(self.command, stdout=self.stdout, stderr=self.stderr)
        except Exception as e:
            # broken process
            process = BrokenProcess(e)
            self.returncode = process.returncode
            self.broken = True
            self.process = process
            return

        # process successfully started to wait for result
        # call wait on Popen process
        self.broken = False
        self.process.process.wait()
        self.returncode = getattr(self.process, 'returncode', None)


class BrokenProcess(object):
    def __init__(self, exception=None):
        self.exception = exception
        self.pid = -1
        self.returncode = 666

    @staticmethod
    def is_running():
        return False


class MultiThreads(ExtendedThread):
    """
    :type threads: list[scripts.core.threads.ExtendedThread]
    """
    def __init__(self, name, progress=False):
        super(MultiThreads, self).__init__(name)
        self.threads = list()
        self.returncodes = dict()
        self.running = 0
        self.stop_on_error = False
        self.counter = None
        self.progress = progress
        self.counter = None
        self.index = 0
        self.stopped = False
        self.separate = False

    def run_next(self):
        if self.stopped:
            return False
        if self.index >= self.total:
            return False

        self.index += 1

        if self.counter:
            if self.separate:
                Printer.separator()
            self.counter.next(locals())

        self.threads[self.index - 1].start()
        return True

    def add(self, thread):
        """
        :type thread: scripts.core.threads.ExtendedThread
        """
        self.threads.append(thread)
        self.returncodes[thread] = None
        thread.on_start += self.on_thread_start
        thread.on_complete += self.on_thread_complete

    @property
    def current_thread(self):
        return self.threads[self.index -1]

    @property
    def returncode(self):
        return max(self.returncodes.values()) if self.returncodes else None

    @property
    def total(self):
        return len(self.threads)

    def on_thread_start(self, thread):
        """
        :type thread: scripts.core.threads.ExtendedThread
        """
        self.running += 1

    def on_thread_complete(self, thread):
        self.returncodes[thread] = thread.returncode
        self.running -= 1

    # aliases
    __len__ = total
    __iadd__ = add
    append = add


class SequentialThreads(MultiThreads):
    def __init__(self, name, progress=True, indent=False):
        super(SequentialThreads, self).__init__(name, progress)
        self.thread_name_property = False
        self.indent = indent

    def _run(self):
        if self.progress:
            if self.thread_name_property:
                self.counter = ProgressCounter('{self.name}: {:02d} of {self.total:02d} | {self.current_thread.name}')
            else:
                self.counter = ProgressCounter('{self.name}: {:02d} of {self.total:02d}')

        if self.indent:
            Printer.open()

        while True:
            if not self.run_next():
                break
            self.current_thread.join()

            if self.stop_on_error and self.current_thread.returncode > 0:
                self.stopped = True
                break

        if self.indent:
            Printer.close()


class ParallelThreads(MultiThreads):
    def __init__(self, n=4, name='runner', progress=True):
        super(ParallelThreads, self).__init__(name, progress)
        self.n = n if type(n) is int else 1
        self.counter = ProgressCounter('Case {:02d} of {self.total:02d}')
        self.stop_on_error = True
        self.separate = True

    def ensure_run_count(self):
        if self.stopped:
            return False
        if self.index >= self.total:
            return False

        if self.running < self.n:
            self.run_next()
        return True

    def _run(self):
        # determine how many parallel batches will be
        # later on take cpu into account
        steps = int(math.ceil(float(self.total)/self.n))

        # run each batch
        for i in range(steps):
            batch = [x for x in range(i*self.n, i*self.n+self.n) if x < len(self.threads)]
            for t in batch:
                self.run_next()

            for t in batch:
                self.threads[t].join()

            if self.stop_on_error and self.returncode > 0:
                self.stopped = True
                # no need to stop processes, just break
                break


class PyPy(ExtendedThread):
    """
    :type executor : scripts.core.threads.BinExecutor
    :type case     : scripts.core.prescriptions.TestPrescription
    """

    returncode_map = {
        '0': 'SUCCESS',
        '1': 'ERROR',
        'None': 'SKIPPED',
        '-1': 'SKIPPED',
    }

    def __init__(self, executor, progress=False, period=.5):
        super(PyPy, self).__init__(name='pypy')
        self.executor = executor
        self.period = period
        self.case = None
        self._progress = None

        self.on_process_start = Event()
        self.on_process_complete = Event()
        self.on_process_update = Event()

        # register monitors
        self.info_monitor = monitors.InfoMonitor(self)
        self.limit_monitor = monitors.LimitMonitor(self)
        self.progress_monitor = monitors.ProgressMonitor(self)
        self.error_monitor = monitors.ErrorMonitor(self)

        self.log = False
        self.output_file = None
        self.output_fp = None

        # default value is to "hide" all output
        self.stdout_stderr = subprocess.PIPE

        # different settings in batch mode
        self.progress = progress

        # dynamic sleeper
        self.sleeper = DynamicSleep()

    @property
    def progress(self):
        return self._progress

    @progress.setter
    def progress(self, value):
        self._progress = value
        self.progress_monitor.active = value
        self.stdout_stderr = subprocess.PIPE if value else Paths.temp_file()

    def _run(self):
        fp = self.stdout_stderr
        if self.stdout_stderr not in (subprocess.PIPE, None):
            self.output_file = self.stdout_stderr

            Paths.ensure_path(self.output_file)
            self.output_fp = open(self.output_file, 'w')
            fp = self.output_fp

        # set outputs
        self.executor.stdout = fp
        self.executor.stderr = subprocess.STDOUT

        # start executor
        self.executor.start()
        wait_for(self.executor, 'process')

        if self.executor.broken:
            Printer.err('Could not start command {}: {}',
                             Command.to_string(self.executor.command),
                             getattr(self.executor, 'exception', 'Unknown error'))
            self.returncode = self.executor.returncode

        # if process is not broken, propagate start event
        self.on_process_start(self)

        while self.executor.process.is_running():
            self.on_process_update(self)
            self.sleeper.sleep()

        # get return code
        self.returncode = getattr(self.executor, 'returncode', None)

        # propagate on_complete event
        self.on_process_complete(self)
