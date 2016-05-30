#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from __future__ import absolute_import

import threading, time, sys, psutil, subprocess
from subprocess import PIPE
from psutil import NoSuchProcess
from scripts.core.base import Printer, Paths
from scripts.execs.test_executor import ProcessUtils
from utils.globals import wait_for
from progressbar import ProgressBar, Timer, AnimatedMarker


class PyProcess(threading.Thread):
    """
    :type monitors: list[AbstractMonitor]
    """
    def __init__(self, executor, printer=None, batch_mode=False, period=.1,):
        """
        :type printer: scripts.core.base.Printer
        :type executor: scripts.execs.test_executor.BinExecutor
        """
        super(PyProcess, self).__init__(name='pypy')
        self.executor = executor
        self.period = period
        self.batch_mode = batch_mode
        self.monitors = []
        self.printer = Printer(Printer.LEVEL_DBG) if not printer else printer

        self.limit_monitor = LimitMonitor(self)
        self.info_monitor = InfoMonitor(self)
        self.progress_monitor = ProgressMonitor(self)

        self.add_monitor(self.info_monitor)
        self.add_monitor(self.limit_monitor)
        self.add_monitor(self.progress_monitor)

        if batch_mode:
            self.progress_monitor.active = False

        self.broken_process = False
        self.name = None
        self.returncode = None

    def add_monitor(self, monitor):
        self.monitors.append(monitor)

    def run(self):
        # alter streams
        self.executor.stdout = self.info_monitor.stdout_stderr
        self.executor.stderr = subprocess.STDOUT

        # start executions
        self.executor.start()

        # start all monitors
        try:
            [m.start() for m in self.monitors if m.active]
        except Exception as e: pass

        # and wait for process to be created
        wait_for(self.executor, 'process')
        from scripts.execs.test_executor import BrokenProcess

        if type(self.executor.process) is BrokenProcess:
            self.broken_process = True

            self.info_monitor.start()
            self.printer.err('Error!: Process is broken: {}'.format(self.executor.process.exception))
            self.info_monitor.stop()
            self.returncode = self.executor.process.returncode
            return

        # wait a bit so other threads can start
        self.sleep()
        # wait for the end
        while True:
            if not self.executor.process.is_running():
                break

            try:
                [m.update() for m in self.monitors if m.active]
            except: pass
            self.sleep()
        self.returncode = getattr(self.executor, 'returncode', None)

        # stop all monitors
        try:
            [m.stop() for m in self.monitors if m.active]
        except: pass

    def sleep(self):
        time.sleep(self.period)


class AbstractMonitor(object):
    def __init__(self, pp):
        """
        :type pp: PyProcess
        """
        self.pp = pp
        self.printer = self.pp.printer.copy()
        self.active = True

    def start(self):
        pass

    def stop(self):
        pass

    def update(self):
        pass


class ExtendedTimer(Timer):
    def __init__(self, limit_monitor):
        """
        :type limit_monitor: scripts.execs.monitor.LimitMonitor
        """
        super(ExtendedTimer, self).__init__()
        self.limit_monitor = limit_monitor

    def update(self, pbar):
        runtime = self.limit_monitor.get_runtime() or pbar.seconds_elapsed
        if not pbar.finished:
            return 'Elapsed Time: {}'.format(Timer.format_time(int(runtime)))

        base = Timer.format_time(int(runtime))
        millis = '{:0.3f}'.format(runtime - int(runtime))
        return 'Elapsed Time: {}:{}'.format(base, millis[2:])


class ProgressMonitor(AbstractMonitor):
    def __init__(self, pp):
        """
        :type pp: scripts.execs.monitor.PyProcess
        """
        widgets = ['Running ', AnimatedMarker(), ' (', ExtendedTimer(pp.limit_monitor), ')']
        super(ProgressMonitor, self).__init__(pp)
        self.progress_bar = ProgressBar(widgets=widgets, maxval=100, fd=sys.stdout)
        self.progress_bar.signal_set = False

    def start(self):
        super(ProgressMonitor, self).start()
        self.progress_bar.start()

    def stop(self):
        super(ProgressMonitor, self).stop()
        self.progress_bar.finish()

    def update(self):
        super(ProgressMonitor, self).update()
        self.progress_bar.update(0)


class InfoMonitor(AbstractMonitor):
    def __init__(self, pp):
        super(InfoMonitor, self).__init__(pp)
        self.n_lines = 30
        self._stdout_stderr = PIPE
        self.fp = None
        self.details = dict(
            command_str=' '.join(self.pp.executor.command),
        )
        self.details['self'] = self
        self.start_fmt = 'Executing {command_str}'
        self.end_fmt = 'Command ({self.pp.executor.process.pid}) ended with {self.pp.executor.process.returncode}'

    @property
    def stdout_stderr(self):
        if self._stdout_stderr in (None, PIPE):
            return self._stdout_stderr
        # it is file
        if not self.fp:
            Paths.ensure_path(self._stdout_stderr)
            self.fp = open(self._stdout_stderr, 'a+')
        return self.fp

    @stdout_stderr.setter
    def stdout_stderr(self, value):
        self._stdout_stderr = value

    def start(self):
        if self.start_fmt:
            self.printer.key(self.start_fmt.format(**self.details))

    def stop(self):
        wait_for(self.pp.executor.process, 'returncode')
        if self.end_fmt:
            self.printer.key(self.end_fmt.format(**self.details))

        if self.pp.executor.process.returncode > 0:
            self.printer.err('Error! Command ({process.pid}) ended with {process.returncode}'.format(
                process=self.pp.executor.process))
            self.printer.err('Command string: {}\nRaw Command:    {}'.format(
                ' '.join(self.pp.executor.command),
                self.pp.executor.command
                ))

            # if file pointer exist try to read errors and outputs
            if self.fp:
                with open(self._stdout_stderr, 'r') as read_fp:
                    lines = read_fp.read().splitlines()[-self.n_lines:]
                    if lines:
                        self.printer.err("## Command's last 10 lines (rest in {})".format(
                            Paths.abspath(self._stdout_stderr)))
                        self.printer.err('#' * 60)
                        for l in lines:
                            self.printer.err('## ' + str(l)[:255])
                        self.printer.err('#' * 60)
                    else:
                        self.printer.err('#' * 60)
                        self.printer.err("## Both stdout and stderr are empty!")
                        self.printer.err("## Could not extract any information from in {}".format(
                            Paths.abspath(self._stdout_stderr)))
                        self.printer.err('#' * 60)

        if self.end_fmt:
            self.printer.key('-' * 60)

        # close file pointer is exists
        if self.fp:
            self.fp.close()


class Limits(object):
    def __init__(self, time_limit=None, memory_limit=None):
        self.time_limit = time_limit
        self.memory_limit = memory_limit


class LimitMonitor(AbstractMonitor):
    """
    :type process: psutil.Process
    """
    def __init__(self, pp):
        super(LimitMonitor, self).__init__(pp)
        self.process = None
        self.memory_limit = None
        self.time_limit = None
        self.monitor_thread = None
        self.terminated = False
        self.terminated_cause = None

    def set_limits(self, limits):
        """
        :type limits: Limits
        """
        # empty Limits object
        if not limits:
            limits = Limits()

        self.memory_limit = limits.memory_limit
        self.time_limit = limits.time_limit

    def start(self):
        # self.printer.dbg('Limits: time: {self.time_limit_str} | memory: {self.memory_limit_str}'.format(self=self))
        wait_for(self.pp.executor, 'process')
        wait_for(self.pp.executor.process, 'pid')
        self.process = psutil.Process(self.pp.executor.process.pid)

    def update(self):
        self.check_limits()

    def get_runtime(self):
        try:
            return time.time() - self.process.create_time()
        except psutil.NoSuchProcess as e1:
            # process has ended
            return 0
        except AttributeError as e2:
            # process did not start
            return 0

    def get_memory_usage(self):
        memory = ProcessUtils.get_memory_info(self.process)
        return memory

    def terminate(self):
        ProcessUtils.secure_kill(self.process)

    def check_limits(self):
        if self.terminated:
            return

        if self.time_limit:
            runtime = self.get_runtime()
            if runtime > self.time_limit:
                self.printer.err('Error: Process running longer then expected! {:1.2f}s of runtime, {:1.2f}s allowed'.format(
                    runtime, self.time_limit
                    )
                )
                self.terminate()
                self.terminated_cause = 'TIME_LIMIT'
                self.terminated = True

        if self.memory_limit:
            memory_usage = self.get_memory_usage()
            if memory_usage > self.memory_limit:
                self.printer.err('Error: Memory usage exceeded limit! {:1.2f}MB used, {:1.2f}MB allowed'.format(
                    memory_usage, self.memory_limit
                    )
                )
                self.terminate()
                self.terminated_cause = 'MEMORY_LIMIT'
                self.terminated = True

    @property
    def time_limit_str(self):
        if self.time_limit:
            return '{:1.2f}s'.format(self.time_limit)

    @property
    def memory_limit_str(self):
        if self.memory_limit:
            return '{:1.2f}MB'.format(self.memory_limit)