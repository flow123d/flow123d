#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import subprocess

import time
from scripts.core.base import Printer, Command, Paths, IO
from utils.counter import ProgressTime


def ensure_active(f):
    def wrapper(self, *args, **kwargs):
        if self.active:
            return f(self, *args, **kwargs)
    return wrapper


class ThreadMonitor(object):
    def __init__(self, pypy):
        """
        :type pypy: scripts.core.threads.PyPy
        """
        self.pypy = pypy
        self.active = True
        self.printer = Printer(Printer.LEVEL_KEY)

        # add listeners
        self.pypy.on_process_start += self.on_start
        self.pypy.on_process_update += self.on_update
        self.pypy.on_process_complete += self.on_complete

    def deactivate(self):
        self.active = False

    def activate(self):
        self.active = True

    @ensure_active
    def on_start(self, pypy=None):
        pass

    @ensure_active
    def on_update(self, pypy=None):
        pass

    @ensure_active
    def on_complete(self, pypy=None):
        pass


class ProgressMonitor(ThreadMonitor):
    def __init__(self, pypy):
        super(ProgressMonitor, self).__init__(pypy)
        self.timer = ProgressTime('Running | elapsed time {}')
        # increase priority
        self.pypy.on_process_complete.set_priority(self.on_complete, 5)

    @ensure_active
    def on_update(self, pypy=None):
        self.timer.update()

    @ensure_active
    def on_complete(self, pypy=None):
        self.timer.stop()


class InfoMonitor(ThreadMonitor):
    def __init__(self, pypy):
        super(InfoMonitor, self).__init__(pypy)
        self.command_str = Command.to_string(self.pypy.executor.command)
        self.start_fmt = 'Executing {self.command_str}'
        self.end_fmt = 'Command ({self.pypy.executor.process.pid}) ended with {self.pypy.returncode}'

    @ensure_active
    def on_start(self, pypy=None):
        if self.start_fmt:
            self.printer.key(self.start_fmt.format(**dict(self=self)))

    @ensure_active
    def on_complete(self, pypy=None):

        # print either error that command failed or on_complete info id exists
        if self.pypy.returncode > 0:
            self.printer.err('Error! Command ({process.pid}) ended with {process.returncode}'.
                             format(process=self.pypy.executor.process))
            self.printer.err(Command.to_string(self.pypy.executor.command))
        elif self.end_fmt:
            self.printer.key(self.end_fmt.format(**dict(self=self)))

        if not self.pypy.progress:
            self.printer.line()
            output = IO.read(self.pypy.output_file)
            if output:
                self.printer.out(output)


class Limits(object):
    def __init__(self, time_limit=None, memory_limit=None):
        self.time_limit = time_limit
        self.memory_limit = memory_limit


class LimitMonitor(ThreadMonitor):
    """
    :type process: scripts.psutils.Process
    """
    def __init__(self, pypy):
        super(LimitMonitor, self).__init__(pypy)
        self.process = None
        self.memory_limit = None
        self.time_limit = None
        self.monitor_thread = None
        self.terminated = False
        self.terminated_cause = None

    def set_limits(self, case):
        """
        :type case: scripts.config.yaml_config.YamlConfigCase
        """
        # empty Limits object
        if not case:
            self.memory_limit = None
            self.time_limit = None
            return

        self.memory_limit = case.memory_limit
        self.time_limit = case.time_limit

    @ensure_active
    def on_start(self, pypy=None):
        self.process = self.pypy.executor.process

    @ensure_active
    def on_update(self, pypy=None):
        if self.terminated:
            return

        if self.time_limit:
            try:
                runtime = self.process.runtime()
                if runtime > self.time_limit:
                    self.printer.err(
                        'Error: Time limit exceeded! {:1.2f}s of runtime, {:1.2f}s allowed'.format(
                            runtime, self.time_limit
                        )
                    )
                    self.terminated_cause = 'TIME_LIMIT'
                    self.terminated = True
                    self.process.secure_kill()
            # except NoSuchProcess as e1:
            #     pass
            except AttributeError as e2:
                pass

        if self.memory_limit:
            try:
                memory_usage = self.process.memory_usage()
                if memory_usage > self.memory_limit:
                    self.printer.err('Error: Memory limit exceeded! {:1.2f}MB used, {:1.2f}MB allowed'.format(
                        memory_usage, self.memory_limit
                        )
                    )
                    self.terminated_cause = 'MEMORY_LIMIT'
                    self.terminated = True
                    self.process.secure_kill()
            # except NoSuchProcess as e1:
            #     pass
            except AttributeError as e2:
                pass


class ErrorMonitor(ThreadMonitor):
    def __init__(self, pypy):
        super(ErrorMonitor, self).__init__(pypy)
        self.message = 'Command failed'
        self.tail = 20

    @ensure_active
    def on_complete(self, pypy=None):
        if self.pypy.returncode > 0:
            if self.message:
                self.printer.line()
                self.printer.err(self.message)

            # if file pointer exist try to read errors and outputs
            if self.pypy.output_file:
                lines = (IO.read(self.pypy.output_file)).splitlines()[-self.tail:]

                if lines:
                    self.printer.err("## Command's last 20 lines (rest in {})".format(
                        Paths.abspath(self.pypy.output_file)))
                    self.printer.err('#' * 60)
                    for l in lines:
                        self.printer.err('## ' + str(l))
                    self.printer.err('#' * 60)
                else:
                    self.printer.err('#' * 60)
                    self.printer.err("## Both stdout and stderr are empty!")
                    self.printer.err("## Could not extract any information from in {}".format(
                        Paths.abspath(self.pypy.output_file)))
                    self.printer.err('#' * 60)
