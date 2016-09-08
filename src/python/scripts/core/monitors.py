#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from scripts.core.base import Printer, Command, Paths
from utils.counter import ProgressTime
# ----------------------------------------------
from utils.strings import format_n_lines
# ----------------------------------------------


def ensure_active(f):
    """
    Wrapper which executes its action if class is active
    :param f:
    :return:
    """
    def wrapper(self, *args, **kwargs):
        if self.active:
            return f(self, *args, **kwargs)
    return wrapper


class ThreadMonitor(object):
    """
    Class ThreadMonitor is abstract class for monitoring other threads.
    Event system is used.
    """

    def __init__(self, pypy):
        """
        :type pypy: scripts.core.threads.PyPy
        """
        self.pypy = pypy
        self.active = True

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
    """
    Class ProgressMonitor monitors PyPy object and reports progress
    """
    
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
        self.timer.format = 'Done    | elapsed time {}'
        self.timer.stop()


class InfoMonitor(ThreadMonitor):
    """
    Class InfoMonitor monitors PyPy object and prints info at the beginning
    and at the end. On error prints details.
    """
    
    def __init__(self, pypy):
        super(InfoMonitor, self).__init__(pypy)
        self.command_str = Command.to_string(self.pypy.executor.command)
        self.start_fmt = 'Executing {self.command_str}'
        self.end_fmt = 'Command ({self.pypy.executor.process.pid}) ended with {self.pypy.returncode}'

    @ensure_active
    def on_start(self, pypy=None):
        if self.start_fmt:
            Printer.out(self.start_fmt.format(**dict(self=self)))

    @ensure_active
    def on_complete(self, pypy=None):

        # print either error that command failed or on_complete info id exists
        if self.pypy.returncode > 0:
            Printer.err('Error! Command ({process.pid}) ended with {process.returncode}'.
                             format(process=self.pypy.executor.process))
            Printer.err(Command.to_string(self.pypy.executor.command))
        elif self.end_fmt:
            Printer.out(self.end_fmt.format(**dict(self=self)))

        if not self.pypy.progress:
            Printer.separator()
            output = self.pypy.executor.output.read()
            if output:
                Printer.out(output)


class LimitMonitor(ThreadMonitor):
    """
    Class LimitMonitor monitors resources of PyPy process and terminates
    PyPy if limits are not withheld.
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
        :type case: scripts.config.yaml_config.ConfigCase
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
                    Printer.out()
                    Printer.err(
                        'Error: Time limit exceeded! {:1.2f}s of runtime, {:1.2f}s allowed'.format(
                            runtime, self.time_limit
                        )
                    )
                    self.terminated_cause = 'TIME_LIMIT'
                    self.terminated = True
                    self.process.secure_kill()
                    return
            except AttributeError as e2:
                pass

        if self.memory_limit:
            try:
                memory_usage = self.process.memory_usage()
                if memory_usage > self.memory_limit:
                    Printer.out()
                    Printer.err('Error: Memory limit exceeded! {:1.2f}MB used, {:1.2f}MB allowed'.format(
                        memory_usage, self.memory_limit
                        )
                    )
                    self.terminated_cause = 'MEMORY_LIMIT'
                    self.terminated = True
                    self.process.secure_kill()
                    return
            # except NoSuchProcess as e1:
            #     pass
            except AttributeError as e2:
                pass


class ErrorMonitor(ThreadMonitor):
    """
    Class ErrorMonitor monitors PyPy object and prints detailed message if error occurs
    """

    def __init__(self, pypy):
        super(ErrorMonitor, self).__init__(pypy)
        self.message = 'Command failed'
        self.tail = 10

    @ensure_active
    def on_complete(self, pypy=None):
        if self.pypy.returncode > 0:
            if self.message:
                Printer.separator()
                Printer.open()
                Printer.out(self.message)
            else:
                Printer.open()

            # if file pointer exist try to read errors and outputs
            output = self.pypy.executor.output.read()
            if output:
                if self.pypy.full_output:
                    Printer.out('Output (last {} lines, rest in {}): ', self.tail, Paths.abspath(self.pypy.full_output))
                else:
                    Printer.out('Output (last {} lines): ', self.tail)
                Printer.wrn(format_n_lines(output, -self.tail, indent=Printer.indent * '    '))
            Printer.close()
