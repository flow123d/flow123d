#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from scripts.core.base import Printer, IO
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


class StartInfoMonitor(ThreadMonitor):
    """
    Class StartInfoMonitor monitors process at the beginning
    """

    def __init__(self, pypy):
        super(StartInfoMonitor, self).__init__(pypy)
        self.format = 'Executing "{self.pypy.escaped_command}"'

    @ensure_active
    def on_start(self, pypy=None):
        # if program did not start at all
        # nothing to do
        if pypy.executor.broken:
            return

        if self.format:
            Printer.all.out(self.format.format(**locals()))


class EndInfoMonitor(ThreadMonitor):
    """
    Class EndInfoMonitor monitors process at the end
    """
    def __init__(self, pypy):
        super(EndInfoMonitor, self).__init__(pypy)
        self.format = 'Command ({self.pypy.executor.process.pid}) ended with {self.pypy.returncode}'

    @ensure_active
    def on_complete(self, pypy=None):
        if self.format:
            Printer.all.out(self.format.format(**locals()))


class ProgressMonitor(ThreadMonitor):
    """
    Class ProgressMonitor monitors pypy process while it is running
    and reports progress
    """
    def __init__(self, pypy):
        super(ProgressMonitor, self).__init__(pypy)
        # set priority so that output is right after is finished and overrides
        # previous progress bar line
        self.pypy.on_process_complete.set_priority(self.on_complete, 5)
        self.timer = ProgressTime('Running | elapsed time {}')

    @ensure_active
    def on_update(self, pypy=None):
        if pypy.executor.broken:
            return

        self.timer.update()

    @ensure_active
    def on_complete(self, pypy=None):
        if pypy.executor.broken:
            return

        self.timer.format = 'Done    | elapsed time {}'
        self.timer.stop()


class LimitMonitor(ThreadMonitor):
    """
    Class LimitMonitor monitors resources of PyPy process and terminates
    PyPy if limits are not withheld.
    :type process: scripts.psutils.Process
    """

    CAUSE_TIME_LIMIT = 1
    CAUSE_MEMORY_LIMIT = 2

    def __init__(self, pypy):
        super(LimitMonitor, self).__init__(pypy)
        self.process = None
        self.memory_limit = None
        self.time_limit = None
        self.monitor_thread = None
        self.terminated = False
        self.terminated_cause = 0

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
                    Printer.console.newline()
                    Printer.all.err(
                        'Time limit exceeded! {:1.2f}s of runtime, {:1.2f}s allowed'.format(
                            runtime, self.time_limit
                        )
                    )
                    self.terminated_cause = self.CAUSE_TIME_LIMIT
                    self.terminated = True
                    self.process.secure_kill()
                    return
            except AttributeError as e2:
                pass

        if self.memory_limit:
            try:
                memory_usage = self.process.memory_usage()
                if memory_usage > self.memory_limit:
                    Printer.console.newline()
                    Printer.all.err('Memory limit exceeded! {:1.2f}MB used, {:1.2f}MB allowed'.format(
                        memory_usage, self.memory_limit
                    )
                    )
                    self.terminated_cause = self.CAUSE_MEMORY_LIMIT
                    self.terminated = True
                    self.process.secure_kill()
                    return
            # except NoSuchProcess as e1:
            #     pass
            except AttributeError as e2:
                pass


class ErrorMonitor(ThreadMonitor):
    def __init__(self, pypy):
        super(ErrorMonitor, self).__init__(pypy)
        self.pypy.on_process_complete.set_priority(self.on_complete, 4)
        self.message = None
        self.indent = 0

    @ensure_active
    def on_complete(self, pypy=None):
        if pypy.was_successful():
            return

        with Printer.all.with_level(self.indent):
            # print custom message on demand
            if self.message:
                Printer.all.err(self.message.format(**locals()))

            if pypy.executor.broken:
                Printer.all.err('Could not execute command! Possible cause is missing binary.')
                with Printer.all.with_level(3):
                    Printer.all.out('Full command: ')
                    Printer.all.out('{self.pypy.escaped_command}'.format(**locals()))

            elif pypy.with_error():
                Printer.all.err('Command ended with {self.pypy.returncode}! (pid={self.pypy.executor.process.pid})'.format(**locals()))
                with Printer.all.with_level(3):
                    Printer.all.out('Full command: ')
                    Printer.all.out('{self.pypy.escaped_command}'.format(**locals()))


class OutputMonitor(ThreadMonitor):
    POLICY_ALWAYS = 1
    POLICY_BATCH_OR_ERROR = 2
    POLICY_ERROR_ONLY = 3

    def __init__(self, pypy):
        super(OutputMonitor, self).__init__(pypy)
        self.log_file = None
        self._content = None
        self.policy = self.POLICY_BATCH_OR_ERROR

    @property
    def content(self):
        if self._content is None:
            self._content = self.pypy.executor.output.read()
        return self._content

    @ensure_active
    def on_complete(self, pypy=None):
        if pypy.executor.broken:
            return

        if self.log_file:
            IO.write(self.log_file, self.content)

        enabled = False
        if self.policy == self.POLICY_ALWAYS:
            enabled = True
        elif self.policy == self.POLICY_BATCH_OR_ERROR:
            enabled = pypy.with_error() or (not Printer.batched.is_muted())
        elif self.policy == self.POLICY_ERROR_ONLY:
            enabled = pypy.with_error()

        if enabled:
            with Printer.all.with_level():

                if self.pypy.full_output:
                    Printer.all.sep()
                    Printer.all.out('Output from file {self.pypy.full_output}'.format(**locals()))

                if not Printer.batched.is_muted():
                    Printer.batched.raw(format_n_lines(self.content, pypy.was_successful()))

                if not Printer.console.is_muted():
                    Printer.console.raw(format_n_lines(self.content, pypy.was_successful()))


# #
# # class InfoMonitor(ThreadMonitor):
# #     """
# #     Class InfoMonitor monitors PyPy object and prints info at the beginning
# #     and at the end. On error prints details.
# #     """
# #
# #     def __init__(self, pypy):
# #         """
# #         :type pypy: scripts.core.threads.PyPy
# #         """
# #         super(InfoMonitor, self).__init__(pypy)
# #         self.command_str = self.pypy.escaped_command
# #         self.start_fmt = 'Executing {self.command_str}'
# #         self.end_fmt = 'Command ({self.pypy.executor.process.pid}) ended with {self.pypy.returncode}'
# #
# #     @ensure_active
# #     def on_start(self, pypy=None):
# #         if pypy.executor.broken:
# #             return
# #
# #         if self.start_fmt:
# #             Printer.all.out(self.start_fmt.format(**dict(self=self)))
# #
# #     @ensure_active
# #     def on_complete(self, pypy=None):
# #
# #         # print either error that command failed or on_complete info id exists
# #         if self.pypy.returncode > 0:
# #             Printer.all.err('Command ({process.pid}) ended with {process.returncode}'.
# #                         format(rc=self.pypy.returncode, process=self.pypy.executor.process))
# #             Printer.all.out('Command: {}', self.pypy.escaped_command)
# #         elif self.end_fmt:
# #             Printer.all.out(self.end_fmt.format(**dict(self=self)))
# #
# #         if not self.pypy.progress:
# #             Printer.all.sep()
# #             output = self.pypy.executor.output.read()
# #             if output:
# #                 Printer.all.raw(output)
# #
# #
# # class LimitMonitor(ThreadMonitor):
# #     """
# #     Class LimitMonitor monitors resources of PyPy process and terminates
# #     PyPy if limits are not withheld.
# #     :type process: scripts.psutils.Process
# #     """
# #
# #     def __init__(self, pypy):
# #         super(LimitMonitor, self).__init__(pypy)
# #         self.process = None
# #         self.memory_limit = None
# #         self.time_limit = None
# #         self.monitor_thread = None
# #         self.terminated = False
# #         self.terminated_cause = None
# #
# #     def set_limits(self, case):
# #         """
# #         :type case: scripts.config.yaml_config.ConfigCase
# #         """
# #         # empty Limits object
# #         if not case:
# #             self.memory_limit = None
# #             self.time_limit = None
# #             return
# #
# #         self.memory_limit = case.memory_limit
# #         self.time_limit = case.time_limit
# #
# #     @ensure_active
# #     def on_start(self, pypy=None):
# #         self.process = self.pypy.executor.process
# #
# #     @ensure_active
# #     def on_update(self, pypy=None):
# #         if self.terminated:
# #             return
# #
# #         if self.time_limit:
# #             try:
# #                 runtime = self.process.runtime()
# #                 if runtime > self.time_limit:
# #                     Printer.all.out()
# #                     Printer.all.err(
# #                         'Error: Time limit exceeded! {:1.2f}s of runtime, {:1.2f}s allowed'.format(
# #                             runtime, self.time_limit
# #                         )
# #                     )
# #                     self.terminated_cause = 'TIME_LIMIT'
# #                     self.terminated = True
# #                     self.process.secure_kill()
# #                     return
# #             except AttributeError as e2:
# #                 pass
# #
# #         if self.memory_limit:
# #             try:
# #                 memory_usage = self.process.memory_usage()
# #                 if memory_usage > self.memory_limit:
# #                     Printer.all.out()
# #                     Printer.all.err('Error: Memory limit exceeded! {:1.2f}MB used, {:1.2f}MB allowed'.format(
# #                         memory_usage, self.memory_limit
# #                     )
# #                     )
# #                     self.terminated_cause = 'MEMORY_LIMIT'
# #                     self.terminated = True
# #                     self.process.secure_kill()
# #                     return
# #             # except NoSuchProcess as e1:
# #             #     pass
# #             except AttributeError as e2:
# #                 pass
# #
# #
# # class ErrorMonitor(ThreadMonitor):
# #     """
# #     Class ErrorMonitor monitors PyPy object and prints detailed message if error occurs
# #     """
# #
# #     def __init__(self, pypy):
# #         super(ErrorMonitor, self).__init__(pypy)
# #         self.message = 'Command failed'
# #         self.tail = 10
# #
# #     @ensure_active
# #     def on_complete(self, pypy=None):
# #         if self.pypy.returncode > 0:
# #             if self.message:
# #                 Printer.all.sep()
# #                 Printer.all.open()
# #                 Printer.all.out(self.message)
# #             else:
# #                 Printer.all.open()
# #
# #             # if file pointer exist try to read errors and outputs
# #             output = self.pypy.executor.output.read()
# #             if output:
# #                 if self.pypy.full_output:
# #                     Printer.all.out('Output (last {} lines, rest in {}): ', self.tail, Paths.abspath(self.pypy.full_output))
# #                 else:
# #                     Printer.all.out('Output (last {} lines): ', self.tail)
# #                 Printer.all.raw(format_n_lines(output, -self.tail, indent=Printer.all.indent))
# #             Printer.all.close()
