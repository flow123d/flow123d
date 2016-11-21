#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import os
import subprocess
import sys
import tempfile
# ----------------------------------------------
from scripts import psutils
from scripts.core.base import IO, Paths, Command
from scripts.core.threads import ExtendedThread, BrokenProcess
from utils.globals import ensure_iterable
# ----------------------------------------------


class OutputMode(object):
    """
    Class OutputMode helper class for redirecting output from a command
    """

    WRITE, APPEND, SHOW, HIDE, VARIABLE = range(5)

    def __init__(self, mode, filename=None, fp=None):
        self.mode = mode
        self.filename = filename
        self.fp = fp
        self.content = None

    def open(self):
        if self.mode in {self.SHOW, self.HIDE}:
            return {self.SHOW: None, self.HIDE: subprocess.PIPE}.get(self.mode)

        if self.mode in {self.WRITE, self.APPEND, self.VARIABLE}:

            # open file manually when append or write
            if self.mode in {self.WRITE, self.APPEND}:
                Paths.ensure_path(self.filename)
                self.fp = open(self.filename, 'w+' if self.mode is self.WRITE else 'a+')

            # create temp file otherwise
            if self.mode is self.VARIABLE:
                self.fp, self.filename = tempfile.mkstemp()
            return self.fp

    def close(self):
        if self.mode in {self.WRITE, self.APPEND, self.VARIABLE}:
            if self.fp is not None:
                if type(self.fp) is int:
                    os.close(self.fp)
                else:
                    self.fp.close()
                self.fp = None

        # remove temp file
        if self.mode in {self.VARIABLE}:
            self.content = IO.read(self.filename)
            os.unlink(self.filename)

    def read(self):
        if self.filename:
            if self.content is None:
                self.content = IO.read(self.filename)
            return self.content

    @classmethod
    def file_write(cls, filename):
        return OutputMode(cls.WRITE, filename)

    @classmethod
    def file_append(cls, filename):
        return OutputMode(cls.APPEND, filename)

    @classmethod
    def null_output(cls):
        return OutputMode(cls.SHOW)

    @classmethod
    def hidden_output(cls):
        return OutputMode(cls.HIDE)

    @classmethod
    def variable_output(cls):
        return OutputMode(cls.VARIABLE)


class BinExecutor(ExtendedThread):
    """
    Class which executes command and saves returncode
    :type process: scripts.psutils.Process
    :type threads: list[scripts.core.execution.BinExecutor]
    :type output : OutputMode
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
        self.output = OutputMode.hidden_output()
        self.exception = "Unknown error"

    def _run(self):
        if self.stopped:
            process = BrokenProcess(Exception('Application terminating'))
            self.returncode = process.returncode
            self.broken = True
            self.process = process
            return

        # run command and block current thread
        try:
            self.process = psutils.Process.popen(
                self.command,
                stdout=self.output.open(),
                stderr=subprocess.STDOUT
            )
        except Exception as e:
            # close file if necessary
            self.output.close()
            # broken process
            process = BrokenProcess(e)
            self.returncode = process.returncode
            self.broken = True
            self.process = process
            self.exception = str(e)
            return

        # process successfully started to wait for result
        # call wait on Popen process
        self.broken = False
        code = self.process.process.wait()
        self.output.close()

        if self.process.terminated:
            self.returncode = 5
            return
        self.returncode = getattr(self.process, 'returncode', None)

    @property
    def escaped_command(self):
        return Command.to_string(self.command)