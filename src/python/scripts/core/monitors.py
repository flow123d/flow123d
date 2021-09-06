#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from loggers import printf
from scripts.core.base import IO
from utils.globals import ensure_iterable
from utils.strings import format_n_lines
from utils.cache import ttl_cache

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
    :type pypy: scripts.core.pypy.PyPy
    """

    def __init__(self, pypy):
        """
        :type pypy: scripts.core.pypy.PyPy
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


class AsciiAnimation(object):
    def __init__(self, frames=None):
        self._frames = frames or ['····', '»···', '·»··', '··»·', '···»']  # cli spinner
        self._frames_len = len(self._frames)
        self._index = 0

    @property
    def animation(self):
        frame = self._frames[self._index]
        self._index = (self._index + 1) % self._frames_len
        return frame


class MainMonitor(ThreadMonitor):
    """
    :type process: scripts.psutils.linux_psutil.Process
    """
    CAUSE_TIME_LIMIT = 1
    CAUSE_MEMORY_LIMIT = 2

    def __init__(self, pypy):
        """
        :type pypy: scripts.core.pypy.PyPy
        """
        super(MainMonitor, self).__init__(pypy)
        self.start_format = None
        self.complete_format = None
        self.color_complete_format = None
        self.error_complete_format = None
        self.update_format = None

        self.memory_limit = None
        self.time_limit = None
        self.proc = 1
        self.terminated = False
        self.terminated_cause = 0
        self.memory_usage = 0
        self.memory_used = 0
        self.runtime = 0
        self.process = None
        self._content = None
        self.log_file = None
        self._animation = None

    @property
    def animation(self):
        if not self._animation:
            self._animation = AsciiAnimation()
        return self._animation.animation

    @ensure_active
    def on_start(self, pypy=None):
        self.process = self.pypy.executor.process

        for fmt in ensure_iterable(self.start_format):
            printf.out(fmt, monitor=self)

    @ensure_active
    def on_update(self, pypy=None):
        self._check_limits()

        if self.terminated:
            return

        for fmt in ensure_iterable(self.update_format):
            if printf.tty():
                printf.rewrite(fmt, monitor=self)
            else:
                printf.out(fmt, monitor=self)

    @property
    def content(self):
        if self._content is None:
            self._content = self.pypy.executor.output.read()
        return self._content

    @ensure_active
    def on_complete(self, pypy=None):

        if self.log_file:
            IO.write(self.log_file, self.content)

        # finish after update
        if self.update_format and (not self.complete_format and not self.color_complete_format):
            printf.finish_rewrite()

        # print regular messages
        for fmt in ensure_iterable(self.complete_format):
            printf.out(fmt, monitor=self)

        # print messages if error
        if self.pypy.with_error():
            for fmt in ensure_iterable(self.error_complete_format):
                printf.error(fmt, monitor=self)

        # print regular color messages based on result
        for fmt in ensure_iterable(self.color_complete_format):
            rc = self.pypy.returncode()
            if rc == 0:
                printf.success(fmt, monitor=self)
            elif rc is None:
                printf.warning(fmt, monitor=self)
            else:
                printf.error(fmt, monitor=self)

        with printf:
            if printf.verbosity() is printf.OutputVerbosity.FULL:
                printf.sep()
                printf.out('Output from file {self.pypy.full_output}'.format(**locals()))
                printf.opt(raw=True).stream(
                    format_n_lines(self.content, pypy.was_successful())
                )
            elif printf.verbosity() is printf.OutputVerbosity.SMART and pypy.with_error():
                printf.sep()
                printf.out('Last 50 lines from file {self.pypy.full_output}'.format(**locals()))
                msg = format_n_lines(self.content, success=False, n_lines=-50)
                printf.opt(raw=True).stream(msg)
            elif printf.verbosity() is printf.OutputVerbosity.MINIMAL and pypy.with_error():
                printf.sep()
                printf.out('Last 50 lines from file {self.pypy.full_output}'.format(**locals()))
                printf.opt(raw=True).stream(
                    format_n_lines(self.content, success=False, n_lines=-50)
                )

    def process_runtime(self):
        try:
            return self.process.runtime()
        except:
            return 0

    @property
    def process_runtime_formatted(self):
        duration = self.process_runtime()
        msec = int((duration - int(duration)) * 1000)
        mins = int(duration) // 60
        secs = int(duration) % 60

        return '{:02d}:{:02d}.{:03d}'.format(mins, secs, msec)

    @property
    def process_memory_usage_formatted(self):
        usage = int(self.process_memory_usage())
        return '{:4d} MiB'.format(usage)

    @property
    def process_memory_used_formatted(self):
        usage = int(self.memory_used)
        return '{:4d} MiB'.format(usage)

    @ttl_cache(ttl=5.0)
    def process_memory_usage(self):
        try:
            return self.process.memory_usage()
        except:
            return self.memory_used

    def _check_limits(self):
        if self.terminated:
            return

        if self.time_limit:
            try:
                self.runtime = self.process_runtime()
                if self.runtime > self.time_limit:
                    printf.finish_rewrite()
                    printf.error(
                        'Time limit exceeded! {:1.2f}s of runtime, {:1.2f}s allowed',
                        self.runtime, self.time_limit
                    )
                    printf.error('Process will now be terminated...')
                    self.terminated_cause = self.CAUSE_TIME_LIMIT
                    self.terminated = True
                    self.process.secure_kill()
                    self.pypy.wakeup()
                    return
            except AttributeError as e2:
                pass

        if self.memory_limit:
            try:
                self.memory_usage = self.process_memory_usage()
                self.memory_used = max(self.memory_used, self.memory_usage)
                if self.memory_usage > self.memory_limit:
                    printf.finish_rewrite()
                    printf.error(
                        'Memory limit exceeded! {:1.2f}MB used, {:1.2f}MB allowed ({:1.2f}MB per cpu)',
                        self.memory_usage, self.memory_limit, self.memory_limit / self.proc
                    )
                    printf.error('Process will now be terminated...')
                    self.terminated_cause = self.CAUSE_MEMORY_LIMIT
                    self.terminated = True
                    self.process.secure_kill()
                    self.pypy.wakeup()
                    return
            # except NoSuchProcess as e1:
            #     pass
            except AttributeError as e2:
                pass

    def set_limits(self, case):
        """
        :type case: scripts.yamlc.yaml_config.ConfigCase
        """

        # empty Limits object
        if not case:
            self.memory_limit = None
            self.time_limit = None
            return

        self.memory_limit = case.memory_limit
        self.time_limit = case.time_limit

        if case.proc > 1:
            self.proc = case.proc
            self.memory_limit = case.memory_limit * self.proc
