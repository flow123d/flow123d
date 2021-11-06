#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import math
import time
import threading
# ----------------------------------------------
from loggers import printf
from scripts.core.base import IO
from scripts.serialization import ResultHolderResult, RuntestTripletResult, ResultParallelThreads
from utils.counter import ProgressCounter
from utils.events import Event
from scripts.core.returncode import RC, RC_NONE, RC_OK, RC_BROKEN
# ----------------------------------------------


class ProcessState(object):
    NOT_STARTED = 0
    STARTED = 1
    FINISHED = 2


class ExtendedThread(threading.Thread):
    """
    Class ExtendedThread is Thread class with extra functionality added
    """

    def __init__(self, name, target=None):
        super(ExtendedThread, self).__init__(name=name)
        self.started = None
        self.target = target

        self.state = ProcessState.NOT_STARTED
        self._returncode = RC_NONE
        self.exception = None

        self.start_time = None
        self.end_time = None

        # create event objects
        self.on_start = Event()
        self.on_complete = Event()
        self.on_update = Event()


    def _run(self):
        if self.target:
            self.exception = None
            try:
                self.target()
                self.returncode = RC_OK
            except BaseException as e:
                self.exception = e
                self.returncode = RC_BROKEN


    @property
    def returncode(self):
        """
        :rtype: scripts.core.returncode.RC
        """
        return self._returncode

    @property
    def duration(self):
        if self.start_time is None:
            return 0.0
        if self.end_time is None:
            return time.time() - self.start_time

        return self.end_time - self.start_time

    @returncode.setter
    def returncode(self, value):
        self._returncode = RC(value)

    def start(self):
        self.state = ProcessState.STARTED
        super(ExtendedThread, self).start()

    def run(self):
        self.state = ProcessState.STARTED
        self.on_start(self)
        self.start_time = time.time()
        self._run()
        self.end_time = time.time()
        self.on_complete(self)
        self.state = ProcessState.FINISHED

    def is_over(self):
        return self.state == ProcessState.FINISHED

    def is_running(self):
        return self.state == ProcessState.STARTED

    def __repr__(self):
        return "<{self.__class__.__name__}:{running} E:{self.returncode}>".format(
            self=self, running="RUNNING" if self.is_alive() else "STOPPED")

    def to_json(self):
        return dict(
            returncode=self.returncode(),
            name=self.name
        )

    def __bool__(self):
        return not self.with_error()

    def was_successful(self):
        """
        Return True if and only if returncode is 0
        :rtype: bool
        """
        return self.returncode == RC_OK

    def with_error(self):
        """
        Return True if returncode is not 0 and is not None
        :rtype: bool
        """
        return self.returncode not in (0, None)


class BrokenProcess(object):
    """
    Class BrokenProcess Dummy object when process execution fails
    """

    def __init__(self, exception=None):
        self.exception = exception
        self.pid = -1
        self.returncode = RC_BROKEN

    @staticmethod
    def is_running():
        return False


class MultiThreads(ExtendedThread):
    """
    Class MultiThreads is base class when executing threads in group
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
                printf.sep()
            self.counter.next(locals())

        self.current_thread.start()
        return True

    def add(self, thread):
        """
        :type thread: scripts.core.threads.ExtendedThread
        """
        self.threads.append(thread)
        self.returncodes[thread] = None
        try:
            thread.on_start += self.on_thread_start
            thread.on_complete += self.on_thread_complete
        except: pass

    @property
    def current_thread(self):
        return self.threads[self.index - 1]

    @property
    def returncode(self):
        # SKIPPED
        if not self.returncodes or set(self.returncodes.values()) == {None}:
            return RC_NONE

        # OK returncode
        if (set(self.returncodes.values()) - {None}) == {0}:
            return RC_OK

        # ERROR returncode
        return RC(max(set(self.returncodes.values()) - {None, 0}))

    @property
    def total(self):
        return len(self.threads)

    def on_thread_start(self, thread):
        """
        :type thread: scripts.core.threads.ExtendedThread
        """
        self.running += 1

    def on_thread_complete(self, thread):
        self.returncodes[thread] = RC(thread.returncode)
        self.running -= 1

    # aliases
    __len__ = total
    __iadd__ = add
    append = add


class SequentialThreads(MultiThreads):
    """
    Class SequentialThreads runs multiple threads in sequential fashion
    """

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

        with printf:
            while True:
                if not self.run_next():
                    break
                self.current_thread.join()

                if self.stop_on_error and self.current_thread.returncode > 0:
                    self.stopped = True
                    break

    def to_json(self):
        items = []
        for thread in self.threads:
            items.append(thread)

        return dict(
            returncode=self.returncode(),
            name=self.name,
            items=items
        )


class ParallelThreads(MultiThreads):
    """
    Class ParallelThreads run multiple threads in parallel fashion
    """

    def __init__(self, n=4, name='runner', progress=True, case_format=False):
        super(ParallelThreads, self).__init__(name, progress)
        self.n = n if type(n) is int else 1
        if case_format:
            self.counter = ProgressCounter('[{:02d} of {self.total:02d}] {self.current_thread.pypy}')
        else:
            self.counter = ProgressCounter(printer=None)
        self.stop_on_error = True
        self.i = 0

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
        steps = int(math.ceil(float(self.total) / self.n))

        # run each batched
        for i in range(steps):
            batch = [x for x in range(i * self.n, i * self.n + self.n) if x < len(self.threads)]
            for t in batch:
                self.run_next()

            for t in batch:
                self.threads[t].join()

            if self.stop_on_error and self.returncode > 0:
                self.stopped = True
                # no need to stop processes, just break
                break

    def dump(self):
        return ResultParallelThreads(self)


class ComparisonMultiThread(SequentialThreads):
    """
    Class ComparisonMultiThread hold comparison results and writes them to a file
    """

    def __init__(self, output, name='Comparison', progress=True, indent=True):
        super(ComparisonMultiThread, self).__init__(name, progress, indent)
        self.output = output

    def on_thread_complete(self, thread):
        """
        :type thread: scripts.core.pypy.PyPy
        """
        super(ComparisonMultiThread, self).on_thread_complete(thread)

        # append ndiff to file
        content = list()
        content.append('-' * 60 + '\n')
        content.append(thread.name + '\n')
        content.append('-' * 60 + '\n')
        content.append(thread.executor.output.read())
        content.append('\n' * 3)
        IO.append(self.output, '\n'.join(content) or '')

    def to_json(self):
        items = []
        for thread in self.threads:
            items.append(dict(
                name=thread.name,
                returncode=thread.returncode(),
                log=self.output,
            ))

        return dict(
            returncode=self.returncode(),
            name=self.name,
            log=self.output,
            tests=items,
        )


class RuntestMultiThread(SequentialThreads):
    """
    Class RuntestMultiThread is simple triplet holding single ConfigCase results
    (CleanThread, PyPy, ComparisonMultiThread)
    :type clean  : scripts.prescriptions.local_run.CleanThread
    :type pypy   : scripts.core.pypy.PyPy
    :type comp   : ComparisonMultiThread
    """

    def __init__(self, clean, pypy, comp):
        super(RuntestMultiThread, self).__init__('test-case', progress=False, indent=False)
        self.clean = clean
        self.pypy = pypy
        self.comp = comp

        self.add(clean)
        self.add(pypy)
        self.add(comp)

    def to_json(self):
        return dict(
            type="test-case",
            returncode=self.returncode(),
            clean=self.clean,
            execution=self.pypy,
            compare=self.comp
        )

    def dump(self):
        return RuntestTripletResult(self)


class ResultHolder(object):
    """
    Class ResultHolder stores object with property returncode
    """

    def __init__(self):
        self.items = list()

    def add(self, item):
        self.items.append(item)

    @property
    def returncode(self):
        rc = [] # TODO Python3 conversion
        for i in self.items:
            rc.append(i.returncode())
        return RC_NONE if not rc else RC(max(rc))

    def singlify(self):
        return self.items[0] if len(self.items) == 1 else self

    def dump(self):
        return ResultHolderResult(self)
