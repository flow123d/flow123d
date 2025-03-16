#!/bin/python3
# author: Jan Hybs
import json
import time

from py123d.loggers import printf
from py123d.scripts.core import monitors
from py123d.scripts.core.base import DynamicSleep, Command, Paths, IO, PathFilters
from py123d.scripts.core.returncode import RC
from py123d.scripts.core.threads import ExtendedThread
from py123d.scripts.serialization import PyPyResult
from py123d.utils.events import Event
from py123d.utils.globals import wait_for


class PyPy(ExtendedThread):
    """
    Class PyPy is main class which executes command having multiple monitors registered
    PyPy = BinExecutor + monitors

    :type executor : scripts.core.execution.BinExecutor
    :type case     : ConfigCase
    """

    returncode_map = {
        '0': 'SUCCESS',
        '1': 'ERROR',
        '5': 'TERM',
        'None': 'SKIP',
        '-1': 'SKIP',
    }

    def __init__(self, executor, period=.5):
        super(PyPy, self).__init__(name='pypy')
        self.executor = executor
        self.period = period
        self.case = None
        self.status_file = None
        self.extra = dict()
        self._short_sleep = False

        self.on_process_start = Event()
        self.on_process_complete = Event()
        self.on_process_update = Event()

        self.on_process_complete += self.generate_status_file
        self.monitor = monitors.MainMonitor(self)

        self.log = False
        self.custom_error = None

        # dynamic sleeper
        self.sleeper = DynamicSleep()

        # path to full output
        self.full_output = None

    def wakeup(self):
        try:
            self._short_sleep = True
            self.sleeper.event.set()
        except:
            pass

    @property
    def escaped_command(self):
        return self.executor.escaped_command

    def _run(self):
        # start executor
        self.executor.start()
        wait_for(self.executor, 'process')

        if self.executor.broken:
            printf.error('Could not start command "{}": {}',
                            Command.to_string(self.executor.command),
                            self.executor.exception)
            self.returncode = RC(self.executor.returncode)

        # if process is not broken, propagate start event
        self.on_process_start(self)

        while self.executor.is_running():
            self.on_process_update(self)
            if self._short_sleep:
                time.sleep(0.01)
            else:
                self.sleeper.sleep()

        # get return code
        rc = getattr(self.executor, 'returncode', None)
        self.returncode = RC(rc if self.custom_error is None or str(rc) == "0" else self.custom_error)

        # reverse return code if death test is set
        if self.case and self.case.death_test is not None:
            self.returncode = RC(self.case.death_test.reverse_return_code(self.returncode))
            self.returncode.reversed = self.case.death_test.value != self.case.death_test.FALSE

        # propagate on_complete event
        self.on_process_complete(self)

    def status(self):
        import getpass
        import platform
        result = dict(
            returncode=self.returncode(),
            duration=self.duration,
            username=getpass.getuser(),
            hostname=platform.node(),
            nodename=platform.node().split('.')[0].strip('0123456789'),
            commit=self.get_commit(),
        )

        if self.case:
            result.update(self.case.info)
        return result

    def to_json(self):
        if self.case:
            return dict(
                returncode=self.returncode(),
                name=self.case.as_string,
                case=self.case,
                log=self.full_output
            )
        json = super(PyPy, self).to_json()
        json['log'] = self.full_output
        json['type'] = 'exec'
        return json

    def dump(self):
        return PyPyResult(self)

    @classmethod
    def get_commit(cls):
        """
        Calls git show on git root to determine unix timestamp of the current commit (HEAD)
        :return:
        """
        import subprocess
        try:
            root = Paths.flow123d_root()
            # get current hash(%H) and date(%ct) from git repo
            result = subprocess.check_output('git show -s --format=%H,%ct HEAD'.split(), cwd=root).decode()
            sha, date = str(result).strip().split(',')
            return dict(
                hash=sha,
                date=int(date)
            )
        except:
            return None

    @classmethod
    def generate_status_file(cls, target):
        """
        Will generate status file if target has option turned on
        :type target: PyPy
        """
        if target.status_file:
            IO.write(
                target.status_file,
                json.dumps(target.status(), indent=4)
            )
            output_dir = Paths.dirname(target.status_file)
            files = Paths.browse(
                output_dir,
                [PathFilters.filter_wildcards('*/profiler_info_*.log.json')]
            )
            # profiler json is missing?
            if not files:
                IO.write(Paths.join(output_dir, 'profiler_info_dummy.log.json'), '{}')

    def __repr__(self):
        if self.case:
            return str(self.case)
        return super().__repr__()
