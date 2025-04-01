#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import random
import sys
import os
import re
import platform
import datetime
import math
import json
import threading
from pathlib import Path
# ----------------------------------------------
from py123d.loggers import printf
from simplejson import JSONEncoder
from py123d.__init__ import py123d_package_dir 
# ----------------------------------------------


is_linux = platform.system().lower().startswith('linux')

flow123d_name = "flow123d" if is_linux else "flow123d.exe"
mpiexec_name = "mpiexec" if is_linux else "mpiexec.hydra"


def find_base_dir():
    """
    Attempts to find base directory from pathfix's module's location
    """
    import os
    path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..', '../..'))
    return path


class GlobalResult(object):
    """
    Class GlobalResult stores global result which can be dumped to json
    """

    items = []
    returncode = None
    error = None
    add = items.append

    @classmethod
    def to_json(cls, f=None):
        obj = dict(
            tests=cls.items,
            returncode=cls.returncode,
            error=cls.error
        )
        content = json.dumps(obj, indent=4, cls=MyEncoder)
        if f:
            with open(f, 'w') as fp:
                fp.write(content)
                print ('\n' * 10)
                print (content)
        return content


class MyEncoder(JSONEncoder):
    """
    Class MyEncoder tries to call to_json on object
     method before converting to string if there is no encoder available
    """

    def default(self, o):
        try:
            return o.to_json()
        except:
            return str(o)


class _Printer(object):
    """
    Class _Printer server as wrapper for all printing and echoing in the program
    """

    LEVEL_BATCH = 1
    LEVEL_CONSOLE = 2
    LEVEL_ALL = LEVEL_BATCH | LEVEL_CONSOLE
    SEPARATOR = '-' * 60
    log_file = None

    # default level
    level = LEVEL_ALL
    _indent = 0
    _indents = {}
    _depth = 0

    def __init__(self, level):
        self.level = level
        self._with_level = 1

    @classmethod
    def indent(cls):
        return '    ' * cls._indent

    def open(self, level=1):
        self.__class__._indent += level

    def close(self, level=1):
        self.__class__._indent -= level

    def with_level(self, level=1):
        self._with_level = level
        return self

    def __enter__(self):
        self.open(self._with_level)
        self.__class__._indents[self.__class__._depth] = self._with_level
        self.__class__._depth += 1

    def __exit__(self, type, value, traceback):
        self.__class__._depth -= 1
        i = self.__class__._indents[self.__class__._depth]

        self.close(i)
        return False

    def is_muted(self):
        return not self.level & self.__class__.level

    def _write(self, str):
        try:
            sys.stdout.write(str)
        except UnicodeEncodeError as e:
            print(str.encode('utf-8'))

        if self.__class__.log_file:
            with open(self.__class__.log_file, "a+") as fp:
                fp.write(str)

    # ------------------------------------------------------
    def out(self, msg, *args, **kwargs):
        if self.is_muted():
            return

        if self.__class__._indent:
            self._write(self.indent())
        if not args and not kwargs:
            self._write(msg)
        else:
            self._write(msg.format(*args, **kwargs))
        self._write('\n')

    def raw(self, msg):
        if self.is_muted():
            return

        self._write(msg)
        self._write('\n')

    def dyn(self, msg, *args, **kwargs):
        if self.is_muted():
            return

        self._write('\r' + ' ' * 80)
        if self.__class__._indent:
            self._write('\r' + self.indent() + msg.format(*args, **kwargs))
        else:
            self._write('\r' + msg.format(*args, **kwargs))
        sys.stdout.flush()

    def sep(self):
        if self.is_muted():
            return

        self._write(self.indent())
        self._write(self.SEPARATOR)
        self._write('\n')

    def suc(self, msg, *args, **kwargs):
        if self.is_muted():
            return

        self._write('[ OK ] |    ')
        self._write(msg.format(*args, **kwargs))
        self._write('\n')

    def err(self, msg, *args, **kwargs):
        if self.is_muted():
            return

        self._write(self.indent())
        self._write('[ERROR] |   ')
        self._write(msg.format(*args, **kwargs))
        self._write('\n')

    def wrn(self, msg, *args, **kwargs):
        if self.is_muted():
            return

        self._write(self.indent())
        self._write('[WARNING] | ')
        self._write(msg.format(*args, **kwargs))
        self._write('\n')

        # import traceback
        # traceback.print_stack(file=sys.stdout)

    def newline(self):
        if self.is_muted():
            return

        self._write('\n')


class Printer(object):
    """
    Class Printer is yet another wrapper for _Printer
    This class offers convenient fields for different print level
    """

    LEVEL_BATCH = _Printer.LEVEL_BATCH
    LEVEL_CONSOLE = _Printer.LEVEL_CONSOLE
    LEVEL_ALL = _Printer.LEVEL_ALL

    all = _Printer(_Printer.LEVEL_ALL)
    console = _Printer(_Printer.LEVEL_CONSOLE)
    batched = _Printer(_Printer.LEVEL_BATCH)

    @classmethod
    def indent(cls):
        return _Printer.indent()

    @classmethod
    def get_level(cls):
        return _Printer.level

    @classmethod
    def set_level(cls, level):
        _Printer.level = level

    @classmethod
    def setup_printer(cls, parser):
        """
        :type parser: utils.argparser.RuntestArgs
        """
        try:
            if parser.batch:
                cls.set_level(cls.LEVEL_BATCH)
            else:
                cls.set_level(cls.LEVEL_CONSOLE)
        except:
            cls.set_level(cls.LEVEL_CONSOLE)


class PathFilters(object):
    """
    Class PathFilters serves as filter library for filtering files and folders
    """

    @staticmethod
    def filter_name(name):
        return lambda x: Paths.basename(x) == name

    @staticmethod
    def filter_ext(ext):
        return lambda x: Paths.basename(x).endswith(ext)

    @staticmethod
    def filter_not(f):
        return lambda x: not f(x)

    @staticmethod
    def filter_type_is_file():
        return lambda x: os.path.isfile(x)

    @staticmethod
    def filter_type_is_dir():
        return lambda x: os.path.isdir(x)

    @staticmethod
    def filter_exists():
        return lambda x: os.path

    @staticmethod
    def filter_wildcards(fmt=""):
        fmt = fmt\
            .replace('.', r'\.')\
            .replace('*', r'.*')\
            .replace('/', r'\/')
        patt = re.compile(fmt)
        return lambda x: patt.match(x)

    @staticmethod
    def filter_dir_contains_file(required_file):
        def filter(file):
            files = Paths.browse(Paths.dirname(file), [PathFilters.filter_name(required_file)])
            return bool(files)

        return filter

    @staticmethod
    def filter_ignore_dirs(dirs):
        def filter(file):
            for d in dirs:
                if file.find(d) > 0:
                    return False
            return True

        return filter

    @staticmethod
    def filter_endswith(suffix=""):
        return lambda x: x.endswith(suffix)


class Paths(object):
    """
    Class Paths is helper class when dealing with files and folders
    """

    _base_dir = find_base_dir()
    cur_dir = os.getcwd()

    @classmethod
    def init(cls, v=None):
        if not v:
            return cls._base_dir

        if os.path.isfile(v):
            # if file is given, we assume file in bin/python was given
            cls._base_dir = os.path.dirname(os.path.dirname(os.path.realpath(v)))
        else:
            # if dir was given we just convert it to real path and use it
            cls._base_dir = os.path.realpath(v)
        return cls._base_dir

    @classmethod
    def current_dir(cls):
        """
        Returns path to current dir, where python was executed
        """
        return cls.cur_dir

    @classmethod
    def flow123d_root(cls):
        """
        Returns path to flow123d root
        """
        return cls._base_dir

    @classmethod
    def python_script_dir(cls):
        """
        Returns path to directory containing python scripts
        """
        return py123d_package_dir

    @classmethod
    def test_paths(cls, *paths):
        status = True
        for path in paths:
            filename = getattr(cls, path)()
            if not cls.exists(filename):
                printf.error('Error: file {:10s} ({}) does not exists!', path, filename)
                status = False

        return status

    @classmethod
    def temp_file(cls, name='{date}-{time}-{rnd}.log'):
        return str(Path(cls.current_dir(), cls.temp_name(name)).absolute())

    @classmethod
    def temp_name(cls, name='{date}-{time}-{rnd}.log'):
        today = datetime.datetime.today()
        time = '{:%H-%M-%S}'.format(today)
        date = '{:%Y-%m-%d}'.format(today)
        dt = '{:%Y-%m-%d_%H-%M-%S}'.format(today)
        rnd = '{:04d}'.format(random.randint(0, 9999))

        return name.format(date=date, time=time, today=today, datetime=dt, random=rnd, rnd=rnd)

    # -----------------------------------

    @classmethod
    def flow123d_dir(cls):
        """
        Returns path to flow123d root dir
        TODO: Simplify, remove try block and use system variable
        """
        try:
            import pathfix
            return Path().joinpath(py123d_package_dir, "../../")
        except ModuleNotFoundError:
            pass
        return Path("/opt/flow123d")

    @classmethod
    def flow123d_bin_dir(cls):
        """
        Returns path to flow123d bin dir containing Flow123d, mpiexec and ndiff
        """
        return Path(cls.flow123d_dir(), "bin")

    @classmethod
    def ndiff(cls):
        return Path(cls.flow123d_bin_dir(), 'ndiff', 'ndiff.pl')

    @classmethod
    def flow123d(cls):
        return Path(cls.flow123d_bin_dir(), flow123d_name)

    @classmethod
    def mpiexec(cls):
        return Path(cls.flow123d_bin_dir(), mpiexec_name)

    # -----------------------------------

    @classmethod
    def join(cls, path, *paths):
        return str(Path(path, *paths).absolute())

    @classmethod
    def dirname(cls, path):
        return str(Path(path).parent.absolute())

    @classmethod
    def without_ext(cls, path):
        return str(Path(path).with_suffix("").absolute())

    @classmethod
    def browse(cls, path, filters=()):
        """
        :rtype: list[str]
        """
        paths = [cls.join(path, p) for p in os.listdir(path)]
        return sorted(cls.filter(paths, filters))


    @classmethod
    def walk(cls, path, filters=()):
        paths = list()
        for root, dirs, files in os.walk(path):
            for name in files:
                paths.append(cls.join(root, name))
            for name in dirs:
                paths.append(cls.join(root, name))

        return sorted(cls.filter(paths, filters))

    @classmethod
    def filter(cls, paths, filters=()):
        for f in filters:
            paths = [p for p in paths if f(p)]
        return paths

    @classmethod
    def match(cls, paths, filters):
        """ Paths that match any of the given filters."""
        result = list()
        for p in paths:
            for f in filters:
                if f(p):
                    result.append(p)
                    break
        return result

    @classmethod
    def ensure_path(cls, f, is_file=True):
        if not f:
            return
        p = os.path.dirname(f) if is_file else f
        if p and not os.path.exists(p):
            os.makedirs(p)

    @classmethod
    def filesize(cls, path, as_string=False):
        size = os.path.getsize(path)
        if not as_string:
            return size

        units = ['B', 'kB', 'MB', 'GB', 'TB', 'PB']
        s = size
        for u in units:
            if s < 100:
                return '{:1.2f}{}'.format(s, u)
            s /= 1000.
        return '[Huge file]'

    @classmethod
    def path_end(cls, path, level=3):
        p = path
        for i in range(level):
            p = cls.dirname(p)
        return cls.relpath(path, p)

    @classmethod
    def path_end_until(cls, path, endswith, level=10):
        p = path
        for i in range(level):
            p = cls.dirname(p)
            if p.endswith(endswith):
                break
        return cls.relpath(path, p)

    @classmethod
    def split(cls, path):
        """
        :rtype: list[str]
        """
        path = cls.abspath(path)
        folders = []
        while 1:
            path, folder = os.path.split(path)
            if folder != "":
                folders.append(folder)
            else:
                if path != "":
                    folders.append(path)

                break
        folders.reverse()
        return folders

    # -----------------------------------

    @staticmethod
    def is_file(*args, **kwargs):
        return os.path.isfile(*args, **kwargs)

    @staticmethod
    def is_dir(*args, **kwargs):
        return os.path.isdir(*args, **kwargs)

    @staticmethod
    def exists(*args, **kwargs):
        return os.path.exists(*args, **kwargs)

    @staticmethod
    def abspath(*args, **kwargs):
        return os.path.abspath(*args, **kwargs)

    @staticmethod
    def relpath(*args, **kwargs):
        return os.path.relpath(*args, **kwargs)

    @staticmethod
    def realpath(*args, **kwargs):
        return os.path.realpath(*args, **kwargs)

    @staticmethod
    def basename(*args, **kwargs):
        return os.path.basename(*args, **kwargs)

    @staticmethod
    def unlink(*args, **kwargs):
        return os.unlink(*args, **kwargs)

    @classmethod
    def rename(cls, path, new_name):
        return cls.join(cls.dirname(path), new_name)


class Command(object):
    """
    Class Command help quote arguments passed to command
    """

    @classmethod
    def escape_command(cls, command):
        """
        :rtype : list[str]
        :type command: list[str]
        """
        import pipes
        return [pipes.quote(x) for x in command]

    @classmethod
    def to_string(cls, command):
        return ' '.join(cls.escape_command(command))


class IO(object):
    """
    Class IO is helper class fo IO operations
    """

    @classmethod
    def read(cls, name, mode='r'):
        """
        :rtype : str or None
        """
        if name and Paths.exists(name):
            with open(name, mode, encoding='utf-8', errors='replace') as fp:
                return fp.read()

    @classmethod
    def write(cls, name, string, mode='w'):
        Paths.ensure_path(name)
        with open(name, mode, encoding='utf-8', errors='replace') as fp:
            fp.write(string)
        return True

    @classmethod
    def append(cls, name, string, mode='a'):
        return cls.write(name, string, mode)

    @classmethod
    def delete(cls, name):
        if Paths.exists(name):
            Paths.unlink(name)
            return True
        return False

    @classmethod
    def delete_all(cls, folder):
        import shutil
        return shutil.rmtree(folder, ignore_errors=True)


class DynamicSleep(object):
    """
    Class DynamicSleep extends sleep duration each time sleep method
    is called. This is useful when we want to have finer resolution at the
    beginning of an operation.
    """

    def __init__(self, min=250, max=1*3000, steps=50):
        # -c * Math.cos(t/d * (Math.PI/2)) + c + b;
        # t: current time, b: begInnIng value, c: change In value, d: duration
        c = float(max - min)
        d = float(steps)
        b = float(min)
        self.steps = list()
        for t in range(steps + 1):
            self.steps.append(float(int(
                -c * math.cos(t / d * (math.pi / 2)) + c + b
            )) / 1000)
        self.current = -1
        self.total = len(self.steps)
        self.event = threading.Event()

    def sleep(self):
        self.event.clear()
        sleep_duration = self.next()
        self.event.wait(sleep_duration)

    def next(self):
        self.current += 1

        if self.current >= self.total:
            return self.steps[-1]

        return self.steps[self.current]


class TestPrinterStatus(object):
    template = '{status_name:11s} | {case_name:45s} [{thread.duration:5.2f} sec] {detail}'
    default = 'failed'

    statuses = {
        '0':    'passed',
        'None': 'skipped',
        '-1':   'skipped',
        '100':  'death',
    }

    errors = {
        'clean': '{sub_detail}',
        'pypy':  '{sub_detail}',
        'comp':  '{sub_detail}'
    }

    @classmethod
    def get(cls, status):
        return cls.statuses.get(status, cls.default)

    @classmethod
    def detail_comp(cls, thread):
        """
        :type thread: scripts.core.threads.RuntestMultiThread
        """
        if thread.comp.returncode.succeeded:
            return ''

        result = '| wrong result: '
        for t in thread.comp.threads:
            if t.returncode != 0:
                result += '\n{:18s}FAILED comparison in {}'.format('', t.name)
        return result

    @classmethod
    def detail_pypy(cls, thread):
        """
        :type thread: scripts.core.threads.RuntestMultiThread
        """
        if thread.pypy.returncode.succeeded:
            if thread.pypy.returncode.reversed:
                return '| OK death test failed'
            return ''
        elif thread.pypy.returncode.failed:
            if thread.pypy.returncode.reversed:
                return '| test did NOT failed but should have failed'
            return '| error while execution'

        # skipped
        return ''

    @classmethod
    def detail_clean(cls, thread):
        """
        :type thread: scripts.core.threads.RuntestMultiThread
        """
        if thread.clean.returncode.failed:
            return '| could not clean directory'
        return ''


class RunnerFormatter(object):
    template = '{status_name:11s} | passed={passed}, failed={failed}, skipped={skipped} in              [{runner.duration:5.2f} sec]'


class StatusPrinter(object):

    @classmethod
    def print_test_result(cls, thread, formatter=TestPrinterStatus):
        """
        :type formatter: TestPrinterStatus
        :type thread: scripts.core.threads.RuntestMultiThread
        """
        status_name = '[ {} ]'.format(formatter.get(str(thread.returncode))).upper()
        case_name = thread.pypy.case.as_string

        detail = ''
        thread_names = ['clean', 'pypy', 'comp']
        thread_rcs = dict()
        thread_msg = dict()

        for ti in thread_names:
            subthread = getattr(thread, ti)
            thread_rcs[ti] = subthread.returncode
            thread_msg[ti] = getattr(formatter, 'detail_{}'.format(ti))(thread)

        # first non zero rc will be formatted
        # otherwise first non empty detail
        for ti in thread_names:
            subthread = getattr(thread, ti)
            if subthread.returncode.failed:
                sub_detail = thread_msg[ti]
                detail = formatter.errors[ti].format(**locals())
                printf.error(formatter.template.format(**locals()))
                return

        # first non empty detail
        for ti in thread_names:
            subthread = getattr(thread, ti)
            sub_detail = thread_msg[ti]
            if sub_detail:
                detail = formatter.errors[ti].format(**locals())
                printf.out(formatter.template.format(**locals()))
                return

        if str(thread.returncode).upper() == 'NONE':
            printf.warning(formatter.template.format(**locals()))
        else:
            printf.success(formatter.template.format(**locals()))

    @classmethod
    def print_runner_stat(cls, runner, formatter=RunnerFormatter):
        """
        :type runner: scripts.core.threads.ParallelThreads
        :type formatter: RunnerFormatter
        """
        returncodes = [t.returncode for t in runner.threads]

        skipped = returncodes.count(None) + returncodes.count(-1)
        passed = returncodes.count(0)
        failed = len(returncodes) - (skipped + passed)

        result = 0 if len(returncodes) == passed else 1
        status_name = '[ {} ]'.format(TestPrinterStatus.get(str(result))).upper()
        if result == 0:
            printf.success(formatter.template.format(**locals()))
        else:
            printf.error(formatter.template.format(**locals()))


class System(object):
    import uuid
    import time as _time

    time = _time.strftime("%H-%M-%S")
    date = _time.strftime("%y.%m.%d")
    datetime = _time.strftime("%y.%m.%d_%H-%M-%S")
    started = _time.time()

    rnd8 = ''.join(random.choice('0123456789ABCDEF') for i in range(8))
    rnd16 = ''.join(random.choice('0123456789ABCDEF') for i in range(16))
    rnd32 = ''.join(random.choice('0123456789ABCDEF') for i in range(32))
    rnd = str(uuid.uuid4().hex)
