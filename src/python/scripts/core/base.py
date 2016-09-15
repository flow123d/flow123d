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
import time
import json
# ----------------------------------------------
from simplejson import JSONEncoder

is_linux = platform.system().lower().startswith('linux')

flow123d_name = "flow123d" if is_linux else "flow123d.exe"
mpiexec_name = "mpiexec" if is_linux else "mpiexec.hydra"


def find_base_dir():
    """
    Attempts to find base directory from pathfix's module's location
    """
    import os
    import pathfix
    path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(pathfix.__file__)), '..', '..'))
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


class Printer(object):
    """
    Class Printer unifies output operation
    """

    indent = 0
    batch_output = True
    dynamic_output = not batch_output

    @classmethod
    def ind(cls):
        return cls.indent * '    '

    @classmethod
    def style(cls, msg='', *args, **kwargs):
        sys.stdout.write(msg.format(*args, **kwargs))
        sys.stdout.write('\n')

    @classmethod
    def separator(cls):
        cls.out('-' * 60)

    @classmethod
    def wrn(cls, msg='', *args, **kwargs):
        sys.stdout.write(msg.format(*args, **kwargs))
        sys.stdout.write('\n')

    @classmethod
    def raw(cls, msg):
        sys.stdout.write(msg)
        sys.stdout.write('\n')

    @classmethod
    def err(cls, msg='', *args, **kwargs):
        if cls.indent:
            sys.stdout.write(cls.ind())
        sys.stdout.write(msg.format(*args, **kwargs))
        sys.stdout.write('\n')

    # ----------------------------------------------

    @classmethod
    def out(cls, msg='', *args, **kwargs):
        if cls.indent:
            sys.stdout.write(cls.ind())
        if not args and not kwargs:
            sys.stdout.write(msg)
        else:
            sys.stdout.write(msg.format(*args, **kwargs))
        sys.stdout.write('\n')

    @classmethod
    def dyn(cls, msg, *args, **kwargs):
        if cls.dynamic_output:
            sys.stdout.write('\r' + ' ' * 80)
            if cls.indent:
                sys.stdout.write('\r' + cls.ind() + msg.format(*args, **kwargs))
            else:
                sys.stdout.write('\r' + msg.format(*args, **kwargs))
            sys.stdout.flush()

    # ----------------------------------------------

    @classmethod
    def open(cls, l=1):
        cls.indent += l

    @classmethod
    def close(cls, l=1):
        cls.indent -= l


def make_relative(f):
    """
    Wrapper which return value as relative absolute or non-changed base on
    value Paths.format
    :param f:
    """
    def wrapper(*args, **kwargs):
        path = f(*args, **kwargs)
        if Paths.format == PathFormat.RELATIVE:
            return os.path.relpath(os.path.abspath(path), Paths.flow123d_root())
        elif Paths.format == PathFormat.ABSOLUTE:
            return os.path.abspath(path)
        return path
    return wrapper


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
        return lambda x: patt.match(x)\


    @staticmethod
    def filter_endswith(suffix=""):
        return lambda x: x.endswith(suffix)


class PathFormat(object):
    """
    Class PathFormat is enum class for different path formats
    """

    CUSTOM = 0
    RELATIVE = 1
    ABSOLUTE = 2


class Paths(object):
    """
    Class Paths is helper class when dealng with files and folders
    """

    _base_dir = find_base_dir()
    format = PathFormat.ABSOLUTE
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
    def test_paths(cls, *paths):
        status = True
        for path in paths:
            filename = getattr(cls, path)()
            if not cls.exists(filename):
                Printer.err('Error: file {:10s} ({}) does not exists!', path, filename)
                status = False

        return status

    @classmethod
    @make_relative
    def temp_file(cls, name='{date}-{time}-{rnd}.log'):
        return cls.path_to(cls.temp_name(name))

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
    @make_relative
    def bin_dir(cls):
        return cls.join(cls.flow123d_root(), 'bin')

    @classmethod
    @make_relative
    def ndiff(cls):
        return cls.path_to(cls.bin_dir(), 'ndiff', 'ndiff.pl')

    @classmethod
    @make_relative
    def flow123d(cls):
        return cls.path_to(cls.bin_dir(), flow123d_name)

    @classmethod
    @make_relative
    def mpiexec(cls):
        return cls.path_to(cls.bin_dir(), mpiexec_name)

    # -----------------------------------

    @classmethod
    @make_relative
    def path_to(cls, *args):
        return os.path.join(cls.current_dir(), *args)

    @classmethod
    @make_relative
    def join(cls, path, *paths):
        return os.path.join(path, *paths)

    @classmethod
    @make_relative
    def dirname(cls, path):
        return os.path.dirname(path)

    @classmethod
    @make_relative
    def without_ext(cls, path):
        return os.path.splitext(path)[0]

    @classmethod
    def browse(cls, path, filters=()):
        """
        :rtype: list[str]
        """
        paths = [cls.join(path, p) for p in os.listdir(path)]
        return cls.filter(paths, filters)

    @classmethod
    def walk(cls, path, filters=()):
        paths = list()
        for root, dirs, files in os.walk(path):
            for name in files:
                paths.append(cls.join(root, name))
            for name in dirs:
                paths.append(cls.join(root, name))

        return cls.filter(paths, filters)

    @classmethod
    def filter(cls, paths, filters=()):
        for f in filters:
            paths = [p for p in paths if f(p)]
        return paths

    @classmethod
    def match(cls, paths, filters):
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
        if not os.path.exists(p):
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
            with open(name, mode) as fp:
                return fp.read()

    @classmethod
    def write(cls, name, string, mode='w'):
        Paths.ensure_path(name)
        with open(name, mode) as fp:
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

    def __init__(self, min=100, max=5000, steps=13):
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

    def sleep(self):
        sleep_duration = self.next()
        time.sleep(sleep_duration)

    def next(self):
        self.current += 1

        if self.current >= self.total:
            return self.steps[-1]

        return self.steps[self.current]


class TestPrinterStatus(object):
    template = '{status_name:11s} | {case_name:40s} [{thread.duration:1.2f} sec] {detail}'
    default = 'failed'

    statuses = {
        '0':    'passed',
        'None': 'skipped',
        '-1': 'skipped',
    }

    errors = {
        'clean': '| could not clean directory',
        'pypy':  '| error while running',
        'comp':  '| wrong result: {sub_detail}'
    }

    @classmethod
    def get(cls, status):
        return cls.statuses.get(status, cls.default)

    @classmethod
    def detail_comp(cls, thread):
        """
        :type thread: scripts.core.threads.RuntestMultiThread
        """
        Printer.open(3)
        result = '\n'
        for t in thread.comp.threads:
            if t.returncode != 0:
                result += '{}ERROR in {}\n'.format(Printer.ind(), t.name)
        Printer.close(3)
        return result


class RunnerFormatter(object):
    template = '{status_name:11s} | passed={passed}, failed={failed}, skipped={skipped} in [{runner.duration:1.2f} sec]'


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
        for ti in ['clean', 'pypy', 'comp']:
            subthread = getattr(thread, ti)
            if subthread.returncode != 0:
                if hasattr(formatter, 'detail_{}'.format(ti)):
                    sub_detail = getattr(formatter, 'detail_{}'.format(ti))(thread)
                detail = formatter.errors[ti].format(**locals())
                break

        Printer.out(formatter.template.format(**locals()))
    
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
        Printer.out(formatter.template.format(**locals()))