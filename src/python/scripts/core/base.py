#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from __future__ import absolute_import
import random
import sys
import os
import re
import platform

is_linux = platform.system().lower().startswith('linux')

flow123d_name = "flow123d" if is_linux else "flow123d.exe"
mpiexec_name = "mpiexec" if is_linux else "mpiexec.hydra"

import datetime


class Printer(object):
    LEVEL_DBG = 00
    LEVEL_KEY = 10
    LEVEL_WRN = 20
    LEVEL_ERR = 30
    indent = 0

    def __init__(self, level=LEVEL_DBG):
        self.level = level

    def dbg(self, *args, **kwargs):
        if self.level <= self.LEVEL_DBG:
            self.out(*args, **kwargs)

    def key(self, *args, **kwargs):
        if self.level <= self.LEVEL_KEY:
            self.out(*args, **kwargs)

    def line(self):
        self.key('-' * 60)

    def wrn(self, *args, **kwargs):
        if self.level <= self.LEVEL_WRN:
            self.out(*args, **kwargs)

    def err(self, *args, **kwargs):
        if self.level <= self.LEVEL_ERR:
            self.out(*args, **kwargs)

    def copy(self):
        return Printer(self.level)
    # ----------------------------------------------

    def out(self, msg='', *args, **kwargs):
        if self.indent:
            sys.stdout.write('    ' * self.indent)
        sys.stdout.write(msg.format(*args, **kwargs))
        sys.stdout.write('\n')

    def out_r(self, msg, *args, **kwargs):
        sys.stdout.write(msg.format(*args, **kwargs))

    def out_rr(self, msg, *args, **kwargs):
        sys.stdout.write(msg.format(*args, **kwargs))
        sys.stderr.write('\r')
        sys.stdout.flush()

    def dyn(self, msg, *args, **kwargs):
        sys.stdout.write('\r' + ' ' * 60)
        sys.stdout.write('\r' + msg.format(*args, **kwargs))
        sys.stdout.flush()
    # ----------------------------------------------

    @classmethod
    def open(cls):
        cls.indent += 1

    @classmethod
    def close(cls):
        cls.indent -= 1


def make_relative(f):
    def wrapper(*args, **kwargs):
        path = f(*args, **kwargs)
        if Paths.format == PathFormat.RELATIVE:
            return os.path.relpath(os.path.abspath(path), Paths.base_dir())
        elif Paths.format == PathFormat.ABSOLUTE:
            return os.path.abspath(path)
        return path
    return wrapper


class PathFilters(object):
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


class PathFormat(object):
    CUSTOM = 0
    RELATIVE = 1
    ABSOLUTE = 2


class Paths(object):
    _base_dir = ''
    format = PathFormat.ABSOLUTE
    printer = Printer(Printer.LEVEL_WRN)

    @classmethod
    def base_dir(cls, v=None):
        if v is None:
            return cls._base_dir
        cls._base_dir = os.path.dirname(os.path.realpath(v))

    @classmethod
    def source_dir(cls):
        return cls.join(cls.dirname(__file__), '..', '..')

    @classmethod
    def test_paths(cls, *paths):
        status = True
        for path in paths:
            filename = getattr(cls, path)()
            if not cls.exists(filename):
                cls.printer.err('Error: file {:10s} ({}) does not exists!', path, filename)
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
        return cls.join(cls.base_dir(), '..')

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
        return os.path.join(cls.base_dir(), *args)

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
    @classmethod
    def read(cls, name, mode='r'):
        """
        :rtype : str or None
        """
        if Paths.exists(name):
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
