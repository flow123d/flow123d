#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from scripts.core.base import Printer, IO
from scripts.core.threads import PyPy
from utils.dotdict import Map
from utils.strings import format_n_lines


class RuntestParser(object):
    def __init__(self, obj, batch):
        self.n_lines = 0 if batch else 15
        self.batch = batch
        self.pypy = self.create_from_fields(
            obj.get('execution', {}),
            returncode=int,
            log=str,
            name=str,
            case=lambda x: self.create_from_fields(x, cpu=int, test_case=str)
        )
        self.clean = self.create_from_fields(
            obj.get('clean', {}),
            returncode=int,
            name=str,
            dir=str,
            error=str
        )
        self.comp = self.create_from_fields(
            obj.get('compare', {}),
            returncode=int,
            name=str,
            tests=lambda x: [Map(y) for y in x]
        )
        self.returncode = int(obj.get('returncode'))

    def get_result(self):
        if self.clean.returncode != 0:
            Printer.out(
                "{} Could not clean directory '{c[dir]}': {c[error]}",
                self.get_status_line(self.clean),
                c=self.clean)
            return

        if self.pypy.returncode != 0:
            Printer.out("{} Run error, case: {p[name]}", self.get_status_line(self.pypy), p=self.pypy)
            self.print_log_file(self.pypy.log, self.n_lines)
            return
        elif self.batch:
            self.print_log_file(self.pypy.log, self.n_lines)

        if self.comp.returncode not in (0, None):
            Printer.out("{} Compare error, case: {p[name]}, Details: ", self.get_status_line(self.comp), p=self.pypy)
            self.print_log_file(self.pypy.log, self.n_lines)
            Printer.open(2)
            for c in self.comp.tests:
                rc = c.returncode
                if rc == 0:
                    Printer.out('[{:^6}]: {}', 'OK', c.name)
                else:
                    Printer.out('[{:^6}]: {}', 'FAILED', c.name)
            Printer.close(2)
            return
        elif self.batch:
            self.print_log_file(self.comp.log, self.n_lines)

    @classmethod
    def create_from_fields(cls, json, **fields):
        obj = Map()
        for k, v in fields.items():
            value = json.get(k, None)
            obj[k] = v(value) if value is not None else None
        return obj

    @classmethod
    def get_status_line(cls, o, map=False):
        if not map:
            return '[{:^6}]:{o.returncode:3} |'.format('ERROR', o=o)
        return '[{:^6}]:{o.returncode:3} |'.format(
            PyPy.returncode_map.get(str(o['returncode']), 'ERROR'), o=o)

    @classmethod
    def print_log_file(cls, f, n_lines):
        log_file = IO.read(f)
        if log_file:
            if n_lines == 0:
                Printer.out('Full log from file {}:', f)
            else:
                Printer.out('Last {} lines from file {}:', abs(n_lines), f)

            Printer.wrn(format_n_lines(log_file.rstrip(), -n_lines, indent=Printer.indent * '    '))


class ExecParser(object):
    def __init__(self, obj, batch):
        self.returncode = obj.get('returncode', None)
        self.name = obj.get('name', None)
        self.log = obj.get('log', None)
        self.n_lines = 0 if batch else 15
        self.batch = batch

    def get_result(self):
        if self.returncode != 0:
            if self.returncode != 0:
                Printer.out("{} Run error, case: {p.name}", RuntestParser.get_status_line(self), p=self)
                RuntestParser.print_log_file(self.log, self.n_lines)
                return
        elif self.batch:
            RuntestParser.print_log_file(self.log, self.n_lines)


class JsonParser(object):
    """
    :type tests : list[RuntestParser]
    """
    parse_map = {
        'test-case': RuntestParser,
        'exec': ExecParser
    }

    def __init__(self, obj, batch):
        self.batch = batch
        self.returncode = int(obj.get('returncode'))
        self.error = obj.get('error')
        self.tests = list()

        for t in obj.get('tests', []):
            self.tests.append(self.parse_map.get(t['type'])(t, batch))
