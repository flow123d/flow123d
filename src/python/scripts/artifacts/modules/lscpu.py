#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import json
import re

from scripts.artifacts.command import Command
from scripts.core.base import Printer


class CommandLSCPU(Command):
    """
    Class CommandLSCPU runs lscpu command and export output to json
    """

    yaml_tag = u'!Command.lscpu'
    fields = dict(
        x64='Architecture',
        modes='CPU op-mode(s)',
        nproc='CPU(s)',
        physical='Core(s) per socket',
        logical='Thread(s) per core',
        sockets='Socket(s)',
        vendor='Vendor ID',
        frequency='CPU MHz',
        virtual='Virtualization',
        l1i='L1i',
        l1d='L1d',
        l2='L2',
        l3='L3',
    )

    list_split = re.compile(r',\s*')

    def __init__(self, output=None, **kwargs):
        if not output:
            raise Exception('Required filed "output" when using !Command.lscpu')

        self.output_file = output
        kwargs['output'] = None
        kwargs['command'] = ['lscpu']
        super(CommandLSCPU, self).__init__(**kwargs)

    def run(self):
        super(CommandLSCPU, self).run()
        lines = self.output.read().splitlines()
        result = dict()

        for line in lines:
            for field, phrase in self.fields.items():
                if line.startswith(phrase):
                    try:
                        result[field] = line.split(':')[1].strip()
                    except: pass

        self.try_expand(result, 'l1i', 'l1d', 'l2', 'l3')
        self.try_convert(result, int, 'physical', 'logical', 'nproc', 'sockets')
        self.try_convert(result, float, 'frequency')
        self.try_convert(result, lambda x: self.list_split.split(x), 'modes')
        self.try_convert(result, lambda x: str(x).find('64') != -1 , 'x64')

        # write file
        with open(self.output_file, 'w') as fp:
            json.dump(result, fp, indent=4)

        Printer.all.out('Written json file {}', self.output_file)

    @staticmethod
    def try_convert(result, method, *fields):
        for field in fields:
            value = result.get(field, None)
            if not value:
                continue

            try:
                value = method(value)
                result[field] = value
                continue
            except:
                pass

    @staticmethod
    def try_expand(result, *fields):
        for field in fields:
            value = result.get(field, '')
            if not value:
                continue

            try:
                if value.endswith('K'):
                    result[field] = int(value.strip('K')) * 1024
                    continue
            except:
                pass