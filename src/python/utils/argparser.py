#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import argparse
# ----------------------------------------------
from scripts.core.base import System
from utils.events import Event


def try_split(lst, sep='--'):
    """
    Method wil try to split list of items by certain value
    given list ['foo', '--', 'bar']
    will produce tuple of two lists (['foor'], ['bar'])
    :rtype: (list[str], list[str])
    :type sep: str
    :type lst: list[str]
    """
    if sep not in lst:
        return lst, []

    i = lst.index(sep)
    return lst[:i], lst[i+1:]


class ParserArgs(object):
    """
    Class ParserArgs is abstract arguments wrapper

     :type rest              : list[str]
     :type args              : list[str]
    """
    def __init__(self, args, rest):
        """
        :type args: argparse.Namespace
        :type rest: list[str]
        """
        self.args, self.rest = try_split(rest)

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return self.__dict__.get(attr)
        return None

    def get(self, attr, default=None):
        if attr in self.__dict__:
            return self.__dict__.get(attr)
        return default

    def __repr__(self):
        from utils.strings import format_dict

        return '<{self.__class__.__name__} \n{y}\n    >'.format(
            self=self,
            y=format_dict(self.__dict__, indent=1)
        )

class ExecWithLimitArgs(ParserArgs):
    """
    Class ExecWithLimitArgs is wrapper for argparse parse exec_with_limit

     :type time_limit        : float
     :type memory_limit      : float
     :type batch             : bool
    """

    def __init__(self, args, rest):
        """
        :type args: argparse.Namespace
        :type rest: list[str]
        """
        super(ExecWithLimitArgs, self).__init__(args, rest)
        self.time_limit = None
        self.memory_limit = None
        self.batch = None
        self.death = None

        for k, v in args._get_kwargs():
            setattr(self, k, v)


class ExecParallelArgs(ParserArgs):
    """
    Class RuntestArgs is wrapper for argparse parse result

     :type valgrind          : bool|str
     :type parallel          : int
     :type batch             : bool

     :type cpu               : list[int]
     :type queue             : bool|str
     :type host              : str

     :type time_limit        : float
     :type memory_limit      : float

     :type root              : str
     :type json              : str
     :type list              : bool
     :type dump              : str
     :type log               : str
    """

    def __init__(self, args, rest):
        """
        :type args: argparse.Namespace
        :type rest: list[str]
        """
        super(ExecParallelArgs, self).__init__(args, rest)
        self.valgrind = None
        self.parallel = None
        self.batch = None

        self.cpu = None
        self.queue = None
        self.host = None

        self.time_limit = None
        self.memory_limit = None

        self.root = None
        self.json = None
        self.list = None
        self.dump = None
        self.log = None

        for k, v in args._get_kwargs():
            setattr(self, k, v)

        # fix more complex args
        self.queue = True if self.queue is None else self.queue
        self.valgrind = True if self.valgrind is None else self.valgrind
        self.cpu = [x for y in self.cpu for x in y] if self.cpu else self.cpu

        # default value for CPU
        if not self.cpu:
            self.cpu = [1]


class RuntestArgs(ParserArgs):
    """
    Class RuntestArgs is wrapper for argparse parse result

     :type keep_going        : bool
     :type valgrind          : bool|str
     :type parallel          : int
     :type batch             : bool
     :type include           : list
     :type exclude           : list

     :type cpu               : list[int]
     :type queue             : bool|str
     :type host              : str

     :type time_limit        : float
     :type memory_limit      : float

     :type root              : str
     :type json              : str
     :type list              : bool
     :type dump              : str
     :type log               : str

     :type export            : bool
     :type random_output_dir : str
     :type no_clean          : bool
     :type no_compare        : bool
     :type death_test        : bool
     :type status_file       : bool
    """

    def __init__(self, args, rest):
        """
        :type args: argparse.Namespace
        :type rest: list[str]
        """
        super(RuntestArgs, self).__init__(args, rest)
        self.keep_going = None
        self.valgrind = None
        self.parallel = None
        self.batch = None
        self.include = None
        self.exclude = None

        self.cpu = None
        self.queue = None
        self.host = None

        self.time_limit = None
        self.memory_limit = None

        self.root = None
        self.json = None
        self.list = None
        self.dump = None
        self.log = None

        self.export = None
        self.random_output_dir = None
        self.no_clean = None
        self.no_compare = None
        self.death_test = None
        self.status_file = None

        for k, v in args._get_kwargs():
            setattr(self, k, v)

        # fix more complex args
        self.queue = True if self.queue is None else self.queue
        self.valgrind = True if self.valgrind is None else self.valgrind
        self.random_output_dir = System.rnd8 if self.random_output_dir is None else self.random_output_dir
        self.cpu = [x for y in self.cpu for x in y] if self.cpu else self.cpu


class SmartFormatter(argparse.HelpFormatter):
    """
    Class SmartFormatter prints help messages without any formatting
        or unwanted line breaks
    """

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return [l.lstrip() for l in text[2:].lstrip().splitlines()]
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


class Parser(object):
    """
    Class Parser is static class which wraps not so well documented
    ArgumentParser.add_argument method
    """

    NONE = '-no-value-provided-'
    IGNORED = ('self', 'name', 'cls', 'parser')

    STORE_TRUE  = 'store_true'
    STORE_FALSE = 'store_false'
    STORE_CONST = 'store_const'
    APPEND = 'append'
    on_parse = Event()

    @classmethod
    def create(cls, *args, **kwargs):
        return argparse.ArgumentParser(*args, formatter_class=SmartFormatter, **kwargs)

    @classmethod
    def add(cls, parser, name, action=NONE, nargs=NONE,
            const=NONE, default=NONE, type=NONE, choices=NONE,
            required=NONE, help=NONE, metavar=NONE, dest=NONE):
        """
        :type parser: argparse.ArgumentParser
        :type dest: str
        :type metavar: str
        :type help: str
        :type required: bool
        :type choices: list
        :type type: class
        :type default: str|int|float|object
        :type const: str|int|float|object
        :type nargs: str|int
        :type action: str
        :type name: str
        """
        items = {k:v for k, v in locals().items() if v != cls.NONE and k not in cls.IGNORED}
        names = [x.strip() for x in name.split(',')]
        names.reverse()
        parser.add_argument(*names, **items)

    @classmethod
    def parse_runtest(cls, parser, args=None):
        """
        Will perform standard parse_args call and
        return RuntestArgs instance
        :type parser: argparse.ArgumentParser
        :rtype: RuntestArgs
        """
        result = RuntestArgs(*parser.parse_known_args(args))
        if not result.args:
            parser.error('No yaml files or folder given')

        cls.on_parse(result)
        return result

    @classmethod
    def parse_exec_with_limit(cls, parser, args=None):
        """
        Will perform standard parse_args call and
        return ExecWithLimitArgs instance
        :type parser: argparse.ArgumentParser
        :rtype: ExecWithLimitArgs
        """
        result = ExecWithLimitArgs(*parser.parse_known_args(args))
        if not result.memory_limit and not result.time_limit:
            parser.error('At least one limit required (either time limit or memory limit)')

        if not result.rest:
            parser.error('Command to execute is missing, specify command you want to run after -- ')

        cls.on_parse(result)
        return result

    @classmethod
    def parse_exec_parallel(cls, parser, args=None):
        """
        Will perform standard parse_args call and
        return ExecParallelArgs instance
        :type parser: argparse.ArgumentParser
        :rtype: ExecParallelArgs
        """
        result = ExecParallelArgs(*parser.parse_known_args(args))
        if not result.rest:
            parser.error('Parallel command and binary command are missing, specify them after -- section')

        if len(result.rest) < 2:
            parser.error('Command to execute is missing, specify command you want to run after -- ')

        cls.on_parse(result)
        return result

#
# import sys
# import os
# import re
# from scripts.core.base import Printer
# from scripts.core.exceptions import ArgumentException
# from utils.events import Event
#
# _long_format = re.compile(r'--[a-z0-9_-]+=')
# _short_eq_format = re.compile(r'-[a-z]=')
# _short_colon_format = re.compile(r'-[a-z]:')
# _short_nospace_format = re.compile(r'-[a-z].')
#
# # python list of int or floats (also empty list): [1, 2, 3.5]
# _list_format = re.compile(r'^\[(\d+\.?\d*)(\s*,\s*(\d+\.?\d*))*\]|\[\]')
# _list_format_convert = eval
#
# # python list of int or floats (also empty list) WITHOUT braces: 1,2,3.5
# _list_nobrace_format = re.compile(r'^(\d+\.?\d*)(\s*,\s*(\d+\.?\d*))*$')
# _list_nobrace_format_convert = lambda x: eval('[{}]'.format(x))
#
# # python list of int or floats WITHOUT commas: [1 2 3.5]
# _list_space_format = re.compile(r'^\[(\d+\.?\d*)(\s+(\d+\.?\d*))*\]|\[\]')
# _list_space_format_convert = lambda x: eval(re.sub('\s+', ',', x))
#
# # python single int or float: 1 or 2.5
# _list_single_digit = re.compile(r'^\d+\.?\d*$')
# _list_single_digit_convert = lambda x: [eval(x)]
#
# # range format such as 1:3
# _list_range_short = re.compile(r'^(\d+):(\d+)$')
#
#
# def _list_range_short_convert(x):
#     """
#     Converts input to string from short list format which is defined as
#     from:to
#     :param x:
#     """
#     args = [int(y) for y in _list_range_short.match(x).groups()]
#     return list(range(args[0], args[1] + 1))
#
# # range format such as 1:10:2
# _list_range_long = re.compile(r'^(\d+):(\d+):(\d+)$')
#
#
# def _list_range_long_convert(x):
#     """
#     Converts input to string from long list format which is defined as
#     from:to:step
#     :param x:
#     """
#     args = [int(y) for y in _list_range_long.match(x).groups()]
#     return list(range(args[0], args[1] + 1, args[2]))
#
# # all list supported format
# _list_formats = [
#     [_list_format, _list_format_convert],
#     [_list_space_format, _list_space_format_convert],
#     [_list_single_digit, _list_single_digit_convert],
#     [_list_range_short, _list_range_short_convert],
#     [_list_range_long, _list_range_long_convert],
#     [_list_nobrace_format, _list_nobrace_format_convert]
# ]
#
# # regex for getting arg name
# _parse_arg_name = re.compile(r'^(--[a-zA-Z0-9_-]+|-[a-zA-Z0-9_-])')
#
#
# class ArgOption(object):
#     """
#     Class ArgOption is simple container for single argument option
#     """
#
#     def __init__(self, short, long, type=str, default=None, name=None, subtype=str, docs='', placeholder=None, hidden=False):
#         self.short = short
#         self.long = long
#         self.type = type
#         self.subtype = subtype
#         self.default = default
#         self.docs = docs
#         self.hidden = hidden
#         self.name = name or self.long[2:] or self.short[1:]
#         self.placeholder = placeholder or self.name
#
#         self.reset()
#
#     def reset(self):
#         if self.type is True or self.type is False:
#             self.value = not self.type
#         elif self.type is list:
#             self.value = list()
#         else:
#             self.value = self.default
#
#     def is_primitive(self):
#         return self.type in (True, False)
#
#     def parse_list(self, value):
#         for fmt, conv in _list_formats:
#             if fmt.match(value):
#                 try:
#                     lst = conv(value)
#                     if type(lst) is not list:
#                         raise Exception('Invalid format {}'.format(value))
#                     lst = [self.subtype(x) for x in lst]
#                     return lst
#                 except:
#                     raise Exception('Invalid format {}'.format(value))
#         # return values as list
#         return [value]
#
#     def usage(self):
#         lsn = ''
#         if self.long and self.short:
#             if self.is_primitive():
#                 lsn = '{self.short}, {self.long}'.format(self=self)
#             else:
#                 lsn = '{self.short}, {self.long} {name}'.format(self=self, name=self.placeholder.upper())
#         else:
#             value = self.short or self.long
#             if self.is_primitive():
#                 lsn = '{value}'.format(value=value)
#             else:
#                 lsn = '{value} {name}'.format(value=value, name=self.placeholder.upper())
#
#         lsn = '  {:30s}'.format(lsn)
#         blank = 32 * ' '
#         if self.docs:
#             if type(self.docs) is str:
#                 return '{} {self.docs}'.format(lsn, self=self)
#             else:
#                 result = ''
#                 for l in self.docs:
#                     result += '{} {}\n'.format(lsn, l)
#                     lsn = blank
#                 return result.rstrip()
#
#     def __repr__(self):
#         return '{self.name}: {self.value}'.format(self=self)
#
#
# class ArgOptions(dict):
#     """
#     Class ArgOptions is dictionary with dot access available
#     """
#
#     def __getattr__(self, attr):
#         return self.get(attr)
#
#     def __setattr__(self, key, value):
#         self.__setitem__(key, value)
#
#     def __setitem__(self, key, value):
#         super(ArgOptions, self).__setitem__(key, value)
#         self.__dict__.update({key: value})
#
#     def __delattr__(self, item):
#         self.__delitem__(item)
#
#     def __delitem__(self, key):
#         super(ArgOptions, self).__delitem__(key)
#         del self.__dict__[key]
#
#
# class ArgParser(object):
#     """
#     Class ArgParser is command-line parsing class.
#
#     Parser support short and long flags with custom type definition.
#     It also return rest of arguments located after double dash --
#     """
#
#     def __init__(self, usage):
#         self._args = [str(x) for x in sys.argv[1:]]
#         self.args = list()
#         self.options = ArgOptions()
#         self.options_map = {}
#         self.others = []
#         self.rest = []
#         self.source = None
#         self.i = None
#         self.keys = None
#         self._usage = usage
#         self.all_options = list()
#         self.add('-h', '--help', type=True, name='help', docs='Display this help and exit')
#         self.on_parse = Event()
#
#     def add(self, short='', long='', type=str, default=None, name=None, subtype=str, docs='', placeholder='', hidden=False):
#         ao = ArgOption(short, long, type, default, name, subtype, docs, placeholder, hidden)
#
#         self.all_options.append(ao)
#         if name:
#             self.options[name] = ao
#         if short:
#             self.options_map[short] = ao
#         if long:
#             self.options_map[long] = ao
#
#     def add_section(self, name):
#         self.all_options.append('\n{}:'.format(name))
#
#     def usage(self):
#         usage_lst = ['Usage: {}'.format(self._usage)]
#         for option in self.all_options:
#             if type(option) is str:
#                 usage_lst.append('{option}\n'.format(option=option))
#             elif not option.hidden:
#                 usage_lst.append('{option}\n'.format(option=option.usage()))
#         return '\n'.join(usage_lst)
#
#     def check_help(self):
#         if self.simple_options.get('help'):
#             self.exit_usage(exit_code=0)
#
#     def exit_usage(self, msg=None, exit_code=1, *args, **kwargs):
#         if msg:
#             Printer.all.raw('Error: {}'.format(msg), *args, **kwargs)
#
#         Printer.all.raw(self.usage())
#
#         if exit_code is not None:
#             raise ArgumentException(exit_code, msg)
#
#     def current(self):
#         """
#         :rtype : str
#         """
#         return self.source[self.i]
#
#     def next(self):
#         """
#         :rtype : str
#         """
#         return None if self.i + 1 >= len(self.source) else self.source[self.i + 1]
#
#     def move_on(self):
#         self.i += 1
#
#     def process_option(self, option):
#         """
#         :type option: ArgOption
#         """
#         if option.type is True:
#             option.value = True
#         elif option.type is False:
#             option.value = False
#         elif option.type in (int, float, long, str):
#             option.value = option.type(self.next())
#             self.move_on()
#         elif option.type is list:
#             option.value.extend(option.parse_list(self.next()))
#             self.move_on()
#         elif type(option.type) is list:
#             # if next arg is -- or if next arg is registered set value to True
#             # otherwise save next value
#             if not self.next() or (self.next() == '--' or self.next_is_registered()):
#                 option.value = option.type[0]
#             else:
#                 option.value = self.next()
#                 self.move_on()
#         else:
#             option.value = option.type(self.next())
#             self.move_on()
#
#         return option.value
#
#     def next_is_registered(self):
#         arg = self.next()
#         match = _parse_arg_name.match(arg)
#         if match:
#             return match.group(1) in self.options_map
#
#     def split_current(self):
#         arg = self.current()
#         result = list()
#         # double dash format
#         if arg.startswith('--'):
#             # format is --NAME=VALUE
#             if _long_format.match(arg):
#                 result.extend(arg.split('=', 1))
#             # format does not contain = sign
#             else:
#                 raise Exception("Invalid input format {}".format(arg))
#
#         # single dash format
#         elif arg.startswith('-'):
#             # format is -f=VALUE
#             if _short_eq_format.match(arg):
#                 result.extend(arg.split('=', 1))
#             # format is -f:VALUE
#             elif _short_colon_format.match(arg):
#                 result.extend(arg.split(':', 1))
#             elif _short_nospace_format.match(arg):
#                 result.append(arg[0:2])
#                 result.append(arg[2:])
#             else:
#                 raise Exception("Invalid input format {}".format(arg))
#         else:
#             raise Exception("Invalid input format {}".format(arg))
#
#         # extend source arg list
#         pre = self.source[0:self.i]
#         curr = result
#         post = self.source[self.i + 1:]
#         self.source = pre + curr + post
#
#     @property
#     def simple_options(self):
#         simple_options = ArgOptions()
#         for k, v in self.options.items():
#             simple_options[k] = v.value
#         return simple_options
#
#     def parse(self, args=None):
#         self.args = []
#         self.i = 0
#         self.keys = sorted(self.options_map.keys(), reverse=True)
#         self.source = args if args is not None else self._args
#         self.others = []
#         self.rest = []
#
#         for opt in self.all_options:
#             if type(opt) is not str:
#                 opt.reset()
#
#         if not self.source:
#             self.exit_usage()
#
#         while self.i < len(self.source):
#             find = False
#             if self.current() in self.options_map:
#                 find = True
#                 self.process_option(self.options_map[self.current()])
#             else:
#                 for k in self.keys:
#                     if k.startswith('--'):
#                         if self.current().startswith(k + "=") or self.current().startswith(k + ":"):
#                             self.split_current()
#                             self.process_option(self.options_map[self.current()])
#                             find = True
#                             break
#                     elif k.startswith('-'):
#                         if self.current().startswith(k):
#                             self.split_current()
#                             self.process_option(self.options_map[self.current()])
#                             find = True
#                             break
#
#             # add to others if not found
#             if not find:
#                 # end of parsing section
#                 if self.current() == '--':
#                     self.rest = self.source[self.i + 1:]
#                     self.check_help()
#                     self.on_parse(self)
#                     return self.simple_options, self.others, self.rest
#                 # just add to others
#                 else:
#                     self.others.append(self.current())
#             self.i += 1
#
#         self.check_help()
#         self.on_parse(self)
#         return self.simple_options, self.others, self.rest
#
#     def __getattr__(self, item):
#         for i in self.options.values():
#             if item == i.name:
#                 return i.value
#
#
# # args = ['-p', '5', '456', '--fff', '45', '-l4456', '-j=6', '--foo-56=789', '--foo=456=9']
# # args = ['-p', '-foo', '--foo', '456', '--true']
# # args = ['--foo=56', '-k', 'cas', '--fo=789', '--', 'foo bar', '-f', 'vsdei']
# # args = ['--ll', '1:5:2', '--', 'fooooooo']
# # ap = ArgParser()
# # # ap.add(long='--foo', tp=int)
# # # ap.add(long='--fo', tp=int)
# # # ap.add('-l', tp=False)
# # # ap.add('-l', tp=False)
# # ap.add(long='--ll', type=list, subtype=int, name="foo")
# # options, others, rest = ap.parse(args)
# # print options.foo
