#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from optparse import OptionParser
import sys, os, re, json

usage="""
Simple script for getting information about python paths and python in general
$> python_path_script ACTION [--format FORMAT]
ACTION must be one of the following:
    sys_prefix                  - python prefix
    sys_exec_prefix             - python exec prefix
    sys_path                    - python path list
    sys_path_without_prefix     - python path list without sys.prefix
    sys_path_subdir             - list of python path first folder from python path list with removed python prefix
FORMAT must be one of the following:
    cmake                       - colon separated joined lists
    json                        - valid json format
    bash                        - '\\n' new line joined format
"""


def get_formatter(f):
    if f.lower() == 'cmake':
        return lambda x: ';'.join(x)
    if f.lower() == 'json':
        return lambda x: json.dumps(x)
    if f.lower() in ('human', 'bash'):
        return lambda x: '\n'.join(x)
    raise "Unsupported formatter " + str(f)


def action_sys_prefix():
    return [sys.prefix]


def action_sys_exec_prefix():
    return [sys.exec_prefix]


def action_sys_path():
    return [p for p in sys.path if p not in ('', os.getcwd())]


def action_sys_path_without_prefix():
    paths = [p for p in sys.path if p not in ('', os.getcwd())]
    paths = [x.replace(sys.prefix, '') for x in paths if x.find(sys.prefix) != -1]
    return paths


def action_sys_path_subdir():
    # remove current dir and loop thru only items containing sys.prefix
    paths = [p for p in sys.path if p not in ('', os.getcwd())]
    paths = [x.replace(sys.prefix, '') for x in paths if x.find(sys.prefix) != -1]
    sub_paths = []
    # store first valid part of each path
    for p in paths:
        if p:
            for x in re.split(r'\\|/', p):
                if x:
                    sub_paths.append(x)
                    break
    # convert result to set to remove duplicates and back to list again
    return list(set(sub_paths))


parser = OptionParser()
parser.usage = usage
parser.add_option("-f", "--format", dest="format", default="cmake",
                  help="output format of the action (cmake, json)", metavar="FORMAT")
options, args = parser.parse_args()
action = args[0] if len(args) else None

actions = dict(
    sys_prefix=action_sys_prefix,
    sys_exec_prefix=action_sys_exec_prefix,
    sys_path=action_sys_path,
    sys_path_without_prefix=action_sys_path_without_prefix,
    sys_path_subdir=action_sys_path_subdir
)

if action is None or actions.get(action, None) is None:
    parser.print_help()
    exit(1)

print get_formatter(options.format)(actions.get(action)())
