#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from optparse import OptionParser


def get_replacement_key(l, placeholders):
    """
    Return replacement if present
    :param l:
    :param placeholders:
    :return:
    """
    ll = l[1:].strip()
    if l.startswith('#') and ll in placeholders:
        return ll


def replace_multiline(lines, placeholders):
    """
    Replace multiple line placeholder
    :param lines:
    :param placeholders:
    :return:
    """
    result = list()
    opened = False
    opened_content = list()

    for line in lines:
        repl_i = get_replacement_key(line, placeholders)
        if repl_i is not None:
            if opened:
                for repl in placeholders[repl_i]:
                    for ll in opened_content:
                        result.append(ll.replace(repl_i, repl))
                opened = False
            else:
                opened = True
                opened_content = list()
            continue

        if opened:
            opened_content.append(line)
        else:
            result.append(line)
    return result


def replace_singleline(lines, placeholders):
    """
    Replace single line placeholder
    :param lines:
    :param placeholders:
    :return:
    """
    result = list()
    for line in lines:
        found = False
        for key in placeholders.keys():
            if line.find(key) != -1:
                found = True
                for repl in placeholders[key]:
                    result.append(line.replace(key, repl))

        if not found:
            result.append(line)
    return result


def load(f):
    """
    Loads a file
    :param f:
    :return:
    """
    with open(f, 'r') as fp:
        return fp.read().splitlines()


parser = OptionParser('jenkins_jobs_script SOURCE [[PLACEHOLDER PLACEHOLDER_FILE], ... ]')
parser.add_option('--prefix', default='$PLACEHOLDER_', dest="prefix")
parser.add_option('--suffix', default='$', dest="suffix")

values, args = parser.parse_args()

# get source file
source = args[0]
args = args[1:]

placeholders = {}
for i in range(0, len(args), 2):
    placeholder_name = values.prefix + args[i] + values.suffix
    placeholder_values = [x.strip() for x in load(args[i + 1]) if not x.startswith('#') and x.strip()]

    placeholders[placeholder_name] = placeholder_values


lines = load(source)
lines = replace_multiline(lines, placeholders)
lines = replace_singleline(lines, placeholders)

print '\n'.join(lines)
