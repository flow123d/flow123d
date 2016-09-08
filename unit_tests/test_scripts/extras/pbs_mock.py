#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# Flow123d mock file
# ----------------------------------------------
import random
import subprocess
import time
import os
# ----------------------------------------------
from optparse import OptionParser
# ----------------------------------------------

prefix = '.qstat_'
parser = OptionParser()
parser.add_option('-o', dest='output')
options, args = parser.parse_args()

actions = ['qsub', 'qstat']
if len(args) == 0:
    raise Exception('Action required')

action = args[0]
if action not in actions:
    raise Exception("Invalid action, action '{}' not in {}".format(action, actions))


def get_qstat_files():
    return [f for f in os.listdir('.') if f.startswith(prefix)]


def get_qstat_status(f):
    qid = int(f[len(prefix):])
    with open(f, 'r') as fp:
        status = fp.read().strip()
    return qid, status


def get_next_id():
    maximum = 0
    for f in get_qstat_files():
        qid, status = get_qstat_status(f)
        if qid > maximum:
            maximum = qid
    return maximum + 1


def write_all(f, c):
    with open(f, 'w') as fp:
        fp.write(c)


def read_all(f):
    with open(f, 'r') as fp:
        return fp.read()


def qsub(pbs_file, queue=0, delay=1, rnd=0, output=None):
    qid = get_next_id()
    pbs_file = os.path.abspath(pbs_file)

    # insert into queue
    output = os.path.abspath(output) if output else output
    qstat_file = os.path.abspath(prefix + str(qid))
    time.sleep(0.001 + queue + random.random() * rnd)
    write_all(qstat_file, 'q')
    print ('Job inserted into queue (ID={qid}, queue=default)'.format(**locals()))

    # run
    time.sleep(0.001 + delay + random.random() * rnd)
    write_all(qstat_file, 'r')
    pbs_file_copy = pbs_file + '.copy'
    content = read_all(pbs_file)
    content += "\n\necho c > {}\n\n".format(qstat_file)
    write_all(pbs_file_copy, content)

    if output:
        command = '/bin/bash {pbs_file_copy} > {output} &'.format(**locals())
        process = subprocess.Popen(command, shell=True)
    else:
        command = '/bin/bash {pbs_file_copy}  &'.format(**locals())
        process = subprocess.Popen(command, shell=True)


def qstat(*args):
    print ('{:8s} {:8s} {:8s}'.format('ID', 'status', 'queue'))
    for f in get_qstat_files():
        qid, status = get_qstat_status(f)
        print ('{:8s} {:8s} {:8s}'.format(str(qid), status, 'default'))


if action == 'qstat':
    qstat(*args[1:])

if action == 'qsub':
    qsub(*args[1:], output=options.output)
