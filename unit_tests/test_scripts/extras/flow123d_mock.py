#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# Flow123d mock file
# ----------------------------------------------
import os
import random
import shutil
from optparse import OptionParser
# ----------------------------------------------
import time
import sys
# ----------------------------------------------

parser = OptionParser()
parser.add_option('-e', '--error',  dest="error",  help="Raise error", action='store_true')
parser.add_option('-s', '--solve',  dest="solve",  help="Solve file")
parser.add_option('-t', '--time',   dest="time",   help="Computation duration, default %default s", default=0.0001)
parser.add_option('-i', '--input',  dest="input",  help="output file")
parser.add_option('-o', '--output', dest="output", help="input file")
parser.add_option('--copy',   dest="copy",   help="Copy reference", action='store_true')
parser.add_option('--clean',  dest="clean",  help="Clean output", action='store_true')
parser.add_option('--random', dest="random", help="Randomize output", action='store_true')

options, args = parser.parse_args()


def randomize(top):
    for root, dirs, files in os.walk(top):
        for name in files:
            f = os.path.join(root, name)
            with open(f, 'w') as fp:
                fp.write(str(random.random()))
                fp.write(' ')
                fp.write(str(10000 * random.randint(1, 200)))
                fp.write(' ')


time.sleep(float(options.time))
if options.solve:
    solve_name = os.path.basename(options.solve)
    solve_dir = os.path.dirname(options.solve)
    ref = os.path.join(solve_dir, 'ref_out', solve_name[:solve_name.find('.')])
    out = os.path.join(solve_dir, 'test_results', solve_name[:solve_name.find('.')])

    # clean on demand
    if options.clean:
        [shutil.rmtree(out + '.' + str(i), ignore_errors=True) for i in range(1, 4)]

    # copy ref on demand
    if options.copy:
        [shutil.copytree(ref, out + '.' + str(i)) for i in range(1, 4)]

    # copy random on demand
    if options.random:
        shutil.copytree(ref, out + '.1')
        randomize(out + '.1')
        [shutil.copytree(out + '.1', out + '.' + str(i)) for i in range(2, 4)]


if options.error:
    raise Exception('Raising exception on demand')
