#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import os
import sys
import argparse
import time
import subprocess
import signal

default_location = os.path.join(os.getcwd(), '.pbs')
default_location = '/tmp/.pbs'

parser = argparse.ArgumentParser()
parser.add_argument('action', help='Specify action (start, qsub, qstat, qdel)')
parser.add_argument('-d', '--dir', dest='location', default=default_location, help='PBS folder')
parser.add_argument('-o', '--output', default='-', help='Script output location')


class Qstat(object):
    def __init__(self, filename):
        self.filename = filename
        content = open(self.filename, 'r').read().strip()
        self.status, self.script, self.output = content.split()
        self.cid = os.path.basename(self.filename).split('.')[1]

    def execute(self):
        if self.output == '-':
            os.chmod(self.script, 0o777)
            process = subprocess.Popen([self.script])
            self.set_status('r')
            process.wait()
            self.set_status('f')
        else:
            with open(self.output, 'w') as fp:
                print(open(self.script, 'r').read())
                os.chmod(self.script, 0o777)
                process = subprocess.Popen([self.script], stderr=subprocess.STDOUT, stdout=fp)
                self.set_status('r')
                process.wait()
                self.set_status('f')

    def set_status(self, status):
        with open(self.filename, 'w') as fp:
            fp.write(status)
            fp.write('\n')
            fp.write(self.script)
            fp.write('\n')
            fp.write(self.output)
        self.status = status

    def __repr__(self):
        return '{}:{} ({})'.format(self.cid, self.status, self.script)


class PBS(object):
    cid_file = None
    args = None

    @classmethod
    def get_qstats(cls, loc):
        return [Qstat(os.path.join(loc, x)) for x in sorted(os.listdir(loc)) if str(x).startswith('qstat')]

    @classmethod
    def next(cls, value=None):
        next_cid = cls.cid() + 1 if value is None else value
        with open(cls.cid_file, 'w+') as fp:
            fp.write(str(next_cid))
            fp.flush()
            fp.close()
        return next_cid

    @classmethod
    def init(cls, args):
        cls.args = args
        cls.cid_file = os.path.join(cls.args.location, '.cid')

    @classmethod
    def cid(cls):
        return int(open(cls.cid_file, 'r+').read().strip())

    @classmethod
    def start(cls, *args):
        print('Starting server in', cls.args.location)
        cls.next(0)

        def signal_handler(signal, frame):
            import shutil
            shutil.rmtree(cls.args.location)
            sys.exit(0)

        signal.signal(signal.SIGINT, signal_handler)

        while True:
            qstats = cls.get_qstats(cls.args.location)
            queued = [x for x in qstats if x.status == 'q']
            for q in shuffled(queued):
                print(q)
                q.execute()
                print(q)
            time.sleep(0.5)

    @classmethod
    def qstat(cls, *args):
        print('{:8s} {:8s} {:8s}'.format('ID', 'status', 'queue'))

        for q in cls.get_qstats(cls.args.location):
            print ('{:8s} {:8s} {:8s}'.format(str(q.cid), q.status, 'default'))

    @classmethod
    def qsub(cls, *args):
        filename = os.path.abspath(args[0])
        cid = cls.cid()
        with open(os.path.join(cls.args.location, 'qstat.'+str(cid)), 'w') as fp:
            fp.write('q')
            fp.write('\n')
            fp.write(filename)
            fp.write('\n')
            fp.write(cls.args.output)
        cls.next()
        print('Job inserted,', cid, 'into queue')


def shuffled(array):
    import random
    c = array.copy()
    random.shuffle(c)
    return c

if __name__ == '__main__':
    args, rest = parser.parse_known_args()

    # create folder if necessary
    if not os.path.exists(args.location):
        os.makedirs(args.location)

    PBS.init(args)

    getattr(PBS, args.action)(*rest)