#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import json
import threading
import formic
import getpass
import platform

from scripts.artifacts.artifacts import ArtifactStep
from scripts.core.base import Paths, Printer
from utils.timer import Timer

TIMEOUT = 5

profiler_int_fields = (
    'call-count-max', 'call-count-min', 'call-count-sum',
    'memory-alloc-min', 'memory-alloc-max', 'memory-alloc-sum',
    'memory-dealloc-min', 'memory-dealloc-max', 'memory-dealloc-sum',
    'memory-peak-min', 'memory-peak-max', 'memory-peak-sum',
    'memory-alloc-called-min', 'memory-alloc-called-max', 'memory-alloc-called-sum',
    'memory-dealloc-called-min', 'memory-dealloc-called-max', 'memory-dealloc-called-sum',
    'memory-petsc-diff-min', 'memory-petsc-diff-max', 'memory-petsc-diff-sum',
    'memory-petsc-peak-min', 'memory-petsc-peak-max', 'memory-petsc-peak-sum',
)
profiler_float_fields = (
    'cumul-time-max',
    'cumul-time-min',
    'cumul-time-sum',
    'percent'
)

class DatabaseMongo(ArtifactStep):
    """
    Class DatabaseMongo collects json files from
    source directory, and publish them into remote mongo db
    """

    yaml_tag = u'!Database.mongodb'

    def __init__(self, server, username, password, source):
        self.server = server
        self.username = username
        self.password = password
        self.source = source

        from pymongo import MongoClient
        self.client = MongoClient(self.server,
                                  connectTimeoutMS=TIMEOUT*1000,
                                  socketTimeoutMS=TIMEOUT*1000,
                                  serverSelectionTimeoutMS=TIMEOUT*1000,
                                  connect=False
                                  )
        self.db = self.client.get_database('bench')
        self.connected = False

    def auth(self):
        def auth_raw():
            try:
                self.connected = self.db.authenticate(
                    name=self.username,
                    password=self.password,
                    source='admin'
                )
                return True
            except Exception as e:
                Printer.all.err(str(e))
                return False

        thread = threading.Thread(target=auth_raw)
        # thread.setDaemon(True)
        thread.start()
        thread.join(TIMEOUT)
        return self.connected

    def run(self):
        # check connection to db first
        if not self.auth():
            Printer.all.err("Could not connect to db")
            return

        profilers = list(formic.FileSet("profiler_info_*.log.json", directory=self.source))
        benchmark = list(formic.FileSet("benchmark*.json", directory=self.source))
        lscpu = list(formic.FileSet("lscpu*.json", directory=self.source))

        nodes_result = None
        if len(benchmark) >= 1 and len(lscpu) >= 1:
            benchmark_data = self._read_json(benchmark[0])
            lscpu_data = self._read_json(lscpu[0])
            lscpu_data['benchmark'] = benchmark_data
            lscpu_data['username'] = getpass.getuser()
            lscpu_data['hostname'] = platform.node()
            lscpu_data['duration'] = Timer.app_timer.duration

            # insert
            nodes_result = self.db.get_collection('nodes').insert_one(lscpu_data)
            Printer.all.out("Inserted 1 document into nodes collection, {}", nodes_result.acknowledged)

        documents = list()
        for filename in profilers:
            info = self._extract_parts(filename)

            # lint to node if present
            if nodes_result:
                info['node_id'] = nodes_result.inserted_id

            # read profiler json and convert fields
            profiler = self._read_json(filename, conversions=[
                (profiler_int_fields, int),
                (profiler_float_fields, float)
            ])

            # merge documents and append to list
            profiler.update(info)
            documents.append(profiler)

        # insert to db
        profilers_result = self.db.get_collection('bench').insert_many(documents)
        Printer.all.out("Inserted {} document(s) into bench collection, {}", len(documents), profilers_result.acknowledged)

    @staticmethod
    def _extract_parts(filename):
        parts = Paths.split(filename)
        test = parts[-3]
        case = parts[-2].split('.')[0]
        nproc = int(parts[-2].split('.')[1])

        return {
            "test-name": test,
            "case-name": case,
            "nproc": nproc
        }

    @classmethod
    def _read_json(cls, filename, conversions=None):
        """
        :rtype: dict
        """
        # read json
        with open(filename, 'r') as fp:
            json_data = json.load(fp)

        # convert fields on demand
        if conversions:
            for fields, method in conversions:
                cls._convert_fields(json_data, method, *fields)

        return json_data

    @classmethod
    def _convert_fields(cls, o, method, *fields):
        """
        :type o: dict
        """
        for k, v in o.items():
            if k in fields:
                o[k] = method(v)

            if type(v) is list:
                for i in v:
                    if type(i) is dict:
                        cls._convert_fields(i, method, *fields)
            if type(v) is dict:
                cls._convert_fields(v, method, *fields)
