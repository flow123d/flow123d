#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import os
import re
import json
import hashlib
import datetime
import uuid

from py123d.scripts.artifacts.collect.modules import ICollectTool, CollectResult, LogPolicy


class Flow123dProfiler(ICollectTool):
    """
    Class Flow123dProfiler is module for project Flow123d
    This module processes json profiler data and runtest.status.json files
    """

    include = 'profiler_info_*.log.json'
    exclude = None

    _floats = [
        re.compile('cumul-time-.+'),
        re.compile('percent'),
        re.compile('timer-resolution'),
    ]

    _ints = [
        re.compile('call-count-.+'),
        re.compile('memory-.+'),
        re.compile('file-line'),
        re.compile('task-size'),
        re.compile('run-process-count'),
    ]

    _dates = [
        re.compile('run-started-at'),
        re.compile('run-finished-at'),
    ]

    _children = 'children'

    @staticmethod
    def _parse_date(s:str) -> datetime.datetime:
        return datetime.datetime.strptime(s, '%m/%d/%y %H:%M:%S').timestamp()

    def try2extract_from_dummy(self, d: str, status):
        pcs = d.split('.')
        if len(pcs) == 3:
            name, cpu, el = pcs
            pcs = el.split('_')
            el = int(pcs[1]) if len(pcs) == 3 else int(pcs[0])
        elif len(pcs) == 2:
            name, cpu = pcs
            el = 0
        else:
            name = d
            cpu = 0
            el = 0
        return {
            'task-size': el,
            'run-process-count': int(cpu),
            'children': [
                {
                    'tag': 'Whole Program'
                }
            ]
        }

    def process_file(self, f: str) -> CollectResult:
        """
        Method processes single file and return collect result
        :param f: file location
        """
        with open(f, 'r') as fp:
            obj = json.load(fp)

        status_file = os.path.join(os.path.dirname(f), 'runtest.status.json')
        status = {}
        if os.path.exists(status_file):
            with open(status_file, 'r') as fp:
                status = json.load(fp)

        # is this just dummy profiler:
        if not obj:
            obj = self.try2extract_from_dummy(os.path.basename(os.path.dirname(f)), status)

        # convert unix timestamp to datetime class to mongo interprets is correctly
        if 'commit' in status and status['commit'] and 'date' in status['commit']:
            status['commit']['date'] = datetime.datetime.fromtimestamp(status['commit']['date'])

        # convert fields to ints and floats
        self._convert_fields(obj, self._ints,   int)
        self._convert_fields(obj, self._floats, float)
        self._convert_fields(obj, self._dates, self._parse_date)

        # unwind
        parts = f.split('/')
        base, start = self._get_base(obj)
        base['test-name'] = parts[-4]
        base['case-name'] = parts[-2].split('.')[0]
        base['uuid'] = uuid.uuid4().hex
        base.update(status)

        items = self._unwind(start, list(), base)
        item = CollectResult(
            items=items,
            log_policy=LogPolicy.ON_ERROR,
            log_folder=os.path.dirname(f),
            logs=[
                os.path.join(os.path.dirname(f), 'job_output.log'),
                os.path.join(os.path.dirname(f), 'flow123.0.log')
            ],
        )
        return item

    def _get_base(self, obj: dict) -> (dict, dict):
        if not obj or self._children not in obj:
            return {}, {}

        base = obj.copy()
        if self._children in base:
            del base[self._children]
            return base, obj[self._children][0]

    def _convert_fields(self, obj, fields, method):
        for key in obj:
            for f in fields:
                if f.match(key):
                    obj[key] = method(obj[key])

        if self._children in obj:
            for child in obj[self._children]:
                self._convert_fields(child, fields, method)

    def _unwind(self, obj: dict, result: list, base: dict = None, path: str = ''):
        # empty obj
        if not obj:
            return result
            
        item = obj.copy()
        if self._children in item:
            del item[self._children]

        tag = item['tag']
        new_path = path + '/' + tag
        indices = dict()
        indices['tag'] = tag
        indices['tag_hash'] = self._md5(tag)
        indices['path'] = new_path
        indices['path_hash'] = self._md5(new_path)
        indices['parent'] = path
        indices['parent_hash'] = self._md5(path)
        item['indices'] = indices
        item['base'] = base.copy()

        result.append(item)

        if self._children in obj:
            for o in obj[self._children]:
                self._unwind(o, result, base, new_path)
        return result

    @classmethod
    def _md5(cls, s):
        return hashlib.md5(s.encode()).hexdigest()
