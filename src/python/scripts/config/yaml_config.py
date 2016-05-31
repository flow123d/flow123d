#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from __future__ import absolute_import

import yaml
import copy
import itertools

from scripts.core.base import Paths
from utils.globals import ensure_iterable


default_values = dict(
    proc=[1],
    time_limit=30,
    memory_limit=400,
    tags=[],
    check_rules=[
        {
            'ndiff': {
                'files': ['*']
            }
        }
    ]
)


class ConfigCaseBase(object):
    def __init__(self, config):
        """
        :type config: scripts.config.yaml_config.YamlConfig
        """
        self.config = config

        self.proc = default_values.get('proc')
        self.time_limit = default_values.get('time_limit')
        self.memory_limit = default_values.get('memory_limit')
        self.check_rules = default_values.get('check_rules')
        self.tags = set(default_values.get('tags'))


    @classmethod
    def _get(cls, o, prop):
        return o.get(prop, default_values.get(prop))

    def __repr__(self):
        return '<{self.__class__.__name__} {self.files}>'.format(self=self)


class DummyConfigCase(ConfigCaseBase):
    def __init__(self, config, yaml_file):
        """
        :type config: scripts.config.yaml_config.YamlConfig
        """
        super(DummyConfigCase, self).__init__(config)
        self.files = ensure_iterable(yaml_file)
        self.proc = [1]


class YamlConfigCase(ConfigCaseBase):
    def __init__(self,  config, o={}):
        """
        :type config: scripts.config.yaml_config.YamlConfig
        """
        super(YamlConfigCase, self).__init__(config)
        self.proc = self._get(o, 'proc')
        self.time_limit = self._get(o, 'time_limit')
        self.memory_limit = self._get(o, 'memory_limit')
        self.check_rules = self._get(o, 'check_rules')
        self.files = ensure_iterable(self._get(o, 'file'))
        self.tags = set(self._get(o, 'tags'))

        for i in range(len(self.files)):
            self.files[i] = Paths.join(config.root, self.files[i])


class YamlConfig(object):
    """
    :type test_cases: list[scripts.config.yaml_config.YamlConfigCase]
    :type common_config: dict
    """
    def __init__(self, filename):
        # prepare paths
        self.filename = filename
        self.root = Paths.dirname(self.filename)
        self.test_results = Paths.join(self.root, 'test_results')
        self.ref_output = Paths.join(self.root, 'ref_output')
        self.input = Paths.join(self.root, 'input')
        self.test_cases = list()

        # read config or use default mini config
        if Paths.exists(self.filename):
            with open(self.filename, 'r') as fp:
                self._yaml = yaml.load(fp)
        else:
            self._yaml = dict(
                common_config=default_values.copy()
            )
        self._iter_index = 0
        self.common_config = None
        self.test_cases = None
        self.include = []
        self.exclude = []

    def update(self, **kwargs):
        """
        Updates all test_case values
        :param kwargs:
        :return:
        """
        for k,v in kwargs.items():
            if v:
                for test_case in self.test_cases:
                    setattr(test_case, k, v)

    def parse(self):
        # update common config using global values
        self.common_config = self._get(('common_config', 'commons'), {})
        self.common_config = self.merge(default_values, self.common_config)

        # update test_cases using common config values
        self.test_cases = list()
        test_cases = self._get('test_cases', [])
        for test_case in test_cases:
            test_case = self.merge(self.common_config, test_case)
            if self.check_tags(test_case):
                self.test_cases.append(YamlConfigCase(self, test_case))

    def check_tags(self, test_case):
        tags = set(test_case['tags'])
        inn = set(self.include)
        exc = set(self.exclude)
        result = True

        if inn:
            # if intersection between tags and include tags is empty
            # do not include this case
            if not tags.intersection(inn):
                return False
        if exc:
            # if intersection between tags and exclude tags is not empty
            # do not include this case
            if tags.intersection(exc):
                return False

        return result

    def get(self, index):
        """
        :rtype : scripts.config.yaml_config.YamlConfigCase
        """
        return self.test_cases[index]

    def _get(self, names, default=None):
        if type(names) in (list, tuple):
            for name in names:
                result = self._yaml.get(name, None)
                if result:
                    return result
            return default
        else:
            return self._yaml.get(names, default)

    def get_cases_for_file(self, prescription_class, yaml_file):
        """
        :type yaml_file: str
        :type prescription_class: class
        :rtype : list[scripts.core.prescriptions.PBSModule]
        """
        tmp_result = list()
        # prepare product of all possible combinations of input arguments for specified file
        # for now we use only proc (cpu list) and files (file)
        for test_case in self.test_cases:
            if yaml_file in test_case.files:
                tmp_result.append(list(itertools.product(
                    ensure_iterable(test_case),
                    ensure_iterable(test_case.proc),
                    ensure_iterable(test_case.files),
                )))

        # if no results exists for this particular config
        # we add dummy case which is basically default values
        if not tmp_result:
            dummy_case = DummyConfigCase(self, yaml_file)
            tmp_result.append(list(itertools.product(
                ensure_iterable(dummy_case),
                ensure_iterable(dummy_case.proc),
                ensure_iterable(dummy_case.files)
            )))

        result = list()
        for lst in tmp_result:
            result.extend([prescription_class(*x) for x in lst])

        # no config was given yaml file was declared
        if not result:
            result.append(prescription_class(default_values, 1, yaml_file))
        return result

    def get_all_cases(self, prescription_class):
        """
        :type prescription_class: class
        :rtype : list[scripts.core.prescriptions.MPIPrescription] or list[scripts.core.prescriptions.PBSModule]
        """
        tmp_result = list()
        # prepare product of all possible combinations of input arguments
        # for now we use only proc (cpu list) and files (file)
        for test_case in self.test_cases:
            tmp_result.append(list(itertools.product(
                ensure_iterable(test_case),
                ensure_iterable(test_case.proc),
                ensure_iterable(test_case.files),
            )))

        result = list()
        for lst in tmp_result:
            result.extend([prescription_class(*x) for x in lst])
        return result

    @classmethod
    def merge(cls, parent, children):
        """
        :type parent: dict
        """
        parent_copy = copy.deepcopy(parent)
        children = ensure_iterable(children)
        for child in children:
            parent_copy.update(child)
        return parent_copy

    def __iter__(self):
        self.iter_index = 0
        return self

    def next(self):
        """
        :rtype : scripts.config.yaml_config.YamlConfigCase
        """
        if self._iter_index >= len(self.test_cases):
            raise StopIteration
        else:
            self._iter_index += 1
            return self.test_cases[self._iter_index - 1]

    __next__ = next