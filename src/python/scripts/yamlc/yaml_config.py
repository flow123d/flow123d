#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import itertools
import yaml
import re
# ----------------------------------------------
from copy import deepcopy
# ----------------------------------------------
from scripts import yamlc
from scripts.core.base import Paths
from scripts.core.base import PathFilters
# ----------------------------------------------
from utils.globals import ensure_iterable


class ConfigCase(object):
    """
    Class ConfigCase represents single ConfigCase
    :type config   : ConfigBase
    """

    sanitise_regex = re.compile(r'[^a-zA-Z0-9_-]+]')

    def __init__(self, o, config):
        o = ConfigBase.merge(yamlc.DEFAULTS, deepcopy(o))

        self.file = o.get(yamlc.TAG_FILES, None)
        self.proc = int(o.get(yamlc.TAG_PROC, None))
        self.time_limit = float(o.get(yamlc.TAG_TIME_LIMIT, None))
        self.memory_limit = float(o.get(yamlc.TAG_MEMORY_LIMIT, None))
        self.tags = set(o.get(yamlc.TAG_TAGS, None))
        self.args = list(o.get(yamlc.TAG_ARGS, []))
        self.check_rules = o.get(yamlc.TAG_CHECK_RULES, None)
        self.input = o.get(yamlc.TAG_INPUTS)
        self.death_test = yamlc.YamlDeathTest(o.get(yamlc.TAG_DEATH_TEST))
        self.config = config

        if self.config:
            self.file = Paths.join(self.config.root, Paths.basename(self.file))
            self.without_ext = Paths.basename(Paths.without_ext(self.file))

            # if input file is non standard we will have slightly different
            # output folder to ensure, we will not have all results in the same folder
            if self.input == 'input':
                self.shortname = '{name}.{proc}'.format(
                    name=self.without_ext,
                    proc=self.proc)
            else:
                self.shortname = '{name}.{proc}.{input}'.format(
                    name=self.without_ext,
                    proc=self.proc,
                    input=self.sanitise_regex.sub('', self.input))

            self.fs = yamlc.ConfigCaseFiles(
                root=self.config.root,
                ref_output=Paths.join(
                    self.config.root, yamlc.REF_OUTPUT_DIR, self.without_ext),
                output=Paths.join(
                    self.config.root,
                    yamlc.TEST_RESULTS,
                    self.shortname
                ),
                input=self.input
            )
        else:
            # create temp folder where files will be
            tmp_folder = Paths.temp_file(o.get('tmp') + '-{date}-{time}-{rnd}')
            Paths.ensure_path(tmp_folder, is_file=False)

            self.fs = yamlc.ConfigCaseFiles(
                root=tmp_folder,
                ref_output=tmp_folder,
                output=tmp_folder,
                input='input'
            )

    @property
    def as_string(self):
        if self.file:
            if self.input == 'input':
                return '{} x {}'.format(
                    self.proc,
                    Paths.path_end(Paths.without_ext(self.file), level=2)
                )
            else:
                return '{} x {} ({})'.format(
                    self.proc,
                    Paths.path_end(Paths.without_ext(self.file), level=2),
                    self.sanitise_regex.sub('_', self.input)
                )
        return 'process'

    @property
    def info(self):
        """
        Will try to gen test name and case name
        from yaml file
        :return: dict
        """
        if self.file:
            parts = str(self.file).split('/')
            return {
                'test-name': parts[-2],
                'case-name': parts[-1].split('.')[0],
            }
        return {}

    def to_json(self):
        return dict(
            cpu=self.proc,
            test_case=self.file,
        )

    def repr(self):
        return self.as_string

    __repr__ = repr


class ConfigBase(object):
    """
    Class ConfigBase represents single configuration yaml file
    """
    # dummy cases will be generated
    MISSING_POLICY_CREATE_DEFAULT   = 1
    # files not specified in config.yaml will not be executed
    MISSING_POLICY_IGNORE           = 2
    
    def __init__(self, yaml_config_file, missing_policy=MISSING_POLICY_CREATE_DEFAULT):
        self.yaml_config_file = yaml_config_file
        self.root = Paths.dirname(self.yaml_config_file)
        self.yamls = self._get_all_yamls()
        self.cases = list()
        self.common_config = None
        self.missing_policy = missing_policy

        # create dummy case for every yaml file in folder
        if not Paths.exists(self.yaml_config_file):
            self.common_config = deepcopy(yamlc.DEFAULTS)
            for y in self.yamls:
                dummy_case = deepcopy(yamlc.DEFAULTS)
                dummy_case['files'] = [y]
                self.cases.append(dummy_case)
        else:
            # setup common config values
            self.yaml_config = self._read_yaml()
            self.common_config = self.merge(
                yamlc.DEFAULTS, self.yaml_config.get('common_config', {}))

            # first process files which are specified in test_cases
            missing = [Paths.basename(y) for y in self.yamls]
            for case in self.yaml_config.get(yamlc.TAG_TEST_CASES, []):
                case_config = self.merge(self.common_config, case)

                # ensure that value are array
                case_config[yamlc.TAG_FILES] = ensure_iterable(
                    case_config.get(yamlc.TAG_FILES, []))
                case_config[yamlc.TAG_INPUTS] = ensure_iterable(
                    case_config.get(yamlc.TAG_INPUTS, []))
                case_config[yamlc.TAG_INPUTS] = [str(x) for x in case_config[yamlc.TAG_INPUTS]]

                # keep correct order
                self.cases.append(case_config)
                for f in case_config[yamlc.TAG_FILES]:
                    if f in missing:
                        missing.remove(f)

            # process rest (dummy case)
            if missing_policy == self.MISSING_POLICY_CREATE_DEFAULT:
                for y in missing:
                    dummy_case = deepcopy(self.common_config)
                    dummy_case[yamlc.TAG_FILES] = [y]
                    self.cases.append(dummy_case)

    def get_all(self):
        """
        :rtype: list[ConfigCase]
        """
        result = list()
        for case in self.cases:
            result.extend(self._get_all_for_case(case))
        return [ConfigCase(r, self) for r in result]

    def get_one(self, yaml_case_file):
        """
        :rtype: list[ConfigCase]
        """
        result = list()
        for case in self.cases:
            for f in case[yamlc.TAG_FILES]:
                if Paths.basename(f) == Paths.basename(yaml_case_file):
                    dummy_case = deepcopy(case)
                    dummy_case[yamlc.TAG_FILES] = [yaml_case_file]
                    result.extend(self._get_all_for_case(dummy_case))
        return [ConfigCase(r, self) for r in result]

    def _read_yaml(self):
        with open(self.yaml_config_file, 'r') as fp:
            result = yaml.load(fp)
        return result or dict()

    def _get_all_yamls(self):
        yamls = Paths.browse(
            self.root, (
                PathFilters.filter_endswith(yamlc.YAML),
                PathFilters.filter_not(
                    PathFilters.filter_endswith(yamlc.CONFIG_YAML))
            ))
        return yamls

    def update(self, proc, time_limit, memory_limit, input, **kwargs):
        for case in self.cases:
            if proc:
                case[yamlc.TAG_PROC] = set(proc)
            if time_limit:
                case[yamlc.TAG_TIME_LIMIT] = time_limit
            if memory_limit:
                case[yamlc.TAG_MEMORY_LIMIT] = memory_limit
            if input:
                case[yamlc.TAG_INPUTS] = ensure_iterable(input)

    def filter_tags(self, include, exclude):
        """
        :type exclude: set
        :type include: set
        """
        # skip tests if no --include and --exclude are set
        if not include and not exclude:
            return

        result = []
        for case in self.cases:
            tags = set(case['tags'])
            match = True

            # every tag in include must be present in tags
            match &= not include or include.issubset(tags)
            # every tag in exclude must not br present in tags
            match &= not exclude or not exclude.intersection(tags)

            if match:
                result.append(case)

        self.cases = result

    @classmethod
    def _get_all_for_case(cls, case):
        result = list()
        # product of all files, cpus and inputs
        changes = list(itertools.product(
            case[yamlc.TAG_FILES],
            case[yamlc.TAG_PROC],
            case[yamlc.TAG_INPUTS],
        ))

        # copy original and edit changes
        for f, p, i in changes:
            dummy_case = deepcopy(case)
            dummy_case[yamlc.TAG_FILES] = f
            dummy_case[yamlc.TAG_PROC] = p
            dummy_case[yamlc.TAG_INPUTS] = i
            result.append(dummy_case)
        return result

    @classmethod
    def merge(cls, parent, *children):
        """
        :type parent: dict
        """
        parent_copy = deepcopy(parent)
        for child in children:
            parent_copy.update(child)
        return parent_copy


class ConfigPool(object):
    """
    Class ConfigPool collects configs from multiple yaml files
    :type configs : dict[str, ConfigBase]
    :type files : dict[str, ConfigBase]
    """

    # match only yaml files not in ref_out or test_results having config.yaml
    #  in the same directory (excluding config.yaml itself)
    yaml_filters= [
        PathFilters.filter_type_is_file(),
        PathFilters.filter_ignore_dirs([yamlc.TEST_RESULTS, yamlc.TEST_RESULTS]),
        PathFilters.filter_not(PathFilters.filter_name('config.yaml')),
        PathFilters.filter_endswith('.yaml'),
        PathFilters.filter_dir_contains_file('config.yaml'),
    ]

    def __init__(self):
        self.configs = dict()
        self.files = dict()

    def add_config(self, yaml_config_file):
        self.configs[yaml_config_file] = None
        return self

    def add_case(self, yaml_case_file):
        config = Paths.join(Paths.dirname(yaml_case_file), yamlc.CONFIG_YAML)
        self.configs[config] = None
        self.files[yaml_case_file] = None
        return self

    def add(self, yaml_file):
        """
        :type yaml_file: str
        """
        if yaml_file.endswith(yamlc.CONFIG_YAML):
            return self.add_config(yaml_file)
        return self.add_case(yaml_file)

    def parse(self, missing_policy=ConfigBase.MISSING_POLICY_CREATE_DEFAULT):
        for k, v in list(self.configs.items()):
            self.configs[k] = ConfigBase(k, missing_policy)

        for k, v in list(self.files.items()):
            config = Paths.join(Paths.dirname(k), yamlc.CONFIG_YAML)
            self.files[k] = self.configs[config]

    __iadd__ = add

    def update(self, proc, time_limit, memory_limit, input, **kwargs):
        for config in list(self.configs.values()):
            config.update(proc, time_limit, memory_limit, input, **kwargs)

    def filter_tags(self, include, exclude):
        include = set(include) if include else None
        exclude = set(exclude) if exclude else None

        for config in list(self.configs.values()):
            config.filter_tags(include, exclude)
