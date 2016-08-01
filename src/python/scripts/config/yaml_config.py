#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import yaml
import itertools
# ----------------------------------------------
from copy import deepcopy
# ----------------------------------------------
from scripts.core.base import Paths
from scripts.core.base import PathFilters
# ----------------------------------------------
from utils.globals import ensure_iterable

YAML = '.yaml'
CONFIG_YAML = 'config.yaml'

DEFAULTS = dict(
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

TAG_FILES = 'files'
TAG_PROC = 'proc'
TAG_TIME_LIMIT = 'time_limit'
TAG_MEMORY_LIMIT = 'memory_limit'
TAG_TEST_CASES = 'test_cases'
TAG_CHECK_RULES = 'check_rules'
TAG_TAGS = 'tags'

REF_OUTPUT_DIR = 'ref_out'

class ConfigPool(object):
    """
    :type configs : dict[str, ConfigBase]
    :type files : dict[str, ConfigBase]
    """
    def __init__(self):
        self.configs = dict()
        self.files = dict()

    def add_config(self, yaml_config_file):
        self.configs[yaml_config_file] = None
        return self

    def add_case(self, yaml_case_file):
        config = Paths.join(Paths.dirname(yaml_case_file), CONFIG_YAML)
        self.configs[config] = None
        self.files[yaml_case_file] = None
        return self

    def add(self, yaml_file):
        """
        :type yaml_file: str
        """
        if yaml_file.endswith(CONFIG_YAML):
            return self.add_config(yaml_file)
        return self.add_case(yaml_file)

    def parse(self):
        for k, v in self.configs.items():
            self.configs[k] = ConfigBase(k)

        for k, v in self.files.items():
            config = Paths.join(Paths.dirname(k), CONFIG_YAML)
            self.files[k] = self.configs[config]

    __iadd__ = add

    def update(self, proc, time_limit, memory_limit, **kwargs):
        for config in self.configs.values():
            config.update(proc, time_limit, memory_limit, **kwargs)


class ConfigCaseFiles(object):
    def __init__(self, root, ref_output, output):
        """
        :type ref_output: str
        :type output: str
        :type root: str
        """
        self.root = root
        self.output = output
        self.ndiff_log = self.in_output('ndiff.log')

        self.pbs_script = self.in_output('pbs_script.qsub')
        self.pbs_output = self.in_output('pbs_output.log')

        self.job_output = self.in_output('job_output.log')
        self.json_output = self.in_output('result.json')

        self.input = self.in_root('input')
        self.ref_output = ref_output

    def in_root(self, *names):
        """
        :rtype: str
        """
        return Paths.join(self.root, *names)

    def in_output(self, *names):
        """
        :rtype: str
        """
        return Paths.join(self.output, *names)


class ConfigCase(object):
    """
    :type config   : scripts.config.base.ConfigBase
    """
    def __init__(self, o, config):
        o = ConfigBase.merge(DEFAULTS, deepcopy(o))

        self.file = o.get(TAG_FILES, None)
        self.proc = int(o.get(TAG_PROC, None))
        self.time_limit = float(o.get(TAG_TIME_LIMIT, None))
        self.memory_limit = float(o.get(TAG_MEMORY_LIMIT, None))
        self.tags = set(o.get(TAG_TAGS, None))
        self.check_rules = o.get(TAG_CHECK_RULES, None)
        self.config = config

        if self.config:
            self.file = Paths.join(self.config.root, self.file)
            self.without_ext = Paths.basename(Paths.without_ext(self.file))
            self.shortname = '{name}.{proc}'.format(name=self.without_ext, proc=self.proc)

            self.fs = ConfigCaseFiles(
                root=self.config.root,
                ref_output=Paths.join(self.config.root, REF_OUTPUT_DIR, self.without_ext),
                output=Paths.join(
                    self.config.root,
                    'test_results',
                    self.shortname
                ))
        else:
            # create temp folder where files will be
            tmp_folder = Paths.temp_file(o.get('tmp') + '-{date}-{time}-{rnd}')
            Paths.ensure_path(tmp_folder, is_file=False)

            self.fs = ConfigCaseFiles(
                root=tmp_folder,
                ref_output=tmp_folder,
                output=tmp_folder
            )

    def to_string(self):
        if self.file:
            return '{} x {}'.format(
                self.proc,
                Paths.path_end(Paths.without_ext(self.file))
            )
        return 'process'

    def to_json(self):
        return dict(
            cpu=self.proc,
            test_case=self.file,
        )
    __repr__ = to_string


class ConfigBase(object):
    def __init__(self, yaml_config_file):
        self.yaml_config_file = yaml_config_file
        self.root = Paths.dirname(self.yaml_config_file)
        self.yamls = self._get_all_yamls()
        self.cases = list()
        self.common_config = None

        # create dummy case for every yaml file in folder
        if not Paths.exists(self.yaml_config_file):
            self.common_config = deepcopy(DEFAULTS)
            for y in self.yamls:
                dummy_case = deepcopy(DEFAULTS)
                dummy_case['file'] = [y]
                self.cases.append(dummy_case)
        else:
            # setup common config values
            self.yaml_config = self._read_yaml()
            self.common_config = self.merge(DEFAULTS, self.yaml_config.get('common_config', {}))

            # first process files which are specified in test_cases
            missing = [Paths.basename(y) for y in self.yamls]
            for case in self.yaml_config.get(TAG_TEST_CASES, []):
                case_config = self.merge(self.common_config, case)

                # ensure that value is array
                case_config[TAG_FILES] = ensure_iterable(case_config.get(TAG_FILES, []))
                self.cases.append(case_config)
                for f in case_config[TAG_FILES]:
                    if f in missing:
                        missing.remove(f)

            # process rest (dummy case)
            for y in missing:
                dummy_case = deepcopy(self.common_config)
                dummy_case[TAG_FILES] = [y]
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
            for f in case[TAG_FILES]:
                if Paths.basename(f) == Paths.basename(yaml_case_file):
                    dummy_case = deepcopy(case)
                    dummy_case[TAG_FILES] = [yaml_case_file]
                    result.extend(self._get_all_for_case(dummy_case))
        return [ConfigCase(r, self) for r in result]

    def _read_yaml(self):
        with open(self.yaml_config_file, 'r') as fp:
            return yaml.load(fp)

    def _get_all_yamls(self):
        yamls = Paths.browse(
            self.root,(
                PathFilters.filter_endswith(YAML),
                PathFilters.filter_not(PathFilters.filter_endswith(CONFIG_YAML))
            ))
        return yamls

    def update(self, proc, time_limit, memory_limit, **kwargs):
        for case in self.cases:
            if proc:
                case[TAG_PROC] = set(proc)
            if time_limit:
                case[TAG_TIME_LIMIT] = time_limit
            if memory_limit:
                case[TAG_MEMORY_LIMIT] = memory_limit

    @classmethod
    def _get_all_for_case(cls, case):
        result = list()
        changes = list(itertools.product(
            case[TAG_FILES],
            case[TAG_PROC]))

        for f, p in changes:
            dummy_case = deepcopy(case)
            dummy_case[TAG_FILES] = f
            dummy_case[TAG_PROC] = p
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
