#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import os
import time
import yaml

from scripts.core.base import Paths, Printer
from utils import strings
from utils.counter import ProgressCounter


class YAMLTag(yaml.YAMLObject):
    """
    Class YAMLTag is simple class which calls constructor
    after parsing
    """

    def __setstate__(self, state):
        """ :type state: dict """
        self.__init__(**state)


class ArtifactStep(YAMLTag):
    """
    Class ArtifactStep represents single process while
    processing artifacts.yml file (can be copy/move files or executing external process)
    """

    def run(self):
        raise Exception('Method not implemented')

    def __repr__(self):
        return self.__class__.__name__


class ArtifactProcessor(object):
    """
    Class ArtifactProcessor will parse yaml artifact file and
    run specific action such as copy files, or publish to database

    :type configuration : dict
    """

    def __init__(self, yaml_file):
        self.yaml_file = Paths.abspath(yaml_file)
        self.configuration = None

        # set default language as english
        os.environ['LC_ALL'] = 'C'

    @property
    def steps(self):
        """
        :rtype: list[ArtifactStep]
        """
        if not self.configuration:
            self.parse_yaml()
        return self.configuration.get('collectors', []) or []

    def run(self):
        counter = ProgressCounter('Artifact step {:02d} / {total:02d}: {step}')
        total = len(self.steps)

        for step in self.steps:
            counter.next(locals())
            with Printer.all.with_level(1):
                step.run()

    def parse_yaml(self):
        # register yaml parser tags
        from scripts.artifacts.collector import Collector
        from scripts.artifacts.command import Command
        from scripts.artifacts.modules.mongodb import DatabaseMongo
        from scripts.artifacts.modules.lscpu import CommandLSCPU

        with open(self.yaml_file, 'r') as fp:
            yaml_data = fp.read()

        yaml_data = strings.replace_placeholders(
            yaml_data,
            _format_ = '<{}>',

            root=Paths.flow123d_root(),
            time=time.strftime("%H-%M-%S"),
            date=time.strftime("%y.%m.%d"),
            datetime=time.strftime("%y.%m.%d_%H-%M-%S"),
        )
        self.configuration = yaml.load(yaml_data) or {}
