#!/usr/bin/python
# -*- coding: utf-8 -*-
# author: Jan Hybs

import os
import yaml


class Config(object):
    _instance = None
    _cfg = dict()

    @classmethod
    def _load_config(cls):
        if cls._instance:
            return cls._instance

        current = 0
        max_limit = 10
        root = os.path.dirname(__file__)

        config_file = os.path.join(root, '.config.yaml')
        while current < max_limit:
            config_file = os.path.join(root, '.config.yaml')
            current += 1
            if os.path.exists(config_file):
                break
            root = os.path.dirname(root)
        # load config and extract username and password

        cls._instance = Config()
        cls._cfg = yaml.load(open(config_file, 'r').read()) or dict()
        return cls._instance

    def get(self, key=None, default=None):
        """
        :rtype: dict or list or string or int or bool
        """
        self.__class__._load_config()
        if key is None:
            return self.__class__._cfg
        return self.__class__._cfg.get(key, default)

    def set(self, key, value):
        """
        :rtype: dict or list or string or int or bool
        """
        self.__class__._load_config()
        self.__class__._cfg[key] = value
        return value

    def __getitem__(self, item):
        return self.get(item)

cfg = Config()
