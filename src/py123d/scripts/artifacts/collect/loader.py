#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import os
import fnmatch

from py123d.scripts.artifacts.collect.modules import CollectResult
from py123d.scripts.artifacts.db.mongo import Mongo


def load_data(location, module, rules=None):
    """
    Method will collect and parse results from given location using specific module
    :rtype: list[CollectResult]
    :type module: scripts.artifacts.collect.modules.flow123d_profiler.Flow123dProfiler
    :param location: path to folder
    :param module: module to use (such as flow123d_profiler)
    :param rules: file filters
    :return:
    """
    if not os.path.exists(location):
        print('[ERROR] Invalid location ', location)
        return list()

    if not os.path.isdir(location):
        location = os.path.dirname(location)

    def recursive_glob(root, pattern):
        """
        Recursive way to go through file structure matching pattern rules
        :param root:
        :param pattern:
        :return:
        """
        results = []
        for base, dirs, files in os.walk(root):
            match = fnmatch.filter(files, pattern)
            results.extend([os.path.join(base, f) for f in match])
        return results

    files = recursive_glob(location, module.include)
    if rules:
        for rule in rules:
            files = [f for f in files if rule(f)]

    result = []
    for file in files:
        result.append(module.process_file(file))

    return result


def save_to_database(data):
    """
    :type data: list[CollectResult]
    """
    mongo = Mongo()
    results = list()
    for item in data:

        # save logs first
        if item.logs and item.items:
            log_ids = mongo.fs.insert_many(item.logs).inserted_ids
            item.update(log_ids)

        # insert rest to db
        if item.items:
            results.append(mongo.flow_long.insert_many(item.items))
    return results
