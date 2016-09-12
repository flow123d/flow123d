#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import pickle


class ConfigCaseResult(object):
    def __init__(self, case):
        """
        :type case: scripts.config.yaml_config.ConfigCase
        """
        self.parent = case.config.yaml_config_file
        self.file = case.file
        self.proc = case.proc
        self.time_limit = case.time_limit
        self.memory_limit = case.memory_limit
        self.check_rules = case.check_rules
        self.repr = case.to_string()


class PyPyResult(object):
    def __init__(self, pypy):
        """
        :type pypy: scripts.core.threads.PyPy
        """
        self.returncode = pypy.returncode
        self.command = pypy.executor.command
        self.name = pypy.name
        self.output = pypy.full_output

        if pypy.case:
            self.case = ConfigCaseResult(pypy.case)


class CleanResult(object):
    def __init__(self, clean):
        """
        :type clean: scripts.prescriptions.local_run.CleanThread
        """
        self.returncode = clean.returncode
        self.dir = clean.dir
        self.name = clean.name
        self.error = clean.error


class ComparisonResult(object):
    def __init__(self, comp):
        """
        :type comp: scripts.core.threads.ComparisonMultiThread
        """
        self.returncode = comp.returncode
        self.name = comp.name
        self.output = comp.output
        self.items = list()
        for thread in comp.threads:
            self.items.append(dict(
                name=thread.name,
                returncode=thread.returncode,
                log=self.output,
            ))


class RuntestTripletResult(object):
    """
    Pickle container for RuntestMultiThread
    """
    def __init__(self, pypy, clean, comp):
        self.pypy = PyPyResult(pypy)
        self.clean = CleanResult(clean)
        self.comp = ComparisonResult(comp)


class ResultHolderResult(object):
    """
    Class ResultHolderResult is pickle container for ResultHolder
    """
    def __init__(self, result):
        """
        :type result: scripts.core.threads.ResultHolder
        """
        self.returncode = result.returncode
        self.items = list()
        for item in result.items:
            self.items.append(item.dump() if hasattr(item, 'dump') else item)


class ResultParallelThreads(object):
    """
    Class ResultParallelThreads is pickle container for ParallelThreads
    :type items : list[RuntestTripletResult]
    """
    def __init__(self, runner):
        """
        :type runner: scripts.core.threads.ParallelThreads
        """
        self.returncode = runner.returncode
        self.items = list()
        for thread in runner.threads:
            self.items.append(thread.dump())



def load_pypy(file):
    """
    :rtype: PyPyResult
    """
    return pickle.load(open(file, 'rb'))


def load_triplet(file):
    """
    :rtype: RuntestTripletResult
    """
    return pickle.load(open(file, 'rb'))


def load_runtest(file):
    """
    :rtype: ResultParallelThreads
    """
    return pickle.load(open(file, 'rb'))