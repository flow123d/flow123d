#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import os
import shutil

import formic

from scripts.artifacts.artifacts import ArtifactStep
from scripts.core.base import Paths, Printer


class Collector(ArtifactStep):
    """
    Class Collector is simple class which collects files from certain
    directory and copies them to specific output directory
    """

    yaml_tag = u'!Collector'

    def __init__(self, source=None, target=None, includes='*', excludes=None, flat=False, name=None, remove_original=False, **kwargs):
        self.source = source
        self.target = target
        self.includes = includes
        self.excludes = excludes
        self.flat = flat
        self.name = name
        self.remove_original = remove_original
        self.__dict__.update(kwargs)

    @staticmethod
    def create_path_dict(filename):
        path = Paths.split(filename)[1:-1]
        result = dict()
        for i in range(len(path)):
            result['+' + str(i)] = path[i]
            result['-' + str(i)] = path[-i]
        return result

    def __iter__(self):
        fileset = formic.FileSet(self.includes, directory=self.source)
        for filename in fileset:

            name = Paths.basename(filename) if not self.name else self.name.format(
                path=self.create_path_dict(filename),
                name=Paths.basename(filename),
            )
            if self.flat:
                root = self.target
            else:
                rel_path = Paths.relpath(Paths.dirname(filename), Paths.abspath(self.source))
                root = Paths.abspath(Paths.join(self.target, rel_path))

            yield CopyRule(filename, Paths.join(root, name), self.remove_original)

    def run(self):
        total = 0
        for p in self:
            total += 1
            p.copy()
        Printer.all.out("Copied out {} files", total)


class CopyRule(object):
    """
    Class CopyRule represents action to copy single file to certain destination
    """

    def __init__(self, source, target, remove_original):
        self.source = source
        self.target = target
        self.remove_original = remove_original

    def __repr__(self):
        return 'CopyRule({self.source} => {self.target})'.format(self=self)

    def copy(self):
        # create dirs for target file
        Paths.ensure_path(self.target)

        # copy file
        shutil.copy(self.source, self.target)

        # remove file if set
        if self.remove_original:
            os.unlink(self.source)