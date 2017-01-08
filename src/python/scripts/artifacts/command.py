#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from scripts.artifacts.artifacts import ArtifactStep
from scripts.core.base import Printer
from scripts.core.execution import OutputMode


class Command(ArtifactStep):
    """
    Class Command is abstract class for artifact step
    Each step does certain action (such as copy files)
    """
    yaml_tag = u'!Command'

    def __init__(self, command, output=None, **kwargs):
        self.command = [str(c) for c in command] if command else None
        self.output = OutputMode.file_write(output) if output else OutputMode.variable_output()
        self.__dict__.update(kwargs)

    def run(self):
        import subprocess
        process = subprocess.Popen(self.command, stdout=self.output.open(), stderr=subprocess.STDOUT)
        returncode = process.wait()
        self.output.close()

        Printer.all.out('Command ended with {}', returncode)
