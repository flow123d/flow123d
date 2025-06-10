#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import re
# ----------------------------------------------
from py123d.scripts.comparisons.modules import InPlaceComparison
from py123d.scripts.core.base import Paths, IO
# ----------------------------------------------

"""
To use this checker, create check_rules section in config.yaml file
Files are read and per-line check is initiated.

Specify either regex keyword where regular expression is set
and/or
substr keyword which performs fast and simple sub string search.


On first match reading is interrupted.

e.g.:

common_config:
  proc: [1, 2]
  memory_limit: 1000
  check_rules:
    - regex:
        files: ["*.*"]
        regex: "[Ee]rror"           # matches the word error or Error in line
        substr: "error"             # looks for the word error in a line
"""


class Regex(InPlaceComparison):
    """
    Class Regex is template for InPlaceComparison
    """

    def compare(self, reference_filepath, other_filepath, **kwargs):
        """
        Method can do anything as long as int value is returned
        :param reference_filepath:
        :param other_filepath:
        :param kwargs:
        :return: 1 - failed, 0 - success
        """

        regex = kwargs.get('regex', None)
        substr = kwargs.get('substr', None)

        if regex is None and substr is None:
            self.output.write("Invalid check rule! Specify either 'regex' or 'substr' keyword.")
            return 1

        # convert to str or compile to regex
        if substr is not None:
            substr = str(substr)
        if regex is not None:
            regex = re.compile(str(regex))

        # for regex comparison we ignore reference file
        # read job_output.log
        content = IO.read(
            Paths.abspath(other_filepath)
        )
        if content is None:
            return 1

        # substr find
        if substr:
            found = False
            for line in content.splitlines():
                if line.find(substr) != -1:
                    found = True
                    self.output.write('[OK] Found substr "{substr}" in line:'.format(**locals()))
                    self.output.write('    ' + line)
                    break
            if not found:
                self.output.write('[ERROR] Could not find substr "{substr}":'.format(**locals()))
                return 1

        # regex match
        if regex:
            found = False
            for line in content.splitlines():
                if regex.findall(line):
                    found = True
                    self.output.write('[OK] Found regex "{regex.pattern}" in line:'.format(**locals()))
                    self.output.write('    ' + line)
                    break
            if not found:
                self.output.write('[ERROR] Could not find regex {regex.pattern}:'.format(**locals()))
                return 1

        # self.output.write("In case of emergency,")
        # self.output.write("    you can provide details on what went wrong")
        # self.output.write("    using self.output.write method")
        # self.output.write("")
        # self.output.write("Error while comparing files \n{} \n{}"
        #                   .format(reference_filepath, other_filepath))

        # must return return-code!
        return 0

# shortcut to Regex so Python naming convention is not violated
regex = Regex
