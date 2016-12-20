#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from scripts.comparisons.modules import InPlaceComparison
from scripts.core.base import Paths
from scripts.core.base import IO
# ----------------------------------------------

"""
To use this checker, create check_rules section in config.yaml file

e.g.:

common_config:
  proc: [1, 2]
  memory_limit: 1000
  check_rules:
    - regex:
        files:
          - "*.*"
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
        :return:
        """

        reference_content = IO.read(
            Paths.abspath(reference_filepath)
        )
        other_content = IO.read(
            Paths.abspath(other_filepath)
        )

        self.output.write("In case of emergency,")
        self.output.write("    you can provide details on what went wrong")
        self.output.write("    using self.output.write method")
        self.output.write("")
        self.output.write("Error while comparing files \n{} \n{}"
                          .format(reference_filepath, other_filepath))

        # must return return-code!
        return 1

# shortcut to Regex so Python naming convention is not violated
regex = Regex
