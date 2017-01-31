#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
from scripts.comparisons.modules import ExecComparison
from scripts.core.base import Paths
# ----------------------------------------------


class Ndiff(ExecComparison):
    """
    Class Ndiff is exec class for comparing two files using external call
    ndiff
    """

    default_r_tol = 0.01
    default_a_tol = 0.0001

    def get_command(self, reference_filepath, other_filepath, **kwargs):
        """
        Method will call ndiff and pass location of two files to be compared
        other arguments are from config.yaml
           - r_tol
           - a_tol
        Arguments are optional so default value is set otherwise
        """
        return [
            Paths.ndiff(),
            '-r', str(kwargs.get('r_tol', self.default_r_tol)),
            '-a', str(kwargs.get('a_tol', self.default_a_tol)),
            Paths.abspath(reference_filepath),
            Paths.abspath(other_filepath)
        ]

# shortcut to Ndiff so Python naming convention is not violated
ndiff = Ndiff
