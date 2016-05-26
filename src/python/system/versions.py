#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

"""
Simple module which provides method for testing running python version
"""

import sys


def require_version_2 ():
    """requires version higher than 3"""
    if sys.version_info > (3, 0):
        print ('Error: Python 2 is required')
        sys.exit (1)


def require_version_3 ():
    """requires version lower than 3"""
    if sys.version_info < (3, 0):
        print ('Error: Python 3 is required')
        sys.exit (1)
