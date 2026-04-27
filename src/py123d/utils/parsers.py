#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import re
from datetime import datetime, timedelta


def parse_float(value):
    """
    Simple method for converting duration value
    Supported format are either number value or HH:MM:SS format

    :type value: str
    :rtype: float
    """
    value = str(value)
    try:
        return float(value)
    except ValueError:
        time = datetime.strptime(value, "%H:%M:%S")
        delta = timedelta(hours=time.hour, minutes=time.minute, seconds=time.second)
        return delta.total_seconds()


def parse_int_list(value):
    """
    Simple method for converting integers in multiple formats
      The value can be:
         - single number
         - set             "[1,3,4]" or "[1 2 4]"
         - range           "1:4"   = "[1,2,3,4]"
         - range with step "1:7:2" = "[1,3,5,7]"
      Method will always return list of integers even if single
      number is entered

    :rtype: list[int]
    :type value: str
    """
    if value.startswith('[') and value.endswith(']'):
        return eval(re.sub(r' +', ',', value))

    if value.find(':') != -1:
        parts = [int(x.strip()) for x in value.split(':')]
        args = [1, 1, 1]
        args[0:len(parts)] = parts
        args[1] += 1
        return list(range(*args))

    return [int(value.strip())]