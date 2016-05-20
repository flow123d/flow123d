#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from datetime import datetime, timedelta


class Duration(object):
    """
    Simple class for converting duration value
    Supported format are either number value or HH:MM:SS format
    """
    @staticmethod
    def parse(value):
        value = str(value)
        try:
            return float(value)
        except ValueError:
            time = datetime.strptime(value,"%H:%M:%S")
            delta = timedelta(hours=time.hour, minutes=time.minute, seconds=time.second)
            return delta.total_seconds()
