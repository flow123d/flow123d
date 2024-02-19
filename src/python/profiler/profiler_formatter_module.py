#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

"""
Module for converting benchmark json reports from project Flow123d into other formats
Package contains:
    class ProfilerJSONDecoder for json parsing
    class ProfilerFormatter setting and conversion
    classes CSVFormatter and SimpleTableFormatter which upon receiving json object returns formatted string
    iface method convert

@url https://github.com/flow123d/flow123d
"""

import os
import io
import json
import datetime
import time
import importlib
from utils.logger import Logger
import collections


class ProfilerJSONDecoder(json.JSONDecoder):
    """Class overriding JSONDecoder which possess default python json decoding method.
    This class decodes given json string and converts values to proper type since everything is represented as string
    Basic conversion are:
        integer
        float
        datetime

    returned object has all values properly typed so
    formatters can make mathematical or other operation without worries
    """
    pass

    def decode(self, json_string):
        """Decodes json_string which is string that is given to json.loads method"""
        default_obj = super(ProfilerJSONDecoder, self).decode(json_string)

        self.intFields = ["file-line", "call-count", "call-count-min", "call-count-max", "call-count-sum"]
        self.floatFields = ["cumul-time", "cumul-time-min", "cumul-time-max", "cumul-time-sum", "percent",
                            "run-duration"]
        self.intFieldsRoot = ["task-size", "run-process-count"]
        self.floatFieldsRoot = ["timer-resolution"]
        self.dateFields = ["run-started-at", "run-finished-at"]

        self.convert_fields(default_obj, self.intFields, int)
        self.convert_fields(default_obj, self.floatFields, float)
        self.convert_fields(default_obj, self.intFieldsRoot, int, False)
        self.convert_fields(default_obj, self.floatFieldsRoot, float, False)
        self.convert_fields(default_obj, self.dateFields, self.parse_date, False)

        return default_obj

    def default_serializer(self, obj):
        """Default JSON serializer."""
        if isinstance(obj, datetime.datetime):
            return obj.strftime("%m/%d/%y %H:%M:%S")
        return str(obj)

    def parse_date(self, str):
        """Default parsing method for date"""
        return datetime.datetime.strptime(str, "%m/%d/%y %H:%M:%S")

    def convert_fields(self, obj, fields, fun, rec=True):
        """Recursive value type conversion"""
        for field in fields:
            for prop in obj:
                if prop == field:
                    obj[prop] = fun(obj[prop])
        if rec:
            try:
                for child in obj["children"]:
                    self.convert_fields(child, fields, fun)
            except:
                pass
            

def get_formater_instance(cls, styles=[]):
    """Method returns class instance upon given name in profiler.formatters.* ns"""
    module = importlib.import_module("profiler.formatters." + cls)
    class_ = getattr(module, cls)
    instance = class_()

    styles = [value.replace('\\n', '\n').replace('\\t', '\t').replace('\\r', '\r') for value in styles]
    styles = dict(item.split(":", 1) for item in styles)
    instance.set_styles(styles)
    return instance


def list_formatters():
    """Method return list of all available formatters.
    Formatter is every class in profiler.formatters.* ns which possesses method format
    """
    result = []
    import pkgutil

    for module_loader, name, ispkg in pkgutil.iter_modules(['formatters']):
        try:
            module = importlib.import_module("formatters." + name)
            class_ = getattr(module, name)
            if getattr(class_, 'format') is not None and isinstance(getattr(class_, 'format'), collections.Callable):
                result.append(name)
        except:
            pass
    return result



def convert(json_location, output_file, formatter, styles=[]):
    """
    Simple iface method for c api
    - it never fails, just tries to log operations and exceptions
    - retries in order to deal with some communication delays
    """
    fmt = get_formater_instance(formatter, styles)
    
    # 
    timeout = 2
    end_time = time.time() + timeout
    n_tries = 0
    while time.time() < end_time and n_tries < 2:
        try:
            with open(json_location, 'r') as f_in:
                json_data = json.load(f_in, encoding="utf-8", cls=ProfilerJSONDecoder)
            Logger.instance().info('File "%s" read', json_location)
            output = fmt.format(json_data)
            with open(output_file, "w") as fp:
                fp.write(output)
            Logger.instance().info('File "%s" generated', output_file)
            return True
        except FileExistsError as e:
            Logger.instance().exception(r"Can not open the file: {json_location}. Waiting.")
            time.sleep(0.2)
            continue
        except Exception as e:
            Logger.instance().exception("Error in profiler JSON. Retrying.", e)
            n_tries += 1
            continue
        break
    Logger.instance().error("Unable to format the profiler output.")
    print(Logger.instance().get_memory_stream())
    return False
