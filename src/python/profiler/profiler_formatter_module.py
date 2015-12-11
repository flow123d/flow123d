# encoding: utf-8
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

import system.versions

system.versions.require_version_2 ()

import json
import datetime
import importlib


class ProfilerJSONDecoder (json.JSONDecoder):
    """Class overriding JSONDecoder which possess default python json decoding method.
    This class decodes given json string and converts values to proper type since everything is represented as string
    Basic conversion are:
        integer
        float
        datetime

    returned object has all values properly typed so
    formatters can make mathematical or other operation without worries
    """


    def decode (self, json_string):
        """Decodes json_string which is string that is given to json.loads method"""
        default_obj = super (ProfilerJSONDecoder, self).decode (json_string)

        self.intFields          = ["file-line", "call-count", "call-count-min", "call-count-max", "call-count-sum"]
        self.floatFields        = ["cumul-time", "cumul-time-min", "cumul-time-max", "cumul-time-sum", "percent", "run-duration"]
        self.intFieldsRoot      = ["task-size", "run-process-count"]
        self.floatFieldsRoot    = ["timer-resolution"]
        self.dateFields         = ["run-started-at", "run-finished-at"]


        self.convert_fields (default_obj, self.intFields, int)
        self.convert_fields (default_obj, self.floatFields, float)
        self.convert_fields (default_obj, self.intFieldsRoot, int, False)
        self.convert_fields (default_obj, self.floatFieldsRoot, float, False)
        self.convert_fields (default_obj, self.dateFields, self.parse_date, False)

        return default_obj


    def default_serializer (self, obj):
        """Default JSON serializer."""
        if isinstance (obj, datetime.datetime):
            return obj.strftime ("%m/%d/%y %H:%M:%S")
        return str (obj)

    def parse_date (self, str):
        """Default parsing method for date"""
        return datetime.datetime.strptime (str, "%m/%d/%y %H:%M:%S")

    def convert_fields (self, obj, fields, fun, rec=True):
        """Recursive value type conversion"""
        for field in fields:
            for prop in obj:
                if prop == field:
                    obj[prop] = fun (obj[prop])
        if rec:
            try:
                for child in obj["children"]:
                    self.convert_fields (child, fields, fun)
            except:
                pass


class ProfilerFormatter (object):
    """
    Class which dynamically loads formatter and perform conversion
    """


    @staticmethod
    def get_class_instance (cls):
        """Method returns class instance upon given name in profiler.formatters.* ns"""
        module = importlib.import_module ("profiler.formatters." + cls)
        class_ = getattr (module, cls)
        instance = class_ ()
        return instance


    @staticmethod
    def list_formatters ():
        """Method return list of all available formatters.
        Formatter is every class in profiler.formatters.* ns which possesses method format
        """
        result = []
        import pkgutil

        for module_loader, name, ispkg in pkgutil.iter_modules (['formatters']):
            try:
                module = importlib.import_module ("formatters." + name)
                class_ = getattr (module, name)
                if getattr (class_, 'format') is not None and callable (getattr (class_, 'format')):
                    result.append (name)
            except:
                pass
        return result

    def convert (self, json_location, output_file=None, formatter="SimpleTableFormatter", styles=[]):
        """Converts file @ json_location to output_file (if set) using given formatter name"""
        # read file to JSON
        #try:
        with open (json_location, 'r') as fp:
            jsonObj = json.load (fp, encoding="utf-8", cls=ProfilerJSONDecoder)
        #except Exception as exception:
        #    # return string with message on error
        #    return str (exception)

        #try:
        # split styles fields declaration
        styles = [value.replace ('\\n', '\n').replace ('\\t', '\t').replace ('\\r', '\r') for value in styles]
        styles = dict (item.split (":", 1) for item in styles)
        # grab instance and hand over styles
        instance = ProfilerFormatter.get_class_instance (formatter)
        instance.set_styles (styles)
        # format json object
        output = instance.format (jsonObj)
        #except Exception as exception:
        #    # return string with message on error
        #    return str (exception)

        #try:
        # if output file is specified write result there
        if output_file is not None:
            with open (output_file, "w") as fp:
                fp.write (output)
            print '{} file generated'.format (output_file)
        # otherwise just print result to stdout
        else:
            print output
        #except Exception as exception:
        #    # return string with message on error
        #    return str (exception)


        # return True on success
        #return True


def convert (json_location, output_file, formatter):
  """Simple iface method for c api"""
  fmt = ProfilerFormatter ()
  return fmt.convert (json_location, output_file, formatter)
