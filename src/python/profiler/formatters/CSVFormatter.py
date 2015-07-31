# encoding: utf-8
# author:   Jan Hybs

from utils.strings import join_iterable
import re, os


class CSVFormatter (object):
    """
    Class which takes json object from flow123d benchmark profiler report
     and returns csv string
    """

    def __init__ (self):
        self.header = []
        self.body = []
        self.headerFields = (
            "percentage", "level", "tag", "call count",
            "max time", "max/min time", "avg time",
            "total", "function", "location"
        )
        self.styles = { "linesep": os.linesep, 'separator': ',', 'prefix': '"', 'suffix': '"' }
        self.linesep = os.linesep

    def format (self, json):
        """"format given json object"""
        self.json = json
        # self.process_header (json)
        self.process_body (json, 0)

        # tmpLst = []
        # tmpLst.append(self.headerFields)
        # tmpLst.extend(self.body)
        # self.maxWidth = self.fix_width (tmpLst)

        result = ""
        result += join_iterable (
            self.headerFields,
            separator=self.styles["separator"],
            prefix=self.styles["prefix"],
            suffix=self.styles["suffix"]) + self.linesep

        for row in self.body:
            result += join_iterable (
                row,
                separator=self.styles["separator"],
                prefix=self.styles["prefix"],
                suffix=self.styles["suffix"]) + self.linesep

        return result

    def set_styles (self, styles):
        """override default styles"""
        self.styles.update (styles)
        self.linesep = self.styles["linesep"]

    def append_to_body (self, values):
        """Appends entry to body row list.
        value is tupple of strings each string is represents cell value
        """
        self.body.append (values)
        pass

    def process_body (self, json, level):
        """Recursive body processing"""
        if level > 0:
            self.append_to_body ((
                "{:1.2f}".format (json["percent"]),
                "{:d}".format (level),
                "{:s}".format (json["tag"]),
                "{:d}".format (json["call-count"]),
                "{:1.4f}".format (json["cumul-time-max"]),
                "{:1.4f}".format (json["cumul-time-max"] / json["cumul-time-min"]),
                "{:1.4f}".format (json["cumul-time-sum"] / json["call-count-sum"]),
                "{:1.4f}".format (json["cumul-time-sum"]),
                "{:s}():{:d}".format (json["function"], json["file-line"]),
                "{:s}".format (json["file-path"])
            ))

        try:
            for child in json["children"]:
                self.process_body (child, level + 1)
        except:
            pass

    def fix_width (self, lst):
        """Gets maximum width for each element of iterable object property"""
        size = len (lst[0])
        maxWidth = [5] * size

        for values in lst:
            maxWidth = [max (maxWidth[i], len (str (values[i]))) for i in range (size)]
        return maxWidth