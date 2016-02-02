# encoding: utf-8
# author:   Jan Hybs

import re, os, sys
from utils.dotdict import DotDict


class SimpleTableFormatter (object):
    """
    Class which takes json object from flow123d benchmark profiler report
     and returns simple table-like text string
    """

    def __init__ (self):
        self.json = None
        self.output = ""

        self.headerCols = []
        self.maxNameSize = 11

        self.bodyRows = []
        self.maxBodySize = None
        self.headerFields = ("tag", "call count", "max time", "max/min t", "avg time", "total", "source", "line")
        self.styles = {
            "linesep": os.linesep, "padding": 0,
            "min_width": 9, "colsep": '',
            "rowsep": '', "space_header": 3,
            "leading_char": " . ",
            "remove_prefix": "/src/"
        }

        self.totalDuration = None
        self.totalDurationMeasured = None

    def set_styles (self, styles):
        """Overrides default styles"""
        self.styles.update (styles)
        # make sure some values are actually ints
        self.styles["min_width"] = int (self.styles["min_width"])
        self.styles["padding"] = int (self.styles["padding"])
        self.styles["space_header"] = int (self.styles["space_header"])
        self.styles["colsep_start"] = self.styles["colsep"] + " "
        self.styles["colsep_end"] = " " + self.styles["colsep"]

    def convert_style (self):
        self.set_styles (self.styles)
        self.styles = DotDict (self.styles)

    def format (self, json):
        """"Formats given json object"""
        self.convert_style ()

        self.json = json
        self.process_header (json)
        self.process_body (json, 0)
        self.maxBodySize = [n + self.styles.padding for n in self.maxBodySize]

        self.maxNameSize = self.maxNameSize + self.styles.space_header

        lineDivider = (sum (self.maxBodySize) + 2 + len (self.maxBodySize) * 2) * self.styles.rowsep
        fmtHead = "{{:{self.maxNameSize}s}}{{}}{self.styles.linesep}".format (self=self)

        for pair in self.headerCols:
            self.output += fmtHead.format (*pair)

        self.output += lineDivider
        self.output += self.styles.linesep
        self.output += self.styles.colsep_start
        for i in range (len (self.headerFields)):
            fmt = "{{:^{maxBodySize}s}}{colsep}".format (maxBodySize=self.maxBodySize[i], colsep=self.styles.colsep_end)
            self.output += fmt.format (self.headerFields[i])
        self.output += self.styles.linesep
        self.output += lineDivider
        self.output += self.styles.linesep

        for tup in self.bodyRows:
            self.output += self.styles.colsep_start
            fields = []
            for i in range (len (self.maxBodySize)):
                fields.append (("{:" + tup[i][0] + "" + str (self.maxBodySize[i]) + "s}").format (tup[i][1]))
            self.output += self.styles.colsep_end.join (fields)
            self.output += self.styles.colsep_end + self.styles.linesep
            # self.output += fmtBody.format (*tup)

        self.output += lineDivider
        return self.output

    def append_to_header (self, name, value=None, linebreak=False):
        """Appends entry to header column list, if no value was given
        value from json object by given name will be taken
        """
        value = value if value is not None else self.json[name.lower ().replace (" ", "-")]
        self.headerCols.append ((name, str(value) + (linebreak if linebreak else "")))

        if self.maxNameSize < len (str (name)):
            self.maxNameSize = len (str (name))

    def append_to_body (self, values):
        """Appends entry to body row list.
        value is tupple of tupples, where inner tupper has two elements, first formatting character, and second value,
        formatting character is used in string format() method
        designating alignment
        < for left
        > for right
        ^ for center
        """
        self.bodyRows.append (values)

        # default empty array
        if self.maxBodySize is None:
            self.maxBodySize = [self.styles.min_width] * len (values)

        # update max length
        for i in range (len (self.maxBodySize)):
            self.maxBodySize[i] = max (self.maxBodySize[i], len (str (values[i][1])))

    def process_header (self, json):
        """Appends header information"""

        self.append_to_header ("Program name")
        self.append_to_header ("Program version")
        self.append_to_header ("Program branch")
        self.append_to_header ("Program revision")
        self.append_to_header ("Program build")

        if 'source-dir' in json:
            self.append_to_header ("Source dir")

        self.append_to_header ("Timer resolution", linebreak=self.styles.linesep)


        desc = re.sub ("\s+", " ", json["task-description"], re.M)
        self.append_to_header ("Task description", desc)
        self.append_to_header ("Task size", linebreak=self.styles.linesep)

        self.append_to_header ("Run process count")
        self.append_to_header ("Run started", json["run-started-at"])
        self.append_to_header ("Run ended", json["run-finished-at"])


    def process_body (self, json, level):
        """Recursive body processing"""

        # first occurrence of cumul-time-sum is whole-program's measured time

        if self.totalDurationMeasured is None and "cumul-time-sum" in json:
            self.totalDurationMeasured = json['cumul-time-sum']

        if level > 0:

            abs_prc = (json["cumul-time-sum"] / self.totalDurationMeasured) * 100
            rel_prc = json["percent"]
            path = json["file-path"]
            if str(json["file-path"]).startswith (self.styles.remove_prefix):
                path = path[len (self.styles.remove_prefix):]

            self.append_to_body ((
                ("<", "{abs_prc:6.2f} {leading} {rel_prc:5.2f} {tag}".format (
                    abs_prc=abs_prc, leading=self.styles.leading_char * (level - 1), rel_prc=rel_prc, tag=json["tag"])),
                ("^", "{:d}".format (json["call-count"])),
                ("^", "{:1.4f}".format (json["cumul-time-max"])),
                ("^", "{:1.4f}".format (json["cumul-time-max"] / json["cumul-time-min"])),
                ("^", "{:1.4f}".format (json["cumul-time-sum"] / json["call-count-sum"])),
                ("^", "{:1.4f}".format (json["cumul-time-sum"])),
                ("<", "{path:s}, {function:s}()".format (function=json["function"], path=path)),
                ("^", "{line:5d}".format (line=json["file-line"]))
            ))

        try:
            for child in json["children"]:
                self.process_body (child, level + 1)
        except Exception as e:
            pass

    def timedelta_milliseconds (self, td):
        return td.days * 86400000 + td.seconds * 1000 + td.microseconds / 1000

