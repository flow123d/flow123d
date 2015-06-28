# encoding: utf-8
# author:   Jan Hybs

import re, os


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
        self.headerFields = ("tag", "call count", "max time", "max/min time", "avg time", "total", "source")
        self.styles = { "linesep": os.linesep, "padding": 0, "minwidth": 11 }
        self.linesep = os.linesep

    def set_styles (self, styles):
        """Overrides default styles"""
        self.styles.update (styles)
        self.linesep = self.styles["linesep"]
        self.styles["minwidth"] = int (self.styles["minwidth"])

    def format (self, json):
        """"Formats given json object"""
        self.json = json
        self.processHeader (json)
        self.processBody (json, 0)
        self.maxBodySize = [n + self.styles['padding'] for n in self.maxBodySize]

        lineDivider = (sum (self.maxBodySize) + 2 + len (self.maxBodySize) * 2) * "-"
        fmtHead = "{:" + str (self.maxNameSize + 2) + "s}{}" + self.linesep

        for pair in self.headerCols:
            self.output += fmtHead.format (*pair)

        self.output += lineDivider
        self.output += self.linesep
        self.output += "| "
        for i in range (len (self.headerFields)):
            self.output += ("{:^" + str (self.maxBodySize[i] + 1) + "s}|").format (self.headerFields[i])
        self.output += self.linesep
        self.output += lineDivider
        self.output += self.linesep

        for tup in self.bodyRows:
            self.output += "| "
            fields = []
            for i in range (len (self.maxBodySize)):
                fields.append (("{:" + tup[i][0] + "" + str (self.maxBodySize[i]) + "s}").format (tup[i][1]))
            self.output += " |".join (fields)
            self.output += " |" + self.linesep
            # self.output += fmtBody.format (*tup)

        self.output += lineDivider
        return self.output

    def appendToHeader (self, name, value=None):
        """Appends entry to header column list, if no value was given
        value from json object by given name will be taken
        """
        value = value if value is not None else self.json[name.lower ().replace (" ", "-")]
        self.headerCols.append ((name, value))

        if self.maxNameSize < len (str (name)):
            self.maxNameSize = len (str (name))

    def appendToBody (self, values):
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
            self.maxBodySize = [self.styles["minwidth"]] * len (values)

        # update max length
        for i in range (len (self.maxBodySize)):
            self.maxBodySize[i] = max (self.maxBodySize[i], len (str (values[i][1])))
        pass

    def processHeader (self, json):
        """Appends header information"""
        self.appendToHeader ("Program name")
        self.appendToHeader ("Program version")
        self.appendToHeader ("Program branch")
        self.appendToHeader ("Program revision")
        self.appendToHeader ("Program build")
        self.appendToHeader ("Timer resolution")

        if 'source-dir' in json:
            self.appendToHeader ("Source dir")

        desc = re.sub ("\s+", " ", json["task-description"], re.M)
        self.appendToHeader ("Task description", desc)
        self.appendToHeader ("Task size")

        self.appendToHeader ("Run process count")
        self.appendToHeader ("Run started", json["run-started-at"])
        self.appendToHeader ("Run ended", json["run-finished-at"])
        self.appendToHeader ("Run duration", json["run-finished-at"] - json["run-started-at"])

    def processBody (self, json, level):
        """Recursive body processing"""
        if level > 0:
            self.appendToBody ((
                ("<", "{:6.2f} {:s} {:s}".format (json["percent"], " >  " * (level - 1) * 1, json["tag"])),
                ("^", "{:d}".format (json["call-count"])),
                (">", "{:1.4f}".format (json["cumul-time-max"])),
                (">", "{:1.4f}".format (json["cumul-time-max"] / json["cumul-time-min"])),
                (">", "{:1.4f}".format (json["cumul-time-sum"] / json["call-count-sum"])),
                (">", "{:1.4f}".format (json["cumul-time-sum"])),
                ("<", "{:>24s} : {:<5d} {:s}".format (json["function"], json["file-line"], json["file-path"]))
            ))

        try:
            for child in json["children"]:
                self.processBody (child, level + 1)
        except:
            pass