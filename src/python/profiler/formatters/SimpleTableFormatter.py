__author__ = 'jan-hybs'
import re
import os

class SimpleTableFormatter (object) :


    def __init__(self):
        self.json = None
        self.output = ""

        self.headerCols = []
        self.maxNameSize = 12

        self.bodyRows = []
        self.maxBodySize = None
        self.headerFields = ("tag", "call count", "max time", "max/min time", "avg time", "total", "source")
        self.styles = {"linesep": os.linesep, "padding": 0}
        self.linesep = os.linesep

    def set_styles (self, styles):
        self.styles.update (styles)
        self.linesep = self.styles["linesep"]

    def format (self, json) :
        self.json = json
        self.processHeader (json)
        self.processBody (json, 0)
        self.maxBodySize = [n + self.styles['padding'] for n in self.maxBodySize]

        lineDivider = (sum (self.maxBodySize) + 2 + len (self.maxBodySize) * 2) * "-"
        fmtHead = "{:" + str (self.maxNameSize + 2) + "s}{}" + self.linesep

        for pair in self.headerCols : self.output += fmtHead.format (*pair)

        self.output += lineDivider
        self.output += self.linesep
        self.output += "| "
        for i in range (len (self.headerFields)) :
            self.output += ("{:^" + str (self.maxBodySize[i]+1) + "s}|").format (self.headerFields[i])
        self.output += self.linesep
        self.output += lineDivider
        self.output += self.linesep

        for tup in self.bodyRows :
            self.output += "| "
            fields = []
            for i in range (len (self.maxBodySize)) :
                fields.append (("{:" + tup[i][0] + "" + str (self.maxBodySize[i]) + "s}").format (tup[i][1]))
            self.output += " |".join(fields)
            self.output += " |" + self.linesep
            # self.output += fmtBody.format (*tup)

        self.output += lineDivider
        return self.output


    def appendToHeader (self, name, value=None) :
        value = value if value is not None else self.json[name.lower ().replace (" ", "-")]
        self.headerCols.append ((name, value))

        if self.maxNameSize < len (str (name)) :
            self.maxNameSize = len (str (name))

    def appendToBody (self, values) :
        self.bodyRows.append (values)

        if self.maxBodySize == None :
            self.maxBodySize = [11] * len (values)

        # update max length
        self.maxBodySize = [max (self.maxBodySize[i], len (str (values[i][1]))) for i in range (len (self.maxBodySize))]
        pass

    def processHeader (self, json) :
        self.appendToHeader ("Program name")
        self.appendToHeader ("Program version")
        self.appendToHeader ("Program branch")
        self.appendToHeader ("Program revision")
        self.appendToHeader ("Program build")

        self.appendToHeader ("Timer resolution")

        desc = re.sub ("\s+", " ", json["task-description"], re.M)
        self.appendToHeader ("Task description", desc)
        self.appendToHeader ("Task size")

        self.appendToHeader ("Run process count")
        self.appendToHeader ("Run started", json["run-started-at"])
        self.appendToHeader ("Run ended", json["run-finished-at"])
        self.appendToHeader ("Run duration", json["run-finished-at"] - json["run-started-at"])


    def processBody (self, json, level) :

        if level > 0 :
            self.appendToBody ((
                ("<", "{:6.2f} {:s} {:s}".format (json["percent"], " >  " * (level - 1) * 1, json["tag"])),
                ("^", "{:d}".format (json["call-count"])),
                (">", "{:1.4f}".format (json["cumul-time-max"])),
                (">", "{:1.4f}".format (json["cumul-time-max"] / json["cumul-time-min"])),
                (">", "{:1.4f}".format (json["cumul-time-sum"] / json["call-count-sum"])),
                (">", "{:1.4f}".format (json["cumul-time-sum"])),
                ("<", "{:>24s} : {:<5d} {:s}".format (json["function"], json["file-line"], json["file-path"]))
            ))

        try :
            for child in json["children"] :
                self.processBody (child, level + 1)
        except :
            pass