#!/usr/bin/env python

import os
import sys
import re

class TagInfo:
    def __init__(self, tag, time, taskSize):
        self.taskSize = taskSize
        self.tag = tag
        self.time = time



#process the output of timers with following tags
tagsToProcess = ["WHOLE PROGRAM","SOLVING MH SYSTEM"]

def main():
    #check whether the user entered an existing directory
    try:
    	dir = sys.argv[1]
    except IndexError:
    	sys.exit("No directory specified")

    dir = os.path.abspath(dir)
    if os.path.exists(dir) == 0 or os.path.isdir(dir) == 0:
    	sys.exit("Specified directory doesn't exist")


    includeSubdirs = 1

    #dictionary containing the information gained from the profiler output files,
    #the key of the dictionary is the number of processors, the value is another dictionary
    #containing timer tags as keys and TagInfo objects as values
    fileInfo = {}

    processDir(dir, includeSubdirs, fileInfo)

    #find the lowest processor count (we will need it when computing effectivities)
    minProcCount = -1
    for proc in fileInfo.keys():
        if minProcCount < proc:
            minProcCount = proc

    #find the longest tag name (for formatting the output)
    longestTag = -1
    for timerTag in tagsToProcess:
        length = len(timerTag)
        if longestTag < length:
            longestTag = length

    #generate report
    #create the output file
    outFile = file(os.path.join(dir, "result"), "w")

    #the first line
    outFile.write("n. proc".ljust(longestTag + 2))
    for proc in sorted(fileInfo.keys()):
        outFile.write(str(proc).rjust(7))

    outFile.write("\n")

    if minProcCount > 0:
        for timerTag in tagsToProcess:
            #write info for each tag we wanted to process
            outFile.write(timerTag.ljust(longestTag + 2))
            if fileInfo[minProcCount].has_key(timerTag):
                minProcInfo = fileInfo[minProcCount][timerTag]

                #write times on the first line
                for proc in fileInfo.keys():
                    info = fileInfo[proc][timerTag]
                    outFile.write(str(info.time).rjust(7))

                outFile.write("\n")
                outFile.write("".ljust(longestTag + 2))

                #compute and write effectivities on the second line
                for proc in fileInfo.keys():
                    info = fileInfo[proc][timerTag]
                    effectivity = (minProcInfo.time / (info.time * proc)) * (info.taskSize / minProcInfo.taskSize)
                    outFile.write(str("%.2f" % effectivity).rjust(7))

                outFile.write("\n")

    outFile.close()

def processDir(dir, includeSubdirs, fileInfo):

    patternNumproc = re.compile("No\. of processors: (?P<num>[0-9]+)")
    patternTaskSize = re.compile("Task size: (?P<num>[0-9]+)")


    fileNames = os.listdir(dir)
    path = ""
    for fileName in fileNames:
        filePath = os.path.join(dir, fileName)

        if os.path.isdir(filePath) and includeSubdirs:
            processDir(filePath, includeSubdirs, fileInfo)
        elif os.path.isfile(filePath) and filePath.endswith("out"):
            #read content of the file
            f = file(filePath, "r")
            fileContent = f.read()
            f.close()

            #find the number of processors and task size using the regular expressions
            numprocmatch = patternNumproc.search(fileContent)
            tagsizematch = patternTaskSize.search(fileContent)

            if tagsizematch:
                size = int(tagsizematch.group("num"))
                for timerTag in tagsToProcess:
                    #for each tag we want to process, create the regular expression
                    patternTag = re.compile("\s*"+timerTag+"\s+[0-9]+\s+(?P<num>[0-9]+(\.[0-9]+)?)")
                    tagmatch = patternTag.search(fileContent)

                    if numprocmatch and tagmatch:
                        time = float(tagmatch.group("num"))
                        proc = int(numprocmatch.group("num"))

                        #insert the parsed time, taskSize and no. of processors into the fileInfo dictionar
                        if fileInfo.has_key(proc):
                            if fileInfo[proc].has_key(timerTag):
                                #choose the minimial time
                                if fileInfo[proc][timerTag].time > time:
                                    fileInfo[proc][timerTag].time = time
                            else:
                                fileInfo[proc][timerTag] = TagInfo(timerTag, time, size)
                        else:
                            fileInfo[proc] = {}
                            fileInfo[proc][timerTag] = TagInfo(timerTag, time, size)

if __name__ == '__main__':
    main()