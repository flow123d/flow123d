#!/usr/bin/python

# Script for submitting flow123 tasks to PBS in background.
# It accepts path to a directory with esired task data to be run and path to flow_run.sh script.
# It creates a special directory  for each number of processors the task is about to be run on,
# copies the content of the passed directory to all the newly created directories and then
# submits each job to PBS. Then it periodically executes check if all the jobs are finished,
# after all a report is generated.

import sys
import os
import datetime
import shutil
import time
import commands
import re
import fcntl
import signal

def main():

    #check if the argument were passed correctly
    try:
    	dir = sys.argv[1]
    except IndexError:
        sys.exit()
        pass

    #create base directory
    now = datetime.datetime.now()
    currentDir = os.getcwd()
    newDir = os.path.join(currentDir, now.strftime("%Y_%m_%d-%H.%M.%S"))

    if os.path.exists(newDir) == 0:
        os.mkdir(newDir)

    measurements = 3
    proc = [1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 32]
    iniFiles = {}

    #create directory for each number of processors the task is about to run on
    for nproc in proc:
        measurement = 0

        while measurement < measurements:
            measurement += 1

            targetDir = os.path.join(newDir, str(nproc)+"_"+str(measurement))
            if os.path.exists(targetDir) == 0:
                os.mkdir(targetDir)

            #copy files from the source directory to the newly created one
            #also try to find either flow.ini or trans.ini file

            fileNames = os.listdir(dir)
            for fileName in fileNames:
                filePath = os.path.join(dir, fileName)

                #copy everything except the .pos files
                if os.path.isfile(filePath) and os.path.splitext(fileName)[1] != ".pos":
                    shutil.copy(filePath, targetDir)

                if fileName == "flow.ini" or fileName == "trans.ini":
                    iniFileFullPath = os.path.join(targetDir, fileName)
                    iniFiles[iniFileFullPath] = nproc

    #submit tasks
    for iniFile in iniFiles.iterkeys():
        os.system("./run_flow.sh -np " + str(iniFiles[iniFile]) + " -s " + iniFile + " -m hydra >> benchmark.log")

    completed = 0

    #wait until all the tasks are finished
    while completed == 0:
        time.sleep(60)
        completed = 1

        for nproc in proc:
            measurement = 0
            while measurement < measurements:
                measurement += 1

                targetDir = os.path.join(newDir, str(nproc)+"_"+str(measurement))
                lockFile = os.path.join(targetDir, "lock")
                if os.path.exists(lockFile) == 1:
                    completed =0

    os.system("./make_report.py \"" + str(newDir) + "\"")
    sys.exit()

if __name__ == '__main__':
    main()